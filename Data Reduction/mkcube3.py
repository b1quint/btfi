#!/usr/bin/python
# -*- coding: utf8 -*-
"""
    mkCube 3.0.0

    This file is part of the Butterfly Software System. It is used to organize
    images obtained through a BTFI scan and use them to create a data-cube.

    by Bruno Quint <bquint at astro.iag.usp.br>
"""

from __future__ import division, print_function

try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
    from sys import exit

    print("\n\tmkcube3 --")
    print("\tAstroPy module was not found in the system.")
    print("\tIt is required to run this program.")
    print("\tPlease, take a look at http://www.astropy.org/ for information")
    print("\tabout installing it on your system.")
    print("\n\tLeaving now.\n")
    exit()


def error(message):
    """
    Prints an error message in the terminal and leaves the program.
    :param message: The message that will be given for the user.
    :return: None
    """
    from sys import exit

    message = '\033[1;38 [!] %s \033[1;m' % message
    print(message)
    exit()


def warning(message):
    """
    Print a colour warning message in the screen.
    """
    message = ' \033[1;33m%s\033[1;m' % message
    print(message)


class MakeCube(object):
    def __init__(self):
        """
        Class constructor. It is used to initialize all the variables used in
        the program silently.
        """
        self.depth = None
        self.height = None
        self.input_files = None
        self.input_table = None
        self.mode = None
        self.output_file = None
        self.verbose = False
        self.width = None
        self.header = None
        self.cube = None

    @property
    def create_input_table(self):
        """
        :return: None
        """
        from numpy import array
        from sys import stdout

        table = []
        header_field = None
        if self.verbose:
            print(" Reading files' headers to organize the data.")
        if self.mode.lower() in ['fp', 'fabry-perot']:
            header_field = 'FPZCURR'
        elif self.mode.lower() in ['ibtf']:
            header_field = 'IBCURANG'
        else:
            error(" I am sorry but I could not find the BTFI mode" +
                  " in the header.")

        i, number_of_files = 0.0, len(self.input_files)
        header_value = None
        for input_file in self.input_files:
            if self.verbose:
                i += 1.0
                stdout.write("\r  %.2f%%" % (i * 100. / number_of_files))
                stdout.flush()
            try:
                header_value = pyfits.getheader(input_file)[header_field]
            except (IOError, KeyError):
                pass

            table.append([input_file, float(header_value)])

        print("")
        table = array(table, dtype=str)
        files = table[:, 0]
        poses = array(table[:, 1], dtype=float)
        indexes = poses.argsort()
        files = files[indexes]
        poses = poses[indexes]
        table = array([files, poses], dtype=str).transpose()

        return table


    def get_cube_dimensions(self):
        """
        Determines what are the data-cube's dimensions.
        :return: None
        """
        self.width = self.get_width
        self.height = self.get_height
        self.depth = self.get_depth
        if self.verbose:
            print(" The data-cube will have the following dimensions:")
            print("  [%d, %d, %d]" % (self.width, self.height, self.depth))


    @property
    def get_depth(self):
        """
        :return: Returns data-cube depth.
        """
        from numpy import array, float32, unique
        pos = self.input_table[:, 1]
        pos = array(unique(pos), dtype=float32)
        return pos.size


    @property
    def get_height(self):
        """
        :return: the data-cube width based on the FITS header.
        """
        height = int(pyfits.getheader(self.input_files[0])['NAXIS2'])
        return height


    @property
    def get_instrument_mode(self):
        """
        Read the instrument mode from the header of a random input file.
        """
        filename = self.input_files[0]
        if self.verbose:
            print(" Reading instrument mode from file:\n  %s" % filename)
        try:
            instrument_mode = pyfits.getheader(filename)['INSTRMOD']
        except KeyError:
            warning("Could not find instrument mode!")
            instrument_mode = 'unknown'
        if self.verbose:
            print(" Instrument mode: %s" % instrument_mode)
        return instrument_mode


    def get_output_filename(self, overwrite=None, verbose=False):
        """
        This is a generic method used to check if a file called 'name' already
        exists. If so, it starts some interaction with the user.

        @keyword overwrite: if False, this method will interact with the user to
        ask if 'name' file shall be overwritten or if a new name will be given. If
        True, 'name' file is automatically overwritten.

        @keyword verbose: force verbose mode on even when overwrite is automatic.

        v1.0.1 - added 'overwrite' keyword.
               - added 'verbose' keyword.
        """
        import os
        import sys

        name = self.output_file
        v = False if (overwrite is True) else True
        if self.verbose:
            print("\n Writing to output file %s" % name)

        while os.path.exists(name):

            if overwrite in ['y', 'Y', True]:
                if v or verbose:
                    print(" Deleting '%s' file now." % name)
                os.remove(name)

            elif overwrite in ['', 'n', 'N', False]:
                name = raw_input("   Please, enter a new filename:\n   > ")

            elif overwrite in ['q']:
                if v:
                    print(" Exiting program.")
                sys.exit()

            else:
                overwrite = \
                    raw_input("  '%s' file exist. Overwrite? (y/[n])" % name)
                if v:
                    print(" Writing data-cube to %s" % name)

        return name


    @property
    def get_width(self):
        """
        :return: the data-cube width based on the FITS header.
        """
        width = int(pyfits.getheader(self.input_files[0])['NAXIS1'])
        return width


    @property
    def make_cube(self):
        from numpy import array, empty, unique
        from sys import stdout

        cube = None
        table = self.input_table
        files = table[:, 0]
        poses = array(table[:, 1], dtype=float)
        unique_poses = unique(poses)

        if self.verbose:
            print(" Creating data-cube...")

        if self.verbose:
            print(" Allocating memory...")
        try:
            cube = empty((self.depth, self.height, self.width), dtype=float)
        except MemoryError:
            error(" Ops! The cube was too big for your computer.\n"
                  " Try making a smaller one.\n"
                  " Leaving now.")
        if self.verbose:
            print(" Ok.")
            print(" Filling data-cube:")

        i, number_of_positions = 0, unique_poses.size
        for pos in unique_poses:
            temp_files = files[poses == pos]
            frame = self.make_frame(temp_files)
            cube[i] = frame
            if self.verbose:
                i += 1
                stdout.write("\r  %3.2f%%" % (i * 100 / number_of_positions))
                stdout.flush()
        print("\n Done.")
        return cube


    def make_frame(self, files):

        from numpy import empty

        dummy_cube = empty((len(files), self.height, self.width))
        for i in range(len(files)):
            try:
                dummy_cube[i] = pyfits.getdata(files[i])
            except IOError:
                warning(" File %s may be corrupted." % files[i])

        return dummy_cube.mean(axis=0)


    @property
    def make_header(self):
        """
        Creates the cube's header.
        :return: header
        """
        header = pyfits.getheader(self.input_files[0])
        header[''] = "--- mkCube3 Calibration ---"

        if self.mode.lower() in ['fp', 'fabry-perot']:
            header['DISPAXIS'] = 3
            header['CRPIX3'] = 1
            header['CRVAL3'] = float(header['FPZINIT'])
            header['CDELT3'] = float(header['FPZDELT'])
            header['CTYPE3'] = 'LINEAR'
            header['CUNIT3'] = 'BNV'
            header['C3_3'] = float(header['FPZDELT'])

        elif header['INSTRMOD'].lower() in ['ibtf']:
            header['DISPAXIS'] = 3
            header['CRPIX3'] = 1
            header['CRVAL3'] = float(header['TFAINIT'])
            header['CDELT3'] = float(header['TFADELT'])
            header['CTYPE3'] = 'LINEAR'
            header['CUNIT3'] = 'degrees'
            header['C3_3'] = float(header['TFADELT'])

        else:
            warning("Invalid BTFI Instrument Mode.")
            warning("Dummy calibration will be added to the data-cube.")
            header['DISPAXIS'] = 3
            header['CRPIX3'] = 1
            header['CRVAL3'] = 1
            header['CDELT3'] = 1
            header['CTYPE3'] = 'LINEAR'
            header['CUNIT3'] = 'channel'
            header['C3_3'] = 1

        return header


    def parse_arguments(self):
        """
        Parses the arguments given by the user and stores the information in
        the class properties.
        :rtype : none
        """
        import argparse
        import glob

        parser = argparse.ArgumentParser(
            description="Build a data-cube from image files.")

        parser.add_argument('-o', '--output', metavar='output', type=str,
                            default="cube.fits", help="Name of the output cube.")

        parser.add_argument('-q', '--quiet', action='store_true',
                            help="Run quietly.")

        parser.add_argument('files', metavar='files', type=str, nargs='+',
                            help="input filenames.")

        args = parser.parse_args()

        self.input_files = args.files
        if len(self.input_files) == 1:
            self.input_files = glob.glob(self.input_files[0])
        self.output_file = args.output
        self.verbose = not args.quiet

    def print_header(self):
        """
        Print this program's header.
        :rtype : None
        """
        if not self.verbose:
            pass
        else:
            print(__doc__)
        return None

    def run(self):
        """
        The main thread.
        """
        self.parse_arguments()
        self.print_header()
        self.mode = self.get_instrument_mode
        self.output_file = self.get_output_filename()
        self.input_table = self.create_input_table
        self.get_cube_dimensions()
        self.header = self.make_header
        self.cube = self.make_cube
        self.save_datacube()
        if self.verbose:
            print(" All done!\n")

    def save_datacube(self):
        """
        Just save the data-cube in the disk.
        """
        if self.verbose:
            print(" Saving data-cube on the disk...")
        pyfits.writeto(self.output_file, self.cube, self.header, output_verify="fix")
        if self.verbose:
            print(" Done!")
        del self.cube
        return


if __name__ == '__main__':
    main = MakeCube()
    main.run()