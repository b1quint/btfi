#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import division, print_function

import argparse
import glob
import os
import pyfits
import time
import watchdog.observers
import watchdog.events

global pattern
global n


class EventHandler(watchdog.events.FileSystemEventHandler):

    def __init__(self, ref_name, **kwargs):
        watchdog.events.FileSystemEventHandler.__init__(self)
        self.ref_name = ref_name
        self.n = 0

    def on_created(self, event):

        print(event.src_path)
        ref_name = self.ref_name
        if os.path.split(event.src_path)[-1] == ref_name:
            return

        path = os.path.split(event.src_path)[0]
        fits_files = glob.glob(os.path.join(path, pattern))
        n = self.n

        if len(fits_files) == 1:
            print(" Creating reference file.")
            time.sleep(0.5)
            fits_file = pyfits.open(event.src_path)[0]
            n = self.n
            fits_file.writeto(ref_name, clobber=True)

        else:
            time.sleep(0.5)
            try:
                fits_file = pyfits.open(event.src_path)[0]
            except IOError:
                print(" Could not read %s" % event.src_path)
                return
            try:
                ref_file = pyfits.open(ref_name)[0]
            except IOError:
                print(" Ref-file does not exist. Creating one")
                ref_file = fits_file
                n = 1
                ref_file.writeto(ref_name, clobber=True)
                return

            n = n + 1
            ref_file.data = ((n - 1) * ref_file.data + fits_file.data) / n
            ref_file.writeto(ref_name, clobber=True)

        self.n = n
        return


if __name__ == '__main__':

    # Parse arguments ---------------------------------------------------------
    parser = argparse.ArgumentParser(description="Try to add fits files as" +
         " they are being acquired in a folder.")
    parser.add_argument('-p', '--path', type=str, default='.',
                        help="Folder to be watched.")
    parser.add_argument('-r', '--ref_name', type=str, default='ref.fits',
                        help="Reference file")
    parser.add_argument('-t', '--template', type=str, default='*.fits')
    args = parser.parse_args()

    # Starting program --------------------------------------------------------
    path = os.path.abspath(args.path)
    pattern = args.template

    print("")
    print(" RealTime FITS Integrator")
    print(" Folder name: %s" % path)
    print(" Running WatchDog. Press Ctrl+C to stop.")
    print("")

    event_handler = EventHandler(args.ref_name)
    observer = watchdog.observers.Observer()
    observer.schedule(event_handler, path)
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()