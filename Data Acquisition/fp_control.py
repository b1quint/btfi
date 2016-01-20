#!/usr/bin/python
# -*- coding: utf8 -*-
import time
import sys
import socket
import argparse
import logging

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)6s ' + \
                              '- %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

handler = logging.StreamHandler()
handler.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(handler)

if __name__ == '__main__':

    ## Parse arguments from command line ---
    parser = argparse.ArgumentParser(description="Control remotely the FP.")
    parser.add_argument('-a','--address', type=str, default='139.229.18.226', help="FP server IP (String, Default: 139.229.18.226)")
    parser.add_argument('-p','--port', type=int, default=6342, help="FP server PORT (Integer, Default: 6342)")
    parser.add_argument('-b','--buffer', type=int, default=1024, help="Buffer size (Integer, Default: 1024)")
    parser.add_argument('-c','--counter', type=int, default=0, help="Counter number (Integer, Default: 0)")
    parser.add_argument('start', type=int, help="Start Z value (Integer)")
    parser.add_argument('step', type=int, help="Step Z value (Integer)")
    parser.add_argument('stop', type=int, help="Stop Z value (Integer)")
    args = parser.parse_args()

    ## Setting up the communication ---
    TCP_IP = args.address
    TCP_PORT = args.port
    BUFFER_SIZE = args.buffer

    logger.debug("Socket AF_INET {}".format(socket.AF_INET))
    logger.debug("Socket SOCK_STREAM {}".format(socket.SOCK_STREAM))
    # s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    logger.info("Connecting to FP server.")
    try:
        # s.connect((TCP_IP, TCP_PORT))
        s = socket.create_connection((TCP_IP, TCP_PORT), timeout=10)
        logger.info("Connection stablished successfully!")
        s.close()
    except socket.error as v:
        logger.error("Connection failed - error code: {0}".format(v))
        sys.exit(1)

    ## Setting up the scan ---
    counter = args.counter
    Z_Start = args.start
    Z_Step = args.step
    Z_Stop = args.stop

    ## Start scan ---
    logger.info("Start scan:")
    for z in range(Z_Start, Z_Stop, Z_Step):
        logger.info("Counter = {}, Z = {} bcv".format(counter, z))

        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((TCP_IP, TCP_PORT))
        s.send('moveabs %d' % z)
        data = s.recv(BUFFER_SIZE)
        logger.debug(data)
        s.close()

        if data[0] != 'D':
            sys.exit()

        # try:
        #     os.system("gphoto2 --capture-image ")
        #     os.system("gphoto2 --get-file=%d"%counter)
        #     counter = counter + 1
        # except:
        #     pass

        time.sleep(2)
        counter += 1