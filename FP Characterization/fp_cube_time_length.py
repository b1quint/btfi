#!/usr/bin/env python3
# -*- coding: utf8 -*-

from datetime import timedelta

# Input ---
exposure_time = 30.0 # [seconds]
number_of_channels = 44 # [--]
readout_time = 2.6 # [seconds]
fp_overhead_time_per_step = 4.15 # [seconds]

total_time = exposure_time * number_of_channels + \
        readout_time * number_of_channels + \
        fp_overhead_time_per_step * (number_of_channels - 1)
total_time = timedelta(seconds=total_time)
print(" Total time: %s" % total_time)
print(" or\n %.2f" % total_time.seconds)
