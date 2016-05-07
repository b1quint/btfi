#!/bin/csh -f 
setenv panport 6342 # service port, must be the same as the one defined on SAM-FP PlugIn GUI 
setenv panhost soarhrc # host where the PlugIn is. It can use machine names or ip addresses
sendsockcmd -h $panhost -p $panport "$argv" -t 40000 # sends the actual command and print the response
