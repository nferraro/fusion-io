#!/bin/bash

time $FIO_ROOT/trace/_$FIO_ARCH/trace -m3dc1 data/m3dc1/C1.h5 -1 -dR 0.1 -p 8 -t 1000 -s 300
