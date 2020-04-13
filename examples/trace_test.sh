#!/bin/bash

if [[ "$1" == "m3dc1" ]]; then
    FILE="-m3dc1 data/m3dc1/C1.h5 -1 -m3dc1 data/m3dc1/C1.h5 1 5."
elif [[ "$1" == "geqdsk" ]]; then
    FILE="-geqdsk data/geqdsk/g158115.04701"
fi

if [ -z "$FILE" ]; then
    echo "Please choose one of:"
    echo " m3dc1"
    echo " geqdsk"    
else
    time $FIO_ROOT/trace/_$FIO_ARCH/trace $FILE -dR 0.05 -p 12 -t 1000 -s 300
fi