#!/bin/bash

mkdir -p snapshot

for angle in {0..359..10}
do
    sed "/^start_angle/s/60/${angle}/" logo.py >| tmp.py
    python tmp.py
    cp logo.png snapshot/t_${angle}.png
done

convert -loop 0 -delay 10 snapshot/t_{0..359..10}.png  namd.gif
