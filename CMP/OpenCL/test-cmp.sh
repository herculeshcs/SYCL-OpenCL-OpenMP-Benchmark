#!/bin/bash
DATA=${DATA:-simple-synthetic.su}
DEVICE=2
./cmp.opencl.x 1.98e-7 1.77e-6 101 600 0.002 $DATA 5 64 $DEVICE

