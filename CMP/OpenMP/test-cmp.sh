#!/bin/bash
TIME=${TIME:-`which time` -f%e}
DATA=${DATA:-simple-syntetic.su}
ARGS=${ARGS:-"1.98e-7 1.77e-6 101 600 0.002"}
PROG=${PROG:-cmp.opencl.x}
$TIME $PROG $ARGS $DATA
