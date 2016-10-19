#!/bin/bash
DATA=${DATA:-simple-syntetic.su}
ARGS=${ARGS:-"1.98e-7 1.77e-6 101 600 0.002"}
PROG=${PROG:-cmp.sycl.x}
$TIME $PROG $ARGS $DATA

