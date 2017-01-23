#!/bin/sh

gcc tbt.cc -o tbt \
    -I$TRACY_LIB/tracy/inc \
    -L$TRACY_LIB/tracy/lib -ltracy \
    -I$NUM_REC/inc -L$NUM_REC/lib -lnum_rec \
    -lstdc++ \
    -lm
