#!/bin/sh
root -L -b << EOF
    .L TwoPCorr.C+
    .x Run_TwoPCorr.C()
EOF
