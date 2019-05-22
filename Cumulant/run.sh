#!/bin/sh
root -L -b << EOF
    .L QCumulant.C+
    .x Run_QCumulant.C()
EOF
