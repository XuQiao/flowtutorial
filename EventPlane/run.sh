#!/bin/sh
root -L -b << EOF
    .L EventPlaneAna.C+
    .x Run_EventPlaneAna.C()
EOF
