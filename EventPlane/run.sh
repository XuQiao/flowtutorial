#!/bin/sh
root -L -b << EOF
    .L EventPlaneAna3sub.C+
    .x Run_EventPlaneAna3sub.C()
EOF
