#!/bin/bash
gcc evec_evol_gw.c -lm -lgsl -lgslcblas
./a.out IC.txt
make && ./rebound
python3 plotcmp.py
