# ganEvection
secular code for evection resonance

## Prerequisites:
The code needs GSL, python3, numpy and matplotlib

## Usage
The initial conditions need to be given in IC.txt

To run only secular code, run:
```
gcc evec_evol_gw.c -lm -lgsl -lgslcblas && ./a.out IC.txt
```
To run N-body code, copy the two folders to examples folder in rebound installation i.e at "rebound/examples". Then run:
```
bash run.sh
```
