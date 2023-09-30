# ibsas_lattice

We implemented ID base lattice aggregate signature by using Crystal-dilithium.

matrix.h and matrix.c
matrix operation functions

sign_agg.h and sign_agg.c include functions
key setup, key derivation, aggregate signature generation, aggregate signture verificaiton

packing_agg.h and packing_agg.c include functions
pack and unpack keys and a signature into/from binary array 

compile
1: prepare source codes from Crystal-dilithium.
These are 
sign.c rounding.c reduce.c randombytes.c poly.c polyvec.c packing.c fips202.c ntt.c symmetric-shake.c and header files 

2: put these files in the same directory

3: gcc -o test_ibsas_lattice sign.c rounding.c reduce.c randombytes.c poly.c polyvec.c packing.c fips202.c ntt.c symmetric-shake.c sign_agg.c matrix.c packing_agg.c test_ibsas_lattice.c

4: run the generated binary file.

