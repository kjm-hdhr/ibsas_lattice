#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

void print_mat33(int32_t mat33[9]);
void print_mat44(int32_t mat44[16]);
uint64_t invmod (uint64_t a, uint64_t p);
// mat33[9]={11,12,13,21,22,23,31,32,33}
// (r-1)*3+(c-1)
// r1c1=0 r1c2=1 r1c3=2
// r2c1=3 r2c2=4 r2c3=5
// r3c1=6 r3c2=7 r3c3=8
//
// mat44[16]={11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44}
// (r-1)*4+(c-1)
// r1c1= 0 r1c2= 1 r1c3= 2 r1c4= 3
// r2c1= 4 r2c2= 5 r2c3= 6 r2c4= 7
// r3c1= 8 r3c2= 9 r3c3=10 r3c4=11
// r4c1=12 r4c2=13 r4c3=14 r4c4=15

void put_value33(int32_t mat33[9], int r, int c, int32_t v);
int32_t get_value33(int32_t mat33[9], int r, int c);
void put_value44(int32_t mat44[16], int r, int c, int32_t v);
int32_t get_value44(int32_t mat44[16], int r, int c);

int32_t mat33_det_elem(int32_t a, int32_t b, int32_t c, int32_t p);
int32_t mat33_det(int32_t mat33[9], int32_t p);

void make_cofactor_matrix(
	int32_t mat44[16], int32_t mat33[9], int r, int c);

int32_t make_adjugate_matrix_elem(
	int r, int c, int32_t mat44[16], int p);

void make_adjugate_mat44(int32_t a[16], int32_t mat44[16],int p);
int32_t mat44_det(int32_t mat44[16], int32_t p);

int32_t multiply(int32_t a, int32_t b, int p);
void multipy_elem(int row, int col, int32_t a[16], int32_t b[16], int32_t c[16], int p);
void multiply_mat(int32_t a[16], int32_t b[16], int32_t c[16], int p);

#endif