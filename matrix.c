#include "matrix.h"
void print_mat33(int32_t mat33[9]){
	int i;
	printf("mat33=[\n");
	for(i=0;i<3;i++){
		printf("%d,",mat33[i*3]);
		printf("%d,",mat33[i*3+1]);
		printf("%d\n",mat33[i*3+2]);
	}
	printf("]\n");
}
void print_mat44(int32_t mat44[16]){
	int i;
	printf("mat44=[\n");
	for(i=0;i<4;i++){
		printf("%d,",mat44[i*4]);
		printf("%d,",mat44[i*4+1]);
		printf("%d,",mat44[i*4+2]);
		printf("%d\n",mat44[i*4+3]);
	}
	printf("]\n");
}

uint64_t invmod (uint64_t a, uint64_t p)
{
  uint64_t j = 1, i = 0, b = p, c = a, x, y;
  while (c != 0)
  {
    x = b / c;
    y = b - x * c;
    b = c; 
    c = y;
    y = j;
    j = i - j * x;
    i = y;
  }
  if ((int64_t)i < 0)
    i += p;
  return i;
}
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

void put_value33(int32_t mat33[9], int r, int c, int32_t v){
	mat33[(r-1)*3+(c-1)]=v;
}
int32_t get_value33(int32_t mat33[9], int r, int c){
	return mat33[(r-1)*3+(c-1)];
}
void put_value44(int32_t mat44[16], int r, int c, int32_t v){
	mat44[(r-1)*4+(c-1)]=v;
}
int32_t get_value44(int32_t mat44[16], int r, int c){
	return mat44[(r-1)*4+(c-1)];
}

int32_t mat33_det_elem(
	int32_t a, int32_t b, int32_t c, int32_t p
){
	//printf("%d x %d x %d = ",a,b,c);
	int64_t ret=((int64_t)a*b)%p;
	ret=(ret*c)%p;
	//printf("%d\n",ret);
	return ret;
}
int32_t mat33_det(int32_t mat33[9], int32_t p){
	int64_t det=0;
	// r1c1*r2c2*r3c3
	det=mat33_det_elem(mat33[0],mat33[4],mat33[8],p);
	// r1c2*r2c3*r3c1
	det=(det+mat33_det_elem(mat33[1],mat33[5],mat33[6],p))%p;
	// r1c3*r2c1*r3c2
	det=(det+mat33_det_elem(mat33[2],mat33[3],mat33[7],p))%p;
	// r1c3*r2c2*r3c1
	det=(det-mat33_det_elem(mat33[2],mat33[4],mat33[6],p))%p;
	// r1c2*r2c1*r3c3
	det=(det-mat33_det_elem(mat33[1],mat33[3],mat33[8],p))%p;
	// r1c1*r2c3*r3c2
	det=(det-mat33_det_elem(mat33[0],mat33[5],mat33[7],p))%p;
	return det%p;
}

void make_cofactor_matrix(int32_t mat44[16], int32_t mat33[9], int r, int c){
	// (r-1)*4 ~ (r-1)*4+4
	// 
	int i,j;
	int index;
	int index33=0;
	for(i=1;i<=4;i++){//row
		if(i==r){
			continue;
		}
		for(j=1;j<=4;j++){//col
			if(j==c){
				continue;
			}
			index=(i-1)*4+(j-1);
			mat33[index33]=mat44[index];
			index33++;
		}
	}
}
int32_t make_adjugate_matrix_elem(int r, int c, int32_t mat44[16], int p){
	int32_t mat33[9];
	make_cofactor_matrix(mat44,mat33,r,c);
	//printf("r=%d, c=%d\n",r,c);
	//print_mat33(mat33);
	return mat33_det(mat33,p);
}
void make_adjugate_mat44(int32_t a[16], int32_t mat44[16],int p){
	int32_t mat44_element;
	int i,j;
	for(i=1;i<=4;i++){
		for(j=1;j<=4;j++){
			mat44_element=make_adjugate_matrix_elem(i,j,mat44,p);
			mat44_element=((i+j)%2==0)?mat44_element:mat44_element*(-1);
			put_value44(a,j,i,mat44_element);
		}
	}
}
int32_t mat44_det(int32_t mat44[16], int32_t p){
	int32_t mat33[9];
	int64_t a;
	int64_t ret;
	make_cofactor_matrix(mat44,mat33,1,1);
	a=get_value44(mat44,1,1);
	a=(a*mat33_det(mat33,p))%p;
	ret=a;

	make_cofactor_matrix(mat44,mat33,1,2);
	a=get_value44(mat44,1,2);
	a=(a*mat33_det(mat33,p))%p;
	ret=(ret-a)%p;

	make_cofactor_matrix(mat44,mat33,1,3);
	a=get_value44(mat44,1,3);
	a=(a*mat33_det(mat33,p))%p;
	ret=(ret+a)%p;

	make_cofactor_matrix(mat44,mat33,1,4);
	a=get_value44(mat44,1,4);
	a=(a*mat33_det(mat33,p))%p;
	ret=(ret-a)%p;

	return ret;
}

int32_t multiply(int32_t a, int32_t b, int p){
	int64_t ret=a;
	ret=(ret*b)%p;
	return ret;
}
void multipy_elem(int row, int col, int32_t a[16], int32_t b[16], int32_t c[16], int p){
	//C=A*B
	c[(row-1)*4+(col-1)]=
		(multiply(a[(row-1)*4],b[(col-1)],p)+
		multiply(a[(row-1)*4+1],b[(col-1)+4],p)+
		multiply(a[(row-1)*4+2],b[(col-1)+8],p)+
		multiply(a[(row-1)*4+3],b[(col-1)+12],p))%p;
}
void multiply_mat(int32_t a[16], int32_t b[16], int32_t c[16], int p){
//C=A*B
// mat44[16]={11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44}
// (r-1)*4+(c-1)
// r1c1= 0 r1c2= 1 r1c3= 2 r1c4= 3
// r2c1= 4 r2c2= 5 r2c3= 6 r2c4= 7
// r3c1= 8 r3c2= 9 r3c3=10 r3c4=11
// r4c1=12 r4c2=13 r4c3=14 r4c4=15
	int i,j;
	for(i=1;i<=4;i++){
		for(j=1;j<=4;j++){
			multipy_elem(i,j,a,b,c,p);
		}
	}
}

