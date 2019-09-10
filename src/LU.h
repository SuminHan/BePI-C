/* This code is modified version of: https://rosettacode.org/wiki/LU_decomposition */
typedef double **mat;
 
void mat_zero(mat, int);
mat mat_new(int);
mat mat_copy(void*, int);
void mat_del(mat);
void mat_show(mat, char *, int);
mat mat_mul(mat, mat, int);
void mat_pivot(mat, mat, int);
void mat_LU(mat, mat, mat, mat, int);
mat mat_invL(mat, int);
mat mat_invU(mat, int);
mat mat_invP(mat, int);
