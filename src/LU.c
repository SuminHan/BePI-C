/* This code is modified version of: https://rosettacode.org/wiki/LU_decomposition */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#define foreach(a, b, c) for (int a = b; a < c; a++)
#define for_i foreach(i, 0, n)
#define for_j foreach(j, 0, n)
#define for_k foreach(k, 0, n)
#define for_ij for_i for_j
#define for_ijk for_ij for_k
#define _dim int n
#define _swap(x, y) { typeof(x) tmp = x; x = y; y = tmp; }
#define _sum_k(a, b, c, s) { s = 0; foreach(k, a, b) s+= c; }
 
typedef double **mat;
 
#define _zero(a) mat_zero(a, n)
void mat_zero(mat x, int n) { for_ij x[i][j] = 0; }
 
#define _new(a) a = mat_new(n)
mat mat_new(_dim)
{
	mat x = (double**) malloc(sizeof(double*) * n);
	x[0]  = (double*) malloc(sizeof(double) * n * n);
 
	for_i x[i] = x[0] + n * i;
	_zero(x);
 
	return x;
}
 
#define _copy(a) mat_copy(a, n)
mat mat_copy(void *s, _dim)
{
	mat x = mat_new(n);
	for_ij x[i][j] = ((double (*)[n])s)[i][j];
	return x;
}
 
#define _del(x) mat_del(x)
void mat_del(mat x) { free(x[0]); free(x); }
 
#define _QUOT(x) #x
#define QUOTE(x) _QUOT(x)
#define _show(a) printf(QUOTE(a)" =");mat_show(a, 0, n)
void mat_show(mat x, char *fmt, _dim)
{
	if (!fmt) fmt = "%8.4g";
	for_i {
		printf(i ? "      " : " [ ");
		for_j {
			printf(fmt, x[i][j]);
			printf(j < n - 1 ? "  " : i == n - 1 ? " ]\n" : "\n");
		}
	}
}
 
#define _mul(a, b) mat_mul(a, b, n)
mat mat_mul(mat a, mat b, _dim)
{
	mat c = _new(c);
	for_ijk c[i][j] += a[i][k] * b[k][j];
	return c;
}

#define _pivot(a, b) mat_pivot(a, b, n)
void mat_pivot(mat a, mat p, _dim)
{
	for_ij { p[i][j] = (i == j); }
	for_i  {
		int max_j = i;
		foreach(j, i, n)
			if (fabs(a[j][i]) > fabs(a[max_j][i])) max_j = j;
 
		if (max_j != i)
			for_k { _swap(p[i][k], p[max_j][k]); }
	}
}
 
#define _LU(a, l, u, p) mat_LU(a, l, u, p, n)
void mat_LU(mat A, mat L, mat U, mat P, _dim)
{
	_zero(L); _zero(U);
	_pivot(A, P);
 
	mat Aprime = _mul(P, A);
 
	for_i  { L[i][i] = 1; }
	for_ij {
		double s;
		if (j <= i) {
			_sum_k(0, j, L[j][k] * U[k][i], s)
			U[j][i] = Aprime[j][i] - s;
		}
		if (j >= i) {
			_sum_k(0, i, L[j][k] * U[k][i], s);
			L[j][i] = (Aprime[j][i] - s) / U[i][i];
		}
	}
 
	_del(Aprime);
}
 
/* Reference: http://www.mymathlib.com/c_source/matrices/linearsystems/upper_triangular.c */
#define _invL(l) mat_invL(l, n)
mat mat_invL(mat L, _dim)
{
	mat X = _new(X);
	/*
	1. for k = 1 to n
	2.   X[k,k] = l/L[k,k]
	3.   for i = k+1 to n
	4.     X[i,k] = -L[i, k:i-1]*X[k:i-1,k]/L[i,i]
	5.   end for i
	6. end for k
	*/
	for_k {
		X[k][k] = 1/L[k][k];
		for (int i = k+1; i < n; i++){
			X[i][k] = 0;
			for (int j = k; j < i; j++){
				X[i][k] += -L[i][j]*X[j][k];
			}
			X[i][k] /= L[i][i];
		}
	}
	return X;
}

#define _invU(u) mat_invU(u, n)
mat mat_invU(mat U, _dim)
{
	mat X = _new(X);

	int i, j, k;
	X[n-1][n-1] = 1/U[n-1][n-1];
	for (i = n-2; i >= 0; i--){
		X[i][i] = 1/U[i][i];
		for (j = n-1; j > i; j--){
			X[i][j] = 0.0;
			for (k = i+1; k <= j; k++){
				X[i][j] += -U[i][k] * X[k][j];
			}
			X[i][j] /= U[i][i];
		}
	}

	return X;
}

int Upper_Triangular_Inverse(double *U, int n)
{
   int i, j, k;
   double *p_i, *p_j, *p_k;
   double sum;

//         Invert the diagonal elements of the upper triangular matrix U.

   for (k = 0, p_k = U; k < n; p_k += (n + 1), k++) {
      if (*p_k == 0.0) return -1;
      else *p_k = 1.0 / *p_k;
   }

//         Invert the remaining upper triangular matrix U.

   for (i = n - 2, p_i = U + n * (n - 2); i >=0; p_i -= n, i-- ) {
      for (j = n - 1; j > i; j--) {
         sum = 0.0;
         for (k = i + 1, p_k = p_i + n; k <= j; p_k += n, k++ ) {
            sum += *(p_i + k) * *(p_k + j);
         }
         *(p_i + j) = - *(p_i + i) * sum;
      }
   }
  
   return 0;
}

#define _invP(p) mat_invP(p, n)
mat mat_invP(mat P, _dim)
{
	mat P_inv = _new(P_inv);

	for_ij {
		P_inv[i][j] = P[j][i];
	}

	return P_inv;
}
 
/*
double A4[][4] = {{1, 2, 3, 4}, {1, 2 ,5, 3}, {3, 5 ,2, 6}, {3, 4, 2 ,1}};
int main()
{
	int n = 4;
	mat A, L, P, U, L_inv, U_inv, M, PA, T, LU, LLi, UUi, P_inv;
 
	n = 4;
 
	_new(L); _new(P); _new(U);
	A = _copy(A4);
	_LU(A, L, U, P);
	P_inv = _invP(P);
	_show(A); _show(L); _show(U); _show(P); _show(P_inv);

	L_inv = _invL(L);
	U_inv = _invU(U);
	_show(L_inv);
	_show(U_inv);

	PA = _mul(P, A);
	_show(_mul(P_inv, PA));
	LU = _mul(L, U);
	
	LLi = _mul(L, L_inv);
	UUi = _mul(U, U_inv);
	_show(LU);
	_show(LLi);
	_show(UUi);

	_show(_mul(U_inv, _mul(L_inv, PA)));

	return 0;
}
*/
