#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

// Initialization
float **matrix_init(int m, int n) {
	int i, j;
	float **A = (float **) malloc (m*sizeof(float *));
	for (i = 0; i < m; i++)
		 *(A+i) = (float *) malloc (n*sizeof(float));
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			*(*(A+i)+j) = 1.;
	return A;
}
void matrix_free(float **A, int m, int n) {
	int i = 0;
	for (i = 0 ; i < m; i++)
		 free(*(A+i));
	free(A);
}

// Scanning & Printing
void matrix_scan(float **A, int m, int n) {
	int i = 0, j = 0;	
	for (i = 0 ; i < m; i++) {
		for (j = 0; j < n; j++)
			scanf("%f", (*(A+i)+j));
	}
}
void matrix_print(float **A, int m, int n) {
	if (!A) {printf("Empty!\n"); return;}
	int i = 0, j = 0;
	printf("\n");
	for  (i = 0 ; i < m; i++) {
		printf(" |");
		for (j = 0; j < n; j++) 
			printf( "%7lg" , *(*(A+i)+j) );
		printf("%6c|",' ');
		printf("\n"); 
	}	printf("\n");
}
// FILE Scanning & Printing
void matrix_fscan(FILE *fp, float **A, int m, int n) {
	int i = 0, j = 0; 
	for (i = 0 ; i < m; i++) {
		for (j = 0; j < n; j++) 
			fscanf(fp, "%f", (*(A+i)+j));
	}
}
void matrix_fprint(FILE *fp, float **A, int m, int n) {
	if (!A) {return;}
	int i = 0, j = 0;
	fprintf(fp, "\n");
	for  (i = 0 ; i < m; i++) {
		fprintf(fp, "  |");
		for (j = 0; j < n; j++) 
			fprintf(fp, "%8lg" , *(*(A+i)+j) );
		fprintf(fp, "%7c|",' ');
		fprintf(fp, "\n"); 
	}	fprintf(fp, "\n");
}

// Copying
void matrix_copy(float **A, int m, int n, float **B) {
	int i = 0, j = 0;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			B[i][j] = A[i][j];
}

// Swap rows and columns
void matrix_swaprow(float **A, int m, int n, int x, int y) {
	int j = 0;
	float s;
	for (j = 0; j < n; j++) {
		s = A[x][j];
		A[x][j] = A[y][j];
		A[y][j] = s;
	}
}
void matrix_swapcol(float **A, int m, int n, int x, int y) {
	int i = 0;
	float s;
	for (i = 0; i < m; i++) {
		s = A[i][x];
		A[i][x] = A[i][y];
		A[i][y] = s;
	}
}

// Extracting submatrix
void matrix_sub(float **A, int M, int N, int x, int y, float **B) {
	int i = 0, ii = 0, j, jj;
	while (i < M) {
		if (i == x) { i++; continue; }
		j = 0; jj = 0;
		while (j < N) {
			if (j == y) { j++; continue; }
			B[ii][jj] = A[i][j];
			jj++; j++;
		}
		ii++; i++;
	}
}	
	
// Transposition
void matrix_trans(float **A, int m, int n, float **B) {
	int i = 0, j = 0;
	for (i = 0; i < m; i++ )
		for (j = 0; j < n; j++)
			B[j][i] = A[i][j];
}

// Arithmetic
void matrix_mult(float **A, int m, int q, float **B, int p, int n, float **C) {
	if (q != p) {C = NULL; return;}
	int i = 0, j = 0, r = 0; 
	float s = 0.0;
	for (i = 0; i < m; i++) 
		for (j = 0; j < n; j++) {
			s = 0.0;
			for (r = 0; r < q; r++) 
				s += (*(*(A+i)+r))*(*(*(B+r)+j));
			*(*(C+i)+j) = s;
		}
}	
void matrix_plus(float **A, int m, int n, float **B, int p, int q, float **C) {
	m = min(m, p);	n = min(n, q);
	int i = 0, j = 0;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			*(*(C+i)+j) = *(*(A+i)+j) + *(*(B+i)+j);
		}
	}
}

// Inverse (square only)
void matrix_inverse(float **A, int m, int n, float **B) {
	if (!A || !B || (m != n)) {B = NULL; return;}
	float D = matrix_det(A, m, m);
	if (!D) {B = NULL; return;}  
	if (m == 1) {
		if (A[0][0] == 0) {B = NULL; return;}
		B[0][0] = 1/A[0][0]; 
		return;
	}
	int i = 0, j = 0;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			B[j][i] = matrix_cofactor(A, m, n, i, j)/D;
	matrix_realzero(B, n, m);
}
float matrix_cofactor(float **A, int m, int n, int x, int y) {
	return odd(x+y)*matrix_minor(A, m, n, x, y);
}
float matrix_minor(float **A, int m, int n, int x, int y) {
	float **B = matrix_init(m-1, n-1);
	matrix_sub(A, m, n, x, y, B);
	float D = matrix_det(B, m-1, n-1);
	matrix_free(B, m-1, n-1);
	return D;
}

// System
float matrix_det(float **A, int m, int n) {
	if ( m != n) {return (-0);}	
	if (m == 1) {return A[0][0];} 
	int k;	
	float s = 0;
	for (k = 0; k < m; k++) {
		s += A[0][k]*matrix_cofactor(A, m, n, 0, k);
	}
	return s;
	
}
void matrix_kramer(float **A, int m, int n, float **X) {
	if (!A || (m+1 != n)) {return;}
	float D = matrix_det(A, m, m);
	if ((D == 0) || (D == -0)) {X = NULL; return;}
	int i = 0, j = 0;
	float **B = matrix_init(m, n);
	for (j = 0; j < m; j++) {
		matrix_copy(A, m, n, B);		
		for (i = 0; i < m; i++) {
			B[i][j] = A[i][n-1];
		}
		X[j][0] = matrix_det(B, m, m)/D;
	}
	matrix_free(B, m, n);
	matrix_realzero(X, n - 1, 1);
}

// Other
void matrix_realzero(float **A, int m, int n) {
	if (!A) {return;}
	int i = 0, j = 0;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			if (A[i][j] == -0) {A[i][j] = 0.0;} 
}
int odd(int i) {
	switch (i%2) {
	case 0:	 return 1;  break;
	case 1:	 return -1; break;
	default: return 0;  break;
	}
}

