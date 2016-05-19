#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#define abs(x) ( (x > 0) ? (x) : (-x) )
#define M 32*64
#define N 16*64
#define L 128*128
#define l 128

float *matrix_init(int m, int n);
void matrix_print(float *A, int m, int n);
void matrix_mult(float *A, int m, int q, float *B, int p, int n, float *C);

__global__ void dot(float *A, float *B, float *C) {
	int i, j;
	int I = blockDim.x*blockIdx.x + threadIdx.x, J = threadIdx.x, k;
	__shared__ float AB[l];
	float s = 0.;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++) {
			AB[J] = A[i*L+I]*B[I*N+j];
			__syncthreads();
			if (!J) {
				s = 0.;
				for (k = 0; k < l; k++)
					s += AB[k];
				atomicAdd((C+i*N+j), s);
			}
		}
}

int main(int argc, char *argv[]) {
	int i, j, k = 0;	time_t dt, ht;
	float *A, *B, *C, *D, *_A, *_B, *_C;
	A = matrix_init(M, L);
	B = matrix_init(L, N);
	C = matrix_init(M, N);
	D = matrix_init(M, N);
	cudaMalloc((void **) &_A, M*L*sizeof(float));
	cudaMalloc((void **) &_B, L*N*sizeof(float));
	cudaMalloc((void **) &_C, M*N*sizeof(float));
	srand(time(NULL));
	for (i = 0; i < M; i++)
		for (j = 0; j < L; j++) 
			*(A+i*L+j) = 1.;
	for (i = 0; i < L; i++)
		for (j = 0; j < N; j++) 
			*(B+i*N+j) = 1.;

	dt = time(NULL);
	cudaMemcpy(_A, A, M*L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(_B, B, L*N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(_C, C, M*N*sizeof(float), cudaMemcpyHostToDevice);
	dot<<<L/l, l>>>(_A, _B, _C);
	cudaMemcpy(C, _C, M*N*sizeof(float), cudaMemcpyDeviceToHost);
	dt = time(NULL) - dt;
	printf("device: %d sec\n", (int) dt);
	fflush(stdout);
	
	ht = time(NULL); 
	matrix_mult(A, M, L, B, L, N, D);
	ht = time(NULL) - ht; 
	printf("host: %d sec\n", (int) ht );		
	printf("acceleration: %.0lf\n", (((double)ht)/((double)dt)) );		

	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++) 
			if(abs(*(C+i*N+j) - *(D+i*N+j)) > .0000001) 
				k++;
	printf("error: %d\n", k);

	cudaFree(_C); cudaFree(_B); cudaFree(_A);
	free(C); free(B); free(A);
	return 0;
}

// Matrix functions
float *matrix_init(int m, int n) {
	int i, j;
	float *A = (float *) malloc (m*n*sizeof(float *));
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			*(A+i*n+j) = 0.;
	return A;
}
void matrix_print(float *A, int m, int n) {
	if (!A) {printf("Empty!\n"); return;}
	int i = 0, j = 0;
	printf("\n");
	for  (i = 0 ; i < m; i++) {
		printf(" |");
		for (j = 0; j < n; j++) 
			printf( "%7lg" , *(A+i*n+j) );
		printf("%6c|",' ');
		printf("\n"); 
	}	printf("\n");
}
void matrix_mult(float *A, int m, int q, float *B, int p, int n, float *C) {
	if (q != p) {C = NULL; return;}
	int i = 0, j = 0, r = 0; 
	float s = 0.0;
	for (i = 0; i < m; i++) 
		for (j = 0; j < n; j++) {
			s = 0.0;
			for (r = 0; r < q; r++) 
				s += (*(A+i*q+r))*(*(B+r*n+j));
			*(C+i*n+j) = s;
		}
}	
