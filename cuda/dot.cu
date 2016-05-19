#include<stdio.h>
#define N (1024*1024)
#define M 1024

__global__ void dot(float *a, float *b, float *c) {
	int i = blockDim.x*blockIdx.x + threadIdx.x, j = threadIdx.x;
	__shared__ float ab[M];
	ab[j] = a[i]*b[i];
	__syncthreads();
	if (!j) {	
		float s = 0.;
		for (i = 0; i < M; i++)
			s += ab[i];
		atomicAdd(c, s);
	}
}

int main(int argc, char *argv[]) {
	int i = 0, size = N*sizeof(float);
	float *a, *b, *c, *dev_a, *dev_b, *dev_c;
	a = (float *) malloc(size);
	b = (float *) malloc(size);
	c = (float *) malloc(sizeof(float));
	cudaMalloc((void **) &dev_a, size);
	cudaMalloc((void **) &dev_b, size);
	cudaMalloc((void **) &dev_c, sizeof(float));

	for (i = 0; i < N; i++) {a[i] = 1.; b[i] = 1.;} *c = 0.;
	cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b, size, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_c, c, sizeof(float), cudaMemcpyHostToDevice);
	dot<<<N/M, M>>>(dev_a, dev_b, dev_c);
	cudaMemcpy(c, dev_c, sizeof(float), cudaMemcpyDeviceToHost);
	printf("%f\n", *c);

	cudaFree(dev_c); cudaFree(dev_b); cudaFree(dev_a);
	free(c); free(b); free(a);
	return 0;
}
