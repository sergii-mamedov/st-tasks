#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

const double pi = 3.14159265358979323846;

void X_init(double *a, int numTicks, double delta) {
	for (int i = 0; i <= numTicks; ++i)
	for (int j = 0; j <= numTicks; ++j)
		a[i*(numTicks+1) + j] = (double)j * delta;
}

void Y_init(double *a, int numTicks, double delta) {
	for (int i = 0; i <= numTicks; ++i)
	for (int j = 0; j <= numTicks; ++j)
		a[i*(numTicks+1) + j] = (double)i * delta;
}

void array_show(double *a, int numTicks) {
	for (int i = 0; i <= numTicks; ++i) {
	for (int j = 0; j <= numTicks; ++j)
		printf("%3d ", (int)a[i*(numTicks+1) + j]);

	printf("\n");
	}
}

__global__ void kernel(double *x, double *y, double *z, double mux, double muy, double b4, int size)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= size)
    	return;

    // initial values
    double x0 = (double)x[tid];
    double y0 = (double)y[tid];
    double px0 = 0.0, py0 = 0.0;

    // проверка вхождение в круг, до первой итерации
    if (x0*x0 + y0*y0 >= 1.0) {
    	z[tid] = 0.0;
    	return;
    }

    int n;
    double x1, px1, px2, y1, py1, py2;
    for (n = 0; n <= 100001; n++) {
        if (x0*x0 + y0*y0 < 1.0) {
            x1  =  x0 * cos(2.0*pi*mux) + px0 * sin(2.0*pi*mux);
            px1 = -x0 * sin(2.0*pi*mux) + px0 * cos(2.0*pi*mux);
            y1  =  y0 * cos(2.0*pi*muy) + py0 * sin(2.0*pi*muy);
            py1 = -y0 * sin(2.0*pi*muy) + py0 * cos(2.0*pi*muy);
            px2 = px1 + b4 * (x1*x1*x1 - 3.0*x1*y1*y1);
            py2 = py1 - b4 * (y1*y1*y1 - 3.0*x1*x1*y1);
            x0  =  x1;  y0 =  y1;
            px0 = px2; py0 = py2;
        } else {
        	break;
        }
    }

    n--;
    z[tid] = (double)(n);
}

int main(void)
{
	const double mux = 0.32;
	const double muy = 0.32;
	const double b4  = 0.50;
	int numTicks = 10;
	double delta = 1.0 / numTicks;
	int arraySize =(numTicks + 1) * (numTicks + 1);
	int numBytes = arraySize * sizeof(double);
	double *x, *y, *z, *x_dev, *y_dev, *z_dev;

	// allocate host memory
	x = (double *) malloc(numBytes);
	y = (double *) malloc(numBytes);
	z = (double *) malloc(numBytes);

	// allocate X, Y array
	X_init(x, numTicks, delta);
	Y_init(y, numTicks, delta);

	// allocate device memory
	cudaMalloc( (void**) &x_dev, numBytes );
	cudaMalloc( (void**) &y_dev, numBytes );
	cudaMalloc( (void**) &z_dev, numBytes );

	// copy X, Y from host to device
	cudaMemcpy( x_dev, x, numBytes, cudaMemcpyHostToDevice);
	cudaMemcpy( y_dev, y, numBytes, cudaMemcpyHostToDevice);
	cudaMemcpy( z_dev, z, numBytes, cudaMemcpyHostToDevice);

	// GPU kernel
	int threadNum = 512;
	kernel <<< arraySize/threadNum + 1, threadNum >>> (x_dev, y_dev, z_dev, mux, muy, b4, arraySize);

	// copy Z from device to host
	cudaMemcpy( z, z_dev, numBytes, cudaMemcpyDeviceToHost);

	// show result
	array_show(z, numTicks);

	// memory free
	cudaFree(x_dev); free(x);
	cudaFree(y_dev); free(y);
	cudaFree(z_dev); free(z);

	return 0;
}