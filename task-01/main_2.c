#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

const double pi = 3.14159265358979323846;

void array_show(int32_t *a, int32_t numPoints) {
    for (int32_t i = 0; i < numPoints; ++i) {
    for (int32_t j = 0; j < numPoints; ++j)
        printf("%6d ", a[i*(numPoints) + j]);

    printf("\n");
    }
}

int32_t set_XY_arrays(int32_t *X, int32_t *Y, int32_t numPoints, double delta) {

    int32_t i, j, size = 0;
    double x, y;
    for (i = 0; i < numPoints; ++i) {
        for (j = 0; j < numPoints; ++j) {
            x = i * delta;
            y = j * delta;
            if (x*x + y*y < 1.0) {
                X[size] = j;
                Y[size] = i;
                size++;
            }
        }
    }

    return size;
}

void calculate(int32_t *X, int32_t *Y, int32_t *Z, int32_t numPoints, int32_t size_XY, double delta, double mux_, double muy_, double b4, int32_t numThreads) {


    int32_t s, n, index;
    double x0, x1, px0, px1, px2, mux;
    double y0, y1, py0, py1, py2, muy;

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (s = 0; s < size_XY; s++) {
        x0  = X[s] * delta;
        y0  = Y[s] * delta;
        mux = 2.0 * pi * mux_;
        muy = 2.0 * pi * muy_;
        px0 = 0.0; py0 = 0.0;
        for (n = 0; n <= 100000; n++) {

            if (x0*x0 + y0*y0 < 1.0) {
                x1  =  x0 * cos(mux) + px0 * sin(mux);
                px1 = -x0 * sin(mux) + px0 * cos(mux);
                y1  =  y0 * cos(muy) + py0 * sin(muy);
                py1 = -y0 * sin(muy) + py0 * cos(muy);
                px2 = px1 + b4 * (x1*x1*x1 - 3.0*x1*y1*y1);
                py2 = py1 - b4 * (y1*y1*y1 - 3.0*x1*x1*y1);
                x0  =  x1;  y0 =  y1;
                px0 = px2; py0 = py2;
            } else break;
        }
        index = numPoints * Y[s] + X[s];
        Z[index] = n - 1;
    }
}


int main(int argc, char const *argv[])
{

    char *p1, *p2;
    int32_t numThreads = (int32_t) strtol(argv[1], &p1, 10);  // first arg
    int32_t numPoints  = (int32_t) strtol(argv[2], &p2, 10);  // second arg

    double delta = 1.0 / (numPoints - 1);
    const double  b4 = 0.5;
    const double  mux = 0.32;
    const double  muy = 0.32;

    int32_t *X, *Y, *Z;
    X = (int32_t *) malloc( numPoints * numPoints * sizeof(int32_t) );
    Y = (int32_t *) malloc( numPoints * numPoints * sizeof(int32_t) );
    Z = (int32_t *) calloc( numPoints * numPoints,  sizeof(int32_t) );  // set zero

    int32_t size_XY;
    size_XY = set_XY_arrays(X, Y, numPoints, delta);

    calculate(X, Y, Z, numPoints, size_XY, delta, mux, muy, b4, numThreads);
    //array_show(Z, numPoints);

    return 0;
}