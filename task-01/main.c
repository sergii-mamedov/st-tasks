#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double pi = 3.14159265358979323846;


void array_show(int *a, int numTicks) {
    for (int i = 0; i <= numTicks; ++i) {
    for (int j = 0; j <= numTicks; ++j)
        printf("%6d ", a[i*(numTicks+1) + j]);

    printf("\n");
    }
}

/* calculate one particle stability */
int one_particle(double x0, double y0, double mux, double muy, double b4) {

    double px0 = 0.0, py0 = 0.0;
    double x1, px1, px2, y1, py1, py2;

    int n;
    for (n = 0; n <= 100000; n++) {
        //printf("%3i %+19.16f %+19.16f %19.16f\n", n, x0, y0, x0*x0 + y0*y0);
        if (x0*x0 + y0*y0 < 1.0) {
            x1  =  x0 * cos(2.0*pi*mux) + px0 * sin(2.0*pi*mux);
            px1 = -x0 * sin(2.0*pi*mux) + px0 * cos(2.0*pi*mux);
            y1  =  y0 * cos(2.0*pi*muy) + py0 * sin(2.0*pi*muy);
            py1 = -y0 * sin(2.0*pi*muy) + py0 * cos(2.0*pi*muy);
            px2 = px1 + b4 * (x1*x1*x1 - 3.0*x1*y1*y1);
            py2 = py1 - b4 * (y1*y1*y1 - 3.0*x1*x1*y1);
            //printf("    %+f %+f\n", px2, py2);
            x0  =  x1;  y0 =  y1;
            px0 = px2; py0 = py2;
        } else break;
    }

    return n-1;
}


/* calculate in a I Quadrant */
void calculate(int numTicks, double delta, double mux, double muy, double b4, int *Z) {

    int n, in_circle;
    double x0, y0;

    for (int y = 0; y <= numTicks; y++)
        for (int x = 0; x <= numTicks; x++)
        {
            x0 = delta * x;
            y0 = delta * y;
            in_circle = (x0*x0 + y0*y0 < 1) ? 1 : 0;

            if (in_circle == 1) {
                n = one_particle(x0, y0, mux, muy, b4);
                Z[y*(numTicks+1) + x] = n;
                //printf("%+5.3f %+5.3f  %5i\n", x0, y0, n);
            } else {
                Z[y*(numTicks+1) + x] = 0;
                //printf("%+5.3f %+5.3f      0\n", x0, y0);
            }
        }   
}


int main(int argc, char const *argv[])
{

    const double b4 = 0.5;
    const double mux = 0.32;
    const double muy = 0.32;
    const int    numTicks = 100;
    double delta = 1.0 / numTicks;

    int *Z;
    Z = (int *) malloc( (numTicks + 1) * (numTicks + 1) * sizeof(int) );

    calculate(numTicks, delta, mux, muy, b4, Z);
    //array_show(Z, numTicks);

    return 0;
}