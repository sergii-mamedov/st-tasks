#include <stdio.h>
#include <math.h>

const double pi = 3.14159265358979323846;

int calculate(double x0, double y0, double mux, double muy, double b4) {

    double px0 = 0.0, py0 = 0.0;
    double x1, px1, px2, y1, py1, py2;

    int n;
    for (n = 0; n <= 10000; n++) {
        //printf("%3i %+5.3f %+5.3f %16.14f\n", n, x0, y0, x0*x0 + y0*y0);
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


int main(int argc, char const *argv[])
{

    double b4;                      // need create array
    double mux = 0.31;
    double muy = 0.32;
    int lower = -100, upper = 100;
    double delta = 0.01;

    int n, in_circle;
    double x0, y0;
    for (int i = 0; i < 1; i++)
    {
        b4 = 0.5;
        for (int x = lower; x <= upper; x++)
            for (int y = lower; y <= upper; y++)
            {
                x0 = delta*x;
                y0 = delta*y;
                in_circle = (x0*x0 + y0*y0 < 1) ? 1 : 0;

                if (in_circle == 1) {
                    n = calculate(x0, y0, mux, muy, b4);
                    printf("%5.3e %+5.3f %+5.3f  %5i\n", b4, x0, y0, n);
                } //else printf("%5.3e %+5.3f %+5.3f      0\n", b4, x0, y0);
            }
    }

    return 0;
}