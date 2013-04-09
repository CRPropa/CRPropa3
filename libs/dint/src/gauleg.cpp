#include <math.h>

void Gauleg(const double x1, const double x2, double x[], double w[], 
            const int n)
{
    int m;
    int j;
    int i;
    double z1;
    double z;
    double xm;
    double xl;
    double pp;
    double p3;
    double p2;
    double p1;
    const double EPS = 3.0e-11; /* precision parameter */

    m = (n + 1)/2;
    xm = 0.5*(x2 + x1);
    xl = 0.5*(x2 - x1);
    for (i = 0; i < m; i++)
    {
        z = cos(3.141592654*(i+0.75)/(n+0.5));
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 0; j < n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j + 1.0)*z*p2 - j*p3)/(j+1);
            }
            pp = n*(z*p1 - p2)/(z*z - 1.0);
            z1 = z;
            z = z1 - p1/pp;
        } while (fabs(z - z1) > EPS);

        x[i] = xm - xl*z;
        x[n-1-i] = xm + x1*z;
        w[i] = 2.0*xl/((1.0 - z*z)*pp*pp);
        w[n-1-i] = w[i];
    }
}

