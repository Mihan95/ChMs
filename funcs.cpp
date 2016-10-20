#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MY 1

void parse_command_line(double &T, double &X, double &h_t, double &h_x)
{
    if (argc == 1)
        print("Wrong parametrs of cmd");
    else
        if (argc > 1)
            T = atoi(arg[1]);
        if (argc > 2)
            X = atoi(arg[2]);
        if (argc > 3)
            h_t = atoi(arg[3]);
        if (argc > 4)
            h_x = atoi(arg[4]);

}

double ro(double t, double x)
{
    return exp(t) * (cos(M_PI*x*0.1) + 1.5);
}

double u(double t, double x)
{
    return cos(M_TWOPI*t) * sin(M_PI*x*x*0.01);
}

double ro(double x)
{
    return cos(M_PI*x*0.1) + 1.5;
}

double u(double x)
{
    return sin(M_PI*x*x*0.01);
}

void ThreeDiagSolve( double *b, double *a, double *c, double *d, int n )
{
    c[0] = c[0] / a[0];
    for( int i = 1; i < n - 1; i++ )
    {
        a[i] = a[i] - d[i - 1] * c[i - 1];
        c[i] = c[i] / a[i];
    }
    a[n - 1] = a[n - 1] - d[n - 2] * c[n - 2];

    b[0] = b[0] / a[0];
    for( int i = 1; i < n; i++ )
        b[i] = (b[i] - d[i - 1] * b[i - 1]) / a[i];

    for( int i = n - 2; i > -1; i-- )
        b[i] = b[i] - c[i] * b[i + 1];
}

double f0(double t, double x)
{
    double pi_10 = M_PI / 10.0;
    return 1.0 - cos(M_TWOPI*t) * sin(M_PI*x*x*0.01) * pi_10 * sin(pi_10*x)
               / (cos(pi_10*x) + 1.5)
         + cos(M_PI*x*x*0.01) * p_10 * 0.2 * x * cos(M_TWOPI*t);
}

double f(double t, double x)
{
    double pi_10 = M_PI / 10.0;
    return M_TWOPI * sin(M_PI*x*x*0.01) * (cos(M_TWOPI*t) * cos(M_PI*x*x*0.01) * x * 0.01
                                         - sin(M_TWOPI*t))
         - (sin(pi_10*x) * pi_10 / (cos(pi_10*x) + 1.5)
         + MY * cos(M_TWOPI*t) * pi_10 * 0.2 * (cos(M_PI*x*x*0.01) - x*x*pi_10*0.02*sin(M_PI*x*x*0.01)))
           / (t + log(cos(pi_10*x) + 1.5));
}

void filling_H_0(double *HL, double *HR, double *H, double *HB, double tau, double h, int n)
{
    HR[0]  = tau * 0.5 * u(h) / h;
    H[0] = 1.0 - tau * 0.5 * u(0) / h;
    HB[0]  = f0(0, 0) + ro(0)
            - (tau * 0.5 / h)
            * (ro(0) * (u(h) - u(0)) + 0.5 * (2.0 * ro(2*h) * u(2*h)
                                                                   - 2.5 * ro(h) * u(h)
                                                                   + 2.0 * ro(0) * u(0)
                                                                 - 0.5 * ro(3*h) * u(3*h)
                                                          + ro(0) * (2 * u(2*h) - 2.5 * u(h) - 0.5 * ro(3*h))));

    for (int j = 1; j < n - 1; j++)
    {
        H[j]    = 1.0;
        HR[j]   = (0.25 * tau / h) * (-1) * (u(j*h) + u((j-1)*h));
        HL[j-1] = (0.25 * tau / h) * (u(j*h) + u((j+1)*h));
        HB[j]    = f0(0, j*h) * ro(j*h) * (1 - 0.25 * tau * u((j+1)*h) - u((j-1)*h)) / h;
    }

    H[n-1]  = 1 + tau * 0.5 * u((n-1)*h) / h;
    HL[n-2] = (-1) * 0.5 * tau  * u((n-2)*h) / h;
    HB[n-1]  = f0(0, (n-1)*h) * ro((n-1)*h) - (0.5 * tau / h) * (ro((n-1)*h) * (u((n-1)*h) - u((n-2)*h)))
            - (0.25 * tau / h) * (2.0 * ro((n-1)*h) * u((n-1)*h) - 2.5 * ro((n-2)*h) * u((n-2)*h)
                                + 2.0 * ro((n-3)*h) * u((n-3)*h) - 0.5 * ro((n-4)*h) * u((n-4)*h)
                                + ro((n-1)*h) * ((-2.5) * u((n-2)*h) + 0.5 * u((n-4)*h) + u((n-3)*h)));
}

void filling_H(double *HL, double *HR, double *H, double *HB, double tau, double h, int n, int i)
{
    int i_tau = i * tau;

    HR[0]  = tau * 0.5 * u(i_tau, h) / h;
    H[0] = 1.0 - tau * 0.5 * u(i_tau, 0) / h;
    HB[0]  = f0(i_tau, 0) + HB[0]
            - (tau * 0.5 / h)
            * (HB[0] * (u(i_tau, h) - u(i_tau, 0)) + 0.5 * (2.0 * HB[2] * u(i_tau, 2*h)
                                                                   - 2.5 * HB[1] * u(i_tau, h)
                                                                   + 2.0 * HB[0] * u(i_tau, 0)
                                                                 - 0.5 * HB[3] * u(i_tau, 3*h)
                                                          + HB[0] * (2 * u(i_tau, 2*h) - 2.5 * u(i_tau, h) - 0.5 * ro(i_tau, 3*h))));

    double HM = HB[n-1], HM_1 = HB[n-2], HM_2 = HB[n-3], HM_3 = HB[n-4];

    for (int j = 0; i < n - 1; j++)
    {

    }

}






























