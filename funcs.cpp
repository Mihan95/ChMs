#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MY        0.01
#define M_TWOPI  (M_PI * 2.0)

using namespace std;

void parse_command_line(int argc, char* argv[], double &T, double &X, int &h_t, int &h_x)
{
    if (argc == 1)
        printf("Wrong parametrs of cmd");
    else
        if (argc > 1)
            T = atoi(argv[1]);
        if (argc > 2)
            X = atoi(argv[2]);
        if (argc > 3)
            h_t = atoi(argv[3]);
        if (argc > 4)
            h_x = atoi(argv[4]);
}

double ro(double t, double x)
{
    return exp(t) * (cos(M_PI*x*0.1) + 1.5);
}

double u(double t, double x)
{
    return cos(M_TWOPI * t) * sin(M_PI*x*x*0.01);
}

double ro(double x)
{
    return cos(M_PI*x*0.1) + 1.5;
}

double u(double x)
{
    return sin(M_PI*x*x*0.01);
}

void matrix_mult(double *b, double *a, double *c, double *d, double *x, int n)
{
    x[0] = a[0] * b[0] + c[0] * b[1];
    cout << x[0] << endl;
    for (int i = 1; i < n-1; i++)
    {
        x[i] = d[i-1] * b[i-1] + a[i] * b[i] + c[i] * b[i+1];
    }
    x[n-1] = d[n-2] * b[n-2] + a[n-1] * b[n-1];
}

void ThreeDiagSolve( double *b, double *a, double *c, double *d, int n )
{
    /*double *cb, *ca, *cd, *cc, *x;
    cb = new double [n];
    ca = new double [n];
    cd = new double [n-1];
    cc = new double [n-1];
    x = new double [n];
    for(int i = 0; i < n-1; i++)
    {
        ca[i] = a[i];
        cb[i] = b[i];
        cd[i] = d[i];
        cc[i] = c[i];
    }
    ca[n-1] = a[n-1];
    cb[n-1] = b[n-1];*/
/**************************************************************/
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
/*************************************************************/
    /*matrix_mult(b, ca, cc, cd, x, n);
    cout << "!!!!!!!!" << endl;
    for(int i = 0; i < n; i++)
    {
        cb[i] -= x[i];
        cout << cb[i] << endl;
    }
    cout << "!!!!!!!!" << endl;*/
}

double f0(double t, double x)
{
    double pi_10 = M_PI / 10.0;
    return exp(t) * (1.5 + cos(pi_10*x))
         + pi_10 * exp(t) * cos(M_TWOPI*t)
            * (x * 0.2 * (1.5 + cos(pi_10*x)) * cos(M_PI*x*x*0.01) - sin(pi_10*x) * sin(M_PI*x*x*0.01));
}

double f(double t, double x)
{
    double pi_10 = M_PI / 10.0;
    return M_PI*exp(t)*(cos(pi_10*x) + 1.5) * (cos(M_TWOPI*t) * sin(M_PI*x*x*0.01) * x*0.02* cos(M_TWOPI*t) * cos(M_PI*x*x*0.01)
                                              - 1.4 * sin(pi_10 * x) * exp(0.4 * t)
                                              - 2.0 * sin(M_TWOPI*t) * sin(M_PI*x*x*0.01))
         - MY * cos(M_TWOPI*t) * pi_10*0.2 * (cos(M_PI*x*x*0.01) - x*x*M_PI*0.02 * sin(M_PI*x*x*0.01));
}

void filling_H_0(double *HL, double *HR, double *H, double *HB, double tau, double h, int n)
{
    HR[0]  = tau * 0.5 * u(h) / h;
    H[0] = 1.0 - tau * 0.5 * u(0) / h;
    HB[0]  = f0(0, 0) + ro(0)
            - (tau * 0.5 / h)
            * ro(0) * (u(h) - u(0)) + (tau * 0.25 / h) * (2.0 * ro(2*h) * u(2*h)
                                                                   - 2.5 * ro(h) * u(h)
                                                                   + 2.0 * ro(0) * u(0)
                                                                 - 0.5 * ro(3*h) * u(3*h)
                                                          + ro(0) * (2.0 * u(2*h) - 2.5 * u(h) - 0.5 * ro(3*h)));
    for (int j = 1; j < n - 1; j++)
    {
        H[j]    = 1.0;
        HR[j]   = - (0.25 * tau / h) * (u(j*h) + u((j-1)*h));
        HL[j-1] = (0.25 * tau / h) * (u(j*h) + u((j+1)*h));
        HB[j]    = f0(0, j*h) + ro(j*h) * (1.0 - (0.25 * tau / h) * u((j+1)*h) - u((j-1)*h));
    }

    H[n-1]  = 1.0 + (tau * 0.5 / h) * u((n-1)*h);
    HL[n-2] = - (0.5 * tau / h) * u((n-2)*h);
    HB[n-1]  = f0(0, (n-1)*h) + ro((n-1)*h) - (0.5 * tau / h) * ro((n-1)*h) * (u((n-1)*h) - u((n-2)*h))
            - (0.25 * tau / h) * (2.0 * ro((n-1)*h) * u((n-1)*h) - 2.5 * ro((n-2)*h) * u((n-2)*h)
                                + 2.0 * ro((n-3)*h) * u((n-3)*h) - 0.5 * ro((n-4)*h) * u((n-4)*h)
                                + ro((n-1)*h) * (u((n-3)*h) - 2.5 * u((n-2)*h) - 0.5 * u((n-4)*h)));
}

void filling_H(double *HL, double *HR, double *H, double *HB, double tau, double h, int n, int i)
{
    int i_tau = i * tau;

    HR[0]  = tau * 0.5 * u(i_tau, h) / h;
    H[0] = 1.0 - tau * 0.5 * u(i_tau, 0) / h;
    HB[0]  = f0(i_tau, 0) + HB[0]
            - (tau * 0.5 / h)
            * HB[0] * (u(i_tau, h) - u(i_tau, 0)) + (tau * 0.25 / h) * (2.0 * HB[2] * u(i_tau, 2*h)
                                                                   - 2.5 * HB[1] * u(i_tau, h)
                                                                   + 2.0 * HB[0] * u(i_tau, 0)
                                                                 - 0.5 * HB[3] * u(i_tau, 3*h)
                                                          + HB[0] * (2 * u(i_tau, 2*h) - 2.5 * u(i_tau, h) - 0.5 * ro(i_tau, 3*h)));

    double HM = HB[n-1], HM_1 = HB[n-2], HM_2 = HB[n-3], HM_3 = HB[n-4];

    for (int j = 1; j < n - 1; j++)
    {
        HB[j] = f0(i_tau, j*h) + HB[j] * (1.0 - (0.25 * tau / h) * (u(i_tau, (j+1)*h) - u(i_tau, (j-1)*h)));
        H[j]  = 1.0;
        HR[j] = - (0.25 * tau / h) * (u(i_tau, j*h) + u(i_tau, (j-1)*h));
        HL[j-1] = (0.25 * tau / h) * (u(i_tau, j*h) + u(i_tau, (j+1)*h));
    }

    HB[n-1] = f0(i_tau, (n-1)*h) + HM - (0.5 * tau / h) * HM * (u(i_tau, (n-1)*h) - u(i_tau, (n-2)*h))
            - (0.25 * tau / h) * (2.0 * HM * u(i_tau, (n-1)*h) - 2.5 * HM_1 * u(i_tau, (n-2)*h)
                                  + 2.0 * HM_2 * u(i_tau, (n-3)*h) - 0.5 * HM_3 * u(i_tau, (n-4)*h)
                                  + HM * (-0.5 * u(i_tau, (n-4)*h) - 2.5 * u(i_tau, (n-2)*h) + u(i_tau, (n-3)*h)));
    H[n-1]  = 1.0 + (0.5 * tau / h) * u(i_tau, (n-1)*h);
    HL[n-2] = - (0.5 * tau / h) * u(i_tau, (n-2)*h);
}

void calculate(double *H, double *HB, double *HL, double *HR, int n, int m, double h, double tau)
{
    filling_H_0(HL, HR, H, HB, tau, h, n);
    ThreeDiagSolve(HB, H, HR, HL, n);

    for (int i = 1; i < m-1; i++)
    {
        filling_H(HL, HR, H, HB, tau, h, n, i);
        ThreeDiagSolve(HB, H, HR, HL, n);
    }

    double res = -1.0;
    double diff = 0.0;
    for (int i = 0; i < n; i++)
    {
        diff = fabs(HB[i] - ro((m-1)*tau, i*h));

        if (diff > res)
            res = diff;
    }
    cout << "\nResidual is " << res << endl << endl;
}




























