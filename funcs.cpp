#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MY        0.01
#define M_TWOPI  (M_PI * 2.0)
#define GAMMA 1.4

using namespace std;

void parse_command_line(int argc, char* argv[], double &T, double &X, int &h_x, int &h_t)
{
    if (argc == 1)
        printf("Wrong parametrs of cmd");
    else
        if (argc > 1)
            T = atoi(argv[1]);
        if (argc > 2)
            X = atoi(argv[2]);
        if (argc > 3)
            h_x = atoi(argv[3]);
        if (argc > 4)
            h_t = atoi(argv[4]);
}


/********Непрерывное решение*****************/
/*double ro (double t, double x)
{
  return exp (t) * (cos (M_PI * x * 0.1) + 1.5);
}

double u (double t, double x)
{
  return cos (2. * M_PI * t) * sin (M_PI * 0.01 * x * x);
}

double f0 (double t, double x)
{
  double ro_ = ro (t, x);
  double u_ = u (t, x);

  return   cos(M_TWOPI*t) * ro(t, x) * cos(M_PI*x*x*0.01) * M_PI * 0.02 * x
         - exp(t) * u(t, x) * sin(M_PI*x*0.1) * M_PI * 0.1
         + ro(t, x);
}

double f (double t, double x)
{
    double ro_ = ro(t, x);

  return (- ro_ * sin(M_PI*x*x*0.01) * sin(M_TWOPI*t) * M_TWOPI
         + ro_ * u(t, x) * cos(M_TWOPI*t) * cos(M_PI*x*x*0.01) * M_PI * x / 50.0
         - 1.4 * M_PI * 0.1 * exp(1.4*t) * sin(M_PI*x*0.1) * pow((cos(M_PI*x*0.1) + 1.5), 0.4)
         - MY * cos(M_TWOPI*t) * M_PI * 0.02 * (cos(M_PI*x*x*0.01) - M_PI*x*x*0.02 * sin(M_PI*x*x*0.01)));
}

double ro(double x)
{
    return cos(M_PI*x*0.1) + 1.5;
}

double u(double x)
{
    return sin(M_PI*x*x*0.01);
}*/
/********************************************/

/********Первое разрывное решение************/
/*double u(double t, double x)
{
    if ((fabs(x - 10.0) < 1e-15) || (fabs(x) < 1e-15))
            return 0.0;
}

double ro(double x)
{
    if (x < 4.5 || x > 5.5)
        return 1.0;
    else return 2.0;
}

double u(double x)
{
    return 0.0;
}

double f(double t, double x)
{
    return 0.0;
}

double f0(double t, double x)
{
    return 0.0;
}*/

/*********************************************/

/********Второе разрывное решение*************/

double f(double t, double x)
{
    return 0.0;
}

double f0(double t, double x)
{
    return 0.0;
}

double u(double x)
{
    if (x < 4.5 || x > 5.5)
        return 0.0;
    else return 1.0;
}

double u(double t, double x)
{
    if ((fabs(x - 10.0) < 1e-15) || (fabs(x) < 1e-15))
            return 0.0;
}

double ro(double x)
{
    return 1.0;
}



/*********************************************/

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



void filling_H_0(double *HL, double *HR, double *H, double *HB, double tau, double h, int n)
{
    HR[0]  = tau * 0.5 * u(h) / h;
    H[0] = 1.0;
    HB[0]  = f0(0, 0) * tau + ro(0) - (tau * 0.5 / h) * ro(0) * u(h)
            + (tau * 0.25 / h) * (- 2.5 * ro(h) * u(h) + 2.0 * ro(2*h) * u(2*h) - 0.5 * ro(3*h) * u(3*h)
                                                          + ro(0) * (2.0 * u(2*h) - 2.5 * u(h) - 0.5 * u(3*h)));
    for (int j = 1; j < n - 1; j++)
    {
        H[j]    = 1.0;
        HR[j]   = (0.25 * tau / h) * (u(j*h) + u((j+1)*h));
        HL[j-1] = - (0.25 * tau / h) * (u(j*h) + u((j-1)*h));
        HB[j]    = f0(0, j*h) * tau + ro(j*h) * (1.0 - (0.25 * tau / h) * (u((j+1)*h) - u((j-1)*h)));
    }

    H[n-1]  = 1.0;
    HL[n-2] = - (0.5 * tau / h) * u((n-2)*h);
    HB[n-1]  = f0(0, (n-1)*h) * tau + ro((n-1)*h) - (0.5 * tau / h) * ro((n-1)*h) * (- u((n-2)*h))
            - (0.25 * tau / h) * (- 2.5 * ro((n-2)*h) * u((n-2)*h)
                                + 2.0 * ro((n-3)*h) * u((n-3)*h) - 0.5 * ro((n-4)*h) * u((n-4)*h)
                                + ro((n-1)*h) * (2.0 * u((n-3)*h) - 2.5 * u((n-2)*h) - 0.5 * u((n-4)*h)));
}

void filling_H(double *HL, double *HR, double *H, double *HB, double tau, double h, int n, int i, double *V)
{
    double i_tau = i * tau;

    HR[0]  = tau * 0.5 * V[1] / h;
    H[0] = 1.0;

    HB[0]  = f0(i_tau, 0) * tau + HB[0] - (0.5 * tau / h) * HB[0] * V[1]
            + (0.25 * tau / h) * (- 2.5 * HB[1] * V[1]
                                  + 2.0 * HB[2] * V[2] - 0.5 * HB[3] * V[3]
                                  + HB[0] * (-0.5 * V[3] - 2.5 * V[1] + 2.0 * V[2]));

    double HM = HB[n-1], HM_1 = HB[n-2], HM_2 = HB[n-3], HM_3 = HB[n-4];

    for (int j = 1; j < n - 1; j++)
    {
        HB[j] = f0(i_tau, j*h) * tau + HB[j] * (1.0 - (0.25 * tau / h) * (V[j+1] - V[j-1]));
        H[j]  = 1.0;
        HR[j] = (0.25 * tau / h) * (V[j] + V[j+1]);
        HL[j-1] = - (0.25 * tau / h) * (V[j] + V[j-1]);
    }

    HB[n-1] = f0(i_tau, (n-1)*h) * tau + HM - (0.5 * tau / h) * HM * (- V[n-2])
            - (0.25 * tau / h) * (- 2.5 * HM_1 * V[n-2]
                                  + 2.0 * HM_2 * V[n-3] - 0.5 * HM_3 * V[n-4]
                                  + HM * (-0.5 * V[n-4] - 2.5 * V[n-2] + 2.0 * V[n-3]));
    H[n-1]  = 1.0;
    HL[n-2] = - (0.5 * tau / h) * V[n-2];
}

void filling_V_0(double *VL, double *V, double *VR, double *VB, double tau, double h, int n, double *H)
{
    double mu = 1.0 / H[0];
    double r = 0;
    for (int i = 1; i < n; i++)
    {
        r = 1.0 / H[i];
        if (r > mu)
            mu = r;
    }
    mu = mu * MY;

    VB[0]   = 0.0;
    VB[n-1] = 0.0;

    VR[0] = (tau / (6.0 * h)) * (u(h) + u(2*h)) - mu * tau / (h * h);
    V[0]  = 1.0 + tau * 2.0 * mu / (h * h);
    VB[1] = u(h) - tau * (pow(H[2], 0.4) - pow(H[0], 0.4)) / (2.0 * h)
            - tau * (mu - MY/H[1]) * (u(2*h) - 2 * u(h) + u(0)) / (h * h) + f(0, h)*tau;

    for (int i = 1; i < n-3; i++)
    {
        V[i]    = 1.0 + tau * 2.0 * mu / (h * h);
        VR[i]   =   (tau / (6.0 * h)) * (u((i+1) * h) + u((i+2) * h)) - mu * tau / (h * h);
        VL[i-1] = - (tau / (6.0 * h)) * (u((i-1) * h) + u(i * h)) - mu * tau / (h * h);
        VB[i+1] = u((i+1) * h) - tau * (pow(H[i+2], 0.4) - pow(H[i], 0.4)) / (2.0 * h)
                - tau * (mu - MY/H[i+1]) * (u((i+2)*h) - 2 * u((i+1)*h) + u(i*h)) / (h * h) + f(0, (i+1)*h)*tau;
    }

     V [n-3] = 1.0 + tau * 2.0 * mu / (h * h);
     VL[n-4] = - (tau / (6.0 * h)) * (u((n-4)*h) + u((n-3)*h)) - mu * tau / (h * h);
     VB[n-2] = u((n-2) * h) - tau * (pow(H[n-1], 0.4) - pow(H[n-3], 0.4)) / (2.0 * h)
             - tau * (mu - MY/H[n-2]) * (u((n-1)*h) - 2 * u((n-2)*h) + u((n-3)*h)) / (h * h) + f(0, (n-2)*h)*tau;
}

void filling_V(double *VL, double *V, double *VR, double *VB, double tau, double h, int n, int j, double *H)
{

       double r = 0;
       double j_tau = j * tau;
       double mu = 1.0 / H[0];

       for (int i = 1; i < n; i++)
       {
           r = 1.0 / H[i];
           if (r > mu)
               mu = r;
       }
       mu = MY * mu;


       VB[0]   = 0.0;
       VB[n-1] = 0.0;

       VR[0] = (tau / (6.0 * h)) * (VB[1] + VB[2]) - mu * tau / (h * h);
       V[0]  = 1.0 + tau * 2.0 * mu / (h * h);

       double prevVB = 0.0, tmpVB = 0.0;
       double prev_prevVB = 0.0;

       prevVB = VB[1];
       tmpVB = 0;
       VB[1] = VB[1] - tau * (pow(H[2], 0.4) - pow(H[0], 0.4)) / (2.0 * h)
               - tau * (mu - MY/H[1]) * (VB[2] - 2 * VB[1]) / (h * h) + f(j_tau, h)*tau;

       for (int i = 1; i < n-3; i++)
       {
           V[i]    = 1.0 + tau * 2.0 * mu / (h * h);
           VR[i]   =   (tau / (6.0 * h)) * (VB[i+1] + VB[i+2]) - mu * tau / (h * h);
           VL[i-1] = - (tau / (6.0 * h)) * (prev_prevVB + prevVB) - mu * tau / (h * h);

           tmpVB = VB[i+1];
           VB[i+1] = VB[i+1] - tau * (pow(H[i+2], 0.4) - pow(H[i], 0.4)) / (2.0 * h)
                   - tau * (mu - MY/H[i+1]) * (VB[i+2] - 2 * VB[i+1] + prevVB) / (h * h) + f(j_tau, (i+1)*h)*tau;
           prev_prevVB = prevVB;
           prevVB = tmpVB;
       }

        V [n-3] = 1.0 + tau * 2.0 * mu / (h * h);
        VL[n-4] = - (tau / (6.0 * h)) * (prev_prevVB + prevVB) - mu * tau / (h * h);
        VB[n-2] = VB[n-2] - tau * (pow(H[n-1], 0.4) - pow(H[n-3], 0.4)) / (2.0 * h)
                - tau * (mu - MY/H[n-2]) * ( - 2 * VB[n-2] + prevVB) / (h * h) + f(j_tau, (n-2)*h)*tau;
}

double integral(double *f, int n, double h)
{
    double integral = 0.0;

    for(int i = 1; i < n-1; i++)
    {
        integral += f[i];
    }
    integral = integral + (f[0] + f[n-1]) / 2;
    integral *= h;
}

void calculate_barsting(double *H, double *HB, double *HL, double *HR, int n, int m, double h, double tau,
                        double *V, double *VB, double *VL, double *VR)
{
    cout << "MY = " << MY << endl;
    FILE *fh = fopen("h.txt", "w");
    FILE *fv = fopen("v.txt", "w");

    fprintf(fh, "%d %d %3.15f %3.15f", n, m, h, tau); fprintf(fh, "\n");
    fprintf(fv, "%d %d %3.15f %3.15f", n, m, h, tau); fprintf(fv, "\n");

    filling_H_0(HL, HR, H, HB, tau, h, n);
    ThreeDiagSolve(HB, H, HR, HL, n);

    filling_V_0(VL, V, VR, VB, tau, h, n, HB);
    ThreeDiagSolve(VB+1, V, VR, VL, n-2);

    for (int i = 0; i < n; i++)
    {
        fprintf(fh, "%3.15f\n", HB[i]);
    }
    fprintf(fh, "\n");

    for (int i = 0; i < n; i++)
    {
        fprintf(fv, "%3.15f\n", VB[i]);
    }
    fprintf(fv, "\n");

    for (int i = 1; i < m-1; i++)
    {
        filling_H(HL, HR, H, HB, tau, h, n, i, VB);
        ThreeDiagSolve(HB, H, HR, HL, n);

        filling_V(VL, V, VR, VB, tau, h, n, i, HB);
        ThreeDiagSolve(VB+1, V, VR, VL, n-2);

        /*for (int i = 0; i < n; i++)
        {
            fprintf(fh, "%3.15f\n", HB[i]);
        }
        fprintf(fh, "\n");

        for (int i = 0; i < n; i++)
        {
            fprintf(fv, "%3.15f\n", VB[i]);
        }
        fprintf(fv, "\n");*/
    }

    double integr_ = integral(HB, n, h);
    cout << "Mass is " << integr_ << endl;

    /***Первое разрывное решение***/
    cout << "Mass residual is " << fabs(11 - integr_) << endl;
    /******************************/

    /***Второе разрывное решение***/
    cout << "Mass residual is " << fabs(10 - integr_) << endl;
    /******************************/

}

void calculate(double *H, double *HB, double *HL, double *HR, int n, int m, double h, double tau,
               double *V, double *VB, double *VL, double *VR)
{

    cout << "MY = " << MY << endl;
    FILE *fh = fopen("h.txt", "w");
    FILE *fv = fopen("v.txt", "w");

    fprintf(fh, "%d %d %3.15f %3.15f", n, m, h, tau); fprintf(fh, "\n");
    fprintf(fv, "%d %d %3.15f %3.15f", n, m, h, tau); fprintf(fv, "\n");

    filling_H_0(HL, HR, H, HB, tau, h, n);
    ThreeDiagSolve(HB, H, HR, HL, n);

    filling_V_0(VL, V, VR, VB, tau, h, n, HB);
    ThreeDiagSolve(VB+1, V, VR, VL, n-2);

    for (int i = 0; i < n; i++)
    {
        fprintf(fh, "%3.15f\n", HB[i]);
    }
    fprintf(fh, "\n");

    for (int i = 0; i < n; i++)
    {
        fprintf(fv, "%3.15f\n", VB[i]);
    }
    fprintf(fv, "\n");

    for (int i = 1; i < m-1; i++)
    {
        filling_H(HL, HR, H, HB, tau, h, n, i, VB);
        ThreeDiagSolve(HB, H, HR, HL, n);

        filling_V(VL, V, VR, VB, tau, h, n, i, HB);
        ThreeDiagSolve(VB+1, V, VR, VL, n-2);

        /*for (int i = 0; i < n; i++)
        {
            fprintf(fh, "%3.15f\n", HB[i]);
        }
        fprintf(fh, "\n");

        for (int i = 0; i < n; i++)
        {
            fprintf(fv, "%3.15f\n", VB[i]);
        }
        fprintf(fv, "\n");*/
    }

    /*filling_H(HL, HR, H, HB, tau, h, n, m-2);
    ThreeDiagSolve(HB, H, HR, HL, n);*/

    double res = -1.0;
    double diff = 0.0;
    for (int i = 0; i < n; i++)
    {
       //diff = fabs(HB[i] - ro((m-1)*tau, i*h));
        if (diff > res)
            res = diff;
    }
    cout << "\nResidual H is " << res << endl << endl;

    /*filling_V(VL, V, VR, VB, tau, h, n, m-2);
    ThreeDiagSolve(VB+1, V, VR, VL, n-2);*/

    res = -1.0;

    for (int i = 0; i < n; i++)
    {
        diff = fabs(VB[i] - u((m-1)*tau, i*h));
        if (diff > res)
        res = diff;
    }
    cout << "\nResidual V is " << res << endl << endl;
}
