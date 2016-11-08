#include <stdio.h>
#include "funcs.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{

    double T = 1,  X = 10;
    int    m = 10, n = 10;

    parse_command_line(argc, argv, T, X, n, m);

    double *H  = new double [n];
    double *HL = new double [n-1];
    double *HR = new double [n-1];
    double *HB = new double [n];

    double *V  = new double [n-2];
    double *VL = new double [n-3];
    double *VR = new double [n-3];
    double *VB = new double [n];

    double h   = X / (n - 1);
    double tau = T / (m - 1);
    cout << "h = " << h << " tau =  " << tau << endl;

    calculate(H, HB, HL, HR, n, m, h, tau, V, VB, VL, VR);

    delete [] H;
    delete [] HL;
    delete [] HR;
    delete [] HB;
    delete [] V;
    delete [] VL;
    delete [] VR;
    delete [] VB;
}
