void parse_command_line(int argc, char *argv[], double &T, double &X, int &h_t, int &h_x);
void calculate(double *H, double *HB, double *HL, double *HR, int n, int m, double h, double tau,
               double *V, double *VB, double *VL, double *VR);
void calculate_barsting(double *H, double *HB, double *HL, double *HR, int n, int m, double h, double tau,
                        double *V, double *VB, double *VL, double *VR);
