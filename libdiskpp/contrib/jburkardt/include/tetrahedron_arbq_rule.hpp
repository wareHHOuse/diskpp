#pragma once

double *kjacoypols3 ( double x, double y, double a, double b, int n );
double *klegeypols ( double x, double y, int n );
double *ortho3eva ( int degree, double xyz[] );
void r8mat_row_copy ( int m, int n, int i, double v[], double a[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_sum ( int n, double a[] );
double *ref_to_koorn ( double r[] );
void rule01 ( int n, double x[], double w[] );
void rule02 ( int n, double x[], double w[] );
void rule03 ( int n, double x[], double w[] );
void rule04 ( int n, double x[], double w[] );
void rule05 ( int n, double x[], double w[] );
void rule06 ( int n, double x[], double w[] );
void rule07 ( int n, double x[], double w[] );
void rule08 ( int n, double x[], double w[] );
void rule09 ( int n, double x[], double w[] );
void rule10 ( int n, double x[], double w[] );
void rule11 ( int n, double x[], double w[] );
void rule12 ( int n, double x[], double w[] );
void rule13 ( int n, double x[], double w[] );
void rule14 ( int n, double x[], double w[] );
void rule15 ( int n, double x[], double w[] );
void tetrahedron_arbq ( int degree, int n, double x[], double w[] );
//void tetrahedron_arbq_gnuplot ( int n, double x[], string header );
int tetrahedron_arbq_size ( int degree );
void tetrahedron_ref ( double v1[], double v2[], double v3[], double v4[] );
void timestamp ( );
