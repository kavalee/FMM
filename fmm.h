//
// Created by kavalee on 4/6/20.
//
typedef struct CauchyMultiplier{
    double* s;
    double* t;
    double* g;
    int n;
    int p;
    double sMin;
    double sMax;
    double tMin;
    double tMax;
    double sH;
    double tH;
    double** BINOMIAL_CACHE;
    unsigned long adds;
    unsigned long muls;
} CauchyMultiplier;

double* slowMultiply(CauchyMultiplier* cm);
double* fastMultiply(CauchyMultiplier* cm);
CauchyMultiplier* newCauchyMultiplier(double* s, double* t, double* g, int n, int p);
void freeCauchyMultiplier(CauchyMultiplier* cm);
