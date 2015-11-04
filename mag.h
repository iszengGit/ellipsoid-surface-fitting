
typedef enum bool{
    false,
    true
}bool;

typedef struct magcoord{
    short int x;
    short int y;
    short int z;
} magcoord;

magcoord rawCoo[8];
#define poolCooLen 300
magcoord poolCoo[poolCooLen];

typedef struct magresult {
    double x0;
    double y0;
    double z0;
    double A;
    double B;
    double R;
}magresult;

typedef struct fitresult {
    double mini;        // The sum of pow of difference of two side of equation
    magresult u;        // The initial value
    magresult delt;     // The margin of error
}fitresult;

#define N 6

void OriginalValueCalculate(void);
double calcDet(double a[][N], int n);
void calcAdj(double ori[][N], int n, double ans[][N]);
bool calcInv(double src[][N], int n, double des[][N]);
