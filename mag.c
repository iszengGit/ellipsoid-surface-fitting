#include "stdio.h"
#include "mag.h"
#include "math.h"

int usedSub[poolCooLen];
int subFlag = 0;

void updateRawCoo()
{
    int i, j, k, iii, jjj;
    int s;
    magcoord mtemp;

    for (i = 0; i < 8; i++) {
        for (iii = 0; iii < poolCooLen; iii++) {
            for (jjj = 0; jjj < subFlag; jjj++) {
                if (iii == usedSub[jjj]) {  
                   goto NEXTiii;
                }
            }
            mtemp = poolCoo[iii];       // New reference point
            usedSub[subFlag++] = iii;   // Record it as used subscript
            break;
NEXTiii:
            continue;
        }
        for (j = 0, s = poolCooLen + 1; j < poolCooLen && i < 6; j++) {
            for (k = 0; k < subFlag; k++) {
                if (j == usedSub[k]) {
                    goto NEXTj;
                }
            }
            if (i == 0) {
                if (poolCoo[j].x < mtemp.x) {
                    mtemp.x = poolCoo[j].x;
                    s = j;
                }
            }
            else if (i == 1) {
                if (poolCoo[j].x > mtemp.x) {
                    mtemp.x = poolCoo[j].x;
                    s = j;
                }
            }
            else if (i == 2) {
                if (poolCoo[j].y < mtemp.y) {
                    mtemp.y = poolCoo[j].y;
                    s = j;
                }
            }
            else if (i == 3) {
                if (poolCoo[j].y > mtemp.y) {
                    mtemp.y = poolCoo[j].y;
                    s = j;
                } 
            }
            else if (i == 4) {
                if (poolCoo[j].z < mtemp.z) {
                    mtemp.z = poolCoo[j].z;
                    s = j;
                }
            }
            else if (i == 5) {
                if (poolCoo[j].z > mtemp.z) {
                    mtemp.z = poolCoo[j].z;
                    s = j;
                } 
            }
NEXTj:
            continue;            
        }
        if (i < 6) {
            if (s == poolCooLen + 1) {
                rawCoo[i] = mtemp;
            } else {
                rawCoo[i] = poolCoo[s];
                usedSub[subFlag++] = s;
            }
        } else {
            rawCoo[i] = mtemp;        
        }
    }
    
}




double matd[8][6];
double transMatd[6][8];
double matc[8][1];
double matw[6][1];

magresult initu;            // INITIAL VALUE

/*
 * Description: Calculate initial value 
 *              d * w = c
 *              w = ?
 */
void initMatd()
{
    int i, j;

    // Update rawCoo
    updateRawCoo(); 

    for (i = 0; i < 8; i++) {
        matd[i][0] = -pow(rawCoo[i].y, 2);
        matd[i][1] = -pow(rawCoo[i].z, 2);
        matd[i][2] = 2*rawCoo[i].x;
        matd[i][3] = 2*rawCoo[i].y;
        matd[i][4] = 2*rawCoo[i].z;
        matd[i][5] = 1;
    }

    for (i = 0; i < 8; i++)
        for (j = 0; j < 6; j++){
            transMatd[j][i] = matd[i][j];
        }

    for (i = 0; i < 8; i++)
        matc[i][0] = pow(rawCoo[i].x, 2);
}

void calcInitialValue()
{
    double prod1[6][6] = {{0.0}}, temp = 0;
    double invProd1[6][6] = {{0.0}};
    double prod2[6][6] = {{0.0}};
    int i, j, k;

    // Calculate the prod1
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++) {
            for (k = 0, temp = 0; k < 8; k++){
                temp += transMatd[j][k]*matd[k][i];
            }
            prod1[j][i] = temp;
        }
    
    // Calculate the inverse
    calcInv(prod1, 6, invProd1);

    // Calculate the prod2
    for (i = 0; i < 8; i++)
        for (j = 0; j < 6; j++){
            for (k = 0, temp = 0; k < 6; k++){
                temp += invProd1[j][k]*transMatd[k][i];            
            }
            prod2[j][i] = temp;
        }
    
    // Calculate matw
    for (i = 0; i < 6; i++) {
        for (j = 0, temp = 0; j < 8; j++) {
            temp += prod2[i][j]*matc[j][0];
        }
        matw[i][0] = temp;
    }

    // Get initial value
    initu.x0 = matw[2][0];
    initu.y0 = matw[3][0]/matw[0][0];
    initu.z0 = matw[4][0]/matw[1][0];
    initu.A = sqrt(matw[0][0]);
    initu.B = sqrt(matw[1][0]); 
    initu.R = sqrt(matw[5][0] + matw[3][0]*initu.y0 + matw[4][0]*initu.z0 + pow(matw[2][0], 2));
}

double math[8][6];
double transMath[6][8];
double maty[8][1];
double matDeltau[6][1];

double miniSum;
magresult deltau = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};       // MARGIN OF ERROR

/*
 * Description: Calculate margin of error
 *              h * Deltau = y
 *              Deltau = ?
 */
void initMath()
{
    int i, j;

    // Update rawCoo
    updateRawCoo();

    // Update INITIAL VALUE
    initu.x0 += deltau.x0;
    initu.y0 += deltau.y0;
    initu.z0 += deltau.z0;
    initu.A += deltau.A;
    initu.B += deltau.B;
    initu.R += deltau.R;
    // Initialize math[8][6]
    for (i = 0; i < 8; i++) {
        math[i][0] = 2*(rawCoo[i].x - initu.x0);
        math[i][1] = 2*pow(initu.A, 2)*(rawCoo[i].y - initu.y0);
        math[i][2] = 2*pow(initu.B, 2)*(rawCoo[i].z - initu.z0);
        math[i][3] = (-2)*initu.A*pow(rawCoo[i].y - initu.y0, 2);
        math[i][4] = (-2)*initu.B*pow(rawCoo[i].z - initu.z0, 2);
        math[i][5] = 2*initu.R;
    }

    // Get transMath[6][8]
    for (i = 0; i < 8; i++)
        for (j = 0; j < 6; j++){
            transMath[j][i] = math[i][j];
        }

    // Initialize maty[8][1]
    for (i = 0; i < 8; i++) {
        maty[i][0] = pow(rawCoo[i].x - initu.x0, 2) +
            pow(initu.A, 2)*pow(rawCoo[i].y - initu.y0, 2) +
            pow(initu.B, 2)*pow(rawCoo[i].z - initu.z0, 2) - 
            pow(initu.R, 2);
    }
}

void calcMarginOfError()
{
    int i, j, k;
    double prod1[6][6] = {{0.0}};
    double inv[6][6] = {{0.0}};
    double prod2[6][8] = {{0.0}};
    double temp;

    // Calc prod1[6][6] 
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++){
            for (k = 0, temp = 0; k < 8; k++) {
                temp += transMath[j][k]*math[k][i];
            }
            prod1[j][i] = temp;
        }

    // The inverse of prod1
    calcInv(prod1, 6, inv);

    // Calc prod2[6][8]
    for (i = 0; i < 8; i++)
        for (j = 0; j < 6; j++){
            for (k = 0, temp = 0; k < 6; k++){
                temp += inv[j][k]*transMath[k][i];
            }
            prod2[j][i] = temp;
        }

    // Get result
    for (i = 0; i < 6; i++) {
        for (j = 0, temp = 0; j < 8; j++){
            temp += prod2[i][j]*maty[j][0];
        }
        matDeltau[i][0] = temp;
    }

    // Update MARGIN OF ERROR
    deltau.x0 = matDeltau[0][0];
    deltau.y0 = matDeltau[1][0];
    deltau.z0 = matDeltau[2][0];
    deltau.A = matDeltau[3][0];
    deltau.B = matDeltau[4][0];
    deltau.R = matDeltau[5][0];
}


/*
 * Description: P(Deltau) = pow(h * Deltau - y, 2)
 *
 *
 */
void powerOfDifference()
{
    int i, j;
    double temp;
    double mleft[8][1], mright[8][1];
    double mdiff[8][1];

    for (i = 0; i < 8; i++) {
        for (j = 0, temp = 0; j < 6; j++) {
            temp += math[i][j]*matDeltau[j][0];
        }
        mleft[i][0] = temp;
        mright[i][0] = maty[i][0];
    }

    for (i = 0, temp = 0; i < 8; i++) {
        mdiff[i][0] = mleft[i][0] - mright[i][0];
        mdiff[i][0] = pow(mdiff[i][0], 2);
        temp += mdiff[i][0];
    }

    miniSum = temp;    
}

fitresult owo[10];     // To print them for watching the fitting process

/*
 * Description: Update intial value
 *
 *
 */

int main()
{
   int i;

   initMatd(); 
   calcInitialValue(); 

   for (i = 0; i < 10; i++) {   
       initMath();
       calcMarginOfError();
       powerOfDifference();
       owo[i].mini = miniSum;
       owo[i].u = initu;
       owo[i].delt = deltau;
   }
   
    
    
    return 0;
}

/*
 * Calculate the determinant of the matrix
 *
 */
double calcDet(double a[][N], int n)
{
	double ans = 0;
	double temp[N][N]={{0.0}};
	double t;
	int i,j,k;

	if (n == 1) {
		return a[0][0];
	}

	for (i = 0; i < n; i++) {
		for(j = 0; j < n-1; j++) {
			for(k = 0; k < n-1; k++) {                             
				temp[j][k] = a[j+1][(k>=i)?k+1:k];    
			}
		}
		t = calcDet(temp, n-1); 
		if (i%2 == 0) {     // + - + - ...                   
			ans += a[0][i]*t;
		} else {
			ans -= a[0][i]*t;
		}
	}      
	return ans;    
}

/*
 * Calculate the adjoint of the matrix
 *
 */
void  calcAdj(double ori[][N], int n, double ans[][N])
{
	int i,j,k,t;
	double temp[N][N]={{0.0}};
	
	if (n == 1) {
		ans[0][0] = 1;
		return;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n-1; k++) {
				for (t = 0; t < n-1; t++) {
					temp[k][t] = ori[k>=i?k+1:k][t>=j?t+1:t];
				}
			}
			ans[j][i] = calcDet(temp,n-1);
			if ((i+j)%2 == 1) {
				ans[j][i] = - ans[j][i];
			}
		}
	}
}

/*
 * Calculate the inverse of the matrix
 *
 */
bool calcInv(double src[][N], int n, double des[][N])
{
	double det = calcDet(src,n);
	double t[N][N]={{0.0}};
	int i, j;
	
	if(det == 0) {
		return false;		// The matrix is singular, haven't multiplicative inverse.
	} else {
		calcAdj(src,n,t);
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				des[i][j]=t[i][j]/det;
			}
		}
	}
	return true;
}
