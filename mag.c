/*
 *  D * w = c
 *  w = (D(t) * D)(-i) * D(t) c
 *  reference from http://www.cnblogs.com/tianya2543/p/3999118.html
 */

int32_t Dmat[8][6];
int32_t transposeDmat[6][8];
int32_t product1[6][6];
int32_t inverseProduct1[6][6];
int32_t product2[6][8];
int32_t c[8][1];
int32_t w[6][1];



void OriginalValueCalculate()
{
	uint8_t i, j, k;
	int32_t temp;
	int32_t tempMartrix[6][12];

	for (i = 0; i < 8; i++)
		for (j = 0; j < 6; j++) {
			switch(j){
				case 0:
					Dmat[i][j] = -pow(pointForOriginalCalculate[i].y, 2);
					break;
				case 1:
					Dmat[i][j] = -pow(pointForOriginalCalculate[i].z, 2);
					break;
				case 2:
					Dmat[i][j] = 2*pointForOriginalCalculate.x;
					break;
				case 3:
					Dmat[i][j] = 2*pointForOriginalCalculate.y;
					break;
				case 4:
					Dmat[i][j] = 2*pointForOriginalCalculate.z;
					break;	
				case 5:
					Dmat[i][j] = 1;
					break;					
			}
			transposeDmat[j][i] = Dmat[i][j];
		}

	// Calculate product between D and transposed D
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++){
			for (k = 0, temp = 0; k < 8; k++) {
				temp += transposeDmat[j][k] * Dmat[k][i];
			}
			product1[j][i] = temp;
		}
	}

	// Calculate inverse matrix of the product
	calcInv(product1, 6, inverseProduct1);

	// Calculate product2
	for (i = 0; i < 8; i++)
		for (j = 0; j < 6; j++) {
			for (k = 0; k < 6; k++) {
				temp += inverseProduct1[j][k]*transposeDmat[k][i];
			}
			product2[j][i] = temp;
		}

	// Get result
	for (i = 0; i < 1; i++)
		for (j = 0; j < 6; j++) {
			for (k = 0; k < 8; k++) {
				temp += product2[j][k]*c[k][i];
			}
			w[j][i] = temp;
		}
		

}

#define N	8//6


/*
 * Calculate the determinant of the matrix
 *
 */
double calcDet(double a[][N], int n)
{
	double ans = 0;
	double temp[N][N]={0.0};
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
		if (i%2 == 0) {                   
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
	double temp[N][N];
	
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
	double flag = calcDet(src,n);
	double t[N][N];
	int i, j;
	
	if(flag == 0) {
		return false;		// The matrix is singular, haven't multiplicative inverse.
	} else {
		calcAdj(src,n,t);
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				des[i][j]=t[i][j]/flag;
			}
		}
	}
	return true;
}
