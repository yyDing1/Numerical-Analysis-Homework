#include<iostream>
#include<stdio.h>
#include <math.h>
using namespace std;

const int N = 5e2 + 5;
double mat[N][N];

struct Matrix {
	double mat[N][N];
	int n;

	void read(void) {
		scanf("%d", &n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n + 1; j++) {
				scanf("%lf", &mat[i][j]);
			}
		}
	}
	void printMat(double mat[N][N], int n, int m) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				printf("%f ", mat[i][j]);
			}
			printf("\n");
		}
	}
	void printVec(double vec[], int n) {
		for (int i = 0; i < n; i++) {
			printf("%f ", vec[i]);
		}
		printf("\n");
	}

	// pivot_choice_method must be in ["none", "column", "all"]
	void GaussianElimination(bool GaussionJordan=true, string pivot_choice_method="none") {
		double ans[N];
		int row, col, pos[N];
		for (int i = 0; i < n; i++) {
			pos[i] = i;
		}
		for (row = 0, col = 0; row < n && col < n; row++, col++) {
			printf("Stage %d\n", row);
			if (pivot_choice_method == "column") {
				int maxrow = row;
				for (int i = row; i < n; i++) {
					if (fabs(mat[i][col]) > fabs(mat[maxrow][col])) {
						maxrow = i;
					}
				}
				for (int j = 0; j < n + 1; j++){
					swap(mat[row][j], mat[maxrow][j]);
				}
			}
			else if (pivot_choice_method == "all") {
				int maxrow = row, maxcol = col;
				for (int i = row; i < n; i++) {
					for (int j = col; j < n; j++) {
						if (fabs(mat[i][j]) > fabs(mat[maxrow][maxcol])) {
							maxrow = i; maxcol = j;
						}
					}
				}
				for (int j = 0; j < n + 1; j++) {
					swap(mat[row][j], mat[maxrow][j]);
				}
				for (int i = 0; i < n; i++) {
					swap(mat[i][col], mat[i][maxcol]);
				}
				swap(pos[col], pos[maxcol]);
			}
			if (mat[row][col] == 0) {
				printf("Fail!\n");
				return;
			}
			double div = mat[row][col];
			for (int j = col; j < n + 1; j++) {
				mat[row][j] /= div;
			}
			for (int i = GaussionJordan? 0: row + 1; i < n; i++) {
				if (i == row) continue;
				double temp = mat[i][col];
				for (int j = col; j < n + 1; j++) {
					mat[i][j] -= mat[row][j] * temp;
				}
				mat[i][col] = 0;
			}
			printMat(mat, n, n + 1);
		}
		for (int i = n - 1; i >= 0; i--) {
			ans[i] = mat[i][n];
			for (int j = i + 1; j < n; j++) {
				ans[i] -= ans[j] * mat[i][j];
			}
		}
		printf("Result: \n");
		for (int i = 0; i < n; i++) {
			printf("ans%d = ans'%d = %f\n", pos[i], i, ans[i]);
		}
	}
	void CroutSplit(void) {
		double l[N][N], u[N][N], x[N], y[N];
		for (int i = 0; i < n; i++) {
			double sum;
			for (int j = 0; j <= i; j++) {
				sum = 0;
				for (int k = 0; k < j; k++) {
					sum += l[i][k] * u[k][j];
				}
				l[i][j] = mat[i][j] - sum;
			}
			for (int j = i + 1; j < n; j++) {
				sum = 0;
				for (int k = 0; k < i; k++) {
					sum += l[i][k] * u[k][j];
				}
				u[i][j] = (mat[i][j] - sum) / l[i][i];
			}
			u[i][i] = 1;
		}
		printf("L = \n"); printMat(l, n, n);
		printf("U = \n"); printMat(u, n, n);
		for (int k = 0; k < n; k++) {
			double sum = 0;
			for (int i = 0; i < k; i++) {
				sum += l[k][i] * y[i];
			}
			y[k] = (mat[k][n] - sum) / l[k][k];
		}
		printf("y = \n"); printVec(y, n);
		for (int k = n - 1; k >= 0; k--) {
			double sum = 0;
			for (int i = k + 1; i < n; i++) {
				sum += u[k][i] * x[i];
			}
			x[k] = (y[k] - sum) / u[k][k];
		}
		printf("x = \n"); printVec(x, n);
	}

	void JacobiMethod(double x[], int iter) {
		double newx[N];
		for (int iter_num = 1; iter_num <= iter; iter_num++) {
			for (int i = 0; i < n; i++) {
				double sum = 0;
				for (int j = 0; j < n; j++) {
					if (j == i) continue;
					sum += mat[i][j] * x[j];
				}
				newx[i] = (mat[i][n] - sum) / mat[i][i];
			}
			printf("Iter %d: \n", iter_num);
			for (int i = 0; i < n; i++) {
				x[i] = newx[i];
			}
			printVec(x, n);
		}
	}

	void GaussSeidelMethod(double x[], int iter) {
		double newx[N];
		for (int iter_num = 1; iter_num <= iter; iter_num++) {
			for (int i = 0; i < n; i++) {
				double sum = 0;
				for (int j = 0; j < i; j++) {
					sum += mat[i][j] * newx[j];
				}
				for (int j = i + 1; j < n; j++) {
					sum += mat[i][j] * x[j];
				}
				newx[i] = (mat[i][n] - sum) / mat[i][i];
			}
			printf("Iter %d: \n", iter_num);
			for (int i = 0; i < n; i++) {
				x[i] = newx[i];
			}
			printVec(x, n);
		}
	}
}now;

int main(void) {
	now.read();
	now.GaussianElimination(true, "all");
	// now.CroutSplit();
	// double x[N] = {0, 0, 0};
	// now.JacobiMethod(x, 20);
	// now.GaussSeidelMethod(x, 20);
	return 0;
}