//============================================================================
// Name        : Code.cpp
// Author      : Yahya Milani
// Version     :
// Copyright   : Use as much as you like with author's name for noncommercial cases ONLY
// Description : C++, Ansi-style
//============================================================================

//#include <stdafx.h>
#include <stdio.h>
//#include <tchar.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

//Global Vars
const int ele_num = 96;
const int node_num = 35;
const int bound_num = 4;
const int boundary[4] = { 0, 3, 6, 7 }; //{ 2,3,5,6 };
const int DOF = 3;
const int length = DOF * (node_num - bound_num);
const double Ro = 1000;
const double m = 782.7586;
const double la = 7866.7;
const double tF = 100;
const double dt = 0.01;
const double dt_save = 1;
const int sparse_data = 2619;

// Global Arrays

//static double V[ele_num];

//


//
double Volume(double x[4], double y[4], double z[4]);
void Jacobi(double x[4], double y[4], double z[4], double Vol, double J[16]);
void shape_function(double B[ele_num * 6 * 12], double J[ele_num][4][4]);
void transpose(double A[], int row_A, int col_A, double A_T[]);
void Mass(double M[length * length], double V[ele_num],
		int ID_ele[ele_num * 4 * DOF]);
void Stiffness(double B_ele[ele_num * 6 * 12], double B_eleT[ele_num * 12 * 6],
		double C[6][6], int ID_ele[4 * DOF * ele_num], double V[ele_num],
		double K[length][length]);
void Multi(double **x, int Ix, int Jx, double **y, int Iy, int Jy);
void SpVec(double *A_sp, int *RA_sp, int *CA_sp, double *b_vec, int length_vec);
double VecVec(double *a, double *b, int n);
void Solve(double A_sp[],int Col_A_sp[],int Row_A_sp[],double B[],double U[]);
//double *Solve(double *A_sp, int *RA_sp, int *CA_sp, double *B, int length);
//double **Material()

int main() {

	//Read Nodes=================================================/
	ifstream nodes;
	nodes.open("nodes.txt");
	double node[node_num * DOF] = { };
	//int row_node = node_num;
	int col_node = DOF;
	int i = 0;
	double data1, data2, data3, data0;
	while (nodes >> data0 >> data1 >> data2 >> data3) {
		node[i * col_node + 0] = data1;
		node[i * col_node + 1] = data2;
		node[i * col_node + 2] = data3;
		i += 1;
	}
	nodes.close();
	//Create Nodes DOF =================================================/
	int ID_node[node_num * DOF] = { };
	int row_ID_node = node_num;
	int col_ID_node = DOF;
	int counter = 1;
	bool bound;
	for (int i = 0; i < row_ID_node; i++) {
		bound = 1;
		for (int j = 0; j < bound_num; j++) {
			if (i == boundary[j]) {
				bound = 0;
			}
		}
		if (bound == 1) {
			for (int k = 0; k < DOF; k++)
				ID_node[i * col_ID_node + k] = counter + k;
			counter += DOF;
		} else {
			for (int k = 0; k < DOF; k++)
				ID_node[i * col_ID_node + k] = 0;
		}
	}

	//Read Elements=================================================/
	ifstream element;
	element.open("elements.txt");
	i = 0;
	int ele[ele_num * 4] = { };
	//int row_ele = ele_num;
	int col_ele = 4;
	int data4, data5, data6, data7, data00;
	while (element >> data00 >> data4 >> data5 >> data6 >> data7) {
		ele[i * col_ele + 0] = data4 - 1;
		ele[i * col_ele + 1] = data5 - 1;
		ele[i * col_ele + 2] = data6 - 1;
		ele[i * col_ele + 3] = data7 - 1;
		i += 1;
	}
	element.close();

	//Create Elemental DOF =================================================/
	int ID_ele[ele_num * 4 * DOF] = { };
	int row_ID_ele = ele_num;
	int col_ID_ele = 4 * DOF;
	for (int i = 0; i < row_ID_ele; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < DOF; k++) {
				ID_ele[i * col_ID_ele + j * DOF + k] = ID_node[(ele[i * col_ele
						+ j]) * col_ID_node + k];
			}
		}

	}

	//=================================================/
	double x[4], y[4], z[4];
	double V[ele_num] = { };
	double B[ele_num][4][4] = { };
	//set B matrix and voulme for every element
	for (int i = 0; i < ele_num; i++) {
		for (int j = 0; j < 4; j++) {
			x[j] = node[(ele[i * col_ele + j]) * col_node + 0];
			y[j] = node[(ele[i * col_ele + j]) * col_node + 1];
			z[j] = node[(ele[i * col_ele + j]) * col_node + 2];
		}

		// Compute Voulme of each element so that none has negative volume and change indexes if there is one
		V[i] = Volume(x, y, z);

		if (V[i] < 0) {
			double temp = x[2];
			x[2] = x[3];
			x[3] = temp;
			temp = y[2];
			y[2] = y[3];
			y[3] = temp;
			temp = z[2];
			z[2] = z[3];
			z[3] = temp;
			V[i] = Volume(x, y, z);

			// Change ID after changing elemnt indexes
			for (int i2 = 2 * DOF; i2 < 3 * DOF; i2++) {
				int temp2 = ID_ele[i * col_ID_ele + i2];
				ID_ele[i * col_ID_ele + i2] = ID_ele[i * col_ID_ele + i2 + DOF];
				ID_ele[i * col_ID_ele + i2 + DOF] = temp2;
			}
		}
		double J[16];
		Jacobi(x, y, z, V[i], J);
		for (int var = 0; var < 4; ++var) {
			for (int var2 = 0; var2 < 4; ++var2) {
				B[i][var][var2] = J[var * 4 + var2];
			}

		}

	}
	//For end of Geometrical issues========================================================

	//=====================================================================================
	//===========================ALL Shit Happens Here=====================================
	//=====================================================================================
	double B_ele[ele_num * 6 * 12] = { };
	int row_B_ele = 6;
	int col_B_ele = 12;
	shape_function(B_ele, B);

	double B_eleT[ele_num * 12 * 6] = { };
	//int row_B_eleT = 12;
	//int col_B_eleT = 6;
	transpose(B_ele, row_B_ele, col_B_ele, B_eleT);

	double M[length * length] = { };
	Mass(M, V, ID_ele);

	clock_t start;
	start = clock();
	double K[length][length];
	for (int n = 0; n < DOF * (node_num - bound_num); n++) {
		for (int n1 = 0; n1 < DOF * (node_num - bound_num); n1++) {
			K[n][n1] = 0;
		}
	}
	double C[6][6] = { { la + 2 * m, la, la, 0, 0, 0 }, { la, la + 2 * m, la, 0,
			0, 0 }, { la, la, la + 2 * m, 0, 0, 0 }, { 0, 0, 0, m, 0, 0 }, { 0,
			0, 0, 0, m, 0 }, { 0, 0, 0, 0, 0, m } };
	Stiffness(B_ele, B_eleT, C,ID_ele,V, K);

	//=====================================================================================
	//=====================================================================================
	//=====================================================================================

	// how to implement counter2??????????????????? VERY IMPORTANT where should I define sparse matrixes
	double trace = 0;
	for (int i = 0; i < DOF * (node_num - bound_num); i++) {
		trace += K[i][i];
	}

	//Sparsize K ============================================================/
	/*static double K_sp[2619]; //[51484];
	 int RK_sp[DOF * (node_num - bound_num) + 1];
	 static int CK_sp[2619]; // [51484];*/
	double K_sp[sparse_data]={}; //[51484];
	int RK_sp[DOF * (node_num - bound_num) + 1]={};
	int CK_sp[sparse_data]={}; // [51484];
	double F[DOF * (node_num - bound_num)]={};
	RK_sp[0] = 0;
	counter = 0;
	for (int i = 0; i < DOF * (node_num - bound_num); i++) {
		for (int j = 0; j < DOF * (node_num - bound_num); j++) {
			if (K[i][j] != 0) {
				K_sp[counter] = K[i][j];
				CK_sp[counter] = j;
				counter += 1;
			}
			RK_sp[i + 1] = counter;
		}
	}
	cout << counter << "\n";
	//Force =================================

	F[56] = -1000; //	F[686] = -1000;

	//Solve =================================


//	for (int j = 0; j < length; j++) {
	double U[length]={};
		Solve(K_sp,CK_sp,RK_sp,F,U);
//	}

	cout << "Bye Bye" << "\n";
	double duration;
	duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << "time: " << duration << endl;
	return 0;
}



//===========================================================================================//
double Volume(double x[4], double y[4], double z[4]) {
	double V;
	V = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
			- x[1] * y[2] * z[0] - x[2] * y[0] * z[1] + x[2] * y[1] * z[0]
			+ x[0] * y[1] * z[3] - x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
			+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1] - x[3] * y[1] * z[0]
			- x[0] * y[2] * z[3] + x[0] * y[3] * z[2] + x[2] * y[0] * z[3]
			- x[2] * y[3] * z[0] - x[3] * y[0] * z[2] + x[3] * y[2] * z[0]
			+ x[1] * y[2] * z[3] - x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
			+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[3] * y[2] * z[1]) / 6;
	return V;

}

void Jacobi(double x[4], double y[4], double z[4], double Vol, double J[16]) {
	J[0] = (x[1] * y[2] * z[3] - x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
			+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[3] * y[2] * z[1])
			/ (6 * Vol);
	J[1] = (y[2] * z[1] - y[1] * z[2] + y[1] * z[3] - y[3] * z[1] - y[2] * z[3]
			+ y[3] * z[2]) / (6 * Vol);
	J[2] = (x[1] * z[2] - x[2] * z[1] - x[1] * z[3] + x[3] * z[1] + x[2] * z[3]
			- x[3] * z[2]) / (6 * Vol);
	J[3] = (x[2] * y[1] - x[1] * y[2] + x[1] * y[3] - x[3] * y[1] - x[2] * y[3]
			+ x[3] * y[2]) / (6 * Vol);
	J[4] = (x[0] * y[3] * z[2] - x[0] * y[2] * z[3] + x[2] * y[0] * z[3]
			- x[2] * y[3] * z[0] - x[3] * y[0] * z[2] + x[3] * y[2] * z[0])
			/ (6 * Vol);
	J[5] = (y[0] * z[2] - y[2] * z[0] - y[0] * z[3] + y[3] * z[0] + y[2] * z[3]
			- y[3] * z[2]) / (6 * Vol);
	J[6] = (x[2] * z[0] - x[0] * z[2] + x[0] * z[3] - x[3] * z[0] - x[2] * z[3]
			+ x[3] * z[2]) / (6 * Vol);
	J[7] = (x[0] * y[2] - x[2] * y[0] - x[0] * y[3] + x[3] * y[0] + x[2] * y[3]
			- x[3] * y[2]) / (6 * Vol);
	J[8] = (x[0] * y[1] * z[3] - x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
			+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1] - x[3] * y[1] * z[0])
			/ (6 * Vol);
	J[9] = (y[1] * z[0] - y[0] * z[1] + y[0] * z[3] - y[3] * z[0] - y[1] * z[3]
			+ y[3] * z[1]) / (6 * Vol);
	J[10] = (x[0] * z[1] - x[1] * z[0] - x[0] * z[3] + x[3] * z[0] + x[1] * z[3]
			- x[3] * z[1]) / (6 * Vol);
	J[11] = (x[1] * y[0] - x[0] * y[1] + x[0] * y[3] - x[3] * y[0] - x[1] * y[3]
			+ x[3] * y[1]) / (6 * Vol);
	J[12] = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
			- x[1] * y[2] * z[0] - x[2] * y[0] * z[1] + x[2] * y[1] * z[0])
			/ (6 * Vol);
	J[13] = (y[0] * z[1] - y[1] * z[0] - y[0] * z[2] + y[2] * z[0] + y[1] * z[2]
			- y[2] * z[1]) / (6 * Vol);
	J[14] = (x[1] * z[0] - x[0] * z[1] + x[0] * z[2] - x[2] * z[0] - x[1] * z[2]
			+ x[2] * z[1]) / (6 * Vol);
	J[15] = (x[0] * y[1] - x[1] * y[0] - x[0] * y[2] + x[2] * y[0] + x[1] * y[2]
			- x[2] * y[1]) / (6 * Vol);
}

void shape_function(double B[ele_num * 6 * 12], double J[ele_num][4][4]) {
	int row_B = 6;
	int col_B = 12;
	for (int k = 0; k < ele_num; k++) {
		//B[col_B*row_B*k+col_B*i+j] [k][i][j]
		B[col_B * row_B * k + col_B * 0 + 0] = J[k][0][1];
		B[col_B * row_B * k + col_B * 0 + 3] = J[k][1][1];
		B[col_B * row_B * k + col_B * 0 + 6] = J[k][2][1];
		B[col_B * row_B * k + col_B * 0 + 9] = J[k][3][1];

		B[col_B * row_B * k + col_B * 1 + 1] = J[k][0][2];
		B[col_B * row_B * k + col_B * 1 + 4] = J[k][1][2];
		B[col_B * row_B * k + col_B * 1 + 7] = J[k][2][2];
		B[col_B * row_B * k + col_B * 1 + 10] = J[k][3][2];

		B[col_B * row_B * k + col_B * 2 + 2] = J[k][0][3];
		B[col_B * row_B * k + col_B * 2 + 5] = J[k][1][3];
		B[col_B * row_B * k + col_B * 2 + 8] = J[k][2][3];
		B[col_B * row_B * k + col_B * 2 + 11] = J[k][3][3];

		B[col_B * row_B * k + col_B * 3 + 0] = J[k][0][2];
		B[col_B * row_B * k + col_B * 3 + 3] = J[k][1][2];
		B[col_B * row_B * k + col_B * 3 + 6] = J[k][2][2];
		B[col_B * row_B * k + col_B * 3 + 9] = J[k][3][2];

		B[col_B * row_B * k + col_B * 3 + 1] = J[k][0][1];
		B[col_B * row_B * k + col_B * 3 + 4] = J[k][1][1];
		B[col_B * row_B * k + col_B * 3 + 7] = J[k][2][1];
		B[col_B * row_B * k + col_B * 3 + 10] = J[k][3][1];

		B[col_B * row_B * k + col_B * 4 + 1] = J[k][0][3];
		B[col_B * row_B * k + col_B * 4 + 4] = J[k][1][3];
		B[col_B * row_B * k + col_B * 4 + 7] = J[k][2][3];
		B[col_B * row_B * k + col_B * 4 + 10] = J[k][3][3];

		B[col_B * row_B * k + col_B * 4 + 2] = J[k][0][2];
		B[col_B * row_B * k + col_B * 4 + 5] = J[k][1][2];
		B[col_B * row_B * k + col_B * 4 + 8] = J[k][2][2];
		B[col_B * row_B * k + col_B * 4 + 11] = J[k][3][2];

		B[col_B * row_B * k + col_B * 5 + 0] = J[k][0][3];
		B[col_B * row_B * k + col_B * 5 + 3] = J[k][1][3];
		B[col_B * row_B * k + col_B * 5 + 6] = J[k][2][3];
		B[col_B * row_B * k + col_B * 5 + 9] = J[k][3][3];

		B[col_B * row_B * k + col_B * 5 + 2] = J[k][0][1];
		B[col_B * row_B * k + col_B * 5 + 5] = J[k][1][1];
		B[col_B * row_B * k + col_B * 5 + 8] = J[k][2][1];
		B[col_B * row_B * k + col_B * 5 + 11] = J[k][3][1];
	}
}

void transpose(double A[], int row_A, int col_A, double A_T[]) {
	for (int k = 0; k < ele_num; k++) {
		for (int i = 0; i < col_A; i++) {
			for (int j = 0; j < row_A; j++) {
				A_T[row_A * col_A * k + row_A * i + j] = A[row_A * col_A * k
						+ col_A * j + i];
			}
		}
	}
}

void Mass(double M[length * length], double V[ele_num],
		int ID_ele[ele_num * 4 * DOF]) {
	int col_ID_ele = 4 * DOF;
	for (int i = 0; i < ele_num; i++) {
		double m2[12 * 12] = { };
		for (int j = 0; j < 12; j++) {
			m2[j * 12 + j] = Ro * V[i] / 4;
		}
		for (int j = 0; j < DOF * 4; j++) {
			for (int k = 0; k < DOF * 4; k++) {
				if (ID_ele[i * col_ID_ele + j] != 0
						&& ID_ele[i * col_ID_ele + k] != 0)
					M[(ID_ele[i * col_ID_ele + j] - 1) * length
							+ (ID_ele[i * col_ID_ele + k] - 1)] +=
							m2[j * 12 + k];
			}
		}
	}
}

void Stiffness(double B_ele[ele_num * 6 * 12], double B_eleT[ele_num * 12 * 6],
		double C[6][6], int ID_ele[4 * DOF * ele_num], double V[ele_num],
		double K[length][length]) {
	int row_B_ele = 6;
	int col_B_ele = 12;
	int row_B_eleT = 12;
	int col_B_eleT = 6;
	int col_ID_ele = 4 * DOF;

	for (int i = 0; i < ele_num; i++) {
		double K1[12][6];
		double sum = 0;
		for (int i2 = 0; i2 < 12; i2++) {
			for (int j2 = 0; j2 < 6; j2++) {
				for (int k2 = 0; k2 < 6; k2++) {
					sum += B_eleT[col_B_eleT * row_B_eleT * i + col_B_eleT * i2
							+ k2] * C[k2][j2]; // *(*(*(B_minT + i2) + k2))*(*(*(C + k2) + j2));
				}
				K1[i2][j2] = sum;
				sum = 0;
			}
		}
		double K2[12][12];
		sum = 0;
		for (int i2 = 0; i2 < 12; i2++) {
			for (int j2 = 0; j2 < 12; j2++) {
				for (int k2 = 0; k2 < 6; k2++) {
					sum += K1[i2][k2]
							* B_ele[col_B_ele * row_B_ele * i + col_B_ele * k2
									+ j2] * V[i]; //(*(*(K1 + i2) + k2))*(*(*(B_min + k2) + j2))*V[i];
				}
				K2[i2][j2] = sum;
				sum = 0;
			}
		}

		for (int n = 0; n < DOF * 4; n++) {
			for (int n1 = 0; n1 < DOF * 4; n1++) {
				if (ID_ele[i * col_ID_ele + n] != 0
						&& ID_ele[i * col_ID_ele + n1] != 0)
					K[(ID_ele[i * col_ID_ele + n] - 1)][(ID_ele[i * col_ID_ele
							+ n1] - 1)] += K2[n][n1];
			}
		}
	}

	int counter2 = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			if (K[i][j] != 0) {
				counter2 += 1;
			}
		}
	}
}

void Multi2(double **x, int Ix, int Jx, double **y, int Iy, int Jy,
		double **z) {
	//Rerurns x*y & Jx=Iy
	double sum = 0;

	for (int i = 0; i < Ix; i++) {
		for (int j = 0; j < Jy; j++) {
			for (int k = 0; k < Iy; k++) {
				sum += (*(*(x + i) + k)) * (*(*(y + k) + j));
			}
			// cout<<sum<<"\n";
			*(*(z + i) + j) = sum;
			sum = 0;
		}
	}
}

void Multi(double **x, int Ix, int Jx, double **y, int Iy, int Jy) {
	//Rerurns x*y & Jx=Iy
	double sum = 0;
	double **R;
	R = new double *[Ix];
	for (int i = 0; i < Ix; i++)
		R[i] = new double[Jy];
	for (int i = 0; i < Ix; i++) {
		for (int j = 0; j < Jy; j++) {
			for (int k = 0; k < Iy; k++) {
				sum += (*(*(x + i) + k)) * (*(*(y + k) + j));
			}
			// cout<<sum<<"\n";
			*(*(R + i) + j) = sum;
			sum = 0;
		}
	}
}

void SpVec(double *A_sp, int *RA_sp, int *CA_sp, double *b_vec,
		int length_vec) {
	//returns A*b
	//static double C[1416]; //define it cause it's not a good idea to return local variable adress
	double sum = 0;
	for (int i = 0; i < length_vec; i++) {
		for (int j = *(RA_sp + i); j < *(RA_sp + i + 1); j++) {
			sum += *(b_vec + *(CA_sp + j)) * (*(A_sp + j));
		}
		//cout<<sum<<"\n";
	//	C[i] = sum;
		sum = 0;
	}
	//return C;
}

double VecVec(double *a, double *b, int n) {
	static double L;
	L = 0;
	for (int i = 0; i < n; i++)
		L += *(a + i) * (*(b + i));
	return L;
}

//====================================================================
void Solve(double A_sp[],int Col_A_sp[],int Row_A_sp[],double B[],double U[]) {
	double p[length];
	double r[length];
//	double U[length];
	double t[length];
	double rho = 0;
	double rhos;
	double alpha;
	double error=0.2;
	bool solved = 0;
	for (int i = 0; i < length; i++) {
		p[i] = B[i];
		r[i] = B[i];
//		U[i] = 0;
	}
	for (int i = 0; i < length; i++) {
		rho += r[i] * r[i];
	}
	cout << "rho    " << rho << "\n";
	for (int j = 0; j < length; j++) {
		if (solved == 0) {
			double sum = 0;
			for (int i = 0; i < length; i++) {
				for (int k = Row_A_sp[i]; k < Row_A_sp[i + 1]; k++) {
					sum += (p[Col_A_sp[k]]) * (A_sp[k]);
				}
				t[i] = sum;
				sum = 0;
			}

			double PT = 0;
			for (int i = 0; i < length; i++) {
				PT += p[i] * t[i];
			}
			alpha = rho / PT;
			for (int i = 0; i < length; i++) {
				U[i] += alpha * p[i];
				r[i] -= alpha * t[i];
			}

			rhos = rho;
			rho = 0;
			for (int i = 0; i < length; i++) {
				rho += r[i] * r[i];
			}
			if ((rho / rhos) < error) {
				solved = 1;
				cout<<endl << "Solved in "<< j<<" Steps!" << "\n";
			}
			for (int i = 0; i < length; i++) {
				p[i] = r[i] + (rho / rhos) * p[i];
			}

		}
	}
	cout << "U(56)= " << *(U + 56) << "\n" << "\n";
}

