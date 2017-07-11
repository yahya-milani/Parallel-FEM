//============================================================================
// Name        : Code.cpp
// Author      : Yahya Milani
// Version     :
// Copyright   : Use as much as you like with or without author's name for noncommercial cases ONLY
// Description : Hello World in C++, Ansi-style
//============================================================================

//#include <stdafx.h>
#include <stdio.h>
//#include <tchar.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

//Globala Vars
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

// Global Arrays
const double C[6][6] = { { la + 2 * m, la, la, 0, 0, 0 }, { la, la + 2 * m, la,
		0, 0, 0 }, { la, la, la + 2 * m, 0, 0, 0 }, { 0, 0, 0, m, 0, 0 }, { 0,
		0, 0, 0, m, 0 }, { 0, 0, 0, 0, 0, m } };
static double B_ele[ele_num][6][12];
static double V[ele_num];
static double K[length][length];
static double M[length][length];
//
static double K_sp[2619]; //[51484];
int RK_sp[DOF * (node_num - bound_num) + 1];
static int CK_sp[2619]; // [51484];
static double F[DOF * (node_num - bound_num)];

//
void Multi(double **x, int Ix, int Jx, double **y, int Iy, int Jy);
void SpVec(double *A_sp, int *RA_sp, int *CA_sp, double *b_vec, int length_vec);
double VecVec(double *a, double *b, int n);
void Solve();
//double *Solve(double *A_sp, int *RA_sp, int *CA_sp, double *B, int length);
//double **Material()

int main() {

	//Read Nodes***********************************************/
	ifstream nodes;
	nodes.open("nodes.txt");
	int i = 0;
	double node[node_num][3];
	double data1, data2, data3, data0;
	while (nodes >> data0 >> data1 >> data2 >> data3) {
		node[i][0] = data1;
		node[i][1] = data2;
		node[i][2] = data3;
		i += 1;
	}
	int node_size = i;
	nodes.close();
	//Create Nodes DOF ***********************************************/
	static int ID_node[node_num][DOF];
	int counter = 1;
	bool bound;
	for (int i = 0; i < node_num; i++) {
		bound = 1;
		for (int j = 0; j < bound_num; j++) {
			if (i == boundary[j]) {
				bound = 0;
			}
		}
		if (bound == 1) {
			for (int k = 0; k < DOF; k++)
				ID_node[i][k] = counter + k;
			counter += DOF;
		} else {
			for (int k = 0; k < DOF; k++)
				ID_node[i][k] = 0;
		}
	}

	//Read Elements***********************************************/
	ifstream element;
	element.open("elements.txt");
	i = 0;
	int ele[ele_num][4];
	int data4, data5, data6, data7, data00;
	while (element >> data00 >> data4 >> data5 >> data6 >> data7) {
		ele[i][0] = data4 - 1;
		ele[i][1] = data5 - 1;
		ele[i][2] = data6 - 1;
		ele[i][3] = data7 - 1;
		i += 1;
	}
	int ele_size = i;
	element.close();

	//Create Elemental DOF ***********************************************/
	static int ID_ele[ele_num][4 * DOF];
	for (int i = 0; i < ele_num; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < DOF; k++) {
				ID_ele[i][j * DOF + k] = ID_node[ele[i][j]][k];
			}
		}

	}

	//************************************************/
	double x[4], y[4], z[4];
	static double B[ele_num][4][4];
	//set B matrix and voulme for every element
	for (int i = 0; i < ele_num; i++) {
		for (int j = 0; j < 4; j++) {
			x[j] = node[(ele[i][j])][0];
			y[j] = node[ele[i][j]][1];
			z[j] = node[ele[i][j]][2];
		}

		// Compute Voulme of each element so that none has negative volume and change indexes if there is one
		V[i] = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
				- x[1] * y[2] * z[0] - x[2] * y[0] * z[1] + x[2] * y[1] * z[0]
				+ x[0] * y[1] * z[3] - x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
				+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1] - x[3] * y[1] * z[0]
				- x[0] * y[2] * z[3] + x[0] * y[3] * z[2] + x[2] * y[0] * z[3]
				- x[2] * y[3] * z[0] - x[3] * y[0] * z[2] + x[3] * y[2] * z[0]
				+ x[1] * y[2] * z[3] - x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
				+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[3] * y[2] * z[1])
				/ 6;
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
			V[i] = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
					- x[1] * y[2] * z[0] - x[2] * y[0] * z[1]
					+ x[2] * y[1] * z[0] + x[0] * y[1] * z[3]
					- x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
					+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1]
					- x[3] * y[1] * z[0] - x[0] * y[2] * z[3]
					+ x[0] * y[3] * z[2] + x[2] * y[0] * z[3]
					- x[2] * y[3] * z[0] - x[3] * y[0] * z[2]
					+ x[3] * y[2] * z[0] + x[1] * y[2] * z[3]
					- x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
					+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2]
					- x[3] * y[2] * z[1]) / 6;

			// Change ID after changing elemnt indexes
			for (int i2 = 2 * DOF; i2 < 3 * DOF; i2++) {
				int temp2 = ID_ele[i][i2];
				ID_ele[i][i2] = ID_ele[i][i2 + DOF];
				ID_ele[i][i2 + DOF] = temp2;
			}
		}
		B[i][0][0] = (x[1] * y[2] * z[3] - x[1] * y[3] * z[2]
				- x[2] * y[1] * z[3] + x[2] * y[3] * z[1] + x[3] * y[1] * z[2]
				- x[3] * y[2] * z[1]) / (6 * V[i]);
		B[i][0][1] = (y[2] * z[1] - y[1] * z[2] + y[1] * z[3] - y[3] * z[1]
				- y[2] * z[3] + y[3] * z[2]) / (6 * V[i]);
		B[i][0][2] = (x[1] * z[2] - x[2] * z[1] - x[1] * z[3] + x[3] * z[1]
				+ x[2] * z[3] - x[3] * z[2]) / (6 * V[i]);
		B[i][0][3] = (x[2] * y[1] - x[1] * y[2] + x[1] * y[3] - x[3] * y[1]
				- x[2] * y[3] + x[3] * y[2]) / (6 * V[i]);
		B[i][1][0] = (x[0] * y[3] * z[2] - x[0] * y[2] * z[3]
				+ x[2] * y[0] * z[3] - x[2] * y[3] * z[0] - x[3] * y[0] * z[2]
				+ x[3] * y[2] * z[0]) / (6 * V[i]);
		B[i][1][1] = (y[0] * z[2] - y[2] * z[0] - y[0] * z[3] + y[3] * z[0]
				+ y[2] * z[3] - y[3] * z[2]) / (6 * V[i]);
		B[i][1][2] = (x[2] * z[0] - x[0] * z[2] + x[0] * z[3] - x[3] * z[0]
				- x[2] * z[3] + x[3] * z[2]) / (6 * V[i]);
		B[i][1][3] = (x[0] * y[2] - x[2] * y[0] - x[0] * y[3] + x[3] * y[0]
				+ x[2] * y[3] - x[3] * y[2]) / (6 * V[i]);
		B[i][2][0] = (x[0] * y[1] * z[3] - x[0] * y[3] * z[1]
				- x[1] * y[0] * z[3] + x[1] * y[3] * z[0] + x[3] * y[0] * z[1]
				- x[3] * y[1] * z[0]) / (6 * V[i]);
		B[i][2][1] = (y[1] * z[0] - y[0] * z[1] + y[0] * z[3] - y[3] * z[0]
				- y[1] * z[3] + y[3] * z[1]) / (6 * V[i]);
		B[i][2][2] = (x[0] * z[1] - x[1] * z[0] - x[0] * z[3] + x[3] * z[0]
				+ x[1] * z[3] - x[3] * z[1]) / (6 * V[i]);
		B[i][2][3] = (x[1] * y[0] - x[0] * y[1] + x[0] * y[3] - x[3] * y[0]
				- x[1] * y[3] + x[3] * y[1]) / (6 * V[i]);
		B[i][3][0] = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2]
				+ x[1] * y[0] * z[2] - x[1] * y[2] * z[0] - x[2] * y[0] * z[1]
				+ x[2] * y[1] * z[0]) / (6 * V[i]);
		B[i][3][1] = (y[0] * z[1] - y[1] * z[0] - y[0] * z[2] + y[2] * z[0]
				+ y[1] * z[2] - y[2] * z[1]) / (6 * V[i]);
		B[i][3][2] = (x[1] * z[0] - x[0] * z[1] + x[0] * z[2] - x[2] * z[0]
				- x[1] * z[2] + x[2] * z[1]) / (6 * V[i]);
		B[i][3][3] = (x[0] * y[1] - x[1] * y[0] - x[0] * y[2] + x[2] * y[0]
				+ x[1] * y[2] - x[2] * y[1]) / (6 * V[i]);
	}
	//For end of Geometrical issues

	for (int k = 0; k < ele_num; k++) {
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 12; j++) {
				B_ele[k][i][j] = 0;
			}
		}
	}

	for (int k = 0; k < ele_num; k++) {
		B_ele[k][0][0] = B[k][0][1];
		B_ele[k][0][3] = B[k][1][1];
		B_ele[k][0][6] = B[k][2][1];
		B_ele[k][0][9] = B[k][3][1];

		B_ele[k][1][1] = B[k][0][2];
		B_ele[k][1][4] = B[k][1][2];
		B_ele[k][1][7] = B[k][2][2];
		B_ele[k][1][10] = B[k][3][2];

		B_ele[k][2][2] = B[k][0][3];
		B_ele[k][2][5] = B[k][1][3];
		B_ele[k][2][8] = B[k][2][3];
		B_ele[k][2][11] = B[k][3][3];

		B_ele[k][3][0] = B[k][0][2];
		B_ele[k][3][3] = B[k][1][2];
		B_ele[k][3][6] = B[k][2][2];
		B_ele[k][3][9] = B[k][3][2];

		B_ele[k][3][1] = B[k][0][1];
		B_ele[k][3][4] = B[k][1][1];
		B_ele[k][3][7] = B[k][2][1];
		B_ele[k][3][10] = B[k][3][1];

		B_ele[k][4][1] = B[k][0][3];
		B_ele[k][4][4] = B[k][1][3];
		B_ele[k][4][7] = B[k][2][3];
		B_ele[k][4][10] = B[k][3][3];

		B_ele[k][4][2] = B[k][0][2];
		B_ele[k][4][5] = B[k][1][2];
		B_ele[k][4][8] = B[k][2][2];
		B_ele[k][4][11] = B[k][3][2];

		B_ele[k][5][0] = B[k][0][3];
		B_ele[k][5][3] = B[k][1][3];
		B_ele[k][5][6] = B[k][2][3];
		B_ele[k][5][9] = B[k][3][3];
		B_ele[k][5][2] = B[k][0][1];
		B_ele[k][5][5] = B[k][1][1];
		B_ele[k][5][8] = B[k][2][1];
		B_ele[k][5][11] = B[k][3][1];
	}

	static double B_eleT[ele_num][12][6];
	for (int k = 0; k < ele_num; k++) {
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 6; j++) {
				B_eleT[k][i][j] = B_ele[k][j][i];
			}
		}
	}

	for (int i = 0; i < ele_num; i++) {
		double m2[12][12];
		for (int i2 = 0; i2 < 12; i2++) {
			m2[i2][i2] = Ro * V[i] / 4;
		}
		for (int n = 0; n < DOF * 4; n++) {
			for (int n1 = 0; n1 < DOF * 4; n1++) {
				if (ID_ele[i][n] != 0 && ID_ele[i][n1] != 0)
					M[(ID_ele[i][n] - 1)][(ID_ele[i][n1] - 1)] += m2[n][n1];
			}
		}
	}

	for (int n = 0; n < DOF * (node_num - bound_num); n++) {
		for (int n1 = 0; n1 < DOF * (node_num - bound_num); n1++) {
			K[n][n1] = 0;
		}
	}

	for (int i = 0; i < ele_num; i++) {
		double K1[12][6];
		double sum = 0;
		for (int i2 = 0; i2 < 12; i2++) {
			for (int j2 = 0; j2 < 6; j2++) {
				for (int k2 = 0; k2 < 6; k2++) {
					sum += B_eleT[i][i2][k2] * C[k2][j2]; // *(*(*(B_minT + i2) + k2))*(*(*(C + k2) + j2));
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
					sum += K1[i2][k2] * B_ele[i][k2][j2] * V[i]; //(*(*(K1 + i2) + k2))*(*(*(B_min + k2) + j2))*V[i];
				}
				K2[i2][j2] = sum;
				sum = 0;
			}
		}

		for (int n = 0; n < DOF * 4; n++) {
			for (int n1 = 0; n1 < DOF * 4; n1++) {
				if (ID_ele[i][n] != 0 && ID_ele[i][n1] != 0)
					K[(ID_ele[i][n] - 1)][(ID_ele[i][n1] - 1)] += K2[n][n1];
			}
		}
	}

	int counter2 = 0;
	for (int i = 0; i < DOF * (node_num - bound_num); i++) {
		for (int j = 0; j < DOF * (node_num - bound_num); j++) {
			if (K[i][j] != 0) {
				counter2 += 1;
			}
		}
	}

	// how to implement counter2??????????????????? VERY IMPORTANT where should I define sparse matrixes
	double trace = 0;
	for (int i = 0; i < DOF * (node_num - bound_num); i++) {
		trace += K[i][i];
	}

	//Sparsize K ********************************************/
	/*static double K_sp[2619]; //[51484];
	int RK_sp[DOF * (node_num - bound_num) + 1];
	static int CK_sp[2619]; // [51484];*/
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
	//Force **********************************/

	//double F[DOF * (node_num - bound_num)];
	for (int i = 0; i < DOF * (node_num - bound_num); i++) {
		F[i] = 0;
	}
	F[56] = -1000; //	F[686] = -1000;

	//Solve********************************

	static double p[length];
	static double r[length];
	static double U[length];
	static double t[length];
	double rho = 0;
	double rhos;
	double rho0;
	double alpha;
	int solved = 0;
	for (int i = 0; i < length; i++) {
		p[i] = F[i];
		r[i] = F[i];
		U[i] = 0;
	}
	for (int i = 0; i < length; i++) {
		rho += r[i] * r[i];
	}
	cout << "rho    " << rho << "\n";
	for (int j = 0; j < length; j++) {
		if (solved == 0) {
			double sum = 0;
			for (int i = 0; i < length; i++) {
				for (int k = RK_sp[i]; k < RK_sp[i + 1]; k++) {
					sum += (p[CK_sp[k]]) * (K_sp[k]);
				}
				t[i] = sum;
				sum = 0;
			}
			cout << "t   " << (t[0]) << "    " << t[1] << "    " << t[2]
					<< "    " << t[3] << "\n";
			double PT = 0;
			for (int i = 0; i < length; i++) {
				PT += p[i] * t[i];
			}
			cout << "PT   " << PT << "\n";
			alpha = rho / PT;
			cout << "alpha       " << alpha << "\n";
			for (int i = 0; i < length; i++) {
				U[i] += alpha * p[i];
				r[i] -= alpha * t[i];
			}
			cout << "U   " << *(U + 56) << "\n" << "\n";
			//cout << "U   " << *(U) << "   " << *(U + 1) << "   " << *(U + 2) << "   " << *(U + 3) << "\n"<<"\n";
			//cout << "r   " << *(r) << "   " << *(r + 1) << "   " << *(r + 2) << "   " << *(r + 3) << "\n";
			//cout << *(U + 1) << "\n";
			rhos = rho;
			rho = 0;
			for (int i = 0; i < length; i++) {
				rho += r[i] * r[i];
			}
			if ((rho / rhos) < 0.2) {
				solved = 1;
				//cout << "HEEEEEEEEEEY" << "\n";
			}
			cout << "rho    " << rho << "\n";
			for (int i = 0; i < length; i++) {
				p[i] = r[i] + (rho / rhos) * p[i];
			}
		}
		Solve();
	}

	cout << "Bye Bye" << "\n";
	//system("pause");
	return 0;
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
	static double C[1416];//define it cause it's not a good idea to return local variable adress
	double sum = 0;
	for (int i = 0; i < length_vec; i++) {
		for (int j = *(RA_sp + i); j < *(RA_sp + i + 1); j++) {
			sum += *(b_vec + *(CA_sp + j)) * (*(A_sp + j));
		}
		//cout<<sum<<"\n";
		C[i] = sum;
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

//***************/
void Solve() {
	static double p[length];
	static double r[length];
	static double U[length];
	static double t[length];
	double rho = 0;
	double rhos;
	double rho0;
	double alpha;
	int solved = 0;
	for (int i = 0; i < length; i++) {
		p[i] = F[i];
		r[i] = F[i];
		U[i] = 0;
	}
	for (int i = 0; i < length; i++) {
		rho += r[i] * r[i];
	}
	cout << "rho    " << rho << "\n";
	for (int j = 0; j < length; j++) {
		if (solved == 0) {
			double sum = 0;
			for (int i = 0; i < length; i++) {
				for (int k = RK_sp[i]; k < RK_sp[i + 1]; k++) {
					sum += (p[CK_sp[k]]) * (K_sp[k]);
				}
				t[i] = sum;
				sum = 0;
			}
			cout << "t   " << (t[0]) << "    " << t[1] << "    " << t[2]
					<< "    " << t[3] << "\n";
			double PT = 0;
			for (int i = 0; i < length; i++) {
				PT += p[i] * t[i];
			}
			cout << "PT   " << PT << "\n";
			alpha = rho / PT;
			cout << "alpha       " << alpha << "\n";
			for (int i = 0; i < length; i++) {
				U[i] += alpha * p[i];
				r[i] -= alpha * t[i];
			}
			cout << "U   " << *(U + 56) << "\n" << "\n";
			//cout << "U   " << *(U) << "   " << *(U + 1) << "   " << *(U + 2) << "   " << *(U + 3) << "\n"<<"\n";
			//cout << "r   " << *(r) << "   " << *(r + 1) << "   " << *(r + 2) << "   " << *(r + 3) << "\n";
			//cout << *(U + 1) << "\n";
			rhos = rho;
			rho = 0;
			for (int i = 0; i < length; i++) {
				rho += r[i] * r[i];
			}
			if ((rho / rhos) < 0.2) {
				solved = 1;
				cout << "HEEEEEEEEEEY" << "\n";
			}
			cout << "rho    " << rho << "\n";
			for (int i = 0; i < length; i++) {
				p[i] = r[i] + (rho / rhos) * p[i];
			}
		}
	}
	/*const int length = 1416;
	 double p[length];
	 double r[length];
	 static double x[length];
	 double *t;
	 double rho;
	 double rhos;
	 double rho0;
	 double alpha;
	 /*
	 for (int i = 0; i<length; i++){
	 p[i] = B[i];
	 r[i] = B[i];
	 }
	 rho = VecVec(r, r, length);
	 cout << "rho    " << rho << "\n";
	 for (int j = 0; j<length; j++)
	 {
	 //t = SpVec(A_sp, RA_sp, CA_sp, p, 3);
	 cout << "t   "<< *(t+1) << "\n";
	 alpha = rho / VecVec(p, t, length);
	 cout<<"alpha       "<<alpha<<"\n";
	 for (int i = 0; i<length; i++)
	 {
	 x[i] += alpha*p[i];
	 r[i] -= alpha*(*(t + i));
	 }
	 cout << *x << "\n";
	 rhos = rho;
	 rho = VecVec(r, r, length);
	 if ((rho / rhos) < 0.02)
	 break;
	 cout << "rho    " << rho << "\n";
	 for (int i = 0; i<length; i++)
	 p[i] = r[i] + (rho / rhos)*p[i];
	 }
	 return x;
	 */
	//return 0;
}
