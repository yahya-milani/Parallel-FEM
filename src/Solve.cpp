/*
 * Solve.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: yahya
 */

void Solve(double A_sp[],double Col_A_sp[],double Row_A_sp[],double B[],double U[]) {
	double p[length];
	double r[length];
//	double U[length];
	double t[length];
	double rho = 0;
	double rhos;
	double alpha;
	double error=2;
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
//			cout << "t   " << (t[0]) << "    " << t[1] << "    " << t[2]
//					<< "    " << t[3] << "\n";
			double PT = 0;
			for (int i = 0; i < length; i++) {
				PT += p[i] * t[i];
			}
//			cout << "PT   " << PT << "\n";
			alpha = rho / PT;
//			cout << "alpha       " << alpha << "\n";
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
			if ((rho / rhos) < error) {
				solved = 1;
				cout << "HEEEEEEEEEEY YOU" << "\n";
			}
			cout << "rho    " << rho << "\n";
			for (int i = 0; i < length; i++) {
				p[i] = r[i] + (rho / rhos) * p[i];
			}
		}
	}
}


