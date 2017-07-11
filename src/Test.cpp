/*
 * Test.cpp
 *
 *  Created on: Jul 1, 2017
 *      Author: yahya
 */
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

//typedef vector<float> matrix_f;

void printer(int X[8]);
int main() {
	//cout<<"Hello world!"<<endl;
	//matrix_f A;
	int B[8];
	for (int i = 0; i < 6; i++) {
		cout << B[i] << "  ";
	}
	cout << endl;
	printer(B);
	for (int i = 0; i < 6; i++) {
		cout << B[i] << "  ";
	}
	return 0;
}

void printer(int X[8]) {
	X[5] = 6;
	for (int i = 0; i < 6; i++) {
		cout << X[i] << "  ";
	}
	cout << endl;
}

/*
 void mult_1 (const matrix_f &  matrixOne, const matrix_f & matrixTwo, matrix_f & result) {
 const int matrixSize = (int)result.size();
 //#pragma omp parallel for simd
 for (int rowResult = 0; rowResult < matrixSize; ++rowResult) {
 for (int colResult = 0; colResult < matrixSize; ++colResult) {
 for (int k = 0; k < matrixSize; ++k) {
 result[rowResult][colResult] += matrixOne[rowResult][k] * matrixTwo[k][colResult];
 }
 }
 }
 }


 void mult (const std::vector<float> &  matrixOne, const std::vector<float> & matrixTwo, std::vector<float> & result, int size) {
 for (int row = 0; row < size; ++row) {
 for (int col = 0; col < size; ++col) {
 for (int k = 0; k <size; ++k) {
 result[(size*row)+col] += matrixOne[(size*row)+k] * matrixTwo[(size*k)+col];
 }
 }
 }
 }
 */
