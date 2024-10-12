#pragma once
#include <vector>
#include <iostream>

using namespace std;
class Iter_method
{
private:
	int n;
	vector<vector<double>> A;
	vector<double> B;
	vector<double> X, new_X;

public:
	Iter_method() {};
	Iter_method(vector<vector<double>>& A_,
		vector<double>& B_) {
		A = A_; B = B_; X = vector<double>(B_.size(), 1); new_X = vector<double>(B_.size(), 1);
	};
	Iter_method(int n_);
	int Yakoby(double &d);
	int Zeidel(double& d);
	int UpRelax(double& d, double o);
	vector<double> ReturnAnswer() { return X; };
	double tochnost();
	friend ostream& operator<<(ostream& out, Iter_method& m);
};

