#include "Iter_method.h"
#include <time.h>
#include <random>
#include <cmath>

double Iter_method::tochnost()
{
	double n_f = 0;
	for (int i = 0; i < B.size(); ++i)
		if (abs(B[i]) > n_f)
			n_f = abs(B[i]);
	double n_fAx = 0;
	for (int i = 0; i < B.size(); ++i)
	{
		double Ax = 0;
		for (int j = 0; j < B.size(); ++j)
			Ax += A[i][j] * new_X[j];

		if (abs(B[i] - Ax) > n_fAx)
			n_fAx = abs(B[i] - Ax);
	}
	return n_fAx / n_f;
}


int Iter_method::Yakoby(double &d)
{
	X = vector<double>(B.size(), 1);
	new_X = vector<double>(B.size(), 1);
	int c = 0;
	int maxItr = 100000;
	while (tochnost() > d && c <= maxItr)
	{
		//cout << "{" << c << "," << tochnost() << "}, ";
		X = new_X;
		for (int i = 0; i < X.size(); ++i)
		{
			new_X[i] = 0;
			for (int j = 0; j < i; ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * X[j];
			}
			for (int j = i + 1; j < X.size(); ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * X[j];
			}
			new_X[i] += B[i] / A[i][i];
		}
		++c;
	}
	return c;
}

int Iter_method::Zeidel(double& d)
{
	int c = 0;
	X = vector<double>(B.size(), 1);
	new_X = vector<double>(B.size(), 1);

	int maxItr = 100000;
	while (tochnost() > d && c <= maxItr)
	{
		//cout << "{" << c << "," << tochnost() << "}, ";
		X = new_X;
		for (int i = 0; i < X.size(); ++i)
		{
			new_X[i] = 0;
			for (int j = 0; j < i; ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * new_X[j];
			}
			for (int j = i + 1; j < X.size(); ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * X[j];
			}
			new_X[i] += B[i] / A[i][i];
		}
		++c;
	}
	return c;
}

int Iter_method::UpRelax(double& d, double o)
{
	int c = 0;
	X = vector<double>(B.size(), 1);
	new_X = vector<double>(B.size(), 1);

	int maxItr = 100000;
	while (tochnost() > d && c <= maxItr)
	{
		//cout << "{" << c << "," << tochnost() << "}, ";
		X = new_X;
		for (int i = 0; i < X.size(); ++i)
		{
			new_X[i] = 0;
			for (int j = 0; j < i; ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * new_X[j];
			}
			for (int j = i + 1; j < X.size(); ++j)
			{
				new_X[i] -= A[i][j] / A[i][i] * X[j];
			}
			new_X[i] += B[i] / A[i][i];

			new_X[i] *= o;
			new_X[i] += (1 - o) * X[i];
		}
		++c;
	}
	return c;
}

Iter_method::Iter_method(int n_)
{
	n = n_;
	srand(time(NULL));
	A = vector<vector<double>>(n_, vector<double>(n_));
	B = vector<double>(n_);
	for (int i = 0; i < n_; ++i)
		for (int j = 0; j < n_; ++j)
			A[i][j] = (rand() % 10);
	for (int i = 0; i < n_; ++i)
		B[i] = (rand() % 10);
}

ostream& operator<<(ostream& out, Iter_method& m)
{
	for (int i = 0; i < m.X.size(); ++i)
		out << "x"<<i<<" = " << m.X[i] << "; ";
	return out;
}
