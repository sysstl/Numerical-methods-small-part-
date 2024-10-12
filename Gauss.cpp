#include "Gauss.h"
#include <random>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
Gauss::Gauss(int n_)
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

Gauss::Gauss(vector< vector<double> >& A_, vector<double>& B_)
{
	if (A_.size() != A_[0].size())
	{
		errFlag = true;
		err = "A is not square\n";
	}
	else if (A_.size() != B_.size())
	{
		errFlag = true;
		err = "A and B have different scale\n";
	}
	else
	{
		errFlag = false;
		n = A_.size();
		A = A_;
		B = B_;
	}
}

ostream & operator <<(ostream& out, const Gauss& sistem)
{
	if (sistem.errFlag)
		out << sistem.err;
	else if (sistem.solved)
	{
		for (int i = 0; i < sistem.n; ++i)
		{
			out << 'x' << i + 1 << " = " << setprecision(20) << sistem.B[i] << "\n";
		}
	}
	else
	{
		out << setw(7) << 'A';
		for (int i = 0; i < sistem.n - 1; ++i)
			out << setw(7) << ' ';
		out << setw(7) << 'B' << '\n';
		for (int i = 0; i < sistem.n; ++i)
		{
			for (int j = 0; j < sistem.n; ++j)
				out << setw(7) << sistem.A[i][j];
			out << setw(7) << sistem.B[i];
			out << '\n';
		}
	}
	return out;
}

void Gauss::PrintA(ostream& out)
{
	if (errFlag)
		out << err;
	else
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				out << setw(7) << setprecision(3) << A[i][j];
			out << '\n';
		}
	}
}

bool Gauss::SolveSLAEGeneral()
{
	if (errFlag)
	{
		return false;
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
			A[i][j] /= A[i][i];
		B[i] /= A[i][i];

		for (int j = i + 1; j < n; ++j)
		{
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[i][k] * A[j][i];
			B[j] -= B[i] * A[j][i];
		}
	}

	for (int i = n - 2; i >= 0; --i)
		for (int j = i + 1; j < n; ++j)
			B[i] -= B[j] * A[i][j];
	solved = true;
	return true;
}
bool Gauss::SolveSLAEMainElement()
{
	if (errFlag)
	{
		return false;
	}
	for (int i = 0; i < n; ++i)
	{
		int m = i;
		for (int j = i; j < n; ++j)
		{
			if (abs(A[m][i]) < abs(A[j][i]))
				m = j;
		}
		A[i].swap(A[m]);
		swap(B[i], B[m]);
		if (abs(A[i][i]) <= 0.00000001)
		{
			errFlag = true;
			err = "Determinant equals 0";
			return false;
		}

		for (int j = i + 1; j < n; ++j)
			A[i][j] /= A[i][i];
		B[i] /= A[i][i];

		for (int j = i + 1; j < n; ++j)
		{
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[i][k] * A[j][i];
			B[j] -= B[i] * A[j][i];
		}
	}

	for (int i = n - 2; i >= 0; --i)
		for (int j = i + 1; j < n; ++j)
			B[i] -= B[j] * A[i][j];
	solved = true;
	return true;
}

Progonka::Progonka(int n_)
{
	n = n_;
	srand(time(NULL));
	a = vector<double>(n_);
	b = vector<double>(n_);
	c = vector<double>(n_);
	d = vector<double>(n_);
	//e = vector<double>(n_-1);
	e = vector<double>(n_);
	f = vector<double>(n_);

	for (int i = 0; i < n_; i++)
	{
		a[i] = rand() % 10 + 1; // dt3 * Kv / (dz3 * dz3);
		b[i] = rand() % 10 + 1; // -(1 + 2 * Kv * dt3 / (dz3 * dz3));
		c[i] = rand() % 10 + 1; //  dt3 * Kv / (dz3 * dz3);
		d[i] = rand() % 10 + 1;
		e[i] = rand() % 10 + 1;
		f[i] = rand() % 10 + 1; // - Tv1[i];
	}

	e[3] = 0;
	d[4] = 0;
	e[4] = 0;

	/*cout << endl << endl;

	cout << c[0] << " x0- " << d[0] << " x1 + " << e[0] << " x2 = " << f[0] << endl;
	cout << -b[1] << " x0+ " << c[1] << " x1 - " << d[1] << " x2 + " << e[1] << " x3 = " << f[1] << endl;
	cout << a[2] << " x0- " << b[2] << " x1 + " << c[2] << " x2 - " <<d[2] <<" x3 + "<< e[2] << " x4 = " << f[2] << endl;
	cout << a[3] << " x1- " << b[3] << " x2 + " << c[3] << " x3 - " <<d[3] << " x4 = " << f[3] << endl;
	cout << a[4] << " x2- " << b[4] << " x3 + " << c[4] << " x4 = " << f[4] << endl;
*/
	//b[0] = 1.;
	//c[0] = 0.;
	//f[0] = 5.; // левая граница

	//a[n_ - 1] = 0.;
	//b[n_ - 1] = 1.;
	//c[n_ - 1] = 0.;
	//f[n_ - 1] = 1.;//  правая гарница;


}

Progonka::Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& f_)
{
	n = b_.size();
	a = a_;
	b = b_;
	c = c_;
}

Progonka::Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& d_, vector<double>& e_, vector<double>& f_)
{

}

vector<double> Progonka::SolveThreeDiagonal()
{
	alpha = vector<double>(n);
	beta = vector<double>(n);
	x = vector<double>(n);

	alpha[1] = (-1) * c[0] / b[0];
	beta[1] = f[0] / b[0];

	for (int i = 1; i <= n - 2; i++)
	{
		alpha[i + 1] = (-1) * c[i] / (b[i] + a[i] * alpha[i]);
		beta[i + 1] = ((-1) * a[i] * beta[i] + f[i]) / (b[i] + a[i] * alpha[i]);
	}

	x[n - 1] = ((-1) * a[n - 1] * beta[n - 1] + f[n - 1]) / (b[n - 1] + a[n - 1] * alpha[n - 1]); // надо b[n-1]

	for (int i = n - 1; i >= 1; i--) {
		x[i - 1] = alpha[i] * x[i] + beta[i];
	}

	cout << endl;
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}

	cout << endl << endl;
	//int ij = 6;
	for (int ij = n - 1; ij >= 0; ij--) {
		cout << x[ij - 1] * a[ij] + x[ij] * b[ij] + x[ij + 1] * c[ij] << " =  " << f[ij] << endl;
	}
	return x;
}

vector<double> Progonka::SolveFiveDiagonal()
{
	alpha = vector<double>(n + 3); // n=5
	beta = vector<double>(n + 3);
	gama = vector<double>(n + 3);
	x = vector<double>(n);

	alpha[1] = d[0] / c[0];
	beta[1] = e[0] / c[0];
	gama[1] = f[0] / c[0];
	alpha[2] = (d[1] - beta[1] * b[1]) / (c[1] - b[1] * alpha[1]);
	beta[2] = e[1] / (c[1] - b[1] * alpha[1]);
	gama[2] = (f[1] + b[1] * gama[1]) / (c[1] - b[1] * alpha[1]);

	for (int i = 2; i <= 4; i++)
	{
		alpha[i + 1] = (d[i] + beta[i] * (a[i] * alpha[i - 1] - b[i])) / (c[i] - b[i] * alpha[i]);
		gama[i + 1] = (f[i] - a[i] * gama[i - 1] - gama[i] * (a[i] * alpha[i - 1] - b[i])) / (c[i] - b[i] * alpha[i]);
		beta[i + 1] = e[i] / (c[i] - b[i] * alpha[i]);
	}

	cout << " Table chisel " << endl;
	cout << 1 << "   " << alpha[1] << "   " << beta[1] << "   " << gama[1] << endl;
	cout << 2 << "   " << alpha[2] << "   " << beta[2] << "   " << gama[2] << endl;
	cout << 3 << "   " << alpha[3] << "   " << beta[3] << "   " << gama[3] << endl;
	cout << 4 << "   " << alpha[4] << "   " << beta[4] << "   " << gama[4] << endl;
	cout << 5 << "   " << alpha[5] << "   " << beta[5] << "   " << gama[5] << endl;
	cout << 6 << "   " << alpha[6] << "   " << beta[6] << "   " << gama[6] << endl;

	//for (int i = 1; i < n - 2; i++)
	//{
	//	beta[i + 1] = e[i] / (c[i] - b[i] * alpha[i]);
	//	alpha[i + 2] = (d[i + 1] + beta[i + 1] * (a[i + 1] * alpha[i] - b[i + 1])) / (c[i + 1] - b[i + 1] * alpha[i + 1]);
	//	gama[i + 2] = (f[i + 1] - a[i + 1] * gama[i] - gama[i + 1] * (a[i + 1] * alpha[i] - b[i + 1])) / (c[i + 1] - b[i + 1] * alpha[i + 1]);
	//}

	//for (int i = n - 1; i < n; i++)//i=4
	//{
	//	alpha[i + 1] = (d[i] + beta[i] * (a[i] * alpha[i - 1] - b[i])) / (c[i] - b[i] * alpha[i]);
	//	gama[i + 1] = (f[i] - a[i] * gama[i - 1] - gama[i] * (a[i] * alpha[i - 1] - b[i])) / (c[i] - b[i] * alpha[i]);
	//}

	// обратный ход

	//x[n - 1] = gama[n];
	//x[n - 2] = alpha[n - 1] * x[n - 1] + gama[n - 1];
	x[4] = gama[5];
	x[3] = alpha[4] * x[4] + gama[4];
	//for (int i = n - 3; i >= 0; i--)
	for (int i = 2; i >= 0; i--)
	{
		x[i] = alpha[i + 1] * x[i + 1] - beta[i + 1] * x[i + 2] + gama[i + 1];
	}

	/*x[n - 1] = gama[n];
	x[n - 2] = alpha[n - 1] * x[n - 1] + gama[n - 1];
	for (int i = n - 3; i >= 0; i--)
	{
		x[i] = alpha[i + 1] * x[i + 1] - beta[i + 1] * x[i + 2] + gama[i + 1];
	}*/

	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}

	cout << endl << endl;
	cout << c[0] * x[0] -  d[0] * x[1] + e[0] * x[2] <<"  =  " << f[0] << endl;
	cout << -b[1] * x[0] + c[1] * x[1] - d[1] * x[2] +  e[1] * x[3] <<"  = " << f[1] << endl;
	cout << a[2] * x[0]-  b[2] * x[1] +  c[2] * x[2] -  d[2] * x[3] + e[2] * x[4] <<" = " << f[2] << endl;
	cout << a[3] * x[1]-  b[3] * x[2] +  c[3] * x[3] -  d[3] * x[4] <<"  = " << f[3] << endl;
	cout << a[4] * x[2]-  b[4] * x[3] +  c[4] * x[4] <<"  = " << f[4] << endl;

	/*cout << c[0] * x[0] - d[0] * x[1] + e[0] * x[2] << " =  " << f[0] << endl;
	cout << -b[1] * x[0] + c[1] * x[1] - d[1] * x[2] + e[1] * x[3]<< " =  " << f[1] << endl;
	cout << a[2] * x[0] - b[2] * x[1] + c[2] * x[2] - d[2] * x[3] + e[2] * x[4]<< " =  " << f[2] << endl;
	cout << a[4] * x[2] - b[4] * x[3]  + c[4] * x[4] << " =  " << f[4] << endl;
*/
	cout << endl << endl;
	/*for (int ij = 2; ij <= n-3; ij++) {
			cout << x[ij - 2] * a[ij] - x[ij-1] * b[ij] + x[ij] * c[ij] - x[ij + 1] * d[ij]+ x[ij + 2] * e[ij] << " =  " << f[ij] << endl;
		}

	cout << a[n - 2] * x[n - 4] - b[n - 2] * x[n - 3] + c[n - 2] * x[n - 2] - d[n - 2] * x[n - 1] << " =  " << f[n - 2] << endl;
	cout << a[n - 1] * x[n - 2] - b[n - 1] * x[n - 1] + c[n - 1] * x[n] << " =  " << f[n - 1] << endl;
*/
	return x;
}