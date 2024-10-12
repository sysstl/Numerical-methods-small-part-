#include "SplineInterpolation.h"
#include "Gauss.h"
#include <iomanip>
#include <fstream>

struct Point
{
	Point(double x, double y) : x{ x }, y{ y } {}
	double x;
	double y;
};

SplineInterpolation::SplineInterpolation(double(*f_)(double))
{
	f = f_;
}

SplineInterpolation::SplineInterpolation(string fiename, double x, double y)
{
	ifstream in(fiename); // окрываем файл для чтения

	if (in.is_open())
	{
		double xx, yy;
		while (in >> xx >> yy)
		{
			experim_data.push_back(Point{ xx, yy });
		}
	}
	in.close();

	for (int i = 0; i < experim_data.size(); i++)
	{
		experim_data[i].x *= x;
		experim_data[i].y *= y;
	}
}

bool SplineInterpolation::InterpolateCubic()
{
	int number_of_knots = experim_data.size();
	double min_x = experim_data[0].x;
	double max_x = experim_data[experim_data.size() - 1].x;
	double d_ = (max_x - min_x) / (number_of_knots - 1);
	X = vector<double>(number_of_knots);
	for (int i = 0; i < number_of_knots; ++i)
	{
		X[i] = min_x;
		min_x += d_;
	}

	Fx = vector<double>(number_of_knots);
	for (int i = 0; i < number_of_knots; ++i)
		Fx[i] = experim_data[i].y; // f(X[i]);

	a = vector<double>(number_of_knots - 1);
	for (int i = 0; i < number_of_knots - 1; ++i)
		a[i] = Fx[i];

	vector<vector<double>> A(number_of_knots, vector<double>(number_of_knots, 0));
	vector<double> B(number_of_knots, 0);
	A[0][0] = -(X[2] - X[1]);
	A[0][1] = X[2] - X[0];
	A[0][2] = -(X[1] - X[0]);
	for (int i = 1; i < number_of_knots - 1; ++i)
	{
		A[i][i - 1] = X[i] - X[i - 1];
		A[i][i] = 2 * (X[i + 1] - X[i - 1]);
		A[i][i + 1] = X[i + 1] - X[i];
		B[i] = 3 * ((Fx[i + 1] - Fx[i]) / (X[i + 1] - X[i]) - (Fx[i] - Fx[i - 1]) / (X[i] - X[i - 1]));
	}
	A[number_of_knots - 1][number_of_knots - 3] = X[number_of_knots - 1] - X[number_of_knots - 2];
	A[number_of_knots - 1][number_of_knots - 2] = -(X[number_of_knots - 1] - X[number_of_knots - 3]);
	A[number_of_knots - 1][number_of_knots - 1] = X[number_of_knots - 2] - X[number_of_knots - 3];
	Gauss g(A, B);
	//cout << g << '\n';
	g.SolveSLAEMainElement();
	//cout << g<<'\n';
	c = g.ReturnAnswer();

	b = vector<double>(number_of_knots - 1);
	for (int i = 0; i < number_of_knots - 1; ++i)
		b[i] = (Fx[i + 1] - Fx[i]) / (X[i + 1] - X[i]) - (c[i + 1] + 2 * c[i]) * (X[i + 1] - X[i]) / 3;

	d = vector<double>(number_of_knots - 1);
	for (int i = 0; i < number_of_knots - 1; ++i)
		d[i] = (c[i + 1] - c[i]) / (3 * (X[i + 1] - X[i]));

	return true;
}

bool SplineInterpolation::InterpolateSquare()
{
	int number_of_knots = experim_data.size();
	double min_x = experim_data[0].x;
	double max_x = experim_data[experim_data.size() - 1].x;
	double d_ = (max_x - min_x) / (number_of_knots - 1);
	X = vector<double>(number_of_knots);
	for (int i = 0; i < number_of_knots; ++i)
	{
		X[i] = min_x;
		min_x += d_;
	}

	Fx = vector<double>(number_of_knots);
	for (int i = 0; i < number_of_knots; ++i)
		Fx[i] = experim_data[i].y;// f(X[i]);

	a = vector<double>(number_of_knots - 1);
	for (int i = 0; i < number_of_knots - 1; ++i)
		a[i] = Fx[i];

	vector<vector<double>> A(number_of_knots, vector<double>(number_of_knots, 0));
	vector<double> B(number_of_knots, 0);
	for (int i = 0; i < number_of_knots - 1; ++i)
	{
		A[i][i] = X[i + 1] - X[i];
		A[i][i + 1] = X[i + 1] - X[i];
		B[i] = 2 * (Fx[i + 1] - Fx[i]);
	}
	A[number_of_knots - 1][number_of_knots - 3] = X[number_of_knots - 1] - X[number_of_knots - 2];
	A[number_of_knots - 1][number_of_knots - 2] = -(X[number_of_knots - 1] - X[number_of_knots - 3]);
	A[number_of_knots - 1][number_of_knots - 1] = X[number_of_knots - 2] - X[number_of_knots - 3];
	Gauss g(A, B);
	//cout << g << '\n';
	g.SolveSLAEMainElement();
	//cout << g<<'\n';
	b = g.ReturnAnswer();

	c = vector<double>(number_of_knots - 1);
	for (int i = 0; i < number_of_knots - 1; ++i)
		c[i] = (b[i + 1] - b[i]) / (2 * (X[i + 1] - X[i]));

	d = vector<double>(number_of_knots - 1, 0);
	return true;
}

double SplineInterpolation::GetValueInPoint(double x)
{
	int i = 0;
	for (int j = 1; j < X.size(); ++j)
		if (x > X[j] && i != X.size() - 2)
			i++;
	double res = 0;
	res = a[i] + b[i] * (x - X[i]) + c[i] * (x - X[i])* (x - X[i]) + d[i] * (x - X[i])* (x - X[i])* (x - X[i]);
	return res;
}

double SplineInterpolation::GetMaxMistake(int n)
{
	double d_ = (X[X.size() - 1] - X[0]) / (n - 1);
	double x = X[0];
	double max_mist = 0;
	for (int i = 0; i < n; ++i)
	{
		double mist = abs(GetValueInPoint(x) - f(x));
		if (mist > max_mist)
			max_mist = mist;
		x += d_;
	}
	return max_mist;
}

ostream& operator<<(ostream& out, const SplineInterpolation& spln)
{
	out << setprecision(10) << "Plot[Piecewise[{";
	for (int i = 0; i < spln.X.size() - 1; ++i)
	{
		out << "{" << spln.a[i];
		if (spln.b[i] != 0)
		{
			if (spln.b[i] < 0)
				out << spln.b[i] << "(x";
			else
				out << "+" << spln.b[i] << "(x";
			if (spln.X[i] < 0)
				out << "+" << -spln.X[i] << ")";
			else
				out << -spln.X[i] << ")";
		}
		if (spln.c[i] != 0)
		{
			if (spln.c[i] < 0)
				out << spln.c[i] << "(x";
			else
				out << "+" << spln.c[i] << "(x";
			if (spln.X[i] < 0)
				out << "+" << -spln.X[i] << ")^2";
			else
				out << -spln.X[i] << ")^2";
		}
		if (spln.d[i] != 0 && !spln.d.empty())
		{
			if (spln.d[i] < 0)
				out << spln.d[i] << "(x";
			else
				out << "+" << spln.d[i] << "(x";
			if (spln.X[i] < 0)
				out << "+" << -spln.X[i] << ")^3";
			else
				out << -spln.X[i] << ")^3";
		}
		out << "," << spln.X[i] << "<=x<=" << spln.X[i + 1] << "}";
		if (i != spln.X.size() - 2)
			out << ",";
	}
	out << "}],{x," << spln.X[0] << "," << spln.X[spln.X.size() - 1] << "}]";
	return out;
}
