#pragma once
#include <vector>
#include <iostream>

using namespace std;

struct Point;

class SplineInterpolation
{
private:
	double(*f)(double);
	vector<Point> experim_data;
	vector<double> a;
	vector<double> b;
	vector<double> c;
	vector<double> d;
	vector<double> X;
	vector<double> Fx;
	bool interpolated = false;

public:
	SplineInterpolation(double(*f_)(double));
	SplineInterpolation(string fiename, double x, double y);
	bool InterpolateCubic();
	bool InterpolateSquare();
	friend ostream& operator<<(ostream& out, const SplineInterpolation& spln);
	double GetValueInPoint(double x);
	double GetMaxMistake(int n);
};
