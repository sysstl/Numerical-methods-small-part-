#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

typedef std::vector<double> VecD;
typedef std::vector<std::vector<double>> VecVecD;

VecD operator-(const VecD& vec1, const VecD& vec2);
VecD operator*(const VecD& vec1, double value);

const double EPS = ((double)(pow(10, -13)));
const double PI = (acos(-1));
const double INF = pow(10, 309);

class Timer;
class SLAE;
//class Matrix;
class Polinomial;
class Splayn;
class NonlinearEquation;
class FunctionOfManyArg;
class SDE;

class Timer
{
public:

	void Start();
	double Finish();

private:

	std::chrono::time_point<std::chrono::steady_clock> start, end;
};

class Matrix
{
public:
	Matrix(size_t Rows = 0, size_t Cols = 0);
	Matrix(VecVecD& values);
	Matrix(VecD values, bool gorizontal = true);
	void Resize(size_t Rows, size_t Cols);
	size_t GetRows();
	size_t GetCols();
	double GetDeterminant();
	Matrix GetTransposed();
	Matrix GetReverse();
	void SwitchRows(size_t Which, size_t Where);
	size_t GetIndexOfBiggestValue(size_t Row, size_t Col);
	double GetNorma();
	void PrintMatrix(ostream& out);

	//---------------------------//

	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix operator*(double Value);
	VecD& operator[](size_t Index);

private:

	size_t _MakeItTriangle();
	//---------------------------//
	VecVecD _Values;
	size_t _Rows;
	size_t _Cols;
};

class SLAE
{
public:

	SLAE(size_t Rows, size_t Cols);
	SLAE(VecVecD& Values, VecD& FreeCol);
	SLAE(Matrix& Matrix, VecD& FreeCol);
	VecD GetSolution();
	VecD SolveProgon();
	Matrix GetMatrix();
	VecD GetFreeCol();

	VecD GetSolutionIteratively(double accuracy);
	VecD GetSolutionRelaxation(double accuracy, double omega);
	VecD GetSolutionGradientDescent(VecD Start, double Accuracy, double Lambda = 0.01);

private:

	VecD _GetInterSolIteratively(VecD& Result);
	VecD _GetInterSolRelaxation(VecD& Result, double omega);

	void _SwitchRows(size_t Which, size_t Where);
	void _MakeItTriangle();

	Matrix _Matrix;
	VecD _FreeCol;
};

Polinomial pow(Polinomial pol, size_t power);

class Polinomial
{
public:
	Polinomial();
	Polinomial(VecD Coefficients);
	void Interpolate(VecD X, VecD Y);
	double GetY(double X);
	size_t GetMaxPower();
	Polinomial GetDerivative(size_t power);

	Polinomial operator*(Polinomial& other);
	Polinomial operator*=(Polinomial other);
	Polinomial operator+(Polinomial& other);
	Polinomial operator+=(Polinomial other);
	Polinomial operator-(Polinomial& other);
	Polinomial operator*=(double value);
	Polinomial operator*(double value);
	Polinomial operator/=(double value);
	double operator()(double X);
	VecD& operator[](int Index);

private:

	double _GetCoefOfDeriv(size_t Power, size_t Derivative);

	VecD _Coefficients;
};

class Splayn
{
public:
	Splayn();
	void Interpolate(VecD X, VecD Y, size_t power);
	double GetY(double X);

private:

	size_t _GetCoeffDeriv(size_t value, size_t power);

	VecD _X;
	std::vector<Polinomial> _Polinoms;
};

class NonlinearEquation
{
public:
	NonlinearEquation(std::function<double(double)> Function);
	double Integrate(double range1, double range2, double accuracy);
	double GetDerivative(double X, size_t power);

	double GetRootNewton(double Point, double Accuracy);
	double GetExtremumGoldenRatio(double Left, double Right, double Accuracy, bool Min = true);
	double GetExtremumNewton(double X, double Accuracy);
	double GetMinimumCauchy(double X, double Start_T, double Min_T);

private:
	double _GetIntermediateResult(double range1, double range2, double h);

	std::function<double(double)> _Function;
};

class FunctionOfManyArg
{
public:

	FunctionOfManyArg(size_t NumberOfArsg, std::function<double(VecD)> Function);
	VecD FindMinGradientDescent(VecD Start, double Accuracy, double Lambda = 0.001);
	VecD Gradient(VecD Point);

private:

	size_t _NumberOfArgs;
	std::function<double(VecD)> _Function;
};

VecD SolveBoundaryProblem(double Left, double Right, double Y_Left, double Y_Right, const std::vector<std::function<double(double)>>& Coefs, double dx);

//class SDE
//{
//public:
//
//	SDE(std::vector<std::function<double(VecD&, double)>> Functions);
//	VecVecD SolveEuler(VecD Start, double Start_T, double Finish_T, double dT);
//	VecVecD SolveRungeKutt(VecD Start, double Start_T, double Finish_T, double dT, bool ClassicType = false);
//	VecVecD SolveRichardson(VecD Start, double Start_T, double Finish_T, double dT, double Accuracy);
//
//private:
//
//	std::vector<std::function<double(VecD&, double)>> _Functions;
//
//};
