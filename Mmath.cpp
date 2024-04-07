#include "Mmath.h"

VecD operator-(const VecD& vec1, const VecD& vec2)
{
	if (vec1.size() != vec2.size())	throw std::exception("Invalid sizes");

	VecD Result(vec1);

	for (size_t i = 0; i < vec2.size(); ++i) {
		Result[i] -= vec2[i];
	}

	return Result;
}

VecD operator*(const VecD& vec1, double value)
{
	VecD Result(vec1.size());
	for (size_t i = 0; i < vec1.size(); ++i) {
		Result[i] = vec1[i] * value;
	}
	return Result;
}

//--------Timer-----------------//

void Timer::Start() {
	start = std::chrono::high_resolution_clock::now();
}

double Timer::Finish()
{
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	return duration.count();
}

//---------Matrix--------------//

Matrix::Matrix(size_t Rows, size_t Cols) : _Rows(Rows), _Cols(Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		this->_Values.push_back(VecD(Cols, 0.0));
	}
}

Matrix::Matrix(VecVecD& values) : _Values(values)
{
	this->_Rows = values.size();
	if (this->_Rows != 0) {
		this->_Cols = values[0].size();
	}
	else {
		this->_Cols = 0;
	}
}

Matrix::Matrix(VecD values, bool gorizontal)
{
	if (gorizontal)
	{
		this->_Rows = 1;
		this->_Cols = values.size();
		this->_Values.push_back(values);
	}
	else
	{
		this->_Cols = 1;
		this->_Rows = values.size();
		for (size_t i = 0; i < this->_Rows; ++i) {
			this->_Values.push_back(VecD(1, values[i]));
		}
	}
}

void Matrix::Resize(size_t Rows, size_t Cols)
{
	this->_Rows = Rows;
	this->_Cols = Cols;
	for (size_t i = 0; i < Rows; ++i) {
		this->_Values.push_back(VecD(Cols, 0.0));
	}
}

size_t Matrix::GetRows()
{
	return this->_Rows;
}

size_t Matrix::GetCols()
{
	return this->_Cols;
}

double Matrix::GetDeterminant()
{
	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
	//double ValueToMulty;
	Matrix temp(*this);

	if (Rows != Cols)	throw std::exception("Imposible to calc determinant");

	try
	{
		NumberOfChanges = temp._MakeItTriangle();

		if (abs(temp._Values[Rows - 1][Cols - 1]) <= EPS)	return 0;

		double Result = 1;

		for (int i = 0; i < Rows; ++i) {
			Result *= temp._Values[i][i];
		}

		return (Result * pow(-1, static_cast<double>(NumberOfChanges)));
	}
	catch (std::exception&)
	{
		return 0;
	}
}

Matrix Matrix::GetTransposed()
{
	Matrix Result(this->_Cols, this->_Rows);

	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[j][i] = this->_Values[i][j];
		}
	}

	return Result;
}

Matrix Matrix::GetReverse()
{
	size_t size = this->_Values.size();
	VecVecD mat(size, VecD(2 * size));

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			mat[i][j] = this->_Values[i][j];
		}
		mat[i][i + size] = 1.0;
	}

	double Divisor, Multiplier;

	for (size_t i = 0; i < size; ++i)
	{
		Divisor = mat[i][i];

		for (size_t j = 0; j < 2 * size; ++j) {
			mat[i][j] /= Divisor;
		}

		for (size_t k = 0; k < size; ++k)
		{
			if (k != i)
			{
				Multiplier = mat[k][i];
				for (size_t j = 0; j < 2 * size; ++j) {
					mat[k][j] -= Multiplier * mat[i][j];
				}
			}
		}
	}

	VecVecD Reveresed(size, VecD(size));
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			Reveresed[i][j] = mat[i][j + size];
		}
	}
	return Matrix(Reveresed);
}

void Matrix::SwitchRows(size_t Which, size_t Where)
{
	if (Which >= this->_Values.size() || Where >= this->_Values.size())	throw "Invalid Index";

	VecD temp = this->_Values[Where];
	this->_Values[Where] = this->_Values[Which];
	this->_Values[Which] = temp;
}

size_t Matrix::GetIndexOfBiggestValue(size_t Row, size_t Col)
{
	double BiggestValue = 0;
	size_t IndexOfRow = Row;

	for (size_t i = Row + 1; i < this->_Rows; ++i)
	{
		if (abs(this->_Values[i][Col]) > BiggestValue)
		{
			BiggestValue = abs(this->_Values[i][Col]);
			IndexOfRow = i;
		}
	}

	return IndexOfRow;
}

double Matrix::GetNorma()
{
	double result = 0;

	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			result += (this->_Values[i][j] * this->_Values[i][j]);
		}
	}

	return sqrt(result);
}

Matrix Matrix::operator-(const Matrix& other)
{
	if (this->_Rows != other._Rows || this->_Cols != other._Cols)	throw std::exception("Imposible to subtract these matrices");

	Matrix Result(this->_Rows, this->_Cols);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[i][j] = this->_Values[i][j] - other._Values[i][j];
		}
	}

	return Result;
}

Matrix Matrix::operator*(const Matrix& other)
{
	if (this->_Cols != other._Rows)		throw std::exception("Imposible to multiply these matrices");

	Matrix Result(this->_Rows, other._Cols);
	for (size_t i = 0; i < Result._Rows; ++i) {
		for (size_t j = 0; j < Result._Cols; ++j) {
			for (size_t k = 0; k < this->_Cols; ++k) {
				Result._Values[i][j] += this->_Values[i][k] * other._Values[k][j];
			}
		}
	}

	return Result;
}

Matrix Matrix::operator*(double Value)
{
	Matrix Result(*this);

	for (size_t i = 0; i < Result.GetRows(); ++i) {
		for (size_t j = 0; j < Result.GetCols(); ++j) {
			Result[i][j] *= Value;
		}
	}

	return Result;
}

VecD& Matrix::operator[](size_t Index)
{
	if (Index >= this->_Rows)	throw std::exception("Invalid index");

	return this->_Values[Index];
}

size_t Matrix::_MakeItTriangle()
{
	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
	double ValueToMulty;

	for (size_t i = 0; i < Rows; ++i)
	{
		if (abs(this->_Values[i][i]) <= EPS)
		{
			this->SwitchRows(i, this->GetIndexOfBiggestValue(i, i));
			if (abs(this->_Values[i][i]) <= EPS)	throw std::exception("Det = 0");
			NumberOfChanges++;
		}
		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Values[j][i] / this->_Values[i][i];
			for (size_t k = i; k < Cols; ++k) {
				this->_Values[j][k] -= (ValueToMulty * this->_Values[i][k]);
			}
		}
	}

	return NumberOfChanges;
}

ostream& operator <<(ostream& out,  Matrix& mtr)
{
	//int Row = mtr.GetRows();
	//int Col = mtr.GetCols();
	for (int i = 0; i < mtr.GetRows(); i++)
	{
		for (int j = 0; j < mtr.GetCols(); j++)
		{
			//mtr._Values;
			out << mtr[i][j] << "\t";
		}
		out << "\n";
	}
	return out;
}

//------------SLAE--------------//

SLAE::SLAE(size_t Rows, size_t Cols) : _Matrix(Rows, Cols), _FreeCol(VecD(Cols, 0.0))
{}

SLAE::SLAE(VecVecD& Values, VecD& FreeCol) : _Matrix(Values), _FreeCol(FreeCol)
{}

SLAE::SLAE(Matrix& Matrix, VecD& FreeCol) : _Matrix(Matrix), _FreeCol(FreeCol)
{}

VecD SLAE::GetSolution()
{
	SLAE temp = *this;
	size_t size = temp._Matrix.GetRows();

	if ((size != temp._Matrix.GetCols()) || size == 0 || temp._Matrix.GetCols() == 0)	throw std::exception("Imposible to solve this slae");

	VecD result(size, 0.0);
	temp._MakeItTriangle();
	double Solution;
	for (int64_t i = size - 1; i >= 0; --i)
	{
		Solution = temp._FreeCol[i];
		for (int64_t j = size - 1; j > i; --j) {
			Solution -= result[j] * temp._Matrix[i][j];
		}

		Solution /= temp._Matrix[i][i];

		if (abs(Solution) < EPS) {
			result[i] = 0;
		}
		else {
			result[i] = Solution;
		}
	}
	return result;
}

VecD SLAE::SolveProgon()
{
	size_t Size = this->_Matrix.GetRows();
	VecD alpha(Size, 0), beta(Size, 0);

	alpha[1] = (-this->_Matrix[0][1]) / this->_Matrix[0][0];
	beta[1] = this->_FreeCol[0] / this->_Matrix[0][0];

	for (size_t i = 2; i < Size; ++i)
	{
		alpha[i] = (-this->_Matrix[i - 1][i]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
		beta[i] = (this->_FreeCol[i - 1] - this->_Matrix[i - 1][i - 2] * beta[i - 1]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
	}

	VecD Result(Size, 0);

	Result[Size - 1] = (this->_FreeCol[Size - 1] - this->_Matrix[Size - 1][Size - 2] * beta[Size - 1]) / (this->_Matrix[Size - 1][Size - 1] + this->_Matrix[Size - 1][Size - 2] * alpha[Size - 1]);

	for (size_t i = Size - 1; i > 0; --i) {
		Result[i - 1] = alpha[i] * Result[i] + beta[i];
	}

	return Result;
}

Matrix SLAE::GetMatrix()
{
	return this->_Matrix;
}

VecD SLAE::GetFreeCol()
{
	return this->_FreeCol;
}

void SLAE::_SwitchRows(size_t Which, size_t Where)
{
	this->_Matrix.SwitchRows(Which, Where);
	double temp = this->_FreeCol[Where];
	this->_FreeCol[Where] = this->_FreeCol[Which];
	this->_FreeCol[Which] = temp;
}

void SLAE::_MakeItTriangle()
{
	size_t Rows = this->_FreeCol.size();
	double ValueToMulty;
	for (size_t i = 0; i < Rows - 1; ++i)
	{
		if (abs(this->_Matrix[i][i]) <= EPS)
		{
			this->_SwitchRows(i, this->_Matrix.GetIndexOfBiggestValue(i, i));
			if (abs(this->_Matrix[i][i]) <= EPS)	throw std::exception("Determinant = 0");
		}

		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Matrix[j][i] / this->_Matrix[i][i];
			for (size_t k = i; k < Rows; ++k) {
				this->_Matrix[j][k] -= (ValueToMulty * this->_Matrix[i][k]);
			}
			this->_FreeCol[j] -= (ValueToMulty * this->_FreeCol[i]);
		}
	}
}

VecD SLAE::GetSolutionIteratively(double accuracy)
{
	SLAE Temp = *this;

	size_t size = Temp._Matrix.GetRows();
	if ((size != Temp._Matrix.GetCols()) || size == 0)	throw std::exception("imposible to solve this slae");

	double ValueToDiv;
	for (size_t i = 0; i < size; ++i)
	{
		ValueToDiv = Temp._Matrix[i][i];
		for (size_t j = 0; j < size; ++j) {
			Temp._Matrix[i][j] /= ValueToDiv;
		}
		Temp._FreeCol[i] /= ValueToDiv;
	}

	VecD Result = Temp._FreeCol;
	double CurrentAccuracy;
	do
	{
		Result = Temp._GetInterSolIteratively(Result);
		CurrentAccuracy = (Matrix(Temp._FreeCol, false) - Temp._Matrix * Matrix(Result, false)).GetNorma();
		CurrentAccuracy /= Matrix(Temp._FreeCol, false).GetNorma();
	} while (CurrentAccuracy >= accuracy);
	return Result;
}

VecD SLAE::GetSolutionRelaxation(double accuracy, double omega)
{
	SLAE Temp = *this;

	size_t size = Temp._Matrix.GetRows();
	if ((size != Temp._Matrix.GetCols()) || size == 0)	throw std::exception("imposible to solve this slae");

	double ValueToDiv;
	for (size_t i = 0; i < size; ++i)
	{
		ValueToDiv = Temp._Matrix[i][i];
		for (size_t j = 0; j < size; ++j) {
			Temp._Matrix[i][j] /= ValueToDiv;
		}
		Temp._FreeCol[i] /= ValueToDiv;
	}

	VecD Result = Temp._FreeCol;
	double CurrentAccuracy;
	do
	{
		Result = Temp._GetInterSolRelaxation(Result, omega);
		CurrentAccuracy = (Matrix(Temp._FreeCol, false) - Temp._Matrix * Matrix(Result, false)).GetNorma();
		CurrentAccuracy /= Matrix(Temp._FreeCol, false).GetNorma();
	} while (CurrentAccuracy >= accuracy);
	return Result;
}

VecD SLAE::_GetInterSolIteratively(VecD& Result)
{
	size_t size = this->_Matrix.GetCols();
	VecD NewResult = this->_FreeCol;

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			if (i != j)
			{
				if (j < i) {
					NewResult[i] -= this->_Matrix[i][j] * NewResult[j];
				}
				else {
					NewResult[i] -= this->_Matrix[i][j] * Result[j];
				}
			}
		}
	}
	return NewResult;
}

VecD SLAE::_GetInterSolRelaxation(VecD& Result, double omega)
{
	size_t size = this->_Matrix.GetCols();
	VecD NewResult = this->_FreeCol;

	for (size_t i = 0; i < size; ++i) {
		NewResult[i] *= omega;
	}

	for (size_t i = 0; i < size; ++i)
	{
		NewResult[i] += (1 - omega) * Result[i];
		for (size_t j = 0; j < size; ++j) {
			if (i != j)
			{
				if (j < i) {
					NewResult[i] -= this->_Matrix[i][j] * NewResult[j] * omega;
				}
				else {
					NewResult[i] -= this->_Matrix[i][j] * Result[j] * omega;
				}
			}
		}
	}
	return NewResult;
}

VecD SLAE::GetSolutionGradientDescent(VecD Start, double Accuracy, double Lambda)
{

	Matrix Result(Start, false);
	do
	{
		Result = Result - (this->_Matrix.GetTransposed() * (this->_Matrix * Result - Matrix(this->_FreeCol, false))) * Lambda;
	} while ((this->_Matrix * Result - Matrix(this->_FreeCol, false)).GetNorma() > Accuracy);

	VecD Res(Result.GetRows());
	for (size_t i = 0; i < Res.size(); ++i) {
		Res[i] = Result[i][0];
	}

	return Res;
}

//-------------Polinomial-----------------//

Polinomial pow(Polinomial pol, size_t power)
{
	Polinomial Result = pol;
	for (size_t i = 1; i < power; ++i) {
		Result *= pol;
	}
	return Result;
}

Polinomial::Polinomial()
{
	this->_Coefficients = VecD(1, 0.0);
}

Polinomial::Polinomial(VecD Coefficients) : _Coefficients(Coefficients)
{}

void Polinomial::Interpolate(VecD X, VecD Y)
{
	size_t size = X.size();
	VecVecD Matrix(size, VecD(size));

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			Matrix[i][j] = pow(X[i], j);
		}
	}

	SLAE slae(Matrix, Y);
	this->_Coefficients = slae.GetSolution();
}

double Polinomial::GetY(double X)
{
	double Result = 0;
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		Result += this->_Coefficients[i] * pow(X, i);
	}
	return Result;
}

size_t Polinomial::GetMaxPower()
{
	return (this->_Coefficients.size() - 1);
}

Polinomial Polinomial::GetDerivative(size_t power)
{
	Polinomial Result(VecD(this->_Coefficients.size() - power, 1.0));
	for (size_t i = 0; i < this->_Coefficients.size() - power; ++i) {
		Result._Coefficients[i] = this->_Coefficients[i + power] * this->_GetCoefOfDeriv(i + power, power);
	}
	return Result;
}

Polinomial Polinomial::operator*(Polinomial& other)
{
	size_t size = ((this->_Coefficients.size() - 1) + (other._Coefficients.size() - 1) + 1);
	VecD Result(size, 0.0);
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		for (size_t j = 0; j < other._Coefficients.size(); ++j) {
			Result[i + j] += this->_Coefficients[i] * other._Coefficients[j];
		}
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator*=(Polinomial other)
{
	size_t size = ((this->_Coefficients.size() - 1) + (other._Coefficients.size() - 1) + 1);
	VecD Result(size, 0.0);
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		for (size_t j = 0; j < other._Coefficients.size(); ++j) {
			Result[i + j] += this->_Coefficients[i] * other._Coefficients[j];
		}
	}
	this->_Coefficients = Result;
	return *this;
}

Polinomial Polinomial::operator+(Polinomial& other)
{
	VecD Result(this->_Coefficients), ToSum(other._Coefficients);
	(this->_Coefficients.size() > other._Coefficients.size()) ? ToSum.resize(Result.size(), 0.0) : Result.resize(ToSum.size(), 0.0);
	for (size_t i = 0; i < ToSum.size(); ++i) {
		Result[i] += ToSum[i];
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator+=(Polinomial other)
{
	*this = *this + other;
	return *this;
}

Polinomial Polinomial::operator-(Polinomial& other)
{
	VecD Result(this->_Coefficients), ToSub(other._Coefficients);
	(this->_Coefficients.size() > other._Coefficients.size()) ? ToSub.resize(Result.size(), 0.0) : Result.resize(ToSub.size(), 0.0);
	for (size_t i = 0; i < ToSub.size(); ++i) {
		Result[i] -= ToSub[i];
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator*=(double value)
{
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		this->_Coefficients[i] *= value;
	}
	return *this;
}

Polinomial Polinomial::operator*(double value)
{
	Polinomial Result = *this;
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		Result._Coefficients[i] = this->_Coefficients[i] * value;
	}
	return Result;
}

Polinomial Polinomial::operator/=(double value)
{
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		this->_Coefficients[i] /= value;
	}
	return *this;
}

double Polinomial::operator()(double X)
{
	return this->GetY(X);
}

VecD& Polinomial::operator[](int Index)
{
	return this->_Coefficients;
}

double Polinomial::_GetCoefOfDeriv(size_t Power, size_t Derivative)
{
	double Result = 1;
	for (size_t i = 0; i < Derivative; i++) {
		Result *= (Power - i);
	}
	return Result;
}

//--------------Splayn----------------//

Splayn::Splayn()
{}

void Splayn::Interpolate(VecD X, VecD Y, size_t power)
{
	this->_X = X;
	this->_Polinoms.clear();
	size_t NumberOfPoints = X.size();
	size_t size = (power + 1) * (NumberOfPoints - 1);
	Matrix mat(size, size);
	VecD Free(size, 0.0);
	int num = 0;
	for (size_t i = 0; i < (NumberOfPoints - 1); ++i)
	{
		for (size_t j = 0; j < 2; ++j)
		{
			for (size_t k = 0; k < (power + 1); ++k)
			{
				mat[(2 * i) + j][(i * (power + 1)) + k] = pow((X[i + j] - X[i]), k);
			}
			Free[(2 * i) + j] = Y[i + j];
			num++;
		}
	}

	for (size_t t = 1; t < power; ++t)
	{
		for (size_t i = 0; i < (NumberOfPoints - 2); ++i)
		{
			for (size_t j = 0; j < 2; ++j)
			{
				for (size_t k = t; k < (power + 1); ++k)
				{
					mat[2 * (NumberOfPoints - 1) - (NumberOfPoints - 2) + (t * (NumberOfPoints - 2)) + i][k + (j * (power + 1)) + (i * (power + 1))] = this->_GetCoeffDeriv(k, t) * pow(X[i + 1] - X[i + j], (k - t)) * pow(-1, j);
				}
			}
			num++;
		}
	}

	for (size_t t = 2; t <= power / 2 + 1; ++t)
	{
		for (size_t k = t; k < (power + 1); ++k)
		{
			mat[num][k] = this->_GetCoeffDeriv(k, t) * pow(X[0] - X[0], (k - t));
		}
		num++;
	}

	for (size_t t = 2; t <= (power)-(power / 2); ++t)
	{
		for (size_t k = power; k >= t; --k)
		{
			mat[num][size - 1 - power + k] = this->_GetCoeffDeriv(k, t) * pow(X[NumberOfPoints - 1] - X[NumberOfPoints - 2], (k - t));
		}
		num++;
	}

	if (power == 2) {
		mat[size - 1][size - 1] = 2;
	}

	SLAE slae(mat, Free);
	VecD coeffs = slae.GetSolution();
	Polinomial temp;
	for (size_t i = 0; i < NumberOfPoints - 1; ++i)
	{
		temp = Polinomial(VecD{ coeffs[(power + 1) * i] });
		for (size_t j = 1; j < (power + 1); ++j) {
			temp += pow(Polinomial(VecD{ -X[i], 1 }), j) * coeffs[(power + 1) * i + j];
		}
		this->_Polinoms.push_back(temp);
	}
}

double Splayn::GetY(double X)
{
	for (size_t i = 0; i < this->_X.size(); ++i)
	{
		if (X <= this->_X[i])
		{
			if (i != 0) {
				return (this->_Polinoms[i - 1].GetY(X));
			}
			else {
				return (this->_Polinoms[i].GetY(X));
			}
		}
	}
	return (this->_Polinoms[_Polinoms.size() - 1].GetY(X));
}

size_t Splayn::_GetCoeffDeriv(size_t value, size_t power)
{
	size_t Result = value;
	for (size_t i = 1; i < power; ++i) {
		Result *= (--value);
	}
	return Result;
}

//-------------------NonlinearEquation---------------//

NonlinearEquation::NonlinearEquation(std::function<double(double)> Function) : _Function(Function)
{}

double NonlinearEquation::Integrate(double range1, double range2, double accuracy)
{
	double h = (range2 - range1) / 1.0;

	double NewResult = this->_GetIntermediateResult(range1, range2, h);
	double Result, Fault, NewFault = INF;

	do
	{
		Fault = NewFault;
		Result = NewResult;
		h /= 2.0;
		NewResult = this->_GetIntermediateResult(range1, range2, h);
		NewFault = ((NewResult - Result) / (pow(2.0, 4.0) - 1));

		if (abs(NewFault) >= abs(Fault))	throw std::exception("accuracy is to low");

	} while (abs(NewFault) >= accuracy);

	return NewResult + NewFault;
}

double NonlinearEquation::GetDerivative(double X, size_t power)
{
	VecD Xs, Ys;
	for (double x = (X - 0.2); x < (X + 0.25); x += 0.05)
	{
		Xs.push_back(x);
		Ys.push_back(this->_Function(x));
	}
	Polinomial pol;
	pol.Interpolate(Xs, Ys);
	return pol.GetDerivative(power).GetY(X);
}

double NonlinearEquation::GetRootNewton(double Point, double Accuracy)
{
	double Prev;
	double dx = 1e-5;
	//int Iters = 0;

	do
	{
		//++Iters;
		Prev = Point;
		Point -= this->_Function(Point) / ((this->_Function(Point + dx) - this->_Function(Point - dx)) / (2 * dx));
	} while (abs(Prev - Point) > Accuracy);

	//-------------//

	/*size_t Multiplicity = 1;
	while (abs(this->GetDerivative(Point, Multiplicity)) < (Accuracy * pow(10, Multiplicity))){
		++Multiplicity;
	}
	std::cout << "Multiplicity : " << Multiplicity << '\n';*/

	//-------------//

	//std::cout << "Newton it : " << Iters << " Calls : " << Iters * 3 << '\n';

	return Point;
}

double NonlinearEquation::GetExtremumGoldenRatio(double Left, double Right, double Accuracy, bool Min)
{
	double GoldenR = (1 + sqrt(5)) / 2.0;
	double X_Left, X_Right;

	do
	{
		X_Right = Left + (Right - Left) / GoldenR;
		X_Left = Right - (Right - Left) / GoldenR;

		if (this->_Function(X_Left) > this->_Function(X_Right))
		{
			if (Min) {
				Left = X_Left;
			}
			else {
				Right = X_Right;
			}
		}
		else
		{
			if (Min) {
				Right = X_Right;
			}
			else {
				Left = X_Left;
			}
		}

	} while (abs(Left - Right) > Accuracy);

	return (Left + Right) / 2;
}

double NonlinearEquation::GetExtremumNewton(double X, double Accuracy)
{
	double Prev;
	double dx = 1e-5;

	do
	{
		Prev = X;
		X -= ((this->_Function(X + dx) - this->_Function(X - dx)) / (2 * dx)) /
			((this->_Function(X + dx) - 2 * this->_Function(X) + this->_Function(X - dx)) / (dx * dx));
	} while (abs(Prev - X) > Accuracy);

	return X;
}

double NonlinearEquation::GetMinimumCauchy(double X, double Start_T, double Min_T)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double Min = X, New_Min, Res = X;
	double T;
	double E_Old, E_New;

	for (size_t k = 1; k < (Start_T / Min_T) + 1; ++k)
	{
		T = Start_T / k;
		while (true)
		{
			New_Min = Min + dist(gen) * 4 - 2;
			E_Old = this->_Function(Min);
			E_New = this->_Function(New_Min);

			if ((exp(-(E_New - E_Old) / T) >= dist(gen)))
			{
				Min = New_Min;
				if (this->_Function(New_Min) < this->_Function(Res)) {
					Res = New_Min;
				}
				break;
			}
		}
	}

	return Res;
}

double NonlinearEquation::_GetIntermediateResult(double range1, double range2, double h)
{
	double Result = 0;
	for (double x = range1; x <= (range2 - h / 2.0); x += h) {
		Result += ((this->_Function(x) + 4.0 * this->_Function(x + h / 2.0) + this->_Function(x + h)) / 6.0) * h;
	}
	return Result;
}

//---------------------------FunctionOfManyArg---------//

FunctionOfManyArg::FunctionOfManyArg(size_t NumberOfArsg, std::function<double(VecD)> Function) : _NumberOfArgs(NumberOfArsg), _Function(Function)
{}

VecD FunctionOfManyArg::FindMinGradientDescent(VecD Start, double Accuracy, double Lambda)
{
	VecD CurrentPoint(Start);
	//size_t k = 1;

	do
	{
		Start = CurrentPoint;
		CurrentPoint = Start - this->Gradient(Start) * Lambda;
		//++k;
	} while (Matrix(CurrentPoint - Start).GetNorma() > Accuracy);

	return CurrentPoint;
}

VecD FunctionOfManyArg::Gradient(VecD Point)
{
	VecD Result(this->_NumberOfArgs, 0);
	VecD DeltaPoint(Point);
	double Delta = 1e-5;

	for (size_t i = 0; i < this->_NumberOfArgs; ++i)
	{
		DeltaPoint[i] += Delta;
		Point[i] -= Delta;

		Result[i] = (this->_Function(DeltaPoint) - this->_Function(Point)) / (2 * Delta);

		DeltaPoint[i] -= Delta;
		Point[i] += Delta;
	}

	return Result;
}

//--------------------SDE------------------//
//
//SDE::SDE(std::vector<std::function<double(VecD&, double)>> Functions) : _Functions(Functions)
//{}
//
//VecVecD SDE::SolveEuler(VecD Start, double Start_T, double Finish_T, double dT)
//{
//	VecVecD Result(1, VecD(Start.size() + 1));
//
//	for (size_t i = 0; i < Start.size(); ++i) {
//		Result[0][i] = Start[i];
//	}
//	Result[0][Start.size()] = Start_T;
//
//	double Temp, IntermediateValue;
//	VecD NextRes(Start), TempRes(Start);
//	NextRes.push_back(Start_T);
//
//	for (double t = Start_T + dT; t <= Finish_T + (dT / 2); t += dT)
//	{
//		for (size_t i = 0; i < Start.size(); i++)
//		{
//			TempRes[i] += dT * (this->_Functions[i](TempRes, t - dT));
//			NextRes[i] += dT * (this->_Functions[i](NextRes, t - dT) + this->_Functions[i](TempRes, t)) / 2;
//		}
//		NextRes[NextRes.size() - 1] = t;
//		Result.push_back(NextRes);
//	}
//
//	return Result;
//}
//
//VecVecD SDE::SolveRungeKutt(VecD Start, double Start_T, double Finish_T, double dT, bool ClassicType)
//{
//	VecD C(4, 0), B(4, 0);
//	VecVecD A(4, VecD(4, 0));
//
//	if (ClassicType)
//	{
//		C[1] = C[2] = 0.5;
//		A[1][0] = A[2][1] = 0.5;
//		B[0] = B[3] = 1.0 / 6.0;
//		B[1] = B[2] = 1.0 / 3.0;
//	}
//	else
//	{
//		C[1] = C[2] = 1.0 / 3.0;
//		A[1][0] = 1.0 / 3.0;
//		A[2][0] = -1.0 / 3.0;
//		A[2][1] = A[3][0] = 1;
//		A[3][1] = -1;
//		B[0] = B[3] = 1.0 / 8.0;
//		B[1] = B[2] = 3.0 / 8.0;
//	}
//	C[3] = 1;
//	A[3][2] = 1;
//
//	VecVecD Result(1, VecD(Start.size() + 1));
//
//	for (size_t i = 0; i < Start.size(); ++i) {
//		Result[0][i] = Start[i];
//	}
//	Result[0][Start.size()] = Start_T;
//
//	VecD Ki(4, 0);
//	VecD NextRes(Start), Temp;
//	double PrevValue;
//	NextRes.push_back(Start_T);
//
//	for (double t = Start_T + dT; t <= Finish_T + (dT / 2); t += dT)
//	{
//		Temp = NextRes;
//		for (size_t i = 0; i < Start.size(); i++)
//		{
//			Ki[0] = this->_Functions[i](Temp, t - dT);
//
//			PrevValue = Temp[i];
//			Temp[i] += A[1][0] * dT * Ki[0];
//			Ki[1] = this->_Functions[i](Temp, t - dT + C[1] * dT);
//			Temp[i] = PrevValue;
//
//			Temp[i] += A[2][0] * dT * Ki[0] + A[2][1] * dT * Ki[1];
//			Ki[2] = this->_Functions[i](Temp, t - dT + C[2] * dT);
//			Temp[i] = PrevValue;
//
//			Temp[i] += A[3][0] * dT * Ki[0] + A[3][1] * dT * Ki[1] + A[3][2] * dT * Ki[2];
//			Ki[3] = this->_Functions[i](Temp, t - dT + C[3] * dT);
//			Temp[i] = PrevValue;
//
//			NextRes[i] += dT * (B[0] * Ki[0] + B[1] * Ki[1] + B[2] * Ki[2] + B[3] * Ki[3]);
//		}
//		NextRes[NextRes.size() - 1] = t;
//		Result.push_back(NextRes);
//	}
//
//	return Result;
//}
//
//VecVecD SDE::SolveRichardson(VecD Start, double Start_T, double Finish_T, double dT, double Accuracy)
//{
//	VecVecD Result(1, VecD(Start.size() + 1));
//
//	for (size_t i = 0; i < Start.size(); ++i) {
//		Result[0][i] = Start[i];
//	}
//	Result[0][Start.size()] = Start_T;
//
//	double p = 4;
//	double k;
//	double Norma;
//	VecD N(Start), Res(Start);
//	Res.push_back(Start_T);
//	VecVecD Temp1, Temp2;
//
//	for (double t = Start_T + dT; t < Finish_T + (dT / 2); t += dT)
//	{
//		k = 1;
//		do
//		{
//			k *= 2;
//			Norma = 0;
//
//			Temp1 = this->SolveRungeKutt(N, t - dT, t, dT / k);
//			Temp2 = this->SolveRungeKutt(N, t - dT, t, dT / k / 2);
//
//			for (size_t i = 0; i < Temp1[0].size() - 1; ++i) {
//				Norma += pow(Temp1[Temp1.size() - 1][i] - Temp2[Temp2.size() - 1][i], 2);
//			}
//			Norma = sqrt(Norma);
//
//		} while (Norma / (pow(2, p) - 1) > Accuracy);
//
//		for (size_t i = 0; i < Temp2[0].size() - 1; ++i)
//		{
//			N[i] = Temp2[Temp2.size() - 1][i];
//			Res[i] = N[i];
//		}
//		Res[Res.size() - 1] = t;
//		Result.push_back(Res);
//	}
//
//	return Result;
//}
//
//VecD SolveBoundaryProblem(double Left, double Right, double Y_Left, double Y_Right, const std::vector<std::function<double(double)>>& Coefs, double dx)
//{
//	size_t Size = static_cast<size_t>(((Right - Left) / dx));
//
//	Matrix Mat(Size, Size);
//	VecD Free(Size, 0);
//
//	Free[0] = Y_Left;
//	Free[Size - 1] = Y_Right;
//
//	Mat[0][0] = Coefs[3](Left);
//	Mat[0][1] = Coefs[4](Left);
//
//	Mat[Size - 1][Size - 2] = Coefs[5](Right);
//	Mat[Size - 1][Size - 1] = Coefs[6](Right);
//
//	for (size_t i = 1; i < Size - 1; ++i)
//	{
//		Mat[i][i - 1] = Coefs[0](Left + i * dx);
//		Mat[i][i] = Coefs[1](Left + i * dx);
//		Mat[i][i + 1] = Coefs[2](Left + i * dx);
//	}
//
//	SLAE slae(Mat, Free);
//	return slae.SolveProgon();
//}
