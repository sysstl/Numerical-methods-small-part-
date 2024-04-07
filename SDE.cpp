#include "SDE.h"

SDE::SDE(std::vector<std::function<double(VecD&, double)>> Functions) : _Functions(Functions)
{}

VecVecD SDE::SolveEuler(VecD Start, double Start_T, double Finish_T, double dT)
{
	VecVecD Result(1, VecD(Start.size() + 1));

	for (size_t i = 0; i < Start.size(); ++i) {
		Result[0][i] = Start[i];
	}
	Result[0][Start.size()] = Start_T;

	double Temp, IntermediateValue;
	VecD NextRes(Start), TempRes(Start);
	NextRes.push_back(Start_T);

	for (double t = Start_T + dT; t <= Finish_T + (dT / 2); t += dT)
	{
		for (size_t i = 0; i < Start.size(); i++)
		{
			TempRes[i] += dT * (this->_Functions[i](TempRes, t - dT));
			NextRes[i] += dT * (this->_Functions[i](NextRes, t - dT) + this->_Functions[i](TempRes, t)) / 2;
		}
		NextRes[NextRes.size() - 1] = t;
		Result.push_back(NextRes);
	}

	return Result;
}

VecVecD SDE::SolveRungeKutt(VecD Start, double Start_T, double Finish_T, double dT, bool ClassicType)
{
	VecD C(4, 0), B(4, 0);
	VecVecD A(4, VecD(4, 0));

	if (ClassicType)
	{
		C[1] = C[2] = 0.5;
		A[1][0] = A[2][1] = 0.5;
		B[0] = B[3] = 1.0 / 6.0;
		B[1] = B[2] = 1.0 / 3.0;
	}
	else
	{
		C[1] = C[2] = 1.0 / 3.0;
		A[1][0] = 1.0 / 3.0;
		A[2][0] = -1.0 / 3.0;
		A[2][1] = A[3][0] = 1;
		A[3][1] = -1;
		B[0] = B[3] = 1.0 / 8.0;
		B[1] = B[2] = 3.0 / 8.0;
	}
	C[3] = 1;
	A[3][2] = 1;

	VecVecD Result(1, VecD(Start.size() + 1));

	for (size_t i = 0; i < Start.size(); ++i) {
		Result[0][i] = Start[i];
	}
	Result[0][Start.size()] = Start_T;

	VecD Ki(4, 0);
	VecD NextRes(Start), Temp;
	double PrevValue;
	NextRes.push_back(Start_T);

	for (double t = Start_T + dT; t <= Finish_T + (dT / 2); t += dT)
	{
		Temp = NextRes;
		for (size_t i = 0; i < Start.size(); i++)
		{
			Ki[0] = this->_Functions[i](Temp, t - dT);

			PrevValue = Temp[i];
			Temp[i] += A[1][0] * dT * Ki[0];
			Ki[1] = this->_Functions[i](Temp, t - dT + C[1] * dT);
			Temp[i] = PrevValue;

			Temp[i] += A[2][0] * dT * Ki[0] + A[2][1] * dT * Ki[1];
			Ki[2] = this->_Functions[i](Temp, t - dT + C[2] * dT);
			Temp[i] = PrevValue;

			Temp[i] += A[3][0] * dT * Ki[0] + A[3][1] * dT * Ki[1] + A[3][2] * dT * Ki[2];
			Ki[3] = this->_Functions[i](Temp, t - dT + C[3] * dT);
			Temp[i] = PrevValue;

			NextRes[i] += dT * (B[0] * Ki[0] + B[1] * Ki[1] + B[2] * Ki[2] + B[3] * Ki[3]);
		}
		NextRes[NextRes.size() - 1] = t;
		Result.push_back(NextRes);
	}

	return Result;
}

VecVecD SDE::SolveRichardson(VecD Start, double Start_T, double Finish_T, double dT, double Accuracy)
{
	VecVecD Result(1, VecD(Start.size() + 1));

	for (size_t i = 0; i < Start.size(); ++i) {
		Result[0][i] = Start[i];
	}
	Result[0][Start.size()] = Start_T;

	double p = 4;
	double k;
	double Norma;
	VecD N(Start), Res(Start);
	Res.push_back(Start_T);
	VecVecD Temp1, Temp2;

	for (double t = Start_T + dT; t < Finish_T + (dT / 2); t += dT)
	{
		k = 1;
		do
		{
			k *= 2;
			Norma = 0;

			Temp1 = this->SolveRungeKutt(N, t - dT, t, dT / k);
			Temp2 = this->SolveRungeKutt(N, t - dT, t, dT / k / 2);

			for (size_t i = 0; i < Temp1[0].size() - 1; ++i) {
				Norma += pow(Temp1[Temp1.size() - 1][i] - Temp2[Temp2.size() - 1][i], 2);
			}
			Norma = sqrt(Norma);

		} while (Norma / (pow(2, p) - 1) > Accuracy);

		for (size_t i = 0; i < Temp2[0].size() - 1; ++i)
		{
			N[i] = Temp2[Temp2.size() - 1][i];
			Res[i] = N[i];
		}
		Res[Res.size() - 1] = t;
		Result.push_back(Res);
	}

	return Result;
}
