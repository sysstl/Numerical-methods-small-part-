#pragma once
#include "Mmath.h"

class SDE
{
public:

	SDE(std::vector<std::function<double(VecD&, double)>> Functions);
	VecVecD SolveEuler(VecD Start, double Start_T, double Finish_T, double dT);
	VecVecD SolveRungeKutt(VecD Start, double Start_T, double Finish_T, double dT, bool ClassicType = false);
	VecVecD SolveRichardson(VecD Start, double Start_T, double Finish_T, double dT, double Accuracy);

private:

	std::vector<std::function<double(VecD&, double)>> _Functions;

};
