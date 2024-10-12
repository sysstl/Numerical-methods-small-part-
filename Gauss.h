#pragma once
#include <vector>
#include <string>
#include <iostream>

using namespace std;
class Gauss
{
private:
	int n;
	vector< vector<double> > A;
	vector<double> B;
	string err = "";
	bool solved = false;

public:
	bool errFlag = false;
	Gauss(int n_);
	Gauss(vector< vector<double> >& A_, vector<double>& B_);
	bool SolveSLAEGeneral();
	bool SolveSLAEMainElement();
	vector<double> ReturnAnswer() { return B; };
	friend ostream & operator <<(ostream&, const Gauss&);
	void PrintA(ostream& out);
	string GetError() { return err; };
};

class Progonka
{
private:
	int n;
	vector<double> a;
	vector<double> b;
	vector<double> c;
	vector<double> d;
	vector<double> e;
	vector<double> f;
	vector<double> alpha;
	vector<double> beta;
	vector<double> gama;
	vector<double> x;

public:

	Progonka(int n_);// дл проверки
	Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& f_);
	Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& d_, vector<double>& e_, vector<double>& f_);
	vector<double> SolveThreeDiagonal();
	vector<double> SolveFiveDiagonal();
};
