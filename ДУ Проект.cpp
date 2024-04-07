// ДУ Пекарского.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <functional>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include <exception>
#include "Mmath.h"
#include "SDE.h"
#include "Gnuplot.h"

using namespace std;

int iter = 0;

struct Point
{
	Point(double x, double y) : x{ x }, y{ y } {}
	double x;
	double y;

	bool operator <(Point& A)
	{
		return y < A.y;
	}

};

template<typename K, typename V>
std::multimap<V, K> invertMap(std::map<K, V> const &map)
{
	std::multimap<V, K> multimap;

	for (auto const &pair : map) {
		multimap.insert(std::make_pair(pair.second, pair.first));
	}

	return multimap;
}

double Func1(VecD& start,double t) {
	return start[1];
	//return  start[1];
}

double Func2(VecD& start, double t) {
	//return -0.2*start[0];
	return -start[0] - 10.*(start[0] * start[0] - 1) * start[1];
}

Matrix SolveWaveEquation_2(std::vector<std::function<double(double, double)>>& Funcs, double alpha, double dx, double dt, double weight, std::pair<double, double> X_Range, std::pair<double, double> T_Range)
{
	size_t XN = static_cast<size_t>((X_Range.second - X_Range.first) / dx);
	size_t TN = static_cast<size_t>((T_Range.second - T_Range.first) / dt);

	Matrix Res(TN + 1, XN + 1);

	for (size_t i = 0; i <= TN; ++i)
	{
		Res[i][0] = Funcs[0](X_Range.first, i* dt);
		Res[i][XN] = Funcs[1](X_Range.second, i* dt);
	}

	for (size_t i = 1; i <= XN; ++i)
	{
		Res[0][i] = Funcs[2](i* dx, T_Range.first);
		Res[1][i] = dt * Funcs[3](i* dx, T_Range.first) + Res[0][i];
	}

	Matrix Mat(XN + 1, XN + 1);
	VecD FreeCol(XN + 1);
	VecD SlaeRes;

	double a = pow(alpha * dt / dx, 2);

	for (size_t j = 1; j < TN; j++)
	{
		Mat[0][0] = 1;
		Mat[XN][XN] = 1;
		FreeCol[0] = Funcs[0](X_Range.first, dt * j);
		FreeCol[XN] = Funcs[1](X_Range.second, dt * j);
		for (size_t i = 1; i < XN; ++i)
		{
			Mat[i][i - 1] = -weight * a;
			Mat[i][i] = (1 + 2 * weight * a);
			Mat[i][i + 1] = -weight * a;
			FreeCol[i] = 2 * Res[j][i] - Res[j - 1][i] +
				a * (1 - 2 * weight) * (Res[j][i + 1] - 2 * Res[j][i] + Res[j][i - 1]) +
				a * weight * (Res[j - 1][i + 1] - 2 * Res[j - 1][i] + Res[j - 1][i - 1]) +
				Funcs[4](dx* i, dt* j) * pow(dt, 2);
		}

		SlaeRes = SLAE(Mat, FreeCol).SolveProgon();

		for (size_t i = 1; i < XN; ++i) {
			Res[j + 1][i] = SlaeRes[i];
		}
	}

	return Res;
}

Matrix SolveHeatEquation_2(std::vector<std::function<double(double, double)>>& Funcs, double alpha, double dx, double dt, std::pair<double, double> X_Range, std::pair<double, double> T_Range)
{
	size_t XN = static_cast<size_t>((X_Range.second - X_Range.first) / dx);
	size_t TN = static_cast<size_t>((T_Range.second - T_Range.first) / dt);

	Matrix Res(TN + 1, XN + 1);

	for (size_t i = 0; i <= XN; ++i) {
		Res[0][i] = Funcs[2](i* dx, 0);
	}

	for (size_t i = 0; i <= TN; ++i)
	{
		Res[i][0] = Funcs[0](0, i* dt);
		Res[i][XN] = Funcs[1](1, i* dt);
	}

	Matrix Mat(XN + 1, XN + 1);
	VecD FreeCol(XN + 1);
	VecD SlaeRes;

	//VecD k(XN), l(XN);

	alpha = dt * pow(alpha / dx, 2);

	for (size_t j = 1; j < TN; ++j)
	{
		/*k[0] = 0;
		l[0] = Funcs[0](0, dt* j);
		for (size_t kN = 1; kN < XN; ++kN)
		{
			k[kN] = -alpha / (alpha * k[kN - 1] - (1 + 2 * alpha));
			l[kN] = -(Res[j - 1][kN] + dt * Funcs[3](kN* dx, (j - 1)* dt) + alpha * l[kN - 1]) / (alpha * k[kN - 1] - (1 + 2 * alpha));
		}*/

		Mat[0][0] = 1;
		Mat[XN][XN] = 1;
		FreeCol[0] = Funcs[0](0, dt * j);
		FreeCol[XN] = Funcs[1](1, dt* j);

		for (size_t i = 1; i < XN; ++i)
		{
			Mat[i][i - 1] = alpha;
			Mat[i][i] = -(1 + 2 * alpha);
			Mat[i][i + 1] = alpha;
			FreeCol[i] = -Res[j - 1][i] - dt * Funcs[3](dx* i, dt* (j - 1));
		}

		/*for (size_t i = 1; i < XN; ++i)
		{
			Mat[i][i] = 1;
			Mat[i][i + 1] = -k[i];
			FreeCol[i] = l[i];
		}*/

		SlaeRes = SLAE(Mat, FreeCol).SolveProgon();

		for (size_t i = 1; i < XN; ++i) {
			Res[j][i] = SlaeRes[i];
		}
	}

	/*for (size_t j = 1; j < TN; j++)
	{
		Mat[0][0] = 1;
		Mat[XN][XN] = 1;
		FreeCol[0] = Funcs[0](X_Range.first, dt* j);
		FreeCol[XN] = Funcs[1](X_Range.second, dt* j);
		for (size_t i = 1; i < XN; ++i)
		{
			Mat[i][i - 1] = -weight * a;
			Mat[i][i] = (1 + 2 * weight * a);
			Mat[i][i + 1] = -weight * a;
			FreeCol[i] = 2 * Res[j][i] - Res[j - 1][i] +
				a * (1 - 2 * weight) * (Res[j][i + 1] - 2 * Res[j][i] + Res[j][i - 1]) +
				a * weight * (Res[j - 1][i + 1] - 2 * Res[j - 1][i] + Res[j - 1][i - 1]) +
				Funcs[4](dx* i, dt* j) * pow(dt, 2);
		}

		SlaeRes = SLAE(Mat, FreeCol).SolveProgon();

		for (size_t i = 1; i < XN; ++i) {
			Res[j + 1][i] = SlaeRes[i];
		}
	}*/

	return Res;
}


Matrix SolveDEPD_Zeidel(std::vector<std::function<double(double)>> Funcs, size_t N, double accuracy)
{
	double h = 1.0 / (N - 1);
	//size_t iter = 0;

	Matrix Res(N, N), Prev;

	for (size_t i = 0; i < N; ++i)
	{
		Res[i][0] = Funcs[0](i* h + 1);  //U(x, 0) = 0
		Res[i][N - 1] = Funcs[1]((i* h + 1));  //U(x, 1) = atan(1/x)
	}

	for (size_t j = 0; j < N; ++j)
	{
		Res[0][j] = Funcs[2](j * h);  //U(1, y) = atan(y)
		Res[N - 1][j] = Funcs[3](j * h);  //U(2, y) = atan(y/2)
	}

	do
	{
		Prev = Res;
		//++iter;
		for (size_t i = 1; i < N - 1; ++i) {
			for (size_t j = 1; j < N - 1; ++j) {
				Res[i][j] = 0.25 * (Res[i + 1][j] + Res[i - 1][j] + Res[i][j + 1] + Res[i][j - 1]);
			}
		}
	} while ((Res - Prev).GetNorma() > accuracy);

	//for (int i = 0; i < N; ++i) {
	//	for (int j = 0; j < N; ++j) {
	//		cout << Res[i][j] << "   ";
	//	}
	//	cout << endl;
	//}

	return Res;
}

Matrix SolveDEPD_Relaxation(std::vector<std::function<double(double)>> Funcs, size_t Nx, size_t Ny, double Right_x, double Left_x, double Right_y, double Left_y, double accuracy, double omega, int& iter)
{
	double hx = (Right_x - Left_x) / (Nx - 1);
	double hy = (Right_y - Left_y) / (Ny - 1);
	//size_t iter = 0;

	Matrix Res(Nx, Ny);
	Matrix Res_new(Nx, Ny);
	Matrix TMP(Nx, Ny);

	for (int i = 0; i < Nx; ++i)
	{
		Res_new[i][0] = Funcs[0](i * hx + Left_x);  //U(x, 0) = 0
		Res_new[i][Ny - 1] = Funcs[1](i * hx + Left_x);  //U(x, 1) = atan(1/x)
	}
	for (int j = 0; j < Ny; ++j)
	{
		Res_new[0][j] = Funcs[2](j * hy + Left_y);  //U(1, y) = atan(y)
		Res_new[Nx - 1][j] = Funcs[3](j * hy + Left_y);  //U(2, y) = atan(y/2)
	}

	TMP = Res_new;

	double t = pow(hx / hy, 2);

	do
	{
		Res = Res_new;
		TMP = Res;
		++iter;
		for (int i = 1; i < Nx - 1; ++i)
		{
			for (int j = 1; j < Ny - 1; ++j)
			{
				Res_new[i][j] = (Res[i + 1][j] + Res[i - 1][j] + t * (Res[i][j + 1] + Res[i][j - 1])) / (2 + 2 * t);
				Res[i][j] = omega * Res_new[i][j] + (1 - omega) * Res[i][j];
			}
		}

	} while ((Res - TMP).GetNorma() > accuracy);

	//cout << iter << endl;

	//cout << endl;
	//TMP = (Res - Res_new);
	//for (int i = 0; i < TMP.GetRows(); i++)
	//{
	//	for (int j = 0; j < TMP.GetCols(); j++)
	//	{
	//		//mtr._Values;
	//		cout << Res[i][j] << "\t";
	//	}
	//	cout << "\n";
	//}

	//cout << endl;

	//for (int i = 0; i < TMP.GetRows(); i++)
	//{
	//	for (int j = 0; j < TMP.GetCols(); j++)
	//	{
	//		//mtr._Values;
	//		cout << TMP[i][j] << "\t";
	//	}
	//	cout << "\n";
	//}


	/*cout << endl << endl;
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			cout << Res[i][j] << "   ";
		}
		cout << endl;
	}*/

	return Res;
}

int main()
{
	VecD start = {1, 0}; // Нач услов x v
	//VecD start = { 1, 1 };
	vector<std::function<double(VecD&, double)>> Functions = {Func1, Func2};
	SDE sd(Functions);
	VecVecD Result = sd.SolveRichardson(start, 0., 50., 0.001, 0.0001);
	//VecVecD Result = sd.SolveRungeKutt(start, 0.1, 10., 0.01);

	ofstream fout("Y(x).txt");

	for (int i = 0; i < Result.size(); i++)
	{
		for (int j = 0; j < Result[i].size(); j++)
		//for (int j = Result[i].size(); j >= 0; j--)
		{
			//fout << Result[i][j] << "   ";
			//cout << Result[i][j] << "   ";
		}
		fout << Result[i][0] << "   "<< Result[i][1];
		//cout << endl;
		fout << endl;
	}
    std::cout << "Hello World!\n";


	vector<string> names_experim = { "Y(x).txt" };
	GnuPlot plt(1, names_experim);
	plt.SetParametrs2D(0, 1, 3, "y(x)", "time", "x,v");
	vector<string> legenda;
	legenda = { "x" };
	plt.ShowDataOnPlot2D(0, 1, legenda, "Y(x)", false);


	//  эллипт  2


	//vector<function<double(double)>> f{ p1, p2, p3, p4 };

	///*for (double i = 0.1; i > 1E-8; i = i / 2) {
	//	Zeydel(f, 100, i);
	//	cout << i << " " << iter << "\n";
	//	iter = 0;
	//}*/

	///* for (double i = 10; i < 200; i++) {
	//	 Zeydel(f, i, 0.0001);
	//	 cout << i << " " << iter << "\n";
	//	 iter = 0;
	// }*/

	//for (double i = 0.1; i < 2; i = i + 0.01) {
	//	Relaxation(f, 100, 0.001, i);
	//	cout << i << " " << iter << "\n";
	//	iter = 0;
	//}
	///*Relaxation(f, 50, 0.0001, 1);
	//cout << iter << "\n";*/

	//system("pause");
	

	// эллиптич

	int iter = 0;
	vector<function<double(double)>> Funcs
	{
		/*[](double x) {return log(x); },
		[](double x) {return 0.5*log(x*x + 1); },
		[](double x) {return 0.5*log(1 + x * x); },
		[](double x) {return 0.5* log(4 + x * x); }*/

		[](double x) {return 0.0; },
		[](double x) {return atan(1.0 / x); },
		[](double x) {return atan(x); },
		[](double x) {return atan(x / 2); }
	};

	cout << endl << endl;
	//SolveDEPD_Zeidel(Funcs, 10, 1e-5);
	//cout << endl << endl;
	//SolveDEPD_Relaxation(Funcs, 10, 10, 2, 1, 1, 0, 1e-5, 1.67, iter);
	system("pause");

	map<double, int> new_points;
	double x, y;

	//SolveDEPD_Zeidel(Funcs, 50, 1e-5);
	/////////////////////////////////////////////
	//for (int i = 10; i <= 50; i++)
	//{
	//	for (int j = 10; j <= 50; j++)
	//	{
	//		for (double omega = 0.999; omega < 1.99; omega += 0.001)
	//		{
	//			iter = 0;
	//			SolveDEPD_Relaxation(Funcs, i, j, 2, 1, 1, 0, 1e-8, omega, iter);
	//			//SolveDEPD_Relaxation(Funcs, 50, 50, 2, 1, 1, 0, 1e-10, omega, iter);
	//			new_points[omega] = iter;
	//			//cout << omega << "   " << new_points[omega] << endl;

	//		}

	//		multimap<int, double> multimap1 = invertMap(new_points);
	//		multimap<int, double> ::iterator it;
	//		it = multimap1.begin();
	//		cout << it->second << "   " << it->first << endl;
	//	}
	//}
	system("pause");

	iter = 0;
	Matrix TMP(10, 10);
	TMP = SolveDEPD_Relaxation(Funcs, 10, 10, 2, 1, 1, 0, 1e-5, 1.5, iter);

	for (int i = 0; i < TMP.GetRows(); i++)
	{
		for (int j = 0; j < TMP.GetCols(); j++)
		{
			//mtr._Values;
			cout << TMP[i][j] << "\t";
		}
		cout << "\n";
	}

	/*iter = 0;
	SolveDEPD_Relaxation(Funcs, 10, 10, 2, 1, 1, 0, 1e-5, 1, iter);
*/

	/*iter = 0;
	SolveDEPD_Relaxation(Funcs, 10, 10, 2, 1, 1, 0, 1e-5, 1, iter);
	cout << endl << iter << endl;
	iter = 0;
	SolveDEPD_Relaxation(Funcs, 10, 10, 2, 1, 1, 0, 1e-5, 0.99, iter);
	cout << endl << iter << endl;*/

	//std::vector<std::function<double(double)>> Funcs
	//{
	//	[](double x) {return 0.0; },
	//	[](double x) {return atan(1.0 / x); },
	//	[](double x) {return atan(x); },
	//	[](double x) {return atan(x / 2); }
	//};

	////SolveDEPD_Zeidel(Funcs, 50, 1e-5);
	//SolveDEPD_Relaxation(Funcs, 50, 1e-5);


	multimap<int, double> multimap1 = invertMap(new_points);
	multimap<int, double> ::iterator it;
	it = multimap1.begin();
	cout << it->second << "   " << it->first<<endl;

	/*for (auto const &pair : multimap1) {
		std::cout << '{' << pair.second << "," << pair.first << '}' << std::endl;
	}*/

	/*Matrix res = SolveDEPD_Relaxation(Funcs, 10, 1e-5, 1, 0);
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			cout << res[i][j] << "   ";
		}
		cout << endl;
	}*/


	std::vector<std::function<double(double, double)>> Funcs1
	{
		[](double x, double t) {return cos(t); },
		[](double x, double t) {return cos(t); },
		[](double x, double t) {return 1 + sin(PI * x); },
		[](double x, double t) {return -sin(x); }
	};

	double dx = 0.1;
	double dt = 0.001;
	double alpha = 1.0 / 2.0;

	Timer t;
	double sum = 0;
	Matrix res;
	for (size_t i = 0; i < 1; i++)
	{
		t.Start();
		res = SolveHeatEquation_2(Funcs1, alpha, dx, dt, { 0, 1 }, { 0, 1 });
		sum += t.Finish();
	}


	system("pause");
	return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
