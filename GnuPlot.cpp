#include "GnuPlot.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

GnuPlot::GnuPlot(int number_of_plots_)
{
	number_of_plots = number_of_plots_;
	file.reserve(number_of_plots);
	Matrix matr(0, 0);

	string text;
	FILE * gnuplotPipe_tmp;
	for (int i = 0; i < number_of_plots; ++i)
	{
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			text.push_back(str[j]);
		}
		string file_name = "gnuplot" + text + ".txt";
		text.clear();
		filename.push_back(file_name);
		file.emplace_back(ofstream{ file_name }); // ������� ������ ������ � ���� � ������ ����� �� ������ ���� ������ � �� �������� file.open(filename);
		gnuplotPipe_tmp = new FILE; // �������� ����� �� �������� ����� ���������� ���� ������ �� ��������
		gnuplotPipe.push_back(gnuplotPipe_tmp); // ������ ������
		gnuplotPipe[i] = _popen(GetGPPath().c_str(), "w"); // ��������� ������� ���� � �������� ��� � ��� �� ������
		P_xyzt.push_back(matr);
	}
}

GnuPlot::GnuPlot(int number_of_plots_, vector<string> filename_)
{
	number_of_plots = number_of_plots_;
	FILE * gnuplotPipe_tmp;
	for (int i = 0; i < number_of_plots; ++i)
	{
		gnuplotPipe_tmp = new FILE; // �������� ����� �� �������� ����� ���������� ���� ������ �� ��������
		gnuplotPipe.push_back(gnuplotPipe_tmp); // ������ ������
		gnuplotPipe[i] = _popen(GetGPPath().c_str(), "w"); // ��������� ������� ���� � �������� ��� � ��� �� ������
		filename = filename_;
	}
}

void GnuPlot::SetParametrs2D(int Current_number_plot, int number_of_lines, int width_line, string title_plot, string xlabel, string ylabel) // Current_number_plot=0,1,2, ....
{
	//gnuplotPipe[Current_number_plot] = _popen(GetGPPath().c_str(), "w"); // ��������� ������� ����
	//fprintf(gnuplotPipe[Current_number_plot], "set xrange [%lf:%lf]\n", 0., 50.); 
	//fprintf(gnuplotPipe[Current_number_plot], "set yrange [%lf:%lf]\n", 0., 50.);

	string tmp_width_line;
	int length = snprintf(NULL, 0, "%i", width_line);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%i", width_line);
	for (int j = 0; j < length; j++)
	{
		tmp_width_line.push_back(str[j]);
	}
	string color;

	for (int i = 1; i <= number_of_lines; i++) // number_of_lines - ����� ����� �� �����, 
	{
		string tmp_num_line;
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			tmp_num_line.push_back(str[j]); // ����� ����� ��� ��� ������
		}

		if (i == 1)
		{
			color = "red";
		}

		if (i == 2)
		{
			color = "yellow";
		}

		if (i == 3)
		{
			color = "orange";
		}

		if (i == 4)
		{
			color = "green";
		}

		if (i == 5)
		{
			color = "cyan";
		}

		if (i == 6)
		{
			color = "blue";
		}

		if (i == 7)
		{
			color = "violet";
		}

		if (i == 8)
		{
			color = "pink";
		}

		if (i == 9)
		{
			color = "dark-red";
		}

		if (i == 10)
		{
			color = "dark-yellow";
		}

		if (i == 11)
		{
			color = "dark-orange";
		}

		if (i == 12)
		{
			color = "dark-green";
		}

		if (i == 13)
		{
			color = "dark-cyan";
		}

		if (i == 14)
		{
			color = "dark-blue";
		}

		if (i == 15)
		{
			color = "dark-violet";
		}

		if (i == 16)
		{
			color = "dark-pink";
		}

		string str_str = "set linetype " + tmp_num_line + " lw " + tmp_width_line + " lc rgb " + "\"" + color + "\"\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str());
	}

	string title = "set title \"" + title_plot + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str()); //��������� ������� �������� ����� (����� Te Ti)
	title.clear();
	title = "set xlabel \"" + xlabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str()); // �������� ������� �������� ���� (���������� �����������)
	title.clear();
	title = "set ylabel \"" + ylabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());

	fprintf(gnuplotPipe[Current_number_plot], ("set term " + GetGPTerminal() + " position 0,0 size %zu,%zu\n").c_str(), GetTermnalWidth(), 415);
	fprintf(gnuplotPipe[Current_number_plot], "plot [][0:1] 2\n");
	printf("Press enter when window will appear");
	fflush(gnuplotPipe[Current_number_plot]); // ����������� ������ ����
}

void GnuPlot::SetParametrsOnPlotColor(int Current_number_plot, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y) // Current_number_plot=0,1,2, ....
{
	//gnuplotPipe[Current_number_plot] = _popen(GetGPPath().c_str(), "w"); // ��������� ������� ����
	fprintf(gnuplotPipe[Current_number_plot], "set xrange [%lf:%lf]\n", 0., right_bondary_x);
	fprintf(gnuplotPipe[Current_number_plot], "set yrange [%lf:%lf]\n", 0., top_bondary_y);

	fprintf(gnuplotPipe[Current_number_plot], "set palette defined ( 0 \"dark-violet\", 1 \"blue\", 2 \"cyan\", 3 \"green\", 4 \"yellow\", 5 \"orange\", 6 \"red\")\n");

	string title = "set title \"" + title_plot + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set xlabel \"" + xlabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set ylabel \"" + ylabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	fprintf(gnuplotPipe[Current_number_plot], ("set term " + GetGPTerminal() + " position 0,0 size %zu,%zu\n").c_str(), GetTermnalWidth(), 415);
	fprintf(gnuplotPipe[Current_number_plot], "plot [][0:1] 2\n");
	printf("Press enter when window will appear");
	fflush(gnuplotPipe[Current_number_plot]); // ����������� ������ ����
}

void GnuPlot::SetGridOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, int number_of_lines, Profile prof)
{
	if (prof == fun_on_x)
	{
		//Matrix mtr1(Nx, number_of_lines + 1);
		//P_xyzt[Current_number_plot] = mtr1; // ������� �������
		P_xyzt[Current_number_plot].Resize(Nx, number_of_lines + 1);// ���� out of range  ����� �� ��� �� 1 �������

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dx;
		}
	}

	if (prof == fun_on_y)
	{
		P_xyzt[Current_number_plot].Resize(Ny, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dy;
		}
	}
}

void GnuPlot::SetGridOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, int number_of_lines, Profile prof)
{
	if (prof == fun_on_x)
	{
		//Matrix mtr1(Nx, number_of_lines + 1);
		//P_xyzt[Current_number_plot] = mtr1; // ������� �������
		P_xyzt[Current_number_plot].Resize(Nx, number_of_lines + 1);// ���� out of range  ����� �� ��� �� 1 �������

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dx;
		}
	}

	if (prof == fun_on_y)
	{
		P_xyzt[Current_number_plot].Resize(Ny, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dy;
		}
	}

	if (prof == fun_on_z)
	{
		P_xyzt[Current_number_plot].Resize(Nz, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dz;
		}
	}
}

void GnuPlot::SetDataOnPlot1D(int Current_number_plot, int Nx, double dx, double** fun_dimensionless, double parametr_for_dimension, int moment_of_time, double dt, Profile prof)
{
	if (prof == fun_on_x) // current_number_of_line ����� �������� � 1 
	{ // ������ ����� ������� ����� ��������� �� �������� �� ������ ������ �.�. �� ��������� �������� ������ ���������� ������ 
		P_xyzt[Current_number_plot].Resize(Nx, moment_of_time + 1);// ���� out of range  ����� �� ��� �� 1 �������

		for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++)
		{
			P_xyzt[Current_number_plot][j][0] = j * dx;
			for (size_t k = 1; k < moment_of_time; k++) {
				P_xyzt[Current_number_plot][j][k] = parametr_for_dimension * fun_dimensionless[j][k];//????
			}
		}

		for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			for (int j = 0; j < moment_of_time + 1; j++)
			{
				file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
			}
			file[Current_number_plot] << endl;
		}
	}

	if (prof == fun_on_t)
	{
		P_xyzt[Current_number_plot].Resize(moment_of_time + 1, Nx);// ���� out of range  ����� �� ��� �� 1 �������

		for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++)
		{
			P_xyzt[Current_number_plot][j][0] = j * dt;
			for (size_t k = 1; k < Nx; k++) {
				P_xyzt[Current_number_plot][j][k] = parametr_for_dimension * fun_dimensionless[j][k];
			}
		}

		for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			for (int j = 0; j < Nx; j++)
			{
				file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
			}
			file[Current_number_plot] << endl;
		}
	}
}

void GnuPlot::SetDataOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, double** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double**> vec, Profile prof)
{
	if (prof == fun_on_x) // current_number_of_line ����� �������� � 1 
	{ // ������ ����� ������� ����� ��������� �� �������� �� ������ ������ �.�. �� ��������� �������� ������ ���������� ������ 
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[j][fixed_point_on_axis_y];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_y)
	{
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_y][j];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_t)
	{
		//  moment_of_time - �������� ��������
		file[Current_number_plot] << moment_of_time << "   ";
		for (int i = 0; i < vec.size(); i++)
		{
			file[Current_number_plot] << parametr_for_dimension * vec[i][fixed_point_on_axis_x][fixed_point_on_axis_y] << "   ";
		}
		file[Current_number_plot] << endl;
	}
}

void GnuPlot::SetDataOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int fixed_point_on_axis_z, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double***> vec, Profile prof)
{//��� ���������� �������� ����� (P �� x ���� ����) (P �� z ���� ����) (T �� �������)

	if (prof == fun_on_x) // current_number_of_line ����� �������� � 1 
	{ // ������ ����� ������� ����� ��������� �� �������� �� ������ ������ �.�. �� ��������� �������� ������ ���������� ������ 
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[j][fixed_point_on_axis_y][fixed_point_on_axis_z];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_y)
	{
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_y][j][fixed_point_on_axis_z];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_z)
	{
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_y][fixed_point_on_axis_y][j];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_t)
	{
		//  moment_of_time - �������� ��������
		file[Current_number_plot] << moment_of_time << "   ";
		for (int i = 0; i < vec.size(); i++)
		{
			file[Current_number_plot] << parametr_for_dimension * vec[i][fixed_point_on_axis_x][fixed_point_on_axis_y][fixed_point_on_axis_z] << "   ";
		}
		file[Current_number_plot] << endl;
	}
}

void GnuPlot::SetDataOnPlotColor3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis, Plane plane)
{
	//  ������� ������� ������� � ������� ������� � �������� ���� (����������: fun � SetDataOnPlotColor3D)!!!!!
	if (plane == xy)
	{
		for (size_t i = 0; i < Nx; i++)
			for (size_t j = 0; j < Ny; j++) {
				file[Current_number_plot] << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][fixed_point_on_axis] << endl;
			}
	}

	if (plane == xz)
	{
		for (size_t i = 0; i < Nx; i++)
			for (size_t j = 0; j < Nz; j++) {
				file[Current_number_plot] << j * dz << "   " << i * dx << "   " << parametr_for_dimension * fun_dimensionless[i][fixed_point_on_axis][j] << endl;
			}
	}

	if (plane == yz)
	{
		for (size_t i = 0; i < Ny; i++)
			for (size_t j = 0; j < Nz; j++) {
				file[Current_number_plot] << j * dz << "   " << i * dy << "   " << parametr_for_dimension * fun_dimensionless[fixed_point_on_axis][i][j] << endl;
			}
	}
}

void GnuPlot::SetDataOnPlotColor2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, double** fun, double parametr_for_dimension)
{
	for (size_t i = 0; i < Nx; i++)
		for (size_t j = 0; j < Ny; j++) {
			file[Current_number_plot] << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun[i][j] << endl;
		}
}

void GnuPlot::ShowDataOnPlot2D(int Current_number_plot, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_)
{
	if (png_)
	{
		//fprintf(gnuplotPipe[Current_number_plot], "set terminal png\n", filename[Current_number_plot].c_str
		fprintf(gnuplotPipe[Current_number_plot], "set terminal png\n");
		// name_of_file - �� ����� png
		string str_str_str = "set output \"" + name_of_file + ".png\" \n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str());
		str_str_str.clear();

		string str_str = "plot ";
		for (int i = 1; i <= number_of_lines; i++) // number_of_lines - ����� ����� �� �����, 
		{
			string tmp_num_line;

			int length = snprintf(NULL, 0, "%i", i + 1);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", i + 1);
			for (int j = 0; j < length; j++)
			{
				tmp_num_line.push_back(str[j]); // ����� �������
			}
			// string name_line - ������� ����� � ������� (�������� ����� ������ �������� ������ �� ���� � ��������� ������� ������� )
			string name_line = list_name_line[i - 1]; // ���� �� ����������, �� ���� ������ � ������ ���������� .c_str(), � ������ ��� C++
			str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
		}
		str_str.pop_back();
		str_str.pop_back();
		str_str += "\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str(), filename[Current_number_plot].c_str());
	}
	else
	{
		// ����� ����� ������ ���� �� ����� ����� ��� ������ 2-�� ������� ���� �� ����
		string str_str = "plot ";
		for (int i = 1; i <= number_of_lines; i++) // number_of_lines - ����� ����� �� �����, 
		{
			string tmp_num_line;

			int length = snprintf(NULL, 0, "%i", i + 1);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", i + 1);
			for (int j = 0; j < length; j++)
			{
				tmp_num_line.push_back(str[j]); // ����� �������
			}
			// string name_line - ������� ����� � ������� 
			string name_line = list_name_line[i - 1];
			str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
		}
		str_str.pop_back();
		str_str.pop_back();
		str_str += "\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str(), filename[Current_number_plot].c_str());
	}
	fflush(gnuplotPipe[Current_number_plot]);
}

void GnuPlot::ShowDataOnPlotColor(int Current_number_plot, string name_of_file, bool png_)
{
	if (png_)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set terminal png\n");
		string str_str_str = "set output \"" + name_of_file + ".png\" \n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str());
		str_str_str.clear();
		str_str_str = "plot '" + filename[Current_number_plot] + "' u 1:2:3 with image\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str(), filename[Current_number_plot].c_str());
	}
	else
	{
		string str_str_str = "plot '" + filename[Current_number_plot] + "' u 1:2:3 with image\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str(), filename[Current_number_plot].c_str());
	}
	fflush(gnuplotPipe[Current_number_plot]);
}

void GnuPlot::Close_and_open_files_for_replot(vector<int> Current_number_of_plots_)
{
	for (int i = 0; i < Current_number_of_plots_.size(); ++i)
	{
		file[Current_number_of_plots_[i]].close();
		file[Current_number_of_plots_[i]].open(filename[Current_number_of_plots_[i]]);
	}
}

void GnuPlot::Close_all_files_and_plots(int number_of_plots_)
{
	number_of_plots = number_of_plots_;
	for (int i = 0; i < number_of_plots; ++i)
	{
		file[i].close();
		fprintf(gnuplotPipe[i], "exit\n");
		_pclose(gnuplotPipe[i]);
	}
}

string GnuPlot::GetGPPath()
{
	return "\"C:\\Program Files (x86)\\gnuplot\\bin\\gnuplot\""; // ��� ����������� ���������
}

string GnuPlot::GetGPTerminal()
{
	return "wxt";
}

size_t GnuPlot::GetTermnalWidth() //������ ������ 
{
	return 650;
}