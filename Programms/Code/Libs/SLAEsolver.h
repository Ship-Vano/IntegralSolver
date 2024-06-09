#pragma once
#include<iostream>
#include<utility>
#include<vector> 
#include <iterator>
#include<fstream>
#include <sstream>
#include<string>
#include <iomanip>
#include <tuple>
using namespace std;

//класс системы
template<typename LT>
class System
{
public:
	System();
	System(vector<vector<LT>> A, vector<LT> b);
	System(string file_loc);
	void output();
	void sysoutput();
	void soloutput();
	void erroutput();
	void QRoutput();
	bool GaussianPartChoiceSolve(LT EPS);
	bool QRDecSolve(LT EPS);
	tuple<vector<vector<LT>>, vector<vector<LT>>> QR_decomposition(LT EPS);
	int n = 0;
	vector<vector<LT>> MatrixA;
	vector<LT> RightVect;
	vector<LT> SolutionX;
	vector<vector<LT>> MatrixR;
	vector<vector<LT>> MatrixQ;
protected:
	vector<vector<LT>> SystemMatrix;
	vector<LT> ApproxRightVect;
	vector<LT> ErrorVect;
	LT error_norm = 0;
protected:
	void file_init(string file_loc); // функция инициализации через файл
	int max_i(int index, LT EPS);
	void coef_reduce(int index, LT coef);
	void str_reduce(int current_ind, int upper_ind);
	//void matr_mult();
	void appr_vec_find();
};

/*КОНСТРУКТОРЫ*/

template<typename LT>
System<LT>::System()  // конструктор из файла
{

}

template<typename LT>
System<LT>::System(vector<vector<LT>> A, vector<LT> b)  // конструктор из файла
{
	MatrixA = A;
	RightVect = b;
	n = MatrixA.size();
	for (int i = 0; i < n; i++)
	{
		SystemMatrix.push_back(MatrixA[i]);
		SystemMatrix[i].push_back(RightVect[i]);
	}
	
}

template<typename LT>
System<LT>::System(string file_loc)  // конструктор из файла
{
	file_init(file_loc);
}

template<typename LT>
void System<LT>::file_init(string file_loc)
{
	ifstream file(file_loc); // открываем файл
	if (file.is_open()) // вызов метода is_open()
	{
		string tmp_line;
		//file >> tmp_line;
		//file >> tmp_line;

		//cout << "Тест №" << tmp_line << endl;
		while (getline(file, tmp_line))
		{
			istringstream ss(tmp_line);
			SystemMatrix.emplace_back(istream_iterator<double>(ss), \
				istream_iterator<double>());
		}
		SystemMatrix.erase(SystemMatrix.begin());
		//cout << "Успешно создали систему!!!" << endl;
		file.close();
		n = SystemMatrix[0].size() - 1;
		//cout << "Размерность матрицы А: " << n << endl;
		//заполнение матрицы системы и вектора правой части:
		for (int i = 0; i < n; i++)
		{
			MatrixA.push_back({});
			for (int j = 0; j <= n; j++)
			{
				if (j == n)
				{
					RightVect.push_back(SystemMatrix[i][j]);
					continue;
				}
				MatrixA[i].push_back(SystemMatrix[i][j]);
			}
		}
	}
	else
	{
		cout << "Файл не открыт!\n\n" << endl;
	}
}

/*ВЫВОДЫ*/

template<typename LT>
void System<LT>::output()
{
	string sep1 = "_____________________________________________";
	string sep2 = "---------------------";
	cout << sep1 << endl;
	cout << "\n Система выглядит следующим образом: " << endl;
	cout << sep2 << endl;
	cout << "Матрица А: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << fixed << setprecision(6) << setw(12) << setfill(' ')<< MatrixA[i][j] << "  ";
		}
		cout << endl;
	}
	cout << sep2 << endl;

	cout << "Вектор правой части: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << fixed << setprecision(6) << setw(12) << setfill(' ') << RightVect[i] << endl;
	}

	cout << sep1 << endl;
}

template<typename LT>
void System<LT>::sysoutput()
{
	cout << "Матрица системы: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n+1; j++)
		{
			cout << fixed << setprecision(6) << setw(12) << setfill(' ') << SystemMatrix[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "-----------" << endl;
}

template<typename LT>
void System<LT>::soloutput()
{
	cout << "Вектор-решение Х: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << fixed << setprecision(6) << setw(12) << setfill(' ') << SolutionX[i] << endl;
	}
}

template<typename LT>
void System<LT>::erroutput()
{
	appr_vec_find();
	cout << "-------------------- \n" << endl;
	cout << "Вектор, полученный подстановкой решения: \n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << fixed << setprecision(6) << setw(12) << setfill(' ') << ApproxRightVect[i] << endl;
	}
	cout << "-------------------- \n" << endl;
	cout << "Вектор-невязка:\n" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << scientific << setprecision(6) << setw(16) << setfill(' ') << ErrorVect[i] << endl;
	}
	cout << "-------------------- \n" << endl;
	cout << "Норма вектора невязки: \n" << endl;
	cout << scientific << setprecision(6) << setw(16) << error_norm << endl;
	

}

template<typename LT>
void System<LT>::QRoutput()
{
	cout << "-------------------- \n" << endl;
	cout << "Матрица Q:" << endl;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			cout << fixed << setprecision(6) << setw(12) << setfill(' ') << MatrixQ[i][j] << " ";
		cout << " " << endl;
	}
	cout << "-------------------- \n" << endl;
	cout << "Матрица R:" << endl;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			//cout << R[i][j] << " "; //чтобы показать, что ровно 0 задаётся под диагональю, а не приблизительно в рез-те расчёта
			cout << fixed << setprecision(6) << setw(12) << setfill(' ') << MatrixR[i][j] << " ";
		cout << " " << endl;
	}
	cout << "-------------------- \n" << endl;
}

/*МЕТОД ГАУССА*/

template<typename LT>
int System<LT>::max_i(int index, LT EPS)
{
	int max_i = index;
	LT max_ai = fabs(SystemMatrix[index][index]);
	for (int i = index; i < n; i++)
		if (fabs(SystemMatrix[i][index]) > max_ai)
		{
			max_i = i;
			max_ai = fabs(SystemMatrix[i][index]);
		}
	if (max_ai < EPS) return -1;
	return max_i;
}

template<typename LT>
void System<LT>::coef_reduce(int index, LT coef)
{
	for (int j = 0; j < n+1; j++)
		SystemMatrix[index][j] /= coef;
}

template<typename LT>
void System<LT>::str_reduce(int current_ind, int upper_ind)
{
	LT coef = SystemMatrix[current_ind][upper_ind];
	for (int i = 0; i < n + 1; i++)
		SystemMatrix[current_ind][i] -= coef * SystemMatrix[upper_ind][i];
}

template<typename LT>
bool System<LT>::GaussianPartChoiceSolve(LT EPS) 
{
	//ПРЯМОЙ ХОД:
	for (int i = 0; i < n; i++)
	{
		int max_pos = max_i(i, EPS);
		if (max_pos != -1)
		{
			SystemMatrix[i].swap(SystemMatrix[max_pos]);
			coef_reduce(i, SystemMatrix[i][i]);
			if (i == n) break;

			for (int j = i + 1; j < n; j++)
			{
					str_reduce(j, i);
			}

		}
		else {
			cout << "Матрица вырождена" << endl;
			return false;
		}
	}
	//ОБРАТНЫЙ ХОД:
	SolutionX.push_back(SystemMatrix[n - 1][n]);
	for (int i = n - 2; i > -1; i--)
	{
		LT xi = SystemMatrix[i][n];
		for (int j = i + 1; j < n; j++)
		{
			xi -= SystemMatrix[i][j] * SolutionX[n - j - 1];
		}
		SolutionX.push_back(xi);
	}
	reverse(SolutionX.begin(), SolutionX.end());
	return true;
}

/*НАХОЖДЕНИЕ ВЕКТОРА-НЕВЯЗКИ*/

template<typename LT>
void System<LT>::appr_vec_find()
{
	for (int i = 0; i < n; i++)
	{
		LT tmp_val = 0;
		for (int j = 0; j < n; j++)
		{
			tmp_val += MatrixA[i][j] * SolutionX[j];
		}
		ErrorVect.push_back(RightVect[i] - tmp_val);
		error_norm += (RightVect[i] - tmp_val)*(RightVect[i] - tmp_val);
		ApproxRightVect.push_back(tmp_val);
	}
	error_norm = sqrt(error_norm);
}

/*QR-факторизация*/

template<typename LT>
tuple<vector<vector<LT>>, vector<vector<LT>>> System<LT>::QR_decomposition(LT EPS)
{
	vector<vector<LT>> R(MatrixA);
	vector<vector<LT>> Q;
	//vector<LT> b(RightVect);
	LT c;
	LT s;
	for (int i = 0; i < n; ++i)
	{
		Q.push_back({});
		for (int j = 0; j < n; ++j)
		{
			LT elem = (j == i) ? 1 : 0;
			Q[i].push_back(elem);
		}
	}
	//предотвращаем случай нулевого ведущего элемента в строке:
	vector<int> swaps;
	bool flag_to_swap = false;
	for (int k = 0; k < n; ++k)
	{
		if (R[k][k] < EPS)
		{
			for (int i = 0; i < n; ++i)
			{
				int max_pos = max_i(i, EPS);
				if (max_pos != -1)
				{
					R[i].swap(R[max_pos]);
					swaps.push_back(i);
					swaps.push_back(max_pos);
					if (i == n) break;

				}
				else {
					cout << "Матрица вырождена" << endl;
					break;
				}
			}
		}
		flag_to_swap = true;
		break;
	}
	for (int i = 0; i < n - 1; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			c = R[i][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
			s = R[j][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
			vector<LT> R_i_str(R[i]);
			vector<LT> R_j_str(R[j]);
			vector<LT> Q_i_str(Q[i]);
			vector<LT> Q_j_str(Q[j]);
			//LT b_i = b[i];
			for (int k = 0; k < n; ++k)
			{
				R[i][k] = c * R_i_str[k] + s * R_j_str[k];
				R[j][k] = c * R_j_str[k] - s * R_i_str[k];

				Q[i][k] = c * Q_i_str[k] + s * Q_j_str[k];
				Q[j][k] = c * Q_j_str[k] - s * Q_i_str[k];
			}
			R[j][i] = 0;
			//b[i] = c * b_i + s * b[j];
			//b[j] = c * b[j] - s * b_i;
		}
	}
	LT temp_el;
	for (int i = 0; i < n; ++i)
		for (int j = i+1; j < n; ++ j)
			if (i != j)
			{
				temp_el = Q[i][j];
				Q[i][j] = Q[j][i];
				Q[j][i] = temp_el;
			}
	if (flag_to_swap)
	{
		for (int i = 0; i < swaps.size()/2 ; ++i)
		{
			Q[swaps[i]].swap(Q[swaps[i++]]);
		}
	}
	/*cout << "Вектор b*" << endl;
	for (int i = 0; i < n; ++i)
		cout << b[i] << endl;*/
	MatrixQ = Q;
	MatrixR = R;
	return { Q, R };
}

template<typename LT>
bool System<LT>::QRDecSolve(LT EPS)
{
	auto At = QR_decomposition(EPS);
	vector<vector<LT>> Q = get<0>(At);
	vector<vector<LT>> R = get<1>(At);
	MatrixQ = Q;
	MatrixR = R;
	//QRoutput();
	vector<LT> b_star;
	for (int i = 0; i < n; ++i)
	{
		LT tmp_val = 0;
		for (int j = 0; j < n; ++j)
		{
			tmp_val += RightVect[j] * Q[j][i];
		}
		b_star.push_back(tmp_val);
	}
	//cout << "Вектор b*" << endl;
	//for (int i = 0; i < n; ++i)
	//cout << b_star[i] << endl;
	
	//System<LT> SysQBB(Q, RightVect);
	//SysQBB.GaussianPartChoiceSolve(EPS); //показать, что b*, найденный разными способами один и тот же
	//SysQBB.soloutput();

	System<LT> SysRXB(R, b_star);
	if (SysRXB.GaussianPartChoiceSolve(EPS))
	{
		SolutionX = SysRXB.SolutionX;
		return true;
	}

	return false;
}