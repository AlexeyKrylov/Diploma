#include "my_matrix.h"
#include "my_vector.h"
#include "functions.h"
#include <fstream>
#include <string>


//2.0.1-3.2
	Matrix::Matrix()
	{
		nC = 0;
		nS = 0;
		data = std::vector<std::vector<double>>();
	}

	Matrix::Matrix(size_t d1, size_t d2) : nS(d1), nC(d2)
	{
		for (size_t i = 0; i < nS; ++i)
		{
			std::vector<double> vec;
			for (size_t j = 0; j < nC; ++j)
				vec.push_back(0);
			data.push_back(vec);
		}
	}

	Matrix::Matrix(std::vector<std::vector<double>> d) : data(d)
	{
		nS = data.size();
		nC = data[0].size();
	}


	Matrix Matrix::operator+(Matrix second) const
	{
		if (nS != second.nS || nC != second.nC)
		{
			std::cout << "Error +Matrix" << std::endl;
			exit(1);
		}

		Matrix res(nS, nC);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				res.data[i][j] = data[i][j] + second.data[i][j];
		return res;
	}

	Matrix Matrix::operator+(Vector vec) const
	{
		if (nS != vec.RetDim())
		{
			std::cout << "Error ++Vector" << std::endl;
			exit(1);
		}

		Matrix res(nS, nC);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				res.data[i][j] = data[i][j] + vec.RetVector()[i];
		return res;
	}

	Matrix Matrix::operator*(Matrix second) const
	{
		if (nC != second.nS)
		{
			std::cout << "Error Mul" << std::endl;
			exit(1);
		}

		Matrix res(nS, second.nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				for (size_t m = 0; m < nC; ++m)
					res.data[i][j] += data[i][m] * second.data[m][j];
		return res;
	}

	Matrix Matrix::operator*(double alpha) const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = data[i][j] * alpha;
		return res;
	}

	Matrix Matrix::T() const
	{
		Matrix res(nC, nS);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = data[j][i];
		return res;
	}

	Matrix Matrix::Mul_by_El(Matrix second) const
	{
		if (nC != second.nC || nS != second.nS)
		{
			std::cout << "Error Mul_by_el" << std::endl;
			exit(1);
		}

		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = data[i][j] * second.data[i][j];
		return res;
	}

	Matrix Matrix::Comprs_by_Col() const
	{
		Matrix v(nS, 1);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				v.data[i][0] += data[i][j];
		return v;
	}

	void Matrix::PrintMatrix() const
	{
		for (size_t i = 0; i < nS; ++i)
		{
			std::cout << "| ";
			for (size_t j = 0; j < nC; ++j)
				std::cout << data[i][j] << ", ";
			std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	Matrix Matrix::sigmoidM() const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = sigmoid(data[i][j]);
		return res;
	}

	Matrix Matrix::sigmoidM_p1() const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = sigmoid_p1(data[i][j]);
		return res;
	}

	Matrix Matrix::sigmoidM_p2() const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = sigmoid_p2(data[i][j]);
		return res;
	}

	Matrix Matrix::sigmoidM_p3() const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i][j] = sigmoid_p3(data[i][j]);
		return res;
	}

	void Matrix::RandomInit()
	{
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i][j] = (rand() % 200) / 100.0 - 1.0;
	}

	void Matrix::InitWithValue(double value)
	{
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i][j] = value;
	}


	Matrix Matrix::sign(double p, double m, double z) const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				if (data[i][j] > 0) res.data[i][j] = p;
				else if (data[i][j] < 0) res.data[i][j] = m;
				else res.data[i][j] = z;
		return res;
	}

	Vector Matrix::toVector() const
	{
		if (nC != 1)
		{
			std::cout << "Error toVector" << std::endl;
			exit(1);
		}
		std::vector<double> vec;
		for (size_t i = 0; i < nS; ++i) vec.push_back(data[i][0]);
		return Vector(vec);
	}

	size_t Matrix::Ret_nS() const { return nS; }

	size_t Matrix::Ret_nC() const { return nC; }

	double Matrix::MaxElementMatrix() const
	{
		double result = 0;
		for(size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
			{
				if (abs(data[i][j]) > result) result = abs(data[i][j]);
			}
		return result;
	}

	void Matrix::change11elem(double value)
	{
		data[0][0] = value;
	}

	std::vector<std::vector<double>> Matrix::RetData() const
	{
		return data;
	}

	void Matrix::from_file_txt(std::string filename)
	{
		std::ifstream fin(filename);

		if (fin.is_open())
		{
			for (size_t i = 0; i < nS; ++i)
				for (size_t j = 0; j < nC; ++j)
					fin >> data[i][j];
		}
		fin.close();
	}

	void Matrix::to_file_txt(std::string filename) const
	{
		std::ofstream fout(filename);

		if (fout.is_open())
		{
			for (size_t i = 0; i < nS; ++i) {
				for (size_t j = 0; j < nC; ++j)
					fout << data[i][j] << ' ';
				fout << std::endl;
			}
		}
		fout.close();
	}

	Matrix Matrix::get_string(size_t index) const
	{
		std::vector<std::vector<double>> vec;
		vec.push_back(data[index]);
		return Matrix(vec);
	}

	void Matrix::write_string(Matrix str, size_t index)
	{
		data[index] = str.data[0];
	}

	Matrix::~Matrix(){}
