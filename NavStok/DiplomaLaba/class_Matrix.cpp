#include "my_matrix.h"
#include "my_vector.h"
#include "functions.h"
#include <fstream>
#include <string>


//2.0.1-3.2
	Matrix::Matrix()
	{
		data = 0;
		nC = 0;
		nS = 0;
	}

	Matrix::Matrix(size_t d1, size_t d2) : nS(d1), nC(d2)
	{
		data = new double[d1 * d2];

		for (size_t i = 0; i < d1 * d2; ++i)
			data[i] = 0.0;
	}

	Matrix::Matrix(std::vector<std::vector<double>> d)
	{
		nS = d.size();
		nC = d[0].size();
		data = new double[nC * nS];

		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i + j * nS] = d[i][j];
	}

	Matrix::Matrix(Matrix const& a) : 
		nC(a.nC), nS(a.nS){
		data = new double[a.nS * a.nC];
		for (int i = 0; i < a.nC * a.nS; ++i) {
			data[i] = a.data[i];
		}
	}

	Matrix& Matrix::operator=(Matrix const& a) {
		if (this != &a)
			Matrix(a).swap(*this);
		return *this;
	}

	void Matrix::swap(Matrix& a) {
		size_t const t1 = nC;
		size_t const s1 = nS;
		nC = a.nC;
		nS = a.nS;
		a.nC = t1;
		a.nS = s1;

		double* const t2 = data;
		data = a.data;
		a.data = t2;
	}

	Matrix Matrix::operator+(Matrix second) const
	{
		if (nS != second.nS || nC != second.nC)
		{
			std::cout << "Error +Matrix" << std::endl;
			exit(1);
		}
		
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		const double* point_left = data;
		const double* point_right = second.data;
		double* current = res.data;
		for (; i < len_res; ++i)
		{
			*current = *point_left + *point_right;
			current++;
			point_left++;
			point_right++;
		}

		return res;
		/*
		size_t len_res = nS * nC;
		size_t i = 0;
		double* point_left = data;
		const double* point_right = second.data;
		for (; i < len_res; ++i)
		{
			*point_left += *point_right;
			point_left++;
			point_right++;
		}

		return *this;
		*/
	}

	Matrix Matrix::operator+(Vector vec) const
	{
		if (nS != vec.RetDim())
		{
			std::cout << "Error ++Vector" << std::endl;
			exit(1);
		}

		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		const double* point_left = data;
		double* current = res.data;
		std::vector<double> vvvvv = vec.RetVector();
		for (; i < len_res; ++i) {
			*current = *point_left + vvvvv[i % nS];
			current++;
			point_left++;
		}
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


		/*
		for (size_t i = 0; i < res1.nS; ++i)
		{
			for (size_t j = 0; j < res1.nC; ++j)
			{
				for (size_t m = 0; m < nC; ++m)
				{
					res1.data[i + res.nS * j] += data[i + nS * m] * second.data[m + second.nS * j];
				}
			}
		}
		*/
		/*
		const double* point_left = data;
		const double* point_right = second.data;
		double* current = res.data;

		size_t i = 0;
		size_t j = 0;
		size_t k = 0;
		for (; i < second.nC; ++i) {
			for (; j < nC; ++j) {
				for (; k < nS; ++k) {
					*current += *point_left * *point_right;
					point_left++;
					current++;
				}
				point_right++;
				current -= nS;
				k = 0;
			}
			j = 0;
			current += nS;
			point_left = data;
		}
		*/

		cblas_dgemm(CblasColMajor, // 1
			CblasNoTrans, // 2
			CblasNoTrans, // 3
			nS, // m
			second.nC, // n 
			nC, // k
			1.0, // alpha
			data, // A
			nS, // m
			second.data, // B
			second.nS, // k
			0.0, // beta
			res.data, // C
			nS); // m

		return res;
	}

	Matrix Matrix::operator*(double alpha) const
	{
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		double* current = res.data;
		double* value = data;
		for (; i < len_res; ++i) {
			*current = (*value) * alpha;
			current++;
			value++;
		}
		return res;
	}

	Matrix Matrix::T() const
	{
		Matrix res(nC, nS);
		for (size_t i = 0; i < res.nS; ++i)
			for (size_t j = 0; j < res.nC; ++j)
				res.data[i+res.nS*j] = data[j+nS*i];
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
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		const double* point_left = data;
		const double* point_right = second.data;
		double* current = res.data;
		for (; i < len_res; ++i)
		{
			*current = *point_left * *point_right;
			current++;
			point_left++;
			point_right++;
		}
		return res;
	}

	Matrix Matrix::Comprs_by_Col() const
	{
		Matrix v(nS, 1);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				v.data[i] += data[i+nS*j];
		return v;
	}

	void Matrix::PrintMatrix() const
	{
		for (size_t i = 0; i < nS; ++i)
		{
			std::cout << "| ";
			for (size_t j = 0; j < nC; ++j)
				std::cout << data[i+nS*j] << ", ";
			std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	Matrix Matrix::sigmoidM() const
	{
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		double* current = res.data;
		double* value = data;
		for (; i < len_res; ++i) {
			*current = sigmoid(*value);
			current++;
			value++;
		}
		return res;
	}

	Matrix Matrix::sigmoidM_p1() const
	{
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		double* current = res.data;
		double* value = data;
		for (; i < len_res; ++i) {
			*current = sigmoid_p1(*value);
			current++;
			value++;
		}
		return res;
	}

	Matrix Matrix::sigmoidM_p2() const
	{
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		double* current = res.data;
		double* value = data;
		for (; i < len_res; ++i){
			*current = sigmoid_p2(*value);
			current++;
			value++;
		}
		return res;
	}

	Matrix Matrix::sigmoidM_p3() const
	{
		Matrix res(nS, nC);
		size_t len_res = res.nS * res.nC;
		size_t i = 0;
		double* current = res.data;
		double* value = data;
		for (; i < len_res; ++i) {
			*current = sigmoid_p3(*value);
			current++;
			value++;
		}
		return res;
	}

	void Matrix::RandomInit()
	{
		double scale = 1.0;
		if (nS == L && nC == L || nS == 3 && nC == L) scale /= sqrt(L);

		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i+nS*j] = ((rand() % 200) / 100.0 - 1.0) * scale;
	}

	void Matrix::InitWithValue(double value)
	{
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i+nS*j] = value;
	}


	Matrix Matrix::sign(double p, double m, double z) const
	{
		Matrix res(nS, nC);
		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				if (data[i+nS*j] > 0) res.data[i+res.nS*j] = p;
				else if (data[i+nS*j] < 0) res.data[i+res.nS*j] = m;
				else res.data[i+res.nS*j] = z;
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
		for (size_t i = 0; i < nS; ++i) vec.push_back(data[i]);
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
				if (abs(data[i+nS*j]) > result) result = abs(data[i+nS*j]);
			}
		return result;
	}

	void Matrix::change11elem(double value)
	{
		data[0] = value;
	}

	std::vector<std::vector<double>> Matrix::RetData() const
	{
		std::vector<std::vector<double>> d(nS, std::vector<double>(nC));
		for (size_t i = 0; i < nS; ++i)
			for(size_t j = 0; j < nC; ++j)
				d[i][j] = data[i + nS * j];
		return d;
	}

	void Matrix::from_file_txt(std::string filename)
	{
		std::ifstream fin(filename);

		if (fin.is_open())
		{
			for (size_t i = 0; i < nS; ++i)
				for (size_t j = 0; j < nC; ++j)
					fin >> data[i+nS*j];
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
					fout << data[i+nS*j] << ' ';
				fout << std::endl;
			}
		}
		fout.close();
	}

	Matrix Matrix::get_string(size_t index) const
	{
		std::vector<std::vector<double>> vec;
		std::vector<double> v;
		for (size_t i = 0; i < nC; ++i)
			v.push_back(data[index + nS * i]);
		vec.push_back(v);
		return Matrix(vec);
	}

	void Matrix::write_string(Matrix str, size_t index)
	{
		for(size_t i = 0; i < nC; ++i)
			data[index + nS*i] = str.data[i];
	}

	void Matrix::change_data(std::vector<std::vector<double>> d)
	{
		nS = d.size();
		nC = d[0].size();	
		data = new double[nC * nS];

		for (size_t i = 0; i < nS; ++i)
			for (size_t j = 0; j < nC; ++j)
				data[i + j * nS] = d[i][j];
	}

	void Matrix::norm2()
	{
		size_t len_res = nS * nC;
		size_t i = 0;
		double* value = data;
		double res = 0;
		for (; i < len_res; ++i) {
			res += *value * *value;
			value++;
		}
		i = 0;
		res = sqrt(res);
		value = data;
		for (; i < len_res; ++i) {
			*value /= res;
			value++;
		}
	}

	void Matrix::del_large_values()
	{
		for (size_t i = 0; i < nS * nC; ++i)
			if (data[i] > 10) data[i] = 10;
			else if (data[i] < -10) data[i] = -10;
	}

	Matrix::~Matrix() {
		delete[] data;
	}