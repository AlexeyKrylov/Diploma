#pragma once
#include "libs_consts.h"
#include "my_vector.h"


//2.0.1-3.2
class Matrix {
private:
	std::vector<std::vector<double>> data;
	size_t nS;
	size_t nC;

public:
	friend void Verification(Matrix& W1, Matrix& W2, Matrix& W3,
		Vector& t1, Vector& t2, Vector& t3,
		std::vector<double>& vec_error);

	friend void BoundaryPoints_2_3(Matrix& W1, Matrix& W2, Matrix& W3,
		Vector& t1, Vector& t2, Vector& t3,
		Matrix& Xgr, Matrix& dXgr_dx, Matrix& dXgr_dy,
		Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
		Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3,
		double& Egr, size_t N, size_t M);

	Matrix();

	Matrix(size_t d1, size_t d2);

	Matrix(std::vector<std::vector<double>> d);

	Matrix operator+(Matrix second) const;

	Matrix operator+(Vector vec) const;

	Matrix operator*(Matrix second) const;

	Matrix operator*(double alpha) const;

	Matrix T() const;

	Matrix Mul_by_El(Matrix second) const;

	Matrix Comprs_by_Col() const;

	void PrintMatrix() const;

	Matrix sigmoidM() const;

	Matrix sigmoidM_p1() const;

	Matrix sigmoidM_p2() const;

	Matrix sigmoidM_p3() const;

	void RandomInit();

	void InitWithValue(double value);

	Matrix sign(double p, double m, double z) const;

	Vector toVector() const;

	double MaxElementMatrix() const;

	size_t Ret_nS() const;

	size_t Ret_nC() const;

	void change11elem(double value);

	std::vector<std::vector<double>> RetData() const;

	void from_file_txt(std::string filename);

	~Matrix();
};