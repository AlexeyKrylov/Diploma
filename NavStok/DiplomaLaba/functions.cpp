#include "my_matrix.h"
#include "functions.h"


double sigmoid(double x) {
	return 1. / (1. + exp(-x));
}

double sigmoid_p1(double x) {
	return sigmoid(x)*(1 - sigmoid(x));
}

double sigmoid_p2(double x) {
	return sigmoid_p1(x) * (1 - 2 * sigmoid(x));
}

double sigmoid_p3(double x) {
	return sigmoid_p2(x) * (1 - 2 * sigmoid(x)) - 2 * sigmoid_p1(x) * sigmoid_p1(x);
}


void changedelta(Matrix & dE_dW, Matrix & dE_dW_k, Matrix & delta_W_k) {
	Matrix sign(dE_dW_k.Ret_nS(), dE_dW_k.Ret_nC());
	sign = dE_dW.Mul_by_El(dE_dW_k).sign(eta_plus, eta_minus, 1);
	delta_W_k = delta_W_k.Mul_by_El(sign);
	dE_dW_k = dE_dW;
}

void changedelta(Matrix & dE_dt, Matrix & dE_dt_k, Vector & delta_t_k) {
	Matrix sign(dE_dt_k.Ret_nS(), dE_dt_k.Ret_nC());
	sign = dE_dt.Mul_by_El(dE_dt_k).sign(eta_plus, eta_minus, 1);
	delta_t_k = Matrix({ delta_t_k.RetVector() }).T().Mul_by_El(sign).toVector();
	dE_dt_k = dE_dt;
}

void parsing_weights(Matrix& W1, Vector& t1,
	Matrix& W2, Vector& t2,
	Matrix& W3, Vector& t3,
	std::string filename) {
	std::ifstream fin(filename, std::ios::binary);

	int N;
	char ch;
	int s1;
	int s2;
	int s3;
	int s4;

	fin.read((char*)&N, sizeof(int));
	std::cout << N << std::endl;
	fin.read((char*)&s1, sizeof(int));
	fin.read((char*)&s2, sizeof(int));
	fin.read((char*)&s3, sizeof(int));
	fin.read((char*)&s4, sizeof(int));
	std::cout << s1 << ' ' << s2 << ' ' << s3 << ' ' << s4 << std::endl;
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;
	fin.read((char*)&N, sizeof(int));
	std::cout << N << std::endl;
	fin.read((char*)&s1, sizeof(int));
	fin.read((char*)&s2, sizeof(int));
	fin.read((char*)&s3, sizeof(int));
	fin.read((char*)&s4, sizeof(int));
	std::cout << s1 << ' ' << s2 << ' ' << s3 << ' ' << s4 << std::endl;
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;

	std::vector<std::vector<double>> w1 = W1.RetData();

	for (size_t j = 0; j < w1[0].size(); ++j)
		for (size_t i = 0; i < w1.size(); ++i) {
			fin.read((char*)&w1[i][j], sizeof(double));
	}
	W1.PrintMatrix();
	W1.change_data(w1);
	W1.PrintMatrix();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;

	std::vector<double> t_1 = t1.RetVector();
	for (size_t i = 0; i < t1.RetDim(); ++i)
		fin.read((char*)&t_1[i], sizeof(double));
	t1.change_data(t_1);
	t1.PrintVector();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;


	std::vector<std::vector<double>> w2 = W2.RetData();
	for (size_t j = 0; j < w2[0].size(); ++j)
		for (size_t i = 0; i < w2.size(); ++i) {
			fin.read((char*)&w2[i][j], sizeof(double));
	}
	W2.change_data(w2);
	W2.PrintMatrix();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;

	std::vector<double> t_2 = t2.RetVector();
	for (size_t i = 0; i < t2.RetDim(); ++i)
		fin.read((char*)&t_2[i], sizeof(double));
	t2.change_data(t_2);
	t2.PrintVector();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;


	std::vector<std::vector<double>> w3 = W3.RetData();
	for (size_t j = 0; j < w3[0].size(); ++j)
		for (size_t i = 0; i < w3.size(); ++i) {
			fin.read((char*)&w3[i][j], sizeof(double));
	}
	W3.change_data(w3);
	W3.PrintMatrix();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;

	std::vector<double> t_3 = t3.RetVector();
	for (size_t i = 0; i < t3.RetDim(); ++i)
		fin.read((char*)&t_3[i], sizeof(double));
	t3.change_data(t_3);
	t3.PrintVector();
	fin.read((char*)&ch, sizeof(char));
	std::cout << ch << std::endl;
	fin.close();
}

void saving_weights(Matrix& W1, Vector& t1,
	Matrix& W2, Vector& t2,
	Matrix& W3, Vector& t3,
	std::string filename) {

	std::ofstream fout(filename, std::ios::binary);

	int N = 4;
	char ch = '#';
	int s1 = W1.Ret_nC();
	int s2 = W2.Ret_nS();
	int s3 = W2.Ret_nC();
	int s4 = W3.Ret_nS();

	fout.write((char*) &N, sizeof(int));

	fout.write((char*)&s1, sizeof(int));
	fout.write((char*)&s2, sizeof(int));
	fout.write((char*)&s3, sizeof(int));
	fout.write((char*)&s4, sizeof(int));
	fout.write((char*)&ch, sizeof(char));

	fout.write((char*)&N, sizeof(int));

	s1 = 1;
	s2 = 0;
	s3 = 0;
	s4 = 1;

	fout.write((char*)&s1, sizeof(int));
	fout.write((char*)&s2, sizeof(int));
	fout.write((char*)&s3, sizeof(int));
	fout.write((char*)&s4, sizeof(int));

	fout.write((char*)&ch, sizeof(char));

	std::vector<std::vector<double>> w1 = W1.RetData();
	for (size_t j = 0; j < w1[0].size(); ++j)
		for (size_t i = 0; i < w1.size(); ++i) {
			fout.write((char*)&w1[i][j], sizeof(double));
	}
	fout.write((char*)&ch, sizeof(char));

	std::vector<double> t_1 = t1.RetVector();
	for(size_t i = 0; i < t1.RetDim(); ++i)
		fout.write((char*)&t_1[i], sizeof(double));
	fout.write((char*)&ch, sizeof(char));


	std::vector<std::vector<double>> w2 = W2.RetData();
	for (size_t j = 0; j < w2[0].size(); ++j)
		for (size_t i = 0; i < w2.size(); ++i) {
			fout.write((char*)&w2[i][j], sizeof(double));
	}
	fout.write((char*)&ch, sizeof(char));

	std::vector<double> t_2 = t2.RetVector();
	for (size_t i = 0; i < t2.RetDim(); ++i)
		fout.write((char*)&t_2[i], sizeof(double));
	fout.write((char*)&ch, sizeof(char));


	std::vector<std::vector<double>> w3 = W3.RetData();
	for (size_t j = 0; j < w3[0].size(); ++j)
		for (size_t i = 0; i < w3.size(); ++i) {
			fout.write((char*)&w3[i][j], sizeof(double));
	}
	fout.write((char*)&ch, sizeof(char));

	std::vector<double> t_3 = t3.RetVector();
	for (size_t i = 0; i < t3.RetDim(); ++i)
		fout.write((char*)&t_3[i], sizeof(double));
	fout.write((char*)&ch, sizeof(char));
	fout.close();
}