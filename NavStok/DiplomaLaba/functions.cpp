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