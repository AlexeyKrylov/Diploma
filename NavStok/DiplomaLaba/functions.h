#pragma once
//2.0
//2.0.4
double sigmoid(double x);
double sigmoid_p1(double x);
double sigmoid_p2(double x);
double sigmoid_p3(double x);

void changedelta(Matrix& dE_dW, Matrix& dE_dW_k, Matrix& delta_W_k);

void changedelta(Matrix& dE_dt, Matrix& dE_dt_k, Vector& delta_t_k);