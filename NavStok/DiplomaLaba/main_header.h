#pragma once
#include "my_matrix.h"
#include "functions.h"
#include <limits>
#include <fstream>


void Initialization_2_1(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xv, size_t& N,
	Matrix& Xgr, size_t& M,
	Matrix& Xgr1, size_t& M1,
	Matrix& Xgr2, size_t& M2,
	Matrix& X_body, Matrix& tau, size_t& M3, size_t& N_add, double& delta_alpha);


void InternalPoints_2_2(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xv, Matrix& dXv_dx, Matrix& dXv_dy,
	Matrix& dEv_dW1, Matrix& dEv_dW2, Matrix& dEv_dW3,
	Matrix& dEv_dt1, Matrix& dEv_dt2, Matrix& dEv_dt3,
	double& Ev, size_t N, size_t M);


void BoundaryPoints_2_3(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xgr, Matrix& dXgr_dx, Matrix& dXgr_dy,
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3,
	double& Egr, size_t N, size_t M, double& delta_alpha);

void BoundaryPoints_2_31(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3, 
	Matrix& Xgr, Matrix& dXgr_dx, Matrix& dXgr_dy,
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3, 
	double& Egr, size_t N, size_t M, double& delta_alpha);

void BoundaryPoints_2_32(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xgr, Matrix& dXgr_dx, Matrix& dXgr_dy,
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3,
	double& Egr, size_t N, size_t M, double& delta_alpha);

void BodyPoints(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xgr, Matrix& tau, Matrix& dXgr_dx, Matrix& dXgr_dy,
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3,
	double& Egr, size_t N, size_t M, double& delta_alpha);

void UpdateWeights_2_4(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& dE_dW1, Matrix& dE_dW2, Matrix& dE_dW3, 
	Matrix& dE_dt1, Matrix& dE_dt2, Matrix& dE_dt3, 
	Matrix& dE_dW1_k, Matrix& dE_dW2_k, Matrix& dE_dW3_k, 
	Matrix& dE_dt1_k, Matrix& dE_dt2_k, Matrix& dE_dt3_k,
	Matrix& dEv_dW1, Matrix& dEv_dW2, Matrix& dEv_dW3, 
	Matrix& dEv_dt1, Matrix& dEv_dt2, Matrix& dEv_dt3, 
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3, 
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3, 
	Matrix& dEgr1_dW1, Matrix& dEgr1_dW2, Matrix& dEgr1_dW3, 
	Matrix& dEgr1_dt1, Matrix& dEgr1_dt2, Matrix& dEgr1_dt3, 
	Matrix& dEgr2_dW1, Matrix& dEgr2_dW2, Matrix& dEgr2_dW3, 
	Matrix& dEgr2_dt1, Matrix& dEgr2_dt2, Matrix& dEgr2_dt3, 
	Matrix& dE_body_dW1, Matrix& dE_body_dW2, Matrix& dE_body_dW3, 
	Matrix& dE_body_dt1, Matrix& dE_body_dt2, Matrix& dE_body_dt3, 
	Matrix& delta_W1_k, Matrix& delta_W2_k, Matrix& delta_W3_k, 
	Vector& delta_t1_k, Vector& delta_t2_k, Vector& delta_t3_k, 
	std::vector<double>& vec_error, double Ev, double Egr, double Egr1, double Egr2, double E_body, size_t N, size_t M, size_t M1, size_t M2, size_t M3, double& delta_alpha);

void Verification(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	std::vector<double>& vec_error);