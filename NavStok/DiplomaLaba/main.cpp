#include "main_header.h"


int main()
{	
	std::cout << "Start!" << std::endl;

	std::vector<double> vec_error;

	Matrix delta_W1_k(L, 2);
	Matrix delta_W2_k(L, L);
	Matrix delta_W3_k(3, L);
	Vector delta_t1_k(L);
	Vector delta_t2_k(L);
	Vector delta_t3_k(3);

	delta_W1_k.InitWithValue(0.0001);
	delta_W2_k.InitWithValue(0.0001);
	delta_W3_k.InitWithValue(0.0001);
	delta_t1_k.InitWithValue(0.0001);
	delta_t2_k.InitWithValue(0.0001);
	delta_t3_k.InitWithValue(0.0001);

	//2.1.2
	srand(4541);
	Matrix W1(L, 2);
	//W1.from_file_txt("W1.txt");
	//W1.PrintMatrix();
	Matrix W2(L, L);
	//W2.from_file_txt("W2.txt");
	//W2.PrintMatrix();
	Matrix W3(3, L);
	//W3.from_file_txt("W3.txt");
	//W3.PrintMatrix();
	Vector t1(L);
	//t1.from_file_txt("T1.txt");
	//t1.PrintVector();
	Vector t2(L);
	//t2.from_file_txt("T2.txt");
	//t2.PrintVector();
	Vector t3(3);
	//t3.from_file_txt("T3.txt");
	//t3.PrintVector();


	Matrix Xv;
	size_t N = 0;
	
	Matrix Xgr;
	size_t M = 0;

	Matrix Xgr1;
	size_t M1 = 0;

	Matrix Xgr2;
	size_t M2 = 0;

	Initialization_2_1(W1, W2, W3, t1, t2, t3, Xv, N, Xgr, M, Xgr1, M1, Xgr2, M2);

	Matrix dE_dW1(L, 2);
	Matrix dE_dt1(L, 1);
	Matrix dE_dW2(L, L);
	Matrix dE_dt2(L, 1);
	Matrix dE_dW3(3, L);
	Matrix dE_dt3(3, 1);

	Matrix dE_dW1_k(L, 2);
	Matrix dE_dt1_k(L, 1);
	Matrix dE_dW2_k(L, L);
	Matrix dE_dt2_k(L, 1);
	Matrix dE_dW3_k(3, L);
	Matrix dE_dt3_k(3, 1);

	// dXv_dx and dXv_dy
	std::vector<double> zero(N, 0);
	std::vector<double> one(N, 1);

	Matrix dXv_dx({ one, zero }); // (12)
	//std::cout << "dXv_dx: " << std::endl;
	//dXv_dx.PrintMatrix();

	Matrix dXv_dy({ zero, one }); // (17)
	//std::cout << "dXv_dy: " << std::endl;
	//dXv_dy.PrintMatrix();

	// dXgr_dx and dXgr_dy
	std::vector<double> zero1(M, 0);
	std::vector<double> one1(M, 1);

	Matrix dXgr_dx({one1, zero1}); // (41)
	//std::cout << "dXgr_dx: " << std::endl;
	//dXgr_dx.PrintMatrix();

	Matrix dXgr_dy({zero1, one1}); // (41)
	//std::cout << "dXgr_dy: " << std::endl;
	//dXgr_dy.PrintMatrix();

	// dXgr1_dx and dXgr1_dy
	std::vector<double> zero11(M1, 0);
	std::vector<double> one11(M1, 1);

	Matrix dXgr1_dx({ one11, zero11 }); // (41)
	//std::cout << "dXgr1_dx: " << std::endl;
	//dXgr1_dx.PrintMatrix();

	Matrix dXgr1_dy({ zero11, one11 }); // (41)
	//std::cout << "dXgr1_dy: " << std::endl;
	//dXgr1_dy.PrintMatrix();

	// dXgr2_dx and dXgr2_dy
	std::vector<double> zero12(M2, 0);
	std::vector<double> one12(M2, 1);

	Matrix dXgr2_dx({ one12, zero12 }); // (41)
	//std::cout << "dXgr2_dx: " << std::endl;
	//dXgr2_dx.PrintMatrix();

	Matrix dXgr2_dy({ zero12, one12 }); // (41)
	//std::cout << "dXgr2_dy: " << std::endl;
	//dXgr2_dy.PrintMatrix();


	std::cout << K << std::endl;

	Matrix dEv_dW1(L, 2);
	Matrix dEv_dt1(L, 1);
	Matrix dEv_dW2(L, L);
	Matrix dEv_dt2(L, 1);
	Matrix dEv_dW3(3, L);
	Matrix dEv_dt3(3, 1);
	double Ev = 0;
	Matrix dEgr_dW1(L, 2);
	Matrix dEgr_dt1(L, 1);
	Matrix dEgr_dW2(L, L);
	Matrix dEgr_dt2(L, 1);
	Matrix dEgr_dW3(3, L);
	Matrix dEgr_dt3(3, 1);
	double Egr = 0;
	Matrix dEgr1_dW1(L, 2);
	Matrix dEgr1_dt1(L, 1);
	Matrix dEgr1_dW2(L, L);
	Matrix dEgr1_dt2(L, 1);
	Matrix dEgr1_dW3(3, L);
	Matrix dEgr1_dt3(3, 1);
	double Egr1 = 0;
	Matrix dEgr2_dW1(L, 2);
	Matrix dEgr2_dt1(L, 1);
	Matrix dEgr2_dW2(L, L);
	Matrix dEgr2_dt2(L, 1);
	Matrix dEgr2_dW3(3, L);
	Matrix dEgr2_dt3(3, 1);
	double Egr2 = 0;


	for (size_t i = 0; i < K; ++i)
	{
		std::cout << "Number of iteration: " << i << std::endl;

		
		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);


		UpdateWeights_2_4(W1, W2, W3, t1, t2, t3, 
						  dE_dW1, dE_dW2, dE_dW3, dE_dt1, dE_dt2, dE_dt3,
						  dE_dW1_k, dE_dW2_k, dE_dW3_k, dE_dt1_k, dE_dt2_k, dE_dt3_k,
						  dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
						  dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
						  dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
						  dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
						  delta_W1_k, delta_W2_k, delta_W3_k, delta_t1_k, delta_t2_k, delta_t3_k,
						  vec_error, Ev, Egr, Egr1, Egr2, N, M, M1, M2);
	
		std::cout << std::endl;
	}

	Verification(W1, W2, W3, t1, t2, t3, vec_error);
	
	std::cout << std::endl;
	InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
		Xv, dXv_dx, dXv_dy,
		dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
		Ev, N, M);

	BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
		Xgr, dXgr_dx, dXgr_dy,
		dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
		Egr, N, M);

	BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
		Xgr1, dXgr1_dx, dXgr1_dy,
		dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
		Egr1, N, M1);

	BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
		Xgr2, dXgr2_dx, dXgr2_dy,
		dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
		Egr2, N, M2);

	/// /////////////////////////////////////////////////
	double pow = 0.001;
	double pow1 = 1000;

		Matrix dddW1(L, 2);
		dddW1.InitWithValue(0.0);
		dddW1.change11elem(pow);

		std::cout << "W1" << std::endl;
		W1.PrintMatrix();

		W1 = W1 + dddW1;

		double error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);
		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);
		double error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddoW1(L, 2);
		dddoW1.InitWithValue(0.0);
		dddoW1.change11elem(pow1);

		Matrix err1(L, 2);
		err1.InitWithValue(0.0);
		err1.change11elem(error777 - error666);

		dE_dW1 = dEv_dW1 + dEgr_dW1 + dEgr1_dW1 + dEgr2_dW1;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "err1.Mul_by_El(dddoW1)" << std::endl;
		err1.Mul_by_El(dddoW1).PrintMatrix();
		std::cout << "dE_dW1" << std::endl;
		dE_dW1.PrintMatrix();
		std::cout << "\n\n\n";
		/// ////////////////////////////////////

		Matrix dddW2(L, L);
		dddW2.InitWithValue(0.0);
		dddW2.change11elem(pow);

		std::cout << "W2" << std::endl;
		W2.PrintMatrix();

		W2 = W2 + dddW2;

		error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddoW2(L, L);
		dddoW2.InitWithValue(0.0);
		dddoW2.change11elem(pow1);

		Matrix err2(L, L);
		err2.InitWithValue(0.0);
		err2.change11elem(error777 - error666);

		dE_dW2 = dEv_dW2 + dEgr_dW2 + dEgr1_dW2 + dEgr2_dW2;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "err2.Mul_by_El(dddoW2)" << std::endl;
		err2.Mul_by_El(dddoW2).PrintMatrix();
		std::cout << "dE_dW2" << std::endl;
		dE_dW2.PrintMatrix();
		std::cout << "\n\n\n";
		///////////////////////////////////////////


		Matrix dddW3(3, L);
		dddW3.InitWithValue(0.0);
		dddW3.change11elem(pow);

		std::cout << "W3" << std::endl;
		W3.PrintMatrix();

		W3 = W3 + dddW3;

		error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);

		error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddoW3(3, L);
		dddoW3.InitWithValue(0.0);
		dddoW3.change11elem(pow1);

		Matrix err(3, L);
		err.InitWithValue(0.0);
		err.change11elem(error777 - error666);

		dE_dW3 = dEv_dW3 + dEgr_dW3 + dEgr1_dW3 + dEgr2_dW3;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "err.Mul_by_El(dddoW3)" << std::endl;
		err.Mul_by_El(dddoW3).PrintMatrix();
		std::cout << "dE_dW3" << std::endl;
		dE_dW3.PrintMatrix();
		std::cout << "\n\n\n";
		//////////////////////////////////////////

		Matrix dddt1(L, 1);
		dddt1.InitWithValue(0.0);
		dddt1.change11elem(pow);

		std::cout << "t1" << std::endl;
		t1.PrintVector();
		t1 = t1 + dddt1.toVector();
		error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);

		error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddot1(L, 1);
		dddot1.InitWithValue(0.0);
		dddot1.change11elem(pow1);

		Matrix errt1(L, 1);
		errt1.InitWithValue(0.0);
		errt1.change11elem(error777 - error666);

		dE_dt1 = dEv_dt1 + dEgr_dt1 + dEgr1_dt1 + dEgr2_dt1;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "errt1.Mul_by_El(dddot1)" << std::endl;
		errt1.Mul_by_El(dddot1).PrintMatrix();
		std::cout << "dE_dt1" << std::endl;
		dE_dt1.PrintMatrix();
		std::cout << "\n\n\n";

		/////////////////////////////////////

		Matrix dddt2(L, 1);
		dddt2.InitWithValue(0.0);
		dddt2.change11elem(pow);

		std::cout << "t2" << std::endl;
		t2.PrintVector();
		t2 = t2 + dddt2.toVector();
		error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);

		error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddot2(L, 1);
		dddot2.InitWithValue(0.0);
		dddot2.change11elem(pow1);

		Matrix errt(L, 1);
		errt.InitWithValue(0.0);
		errt.change11elem(error777 - error666);

		dE_dt2 = dEv_dt2 + dEgr_dt2 + dEgr1_dt2 + dEgr2_dt2;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "errt.Mul_by_El(dddot2)" << std::endl;
		errt.Mul_by_El(dddot2).PrintMatrix();
		std::cout << "dE_dt2" << std::endl;
		dE_dt2.PrintMatrix();
		std::cout << "\n\n\n";
		///////////////////////////////////////////////
		Matrix dddt3(3, 1);
		dddt3.InitWithValue(0.0);
		dddt3.change11elem(pow);

		std::cout << "t3" << std::endl;
		t3.PrintVector();
		t3 = t3 + dddt3.toVector();
		error666 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error666 = " << error666 << std::endl;

		InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
			Xv, dXv_dx, dXv_dy,
			dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
			Ev, N, M);

		BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
			Xgr, dXgr_dx, dXgr_dy,
			dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
			Egr, N, M);

		BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
			Xgr1, dXgr1_dx, dXgr1_dy,
			dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
			Egr1, N, M1);

		BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
			Xgr2, dXgr2_dx, dXgr2_dy,
			dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
			Egr2, N, M2);

		error777 = sqrt((Ev + Egr + Egr1 + Egr2) / (double(N) + double(M) + double(M1) + double(M2)));
		std::cout << "error777 = " << error777 << std::endl;

		Matrix dddot3(3, 1);
		dddot3.InitWithValue(0.0);
		dddot3.change11elem(pow1);

		Matrix errt3(3, 1);
		errt3.InitWithValue(0.0);
		errt3.change11elem(error777 - error666);

		dE_dt3 = dEv_dt3 + dEgr_dt3 + dEgr1_dt3 + dEgr2_dt3;
		std::cout << "pow = " << pow << std::endl;
		std::cout << "errt3.Mul_by_El(dddot3)" << std::endl;
		errt3.Mul_by_El(dddot3).PrintMatrix();
		std::cout << "dE_dt3" << std::endl;
		dE_dt3.PrintMatrix();
		std::cout << "\n\n\n";
		/////////////////////////////////////


	std::cout << "The End!" << std::endl;
	
	return 0;
}