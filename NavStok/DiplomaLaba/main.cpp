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

	double delta_alpha = lr_scheduler;

	Matrix Xv;
	size_t N = 0;
	
	Matrix Xgr;
	size_t M = 0;

	Matrix Xgr1;
	size_t M1 = 0;

	Matrix Xgr2;
	size_t M2 = 0;

	Matrix X_body;
	Matrix tau;
	size_t M3 = 0;
	size_t N_add = 0;

	Initialization_2_1(W1, W2, W3, t1, t2, t3, Xv, N, Xgr, M, Xgr1, M1, Xgr2, M2, X_body, tau, M3, N_add, delta_alpha);

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

	std::vector<double> zero_body(M3, 0);
	std::vector<double> one_body(M3, 1);

	Matrix dX_body_dx({ one_body, zero_body }); // (41)
	//std::cout << "dX_body_dx: " << std::endl;
	//dX_body_dx.PrintMatrix();

	Matrix dX_body_dy({ zero_body, one_body }); // (41)
	//std::cout << "dX_body_dy: " << std::endl;
	//dX_body_dy.PrintMatrix();

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
	Matrix dE_body_dW1(L, 2);
	Matrix dE_body_dt1(L, 1);
	Matrix dE_body_dW2(L, L);
	Matrix dE_body_dt2(L, 1);
	Matrix dE_body_dW3(3, L);
	Matrix dE_body_dt3(3, 1);
	double E_body = 0;
	double error1;
	double error0;
	
	for (size_t i = 0; i < K; ++i)
	{
		std::cout << "Number of iteration: " << i << std::endl;
		if (K >= checkpoint - 1) delta_alpha = 1.0;
		std::cout << "delta alpha = " << delta_alpha << std::endl;
		
			InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
				Xv, dXv_dx, dXv_dy,
				dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
				Ev, N, M);


			BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
				Xgr, dXgr_dx, dXgr_dy,
				dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
				Egr, N, M, delta_alpha);

			BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
				Xgr1, dXgr1_dx, dXgr1_dy,
				dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
				Egr1, N, M1, delta_alpha);

			BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
				Xgr2, dXgr2_dx, dXgr2_dy,
				dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
				Egr2, N, M2, delta_alpha);
		
		
			BodyPoints(W1, W2, W3, t1, t2, t3,
				X_body, tau, dX_body_dx, dX_body_dy,
				dE_body_dW1, dE_body_dW2, dE_body_dW3,
				dE_body_dt1, dE_body_dt2, dE_body_dt3,
				E_body, N, M3, delta_alpha);
		
		UpdateWeights_2_4(W1, W2, W3, t1, t2, t3, 
						  dE_dW1, dE_dW2, dE_dW3, dE_dt1, dE_dt2, dE_dt3,
						  dE_dW1_k, dE_dW2_k, dE_dW3_k, dE_dt1_k, dE_dt2_k, dE_dt3_k,
						  dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
						  dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
						  dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
						  dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
						  dE_body_dW1, dE_body_dW2, dE_body_dW3, dE_body_dt1, dE_body_dt2, dE_body_dt3,
						  delta_W1_k, delta_W2_k, delta_W3_k, delta_t1_k, delta_t2_k, delta_t3_k,
						  vec_error, Ev, Egr, Egr1, Egr2, E_body, N, M, M1, M2, M3, delta_alpha);
	
		std::cout << std::endl;
	}
	
	Verification(W1, W2, W3, t1, t2, t3, vec_error);
	
	std::cout << std::endl;
	
	
	/// /////////////////////////////////////////////////
	/*
	for (int index = 0; index < 5; ++index) {
		switch (index)
		{
		case 0:
			InternalPoints_2_2(W1, W2, W3, t1, t2, t3,
				Xv, dXv_dx, dXv_dy,
				dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
				Ev, N, M);
			error0 = Ev;
			std::cout << "dEv_dW1: \n" << std::endl;
			dEv_dW1.PrintMatrix();
			std::cout << "dEv_dW2: \n" << std::endl;
			dEv_dW2.PrintMatrix();
			std::cout << "dEv_dW3: \n" << std::endl;
			dEv_dW3.PrintMatrix();
			std::cout << "dEv_dt1: \n" << std::endl;
			dEv_dt1.PrintMatrix();
			std::cout << "dEv_dt2: \n" << std::endl;
			dEv_dt2.PrintMatrix();
			std::cout << "dEv_dt3: \n" << std::endl;
			dEv_dt3.PrintMatrix();
			break;
		case 1:
			BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3,
				Xgr, dXgr_dx, dXgr_dy,
				dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
				Egr, N, M, delta_alpha);
			error0 = Egr;
			std::cout << "dEgr_dW1: \n" << std::endl;
			dEgr_dW1.PrintMatrix();
			std::cout << "dEgr_dW2: \n" << std::endl;
			dEgr_dW2.PrintMatrix();
			std::cout << "dEgr_dW3: \n" << std::endl;
			dEgr_dW3.PrintMatrix();
			std::cout << "dEgr_dt1: \n" << std::endl;
			dEgr_dt1.PrintMatrix();
			std::cout << "dEgr_dt2: \n" << std::endl;
			dEgr_dt2.PrintMatrix();
			std::cout << "dEgr_dt3: \n" << std::endl;
			dEgr_dt3.PrintMatrix();
			break;
		case 2:
			BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3,
				Xgr1, dXgr1_dx, dXgr1_dy,
				dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
				Egr1, N, M1, delta_alpha);
			error0 = Egr1;
			std::cout << "dEgr1_dW1: \n" << std::endl;
			dEgr1_dW1.PrintMatrix();
			std::cout << "dEgr1_dW2: \n" << std::endl;
			dEgr1_dW2.PrintMatrix();
			std::cout << "dEgr1_dW3: \n" << std::endl;
			dEgr1_dW3.PrintMatrix();
			std::cout << "dEgr1_dt1: \n" << std::endl;
			dEgr1_dt1.PrintMatrix();
			std::cout << "dEgr1_dt2: \n" << std::endl;
			dEgr1_dt2.PrintMatrix();
			std::cout << "dEgr1_dt3: \n" << std::endl;
			dEgr1_dt3.PrintMatrix();
			break;
		case 3:
			BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3,
				Xgr2, dXgr2_dx, dXgr2_dy,
				dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
				Egr2, N, M2, delta_alpha);
			error0 = Egr2;
			std::cout << "dEgr2_dW1: \n" << std::endl;
			dEgr2_dW1.PrintMatrix();
			std::cout << "dEgr2_dW2: \n" << std::endl;
			dEgr2_dW2.PrintMatrix();
			std::cout << "dEgr2_dW3: \n" << std::endl;
			dEgr2_dW3.PrintMatrix();
			std::cout << "dEgr2_dt1: \n" << std::endl;
			dEgr2_dt1.PrintMatrix();
			std::cout << "dEgr2_dt2: \n" << std::endl;
			dEgr2_dt2.PrintMatrix();
			std::cout << "dEgr2_dt3: \n" << std::endl;
			dEgr2_dt3.PrintMatrix();
			break;
		case 4:
			BodyPoints(W1, W2, W3, t1, t2, t3,
				X_body, tau, dX_body_dx, dX_body_dy,
				dE_body_dW1, dE_body_dW2, dE_body_dW3,
				dE_body_dt1, dE_body_dt2, dE_body_dt3,
				E_body, N, M3, delta_alpha);
			error0 = E_body;
			std::cout << "dE_body_dW1: \n" << std::endl;
			dE_body_dW1.PrintMatrix();
			std::cout << "dE_body_dW2: \n" << std::endl;
			dE_body_dW2.PrintMatrix();
			std::cout << "dE_body_dW3: \n" << std::endl;
			dE_body_dW3.PrintMatrix();
			std::cout << "dE_body_dt1: \n" << std::endl;
			dE_body_dt1.PrintMatrix();
			std::cout << "dE_body_dt2: \n" << std::endl;
			dE_body_dt2.PrintMatrix();
			std::cout << "dE_body_dt3: \n" << std::endl;
			dE_body_dt3.PrintMatrix();
			break;
		}

		for (double pow = 0.1; pow > 10e-2; pow /= 10) {
			
			Matrix dddW1(L, 2);
			Matrix W1_new(L, 2);
			dddW1.InitWithValue(0.0);
			dddW1.change11elem(pow);
			W1_new = W1 + dddW1;
			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1_new, W2, W3, t1, t2, t3,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1_new, W2, W3, t1, t2, t3,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1_new, W2, W3, t1, t2, t3,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1_new, W2, W3, t1, t2, t3,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1_new, W2, W3, t1, t2, t3,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dW1:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";
			/// ////////////////////////////////////

			Matrix dddW2(L, L);
			Matrix W2_new(L, L);
			dddW2.InitWithValue(0.0);
			dddW2.change11elem(pow);
			W2_new = W2 + dddW2;

			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1, W2_new, W3, t1, t2, t3,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1, W2_new, W3, t1, t2, t3,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1, W2_new, W3, t1, t2, t3,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1, W2_new, W3, t1, t2, t3,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1, W2_new, W3, t1, t2, t3,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dW2:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";

			///
			Matrix dddW3(3, L);
			Matrix W3_new(3, L);
			dddW3.InitWithValue(0.0);
			dddW3.change11elem(pow);
			W3_new = W3 + dddW3;

			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1, W2, W3_new, t1, t2, t3,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1, W2, W3_new, t1, t2, t3,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1, W2, W3_new, t1, t2, t3,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1, W2, W3_new, t1, t2, t3,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1, W2, W3_new, t1, t2, t3,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dW3:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";
			/////////////////////////////////////////

			Matrix dddt1(L, 1);
			Vector t1_new(L);
			dddt1.InitWithValue(0.0);
			dddt1.change11elem(pow);
			t1_new = t1 + dddt1.toVector();

			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1, W2, W3, t1_new, t2, t3,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1, W2, W3, t1_new, t2, t3,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1, W2, W3, t1_new, t2, t3,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1, W2, W3, t1_new, t2, t3,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1, W2, W3, t1_new, t2, t3,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dt1:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";

			/////////////////////////////////////

			Matrix dddt2(L, 1);
			Vector t2_new(L);
			dddt2.InitWithValue(0.0);
			dddt2.change11elem(pow);
			t2_new = t2 + dddt2.toVector();
			
			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1, W2, W3, t1, t2_new, t3,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1, W2, W3, t1, t2_new, t3,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1, W2, W3, t1, t2_new, t3,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1, W2, W3, t1, t2_new, t3,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1, W2, W3, t1, t2_new, t3,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dt2:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";

			///////////////////////////////////////////////
			Matrix dddt3(3, 1);
			Vector t3_new(3);
			dddt3.InitWithValue(0.0);
			dddt3.change11elem(pow);
			t3_new = t3 + dddt3.toVector();

			switch (index)
			{
			case 0:
				InternalPoints_2_2(W1, W2, W3, t1, t2, t3_new,
					Xv, dXv_dx, dXv_dy,
					dEv_dW1, dEv_dW2, dEv_dW3, dEv_dt1, dEv_dt2, dEv_dt3,
					Ev, N, M);
				error1 = Ev;
				break;
			case 1:
				BoundaryPoints_2_3(W1, W2, W3, t1, t2, t3_new,
					Xgr, dXgr_dx, dXgr_dy,
					dEgr_dW1, dEgr_dW2, dEgr_dW3, dEgr_dt1, dEgr_dt2, dEgr_dt3,
					Egr, N, M, delta_alpha);
				error1 = Egr;
				break;
			case 2:
				BoundaryPoints_2_31(W1, W2, W3, t1, t2, t3_new,
					Xgr1, dXgr1_dx, dXgr1_dy,
					dEgr1_dW1, dEgr1_dW2, dEgr1_dW3, dEgr1_dt1, dEgr1_dt2, dEgr1_dt3,
					Egr1, N, M1, delta_alpha);
				error1 = Egr1;
				break;
			case 3:
				BoundaryPoints_2_32(W1, W2, W3, t1, t2, t3_new,
					Xgr2, dXgr2_dx, dXgr2_dy,
					dEgr2_dW1, dEgr2_dW2, dEgr2_dW3, dEgr2_dt1, dEgr2_dt2, dEgr2_dt3,
					Egr2, N, M2, delta_alpha);
				error1 = Egr2;
				break;
			case 4:
				BodyPoints(W1, W2, W3, t1, t2, t3_new,
					X_body, tau, dX_body_dx, dX_body_dy,
					dE_body_dW1, dE_body_dW2, dE_body_dW3,
					dE_body_dt1, dE_body_dt2, dE_body_dt3,
					E_body, N, M3, delta_alpha);
				error1 = E_body;
				break;
			}

			std::cout << "pow = " << pow << std::endl;
			std::cout << "dE_dt3:" << std::endl;
			std::cout << (error1 - error0) / pow << std::endl;
			std::cout << "\n\n\n";

		}
	}
	*/
	std::cout << "The End!" << std::endl;
	
	return 0;
}