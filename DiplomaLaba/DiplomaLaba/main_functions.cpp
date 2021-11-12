#include "main_header.h"


void Initialization_2_1(Matrix& W1, Matrix& W2, Matrix& W3,
						Vector& t1, Vector& t2, Vector& t3,
						Matrix& Xv, size_t& N,
						Matrix& Xgr, size_t& M)
{
	std::cout << std::endl;
	//std::cout << "Start of Initialization_2_1" << std::endl;

	W1.RandomInit();
	//std::cout << "W1: " << std::endl;
	//W1.PrintMatrix();

	W2.RandomInit();
	//std::cout << "W2: " << std::endl;
	//W2.PrintMatrix();

	W3.RandomInit();
	//std::cout << "W3: " << std::endl;
	//W3.PrintMatrix();

	t1.RandomInit();
	//std::cout << "t1: " << std::endl;
	//t1.PrintVector();

	t2.RandomInit();
	//std::cout << "t2: " << std::endl;
	//t2.PrintVector();

	t3.RandomInit();
	//std::cout << "t3: " << std::endl;
	//t3.PrintVector();


	//2.1.3
	std::vector<std::vector<std::pair<double, double>>> grid;
	for (double i = -1.0; i <= 1.0 + std::numeric_limits<double>::epsilon(); i += lambda)
	{
		std::vector<std::pair<double, double>> vec;
		for (double j = -1.0; j <= 1.0 + std::numeric_limits<double>::epsilon(); j += lambda) vec.push_back(std::make_pair(i, j));
		grid.push_back(vec);
	}

	// x^2 + y^2 <= 1

	for (size_t i = 0; i < grid.size(); ++i)
		for (size_t j = 0; j < grid[0].size(); ++j) {
			if (grid[i][j].first * grid[i][j].first + grid[i][j].second * grid[i][j].second > 1.0 + std::numeric_limits<double>::epsilon()) grid[i][j] = std::make_pair(777, 777);
			else ++N;
		}

	std::cout << "N = " << N << std::endl;

	std::vector<std::vector<double>> vecxv;
	std::vector<double> x;
	std::vector<double> y;

	for (size_t i = 0; i < grid.size(); ++i)
		for (size_t j = 0; j < grid[0].size(); ++j)
			if (grid[i][j].first != 777)
			{
				x.push_back(grid[i][j].first);
				y.push_back(grid[i][j].second);
			}

	vecxv = { x, y };
	Xv = Matrix(vecxv);
	//std::cout << "Xv: " << std::endl;
	//Xv.PrintMatrix();

	//2.1.4
	std::vector<std::vector<double>> vecxgr;
	std::vector<double> x_phi;
	std::vector<double> y_phi;

	for (double phi = 0; phi < 2 * M_PI + std::numeric_limits<double>::epsilon(); phi += lambda)
	{
		x_phi.push_back(cos(phi));
		y_phi.push_back(sin(phi));
	}

	M = x_phi.size();
	std::cout << "M = " << M << std::endl;

	vecxgr = { x_phi, y_phi };
	Xgr = Matrix(vecxgr); // (40)

	//std::cout << "Xgr: " << std::endl;
	//Xgr.PrintMatrix();	

	//std::cout << "End of Initialization_2_1" << std::endl;
}
	
void InternalPoints_2_2(Matrix& W1, Matrix& W2, Matrix& W3,
						Vector& t1, Vector& t2, Vector& t3,
						Matrix& Xv, Matrix& dXv_dx, Matrix& dXv_dy,
						Matrix& dEv_dW1, Matrix& dEv_dW2, Matrix& dEv_dW3,
						Matrix& dEv_dt1, Matrix& dEv_dt2, Matrix& dEv_dt3,
						double& Ev, size_t N, size_t M)
{
	//std::cout << "Start of InternalPoints_2_2" << std::endl;
	//2.2 Internal points
	
	//2.1.5
	Matrix A(L, N);
	Matrix B(L, N);
	Matrix u(1, N);
	Matrix Ax(L, N);
	Matrix Ax1(L, N);
	Matrix Bx(L, N);
	Matrix Bx1(L, N);
	Matrix ux(1, N);
	Matrix Axx1(L, N);
	Matrix Bxx(L, N);
	Matrix Bxx1(L, N);
	Matrix uxx(1, N);
	Matrix Ax2(L, N);
	Matrix Bx2(L, N);
	Matrix Axx2(L, N);
	Matrix Bxx2(L, N);
	Matrix Ay(L, N);
	Matrix Ay1(L, N);
	Matrix By(L, N);
	Matrix By1(L, N);
	Matrix uy(1, N);
	Matrix Ayy1(L, N);
	Matrix Byy(L, N);
	Matrix Byy1(L, N);
	Matrix uyy(1, N);
	Matrix Ay2(L, N);
	Matrix By2(L, N);
	Matrix Ayy2(L, N);
	Matrix Byy2(L, N);

	std::vector<double> tmp1(1, -1.0);

	Matrix dEv_duxx(1, N);
	Matrix dEv_duyy(1, N);
	Matrix dEv_dux(1, N); //        (25)
	Matrix dEv_duy(1, N); //        (25)
	Matrix dEv_du(1, N); //         (25)

	Matrix dEv_dB(L, N);
	Matrix dEv_dBx(L, N);
	Matrix dEv_dBxx(L, N);
	Matrix dEv_dBy(L, N);
	Matrix dEv_dByy(L, N);

	Matrix dEv_dA(L, N);
	Matrix dEv_dAx(L, N);
	Matrix dEv_dAy(L, N);

//2.2.1
	A = W1 * Xv + t1; //           (8)
	//std::cout << "A: " << std::endl;
	//A.PrintMatrix();

	B = W2 * A.sigmoidM() + t2;//  (9)
	//std::cout << "B: " << std::endl;
	//B.PrintMatrix();

	u = W3 * B.sigmoidM() + t3;//  (10)
	//std::cout << "u: " << std::endl;
	//u.PrintMatrix();

//2.2.2
	Ax = W1 * dXv_dx; //           (11)
	//std::cout << "Ax: " << std::endl;
	//Ax.PrintMatrix();
	Ax1 = A.sigmoidM_p1().Mul_by_El(Ax);
	//std::cout << "Ax1: " << std::endl;
	//Ax1.PrintMatrix();
	Bx = W2 * Ax1; //              (13)
	//std::cout << "Bx: " << std::endl;
	//Bx.PrintMatrix();
	Bx1 = B.sigmoidM_p1().Mul_by_El(Bx);
	//std::cout << "Bx1: " << std::endl;
	//Bx1.PrintMatrix();
	ux = W3 * Bx1; //              (15)
	//std::cout << "ux: " << std::endl;
	//ux.PrintMatrix();
	Axx1 = A.sigmoidM_p2().Mul_by_El(Ax).Mul_by_El(Ax);
	//std::cout << "Axx1: " << std::endl;
	//Axx1.PrintMatrix();
	Bxx = W2 * Axx1; //            (18)
	//std::cout << "Bxx: " << std::endl;
	//Bxx.PrintMatrix();
	Bxx1 = B.sigmoidM_p2().Mul_by_El(Bx).Mul_by_El(Bx) + B.sigmoidM_p1().Mul_by_El(Bxx);
	//std::cout << "Bxx1: " << std::endl;
	//Bxx1.PrintMatrix();
	uxx = W3 * Bxx1; //            (20)
	//std::cout << "uxx: " << std::endl;
	//uxx.PrintMatrix();
	Ax2 = A.sigmoidM_p2().Mul_by_El(Ax); // (14)
	//std::cout << "Ax2: " << std::endl;
	//Ax2.PrintMatrix();
	Bx2 = B.sigmoidM_p2().Mul_by_El(Bx);//  (16)
	//std::cout << "Bx2: " << std::endl;
	//Bx2.PrintMatrix();
	Axx2 = A.sigmoidM_p3().Mul_by_El(Ax).Mul_by_El(Ax); // (19)
	//std::cout << "Axx2: " << std::endl;
	//Axx2.PrintMatrix();
	Bxx2 = B.sigmoidM_p3().Mul_by_El(Bx).Mul_by_El(Bx) + B.sigmoidM_p2().Mul_by_El(Bxx); // (21)
	//std::cout << "Bxx2: " << std::endl;
	//Bxx2.PrintMatrix();

//2.2.3
	Ay = W1 * dXv_dy; //           (11)
	//std::cout << "Ay: " << std::endl;
	//Ay.PrintMatrix();
	Ay1 = A.sigmoidM_p1().Mul_by_El(Ay);
	//std::cout << "Ay1: " << std::endl;
	//Ay1.PrintMatrix();
	By = W2 * Ay1; //              (13)
	//std::cout << "By: " << std::endl;
	//By.PrintMatrix();
	By1 = B.sigmoidM_p1().Mul_by_El(By);
	//std::cout << "By1: " << std::endl;
	//By1.PrintMatrix();
	uy = W3 * By1; //              (15)
	//std::cout << "uy: " << std::endl;
	//uy.PrintMatrix();
	Ayy1 = A.sigmoidM_p2().Mul_by_El(Ay).Mul_by_El(Ay);
	//std::cout << "Ayy1: " << std::endl;
	//Ayy1.PrintMatrix();
	Byy = W2 * Ayy1; //            (18)
	//std::cout << "Byy: " << std::endl;
	//Byy.PrintMatrix();
	Byy1 = B.sigmoidM_p2().Mul_by_El(By).Mul_by_El(By) + B.sigmoidM_p1().Mul_by_El(Byy);
	//std::cout << "Byy1: " << std::endl;
	//Byy1.PrintMatrix();
	uyy = W3 * Byy1; //            (20)
	//std::cout << "uyy: " << std::endl;
	//uyy.PrintMatrix();
	Ay2 = A.sigmoidM_p2().Mul_by_El(Ay); // (14)
	//std::cout << "Ay2: " << std::endl;
	//Ay2.PrintMatrix();
	By2 = B.sigmoidM_p2().Mul_by_El(By);//  (16)
	//std::cout << "By2: " << std::endl;
	//By2.PrintMatrix();
	Ayy2 = A.sigmoidM_p3().Mul_by_El(Ay).Mul_by_El(Ay); // (19)
	//std::cout << "Ayy2: " << std::endl;
	//Ayy2.PrintMatrix();
	Byy2 = B.sigmoidM_p3().Mul_by_El(By).Mul_by_El(By) + B.sigmoidM_p2().Mul_by_El(Byy); // (21)
	//std::cout << "Byy2: " << std::endl;
	//Byy2.PrintMatrix();
	//std::cout << "2.2.2" << std::endl;
//2.2.4
	Ev = (uxx + uyy + tmp1).Mul_by_El(uxx + uyy + tmp1).Comprs_by_Col().RetData()[0][0] * 0.5; // (22)
	std::cout << "Ev = " << Ev << std::endl;

//2.2.5
	dEv_duxx = uxx + uyy + tmp1; // (23)
	dEv_duyy = uxx + uyy + tmp1; // (24)

//2.2.6
	dEv_dW3 = dEv_du * B.sigmoidM().T() + dEv_dux * Bx1.T() + dEv_duy * By1.T() + dEv_duxx * Bxx1.T() + dEv_duyy * Byy1.T(); // (26)
	//std::cout << "dEv_dW3: " << std::endl;
	//dEv_dW3.PrintMatrix();

	dEv_dt3 = dEv_du.Comprs_by_Col(); // (27)
	//std::cout << "dEv_dt3: " << std::endl;
	//dEv_dt3.PrintMatrix();

//2.2.7
	dEv_dB = (W3.T() * dEv_du).Mul_by_El(B.sigmoidM_p1()) +
		(W3.T() * dEv_dux).Mul_by_El(Bx2) +
		(W3.T() * dEv_duy).Mul_by_El(By2) +
		(W3.T() * dEv_duxx).Mul_by_El(Bxx2) +
		(W3.T() * dEv_duyy).Mul_by_El(Byy2); // (29)
	dEv_dBx = (W3.T() * dEv_dux).Mul_by_El(B.sigmoidM_p1()) + (W3.T() * dEv_duxx).Mul_by_El(Bx2) * 2.0; // (30)
	dEv_dBxx = (W3.T() * dEv_duxx).Mul_by_El(B.sigmoidM_p1()); // (31)
	dEv_dBy = (W3.T() * dEv_duy).Mul_by_El(B.sigmoidM_p1()) + (W3.T() * dEv_duyy).Mul_by_El(By2) * 2.0; // (32)
	dEv_dByy = (W3.T() * dEv_duyy).Mul_by_El(B.sigmoidM_p1()); // (33)

//2.2.8
	dEv_dW2 = dEv_dB * A.sigmoidM().T() +
		dEv_dBx * Ax1.T() +
		dEv_dBy * Ay1.T() +
		dEv_dBxx * Axx1.T() +
		dEv_dByy * Ayy1.T(); // (28)
	//std::cout << "dEv_dW2: " << std::endl;
	//dEv_dW2.PrintMatrix();

	dEv_dt2 = dEv_dB.Comprs_by_Col(); // (34)
	//std::cout << "dEv_dt2: " << std::endl;
	//dEv_dt2.PrintMatrix();

//2.2.9
	dEv_dA = (W2.T() * dEv_dB).Mul_by_El(A.sigmoidM_p1()) +
		(W2.T() * dEv_dBx).Mul_by_El(Ax2) +
		(W2.T() * dEv_dBy).Mul_by_El(Ay2) +
		(W2.T() * dEv_dBxx).Mul_by_El(Axx2) +
		(W2.T() * dEv_dByy).Mul_by_El(Ayy2); // (35)
	dEv_dAx = (W2.T() * dEv_dBx).Mul_by_El(A.sigmoidM_p1()) + (W2.T() * dEv_dBxx).Mul_by_El(Ax2) * 2.0; // (36)
	dEv_dAy = (W2.T() * dEv_dBy).Mul_by_El(A.sigmoidM_p1()) + (W2.T() * dEv_dByy).Mul_by_El(Ay2) * 2.0; // (37)

//2.2.10
	dEv_dW1 = dEv_dA * Xv.T() + dEv_dAx * dXv_dx.T() + dEv_dAy * dXv_dy.T(); // (38)
	//std::cout << "dEv_dW1: " << std::endl;
	//dEv_dW1.PrintMatrix();

	dEv_dt1 = dEv_dA.Comprs_by_Col(); // (39)
	//std::cout << "dEv_dt1: " << std::endl;
	//dEv_dt1.PrintMatrix();

	//std::cout << "End of InternalPoints_2_2" << std::endl;
}

void BoundaryPoints_2_3(Matrix& W1, Matrix& W2, Matrix& W3,
	Vector& t1, Vector& t2, Vector& t3,
	Matrix& Xgr, Matrix& dXgr_dx, Matrix& dXgr_dy,
	Matrix& dEgr_dW1, Matrix& dEgr_dW2, Matrix& dEgr_dW3,
	Matrix& dEgr_dt1, Matrix& dEgr_dt2, Matrix& dEgr_dt3,
	double& Egr, size_t N, size_t M)
{
	//std::cout << "Start of BoundaryPoints_2_2" << std::endl;

	Matrix Agr(L, M);
	Matrix Bgr(L, M);
	Matrix ugr(1, M);
	Matrix Agrx(L, M);
	Matrix Agrx1(L, M);
	Matrix Bgrx(L, M);
	Matrix Bgrx1(L, M);
	Matrix ugrx(1, M);
	Matrix Agrxx1(L, M);
	Matrix Bgrxx(L, M);
	Matrix Bgrxx1(L, M);
	Matrix ugrxx(1, M);
	Matrix Agrx2(L, M);
	Matrix Bgrx2(L, M);
	Matrix Agrxx2(L, M);
	Matrix Bgrxx2(L, M);
	Matrix Agry(L, M);
	Matrix Agry1(L, M);
	Matrix Bgry(L, M);
	Matrix Bgry1(L, M);
	Matrix ugry(1, M);
	Matrix Agryy1(L, M);
	Matrix Bgryy(L, M);
	Matrix Bgryy1(L, M);
	Matrix ugryy(1, M);
	Matrix Agry2(L, M);
	Matrix Bgry2(L, M);
	Matrix Agryy2(L, M);
	Matrix Bgryy2(L, M);

	std::vector<double> tmp2(1, -1);
	Matrix dEgr_dugrxx(1, M);
	Matrix dEgr_dugryy(1, M);
	Matrix dEgr_dugrx(1, M); //        (25)
	Matrix dEgr_dugry(1, M); //        (25)
	Matrix dEgr_dugr(1, M); //         (25)

	Matrix dEgr_dBgr(L, M);
	Matrix dEgr_dBgrx(L, M);
	Matrix dEgr_dBgrxx(L, M);
	Matrix dEgr_dBgry(L, M);
	Matrix dEgr_dBgryy(L, M);

	Matrix dEgr_dAgr(L, M);
	Matrix dEgr_dAgrx(L, M);
	Matrix dEgr_dAgry(L, M);

	//2.3 Boundary points
	//2.3.1 
	Agr = W1 * Xgr + t1; // (8)
	//std::cout << "Agr: " << std::endl;
	//Agr.PrintMatrix();

	Bgr = W2 * Agr.sigmoidM() + t2; // (9)
	//std::cout << "Bgr: " << std::endl;
	//Bgr.PrintMatrix();

	ugr = W3 * Bgr.sigmoidM() + t3; // (10)
	//std::cout << "ugr: " << std::endl;
	//ugr.PrintMatrix();


	//2.3.2	
	Agrx = W1 * dXgr_dx; //           (11)
	//std::cout << "Agrx: " << std::endl;
	//Agrx.PrintMatrix();
	Agrx1 = Agr.sigmoidM_p1().Mul_by_El(Agrx);
	//std::cout << "Agrx1: " << std::endl;
	//Agrx1.PrintMatrix();
	Bgrx = W2 * Agrx1; //              (13)
	//std::cout << "Bgrx: " << std::endl;
	//Bgrx.PrintMatrix();
	Bgrx1 = Bgr.sigmoidM_p1().Mul_by_El(Bgrx);
	//std::cout << "Bgrx1: " << std::endl;
	//Bgrx1.PrintMatrix();
	ugrx = W3 * Bgrx1; //              (15)
	//std::cout << "ugrx: " << std::endl;
	//ugrx.PrintMatrix();
	Agrxx1 = Agr.sigmoidM_p2().Mul_by_El(Agrx).Mul_by_El(Agrx);
	//std::cout << "Agrxx1: " << std::endl;
	//Agrxx1.PrintMatrix();
	Bgrxx = W2 * Agrxx1; //            (18)
	//std::cout << "Bgrxx: " << std::endl;
	//Bgrxx.PrintMatrix();
	Bgrxx1 = Bgr.sigmoidM_p2().Mul_by_El(Bgrx).Mul_by_El(Bgrx) + Bgr.sigmoidM_p1().Mul_by_El(Bgrxx);
	//std::cout << "Bgrxx1: " << std::endl;
	//Bgrxx1.PrintMatrix();
	ugrxx = W3 * Bgrxx1; //            (20)
	//std::cout << "ugrxx: " << std::endl;
	//ugrxx.PrintMatrix();
	Agrx2 = Agr.sigmoidM_p2().Mul_by_El(Agrx); // (14)
	//std::cout << "Agrx2: " << std::endl;
	//Agrx2.PrintMatrix();
	Bgrx2 = Bgr.sigmoidM_p2().Mul_by_El(Bgrx);//  (16)
	//std::cout << "Bgrx2: " << std::endl;
	//Bgrx2.PrintMatrix();
	Agrxx2 = Agr.sigmoidM_p3().Mul_by_El(Agrx).Mul_by_El(Agrx); // (19)
	//std::cout << "Agrxx2: " << std::endl;
	//Agrxx2.PrintMatrix();
	Bgrxx2 = Bgr.sigmoidM_p3().Mul_by_El(Bgrx).Mul_by_El(Bgrx) + Bgr.sigmoidM_p2().Mul_by_El(Bgrxx); // (21)
	//std::cout << "Bgrxx2: " << std::endl;
	//Bgrxx2.PrintMatrix();

	//2.3.3
	Agry = W1 * dXgr_dy; //           (11)
	//std::cout << "Agry: " << std::endl;
	//Agry.PrintMatrix();
	Agry1 = Agr.sigmoidM_p1().Mul_by_El(Agry);
	//std::cout << "Agry1: " << std::endl;
	//Agry1.PrintMatrix();
	Bgry = W2 * Agry1; //              (13)
	//std::cout << "Bgry: " << std::endl;
	//Bgry.PrintMatrix();
	Bgry1 = Bgr.sigmoidM_p1().Mul_by_El(Bgry);
	//std::cout << "Bgry1: " << std::endl;
	//Bgry1.PrintMatrix();
	ugry = W3 * Bgry1; //              (15)
	//std::cout << "ugry: " << std::endl;
	//ugry.PrintMatrix();
	Agryy1 = Agr.sigmoidM_p2().Mul_by_El(Agry).Mul_by_El(Agry);
	//std::cout << "Agryy1: " << std::endl;
	//Agryy1.PrintMatrix();
	Bgryy = W2 * Agryy1; //            (18)
	//std::cout << "Bgryy: " << std::endl;
	//Bgryy.PrintMatrix();
	Bgryy1 = Bgr.sigmoidM_p2().Mul_by_El(Bgry).Mul_by_El(Bgry) + Bgr.sigmoidM_p1().Mul_by_El(Bgryy);
	//std::cout << "Bgryy1: " << std::endl;
	//Bgryy1.PrintMatrix();
	ugryy = W3 * Bgryy1; //            (20)
	//std::cout << "ugryy: " << std::endl;
	//ugryy.PrintMatrix();
	Agry2 = Agr.sigmoidM_p2().Mul_by_El(Agry); // (14)
	//std::cout << "Agry2: " << std::endl;
	//Agry2.PrintMatrix();
	Bgry2 = Bgr.sigmoidM_p2().Mul_by_El(Bgry);//  (16)
	//std::cout << "Bgry2: " << std::endl;
	//Bgry2.PrintMatrix();
	Agryy2 = Agr.sigmoidM_p3().Mul_by_El(Agry).Mul_by_El(Agry); // (19)
	//std::cout << "Agryy2: " << std::endl;
	//Agryy2.PrintMatrix();
	Bgryy2 = Bgr.sigmoidM_p3().Mul_by_El(Bgry).Mul_by_El(Bgry) + Bgr.sigmoidM_p2().Mul_by_El(Bgryy); // (21)
	//std::cout << "Bgryy2: " << std::endl;
	//Bgryy2.PrintMatrix();

	//2.3.3
	Egr = (ugr.Mul_by_El(ugr) + (Matrix({ Xgr.data[1]}).Mul_by_El(ugrx) + Matrix({ Xgr.data[0]}).Mul_by_El(ugry) * (-1)).Mul_by_El(Matrix({ Xgr.data[1] }).Mul_by_El(ugrx) + Matrix({ Xgr.data[0] }).Mul_by_El(ugry) * (-1))).Comprs_by_Col().RetData()[0][0] * (double(N) / double(M)) * 0.5; // (43)
	std::cout << "Egr = " << Egr << std::endl;

	//2.3.4
	dEgr_dugrx = Matrix({ Xgr.data[1] }).Mul_by_El(Matrix({ Xgr.data[1] }).Mul_by_El(ugrx) + Matrix({ Xgr.data[0] }).Mul_by_El(ugry) * (-1)) * (double(N) / double(M)); // (45)
	dEgr_dugry = Matrix({ Xgr.data[0] }).Mul_by_El(Matrix({ Xgr.data[1] }).Mul_by_El(ugrx) + Matrix({ Xgr.data[0] }).Mul_by_El(ugry) * (-1)) * (-1) * (double(N) / double(M)); // (46)	
	dEgr_dugr = ugr * (double(N) / double(M)); // (47)

	//2.3.5.6
	dEgr_dW3 = dEgr_dugr * Bgr.sigmoidM().T() + dEgr_dugrx * Bgrx1.T() + dEgr_dugry * Bgry1.T() + dEgr_dugrxx * Bgrxx1.T() + dEgr_dugryy * Bgryy1.T(); // (26)
	//std::cout << "dEgr_dW3: " << std::endl;
	//dEgr_dW3.PrintMatrix();

	dEgr_dt3 = dEgr_dugr.Comprs_by_Col(); // (27)
	//std::cout << "dEgr_dt3: " << std::endl;
	//dEgr_dt3.PrintMatrix();

	//2.3.5.7
	dEgr_dBgr = (W3.T() * dEgr_dugr).Mul_by_El(Bgr.sigmoidM_p1()) +
			(W3.T() * dEgr_dugrx).Mul_by_El(Bgrx2) +
			(W3.T() * dEgr_dugry).Mul_by_El(Bgry2) +
			(W3.T() * dEgr_dugrxx).Mul_by_El(Bgrxx2) +
			(W3.T() * dEgr_dugryy).Mul_by_El(Bgryy2); // (29) 	
	dEgr_dBgrx = (W3.T() * dEgr_dugrx).Mul_by_El(Bgr.sigmoidM_p1()) + (W3.T() * dEgr_dugrxx).Mul_by_El(Bgrx2) * 2.0; // (30)
	dEgr_dBgrxx = (W3.T() * dEgr_dugrxx).Mul_by_El(Bgr.sigmoidM_p1()); // (31)
	dEgr_dBgry = (W3.T() * dEgr_dugry).Mul_by_El(Bgr.sigmoidM_p1()) + (W3.T() * dEgr_dugryy).Mul_by_El(Bgry2) * 2.0; // (32)
	dEgr_dBgryy = (W3.T() * dEgr_dugryy).Mul_by_El(Bgr.sigmoidM_p1()); // (33)

	//2.3.5.8
	dEgr_dW2 = dEgr_dBgr * Agr.sigmoidM().T() +
			dEgr_dBgrx * Agrx1.T() +
			dEgr_dBgry * Agry1.T() +
			dEgr_dBgrxx * Agrxx1.T() +
			dEgr_dBgryy * Agryy1.T(); // (28)
	//std::cout << "dEgr_dW2: " << std::endl;
	//dEgr_dW2.PrintMatrix();

	dEgr_dt2 = dEgr_dBgr.Comprs_by_Col(); // (34)
	//std::cout << "dEgr_dt2: " << std::endl;
	//dEgr_dt2.PrintMatrix();

	//2.3.5.9
	dEgr_dAgr = (W2.T() * dEgr_dBgr).Mul_by_El(Agr.sigmoidM_p1()) +
			(W2.T() * dEgr_dBgrx).Mul_by_El(Agrx2) +
			(W2.T() * dEgr_dBgry).Mul_by_El(Agry2) +
			(W2.T() * dEgr_dBgrxx).Mul_by_El(Agrxx2) +
			(W2.T() * dEgr_dBgryy).Mul_by_El(Agryy2); // (35)
	dEgr_dAgrx = (W2.T() * dEgr_dBgrx).Mul_by_El(Agr.sigmoidM_p1()) + (W2.T() * dEgr_dBgrxx).Mul_by_El(Agrx2) * 2.0; // (36)
	dEgr_dAgry = (W2.T() * dEgr_dBgry).Mul_by_El(Agr.sigmoidM_p1()) + (W2.T() * dEgr_dBgryy).Mul_by_El(Agry2) * 2.0; // (37)	

	//2.3.5.10
	dEgr_dW1 = dEgr_dAgr * Xgr.T() + dEgr_dAgrx * dXgr_dx.T() + dEgr_dAgry * dXgr_dy.T(); // (38)
	//std::cout << "dEgr_dW1: " << std::endl;
	//dEgr_dW1.PrintMatrix();

	dEgr_dt1 = dEgr_dAgr.Comprs_by_Col(); // (39)
	//std::cout << "dEgr_dt1: " << std::endl;
	//dEgr_dt1.PrintMatrix();

	//std::cout << "End of BoundaryPoints_2_2" << std::endl;
}

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
					   Matrix& delta_W1_k, Matrix& delta_W2_k, Matrix& delta_W3_k,
					   Vector& delta_t1_k, Vector& delta_t2_k, Vector& delta_t3_k,
					   std::vector<double>& vec_error, double Ev, double Egr, size_t N, size_t M)
{

//std::cout << "Start of UpdateWeights_2_4" << std::endl;

//2.4 Updating the weights
//2.4.1 
	dE_dW1 = dEv_dW1 + dEgr_dW1;
	//std::cout << "dE_dW1: " << std::endl;
	//dE_dW1.PrintMatrix();

	dE_dt1 = dEv_dt1 + dEgr_dt1;
	//std::cout << "dE_dt1: " << std::endl;
	//dE_dt1.PrintMatrix();

	dE_dW2 = dEv_dW2 + dEgr_dW2;
	//std::cout << "dE_dW2: " << std::endl;
	//dE_dW2.PrintMatrix();

	dE_dt2 = dEv_dt2 + dEgr_dt2;
	//std::cout << "dE_dt2: " << std::endl;
	//dE_dt2.PrintMatrix();

	dE_dW3 = dEv_dW3 + dEgr_dW3;
	//std::cout << "dE_dW3: " << std::endl;
	//dE_dW3.PrintMatrix();

	dE_dt3 = dEv_dt3 + dEgr_dt3;
	//std::cout << "dE_dt3: " << std::endl;
	//dE_dt3.PrintMatrix();

	//2.4.2
	changedelta(dE_dW1, dE_dW1_k, delta_W1_k);
	//std::cout << "changedelta W1 completed" << std::endl;

	changedelta(dE_dW2, dE_dW2_k, delta_W2_k);
	//std::cout << "changedelta W2 completed" << std::endl;

	changedelta(dE_dW3, dE_dW3_k, delta_W3_k);
	//std::cout << "changedelta W3 completed" << std::endl;

	changedelta(dE_dt1, dE_dt1_k, delta_t1_k);
	//std::cout << "changedelta t1 completed" << std::endl;

	changedelta(dE_dt2, dE_dt2_k, delta_t2_k);
	//std::cout << "changedelta t2 completed" << std::endl;

	changedelta(dE_dt3, dE_dt3_k, delta_t3_k);
	//std::cout << "changedelta t3 completed" << std::endl;

	//2.4.3
	W1 = W1 + delta_W1_k.Mul_by_El(dE_dW1.sign(-1, 1, 0));
	//std::cout << "W1: " << std::endl;
	//W1.PrintMatrix();

	W2 = W2 + delta_W2_k.Mul_by_El(dE_dW2.sign(-1, 1, 0));
	//std::cout << "W2: " << std::endl;
	//W2.PrintMatrix();

	W3 = W3 + delta_W3_k.Mul_by_El(dE_dW3.sign(-1, 1, 0));
	//std::cout << "W3: " << std::endl;
	//W3.PrintMatrix();

	
	t1 = t1 + Matrix({ delta_t1_k.RetVector() }).T().Mul_by_El(dE_dt1.sign(-1, 1, 0)).toVector();
	//std::cout << "t1: " << std::endl;
	//t1.PrintVector();

	t2 = t2 + Matrix({ delta_t2_k.RetVector() }).T().Mul_by_El(dE_dt2.sign(-1, 1, 0)).toVector();
	//std::cout << "t2: " << std::endl;
	//t2.PrintVector();

	t3 = t3 + Matrix({ delta_t3_k.RetVector() }).T().Mul_by_El(dE_dt3.sign(-1, 1, 0)).toVector();
	//std::cout << "t3: " << std::endl;
	//t3.PrintVector();

	//2.4.4
	double error = sqrt((Ev + Egr) / (double(N) + double(M)));
	std::cout << "Error = " << error << std::endl;
	vec_error.push_back(error);

	//std::cout << "End of UpdateWeights_2_4" << std::endl;
}

void Verification(Matrix& W1, Matrix& W2, Matrix& W3,
				  Vector& t1, Vector& t2, Vector& t3,
				  std::vector<double>& vec_error)
{
	std::cout << "Start of Verification" << std::endl;

	//2.5
//2.5.1
	std::cout << "Vector of errors: " << std::endl;

	for (size_t k = 0; k < K; ++k)
		std::cout << vec_error[k] << ' ';
	std::cout << std::endl;

	//2.5.2
	std::vector<std::vector<std::pair<double, double>>> grid_test;
	for (double i = -1.0; i <= 1.0 + std::numeric_limits<double>::epsilon(); i += 0.02)
	{
		std::vector<std::pair<double, double>> vec_test;
		for (double j = -1.0; j <= 1.0 + std::numeric_limits<double>::epsilon(); j += 0.02) {
			vec_test.push_back(std::make_pair(i, j));
		}
		grid_test.push_back(vec_test);
	}

	// x^2 + y^2 <= 1
	size_t N_test = 0;

	for (size_t i = 0; i < grid_test.size(); ++i)
		for (size_t j = 0; j < grid_test[0].size(); ++j) {
			if (grid_test[i][j].first * grid_test[i][j].first + grid_test[i][j].second * grid_test[i][j].second > 1.0 + std::numeric_limits<double>::epsilon()) grid_test[i][j] = std::make_pair(777, 777);
			else ++N_test;
		}


	std::vector<std::vector<double>> vecx_test;
	std::vector<double> x_test;
	std::vector<double> y_test;

	for (size_t i = 0; i < grid_test.size(); ++i)
		for (size_t j = 0; j < grid_test[0].size(); ++j)
			if (grid_test[i][j].first != 777)
			{
				x_test.push_back(grid_test[i][j].first);
				y_test.push_back(grid_test[i][j].second);
			}
	vecx_test = { x_test, y_test };
	Matrix X_test(vecx_test);

	// dXv_test_dx and dXv_test_dy
	std::vector<double> zero_test(N_test, 0);
	std::vector<double> one_test(N_test, 1);

	Matrix dX_test_dx({ one_test, zero_test }); // (12)
	Matrix dX_test_dy({ zero_test, one_test }); // (17)

//2.5.3
	Matrix A_test(L, N_test);
	Matrix B_test(L, N_test);
	Matrix u_test(1, N_test);
	Matrix Ax_test(L, N_test);
	Matrix Ax1_test(L, N_test);
	Matrix Bx_test(L, N_test);
	Matrix Bx1_test(L, N_test);
	Matrix ux_test(1, N_test);
	Matrix Axx1_test(L, N_test);
	Matrix Bxx_test(L, N_test);
	Matrix Bxx1_test(L, N_test);
	Matrix uxx_test(1, N_test);


	Matrix Ay_test(L, N_test);
	Matrix Ay1_test(L, N_test);
	Matrix By_test(L, N_test);
	Matrix By1_test(L, N_test);
	Matrix uy_test(1, N_test);
	Matrix Ayy1_test(L, N_test);
	Matrix Byy_test(L, N_test);
	Matrix Byy1_test(L, N_test);
	Matrix uyy_test(1, N_test);

	std::vector<double> tmp1(1, -1);

	//2.2.1
	A_test = W1 * X_test + t1; //           (8)
	B_test = W2 * A_test.sigmoidM() + t2;//  (9)
	u_test = W3 * B_test.sigmoidM() + t3;//  (10)
	//std::cout << "u_test: " << std::endl;
	//u_test.PrintMatrix();

	//2.2.2
	Ax_test = W1 * dX_test_dx; //           (11)
	Ax1_test = A_test.sigmoidM_p1().Mul_by_El(Ax_test);
	Bx_test = W2 * Ax1_test; //              (13)
	Bx1_test = B_test.sigmoidM_p1().Mul_by_El(Bx_test);
	ux_test = W3 * Bx1_test; //              (15)
	Axx1_test = A_test.sigmoidM_p2().Mul_by_El(Ax_test).Mul_by_El(Ax_test);
	Bxx_test = W2 * Axx1_test; //            (18)
	Bxx1_test = B_test.sigmoidM_p2().Mul_by_El(Bx_test).Mul_by_El(Bx_test) + B_test.sigmoidM_p1().Mul_by_El(Bxx_test);
	uxx_test = W3 * Bxx1_test; //            (20)
	//2.2.3
	Ay_test = W1 * dX_test_dy; //           (11)
	Ay1_test = A_test.sigmoidM_p1().Mul_by_El(Ay_test);
	By_test = W2 * Ay1_test; //              (13)
	By1_test = B_test.sigmoidM_p1().Mul_by_El(By_test);
	uy_test = W3 * By1_test; //              (15)
	Ayy1_test = A_test.sigmoidM_p2().Mul_by_El(Ay_test).Mul_by_El(Ay_test);
	Byy_test = W2 * Ayy1_test; //            (18)
	Byy1_test = B_test.sigmoidM_p2().Mul_by_El(By_test).Mul_by_El(By_test) + B_test.sigmoidM_p1().Mul_by_El(Byy_test);
	uyy_test = W3 * Byy1_test; //            (20)	
//2.5.4
	double E_test = (uxx_test + uyy_test + tmp1).Mul_by_El(uxx_test + uyy_test + tmp1).Comprs_by_Col().RetData()[0][0] * 0.5; // (22)
	std::cout << "Test error = " << E_test << std::endl;
	std::cout << "Sqrt Test error = " << sqrt(E_test / N_test) << std::endl;
	
	//2.5.5
	std::vector<double> u_true_test;
	for (size_t i = 0; i < X_test.nC; ++i)
	{
		u_true_test.push_back((X_test.data[0][i] * X_test.data[0][i] + X_test.data[1][i] * X_test.data[1][i] - 1) / 4.0);
	}
	Matrix u_tr_test({ u_true_test });

	std::cout << "MSE = " << sqrt((u_test + u_tr_test * (-1)).Mul_by_El(u_test + u_tr_test * (-1)).Comprs_by_Col().RetData()[0][0]) / N_test << std::endl;
	std::cout << "Max Error = " << (u_test + u_tr_test * (-1)).MaxElementMatrix() << std::endl;
	std::cout << std::endl;
	std::cout << "End of Verification" << std::endl;
}
