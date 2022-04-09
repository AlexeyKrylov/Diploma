#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "mkl.h"

const double M_PI = 3.14159265358979323846;  /* pi */
const double eta_plus = 1.2;
const double eta_minus = 0.5;
const unsigned int K = 1000;
const unsigned int L = 128;
const double lambda = 0.6;
const double ro_1 = 1.0;
const double mu = -1.0;
const double p1_const = 0.0;
const double p2_const = 0.0;
const double a_const = 0.0;
const double b_const = 0.0;
const double c_const = 1.0;
const double lr_scheduler = 1000;
const size_t checkpoint = 950;
const bool debug = false;
const double H = 24;