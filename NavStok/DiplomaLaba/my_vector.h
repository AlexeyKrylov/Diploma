#pragma once
#include "libs_consts.h"

//2.0.1-3.1
class Vector {
private:
	std::vector<double> data;

	size_t dim;

public:
	Vector(size_t d);

	Vector(std::vector<double> d);

	std::vector<double> RetVector() const;

	void PrintVector() const;

	void RandomInit();

	void InitWithValue(double value);

	Vector operator+(Vector second) const;

	size_t RetDim() const;

	void from_file_txt(std::string filename);

	void to_file_txt(std::string filename) const;

	~Vector();
};