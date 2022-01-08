#include "my_vector.h"	
#include "libs_consts.h"


Vector::Vector(size_t d) : dim(d)
{
	for (size_t i = 0; i < dim; ++i)
		data.push_back(0);
}

Vector::Vector(std::vector<double> d) : data(d)
{
	dim = d.size();
}

std::vector<double> Vector::RetVector() const
{
	return data;
}

void Vector::PrintVector() const
{
	std::cout << "| ";
	for (size_t i = 0; i < dim; ++i)
		std::cout << data[i] << ", ";
	std::cout << "|" << std::endl;
}

void Vector::RandomInit()
{
	for (size_t i = 0; i < dim; ++i)
		data[i] = (rand() % 200) / 100.0 - 1.0;
}

void Vector::InitWithValue(double value)
{
	for (size_t i = 0; i < dim; ++i)
		data[i] = value;
}

Vector Vector::operator+(Vector second) const
{
	if (dim != second.dim)
	{
		std::cout << "Error vector+Vector" << std::endl;
		exit(1);
	}
	Vector res(dim);
	for (size_t i = 0; i < dim; ++i)
		res.data[i] = data[i] + second.data[i];
	return res;
}

size_t Vector::RetDim() const { return dim; }

void Vector::from_file_txt(std::string filename)
{
	std::ifstream fin(filename);

	if (fin.is_open())
	{
		for (size_t i = 0; i < dim; ++i)
			fin >> data[i];
	}
	fin.close();
}

void Vector::to_file_txt(std::string filename) const
{
	std::ofstream fout(filename);

	if (fout.is_open())
	{
		for (size_t i = 0; i < dim; ++i)
			fout << data[i] << ' ';
		fout << std::endl;
	}
	fout.close();
}

Vector::~Vector() {}