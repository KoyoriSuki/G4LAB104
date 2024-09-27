#ifndef VoltageFieldData_hh
#define VoltageFieldData_hh
#include <stdio.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <random>

std::string MakeKey(double x, double y, double z, double n)
{
	//
	//definition
	//
	std::ostringstream x_strstream, y_strstream, z_strstream;
	std::string x_str, y_str, z_str, key;

	//
	//prevent -0.0 (e.g. -0.04->0.0 instead of -0.0)
	//
	if (fabs(x) < 0.5 * pow(10, -n))x = 0.0;
	if (fabs(y) < 0.5 * pow(10, -n))y = 0.0;
	if (fabs(z) < 0.5 * pow(10, -n))z = 0.0;

	//
	//convert double to string
	//
	x_strstream << std::fixed << std::setprecision(n) << x;
	y_strstream << std::fixed << std::setprecision(n) << y;
	z_strstream << std::fixed << std::setprecision(n) << z;
	key = x_strstream.str() + " " + y_strstream.str() + " " + z_strstream.str();
	return key;
}

double UniformRand(double min, double max) {
	static std::default_random_engine generator;
	static int seedflag = 1;
	if (seedflag) {
		generator.seed(time(NULL));
		seedflag = 0;
	}
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}

double floor_n(double x, double n)
{
	return floor(x * pow(10, n)) / pow(10, n);
}

double ceil_n(double x, double n)
{
	return (floor(x * pow(10, n)) + 1) / pow(10, n);
}

double FieldValueInterpolation(double x, double y, double z, double n, std::unordered_map<std::string, double> &VoltageField)
{
	//
	//weight
	//
	double u = 1, v = 0, w = 0;
	u = (x - floor_n(x, n)) / pow(10, -n);
	v = (y - floor_n(y, n)) / pow(10, -n);
	w = (z - floor_n(z, n)) / pow(10, -n);

	//
	//average by weight
	//
	double value = 0;
	value =
		(1 - u) * (1 - v) * (1 - w) * VoltageField.at(MakeKey(floor_n(x, n), floor_n(y, n), floor_n(z, n), n)) +
		(u) * (1 - v) * (1 - w) * VoltageField.at(MakeKey(ceil_n(x, n), floor_n(y, n), floor_n(z, n), n)) +
		(u) * (v) * (1 - w) * VoltageField.at(MakeKey(ceil_n(x, n), ceil_n(y, n), floor_n(z, n), n)) +
		(u) * (v) * (w)*VoltageField.at(MakeKey(ceil_n(x, n), ceil_n(y, n), ceil_n(z, n), n)) +
		(1 - u) * (v) * (1 - w) * VoltageField.at(MakeKey(floor_n(x, n), ceil_n(y, n), floor_n(z, n), n)) +
		(1 - u) * (v) * (w)*VoltageField.at(MakeKey(floor_n(x, n), ceil_n(y, n), ceil_n(z, n), n)) +
		(1 - u) * (1 - v) * (w)*VoltageField.at(MakeKey(floor_n(x, n), floor_n(y, n), ceil_n(z, n), n)) +
		(u) * (1 - v) * (w)*VoltageField.at(MakeKey(ceil_n(x, n), floor_n(y, n), ceil_n(z, n), n));
	return value;
}

int NearAnode(double x, double y, double z)
{
	if (x > 3.8 && x < 6.2 && y>3.8 && y < 5 && z>3.8 && z < 6.2)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

double GetFieldValue(double x, double y, double z, std::unordered_map<std::string, double>& VoltageField)
{
	int n = 2;
	if (NearAnode(x, y, z))
	{
		n = 2;
	}
	else
	{
		n = 1;
	}
	return FieldValueInterpolation(x, y, z, n, VoltageField);
}

//int VoltageFieldData()
//{
//	double n = 2;
//	std::unordered_map<std::string, double> VoltageField;
//	
//	//
//	//construct a hash table to store the voltage field data
//	//
//	FILE* data = fopen("D:/aNinmu/Aresearch/share/root/CZTVoltageField/input_file/finetest.fld", "r+");
//	double x = 0, y = 0, z = 0, v = 0;
//	int count = 0;
//	while (!feof(data))
//	{
//		fscanf(data, "%lf %lf %lf %lf", &x, &y, &z, &v);
//		x *= 1000;//(unit) m to mm
//		x += 5;
//		y *= 1000;
//		y += 2.5;
//		z *= 1000;
//		z += 5;
//		VoltageField[MakeKey(x, y, z, n)] = v;
//		count++;
//		if (count % 100000 == 0)std::cout << count << std::endl;
//	}
//	fclose(data);
//	n = 1;
//	data = fopen("D:/aNinmu/Aresearch/share/root/CZTVoltageField/input_file/444test.fld", "r+");
//	while (!feof(data))
//	{
//		fscanf(data, "%lf %lf %lf %lf", &x, &y, &z, &v);
//		x *= 1000;//(unit) m to mm
//		x += 5;
//		y *= 1000;
//		y += 2.5;
//		z *= 1000;
//		z += 5;
//		VoltageField[MakeKey(x, y, z, n)] = v;
//	}
//	fclose(data);
//
//	//
//	//hash table test
//	//
//	for (int i = 0; i < 10000; i++)
//	{
//		int n = 2;
//		double x_test = UniformRand(0, 10), y_test = UniformRand(0, 5), z_test = UniformRand(0, 10);
//		if (NearAnode(x_test, y_test, z_test))
//		{
//			n = 2;
//		}
//		else
//		{
//			n = 1;
//		}
//		std::cout << MakeKey(x_test, y_test, z_test, n) << std::endl;
//		//cout << x_test << " " << y_test << " " << z_test << " " << VoltageField.at(MakeKey(x_test, y_test, z_test, n)) << " " << VoltageField.count(MakeKey(x_test, y_test, z_test, n)) << endl;
//		std::cout << GetFieldValue(x_test, y_test, z_test, VoltageField) << std::endl;
//	}
//	return 0;
//}
#endif