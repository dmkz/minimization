#pragma once 

#include <cmath>
#include <fstream>
#include "tmsnet.hpp"

class TMSNetTestElement
{
public:
	virtual int RunTest() = 0;
	~TMSNetTestElement();
	
	
	std::string test_name;
	TMSNet* generator;
	bool write_output;
};

TMSNetTestElement::~TMSNetTestElement()
{
	delete generator;
}

class UniquenessTest : public TMSNetTestElement
{
public:
	int RunTest();
	
	uint32_t dimension;
	uint32_t point_num;
};

int
UniquenessTest::RunTest()
{
	// TODO
	std::cout << "\n" << test_name << " UniquenessTest";
	auto p = generator->GeneratePoint();
	p = generator->GeneratePoint();
	std::cout << p.coordinate[0] << p.coordinate[p.coordinate.size()-1];
	return 0;
}


class IntegrationTest : public TMSNetTestElement
{
public:
	int RunTest();

	int SubcubeTest();
	double SubcubeTestFunction(PointReal point, uint32_t d);
	double SubcubeTestAnalyticValue(uint32_t d);
	
	uint32_t dimension;
	uint32_t point_num;
	std::string function_key;
};

int
IntegrationTest::RunTest()
{	
	if (function_key == "subcube")
	{
		return SubcubeTest();
	}
	std::cout << "\nNo new tests for this function key have been implemented yet or wrong function key!(" << function_key << ")";
	return -1;
}

int
IntegrationTest::SubcubeTest()
{
	std::ofstream output_stream;
	output_stream.open(test_name);
	for (uint32_t d = 1; d <= dimension; d++)
	{
		double result = 0;
		for (uint32_t i = 0; i < point_num; i++)
		{
			auto point = generator->GeneratePoint();
			result += SubcubeTestFunction(point, d);
		}
		result /= point_num;
		std::cout << "\n" << result << " " << point_num << " " << SubcubeTestAnalyticValue(d);
		result = fabs(result - SubcubeTestAnalyticValue(d));
		std::cout << "\nD: " << d << ", delta = " << result;
		if (write_output)
		{
			output_stream << d << result;
		}
		generator->Reset();
	}
	output_stream.close();
	
	return 0;
}

double
IntegrationTest::SubcubeTestFunction(PointReal point, uint32_t d)
{
	for (uint32_t i = 0; i < d; i++)
	{
		auto coord = point.coordinate[i];
		if (coord > 0.5) return 0.0;
	}
	
	return 1.0;
}

double
IntegrationTest::SubcubeTestAnalyticValue(uint32_t d)
{
	return pow(0.5, d);
}