/*
 * Авторы: Бураханова А., Золкин А., Арцишевский Л.
 */
 
#pragma once 

#include <cmath>
#include <fstream>
#include <set>
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
	std::cout << "\n" << test_name << " start";
    std::set<Real> s;
    for (uint32_t i = 0; i < dimension; i++)
    {
        for (uint32_t j = 0; j < point_num; j++)
        {
			auto point = generator->GeneratePoint();
            s.insert(point.coordinate[i]);
        }
        if (s.size() < point_num)
        {
			std::cout << "\nGenerator does not pass the test at dimension " << i << "!"; 
			std::cout << "\n" << test_name << " end";
            return 1;
        }
        s.clear();
		generator->Reset();
    }
	
	std::cout << "\nGenerator passed the test " << test_name << "!"; 
	std::cout << "\n" << test_name << " end";
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
	std::cout << "\n" << test_name << " start";
	if (function_key == "subcube")
	{
		int res = SubcubeTest();
		std::cout << "\n" << test_name << " end";
		return res;
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
		Real result = 0;
		for (uint32_t i = 0; i < point_num; i++)
		{
			auto point = generator->GeneratePoint();
			result += SubcubeTestFunction(point, d);
		}
		result /= point_num;
		result = fabs(result - SubcubeTestAnalyticValue(d));
		if (write_output)
		{
			output_stream << d << " " << result << "\n";
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

class ProjectionTest : public TMSNetTestElement
{
public:
	int RunTest();
	
	uint32_t dimension;
	uint32_t x;
	uint32_t y;
	uint32_t point_num;
};

int
ProjectionTest::RunTest()
{
	std::cout << "\n" << test_name << " start";
	std::ofstream output_stream;
	output_stream.open(test_name);
	for (uint32_t i = 0; i < point_num; i++)
	{
		auto point = generator->GeneratePoint();
		output_stream << point.coordinate[x-1] << " " << point.coordinate[y-1] << "\n";
	}
	output_stream.close();
	std::cout << "\n" << test_name << " end";

	return 0;
}

class OrthogonalityTest : public TMSNetTestElement
{
public:
	int RunTest();
	
	uint32_t dimension;
	uint32_t b;
	uint32_t point_num;
};

int
OrthogonalityTest::RunTest()
{
	std::cout << "\n" << test_name << " start";
	uint32_t r = 2u;
	uint32_t b_pow_r = std::pow(b, r);
	uint32_t lambda = point_num / b_pow_r;
	

	for(uint32_t j1 = 0; j1 < dimension; j1++)
	{
		for(uint32_t j2 = j1 + 1; j2 < dimension; j2++)
		{
			std::vector<uint32_t> appearance(b_pow_r, 0u);
			for(uint32_t i = 0; i < point_num; i++)
			{
				auto point = generator->GeneratePoint();
				uint32_t first_coord = 0;
				uint32_t second_coord = 0;
				for (uint32_t j = 0; j < b; j++)
				{
					if (point.coordinate[j1] >= (double)(j)/b && point.coordinate[j1] < (double)((j+1))/b)
					{
						first_coord = j;
					}
				}
				for (uint32_t j = 0; j < b; j++)
				{
					if (point.coordinate[j2] >= (double)(j)/b && point.coordinate[j2] < (double)((j+1))/b)
					{
						second_coord = j;
					}
				}
				appearance[first_coord + second_coord * b] += 1;
			}
			
			for(uint32_t j = 0; j < appearance.size(); j++)
			{
				if (appearance[j] != lambda)
				{
					std::cerr << "\nOrthogonality test(" << test_name << "):Wrong number of pair " << j1 <<", "<< j2 << "(axis pair) appearance: pair(" << j % b << ", " << j / b << ") appears " << appearance[j] << " times!!!";					
					std::cout << "\n" << test_name << " end";
					return 1;
				}
			}
			generator->Reset();
		}
	}
	
	std::cout << "\nGenerator passed the test " << test_name << "!"; 
	std::cout << "\n" << test_name << " end";
	
	return 0;
}
