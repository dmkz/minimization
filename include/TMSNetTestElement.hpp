#pragma once 

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

	uint32_t dimension;
	uint32_t point_num;
	std::string function_key;
	std::map<std::string, double> function_params;
};

int
IntegrationTest::RunTest()
{
	// TODO: Test
	std::cout << "\n" << test_name << " IntegrationTest";
	return 0;
}