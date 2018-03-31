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
