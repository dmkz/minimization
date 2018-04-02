/*
 * Авторы: Бураханова А., Золкин А., Арцишевский Л.
 */
 
#pragma once 

#include "tmsnet.hpp"
#include <cmath>
#include <ctime>
#include <fstream>
#include <set>

class TMSNetTestElement
{
public:
	~TMSNetTestElement();
	// Прогонка теста
	virtual int Run() = 0;
	
	// Получение даты\времени в формате YYYY-MM-DD.HH:mm:ss
	const std::string CurrentDateTime();
	
	
	std::string test_name;
	TMSNet* generator;
	bool write_output;
};

TMSNetTestElement::~TMSNetTestElement()
{}

const std::string 
TMSNetTestElement::CurrentDateTime() {
    time_t     now = std::time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
    
    return buf;
}
