/*
 * Авторы: Бураханова А., Золкин А., Арцишевский Л.
 */

#include <cmath>
#include <fstream>
#include <set>
#include "TMSNetTestElement.hpp"

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