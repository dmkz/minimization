/*
 * Авторы: Бураханова А., Золкин А.
 */

#include "TestingSuite.hpp"

int main (int argc, char* argv[])
{
	TestingSuite suite;
	if (argc > 1)
	{
		for (int i = 1; i < argc; i++) 
		{
			suite.RegisterCustomConfigurationFile(argv[i]);
		}
	} 
	else 
	{
		suite.RegisterDefaultConfigurationFile();
	}
	
	suite.RunAllTests();
}