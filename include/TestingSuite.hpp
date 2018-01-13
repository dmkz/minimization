#pragma once

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include "tinyxml2.hpp"
#include "TMSNetTestElement.hpp"
#include "sobolseqgenerator.hpp"

#define DEFAULT_CONFIGURATION_FILENAME "default.xml"
#define TESTGROUP_TAGNAME "testgroup"

class TestingSuite
{
public:
	int RegisterDefaultConfigurationFile();
	int RegisterCustomConfigurationFile(std::string config_filename);
	int RunAllTests();
	
	int WriteTestGroupConfiguration(tinyxml2::XMLElement* test_group);
	
	~TestingSuite();
private:
	std::vector<TMSNetTestElement* > tests;
};

TestingSuite::~TestingSuite()
{
	for (auto& test: tests)
	{
		delete test;
	}
}

int 
TestingSuite::RegisterDefaultConfigurationFile()
{
	return RegisterCustomConfigurationFile(std::string(DEFAULT_CONFIGURATION_FILENAME));
}

int
TestingSuite::RegisterCustomConfigurationFile(std::string config_filename)
{
	tinyxml2::XMLDocument doc;
	auto error_code = doc.LoadFile(config_filename.c_str());
    if (error_code != tinyxml2::XMLError::XML_SUCCESS)
	{
		std::cerr << "\nError loading \"" + config_filename + "\" configuration file: " << error_code  << " (tinyxml2::XMLError)";
		return error_code;
	}
	
	tinyxml2::XMLNode* config_root = doc.FirstChild();
	if (config_root == nullptr) return tinyxml2::XMLError::XML_ERROR_FILE_READ_ERROR;
	
	tinyxml2::XMLElement* test_group = config_root->FirstChildElement(TESTGROUP_TAGNAME);
	while (test_group != nullptr)
	{
		WriteTestGroupConfiguration(test_group);
		test_group = test_group->NextSiblingElement(TESTGROUP_TAGNAME);
	}
	
	return 0;
}

int
TestingSuite::WriteTestGroupConfiguration(tinyxml2::XMLElement* test_group)
{
	std::string new_testgroup_name = test_group->Attribute("name");
	if (new_testgroup_name.empty()) 
	{
		new_testgroup_name = "testnoname";
		std::cerr << "\nTest doesn't have a name. Assigning name \"" + new_testgroup_name + "\"";
	}
	
	for(auto tms_net_generator = test_group->FirstChildElement("tmsnetgenerator"); tms_net_generator != nullptr; tms_net_generator = tms_net_generator->NextSiblingElement("tmsnetgenerator"))
	{
		std::string new_generator_name = tms_net_generator->Attribute("name");
		if (new_generator_name.empty()) 
		{
			std::cerr << "\nGenerator in test group \"" + new_testgroup_name + "\" doesn't have a name. Skipping";
			continue;
		}
		
		for(auto uniqueness_test = test_group->FirstChildElement("uniqueness"); uniqueness_test != nullptr; uniqueness_test = uniqueness_test->NextSiblingElement("uniqueness"))
		{
			UniquenessTest* new_test = new UniquenessTest();
			uniqueness_test->QueryBoolAttribute("writeoutput", &new_test->write_output);
			uniqueness_test->QueryUnsignedAttribute("dimension", &new_test->dimension);
			uniqueness_test->QueryUnsignedAttribute("pointnum", &new_test->point_num);
			if (new_generator_name == "joe-kuo")
			{
				std::string filename = tms_net_generator->Attribute("filename");
				std::string new_test_name = new_testgroup_name + "_" + 
					filename + 
					"_uniqueness_" + 
					uniqueness_test->Attribute("dimension") + "_" + 
					uniqueness_test->Attribute("pointnum");
				new_test->test_name = new_test_name;
				new_test->generator = new SobolSeqGenerator();
				new_test->generator->Init(new_test->dimension, new_test->point_num, filename);
				TMSNetTestElement* test_ptr = new_test; 
				tests.push_back(test_ptr);
			}
			else
			{
				std::cerr << "\nNo new tests for this generator have been implemented yet or wrong generator name!(" << new_generator_name << ")";
				delete new_test;
				break;
			}
		}
		
		for(auto integration_test = test_group->FirstChildElement("integration"); integration_test != nullptr; integration_test = integration_test->NextSiblingElement("integration"))
		{
			for(auto integration_params = integration_test->FirstChildElement("parameters"); integration_params != nullptr; integration_params = integration_params->NextSiblingElement("parameters"))
			{
				IntegrationTest* new_test = new IntegrationTest();
				integration_test->QueryBoolAttribute("writeoutput", &new_test->write_output);
				integration_test->QueryUnsignedAttribute("dimension", &new_test->dimension);
				integration_test->QueryUnsignedAttribute("pointnum", &new_test->point_num);
				new_test->function_key = integration_params->Attribute("function_key");
				std::cout << "\nRegistration of integration test: START";
				if (new_generator_name == "joe-kuo")
				{
					std::string filename = tms_net_generator->Attribute("filename");
					
					std::cout << "\nRegistration of integration test: STEP";
					std::string new_test_name = new_testgroup_name + "_" + 
						filename + "_" +
						new_test->function_key + "_" +
						"integration" + "_" + 
						integration_test->Attribute("dimension") + "_" + 
						integration_test->Attribute("pointnum");
						
					std::cout << "\nRegistration of integration test: STEP";
					new_test->test_name = new_test_name;
					
					std::cout << "\nRegistration of integration test: STEP";
					new_test->generator = new SobolSeqGenerator();
					
					std::cout << "\nRegistration of integration test: STEP";
					new_test->generator->Init(new_test->dimension, new_test->point_num, filename);
					
					std::cout << "\nRegistration of integration test: STEP";
					TMSNetTestElement* test_ptr = new_test; 
					
					std::cout << "\nRegistration of integration test: STEP";
					tests.push_back(test_ptr);
					
					std::cout << "\nRegistration of integration test: END";
				}
				else
				{
					std::cerr << "\nNo new tests for this generator have been implemented yet or wrong generator name!(" << new_generator_name << ")";
					delete new_test;
					break;
				}
			}
		}
	}
	
	return 0;
}

int
TestingSuite::RunAllTests()
{
	for (auto& test: tests) 
	{
		test->RunTest();
	}
	
	return 0;
}
