/*
 * Авторы: Бураханова А., Золкин А., Арцишевский Л.
 */
 
#pragma once

#include "TMSNetTestElement.hpp"
#include <set>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip>

class UniquenessTest : public TMSNetTestElement
{
public:
	int Init(uint32_t dimension_, uint32_t point_num_, TMSNet* generator_ptr);
	int Run();

	uint32_t dimension;
	uint32_t point_num;
	uint32_t start;
	uint32_t finish;
};

int 
UniquenessTest::Init(uint32_t dimension_, uint32_t point_num_, TMSNet* generator_ptr) {
	generator = generator_ptr;
	dimension = dimension_;
	point_num = point_num_;
	start = 1;
	finish = dimension;
	
	return 0;
}


int
UniquenessTest::Run()
{
    std::cout << "-- Test Components Unique from dim = " << start << " to " << finish << std::endl;
    for (uint32_t dim_num = start; dim_num <= finish; dim_num++) {
		uint32_t cur_point_num = (point_num < 1ULL << dim_num) && point_num !=0 ? point_num : 1ULL << dim_num;
        std::vector<std::vector<long double>> table(dim_num, std::vector<long double>(cur_point_num));
		if(point_num == 0)
		{
			generator->Init(dim_num);
        } 
		else 
		{
			generator->Init(dim_num, point_num);
		}
		
        for (uint64_t i = 0; i < cur_point_num; i++) {
            auto r = generator->GeneratePoint();
            for (uint64_t j = 0; j < dim_num; ++j) {
                table[j][i] = r.coordinate[j];
            }
        }
        bool is_success = true;

        for (uint32_t i = 0; i < dim_num; ++i) {
            auto prev_size = table[i].size();
            std::sort(table[i].begin(), table[i].end());
            table[i].erase(std::unique(table[i].begin(), table[i].end()), table[i].end());
            auto curr_size = table[i].size();
            if (prev_size > curr_size) {
                is_success = false;
                break;
            }
        }

        std::cout << "\tdim = " << std::setw(2) << dim_num << (is_success ? ": ok, projections for each component are unique!" : ": bad") << std::endl;
    }
    return 0;
}