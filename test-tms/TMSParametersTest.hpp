/*
 * Авторы: Бураханова А., Золкин А.
 */
 
#pragma once

#include "TMSNetTestElement.hpp"
#include <set>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <numeric>

class TMSParametersTest : public TMSNetTestElement
{
public:
	int Init(uint32_t t_, uint32_t m_, uint32_t s_, uint32_t b_, TMSNet* generator_ptr);
	int Run();
	uint32_t CountNumberOfPointsWithinInterval(std::vector<uint32_t> elem_interval_sides_power,
		std::vector<uint32_t> elem_interval_sides_offset,	uint32_t elem_interval_vol_power);

	uint32_t t;
	uint32_t m;
	uint32_t s;
	uint32_t b;
};

int 
TMSParametersTest::Init(uint32_t t_, uint32_t m_, uint32_t s_, uint32_t b_, TMSNet* generator_ptr) {
	t = t_;
	m = m_;
	s = s_;
	b = b_;
	generator = generator_ptr;
	
	uint32_t point_num = 1ULL << m;
	generator->Init(s, point_num);
	
	return 0;
}


int
TMSParametersTest::Run()
{
	uint32_t point_num = 1ULL << m;
	uint32_t elem_interval_vol_power = m - t;
	uint32_t target_point_num = std::pow(b, t);
	uint32_t number_of_intervals = 0;
	
	std::vector<uint32_t> elem_interval_sides_power(s, 0);
	
	for(uint32_t i = 0; i < std::pow(elem_interval_vol_power + 1, s); i++)
	{
		if(elem_interval_sides_power[s-1] == elem_interval_vol_power && elem_interval_vol_power != 0)
		{
			break;
		}

		for(uint32_t j = 0; j < s; j++)
		{
			if(elem_interval_sides_power[j] == elem_interval_vol_power)
			{
				elem_interval_sides_power[j] = 0;
			} else {
				elem_interval_sides_power[j] += 1;
				break;
			}
		}		

		if(std::accumulate(elem_interval_sides_power.begin(), elem_interval_sides_power.end(), 0ULL) == elem_interval_vol_power)
		{
			std::vector<uint32_t> elem_interval_sides_offset(s, 0);
			while(elem_interval_sides_offset[s-1] != std::pow(b, elem_interval_sides_power[s-1]))
			{
				
				number_of_intervals++;
				
				uint32_t elem_interval_point_num = CountNumberOfPointsWithinInterval(elem_interval_sides_power, elem_interval_sides_offset, point_num);
				if(elem_interval_point_num != target_point_num)
				{
					/*std::cout << "\nFail: number of points in elementary interval " << elem_interval_point_num;
					std::cout << "\nElementary interval: ";
					for(auto &side_power : elem_interval_sides_power)
					{
						std::cout << side_power << " ";
					}				
					std::cout << "\nThe interval offest: ";
					for(auto &side_offset : elem_interval_sides_offset)
					{
						std::cout << side_offset << " ";
					}
					*/
					return 1;
				}
				
				for(uint32_t j = 0; j < s; j++)
				{
					if(elem_interval_sides_offset[j] == std::pow(b, elem_interval_sides_power[j])-1)
					{
						if(j == s - 1)
						{
							elem_interval_sides_offset[j] += 1;
							break;
						} else {
							elem_interval_sides_offset[j] = 0;
						}
					} else {
						elem_interval_sides_offset[j] += 1;
						break;
					}
				}
			}
		}
	}
	
	//std::cout << "\nTest is finished: number of intervals = " << number_of_intervals;
    return 0;
}

uint32_t
TMSParametersTest::CountNumberOfPointsWithinInterval(std::vector<uint32_t> elem_interval_sides_power, 
	std::vector<uint32_t> elem_interval_sides_offset, uint32_t point_num)
{
	uint32_t res = 0;
	generator->Reset();
	
	for(uint32_t i = 0; i < point_num; i++)
	{
		auto point = generator->GeneratePoint().coordinate;
		bool is_inside_interval = true;
		
		for(uint32_t j = 0; j < point.size(); j++)
		{
			if( (elem_interval_sides_offset[j]  )/std::pow(b, elem_interval_sides_power[j]) <= point[j] 
			 && (elem_interval_sides_offset[j]+1)/std::pow(b, elem_interval_sides_power[j]) >  point[j])
			{
				continue;
			} else {
				is_inside_interval = false;
				break;
			}
		}
		
		if(is_inside_interval)
		{
			res++;
		}
	}
	
	return res;
}