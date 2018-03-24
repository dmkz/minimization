/*
 * Авторы: Бураханова А., Золкин А.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "sobolseqgenerator.hpp"

int 
SobolSeqGenerator::Init()
{
    return Init(1, 1, "sobol_Net_default.txt");
}

int 
SobolSeqGenerator::Init(uint32_t N_, uint32_t D_, std::string dir_file)
{
    N = N_;
    D = D_;
    current_point_number = -1;

    L = (uint32_t)std::ceil(std::log(Real(N))/std::log(2.0L));

    std::ifstream infile(dir_file, std::ios::in);
    if (!infile) {
		std::cout << "\nDirection numbers file not found!" << "(" << dir_file << ")";
		return -1;
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');

    dir_num_params.push_back(DirectionalNumbersParams());

    for (unsigned j = 1; j < D; j++)
    {
        dir_num_params.push_back(DirectionalNumbersParams());        
        // Чтение параметров из файла
        unsigned D_, s_;
        unsigned a_;
        infile >> D_ >> s_ >> a_;
        dir_num_params[j].d = D_;
        dir_num_params[j].s = s_;
        dir_num_params[j].a = a_;
        dir_num_params[j].m.push_back(1);
        for (unsigned i = 1; i <= s_; i++)
        {
            unsigned m_;
            infile >> m_;
            dir_num_params[j].m.push_back(m_);
        }
    }

    infile.close();

    params_filename = dir_file;
    return 0;
}

PointReal 
SobolSeqGenerator::GeneratePoint()
{
    if (N == 0 || D == 0)
    {
        std::cout << "Sobol sequence generator isnt initialized correctly! N = " << N << ", D = " << D;
        return PointReal();
    }

    if(current_point_number == N - 1)

    {
        std::cout << "Generator has already generated N points!";
        return PointReal();
    }

    current_point_number += 1;

    if(current_point_number == 0)
    {
        last_generated_point = PointUnsigned(D, std::vector<uint32_t>(D, 0));
        return PointReal(D, std::vector<Real>(D, 0));
    }
	
    PointUnsigned result_point = PointUnsigned(D);

    C = 1;

    uint32_t value = current_point_number - 1;
    while (value & 1) {

        value >>= 1;
        C++;
    }
	
    // Вычислить направляюoее число V, умноженное на pow(2,32)
    uint32_t V_first = 1 << (32u-C);

    // Вычислить первую координату, умноженную на pow(2,32)
    result_point.coordinate[0] = last_generated_point.coordinate[0] ^ V_first;
	
    // ----- Вычислить остальные координаты -----
    for (unsigned j = 1; j < D; j++)
    {

        // Чтение параметров из файла
        unsigned s = dir_num_params[j].s;
        unsigned a = dir_num_params[j].a;
        auto m = dir_num_params[j].m;

        // Вычислить направляющие числа V, умноженные на pow(2,32)
        std::vector<uint32_t> V(L + 1);
        if (L <= s)
        {
            for (uint32_t i = 1; i <= L; i++)
            {
                V[i] = m[i] << (32u - i);
            }
        }
        else
        {
            for (uint32_t i = 1; i <= s; i++)
            {
                V[i] = m[i] << (32u - i);
            }

            for (uint32_t i = s + 1; i <= L; i++)
            {
                V[i] = V[i-s] ^ (V[i-s] >> s);
  	            for (uint32_t k = 1; k <= s-1; k++)
                {
                    V[i] ^= (((a >> (s - 1 - k)) & 1u) * V[i-k]);
                }
            }
        }

        // Вычислить j-ю координату, умноженную на pow(2,32)
        result_point.coordinate[j] = last_generated_point.coordinate[j] ^ V[C];
    }

    last_generated_point = result_point;

    PointReal final_result = PointReal(D);
    for (uint32_t j = 0; j < D; j++)
    {
        final_result.coordinate[j] = (Real) result_point.coordinate[j] / std::pow(Real(2.0), 32);
    }

    return final_result;
}

int
SobolSeqGenerator::Reset()
{
	current_point_number = -1;
	return 0;
}
