#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "sobolseqgenerator.h"

int SobolSeqGenerator::Init()
{
    return Init(1, 1, "sobol_net_default.txt");
}

int SobolSeqGenerator::Init(uint32_t _N, uint32_t _D, std::string dir_file)
{
    N = _N;
    D = _D;
    current_point_number = -1;
    L = (uint32_t)std::ceil(std::log(Real(N))/std::log(2.0L));

    std::ifstream infile(dir_file, std::ios::in);
    if (!infile) {
		std::cout << "Не найден файл содержащий направляющие числа!\n";
		return -1;
    }

    infile.close();

    params_filename = dir_file;
    return 0;
}

PointReal SobolSeqGenerator::GeneratePoint()
{
    if (N == 0 || D == 0)
    {
        std::cout << "Генератор сеток Соболева некорректно инициализирован! N = " << N << ", D = " << D;
        return PointReal();
    }
    
	std::ifstream infile(params_filename, std::ios::in);
    if (!infile) {
		std::cout << "Не найден файл содержащий направляющие числа!\n";
		return -1;
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');
	
    if (current_point_number == N - 1)
    {
        std::cout << "Генератор уже сгенерировал все N точек!";
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
    for (uint32_t j = 1; j < D; j++) {

        // Чтение параметров из файла
        uint32_t d, s;
        uint32_t a;
        infile >> d >> s >> a;
        std::vector<uint32_t> m(s+1);
        for (uint32_t i = 1; i <= s; i++)
        {
            infile >> m[i];
        }

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

    infile.close();
	
    last_generated_point = result_point;

    PointReal final_result = PointReal(D);
    for (uint32_t j = 0; j < D; j++)
    {
        final_result.coordinate[j] = (Real) result_point.coordinate[j] / std::pow(Real(2.0), 32);
    }

    return final_result;
}
