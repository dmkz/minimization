#pragma once

#include <vector>
#include <string>
#include "tmsnet.hpp"

typedef Point<uint32_t> PointUnsigned;

struct DirectionalNumbersParams
{
        unsigned d, s;
        unsigned a;
        std::vector<unsigned> m;
};

class SobolSeqGenerator : public TMSNet
{
public:
    int Init();
    int Init(uint32_t N_, uint32_t D_, std::string dir_file);
	int Reset();
	
    PointReal GeneratePoint();

private:
	// Количество точек, которые нужно сгенерировать 
    uint32_t N;
	// Размерность сетки
    uint32_t D;
	
    // Максимальное число нужных битов
    uint32_t L;
	// ???
    uint32_t C;

    // Параметры направляющих чисел
    std::string params_filename;
    // Направляющие числа
    std::vector<DirectionalNumbersParams> dir_num_params;

	// Номер последней сгенерированной точки
    uint32_t current_point_number;
    // Последняя сгенерированная точка
    PointUnsigned last_generated_point;

};
