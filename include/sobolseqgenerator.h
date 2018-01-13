#pragma once

#include <vector>
#include <string>
#include "tmsnet.h"

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
    int Init(uint32_t _N, uint32_t _D, std::string dir_file);
    PointReal GeneratePoint();

private:
    uint32_t N;
    uint32_t D;
    uint32_t current_point_number;

    // L = максимальное число нужных битов
    uint32_t L;
    uint32_t C;

    // Параметры направляющих чисел
    std::string params_filename;
    std::vector<DirectionalNumbersParams> dir_num_params;

    // Последняя сгенерированная одиночная точка
    PointUnsigned last_generated_point;

};
