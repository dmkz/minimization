#pragma once

#include <vector>
#include <string>
#include "tmsnet.h"

typedef Point<uint32_t> PointUnsigned;

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

    std::string params_filename;

    // Последняя сгенерированная одиночная точка
    PointUnsigned last_generated_point;
};
