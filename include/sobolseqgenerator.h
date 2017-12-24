#pragma once

#include <vector>
#include <string>
#include "tmsnet.h"

using namespace std;

typedef Point<unsigned> PointUnsigned;

class SobolSeqGenerator : public TMSNet
{
public:
    int Init();
    int Init(long unsigned int _N, long unsigned int _D, string dir_file);
    PointReal GeneratePoint();

private:
    long unsigned int N;
    long unsigned int D;
    long unsigned int current_point_number;

    // L = максимальное число нужных битов
    unsigned L;
    unsigned C;

    string params_filename;

    // Последняя сгенерированная одиночная точка
    PointUnsigned last_generated_point;
};
