#pragma once

#include <vector>
#include <string>
#include "tmsnet.h"

using namespace std;

typedef Point<unsigned> PointUnsigned;

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
    int Init(long unsigned int _N, long unsigned int _D, string dir_file);
    PointReal GeneratePoint();

private:
    long unsigned int N;
    long unsigned int D;
    long unsigned int current_point_number;

    // L = максимальное число нужных битов
    unsigned L;
    unsigned C;

    // Параметры направляющих чисел
    string params_filename;
    vector<DirectionalNumbersParams> dir_num_params;

    // Последняя сгенерированная одиночная точка
    PointUnsigned last_generated_point;

};
