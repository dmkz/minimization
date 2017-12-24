#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "sobolseqgenerator.h"

using namespace std;

int SobolSeqGenerator::Init()
{
    return Init(1, 1, "sobol_net_default.txt");
}

int SobolSeqGenerator::Init(long unsigned int _N, long unsigned int _D, string dir_file)
{
    N = _N;
    D = _D;
    current_point_number = -1;
    L = (unsigned)ceil(log((double)N)/log(2.0));

    ifstream infile(dir_file,ios::in);
    if (!infile) {
      cout << "Не найден файл содержащий направляющие числа!\n";
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
        cout << "Генератор сеток Соболева некорректно инициализирован! N = " << N << ", D = " << D;
        return PointReal();
    }

    ifstream infile(params_filename, ios::in);
    if (!infile) {
      cout << "Не найден файл содержащий направляющие числа!\n";
      return -1;
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');

    if(current_point_number == N - 1)
    {
        cout << "Генератор уже сгенерировал все N точек!";
        return PointReal();
    }

    current_point_number += 1;

    if(current_point_number == 0)
    {
        last_generated_point = PointUnsigned(D, vector<unsigned>(D, 0));
        return PointReal(D, vector<Real>(D, 0));
    }

    PointUnsigned result_point = PointUnsigned(D);

    C = 1;
    unsigned value = current_point_number - 1;
    while (value & 1) {
        value >>= 1;
        C++;
    }

    // Вычислить направляюoее число V, умноженное на pow(2,32)
    unsigned V_first = 1 << (32-C);

    // Вычислить первую координату, умноженную на pow(2,32)
    result_point.coordinate[0] = last_generated_point.coordinate[0] ^ V_first;

    // ----- Вычислить остальные координаты -----
    for (unsigned j = 1; j < D; j++) {

        // Чтение параметров из файла
        unsigned d, s;
        unsigned a;
        infile >> d >> s >> a;
        unsigned *m = new unsigned [s + 1];
        for (unsigned i = 1; i <= s; i++)
        {
            infile >> m[i];
        }

        // Вычислить направляющие числа V, умноженные на pow(2,32)
        unsigned *V = new unsigned [L + 1];
        if (L <= s)
        {
            for (unsigned i = 1; i <= L; i++)
            {
                V[i] = m[i] << (32 - i);
            }
        }
        else
        {
            for (unsigned i = 1; i <= s; i++)
            {
                V[i] = m[i] << (32 - i);
            }

            for (unsigned i = s + 1; i <= L; i++)
            {
                V[i] = V[i-s] ^ (V[i-s] >> s);
  	            for (unsigned k = 1; k <= s-1; k++)
                {
                    V[i] ^= (((a >> (s - 1 - k)) & 1) * V[i-k]);
                }
            }
        }

        // Вычислить j-ю координату, умноженную на pow(2,32)
        result_point.coordinate[j] = last_generated_point.coordinate[j] ^ V[C];
    }

    infile.close();

    last_generated_point = result_point;

    PointReal final_result = PointReal(D);
    for (unsigned j = 0; j < D; j++)
    {
        final_result.coordinate[j] = (Real) result_point.coordinate[j] / pow(2.0, 32);
    }

    return final_result;
}
