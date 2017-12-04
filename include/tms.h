#ifndef MINIMIZATION_TMS_H
#define MINIMIZATION_TMS_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>

void GENIN2(unsigned long dimen_, unsigned long seqlen_, unsigned dots_, double (*fun)(std::vector < double > &),
            std::vector < std::vector < double > > &bounds);

void SETFLD(int qin);

void PLYMUL(int *pa, int *pb, int *pc);

void CALCV(int *px, int *b, int *v, int maxv);

void CALCC2();

void INLO2(unsigned long dim, int skip);

void GOLO2(double *quasi, std::vector < std::vector < double > > &bounds);

#endif //MINIMIZATION_TMS_H
