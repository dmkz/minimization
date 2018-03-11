/*
 * Авторы изменений: Бураханова А., Золкин А., Хохлов С.
 */

#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include "tmsnet.hpp"

class NiederreiterGenerator : public TMSNet
{
public:
    // Инициализация таблицы сложений, умножений, вычитаний
    void SETFLD(int qin);

    // Метод для умножения полиномов
    void PLYMUL(int *pa, int *pb, int *pc);

    /*
      C помощью этой программы вычисляем значения констант V(J,R)

      px соответствующего неприводимого полинома для  текущего измерения
      На входе в подпрограмму, степень B  определяет параметр J, как deg(B) = Е*(J-1), где Е - степень px.
      MAXV дает размерность массива V.

      результат выполнения - B умноженный на px, его степень на данный момент Е*J.
      V содержит требуемыые значения.
      подпрограмма PLYMUL используется для умножения полиномов
    */
    void CALCV(int *px, int *b, int *v, int maxv);

    /*
    C помощью этой программы вычисляем значения константы  C(I,J,R).
    для каждого значения I, мы сначала рассчитываем все соответствующие значения  C : они хранятся в массиве CI.
    Все  значения C  0 или 1. Далее мы помещаем значения в массив  СJ, таким образом,  CJ(I,R) имеет значения С
    для указанных значений I и R, и для каждого J от 1 до NBITS
    И когда значения CJ(I,R) расчитаны, мы возвращаемсяк этому массиву в вызываемой программе
    */
    void CALCC2();

    /*
    Эта подпрограмма вычисляет значения нидеррейтера C(I,J,R) , вызывая CALCC2 и выполняет другие необходимые инициализации перед вызовом GOLO2
    dim - измерение последовательностей, которые должны быть созданы.
    skip - количество значений, которые должны быть отброшены в начале последовательности
    */
    void INLO2(unsigned long dim, int skip);

    /*
    Эта функция создает новый квази-случайный вектор на каждый вызов программы
    bounds - границы области, в которой генерируются точки.
    */
    void GOLO2(double *quasi, std::vector < std::vector < double > > &bounds);

    /*
    C помощью этой программы тестируем на точность численное интегрирование
    C помощью малых отклонений двоичных последовательностей Нидеррайтера (1988), как это реализовано в INLO2, GOLO2.
    Различные возможные проверки интегралов GENIN2 генерирует только  последовательности с основанием 2.
    dime_ - размерность пространвства
    sqlen_ - колличество точе
    m_ - колличество минимальныйх точек для запуска функций минимизации.
    bounds - границы области
    */
    void GENIN2(unsigned long dimen_, unsigned long seqlen_, unsigned dots_, double (*fun)(std::vector < double > &),
        std::vector < std::vector < double > > &bounds);
    //Функция сравнения точек для сортировки
    static int cmp(const void *a, const void *b);

    //Сгенерировать точки в файл
    int GeneratePointsToFile(int argc, char *argv[]);

    // COMM
    static const int MAX_DIMENSION = 20;
    static const int NUMBER_OF_BITS = 31;
    double recip;
    unsigned long dimen;

    // COMM2
    int cj[MAX_DIMENSION][NUMBER_OF_BITS/*!*/];
    int count2;
    unsigned long dimen2;
    int nextq2[MAX_DIMENSION];

    // FIELD
    static const int MAX_Q_VALUE = 50;
    static const int MAX_POLY_DEGREE = 50;
    static const int CURRENT_DEGREE_VALUE = 0;
    int p, q, add[MAX_Q_VALUE/*!*/][MAX_Q_VALUE/*!*/], mul[MAX_Q_VALUE/*!*/][MAX_Q_VALUE/*!*/],
            sub[MAX_Q_VALUE/*!*/][MAX_Q_VALUE/*!*/];

    std::vector< std::vector<double> > vResult;
    std::vector<double> vector_;
    int ResultCounter = 0;
};
