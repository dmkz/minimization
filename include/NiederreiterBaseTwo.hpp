/*
 * Авторы: Бураханова А., Золкин А., Казарян М., Хохлов С., Малоенко С.
 * Генератор сетоr Нидеррайтера, использующий поле F2
 * файл NiederreiterGenerator2.hpp
 */
 
//  Local Parameters:
//
//    int64_t MAXDEG, the highest degree of polynomial
//    to be handled. 
//
//    int64_t DIM_MAX, the maximum dimension that will be used.
//
//    int64_t NBITS, the number of bits (not counting the sign) in a
//    fixed-point integer.
//
//    long double RECIP, the multiplier which changes the
//    integers in NEXTQ into the required real values in QUASI.
//
# pragma once
# include "tmsnet.hpp"

# define MAXDEG 250
# define DIM_MAX 63
# define NBITS 63
 
class NiederreiterBaseTwo
{
private:
    const long double RECIP = 1.0 / ( long double ) ( uint64_t(1) << NBITS ); // Коэффициент масштабирования
    int MAXE; // Максимальная степень многочленов
    
    int add[2][2], mul[2][2], sub[2][2]; // Таблицы умножения, сложения, вычитания в поле по модулю 2
    
    int64_t dimension;
    int64_t seed;
    
    int64_t dim_save;
    int64_t seed_save;
    
    int64_t cj[DIM_MAX][NBITS];
    
    int64_t nextq[DIM_MAX];
    
    std::vector<std::vector<int64_t>> irred; // Неприводимые многочлены
    std::vector<int> irred_deg; // Степени неприводимых многочленов
    
    
    
private: // Вспомогательные функции
    void setfld2 ();
    void plymul2 ( 
        int64_t pa_deg, int64_t pa[MAXDEG+1], 
        int64_t pb_deg, int64_t pb[MAXDEG+1], 
        int64_t *pc_deg, int64_t pc[MAXDEG+1] ) const;
        
    void calcc2 ();

    void calcv2 ( int64_t maxv, int64_t px_deg, int64_t px[MAXDEG+1], int64_t *b_deg, int64_t b[MAXDEG+1], int64_t v[] ) const;
    
public:
    // Инициализация:
    int Init();
    int Init(uint32_t dimension);
    int Init(uint32_t dimension, int64_t seed);
    // Генерация одной точки:
    std::vector<Real> GeneratePoint();
};
 