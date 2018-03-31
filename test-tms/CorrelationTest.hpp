/*
 * Авторы: Бураханова А., Золкин А., Арцишевский Л.
 */

#pragma once

#include "TMSNetTestElement.hpp"

// Среднее значение вектора:
long double avg(const std::vector<long double>& v) {
    long double sum = 0;
    for (auto & it : v) sum += it;
    return sum / v.size();
}

// Значение линейного коэффициента ккорреляции Пирса:
long double corr(const std::vector<long double>& x, const std::vector<long double>& y) {
    assert(x.size() == y.size());
    auto x_avg = avg(x);
    auto y_avg = avg(y);
    
    long double p  = 0;
    long double q1 = 0;
    long double q2 = 0;
    for (int i = 0; i < (int)x.size(); ++i) {
        p  += (x[i]-x_avg)*(y[i]-y_avg);
        q1 += (x[i]-x_avg)*(x[i]-x_avg);
        q2 += (y[i]-y_avg)*(y[i]-y_avg);
    }
    return p / std::sqrt(q1 * q2);
}

// Максимальная величина корреляции по модулю среди всех компонент:
long double max_abs_corr(const std::vector<std::vector<long double>>& table) {    
    const int dim = (int)table.size();
    long double max = -1;
    
    for (int i = 0; i < dim; ++i) {
        for (int j = i+1; j < dim; ++j) {
            auto c = corr(table[i], table[j]);
            if (max < std::abs(c)) {
                max = std::abs(c);
            }
            max = std::max(max, std::abs(c));
        }
    }
    
    return max;
}

void test_corr_range(int start, int finish) {
// Тест на покомпонентную корреляцию - протестировано до 24 размерности, работает!
    std::cout << "-- Test Components Correlation from dim = " << start << " to " << finish << std::endl;
    for (int dim_num = start; dim_num <= finish; dim_num++) {
        
        std::vector<std::vector<long double>> table(dim_num, std::vector<long double>(1LL << dim_num));
        
        NiederreiterBaseTwo generator;
        generator.Init(dim_num);
        
        for (int64_t i = 0; i < (1LL << dim_num); i++) {
            auto r = generator.GeneratePoint();
            for (int j = 0; j < dim_num; ++j) {
                table[j][i] = r[j];
            }
        }
        
        auto res = max_abs_corr(table);
        
        std::cout << "\tdim = " << std::setw(2) << dim_num << ": max_abs_corr(x,y) = " << std::fixed << std::setprecision(8) << std::setw(12) << res << std::endl;
    }
    return;
}
