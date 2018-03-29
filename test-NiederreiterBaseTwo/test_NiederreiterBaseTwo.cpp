/*
 * Авторы: Бураханова А., Золкин А., Казарян М., Хохлов С., Малоенко С.
 */

#include "NiederreiterBaseTwo.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string CurrentDateTime() {
    time_t     now = std::time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
    
    return buf;
}

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

void test_comp_sum_range (int start, int finish) {
// Тест на покомпонентную сумму с вычитанием среднего значения. Все правильно, если сумма равна нулю
// - протестировано до 24 размерности, работает
    std::cout << "-- Test Components Sum diff Mean equal zero? from dim = " << start << " to " << finish << std::endl;
    for (int dim_num = start; dim_num <= finish; dim_num++) {

        std::vector<bool> matrix(1LL << dim_num);

        std::vector<std::vector<long double>> table(dim_num, std::vector<long double>(1LL << dim_num));
        
        NiederreiterBaseTwo generator;
        generator.Init(dim_num);
        
        for (int64_t i = 0; i < (1LL << dim_num); i++) {
            auto r = generator.GeneratePoint();
            for (int j = 0; j < dim_num; ++j) {
                table[j][i] = r[j];
            }
        }

        long double max_abs_sum = -1;
        for (int i = 0; i < dim_num; ++i) {
            long double x_avg = avg(table[i]);
            long double sum = 0;
            for (const auto & it : table[i]) {
                sum += (it - x_avg);
            }
            max_abs_sum = std::max(max_abs_sum, std::abs(sum));
        }

        std::cout << "\tdim = " << std::setw(2) << dim_num << ": max_abs_comp_sum_diff_mean(x) = " << max_abs_sum << (max_abs_sum < 1e-8 ? " OK, EQUAL ZERO" : "BAD" ) << std::endl;
    }
    return;
}

void test_comp_uniq_range (int start, int finish) {
// Тест на покомпонентную уникальность
// протестированно от 2 до 24 включительно. Вердикт - ОК
    std::cout << "-- Test Components Unique from dim = " << start << " to " << finish << std::endl;
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

        bool is_success = true;

        for (int i = 0; i < dim_num; ++i) {
            auto prev_size = table[i].size();
            std::sort(table[i].begin(), table[i].end());
            table[i].erase(std::unique(table[i].begin(), table[i].end()), table[i].end());
            auto curr_size = table[i].size();
            if (prev_size > curr_size) {
                is_success = false;
                break;
            }
        }

        std::cout << "\tdim = " << std::setw(2) << dim_num << (is_success ? ": ok, projections for each component are unique!" : ": bad") << std::endl;
    }
    return;
}

int main() {
    std::cout << "Runned!" << std::endl;
    test_comp_sum_range(2, 20);  std::cout << std::endl << std::endl;
    test_comp_uniq_range(2, 20); std::cout << std::endl << std::endl;
    test_corr_range(2, 20);      std::cout << std::endl << std::endl;
    return 0;
}