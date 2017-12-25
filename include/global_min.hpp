#pragma once

#include "sobolseqgenerator.h"
#include "math.hpp"
#include "nesterov.hpp"
#include "hessian_free.hpp"
#include "bfgs.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <thread>
#include <mutex>

// Первый этап: вычисление значений функции в точках сетки:
std::vector<std::pair<ld, Vector>>
calc_f_with_threads(Function f, const std::vector<Vector> & inData);

// Второй этап: запуск методов локальной минимизации в попытках улучшить результат: 
std::vector<std::pair<ld, Vector>>
find_local_mins_with_threads(Function f, const std::vector<std::pair<ld, Vector>>& inData);

// Основная функция поиска абсолютных минимумов:
std::vector<std::pair<ld, Vector>>
find_absmin(Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max);