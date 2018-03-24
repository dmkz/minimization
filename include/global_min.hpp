#pragma once

#include "sobolseqgenerator.hpp"
#include "math.hpp"
#include "nesterov.hpp"
#include "hessian_free.hpp"
#include "bfgs.hpp"
#include "dfp.hpp"
#include "powell.hpp"
#include "StopCondition.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <thread>
#include <mutex>

// Первый этап: вычисление значений функции в точках сетки:
// Автор: Козырев Дмитрий
std::vector<std::pair<Real, Vector>>
calc_f_with_threads(Function f, const std::vector<Vector> & inData);

// Второй этап: запуск методов локальной минимизации в попытках улучшить результат: 
// Автор: Козырев Дмитрий
std::vector<std::pair<Real, Vector>>
find_local_mins_with_threads(Function f, const StopCondition& stop_condition, const std::vector<std::pair<Real, Vector>>& inData);

// Основная функция поиска абсолютных минимумов:
// Автор: Козырев Дмитрий
std::vector<std::pair<Real, Vector>>
find_absmin(Function f, const StopCondition& stop_condition, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max);