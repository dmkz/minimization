#pragma once

#include "NiederreiterBaseTwo.hpp"
//#include "sobolseqgenerator.hpp"
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

// Структура данных для вывода информации о лучших методах
// Автор: Козырев Дмитрий, Бураханова Алена, Золкин Артем
struct GlobalTestData
{
	GlobalTestData() : x()
	{
		method = std::string("");
	}
	
	GlobalTestData(const Vector& x_, std::string method_, std::string stop_cause_)
	{
		x = x_;
		method = method_;
		stop_cause = stop_cause_;
	}
	
	Vector x;
	std::string method;
	std::string stop_cause;
};

// Структуры для определения operator() для использования в шаблонах(std::max_element, std::sort, и.т.д.)
struct comparePairRealVector
{
    inline bool operator() (const std::pair<Real, Vector>& i, const std::pair<Real, Vector>& j)
    {
		return i.first < j.first;
    }
};


struct comparePairRealGlobalTestData
{
    inline bool operator() (const std::pair<Real, GlobalTestData>& i, const std::pair<Real, GlobalTestData>& j)
    {
		return i.first < j.first;
    }
};

// Первый этап: вычисление значений функции в точках сетки:
std::vector<std::pair<Real, Vector>>
calc_f_with_threads(Function f, const std::vector<Vector> & inData);

// Второй этап: запуск методов локальной минимизации в попытках улучшить результат: 
std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_local_mins_with_threads(Function f, const StopCondition& stop_condition, const std::vector<std::pair<Real, Vector>>& inData);

// Основная функция поиска абсолютных минимумов:
std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_absmin(Function f, const StopCondition& stop_condition, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max);