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
// Автор: Бураханова Алена
struct GlobalTestData
{
	GlobalTestData() : x()
	{
		method = std::string("");
	}
	
	GlobalTestData(const Vector& x_, std::string method_)
	{
		x = x_;
		method = method_;
	}
	
	Vector x;
	std::string method;
};

// Структуры для определения operator() для использования в шаблонах(std::max_element, std::sort, и.т.д.)
// Автор: Бураханова Алена
struct comparePairRealVector
{
    inline bool operator() (const std::pair<Real, Vector>& i, const std::pair<Real, Vector>& j)
    {
		return i.first < j.first;
    }
};


// Автор: Бураханова Алена
struct comparePairRealGlobalTestData
{
    inline bool operator() (const std::pair<Real, GlobalTestData>& i, const std::pair<Real, GlobalTestData>& j)
    {
		return i.first < j.first;
    }
};

// Первый этап: вычисление значений функции в точках сетки:
// Автор: Козырев Дмитрий
std::vector<std::pair<Real, Vector>>
calc_f_with_threads(Function f, const std::vector<Vector> & inData);

// Второй этап: запуск методов локальной минимизации в попытках улучшить результат: 
// Автор: Козырев Дмитрий, Бураханова Алена
std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_local_mins_with_threads(Function f, const StopCondition& stop_condition, const std::vector<std::pair<Real, Vector>>& inData);

// Основная функция поиска абсолютных минимумов:
// Автор: Козырев Дмитрий, Бураханова Алена
std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_absmin(Function f, const StopCondition& stop_condition, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max);