// Файл StopCondition.hpp
// Авторы: Козырев Дмитрий
#pragma once

#include <functional>

#include "math.hpp"

struct IterationData;
typedef std::function<bool(const IterationData&)> StopCondition;

// Структура, хранящая информацию о текущей итерации:
struct IterationData {
    Vector x_prev;            // Предущыщая точка
    Real   f_prev;            // Значение функции в предыдущей точке
    Vector x_curr;            // Текущая точка    
    Real   f_curr;            // Значение функции в текущей точке
    int    iter_counter;      // Текущее количество итераций
    std::string method_title; // Название метода
    
    IterationData(); // Конструктор по-умолчанию
    
    void next(const Vector& x_next, const Real f_next); // Переход к следующей итерации
};

// Условие остановы итерации, применяемое по-умолчанию:
bool default_stop_condition(const IterationData&);