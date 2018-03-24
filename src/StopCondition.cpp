// Файл StopCondition.cpp
// Авторы: Козырев Дмитрий, Кондратов Юрий

#include "StopCondition.hpp"

// Конструктор по-умолчанию:
IterationData::IterationData() {
    f_prev = f_curr = iter_counter = 0;
}

// Переход к следующей итерации
void IterationData::next(const Vector& x_next, const Real f_next) {
    f_prev = f_curr;
    x_prev = x_curr;
    x_curr = x_next;
    f_curr = f_next;
    iter_counter++;
}

// Условие останова метода по-умолчанию
bool default_stop_condition(const IterationData& data) { 
    return data.iter_counter >= 100 || std::abs(data.f_curr - data.f_prev) < 1e-8;
}


