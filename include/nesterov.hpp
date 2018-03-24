#pragma once

#include "math.hpp"
#include "StopCondition.hpp"

// Метод Нестерова (возвращается результат - точка минимума и количество сделанных итераций)
// Авторы: Петрухина Светлана, Кулага Иван
IterationData nesterov(Function f, Vector startingPoint, const StopCondition& stop_condition = default_stop_condition);
// f - указатель на целевую функцию
// startingPoint - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации