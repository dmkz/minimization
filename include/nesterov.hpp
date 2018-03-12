#pragma once

#include "math.hpp"
#include "iteration_object.hpp"

// Метод Нестерова (возвращается результат - точка минимума и количество сделанных итераций)
// Авторы: Петрухина Светлана, Кулага Иван
void nesterov(Function f, Vector startingPoint, BasicIterationObject* iter_object);
// f - указатель на целевую функцию
// startingPoint - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации