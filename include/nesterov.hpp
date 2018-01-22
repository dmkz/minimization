#pragma once

#include "math.hpp"

// Метод Нестерова (возвращается результат - точка минимума и количество сделанных итераций)
// Авторы: Петрухина Светлана, Кулага Иван
std::pair<Vector, int> nesterov(Function f, Vector startingPoint, int iter_limit = 100);