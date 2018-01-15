#pragma once

#include "math.hpp"
#include "stopCondition.h"

// Метод Нестерова (возвращается результат - точка минимума и количество сделанных итераций)
std::pair<Vector, int> nesterov(Function f, Vector startingPoint, stopCondition condition);
