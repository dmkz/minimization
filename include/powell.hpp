#pragma once

// #include "math.hpp"
// #include <cmath>  
// #include <math.h>
// #include <limits>

// Авторы: Бабичев Денис (теория), Бессонов Трофим (тестирование, теория), Данилов Алексей (реализация), Киселев Николай (реализация)

#include <fstream>
#include <complex>
#include <iostream>
#include <cmath>
#include "math.hpp"
#include "StopCondition.hpp"

Real brent(const Real ax, const Real bx, const Real cx, Real f(const Real, int, Vector *, Vector *, Function),
	const Real tol, int ncom, Vector *pcom_p, Vector *xicom_p, Function);

Real f1dim(const Real x, int ncom, Vector *pcom_p, Vector *xicom_p, Function func);

std::pair<Vector, Vector> linmin(Vector &pInit, Vector &xiInit, Function func);

std::pair<Vector, Vector> mnbrak(Real axInit, Real bxInit, Real func(const Real, int, Vector *, Vector *, Function), int ncom, Vector *pcom_p, Vector *xicom_p, Function);
      
IterationData powell(Function func, Vector p, const StopCondition& stop_condition = default_stop_condition);
// func - указатель на целевую функцию
// p - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации