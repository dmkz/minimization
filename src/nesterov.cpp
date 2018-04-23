#include "nesterov.hpp"

// Метод Нестерова
// Авторы: Петрухина Светлана, Кулага Иван
IterationData nesterov(Function f, Vector startingPoint, const StopCondition& stop_condition) {
// f - указатель на целевую функцию
// startingPoint - начальное приближение
// stop_condition - критерий останова
// Результат работы метода будет лежать в структуре данных о последней итерации
    
    Real ro = 2.0;
	Real teta = 1.1;
	Real alfa = 1.0;

	Vector v = startingPoint;
	Vector x = startingPoint;
	Vector xNext = startingPoint;
	Vector y = startingPoint;
	Real alfaNext = 0.0;
	Real A = 0;
    
    // Инициализируем начальной точкой структуру данных итерации:
    IterationData iter_data;
    iter_data.x_curr = x;
    iter_data.f_curr = f(x);
    iter_data.iter_counter = 0;
    iter_data.method_title = "Nesterov";
    
	do {
		for (int local_iter=0; local_iter < 30; ++local_iter) {
            // alfaNext^2 = 2 * alfa * (A + alfaNext)
            // alfaNext^2 - 2 * alfa * (A + alfaNext) = 0
            // alfaNext^2 - 2 * alfa * alfaNext - 2 * alfa * A = 0
            // D = 4 * alfa * alfa + 8 * alfa * A = 4*(alfa * alfa + 2 * alfa * A)
            // alfaNext1 = alfa + sqrt(alfa*alfa + 2 * alfa * A)
			alfaNext = alfa + std::sqrt(alfa*alfa + 2 * alfa*A);
			y = (A / (A + alfaNext)) * x + (alfaNext / (A + alfaNext)) * v;
			xNext = y - alfa * grad(f, y);
			Real temp1 = f(xNext);
			Real temp2 = (f(y) + dot(grad(f, y), xNext - y) + 1 / (2 * alfa) *  std::pow(norm(xNext - y), 2));
			if (temp1 - temp2 <= COMPARE_EPS)
				break;
			alfa = alfa / ro;
		}
		v = v - alfaNext * (grad(f, xNext));
		A += alfaNext;
		alfaNext = teta * alfa;
        /*
        if (dot(y-xNext,xNext-x) > 0) {
            xNext = x; v = x; A = 0;
        } */
        iter_data.next(xNext, f(xNext));
        
	} while (!stop_condition(iter_data));
    return iter_data;
}
