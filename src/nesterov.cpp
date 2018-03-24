#include "nesterov.hpp"

// Метод Нестерова
// Авторы: Петрухина Светлана, Кулага Иван
IterationData nesterov(Function f, Vector startingPoint, const StopCondition& stop_condition) {
// f - указатель на целевую функцию
// startingPoint - начальное приближение
// stop_condition - критерий остановы
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
		//for (int iter=0; iter < 100; ++iter) {
        while (true) {
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

        iter_data.next(xNext, f(xNext));
	} while (!stop_condition(iter_data));
    return iter_data;
}
