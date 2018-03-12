#include "nesterov.hpp"

// Метод Нестерова
// Авторы: Петрухина Светлана, Кулага Иван
void nesterov(Function f, Vector startingPoint, BasicIterationObject* iter_object) {
// f - указатель на целевую функцию
// startingPoint - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации

    Real ro = 2.0;
	Real teta = 1.1;
	Real alfa = 1.0;

	Vector v = startingPoint;
	Vector x = startingPoint;
	Vector xNext = startingPoint;
	Vector y = startingPoint;
	Real alfaNext = 0.0;
	Real A = 0;
    
    // Инициализируем начальной точкой объект контроля итераций:
    iter_object->set_x_curr(x);
    iter_object->set_f_curr(f(x));
    iter_object->set_iter_counter(0);
    iter_object->set_method_title("Nesterov");
    
	do {
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

        iter_object->next_iteration(xNext, f(xNext));
	} while (!iter_object->is_stopped());
}
