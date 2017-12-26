#include "nesterov.hpp"

// Метод Нестерова (возвращается результат - точка минимума и количество сделанных итераций)
std::pair<Vector, int> nesterov(Function f, Vector startingPoint, int iter_limit) {
    ld ro = 2.0;
	ld teta = 1.1;
	ld alfa = 1.0;

	Vector v = startingPoint;
	Vector x = startingPoint;
	Vector xNext = startingPoint;
	Vector y = startingPoint;
	ld alfaNext = 0.0;
	ld A = 0;

    int iterations = 0;
	while (iterations < iter_limit) {
		while (true) {
			alfaNext = alfa + std::sqrt(alfa*alfa + 2 * alfa*A);
			y = (A / (A + alfaNext)) * x + (alfaNext / (A + alfaNext)) * v;
			xNext = y - alfa * grad(f, y);
			ld temp1 = f(xNext);
			ld temp2 = (f(y) + dot(grad(f, y), xNext - y) + 1 / (2 * alfa) *  std::pow(norm(xNext - y), 2));
			if (temp1 - temp2 <= COMPARE_EPS)
				break;
			alfa = alfa / ro;
		}
		v = v - alfaNext * (grad(f, xNext));
		A += alfaNext;
		alfaNext = teta * alfa;
        ++iterations;
	}
	return {xNext, iterations};
}