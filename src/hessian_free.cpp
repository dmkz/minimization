#include "hessian_free.hpp"

// Авторы: Козырев Дмитрий (реализация), Бадави Полина (теория)
IterationData slow_hessian_free(Function f, Vector x, const StopCondition& stop_condition) {
// f - указатель на целевую функцию
// x - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации

    
// Инициализируем начальной точкой структуру данных итерации:
    IterationData iter_data;
    iter_data.x_curr = x;
    iter_data.f_curr = f(x);
    iter_data.iter_counter = 0;
    iter_data.method_title = "Slow Hessian Free";
    
    const int n = (int)x.size();
    do {
        x += conjugade_gradient(hess(f, x), grad(f, x), Vector(n, 0));
        iter_data.next(x, f(x));
    } while (!stop_condition(iter_data));
    return iter_data;
}

// Авторы: Козырев Дмитрий (реализация), Бадави Полина (теория)
IterationData hessian_free(Function f, Vector x, const StopCondition& stop_condition) {
// f - указатель на целевую функцию
// x - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в объекте итерации

// Инициализируем начальной точкой структуру данных итерации:
    IterationData iter_data;
    iter_data.x_curr = x;
    iter_data.f_curr = f(x);
    iter_data.iter_counter = 0;
    iter_data.method_title = "Hessian Free";
    
    const int n = (int)x.size();
    do {
        Vector dx(n, 0);
        auto grad_f_x = grad(f, x);
        auto d = - (hess_prod_vect(f, x, dx) + grad_f_x);
        for (int j = 0; j < n && !is_zero(d); ++j) {
            auto hess_prod_vect_f_x_d = hess_prod_vect(f, x, d);
            Real alpha = -(dot(d, hess_prod_vect(f, x, dx) + grad_f_x) / dot(d, hess_prod_vect_f_x_d));
            dx += alpha * d;
            auto temp = hess_prod_vect(f, x, dx) + grad_f_x;
            d = -temp + dot(temp, hess_prod_vect_f_x_d) / dot(d, hess_prod_vect_f_x_d) * d;
        }
        x += dx;
        iter_data.next(x, f(x));
    } while (!stop_condition(iter_data));
    return iter_data;
}

// Апроксимация умножения градиента функции f в точке х на вектор dx
// Погрешность O(||h||^2*||dx||^3), где h - выбранный шаг дифференцирования
// Авторы: Козырев Дмитрий (реализация), Бадави Полина (теория)
Real grad_prod_vect(Function f, const Vector& x, const Vector& dx) {
    if (is_zero(dx)) {
        return 0;
    }
    // Подбор шага h таким образом, чтобы погрешность составила 1e-8:
    Real h = std::sqrt(1e-8 / std::pow(norm(dx), 3.0));
	return (1 / (2*h)) * (f(x + h * dx) - f(x - h * dx));
}

// Апроксимация умножения матрицы Гессе в точке x на вектор dx 
// Погрешность O(||h||^2*||dx||^3)
// Авторы: Козырев Дмитрий (реализация), Бадави Полина (теория)
Vector hess_prod_vect(Function f, const Vector& x, const Vector& dx) {
    if (is_zero(dx)) {
        return Vector((int)x.size(), 0);
    }
    // Подбор шага h таким образом, чтобы погрешность составила 1e-8:
    Real h = std::sqrt(1e-8 / std::pow(norm(dx), 3.0));
	return (1 / (2*h)) * (grad(f, x+h*dx) - grad(f, x-h*dx));
}

// Авторы: Козырев Дмитрий (реализация), Бадави Полина (теория)
Vector conjugade_gradient(Matrix A, Vector b, Vector x) {
    auto d = - (A * x + b);
    for (int i = 0; i < (int)x.size() && !is_zero(d); ++i) {
        Real alpha = -(dot(d, A*x+b)) / dot(d, A * d);
        x += alpha * d;
        d = -(A * x + b) + (dot(A * x + b, A * d) / dot(d, A * d)) * d;    
    }
    return x;
}