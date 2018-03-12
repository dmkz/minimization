#include "math.hpp"

/*
    Библиотека со всей математикой. 
    Автор: Дмитрий Козырев
*/

/*
// Математические константы:
const Real M_E = std::exp(1.0L);
const Real M_PI = std::acos(-1.0L);
const Real M_PI_2 = M_PI / 2.0L;
*/

// Операторы вывода в поток:
std::ostream& operator<<(std::ostream& os, const Vector& v) {
	for (auto& it : v) {
		os << std::setprecision(8) << std::fixed << std::setw(16) << it;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
	for (auto& it : m) {
		os << it << std::endl;
	}
	return os;
}

// Операция умножения вектора на число и различные ее вариации:
Vector& operator*=(Vector& v, const Real value) {
	for (auto & it : v) {
		it *= value;
	}
	return v;
}

Vector operator*(const Real value, Vector v) {
	return v *= value;
}

Vector operator*(Vector v, const Real value) {
	return v *= value;
}

// Унарный минус для вектора:
Vector operator-(Vector v) {
    return v *= -1;   
}

// Сложение и вычитание векторов:
Vector& operator+=(Vector & v1, const Vector& v2) {
    assert(v1.size() == v2.size());
	for (int i = 0; i < (int)v1.size(); ++i) {
		v1[i] += v2[i];
	}
	return v1;
}

Vector& operator-=(Vector & v1, const Vector& v2) {
	return v1 += -v2;
}

Vector operator+(Vector v1, const Vector& v2) {
	return v1 += v2;
}

Vector operator-(Vector v1, const Vector& v2) {
	return v1 -= v2;
}

// Скалярное произведение векторов:
Real dot(const Vector& v1, const Vector& v2) {
    assert(v1.size() == v2.size());
	Real sum = 0;
	for (int i = 0; i < (int)v1.size(); ++i) {
		sum += v1[i] * v2[i];
	}
	return sum;
}

// Норма вектора:
Real norm(const Vector& v) {
    return std::sqrt(dot(v, v));
}

// Умножение матрицы на вектор:
Vector operator*(const Matrix& m, const Vector& v) {
	int nRows = (int)m.size();
	Vector ans(nRows);
	for (int i = 0; i < nRows; ++i) {
		ans[i] = dot(m[i], v);
	}
	return ans;
}

// Проверка равенства вектора нулю
bool is_zero(const Vector& v) {
    for (const auto & it : v) {
        if (std::abs(it) > COMPARE_EPS) {
            return false;
        }
    }
    return true;
}

// Базисный вектор в пространстве R^n:
// на i-ом месте 1, на всех остальных нули
Vector id_vect(int size, int i) {
	Vector answer(size);
	answer[i] = 1;
	return answer;
}

// Явное взятие градиента от функции:
Vector grad(Function f, const Vector& point, const Real h) {
	int n = (int)point.size();
	Vector right = point;
	Vector left = point;
	Vector answer(n);
	for (int i = 0; i < n; ++i) {
		// Придаем приращение i-му аргументу:
        right[i] += h;
		left[i]  -= h;
		// По формуле центральных разностей с погрешностью O(h^2):
		answer[i] = (f(right) - f(left)) / (2*h);
		// Убираем приращение обратно:
		right[i] -= h;
		left[i]  += h;
	}
	return answer;
}

// Явное вычисление матрицы Гессе (слишком затратно, годится только для проверки результатов)
// в точке x с шагом h и погрешностью O(h)
Matrix hess(Function f, const Vector& x, const Real h) {
	int n = (int)x.size();
	Matrix ans(n, Vector(n));
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			ans[i][j] = (
				f(x+(id_vect(n, i)+id_vect(n, j))*h)
				-f(x+id_vect(n, i)*h)
				-f(x+id_vect(n, j)*h)
				+f(x)) / h / h;
		}
	return ans;
}
