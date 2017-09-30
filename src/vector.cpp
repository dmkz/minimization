#include "vector.h"

// ----- Методы класса "Вектор" -----

// Оператор присваивания: 
Vector& Vector::operator=(const Vector& other) {
	// использовать оператор присваивания для std::vector
}

// Оператор перемещения временного объекта:
Vector& Vector::operator=(Vector&& other) {
	// использовать реализованную семантику перемещения для std::vector
}
	
// Оператор доступа к компоненте вектора:
double& Vector::operator[](const int index) {
	// использовать соответствующий метод контейнера std::vector
}

// Оператор константного доступа к компоненте вектора:
const double& Vector::operator[](const int index) const {
	// использовать соответствующий метод контейнера std::vector
}
	
// Проверка на пустоту:
bool Vector::empty() const {
	// использовать соответствующий метод контейнера std::vector
}

// Получение итераторов на начало и конец вектора:
std::vector<double>::iterator Vector::begin() {
	// использовать соответствующий метод контейнера std::vector
}

std::vector<double>::iterator Vector::end() {
	// использовать соответствующий метод контейнера std::vector
}

const std::vector<double>::iterator Vector::begin() const {
	// использовать соответствующий метод контейнера std::vector
}

const std::vector<double>::iterator Vector::end() const {
	// использовать соответствующий метод контейнера std::vector
}

// Возвращение размерности вектора
int Vector::size() const {
	// использовать соответствующий метод контейнера std::vector с приведением к int
}
	
// Прибавление к вектору содержимое другого вектора (покомпонентно):
Vector& Vector::operator+=(const Vector& other) {
	// реализовать, не забыть вернуть ссылку на текущий объект
}

// Вычитание из вектора содержимого другого вектора (покомпонентно):
Vector& Vector::operator-=(const Vector& other) {
	// реализовать, не забыть вернуть ссылку на текущий объект
}

// Умножение всего вектора на число с изменением *this:
Vector& Vector::operator*=(const double value) {
	// реализовать, не забыть вернуть ссылку на текущий объект
}

// Деление всего вектора на число с изменением *this:
Vector& Vector::operator/=(const double value) {
	// реализовать, не забыть вернуть ссылку на текущий объект
}

// Унарный минус - получение противоположного вектора:
Vector Vector::operator-() const {
	// использовать умножение на -1
}

// ----- Дружественные функции -----
	
// Умножение числа на вектор:
Vector operator*(const double value, Vector vect) {
	// реализовать с использованием *= и семантикой перемещения
} 

// Умножение вектора на число:
Vector operator*(Vector vect, const double value) {
	// реализовать с использованием *= и семантикой перемещения
}

// Деление вектора на число:
Vector operator/(const Vector vect, const double value) {
	// реализовать с использованием *= на обратное и семантикой премещения
}

// Сложение векторов:
Vector operator+(Vector left, const Vector& right) {
	// реализовать при помощи += и семантикой перемещения
}

// Разность векторов:
Vector operator-(Vector left, const Vector& right) {
	// реализовать при помощи -= и семантикой перемещения
}

// Сравнение двух векторов
bool operator==(const Vector& left, const Vector& right) {
	// реализовать, используя сравнение std::vector
}

// Вывод вектора в поток:
std::ostream& operator<<(std::ostream& os, const Vector& vect) {
	// реализовать
} 

// Скалярное произведение векторов:
double dot(const Vector& left, const Vector& right) {
	// реализовать
}