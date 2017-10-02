#include "matrix.h"
#include <stdexcept> // Для выбрасывания исключений std::runtime_error("Поясняющий текст")

// ----- Методы класса Матриц -----

// Конструктор квадратной матрицы фиксированного размера + конструктор по-умолчанию
Matrix::Matrix(int size, double value)
	: data(size, Vector(size, value)) 
{ }

// Конструктор копирования матрицы
Matrix::Matrix(const Matrix& other)
	: data(other.data)
{ }

// Конструктор перемещения матрицы
Matrix::Matrix(Matrix&& other)
	: data(std::move(other.data))
{ }

// Конструктор копирования стандартного контейнера std::vector
Matrix::Matrix(const std::vector<std::vector<double>>& other)
	: data(other.size())
{ 
	// реализовать копирование каждой вектор-строки в цикле отдельно
	// проверить, что размерность столбцов совпадает с размерностью строк (матрица квадратная)
	// если это не так - выбросить исключение при помощи std::runtime_error()
}

// Конструктор перемещения стандартного контейнера std::vector
Matrix::Matrix(std::vector<std::vector<double>>&& other) 
	: data(other.size())
{ 
	// реализовать перемещение каждой вектор-строки в цикле отдельно
	// проверить, что размерность столбцов совпадает с размерностью строк (матрица квадратная)
	// если это не так - выбросить исключение при помощи std::runtime_error()
}

// Оператор присваивания
Matrix& Matrix::operator=(const Matrix& other) {
	// использовать присвоение у data
	// вернуть ссылку на this
}

// Оператор перемещения
Matrix& Matrix::operator=(Matrix&& other) {
	// использовать перемещение у data
	// вернуть ссылку на this
}

// Размер матрицы и проверка на пустоту:
int Matrix::size() const {
	// использовать получение размера у data
}

bool Matrix::empty() const {
	// использовать проверку на пустоту у data
}

// Доступ к итераторам:
std::vector<Vector>::iterator Matrix::begin() {
	// использовать метод у data
}

std::vector<Vector>::iterator Matrix::end() {
	// использовать метод у data
}
std::vector<Vector>::const_iterator Matrix::begin() const {
	// использовать метод у data
}
std::vector<Vector>::const_iterator Matrix::end() const {
	// использовать метод у data
}

// Доступ к отдельной строке при помощи []:
Vector& Matrix::operator[](const int index) {
	// использовать метод [] у data
}

const Vector& Matrix::operator[](const int index) const {
	// использовать метод [] у data
}

// Доступ к отдельному элементу при помощи ():
double& Matrix::operator()(const int row, const int col) {
	// использовать метод [] у data и у вложенного вектора
}

const double& Matrix::operator()(const int row, const int col) const {
	// использовать метод [] у data и у вложенного вектора
}

// Прибавление, вычитание другой матрицы:
Matrix& Matrix::operator+=(const Matrix& other) {
	// проверить на совпадение размерности - выбросить исключение
	// вызвать += у каждой строки
	// вернуть ссылку на this	
}

Matrix& Matrix::operator-=(const Matrix& other) {
	// проверить на совпадение размерности - выбросить исключение
	// вызвать -= у каждой строки
	// вернуть ссылку на this	
}

// Умножение, деление матрицы на число:
Matrix& Matrix::operator*=(const double value) {
	// вызвать *= у каждой строки
	// вернуть ссылку на this	
}

Matrix& Matrix::operator/=(const double value) {
	// проверить деление на ноль - выбросить исключение
	// вызвать /= у каждой строки
	// вернуть ссылку на this	
}

// Сравнение двух матриц:
bool operator==(const Matrix& left, const Matrix& right) {
	// вызвать сравнение у data
}

// Сложение, вычитание двух матриц:
Matrix operator+(Matrix left, const Matrix& right) {
	// использовать += с семантикой перемещения
}

Matrix operator-(Matrix left, const Matrix& right) {
	// использовать -= с семантикой перемещения
}

// Умножение матрицы на матрицу и на вектор справа и слева:
Matrix operator*(const Matrix& left, const Matrix& right) {
	// проверить на совпадение размерностей - выбросить исключение
	// перемножить матрицы во временный объект обычным умножением "каждая строка на аждый столбец"
	// переместить временный объект в результат
}

Vector operator*(const Matrix& matr, const Vector& vect) {
	// проверить на совпадение размерностей - выбросить исключение
	// перемножить матрицу на вектор во временный объект обычным умножением "каждая строка на столбец"
	// переместить временный объект в результат
}

Vector operator*(const Vector& vect, const Matrix& matr) {
	// проверить на совпадение размерностей - выбросить исключение
	// перемножить вектор на матрицу во временный объект обычным умножением "строка на каждый столбец"
	// переместить временный объект в результат
}

// Сравнение размерностей:
bool same_dimension(const Matrix& left, const Matrix& right) {
	// сравнить размеры
}

bool same_dimension(const Matrix& matr, const Vector& vect) {
	// сравнить размеры
}

bool same_dimension(const Vector& vect, const Matrix& matr) {
	// сравнить размеры
}

// Оператор вывода в поток:
std::ostream& operator<<(std::ostream& os, const Matrix& matr) {
	// вывести каждую строку матрицы отдельно с перемносом строки
	// вернуть ссылку на поток вывода
}
