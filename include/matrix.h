#pragma once

#include "vector.h" 

#include <vector> // для контейнера std::vector
#include <initializer_list> // для списка инициализации std::initializer_list
#include <algorithm> // для семантики перемещения std::move
#include <iostream> // для потока вывода

// Тип матрицы инициализации (например, вида {{a, b},{c,d}}):
typedef std::initializer_list<Vector> Init_Matr;

// Проверка на то, что матрица инициализации является квадратной с выбрасыванием исключения:
Init_Matr throw_if_not_square_matrix(Init_Matr M);

class Matrix;

// Сравнение двух матриц:
bool operator==(const Matrix& left, const Matrix& right);
	
// Сложение, вычитание двух матриц:
Matrix operator+(Matrix left, const Matrix& right);
Matrix operator-(Matrix left, const Matrix& right);
	
// Умножение матрицы на матрицу и на вектор справа и слева:
Matrix operator*(const Matrix& left, const Matrix& right);
Vector operator*(const Matrix& matr, const Vector& vect);
Vector operator*(const Vector& vect, const Matrix& matr);
	
// Сравнение размерностей:
bool same_dimension(const Matrix& left, const Matrix& right);
bool same_dimension(const Matrix& matr, const Vector& vect);
bool same_dimension(const Vector& vect, const Matrix& matr);

// Вывод матрицы в поток:
std::ostream& operator<<(std::ostream& os, const Matrix& matr);

class Matrix { // Класс квадратных матриц
private:
	std::vector<Vector> data;
public:
	// Конструктор квадратной матрицы фиксированного размера + конструктор по-умолчанию
	Matrix(int size = 0, double value = 0);
	
	// Конструктор копирования матрицы
	Matrix(const Matrix& other);
	
	// Конструктор перемещения матрицы
	Matrix(Matrix&& other);
	
	// Конструктор перемещения матрицы инициализации
	Matrix(Init_Matr other);
	
	// Оператор присваивания
	Matrix& operator=(const Matrix& other);
	
	// Оператор перемещения
	Matrix& operator=(Matrix&& other);
	
	// Размер матрицы и проверка на пустоту:
	int size() const;
	bool empty() const;
	
	// Доступ к итераторам:
	std::vector<Vector>::iterator begin();
	std::vector<Vector>::iterator end();
	std::vector<Vector>::const_iterator begin() const;
	std::vector<Vector>::const_iterator end() const;
	
	// Доступ к отдельной строке при помощи []:
	Vector& operator[](const int index);
	const Vector& operator[](const int index) const;
	
	// Доступ к отдельному элементу при помощи ():
	double& operator()(const int row, const int col);
	const double& operator()(const int row, const int col) const;
	
	// Прибавление, вычитание другой матрицы:
	Matrix& operator+=(const Matrix& other);
	Matrix& operator-=(const Matrix& other);
	
	// Умножение, деление матрицы на число:
	Matrix& operator*=(const double value);
	Matrix& operator/=(const double value);
	
	// Сравнение двух матриц:
	friend bool operator==(const Matrix& left, const Matrix& right);
	
	// Сложение, вычитание двух матриц:
	friend Matrix operator+(Matrix left, const Matrix& right);
	friend Matrix operator-(Matrix left, const Matrix& right);
	
	// Умножение матрицы на матрицу и на вектор справа и слева:
	friend Matrix operator*(const Matrix& left, const Matrix& right);
	friend Vector operator*(const Matrix& matr, const Vector& vect);
	friend Vector operator*(const Vector& vect, const Matrix& matr);
	
	// Сравнение размерностей:
	friend bool same_dimension(const Matrix& left, const Matrix& right);
	friend bool same_dimension(const Matrix& matr, const Vector& vect);
	friend bool same_dimension(const Vector& vect, const Matrix& matr);
	
	// Оператор вывода в поток:
	friend std::ostream& operator<<(std::ostream& os, const Matrix& matr);
};