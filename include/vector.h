#pragma once

#include <vector> // для контейнера std::vector
#include <algorithm> // для семантики перемещения std::move

class Vector;

Vector operator*(const double value, Vector vect); // умножение числа на вектор
Vector operator*(Vector vect, const double value); // умножение вектора на число
Vector operator/(Vector vect, const double value); // деление вектора на число
Vector operator+(Vector left, const Vector& right); // сложение векторов
Vector operator-(Vector left, const Vector& right); // разность векторов
bool operator==(const Vector& left, const Vector& right); // сравнение двух векторов
std::ostream& operator<<(std::ostream& os, const Vector& vect); // вывод вектора в поток
double dot(const Vector& left, const Vector& right); // скалярное произведение векторов

class Vector {
private:
	std::vector<double> data;
public:
	// Конструктор вектора фиксированного размера + конструктор по-умолчанию:
	Vector(int size = 0, double value = 0)
		: data(size, value)
	{ };

	// Конструктор копирования "std::vector":
	Vector(const std::vector<double>& data = {}) 
		: data(data)
	{ };
	
	// Конструктор перемещения объекта класса "std::vector":
	Vector(std::vector<double>&& data)
		: data(std::move(data))
	{ };
	
	// Конструктор копирования объекта класса "Вектор":
	Vector(const Vector& other) 
		: data(other.data)
	{ };

	// Конструктор перемещения объекта класса "Вектор":
	Vector(Vector&& other) 
		: data(std::move(other.data)) 
	{ };
	
	// Операторы присваивания и перемещения объектов класса "Вектор": 
	Vector& operator=(const Vector& other);
	Vector& operator=(Vector&& other);
	
	// Оператор доступа к компоненте вектора:
	double& operator[](const int index);
	const double& operator[](const int index) const;
	
	// Получение размера вектора и проверка на пустоту:
	bool empty() const;
	int size() const;
	
	// Получение итераторов на начало и конец вектора:
	std::vector<double>::iterator begin();
	std::vector<double>::iterator end();
	const std::vector<double>::iterator begin() const;
	const std::vector<double>::iterator end() const;
	
	// Основные арифметические операции с изменением вектора *this:
	Vector& operator+=(const Vector& other);
	Vector& operator-=(const Vector& other);
	Vector& operator*=(const double value);
	Vector& operator/=(const double value);
		
	// Унарный минус - получение противоположного вектора:
	Vector operator-() const;
	
	// Дружественные операторы и функции:
	friend Vector operator*(const double value, Vector vect); // умножение числа на вектор
	friend Vector operator*(Vector vect, const double value); // умножение вектора на число
	friend Vector operator/(Vector vect, const double value); // деление вектора на число
	friend Vector operator+(Vector left, const Vector& right); // сложение векторов
	friend Vector operator-(Vector left, const Vector& right); // разность векторов
	friend bool operator==(const Vector& left, const Vector& right); // сравнение двух векторов
	friend std::ostream& operator<<(std::ostream& os, const Vector& vect); // вывод вектора в поток
	friend double dot(const Vector& left, const Vector& right); // скалярное произведение векторов
};
