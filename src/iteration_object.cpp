// Файл iteration_object.cpp
// Авторы: Козырев Дмитрий, Кондратов Юрий

#include "iteration_object.hpp"

// Базовый объект контроля за итерациями: конструктор
BasicIterationObject::BasicIterationObject() : iter_counter(0) { 
    
}

// Деструктор
BasicIterationObject::~BasicIterationObject() { 
    
}

// Создание нового объекта:
BasicIterationObject* BasicIterationObject::new_object() const {
    return new BasicIterationObject();
}

// Условие останова метода
bool BasicIterationObject::is_stopped() const { 
    return iter_counter >= 100 || std::abs(f_curr - f_prev) < 1e-8;
}

// Переход к следующей итерации
void BasicIterationObject::next_iteration(const Vector& x_next, Real f_next) {
    f_prev = f_curr;
    x_prev = x_curr;
    x_curr = x_next;
    f_curr = f_next;
    iter_counter++;
}

// Получение текущего количества итераций
int BasicIterationObject::get_iter_counter() const {
    return iter_counter;
}

// Получение текущего значения функции f
Real BasicIterationObject::get_f_curr() const {
    return f_curr;
}

// Получение предыдущего значения функции f
Real BasicIterationObject::get_f_prev() const {
    return f_prev;
}

// Получение текущего значения точки x
const Vector& BasicIterationObject::get_x_curr() const {
    return x_curr;
}

// Получение предыдущего значения точки x
const Vector& BasicIterationObject::get_x_prev() const {
    return x_prev;
}

// Изменение текущего значения функции f
void BasicIterationObject::set_f_curr(Real f) {
    f_curr = f;
}

// Изменение предыдущего значения функции f
void BasicIterationObject::set_f_prev(Real f) {
    f_prev = f;
}

// Изменение текущего значения точки x
void BasicIterationObject::set_x_curr(const Vector& x) {
    x_curr = x;
}

// Получение предыдущего значения точки x
void BasicIterationObject::set_x_prev(const Vector& x) {
    x_prev = x;
}

// Установка счетчика итераций
void BasicIterationObject::set_iter_counter(int value) {
    iter_counter = value;
}

// Получение названия метода
std::string BasicIterationObject::get_method_title() const {
    return method_title;
}

// Изменение названия метода
void BasicIterationObject::set_method_title(const std::string& s) {
    method_title = s;
}

// Продвинутый объект контроля за итерациями: 
// в нем сохраняются все промежуточные посещенные точки и значения функции в них
AdvancedIterationObject::AdvancedIterationObject() : BasicIterationObject() { }

// Создание нового объекта:
BasicIterationObject* AdvancedIterationObject::new_object() const {
    return new AdvancedIterationObject();
}

// Переход к следующей итерации
void AdvancedIterationObject::next_iteration(const Vector& x_next, Real f_next) {
    x_sequence.push_back(x_curr);
    f_sequence.push_back(f_curr);
    x_prev = x_curr;
    f_prev = f_curr;
    x_curr = x_next;
    f_curr = f_next;
    iter_counter++;
}

// Получение последовательности приближений x:
const std::vector<Vector>& AdvancedIterationObject::get_x_seq() const {
    return x_sequence;
}

// Получение последовательности значений функции f:
const std::vector<Real>& AdvancedIterationObject::get_f_seq() const {
    return f_sequence;
}