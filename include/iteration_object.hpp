// Файл iteration_object.hpp
// Авторы: Козырев Дмитрий
#pragma once

#include "math.hpp"

// Базовый объект контроля за итерациями
class BasicIterationObject {
protected:
    Vector x_prev, x_curr;    // Предыдущая и текущая точки x
    Real f_prev, f_curr;      // Предыдущее и текущее значение функции в этих точках
    int iter_counter;         // Текущее количество итераций
    std::string method_title; // Название метода, в который объект передан
    
public:    
    BasicIterationObject();
    
    virtual ~BasicIterationObject();
    
    virtual BasicIterationObject* new_object() const;
    
    virtual bool is_stopped() const; // Условие останова метода
    
    virtual void next_iteration(const Vector& x_next, Real f_next); // Переход к следующей итерации
    
    virtual int get_iter_counter() const; // Получение текущего количества итераций
    
    virtual Real get_f_curr() const; // Получение текущего значения функции f
    virtual Real get_f_prev() const; // Получение предыдущего значения функции f
    
    virtual const Vector& get_x_curr() const; // Получение текущего значения точки x
    virtual const Vector& get_x_prev() const; // Получение предыдущего значения точки x
    
    virtual void set_f_curr(Real); // Изменение текущего значения функции f
    virtual void set_f_prev(Real); // Изменение предыдущего значения функции f
    
    virtual void set_x_curr(const Vector&); // Изменение текущего значения точки x
    virtual void set_x_prev(const Vector&); // Изменение предыдущего значения точки x
    
    virtual void set_iter_counter(int); // Установка счетчика итераций
    
    virtual std::string get_method_title() const; // Получение названия метода
    virtual void set_method_title(const std::string&); // Изменение названия метода
};

// Продвинутый объект контроля за итерациями: 
// в нем сохраняются все промежуточные посещенные точки и значения функции в них
class AdvancedIterationObject : public BasicIterationObject {
protected:
    std::vector<Vector> x_sequence; // Последовательность всех посещенных точек x
    std::vector<Real> f_sequence;   // Последовательность значений функции в этих точках

public:    
    AdvancedIterationObject();
    
    virtual BasicIterationObject* new_object() const override;
    
    virtual void next_iteration(const Vector& x_next, Real f_next) override; // Переход к следующей итерации
    
    virtual const std::vector<Vector>& get_x_seq() const; // Получение последовательности приближений x
    virtual const std::vector<Real>& get_f_seq() const; // Получение последовательности значений функции f
};