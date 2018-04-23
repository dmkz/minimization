#include "hessian_free.hpp"
#include "nesterov.hpp"
#include "bfgs.hpp"
#include "dfp.hpp"
#include "powell.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <set>

/*
    Тестирование методов отдельно.
    Автор: Юрий Кондратов, Бураханова А., Казарян М.
*/

typedef IterationData(*Method)(Function, Vector, const StopCondition&);

// Тестовые функции
Real f1(const Vector& x) {
    return 1 + x[0] + x[1] - x[0] * x[1] + x[0] * x[0] + x[1] * x[1];
}

Real f2(const Vector& x) {
    return 1 + 7 * x[0] + 5 * x[1] + 0.5 * x[0] * x[1] + 3 * x[0] * x[0] + x[1] * x[1];
}

Real f3(const Vector& x) {
    return 100 + 7 * x[0] + 5 * x[1] - 10 * x[0] * x[1] + 3 * x[0] * x[0] + 10 * x[1] * x[1];
}

Real f4(const Vector& x) {
    return 100 + 7 * x[0] + 5 * x[1] - 10.95 * x[0] * x[1] + 3 * x[0] * x[0] + 10 * x[1] * x[1];
}

Real f5(const Vector& x) {
    return 1+x[0]+x[1]+x[2] + x[0] * x[1] + x[0] * x[2] + x[1] * x[2] + x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

Real f6(const Vector& x) {
    return 10 * std::pow(x[0], 4) + 15 * std::pow(x[1], 4) + 15 * x[0] * x[1];
}

Real f7(const Vector& v) {
    Real x = v[0];
    Real y = v[1];
    return 10 * std::pow(x, 6) + 15 * std::pow(y, 6) - 20 * (std::pow(x, 3) * y + x * std::pow(y, 3));
}

Real f8(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 2 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + x*x + y*y;
}

Real f9(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 3 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + x*x + y*y;
}

Real f10(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 2 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + std::pow(x, 4) + std::pow(y, 4) - x*x - y*y;
}

Real f11(const Vector &v) {
    Real x = v[0], y = v[1];
    return (x*x + y - 11)*(x*x + y - 11) + (y*y + x - 7)*(y*y + x - 7);
}

Real f12(const Vector &v) {
    // (x1+x2+x3+x4-8)^2+(x1^2+x2^2+x3^2+x4^2-18)^2+(x1^3+x2^3+x3^3+x4^3-44)^2+(x1^4+x2^4+x3^4+x4^4-114)^2
    assert(v.size() == 4u);
    
    Vector b = {8, 18, 44, 114};    
    Real fun = 0;
    for (int i = 0; i < 4; ++i) {
        Real temp = -b[i];
        for (auto x : v) {
            temp += std::pow(x, i+1);
        }
        fun += temp * temp;
    }
    return fun;
}

Real f13(const Vector &v) {
    assert(v.size()==2u);
    Real x = v[0], y = v[1];
    return std::pow(y-5.1/(4*M_PI*M_PI)*x*x+5/M_PI*x-6, 2)+10*(1-1/(8*M_PI))*cos(x)+10;
}

Real f14(const Vector &v) {
    return std::sin(v[0] + v[1]) + std::pow(v[0] - v[1], 2) + 1.5 * v[0]*v[0] + 2.5 * v[1]*v[1] + 1;
}

Real f15(const Vector &v) {
    return 0.26 * (std::pow(v[0], 2) + std::pow(v[1], 2)) - 0.48 * v[0] * v[1];
}

Real f16(const Vector &v) {
    return norm(v) < 1e-16 ? 0 : -(std::pow(std::sin(v[0] - v[1]), 2) * std::pow(std::sin(v[0] + v[1]), 2))/norm(v);
}

Real f17(const Vector &v) {
    Real sum_sqr = 0, sum_sin = 0;
    for (auto x : v) {
        sum_sqr += x * x;
        sum_sin += sin(x)*sin(x);
    }
    return sum_sin-exp(-sum_sqr);
}

Real f18(const Vector &v) {
    return v[0] / 4 + std::pow(v[0] * v[0] - 2 * v[0] + v[1] * v[1], 2);
}

Real f19(const Vector &v) {
    Real x = v[0], y = v[1];
    return (1.5-x*(1-y))*(1.5-x*(1-y))+(2.25-x*(1-y*y))*(2.25-x*(1-y*y))+(2.625-x*(1-y*y*y))*(2.625-x*(1-y*y*y));
}

Real f20(const Vector &v) {
    Real x1 = v[0], x2 = v[1], x3 = v[2], x4 = v[3];
    return 100*(x2-x1*x1)*(x2-x1*x1)+(1-x1)*(1-x1)+90*(x4-x3*x3)*(x4-x3*x3)+(1-x3)*(1-x3)+10.1*(x2-1)*(x2-1)+10.1*(x4-1)*(x4-1);
}

Real f21(const Vector &v) {
    return std::pow(v[0] * v[0] - v[1], 2) + std::pow(1 - v[0], 2);
}

Real f22(const Vector &v) {
    return std::pow(v[0] * v[0] - v[1], 2) + 100 * std::pow(1 - v[0], 2);
}
Real f23(const Vector &v) {
    return 100 * std::pow(v[1] - std::pow(v[0], 3), 2) + std::pow(1 - v[0], 2);
}
Real f24(const Vector &v) {
    Real x1 = v[0], x2 = v[1], x3 = v[2], x4 = v[3];
    return std::pow((x1 + 10 * x2),2) + 5 * std::pow((x4-x3),2)+ std::pow((x2 - 2 * x3),4) + std::pow(10 * (x1 - x4),4);
}
Real f25(const Vector &v) {
    Real x1 = v[0], x2 = v[1], x3 = v[2], x4 = v[3];
    return std::pow(x1*x1 - x2+1, 4) + 100 * std::pow(x2-x3, 6) +
           std::pow(std::tan(x3-x4), 4) + std::pow(x1, 8) + std::pow(x4 - 1, 2);
}

Real f26(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=0; i < v.size(); ++i) {
        fun += std::pow(v[i], 2);
    }
    return fun;
}

Real f27(const Vector& v) {
    Real fun = std::pow(v[0] - 1, 2) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::pow(v[i] - 2 * std::pow(v[i-1], 2) + 1, 2);
    }
    return fun;
}

Real f28(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=1; i < v.size(); i += 2) {
        fun += 100 * std::pow(std::pow(v[i-1], 2) - v[i], 2) + std::pow(v[i-1] - 1, 2);
    }
    return fun;
}

Real f29(const Vector &v) {
    Real fun = 0;
    for (int i = 0; i < (int)v.size() - 1; i++) {
        fun += v[i] * v[i] + 2 * v[i+1] * v[i+1] - 0.3 * std::cos(3 * M_PI * v[i]) -
               0.4 * std::cos(4 * M_PI * v[i+1]) + 0.7;
    }
    return fun;
}

Real f30(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=1; i < v.size(); i++) {
        fun += v[i] * v[i];
    }
    fun *= 1000000;
    fun += v[0] * v[0];
    return fun;
}

Real f31(const Vector &x){
	return (1+std::pow((x[0]+x[1]+1),2)*(19-14*x[0]+3*x[0]*x[0] - 14*x[1] + 6*x[0]*x[1]+3*x[1]*x[1]))*
	(30+std::pow((2*x[0]-3*x[1]),2)*(18-32*x[0]+12*x[0]*x[0]+48*x[1]-36*x[0]*x[1]+27*x[1]*x[1]));
}

Real f32(const Vector &x){
	return (100*sqrt(abs(x[1]-0.01*x[0]*x[0]))+0.01*abs(x[0]+10));
}

Real f33(const Vector &v)
{
	Real fun = 0, fun1 = 0, fun2 = 0;
	
	for (unsigned int i = 0; i < v.size(); i++)
	{
		fun1 = fun1 + std::cos(5* M_PI*v[i]);
		fun2 = fun2 + v[i]*v[i];
	}
	fun = (-0.1)*fun1 + fun2;
	return fun;
}

Real f34(const Vector &x){
	return (100*std::pow(x[1]-0.01*x[0]*x[0]+1, 2)+0.01*(x[0]+10)*(x[0]+10));
}

Real f35(const Vector& v){
    assert(v.size() == 2u);
    return abs(v[0]*v[0]+v[1]*v[1]+v[0]*v[1])+abs(sin(v[0]))+abs(sin(v[1]));
}

Real f36(const Vector &x){
	return std::pow((2*std::pow(x[0],3)*x[1]-std::pow(x[1],3)),2)+std::pow((6*x[0]-x[1]*x[1]+x[1]),2);
}

// ThreeHumpCamel
Real f37(const Vector &v) {
    return (2*std::pow(v[0], 2) - 1.05 * std::pow(v[0], 4) + std::pow(v[0], 6)/6 + v[0]*v[1] + std::pow(v[1],2));
}

// SixHumpCamel
Real f38(const Vector &v) {
    return (4 - 2.1*std::pow(v[0], 2)+std::pow(v[0], 4)/3)*std::pow(v[0], 2) + v[0]*v[1] + (-4 + 4*std::pow(v[1], 2))*std::pow(v[1], 2);
}

// Branin01
Real f39(const Vector &v) {
    return std::pow((-1.275*std::pow(v[0],2)/std::pow(M_PI,2)+5*v[0]/M_PI+v[1]-6),2) + (10-5/(4*M_PI))*std::cos(v[0])+10;
}

// Branin02
Real f40(const Vector &v) {
    return std::pow((-1.275*std::pow(v[0],2)/std::pow(M_PI,2)+5*v[0]/M_PI+v[1]-6),2) + (10-5/(4*M_PI))*std::cos(v[0])*std::cos(v[1])+std::log(std::pow(v[0],2) + std::pow(v[1],2) + 1)+10;
}

// RotatedEllipse01 7x*x-6*sqrt(3)*x*y+13y*y
Real f41(const Vector &v) {
    Real x = v[0], y = v[1];
    return 7*x*x-6*std::sqrt(3)*x*y+13*y*y;
}

// EggCrate
Real f42(const Vector &v) {
    Real x = v[0], y = v[1];
    Real sin_x = sin(x), sin_y = sin(y);
    return x*x+y*y+25*(sin_x*sin_x+sin_y*sin_y);
}

// RotatedEllipse02 x*x-x*y+y*y
Real f43(const Vector& v) {
    //-x1*x2(72-2x1-2x2)
    assert(v.size()==2u);
    Real x = v[0], y = v[1];
    return x*x-x*y+y*y;
}

// Bird (x-y)^2+exp((1-sin(x))^2)cos(y)+exp((1-cos(y))^2)*sin(x)
Real f44(const Vector& v) {
    assert(v.size() == 2u);
    Real x = v[0], y = v[1];
    return (x-y)*(x-y)+exp(pow(1-sin(x),2))*cos(y)+exp(pow(1-cos(y),2))*sin(x);
}

// Hosaki (1-8x+7x^2-7.0/3*x^3+1.0/4*x^4)*y^2*exp(-y)
Real f45(const Vector& v) {
    assert(v.size() == 2u);
    Real x = v[0], y = v[1];
    return (1-8*x+7*x*x-7.0/3*x*x*x+1.0/4*x*x*x*x)*y*y*std::exp(-y*y);
}

// El-Attar-Vidyasagar-Dutta
Real f46(const Vector &x){
	return std::pow((x[0]*x[0]+x[1]-10),2)+std::pow((x[0]+x[1]*x[1]-7),2)+std::pow((x[0]*x[0]+std::pow(x[1],3)-1),2);
}

//Ursem01 
Real f47(const Vector &x)
{
	return -sin(2*x[0]-0.5*M_PI)-3*cos(x[1])+0.5*x[0]*x[0];
}

//Alpine 1
Real f48(const Vector &v)
{
	Real fun = 0;
	
	for (unsigned int i = 0; i < v.size(); i++)
	{
		fun = fun + std::abs(v[i]*sin(v[i])+0.1*v[i]);
	}
	return fun;
}

// Levy13
Real f49(const Vector &v) {
    return std::pow((v[0]-1),2)*(std::pow(std::sin(3*M_PI*v[1]),2) + 1) + std::pow((v[1]-1),2)*(std::pow(std::sin(2*M_PI*v[1]),2) + 1) + std::pow(std::sin(3*M_PI*v[0]),2);
}

// Mishra08
Real f50(const Vector &v) {
    return 0.001*std::pow((std::abs(std::pow(v[0],10) - 20*std::pow(v[0],9) + 180*std::pow(v[0],8) - 960*std::pow(v[0],7) + 3360*std::pow(v[0],6) -8064*std::pow(v[0],5) + 13340*std::pow(v[0],4) - 15360*std::pow(v[0],3) + 11520*std::pow(v[0],2) - 5120*v[0] + 2624)*std::abs(std::pow(v[1],4) + 12*std::pow(v[1],3) + 54*std::pow(v[1],2) + 108*v[1] + 81)),2);
}

struct ControlPoint {
    Vector x;
    std::string type;

    ControlPoint(Vector x, std::string type)
        : x(x)
        , type(type)
    { }
};

Matrix gen_start_points(int Dimensions, Real left, Real right) {
    Vector point(Dimensions, 1);
    point[0] = -1;
    Matrix start_points{point};
    start_points.push_back(Vector(Dimensions, left));
    start_points.push_back(Vector(Dimensions, right));
    point = Vector(Dimensions, right);
    for (int i = Dimensions / 2; i < Dimensions; ++i) {
        point[i] = left;
    }
    point = Vector(Dimensions, left);
    for (int i = Dimensions / 2; i < Dimensions; ++i) {
        point[i] = right;
    }
    start_points.push_back(point);
    return start_points;
}

struct Test {
    std::string id;                         // Идентификатор теста (например: "1", "2", "30_2")
    
    Function f;                             // Целевая функция
    std::string description_f;              // Ее символьное описание
    
    StopCondition stop_condition;           // Условие останова
    std::string description_stop_condition; // Его символьное описание
    
    std::vector<ControlPoint> expected;     // Ожидаемые точки
    
    std::vector<Vector> start_points;       // Стартовые точки
    std::vector<std::vector<IterationData>> result;    // Результаты тестирования из стартовых точек
};

std::vector<Test> Tests; // Массив всех тестов

bool example_stop_condition(const IterationData& iter_data) {
    return iter_data.iter_counter >= 100 || std::abs(iter_data.f_curr - iter_data.f_prev) < 1e-8;
}

const std::string descript_ex_stop_cond = "iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001";

void prepare_tests() {
    // Добавление теста 01:
    Tests.push_back(Test{
        "Test 01, dim 02", f1, "Гладкая функция:\n\tf(x,y) = 1+x+y-xy+x^2+y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-1, -1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5)/* стартовые точки */, {}
    });
    // Добавление теста 02:
    Tests.push_back(Test{
        "Test 02, dim 02", f2, "Гладкая функция:\n\tf(x,y) = 1+7x+5y+0.5xy+3x^2+y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{Real(-46) / 47, Real(-106) / 47}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 03:
    Tests.push_back(Test{
        "Test 03, dim 02", f3, "Гладкая функция:\n\tf(x,y) = 100+7x+5y-10xy+3x^2+10y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{Real(-19) / 2, -5}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 04:
    Tests.push_back(Test{
        "Test 04, dim 02", f4, "Гладкая функция:\n\tf(x,y) = 100+7x+5y-10.95xy+3x^2+10y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{Real(-77900) / 39, Real(-14220)/13}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 05:
    Tests.push_back(Test{
        "Test 05, dim 03", f5, "Гладкая функция:\n\tf(x,y,z) = 1+x+y+z+xy+xz+yz+x^2+y^2+z^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-0.25, -0.25, -0.25}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(3, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 06:
    Tests.push_back(Test{
        "Test 06, dim 02", f6, "Гладкая функция:\n\tf(x,y) = 10x^4+15y^4+15xy", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-0.582109, 0.525995}, "Global Min"}, {{0.582109, -0.525995}, "Global Min"}, {{0, 0}, "Saddle point"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 07:
    Tests.push_back(Test{
        "Test 07, dim 02", f7, "Гладкая функция:\n\tf(x,y) = 10x^6+15y^6-20(x^3y+xy^3)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-1.08789, -1.00318}, "Global Min"}, {{1.08789, 1.00318}, "Global Min"}, {{0, 0}, "Saddle point"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 08:
    Tests.push_back(Test{
        "Test 08, dim 02", f8, "Гладкая функция:\n\tf(x,y) = x^6+y^6-2(x^3y+xy^3)+x^2+y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-1, -1}, "Global Min"}, {{0, 0}, "Global Min"}, {{1, 1}, "Global Min"},
            {{-0.57735, -0.57735}, "Saddle point"}, {{0.57735, 0.57735}, "Saddle point"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 09:
    Tests.push_back(Test{
        "Test 09, dim 02", f9, "Гладкая функция:\n\tf(x,y) = x^6+y^6-3(x^3y+xy^3)+x^2+y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-1.34777, -1.34777}, "Global Min"}, {{1.34777, 1.34777}, "Global Min"}, {{0, 0}, " Local Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 10:
    Tests.push_back(Test{
        "Test 10, dim 02", f10, "Гладкая функция:\n\tf(x,y) = x^6+y^6-2(x^3y+xy^3)+x^4+y^4-x^2-y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
            {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 11:
    Tests.push_back(Test{
        "Test 11, dim 02", f11, "Гладкая функция Химмельблау:\n\tf(x,y) = (x^2+y-11)^2+(x+y^2-7)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{3, 2}, "Global Min"}, {{-3.77931, -3.28319}, "Global Min"},
            {{-2.80512, 3.13131}, "Global Min"}, {{3.58443, -1.84813}, "Global Min"},
            {{-0.270845, -0.923039}, " Local Max"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 12:
    Tests.push_back(Test{
        "Test 12, dim 04", f12, "Гладкая функция:\n\tf(x) = (x1+x2+x3+x4-8)^2+(x1^2+x2^2+x3^2+x4^2-18)^2+(x1^3+x2^3+x3^3+x4^3-44)^2+(x1^4+x2^4+x3^4+x4^4-114)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1,2,2,3}, "Global Min"}}, // Ожидаемые точки
        {{2,2,2,2},{-2,2,-2,2},{2,-2,2,-2},{-2,-2,-2,-2}} /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 13:
    Tests.push_back(Test{
        "Test 13, dim 02", f13, "Гладкая функция Бранина:\n\tf(x,y) = (y-5.1*x^2/(4*pi^2)+5*x/pi-6)^2+10(1-1/(8pi)*cos(x))+10", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-M_PI,12.275},  "Global Min"},
            {{M_PI,2.275},    "Global Min"},
            {{9.42478,2.475}, "Global Min"}}, // Ожидаемые точки
        {{-10,-10},{-10,10},{10,-10},{10,10}} /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 14:
    Tests.push_back(Test{
        "Test 14, dim 02", f14, "Модифицированная гладкая функция МакКормика:\n\tf(x,y) = sin(x+y)+(x-y)^2+1.5x^2+2.5y^2+1", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-0.25988392,-0.20213194},  "Global Min"}},
        {{-10,-10},{-10,10},{10,-10},{10,10}} /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 15:
    Tests.push_back(Test{
        "Test 15, dim 02", f15, "Гладкая функция Матиаса:\n\tf(x,y) = 0.26(x^2+y^2)-0.48xy", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -10, 10) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 16:
    Tests.push_back(Test{
        "Test 16, dim 02", f16, "Гладкая функция Кин:\n\tf(x,y) = -sin(x-y)^2*sin(x+y)^2/sqrt(x^2+y^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{ 0.00000,  1.39325}, "Global Min"}, 
            {{ 1.39325,  0.00000}, "Global Min"},
            {{ 0.00000, -1.39325}, "Global Min"},
            {{-1.39325,  0.00000}, "Global Min"}}, // Ожидаемые точки
        {{0,-2},{0, 2},{-2,0},{2,0}} /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_2:
    Tests.push_back(Test{
        "Test 17, dim 02", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}},  // Ожидаемые точки
        {{1,1},{1.2,1.2},{1.4,1.4},{1.6,1.6},{1.8,1.8},{2,2}} /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_4:
    Tests.push_back(Test{
        "Test 17, dim 04", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0, 0, 0}, "Global Min"}},  // Ожидаемые точки
        {
            {1.0,1.0,1.0,1.0},
            {1.2,1.2,1.2,1.2},
            {1.4,1.4,1.4,1.4},
            {1.6,1.6,1.6,1.6},
            {1.8,1.8,1.8,1.8},
            {2.0,2.0,2.0,2.0}
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_8:
    Tests.push_back(Test{
        "Test 17, dim 08", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(8,0), "Global Min"}},  // Ожидаемые точки
        {
            Vector(8, 1.0),
            Vector(8, 1.2),
            Vector(8, 1.4),
            Vector(8, 1.6),
            Vector(8, 1.8),
            Vector(8, 2.0)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_16:
    Tests.push_back(Test{
        "Test 17, dim 16", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(16,0), "Global Min"}},  // Ожидаемые точки
        {
            Vector(16, 1.0),
            Vector(16, 1.2),
            Vector(16, 1.4),
            Vector(16, 1.6),
            Vector(16, 1.8),
            Vector(16, 2.0)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_32:
    Tests.push_back(Test{
        "Test 17, dim 32", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(32,0), "Global Min"}},  // Ожидаемые точки
        {
            Vector(32, 1.0),
            Vector(32, 1.2),
            Vector(32, 1.4),
            Vector(32, 1.6),
            Vector(32, 1.8),
            Vector(32, 2.0)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 17_64:
    Tests.push_back(Test{
        "Test 17, dim 64", f17, "Модифицированная гладкая функция XinSheYang04:\n\tf(x) = sum(sin(x(i))^2)-exp(-sum(x(i)^2))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(64,0), "Global Min"}},  // Ожидаемые точки
        {
            Vector(64, 1.0),
            Vector(64, 1.2),
            Vector(64, 1.4),
            Vector(64, 1.6),
            Vector(64, 1.8),
            Vector(64, 2.0)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 18:
    Tests.push_back(Test{
        "Test 18, dim 02", f18, "Гладкая функция Зеттла:\n\tf(x,y) = x/4+(x^2-2x+y^2)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-0.029896, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -1, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 19:
    Tests.push_back(Test{
        "Test 19, dim 02", f19, "Гладкая функция Биля:\n\tf(x,y) = (xy-x+1.5)^2+(xy^2-x+2.25)^2+(x*y^3-x+2.625)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{3, 0.5}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -0.5, 0.5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
     
    // Добавление теста 20:
    Tests.push_back(Test{
        "Test 20, dim 04", f20, "Модифицированная гладкая функция Колвилля:\n\tf(x) = 100*(x2-x1*x1)^2+(1-x1)^2+90*(x4-x3*x3)^2+(1-x3)^2+10.1*(x2-1)^2+10.1*(x4-1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(4,1), "Global Min"}}, // Ожидаемые точки
        {{-3,-1,-3,-1},{-0.5, 0.5, -0.5, 0.5}, Vector(4, 30), Vector(4, 100)} /* стартовые точки */, {} /* пустой вектор результатов */
    }); 
    
    // Добавление теста 21:
    Tests.push_back(Test{
        "Test 21, dim 02", f21, "Гладкая функция:\n\tf(x,y) = (y-x^2)^2+(1-x)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 22:
    Tests.push_back(Test{
        "Test 22, dim 02", f22, "Гладкая функция:\n\tf(x,y) = (y-x^2)^2+100(1-x)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 23:
    Tests.push_back(Test{
        "Test 23, dim 02", f23, "Гладкая функция:\n\tf(x,y) = 100(y-x^3)^2+(1-x)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 24:
    Tests.push_back(Test{
        "Test 24, dim 04", f24, "Гладкая функция:\n\tf(x,y,z,t) = (x+10y)^2+5(z-t)^2+(y-2z)^4+10(x-t)^4", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0, 0, 0}, "Global Min"}}, // Ожидаемые точки
        {
            {-3,-1, 0, 1},
            {-1, 1, 1, 1},
            {-5,-5,-5,-5},
            { 5, 5, 5, 5},
            {-5,-5, 5, 5}
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    
	// Добавление теста 25:
    Tests.push_back(Test{
      "Test 25, dim 04", f25, "Гладкая функция:\n\tf(x,y,z,t) = (x^2-y+1)^4+100(y-z)^6+tg^4(z-t)+x^8+(t-1)^2", // Номер теста, функция, ее описание
       example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        {
			{ 1, 2, 2, 2},
            {-1, 1, 1, 1},
            {-2,-2,-2,-2},
            { 2, 1, 0,-1},
            {-2,-2, 2, 2}
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 26_2:
    Tests.push_back(Test{
        "Test 26, dim 02", f26, "Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 26_4:
    Tests.push_back(Test{
        "Test 26, dim 04", f26, "Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0, 0, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(4, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 26_8:
    Tests.push_back(Test{
        "Test 26, dim 08", f26, "Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(8, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 26_12:
    Tests.push_back(Test{
        "Test 26, dim 12", f26, "Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(12, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 27_2:
    Tests.push_back(Test{
        "Test 27, dim 02", f27, "Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 27_4:
    Tests.push_back(Test{
        "Test 27, dim 04", f27, "Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(4, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 27_8:
    Tests.push_back(Test{
        "Test 27, dim 08", f27, "Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(8, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 27_12:
    Tests.push_back(Test{
        "Test 27, dim 12", f27, "Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(12, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 28_2:
    Tests.push_back(Test{
        "Test 28, dim 02", f28, "Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 28_4:
    Tests.push_back(Test{
        "Test 28, dim 04", f28, "Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(4, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 28_8:
    Tests.push_back(Test{
        "Test 28, dim 08", f28, "Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(8, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 28_12:
    Tests.push_back(Test{
        "Test 28, dim 12", f28, "Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(12, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	
    // Добавление теста 29_2:
    Tests.push_back(Test{
        "Test 29, dim 02", f29, "Гладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(2, 0), "Global Min"}}, // Ожидаемые точки
        {
            Vector(2,1.0),
            Vector(2,1.5),
            Vector(2,2.0),
            Vector(2,2.5)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 29_4:
    Tests.push_back(Test{
        "Test 29, dim 04", f29, "Гладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(4, 0), "Global Min"}}, // Ожидаемые точки
        {
            Vector(4,1.0),
            Vector(4,1.5),
            Vector(4,2.0),
            Vector(4,2.5)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 29_08:
    Tests.push_back(Test{
        "Test 29, dim 08", f29, "Гладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(8, 0), "Global Min"}}, // Ожидаемые точки
        {
            Vector(8,1.0),
            Vector(8,1.5),
            Vector(8,2.0),
            Vector(8,2.5)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    // Добавление теста 29_16:
    Tests.push_back(Test{
        "Test 29, dim 16", f29, "Гладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{Vector(16, 0), "Global Min"}}, // Ожидаемые точки
        {
            Vector(16,1.0),
            Vector(16,1.5),
            Vector(16,2.0),
            Vector(16,2.5)
        } /* стартовые точки */, {} /* пустой вектор результатов */
    });
    
    
	// Добавление теста 30_2:
    Tests.push_back(Test{
        "Test 30, dim 02", f30, "Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0,0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(2, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 30_4:
    Tests.push_back(Test{
        "Test 30, dim 04", f30, "Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0,0,0,0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(4, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 30_8:
    Tests.push_back(Test{
        "Test 30, dim 08", f30, "Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0,0,0,0,0,0,0,0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(8, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 30_12:
    Tests.push_back(Test{
        "Test 30, dim 12", f30, "Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{0,0,0,0,0,0,0,0,0,0,0,0}, "Global Min"}}, // Ожидаемые точки
        gen_start_points(12, -5, 5) /* стартовые точки */, {} /* пустой вектор результатов */
    });
	// Добавление теста 31:
    Tests.push_back(Test{
        "Test 31, dim 02", f31, "Гладкая функция Голдштейна-Прайса:\n\tf(x,y) = [1+(x+y+1)^2(19-14x+3x^2-14y+6xy+3y^2)][30+(2x-3y)^2(18-32x+12x^2+48y-36xy+27y^2)]", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
		{{     0,     -1}, "Global Min"},
        {{ 4.0/5,  1.0/5}, "Local Max"},
        {{-3.0/5, -2.0/5}, "Local Min"},
        {{ 6.0/5,  4.0/5}, "Local Min"},
        {{ 9.0/5,  1.0/5}, "Local Min"}
        }, // Ожидаемые точки
        Matrix{
            {1, 1},
            {0.25, -0.125},
            {0.5, -1.5},
            {0, -1.2}
        }, {} // стартовые точки
    });
	// Добавление теста 32:
    Tests.push_back(Test{
        "Test 32, dim 02", f32, "Негладкая функция Bukin06:\n\tf(x,y) = 100*sqrt(|y-0.01*x^2|)+0.01|x+10|", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{{{-10, 1}, "Global Min"}}, // Ожидаемые точки
        Matrix{
            {-15, -3},
            {-15,  3},
            { -5, -3},
            { -5,  3}
        } /* стартовые точки */, {}
    });
    
    // Добавление теста 33 размерности 2:
    Tests.push_back(Test{
        "Test 33, dim 02", f33, "Гладкая функция Cosine Mixture:\n\tf(x1, ..., xn) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i)))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(2, -1, 1)/* стартовые точки */, {}
    });
    // Добавление теста 33 размерности 4:
    Tests.push_back(Test{
        "Test 33, dim 04", f33, "Гладкая функция Cosine Mixture:\n\tf(x1, ..., xn) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i)))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0, 0, 0}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(4, -1, 1)/* стартовые точки */, {}
    });
    // Добавление теста 33 размерности 8:
    Tests.push_back(Test{
        "Test 33, dim 08", f33, "Гладкая функция Cosine Mixture:\n\tf(x1, ..., xn) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i)))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(8, -1, 1)/* стартовые точки */, {}
    });
    // Добавление теста 33 размерности 12:
    Tests.push_back(Test{
        "Test 33, dim 12", f33, "Гладкая функция Cosine Mixture:\n\tf(x1, ..., xn) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i)))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(12, -1, 1)/* стартовые точки */, {}
    });
    // Добавление теста 34:
    Tests.push_back(Test{
        "Test 34, dim 02", f34, "Гладкая функция Bukin02:\n\tf(x,y) = 100*(y-0.01*x^2+1)^2+0.01*(x+10)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-10, 0}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-15, -3},
            {-15, 3},
            {-5, -3},
            {-5, 3},
            {-10.5, 1.5}
        }/* стартовые точки */, {}
    });
    // Добавление теста 35:
    Tests.push_back(Test{
        "Test 35, dim 02", f35, "Негладкая функция BartelsConn:\n\tf(x,y) = |x^2+y^2+x*y|+|sin(x)|+|sin(y)|", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(2, -5, 5)/* стартовые точки */, {}
    });
    
    // Добавление теста 36:
    Tests.push_back(Test{
        "Test 36, dim 02", f36, "Гладкая функция Price04:\n\tf(x,y) = (2*x^3*y-y^3)^2+(6x-y^2+y)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"},
            {{2, 4}, "Global Min"},
            {{1.464, -2.506}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {3, 5.5},
            {1.5, -3},
            {0, -5},
            {-1, 1.5}
        }/* стартовые точки */, {}
    });
    // Добавление теста 37:
    Tests.push_back(Test{
        "Test 37, dim 02", f37, "Гладкая функция:\n\tf(x,y) = (2x^2 -1.05x^4 + x^6/6 + xy) + y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"},
            {{-1.74755, 0.873776}, " Local Min"},
            {{1.74755, -0.873776}, " Local Min"}
        }, // Ожидаемые точки
        gen_start_points(2, -1, 1)/* стартовые точки */, {}
    });
    // Добавление теста 38:
    Tests.push_back(Test{
        "Test 38, dim 02", f38, "Гладкая функция:\n\tf(x,y) = (4-2.1*x1^2+x1^4/3)*x1^2+x1*x2+(-4+4*x2^2)*x2^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-0.0898, 0.7126}, "Global Min"},
            {{0.0898, -0.7126}, "Global Min"},
            {{-1.70361,  0.796084}, " Local Min"},
            {{ -1.6071, -0.568651}, " Local Min"},
            {{  1.6071,  0.568651}, " Local Min"},
            {{ 1.70361, -0.796084}, " Local Min"},
            {{-1.23023, -0.162335}, " Local Max"},
            {{ 1.23023,  0.162335}, " Local Max"}
        }, // Ожидаемые точки
        Matrix{
            {0,0},
            {1,1},
            {-1, -1},
            {-2, 0},
            {0, 2}
        }/* стартовые точки */, {}
    });
    // Добавление теста 39:
    Tests.push_back(Test{
        "Test 39, dim 02", f39, "Гладкая функция:\n\tf(x,y) = (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+10", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-M_PI, 12.275}, "Global Min"},
            {{M_PI, 2.275}, "Global Min"},
            {{9.42478, 2.475}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-5, 0},
            {10, 0},
            {-5, 15},
            {10, 15}
        }/* стартовые точки */, {}
    });
    // Добавление теста 40:
    Tests.push_back(Test{
        "Test 40, dim 02", f40, "Гладкая функция:\n\tf(x,y) = (-1.275*x1^2/M_PI^2+5*x1/M_PI+x2-6)^2 + (10-5/(4*M_PI))*cos(x1)*cos(x2)+log(x1^2+x2^2+1)+10", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-3.2, 12.53}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-2, -2},
            {-2, 11},
            {10, -2},
            {10, 10}
        }/* стартовые точки */, {}
    });
    // Добавление теста 41:
    Tests.push_back(Test{
        "Test 41, dim 02", f41, "Гладкая функция RotatedEllipse01:\n\tf(x,y) = 7x^2-6*sqrt(3)*x*y+13y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-15, -5},
            {-500, 500}
        }/* стартовые точки */, {}
    });
    // Добавление теста 42:
    Tests.push_back(Test{
        "Test 42, dim 02", f42, "Гладкая функция:\n\tf(x,y) = x^2 + y^2 + 25[sin^2(x) + sin^2(y)]", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-1,-1},
            {-1,1},
            {1,-1},
            {1,1}
        }/* стартовые точки */, {}
    });
    // Добавление теста 43:
    Tests.push_back(Test{
        "Test 43, dim 02", f43, "Гладкая функция RotatedEllipse02:\n\tf(x,y) = x^2-x*y+y^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0, 0}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {500, 0},
            {0, 500},
            {500, 500},
            {100, 100}
        }/* стартовые точки */, {}
    });
    // Добавление теста 44:
    Tests.push_back(Test{
        "Test 44, dim 02", f44, "Гладкая функция Bird:\n\tf(x,y) = (x-y)^2+exp((1-sin(x))^2)cos(y)+exp((1-cos(y))^2)*sin(x)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{4.701055751981055, 3.152946019601391}, "Global Min"},
            {{-1.582142172055011, -3.130246799635430}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {-M_PI, -M_PI},
            {-M_PI, M_PI},
            {M_PI, -M_PI},
            {M_PI, M_PI}
        }/* стартовые точки */, {}
    });
    // Добавление теста 45:
    Tests.push_back(Test{
        "Test 45, dim 02", f45, "Гладкая функция Hosaki:\n\tf(x,y) = (1-8x+7x^2-7.0/3*x^3+1.0/4*x^4)*y^2*exp(-y^2)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{4,-1}, "Global Min"},
            {{4, 1}, "Global Min"},
            {{1, 1}, " Local Min"},
            {{1,-1}, " Local Min"}
        }, // Ожидаемые точки
        Matrix{
            {0,0},
            {-3,-1},
            {-3,1},
            {3,-1},
            {3,1}
        }/* стартовые точки */, {}
    });
    // Добавление теста 46:
    Tests.push_back(Test{
        "Test 46, dim 02", f46, "Гладкая функция El-Attar-Vidyasagar-Dutta:\n\tf(x,y) = (x^2+y-10)^2+(x+y^2-7)^2+(x^2+y^3-1)^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{3.40919, -2.17143}, "Global Min"},
            {{-3.62315, -2.38415}, "Local Min"},
            {{-1.52071, 1.41228}, "Local Min"},
            {{2.27617, 0.864777}, "Local Min"}
        }, // Ожидаемые точки
        Matrix{
            {3.5, -2.5},
            {5, 0},
            {4.125, -0.125},
            {3.5, -1.5}
        }/* стартовые точки */, {}
    });
    // Добавление теста 47:
    Tests.push_back(Test{
        "Test 47, dim 02", f47, "Гладкая функция Ursem01:\n\tf(x,y) = -sin(2x-0.5*pi)-3cos(y)+0.5x^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{-1.23729,0}, "Global Min"},
            {{1.23729,0}, "Global Min"},
        }, // Ожидаемые точки
        Matrix{
            {-2, -2},
            {2, -2},
            {-2, 2},
            {2, 2}
        }/* стартовые точки */, {}
    });
    // Добавление теста 48:
    Tests.push_back(Test{
        "Test 48, dim 02", f48, "Негладкая функция Alpine01:\n\tf(x,y) = sum(abs(x(i)*sin(x(i))+0.1x(i)))", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{0,0}, "Global Min"}
        }, // Ожидаемые точки
        Matrix{
            {0.5, 0.5},
            {1, 1},
            {1.5, 1.5},
            {2, 2},
            {4, 4}
        }/* стартовые точки */, {}
    });
    
    // Добавление теста 49:
    Tests.push_back(Test{
        "Test 49, dim 02", f49, "Гладкая функция Levy13:\n\tf(x,y) = (x-1)^2*[sin^2(3 pi y) + 1] + (y-1)^2*[sin^2(2 pi y) + 1] + sin^2(3 pi x)", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{1,1}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(2, -10, 10) /* стартовые точки */, {}
    });
    
    // Добавление теста 50:
    Tests.push_back(Test{
        "Test 50, dim 02", f50, "Негладкая функция Mishra08:\n\tf(x,y) = 0.001[|x^10 -20x^9 +180x^8 -960x^7 +3360x^6 -8064x^5 +13340x^4 -15360x^3 +11520x^2 -5120x +2624||y^4 +12y^3 +54y^2 +108y +81|]^2", // Номер теста, функция, ее описание
        example_stop_condition, descript_ex_stop_cond,       // Условие останова и его описание
        std::vector<ControlPoint>{
            {{2, -3}, "Global Min"}
        }, // Ожидаемые точки
        gen_start_points(2, -10, 10) /* стартовые точки */, {}
    });
}

int find_point_to_converge(Function f, const Vector& x, const std::vector<ControlPoint>& expected) {
// Поиск экстремума, к которому сходится метод
// f - указатель на функцию, x - точка, полученная в результате рабоыт метода, expected - экстремумы
// возврат номера экстремума в списке экстремумов

    int id_to_converge = 0;
    Real min_diff_f = 1e9;
    Real min_diff_x = 1e9;
    for (int i = 0; i < (int)expected.size(); ++i) {
        Real curr_diff_f = std::abs(f(x) - f(expected[i].x));
        Real curr_diff_x = norm(x - expected[i].x);
        if ((std::abs(curr_diff_f-min_diff_f) <= 2 * COMPARE_EPS && curr_diff_x < min_diff_x) || curr_diff_f < min_diff_f) {
            id_to_converge = i;
            min_diff_f = curr_diff_f;
            min_diff_x = curr_diff_x;
        }
    }
    return id_to_converge;
}

std::map<std::string, std::ofstream> files; // Файлы для вывода для каждого метода отдельно!
std::map<std::string, std::string> method_to_filename; // Получение по названию метода имени файла

void run_tests() { // Запуск тестирования с сохранением результата
    for (auto & t : Tests) { // Проход по всем тестам
        
        // Вывод идентификатора текущего теста в консоль:
        std::cout << "-- Starting " << t.id << std::endl;
        
        const int nPoints = (int)t.start_points.size(); // Количество точек в тесте
        t.result.resize(nPoints);                       // Выделяем память для вектора под ответ
        
        // Проход по всем точкам - запуск всех методов из каждой точки:
        for (int i = 0; i < nPoints; ++i) {
            
            // Проход по всем методам - получение результатов от каждого метода:
            for (auto method : {bfgs, dfp, powell, hessian_free, nesterov}) {
                t.result[i].push_back(method(t.f, t.start_points[i], t.stop_condition));
            }
            
            // Сортировка в порядке возрастания значений функции в найденных точках:
            std::stable_sort(t.result[i].begin(), t.result[i].end(), [](const IterationData& left, const IterationData& right){
                if (2 * COMPARE_EPS > std::abs(left.f_curr - right.f_curr)) {
                    return left.iter_counter < right.iter_counter;
                }
                
                if (norm(left.x_curr) > 1e9) {
                    if (norm(right.x_curr) > 1e9) {
                        return norm(left.x_curr) < norm(right.x_curr);
                    } else {
                        return false;
                    }
                } else {
                    if (norm(right.x_curr) > 1e9) {
                        return true;
                    } else {
                        return left.f_curr < right.f_curr;
                    }
                }
            });
            
            // Проход по результату с целью заранее открыть файлы для вывода результатов в них
            for (const auto & r : t.result[i]) {
                // Открываем файл для вывода на будущее если файл еще был не открыт:
                if (files.find(r.method_title) == files.end()) {
                    char buf[50];
                    sprintf(buf, "test %s.txt", r.method_title.c_str());
                    for (char* it = buf; *it != '\0'; ++it) *it = (*it == ' ') ? ('-') : (*it);
                    method_to_filename[r.method_title] = buf;
                    files[r.method_title].open(buf);
                }
            }
        }
    }
}

void print_results_per_methods() { // Вывод информации о тестировании методов для каждого метода отдельно!
    for (auto & it : files) // В каждый файл выводим название метода
        it.second << it.first << "\n\n";
    
    for (auto & t : Tests) { // Проход по всем тестам
        for (auto & it : files) { // Запись первоначальной информации о тесте в каждый файл
            it.second << std::string(48, '-') << " " << t.id << " " << std::string(48, '-') << std::endl << std::endl;
            it.second << t.description_f << "\nПодробнее в документе \"Тестовые функции\"\n\n";
            it.second << "Условие останова: " << t.description_stop_condition << "\n\n";
        }
        const int nPoints = (int)t.result.size();
        // Проход по вектору-результату:
        for (int i = 0; i < nPoints; ++i) {
            for (auto & r : t.result[i]) {
                auto& fout_txt = files[r.method_title];
                
                // Находим точку, к которой сходится метод:
                int best_id = find_point_to_converge(t.f, r.x_curr, t.expected);
                auto best_point = t.expected[best_id];
                fout_txt << "\tИз точки #" << i << " сходится к " << best_point.type << " #" << best_id
                    << " (итераций: " << std::setw(8) << r.iter_counter << "), отклонение = " 
                    << std::fixed << std::setprecision(24) << std::setw(30) << r.f_curr - t.f(best_point.x) << std::endl;                
                
                fout_txt << "\t\t                Начальная точка: " << t.start_points[i] << std::endl;
                fout_txt << "\t\t  Предполагаемая точка минимума: " << best_point.x << std::endl;
                fout_txt << "\t\t      Полученная точка минимума: " << r.x_curr << std::endl;
                fout_txt << "\t\tПредполагаемое значение функции: " 
                    << std::setprecision(8) << std::fixed << std::setw(16) << t.f(best_point.x) << std::endl;
                fout_txt << "\t\t    Полученное значение функции: " 
                    << std::setprecision(8) << std::fixed << std::setw(16) << r.f_curr << std::endl << std::endl;
            }
        }
    }
    for (auto& it : files) it.second.close(); // Закрытие файлов
    
    // Вывод сопроводительной информации об успешном записи файлов в консоль:
    std::cout << "\nResults of separate testing of minimization methods are located in the files:\n";
    for (auto& it : files) {
        std::cout << "\t* " << method_to_filename[it.first] << std::endl;
    }
    std::cout << std::endl;
    files.clear(); // Очистка контейнера с файлами
}

// Подсчет рейтинга методов:
void calc_score() {
    std::map<std::string, int> score;
    
    // Подсчет рейтинга и вывод на экран
    for (auto & t : Tests)
        for (auto & r : t.result) {
            score[r[0].method_title] += 3;
            score[r[1].method_title] += 2;
            score[r[2].method_title] += 1;
        }
    
    std::set<std::pair<int, std::string>> table;
    for (const auto & it : score) {
        table.insert({it.second, it.first});
    }
    std::cout << std::setw(12) << "Place" << " | " << std::setw(12) << "Total score" << " | " << std::setw(16) << "Method title" << std::endl;
    std::cout << std::string(12,'-') << "-+-" << std::string(12,'-') << "-+-" << std::string(16,'-') << std::endl;
    int i = 1;
    for (auto it = table.rbegin(); it != table.rend(); ++it) {
        std::cout << std::setw(12) << i++ << " | " << std::setw(12) << it->first << " | " << std::setw(16) << it->second << std::endl;
    }
    std::cout << std::endl;
}

void compare_tests() {  
    std::string filename = "compare_tests.txt";
    std::ofstream fout(filename); // Открываем файл
	for (auto & t : Tests) { // Цикл по всем тестам
        // Вычислим глобальный минимум целевой функции:
        Real f_min = t.f(t.expected[0].x);
        
        // Выведем шапку теста:
        fout << std::string(80, '-') << "\n";
        fout << t.id << ": " << t.description_f << " (критерий останова: " << t.description_stop_condition << ")\n\n";
        
        // Для каждой начальной точки выведем лучший метод:
        for (int i = 0; i < (int)t.start_points.size(); ++i) { 
            auto& r = t.result[i];           // Результаты по текущему методу
            Real diff = r[0].f_curr - f_min; // Отклонение относительно глобального минимума
            fout << "\tИз точки (";
			for(auto& coord: t.start_points[i]) {
				fout << std::setprecision(2) << std::fixed << std::setw(7) << coord;
			}
			fout << ") справился лучше всего метод " << std::setw(12) << r[0].method_title 
                << " (отклонение от глобального минимума = " << std::setw(12) << diff << ", число итераций = " << std::setw(4) << r[0].iter_counter << ")" << std::endl;
        }
        fout << std::endl;
    }
    fout.close(); // закрываем файл
    std::cout << "-- Comparison of all methods for each test finished!\n";
    std::cout << "-- Results are available in the file \"" << filename << "\""<< std::endl;
}


int main() {
    prepare_tests();
    run_tests();
    print_results_per_methods();
    calc_score();
	compare_tests();
    
    return 0;
}
