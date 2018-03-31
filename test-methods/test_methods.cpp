#include "hessian_free.hpp"
#include "nesterov.hpp"
#include "bfgs.hpp"
#include "dfp.hpp"
#include "powell.hpp"
#include <fstream>

/*
    Тестирование методов отдельно.
    Автор: Юрий Кондратов, Бураханова А., Казарян М.
*/

typedef IterationData(*Method)(Function, Vector, const StopCondition&);

std::ofstream fout_txt;

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
    const int b[] = {8, 18, 44, 144};
    Real fun = 0;
    Real g_fun;
    for (int i = 0; i < 4; ++i) {
        g_fun = 0;
        for (int j = 0; j < 4; ++j) {
            g_fun += std::pow(v[j], i);
        }
        g_fun -= b[i];
        fun += std::pow(g_fun, 2);
    }
    return fun;
}

Real f13(const Vector &v) {
    return std::pow(v[1] - 5.1 / (4 * M_PI_2) * std::pow(v[0], 2) + 5 * v[0] / M_PI - 6, 2) +
           10 * (1 - 1 / (8 * M_PI)) * std::cos(v[0]) + 10;
}

Real f14(const Vector &v) {
    return std::sin(v[0] + v[1]) + std::pow(v[0] - v[1], 2) - 1.5 * v[0] + 2.5 * v[1] + 1;
}

Real f15(const Vector &v) {
    return 0.26 * (std::pow(v[0], 2) + std::pow(v[1], 2)) - 0.48 * v[0] * v[1];
}

Real f16(const Vector &v) {
    return (std::pow(std::sin(v[0] - v[1]), 2) * std::pow(std::sin(v[0] + v[1]), 2));
}

Real f17(const Vector &v) {
    return 1 / (1 + std::pow(v[0] - v[1], 2)) + std::sin((M_PI * v[1] + v[2]) / 2) +
           std::exp(std::pow((v[0] + v[1]) / v[1] - 2, 2));
}

Real f18(const Vector &v) {
    return v[0] / 4 + std::pow(v[0] * v[0] - 2 * v[0] + v[1] * v[1], 2);
}

Real f19(const Vector &v) {
    Real x = v[0], y = v[1];
    return (1.5-x*(1-y))*(1.5-x*(1-y))+(2.25-x*(1-y*y))*(2.25-x*(1-y*y))+(2.625-x*(1-y*y*y))*(2.625-x*(1-y*y*y));
}

Real f20(const Vector &v) {
// Не имеет глобального минимума - не использовать.
    Real x1 = v[0], x2 = v[1], x3 = v[2], x4 = v[3];
    return 100*(x2-x1*x1)*(x2-x1*x1)+(1-x1)*(1-x1)+90*(x4-x3*x3)*(x4-x3*x3)+(1-x3)*(1-x3)*(1-x3)+10.1*(x2-1)*(x2-1)+(x4-1)*(x4-1)+19.8*(x2-1)*(x4-1);
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
    for (unsigned int i=0; i < v.size() - 1; i++) {
        fun += v[i] * v[i] + 2 * v[i+1] * v[i+1] -0.3 * std::cos(3 * M_PI * v[i]) -
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
	return -sin(2*x[0]-0.5*M_PI)-3*cos(x[1])-0.5*x[0];
}

//Alpine 1
Real f48(const Vector &v)
{
	Real fun = 0;
	
	for (unsigned int i = 1; i < v.size(); i++)
	{
		fun = fun + std::abs(v[i]*sin(v[i])+0.1*v[i]);
	}
	return fun;
}
struct ControlPoint {
    Vector x;
    std::string type;

    ControlPoint(Vector x, std::string type)
        : x(x)
        , type(type)
    { }
};

void test_method(
    Method method,
    Function f,
    const StopCondition& stop_condition,
    const Matrix& start_points,
    const std::vector<ControlPoint>& expected
) {
    for (int i = 0; i < (int)start_points.size(); ++i) {
        const auto & p = start_points[i];
        
        IterationData iter_data = method(f, p, stop_condition);
        const auto x = iter_data.x_curr;
        const auto f_test =  iter_data.f_curr;
        const auto iter_counter = iter_data.iter_counter;

        auto best_eps_f = 1e9;
        auto best_eps_x = 1e9;
        auto best_point = expected.front();
        int best_id = 0;
        for (int curr_id = 0; curr_id < (int)expected.size(); ++curr_id) {
            const auto & expected_point = expected[curr_id];
            auto f_true = f(expected_point.x);
            auto temp_eps_f = f_test-f_true;
            auto temp_eps_x = norm(x-expected_point.x);
            if (temp_eps_x < best_eps_x) {
                best_point = expected_point;
                best_eps_x = temp_eps_x;
                best_eps_f = temp_eps_f;
                best_id = curr_id;
            }
        }
        fout_txt << "\tИз точки #" << i << " сходится к " << best_point.type << " #" << best_id
            << " (итераций: " << std::setw(8) << iter_counter << ") ";
        if (best_eps_f == 1e9) {
            fout_txt << "все очень плохо :(\n";
        } else {
            fout_txt << ", отклонение = " << std::fixed << std::setprecision(24) << std::setw(30) << best_eps_f << std::endl;
        }
        fout_txt << "\t\t                Начальная точка: " << p << std::endl;
        fout_txt << "\t\t  Предполагаемая точка минимума: " << best_point.x << std::endl;
        fout_txt << "\t\t      Полученная точка минимума: " << x << std::endl;
        fout_txt << "\t\tПредполагаемое значение функции: " 
            << std::setprecision(8) << std::fixed << std::setw(16) << f(best_point.x) << std::endl;
        fout_txt << "\t\t    Полученное значение функции: " 
            << std::setprecision(8) << std::fixed << std::setw(16) << f_test << std::endl << std::endl;
    }
}

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
//    Real mid = (left + right) / 2;
//    point = Vector(Dimensions, mid);
//    for (int i = 0; i < Dimensions; i += 2) {
//        point[i] = -mid;
//    }
//    start_points.push_back(point);
    return start_points;
}

bool example_stop_condition(const IterationData& iter_data) {
    return iter_data.iter_counter >= 100 || std::abs(iter_data.f_curr - iter_data.f_prev) < 1e-8;
}

void test1(Method method) {
    auto expected = std::vector<ControlPoint>{{{-1, -1}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 1 -----------------------------------\n\n";
    fout_txt << "01. Гладкая функция: f(x,y) = 1+x+y-xy+x^2+y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f1, example_stop_condition, start_points, expected);
}

void test2(Method method) {
    auto expected = std::vector<ControlPoint>{{{Real(-46) / 47, Real(-106) / 47}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 2 -----------------------------------\n\n";
    fout_txt << "02. Гладкая функция: f(x,y) = 1+7x+5y+0.5xy+3x^2+y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f2, example_stop_condition, start_points, expected);
}

void test3(Method method) {
    auto expected = std::vector<ControlPoint>{{{Real(-19) / 2, -5}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 3 -----------------------------------\n\n";
    fout_txt << "03. Гладкая функция: f(x,y) = 100+7x+5y-10xy+3x^2+10y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f3, example_stop_condition, start_points, expected);
}

void test4(Method method) {
    auto expected = std::vector<ControlPoint>{{{Real(-77900) / 39, Real(-14220)/13}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 4 -----------------------------------\n\n";
    fout_txt << "04. Гладкая функция: f(x,y) = 100+7x+5y-10.95xy+3x^2+10y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f4, example_stop_condition, start_points, expected);
}

void test5(Method method) {
    auto expected = std::vector<ControlPoint>{{{-0.25, -0.25, -0.25}, "Global Min"}};
    auto start_points = gen_start_points(3, -5, 5);
    fout_txt << "----------------------------------- Тест 5 -----------------------------------\n\n";
    fout_txt << "05. Гладкая функция: f(x,y,z) = 1+x+y+z+xy+xz+yz+x^2+y^2+z^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f5, example_stop_condition, start_points, expected);
}

void test6(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-0.582109, 0.525995}, "Global Min"}, {{0.582109, -0.525995}, "Global Min"}, {{0, 0}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 6 -----------------------------------\n\n";
    fout_txt << "06. Гладкая функция: f(x,y) = 10x^4+15y^4+15xy, имеющая два глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f6, example_stop_condition, start_points, expected);
}

void test7(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-1.08789, -1.00318}, "Global Min"}, {{1.08789, 1.00318}, "Global Min"}, {{0, 0}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 7 -----------------------------------\n\n";
    fout_txt << "07. Гладкая функция: f(x,y) = 10x^6+15y^6-20x^3y+xy^3, имеющая два глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f7, example_stop_condition, start_points, expected);
}

void test8(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{0, 0}, "Global Min"}, {{1, 1}, "Global Min"},
        {{-0.57735, -0.57735}, "Saddle point"}, {{0.57735, 0.57735}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 8 -----------------------------------\n\n";
    fout_txt << "08. Гладкая функция: f(x,y) = x^6+y^6-2(x^3y+xy^3)+x^2+y^2, имеющая три глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f8, example_stop_condition, start_points, expected);
}

void test9(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-1.34777, -1.34777}, "Global Min"}, {{1.34777, 1.34777}, "Global Min"}, {{0, 0}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 9 -----------------------------------\n\n";
    fout_txt << "09. Гладкая функция: f(x,y) = x^6+y^6-3(x^3y+xy^3)+x^2+y^2, имеющая два глобальных минимума и один локальный. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f9, example_stop_condition, start_points, expected);
}

void test10(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
        {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 10 ----------------------------------\n\n";
    fout_txt << "10. Гладкая функция: f(x,y) = x^6+y^6-2(x^3y+xy^3)+x^4+y^4-x^2-y^2, имеющая два глобальных минимума и два локальных. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f10, example_stop_condition, start_points, expected);
}

void test11(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{3, 2}, "Global Min"}, {{-3.77931, -3.28319}, "Global Min"},
        {{-2.80512, 3.13131}, "Global Min"}, {{3.58443, -1.84813}, "Global Min"},
        {{-0.270845, -0.923039}, " Local Max"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 11 ----------------------------------\n\n";
    fout_txt << "11. Гладкая функция Химмельблау: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2, имеющая четыре глобальных минимума и один локальный. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f11, example_stop_condition, start_points, expected);
}

void test15(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -10, 10);
    fout_txt << "----------------------------------- Тест 15 ----------------------------------\n\n";
    fout_txt << "15. Гладкая функция Матиаса: f(x,y)=0.26(x^2+y^2)-0.48xy, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f15, example_stop_condition, start_points, expected);
}

void test18(Method method) {
    auto start_points = gen_start_points(2, -1, 5);
    auto expected = std::vector<ControlPoint>{
        {{-0.029896, 0}, "Global Min"}
    };
    fout_txt << "----------------------------------- Тест 18 ----------------------------------\n\n";
    fout_txt << "18. Гладкая функция Зеттла: f(x,y)=x/4+(x^2-2x+y^2)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f18, example_stop_condition, start_points, expected);
}

void test19(Method method) {
    auto start_points = gen_start_points(2, -0.5, 0.5);
    auto expected = std::vector<ControlPoint>{
        {{3, 0.5}, "Global Min"}
    };
    fout_txt << "----------------------------------- Тест 19 ----------------------------------\n\n";
    fout_txt << "19. Гладкая функция Биля: f(x,y)=(xy-x+1.5)^2+(xy^2-x+2.25)^2+(x*y^3-x+2.625)^2, имеющая два глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f19, example_stop_condition, start_points, expected);
}

void test20(Method method) {
    auto start_points = Matrix{{-3,-1,-3,-1}};
    for (auto & it : gen_start_points(4, -0.5, 0.5)) {
        start_points.push_back(it);
    }
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    fout_txt << "----------------------------------- Тест 20 ----------------------------------\n\n";
    fout_txt << "20. Гладкая функция: f(x,y,z,t) = (x-1)^2+100(x^2-y)^2+10.1(y-1)^2+(z-1)^2+90(z^2-t)^2+10.1(t-1)^2+19.8(t-1)/y, имеет минимум на -10 < x_i < 10. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f20, example_stop_condition, start_points, expected);
}

void test21(Method method) {
    auto start_points = gen_start_points(2, -5, 5);
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    fout_txt << "----------------------------------- Тест 21 ----------------------------------\n\n";
    fout_txt << "21. Гладкая функция: f(x,y) = (y-x^2)^2+(1-x)^2, имеет один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f21, example_stop_condition, start_points, expected);
}

void test22(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 22 ----------------------------------\n\n";
    fout_txt << "22. Гладкая функция: f(x,y) = (y-x^2)^2+100(1-x)^2, имеет один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f22, example_stop_condition, start_points, expected);
}

void test23(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 23 ----------------------------------\n\n";
    fout_txt << "23. Гладкая функция: f(x,y) = 100(y-x^3)^2+(1-x)^2, имеет один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f23, example_stop_condition, start_points, expected);
}

void test24(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = Matrix{{-3,-1,0,1}};
    for (auto& it : gen_start_points(4, -5, 5)) {
        start_points.push_back(it);
    }
    fout_txt << "----------------------------------- Тест 24----------------------------------\n\n";
    fout_txt << "24. Гладкая функция: f(x,y,z,t) = (x+10y)^2+5(z-t)^2+(y-2z)^4+10(x-t)^4, имеет один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f24, example_stop_condition, start_points, expected);
}

void test25(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 1, 1, 1}, "Global Min"}
    };
    auto start_points = Matrix{{1, 2, 2, 2}};
    for (auto & it : gen_start_points(4, -5, 5)) {
        start_points.push_back(it);
    }
    
    fout_txt << "----------------------------------- Тест 25 ----------------------------------\n\n";
    fout_txt << "25. Гладкая функция: f(x,y,z,t) = (x^2-y+1)^4+100(y-z)^6+tg^4(z-t)+x^8+(t-1)^2, имеет один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f25, example_stop_condition, start_points, expected);
}

void test26_2(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 26_2 ----------------------------------\n\n";
    fout_txt << "26. Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f26, example_stop_condition, start_points, expected);
}

void test26_4(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "----------------------------------- Тест 26_4 ----------------------------------\n\n";
    fout_txt << "26. Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f26, example_stop_condition, start_points, expected);
}

void test26_8(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "----------------------------------- Тест 26_8 ----------------------------------\n\n";
    fout_txt << "26. Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f26, example_stop_condition, start_points, expected);
}

void test26_12(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "----------------------------------- Тест 26_12 ----------------------------------\n\n";
    fout_txt << "26. Гладкая функция:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f26, example_stop_condition, start_points, expected);
}

void test27_2(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 27_2 ----------------------------------\n\n";
    fout_txt << "27. Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f27, example_stop_condition, start_points, expected);
}

void test27_4(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "----------------------------------- Тест 27_4 ----------------------------------\n\n";
    fout_txt << "27. Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f27, example_stop_condition, start_points, expected);
}

void test27_8(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "----------------------------------- Тест 27_8 ----------------------------------\n\n";
    fout_txt << "27. Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f27, example_stop_condition, start_points, expected);
}

void test27_12(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "----------------------------------- Тест 27_12 ----------------------------------\n\n";
    fout_txt << "27. Гладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f27, example_stop_condition, start_points, expected);
}

void test28_2(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 28_2 ----------------------------------\n\n";
    fout_txt << "28. Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f28, example_stop_condition, start_points, expected);
}

void test28_4(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "----------------------------------- Тест 28_4 ----------------------------------\n\n";
    fout_txt << "28. Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f28, example_stop_condition, start_points, expected);
}

void test28_8(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "----------------------------------- Тест 28_8 ----------------------------------\n\n";
    fout_txt << "28. Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f28, example_stop_condition, start_points, expected);
}

void test28_12(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "----------------------------------- Тест 28_12 ----------------------------------\n\n";
    fout_txt << "28. Гладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f28, example_stop_condition, start_points, expected);
}

void test30_2(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "----------------------------------- Тест 30_2 ----------------------------------\n\n";
    fout_txt << "30. Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f30, example_stop_condition, start_points, expected);
}

void test30_4(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "----------------------------------- Тест 30_4 ----------------------------------\n\n";
    fout_txt << "30. Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f30, example_stop_condition, start_points, expected);
}

void test30_8(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "----------------------------------- Тест 30_8 ----------------------------------\n\n";
    fout_txt << "30. Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f30, example_stop_condition, start_points, expected);
}

void test30_12(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "----------------------------------- Тест 30_12 ----------------------------------\n\n";
    fout_txt << "30. Гладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f30, example_stop_condition, start_points, expected);
}

void test31(Method method)
{	
	auto expected = std::vector<ControlPoint>{
        {{     0,     -1}, "Global Min"},
        {{ 4.0/5,  1.0/5}, "Local Max"},
        {{-3.0/5, -2.0/5}, "Local Min"},
        {{ 6.0/5,  4.0/5}, "Local Min"},
        {{ 9.0/5,  1.0/5}, "Local Min"}
    };
    //auto start_points = gen_start_points(2, -5, 5);
	std::vector<std::vector<Real>> start_points = {{1, 1},{0.25, -0.125},{0.5, -1.5},{0, -1.2}};

    fout_txt << "----------------------------------- Тест 31 -----------------------------------\n\n";
    fout_txt << "31. Гладкая функция Голдштейна-Прайса: f(x,y) = [1+(x+y+1)^2(19-14x+3x^2-14y+6xy+3y^2)][30+(2x-3y)^2(18-32x+12x^2+48y-36xy+27y^2)], имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f31, example_stop_condition, start_points, expected);
}

void test32(Method method)
{
	auto expected = std::vector<ControlPoint>{{{-10, 1}, "Global Min"}};
	//auto start_points = gen_start_points(2, -5, 5);
    
    // -15 <= x <= -5
    // -3 <= y <= 3
	std::vector<std::vector<Real>> start_points = {{-15, -3}, {-15, 3}, {-5, -3}, {-5, 3}, {-10.5, 1.5}};
	
	fout_txt << "----------------------------------- Тест 32 -----------------------------------\n\n";
    fout_txt << "32. Негладкая функция Букин06: f(x,y) = 100*sqrt(|y-0.01*x^2|)+0.01|x+10|, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f32, example_stop_condition, start_points, expected);
	
}

void test33_2(Method method)
{
	auto expected = std::vector<ControlPoint>{{{0, 0}, "Global Min"}};
	auto start_points = gen_start_points(2, -1, 1);
	
	fout_txt << "----------------------------------- Тест 33_2 -----------------------------------\n\n";
    fout_txt << "33. Гладкая функция Cosine Mixture: f(x,y) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i))), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f33, example_stop_condition, start_points, expected);
	
}

void test33_4(Method method)
{
	auto expected = std::vector<ControlPoint>{{{0, 0, 0, 0}, "Global Min"}};
	auto start_points = gen_start_points(4, -1, 1);
	
	fout_txt << "----------------------------------- Тест 33_4 -----------------------------------\n\n";
    fout_txt << "33. Гладкая функция Cosine Mixture: f(x,y) = sum(x(i)^2)-0.1*sum(cos(5*pi*x(i))), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f33, example_stop_condition, start_points, expected);
	
}

void test34(Method method)
{	
	auto expected = std::vector<ControlPoint>{{{-10,0}, "Global Min"}};
	//auto start_points = gen_start_points(2, -5, 5);
	//std::vector<std::vector<Real>> start_points = {{-11, -3},{-10, 0.125},{-12, -1.5},{-15, 3}};
    std::vector<std::vector<Real>> start_points = {{-15, -3}, {-15, 3}, {-5, -3}, {-5, 3}, {-10.5, 1.5}};
    fout_txt << "----------------------------------- Тест 34 -----------------------------------\n\n";
    fout_txt << "34. Гладкая функция Букин02: f(x,y) = 100*(y-0.01*x^2+1)^2+0.01*(x+10)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f34, example_stop_condition, start_points, expected);
}

void test35(Method method) {
    // BartelsConn |x1^2+x2^2+x1*x2|+|sin(x1)|+|sin(x2)|
    auto expected = std::vector<ControlPoint>{{{0, 0}, "Global Min"}};
	auto start_points = gen_start_points(2, -5, 5);
	
	fout_txt << "----------------------------------- Тест 35 -----------------------------------\n\n";
    fout_txt << "35. Негладкая функция BartelsConn: f(x,y) = |x^2+y^2+x*y|+|sin(x)|+|sin(y)|, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f35, example_stop_condition, start_points, expected);
}

void test36(Method method) {
    auto expected = std::vector<ControlPoint>{{{0, 0}, "Global Min"}, {{2, 4}, "Global Min"}, {{1.464, -2.506}, "Global Min"}};
	//auto start_points = gen_start_points(2, -5, 5);
	std::vector<std::vector<Real>> start_points = {{3, 5.5}, {1.5, -3}, {0, -5}, {-1, 1.5}};
	
	fout_txt << "----------------------------------- Тест 36 -----------------------------------\n\n";
    fout_txt << "36. Гладкая функция Price04: f(x,y) = (2*x^3*y-y^3)^2+(6x-y^2+y)^2, имеющая три глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f36, example_stop_condition, start_points, expected);
}


void test37(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"},
        {{-1.74755, 0.873776}, " Local Min"},
        {{1.74755, -0.873776}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -1, 1);
    fout_txt << "----------------------------------- Тест 37 -----------------------------------\n\n";
    fout_txt << "37. Гладкая функция: f(x,y) = (2x^2 -1.05x^4 + x^6/6 + xy) + y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f37, example_stop_condition, start_points, expected);
}

void test38(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-0.0898, 0.7126}, "Global Min"},
        {{0.0898, -0.7126}, "Global Min"},
        {{-1.70361,  0.796084}, " Local Min"},
        {{ -1.6071, -0.568651}, " Local Min"},
        {{  1.6071,  0.568651}, " Local Min"},
        {{ 1.70361, -0.796084}, " Local Min"},
        {{-1.23023, -0.162335}, " Local Max"},
        {{ 1.23023,  0.162335}, " Local Max"}
    };
    std::vector<std::vector<Real>> start_points = {
        {0,0}, {1,1}, {-1, -1}, {-2, 0}, {0, 2}
    };
    fout_txt << "----------------------------------- Тест 38 -----------------------------------\n\n";
    fout_txt << "38. Гладкая функция: f(x,y) = (4-2.1*x1^2+x1^4/3)*x1^2+x1*x2+(-4+4*x2^2)*x2^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f38, example_stop_condition, start_points, expected);
}

void test39(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-M_PI, 12.275}, "Global Min"},
        {{M_PI, 2.275}, "Global Min"},
        {{9.42478, 2.475}, "Global Min"}
    };
    std::vector<std::vector<Real>> start_points = {
        {-5, 0}, {10, 0}, {-5, 15}, {10, 15}
    };
    fout_txt << "----------------------------------- Тест 39 -----------------------------------\n\n";
    fout_txt << "39. Гладкая функция: f(x,y) = (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+10, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f39, example_stop_condition, start_points, expected);
}

void test40(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{-3.2, 12.53}, "Global Min"}
    };
    std::vector<std::vector<Real>> start_points = {
        {-2, -2}, {-2, 11}, {10, -2}, {10, 10}
    };
    fout_txt << "----------------------------------- Тест 40 -----------------------------------\n\n";
    fout_txt << "40. Гладкая функция: f(x,y) = (-1.275*x1^2/M_PI^2+5*x1/M_PI+x2-6)^2 + (10-5/(4*M_PI))*cos(x1)*cos(x2)+log(x1^2+x2^2+1)+10, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f40, example_stop_condition, start_points, expected);
}

void test41(Method method) {
    auto expected = std::vector<ControlPoint>{{{0, 0}, "Global Min"}};
    std::vector<std::vector<Real>> start_points = {
        {-15, -5}, {-500, 500}
    };
    fout_txt << "----------------------------------- Тест 41 -----------------------------------\n\n";
    fout_txt << "41. Гладкая функция RotatedEllipse01: f(x,y) = 7x^2-6*sqrt(3)*x*y+13y^2, имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f41, example_stop_condition, start_points, expected);
}

void test42(Method method) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    std::vector<std::vector<Real>> start_points = {{-1,-1},{-1,1},{1,-1},{1,1}};
    fout_txt << "----------------------------------- Тест 42 -----------------------------------\n\n";
    fout_txt << "42. Гладкая функция: f(x,y) = x^2 + y^2 + 25[sin^2(x) + sin^2(y)], имеющая единственный глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f42, example_stop_condition, start_points, expected);
}

void test43(Method method) {
    // RotatedEllipse02 x*x-x*y+y*y
    auto expected = std::vector<ControlPoint>{{{0, 0}, "Global Min"}};
	std::vector<std::vector<Real>> start_points = {
        {500, 0}, {0, 500}, {500, 500}, {0, 0}
    };
	
	fout_txt << "----------------------------------- Тест 43 -----------------------------------\n\n";
    fout_txt << "43. Гладкая функция RotatedEllipse02: f(x,y) = x^2-x*y+y^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f43, example_stop_condition, start_points, expected);
}

void test44(Method method) {
    // Bird (x-y)^2+exp((1-sin(x))^2)cos(y)+exp((1-cos(y))^2)*sin(x)
    auto expected = std::vector<ControlPoint>{
        {{4.701055751981055, 3.152946019601391}, "Global Min"},
        {{-1.582142172055011, -3.130246799635430}, "Global Min"}
    };
	std::vector<std::vector<Real>> start_points = {
        {-M_PI, -M_PI}, {-M_PI, M_PI}, {M_PI, -M_PI}, {M_PI, M_PI}
    };
	
	fout_txt << "----------------------------------- Тест 44 -----------------------------------\n\n";
    fout_txt << "44. Гладкая функция Bird: f(x,y) = (x-y)^2+exp((1-sin(x))^2)cos(y)+exp((1-cos(y))^2)*sin(x), имеющая два глобальных минимума. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f44, example_stop_condition, start_points, expected);
}

void test45(Method method) {
    // Hosaki (1-8x+7x^2-7.0/3*x^3+1.0/4*x^4)*y^2*exp(-y)
    auto expected = std::vector<ControlPoint>{
        {{4,-1}, "Global Min"},
        {{4, 1}, "Global Min"},
        {{1, 1}, " Local Min"},
        {{1,-1}, " Local Min"}
    };
	std::vector<std::vector<Real>> start_points = {
        {0,0},{-3,-1},{-3,1},{3,-1},{3,1}
    };
	
	fout_txt << "----------------------------------- Тест 45 -----------------------------------\n\n";
    fout_txt << "45. Гладкая функция Hosaki: f(x,y) = (1-8x+7x^2-7.0/3*x^3+1.0/4*x^4)*y^2*exp(-y^2), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f45, example_stop_condition, start_points, expected);
}

void test46(Method method) {
    auto expected = std::vector<ControlPoint>{{{3.40919, -2.17143}, "Global Min"},
	{{-3.62315, -2.38415}, "Local Min"},
	{{-1.52071, 1.41228}, "Local Min"},
	{{2.27617, 0.864777}, "Local Min"}};
	//auto start_points = gen_start_points(2, -5, 5);
	std::vector<std::vector<Real>> start_points = {{3.5, -2.5}, {5, 0}, {4.125, -0.125}, {3.5, -1.5}};
	
	fout_txt << "----------------------------------- Тест 46 -----------------------------------\n\n";
    fout_txt << "46. Гладкая функция El-Attar-Vidyasagar-Dutta: f(x,y) = (x^2+y-10)^2+(x+y^2-7)^2+(x^2+y^3-1)^2, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f46, example_stop_condition, start_points, expected);
}

void test47(Method method) {
    auto expected = std::vector<ControlPoint>{{{1.697136443570341,0}, "Global Min"}};
	auto start_points = gen_start_points(2, -2, 2);
	//std::vector<std::vector<Real>> start_points = {{3.5, -2.5}, {5, 0}, {4.125, -0.125}, {3.5, -1.5}};
	
	fout_txt << "----------------------------------- Тест 47 -----------------------------------\n\n";
    fout_txt << "47. Гладкая функция Ursem01: f(x,y) = -sin(2x-0.5*pi)-3cos(y)-0.5x, имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f47, example_stop_condition, start_points, expected);
}

void test48(Method method) {
    auto expected = std::vector<ControlPoint>{{{0,0}, "Global Min"}};
	auto start_points = gen_start_points(2, -5, 5);
	//std::vector<std::vector<Real>> start_points = {{3.5, -2.5}, {5, 0}, {4.125, -0.125}, {3.5, -1.5}};
	
	fout_txt << "----------------------------------- Тест 48 -----------------------------------\n\n";
    fout_txt << "48. Негладкая функция Alpine01: f(x,y) = sum(abs(x(i)*sin(x(i))+0.1x(i))), имеющая один глобальный минимум. \nПодробнее в документе \"Тестовые функции\"\n\n";
    fout_txt << "Условие остановы: iter_counter >= 100 || |f_i-f_(i-1)| < 0.00000001\n\n";
    test_method(method, f48, example_stop_condition, start_points, expected);
}

void Test(Method method) {
    test1  (method);    std::cout << "1";     std::cout.flush();  fout_txt.flush();
    test2  (method);    std::cout << ", 2";   std::cout.flush();  fout_txt.flush();
    test3  (method);    std::cout << ", 3";   std::cout.flush();  fout_txt.flush();
    test4  (method);    std::cout << ", 4";   std::cout.flush();  fout_txt.flush();
    test5  (method);    std::cout << ", 5";   std::cout.flush();  fout_txt.flush();
    test6  (method);    std::cout << ", 6";   std::cout.flush();  fout_txt.flush();
    test7  (method);    std::cout << ", 7";   std::cout.flush();  fout_txt.flush();
    test8  (method);    std::cout << ", 8";   std::cout.flush();  fout_txt.flush();
    test9  (method);    std::cout << ", 9";   std::cout.flush();  fout_txt.flush();
    test10 (method);    std::cout << ", 10";  std::cout.flush();  fout_txt.flush();
    test11 (method);    std::cout << ", 11";  std::cout.flush();  fout_txt.flush();
//    test12 (method);    std::cout << ", 12";  std::cout.flush();  fout_txt.flush();
//    test13 (method);    std::cout << ", 13";  std::cout.flush();  fout_txt.flush();
//    test14 (method);    std::cout << ", 14";  std::cout.flush();  fout_txt.flush();
    test15 (method);    std::cout << ", 15";  std::cout.flush();  fout_txt.flush();
//    test16 (method);    std::cout << ", 16";  std::cout.flush();  fout_txt.flush();
//    test17 (method);    std::cout << ", 17";  std::cout.flush();  fout_txt.flush();
    test18 (method);    std::cout << ", 18";  std::cout.flush();  fout_txt.flush();
    test19 (method);    std::cout << ", 19";  std::cout.flush();  fout_txt.flush();
    // test20 (method);    std::cout << ", 20";  std::cout.flush();  fout_txt.flush();
    test21 (method);    std::cout << ", 21";  std::cout.flush();  fout_txt.flush();
    test22 (method);    std::cout << ", 22";  std::cout.flush();  fout_txt.flush();
    test23 (method);    std::cout << ", 23";  std::cout.flush();  fout_txt.flush();
    test24 (method);    std::cout << ", 24";  std::cout.flush();  fout_txt.flush();
    test25 (method);    std::cout << ", 25";  std::cout.flush();  fout_txt.flush();
    test26_2(method);   std::cout << ", 26_2";  std::cout.flush();  fout_txt.flush();
    test26_4(method);   std::cout << ", 26_4";  std::cout.flush();  fout_txt.flush();
    test26_8(method);   std::cout << ", 26_8";  std::cout.flush();  fout_txt.flush();
    test26_12(method);  std::cout << ", 26_12"; std::cout.flush();  fout_txt.flush();
    test27_2(method);   std::cout << ", 27_2";  std::cout.flush();  fout_txt.flush();
    test27_4(method);   std::cout << ", 27_4";  std::cout.flush();  fout_txt.flush();
    test27_8(method);   std::cout << ", 27_8";  std::cout.flush();  fout_txt.flush();
    test27_12(method);  std::cout << ", 27_12"; std::cout.flush();  fout_txt.flush();
    test28_2(method);   std::cout << ", 28_2";  std::cout.flush();  fout_txt.flush();
    test28_4(method);   std::cout << ", 28_4";  std::cout.flush();  fout_txt.flush();
    test28_8(method);   std::cout << ", 28_8";  std::cout.flush();  fout_txt.flush();
    test28_12(method);  std::cout << ", 28_12"; std::cout.flush();  fout_txt.flush();
//    test29_2(method);   std::cout << ", 29_2";  std::cout.flush();  fout_txt.flush();
//    test29_4(method);   std::cout << ", 29_4";  std::cout.flush();  fout_txt.flush();
//    test29_8(method);   std::cout << ", 29_8";  std::cout.flush();  fout_txt.flush();
//    test29_12(method);  std::cout << ", 29_12"; std::cout.flush();  fout_txt.flush();
    test30_2(method);   std::cout << ", 30_2";  std::cout.flush();  fout_txt.flush();
    test30_4(method);   std::cout << ", 30_4";  std::cout.flush();  fout_txt.flush();
    test30_8(method);   std::cout << ", 30_8";  std::cout.flush();  fout_txt.flush();
    test30_12(method);  std::cout << ", 30_12"; std::cout.flush();  fout_txt.flush();
	test31(method);  std::cout << ", 31_2"; std::cout.flush();  fout_txt.flush();
	test32(method);  std::cout << ", 32"; std::cout.flush();  fout_txt.flush();
	test33_2(method);  std::cout << ", 33_2"; std::cout.flush();  fout_txt.flush();
	test33_4(method);  std::cout << ", 33_4"; std::cout.flush();  fout_txt.flush();
	test34(method);  std::cout << ", 34"; std::cout.flush();  fout_txt.flush();
    test35(method);  std::cout << ", 35"; std::cout.flush();  fout_txt.flush();
	test36(method);  std::cout << ", 36"; std::cout.flush();  fout_txt.flush();
    test37(method);  std::cout << ", 37"; std::cout.flush();  fout_txt.flush();
    test38(method);  std::cout << ", 38"; std::cout.flush();  fout_txt.flush();
    test39(method);  std::cout << ", 39"; std::cout.flush();  fout_txt.flush();
    test40(method);  std::cout << ", 40"; std::cout.flush();  fout_txt.flush();
    test41(method);  std::cout << ", 41"; std::cout.flush();  fout_txt.flush();
    test42(method);  std::cout << ", 42"; std::cout.flush();  fout_txt.flush();
    test43(method);  std::cout << ", 43"; std::cout.flush();  fout_txt.flush();
    test44(method);  std::cout << ", 44"; std::cout.flush();  fout_txt.flush();
    test45(method);  std::cout << ", 45"; std::cout.flush();  fout_txt.flush();
	test46(method);  std::cout << ", 46"; std::cout.flush();  fout_txt.flush();
	test47(method);  std::cout << ", 47"; std::cout.flush();  fout_txt.flush();
	test48(method);  std::cout << ", 48"; std::cout.flush();  fout_txt.flush();
    std::cout << std::endl;
}

int main() {
    std::cout << std::endl;
    std::cout << "-- Start BFGS Method Tests. Results in test_bfgs.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_bfgs.txt");
    fout_txt << "BFGS method:\n\n";
    Test(bfgs);
    fout_txt.close();

    std::cout << std::endl;
    std::cout << "-- Start DFP Method Tests. Results in test_dfp.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_dfp.txt");
    fout_txt << "DFP method:\n\n";
    Test(dfp);
    fout_txt.close();

    std::cout << std::endl;
    std::cout << "-- Start Powell Method Tests. Results in test_powell.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_powell.txt");
    fout_txt << "Powell method:\n\n";
    Test(powell);
    fout_txt.close();
    
    std::cout << std::endl;
    std::cout << "-- Start Hessian Free Method Tests. Results in test_hessianfree.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_hessianfree.txt");
    fout_txt << "Hessian Free method:\n\n";
	Test(hessian_free);
    fout_txt.close();
    
    std::cout << std::endl;
    std::cout << "-- Start Nesterov Method Tests. Results in test_nesterov.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_nesterov.txt");
    fout_txt << "Nesterov method:\n\n";
    Test(nesterov);
    fout_txt.close();
    
    std::cout << std::endl;
    std::cout << "-- Finish testing... The results are written to a files:\n";
    std::cout << "\t * test_bfgs.txt\n";
    std::cout << "\t * test_dfp.txt\n";
    std::cout << "\t * test_powell.txt\n";
    std::cout << "\t * test_hessianfree.txt\n";
    std::cout << "\t * test_nesterov.txt\n";
    
    return 0;
}
