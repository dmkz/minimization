#include "hessian_free.hpp"
#include "nesterov.hpp"
#include "bfgs.hpp"
#include "dfp.hpp"
#include "powell.hpp"
#include <fstream>

/*
    Тестирование методов отдельно.
    Автор: Юрий Кондратов
*/

typedef void(*Method)(Function, Vector, BasicIterationObject*);

std::ofstream fout_txt, fout_csv;

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
    return std::pow(std::pow(v[0], 2) + v[1] - 11, 2) + std::pow(std::pow(v[1], 2) + v[0] - 7, 2);
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
    Real a = (1.5-x*(1-y));
    Real b = (2.25-x*(1-y*y));
    Real c = (2.625-x*(1-y*y*y));
    return a*a+b*b+c*c;
}

Real f20(const Vector &v) {
    return std::pow(v[0] - 1, 2) + 100 * std::pow(v[0] * v[0] - v[1], 2) +
           10.1 * std::pow(v[2] - 1, 2) + std::pow(v[2] - 1, 2) +
           90 * std::pow(v[2] * v[2] - v[3], 2) + 10.1 * std::pow(v[3] - 1, 2) +
           19.8 * (v[3] - 1) * (v[1] - 1);
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
    return v[0] + 10 * v[1] + 5 * std::pow(v[3] - v[2], 2) +
           std::pow(v[1] - 2 * v[2], 4) + 10 * std::pow(v[0] - v[3], 4);
}
Real f25(const Vector &v) {
    return std::pow(v[0] * v[0] - v[1], 4) + 100 * std::pow(v[1] - v[2], 6) +
           std::pow(std::tan(v[2]-v[3]), 4) + std::pow(v[0], 8) + std::pow(v[3] - 1, 2);
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
    const Matrix& start_points,
    const std::vector<ControlPoint>& expected,
    const int id_test,
    const std::string& method_title
) {
    for (int i = 0; i < (int)start_points.size(); ++i) {
        const auto & p = start_points[i];
        BasicIterationObject iter_object;
        method(f, p, &iter_object);
        const auto x = iter_object.get_x_curr();
        auto f_test =  iter_object.get_f_curr();
        auto iter_counter = iter_object.get_iter_counter();

        auto best_eps_f = 1e9;
        auto best_eps_x = 1e9;
        auto best_point = expected.front();
        int best_id = 0;
        for (int curr_id = 0; curr_id < (int)expected.size(); ++curr_id) {
            const auto & expected_point = expected[curr_id];
            auto f_true = f(expected_point.x);
            auto temp_eps_f = std::abs(f_test-f_true);
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
            fout_txt << "с точностью " << std::fixed << std::setprecision(24) << std::setw(30) << best_eps_f << std::endl;
        }
        fout_txt << "\t\tНачальное приближение: " << p << std::endl;
        fout_txt << "\t\t Предполагаемый ответ: " << best_point.x << std::endl;
        fout_txt << "\t\t     Полученный ответ: " << x << std::endl << std::endl;
        fout_csv << id_test << ", " << p << "," << x << ", " << best_point.x << ", " << method_title << ";" << std::endl;
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

void test1(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{-1, -1}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 1 ------------------\n";
    test_method(method, f1, start_points, expected, 1, method_title);
}

void test2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{Real(-46) / 47, Real(-106) / 47}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 2 ------------------\n";
    test_method(method, f2, start_points, expected, 2, method_title);
}

void test3(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{Real(-19) / 2, -5}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 3 ------------------\n";
    test_method(method, f3, start_points, expected, 3, method_title);
}

void test4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{Real(-77900) / 39, Real(-14220)/13}, "Global Min"}};
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 4 ------------------\n";
    test_method(method, f4, start_points, expected, 4, method_title);
}

void test5(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{-0.25, -0.25, -0.25}, "Global Min"}};
    auto start_points = gen_start_points(3, -5, 5);
    fout_txt << "------------------ Тест 5 ------------------\n";
    test_method(method, f5, start_points, expected, 5, method_title);
}

void test6(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-0.582109, 0.525995}, "Global Min"}, {{0.582109, -0.525995}, "Global Min"}, {{0, 0}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 6 ------------------\n";
    test_method(method, f6, start_points, expected, 6, method_title);
}

void test7(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1.08789, -1.00318}, "Global Min"}, {{1.08789, 1.00318}, "Global Min"}, {{0, 0}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 7 ------------------\n";
    test_method(method, f7, start_points, expected, 7, method_title);
}

void test8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{0, 0}, "Global Min"}, {{1, 1}, "Global Min"},
        {{-0.57735, -0.57735}, "Saddle point"}, {{0.57735, 0.57735}, "Saddle point"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 8 ------------------\n";
    test_method(method, f8, start_points, expected, 8, method_title);
}

void test9(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1.34777, -1.34777}, "Global Min"}, {{1.34777, 1.34777}, "Global Min"}, {{0, 0}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 9 ------------------\n";
    test_method(method, f9, start_points, expected, 9, method_title);
}

void test10(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
        {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 10 -----------------\n";
    test_method(method, f10, start_points, expected, 10, method_title);
}

void test11(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
        {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 11 -----------------\n";
    test_method(method, f11, start_points, expected, 11, method_title);
}

void test15(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -10, 10);
    fout_txt << "------------------ Тест 15 -----------------\n";
    test_method(method, f15, start_points, expected, 15, method_title);
}

void test18(Method method, const std::string& method_title) {
    auto start_points = gen_start_points(2, -1, 5);
    auto expected = std::vector<ControlPoint>{
        {{-0.029896, 0}, "Global Min"}
    };
    fout_txt << "------------------ Тест 18 -----------------\n";
    test_method(method, f18, start_points, expected, 18, method_title);
}

void test19(Method method, const std::string& method_title) {
    auto start_points = gen_start_points(2, -0.5, 0.5);
    auto expected = std::vector<ControlPoint>{
        {{3, 0.5}, "Global Min"}
    };
    fout_txt << "------------------ Тест 19 -----------------\n";
    test_method(method, f19, start_points, expected, 19, method_title);
}

void test20(Method method, const std::string& method_title) {
    auto start_points = gen_start_points(4, -10, 10);
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    fout_txt << "------------------ Тест 20 -----------------\n";
    test_method(method, f20, start_points, expected, 20, method_title);
}

void test21(Method method, const std::string& method_title) {
    auto start_points = gen_start_points(2, -5, 5);
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    fout_txt << "------------------ Тест 21 -----------------\n";
    test_method(method, f21, start_points, expected, 21, method_title);
}

void test22(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 22 -----------------\n";
    test_method(method, f22, start_points, expected, 22, method_title);
}

void test23(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 23 -----------------\n";
    test_method(method, f23, start_points, expected, 23, method_title);
}

void test24(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 24-----------------\n";
    test_method(method, f24, start_points, expected, 24, method_title);
}

void test25(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 1, 1, 1}, "Global Min"}
    };
//    auto start_points = Matrix{{1, 2, 2, 2});
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 25 -----------------\n";
    test_method(method, f25, start_points, expected, 25, method_title);
}

void test26_2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 26_2 -----------------\n";
    test_method(method, f26, start_points, expected, 26, method_title);
}

void test26_4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 26_4 -----------------\n";
    test_method(method, f26, start_points, expected, 26, method_title);
}

void test26_8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "------------------ Тест 26_8 -----------------\n";
    test_method(method, f26, start_points, expected, 26, method_title);
}

void test26_12(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "------------------ Тест 26_12 -----------------\n";
    test_method(method, f26, start_points, expected, 26, method_title);
}

void test27_2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 27_2 -----------------\n";
    test_method(method, f27, start_points, expected, 27, method_title);
}

void test27_4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 27_4 -----------------\n";
    test_method(method, f27, start_points, expected, 27, method_title);
}

void test27_8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "------------------ Тест 27_8 -----------------\n";
    test_method(method, f27, start_points, expected, 27, method_title);
}

void test27_12(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "------------------ Тест 27_12 -----------------\n";
    test_method(method, f27, start_points, expected, 27, method_title);
}

void test28_2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 28_2 -----------------\n";
    test_method(method, f28, start_points, expected, 28, method_title);
}

void test28_4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 28_4 -----------------\n";
    test_method(method, f28, start_points, expected, 28, method_title);
}

void test28_8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "------------------ Тест 28_8 -----------------\n";
    test_method(method, f28, start_points, expected, 28, method_title);
}

void test28_12(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "------------------ Тест 28_12 -----------------\n";
    test_method(method, f28, start_points, expected, 28, method_title);
}

void test30_2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(2, -5, 5);
    fout_txt << "------------------ Тест 30_2 -----------------\n";
    test_method(method, f30, start_points, expected, 30, method_title);
}

void test30_4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(4, -5, 5);
    fout_txt << "------------------ Тест 30_4 -----------------\n";
    test_method(method, f30, start_points, expected, 30, method_title);
}

void test30_8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(8, -5, 5);
    fout_txt << "------------------ Тест 30_8 -----------------\n";
    test_method(method, f30, start_points, expected, 30, method_title);
}

void test30_12(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Global Min"}
    };
    auto start_points = gen_start_points(12, -5, 5);
    fout_txt << "------------------ Тест 30_12 -----------------\n";
    test_method(method, f30, start_points, expected, 30, method_title);
}

void Test(Method method, const std::string& method_title) {
    test1  (method, method_title);    std::cout << "1";     std::cout.flush();  fout_txt.flush();
    test2  (method, method_title);    std::cout << ", 2";   std::cout.flush();  fout_txt.flush();
    test3  (method, method_title);    std::cout << ", 3";   std::cout.flush();  fout_txt.flush();
    test4  (method, method_title);    std::cout << ", 4";   std::cout.flush();  fout_txt.flush();
    test5  (method, method_title);    std::cout << ", 5";   std::cout.flush();  fout_txt.flush();
    test6  (method, method_title);    std::cout << ", 6";   std::cout.flush();  fout_txt.flush();
    test7  (method, method_title);    std::cout << ", 7";   std::cout.flush();  fout_txt.flush();
    test8  (method, method_title);    std::cout << ", 8";   std::cout.flush();  fout_txt.flush();
    test9  (method, method_title);    std::cout << ", 9";   std::cout.flush();  fout_txt.flush();
    test10 (method, method_title);    std::cout << ", 10";  std::cout.flush();  fout_txt.flush();
    test11 (method, method_title);    std::cout << ", 11";  std::cout.flush();  fout_txt.flush();
//    test12 (method, method_title);    std::cout << ", 12";  std::cout.flush();  fout_txt.flush();
//    test13 (method, method_title);    std::cout << ", 13";  std::cout.flush();  fout_txt.flush();
//    test14 (method, method_title);    std::cout << ", 14";  std::cout.flush();  fout_txt.flush();
    test15 (method, method_title);    std::cout << ", 15";  std::cout.flush();  fout_txt.flush();
//    test16 (method, method_title);    std::cout << ", 16";  std::cout.flush();  fout_txt.flush();
//    test17 (method, method_title);    std::cout << ", 17";  std::cout.flush();  fout_txt.flush();
    test18 (method, method_title);    std::cout << ", 18";  std::cout.flush();  fout_txt.flush();
    test19 (method, method_title);    std::cout << ", 19";  std::cout.flush();  fout_txt.flush();
    test20 (method, method_title);    std::cout << ", 20";  std::cout.flush();  fout_txt.flush();
    test21 (method, method_title);    std::cout << ", 21";  std::cout.flush();  fout_txt.flush();
    test22 (method, method_title);    std::cout << ", 22";  std::cout.flush();  fout_txt.flush();
    test23 (method, method_title);    std::cout << ", 23";  std::cout.flush();  fout_txt.flush();
    test24 (method, method_title);    std::cout << ", 24";  std::cout.flush();  fout_txt.flush();
    test25 (method, method_title);    std::cout << ", 25";  std::cout.flush();  fout_txt.flush();
    test26_2(method, method_title);   std::cout << ", 26_2";  std::cout.flush();  fout_txt.flush();
    test26_4(method, method_title);   std::cout << ", 26_4";  std::cout.flush();  fout_txt.flush();
    test26_8(method, method_title);   std::cout << ", 26_8";  std::cout.flush();  fout_txt.flush();
    test26_12(method, method_title);  std::cout << ", 26_12"; std::cout.flush();  fout_txt.flush();
    test27_2(method, method_title);   std::cout << ", 27_2";  std::cout.flush();  fout_txt.flush();
    test27_4(method, method_title);   std::cout << ", 27_4";  std::cout.flush();  fout_txt.flush();
    test27_8(method, method_title);   std::cout << ", 27_8";  std::cout.flush();  fout_txt.flush();
    test27_12(method, method_title);  std::cout << ", 27_12"; std::cout.flush();  fout_txt.flush();
    test28_2(method, method_title);   std::cout << ", 28_2";  std::cout.flush();  fout_txt.flush();
    test28_4(method, method_title);   std::cout << ", 28_4";  std::cout.flush();  fout_txt.flush();
    test28_8(method, method_title);   std::cout << ", 28_8";  std::cout.flush();  fout_txt.flush();
    test28_12(method, method_title);  std::cout << ", 28_12"; std::cout.flush();  fout_txt.flush();
//    test29_2(method, method_title);   std::cout << ", 29_2";  std::cout.flush();  fout_txt.flush();
//    test29_4(method, method_title);   std::cout << ", 29_4";  std::cout.flush();  fout_txt.flush();
//    test29_8(method, method_title);   std::cout << ", 29_8";  std::cout.flush();  fout_txt.flush();
//    test29_12(method, method_title);  std::cout << ", 29_12"; std::cout.flush();  fout_txt.flush();
    test30_2(method, method_title);   std::cout << ", 30_2";  std::cout.flush();  fout_txt.flush();
    test30_4(method, method_title);   std::cout << ", 30_4";  std::cout.flush();  fout_txt.flush();
    test30_8(method, method_title);   std::cout << ", 30_8";  std::cout.flush();  fout_txt.flush();
    test30_12(method, method_title);  std::cout << ", 30_12"; std::cout.flush();  fout_txt.flush();
    std::cout << std::endl;
}

int main() {
    fout_csv.open("test_result.csv");

    std::cout << "-- Start BFGS Method Tests. Results in test_bfgs.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_bfgs.txt");
    fout_txt << "BFGS method:\n\n";
    Test(bfgs, "BFGS");
    fout_txt.close();

    std::cout << "-- Start DFP Method Tests. Results in test_dfp.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_dfp.txt");
    fout_txt << "DFP method:\n\n";
    Test(dfp, "DFP");
    fout_txt.close();

    std::cout << "-- Start Powell Method Tests. Results in test_powell.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_powell.txt");
    fout_txt << "Powell method:\n\n";
    Test(powell, "Powell");
    fout_txt.close();

    std::cout << "-- Start Hessian Free Method Tests. Results in test_hessianfree.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_hessianfree.txt");
    fout_txt << "Hessian Free method:\n\n";
	Test(hessian_free, "Hessian Free");
    fout_txt.close();
    
    std::cout << "-- Start Nesterov Method Tests. Results in test_nesterov.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout_txt.open("test_nesterov.txt");
    fout_txt << "Nesterov method:\n\n";
    Test(nesterov, "Nesterov");
    fout_txt.close();
    
    fout_csv.close();

    return 0;
}
