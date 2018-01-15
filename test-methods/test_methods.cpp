#include "hessian_free.hpp"
#include "nesterov.hpp"
#include "bfgs.hpp"
#include <fstream>

// typedef Vector(*Method)(Function, Vector);

// Тестовые функции
ld f1(const Vector& x) {
    return 1 + x[0] + x[1] - x[0] * x[1] + x[0] * x[0] + x[1] * x[1];
}

ld f2(const Vector& x) {
    return 1 + 7 * x[0] + 5 * x[1] + 0.5 * x[0] * x[1] + 3 * x[0] * x[0] + x[1] * x[1];
}

ld f3(const Vector& x) {
    return 100 + 7 * x[0] + 5 * x[1] - 10 * x[0] * x[1] + 3 * x[0] * x[0] + 10 * x[1] * x[1];
}

ld f4(const Vector& x) {
    return 100 + 7 * x[0] + 5 * x[1] - 10.95 * x[0] * x[1] + 3 * x[0] * x[0] + 10 * x[1] * x[1];
}

ld f5(const Vector& x) {
    return 1+x[0]+x[1]+x[2] + x[0] * x[1] + x[0] * x[2] + x[1] * x[2] + x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

ld f6(const Vector& x) {
    return 10 * std::pow(x[0], 4) + 15 * std::pow(x[1], 4) + 15 * x[0] * x[1];
}

ld f7(const Vector& v) {
    ld x = v[0];
    ld y = v[1];
    return 10 * std::pow(x, 6) + 15 * std::pow(y, 6) - 20 * (std::pow(x, 3) * y + x * std::pow(y, 3));
}

ld f8(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 2 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + x*x + y*y;
}

ld f9(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 3 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + x*x + y*y;
}

ld f10(const Vector& v) {
    auto x = v[0], y = v[1];
    return std::pow(x, 6) + std::pow(y, 6) - 2 * (std::pow(x, 3) * y + x * std::pow(y, 3)) + std::pow(x, 4) + std::pow(y, 4) - x*x - y*y;
}

ld f11(const Vector& v) {
    ld fun = std::pow(v[0] - 1, 2) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::pow(v[i] - 2 * std::pow(v[i-1], 2) + 1, 2);
    }
    return fun;
}

ld f12(const Vector &v) {
    ld fun = std::pow(v[0] - 1, 2) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::abs(v[i] - 2 * std::pow(v[i-1], 2) + 1);
    }
    return fun;
}

ld f13(const Vector &v) {
    ld fun = std::abs(v[0] - 1) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::abs(v[i] - 2 * std::pow(v[i-1], 2) + 1);
    }
    return fun;
}

ld f14(const Vector &v) {
    ld fun = 0;
    for (unsigned int i=1; i < v.size(); i += 2) {
        fun += 100 * std::pow(std::pow(v[i-1], 2) - v[i], 2) + std::pow(v[i-1] - 1, 2);
    }
    return fun;
}

ld f15(const Vector &v) {
    ld fun = 0;
    for (unsigned int i=1; i < v.size(); i += 2) {
        fun += 10 * std::abs(std::pow(v[i-1], 2) - v[i]) + std::abs(v[i-1] - 1);
    }
    return fun;
}

ld f16(const Vector &v) {
    ld fun = 0;
    for (int i=0; i < 20; ++i) {
        fun = std::max(fun, std::pow(v[i], 2));
    }
    return fun;
}

ld f17(const Vector &v) {
    ld fun = 0;
    ld g_fun, t_val;
    for (int i=0; i < 21; ++i) {
        t_val = 0.25 + 0.75 * i / 20;
        g_fun = std::abs(v[3] - std::pow(v[0] *
                std::pow(t_val, 2) + v[1] * t_val + v[2], 2) - std::sqrt(t_val));
        fun = std::max(fun, g_fun);
    }
    return fun;
}

ld f18(const Vector &v) {
    ld fun = 0;
    ld g_fun, t_val;
    for (int i=0; i < 51; ++i) {
        t_val = 0.1 * i;
        g_fun = 0.5 * std::exp(-t_val) - std::exp(-2 * t_val) +
                0.5 * std::exp(-3 * t_val) + 1.5 * std::exp(-1.5 * t_val) * std::sin(7 * t_val) +
                std::exp(-2.5 * t_val) * std::sin(5 * t_val);
        fun += std::abs(v[0] * exp(-v[1] * t_val) * std::cos(v[2] * t_val + v[3]) +
                        v[4] * std::exp(-v[5] * t_val) - g_fun);

    }
    return fun;
}

ld f19(const Vector &v) {
    ld fun = 0;
    ld g_fun;
    for (int i = 0; i < 50; ++i) {
        g_fun = 0;
        for (int j=0; j < 50; ++j) {
            g_fun += 1 / (i + j + 1);
        }
        g_fun *= std::abs(v[i]);
        fun += g_fun;
    }
    return fun;
}

ld f20(const Vector &v) {
    ld g_fun = 1e-8 * std::pow(v[0], 2) + std::pow(v[2], 2) + 4 * std::pow(v[3], 2) +
               std::pow(v[5], 2) + std::pow(v[6], 2) + std::pow(v[7], 2) +
               std::pow(v[8], 2) + std::pow(v[9], 2);
    ld fun = std::exp(g_fun + 2 * std::pow(v[1], 2));
    fun = std::max(fun, std::exp(g_fun - 2 * std::pow(v[1], 2)));
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

void test_method(Method method, Function f, const Matrix& start_points, const std::vector<ControlPoint>& expected, std::ofstream& fout) {
    const int ITERATIONS_LIMIT = 1000;
    for (int i = 0; i < (int)start_points.size(); ++i) {
        const auto & p = start_points[i];
        auto res = method(f, p, ITERATIONS_LIMIT);
        const auto& x = res.first;
        auto f_test = f(x);
        auto best_eps = 1e9;
        auto best_point = expected.front();
        int best_id = 0;
        for (int curr_id = 0; curr_id < (int)expected.size(); ++curr_id) {
            const auto & expected_point = expected[curr_id];
            auto f_true = f(expected_point.x);
            auto temp_eps = std::abs(f_test-f_true);
            if (best_eps > temp_eps) {
                best_point = expected_point;
                best_eps = temp_eps;
                best_id = curr_id;
            }
        }
        fout << "\tИз точки #" << i << " сходится к " << best_point.type << " #" << best_id
            << " (итераций: " << std::setw(8) << res.second << ") ";
        if (best_eps == 1e9) {
            fout << "все очень плохо :(\n";
        } else {
            fout << "с точностью " << std::fixed << std::setprecision(24) << std::setw(30) << best_eps << std::endl;
        }
        fout << "\t\tНачальное приближение: " << p << std::endl;
        fout << "\t\t Предполагаемый ответ: " << best_point.x << std::endl;
        fout << "\t\t     Полученный ответ: " << x << std::endl << std::endl;
    }
}

void test1(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{{-1, -1}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 1 ------------------\n";
    test_method(method, f1, start_points, expected, fout);
}

void test2(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{{ld(-46) / 47, ld(-106) / 47}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 2 ------------------\n";
    test_method(method, f2, start_points, expected, fout);
}

void test3(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{{ld(-19) / 2, -5}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 3 ------------------\n";
    test_method(method, f3, start_points, expected, fout);
}

void test4(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{{ld(-77900) / 39, ld(-14220)/13}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 4 ------------------\n";
    test_method(method, f4, start_points, expected, fout);
}

void test5(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{{-0.25, -0.25, -0.25}, "Global Min"}};
    auto start_points = Matrix{{-20, -20, -20}, {-20, -20, 20}, {-20, 20, -20}, {-20, 20, 20}, {20, -20, -20}, {20, -20, 20}, {20, 20, -20}, {20, 20, 20}, {0, 0, 0}};
    fout << "------------------ Тест 5 ------------------\n";
    test_method(method, f5, start_points, expected, fout);
}

void test6(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{
        {{-0.582109, 0.525995}, "Global Min"}, {{0.582109, -0.525995}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0.5, 0.5}};
    fout << "------------------ Тест 6 ------------------\n";
    test_method(method, f6, start_points, expected, fout);
}

void test7(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{
        {{-1.08789, -1.00318}, "Global Min"}, {{1.08789, 1.00318}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 7 ------------------\n";
    test_method(method, f7, start_points, expected, fout);
}

void test8(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{0, 0}, "Global Min"}, {{1, 1}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 8 ------------------\n";
    test_method(method, f8, start_points, expected, fout);
}

void test9(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{
        {{-1.34777, -1.34777}, "Global Min"}, {{1.34777, 1.34777}, "Global Min"}, {{0, 0}, " Local Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout << "------------------ Тест 9 ------------------\n";
    test_method(method, f9, start_points, expected, fout);
}

void test10(Method method, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
        {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0.5, 0.5}};
    fout << "------------------ Тест 10 -----------------\n";
    test_method(method, f10, start_points, expected, fout);
}

void test_Nesterov_n_Rosenbrock(Method method, int type, int Dimensions, std::ofstream& fout) {
    auto expected = std::vector<ControlPoint>{{Vector(Dimensions, 1), "Global Min"}};
    Vector point(Dimensions, 1);
    point[0] = -1;
    auto start_points = Matrix{point};
    start_points.push_back(Vector(Dimensions, -1));
    start_points.push_back(Vector(Dimensions, 0));
    point = Vector(Dimensions, 10);
    for (int i = Dimensions / 2; i < Dimensions; ++i) {
        point[i] = -10;
    }
    start_points.push_back(point);
    point = Vector(Dimensions, -1);
    for (int i = 0; i < Dimensions; i += 2) {
        point[i] = 0;
    }
    start_points.push_back(point);
    switch (type) {
    case 1:
        test_method(method, f11, start_points, expected, fout);
        break;
    case 2:
        test_method(method, f12, start_points, expected, fout);
        break;
    case 3:
        test_method(method, f13, start_points, expected, fout);
        break;
    case 4:
        test_method(method, f14, start_points, expected, fout);
        break;
    case 5:
        test_method(method, f15, start_points, expected, fout);
    }
}

void test11(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 11 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 2, fout);
}

void test12(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 12 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 3, fout);
}

void test13(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 13 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 4, fout);
}

void test14(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 14 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 5, fout);
}

void test15(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 15 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 6, fout);
}

void test16(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 16 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 2, fout);
}

void test17(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 17 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 5, fout);
}

void test18(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 18 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 10, fout);
}

void test19(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 19 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 100, fout);
}

void test20(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 20 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 1000, fout);
}

void test21(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 21 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 2, fout);
}

void test22(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 22 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 5, fout);
}

void test23(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 23 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 10, fout);
}

void test24(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 24 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 100, fout);
}

void test25(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 25 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 1000, fout);
}

void test26(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 26 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 2, fout);
}

void test27(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 27 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 6, fout);
}

void test28(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 28 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 10, fout);
}

void test29(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 29 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 100, fout);
}

void test30(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 30 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 1000, fout);
}

void test31(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 31 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 2, fout);
}

void test32(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 32 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 6, fout);
}

void test33(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 33 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 10, fout);
}

void test34(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 34 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 100, fout);
}

void test35(Method method, std::ofstream& fout) {
    fout << "------------------ Тест 35 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 1000, fout);
}

void Test(Method method, std::ofstream& fout) {
    test1  (method, fout);    std::cout << "1";     std::cout.flush();  fout.flush();
    test2  (method, fout);    std::cout << ", 2";   std::cout.flush();  fout.flush();
    test3  (method, fout);    std::cout << ", 3";   std::cout.flush();  fout.flush();
    test4  (method, fout);    std::cout << ", 4";   std::cout.flush();  fout.flush();
    test5  (method, fout);    std::cout << ", 5";   std::cout.flush();  fout.flush();
    test6  (method, fout);    std::cout << ", 6";   std::cout.flush();  fout.flush();
    test7  (method, fout);    std::cout << ", 7";   std::cout.flush();  fout.flush();
    test8  (method, fout);    std::cout << ", 8";   std::cout.flush();  fout.flush();
    test9  (method, fout);    std::cout << ", 9";   std::cout.flush();  fout.flush();
    test10 (method, fout);    std::cout << ", 10";  std::cout.flush();  fout.flush();
    test11 (method, fout);    std::cout << ", 11";  std::cout.flush();  fout.flush();
    test12 (method, fout);    std::cout << ", 12";  std::cout.flush();  fout.flush();
    test13 (method, fout);    std::cout << ", 13";  std::cout.flush();  fout.flush();
    test14 (method, fout);    std::cout << ", 14";  std::cout.flush();  fout.flush();
    test15 (method, fout);    std::cout << ", 15";  std::cout.flush();  fout.flush();
    test16 (method, fout);    std::cout << ", 16";  std::cout.flush();  fout.flush();
    test17 (method, fout);    std::cout << ", 17";  std::cout.flush();  fout.flush();
    test18 (method, fout);    std::cout << ", 18";  std::cout.flush();  fout.flush();
    test19 (method, fout);    std::cout << ", 19";  std::cout.flush();  fout.flush();
    // test20 (method);       std::cout << ", 20";  std::cout.flush();  fout.flush();
    test21 (method, fout);    std::cout << ", 21";  std::cout.flush();  fout.flush();
    test22 (method, fout);    std::cout << ", 22";  std::cout.flush();  fout.flush();
    test23 (method, fout);    std::cout << ", 23";  std::cout.flush();  fout.flush();
    test24 (method, fout);    std::cout << ", 24";  std::cout.flush();  fout.flush();
    // test25 (method);       std::cout << ", 25";  std::cout.flush();  fout.flush();
    test26 (method, fout);    std::cout << ", 26";  std::cout.flush();  fout.flush();
    test27 (method, fout);    std::cout << ", 27";  std::cout.flush();  fout.flush();
    test28 (method, fout);    std::cout << ", 28";  std::cout.flush();  fout.flush();
    test29 (method, fout);    std::cout << ", 29";  std::cout.flush();  fout.flush();
    // test30 (method);       std::cout << ", 30";  std::cout.flush();  fout.flush();
    test31 (method, fout);    std::cout << ", 31";  std::cout.flush();  fout.flush();
    test32 (method, fout);    std::cout << ", 32";  std::cout.flush();  fout.flush();
    test33 (method, fout);    std::cout << ", 33";  std::cout.flush();  fout.flush();
    test34 (method, fout);    std::cout << ", 34";  std::cout.flush();  fout.flush();
    // test35 (method);       std::cout << ", 35";  std::cout.flush();  fout.flush();
    std::cout << std::endl;
}

int main() {    
    std::ofstream fout;
    
    std::cout << "-- Start BFGS Method Tests. Results in test_bfgs.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout.open("test_bfgs.txt");
    fout << "BFGS method:\n\n";
    Test(bfgs, fout);
    fout.close();
    
    std::cout << "-- Start Hessian Free Method Tests. Results in test_hessianfree.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout.open("test_hessianfree.txt");
    fout << "Hessian Free method:\n\n";
	Test(hessian_free, fout);
    fout.close();
    
    std::cout << "-- Start Nesterov Method Tests. Results in test_nesterov.txt" << std::endl;
    std::cout << "-- Tests: ";
    fout.open("test_nesterov.txt");
    fout << "Nesterov method:\n\n";
    Test(nesterov, fout);
    fout.close();

    return 0;
}
