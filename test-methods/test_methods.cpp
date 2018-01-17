#include "hessian_free.hpp"
#include "nesterov.hpp"
#include "bfgs.hpp"
#include <fstream>

// typedef Vector(*Method)(Function, Vector);

std::ofstream fout_txt, fout_csv;

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

void test_method(
    Method method, 
    Function f, 
    const Matrix& start_points, 
    const std::vector<ControlPoint>& expected, 
    const int id_test,
    const std::string& method_title
) {
    const int ITERATIONS_LIMIT = 1000;
    for (int i = 0; i < (int)start_points.size(); ++i) {
        const auto & p = start_points[i];
        auto res = method(f, p, ITERATIONS_LIMIT);
        const auto& x = res.first;
        auto f_test = f(x);
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
            << " (итераций: " << std::setw(8) << res.second << ") ";
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

void test1(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{-1, -1}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 1 ------------------\n";
    test_method(method, f1, start_points, expected, 1, method_title);
}

void test2(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{ld(-46) / 47, ld(-106) / 47}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 2 ------------------\n";
    test_method(method, f2, start_points, expected, 2, method_title);
}

void test3(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{ld(-19) / 2, -5}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 3 ------------------\n";
    test_method(method, f3, start_points, expected, 3, method_title);
}

void test4(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{ld(-77900) / 39, ld(-14220)/13}, "Global Min"}};
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 4 ------------------\n";
    test_method(method, f4, start_points, expected, 4, method_title);
}

void test5(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{{{-0.25, -0.25, -0.25}, "Global Min"}};
    auto start_points = Matrix{{-20, -20, -20}, {-20, -20, 20}, {-20, 20, -20}, {-20, 20, 20}, {20, -20, -20}, {20, -20, 20}, {20, 20, -20}, {20, 20, 20}, {0, 0, 0}};
    fout_txt << "------------------ Тест 5 ------------------\n";
    test_method(method, f5, start_points, expected, 5, method_title);
}

void test6(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-0.582109, 0.525995}, "Global Min"}, {{0.582109, -0.525995}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0.5, 0.5}};
    fout_txt << "------------------ Тест 6 ------------------\n";
    test_method(method, f6, start_points, expected, 6, method_title);
}

void test7(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1.08789, -1.00318}, "Global Min"}, {{1.08789, 1.00318}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 7 ------------------\n";
    test_method(method, f7, start_points, expected, 7, method_title);
}

void test8(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{0, 0}, "Global Min"}, {{1, 1}, "Global Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 8 ------------------\n";
    test_method(method, f8, start_points, expected, 8, method_title);
}

void test9(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1.34777, -1.34777}, "Global Min"}, {{1.34777, 1.34777}, "Global Min"}, {{0, 0}, " Local Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0, 0}};
    fout_txt << "------------------ Тест 9 ------------------\n";
    test_method(method, f9, start_points, expected, 9, method_title);
}

void test10(Method method, const std::string& method_title) {
    auto expected = std::vector<ControlPoint>{
        {{-1, -1}, "Global Min"}, {{1, 1}, "Global Min"},
        {{0, 0}, " Local Max"}, {{-0.39332, 0.39332}, " Local Min"}, {{0.39332, -0.39332}, " Local Min"}
    };
    auto start_points = Matrix{{-20, -20}, {-20, 20}, {20, -20}, {20, 20}, {0.5, 0.5}};
    fout_txt << "------------------ Тест 10 -----------------\n";
    test_method(method, f10, start_points, expected, 10, method_title);
}

void test_Nesterov_n_Rosenbrock(Method method, int type, int Dimensions, int id_test, const std::string& method_title) {
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
        test_method(method, f11, start_points, expected, id_test, method_title);
        break;
    case 2:
        test_method(method, f12, start_points, expected, id_test, method_title);
        break;
    case 3:
        test_method(method, f13, start_points, expected, id_test, method_title);
        break;
    case 4:
        test_method(method, f14, start_points, expected, id_test, method_title);
        break;
    case 5:
        test_method(method, f15, start_points, expected, id_test, method_title);
    }
}

void test11(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 11 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 2, 11, method_title);
}

void test12(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 12 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 3, 12, method_title);
}

void test13(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 13 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 4, 13, method_title);
}

void test14(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 14 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 5, 14, method_title);
}

void test15(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 15 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 1, 6, 15, method_title);
}

void test16(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 16 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 2, 16, method_title);
}

void test17(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 17 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 5, 17, method_title);
}

void test18(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 18 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 10, 18, method_title);
}

void test19(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 19 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 100, 19, method_title);
}

void test20(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 20 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 2, 1000, 20, method_title);
}

void test21(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 21 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 2, 21, method_title);
}

void test22(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 22 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 5, 22, method_title);
}

void test23(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 23 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 10, 23, method_title);
}

void test24(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 24 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 100, 24, method_title);
}

void test25(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 25 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 3, 1000, 25, method_title);
}

void test26(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 26 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 2, 26, method_title);
}

void test27(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 27 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 6, 27, method_title);
}

void test28(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 28 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 10, 28, method_title);
}

void test29(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 29 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 100, 29, method_title);
}

void test30(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 30 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 4, 1000, 30, method_title);
}

void test31(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 31 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 2, 31, method_title);
}

void test32(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 32 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 6, 32, method_title);
}

void test33(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 33 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 10, 33, method_title);
}

void test34(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 34 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 100, 34, method_title);
}

void test35(Method method, const std::string& method_title) {
    fout_txt << "------------------ Тест 35 ------------------\n";
    test_Nesterov_n_Rosenbrock(method, 5, 1000, 35, method_title);
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
    test12 (method, method_title);    std::cout << ", 12";  std::cout.flush();  fout_txt.flush();
    test13 (method, method_title);    std::cout << ", 13";  std::cout.flush();  fout_txt.flush();
    test14 (method, method_title);    std::cout << ", 14";  std::cout.flush();  fout_txt.flush();
    test15 (method, method_title);    std::cout << ", 15";  std::cout.flush();  fout_txt.flush();
    test16 (method, method_title);    std::cout << ", 16";  std::cout.flush();  fout_txt.flush();
    test17 (method, method_title);    std::cout << ", 17";  std::cout.flush();  fout_txt.flush();
    test18 (method, method_title);    std::cout << ", 18";  std::cout.flush();  fout_txt.flush();
    test19 (method, method_title);    std::cout << ", 19";  std::cout.flush();  fout_txt.flush();
    // test20 (method);       std::cout << ", 20";  std::cout.flush();  fout_txt.flush();
    test21 (method, method_title);    std::cout << ", 21";  std::cout.flush();  fout_txt.flush();
    test22 (method, method_title);    std::cout << ", 22";  std::cout.flush();  fout_txt.flush();
    test23 (method, method_title);    std::cout << ", 23";  std::cout.flush();  fout_txt.flush();
    test24 (method, method_title);    std::cout << ", 24";  std::cout.flush();  fout_txt.flush();
    // test25 (method);       std::cout << ", 25";  std::cout.flush();  fout_txt.flush();
    test26 (method, method_title);    std::cout << ", 26";  std::cout.flush();  fout_txt.flush();
    test27 (method, method_title);    std::cout << ", 27";  std::cout.flush();  fout_txt.flush();
    test28 (method, method_title);    std::cout << ", 28";  std::cout.flush();  fout_txt.flush();
    test29 (method, method_title);    std::cout << ", 29";  std::cout.flush();  fout_txt.flush();
    // test30 (method);       std::cout << ", 30";  std::cout.flush();  fout_txt.flush();
    test31 (method, method_title);    std::cout << ", 31";  std::cout.flush();  fout_txt.flush();
    test32 (method, method_title);    std::cout << ", 32";  std::cout.flush();  fout_txt.flush();
    test33 (method, method_title);    std::cout << ", 33";  std::cout.flush();  fout_txt.flush();
    test34 (method, method_title);    std::cout << ", 34";  std::cout.flush();  fout_txt.flush();
    // test35 (method);       std::cout << ", 35";  std::cout.flush();  fout_txt.flush();
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
