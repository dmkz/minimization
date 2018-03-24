#include "global_min.hpp"

/*
    Тестирование глобальной оптимизации.
    Автор: Юрий Кондратов.
*/

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>

#define N_max 100
#define M_E 2.71828182845904523536
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923

// Тестовые функции:
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

Real f11(const Vector& v) {
    Real fun = std::pow(v[0] - 1, 2) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::pow(v[i] - 2 * std::pow(v[i-1], 2) + 1, 2);
    }
    return fun;
}

Real f12(const Vector &v) {
    Real fun = std::pow(v[0] - 1, 2) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::abs(v[i] - 2 * std::pow(v[i-1], 2) + 1);
    }
    return fun;
}

Real f13(const Vector &v) {
    Real fun = std::abs(v[0] - 1) / 4;
    for (unsigned int i=1; i < v.size(); ++i) {
        fun += std::abs(v[i] - 2 * std::pow(v[i-1], 2) + 1);
    }
    return fun;
}

Real f14(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=1; i < v.size(); i += 2) {
        fun += 100 * std::pow(std::pow(v[i-1], 2) - v[i], 2) + std::pow(v[i-1] - 1, 2);
    }
    return fun;
}

Real f15(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=1; i < v.size(); i += 2) {
        fun += 10 * std::abs(std::pow(v[i-1], 2) - v[i]) + std::abs(v[i-1] - 1);
    }
    return fun;
}

Real f16(const Vector &v) {
    Real fun = 0;
    for (unsigned int i=0; i < v.size(); ++i) {
        fun += std::pow(v[i], 2);
    }
    return fun;
}

Real f17(const Vector &v) {
    Real fun = 10 * v.size();
    for (unsigned int i=0; i < v.size(); ++i) {
        fun += std::pow(v[i] - 10 * std::cos(2 * M_PI * v[i]), 2);
    }
    return fun;
}

Real f18(const Vector &v) {
    return std::pow(std::pow(v[0], 2) + v[1] - 11, 2) + std::pow(std::pow(v[1], 2) + v[0] - 7, 2);
}

Real f19(const Vector &v) {
    const int b[] = {8, 18, 44, 144};
    Real fun = 0;
    Real g_fun;
    for (int i = 0; i < 4; ++i) {
        g_fun = 0;
        for (int j = 0; j < 4; ++j) {
            g_fun += 1;
            g_fun *= v[i];
        }
        g_fun -= b[i];
        fun += std::pow(g_fun, 2);
    }
    return fun;
}

Real f20(const Vector &v) {
    return std::pow(v[1] - 5.1 / (4 * M_PI_2) * std::pow(v[0], 2) + 5 * v[0] / M_PI - 6, 2) +
           10 * (1 - 1 / (8 * M_PI)) * std::cos(v[0]) + 10;
}

Real f21(const Vector &v) {
    return std::sin(v[0] + v[1]) + std::pow(v[0] - v[1], 2) - 1.5 * v[0] + 2.5 * v[1] + 1;
}

Real f22(const Vector &v) {
    return 0.26 * (std::pow(v[0], 2) + std::pow(v[1], 2)) - 0.48 * v[0] * v[1];
}

void test(std::string title, Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max, std::ofstream& fout) {
    std::cout << "-- " << title << std::endl; 
    std::cout.flush();
	fout << "----------------------------------------- " << title << " -----------------------------------------" << std::endl;
	for (auto & rec : find_absmin(f, default_stop_condition, dim, nBestPoints, nAllPoints, min, max)) {
		fout << "\tf_min = " << std::fixed << std::setprecision(6) << std::setw(12) << rec.first << ", point = {" << rec.second << '}' << std::endl;
	}
    fout.flush();
}

void test1(std::ofstream& fout) {
	test("Test  1", f1, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test2(std::ofstream& fout) {
	test("Test  2", f2, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test3(std::ofstream& fout) {
	test("Test  3", f3, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test4(std::ofstream& fout) {
	test("Test  4", f4, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test5(std::ofstream& fout) {
	test("Test  5", f5, 3, 5, 32, Vector(3, -5), Vector(3, 5), fout);
}

void test6(std::ofstream& fout) {
	test("Test  6", f6, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test7(std::ofstream& fout) {
	test("Test  7", f7, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test8(std::ofstream& fout) {
	test("Test  8", f8, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test9(std::ofstream& fout) {
	test("Test  9", f9, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test10(std::ofstream& fout) {
	test("Test 10", f10, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test11(std::ofstream& fout) {
    test("Test 11", f18, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test12(std::ofstream& fout) {
    test("Test 12", f20, 2, 5, 32, {-5, 0}, {10, 15}, fout);
}

void test13(std::ofstream& fout) {
    test("Test 13", f21, 2, 10, 32, {-1.5, -3}, {4, 4}, fout);
}

void test14(std::ofstream& fout) {
    test("Test 14", f22, 2, 5, 32, Vector(2, -10), Vector(2, 10), fout);
}

void test15(std::ofstream& fout) {
    test("Test 15", f11, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test16(std::ofstream& fout) {
    test("Test 16", f12, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test17(std::ofstream& fout) {
    test("Test 17", f13, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test18(std::ofstream& fout) {
    test("Test 18", f14, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test19(std::ofstream& fout) {
    test("Test 19", f15, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test20(std::ofstream& fout) {
    test("Test 20", f16, 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test21(std::ofstream& fout) {
    test("Test 21", f17, 2, 10, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test22(std::ofstream& fout) {
    test("Test 22", f19, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test23(std::ofstream& fout) {
    test("Test 23", f11, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test24(std::ofstream& fout) {
    test("Test 24", f12, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test25(std::ofstream& fout) {
    test("Test 25", f13, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test26(std::ofstream& fout) {
    test("Test 26", f14, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test27(std::ofstream& fout) {
    test("Test 27", f15, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test28(std::ofstream& fout) {
    test("Test 28", f16, 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test29(std::ofstream& fout) {
    test("Test 29", f17, 4, 20, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test30(std::ofstream& fout) {
    test("Test 30", f11, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test31(std::ofstream& fout) {
    test("Test 31", f12, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test32(std::ofstream& fout) {
    test("Test 32", f13, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test33(std::ofstream& fout) {
    test("Test 33", f14, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test34(std::ofstream& fout) {
    test("Test 34", f15, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test35(std::ofstream& fout) {
    test("Test 35", f16, 10, 5, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test36(std::ofstream& fout) {
    test("Test 36", f17, 10, 50, 2048, Vector(10, -5), Vector(10, 5), fout);
}

void test37(std::ofstream& fout) {
    test("Test 37", f11, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test38(std::ofstream& fout) {
    test("Test 38", f12, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test39(std::ofstream& fout) {
    test("Test 39", f13, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test40(std::ofstream& fout) {
    test("Test 40", f14, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test41(std::ofstream& fout) {
    test("Test 41", f15, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test42(std::ofstream& fout) {
    test("Test 42", f16, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test43(std::ofstream& fout) {
    test("Test 43", f17, 20, 50, 1024u * 1024u, Vector(20, -5), Vector(20, 5), fout);
}

void test44(std::ofstream& fout) {
    test("Test 44", f11, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test45(std::ofstream& fout) {
    test("Test 45", f12, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test46(std::ofstream& fout) {
    test("Test 46", f13, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test47(std::ofstream& fout) {
    test("Test 47", f14, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test48(std::ofstream& fout) {
    test("Test 48", f15, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test49(std::ofstream& fout) {
    test("Test 49", f16, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

void test50(std::ofstream& fout) {
    test("Test 50", f17, N_max, 50, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5), fout);
}

int main() {
	/*
    auto start_d = std::chrono::system_clock::now();
    auto start_t = std::chrono::system_clock::to_time_t(start_d);
	std::string s = std::ctime(&start_t);
	s.pop_back();
	s = "Test " + s + ".txt";
	*/
    std::ofstream fout;
    fout.open("test_results.txt");
	
	std::cout << "-- Number of Cores = " << std::thread::hardware_concurrency() << std::endl;
	std::cout << "-- Start testing... The results will be written to a file test_results.txt" << std::endl;
	test1 (fout);
	test2 (fout);
	test3 (fout);
	test4 (fout);
	test5 (fout);
	test6 (fout);
	test7 (fout);
	test8 (fout);
	test9 (fout);
	test10(fout);
	test11(fout);
	test12(fout);
	test13(fout);
	test14(fout);
	test15(fout);
	test16(fout);
	test17(fout);
	test18(fout);
	test19(fout);
	test20(fout);
	test21(fout);
	test22(fout);
	test23(fout);
	test24(fout);
	test25(fout);
	test26(fout);
	test27(fout);
	test28(fout);
	test29(fout);
	test30(fout); 
	test31(fout);
	test32(fout);
	test33(fout);
	test34(fout);
	test35(fout);
	test36(fout); /*
	test37(fout);
	test38(fout);
	test39(fout);
	test40(fout);
	test41(fout);
	test42(fout); 
	test43(fout);
	test44(fout);
	test45(fout);
	test46(fout);
	test47(fout);
	test48(fout);
	test49(fout);
	test50(fout); */
    fout.close();
    std::cout << "-- Finish testing... The results are written to a file test_results.txt" << std::endl;
	return 0;
}
