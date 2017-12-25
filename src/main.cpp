#include "global_min.hpp"

#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <algorithm>

#define N_max 100

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
    for (unsigned int i=0; i < v.size(); ++i) {
        fun += std::pow(v[i], 2);
    }
    return fun;
}

ld f17(const Vector &v) {
    ld fun = 10 * v.size();
    for (unsigned int i=0; i < v.size(); ++i) {
        fun += std::pow(v[i] - 10 * std::cos(2 * M_PI * v[i]), 2);
    }
    return fun;
}

ld f18(const Vector &v) {
    return std::pow(std::pow(v[0], 2) + v[1] - 11, 2) + std::pow(std::pow(v[1], 2) + v[0] - 7, 2);
}

ld f19(const Vector &v) {
    const int b[] = {8, 18, 44, 144};
    ld fun = 0;
    ld g_fun;
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

ld f20(const Vector &v) {
    return std::pow(v[1] - 5.1 / (4 * M_PI_2) * std::pow(v[0], 2) + 5 * v[0] / M_PI - 6, 2) +
           10 * (1 - 1 / (8 * M_PI)) * std::cos(v[0]) + 10;
}

ld f21(const Vector &v) {
    return std::sin(v[0] + v[1]) + std::pow(v[0] - v[1], 2) - 1.5 * v[0] + 2.5 * v[1] + 1;
}

ld f22(const Vector &v) {
    return 0.26 * (std::pow(v[0], 2) + std::pow(v[1], 2)) - 0.48 * v[0] * v[1];
}

void test(std::string title, Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max) {
	std::cout << "----------------------------------------- " << title << " -----------------------------------------" << std::endl;
	for (auto & rec : find_absmin(f, dim, nBestPoints, nAllPoints, min, max)) {
		std::cout << "\tf_min = " << std::fixed << std::setprecision(6) << std::setw(12) << rec.first << ", point = {" << rec.second << '}' << std::endl;
	}
}

void test1() {
	test("Test  1", f1, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test2() {
	test("Test  2", f2, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test3() {
	test("Test  3", f3, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test4() {
	test("Test  4", f4, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test5() {
	test("Test  5", f5, 3, 5, 32, Vector(3, -5), Vector(3, 5));
}

void test6() {
	test("Test  6", f6, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test7() {
	test("Test  7", f7, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test8() {
	test("Test  8", f8, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test9() {
	test("Test  9", f9, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test10() {
	test("Test 10", f10, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test11() {
    test("Test 11", f18, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test12() {
    test("Test 12", f20, 2, 5, 32, {-5, 0}, {10, 15});
}

void test13() {
    test("Test 13", f21, 2, 5, 32, {-1.5, -3}, {4, 4});
}

void test14() {
    test("Test 14", f22, 2, 5, 32, Vector(2, -10), Vector(2, 10));
}

void test15() {
    test("Test 15", f11, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test16() {
    test("Test 16", f12, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test17() {
    test("Test 17", f13, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test18() {
    test("Test 18", f14, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test19() {
    test("Test 19", f15, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test20() {
    test("Test 20", f16, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test21() {
    test("Test 21", f17, 2, 5, 32, Vector(2, -5), Vector(2, 5));
}

void test22() {
    test("Test 22", f19, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test23() {
    test("Test 23", f11, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test24() {
    test("Test 24", f12, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test25() {
    test("Test 25", f13, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test26() {
    test("Test 26", f14, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test27() {
    test("Test 27", f15, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test28() {
    test("Test 28", f16, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test29() {
    test("Test 29", f17, 4, 5, 128, Vector(4, -5), Vector(4, 5));
}

void test30() {
    test("Test 30", f11, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test31() {
    test("Test 31", f12, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test32() {
    test("Test 32", f13, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test33() {
    test("Test 33", f14, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test34() {
    test("Test 34", f15, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test35() {
    test("Test 35", f16, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test36() {
    test("Test 36", f17, 10, 5, 2048, Vector(10, -5), Vector(10, 5));
}

void test37() {
    test("Test 37", f11, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test38() {
    test("Test 38", f12, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test39() {
    test("Test 39", f13, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test40() {
    test("Test 40", f14, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test41() {
    test("Test 41", f15, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test42() {
    test("Test 42", f16, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test43() {
    test("Test 43", f17, 20, 5, 1024u * 1024u, Vector(20, -5), Vector(20, 5));
}

void test44() {
    test("Test 44", f11, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test45() {
    test("Test 45", f12, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test46() {
    test("Test 46", f13, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test47() {
    test("Test 47", f14, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test48() {
    test("Test 48", f15, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test49() {
    test("Test 49", f16, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

void test50() {
    test("Test 50", f17, N_max, 5, 1024u * 1024u, Vector(N_max, -5), Vector(N_max, 5));
}

int main() {
    freopen("test_result.txt", "wt", stdout);
	test1 ();
	test2 ();
	test3 ();
	test4 ();
	test5 ();
	test6 ();
	test7 ();
	test8 ();
	test9 ();
	test10();
	test21();
	test12();
	test13();
	test14();
	test15();
	test16();
	test17();
	test18();
	test19();
	test20();
	test21();
	test22();
	test23();
	test24();
	test25();
	test26();
	test27();
	test28();
	test29();
	test30();
	test31();
	test32();
	test33();
	test34();
	test35();
	test36();
	test37();
	test38();
	test39();
	test40();
	test41();
	test42();
	test43();
	test44();
	test45();
	test46();
	test47();
	test48();
	test49();
	test50();
	return 0;
}
