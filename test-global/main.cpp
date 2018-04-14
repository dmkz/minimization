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
	return std::sin(v[0] + v[1]) + std::pow(v[0] - v[1], 2) + 1.5 * v[0]*v[0] + 2.5 * v[1]*v[1] + 1;
}

Real f15(const Vector &v) {
    return 0.26 * (std::pow(v[0], 2) + std::pow(v[1], 2)) - 0.48 * v[0] * v[1];
}

Real f16(const Vector &v) {
    //Keane
    return -(std::pow(std::sin(v[0] - v[1]), 2) * std::pow(std::sin(v[0] + v[1]), 2)) / std::sqrt(v[0] * v[0] + v[1] *v[1]);
}

Real f17(const Vector &v) {
    //SchmidtVetters
    return 1 / (1 + std::pow(v[0] - v[1], 2)) + std::sin((M_PI * v[1] + v[2]) / 2) +
           std::exp(std::pow((v[0] + v[1]) / v[1] - 2, 2));
}

Real f18(const Vector &v) {
    //Zettl
    return v[0] / 4 + std::pow(v[0] * v[0] - 2 * v[0] + v[1] * v[1], 2);
}

Real f19(const Vector &v) {
    // Beale
    return std::pow(v[0] * v[1] - v[0] + 1.5, 2) +
           std::pow(v[0] * v[1] * v[1] - v[0] + 2.25, 2) +
           std::pow(v[0] * std::pow(v[2], 3) - v[0] + 2.625, 2);
}

Real f20(const Vector &v) {
    // Colville
    return std::pow(v[0] - 1, 2) + 100 * std::pow(v[0] * v[0] - v[1], 2) +
           10.1 * std::pow(v[2] - 1, 2) + std::pow(v[2] - 1, 2) +
           90 * std::pow(v[2] * v[2] - v[3], 2) + 10.1 * std::pow(v[3] - 1, 2) +
           19.8 * (v[3] - 1) * (v[1] - 1);
}

Real f21(const Vector &v) {
    // II
    return std::pow(v[0] * v[0] - v[1], 2) + std::pow(1 - v[0], 2);
}

Real f22(const Vector &v) {
    // III
    return std::pow(v[0] * v[0] - v[1], 2) + 100 * std::pow(1 - v[0], 2);
}
Real f23(const Vector &v) {
    // IV
    return 100 * std::pow(v[1] - std::pow(v[0], 3), 2) + std::pow(1 - v[0], 2);
}
Real f24(const Vector &v) {
    //  VII
    return std::pow(v[0] + 10 * v[1],2) + 5 * std::pow(v[3] - v[2], 2) +
           std::pow(v[1] - 2 * v[2], 4) + 10 * std::pow(v[0] - v[3], 4);
}
Real f25(const Vector &v) {
    // VIII f_min(0,1,1,1+-pi*n)=0
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
    // Bohachevsky
    Real fun = 0;
    for (unsigned int i=0; i < v.size() - 1; i++) {
        fun += v[i] * v[i] + 2 * v[i+1] * v[i+1] -0.3 * std::cos(3 * M_PI * v[i]) -
               0.4 * std::cos(4 * M_PI * v[i+1]) + 0.7;
    }
    return fun;
}

Real f30(const Vector &v) {
    // Cigar problem
    // Медленный спуск
    Real fun = 0;
    for (unsigned int i=1; i < v.size(); i++) {
        fun += v[i] * v[i];
    }
    fun *= 1000000;
    fun += v[0] * v[0];
    return fun;
}


void test(std::string title, Function f, std::string description_f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max, std::ofstream& fout) {
    std::cout << "-- " << title << std::endl;
    std::cout.flush();
	fout << "\n";
	fout << "----------------------------------------- " << title << " -----------------------------------------" << std::endl;
    fout << description_f << "\n" << std::endl;
	for (auto & rec : find_absmin(f, default_stop_condition, dim, nBestPoints, nAllPoints, min, max)) {
		fout << "\tf_min = " << std::fixed << std::setprecision(6) << std::setw(12) << rec.first << ", point = {" << rec.second << '}' << std::endl;
	}
    fout.flush();
}

void test1(std::ofstream& fout) {
	test("Test  1, dim 02", f1, "\nГладкая функция:\n\tf(x,y) = 1+x+y-xy+x^2+y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test2(std::ofstream& fout) {
	test("Test  2, dim 02", f2, "\nГладкая функция:\n\tf(x,y) = 1+7x+5y+0.5xy+3x^2+y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test3(std::ofstream& fout) {
	test("Test  3, dim 02", f3, "\nГладкая функция:\n\tf(x,y) = 100+7x+5y-10xy+3x^2+10y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test4(std::ofstream& fout) {
	test("Test  4, dim 02", f4, "\nГладкая функция:\n\tf(x,y) = 100+7x+5y-10.95xy+3x^2+10y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test5(std::ofstream& fout) {
	test("Test  5, dim 03", f5, "\nГладкая функция:\n\tf(x,y,z) = 1+x+y+z+xy+xz+yz+x^2+y^2+z^2", 3, 5, 64, Vector(3, -5), Vector(3, 5), fout);
}

void test6(std::ofstream& fout) {
	test("Test  6, dim 02", f6, "\nГладкая функция:\n\tf(x,y) = 10x^4+15y^4+15xy", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test7(std::ofstream& fout) {
	test("Test  7, dim 02", f7, "\nГладкая функция:\n\tf(x,y) = 10x^6+15y^6-20(y*x^3+x*y^3)", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test8(std::ofstream& fout) {
	test("Test  8, dim 02", f8, "\nГладкая функция:\n\tf(x,y) = x^6+y^6-2(y*x^3+x*y^3)+x^2+y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test9(std::ofstream& fout) {
	test("Test  9, dim 02", f9, "\nГладкая функция:\n\tf(x,y) = x^6+y^6-3(y*x^3+x*y^3)+x^2+y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test10(std::ofstream& fout) {
	test("Test 10, dim 02", f10, "\nГладкая функция:\n\tf(x,y) = x^6+y^6-2(y*x^3+x*y^3)+x^4+y^4-x^2-y^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test11(std::ofstream& fout) {
    test("Test 11, dim 02", f11, "\nГладкая функция:\n\tf(x,y) = (x^2+y-11)^2+(y^2+x-7)^2", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test12(std::ofstream& fout) {
    test("Test 12, dim 04", f12, "\nГладкая функция:\n\tf(x1,x2,x3,x4) = (x1+x2+x3+x4-8)^2+(x1^2+x2^2+x3^2+x4^2-18)^2+(x1^3+x2^3+x3^3+x4^3-44)^2+(x1^4+x2^4+x3^4+x4^4-114)^2", 4, 5, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test13(std::ofstream& fout) {
    test("Test 13, dim 02", f13, "\nГладкая функция:\n\tf(x,y) = (y-5.1/(4*M_PI_2)*x^2+5/M_PI*x-6)^2+10(1-1/(8*M_PI))cos(x)+10", 2, 10, 32, Vector(2, -10), Vector(2, 10), fout);
}

void test14(std::ofstream& fout) {
    test("Test 14, dim 02", f14, "\nГладкая функция МакКормика:\n\tf(x,y) = sin(x+y)+(x-y)^2+1.5x^2+2.5y^2+1", 2, 5, 32, Vector(2, -10), Vector(2, 10), fout);
}

void test15(std::ofstream& fout) {
    test("Test 15, dim 02", f15, "\nГладкая функция:\n\tf(x,y) = 0.26(x^2+y^2)-0.48xy", 2, 5, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test16(std::ofstream& fout) {
    test("Test 16, dim 02", f16, "\nГладкая функция Keane:\n\tf(x,y) = -((sin(x-y))^2*(sin(x+y))^2)/(x^2+y^2)^(1/2)", 2, 5, 32, Vector(2, 0), Vector(2, 10), fout);
}

void test17(std::ofstream& fout) {
    test("Test 17, dim 03",  f17, "\nГладкая функция SchmidtVetters:\n\tf(x,y,z) = 1/(1+(x-y)^2)+sin((M_PI*y+z)/2)+M_E^(((x+y)/y-2)^2)", 3, 5, 64, Vector(3, 0), Vector(3, 10), fout);
}

void test18(std::ofstream& fout) {
    test("Test 18, dim 02", f18, "\nГладкая функция Zettl:\n\tf(x,y) = x/4+(x^2-2x+y^2)^2", 2, 5, 32, Vector(2, -1), Vector(2, 5), fout);
}

void test19(std::ofstream& fout) {
    test("Test 19, dim 02", f19, "\nГладкая функция Beale:\n\tf(x,y) = (1.5-x+xy)^2+(2.25-x+xy^2)^2+(2.625-x+xy^3)^2", 2, 5, 32, Vector(2, -10), Vector(2, 10), fout);
}

void test20(std::ofstream& fout) {
    test("Test 20, dim 04", f20, "\nГладкая функция Colville:\n\tf(x,y,z,t) = 100(x-y^2)^2+(1-x)^2+90(t-z^2)^2+(1-z)^2+10.1((y-1)^2+(t-1)^2)+19.8(y-1)(t-1)", 4, 5, 128, Vector(4, -10), Vector(4, 10), fout);
}

void test21(std::ofstream& fout) {
    test("Test 21, dim 02", f21, "\nГладкая функция:\n\tf(x,y) = (x^2-y)^2+(1-x)^2",2, 10, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test22(std::ofstream& fout) {
    test("Test 22, dim 02", f22, "\nГладкая функция:\n\tf(x,y) = (x^2-y)^2+100(1-x)^2", 2, 10, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test23(std::ofstream& fout) {
    test("Test 23, dim 02", f23, "\nГладкая функция:\n\tf(x,y) = 100(y-x^3)^2+(1-x)^2", 2, 10, 32, Vector(2, -5), Vector(2, 5), fout);
}

void test24(std::ofstream& fout) {
    test("Test 24, dim 04", f24, "\nГладкая функция:\n\tf(x,y,z,t) = (x+10y)^2+5(t-z)^2+(y-2z)^4+10(x-t)^4", 4, 10, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test25(std::ofstream& fout) {
    test("Test 25, dim 04", f25, "\nГладкая функция:\n\tf(x,y,z,t) = (x^2-y)^4+100(y-z)^6+(tan(z-t))^4+x^8+(t^3-1)^2)", 4, 10, 128, Vector(4, -5), Vector(4, 5), fout);
}

void test26_2(std::ofstream& fout) {
    test("Test 26, dim 02", f26, "\nГладкая функция:\n\tf(x1,..,xn) = sum_(i=1)^(n-1)(x_i^2)", 2, 10, 32, Vector(2, -5), Vector(2, -5), fout);
}

void test26_4(std::ofstream& fout) {
    test("Test 26, dim 04", f26, "\nГладкая функция:\n\tf(x1,..,xn) = sum_(i=1)^(n-1)(x_i^2)", 4, 10, 128, Vector(4, -5), Vector(4, -5), fout);
}

void test26_8(std::ofstream& fout) {
    test("Test 26, dim 08", f26, "\nГладкая функция:\n\tf(x1,..,xn) = sum_(i=1)^(n-1)(x_i^2)", 8, 10, 128, Vector(8, -5), Vector(8, -5), fout);
}

void test26_12(std::ofstream& fout) {
    test("Test 26, dim 12", f26, "\nГладкая функция:\n\tf(x1,..,xn) = sum_(i=1)^(n-1)(x_i^2)", 12, 10, 128, Vector(12, -5), Vector(12, -5), fout);
}


void test27_2(std::ofstream& fout) {
    test("Test 27, dim 02", f27, "\nГладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", 2, 10, 32, Vector(2, -5), Vector(2, -5), fout);
}

void test27_4(std::ofstream& fout) {
    test("Test 27, dim 04", f27, "\nГладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", 4, 10, 128, Vector(4, -5), Vector(4, -5), fout);
}

void test27_8(std::ofstream& fout) {
    test("Test 27, dim 08", f27, "\nГладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", 8, 10, 128, Vector(8, -5), Vector(8, -5), fout);
}

void test27_12(std::ofstream& fout) {
    test("Test 27, dim 12", f27, "\nГладкая функция Нестерова Чебышева-Розенброка 1:\n\tf(x1, ..., xn) = (x1-1)^2/4+sum_(i=1)^(n-1)(x_(i+1)-2*x_i^2+1)^2", 12, 10, 128, Vector(12, -5), Vector(12, -5), fout);
}

void test28_2(std::ofstream& fout) {
    test("Test 28, dim 02", f28, "\nГладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", 2, 10, 32, Vector(2, -5), Vector(2, -5), fout);
}

void test28_4(std::ofstream& fout) {
    test("Test 28, dim 04", f28, "\nГладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", 4, 10, 128, Vector(4, -5), Vector(4, -5), fout);
}

void test28_8(std::ofstream& fout) {
    test("Test 28, dim 08", f28, "\nГладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", 8, 10, 128, Vector(8, -5), Vector(8, -5), fout);
}

void test28_12(std::ofstream& fout) {
    test("Test 28, dim 12", f28, "\nГладкая функция Розенброка:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (100(x_i^2-x_(i+1))^2+(x_i-1)^2)", 12, 10, 128, Vector(12, -5), Vector(12, -5), fout);
}

void test29_2(std::ofstream& fout) {
    test("Test 29, dim 02", f29, "\nГладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", 2, 10, 32, Vector(2, -5), Vector(2, -5), fout);
}

void test29_4(std::ofstream& fout) {
    test("Test 29, dim 04", f29, "\nГладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", 4, 10, 128, Vector(4, -5), Vector(4, -5), fout);
}

void test29_8(std::ofstream& fout) {
    test("Test 29, dim 08", f29, "\nГладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", 8, 10, 128, Vector(8, -5), Vector(8, -5), fout);
}

void test29_12(std::ofstream& fout) {
    test("Test 29, dim 12", f29, "\nГладкая функция Бохачевсокого:\n\tf(x1, ..., xn) = sum_(i=1)^(n-1) (x(i)^2+2x(i+1)^2-0.3*cos(3*pi*x(i)-0.4*cos(4*pi*x(i+1))+0.7)", 12, 10, 128, Vector(12, -5), Vector(12, -5), fout);
}


void test30_2(std::ofstream& fout) {
    test("Test 30, dim 02", f30, "\nГладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", 2, 10, 32, Vector(2, -5), Vector(2, -5), fout);
}

void test30_4(std::ofstream& fout) {
    test("Test 30, dim 04", f30, "\nГладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", 4, 10, 128, Vector(4, -5), Vector(4, -5), fout);
}

void test30_8(std::ofstream& fout) {
    test("Test 30, dim 08", f30, "\nГладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", 8, 10, 128, Vector(8, -5), Vector(8, -5), fout);
}

void test30_12(std::ofstream& fout) {
    test("Test 30, dim 12", f30, "\nГладкая функция:\n\tf(x1, ..., xn) = x_1^2+10^6*sum_(i=1)^(n-1) (x_i^2)", 12, 10, 128, Vector(12, -5), Vector(12, -5), fout);
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
	test26_2(fout);
	test26_4(fout);
	test26_8(fout);
	test26_12(fout);
	test27_2(fout);
	test27_4(fout);
	test27_8(fout);
	test27_12(fout);
	test28_2(fout);
	test28_4(fout);
	test28_8(fout);
	test28_12(fout);
	test29_2(fout);
	test29_4(fout);
	test29_8(fout);
	test29_12(fout);
	test30_2(fout);
	test30_4(fout);
	test30_8(fout);
	test30_12(fout);

    fout.close();
    std::cout << "-- Finish testing... The results are written to a file test_results.txt" << std::endl;
	return 0;
}
