#include "sobolseqgenerator.h"
#include "math.hpp"
#include "nesterov.hpp"
#include "hessian_free.hpp"
#include "bfgs.hpp"

#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <algorithm>

std::vector<std::pair<ld, Vector>>
find_absmin(Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max) {
	assert(dim > 0u && dim == min.size() && dim == max.size());
	assert(nBestPoints <= nAllPoints && nBestPoints > 0u);
	for (uint32_t i = 0; i < dim; ++i) {
		assert(min[i] <= max[i]);
	}
	
	SobolSeqGenerator net;
	net.Init(nAllPoints, dim, "new-joe-kuo-6.21201.txt");

	// Формирование списка лучших кандидатов
	std::set<std::pair<ld, Vector>> candidates;
	for (uint32_t i = 0; i < nAllPoints; ++i) {
		// Получение текущего узла сетки:
		const auto net_point = net.GeneratePoint().coordinate;
		// Преобразование узла к точке в ограниченной min и max области:
		Vector curr(dim);
		for (uint32_t i = 0; i < dim; ++i) {
			curr[i] = net_point[i] * (max[i] - min[i]) + min[i];
		}
		candidates.insert({f(curr), curr});
		if (candidates.size() > nBestPoints) {
			candidates.erase(std::prev(candidates.end()));
		}
	}
	
	std::vector<std::pair<ld, Vector>> answer(candidates.size());
	
	auto it = candidates.begin();
	for (uint32_t i = 0; i < answer.size(); ++i) {
		auto x1 = bfgs(f, it->second, 10).first;
		auto x2 = hessian_free(f, it->second, 10).first;
		auto x3 = nesterov(f, it->second, 10).first;
		answer[i] = std::min({
			*it,
			std::make_pair(f(x1), x1),
			std::make_pair(f(x2), x2), 
			std::make_pair(f(x3), x3)
		});
		it++;
	}
	
	std::sort(answer.begin(), answer.end());
	
	return answer;
}

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

void test(std::string title, Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max) {
	std::cout << "----------------------------------------- " << title << " -----------------------------------------\n";
	for (auto & rec : find_absmin(f, dim, nBestPoints, nAllPoints, min, max)) {
		std::cout << "\tf_min = " << std::fixed << std::setprecision(6) << std::setw(12) << rec.first << ", point = {" << rec.second << "}\n";
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

int main() {
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
	
	return 0;
}