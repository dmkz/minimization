#include <iostream>
#include "vector.h"

int main() {
	Vector v(5, 1);
	auto w = 2 * v;
	std::cout << w << std::endl;
	std::cout << v << std::endl;
	return 0;
}