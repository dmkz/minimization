#include <iostream>
#include "vector.h"
#include "matrix.h"

int main() {	
	Vector v(5, 1), w(5, 2);
	std::cout << v+w << std::endl << std::endl;

	Matrix m = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	
	std::cout << m << std::endl;
	return 0;
}