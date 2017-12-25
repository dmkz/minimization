#include "sobolseqgenerator.h"
#include "math.hpp"

#include <iostream>
#include <iomanip>

int main() {
	SobolSeqGenerator net_gen;
    net_gen.Init(10, 3, "new-joe-kuo-6.21201.txt");
	for (int i = 0; i < 10; ++i) {
		std::cout << Vector(net_gen.GeneratePoint().coordinate) << "\n";
	}
	return 0;
}