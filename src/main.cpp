#include "sobolseqgenerator.h"

#include <iostream>
#include <iomanip>

int main() {
	SobolSeqGenerator net_gen;
    net_gen.Init(10, 3, "new-joe-kuo-6.21201.txt");
	for (int i = 0; i < 10; ++i) {
		for (auto & it : net_gen.GeneratePoint().coordinate) {
			std::cout << std::fixed << std::setw(8) << std::setprecision(6) << it << " ";
		}
		std::cout << "\n";
	}
	return 0;
}