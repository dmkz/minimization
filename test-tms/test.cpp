/*
 * Авторы: Бураханова А., Золкин А., Казарян М., Хохлов С., Малоенко С.
 */

#include "UniquenessTest.hpp"
#include "sobolseqgenerator.hpp"

int main() {
	SobolSeqGenerator* generator_ptr = new SobolSeqGenerator();
	uint32_t dim = 1024;
	uint32_t point_num = 1ULL << 20;
	UniquenessTest test;
	test.Init(dim, point_num, generator_ptr);
	test.Run();
    return 0;
}