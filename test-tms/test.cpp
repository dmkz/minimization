/*
 * Авторы: Бураханова А., Золкин А., Казарян М., Хохлов С., Малоенко С.
 */

#include "UniquenessTest.hpp"
#include "TMSParametersTest.hpp"
#include "NiederreiterBaseTwo.hpp"
#include "sobolseqgenerator.hpp"

int main() {
	{
		SobolSeqGenerator* generator_ptr = new SobolSeqGenerator();
		uint32_t dim = 22;
		uint32_t point_num = 1ULL << 19;
		UniquenessTest test;
		test.Init(dim, point_num, generator_ptr);
		test.Run();
		delete generator_ptr;
	}
	
	{
		NiederreiterBaseTwo* generator_ptr = new NiederreiterBaseTwo();
		uint32_t dim = 22;
		uint32_t point_num = 1ULL << 19;
		UniquenessTest test;
		test.Init(dim, point_num, generator_ptr);
		test.Run();
		delete generator_ptr;
	}

	{
		SobolSeqGenerator* generator_ptr = new SobolSeqGenerator();
		uint32_t s_start = 1;
		uint32_t s_finish = 8;
		uint32_t m_start = 0;
		uint32_t m_finish = 8;
		
		std::vector<std::vector<int> > results(s_finish-s_start + 1, std::vector<int>(m_finish-m_start + 1));
		
		uint32_t b = 2;
		
		for(uint32_t s = s_start; s <= s_finish; s++) 
		{
			for(uint32_t m = m_start; m <= m_finish; m++) 
			{
				uint32_t t = m + 1; 
				
				TMSParametersTest test;
				std::cout << "\nStarted: " << "m = " << m << ", s = " << s << ". Timestamp: "<< test.CurrentDateTime();
				do  
				{
					t--;
					test.Init(t, m, s, b, generator_ptr);
					int test_result = test.Run();
					if(test_result == 1)
					{
						std::cout << "\nFinished: " << "t = " <<  results[s-s_start][m-m_start] << ", time = " << test.CurrentDateTime();
						break;
					} else {
						results[s-s_start][m-m_start] = t;
						if(t == 0)
						{
							std::cout << "\nFinished: " << "t = " <<  results[s-s_start][m-m_start] << ", time = " << test.CurrentDateTime();
						}
					}
				} while(t != 0);
			}
		}
		
		std::cout << "\nMinimal defect test for SobolSeqGenerator is finished!";

		for(uint32_t m = m_start; m <= m_finish; m++) 
		{
			std::cout << "\n";
			for(uint32_t s = s_start; s <= s_finish; s++) 
			{
				std::cout << results[s-s_start][m-m_start] << " ";
			}
		}
		delete generator_ptr;
	}
	
	{
		NiederreiterBaseTwo* generator_ptr = new NiederreiterBaseTwo();
		uint32_t s_start = 1;
		uint32_t s_finish = 8;
		uint32_t m_start = 0;
		uint32_t m_finish = 8;
		
		std::vector<std::vector<int> > results(s_finish-s_start + 1, std::vector<int>(m_finish-m_start + 1));
		
		uint32_t b = 2;
		
		for(uint32_t s = s_start; s <= s_finish; s++) 
		{
			for(uint32_t m = m_start; m <= m_finish; m++) 
			{
				uint32_t t = m + 1; 
				
				TMSParametersTest test;
				std::cout << "\nStarted: " << "m = " << m << ", s = " << s << ". Timestamp: "<< test.CurrentDateTime();
				do  
				{
					t--;
					test.Init(t, m, s, b, generator_ptr);
					int test_result = test.Run();
					if(test_result == 1)
					{
						std::cout << "\nFinished: " << "t = " <<  results[s-s_start][m-m_start] << ", time = " << test.CurrentDateTime();
						break;
					} else {
						results[s-s_start][m-m_start] = t;
						if(t == 0)
						{
							std::cout << "\nFinished: " << "t = " <<  results[s-s_start][m-m_start] << ", time = " << test.CurrentDateTime();
						}
					}
				} while(t != 0);
			}
		}
		
		std::cout << "\nMinimal defect test for NiederreiterBaseTwo is finished!";

		for(uint32_t m = m_start; m <= m_finish; m++) 
		{
			std::cout << "\n";
			for(uint32_t s = s_start; s <= s_finish; s++) 
			{
				std::cout << results[s-s_start][m-m_start] << " ";
			}
		}
		delete generator_ptr;
	}
    return 0;
}