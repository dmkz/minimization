#include "global_min.hpp"

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