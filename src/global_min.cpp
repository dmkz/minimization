#include "global_min.hpp"
#include <thread>
#include <mutex>

// Автор: Козырев Дмитрий, Бураханова Алена, Золкин Артем
std::vector<std::pair<Real, Vector>>
calc_f_with_threads(Function f, const std::vector<Vector> & inData) {
	// Создаем вектор под ответ:
	std::vector<std::pair<Real, Vector>> outData(inData.size());
	
	// Количество ядер:
	uint32_t nCores = std::max(1u, std::thread::hardware_concurrency());
	
	// Вектор из тредов:
	std::vector<std::thread> t;
	
	// Мьютексы на чтение и запись:
	std::mutex inRead, outWrite;
	
	uint32_t globalIndex = 0;
	
	// Создаем столько тредов, сколько ядер:
	for (uint32_t i = 0; i < nCores; i++) {
		t.push_back(std::thread([&] {
			while (1) {
				inRead.lock();
				uint32_t i = globalIndex;
				if (globalIndex >= inData.size()) {
					inRead.unlock();
					return;
				}
				
				auto x = inData[i];
				globalIndex++;
				inRead.unlock();

				auto res = f(x);

				outWrite.lock();
				outData[i] = {res, x};
				outWrite.unlock();
			}
			return;
		}));
	}
	
	// Присоединяем все треды к материнскому треду, гарантируя то, что все они будут выполнены
	for (uint32_t i = 0; i < nCores; ++i) {
		t[i].join();
	}
	
	assert(globalIndex == inData.size());
	
	return outData;
}

std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_local_mins_with_threads(Function f, const StopCondition& stop_condition, const std::vector<std::pair<Real, Vector>>& inData) {
	// Создаем вектор под ответ:
	std::vector<std::vector<std::pair<Real, GlobalTestData>>> outData(inData.size());
	
	// Количество ядер:
	uint32_t nCores = std::max(1u, std::thread::hardware_concurrency());
	
	// Вектор из тредов:
	std::vector<std::thread> t;
    
	// Мьютексы на чтение и запись:
	std::mutex inRead, outWrite;
	
	uint32_t globalIndex = 0;
	
	// Создаем столько тредов, сколько ядер:
	for (uint32_t thread_id = 0; thread_id < nCores; thread_id++) {
		t.push_back(std::thread([&] {
			while (1) {
				inRead.lock();
				uint32_t i = globalIndex;
				if (globalIndex >= inData.size()) {
					inRead.unlock();
					return;
				}
				
				auto it = inData[i];
				globalIndex++;
				inRead.unlock();
				
				// После чтения вызываем методы минимизации:
				auto iter_data = bfgs(f, it.second, stop_condition);
                auto x1 = iter_data.x_curr;
				auto f1 = iter_data.f_curr;
				auto sc1 = "Кол-во итераций: " + std::to_string(iter_data.iter_counter) 
					+ ", отклонение: " + std::to_string(std::abs(iter_data.f_curr - iter_data.f_prev));
                
                iter_data = hessian_free(f, it.second, stop_condition);
				auto x2 = iter_data.x_curr;
                auto f2 = iter_data.f_curr;
				auto sc2 = "Кол-во итераций: " + std::to_string(iter_data.iter_counter) 
					+ ", отклонение: " + std::to_string(std::abs(iter_data.f_curr - iter_data.f_prev));
                
                iter_data = nesterov(f, it.second, stop_condition);
				auto x3 = iter_data.x_curr;
                auto f3 = iter_data.f_curr;
				auto sc3 = "Кол-во итераций: " + std::to_string(iter_data.iter_counter) 
					+ ", отклонение: " + std::to_string(std::abs(iter_data.f_curr - iter_data.f_prev));
                
                iter_data = dfp(f, it.second, stop_condition);
                auto x4 = iter_data.x_curr;
                auto f4 = iter_data.f_curr;
				auto sc4 = "Кол-во итераций: " + std::to_string(iter_data.iter_counter) 
					+ ", отклонение: " + std::to_string(std::abs(iter_data.f_curr - iter_data.f_prev));
                
                iter_data = powell(f, it.second, stop_condition);
                auto x5 = iter_data.x_curr;
                auto f5 = iter_data.f_curr;
				auto sc5 = "Кол-во итераций: " + std::to_string(iter_data.iter_counter) 
					+ ", отклонение: " + std::to_string(std::abs(iter_data.f_curr - iter_data.f_prev));
                
				// Записываем ответ:
				outWrite.lock();
				outData[i].push_back(std::make_pair(it.first, GlobalTestData(it.second, std::string("Initial point(no method)"), std::string(""))));
				outData[i].push_back(std::make_pair(f1, GlobalTestData(x1, std::string("BFGS"), sc1)));
				outData[i].push_back(std::make_pair(f2, GlobalTestData(x2, std::string("Hessian Free"), sc2)));
				outData[i].push_back(std::make_pair(f3, GlobalTestData(x3, std::string("Nesterov"), sc3)));
				outData[i].push_back(std::make_pair(f4, GlobalTestData(x4, std::string("DFP"), sc4)));
				outData[i].push_back(std::make_pair(f5, GlobalTestData(x5, std::string("Powell"), sc5)));
				
				outWrite.unlock();
			}
			return;
		}));
	}
	
	// Присоединяем все треды к материнскому треду, гарантируя то, что все они будут выполнены
	for (uint32_t i = 0; i < nCores; ++i) {
		t[i].join();
	}
    
	assert(globalIndex == inData.size());
	
	return outData;
	
}

std::vector<std::vector<std::pair<Real, GlobalTestData>>>
find_absmin(Function f, const StopCondition& stop_condition, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max) {
	// Несколько проверок на входные данные:
	assert(dim > 0u && dim == min.size() && dim == max.size());
	assert(nBestPoints <= nAllPoints && nBestPoints > 0u);
	for (uint32_t i = 0; i < dim; ++i) {
		assert(min[i] <= max[i]);
	}
	
	// Объект-генератор сетки:
	NiederreiterBaseTwo net;
	net.Init(dim, nAllPoints);
	
	// ----- Первый этап: вычисление значений функции в узлах сетки с отбором точек-кандидатов -----
	
	// Лимит на группу из одновременно обрабатываемых точек:
	const uint32_t GROUP_SIZE = 1024u;
	
	// Формирование списка лучших кандидатов
	std::set<std::pair<Real, Vector>> candidates;
	
	// Сначала складываем точки группами по GROUP_SIZE:
	std::vector<Vector> group;
	for (uint32_t i = 0; i < nAllPoints; ++i) {
		// Получение текущего узла сетки:
		const auto net_point = net.GeneratePoint().coordinate;
		// Преобразование узла к точке в ограниченной min и max области:
		Vector curr(dim);
		for (uint32_t i = 0; i < dim; ++i) {
			curr[i] = net_point[i] * (max[i] - min[i]) + min[i];
		}
		if (i == nAllPoints - 1 || group.size() == GROUP_SIZE) {
			// Запускаем многопоточное вычисление значений во всех точках
			for (auto& it : calc_f_with_threads(f, group)) {
				candidates.insert(it);
				if (candidates.size() > nBestPoints) {
					auto max_element_it = std::max_element(candidates.begin(), candidates.end(), comparePairRealVector());
					candidates.erase(max_element_it);
				}
			}
			group.clear();
		} else {
			group.push_back(curr);
		}
	}
	
	// ----- Второй этап: запуск алгоритмов поиска локального минимума из выбранных точек -----
	// Подготовка (перекладываем точки из set в vector - возможен рост скорости при последовательном размещении в памяти точек):
	std::vector<std::pair<Real, Vector>> temp;
	for (auto & it : candidates) {
		temp.push_back(it);
	}
	
	// Многопоточная обработка кандидатов:
	auto answer = find_local_mins_with_threads(f, stop_condition, temp);
	
	// Итоговая сортировка всех найденных точек по неубыванию значения функции в них:
	for(uint32_t i = 0; i < candidates.size(); ++i)
	{
		std::sort(answer[i].begin(), answer[i].end(), comparePairRealGlobalTestData());
	}
	
	return answer;
}