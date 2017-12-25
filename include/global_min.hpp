#pragma once

#include "sobolseqgenerator.h"
#include "math.hpp"
#include "nesterov.hpp"
#include "hessian_free.hpp"
#include "bfgs.hpp"

#include <vector>
#include <set>
#include <algorithm>

std::vector<std::pair<ld, Vector>>
find_absmin(Function f, uint32_t dim, uint32_t nBestPoints, uint32_t nAllPoints, Vector min, Vector max);