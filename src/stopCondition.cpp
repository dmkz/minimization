#include "stopCondition.h"
#include "math.hpp"

class stopCondition {

	Vector lastIterationPoint;
	Vector currentIterationPoint;
	ld lastIterationFunctionValue;
	ld currentIterationFunctionValue;
	long int iteration;
	long int iterationLimit = 10000;

	stopCondition() {
	}

	stopCondition(long int iterationLimit) {
		this->iterationLimit = iterationLimit;
	}

	bool stopConditionHappened() {
		return iteration > iterationLimit;
	}
	void update(Vector lastIterationPoint, Vector currentIterationPoint, ld lastIterationFunctionValue, ld currentIterationFunctionValue, long int iteration) {
		this->lastIterationPoint = lastIterationPoint;
		this->currentIterationPoint = currentIterationPoint;
		this->lastIterationFunctionValue = lastIterationFunctionValue;
		this->currentIterationFunctionValue = currentIterationFunctionValue;
		this->iteration = iteration;
	}
};
