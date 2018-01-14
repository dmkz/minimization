#pragma once
#include "math.hpp"

class stopCondition {

private:
	Vector lastIterationPoint;
	Vector currentIterationPoint;
	ld lastIterationFunctionValue;
	ld currentIterationFunctionValue;
	long int iteration;
	long int iterationLimit;

public:
	stopCondition();
	stopCondition(long int iterationLimit);
	bool stopConditionHappened();
	void update(Vector lastIterationPoint, Vector currentIterationPoint, ld lastIterationFunctionValue, ld currentIterationFunctionValue, long int iteration);
};
