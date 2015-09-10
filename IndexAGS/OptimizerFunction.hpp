#ifndef OPTIMIZER_FUNCTION_HPP
#define OPTIMIZER_FUNCTION_HPP

#include "OptimizerCoreGlobal.hpp"
#include "CoreUtils.hpp"

namespace optimizercore	{

	class OptimizerFunction
	{

	public:
		virtual ~OptimizerFunction() {};
		virtual double Calculate(const double*) = 0;
	};

	using OptimizerFunctionPtr = std::shared_ptr<OptimizerFunction>;
}

#endif