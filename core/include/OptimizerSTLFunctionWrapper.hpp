#ifndef OPTIMIZER_STL_FUNCTION_WRAPPER_HPP
#define OPTIMIZER_STL_FUNCTION_WRAPPER_HPP

#include"OptimizerCoreGlobal.hpp"
#include "OptimizerFunction.hpp"

#include <functional>

namespace optimizercore{

	class EXPORT_API OptimizerSTLFunctionWrapper final : public OptimizerFunction
	{

	private:

		std::function<double(const double*)> mFunction;
		int mCalculationsCounter;

	public:

		OptimizerSTLFunctionWrapper(std::function<double(const double*)>);

		void SetFunction(std::function<double(const double*)>);
		double Calculate(const double *x) override;
		int GetCalculationsCounter() const;
		void ResetCalculationsCounter();
	};

}
#endif