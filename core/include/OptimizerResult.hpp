#ifndef OPTIMIZER_RESULT_HPP
#define OPTIMIZER_RESULT_HPP 

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerSolution.hpp"

namespace optimizercore
{
	class EXPORT_API OptimizerResult final
	{
	private:
		OptimizerSolution mSolution;
		SharedIntVector mFunctionalsCalculationStat;
		int mNumberOfFunctionals;
		bool mIsInitialized;

	public:
		OptimizerResult();
		OptimizerResult(const OptimizerSolution& Solution);
		OptimizerResult(const OptimizerSolution& Solution,
			SharedIntVector& functionalsCalculationStat,
			int numberOfFunctionals);

		~OptimizerResult();

		OptimizerSolution GetSolution() const;
		int GetNumberOfCalculations(int fNumber) const;
		int GetNumberOfFunctionals() const;

	private:
		void CheckIsInitialized() const;
	};
}

#endif