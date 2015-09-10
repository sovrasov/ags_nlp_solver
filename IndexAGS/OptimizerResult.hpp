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
		bool mIsInitialized;

	public:
		OptimizerResult();
		OptimizerResult(const OptimizerSolution& Solution);
		~OptimizerResult();

		OptimizerSolution GetSolution() const;

	private:
		void CheckIsInitialized() const;
	};
}

#endif