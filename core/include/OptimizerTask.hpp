#ifndef OPTIMIZER_TASK_HPP
#define	OPTIMIZER_TASK_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerFunction.hpp"
#include "OptimizerSpaceTransformation.hpp"

namespace optimizercore
{

	class EXPORT_API OptimizerTask final
	{

	private:

		std::shared_ptr<OptimizerFunctionPtr> mTaskFunctions;

		unsigned mNumberOfRestrictions;
		unsigned mTaskDimension;
		SharedVector mOptimumPoint;
		OptimizerSpaceTransformation mSpaceTransform;

		bool mIsInitialized;

	public:
		OptimizerTask(std::shared_ptr<OptimizerFunctionPtr> taskFunctions,
			unsigned numberOfRestrictions, unsigned taskDimention,
			SharedVector leftBound, SharedVector rightBound);
		OptimizerTask();
		~OptimizerTask();

		std::shared_ptr<OptimizerFunctionPtr> GetTaskFunctions() const;
		unsigned GetNumberOfRestrictions() const;
		unsigned GetTaskDimension() const;
		PROPERTY(SharedVector, OptimumPoint);
		OptimizerSpaceTransformation GetSpaceTransformation() const;

	private:
		void CheckIsInitialized() const;
	};
}
#endif