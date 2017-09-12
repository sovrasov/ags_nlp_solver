#ifndef OPTIMIZER_SEARCH_INFORMATION_STORAGE_HPP
#define OPTIMIZER_SEARCH_INFORMATION_STORAGE_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerDataStructures.hpp"

namespace optimizercore
{

	class OptimizerSearchInformationStorage
	{

	private:

		size_t mCurrentSize;
		size_t mMaxSize;

		bool mIsInitialized;

	public:

		OptimizerSearchInformationStorage();
		OptimizerSearchInformationStorage(size_t maxSize);
		~OptimizerSearchInformationStorage();

		bool IsEmpty() const;

		bool Insert(OptimizerTrialPoint &point);

		void Clear();
	};

}

#endif