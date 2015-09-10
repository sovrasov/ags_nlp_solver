#ifndef OPTIMIZER_ALGORITHM_HPP
#define OPTIMIZER_ALGORITHM_HPP

#include "OptimizerDataStructures.hpp"
#include "OptimizerFunction.hpp"
#include "OptimizerCoreGlobal.hpp"
#include "OptimizerSearchSequence.hpp"
#include "OptimizerSolution.hpp"
#include "OptimizerResult.hpp"
#include "OptimizerTask.hpp"

#include <set>

namespace optimizercore	{

	enum class StopCriterionType : int { OptimalPoint = 0, Precision = 1 };

	class EXPORT_API OptimizerAlgorithm final
	{

	private:

		size_t mCurrentStorageSize;
		size_t mStorageReallocStep;
		size_t mStorageReallocBarrier;

		bool mLocalMixType;
		bool mIsAlgorithmMemoryAllocated;
		bool mIsParamsInitialized;
		bool mIsTaskInitialized;

		int mNumberOfThreads;
		int mLocalStartIterationNumber;
		int mMaxNumberOfIterations;
		int mRestrictionsNumber;
		int mMapTightness;

		OptimizerTask mTask;
		OptimizerFunction *mTargetFunction, **mRestrictions;

		OptimaizerInterval *mIntervalsForTrials;
		std::set<OptimizerTrialPoint> mSearchInformationStorage;

		int it_count, ins_pos, v_max, m_type, dim, local_percent, alpha;
		OptimizerTrialPoint *sq, result, *next_points;
		double *lip_const, *set_ranks, eps, *r, r_point, l_point, reserves;
		double **mNextPoints, *y;
		IndxSet **v_indexes;

		int GetIndex(OptimizerTrialPoint* oneDimPoint, double* point);
		void AllocMem();
		void InitializeInformationStorages();
		void UpdateLipConsts(IndxSet* set, int setNumber);
		void UpdateRanks(bool is_local);

	public:

		OptimizerAlgorithm();
		~OptimizerAlgorithm();

		void SetTask(OptimizerTask task);
		void SetThreadsNum(int num);
		void SetParameters(OptimizerParameters params);

		OptimizerResult StartOptimization(const double* x_opt, StopCriterionType stopType);

		double GetLipschitzConst(int fNumber) const;
		OptimizerSearchSequence GetSearchSequence() const;
	};

}
#endif