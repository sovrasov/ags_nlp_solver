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

	class EXPORT_API OptimizerAlgorithm final
	{

	private:

		size_t mCurrentStorageSize;
		size_t mStorageReallocStep;

		bool mLocalMixType;
		bool mIsAlgorithmMemoryAllocated;
		bool mIsParamsInitialized;
		bool mIsTaskInitialized;
		bool mNeedLocalVerification;

		int mNumberOfThreads;
		int mLocalStartIterationNumber;
		int mMaxNumberOfIterations;
		int mRestrictionsNumber;
		int mMapTightness;
		int mMethodDimention;
		int mAlpha;
		int mLocalMixParameter;
		int mMapType;

		OptimizerTask mTask;
		OptimizerSpaceTransformation mSpaceTransform;
		OptimizerFunction *mTargetFunction, **mRestrictions;

		OptimaizerInterval *mIntervalsForTrials;
		std::set<OptimizerTrialPoint> mSearchInformationStorage;
		OptimizerTrialPoint mOptimumEvaluation, *mNextTrialsPoints;

		double *lip_const, *set_ranks, eps, *r, reserves;
		double **mNextPoints;
		IndxSet **v_indexes;

		int GetIndex(OptimizerTrialPoint* oneDimPoint, double* point);
		void AllocMem();
		void InitializeInformationStorages();
		void UpdateLipConsts(IndxSet* set, const OptimizerTrialPoint& value);
		int UpdateRanks(bool isLocal);
		bool InsertNewTrials(int trailsNumber);
		OptimizerSolution DoLocalVerification(OptimizerSolution startPoint);

	public:

		OptimizerAlgorithm();
		~OptimizerAlgorithm();

		void SetTask(OptimizerTask task);
		void SetThreadsNum(int num);
		void SetParameters(OptimizerParameters params);

		OptimizerResult StartOptimization(const double* xOpt,
			StopCriterionType stopType);

		double GetLipschitzConst(int fNumber) const;
		OptimizerSearchSequence GetSearchSequence() const;
	};

}
#endif