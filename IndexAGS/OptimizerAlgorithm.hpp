#ifndef OPTIMIZER_ALGORITHM_HPP
#define OPTIMIZER_ALGORITHM_HPP

#include "OptimizerDataStructures.hpp"
#include "OptimizerFunction.hpp"
#include "OptimizerCoreGlobal.hpp"
#include "OptimizerSearchSequence.hpp"
#include "OptimizerSolution.hpp"
#include "OptimizerResult.hpp"
#include "OptimizerTask.hpp"
#include "OptimizerMultiMap.hpp"
#include "OptimizerMap.hpp"
#include "OptimizerQueue.hpp"

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
		bool mNeedQueueRefill;

		int mNumberOfThreads;
		int mNumberOfMaps;
		int mLocalStartIterationNumber;
		int mMaxNumberOfIterations;
		int mRestrictionsNumber;
		int mMapTightness;
		int mMethodDimension;
		int mAlpha;
		int mLocalMixParameter;
		int mMapType;

		OptimizerTask mTask;
		OptimizerMap *mPMap;
		OptimizerSpaceTransformation mSpaceTransform;
		OptimizerFunction *mTargetFunction, **mRestrictions;

		OptimizerInterval *mIntervalsForTrials;
		std::set<OptimizerTrialPoint> mSearchInformationStorage;
		OptimizerTrialPoint mOptimumEvaluation, *mNextTrialsPoints;

		OptimizerQueue mQueue;

		double *lip_const, *set_ranks, eps, *r, reserves;
		double **mNextPoints;
		MultimapIndxSet **v_indexes;

		int GetIndex(OptimizerTrialPoint* oneDimPoint, double* point);
		void AllocMem();
		void InitializeInformationStorages();
		void UpdateLipConsts(MultimapIndxSet* set, const OptimizerTrialPoint& value);
		double CalculateGlobalR(const OptimizerInterval&);
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