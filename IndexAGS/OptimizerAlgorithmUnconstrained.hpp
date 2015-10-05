#ifndef OPTIMIZER_ALGORITHM_UNCONSTRAINED_HPP
#define OPTIMIZER_ALGORITHM_UNCONSTRAINED_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerTask.hpp"
#include "OptimizerSolution.hpp"
#include "OptimizerFunction.hpp"
#include "OptimizerDataStructures.hpp"
#include "OptimizerSolution.hpp"
#include "OptimizerResult.hpp"
#include "OptimizerSearchSequence.hpp"

#include <set>

namespace optimizercore
{
	
	class EXPORT_API OptimizerAlgorithmUnconstrained final
	{

	private:

		bool mLocalMixType;
		bool mIsAlgorithmMemoryAllocated;
		bool mIsParamsInitialized;
		bool mIsTaskInitialized;
		bool mNeedLocalVerification;

		int mNumberOfThreads;
		int mLocalStartIterationNumber;
		int mMaxNumberOfIterations;
		int mMapTightness;
		int mMethodDimension;
		int mAlpha;
		int mLocalMixParameter;
		int mMapType;

		OptimizerSpaceTransformation mSpaceTransform;
		OptimizerFunction *mTargetFunction;
		OptimizerFunctionPtr mTargetFunctionSmartPtr;

		OptimizerInterval *mIntervalsForTrials;
		std::set<OptimizerTrialPoint> mSearchInformationStorage;
		OptimizerTrialPoint mOptimumEvaluation, *mNextTrialsPoints;

		LocalTuningMode mLocalTuningMode;

		double mGlobalM, mZ, eps, r, mMaxIntervalNorm;
		double **mNextPoints;

		void AllocMem();
		void InitializeInformationStorage();
		void UpdateGlobalM(std::set<OptimizerTrialPoint>::iterator&);
		int UpdateRanks(bool isLocal);
		bool InsertNewTrials(int trailsNumber);
		OptimizerSolution DoLocalVerification(OptimizerSolution startPoint);

	public:
		OptimizerAlgorithmUnconstrained();
		~OptimizerAlgorithmUnconstrained();

		void SetTask(OptimizerFunctionPtr function,
			OptimizerSpaceTransformation spaceTransform);
		void SetThreadsNum(int num);
		void SetParameters(OptimizerParameters params);

		OptimizerResult StartOptimization(const double* xOpt,
			StopCriterionType stopType);

		double GetLipschitzConst() const;
		OptimizerSearchSequence GetSearchSequence() const;

	};
}
#endif 