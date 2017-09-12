#ifndef HOOKE_JEEVES_LOCAL_METHOD_HPP
#define HOOKE_JEEVES_LOCAL_METHOD_HPP

#include "LocalMethodCommon.hpp"
#include "OptimizerFunction.hpp"
#include "OptimizerTask.hpp"

using namespace optimizercore;

namespace localoptimizer
{
	class EXPORT_API HookeJeevesLocalMethod final : public LocalMethodCommon
	{
	private:

		int mDimension;
		int mConstraintsNumber;
		OptimizerFunction** mFunctions;

		double *mStartPoint;
		double mEps;
		double mStep;
		double mStepMultiplier;

		double* mCurrentPoint;
		double* mCurrentResearchDirection;
		double* mPreviousResearchDirection;

		double MakeResearch(double*);
		void DoStep();
		double EvaluateTargetFunctiuon(const double*) const;

	public:
		HookeJeevesLocalMethod();
		~HookeJeevesLocalMethod();

		void SetProblem(OptimizerTask);
		void SetStartPoint(const double*, int);
		void SetEps(double);
		void SetInitialStep(double);
		void SetStepMultiplier(double);

		virtual void StartOptimization(double*) override;
	};
}

#endif