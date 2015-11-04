#ifndef OPTIMIZER_SPACE_TRANSFORMATION_HPP
#define OPTIMIZER_SPACE_TRANSFORMATION_HPP

#include "OptimizerCoreGlobal.hpp"
#include"OptimizerFunction.hpp"

namespace optimizercore
{
	class EXPORT_API OptimizerSpaceTransformation final
	{
	private:

		SharedVector mLeftDomainBound;
		SharedVector mRightDomainBound;
		SharedVector mSpaceShiftValues;
		int mDimension;
		double mRho;
		double* mLeftBoundPtr;
		double* mRightBoundPtr;
		double* mSpaceShiftValuesPtr;

		OptimizerFunctionPtr mZeroConstraint;

		bool mNeedZeroConstraint;
		bool mIsInitialized;

	public:

		OptimizerSpaceTransformation();

		OptimizerSpaceTransformation(SharedVector leftBound, SharedVector rightBound,
			int domainDimention);

		void Transform(const double*, double*) const;
		void InvertTransform(const double*, double*) const;
		OptimizerFunctionPtr GetZeroConstraint() const;
		bool IsZeroConstraintActive() const;
		SharedVector GetLeftDomainBound() const;
		SharedVector GetRightDomainBound() const;
		int GetDomainDimension() const;

	};
}

#endif