#ifndef OPTIMIZER_SEARCH_SEQUENCE_HPP
#define OPTIMIZER_SEARCH_SEQUENCE_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerDataStructures.hpp"
#include "Map.hpp"
#include "OptimizerSpaceTransformation.hpp"
#include <set>

namespace optimizercore	{

	class EXPORT_API OptimizerSearchSequence final
	{

	private:

		size_t mSize;
		unsigned mDimension;
		MapType mMapType;
		unsigned mMapDensity;
		double* mPointsMemPtr;
		double* mValuesMemPtr;
		SharedVector mPoints;
		SharedVector mValues;
		OptimizerSpaceTransformation mSpaceTransform;
		bool mIsInitialized;

	public:

		size_t GetSize() const;
		unsigned GetMapDensity() const;
		unsigned GetDimention() const;
		MapType GetMapType() const;

		void GetPoint(int indx, double* x);
		double GetOneDimPoint(int indx);
		double GetValue(int indx);

		OptimizerSearchSequence();
		OptimizerSearchSequence(const std::set<OptimizerTrialPoint>& sequence,
			unsigned dimention,	MapType mapType, unsigned mapDensity,
			OptimizerSpaceTransformation transform);
		~OptimizerSearchSequence();

	private:
		void CheckIsInitialized() const;
	};

}
#endif