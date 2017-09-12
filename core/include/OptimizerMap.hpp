#ifndef OPTIMIZER_MAP_HPP
#define OPTIMIZER_MAP_HPP

#include "OptimizerCoreGlobal.hpp"
#include "Map.hpp"

namespace optimizercore {

#define MAX_PREIMAGES 32

	class EXPORT_API OptimizerMap
	{
	protected:
		int mDimension;
		int mTightness;

		bool mIsInitialized;
	private:
		MapType mMapType;
		int mMapKey;

	public:
		OptimizerMap();
		OptimizerMap(int dimension, int tightness, MapType type);
		virtual ~OptimizerMap();

		virtual void GetImage(double x, double y[]);
		virtual int GetAllPreimages(double* p, double xp[]);
	};
	
}

#endif
