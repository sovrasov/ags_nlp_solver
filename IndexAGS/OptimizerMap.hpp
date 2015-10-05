#ifndef OPTIMIZER_MAP_HPP
#define OPTIMIZER_MAP_HPP

#include "OptimizerCoreGlobal.hpp"
#include "Map.hpp"

#define MAX_PREIMAGES 32

namespace optimizercore {

	class OptimizerMap
	{
	private:
		MapType mMapType;
		int mDimension;
		int mTightness;

	public:
		OptimizerMap();
		virtual ~OptimizerMap();
		OptimizerMap(int dimension, int tightness, MapType type);

		virtual void GetImage(double x, double y[]);
		virtual int GetAllPreimages(double* p, double xp[]);
	};
	
}

#endif
