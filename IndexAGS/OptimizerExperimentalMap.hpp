#ifndef OPTIMIZER_EXPERIMANTAL_MAP_HPP
#define OPTIMIZER_EXPERIMANTAL_MAP_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerMap.hpp"

namespace optimizercore {

	class OptimizerExperimentalMap : public OptimizerMap
	{
	protected:


	public:
		OptimizerExperimentalMap();
		virtual ~OptimizerExperimentalMap();

		virtual void GetImage(double x, double y[]) override;
		virtual int GetAllPreimages(double* p, double xp[]) override;
	};
}

#endif
