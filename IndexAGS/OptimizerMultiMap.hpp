#ifndef OPTIMIZER_MULTIMAP_HPP
#define OPTIMIZER_MULTIMAP_HPP

#include "OptimizerCoreGlobal.hpp"
#include "CoreUtils.hpp"
#include "OptimizerMap.hpp"

namespace optimizercore
{
	enum class MultimapType {Rotated, Set};

	class EXPORT_API OptimizerMultiMap final : public OptimizerMap
	{
	private:
		class MatrixMemHolder final {
		private:
			int **mMem, mMemSize;
		public:
			MatrixMemHolder(int **mem, int size) : mMem(mem), mMemSize(size)
			{	}
			~MatrixMemHolder() {
				utils::DeleteMatrix(mMem, mMemSize);
			}
		};
	private:
		int mNumberOfMaps;
		int mRotatedMapPlanesCount;
		MultimapType mMapType;

		int **mRotationPlanes;
		std::shared_ptr<MatrixMemHolder> mRotationPlanesHolder;
		SharedVector mTmpVector; //smart ptr for p2
		double *p2; 

		void(OptimizerMultiMap:: *mCurrentMap)(double x, double y[]);
		int(OptimizerMultiMap:: *mCurrentInvMap)(double* p, double xp[]);

		void GetSetMapImage(double x, double y[]);
		void GetRotatedMapImage(double x, double y[]);
		int GetAllSetMapPreimages(double* p, double xp[]);
		int GetAllRotatedMapPreimages(double* p, double xp[]);

	public:
		OptimizerMultiMap();
		OptimizerMultiMap(MultimapType mapType, int n, int m, int l);
		~OptimizerMultiMap();

		int GetNumberOfMaps() const;
		void GetImage(double x, double y[]);
		int GetAllPreimages(double* p, double xp[]);

	private:
		void InitRotatedMap();
	};

}
#endif