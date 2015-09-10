#ifndef OPTIMIZER_MULTIMAP_HPP
#define OPTIMIZER_MULTIMAP_HPP

namespace optimizer
{
	enum class MultimapType {Rotated, Set};
	class OptimizerMultiMap final
	{
	protected:
		int mDimension;
		int mTightness;
		int mNumberOfMaps;
		MultimapType mMapType;

		int **mRotationPlanes;
		double *p2;

		bool mIsInitialized;
	public:
		OptimizerMultiMap();
		OptimizerMultiMap(MultimapType mapType, int n, int m, int l);
		~OptimizerMultiMap();

		int GetNumberOfMaps() const;
		void GetImage(double x, double y[]);
		void GetAllPreimages(double* p, double xp[]);

	protected:
		void InitRotatedMap();
	};

}
#endif