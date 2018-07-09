#pragma once
#include <vector>

#define MAX_PREIMAGES 32

enum MapType {
  Simple = 1, Linear = 2, Noninjective = 3
};

class Evolvent
{
protected:
  int mDimension;
  int mTightness;

  std::vector<double> mRho;
  std::vector<double> mShiftScalars;
  std::vector<double> p2;

  void TransformToStandardCube(const double *y, double *z);
  void TransformToSearchDomain(const double *y, double *z);

  bool mIsInitialized;
private:
  MapType mMapType;
  int mMapKey;

public:
  Evolvent();
  Evolvent(int dimension, int tightness, const double* lb, const double* ub, MapType type = Simple);
  ~Evolvent();

  virtual void GetImage(double x, double y[]);
  virtual int GetAllPreimages(const double* p, double xp[]);
};
