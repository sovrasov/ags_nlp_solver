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

  bool mIsInitialized;
private:
  MapType mMapType;
  int mMapKey;

public:
  Evolvent();
  Evolvent(int dimension, int tightness, double* lb, double*ub, MapType type = Simple);
  ~Evolvent();

  void GetImage(double x, double y[]) const;
  int GetAllPreimages(double* p, double xp[]);
};
