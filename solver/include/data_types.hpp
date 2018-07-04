#pragma once

const unsigned solverMaxDim = 5;
const unsigned solverMaxConstraints = 5;

struct Trial
{
  double x;
  double y[solverMaxDim];
  double g[solverMaxConstraints + 1];
  int idx;
  Trial() {}
  Trial(double _x) : x(_x) {}
};

struct Interval
{
  Trial pl;
  Trial pr;
  double R;
  double delta;
  Interval() {}
  Interval(const Trial& _pl, const Trial& _pr) : pl(_pl), pr(_pr) {}
};

struct CompareIntervals
{
  bool operator() (const Interval* i1, const Interval* i2)
  {
    return i1->pl.x < i2->pl.x;
  }
};

class CompareByR
{
public:
  bool operator() (const Interval* i1, const Interval* i2)
  {
    return i1->R < i2->R;
  }
};
