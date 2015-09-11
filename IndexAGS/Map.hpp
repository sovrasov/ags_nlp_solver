#ifndef MAP_HPP
#define MAP_HPP

#include "OptimizerCoreGlobal.hpp"

enum class MapType : int { Simple = 1, Linear = 2, Noninjective = 3, Set = 4, Rotated = 5};

EXPORT_API void mapd(double x, int m, double* y, int n, int key = 1);    /* map x to y         */
EXPORT_API void invmad(int, double *, int, int *, double *, int, int);  /* map y to x         */
EXPORT_API void xyd(double *xx, int m, double* y, int n);        /* get preimage       */

#endif