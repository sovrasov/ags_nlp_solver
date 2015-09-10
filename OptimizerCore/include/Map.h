/* Программы разверток */

#ifndef _MAP
#define _MAP

enum class MapType : int { Simple = 1, Linear = 2, Noninjective = 3};

__declspec(dllexport) void mapd(double x, int m, double* y, int n, int key);    /* map x to y         */
__declspec(dllexport) void invmad(int, double *, int, int *, double *, int, int);  /* map y to x         */
void xyd(double *xx, int m, double* y, int n);        /* get preimage       */

#endif