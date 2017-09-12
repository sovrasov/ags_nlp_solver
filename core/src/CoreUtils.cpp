#include "CoreUtils.hpp"
#include <cmath>

namespace optimizercore{
	namespace utils {

		double NormNDimMax(const double *x, const double* y, int dim)
		{
			double max = std::abs(y[0] - x[0]);
			for (int i = 1; i < dim; i++)
				max = std::fmax(max, std::abs(y[i] - x[i]));
			return max;
		}

		double Pow2(double x)
		{
			return x*x;
		}

		int sgn(double x)
		{
			if (x > 0)
				return 1;
			else 
				return -1;
		}

	}
}
