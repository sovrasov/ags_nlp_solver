#ifndef VAGrisFunction_H
#define VAGrisFunction_H

#include "VAGrisCoreGlobal.h"

namespace vagrisfunction {

#ifndef M_PI
#define M_PI 3.141592654
#endif

	class EXPORT_API  VAGrisFunction final {

	private:
		int mFunctionNumber;
		unsigned char icnf[45];
		double af[7][7], bf[7][7], cf[7][7], df[7][7];

		double rndm20(unsigned char k[]);
		void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);

	public:
		VAGrisFunction();
		~VAGrisFunction();

		PROPERTY(int, FunctionNumber); 

		double EvaluateFunction(double x, double y) const;
		double EvaluateXDerivative(double x, double y) const;
		double EvaluateYDerivative(double x, double y) const;
		double GetMinimumXCoordinate(int fNumber) const;
		double GetMinimumYCoordinate(int fNumber) const;
	};
}

#endif