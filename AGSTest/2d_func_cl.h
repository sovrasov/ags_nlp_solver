#ifndef D_FUNCTION_CL_H
#define D_FUNCTION_CL_H
#include "OptimizerFunction.hpp"
//#include "vagris.h"
#include "gkls.h"

/*class grishagin_func : public func
{
	public:
	void set_func(int num)
	{
		set_random(num);
	}
	double calculate(double *x)
	{
		return random_func(x[0]+0.5, x[1]+0.5);
	}
};*/
class gkls_func : public optimizercore::OptimizerFunction
{
private:
	int m, num;
public:
	gkls_func()
	{}
	void set_func(int num)
	{
		GKLS_arg_generate(num);
	}
	double Calculate(const double *x) override
	{
		double y[5];
		for (unsigned i = 0; i < GKLS_dim; i++)
			y[i] = 2*x[i];
		return GKLS_D_func(y);
	}
	~gkls_func()
	{
	}
};
#endif