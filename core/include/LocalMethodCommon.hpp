#ifndef LOCAL_METHOD_COMMON_HPP
#define LOCAL_METHOD_COMMON_HPP

namespace localoptimizer
{
	class LocalMethodCommon
	{
	public:
		virtual void StartOptimization(double*) = 0;
	};
}

#endif