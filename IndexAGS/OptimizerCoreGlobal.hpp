#ifndef OptimizerCoreClobal_HPP
#define OptimizerCoreClobal_HPP

//Export directives
#if defined(_MSC_VER)
//  Microsoft 
#pragma warning(disable: 4251)
#define EXPORT_API __declspec(dllexport)
#define STD_CALL __stdcall
#elif defined(__GNUG__)
//  GCC
#define EXPORT_API __attribute__((visibility("default")))
#define STD_CALL
#else
//  Unknown
#define EXPORT_API
#define STD_CALL
#pragma warning Unknown dynamic link import/export semantics.
#endif

#define PROPERTY(T, N)     \
	T Get ## N() const;     \
	void Set ## N(T value);

#include <memory>
namespace optimizercore {
	using SharedVector = std::shared_ptr<double>;
	using SharedIntVector = std::shared_ptr<int>;
}
#endif