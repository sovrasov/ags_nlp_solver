#ifndef GKLSCoreGlobal_H
#define GKLSCoreGlobal_H

//Export directives
#if defined(_MSC_VER)
//  Microsoft 
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


#endif