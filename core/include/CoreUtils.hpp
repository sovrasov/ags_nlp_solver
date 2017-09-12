#ifndef CoreUtils_HPP
#define CoreUtils_HPP

#include "OptimizerCoreGlobal.hpp"

namespace optimizercore{
	namespace utils {

		template<typename T>
		struct array_deleter {
			void operator()(T const * p) {
				delete[] p;
			}
		};

		template<typename T>
		T** AllocateMatrix(size_t rowsNumber, size_t columnsNumber)
		{
			T** matrix;
			matrix = new T*[rowsNumber];
			for (unsigned i = 0; i < rowsNumber; i++)
				matrix[i] = new T[columnsNumber];
			return matrix;
		}

		template<typename T>
		void DeleteMatrix(T** matrix, size_t rowsNumber)
		{
			for (unsigned i = 0; i < rowsNumber; i++)
				delete[] matrix[i];
			delete[] matrix;
		}

		template<typename T>
		int reAllocMem(T **ptr, size_t cur_size, size_t step)
		{
			T *tmp = new T[cur_size + step];
			if (tmp != nullptr)
			{
				for (size_t i = 0; i < cur_size; i++)
					tmp[i] = (*ptr)[i];
				delete[] * ptr;
				*ptr = tmp;
				return 1;
			}
			return 0;
		}

		template<typename T>
		int insert(T *mem, T value, int size)
		{
			int left = 0, middle = 0;
			int right = size;
			do {
				middle = (left + right) / 2;
				if (value.x > mem[middle].x)
					left = middle + 1;
				else
					right = middle - 1;
			} while (left <= right);

			for (int j = size; j > left; j--)
				mem[j] = mem[j - 1];
			mem[left] = value;
			return left;
		}

		double EXPORT_API NormNDimMax(const double *x, const double* y, int dim);
		double Pow2(double x);
		int sgn(double x);
	}
}
#endif