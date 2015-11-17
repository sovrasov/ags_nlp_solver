#ifndef OPTIMIZER_QUEUE_HPP
#define OPTIMIZER_QUEUE_HPP

#include "OptimizerCoreGlobal.hpp"
#include "OptimizerDataStructures.hpp"

namespace optimizercore {

	class OptimizerQueue final
	{
	protected:
		int MaxSize;
		int CurSize;
		OptimizerInterval *pMem;
		int GetIndOfMinElem();
		void DeleteMinElem();
		void ReBuild(int Index);

	public:
		// размер очереди должен быть равен 2^k - 1
		OptimizerQueue(int _MaxSize = 131071);
		~OptimizerQueue() { delete[] pMem; }

		int GetSize() const;
		bool IsEmpty() const;
		bool IsFull() const;

		void Push(const OptimizerInterval &value);
		void PushWithPriority(const OptimizerInterval &value);
		OptimizerInterval Pop();

		void Clear();
	};
}
#endif