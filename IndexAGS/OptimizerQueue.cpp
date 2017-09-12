#include "OptimizerQueue.hpp"

#include <cmath>
#include <algorithm>

using namespace optimizercore;

OptimizerQueue::OptimizerQueue(int _MaxSize)
{
	MaxSize = _MaxSize;
	CurSize = 0;
	pMem = new OptimizerInterval[MaxSize];
}

// ------------------------------------------------------------------------------------------------
/*optimizercore::OptimizerQueue::~OptimizerQueue()
{
	delete[] pMem;
}
*/
// ------------------------------------------------------------------------------------------------
bool OptimizerQueue::IsEmpty() const
{
	return CurSize == 0;
}
// ------------------------------------------------------------------------------------------------
int OptimizerQueue::GetSize() const
{
	return CurSize;
}

// ------------------------------------------------------------------------------------------------
bool OptimizerQueue::IsFull() const
{
	return CurSize == MaxSize;
}

// ------------------------------------------------------------------------------------------------
void OptimizerQueue::Push(const OptimizerInterval &value)
{
	if (IsFull())
	{
		int MinInd = GetIndOfMinElem();
		if (value.R > pMem[MinInd].R)
			DeleteMinElem();
		else
			return;
	}
	CurSize++;
	pMem[CurSize - 1] = value;
	if (CurSize > 1)
		ReBuild(CurSize - 1);
}

// ------------------------------------------------------------------------------------------------
void OptimizerQueue::PushWithPriority(const OptimizerInterval &value)
{
	if (IsEmpty())
	{
		CurSize++;
		pMem[CurSize - 1] = value;
	}
	else
	{
		int MinInd = GetIndOfMinElem();

		// В очереди должны быть элементы с одинаковыми характеристиками!!!
		if (value.R >= pMem[MinInd].R)
		{
			if (IsFull())
				DeleteMinElem();
			CurSize++;
			pMem[CurSize - 1] = value;
			if (CurSize > 1)
				ReBuild(CurSize - 1);
		}
	}
}

void optimizercore::OptimizerQueue::DeleteInterval(const OptimizerInterval & value)
{
	for (int i = 0; i < CurSize; i++)
		if (pMem[i].right.x == value.right.x)
		{
			pMem[i] = pMem[CurSize - 1];
			CurSize--;
			if (CurSize > 1)
				ReBuild(i);
			break;
		}
}

// ------------------------------------------------------------------------------------------------
OptimizerInterval OptimizerQueue::Pop()
{
	OptimizerInterval tmp = pMem[0];
	pMem[0] = pMem[CurSize - 1];
	CurSize--;
	if (CurSize > 1)
		ReBuild(0);

	return tmp;
}

// ------------------------------------------------------------------------------------------------
int OptimizerQueue::GetIndOfMinElem()
{
	int i, StartIndex;
	double min = HUGE_VAL;
	int MinIndex = -1;

	// Находим левую границу для поиска минимума. Левая граница - первый лист в куче.
	// Если последний уровень дерева заполнен не до конца, то первый лист в предпоследнем уровне.
	// Если у узла есть потомки, то в этом узле не может быть минимума (по свойству d-кучи)
	if (CurSize % 2) // Элемент CurSize - левый потомок своего предка
		StartIndex = (CurSize - 1) / 2; // предок первого свободного листа в куче - граница поиска
	else // Элемент CurSize - правый потомок своего предка
		StartIndex = (CurSize - 1) / 2 + 1; // граница поиска - следующий предок за предком CurSize

	for (i = StartIndex; i < CurSize; i++)
		if (min > pMem[i].R)
		{
			MinIndex = i;
			min = pMem[i].R;
		}

	return MinIndex;
}

// ------------------------------------------------------------------------------------------------
void OptimizerQueue::DeleteMinElem()
{
	int MinInd = GetIndOfMinElem();
	pMem[MinInd] = pMem[CurSize - 1];
	CurSize--;
	if (CurSize > 1)
		ReBuild(MinInd);
}

// ------------------------------------------------------------------------------------------------
void OptimizerQueue::ReBuild(int Index)
{
	int i, j, k, level = 0;
	if (Index == 0) // восстановление структуры двоичной кучи от корня (погружение)
	{
		i = Index;
		j = 2 * i + 1; // левый потомок корня (всегда существует)
		k = 2 * i + 2; // правый потомок корня (может отсутствовать)
	//    if (k < CurSize - 1) // оба потомка существуют
		if (k < CurSize) // оба потомка существуют
			if (pMem[j].R < pMem[k].R)
				j = k;
		// теперь j - наибольший потомок узла i
		while (true)
		{
			if (pMem[i].R >= pMem[j].R)
				break;
			std::swap(pMem[i], pMem[j]);

			i = j;
			j = 2 * i + 1;
			k = 2 * i + 2;
			if (j > CurSize - 1) // у узла i нет потомков
				break;
			//      if (k < CurSize - 1) // оба потомка существуют
			if (k < CurSize) // оба потомка существуют
				if (pMem[j].R < pMem[k].R)
					j = k;
		}
	}
	else // восстановление структуры двоичной кучи от листа (всплытие)
	{
		i = Index;
		j = (i - 1) / 2; // предок узла i
		while ((i > 0) && (pMem[j].R <= pMem[i].R)) //А можно ли тут поставить <= вместо < ?
		{
			std::swap(pMem[i], pMem[j]);
			i = j;
			j = (i - 1) / 2;
		}
	}
}

void OptimizerQueue::Clear()
{
	CurSize = 0;
}
// - end of file ----------------------------------------------------------------------------------
