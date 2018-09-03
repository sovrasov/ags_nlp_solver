#pragma once

#include <algorithm>
#include <stdexcept>
#include <functional>

  // возвращает двоичный логарифм от val
inline unsigned int log2(unsigned int val)
{
  //if (val == 0)
  //	throw exception("log2(0) не определен\n");

  unsigned int result = 0;

  while (val)
  {
    val >>= 1;
    ++result;
  }

  return result - 1;
}

//	T - тип элементов в куче
//  Container - тип контейнера, используемый для хранения кучи
//  Compare - функция сравнения, чтобы определить как сравнивать элементы в куче

// Считаем, что корень - Max Level
template<class T, class Compare = std::less<T> >
class MinMaxHeap
{
  T* m_heap; // здесь хранится содержимое кучи
  Compare m_compare; // объект для сравнения
  unsigned m_heapsize; // размер кучи
  unsigned m_currentheapsize; // текущее количество элементов в куче


  static inline unsigned int parent(unsigned int index)
  {
    return (index - 1) / 2; // возвращает индекс родителя узла, заданного index-ом
  }

  /*
  inline void swapLinkedElems(T& arg1, T& arg2)
  {
    if (arg1.pLinkedElement != NULL && arg2.pLinkedElement != NULL)
      std::swap(arg1.pLinkedElement->pLinkedElement, arg2.pLinkedElement->pLinkedElement);
    else if (arg1.pLinkedElement != NULL)
      arg1.pLinkedElement->pLinkedElement = &arg2;
    else if (arg2.pLinkedElement != NULL)
      arg2.pLinkedElement->pLinkedElement = &arg1;
  }
  */

  // ------------------------------------------------------------------------------------------------
  static inline unsigned int leftChild(unsigned int index)
  {
    return 2 * index + 1;  // возвращает индекс левого сына узла, заданного index-ом
  }

  // ------------------------------------------------------------------------------------------------
  static inline unsigned int rightChild(unsigned int index)
  {
    return 2 * index + 2;  // возвращает индекс правого сына узла, заданного index-ом
  }

  // ------------------------------------------------------------------------------------------------
  static inline bool isOnMinLevel(unsigned int index)
  {
    return log2(index + 1) % 2 == 1; // true - если вершина, соответствующая index находится на Min Level
  }

  // ------------------------------------------------------------------------------------------------
  static inline bool isOnMaxLevel(unsigned int index)
  {
    return !isOnMinLevel(index); // true - если вершина, соответствующая index находится на Max Level
  }

  // ------------------------------------------------------------------------------------------------
  template<bool MaxLevel>
  int trickleUp_(unsigned int index) // вспомогательный метод для всплытия
  {
    if (index == 0) // не можем всплывать дальше
      return index;

    unsigned int index_grandparent = parent(index); // найдем первый родительский пройденный уровень

    if (index_grandparent == 0) // если такого нет, выходим
      return index;

    index_grandparent = parent(index_grandparent); // найти прародителя

    if (m_compare(m_heap[index], m_heap[index_grandparent]) ^ MaxLevel) // убедимся, что нужно поменяться местами с прародителем
    {
      //swapLinkedElems(m_heap[index_grandparent], m_heap[index]);
      std::swap(m_heap[index_grandparent], m_heap[index]);

      return trickleUp_<MaxLevel>(index_grandparent);
    }
    return index;
  }

  // ------------------------------------------------------------------------------------------------
  int trickleUp(unsigned int index) // размещаем узел на соответствующем уровне (Min или Max)
  {
    if (index == 0)  // не можем всплывать дальше
      return index;

    unsigned int index_parent = parent(index); // найдем первый родительский пройденный уровень

    if (isOnMinLevel(index))
    {
      // убедимся, что нужно поменяться местами с родителем
      if (m_compare(m_heap[index_parent], m_heap[index]))
      {
        //swapLinkedElems(m_heap[index_parent], m_heap[index]);
        std::swap(m_heap[index_parent], m_heap[index]);

        return trickleUp_<true>(index_parent);
      }
      else
        return trickleUp_<false>(index);
    }
    else
    {
      // убедимся, что нужно поменяться местами с родителем
      if (m_compare(m_heap[index], m_heap[index_parent]))
      {
        //swapLinkedElems(m_heap[index_parent], m_heap[index]);
        std::swap(m_heap[index_parent], m_heap[index]);

        return trickleUp_<false>(index_parent);
      }
      else
        return trickleUp_<true>(index);
    }
  }

  // ------------------------------------------------------------------------------------------------
  template<bool MaxLevel>
  void trickleDown_(unsigned int index) // вспомогательный метод для погружения
  {
    //  if ( index >= m_currentheapsize ) // убедимся, что элемент существует
     //     throw exception("Элемент с таким индексом не существует\n");

    unsigned int smallestNode = index; // храним индекс наименьшего узла
    unsigned int left = leftChild(index); // получаем правого сына

    if (left < m_currentheapsize && (m_compare(m_heap[left], m_heap[smallestNode]) ^ MaxLevel)) // проверяем левого и правого сыновей
      smallestNode = left;
    if (left + 1 < m_currentheapsize && (m_compare(m_heap[left + 1], m_heap[smallestNode]) ^ MaxLevel))
      smallestNode = left + 1;

    unsigned int leftGrandchild = leftChild(left); // проверяем внуков
    for (unsigned int i = 0; i < 4 && leftGrandchild + i < (unsigned int)m_currentheapsize; ++i)
      if (m_compare(m_heap[leftGrandchild + i], m_heap[smallestNode]) ^ MaxLevel)
        smallestNode = leftGrandchild + i;

    if (index == smallestNode) // если текущий узел наименьший, ничего не делаем и выходим
      return;

    //swapLinkedElems(m_heap[index], m_heap[smallestNode]);
    std::swap(m_heap[index], m_heap[smallestNode]); // меняем местами текущий узел и наименьший

    if (smallestNode - left > 1)
    {
      // если родитель наименьшего узла больше, чем сам узел, меняем местами
      if (m_compare(m_heap[parent(smallestNode)], m_heap[smallestNode]) ^ MaxLevel) {
        //swapLinkedElems(m_heap[parent(smallestNode)], m_heap[smallestNode]);
        std::swap(m_heap[parent(smallestNode)], m_heap[smallestNode]);
      }

      trickleDown_<MaxLevel>(smallestNode);
    }
  }

  // ------------------------------------------------------------------------------------------------
  void trickleDown(unsigned int index) // погружение
  {
    if (isOnMinLevel(index))
      trickleDown_<false>(index);
    else
      trickleDown_<true>(index);
  }

  // ------------------------------------------------------------------------------------------------
  unsigned int findMinIndex() const // поиск индекса наименьшего узла
  {
    switch (m_currentheapsize)
    {
    case 0:
      // куча пуста
      throw std::runtime_error("Куча пуста\n");
      break;
    case 1:
      // в куче только один элемент
      return 0;
    case 2:
      // в куче 2 элемента => сын должен быть минимумом
      return 1;
    default:
      // в куче больше 2х элементов
      return m_compare(m_heap[1], m_heap[2]) ? 1 : 2;
    }
  }

  // ------------------------------------------------------------------------------------------------
  void deleteElement(unsigned int index) // удаление элемента из кучи
  {
    //     if (index >= (unsigned int)m_currentheapsize) // проверить существование элемента
     //        throw exception("Элемент с таким индексом не существует\n");

         // если мы удаляем последний элемент из кучи
    if (index == m_currentheapsize - 1)
    {
      m_currentheapsize--;
      return;
    }

    //swapLinkedElems(m_heap[index], m_heap[m_currentheapsize - 1]);
    std::swap(m_heap[index], m_heap[m_currentheapsize - 1]); // меняем местами элемент с последним в куче

    m_currentheapsize--; // удаляем последний элемент из кучи

    trickleDown(index); // погружаем тот элемент, который поместили на место index
  }

public:
  MinMaxHeap(unsigned heapsize) : m_heap(NULL)
  {
    m_heapsize = heapsize;
    m_heap = new T[m_heapsize];
    m_currentheapsize = 0;
  }

  ~MinMaxHeap()
  {
    if (m_heap != NULL)
      delete[] m_heap;
  }

  void clear()
  {
    m_currentheapsize = 0;
  }

  bool empty() const
  {
    return m_currentheapsize == 0; // проверка на пустоту кучи
  }

  bool full() const
  {
    return m_currentheapsize == m_heapsize;
  }

  unsigned int size() const
  {
    return (unsigned int)m_currentheapsize; // количество элементов в куче
  }


  T* push(const T & val) // добавление элемента в кучу
  {
    if (m_currentheapsize < m_heapsize)
    {
      m_heap[m_currentheapsize] = val; // добавляем в конец
      m_currentheapsize++;
      return m_heap + trickleUp(m_currentheapsize - 1); // всплываем
    }
    return NULL;
    //	else
    //		throw exception("Вставка в кучу невозможна! Переполнение!\n");
  }

  const T & findMax() const // возвращает максимум за O(1) без удаления из кучи
  {
    //      if (empty()) // проверка на пустоту
     //         throw exception("Куча пуста \n");

    return m_heap[0];
  }

  const T & findMin() const
  {
    return m_heap[findMinIndex()];  // возвращает минимум за O(1) без удаления из кучи
  }

  T popMax() // извлечение (с удалением из кучи) максимального элемента
  {
    //   if (empty()) // проверяем кучу на пустоту
     //      throw exception("Куча пуста \n");

    T temp = m_heap[0]; // сохраняем максимум
    int delIndx = 0;
    deleteElement(delIndx); // удаляем

    return temp;
  }

  T pop()
  {
    return popMax(); // извлечение (с удалением из кучи) максимального элемента
  }


  T popMin() // извлечение (с удалением из кучи) минимального элемента
  {
    //   if (empty()) // проверяем не пуста ли куча
     //      throw exception("Куча пуста \n");

    unsigned int smallest = findMinIndex(); // сохраняем индекс минимума
    T temp = m_heap[smallest]; // сохраняем минимальное значение
    deleteElement(smallest); // удаляем

    return temp;
  }

  void deleteElement(const T* ptr)
  {
    unsigned index = unsigned(ptr - m_heap);
    //if(index >= 0 && index < m_heapsize)
    deleteElement(index);
    //	else
    //		throw std::exception("Куча пуста \n");
  }

  T* getHeapMemPtr() const
  {
    return m_heap;
  }
};
