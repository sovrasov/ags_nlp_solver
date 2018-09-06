/**
 * @file MinMaxHeap.hpp
 * @author John Sullivan (jsull003 at ucr.edu), Vladislav Sovrasov
 * @date December 28, 2010
 *
 * @brief Definition and implementation for class @ref MinMaxHeap.
 *
 * @see Paper Introducing the Min-Max Heap (http://www.cs.otago.ac.nz/staffpriv/mike/Papers/MinMaxHeaps/MinMaxHeaps.pdf)
 * @see Alternative Implementation (http://www.coldbrains.com/code/code/C++/Data_Structures/Min-Max_Heap/MinMaxHeap.C.html,
 *                                  http://www.coldbrains.com/code/code/C++/Data_Structures/Min-Max_Heap/MinMaxHeap.H.html)
 **/
#pragma once

#include <algorithm>
#include <stdexcept>
#include <functional>

namespace
{
  inline unsigned log2(unsigned val)
  {
    unsigned int result = 0;
    while (val)
    {
      val >>= 1;
      ++result;
    }
    return result - 1;
  }
}

template<class T, class Compare = std::less<T>>
class MinMaxHeap
{
  T* m_heap;
  Compare m_compare;
  unsigned m_heapsize;
  unsigned m_currentheapsize;

  static inline unsigned int parent(unsigned int index)
  {
    return (index - 1) / 2;
  }

  static inline unsigned int leftChild(unsigned int index)
  {
    return 2 * index + 1;
  }

  static inline unsigned int rightChild(unsigned int index)
  {
    return 2 * index + 2;
  }

  static inline bool isOnMinLevel(unsigned int index)
  {
    return log2(index + 1) % 2 == 1;
  }

  static inline bool isOnMaxLevel(unsigned int index)
  {
    return !isOnMinLevel(index);
  }

  template<bool MaxLevel>
  int trickleUp_(unsigned int index)
  {
    if (index == 0)
      return index;
    unsigned int index_grandparent = parent(index);
    if (index_grandparent == 0)
      return index;
    index_grandparent = parent(index_grandparent);
    if (m_compare(m_heap[index], m_heap[index_grandparent]) ^ MaxLevel)
    {
      swap(m_heap[index_grandparent], m_heap[index]);
      return trickleUp_<MaxLevel>(index_grandparent);
    }
    return index;
  }

  int trickleUp(unsigned int index)
  {
    if (index == 0)
      return index;
    unsigned int index_parent = parent(index);
    if (isOnMinLevel(index))
    {
      if (m_compare(m_heap[index_parent], m_heap[index]))
      {
        swap(m_heap[index_parent], m_heap[index]);
        return trickleUp_<true>(index_parent);
      }
      else
        return trickleUp_<false>(index);
    }
    else
    {
      if (m_compare(m_heap[index], m_heap[index_parent]))
      {
        swap(m_heap[index_parent], m_heap[index]);
        return trickleUp_<false>(index_parent);
      }
      else
        return trickleUp_<true>(index);
    }
  }

  template<bool MaxLevel>
  void trickleDown_(unsigned int index)
  {
    if (index >= m_currentheapsize)
      throw std::invalid_argument("Element specified by zindex does not exist");

    unsigned int smallestNode = index;
    unsigned int left = leftChild(index);

    if (left < m_currentheapsize && (m_compare(m_heap[left], m_heap[smallestNode]) ^ MaxLevel))
      smallestNode = left;
    if (left + 1 < m_currentheapsize && (m_compare(m_heap[left + 1], m_heap[smallestNode]) ^ MaxLevel))
      smallestNode = left + 1;

    unsigned int leftGrandchild = leftChild(left);
    for (unsigned int i = 0; i < 4 && leftGrandchild + i < (unsigned int)m_currentheapsize; ++i)
      if (m_compare(m_heap[leftGrandchild + i], m_heap[smallestNode]) ^ MaxLevel)
        smallestNode = leftGrandchild + i;

    if (index == smallestNode)
      return;

    swap(m_heap[index], m_heap[smallestNode]);

    if (smallestNode - left > 1)
    {
      if (m_compare(m_heap[parent(smallestNode)], m_heap[smallestNode]) ^ MaxLevel)
        swap(m_heap[parent(smallestNode)], m_heap[smallestNode]);

      trickleDown_<MaxLevel>(smallestNode);
    }
  }

  void trickleDown(unsigned int index)
  {
    if (isOnMinLevel(index))
      trickleDown_<false>(index);
    else
      trickleDown_<true>(index);
  }

  unsigned int findMinIndex() const
  {
    switch (m_currentheapsize)
    {
    case 0:
      throw std::underflow_error("No min element exists because "
                                 "there are no elements in the heap.");
    case 1:
      return 0;
    case 2:
      return 1;
    default:
      return m_compare(m_heap[1], m_heap[2]) ? 1 : 2;
    }
  }

  void deleteElement_(unsigned int index)
  {
    if (index >= (unsigned int)m_currentheapsize)
      throw std::underflow_error("Cannot delete specified element from "
                                "the heap because it does not exist.");
    if (index == m_currentheapsize - 1)
    {
      m_currentheapsize--;
      return;
    }

    swap(m_heap[index], m_heap[m_currentheapsize - 1]);
    m_currentheapsize--;
    trickleDown(index);
  }

public:
  MinMaxHeap(unsigned heapsize) : m_heap(nullptr)
  {
    m_heapsize = heapsize;
    m_heap = new T[m_heapsize];
    m_currentheapsize = 0;
  }

  ~MinMaxHeap()
  {
    if (m_heap != nullptr)
      delete[] m_heap;
  }

  void clear()
  {
    m_currentheapsize = 0;
  }

  bool empty() const
  {
    return m_currentheapsize == 0;
  }

  bool full() const
  {
    return m_heapsize == m_currentheapsize;
  }

  unsigned int size() const
  {
    return (unsigned int)m_currentheapsize;
  }

  T* push(const T & val)
  {
    if (m_currentheapsize < m_heapsize)
    {
      m_heap[m_currentheapsize] = val;
      m_currentheapsize++;
      return m_heap + trickleUp(m_currentheapsize - 1);
    }
    else
      throw std::runtime_error("Minmax heap overflow");
  }

  const T & findMax() const
  {
    if (empty())
      throw std::runtime_error("No max element exists, heap is empty");
    return m_heap[0];
  }

  const T & findMin() const
  {
    return m_heap[findMinIndex()];
  }

  T popMax()
  {
    if (empty())
      throw std::runtime_error("Trying to pop form the empty heap");
    T temp = m_heap[0];
    deleteElement_(0);
    return temp;
  }

  T pop()
  {
    return popMax();
  }

  T popMin()
  {
    if (empty())
      throw std::runtime_error("Trying to pop form the empty heap");
    unsigned int smallest = findMinIndex();
    T temp = m_heap[smallest];
    deleteElement_(smallest);
    return temp;
  }

  void deleteElement(const T* ptr)
  {
    unsigned index = unsigned(ptr - m_heap);
    if(index >= 0 && index < m_heapsize)
      deleteElement_(index);
    else
        throw std::runtime_error("Invelid element pointer to delete");
  }
};
