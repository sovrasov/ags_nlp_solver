/*
Copyright (C) 2018 Sovrasov V. - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license.
 * You should have received a copy of the MIT license with
 * this file. If not visit https://opensource.org/licenses/MIT
*/

#include "intervals_queue.hpp"

using namespace ags;

bool SingleIntervalsQueue::CompareByR::operator() (const Interval* i1, const Interval* i2) const
{
  return i1->R < i2->R;
}

SingleIntervalsQueue::SingleIntervalsQueue() {}

bool SingleIntervalsQueue::empty() const
{
  return mQueue.empty();
}

Interval* SingleIntervalsQueue::pop(bool is_local)
{
  NLP_SOLVER_ASSERT(!is_local, "Trying to get local R from the single queue");
  auto ret_val = mQueue.top();
  mQueue.pop();
  return ret_val;
}

void SingleIntervalsQueue::push(Interval* i)
{
  mQueue.push(i);
}

void SingleIntervalsQueue::clear()
{
  mQueue = std::priority_queue<Interval*, std::vector<Interval*>, CompareByR>();
}

/* Dual queue */

DualIntervalsQueue::DualIntervalsQueue(size_t size)
{
  NLP_SOLVER_ASSERT(size > 0, "Zero-sized queue");
  mPGlobalHeap = std::make_shared<MinMaxHeap<QueueElement, _less_global>>(size);
  mPLocalHeap = std::make_shared<MinMaxHeap<QueueElement, _less_local>>(size);
}

bool DualIntervalsQueue::empty() const
{
  return mPGlobalHeap->empty() || mPLocalHeap->empty();
}

Interval* DualIntervalsQueue::pop(bool is_local)
{
  if (!is_local)
  {
    auto element = mPGlobalHeap->popMax();
    if (element.pLinkedElement != NULL)
      mPLocalHeap->deleteElement(element.pLinkedElement);
    return element.pInterval;
  }
  else
  {
    auto element = mPLocalHeap->popMax();
    if (element.pLinkedElement != NULL)
      mPGlobalHeap->deleteElement(element.pLinkedElement);
    return element.pInterval;
  }
}

void DualIntervalsQueue::push(Interval* i)
{
  QueueElement *pGlobalElem = nullptr, *pLocalElem = nullptr;
  if (!mPGlobalHeap->full())
    pGlobalElem = mPGlobalHeap->push(QueueElement(i));
  else
  {
    if(i->R > mPGlobalHeap->findMin().pInterval->R)
    {
      QueueElement tmp = mPGlobalHeap->popMin();
      if (tmp.pLinkedElement != nullptr)
        tmp.pLinkedElement->pLinkedElement = nullptr;
      pGlobalElem = mPGlobalHeap->push(QueueElement(i));
    }
  }
  if (!mPLocalHeap->full())
    pLocalElem = mPLocalHeap->push(QueueElement(i));
  else
  {
    if (i->local_R > mPLocalHeap->findMin().pInterval->local_R)
    {
      QueueElement tmp = mPLocalHeap->popMin();
      if (tmp.pLinkedElement != nullptr)
        tmp.pLinkedElement->pLinkedElement = nullptr;
      pLocalElem = mPLocalHeap->push(QueueElement(i));
    }
  }
  //link elements
  if (pGlobalElem != nullptr && pLocalElem != nullptr)
  {
    pGlobalElem->pLinkedElement = pLocalElem;
    pLocalElem->pLinkedElement = pGlobalElem;
  }
}

void DualIntervalsQueue::clear()
{
  mPGlobalHeap->clear();
  mPLocalHeap->clear();
}
