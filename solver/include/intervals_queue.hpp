/*
Copyright (C) 2018 Sovrasov V. - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license.
 * You should have received a copy of the MIT license with
 * this file. If not visit https://opensource.org/licenses/MIT
*/
#pragma once

#include "data_types.hpp"
#include "minmaxheap.hpp"
#include <queue>

namespace ags
{

class IIntervalsQueue
{
public:
  virtual bool empty() const = 0;
  virtual Interval* pop(bool is_local) = 0;
  virtual void push(Interval*) = 0;
  virtual void clear() = 0;
};

class SingleIntervalsQueue : public IIntervalsQueue
{
private:
  struct CompareByR
  {
    bool operator() (const Interval* i1, const Interval* i2) const;
  };
  std::priority_queue<Interval*, std::vector<Interval*>, CompareByR> mQueue;
public:
  SingleIntervalsQueue();
  virtual bool empty() const override;
  virtual Interval* pop(bool is_local) override;
  virtual void push(Interval*) override;
  virtual void clear() override;
};

class DualIntervalsQueue : public IIntervalsQueue
{
private:
  struct QueueElement
  {
    QueueElement *pLinkedElement;
    Interval* pInterval;
    QueueElement() {}
    QueueElement(Interval* _pValue) :
      pLinkedElement(nullptr), pInterval(_pValue) {}
    QueueElement(double _Key, Interval* _pValue, QueueElement* _pLinkedElement) :
      pLinkedElement(_pLinkedElement), pInterval(_pValue) {}

    friend void swap(QueueElement& arg1, QueueElement& arg2)
    {
      using std::swap;
      if (arg1.pLinkedElement != nullptr && arg2.pLinkedElement != nullptr)
        std::swap(arg1.pLinkedElement->pLinkedElement, arg2.pLinkedElement->pLinkedElement);
      else if (arg1.pLinkedElement != nullptr)
        arg1.pLinkedElement->pLinkedElement = &arg2;
      else if (arg2.pLinkedElement != nullptr)
        arg2.pLinkedElement->pLinkedElement = &arg1;
      std::swap(arg1, arg2);
    }
  };

  struct _less_local
  {
    bool operator()(const QueueElement& _Left, const QueueElement& _Right) const
    {
      return _Left.pInterval->local_R < _Right.pInterval->local_R;
    }
  };

  struct _less_global
  {
    bool operator()(const QueueElement& _Left, const QueueElement& _Right) const
    {
      return _Left.pInterval->R < _Right.pInterval->R;
    }
  };

  MinMaxHeap<QueueElement, _less_global>* mPGlobalHeap;
  MinMaxHeap<QueueElement, _less_local>* mPLocalHeap;

public:
  DualIntervalsQueue();
  virtual bool empty() const override;
  virtual Interval* pop(bool is_local) override;
  virtual void push(Interval*) override;
  virtual void clear() override;
};

}
