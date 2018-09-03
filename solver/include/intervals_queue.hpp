/*
Copyright (C) 2018 Sovrasov V. - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license.
 * You should have received a copy of the MIT license with
 * this file. If not visit https://opensource.org/licenses/MIT
*/
#pragma once

#include "data_types.hpp"
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
  std::priority_queue<Interval*, std::vector<Interval*>, CompareByR> mQueue;
public:
  SingleIntervalsQueue();
  virtual bool empty() const override;
  virtual Interval* pop(bool is_local) override;
  virtual void push(Interval*) override;
  virtual void clear() override;
};

}
