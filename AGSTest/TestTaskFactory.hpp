#ifndef TEST_TASK_FACTORY_HPP
#define TEST_TASK_FACTORY_HPP

#include "OptimizerTask.hpp"

class TestTaskFactory
{

public:

	static optimizercore::OptimizerTask GetTask(int taskNum);

};

#endif