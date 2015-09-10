#define _USE_MATH_DEFINES
#include "TestTaskFactory.hpp"
#include "OptimizerSTLFunctionWrapper.hpp"
#include <cmath>

using namespace optimizercore;

OptimizerTask TestTaskFactory::GetTask(int taskNum)
{
	OptimizerFunctionPtr *functions;
	OptimizerTask task;
	SharedVector leftBound;
	SharedVector rightBound;

	switch (taskNum)
	{
	case 0:
		leftBound = SharedVector(new double[2]);
		rightBound = SharedVector(new double[2]);
		leftBound.get()[0] = 0;
		leftBound.get()[1] = -1;
		rightBound.get()[0] = 4;
		rightBound.get()[1] = 3;


		functions = new OptimizerFunctionPtr[4];

		functions[0] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 0.01*(pow(x[0] - 2.2, 2) + pow(x[1] - 1.2, 2) - 2.25);
		}
		)));

		functions[1] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 100 * (1 - pow(x[0] - 2, 2) / 1.44 - pow(0.5*x[1], 2));
		}
		)));

		functions[2] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 10 * (x[1] - 1.5 - 1.5*sin(2*M_PI*(x[0] - 1.75)));
		}
		)));
			
		functions[3] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return -1.5*pow(x[0], 2) * exp(1 - pow(x[0], 2)
				- 20.25*pow(x[0] - x[1], 2)) - pow(0.5 * (x[1] - 1)*(x[0]- 1), 4)
				* exp(2 - pow(0.5 * (x[0] - 1), 4) - pow(x[1] - 1, 4));
		}
		)));
		
		task = OptimizerTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
			utils::array_deleter<OptimizerFunctionPtr>()), 3, 2, leftBound, rightBound);

		break;

	case 1:
		leftBound = SharedVector(new double[2]);
		rightBound = SharedVector(new double[2]);
		leftBound.get()[0] = 0.5;
		leftBound.get()[1] = -.5;
		rightBound.get()[0] = .5;
		rightBound.get()[1] = .5;

		functions = new OptimizerFunctionPtr[2];

		functions[1] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return pow(x[0], 2) + pow(x[1], 2);
		}
		)));
		functions[0] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return x[0] + x[1] + .1;
		}
		)));

		task = OptimizerTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
			utils::array_deleter<OptimizerFunctionPtr>()), 1, 2, leftBound, rightBound);

		break;

	case 2:
		leftBound = SharedVector(new double[2]);
		rightBound = SharedVector(new double[2]);
		leftBound.get()[0] = 0;
		leftBound.get()[1] = -1;
		rightBound.get()[0] = 4;
		rightBound.get()[1] = 3;

		functions = new OptimizerFunctionPtr[3];

		functions[0] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return pow(x[0] - 2.2, 2) + pow(x[1] - 1.2, 2) - 1.25;
		}
		)));

		functions[1] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 1.21 - (pow(x[0] - 2.2, 2) + pow(x[1] - 1.2, 2));
		}
		)));

		functions[2] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return -1.5*pow(x[0], 2) * exp(1 - pow(x[0], 2)
				- 20.25*pow(x[0] - x[1], 2)) - pow(0.5 * (x[1] - 1)*(x[0] - 1), 4)
				* exp(2 - pow(0.5 * (x[0] - 1), 4) - pow(x[1] - 1, 4));
		}
		)));

		task = OptimizerTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
			utils::array_deleter<OptimizerFunctionPtr>()), 2, 2, leftBound, rightBound);
		break;

	case 3:
		leftBound = SharedVector(new double[2]);
		rightBound = SharedVector(new double[2]);
		leftBound.get()[0] = 0;
		leftBound.get()[1] = 0;
		rightBound.get()[0] = 80;
		rightBound.get()[1] = 80;

		functions = new OptimizerFunctionPtr[5];
		functions[0] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* y)->double
		{
			return 450 - y[0] * y[1];
		}
		)));
		functions[1] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* y)->double
		{
			return (0.1*y[0] - 1)*(0.1*y[0] - 1) - y[1];
		}
		)));
		functions[2] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* y)->double
		{
			return 8 * (y[0] - 40) - (y[1] - 30)*(y[1] - 55);
		}
		)));
		functions[3] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* y)->double
		{
			return y[1] + (y[0] - 35)*(y[0] - 30) / 125 - 80;
		}
		)));
		functions[4] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* y)->double
		{
			static const double Coeff[20] = { 75.1963666677, -3.8112755343, 0.1269366345, -0.0020567665, 0.000010345,
				-6.8306567631, 0.0302344793, -0.0012813448, 0.0000352559, -0.0000002266,
				0.2564581253, -0.0034604030, 0.0000135139, -28.1064434908, -0.0000052375,
				-0.0000000063, 0.0000000007, 0.0003405462, -0.0000016638, -2.8673112392 };
			return -(Coeff[0] + Coeff[1] * y[0] + Coeff[2] * y[0] * y[0] + Coeff[3] * pow(y[0], 3) + Coeff[4] * pow(y[0], 4) +
				Coeff[5] * y[1] + Coeff[6] * y[0] * y[1] + Coeff[7] * y[0] * y[0] * y[1] + Coeff[8] * pow(y[0], 3)*y[1] + Coeff[9] * pow(y[0], 4)*y[1] +
				Coeff[10] * y[1] * y[1] + Coeff[11] * pow(y[1], 3) + Coeff[12] * pow(y[1], 4) + Coeff[13] / (y[1] + 1) + Coeff[14] * y[0] * y[0] * y[1] * y[1] +
				Coeff[15] * pow(y[0], 3)*y[1] * y[1] + Coeff[16] * pow(y[0], 3)*pow(y[1], 3) + Coeff[17] * y[0] * y[1] * y[1] + Coeff[18] * y[0] * pow(y[1], 3) + Coeff[19] * exp(0.0005*y[0] * y[1]));;
		}
		)));

		task = OptimizerTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
			utils::array_deleter<OptimizerFunctionPtr>()), 4, 2, leftBound, rightBound);
		break;

	case 4:
		leftBound = SharedVector(new double[5]);
		rightBound = SharedVector(new double[5]);
		leftBound.get()[0] = leftBound.get()[1] = leftBound.get()[2] = -3;
		rightBound.get()[0] = rightBound.get()[1] = rightBound.get()[2] = 3;

		leftBound.get()[3] = leftBound.get()[4] = -10;
		rightBound.get()[3] = rightBound.get()[4] = 10;

		functions = new OptimizerFunctionPtr[6];
		functions[0] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return -(x[0] + x[1] + x[2] + x[3] + x[4]);
		}
		)));
		functions[1] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return x[1] * x[1] / 9 + x[3] * x[3]/100 - 1.4;
		}
		)));
		functions[2] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 3 - pow(x[0] + 1, 2) - pow(x[1] + 2, 2) - pow(x[2] - 2, 2) - pow(x[4] + 5, 2);
		}
		)));
		functions[3] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return 4 * x[0] * x[0] * sin(x[0]) + x[1] * x[1] * cos(x[1] + x[3]) +
				x[2] * x[2] * (sin(x[2] + x[4]) + sin(10 * (x[2] - x[3]) / 3)) - 4;
		}
		)));
		functions[4] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return x[0] * x[0] + 
				x[1] * x[1] * pow(sin((x[0] + x[3]) / 3 + 6.6) + sin((x[1] + x[4]) / 2 + 0.9), 2) -
				17 * pow(cos(x[0] + x[2] + 1), 2) + 16;
		}
		)));
		functions[5] = OptimizerFunctionPtr(new OptimizerSTLFunctionWrapper(
			std::function<double(const double *)>(
			[](const double* x)->double
		{
			return sin(x[0] * x[2]) - (x[1] * x[4] + x[2] * x[3])*cos(x[0] * x[1]);
		}	
		)));

		task = OptimizerTask(std::shared_ptr<OptimizerFunctionPtr>(functions,
			utils::array_deleter<OptimizerFunctionPtr>()), 5, 5, leftBound, rightBound);
		break;
	}

	return task;
}