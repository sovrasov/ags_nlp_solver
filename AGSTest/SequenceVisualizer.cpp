#include "SequenceVisualizer.hpp"
#include "OptimizerFunctionWrappers.hpp"
#include <cmath>
#include <discpp.h>

double Get2dPointDifference(const double *x, const double *y)
{
	return fmax(fabs(x[0] - y[0]), fabs(x[1] - y[1]));
}

void VisualizeSequences(OptimizerFunction **task, OptimizerSearchSequence sequence1,
	OptimizerSearchSequence sequence2, char* fileName)
{
	OptimizerFunction *function = task[3];
	int n = 200, i, j, width = 800, height = 1100;
	double  fpi = 1, step, x_left, y_min;
	double arg[2], zlev;

	double xray[200], yray[200], zmat[200][200];

	x_left = y_min = 0;
	fpi = 1;
	step = fpi / (n - 1);
	for (i = 0; i < n; i++)
	{
		xray[i] = i * step;
		yray[i] = i * step;
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			arg[0] = xray[i] - 0.5; arg[1] = yray[j] - 0.5;
			zmat[i][j] = function->Calculate(arg);
		}

	Dislin g;
	g.winsiz(width, height);
	g.page(2200, 2400);
	g.metafl("png");
	g.setfil(fileName);
	g.scrmod("revers");
	g.disini();
	g.complx();
	//g.wintit("Линии уровня");

	//g.titlin(taskName, 3);

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.axspos(240, 2200);
	g.axslen(1900, 1900);

	g.graf(x_left, x_left + fpi, x_left, fpi / 4, y_min, y_min + fpi, y_min, fpi / 4);

	g.height(30);
	for (i = 0; i < 20; i++)
	{
		//zlev = ((VAGRisFunctionWrapper*)function)->GetMinValue() + 0.1 + i;
		zlev = -1.48 +  i*0.09;
		g.contur(xray, n, yray, n, (double*)zmat, zlev);
	}
	
	auto points = &sequence2;
	if (points != nullptr)
	{
		g.color("white");
		double x[2], y[2], minPointDifference;
		g.hsymbl(24);
		int size = points->GetSize();
		for (int i = 0; i < size; i++)
		{
			points->GetPoint(i, y);

			minPointDifference = 10;
			int seq1Size = sequence1.GetSize();
			for (int j = 0; j < seq1Size; j++)
			{
				sequence1.GetPoint(j, x);
				minPointDifference = fmin(minPointDifference, Get2dPointDifference(x, y));
			}

			if (minPointDifference > 0.005)
				g.rlsymb(3, (float)(y[0] + .5), (float)(y[1] + .5));
				
		}
	}


	auto points2 = &sequence1;
	if (points2 != nullptr)
	{
		g.color("white");
		double y[2];
		g.hsymbl(15);
		int size = points2->GetSize();
		for (int i = 0; i < size; i++)
		{
			points2->GetPoint(i, y);
			g.rlsymb(21, (float)(y[0] + .5), (float)(y[1] + .5));
		}
		g.color("red");
		g.hsymbl(30);
	//	g.rlsymb(21, ((VAGRisFunctionWrapper*)function)->GetMinXCoordinate(), ((VAGRisFunctionWrapper*)function)->GetMinYCoordinate());
	}

	g.height(50);
	g.color("fore");
	g.title();
	g.disfin();
}

void VisualizeSolution(OptimizerTask task, 
	OptimizerSearchSequence sequence, OptimizerSolution solution, char* fileName)
{
	std::shared_ptr<OptimizerFunctionPtr> functions = task.GetTaskFunctions();
	OptimizerFunction *targetFunction =
		functions.get()[task.GetNumberOfRestrictions()].get();
	
	int n = 200, i, j, width = 800, height = 1100;
	double  x_step, y_step, x_left, y_left, x_right, y_right;
	double arg[2], zlev;

	double xray[200], yray[200], zmat[200][200];

	x_left = task.GetSpaceTransformation().GetLeftDomainBound().get()[0];
	y_left = task.GetSpaceTransformation().GetLeftDomainBound().get()[1];
	x_right = task.GetSpaceTransformation().GetRightDomainBound().get()[0];
	y_right = task.GetSpaceTransformation().GetRightDomainBound().get()[1];

	x_step = (x_right - x_left) / (n - 1);
	y_step = (y_right - y_left) / (n - 1);

	for (i = 0; i < n; i++)
	{
		xray[i] = x_left + i * x_step;
		yray[i] = y_left + i * y_step;
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			arg[0] = xray[i]; arg[1] = yray[j];
			zmat[i][j] = targetFunction->Calculate(arg);
		}

	Dislin g;
	g.winsiz(width, height);
	g.page(2200, 2400);
	g.metafl("png");
	g.setfil(fileName);
	g.scrmod("revers");
	g.disini();
	g.complx();

	g.name("X-axis", "x");
	g.name("Y-axis", "y");

	g.axspos(240, 2200);
	g.axslen(1900, 1900);

	g.graf(x_left, x_right, x_left, (x_right - x_left) / 4,
		y_left, y_right, y_left, (y_right - y_left) / 4);

	g.height(30);
	for (i = 0; i < 20; i++)
	{
		zlev = solution.GetOptimumValue() + i*0.1;
		g.contur(xray, n, yray, n, (double*)zmat, zlev);
	}

	for (unsigned k = 0; k < task.GetNumberOfRestrictions(); k++)
	{
		targetFunction =
			functions.get()[k].get();
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				arg[0] = xray[i]; arg[1] = yray[j];
				zmat[i][j] = targetFunction->Calculate(arg);
			}
		g.contur(xray, n, yray, n, (double*)zmat, 0);
	}

	auto points2 = &sequence;
	if (points2 != nullptr)
	{
		g.color("white");
		double y[2];
		g.hsymbl(15);
		int size = points2->GetSize();
		for (int i = 0; i < size; i++)
		{
			points2->GetPoint(i, y);
			g.rlsymb(21, (float)(y[0]), (float)(y[1]));
		}
		g.color("red");
		g.hsymbl(30);
		g.rlsymb(21, solution.GetOptimumPoint().get()[0], solution.GetOptimumPoint().get()[1]);
	}

	g.height(50);
	g.color("fore");
	g.title();
	g.disfin();
}