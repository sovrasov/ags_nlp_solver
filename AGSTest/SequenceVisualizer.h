#include "discpp.h"
#include "OptimizerFunction.hpp"
#include "OptimizerSearchSequence.hpp"
#include "OptimizerTask.hpp"
#include "OptimizerSolution.hpp"

using namespace optimizercore;

void VisualizeSequences(OptimizerFunction **functions, OptimizerSearchSequence sequence1,
	OptimizerSearchSequence sequence2, char* fileName);
void VisualizeSolution(OptimizerTask task, 
	OptimizerSearchSequence sequence, OptimizerSolution solution, char* fileName);