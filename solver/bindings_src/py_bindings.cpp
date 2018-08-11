#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "solver.hpp"

namespace py = pybind11;

void call_python_func(py::function func)
{
  std::vector<double> x = {0, .25};
  py::object result_py = func(x);
  double result = result_py.cast<double>();
  std::cout << result << '\n';
}

class AGSPyWrapper
{
private:
  ags::NLPSolver mSolver;

public:
  AGSPyWrapper() {}

  void SetParameters(const ags::SolverParameters& params)
  {
    mSolver.SetParameters(params);
  }

  void SetProblem(const std::vector<py::function> functions,
                  const std::vector<double> leftBound, const std::vector<double> rightBound)
  {

  }

  std::tuple<std::vector<double>, double, int> Solve()
  {
    ags::Trial result = mSolver.Solve();

  }

  std::vector<unsigned> GetCalculationsStatistics() const
  {
    return mSolver.GetCalculationsStatistics();
  }

  std::vector<double> GetHolderConstantsEstimations() const
  {
    return mSolver.GetHolderConstantsEstimations();
  }
};

PYBIND11_MODULE(ags_solver, m)
{
  m.def("call_python_func", &call_python_func);

  py::class_<ags::SolverParameters>(m, "AGSParameters")
    .def(py::init<>())
    .def_readwrite("eps", &ags::SolverParameters::eps)
    .def_readwrite("stopVal", &ags::SolverParameters::stopVal)
    .def_readwrite("r", &ags::SolverParameters::r)
    .def_readwrite("numPoints", &ags::SolverParameters::numPoints)
    .def_readwrite("itersLimit", &ags::SolverParameters::itersLimit)
    .def_readwrite("evolventDensity", &ags::SolverParameters::evolventDensity)
    .def_readwrite("epsR", &ags::SolverParameters::epsR)
    .def_readwrite("refineSolution", &ags::SolverParameters::refineSolution)
    ;

    py::class_<AGSPyWrapper>(m, "Solver")
      .def(py::init<>())
      .def("GetHolderConstantsEstimations", &AGSPyWrapper::GetHolderConstantsEstimations)
      ;
}
