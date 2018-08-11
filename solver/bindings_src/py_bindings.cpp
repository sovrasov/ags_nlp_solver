#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "solver.hpp"

namespace py = pybind11;

class AGSPyWrapper
{
private:
  ags::NLPSolver mSolver;
  int mDimension;

public:
  AGSPyWrapper() : mDimension(0) {}

  void SetParameters(const ags::SolverParameters& params)
  {
    mSolver.SetParameters(params);
  }

  void SetProblem(const std::vector<py::function> py_functions,
                  const std::vector<double> leftBound, const std::vector<double> rightBound)
  {
    std::vector<ags::NLPSolver::FuncPtr> functions;
    size_t dim = leftBound.size();
    for (const auto& py_f : py_functions)
    {
      functions.push_back([py_f, dim](const double* y)
      {
        std::vector<double> y_vec(y, y + dim);
        py::object result_py = py_f(y_vec);
        return result_py.cast<double>();
      }
      );
    }
    mDimension = static_cast<int>(dim);
    mSolver.SetProblem(functions, leftBound, rightBound);
  }

  std::tuple<std::vector<double>, double, int> Solve()
  {
    ags::Trial result = mSolver.Solve();
    std::vector<double> point(result.y, result.y + mDimension);
    return std::make_tuple(point, result.g[result.idx], result.idx);
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
  py::class_<ags::SolverParameters>(m, "Parameters")
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
      .def("GetCalculationsStatistics", &AGSPyWrapper::GetCalculationsStatistics)
      .def("SetParameters", &AGSPyWrapper::SetParameters)
      .def("Solve", &AGSPyWrapper::Solve)
      .def("SetProblem", &AGSPyWrapper::SetProblem)
      ;
}
