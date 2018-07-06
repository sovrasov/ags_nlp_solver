# Global search NLP solver

An implementation of the algorithm AGS to solve constrained nonlinear programming problems with Lipschitzian functions. AGS was introduced by prof. R.G. Strongin (see R. G. Strongin, D. L. Markin, ,Minimization of multiextremal functions under nonconvex constraints, Cybernetics 22(4), 486–493. Translated from Russian. Consultant Bureau. New York, 1986. [[link]][paper]). The method exploits Peano-type curve to reduce dimension of the source bounded multidimensional constrained NLP problem and then solves a univariate one.

AGS is proven to converge to a global optima if all objectives and constraints satisfy Lipschitz condition in a given hyperrectangle, the reliability parameter `r` is large enough and accuracy parameter `eps` is zero.

## Clone & build
- on Linux:
```bash
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake ..
make -j 4
```
- on Windows:
```batch
git clone --recursive https://github.com/sovrasov/glob_search_nlp_solver.git
cd multicriterial-go
mkdir build
cd build
cmake .. -G "NMake Makefiles"
nmake
```
[paper]: https://www.tandfonline.com/doi/abs/10.1080/17442508908833568?journalCode=gssr19
