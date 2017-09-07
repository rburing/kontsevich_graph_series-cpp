# kontsevich_graph_series-cpp
Library for manipulation of Kontsevich graphs and sums and series of them (in C++11).

General features:
- generate Kontsevich graphs,
- manipulate sums and series of graphs with arbitrary (e.g. numeric, symbolic) coefficients,
- reduce sums and series of graphs modulo skew-symmetry,
- reduce graph series modulo the Jacobi identity and its differential consequences,
- evaluate graph series at particular Poisson structures,
- for a graph series that should vanish for all Poisson structures,
  obtain relations between the coefficients by substituting particular Poisson structures,
- substitute relations between coefficients back into a graph series,
- write the polydifferential operator associated to a graph series as a LaTeX formula.

Star product features:
- construct a star product from a list of admissible graphs and (possibly undetermined) weights,
- obtain cyclic weight relations for graphs with undetermined weights,
- compute the associator (f★g)★h - f★(g★h) for a graph star product,
- obtain weight relations from the associativity constraint for particular Poisson structures,
- gauge-transform a star product,
- invert a gauge transformation,
- calculate the weight integrand of a graph as a rational function of Cartesian coordinates.

Poisson cohomology features:
- skew-symmetrize graph series,
- calculate the Schouten bracket of polyvectorfields given by skew Kontsevich graph series.

Dependencies:
- GiNaC and its dependency CLN, for:
  - `kontsevich_graph_operator.hpp`,
  - `kontsevich_graph_weight.hpp`,
  - `util/continued_fraction.hpp`,
  - `util/poisson_structure.hpp`,
  - `util/poisson_structure_examples.hpp`,
  - and the test programs.
- Eigen, for linear algebra in
  - `tests/reduce_mod_jacobi.cpp`,
  - `tests/reduce_mod_jacobi_iterative.cpp`.
