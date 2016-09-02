# kontsevich_graph_series-cpp
Library for manipulation of Kontsevich graphs and sums and series of them (in C++11).

Features:
- generate Kontsevich graphs,
- manipulate sums and series of graphs with arbitrary coefficients,
- reduce sums and series of graphs modulo skew-symmetry,
- construct a star product from a list of admissible graphs and (possibly undetermined) weights,
- obtain cyclic weight relations for graphs with undetermined weights,
- compute the associator (f★g)★h - f★(g★h) for a graph star product,
- evaluate graph series (such as the star product's associator) at particular Poisson structures (e.g. to obtain weight relations),
- substitute weight relations back into a list of graphs and (possibly undetermined) weights,
- reduce graph series modulo the Jacobi identity and its differential consequences.

Dependencies:
- GiNaC and its dependency CLN, for:
  - `kontsevich_graph_operator.hpp`,
  - `util/continued_fraction.hpp`,
  - the test programs.
