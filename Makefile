CC=g++
CFLAGS=-std=c++11 -O3 -pedantic -Wall -Wextra # -Werror -I${HOME}/include
LDFLAGS=
GINAC_LDFLAGS=-lcln -lginac # -L${HOME}/lib 
EIGEN_CFLAGS=# -I${HOME}/src/eigen3

.PHONY: all
all: bin \
  bin/kontsevich_graph_tests \
  bin/poisson_evaluate \
  bin/poisson_make_vanish \
  bin/generate_graphs \
  bin/star_product \
  bin/cyclic_weight_relations \
  bin/substitute_relations \
  bin/reduce_mod_skew \
  bin/invert \
  bin/gauge \
  bin/star_product_associator \
  bin/reduce_mod_jacobi \
  bin/symmetrize \
  bin/skew_symmetrize \
  bin/weight_integrands \
  bin/extract_coefficient \
  bin/gerstenhaber_bracket \
  bin/schouten_bracket \
  bin/operator_latex \
  bin/reduce_mod_jacobi_iterative \
  bin/reduce_mod_coboundary \
  bin/leibniz_expand \
  bin/leibniz_reduce

bin:
	mkdir bin

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/%.o: tests/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%: tests/%.o kontsevich_graph.o
	$(CC) -o $@ $< kontsevich_graph.o $(LDFLAGS) $(GINAC_LDFLAGS)

tests/reduce_mod_jacobi.o: tests/reduce_mod_jacobi.cpp
	$(CC) $(CFLAGS) $(EIGEN_CFLAGS) -c tests/reduce_mod_jacobi.cpp -o tests/reduce_mod_jacobi.o

.PHONY: clean
clean:
	rm -f kontsevich_graph.o
	rm -f tests/*.o
	rm -rf bin
