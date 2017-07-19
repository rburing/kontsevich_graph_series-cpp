CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3
LDFLAGS=
GINAC_LDFLAGS=-lcln -lginac
EIGEN_CFLAGS=-I/usr/include/eigen3

.PHONY: all
all: bin bin/kontsevich_graph_tests bin/poisson_evaluate bin/poisson_make_vanish bin/generate_graphs bin/star_product bin/cyclic_weight_relations bin/substitute_relations bin/reduce bin/invert bin/gauge bin/star_product_associator bin/reduce_mod_jacobi bin/skew_symmetrize bin/weight_integrands bin/extract_coefficient bin/schouten_bracket bin/operator_latex bin/reduce_mod_jacobi_iterative

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

tests/reduce_mod_jacobi_iterative.o: tests/reduce_mod_jacobi_iterative.cpp
	$(CC) $(CFLAGS) $(EIGEN_CFLAGS) -c tests/reduce_mod_jacobi_iterative.cpp -o tests/reduce_mod_jacobi_iterative.o

.PHONY: clean
clean:
	rm -f kontsevich_graph.o
	rm -f tests/*.o
	rm -rf bin
