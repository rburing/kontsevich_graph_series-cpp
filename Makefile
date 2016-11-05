CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3
LDFLAGS=
GINAC_LDFLAGS=-lcln -lginac
EIGEN_CFLAGS=-I/usr/include/eigen3

.PHONY: all
all: bin bin/kontsevich_graph_tests bin/poisson_evaluate bin/generate_graphs bin/star_product bin/cyclic_weight_relations bin/substitute_relations bin/reduce bin/invert bin/gauge bin/star_product_associator bin/reduce_mod_jacobi bin/skew_symmetrize bin/reduce_mod_permutations bin/weight_integrands

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
