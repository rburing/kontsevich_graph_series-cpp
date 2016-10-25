CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3
LDFLAGS=
EIGEN_CFLAGS=-I/usr/include/eigen3

.PHONY: all
all: bin bin/kontsevich_graph_tests bin/poisson_evaluate bin/generate_graphs bin/star_product bin/cyclic_weight_relations bin/substitute_relations bin/reduce bin/invert bin/gauge bin/star_product_associator bin/reduce_mod_jacobi bin/skew_symmetrize bin/reduce_mod_permutations

bin:
	mkdir bin

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/kontsevich_graph_tests.o:
	$(CC) $(CFLAGS) -c tests/kontsevich_graph_tests.cpp -o tests/kontsevich_graph_tests.o

bin/kontsevich_graph_tests: bin tests/kontsevich_graph_tests.o kontsevich_graph.o
	$(CC) -o bin/kontsevich_graph_tests tests/kontsevich_graph_tests.o kontsevich_graph.o $(LDFLAGS)

tests/poisson_evaluate.o:
	$(CC) $(CFLAGS) -c tests/poisson_evaluate.cpp -o tests/poisson_evaluate.o

bin/poisson_evaluate: bin tests/poisson_evaluate.o kontsevich_graph.o
	$(CC) -o bin/poisson_evaluate tests/poisson_evaluate.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/generate_graphs.o:
	$(CC) $(CFLAGS) -c tests/generate_graphs.cpp -o tests/generate_graphs.o

bin/generate_graphs: bin tests/generate_graphs.o kontsevich_graph.o
	$(CC) -o bin/generate_graphs tests/generate_graphs.o kontsevich_graph.o $(LDFLAGS)

tests/star_product.o:
	$(CC) $(CFLAGS) -c tests/star_product.cpp -o tests/star_product.o

bin/star_product: bin tests/star_product.o kontsevich_graph.o
	$(CC) -o bin/star_product tests/star_product.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/cyclic_weight_relations.o:
	$(CC) $(CFLAGS) -c tests/cyclic_weight_relations.cpp -o tests/cyclic_weight_relations.o

bin/cyclic_weight_relations: bin tests/cyclic_weight_relations.o kontsevich_graph.o
	$(CC) -o bin/cyclic_weight_relations tests/cyclic_weight_relations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/substitute_relations.o:
	$(CC) $(CFLAGS) -c tests/substitute_relations.cpp -o tests/substitute_relations.o

bin/substitute_relations: bin tests/substitute_relations.o kontsevich_graph.o
	$(CC) -o bin/substitute_relations tests/substitute_relations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce.o:
	$(CC) $(CFLAGS) -c tests/reduce.cpp -o tests/reduce.o

bin/reduce: bin tests/reduce.o kontsevich_graph.o
	$(CC) -o bin/reduce tests/reduce.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/invert.o:
	$(CC) $(CFLAGS) -c tests/invert.cpp -o tests/invert.o

bin/invert: bin tests/invert.o kontsevich_graph.o
	$(CC) -o bin/invert tests/invert.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/gauge.o:
	$(CC) $(CFLAGS) -c tests/gauge.cpp -o tests/gauge.o

bin/gauge: bin tests/gauge.o kontsevich_graph.o
	$(CC) -o bin/gauge tests/gauge.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/star_product_associator.o:
	$(CC) $(CFLAGS) -c tests/star_product_associator.cpp -o tests/star_product_associator.o

bin/star_product_associator: bin tests/star_product_associator.o kontsevich_graph.o
	$(CC) -o bin/star_product_associator tests/star_product_associator.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce_mod_jacobi.o:
	$(CC) $(CFLAGS) $(EIGEN_CFLAGS) -c tests/reduce_mod_jacobi.cpp -o tests/reduce_mod_jacobi.o

bin/reduce_mod_jacobi: bin tests/reduce_mod_jacobi.o kontsevich_graph.o
	$(CC) -o bin/reduce_mod_jacobi tests/reduce_mod_jacobi.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/skew_symmetrize.o:
	$(CC) $(CFLAGS) -c tests/skew_symmetrize.cpp -o tests/skew_symmetrize.o

bin/skew_symmetrize: bin tests/skew_symmetrize.o kontsevich_graph.o
	$(CC) -o bin/skew_symmetrize tests/skew_symmetrize.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce_mod_permutations.o:
	$(CC) $(CFLAGS) -c tests/reduce_mod_permutations.cpp -o tests/reduce_mod_permutations.o

bin/reduce_mod_permutations: bin tests/reduce_mod_permutations.o kontsevich_graph.o
	$(CC) -o bin/reduce_mod_permutations tests/reduce_mod_permutations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

.PHONY: clean
clean:
	rm -f kontsevich_graph.o
	rm -f tests/*.o
	rm -rf bin
