CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3
LDFLAGS=

.PHONY: all
all: bin bin/kontsevich_graph_tests bin/star_product_associativity bin/relevant_graphs bin/star_product bin/cyclic_weight_relations bin/substitute_weight_relations bin/reduce bin/star_product_associator bin/reduce_mod_jacobi bin/reduce_mod_permutations

bin:
	mkdir bin

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/kontsevich_graph_tests.o:
	$(CC) $(CFLAGS) -c tests/kontsevich_graph_tests.cpp -o tests/kontsevich_graph_tests.o

bin/kontsevich_graph_tests: bin tests/kontsevich_graph_tests.o kontsevich_graph.o
	$(CC) -o bin/kontsevich_graph_tests tests/kontsevich_graph_tests.o kontsevich_graph.o $(LDFLAGS)

tests/star_product_associativity.o:
	$(CC) $(CFLAGS) -c tests/star_product_associativity.cpp -o tests/star_product_associativity.o

bin/star_product_associativity: bin tests/star_product_associativity.o kontsevich_graph.o
	$(CC) -o bin/star_product_associativity tests/star_product_associativity.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/relevant_graphs.o:
	$(CC) $(CFLAGS) -c tests/relevant_graphs.cpp -o tests/relevant_graphs.o

bin/relevant_graphs: bin tests/relevant_graphs.o kontsevich_graph.o
	$(CC) -o bin/relevant_graphs tests/relevant_graphs.o kontsevich_graph.o $(LDFLAGS)

tests/star_product.o:
	$(CC) $(CFLAGS) -c tests/star_product.cpp -o tests/star_product.o

bin/star_product: bin tests/star_product.o kontsevich_graph.o
	$(CC) -o bin/star_product tests/star_product.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/cyclic_weight_relations.o:
	$(CC) $(CFLAGS) -c tests/cyclic_weight_relations.cpp -o tests/cyclic_weight_relations.o

bin/cyclic_weight_relations: bin tests/cyclic_weight_relations.o kontsevich_graph.o
	$(CC) -o bin/cyclic_weight_relations tests/cyclic_weight_relations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/substitute_weight_relations.o:
	$(CC) $(CFLAGS) -c tests/substitute_weight_relations.cpp -o tests/substitute_weight_relations.o

bin/substitute_weight_relations: bin tests/substitute_weight_relations.o kontsevich_graph.o
	$(CC) -o bin/substitute_weight_relations tests/substitute_weight_relations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce.o:
	$(CC) $(CFLAGS) -c tests/reduce.cpp -o tests/reduce.o

bin/reduce: bin tests/reduce.o kontsevich_graph.o
	$(CC) -o bin/reduce tests/reduce.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/star_product_associator.o:
	$(CC) $(CFLAGS) -c tests/star_product_associator.cpp -o tests/star_product_associator.o

bin/star_product_associator: bin tests/star_product_associator.o kontsevich_graph.o
	$(CC) -o bin/star_product_associator tests/star_product_associator.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce_mod_jacobi.o:
	$(CC) $(CFLAGS) -c tests/reduce_mod_jacobi.cpp -o tests/reduce_mod_jacobi.o

bin/reduce_mod_jacobi: bin tests/reduce_mod_jacobi.o kontsevich_graph.o
	$(CC) -o bin/reduce_mod_jacobi tests/reduce_mod_jacobi.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/reduce_mod_permutations.o:
	$(CC) $(CFLAGS) -c tests/reduce_mod_permutations.cpp -o tests/reduce_mod_permutations.o

bin/reduce_mod_permutations: bin tests/reduce_mod_permutations.o kontsevich_graph.o
	$(CC) -o bin/reduce_mod_permutations tests/reduce_mod_permutations.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

.PHONY: clean
clean:
	rm -f kontsevich_graph.o
	rm -f tests/*.o
	rm -rf bin
