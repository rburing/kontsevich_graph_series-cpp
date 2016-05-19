CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3
LDFLAGS=

.PHONY: all
all: tests/kontsevich_graph_tests tests/star_product_associativity tests/relevant_graphs

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/kontsevich_graph_tests.o:
	$(CC) $(CFLAGS) -c tests/kontsevich_graph_tests.cpp -o tests/kontsevich_graph_tests.o

tests/kontsevich_graph_tests: tests/kontsevich_graph_tests.o kontsevich_graph.o
	$(CC) -o tests/kontsevich_graph_tests tests/kontsevich_graph_tests.o kontsevich_graph.o $(LDFLAGS)

tests/star_product_associativity.o:
	$(CC) $(CFLAGS) -c tests/star_product_associativity.cpp -o tests/star_product_associativity.o

tests/star_product_associativity: tests/star_product_associativity.o kontsevich_graph.o
	$(CC) -o tests/star_product_associativity tests/star_product_associativity.o kontsevich_graph.o -lcln -lginac $(LDFLAGS)

tests/relevant_graphs.o:
	$(CC) $(CFLAGS) -c tests/relevant_graphs.cpp -o tests/relevant_graphs.o

tests/relevant_graphs: tests/relevant_graphs.o kontsevich_graph.o
	$(CC) -o tests/relevant_graphs tests/relevant_graphs.o kontsevich_graph.o $(LDFLAGS)

.PHONY: clean
clean:
	rm -f kontsevich_graph.o
	rm -f tests/kontsevich_graph_tests.o tests/kontsevich_graph_tests
	rm -f tests/star_product_associativity.o tests/star_product_associativity
	rm -f tests/relevant_graphs.o tests/relevant_graphs
