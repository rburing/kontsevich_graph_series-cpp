CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3

all: tests/kontsevich_graph_tests tests/star_product_associativity

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/kontsevich_graph_tests.o:
	$(CC) $(CFLAGS) -c tests/kontsevich_graph_tests.cpp -o tests/kontsevich_graph_tests.o

tests/kontsevich_graph_tests: tests/kontsevich_graph_tests.o kontsevich_graph.o
	$(CC) -o tests/kontsevich_graph_tests tests/kontsevich_graph_tests.o kontsevich_graph.o

tests/star_product_associativity.o:
	$(CC) $(CFLAGS) -c tests/star_product_associativity.cpp -o tests/star_product_associativity.o

tests/star_product_associativity: tests/star_product_associativity.o kontsevich_graph.o
	$(CC) -o tests/star_product_associativity tests/star_product_associativity.o kontsevich_graph.o -lcln -lginac

clean:
	rm -f kontsevich_graph.o tests/kontsevich_graph_tests.o tests/star_product_associativity.o
	rm -f tests/kontsevich_graph_tests tests/star_product_associativity
