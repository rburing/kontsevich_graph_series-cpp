CC=g++
CFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -O3

all: tests/kontsevich_graph_tests

kontsevich_graph.o:
	$(CC) $(CFLAGS) -c kontsevich_graph.cpp

tests/kontsevich_graph_tests.o:
	$(CC) $(CFLAGS) -c tests/kontsevich_graph_tests.cpp -o tests/kontsevich_graph_tests.o

tests/kontsevich_graph_tests: tests/kontsevich_graph_tests.o kontsevich_graph.o
	$(CC) -o tests/kontsevich_graph_tests tests/kontsevich_graph_tests.o kontsevich_graph.o

clean:
	rm -f kontsevich_graph.o tests/kontsevich_graph_tests.o
	rm -f tests/kontsevich_graph_tests
