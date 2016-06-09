#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <star-product-filename>\n";
        return 1;
    }

    // Reading in star product:
    string star_product_filename(argv[1]);
    ifstream star_product_file(star_product_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> star_product;
    KontsevichGraphSum<ex> term;
    size_t order = 0;
    for (string line; getline(star_product_file, line); )
    {
        if (line.length() == 0)
            continue;
        if (line[0] == 'h')
        {
            star_product[order] = term;
            term = KontsevichGraphSum<ex>({ });
            order = stoi(line.substr(2));
        }
        else
        {
            KontsevichGraph graph;
            stringstream ss(line);
            ss >> graph;
            string coefficient_str;
            ss >> coefficient_str;
            ex coefficient = coefficient_reader(coefficient_str);
            term += KontsevichGraphSum<ex>({ { coefficient, graph } });
        }
    }
    star_product[order] = term; // the last one

    star_product.reduce();

    KontsevichGraphSeries<ex> arg = { { 0, { { 1, KontsevichGraph(0, 1, {}) } }} };
    KontsevichGraphSeries<ex> assoc = star_product({ star_product, arg }) - star_product({ arg, star_product });

    assoc.reduce();

    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& term : assoc[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
