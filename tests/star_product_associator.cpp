#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
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
    KontsevichGraphSeries<ex> star_product = KontsevichGraphSeries<ex>::from_istream(star_product_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    size_t order = star_product.precision();
    star_product.reduce_mod_skew();

    KontsevichGraphSeries<ex> arg = { { 0, { { 1, KontsevichGraph(0, 1, {}) } }} };
    KontsevichGraphSeries<ex> assoc = star_product({ star_product, arg }) - star_product({ arg, star_product });

    assoc.reduce_mod_skew();

    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& indegree : assoc[n].in_degrees(true))
        {
            cout << "# ";
            for (size_t in : indegree)
                cout << in << " ";
            cout << "\n";
            for (auto& term : assoc[n][indegree])
            {
                cout << term.second.encoding() << "    " << term.first << "\n";
            }
        }
    }
}
