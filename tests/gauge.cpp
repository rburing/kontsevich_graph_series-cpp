#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <star-product-filename> <gauge-series-filename>\n";
        return 1;
    }
    
    parser coefficient_reader;
    // Reading in star product series:
    string star_product_filename(argv[1]);
    ifstream star_product_file(star_product_filename);
    KontsevichGraphSeries<ex> star_product = KontsevichGraphSeries<ex>::from_istream(star_product_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    // Reading in gauge series:
    string gauge_series_filename(argv[2]);
    ifstream gauge_series_file(gauge_series_filename);
    KontsevichGraphSeries<ex> gauge_series = KontsevichGraphSeries<ex>::from_istream(gauge_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    KontsevichGraphSeries<ex> gauged_product = star_product.gauge_transform(gauge_series);

    for (size_t n = 0; n <= gauged_product.precision(); ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& term : gauged_product[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
