#include "../kontsevich_graph_series.hpp"
#include "../kontsevich_graph_operator.hpp"
#include "../util/poisson_structure.hpp" // for poisson_structures
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 3 || poisson_structures.find(argv[2]) == poisson_structures.end())
    {
        cerr << "Usage: " << argv[0] << " <graph-series-filename> <poisson-structure>\n\n"
             << "Poisson structures can be chosen from the following list:\n";
        for (auto const& entry : poisson_structures)
        {
            cerr << "- " << entry.first << "\n";
        }
        return 1;
    }

    PoissonStructure const& poisson = poisson_structures[argv[2]];

    cout << "Coordinates: ";
    for (symbol coordinate : poisson.coordinates)
        cout << coordinate << " ";
    cout << "\n";
    cout << "Poisson structure matrix:\n";
    cout << "[";
    for (auto row = poisson.bivector.begin(); row != poisson.bivector.end(); ++row)
    {
        cout << "[";
        for (auto entry = row->begin(); entry != row->end(); ++entry)
        {
            cout << *entry;
            if (entry + 1 != row->end())
                cout << ", ";
        }
        cout << "]";
        if (row + 1 != poisson.bivector.end())
            cout << "\n";
    }
    cout << "]\n\n";

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    size_t order = graph_series.precision();

    for (size_t n = 0; n <= order; ++n)
    {
        if (graph_series[n] != 0 || n == order)
            cout << "h^" << n << ":\n";
        for (std::vector<size_t> indegrees : graph_series[n].in_degrees(true))
        {
            cout << "# ";
            for (size_t j = 0; j != indegrees.size(); ++j)
                cout << indegrees[j] << " ";
            cout << "\n";
            
            map< multi_indexes, ex > coefficients;
            for (auto& term : graph_series[n][indegrees])
            {
                map_operator_coefficients_from_graph(term.second, poisson, [&coefficients, &term](multi_indexes arg_derivatives, GiNaC::ex summand) {
                    ex result = (term.first * summand).expand();
                    coefficients[arg_derivatives] += result;
                });
            }
            for (auto& entry : coefficients)
            {
                cout << "# ";
                for (auto mindex = entry.first.begin(); mindex != entry.first.end(); mindex++)
                {
                    cout << "[ ";
                    for (auto &index : *mindex)
                    {
                        cout << poisson.coordinates[index] << " ";
                    }
                    cout << "]";
                    if (mindex + 1 != entry.first.end())
                        cout << " ";
                }
                cout << "\n";
                ex result = entry.second;
                cout << result << "\n";
            }
            cout.flush();
        }
    }
}
