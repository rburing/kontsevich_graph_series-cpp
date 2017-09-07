#include "../kontsevich_graph_series.hpp"
#include "../util/partitions.hpp"
#include "../util/factorial.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename>\n\n"
             << "The file should contain a graph series with n-internal-vertex prime graphs modulo mirror images at order n.\n";
        return 1;
    }
    string filename(argv[1]);
    map< size_t, set<KontsevichGraph> > relevants;
    map<KontsevichGraph, ex> weights;
    // Reading in known graphs and their (possibly symbolic) weights:
    ifstream weights_file(filename);
    symtab weights_table;
    parser weights_reader(weights_table);
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(weights_file, [&weights_reader](std::string s) -> ex { return weights_reader(s); });
    for (auto& order : graph_series)
    {
        for (auto& term : order.second)
        {
            weights[term.second] = term.first;
            relevants[order.first].insert(term.second);
        }
    }
    if (relevants.empty())
    {
        cout << "Found no graphs with weights. Aborting mission.\n";
        return 1;
    }
    size_t order = relevants.rbegin()->first;
    // Making a table of primes and weights:
    map< size_t, vector<KontsevichGraph> > primes;
    for (size_t n = 1; n <= order; ++n)
    {
        primes[n].reserve(2*relevants[n].size());
        for (KontsevichGraph g : relevants[n])
        {
            primes[n].push_back(g);
            KontsevichGraph mirror = g.mirror_image();
            if (g.abs() != mirror.abs())
            {
                weights[mirror] = weights[g];
                if (n % 2 == 1)
                    weights[mirror] *= -1;
                primes[n].push_back(mirror);
            }
        }
    }
    // Constructing star product:
    KontsevichGraphSeries<ex> star_product;
    star_product[0] = KontsevichGraphSum<ex> { { 1, KontsevichGraph(0, 2, {}) } };
    for (size_t n = 1; n <= order; ++n)
    {
        ex major_coeff = 1/(ex)factorial(n);
        Partitions partitions(n);
        for (auto partition = partitions.begin(); partition != partitions.end(); ++partition)
        {
            // count multiplicities of parts
            std::map<size_t, std::vector<size_t> > multiplicity;
            for (size_t i = 0; i != (*partition).size(); ++i)
                multiplicity[(*partition)[i]].push_back(i);
            // ignore multiplicity 1
            for (auto part : *partition)
                if (multiplicity[part].size() == 1)
                    multiplicity.erase(part);

            std::vector<size_t> prime_sizes((*partition).size());
            for (size_t i = 0; i != (*partition).size(); ++i)
                prime_sizes[i] = primes[(*partition)[i]].size();
            CartesianProduct decompositions(prime_sizes);
            for (auto decomposition = decompositions.begin(); decomposition != decompositions.end(); ++decomposition)
            {
                bool accept = true;
                for (auto& part : multiplicity)
                {
                    size_t prev = 0;
                    for (auto& idx : part.second) // when drawing from multiple equal sets
                    {
                        size_t current = (*decomposition)[idx];
                        if (current < prev) // accept only non-decreasing indices, to avoid duplicates
                            accept = false;
                        prev = current;
                    }
                }
                if (!accept)
                    continue;

                KontsevichGraph composite(0, 2, {});
                ex coeff = major_coeff;
                for (size_t i = 0; i != (*partition).size(); ++i)
                {
                    composite *= primes[(*partition)[i]][(*decomposition)[i]];
                    coeff *= weights[primes[(*partition)[i]][(*decomposition)[i]]];
                }
                coeff *= composite.multiplicity();
                star_product[n] += KontsevichGraphSum<ex>({ { coeff, composite} });
            }
        }
        // TODO: n is not yet considered a partition of n, so treat it separately for now:
        for (KontsevichGraph prime : primes[n])
        {
            star_product[n] += KontsevichGraphSum<ex>({ { major_coeff * weights[prime] * prime.multiplicity(), prime } });
        }
    }
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& indegree : star_product[n].in_degrees(true))
        {
            cout << "# ";
            for (size_t in : indegree)
                cout << in << " ";
            cout << "\n";
            for (auto& term : star_product[n][indegree])
            {
                cout << term.second.encoding() << "    " << term.first << "\n";
            }
        }
    }
}
