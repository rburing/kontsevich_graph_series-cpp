#include "../kontsevich_graph_series.hpp"
#include "../leibniz_graph.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc == 1 || (argc == 2 && string(argv[1]) == "--help"))
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> [optional-arguments]\n\n"
             << "Available optional arguments:\n"
             << "--max-iterations=i     perform at most i iterations.\n"
             << "--max-jac-indegree=k   restricts the number of arrows falling on Jacobiators to be <= k.\n"
             << "--skew-leibniz         skew-symmetrizes each Leibniz graph before subtracting it with an undetermined coefficient.\n"
             << "--leibniz-in=filename  input graph series already contains Leibniz graphs, with encodings in filename.\n"
             << "--leibniz-out=filename store Leibniz graph encodings in filename (default: standard output).\n"
             << "--linsys-out=filename  store linear system in filename (default: standard output).\n"
             << "--coeff-prefix=c       let the coefficients of leibniz graphs be c_n.\n"
             << "--solve                the undetermined variables in the input are added to the linear system to-be-solved.\n"
             << "--interactive          ask whether to continue to the next iteration.\n";
        return 1;
    }

    bool interactive = false;
    bool solve = false;
    size_t max_jac_indegree = numeric_limits<size_t>::max();
    size_t max_iterations = numeric_limits<size_t>::max();
    bool skew_leibniz = false;
    string leibniz_in_filename = "";
    string leibniz_out_filename = "";
    string linsys_out_filename = "";
    string coefficient_prefix = "c";

    // Process arguments
    for (int idx = 2; idx < argc; ++idx)
    {
        string argument = argv[idx];
        size_t equals_pos = argument.find('=');
        if (equals_pos != string::npos)
        {
            string key = argument.substr(0, equals_pos);
            string value = argument.substr(equals_pos+1);
            if (key == "--max-iterations")
                max_iterations = stoi(value);
            else if (key == "--max-jac-indegree")
                max_jac_indegree = stoi(value);
            else if (key == "--coeff-prefix")
                coefficient_prefix = value;
            else if (key == "--leibniz-in")
                leibniz_in_filename = value;
            else if (key == "--leibniz-out")
                leibniz_out_filename = value;
            else if (key == "--linsys-out")
                linsys_out_filename = value;
            else {
                cout << "Unrecognized option: " << argument << "\n";
                return 1;
            }
        }
        else
        {
            if (argument == "--skew-leibniz")
                skew_leibniz = true;
            else if (argument == "--solve")
                solve = true;
            else if (argument == "--interactive")
                interactive = true;
            else {
                cout << "Unrecognized option: " << argument << "\n";
                return 1;
            }
        }
    }

    cout << "Options: "
         << "max-iterations = ";
    if (max_iterations == numeric_limits<size_t>::max())
        cout << "none";
    else
        cout << max_iterations;
    cout << ", max-jac-indegree = ";
    if (max_jac_indegree == numeric_limits<size_t>::max())
        cout << "none";
    else
        cout << max_jac_indegree;
    cout << ", solve = " << (solve ? "yes" : "no")
         << ", skew-leibniz = " << (skew_leibniz ? "yes" : "no")
         << ", leibniz-in = " << (leibniz_in_filename == "" ? "none" : leibniz_in_filename)
         << ", leibniz-out = " << (leibniz_out_filename == "" ? "stdout" : leibniz_out_filename)
         << ", linsys-out = " << (linsys_out_filename == "" ? "stdout" : linsys_out_filename)
         << ", coeff-prefix = " << coefficient_prefix
         << ", interactive = " << (interactive ? "yes" : "no") << "\n";

    // Reading in Leibniz graphs
    parser coefficient_reader;
    map< LeibnizGraph<ex>, ex> leibniz_graphs;

    if (leibniz_in_filename != "")
    {
        ifstream leibniz_in_file(leibniz_in_filename);
        leibniz_graphs = LeibnizGraph<ex>::map_from_istream(leibniz_in_file,
                                                                [&coefficient_reader](string s) -> ex
                                                                {
                                                                    return coefficient_reader(s);
                                                                });
    }
    std::vector<ex> leibniz_coeffs;
    for (auto const& pair : leibniz_graphs)
        leibniz_coeffs.push_back(pair.second);

    // Reading in graph series
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    map<size_t, set< vector<size_t> > > in_degrees;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file,
        [&coefficient_reader](string s) -> ex { return coefficient_reader(s); },
        [&in_degrees](KontsevichGraph graph, size_t order) -> bool
                                   {
                                       in_degrees[order].insert(graph.in_degrees());
                                       return true;
                                   }
    );
    size_t order = graph_series.precision();

    graph_series.reduce_mod_skew();

    vector<symbol> unknowns_list;
    for (auto const& namevar : coefficient_reader.get_syms())
        if (find(leibniz_coeffs.begin(), leibniz_coeffs.end(), namevar.second) == leibniz_coeffs.end())
            unknowns_list.push_back(ex_to<symbol>(namevar.second));

    if (unknowns_list.size() != 0 && !solve)
    {
        cerr << "Input contains unknowns; remove them or specify --solve to solve for them.\n";
        return 1;
    }

    vector<symbol> coefficient_list = unknowns_list;
    coefficient_list.reserve(leibniz_coeffs.size());
    for (auto& leibniz_coeff : leibniz_coeffs)
        coefficient_list.push_back(ex_to<symbol>(leibniz_coeff));

    KontsevichGraphSeries<ex> leibniz_graph_series = graph_series;

    set<KontsevichGraph> processed_graphs;

    size_t counter = 0;
    bool converged = false;
    size_t step = 0;
    while (!converged && ++step <= max_iterations)
    {
        size_t old_counter = counter;
        for (size_t n = 0; n <= order; ++n)
        {
            size_t termcounter = 0;
            for (auto& term : graph_series[n])
            {
                cerr << "\rProcessing term " << ++termcounter << " / " << graph_series[n].size() << ". ";
                if (skew_leibniz)
                    cerr << "Skew-";
                cerr << "Leibniz graphs so far: " << counter;

                KontsevichGraph& graph = term.second;

                // Check if this graph has already been processed
                // TODO: optimization: take minimum of skew symmetrization, if skew_leibniz
                if (find(processed_graphs.begin(), processed_graphs.end(), graph) != processed_graphs.end())
                    continue;

                // Use that the graph series is reduced (graphs are in normal form with sign +1)
                processed_graphs.insert(graph);

                // Subtract Leibniz graphs that yield this Kontsevich graph
                for (LeibnizGraph<ex> leibniz_graph : LeibnizGraph<ex>::those_yielding_kontsevich_graph(graph, skew_leibniz))
                {
                    leibniz_graph.normalize();
                    leibniz_graph.sign(1);
                    if (leibniz_graphs.find(leibniz_graph) != leibniz_graphs.end())
                        continue;
                    symbol coefficient(coefficient_prefix + "_" + to_string(counter));
                    KontsevichGraphSum<ex> graph_sum = leibniz_graph.expansion(coefficient);
                    if (graph_sum.size() != 0)
                    {
                        coefficient_list.push_back(coefficient);
                        leibniz_graphs[leibniz_graph] = coefficient;
                        ++counter;
                        leibniz_graph_series[n] -= graph_sum;
                    }
                }
            }
        }

        cout << "\nNumber of Leibniz graphs: " << leibniz_graphs.size() << "\n";
        cout << "\nNumber of terms (before reducing): " << leibniz_graph_series[order].size() << "\n";

        cout.flush();

        cerr << "\nReducing...\n";
        leibniz_graph_series.reduce_mod_skew();

        cout << "\nNumber of terms (after reducing): " << leibniz_graph_series[order].size() << "\n";

        lst equations;

        ostream* linsys_out_stream = &cout;
        ofstream linsys_out_fstream;
        if (linsys_out_filename != "")
        {
            cout << "Writing linear system to " << linsys_out_filename << "\n";
            linsys_out_fstream.open(linsys_out_filename);
            linsys_out_stream = &linsys_out_fstream;
        }
        for (size_t n = 0; n <= order; ++n)
            for (auto& term : leibniz_graph_series[n])
            {
                (*linsys_out_stream) << term.second.encoding() << "    " << term.first << "==0\n";
                equations.append(term.first);
            }

        size_t rows = equations.nops();
        size_t cols = coefficient_list.size();

        cout << "Got linear system of size " << rows << " x " << cols << ".\n";

        char iterate = 'Y';
        if (interactive)
        {
            cerr << "Next iteration? (Y/N) ";
            cin >> iterate;
        }
        if (iterate != 'Y')
            break;

        graph_series = leibniz_graph_series;

        converged = counter == old_counter;
    }
    if (converged)
        cout << "\nConverged in " << step << " steps.\n\n";
    if (skew_leibniz)
        cout << "Skew-";
    cout << "Leibniz graphs:";
    ostream* leibniz_out_stream = &cout;
    ofstream leibniz_out_fstream;
    if (leibniz_out_filename != "")
    {
        cout << " writing to " << leibniz_out_filename << "\n";
        leibniz_out_fstream.open(leibniz_out_filename);
        leibniz_out_stream = &leibniz_out_fstream;
    }
    else
        cout << "\n";
    for (auto& pair : leibniz_graphs)
    {
        (*leibniz_out_stream) << pair.first.encoding() << "    " << pair.second << "\n";
    }
    leibniz_out_fstream.close();
}
