#include "../kontsevich_graph_series.hpp"
#include <string>
#include <functional>

using namespace std;

int main(int argc, char* argv[])
{
    std::function<int(void)> usage = [argv]() {
        cerr << "Usage: " << argv[0] << " <internal> [external] [graph-options] [additional-options]\n\n"
             << "Available graph-options: --prime --zero --positive-differential-order\n"
             << "Possible graph-option values: yes, no.\n"
             << "Omitting a graph-option indicates indifference.\n"
             << "Additional options: --with-coefficients (default: no), --modulo-mirror-images (default: no), --normal-forms (default: no)\n"
             << "Example: " << argv[0] << " 3 --prime=yes --normal-forms=yes --modulo-mirror-images=yes\n";
        return 1;
    };

    if (argc < 2)
    {
        return usage();
    }
    size_t internal;
    try {
        internal = stoi(argv[1]);
    }
    catch (std::invalid_argument) {
        return usage();
    }
    size_t external = 2;
    int optional_argument_start = 3;
    if (argc > 2)
    {
        try {
            external = stoi(argv[2]);
        }
        catch (std::invalid_argument)
        {
            optional_argument_start = 2;
        }
    }

    // TODO: factor out parsing of optional arguments

    enum class Answer { Yes, No, Indifferent };
    map< string, Answer > answer_name = { { "yes",         Answer::Yes },
                                          { "no",          Answer::No },
                                          { "indifferent", Answer::Indifferent } };

    map< string, Answer > option_values = { { "prime",                       Answer::Indifferent },
                                            { "zero",                        Answer::Indifferent },
                                            { "positive-differential-order", Answer::Indifferent },
                                            { "with-coefficients",           Answer::No },
                                            { "normal-forms",                Answer::No },
                                            { "modulo-mirror-images",        Answer::No } };

    if (optional_argument_start < argc)
    {
        for (int idx = optional_argument_start; idx != argc; ++idx)
        {
            string argument = argv[idx];
            if (argument.length() < 2 || argument.substr(0,2) != "--" || argument.find("=") == string::npos)
            {
                cerr << "Invalid argument: " << argument << "\n";
                return 1;
            }
            string option = argument.substr(2, argument.find("=") - 2);
            string answer = argument.substr(argument.find("=") + 1, string::npos);
            try {
                option_values.at(option) = answer_name.at(answer);
            }
            catch (std::out_of_range)
            {
                cerr << "Invalid argument: " << argument << "\n";
                return 1;
            }
        }
    }

    size_t counter = 0;
    KontsevichGraph::graphs(internal, external, option_values["normal-forms"] == Answer::Yes, option_values["modulo-mirror-images"] == Answer::Yes,
        [internal, &counter, &option_values](KontsevichGraph& g)
        {
            cout << g.encoding();
            if (option_values["with-coefficients"] == Answer::Yes)
            {
                cout << "    ";
                if (g.is_zero() || !g.positive_differential_order())
                    cout << "0";
                else
                    cout << "w_" << internal << "_" << (++counter);
            }
            cout << "\n";
            cout.flush();
        },
        [&option_values](KontsevichGraph& g) -> bool
        {
            bool answer = true;
            answer &= (option_values["prime"] == Answer::Indifferent) ||
                      (option_values["prime"] == Answer::Yes && g.is_prime()) ||
                      (option_values["prime"] == Answer::No && !g.is_prime());
            answer &= (option_values["zero"] == Answer::Indifferent) ||
                      (option_values["zero"] == Answer::Yes && g.is_zero()) ||
                      (option_values["zero"] == Answer::No && !g.is_zero());
            answer &= (option_values["positive-differential-order"] == Answer::Indifferent) ||
                      (option_values["positive-differential-order"] == Answer::Yes && g.positive_differential_order()) ||
                      (option_values["positive-differential-order"] == Answer::No && !g.positive_differential_order());
            return answer;
        }
    );
}
