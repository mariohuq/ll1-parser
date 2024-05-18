#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>

using Terminal = char;
using Nonterminal = char;

using ParseTable = std::vector<std::vector<std::set<size_t>>>;
using TermSet = std::map<Nonterminal, std::set<Terminal>>;
using Production = std::pair<Nonterminal, std::string>;

constexpr Terminal EPSILON = 'e';
constexpr Terminal END_OF_INPUT = '$';

enum CheckerResult {
    INPUT_INVALID = -1,
    ACCEPTED = 1,
    REJECTED = 0
};

struct Grammar {
    // Order of declaration is important, because functions rely on previously set fields!
    std::vector<Production> productions;
    std::set<Nonterminal> non_terms;
    std::set<Terminal> terms;
    TermSet firsts;
    TermSet follows;
    ParseTable parse_table;

    explicit Grammar(const std::vector<Production> &productions)
        : productions(productions)
        , non_terms{get_non_terms()}
        , terms{get_terms()}
        , firsts{compute_firsts()}
        , follows{compute_follows()}
        , parse_table{build_parse_table()}
        {}

    [[nodiscard]] Nonterminal starter() const {
        return productions.front().first;
    }

    auto operator[](size_t index) const {
        return productions[index];
    }

    static bool is_terminal(char ch) {
        return !isupper(ch);
    }

private:
    void find_follow(TermSet &the_follows, Nonterminal non_term) const {

        // cout<<"Finding follow of "<<non_term<<"\n";

        for (const auto &[lhs, rhs]: productions) {
            // finished when finding FOLLOW from this production is complete
            bool finished = true;
            auto ch = rhs.begin();

            // Skip variables till read non-terminal
            //ch = std::find(ch, rhs.end(), non_term);
            for (; ch != rhs.end(); ++ch) {
                if (*ch == non_term) {
                    finished = false;
                    break;
                }
            }
            ++ch;
            for (; ch != rhs.end() && !finished; ++ch) {
                // If terminal, just append to FOLLOW
                if (is_terminal(*ch)) {
                    the_follows[non_term].insert(*ch);
                    finished = true;
                    break;
                }

                auto firsts_copy = firsts.at(*ch);
                // If char's FIRSTs don't have ϵ FOLLOW search is over
                if (firsts_copy.find(EPSILON) == firsts_copy.end()) {
                    the_follows[non_term].insert(firsts_copy.begin(), firsts_copy.end());
                    finished = true;
                    break;
                }
                // Else next char has to be checked after appending FIRSTs to FOLLOW
                firsts_copy.erase(EPSILON);
                the_follows[non_term].insert(firsts_copy.begin(), firsts_copy.end());
            }

            // If end of production reached, FOLLOW is same as FOLLOW of variable
            if (ch == rhs.end() && !finished) {
                // Find FOLLOW if it doesn't have
                if (the_follows[lhs].empty()) {
                    find_follow(the_follows, lhs);
                }
                the_follows[non_term].insert(the_follows[lhs].begin(), the_follows[lhs].end());
            }
        }
    }

    void find_first(TermSet &the_firsts, Nonterminal non_term) const {

        // cout<<"Finding firsts of "<<non_term<<"\n";

        for (const auto &[lhs, rhs]: productions) {
            // Find productions of the non-terminal
            if (lhs != non_term) {
                continue;
            }
            // cout<<"Processing production "<<lhs<<"->"<<rhs<<"\n";
            // Loop till a non-terminal or no epsilon variable found
            for (auto ch = rhs.begin(); ch != rhs.end(); ++ch) {
                // If first char in production a non term, add it to firsts list
                if (is_terminal(*ch)) {
                    the_firsts[non_term].insert(*ch);
                    break;
                } else {
                    // If char in prod is non-terminal and whose firsts has no yet been found out
                    // Find first for that non-terminal
                    if (the_firsts[*ch].empty()) {
                        find_first(the_firsts, *ch);
                    }
                    // If variable doesn't have epsilon, stop loop
                    if (the_firsts[*ch].find(EPSILON) == the_firsts[*ch].end()) {
                        the_firsts[non_term].insert(the_firsts[*ch].begin(), the_firsts[*ch].end());
                        break;
                    }

                    auto firsts_copy = the_firsts[*ch];

                    // Remove epsilon from firsts if not the last variable
                    if (ch + 1 != rhs.end()) {
                        firsts_copy.erase(EPSILON);
                    }

                    // Append firsts of that variable
                    the_firsts[non_term].insert(firsts_copy.begin(), firsts_copy.end());
                }
            }

        }
    }

    [[nodiscard]] TermSet compute_firsts() const {
        TermSet result;
        for (Nonterminal non_term: non_terms) {
            if (result[non_term].empty()) {
                find_first(result, non_term);
            }
        }
        return result;
    }

    [[nodiscard]] TermSet compute_follows() const {
        TermSet result;
        // Find follow of start variable first
        auto start_var = starter();
        result[start_var].insert(END_OF_INPUT);
        find_follow(result, start_var);
        // Find follows for rest of variables
        for (Nonterminal non_term: non_terms) {
            if (result[non_term].empty()) {
                find_follow(result, non_term);
            }
        }
        return result;
    }

    [[nodiscard]] std::set<Nonterminal> get_non_terms() const {
        std::set<Nonterminal> result;
        for (auto &[lhs, _]: productions) {
            result.insert(lhs);
        }
        return result;
    }

    [[nodiscard]] std::set<Terminal> get_terms() const {
        std::set<Terminal> result;
        for (const auto &[_, rhs]: productions) {
            for (char ch: rhs) {
                if (is_terminal(ch)) {
                    result.insert(ch);
                }
            }
        }
        // Remove epsilon and add end character $
        result.erase(EPSILON);
        result.insert(END_OF_INPUT);
        return result;
    }

    [[nodiscard]] ParseTable build_parse_table() const {
        ParseTable result{non_terms.size(), std::vector<std::set<size_t>>(terms.size())};

        size_t prod_num = 0;
        for (const auto &[lhs, rhs]: productions) {

            std::set<char> next_list;
            bool finished = false;
            for (char rh: rhs) {
                if (is_terminal(rh)) {
                    if (rh != EPSILON) {
                        next_list.insert(rh);
                        finished = true;
                        break;
                    }
                    continue;
                }

                auto firsts_copy = firsts.at(rh);
                if (firsts_copy.find(EPSILON) == firsts_copy.end()) {
                    next_list.insert(firsts_copy.begin(), firsts_copy.end());
                    finished = true;
                    break;
                }
                firsts_copy.erase(EPSILON);
                next_list.insert(firsts_copy.begin(), firsts_copy.end());
            }
            // If the whole rhs can be skipped through epsilon or reaching the end
            // Add follow to next list
            if (!finished) {
                const auto &my_follows = follows.at(lhs);
                next_list.insert(my_follows.begin(), my_follows.end());
            }

            size_t row = distance(non_terms.begin(), non_terms.find(lhs));

            for (char ch: next_list) {
                size_t col = distance(terms.begin(), terms.find(ch));
                if (!result[row][col].empty()) {
                    // cout<<"Collision at ["<<lhs<<"]["<<ch<<"] for production "<<prod_num<<"\n";
                    // continue;
                }
                result[row][col].insert(prod_num);
            }
            prod_num++;
        }
        return result;
    }
};


std::vector<Production> parse_file(std::istream &grammar_file) {
    std::vector<Production> gram;

    while (!grammar_file.eof()) {
        char buffer[20];
        grammar_file.getline(buffer, sizeof(buffer));
        gram.emplace_back(buffer[0], buffer + 3);
    }
    return gram;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &s) {
    os << "{";
    bool is_first = true;
    for (T e: s) {
        if (is_first) {
            is_first = false;
        } else {
            os << ", ";
        }
        os << e;
    }
    return os << "}";
}


std::ostream &operator<<(std::ostream &out, const Grammar &gram) {
    for (int count = 1; const auto &[lhs, rhs]: gram.productions) {
        out << count++ << ". " << lhs << " → " << (rhs[0] == EPSILON ? "ϵ" : rhs) << "\n";
    }
    out << "\n"
        << "The non-terminals in the grammar are: " << gram.non_terms << "\n"
        << "The terminals in the grammar are: " << gram.terms << "\n"
        << "\n"
        << "Firsts list: \n";
    for (const auto &[nonterminal, set]: gram.firsts) {
        out << "FIRST(" << nonterminal << ") = " << set << "\n";
    }
    out << "\n"
        << "Follows list: \n";
    for (const auto &[nonterminal, set]: gram.follows) {
        out << "FOLLOW(" << nonterminal << ") = " << set << "\n";
    }
    out << "\n"
        << "Parsing Table: \n"
        << "\t";
    for (char term: gram.terms) {
        out << term << "\t";
    }
    out << "\n";
    for (size_t row_num = 0; Nonterminal n: gram.non_terms) {
        out << n << "\t";
        for (const auto &el: gram.parse_table[row_num++]) {
            out << el << "\t";
        }
        out << "\n";
    }
    return out;
}

struct Checker {

    const Grammar &gram;

    explicit Checker(const Grammar &gram) : gram(gram) {
    }

    bool is_acc2(std::string input_string, std::stack<char> st, size_t prod) {
        st.pop();
        auto [_, rhs] = gram[prod];
        if (rhs[0] != EPSILON) {
            for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                st.push(*ch);
            }
        }
        return is_acc(std::move(input_string), st);
    }

    bool is_acc(std::string input_string, std::stack<char> st) {
        // cout<<"Processing input string\n";
        while (!st.empty() && !input_string.empty()) {
            // If stack top same as input string char remove it

            if (input_string[0] == st.top()) {
                st.pop();
                input_string.erase(0, 1);
                continue;
            }
            if (Grammar::is_terminal(st.top())) {
                //cout<<"Unmatched terminal found\n";
                return false;
            }
            char stack_top = st.top();
            size_t row = distance(gram.non_terms.begin(), gram.non_terms.find(stack_top));
            size_t col = distance(gram.terms.begin(), gram.terms.find(input_string[0]));
            auto prod_num = gram.parse_table[row][col];

            if (prod_num.empty()) {
                //cout<<"No production found in parse table\n";
                return false;
            }
            if (prod_num.size() == 1) {
                st.pop();
                std::string rhs = gram[*prod_num.begin()].second;
                if (rhs[0] == EPSILON) {
                    continue;
                }
                for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                    st.push(*ch);
                }
            } else {
                bool accepted = false;
                for (size_t num: prod_num) {
                    if (is_acc2(input_string, st, num)) {
                        accepted = true;
                    }
                }
                return accepted;
            }
        }
        return true;
    }

    static CheckerResult is_accepted(std::string my_input_string,
                                     const Grammar &gram) {
        my_input_string.push_back(END_OF_INPUT);
        std::stack<char> my_st;
        my_st.push(END_OF_INPUT);
        my_st.push(gram.starter());

        // Check if input string is valid
        for (char &ch: my_input_string) {
            if (gram.terms.find(ch) == gram.terms.end()) {
                return INPUT_INVALID;
            }
        }

        return Checker{gram}.is_acc(my_input_string, my_st) ? ACCEPTED : REJECTED;
    }
};

void verdict(const std::string &str, CheckerResult accepted) {
    std::cout << '[' << str << "] ";
    switch (accepted) {
        case INPUT_INVALID:
            std::cout << "has unknown symbols";
            break;
        case ACCEPTED:
            std::cout << "accepted";
            break;
        case REJECTED:
            std::cout << "rejected";
            break;
    }
    std::cout << "\n";
}

int main(int argc, char const *argv[]) {
    using std::cout;
    if (argc != 2) {
        cout << "Usage:\n"
             << argv[0] << " <path to grammar file>\n";
        return EXIT_FAILURE;
    }

    // Parsing the grammar file
    std::ifstream grammar_file{argv[1], std::ios::in};
    if (grammar_file.fail()) {
        cout << "Error in opening grammar file\n";
        return EXIT_FAILURE;
    }
    Grammar gram{parse_file(grammar_file)};
    cout << "Grammar parsed: \n" << gram << "\n";

    std::ifstream rights{"right-strings.txt"};
    while (rights) {
        std::string str;
        std::getline(rights, str);
        if (str.empty()) { continue; }
        auto accepted = Checker::is_accepted(str, gram);
        verdict(str, accepted);
    }
    cout << "Press Ctrl+D to finish\n";
    while (true) {
        std::string str;
        cout << "> ";
        if (!std::getline(std::cin, str)) {
            break;
        }
        auto accepted = Checker::is_accepted(str, gram);
        verdict(str, accepted);
    }
    return EXIT_SUCCESS;
}
