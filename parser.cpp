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

const std::set<char> BRACES = {
'(',')',
'[',']',
'{','}',
'<','>' };

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

    static bool is_semantic(char ch) {
        return BRACES.find(ch) != BRACES.end();
    }

private:
    void find_follow(TermSet &the_follows, Nonterminal x) const {
        // cout<<"Finding follow of "<<non_term<<"\n";
        for (const auto &[lhs, rhs]: productions) {
            // Skip variables till read non-terminal
            auto ch = std::find(rhs.begin(), rhs.end(), x);
            if (ch == rhs.end()) { continue; }
            // finished when finding FOLLOW from this production is complete
            bool finished = false;
            for (++ch; ch != rhs.end(); ++ch) {
                if (is_semantic(*ch)) {
                    continue;
                }
                // If terminal, just append to FOLLOW
                if (is_terminal(*ch)) {
                    the_follows[x].insert(*ch);
                    finished = true;
                    break;
                }

                const auto& ch_firsts = firsts.at(*ch);
                // If char's FIRSTs don't have ε FOLLOW search is over
                if (ch_firsts.find(EPSILON) == ch_firsts.end()) {
                    the_follows[x].insert(ch_firsts.begin(), ch_firsts.end());
                    finished = true;
                    break;
                }
                // Else next char has to be checked after appending FIRSTs to FOLLOW
                auto ch_firsts_copy = ch_firsts;
                ch_firsts_copy.erase(EPSILON);
                the_follows[x].insert(ch_firsts_copy.begin(), ch_firsts_copy.end());
            }
            if (finished) { continue; }
            // If end of production reached, FOLLOW ⊃ FOLLOW of lhs
            if (ch == rhs.end()) {
                // Find FOLLOW if it doesn't have
                if (the_follows[lhs].empty()) {
                    find_follow(the_follows, lhs);
                }
                the_follows[x].insert(the_follows[lhs].begin(), the_follows[lhs].end());
            }
        }
    }

    void find_first(TermSet &the_firsts, Nonterminal x) const {
        // cout<<"Finding firsts of "<<non_term<<"\n";
        for (const auto &[lhs, rhs]: productions) {
            // Find productions of the non-terminal
            if (lhs != x) {
                continue;
            }
            // cout<<"Processing production "<<lhs<<"->"<<rhs<<"\n";
            // Loop till a non-terminal or no ε found
            for (auto ch = rhs.begin(); ch != rhs.end(); ++ch) {
                if (is_semantic(*ch)) {
                    continue;
                }
                // If first char in production a terminal, add it to firsts list
                if (is_terminal(*ch)) {
                    the_firsts[x].insert(*ch);
                    break;
                }
                // If char in rhs is non-terminal and whose FIRST has not yet been found out
                // Find FIRST for that non-terminal
                const auto& ch_firsts = the_firsts[*ch];
                if (ch_firsts.empty()) {
                        find_first(the_firsts, *ch);
                }
                // If variable doesn't have ε, go to next production
                if (ch_firsts.find(EPSILON) == ch_firsts.end()) {
                    the_firsts[x].insert(ch_firsts.begin(), ch_firsts.end());
                    break;
                }
                auto ch_firsts_copy = ch_firsts;
                // Remove ε from FIRST if not the last variable
                if (!is_last(ch, rhs.end())) {
                    ch_firsts_copy.erase(EPSILON);
                }
                // Append firsts of that variable
                the_firsts[x].insert(ch_firsts_copy.begin(), ch_firsts_copy.end());
            }
        }
    }

    template<typename It>
    static bool is_last(It it, It end) {
        return std::find_if_not(it, end, is_semantic) != end;
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
                if (is_terminal(ch) && !is_semantic(ch)) {
                    result.insert(ch);
                }
            }
        }
        // Remove ε and add end character $
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
            for (char c: rhs) {
                if (is_semantic(c)) { continue; }
                if (is_terminal(c)) {
                    if (c != EPSILON) {
                        next_list.insert(c);
                        finished = true;
                        break;
                    }
                    continue;
                }

                auto firsts_copy = firsts.at(c);
                if (firsts_copy.find(EPSILON) == firsts_copy.end()) {
                    next_list.insert(firsts_copy.begin(), firsts_copy.end());
                    finished = true;
                    break;
                }
                firsts_copy.erase(EPSILON);
                next_list.insert(firsts_copy.begin(), firsts_copy.end());
            }
            // If the whole rhs can be skipped through ε or reaching the end
            // Add FOLLOW to NEXT list
            if (!finished) {
                const auto &my_follows = follows.at(lhs);
                next_list.insert(my_follows.begin(), my_follows.end());
            }
            size_t row = distance(non_terms.begin(), non_terms.find(lhs));
            for (char c: next_list) {
                size_t col = distance(terms.begin(), terms.find(c));
                //if (!result[row][col].empty()) {
                    // cout<<"Collision at ["<<lhs<<"]["<<ch<<"] for production "<<prod_num<<"\n";
                    // continue;
                //}
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
    if (s.empty()) {
        return os << "∅";
    }
    os << "{";
    bool is_first = true;
    for (T e: s) {
        if (is_first) {
            is_first = false;
        } else {
            os << ", ";
        }
        if (e == EPSILON) {
            os << "ε";
        } else {
            os << e;
        }
    }
    return os << "}";
}


std::ostream &operator<<(std::ostream &out, const Grammar &gram) {
    for (int count = 0; const auto &[lhs, rhs]: gram.productions) {
        out << count++ << ". " << lhs << " → " << (rhs[1] == EPSILON ? (std::string{ rhs[0] } + "ε" + rhs[2]) : rhs) << "\n";
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

    std::pair<bool, std::string> helper(std::string input, std::stack<char> stack) {
        std::string output;
        // cout<<"Processing input string\n";
        while (!stack.empty() && !input.empty()) {
            if (Grammar::is_semantic(stack.top())) {
                output.push_back(stack.top());
                stack.pop();
                continue;
            }

            // If stack top same as input string char remove it
            if (input[0] == stack.top()) {
                if (input[0] != END_OF_INPUT) {
                    output.push_back(input[0]);
                }
                stack.pop();
                input.erase(0, 1);
                continue;
            }
            if (Grammar::is_terminal(stack.top())) {
                //cout<<"Unmatched terminal found\n";
                return {false, {}};
            }
            Nonterminal stack_top = stack.top();
            stack.pop();
            output.push_back(stack_top);
            size_t row = distance(gram.non_terms.begin(), gram.non_terms.find(stack_top));
            size_t col = distance(gram.terms.begin(), gram.terms.find(input[0]));
            auto prod_nums = gram.parse_table[row][col];

            if (prod_nums.empty()) {
                //cout<<"No production found in parse table\n";
                return {false, {}};
            }
            if (prod_nums.size() == 1) {
                auto [_, rhs] = gram[*prod_nums.begin()];
                if (rhs[1] != EPSILON) {
                    for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                        stack.push(*ch);
                    }
                } else {
                    stack.push(rhs[2]);
                    stack.push(rhs[0]);
                }
                continue;
            }

            for (size_t num: prod_nums) {
                auto [_, rhs] = gram[num];
                auto _st = stack;
                if (rhs[1] != EPSILON) {
                    for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                        _st.push(*ch);
                    }
                } else {
                    _st.push(rhs[2]);
                    _st.push(rhs[0]);
                }
                auto [acc, out] = helper(input, std::move(_st));
                if (acc) {
                    return {true, output + out};
                }
            }
            return {false, {}};
        }
        return {true, output};
    }

    static std::pair<CheckerResult, std::string> is_accepted(std::string input,
                                     const Grammar &gram) {
        input.push_back(END_OF_INPUT);
        std::stack<char> stack;
        stack.push(END_OF_INPUT);
        stack.push(gram.starter());

        // Check if input string is valid
        for (char c: input) {
            if (gram.terms.find(c) == gram.terms.end()) {
                return {INPUT_INVALID, {}};
            }
        }

        auto [acc, out] = Checker{gram}.helper(input, stack);
        if (acc) {
            return {ACCEPTED, out};
        } else {
            return {REJECTED, {}};
        }
    }
};

void verdict(const std::string &str, const std::pair<CheckerResult, std::string>& accepted) {
    std::cout << '[' << str << "] ";
    switch (accepted.first) {
        case INPUT_INVALID:
            std::cout << "has unknown symbols";
            break;
        case ACCEPTED:
            std::cout << "accepted. Output: " << accepted.second;
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
        verdict(str, Checker::is_accepted(str, gram));
    }
    cout << "Press Ctrl+D to finish\n";
    while (true) {
        std::string str;
        cout << "> ";
        if (!std::getline(std::cin, str)) {
            break;
        }
        verdict(str, Checker::is_accepted(str, gram));
    }
    return EXIT_SUCCESS;
}
