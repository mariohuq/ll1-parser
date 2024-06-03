#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>
#include <regex>
#include <random>

using Terminal = char;

using Id = size_t;
using NonterminalId = Id;
using TerminalId = Id;
Id EPSILON_ID = -1;
Id END_OF_INPUT_ID = EPSILON_ID - 1;
size_t TERMINALS_OFFSET = EPSILON_ID / 2;

using ProductionIndex = size_t;
using ParseTable = std::vector<std::vector<std::set<ProductionIndex>>>;
using TermSet = std::map<NonterminalId, std::set<TerminalId>>;
using Production = std::pair<NonterminalId, std::vector<Id>>;

constexpr Terminal EPSILON = 'e';
constexpr Terminal END_OF_INPUT = '$';

enum CheckerResult {
    INPUT_INVALID = -1,
    ACCEPTED = 1,
    REJECTED = 0
};


struct BiMap {

    size_t offset;

    explicit BiMap(size_t offset = 0) : offset(offset) {}

    std::map<Id, std::string> fromId;
    std::map<std::string, Id> toId;
    [[nodiscard]] size_t size() const {
        return fromId.size();
    }
    Id addOrGet(std::string s) {
        const auto it = toId.find(s);
        if (!(it == toId.end())) {
            return it->second;
        }
        Id nextId = size() + offset;
        toId.emplace(s, nextId);
        fromId.emplace(nextId, s);
        return nextId;
    }
    [[nodiscard]] Id get(const std::string& s) const {
        return toId.at(s);
    }
    [[nodiscard]] std::string get(Id id) const {
        return fromId.at(id);
    }

    friend std::ostream &operator<<(std::ostream &os, const BiMap &map) {
        if (map.fromId.empty()) {
            return os << "∅";
        }
        os << "{";
        bool is_first = true;
        for (const auto& [id, s]: map.fromId) {
            if (is_first) {
                is_first = false;
            } else {
                os << ", ";
            }
            if (id == EPSILON_ID) {
                os << "ε";
            } else {
                os << s;
            }
        }
        return os << "}";
        }
};

struct Grammar {
    // Order of declaration is important, because functions rely on previously set fields!
    std::vector<Production> productions;
    BiMap non_terms;
    BiMap terms;
    TermSet firsts;
    TermSet follows;
    ParseTable parse_table;

    explicit Grammar(const std::vector<Production> &productions, BiMap nonterminals, BiMap terminals)
        : productions(productions)
        , non_terms{std::move(nonterminals)}
        , terms{std::move(terminals)}
        , firsts{compute_firsts()}
        , follows{compute_follows()}
        , parse_table{build_parse_table()}
        {}

    [[nodiscard]] NonterminalId starter() const {
        return productions.front().first;
    }

    auto operator[](size_t index) const {
        return productions[index];
    }

    static bool is_terminal(Id id) {
        return id >= TERMINALS_OFFSET;
    }

private:
    void find_follow(TermSet &the_follows, NonterminalId x) const {
        // cout<<"Finding follow of "<<non_term<<"\n";
        for (const auto &[lhs, rhs]: productions) {
            // Skip variables till read non-terminal
            auto ch = std::find(rhs.begin(), rhs.end(), x);
            if (ch == rhs.end()) { continue; }
            // finished when finding FOLLOW from this production is complete
            bool finished = false;
            for (++ch; ch != rhs.end(); ++ch) {
                // If terminal, just append to FOLLOW
                if (is_terminal(*ch)) {
                    the_follows[x].insert(*ch);
                    finished = true;
                    break;
                }

                const auto& ch_firsts = firsts.at(*ch);
                // If char's FIRSTs don't have ϵ FOLLOW search is over
                if (ch_firsts.find(EPSILON_ID) == ch_firsts.end()) {
                    the_follows[x].insert(ch_firsts.begin(), ch_firsts.end());
                    finished = true;
                    break;
                }
                // Else next char has to be checked after appending FIRSTs to FOLLOW
                auto ch_firsts_copy = ch_firsts;
                ch_firsts_copy.erase(EPSILON_ID);
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

    void find_first(TermSet &the_firsts, NonterminalId x) const {
        // cout<<"Finding firsts of "<<non_term<<"\n";
        for (const auto &[lhs, rhs]: productions) {
            // Find productions of the non-terminal
            if (lhs != x) {
                continue;
            }
            // cout<<"Processing production "<<lhs<<"->"<<rhs<<"\n";
            // Loop till a non-terminal or no ϵ found
            for (auto ch = rhs.begin(); ch != rhs.end(); ++ch) {
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
                // If variable doesn't have ϵ, go to next production
                if (ch_firsts.find(EPSILON_ID) == ch_firsts.end()) {
                    the_firsts[x].insert(ch_firsts.begin(), ch_firsts.end());
                    break;
                }
                auto ch_firsts_copy = ch_firsts;
                // Remove ϵ from FIRST if not the last variable
                if (!is_last(ch, rhs.end())) {
                    ch_firsts_copy.erase(EPSILON_ID);
                }
                // Append firsts of that variable
                the_firsts[x].insert(ch_firsts_copy.begin(), ch_firsts_copy.end());
            }
        }
    }

    template<typename It>
    static bool is_last(It it, It end) {
        return std::next(it) != end;
    }

    [[nodiscard]] TermSet compute_firsts() const {
        TermSet result;
        for (const auto& [non_term_id, _]: non_terms.fromId) {
            if (result[non_term_id].empty()) {
                find_first(result, non_term_id);
            }
        }
        return result;
    }

    [[nodiscard]] TermSet compute_follows() const {
        TermSet result;
        // Find follow of start variable first
        auto start_var = starter();
        result[start_var].insert(END_OF_INPUT_ID);
        find_follow(result, start_var);
        // Find follows for rest of variables
        for (const auto & [non_term, _] : non_terms.fromId) {
            if (result[non_term].empty()) {
                find_follow(result, non_term);
            }
        }
        return result;
    }

    [[nodiscard]] ParseTable build_parse_table() const {
        ParseTable result{non_terms.size(), std::vector<std::set<size_t>>(terms.size())};

        size_t prod_num = 0;
        for (const auto &[lhs, rhs]: productions) {

            std::set<TerminalId> next_list;
            bool finished = false;
            for (Id c: rhs) {
                if (is_terminal(c)) {
                    if (c != EPSILON_ID) {
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
            // If the whole rhs can be skipped through ϵ or reaching the end
            // Add FOLLOW to NEXT list
            if (!finished) {
                const auto &my_follows = follows.at(lhs);
                next_list.insert(my_follows.begin(), my_follows.end());
            }
            size_t row = lhs;
            for (TerminalId c: next_list) {
                size_t col = c - TERMINALS_OFFSET;
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


std::pair<std::vector<Production>, std::pair<BiMap, BiMap>> parse_file(std::istream &grammar_file) {
    std::vector<Production> gram;
    BiMap nonts;
    BiMap terms{TERMINALS_OFFSET};

    std::regex line_regex{R"(<(\w+?)>\s*::=\s*(.*?)\s*)"};
    std::regex terminal_regex{R"(`(.+?)')"};
    std::regex split_regex{R"(\s*\|\s*)"};
    std::regex nont_regex{R"(<(\w+?)>)"};

    while (!grammar_file.eof()) {
        std::string line;
        std::getline(grammar_file, line);
        if (line.empty()) {
            continue;
        }
        std::smatch sm;
        if (!std::regex_match(line, sm, line_regex)) {
            std::cerr << "Could not match line [" << line << "]\n";
            continue;
        }
        Id lhs = nonts.addOrGet(sm[1]);
        line = sm[2];
        std::sregex_token_iterator iter{line.begin(), line.end(), split_regex, -1};
        std::sregex_token_iterator end;
        for ( ; iter != end; ++iter) {
            if (std::string{*iter}.empty()) {
                continue;
            }
            std::vector<Id> rhs{};
            std::smatch match;
            std::stringstream ss { *iter };
            while (ss) {
                std::string tok;
                ss >> tok;
                if (tok.empty()) {
                    continue;
                }
                if (std::regex_match(tok, match, nont_regex)) {
                    rhs.push_back(nonts.addOrGet(match[1]));
                }
                else if (std::regex_match(tok, match, terminal_regex)) {
                    rhs.push_back(terms.addOrGet(match[1]));
                }
                else if (tok == std::string{EPSILON}) {
                    rhs.push_back(EPSILON_ID);
                }
                else {
                    std::cerr << "Could not parse gram at " << line << ". Error in ["<< tok <<"]\n";
                    exit(EXIT_FAILURE);
                }
            }
            gram.emplace_back(lhs, rhs);
        }
    }
    // add end character $
    terms.fromId[END_OF_INPUT_ID] = {END_OF_INPUT};
    terms.toId[{END_OF_INPUT}] = END_OF_INPUT_ID;
    return {gram, {nonts, terms}};
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &s) {
    if (s.empty()) {
        return os << "∅";
    }
    os << "{";
    bool is_first = true;
    for (const T& e: s) {
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
    for (int count = 0; const auto &[lhs, rhs]: gram.productions) {
        out << count++ << ". <" << gram.non_terms.fromId.at(lhs) << "> ::= ";
        if (rhs[0] == EPSILON_ID) {
            out << "ϵ\n";
            continue;
        }
        for (Id termnt: rhs) {
            if (Grammar::is_terminal(termnt)) {
                out << "`" << gram.terms.fromId.at(termnt) << "' ";
            } else {
                out << "<" << gram.non_terms.fromId.at(termnt) <<"> ";
            }
        }
        out << "\n";
    }
    out << "\n"
        << "The non-terminals in the grammar are: " << gram.non_terms << "\n"
        << "The terminals in the grammar are: " << gram.terms << "\n"
        << "\n"
        << "Firsts list: \n";
    for (const auto &[nonterminal, set]: gram.firsts) {
        std::set<std::string> names;
        for (auto s :set) {
            if (s == EPSILON_ID) {
                names.insert("ε");
            } else {
            names.insert(gram.terms.get(s));
            }
        }
        out << "FIRST(" << gram.non_terms.get(nonterminal) << ") = " << names << "\n";
    }
    out << "\n"
        << "Follows list: \n";
    for (const auto &[nonterminal, set]: gram.follows) {
        std::set<std::string> names;
        for (auto s :set) {
            if (s == END_OF_INPUT_ID) {
                names.insert("$");
            } else {
                names.insert(gram.terms.get(s));
            }
        }
        out << "FOLLOW(" << gram.non_terms.get(nonterminal) << ") = " << names << "\n";
    }
    out << "\n"
        << "Parsing Table: \n"
        << "\t";
    for (const auto& [id, term]: gram.terms.fromId) {
        out << term << "\t";
    }
    out << "\n";
    for (size_t row_num = 0; auto& [_, n]: gram.non_terms.fromId) {
        out << n << "\t";
        for (const auto &el: gram.parse_table[row_num++]) {
            out << el << "\t";
        }
        out << "\n";
    }
    return out;
}

struct Checker {
    static const std::regex split_regex;
    const Grammar &gram;

    explicit Checker(const Grammar &gram) : gram(gram) {
    }

    bool helper(std::string input, std::stack<Id> stack) {
        std::sregex_token_iterator iter{input.begin(), input.end(), split_regex, -1};
        std::sregex_token_iterator end;

        // cout<<"Processing input string\n";
        while (!stack.empty() && iter != end) {
            if (std::string{*iter}.empty()) {
                iter++;
                continue;
            }
            // If stack top same as input string char remove it

            std::string token = *iter;
            if (Grammar::is_terminal(stack.top())) {
                if (gram.terms.get(stack.top()) == token) {
                    stack.pop();
                    iter++;
                    continue;
                }
                //cout<<"Unmatched terminal found\n";
                return false;
            }
            NonterminalId stack_top = stack.top();
            stack.pop();
            size_t row = stack_top;
            size_t col = gram.terms.toId.at(token) - TERMINALS_OFFSET;
            auto prod_nums = gram.parse_table[row][col];

            if (prod_nums.empty()) {
                //cout<<"No production found in parse table\n";
                return false;
            }
            if (prod_nums.size() == 1) {
                auto [_, rhs] = gram[*prod_nums.begin()];
                if (rhs[0] == EPSILON_ID) continue;
                for (auto ch1 = rhs.rbegin(); ch1 != rhs.rend(); ++ch1) {
                    stack.push(*ch1);
                }
                continue;
            }

            for (size_t num: prod_nums) {
                auto [_, rhs] = gram[num];
                auto _st = stack;
                if (rhs[0] != EPSILON_ID) {
                    for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                        _st.push(*ch);
                    }
                }

                auto it = iter->second - input.begin();
                auto acc = helper(token + " " + input.substr(it), std::move(_st));
                if (acc) {
                    return true;
                }
            }
            return false;
        }
        return true;
    }

    static CheckerResult is_accepted(std::string input,
                                     const Grammar &gram) {
        if (input.empty()) {
            return REJECTED;
        }
        input[0] = tolower(input[0]);
        input.push_back(' ');
        input.push_back(END_OF_INPUT);
        std::stack<Id> stack;
        stack.push(END_OF_INPUT_ID);
        stack.push(gram.starter());

        std::sregex_token_iterator iter{input.begin(), input.end(), split_regex, -1};
        std::sregex_token_iterator end;
        // Check if input string is valid
        for (; iter != end; iter++) {
            std::string token{*iter};
            if (token.empty()) {
                continue;
            }
            if (gram.terms.toId.find(token) == gram.terms.toId.end()) {
                std::cerr << "Unknown word: [" << token << "]\n";
                return INPUT_INVALID;
            }
        }

        if (Checker{gram}.helper(input, stack)) {
            return ACCEPTED;
        } else {
            return REJECTED;
        }
    }
};

const std::regex Checker::split_regex{R"(\s+|(?=\.|\?))"};

void verdict(const std::string &str, CheckerResult accepted) {
    std::cout << '[' << str << "] ";
    switch (accepted) {
        case INPUT_INVALID:
            std::cout << "has unknown words";
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

std::string generate(const Grammar& gram) {
    using D = std::uniform_int_distribution<size_t>;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator

    std::stringstream ss;

    std::vector<std::vector<TerminalId>> alternatives{gram.parse_table.size()};
    for (size_t i = 0; i < alternatives.size(); ++i) {
        for (size_t j = 0; j < gram.parse_table[i].size(); ++j) {
            const auto& s = gram.parse_table[i][j];
            if (s.empty()) {
                continue;
            }
            alternatives[i].push_back(j);
        }
    }
    std::stack<NonterminalId> stack;
    stack.push(gram.starter());
    while (!stack.empty()) {
        if (Grammar::is_terminal(stack.top())) {
            ss << gram.terms.get(stack.top()) << " ";
            stack.pop();
            continue;
        }
        NonterminalId nid = stack.top();
        stack.pop();
        size_t n = alternatives[nid].size();
        auto termid = alternatives[nid][D{0, n-1}(gen)];
        auto variants = gram.parse_table[nid][termid];

        const auto& [_, rhs] = gram.productions[*std::next(variants.begin(),D{0, variants.size() - 1}(gen))];
        if (rhs[0] == EPSILON_ID) continue;
        for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
            stack.push(*ch);
        }
    }
    auto result = ss.str();
    result[0] = toupper(result[0]);
    result.erase(result.size() - 3, 1);
    result.erase(result.size() - 1, 1);
    return result;
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
    auto [a, b] = parse_file(grammar_file);
    Grammar gram{a,b.first, b.second};
    cout << "Grammar parsed: \n" << gram << "\n";

//    std::ifstream rights{"right-strings.txt"};
//    while (rights) {
//        std::string str;
//        std::getline(rights, str);
//        if (str.empty()) { continue; }
//        verdict(str, Checker::is_accepted(str, gram));
//    }



    cout << "Press Ctrl+D to finish, ! to generate\n";
    while (true) {
        std::string str;
        cout << "> ";
        if (!std::getline(std::cin, str)) {
            break;
        }
        if (!str.empty() && str[0] == '!') {
            cout << "Random sentence: [" << generate(gram) << "]\n";
        } else {
            verdict(str, Checker::is_accepted(str, gram));
        }
    }
    return EXIT_SUCCESS;
}
