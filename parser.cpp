#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <stack>

constexpr int INPUT_STRING_INVALID = -1;
using namespace std;


using Terminal = char;
using Nonterminal = char;

using ParseTable = vector<vector<set<size_t>>>;
using TermSet = map< Nonterminal, set<Terminal> >;
using Production = pair<Nonterminal, string>;

constexpr Terminal EPSILON = 'e';
constexpr char END_OF_INPUT = '$';

struct Grammar {
    vector< Production > productions;
    ParseTable parse_table;
    TermSet firsts;
    TermSet follows;

    set<Nonterminal> non_terms;
    set<Terminal> terms;

    explicit Grammar(const vector<Production> &productions) : productions(productions) {
        // Gather all non-terminals
	    non_terms = get_non_terms();
        // Gather all terminals
        terms = get_terms();

        firsts = compute_firsts();
        follows = compute_follows();

        parse_table = build_parse_table();
    }

    [[nodiscard]] Nonterminal starter() const {
        return productions.front().first;
    }

    auto operator[](size_t index) const {
        return productions[index];
    }

    void find_follow(Nonterminal non_term) {

        // cout<<"Finding follow of "<<non_term<<"\n";

        for(const auto& [lhs, rhs] : productions) {
            // finished is true when finding follow from this production is complete
            bool finished = true;
            auto ch = rhs.begin();

            // Skip variables till read non-terminal
            for(;ch != rhs.end() ; ++ch) {
                if(*ch == non_term) {
                    finished = false;
                    break;
                }
            }
            ++ch;
            for(;ch != rhs.end() && !finished; ++ch) {
                // If non-terminal, just append to follow
                if(!isupper(*ch)) {
                    follows[non_term].insert(*ch);
                    finished = true;
                    break;
                }

                set<char> firsts_copy(firsts[*ch]);
                // If char's firsts doesn't have epsilon follow search is over
                if(firsts_copy.find(EPSILON) == firsts_copy.end()) {
                    follows[non_term].insert(firsts_copy.begin(), firsts_copy.end());
                    finished = true;
                    break;
                }
                // Else next char has to be checked after appending firsts to follow
                firsts_copy.erase(EPSILON);
                follows[non_term].insert(firsts_copy.begin(), firsts_copy.end());

            }

            // If end of production, follow same as follow of variable
            if(ch == rhs.end() && !finished) {
                // Find follow if it doesn't have
                if(follows[lhs].empty()) {
                    find_follow(lhs);
                }
                follows[non_term].insert(follows[lhs].begin(), follows[lhs].end());
            }
        }
    }

    void find_first(Nonterminal non_term) {

        // cout<<"Finding firsts of "<<non_term<<"\n";

        for(const auto& [lhs, rhs] : productions) {
            // Find productions of the non-terminal
            if(lhs != non_term) {
                continue;
            }
            // cout<<"Processing production "<<lhs<<"->"<<rhs<<"\n";
            // Loop till a non-terminal or no epsilon variable found
            for(auto ch = rhs.begin(); ch != rhs.end(); ++ch) {
                // If first char in production a non term, add it to firsts list
                if(!isupper(*ch)) {
                    firsts[non_term].insert(*ch);
                    break;
                }
                else {
                    // If char in prod is non-terminal and whose firsts has no yet been found out
                    // Find first for that non-terminal
                    if(firsts[*ch].empty()) {
                        find_first(*ch);
                    }
                    // If variable doesn't have epsilon, stop loop
                    if(firsts[*ch].find(EPSILON) == firsts[*ch].end()) {
                        firsts[non_term].insert(firsts[*ch].begin(), firsts[*ch].end());
                        break;
                    }

                    set<char> firsts_copy(firsts[*ch].begin(), firsts[*ch].end());

                    // Remove epsilon from firsts if not the last variable
                    if(ch + 1 != rhs.end()) {
                        firsts_copy.erase(EPSILON);
                    }

                    // Append firsts of that variable
                    firsts[non_term].insert(firsts_copy.begin(), firsts_copy.end());
                }
            }

        }
    }

    TermSet compute_firsts() {
        TermSet result;
        for(Nonterminal non_term : non_terms) {
            if(result[non_term].empty()){
                find_first(non_term);
            }
        }
        return result;
    }

    TermSet compute_follows() {
        TermSet result;
        // Find follow of start variable first
        auto start_var = starter();
        result[start_var].insert(END_OF_INPUT);
        find_follow(start_var);
        // Find follows for rest of variables
        for(Nonterminal non_term : non_terms) {
            if(result[non_term].empty()) {
                find_follow(non_term);
            }
        }
        return result;
    }

    [[nodiscard]] set<Nonterminal> get_non_terms() const {
        set<Nonterminal> result;
        for(auto & [nont, _] : productions) {
            result.insert(nont);
        }
        return result;
    }


    [[nodiscard]] set<Terminal> get_terms() const {
        set<Terminal> result;
        for(const auto & [_, rhs] : productions) {
            for(char ch : rhs) {
                if(!isupper(ch)) {
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
        ParseTable result{non_terms.size(), vector<set<size_t>>(terms.size()) };

        size_t prod_num = 0;
        for(const auto& [lhs, rhs] : productions) {

            set<char> next_list;
            bool finished = false;
            for(char rh : rhs) {
                if(!isupper(rh)) {
                    if(rh != EPSILON) {
                        next_list.insert(rh);
                        finished = true;
                        break;
                    }
                    continue;
                }

                set<char> firsts_copy(firsts.at(rh).begin(), firsts.at(rh).end());
                if(firsts_copy.find(EPSILON) == firsts_copy.end()) {
                    next_list.insert(firsts_copy.begin(), firsts_copy.end());
                    finished = true;
                    break;
                }
                firsts_copy.erase(EPSILON);
                next_list.insert(firsts_copy.begin(), firsts_copy.end());
            }
            // If the whole rhs can be skipped through epsilon or reaching the end
            // Add follow to next list
            if(!finished) {
                const set<Terminal> &my_follows = follows.at(lhs);
                next_list.insert(my_follows.begin(), my_follows.end());
            }

            size_t row = distance(non_terms.begin(), non_terms.find(lhs));

            for(char ch : next_list) {
                size_t col = distance(terms.begin(), terms.find(ch));
                if(!result[row][col].empty()) {
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


vector<Production> parse_file(fstream& grammar_file) {
    vector<Production> gram;

	while(!grammar_file.eof()) {
		char buffer[20];
		grammar_file.getline(buffer, 19);

		char lhs = buffer[0];
		string rhs = buffer+3;
		pair <char, string> prod (lhs, rhs);
		gram.push_back(prod);
	}
    return gram;
}

template<typename T>
ostream& operator<<(ostream& os, const set<T>& s) {
    os << "{";
    bool is_first = true;
    for (T e : s) {
        if (is_first) {
            is_first = false;
        } else {
            os << ", ";
        }
        os << e;
    }
    return os << "}";
}


void show(const Grammar& gram, ostream& out) {
    out<<"Grammar parsed: \n";
    int count = 0;
    for (const auto& [lhs, rhs] : gram.productions) {
        out << count++ << ".  " << lhs << " -> " << rhs << "\n";
    }
    out<<"\n"
	   <<"The non-terminals in the grammar are: ";
	for(Nonterminal non_term : gram.non_terms) {
		out<<non_term<<" ";
	}
	out<<"\n"
 	   <<"The terminals in the grammar are: ";
	for(Terminal term : gram.terms) {
		out<<term<<" ";
	}
	out<<"\n\n"
	   <<"Firsts list: \n";
	for(auto & first : gram.firsts) {
		out<<first.first<<" : ";
		for(char firsts_it : first.second) {
			out<<firsts_it<<" ";
		}
		out<<"\n";
	}
	out<<"\n"
	   <<"Follows list: \n";
	for(auto & follow : gram.follows) {
		out<<follow.first<<" : ";
		for(char follows_it : follow.second) {
			out<<follows_it<<" ";
		}
		out<<"\n";
	}
	out<<"\n"
       <<"Parsing Table: \n"
       <<"   ";
	for(char term : gram.terms) {
		out<<term<<" ";
	}
	out<<"\n";
    size_t row_num = 0;
	for(auto n : gram.non_terms) {
        out << n << "  ";
        for (int col = 0; col < gram.terms.size(); ++col) {
            out << gram.parse_table[row_num][col] << " ";
        }
        out << "\n";
        row_num++;
    }
}

struct Checker {

    const Grammar& gram;

    explicit Checker(const Grammar &gram) : gram(gram) {
    }

    bool is_acc2(string input_string, stack<char> st, size_t prod) {
        st.pop();
        auto [_, rhs] = gram[prod];
        if (rhs[0] != EPSILON) {
            for (auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                st.push(*ch);
            }
        }
        return is_acc(std::move(input_string), st);
    }
    bool is_acc(string input_string, stack<char> st) {
        // cout<<"Processing input string\n";
        while(!st.empty() && !input_string.empty()) {
            // If stack top same as input string char remove it

            if(input_string[0] == st.top()) {
                st.pop();
                input_string.erase(0, 1);
                continue;
            }
            if(!isupper(st.top())) {
                //cout<<"Unmatched terminal found\n";
                return false;
            }
            char stack_top = st.top();
            size_t row = distance(gram.non_terms.begin(), gram.non_terms.find(stack_top));
            size_t col = distance(gram.terms.begin(), gram.terms.find(input_string[0]));
            auto prod_num = gram.parse_table[row][col];

            if(prod_num.empty()) {
                //cout<<"No production found in parse table\n";
                return false;
            }
            if (prod_num.size() == 1) {
                st.pop();
                string rhs = gram[*prod_num.begin()].second;
                if(rhs[0] == EPSILON) {
                    continue;
                }
                for(auto ch = rhs.rbegin(); ch != rhs.rend(); ++ch) {
                    st.push(*ch);
                }
            } else {
                bool accepted = false;
                for (size_t num : prod_num) {
                    if (is_acc2(input_string, st, num) ) {
                        accepted = true;
                    }
                }
                return accepted;
            }
        }
        return true;
    }

static int is_accepted(string my_input_string,
                const Grammar& gram) {
    stack<char> my_st;
    my_input_string.push_back(END_OF_INPUT);
    my_st.push(END_OF_INPUT);
    my_st.push(gram.starter());

    // Check if input string is valid
    for(char & ch : my_input_string) {
        if(gram.terms.find(ch) == gram.terms.end()) {
            return INPUT_STRING_INVALID;
        }
    }

    return Checker{gram}.is_acc(my_input_string, my_st);
}
};

int main(int argc, char const *argv[])
{
	if(argc != 2) {
		cout<<"Arguments should be <grammar file>\n";
		return 1;
	}
	
	// Parsing the grammar file
	fstream grammar_file;
	grammar_file.open(argv[1], ios::in);
	if(grammar_file.fail()) {
		cout<<"Error in opening grammar file\n";
		return 2;
	}
	Grammar gram{parse_file(grammar_file)};
    show(gram, cout);
    cout<<"\n";

    ifstream rights{ "right-strings.txt" };

    while (rights) {
        string str;
        // cout << "> ";
        std::getline(rights, str);
        if (str.empty()) { continue; }
        int accepted = Checker::is_accepted(str, gram);
        if (accepted == INPUT_STRING_INVALID) {
            cout<< '[' << str << "] has unknown symbols\n";
        }
        else if(accepted) {
            cout<< '[' << str << "] accepted\n";
        }
        else {
            cout << '[' << str << "] rejected\n";
        }
    }

    while (true) {
        string str;
        cout << "> ";
        std::getline(cin, str);
        int accepted = Checker::is_accepted(str, gram);
        if (accepted == INPUT_STRING_INVALID) {
            cout<< '[' << str << "] has unknown symbols\n";
        }
        else if(accepted) {
            cout<< '[' << str << "] accepted\n";
        }
        else {
            cout << '[' << str << "] rejected\n";
        }
    }

	return 0;
}




