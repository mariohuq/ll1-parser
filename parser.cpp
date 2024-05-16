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
using Grammar = vector< pair<Nonterminal, string> >;
using ParseTable = vector<vector<set<size_t>>>;
using TermSet = map< Nonterminal, set<Terminal> >;

constexpr Nonterminal EPSILON = 'e';
constexpr char END_OF_INPUT = '$';
constexpr size_t NO_LINK = -1;

void find_first(const Grammar& gram,
	TermSet &firsts,
	Nonterminal non_term);

void find_follow(const Grammar& gram,
	TermSet &follows,
	TermSet firsts,
	Nonterminal non_term);

TermSet compute_firsts(const Grammar& gram, const set<Nonterminal>& non_terms) {
    TermSet result;
    for(Nonterminal non_term : non_terms) {
		if(result[non_term].empty()){
			find_first(gram, result, non_term);
		}
	}
    return result;
}

TermSet compute_follows(const Grammar& gram, const set<Nonterminal>& non_terms, const TermSet& firsts) {
	TermSet result;
	// Find follow of start variable first
	char start_var = gram.begin()->first;
	result[start_var].insert(END_OF_INPUT);
	find_follow(gram, result, firsts, start_var);
	// Find follows for rest of variables
	for(char non_term : non_terms) {
		if(result[non_term].empty()) {
			find_follow(gram, result, firsts, non_term);
		}
	}
    return result;
}

Grammar parse_file(fstream& grammar_file) {
    Grammar gram;

	while(!grammar_file.eof()) {
		char buffer[20];
		grammar_file.getline(buffer, 19);

		char lhs = buffer[0];
		string rhs = buffer+3;
		pair <char, string> prod (lhs, rhs);
		gram.push_back(prod);
	}
	cout<<"\n";
    return gram;
}

void show(const Grammar& gram) {
    cout<<"Grammar parsed: \n";
    int count = 0;
    for (const auto& prod : gram) {
        cout<<count++<<".  "<<prod.first<<" -> "<<prod.second<<"\n";
    }
}

set<Nonterminal> get_non_terms(const Grammar& gram) {
    set<Nonterminal> result;
	for(auto & i : gram) {
		result.insert(i.first);
	}
    return result;
}


set<Terminal> get_terms(const Grammar& gram) {
    set<Terminal> result;
	for(auto & i : gram) {
		for(char ch : i.second) {
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

ParseTable build_parse_table(const Grammar& gram,
                const set<Nonterminal>& non_terms,
                const set<Terminal>& terms, const TermSet& firsts, const TermSet& follows) {
    ParseTable parse_table{ non_terms.size(), vector<set<size_t>>(terms.size()) };

	for(auto prod = gram.begin(); prod != gram.end(); ++prod) {
		string rhs = prod->second;

		set<char> next_list;
		bool finished = false;
		for(char & rh : rhs) {
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
            const set<Terminal> &my_follows = follows.at(prod->first);
            next_list.insert(my_follows.begin(), my_follows.end());
		}

		for(char ch : next_list) {
			size_t row = distance(non_terms.begin(), non_terms.find(prod->first));
			size_t col = distance(terms.begin(), terms.find(ch));
			size_t prod_num = distance(gram.begin(), prod);
			if(!parse_table[row][col].empty()) {
				//cout<<"Collision at ["<<prod->first<<"]["<<ch<<"] for production "<<prod_num<<"\n";
				//continue;
			}
			parse_table[row][col].insert(prod_num);
		}

	}
    return parse_table;
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

struct Checker {

    const Grammar& gram;
    const set<Nonterminal>& non_terms;
    const set<Terminal>& terms;
    const ParseTable& parse_table;

    bool is_acc2(string input_string, stack<char> st, size_t prod) {
        st.pop();
        string rhs = gram[prod].second;
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
            size_t row = distance(non_terms.begin(), non_terms.find(stack_top));
            size_t col = distance(terms.begin(), terms.find(input_string[0]));
            auto prod_num = parse_table[row][col];

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
                const Grammar& gram,
                const set<Nonterminal>& non_terms,
                const set<Terminal>& terms,
                const ParseTable& parse_table) {
    stack<char> my_st;
    my_input_string.push_back(END_OF_INPUT);
    my_st.push(END_OF_INPUT);
    my_st.push(gram.front().first);

    // Check if input string is valid
    for(char & ch : my_input_string) {
        if(terms.find(ch) == terms.end()) {
            return INPUT_STRING_INVALID;
        }
    }

    return Checker{gram, non_terms, terms, parse_table}.is_acc(my_input_string, my_st);
}
};

string generate(const Grammar& gram, const ParseTable& parse_table, const set<Nonterminal>& non_terms, const set<Terminal>& terms) {
    string result = {gram[0].first};
    char random_ch = rand() % 2 ? 'a' : 'b';
    size_t row = distance(non_terms.begin(), non_terms.find(result.back()));
    size_t col = distance(terms.begin(), terms.find(random_ch));
    auto nexts = parse_table[row][col];
}


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
	Grammar gram = parse_file(grammar_file);
    show(gram);

	// Gather all non-terminals
	set<Nonterminal> non_terms = get_non_terms(gram);
	cout<<"The non-terminals in the grammar are: ";
	for(Nonterminal non_term : non_terms) {
		cout<<non_term<<" ";
	}
	cout<<"\n";
	// Gather all terminals
    set<Terminal> terms = get_terms(gram);
 	cout<<"The terminals in the grammar are: ";
	for(Terminal term : terms) {
		cout<<term<<" ";
	}
	cout<<"\n\n";

	TermSet firsts = compute_firsts(gram, non_terms);

	cout<<"Firsts list: \n";
	for(auto & first : firsts) {
		cout<<first.first<<" : ";
		for(char firsts_it : first.second) {
			cout<<firsts_it<<" ";
		}
		cout<<"\n";
	}
	cout<<"\n";


    auto follows = compute_follows(gram, non_terms, firsts);

	cout<<"Follows list: \n";
	for(auto & follow : follows) {
		cout<<follow.first<<" : ";
		for(char follows_it : follow.second) {
			cout<<follows_it<<" ";
		}
		cout<<"\n";
	}
	cout<<"\n";

    auto parse_table = build_parse_table(gram, non_terms, terms, firsts, follows);

	// Print parse table
	cout<<"Parsing Table: \n";
	cout<<"   ";
	for(char term : terms) {
		cout<<term<<" ";
	}
	cout<<"\n";
	for(auto row = non_terms.begin(); row != non_terms.end(); ++row) {
		cout<<*row<<"  ";
		for(int col = 0; col < terms.size(); ++col) {
			size_t row_num = distance(non_terms.begin(), row);
			if(parse_table[row_num][col].empty()) {
				cout<<"- ";
				continue;
			}
			cout<<parse_table[row_num][col]<<" ";
		}
		cout<<"\n";
	}
	cout<<"\n";

    ifstream rights{ "right-strings.txt" };

    while (rights) {
        string str;
//        cout << "> ";
        std::getline(rights, str);
        if (str.empty()) { continue; }
        int accepted = Checker::is_accepted(str, gram, non_terms, terms, (ParseTable)parse_table);
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
        int accepted = Checker::is_accepted(str, gram, non_terms, terms, (ParseTable)parse_table);
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

void find_first(const Grammar& gram,
	map< Nonterminal , set<Terminal> > &firsts,
	Nonterminal non_term) {

	// cout<<"Finding firsts of "<<non_term<<"\n";

	for(auto it = gram.begin(); it != gram.end(); ++it) {
		// Find productions of the non-terminal
		if(it->first != non_term) {
			continue;
		}

		// cout<<"Processing production "<<it->first<<"->"<<it->second<<"\n";

		string rhs = it->second;
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
					find_first(gram, firsts, *ch);
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

void find_follow(const Grammar& gram,
	TermSet &follows,
	TermSet firsts,
	Nonterminal non_term) {

	// cout<<"Finding follow of "<<non_term<<"\n";

	for(auto it = gram.begin(); it != gram.end(); ++it) {

		// finished is true when finding follow from this production is complete
		bool finished = true;
		auto ch = it->second.begin();

		// Skip variables till read non-terminal
		for(;ch != it->second.end() ; ++ch) {
			if(*ch == non_term) {
				finished = false;
				break;
			}
		}
		++ch;

		for(;ch != it->second.end() && !finished; ++ch) {
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
		if(ch == it->second.end() && !finished) {
			// Find follow if it doesn't have
			if(follows[it->first].empty()) {
				find_follow(gram, follows, firsts, it->first);
			}
			follows[non_term].insert(follows[it->first].begin(), follows[it->first].end());
		}
	}

}
