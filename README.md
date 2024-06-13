# LL(1) Parser

An _nondeterministic_ LL(1) parser made in C++ with STL libraries.

Like the original program, it prints out the first, follow, parse table and checks if the given input string is accepted by the given grammar.

Also, it generates random string that accepted by the grammar (enter `!`).

More over, it accepts _semantic_ symbols in grammar (four kinds of braces `{}[]()<>`) which are used to produce brace expression representing parsing tree of input string.

## Usage

Compile the code with g++ and run with arguments `parser <grammar-file>`

Grammar file contains all the productions for the grammar. It follows these rules:

* Each line contains a production for the grammar in the format `LHS -> RHS`
* All upper case alphabets are non-terminals or variables.
* `e` is considered as empty string.
* `$` is used as END character.
* Any other character is taken as a terminal for the grammar.

Sample Usage:

```bash
g++ parser.cpp -o parser
./parser grammar0.txt
```

Output:

```text
Grammar parsed: 
0. S → (aB)
1. S → (bA)
2. A → [aT]
3. A → [baA]
4. B → {bT}
5. B → {aBB}
6. T → <ϵ>
7. T → <S>

The non-terminals in the grammar are: {A, B, S, T}
The terminals in the grammar are: {$, a, b}

Firsts list: 
FIRST(A) = {a, b}
FIRST(B) = {a, b}
FIRST(S) = {a, b}
FIRST(T) = {a, b, ε}

Follows list: 
FOLLOW(A) = {$}
FOLLOW(B) = {$, a, b}
FOLLOW(S) = {$, a, b}
FOLLOW(T) = {$, a, b}

Parsing Table: 
	$	a	b	
A	∅	{2}	{3}	
B	∅	{5}	{4}	
S	∅	{0}	{1}	
T	{6}	{6, 7}	{6, 7}	

> aabaabbabb
[aabaabbabb] accepted. Output: S(aB{aB{bT<>}B{aB{aB{bT<>}B{bT<>}}B{aB{bT<>}B{bT<>}}}})
```

# `parser2`

`parser2` is the same as `parser` except these differences:

1. It accepts multicharacter terminals and nonterminals.
2. Grammar file format is different:
    ```
   <StartingNont> ::= `terminal1' | e | <A>
   <A> ::= `terminal2'
   ```
   
   - uses `::=` for splitting left- and right-hand-side of rules
   - terminals are quoted using `` ` `` and `'` characters
   - nonterminals are quoted using `<` and `>`
   - `e` again means empty string
   - `$` again means "end of input" terminal that is appended to input string before processing 
3. Input string is preprocessed before input and postprocessed to match English sentence rules for capitalizing first letter and use no space before the dot.
4. no semantic symbols in this program

To run use `parser2 grammarFutureCont.txt`.