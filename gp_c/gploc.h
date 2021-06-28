// Initialization & Execution
extern void initLoc();
extern void execTrees(TreeSet* trees, int run_number, int gen);
extern void execSingleTree(Tree* tree);
extern void clockTic();
extern void terminate();

// Reporting
extern void treeResultHeader(short gen);
extern void reportTreeResult(Tree* t, short i, short gen, double avgRawF);

// Data 
extern Control control;
extern Terminal terminal[];
extern Function function[];

// Sample trees
extern Tree* createWinnerSampleTree();

#define TERM_MOVE       0
#define TERM_LEFT       1
#define TERM_RIGHT      2
#define NO_TERMINALS    3

#define FUNC_IF_FOOD_AHEAD 0
#define FUNC_PROG_N2       1
#define FUNC_PROG_N3       2
#define NO_FUNCTIONS       3

