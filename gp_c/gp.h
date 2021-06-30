//#define VERBOSE
//#define TRACE_OUTPUT
//#define CHOICE_LOG
#define RNG_LOG             # Random Number Generator logging

#define TRUE 1
#define FALSE 0
#define DIE(x) die(x, __FILE__, __LINE__);


// SELECTION METHOD
// ONLY PICK ONE OF THE FOLLOWING
#define USE_FITNESS_PROPORTIONATE_BINARY
//#define USE_FITNESS_PROPORTIONATE_LINEAR
//#define USE_TOURNAMENT

//
//#define RUN_TESTS
//#define SHOW_ALL_TREE_RESULTS
#define SHOW_BEST_TREE_RESULTS
#define SHOW_CONTROLS
//#define SHOW_ALL_TREES
//#define SHOW_TREE_CREATIONS

typedef struct {
    short M;        // Number of individuals in each generation
    short G;        // Number of generations to run
    short Di;       // Maximum depth of S Expressions for an initial tree
    short Dc;       // Maximum depth of S Expressions for a created tree
    double Pc;      // Probability of cross over
    double Pr;      // Probability of reproduction
    double Pip;     // Probability of cross over internal point
    short noFunctions;
    short noTerminals;
    double GRc;      // Greedy C Value (Function of M)
                    // n   M     C   GRc 
                    // 0  <1000  -    -
                    // 0   1000  1   .32
                    // 1   2000  2   .16
                    // 2   4000  4   .08
                    // 3   8000  8   .04
                    //                   GRc = 32/2^n
    short noFitnessCases;
} Control;

typedef struct {
    char type;          // 't' for a terminal, 'f' for a function
} Node;

typedef struct {
    char type;          // always 't'
    short tid;
    char* name;
    gptype (*code)();
} Terminal;

struct FunctionNode;
typedef struct {
    short fid;
    char* name;
    short arity;
    gptype (*code)(struct FunctionNode* fn); 
} Function;

typedef struct FunctionNode {
    char type;          // always 'f'
    Function* func;
    short arity;        // copy of func->arity to improve execution speed
    gptype (*code)(struct FunctionNode* this); 
                        // copy of func->code to improve execution speed
    Node *branch[];     // mixed array of Terminal(s) and FunctionNode(s)
                        // passed as args to code call
} FunctionNode;
 
#define DL_SHIFT    1000000000LL
typedef struct BaseFitness {
    long long lng_nfr;
    long long lng_n;
    long long lng_a;
    long long lng_raw;         // note that if loc Fitness might uses short type
                        // and then this just contains a casted copy
} BaseFitness;

typedef struct {
    short tfid;     // This is Tree's zero based index within TreeSet->treeArray
                   // after sorting for fitness (least fit are lower valued)
    short tcid;     // The id of the Tree when first created and put into the array
    Node* root;
    BaseFitness* fitness;
    short nFunctions;
    short nTerminals;
    int hits;
} Tree;
#define NFR(tree,i)   tree[i]->fitness->lng_nfr

typedef struct NodeLocation {
    Tree* tree;
    Node* node;
    FunctionNode* parent;
    short pos;     // arg number in parent (1,2,3...,parent->arity)
    short cur;     // cur - used in traversing the tree
} NodeLocation;

typedef struct {
    short n;            // current size (population)
    short size;         // max size (population)
                        // (After population is created, normally n == size)
    char haveWinner;   // A winner is an individual with a best possible
                       // score
    short winningIndex; //
    double avgRawF;
    Tree** treeArray;
    char*  tagArray;    // used for temporary individual state marking
                        // e.g. used during tournament selection to insure
                        // unique individuals compete
} TreeSet;

extern gptype execNode(Node* node);
extern FunctionNode* newFunctionNode(short fid);
extern Tree* newTree(Node* root);
extern void initBaseFitness(BaseFitness* f, short raw);
extern void printTree(Tree* t);
extern void tp(short i);
extern int rnd(short max);
extern void setArg(FunctionNode *fn, short pos, Node* arg);
extern double intToFloat(long long llval);
extern long long float_to_int(double dblval);
extern int choice_log_count;
