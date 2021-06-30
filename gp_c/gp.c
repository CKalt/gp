#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "gptype.h"
#include "gp.h"
#include "gploc.h"
#include "time.h"
#include "string.h"

double intToFloat(long long llval) {
    return ((double) llval) / ((double) DL_SHIFT);
}

long long float_to_int(double dblval) {
    double lldbl = dblval * ((double) DL_SHIFT);
    double r_lldbl = lldbl > 0.0 ?
        ( lldbl + (double)0.5 ) : ( lldbl - (double)0.5 );

    long long llResult = (long long) r_lldbl;
    return llResult;
}

void die(char* msg, const char* file, int line) {
    printf("die %s(%d): %s\n", file, line, msg);
    exit(0);
}

#ifdef RNG_LOG
    int TRACE_COUNT = 0;
    FILE *rngFp = NULL;
    int log_rand() {
        assert(rngFp != NULL);
        int r = rand();
        fprintf(rngFp, "%d,%d\n", ++TRACE_COUNT, r);
        return r;
    }
    #define GP_RAND(x) log_rand(x)
#else
    #define GP_RAND(x) rand(x)
#endif
#ifdef CHOICE_LOG
FILE *choiceFp = NULL;
int choice_log_count = 0;

char choiceLogBuf[1024];
short get_choice_from_log(short expectedTp) {
    int tp, result;
    if (fgets(choiceLogBuf, sizeof(choiceLogBuf), choiceFp) == NULL) {
        printf("Ran out of choice.log lines.");
#ifdef RNG_LOG
        if (rngFp != NULL) {
            fclose(rngFp);
        }
#endif
        exit(1);
    }

    sscanf(choiceLogBuf, "%d, %d,%d", &choice_log_count, &tp, &result);

    if ((short) tp != expectedTp) {
        printf("tp at log=%d was %d, but expected %d.", choice_log_count, tp,
                expectedTp);
        exit(1);
    }

    return (short) result;
}
#endif

void tp(short i) {
    printf("tp(%d)\n", i);
}

FunctionNode* newFunctionNode(short fid) {
    Function* fnc = &function[fid];
    FunctionNode *fn = (FunctionNode *)malloc(
                        sizeof(FunctionNode) + (fnc->arity)*sizeof(void*)
                        );
    fn->type = 'f';
    fn->func = fnc;
    fn->code = fnc->code;
    fn->arity = fnc->arity;

    return fn;
}

//
// rnd - return random short value between 0 and max (inclusive)
//       i.e. 0 >= return value =< max
int rnd(short max) {
    int r = GP_RAND();
    double d = r/(double)RAND_MAX;
    double f = d*max + 0.5;
    double f2 = (d * ((double)max)) + 0.5;
    int result = floorl(f);
    return result;

}

//
// rnd - return random double value between 0.0 and 1.0 (inclusive)
//   
double rnd_dbl() {
    double result = (double)GP_RAND()/(double)RAND_MAX;
    return result;
}

long long rnd_greedy_dbl_as_ll() {
    double r = rnd_dbl();
    double result;

//#ifdef TRACE_OUTPUT
//    printf("TP005.2:(rnd_greedy_dbl) rlog=%d,r=%lf\n", TRACE_COUNT, r);
//#endif

    if (control.GRc < .0001) {
        result = r;
    }
    else if (rnd(9) > 1) {
        // select from top set
        result = ((double)1.0) - control.GRc*r;
    }
    else {
        // select from lower set
        result = ((double)1.0 - control.GRc)*r;
    }
    long long llResult = float_to_int(result);

    return llResult;
}

Node* newRndFTNode() {
#ifdef CHOICE_LOG
    short r = get_choice_from_log(2);
    if (r < control.noTerminals) {
        return (Node*)&terminal[r];
    }
    else {
        return (Node*)newFunctionNode(r - control.noTerminals);
    }
#else 
    short r = rnd(control.noTerminals + control.noFunctions - 1);
    if (r < control.noTerminals) {
        return (Node*)&terminal[r];
    }
    else {
        return (Node*)newFunctionNode(r - control.noTerminals);
    }
#endif
}

FunctionNode* newRndFNode() {
#ifdef CHOICE_LOG
    short r = get_choice_from_log(4);
#else
    short r = rnd(control.noFunctions - 1);
#endif
    return newFunctionNode(r);
}

Terminal* getRndTNode() {
#ifdef CHOICE_LOG
    short r = get_choice_from_log(3);
#else
    short r = rnd(control.noTerminals - 1);
#endif
    return &terminal[r];
}

void setArg(FunctionNode *fn, short pos, Node* arg) {
    // pos: 1,2,...fn->arity
    assert(fn->type == 'f');
    assert(pos >=1 && pos <= fn->arity);
    fn->branch[pos-1] = arg;
}

Node* getArg(FunctionNode *fn, short pos) {
    // pos: 1,2,...fn->arity
    assert(fn->type == 'f');
    assert(pos >=1 && pos <= fn->arity);
    return fn->branch[pos-1];
}

void genTreeFullMethodR(FunctionNode* fn, short level, short depth) {
    int i;
    if (level < depth) {
        short cDepth = level+1;
        for (i = 1; i <= fn->arity; i++) {
            FunctionNode* cfn = newRndFNode();
            setArg(fn, i, (Node*)cfn);
            genTreeFullMethodR(cfn, cDepth, depth);
        }
    }
    else {
        for (i = 1; i <= fn->arity; i++) {
            Terminal* ct = getRndTNode();
            setArg(fn, i, (Node*)ct);
        }
    }
}

Tree* newTree(Node* root) {
    Tree* tree = malloc(sizeof(Tree));
    tree->root = root;
    tree->fitness = NULL;
    tree->hits = 0;
    tree->nFunctions = -1;
    tree->nTerminals = -1;
    tree->tfid = tree->tcid = -1;
    return tree;
}

void initBaseFitness(BaseFitness* f, short raw) {
    f->lng_n = f->lng_a = f->lng_nfr = -1LL;
    f->lng_raw = ((long long) raw) * DL_SHIFT;
}

Tree* genTreeFullMethod(short depth) {
    FunctionNode* root = newRndFNode();
    genTreeFullMethodR(root, 2, depth);
    Tree *tree = newTree((Node*)root);
    return tree;
}

void genTreeGrowMethodR(FunctionNode* fn, short level, short depth) {
    int i;
    if (level < depth) {
        short cDepth = level+1;
        for (i = 1; i <= fn->arity; i++) {
            Node* cn = (Node*)newRndFTNode();
            setArg(fn, i, cn);
            if (cn->type == 'f') {
                genTreeGrowMethodR((FunctionNode*) cn, cDepth, depth);
            }
        }
    }
    else {
        for (i = 1; i <= fn->arity; i++) {
            Terminal* ct = getRndTNode();
            setArg(fn, i, (Node*)ct);
        }
    }
}

Tree* genTreeGrowMethod(short depth) {
    FunctionNode* root = newRndFNode();
    genTreeGrowMethodR(root, 2, depth);
    Tree *tree = newTree((Node*)root);
    return tree;
}

Node* cloneNode(Node *node) {
    if (node->type == 't') {
        short tid = ((Terminal*) node)->tid;
        return (Node*)&terminal[tid];
    }
    else {
        assert(node->type == 'f');
        FunctionNode* fn = (FunctionNode*)node;
        FunctionNode* newFn = newFunctionNode(fn->func->fid);
        int i;
        for (i = 1; i <= fn->arity; i++) {
            Node* cn = getArg(fn, i);
            Node* newCn = cloneNode(cn);
            setArg(newFn, i, newCn);
        }
        assert(newFn->type == 'f');
        return (Node*)newFn;
    }
}

Tree* cloneTree(Tree* t) {
    Node* newRoot = cloneNode(t->root);
    Tree *tree = newTree(newRoot);
    tree->nFunctions = t->nFunctions;
    tree->nTerminals = t->nTerminals;
    return tree;
}

void tab(short count) {
    while (count--) {
        printf(" ");
    }
}

void printTreeR(Node* node, short depth) {
    if (node->type == 't') {
        Terminal* t = (Terminal*)node;
        printf(" %s", t->name);
    }
    else {
        int i;
        assert(node->type == 'f');
        FunctionNode* fn = (FunctionNode*)node;
        printf("\n");
        tab(depth*2);
        printf("(%s", fn->func->name);
        for (i = 0; i < fn->arity; i++) {
            Node* sub = fn->branch[i];
            printTreeR(sub, depth+1);
        }
        if (depth == 0) {
            for (i = 0; i < fn->arity; i++) {
                Node* sub = fn->branch[i];
                if (sub->type == 'f') {
                    printf("\n");
                    break;      //  if there is one or more functions
                                //  we skip a line for readability
                }
            }
        }

        printf(")");
    }
}

void printNode(Node* node) {
    assert(node->type == 't' || node-> type == 'f');
    printTreeR(node, 0);
}

void printNorNone(short num) {
    if (num == -1) {
        printf("None");
    } else {
        printf("%d", num);
    }
}

void printOptional(short num) {
    if (num == -1) {
        printf("None");
    } else {
        printf("Some(%d)", num);
    }
}

void printTree(Tree* t) {
    char tcidStr[20];

#ifdef TRACE_OUTPUT
    #ifdef CHOICE_LOG
        printf("-log=%-10d------tree (", choice_log_count);
    #else
        printf("-log=%-10d------tree (", 0);
    #endif
    printNorNone(t->tfid);
    printf("/");
    printNorNone(t->tcid);
    printf(")----------------\n");
#else
    printf("----------------tree (");
    printNorNone(t->tfid);
    printf("/");
    printNorNone(t->tcid);
    printf(")----------------\n");
#endif
    printf("nFunctions = ");
    printOptional(t->nFunctions);
    printf("\nnTerminals= ");
    printOptional(t->nTerminals);
    printf("\n");

    printNode(t->root);
    printf("\n");
}

Tree* printTreeFromID(TreeSet* trees, short i) {
    Tree** treeArray = trees->treeArray;
    Tree* t = treeArray[i];
    printTree(t);
    return t;
}

void printTrees(TreeSet* trees) {
    short i;
    short size = trees->n;
    for (i = 0; i < size; i++) {
        printTreeFromID(trees, i);
    }
}

gptype execNode(Node* node) {
    if (node->type == 't') {
        Terminal* t = (Terminal*)node;
        return t->code();
    }
    else {
        assert(node->type == 'f');
        FunctionNode* fn = (FunctionNode*)node;
        return fn->code(fn);
    }
}

void freeTreeNode(Node *node) {
    if (node->type == 'f') {
        int i;
        FunctionNode* fn = (FunctionNode*)node;
        for (i = 0; i < fn->arity; i++) {
            Node* sub = fn->branch[i];
            freeTreeNode(sub);
        }
        free(node);
    }
    else {
        assert(node->type == 't');
    }
}

void freeTree(Tree* t) {
    freeTreeNode(t->root);
    if (t->fitness)
        free(t->fitness);
    free(t);
}

char treeNodeMatch(Node* n1, Node* n2) {
    if (n1->type != n2->type) {
        return FALSE;
    }
    else if (n1->type == 't') {
        return ((Terminal*) n1)->tid == ((Terminal*) n2)->tid;
    }
    else {
        int i;
        assert(n1->type == 'f');
        FunctionNode* fn1 = (FunctionNode*)n1;
        FunctionNode* fn2 = (FunctionNode*)n2;

        Function *f1 = fn1->func;
        Function *f2 = fn2->func;
        if (f1->fid != f2->fid) {
            return FALSE;
        }

        for (i = 0; i < f1->arity; i++) {
            Node* sub1 = fn1->branch[i];
            Node* sub2 = fn2->branch[i];
            if (!treeNodeMatch(sub1, sub2)) {
                return FALSE;
            }
        }
    }
    return TRUE;
}

char treeMatch(Tree* t1, Tree* t2) {
    return treeNodeMatch(t1->root, t2->root);
}

char treeMatchExists(TreeSet* trees, Tree* t) {
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    int i = 0;
    for (i = 0; i < n; i++) {
        if (treeMatch(treeArray[i], t))
            return TRUE;
    }
    return FALSE;
}

Tree* createUniqueTree(TreeSet* trees, short d) {
    int i = 1;
    char fullMethod = (trees->n % 2) == 0;
    Tree* t = fullMethod ? genTreeFullMethod(d) : genTreeGrowMethod(d);
    while (treeMatchExists(trees, t)) {
        freeTree(t);
        if (i++ > 10) {
            d++;
            i = 1;
        }
        t = fullMethod ? genTreeFullMethod(d) : genTreeGrowMethod(d);
    }

    return t;
}

TreeSet* newTreeSet() {
    short size = control.M;
    TreeSet* trees = malloc(sizeof(TreeSet));
    trees->treeArray = calloc(sizeof(Tree*), size);
    trees->tagArray = calloc(sizeof(char), size);
    trees->size = size;
    trees->haveWinner = FALSE;
    trees->n = 0;
    trees->avgRawF = 0.0;
    return trees;
}

void clearTagArray(TreeSet* trees) {
    int i;
    char* tagArray = trees->tagArray;
    for (i = 0; i < trees->n; i++) {
        tagArray[i] = 0;
    }
}

void pushTree(TreeSet* trees, Tree* t) {
    short n = trees->n;
    assert(n < trees->size);
    if (n < trees->size) {
        Tree** treeArray = trees->treeArray;
        treeArray[n] = t;
        t->tfid = -1;
        t->tcid = n;
        trees->n = n+1;
    }
    else DIE("pushTree: array exeeded");
}

void countTreeNodesR(Tree* t, Node *node) {
    if (node->type == 't') {
        t->nTerminals++;
    }
    else {
        assert(node->type == 'f');
        t->nFunctions++;
        FunctionNode* fn = (FunctionNode*)node;
        int i;
        for (i = 1; i <= fn->arity; i++) {
            Node* cn = getArg(fn, i);
            countTreeNodesR(t, cn);
        }
    }
}

Tree* countTreeNodes(Tree* t) {
    Node* root = t->root;
    t->nFunctions = t->nTerminals = 0;
    countTreeNodesR(t, root);
    return t;
}

TreeSet* createInitialPopulation() {
    // Following Koza's recipe, he calls "ramped half-and-half",
    // we will evenly produce population segments starting with 2
    // up to the maxium depth (control.Di) and alternate
    // between Full Method and Grow Method for S Expressions.
    //

#ifdef TRACE_OUTPUT
    Tree *debugTree;
#endif
    
    TreeSet* trees = newTreeSet();
    short size = trees->size;
    short i = trees->n;   // s/b zero of course 'cause we just created array
    assert(i == 0);
    int d;
    double seg = (double) size / ((double)control.Di - 1.0);
    double bdr = 0;
#ifdef TRACE_OUTPUT
    printf("TP001:create_init_pop start\n");
#endif
    for (d = 2; d <= control.Di; d++) {
        bdr += seg;
        while (trees->n < bdr && trees->n < size) {

#ifdef TRACE_OUTPUT
            pushTree(trees, debugTree=createUniqueTree(trees, d));
            countTreeNodes(debugTree);
            printTree(debugTree);
#else
            pushTree(trees, createUniqueTree(trees, d));
#endif

        }
    }

    // fill out to end in case there are "left-overs" due to rounding
    while(trees->n < size) {
#ifdef TRACE_OUTPUT
        pushTree(trees, debugTree=createUniqueTree(trees, control.Di));
        countTreeNodes(debugTree);
        printTree(debugTree);
#else
        pushTree(trees, createUniqueTree(trees, control.Di));
#endif
    }
#ifdef TRACE_OUTPUT
    printf("TP002:create_init_pop done\n");
#endif
    assert(trees->n == trees->size);

    return trees;
}

void computeNormalizedFitness(TreeSet* trees) {
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    short i;
    long long sum_a = 0LL;
    long long sumraw = 0LL;
    for (i = 0; i < n; i++) {
        Tree* t = treeArray[i];
        sum_a += t->fitness->lng_a;
        sumraw += t->fitness->lng_raw;
    }
    long long avgRawF = (long long) (((double)sumraw) / ((double)trees->n));
    trees->avgRawF = intToFloat(avgRawF);
    double flt_sum_a = intToFloat(sum_a);

    for (i = 0; i < n; i++) {
        Tree* t = treeArray[i];
        double dbl_n = intToFloat(t->fitness->lng_a) / flt_sum_a;
        t->fitness->lng_n = float_to_int(dbl_n);

#ifdef TRACE_OUTPUT
        printf("TP1:i=%d, tcid=%d, hits=%d, t.fitness.n=%10.9lf a=%10.9lf sum_a=%10.9lf\n",
                      i, t->tcid, t->hits,
                      intToFloat(t->fitness->lng_n),
                      intToFloat(t->fitness->lng_a),
                      flt_sum_a);

//        if (t->tcid == 2) {
//            printTree(t);
//            execSingleTree(t);
//        }

#endif
    }
}

static int cmpTreesNFitness(const void* p1, const void* p2) {
    const Tree* t1 = (* (Tree * const *) p1);
    const Tree* t2 = (* (Tree * const *) p2);
    BaseFitness* f1 = t1->fitness;
    BaseFitness* f2 = t2->fitness;
    long long n1 = f1->lng_n;
    long long n2 = f2->lng_n;

    if (n1 > n2) {
        return 1;
    } else if (n1 < n2) {
        return -1;
    } else {
        return 0;
    }
}

void sortByNormalizedFitness(TreeSet *trees) {
    short i;
    Tree** treeArray = trees->treeArray;
    qsort(&treeArray[0], trees->n, sizeof(Tree*), cmpTreesNFitness);
    for (i = 0; i< trees->n; i++)
        treeArray[i]->tfid = i;
}

void assignNFRankings(TreeSet *trees) {
    int i;
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    long long nfr = 0LL;
    for (i = 0; i < n; i++) {
        Tree* t = treeArray[i];
        nfr += t->fitness->lng_n;
        t->fitness->lng_nfr = nfr;
    }
}

void countNodes(TreeSet* trees) {
    short i;
    short size = trees->n;
    Tree** treeArray = trees->treeArray;
    for (i = 0; i < size; i++) {
        Tree* t = treeArray[i];
        countTreeNodes(t);
    }
}

short getRndTerminalIndex(Tree* t) {
    return rnd(t->nTerminals - 1);
}

short getRndFunctionIndex(Tree* t) {
    return rnd(t->nFunctions - 1);
}

short getRndTreeIndex(TreeSet* trees) {
    short n = trees->n;
    return rnd(n-1);
}

void freeTreeSet(TreeSet* trees) {
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    short i;
    for (i = 0; i < n; i++) {
        Tree* t = treeArray[i];
        freeTree(t);
    }
    free(treeArray);
    free(trees->tagArray);
    free(trees);
}

char rndReproduction() {
    return rnd_dbl() < control.Pr;
}

#if 1

char useReproduction(short i) {
#ifdef CHOICE_LOG
    return (char) get_choice_from_log(9);
#else
    return ((double)i)/((double)control.M) < control.Pr;
#endif
}
#else
// this is the GpInt way
char useReproduction(short i) {
    long long int_index_ratio = float_to_int((double)i / (double) control.M);
    long long int_control_ratio = float_to_int((double)control.Pr);

    return int_index_ratio < int_control_ratio;
}
#endif

char rndInternalPoint() {
#ifdef CHOICE_LOG
    return (char) get_choice_from_log(1);
#else 
    return rnd_dbl() < control.Pip;
#endif
}

#ifdef USE_FITNESS_PROPORTIONATE_BINARY
short selectIndBinR(TreeSet* trees, long long level, short lo, short hi) {
    short gap = ((hi - lo) / 2);
    short guess = lo + gap + 1;

    if (trees->treeArray[guess-1]->fitness->lng_nfr < level &&
        trees->treeArray[guess]->fitness->lng_nfr >= level) {
        return guess;
    }
    else if (trees->treeArray[guess]->fitness->lng_nfr < level) {
        // new_lo = guess
        return selectIndBinR(trees, level, guess, hi);
    }
    else {
        // new_hi = guess
        return selectIndBinR(trees, level, lo, guess-1);
    }
}

//
// selectIndBin - uses binary search to locate fitness proportional
// individual.  (This requires that trees->treeArray is sorted by
// normalized fitness measure.)
//
short selectIndBin(TreeSet* trees, long long level) {
    assert(trees-> n > 0);
    if (trees->treeArray[0]->fitness->lng_nfr >= level ||
        trees->n == 1) {
        return 0;
    } else if (level > trees->treeArray[trees->n - 2]->fitness->lng_nfr) {
        return trees->n - 1;
    } else {
        return selectIndBinR(trees, level, 0, trees->n - 1);
    }
}

Tree* selectTree(TreeSet* trees) {
#ifdef CHOICE_LOG
    short i = get_choice_from_log(6);
#else
    short i = selectIndBin(trees, rnd_greedy_dbl_as_ll());
#endif
    return trees->treeArray[i];
}
#endif
#ifdef USE_FITNESS_PROPORTIONATE_LINEAR
//
// selectIndLinear - slower but does not require
// sorted array, also does not work with greedy
// selection
short selectIndLinear(TreeSet* trees, double r) {
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    int i;
    double sum = 0.0;
    for (i=0; i<n; i++) {
        sum += treeArray[i]->fitness->n;
        if (sum > r)
            return i;
    }
    assert(0);
    return 0;
}

Tree* selectTree(TreeSet* trees) {
    short i = selectIndLinear(trees, rnd_greedy_dbl_as_ll());
    return trees->treeArray[i];
}
#endif
#ifdef USE_TOURNAMENT
Tree* selectTree(TreeSet* trees) {
    // Return best tree among seven distinct and randomly selected trees
    // tagArray is used to all seven are distinct.
    assert(trees-> n > 50);         // anything lower try other selection

    Tree** treeArray = trees->treeArray;
    short lastTi = trees->n - 1;
    short r;

    Tree* bestTree = treeArray[r = rnd(lastTi)];
    BaseFitness* bestFitness = bestTree->fitness;
    double bestN = bestFitness->n;
    clearTagArray(trees);
    char* tagArray = trees->tagArray;
    tagArray[r] = 1;

    short i;
    for (i = 2; i <= 7; i++ ) {
        short sanity = 100;
        do {
            r = rnd(lastTi);
            if (--sanity < 1) {
                DIE("selectTree: sanity limit reached.");
            }
        } while (tagArray[r]);
        tagArray[r] = 1;

        Tree* trialTree = treeArray[r];

        BaseFitness* trialFitness = trialTree->fitness;
        double trialN = trialFitness->n;
        int diff = (int) (100000.0*bestN - 100000.0*trialN);

        if (diff < 0) {
            bestTree = trialTree;
            bestN = trialN;
        }
    }

    return bestTree;
}
#endif

char findTreeFunctionNodeR(
    NodeLocation *nl, FunctionNode* parent, short pos, FunctionNode* fn, 
    short fi
) {
    assert(fn->type == 'f');
    int i;
    if (fi == nl->cur) {
        nl->node = (Node*)fn;
        nl->parent = parent;
        nl->pos = pos;
        assert(nl->node->type == 'f');
        return TRUE;
    }
    assert(nl->cur < fi);
    for (i = 1; i <= fn->arity; i++) {
        Node* cn = getArg(fn, i);
        if (cn->type == 'f') {
            nl->cur++;
            if (findTreeFunctionNodeR(
                    nl, fn, i, (FunctionNode*)cn, fi)
            ) {
                return TRUE;
            }
        }
    }
    return FALSE;
}

void findTreeFunctionNode(NodeLocation* nl, Tree* t, short fi) {
    char result;
    Node* root = t->root;
    if (root->type == 't') {
        // we never can find this function and we should DIE right here
        printTree(t);
        assert(0);
        DIE("Invalid attempt to find a Function node with a terminal tree.");
    }

    nl->tree = t;
    assert(fi < t->nFunctions);
    nl->cur = 0;  // when fi == nl->cur we've found it
    result = findTreeFunctionNodeR(nl, NULL, 0, (FunctionNode*)root, fi);
    assert(result);
}

void findRndTreeFunctionNode(NodeLocation* nl, Tree* t) {
    findTreeFunctionNode(nl, t, getRndFunctionIndex(t));
}

char findTreeTerminalR(
    NodeLocation *nl, short pos, FunctionNode* fn, 
    short ti
) {
    assert(fn->type == 'f');
    int i;
    for (i = 1; i <= fn->arity; i++) {
        Node* cn = getArg(fn, i);
        if (cn->type == 't') {
            if (ti == nl->cur) {
                nl->node = cn;
                nl->parent = fn;
                nl->pos = i;
                return TRUE;
            }
            assert(nl->cur < ti);
            nl->cur++;
        }
        else {
            if (findTreeTerminalR(nl, i, (FunctionNode*)cn, ti)) {
                return TRUE;
            }
        }
    }
    return FALSE;
}

void findTreeTerminal(NodeLocation* nl, Tree* t, short ti) {
    char result;
    Node* root = t->root;
    nl->tree = t;
    if (root->type == 't') {
        // Rare - but crossover can result in this.  If so
        // we only have a single terminal so this is the one
        nl->node = root;
        nl->parent = NULL;
        nl->pos = 0;
        assert(ti == 0);
        return;
    }
    assert(ti < t->nTerminals);
    nl->cur = 0; // when ti == nl->cur we've found it
    result = findTreeTerminalR(nl, 0, (FunctionNode*)root, ti);
    assert(result);
}

void findRndTreeTerminal(NodeLocation* nl, Tree* t) {
    findTreeTerminal(nl, t, getRndTerminalIndex(t));
}

void performCrossover(Tree* t1, Tree* t2) {
    assert(t1->nTerminals > 0 && t2->nTerminals > 0);
    NodeLocation nl1, nl2;

    if (t1->nFunctions > 0 && rndInternalPoint()) {
#ifdef CHOICE_LOG
        short fi = get_choice_from_log(7);
#else
        short fi = rnd(t1->nFunctions - 1);
#endif
        findTreeFunctionNode(&nl1, t1, fi);
        assert(fi > 0 || nl1.parent == NULL);
        assert(nl1.tree == t1);
    }
    else {
#ifdef CHOICE_LOG
        short ti = get_choice_from_log(8);
#else
        short ti = rnd(t1->nTerminals - 1);
#endif
        findTreeTerminal(&nl1, t1, ti);
        assert(nl1.tree == t1);
    }

    if (t2->nFunctions > 0 && rndInternalPoint()) {
#ifdef CHOICE_LOG
        short fi = get_choice_from_log(7);
#else
        short fi = rnd(t2->nFunctions - 1);
#endif
        findTreeFunctionNode(&nl2, t2, fi);
        assert(fi > 0 || nl2.parent == NULL);
        assert(nl2.tree == t2);
    }
    else {
#ifdef CHOICE_LOG
        short ti = get_choice_from_log(8);
#else
        short ti = rnd(t2->nTerminals - 1);
#endif
        findTreeTerminal(&nl2, t2, ti);
        assert(nl2.tree == t2);
    }

    if (nl1.parent == NULL) {
        nl1.tree->root = nl2.node;
    }
    else {
        setArg(nl1.parent, nl1.pos, nl2.node);
    }

    if (nl2.parent == NULL) {
        nl2.tree->root = nl1.node;
    }
    else {
        setArg(nl2.parent, nl2.pos, nl1.node);
    }
}

#ifdef RUN_TESTS
void runTests(TreeSet* trees) {
    int i;
    int internalPoint = 0;
    countNodes(trees);
    printf("Testing getRndFunctionIndex and getRndTerminalIndex\n");
    for (i = 1; i <= 10; i++) {
        printf("=============\n| Trial %3d |\n=============\n", i);
        NodeLocation nl;
        short ti, fi;
#ifdef CHOICE_LOG
        short id = get_choice_from_log(5);
#else
        short id = getRndTreeIndex(trees);
#endif
        Tree* t = printTreeFromID(trees, id);

#ifdef CHOICE_LOG
        fi = get_choice_from_log(7);
#else
        fi = getRndFunctionIndex(t);
#endif
        printf("\n--------\nRnd fi=%d Subtree -->\n", fi);
        findTreeFunctionNode(&nl, t, fi);
        assert(nl.node->type == 'f');
        printNode(nl.node);

#ifdef CHOICE_LOG
        ti = get_choice_from_log(8);
#else
        ti = getRndTerminalIndex(t);
#endif
        printf("\n--------\nRnd ti=%d Subtree -->\n", ti);
        findTreeTerminal(&nl, t, ti);
        printNode(nl.node);
        printf("\n--------\n");
    }
}
#endif

void reportResults(int gen, TreeSet* trees, char* headerNeed) {
//    int i;
#ifdef SHOW_ALL_TREES
    printf("Generation %d\n", gen);
    printTrees(trees);
#endif

#ifdef SHOW_ALL_TREE_RESULTS
    Tree** treeArray = trees->treeArray;
    short n = trees->n;
    printf("gen %d\n", gen);
    treeResultHeader(-1);
    for (int i = 0; i < n; i++) {
        Tree* t = treeArray[i];
        reportTreeResult(t, i, -1, -1.0);
    }
#endif
#ifdef SHOW_BEST_TREE_RESULTS
    {
        Tree** treeArray = trees->treeArray;
        short n = trees->n;
        short i = n-1;
        if (*headerNeed) {
            treeResultHeader(gen);
            *headerNeed = (char) FALSE;
        }
        Tree* t = treeArray[i];
        reportTreeResult(t, i, gen, trees->avgRawF);
    }
#endif

#ifdef RUN_TESTS
    runTests(trees);
#endif
    fflush(stdout);
}

char nodeDepthGT(Node* node, short d, short sofar) {
    if (sofar > d)
        return TRUE;
    if (node->type == 't')
        return FALSE;

    assert(node->type == 'f');
    FunctionNode* fn = (FunctionNode*)node;
    int i;
    for (i = 1; i <= fn->arity; i++) {
        Node* cn = getArg(fn, i);
        assert(cn->type == 'f' || cn->type == 't');
        if (nodeDepthGT(cn, d, sofar+1))
            return TRUE;
    }
    return FALSE;
}

char treeDepthGT(Tree* t, short d) {
    return nodeDepthGT(t->root, d, 1);
}

char newTreeQualifies(Tree* t) {
    return !treeDepthGT(t, control.Dc);
}

Tree* run(int run_number) {
    int i;
    int gen = 0;
    char haveWinner;
    Tree* winner = NULL;
#ifdef SHOW_CONTROLS
    printf("M = %d, G = %d, D = %d\n", control.M, control.G, control.Di);
#endif

    TreeSet* trees = createInitialPopulation();

    countNodes(trees);

    char headerNeed = TRUE;
    while (gen <= control.G && !trees->haveWinner) { 
        execTrees(trees, run_number, gen);
        if (trees->haveWinner)
            break;

#ifdef TRACE_OUTPUT
    printf("TP0:trace log for gen=%d\n", gen);
#endif
        computeNormalizedFitness(trees);
#ifdef TRACE_OUTPUT
//        if (gen == 9) {
//            DIE("pause");
//        }
//        else {
//            printf("at pause point gen=%d skipping until gen=1 here.\n", gen);
//        }
#endif
        sortByNormalizedFitness(trees);
        assignNFRankings(trees);

        reportResults(gen, trees, &headerNeed);

//if (gen == 2) { 
//    DIE("pause");
//}

        if (gen >= control.G) {
            // We're done
            break;
        }

        TreeSet* trees2 = newTreeSet();
#ifdef TRACE_OUTPUT
        Tree *debugTree;
        printf("TP003:breed start\n");
#endif
        for (i = 0; i < control.M; i++) {
            if (useReproduction(i)) {
#ifdef TRACE_OUTPUT
                printf("TP003.1: breeding reproduction branch\n");
#endif
                // do reproduction
                Tree* t = selectTree(trees);
#ifdef TRACE_OUTPUT
                pushTree(trees2, debugTree=countTreeNodes(cloneTree(t)));
                countTreeNodes(debugTree);
                printTree(debugTree);
#else
                pushTree(trees2, countTreeNodes(cloneTree(t)));
#endif
            }
            else {
                // do crossover
#ifdef TRACE_OUTPUT
                printf("TP003.2: breeding crossover branch\n");
#endif
                Tree *nt1, *nt2;
                for(;;) {
                    Tree* t1 = selectTree(trees);
                    Tree* t2 = selectTree(trees);

                    nt1 = cloneTree(t1);
                    nt2 = cloneTree(t2);
                    performCrossover(nt1, nt2);

                    if (newTreeQualifies(nt1) && newTreeQualifies(nt2))
                        break;

                    freeTree(nt1);
                    freeTree(nt2); 
                }
#ifdef TRACE_OUTPUT
                pushTree(trees2, debugTree=countTreeNodes(nt1));
                countTreeNodes(debugTree);
                printTree(debugTree);
#else
                pushTree(trees2, countTreeNodes(nt1));
#endif

                if (++i < control.M) {
#ifdef TRACE_OUTPUT
                    pushTree(trees2, debugTree=countTreeNodes(nt2));
                    countTreeNodes(debugTree);
                    printTree(debugTree);
#else
                    pushTree(trees2, countTreeNodes(nt2));
#endif
                }
            }
        }

        assert(trees2->n == control.M);
#ifdef TRACE_OUTPUT
        printf("TP004:breed done\n");
#endif
        freeTreeSet(trees);
        trees = trees2;
        gen++;
    }
    if (haveWinner = trees->haveWinner) {
        winner = cloneTree(trees->treeArray[trees->winningIndex]);
    }
    
    freeTreeSet(trees);
    return winner;
}

#if 0
void runSample() {
    Tree* sample = createWinnerSampleTree();
    execSingleTree(sample);
    treeResultHeader(0);
    reportTreeResult(sample, 0, 0);
    freeTree(sample);
}
#endif

int main(int argc, char* argv[]) {
#ifdef CHOICE_LOG
    choiceFp = fopen("choice.log", "r");
    if (choiceFp == NULL) {
        DIE("could not open choice.log. Does it exist?");
    }
#endif
#ifdef RNG_LOG
    rngFp = fopen("rng.log", "w");
    if (rngFp == NULL) {
        DIE("could not open rng.log for writing.");
    }
#endif
    int seed = time(NULL) % 932982389;
//    srand(seed);
//    srand(28);  // 11
    srand(9);   // winner during 8th run
    initLoc();
    int totalRuns = 0;
    Tree* winner = NULL;
    for (;;) {
        printf("Run #%d\n", ++totalRuns);
        if (winner = run(totalRuns)) 
            break;
    }
    printf("total_runs=%d\n", totalRuns);
    printTree(winner);
    execSingleTree(winner);
    freeTree(winner);
#ifdef CHOICE_LOG
    if (choiceFp != NULL) {
        fclose(choiceFp);
    }
#endif
#ifdef RNG_LOG
    if (rngFp != NULL) {
        fclose(rngFp);
    }
#endif
}
