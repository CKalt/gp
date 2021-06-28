#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <assert.h>
#include "gptype.h"
#include "gp.h"
#include "gploc.h"

//#define SANTE_FE_TRAIL
#define LOS_ALTOS_TRAIL

#ifdef SANTE_FE_TRAIL
// Sante Fe Trail
#define MAX_CLOCK 600       // maximum number of operations

Control control = {
    1000,           // M - Number of individuals in each generation
    51,             // G - Number of generations to run
    6,              // Di - Maximum depth of S Expressions for an initial tree
    17,             // Dc - Maximum depth of S Expressions for a created tree
    .90,            // Pc - Probability of cross over
    .10,            // Pr - Probability of reproduction
    .90,            // Pip - Probability of cross over internal point
    NO_FUNCTIONS,
    NO_TERMINALS,
    0
};

#define MAXY 31
#define MAXX 31
char grid[MAXY+1][MAXX+1];
char food[][2] = {
    {1,0}, {2,0}, {3,0}, {3,1}, {3,2}, {3,3}, {3,4}, {3,5}, {4,5}, {5,5}, 
    {6,5}, {8,5}, {9,5}, {10,5}, {11,5}, {12,5}, {12,6}, {12,7}, {12,8},
    {12,9}, {12,11}, {12,12}, {12,13}, {12,14}, {12,17}, {12,18}, {12,19},
    {12,20}, {12,21}, {12,22}, {12,23}, {11,24}, {10,24}, {9,24}, {8,24},
    {7,24}, {4,24}, {3,24}, {1,25}, {1,26}, {1,27}, {1,28}, {2,30}, {3,30},
    {4,30}, {5,30}, {7,29}, {7,28}, {8,27}, {9,27}, {10,27}, {11,27}, {12,27},
    {13,27}, {14,27}, {16,26}, {16,25}, {16,24}, {16,21}, {16,20}, {16,19},
    {16,18}, {17,15}, {20,14}, {20,13}, {20,10}, {20,9}, {20,8}, {20,7}, 
    {21,5}, {22,5}, {24,4}, {24,3}, {25,2}, {26,2}, {27,2}, {29,3}, {29,4},
    {29,6}, {29,9}, {29,12}, {28,14}, {27,14}, {26,14}, {23,15}, {24,18}, 
    {27,19}, {26,22}, {23,23}, {-1,-1},
};
#endif

#ifdef LOS_ALTOS_TRAIL
#define MAX_CLOCK 3000       // maximum number of operations

Control control = {
    2000,           // M - Number of individuals in each generation
    51,                // G - Number of generations to run
    6,              // Di - Maximum depth of S Expressions for an initial tree
    17,             // Dc - Maximum depth of S Expressions for a created tree
    .90,            // Pc - Probability of cross over
    .10,            // Pr - Probability of reproduction
    .90,            // Pip - Probability of cross over internal point
    NO_FUNCTIONS,
    NO_TERMINALS,
    .16           // GRc
};

#define MAXY 66
#define MAXX 43

char grid[MAXY+1][MAXX+1];
char food[][2] = {
    {1,0}, {2,0}, {3,0}, {3,1}, {3,2}, {3,3}, {3,4}, {3,5}, {4,5}, 
    {5,5}, {6,5}, {8,5}, {9,5}, {10,5}, {11,5}, {12,5}, {12,6}, {12,7}, 
    {12,8}, {12,9}, {12,11}, {12,12}, {12,13}, {12,14}, {12,17}, {12,18}, 
    {12,19}, {12,20}, {12,21}, {12,22}, {12,23}, {11,24}, {10,24}, {9,24}, 
    {8,24}, {7,24}, {4,24}, {3,24}, {1,25}, {1,26}, {1,27}, {1,28}, 
    {1,29}, {1,30}, {1,31}, {2,33}, {3,33}, {4,33}, {5,33}, {7,32}, {7,31},
    {8,30}, {9,30}, {10,30}, {11,30}, {12,30}, {17,30}, {18,30}, {20,29},
    {20,28}, {20,27}, {20,26}, {20,25}, {20,24}, {20,21}, {20,20}, {20,19},
    {20,18}, {21,15}, {24,14}, {24,13}, {24,10}, {24,9}, {24,8}, {24,7},
    {25,5}, {26,5}, {27,5}, {28,5}, {29,5}, {30,5}, {31,5}, {32,5}, {33,5}, 
    {34,5}, {35,5}, {36,5}, {38,4}, {38,3}, {39,2}, {40,2}, {41,2}, {43,3}, 
    {43,4}, {43,6}, {43,9}, {43,12}, {42,14}, {41,14}, {40,14}, {37,15}, 
    {38,18}, {41,19}, {40,22}, {37,23}, {37,24}, {37,25}, {37,26}, {37,27}, 
    {37,28}, {37,29}, {37,30}, {37,31}, {37,32}, {37,33}, {37,34}, {35,34}, 
    {35,35}, {35,36}, {35,37}, {35,38}, {35,39}, {35,40}, {35,41}, {35,42}, 
    {34,42}, {33,42}, {32,42}, {31,42}, {30,42}, {30,44}, {30,45}, {30,46}, 
    {30,47}, {30,48}, {30,49}, {28,50}, {28,51}, {28,52}, {28,53}, {28,54}, 
    {28,55}, {28,56}, {27,56}, {26,56}, {25,56}, {24,56}, {23,56}, {22,58}, 
    {22,59}, {22,60}, {22,61}, {22,62}, {22,63}, {22,64}, {22,65}, {22,66}, 
    {-1,-1}
};
#endif

long antX, antY;
short antXD, antYD;
short eatCount;
short clock;
short nPellets;

typedef struct Fitness {
// base values
    long long lng_nfr;
    long long lng_n;
    long long lng_a;
    long long lng_raw;
              
// derived values
    short r;
    short s;
} Fitness;

//
// private functions 
//
void clearGrid() {
    short i,j;
    for (j = 0; j <= MAXY; j++) {
        for (i = 0; i <= MAXY; i++) {
            grid[j][i] = 0;
        }
    }
}

void printGrid(const char* label, int run_number, int gen);
void prepareRun() {
    clearGrid();
    int j = 0;
    short x,y;
    while (food[j][0] > -1) {
        x = food[j][0];
        y = food[j][1];
        grid[y][x] = 1;
        j++;
    }
    nPellets = j;
    antX = 0;
    antY = 0;
    antXD = 1;
    antYD = 0;
    eatCount = 0;
    clock = 0;
}

char computeFitness(Tree* t) {
    if (!t->fitness) 
        t->fitness = malloc(sizeof (Fitness));
    Fitness* f = (Fitness*) t->fitness;
    t->hits = f->r = eatCount;
    initBaseFitness((BaseFitness*)f, f->r);
                                        // average over generation
    f->s = nPellets - eatCount;
    double a = 1.0 / (1.0 + (double)f->s);
    f->lng_a = float_to_int(a);

    if (f->s == 0)
        return TRUE;        // have a winner!
        
    return FALSE;
}

///
/// public functions - called by gp engine (gp.c))
///
void initLoc() { }

void treeResultHeader(short gen) {
    printf("n_pellets=%d\n", nPellets);
    if (gen < 0) {
        printf("%4s %4s %6s %6s %6s %8s %8s %8s\n", 
            "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr");
        printf("%4s %4s %6s %6s %6s %8s %8s %8s\n", 
            "---", "---", "-----", "---", "---", "------", "------", "------");
    }
    else {
        printf("%6s %4s %4s %6s %6s %6s %8s %8s %8s %6s\n", 
            "gen", "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr", "avgRawF");
        printf("%6s %4s %4s %6s %6s %6s %8s %8s %8s\n", 
            "----", "---", "---", "-----", "---", "---", "------", "------", "------");
    }
}

void reportTreeResult(Tree* t, short i, short gen, double avgRawF) {
    assert(i == t->tfid);
    Fitness* f = (Fitness*)t->fitness;
    if (gen < 0) 
        printf("%4d %4d %6d %6d %6d %6.6f %6.6f %6.6f\n", 
           t->tfid, t->tcid, t->hits, f->r, f->s,
           intToFloat(f->lng_a), intToFloat(f->lng_n), intToFloat(f->lng_nfr));
    else
        printf("%6d %4d %4d %6d %6d %6d %6.6f %6.6f %6.6f %6.2f\n", 
           gen, t->tfid, t->tcid, t->hits, f->r, f->s,
           intToFloat(f->lng_a), intToFloat(f->lng_n), intToFloat(f->lng_nfr), 
           avgRawF);
}

void printGrid(const char* label, int run_number, int gen) {
    short x,y;
    if (run_number > 0) {
        printf(label, run_number, gen);
    }
    else {
        printf(label);
    }
    printf("---------------------------------------------------------------------------------------\n");

    for (y = 0; y <= MAXY ; y++) {
        for (x = 0; x <= MAXX ; x++) {
            if (grid[y][x] == 1) 
                printf("X ");           // Food Not Found
            else if (grid[y][x] == 2) 
                printf("@ ");           // Food Found and Eaten
            else if (grid[y][x] == 3) 
                printf("O ");           // Cell Visited but no food
            else
                printf(". ");           // Cell Not Visited and no food
        }
        printf("\n");
    }
    printf("---------------------------------------------------------------------------------------\n");
}

jmp_buf timeout;
void execTrees(TreeSet* trees, int run_number, int gen) {
    Tree** t = trees->treeArray;
    short size = trees->n;
    short i;
    for (i = 0; i < size; i++) {
        if (t[i]->tcid != i)
            printf("tp: tcid=%d\n", t[i]->tcid);

        assert(t[i]->tcid == i);
        prepareRun();
        setjmp(timeout);
        while(clock < MAX_CLOCK) {
            execNode(t[i]->root);
        }
        if (computeFitness(t[i])) {
            Tree* tree = t[i];
            reportTreeResult(tree, tree->tfid, -1, -1.0);
            printGrid("Have Winner! - Run# %d Gen# %d\n", run_number, gen);
            trees->haveWinner = TRUE;
            trees->winningIndex = i;
            break;
        }
    }
}

void execSingleTree(Tree* tree) {
    prepareRun();
    printGrid("Before Run\n", -1, -1);
    setjmp(timeout);
    while(clock < MAX_CLOCK) {
        execNode(tree->root);
    }
    if (computeFitness(tree)) {
        printf("Have Winner!\n");
    }
    printGrid("After Run\n", -1, -1);
}

void clockTic() {
    if (++clock > MAX_CLOCK)
        longjmp(timeout, 0); 
}

void terminate() {
    clock = MAX_CLOCK;
    clockTic();
}

// Terminal Implementations
gptype move() {
    clockTic();
    antX += antXD;
    antY += antYD;
    if (antX < 0 || antX > MAXX || antY < 0 || antY > MAXY) {
        return 0;
    }

    if (grid[antY][antX] == 1) {
        grid[antY][antX] = 2;   // yummy!
        if (++eatCount == nPellets)
            terminate();
        assert(eatCount < nPellets);
    }
    else if (grid[antY][antX] == 0) {
        grid[antY][antX] = 3;   // been here
    }

    return 0;
}

gptype left() {
    short T = antXD;
    clockTic();
    antXD = antYD;
    antYD = -1 * T;
    return 0;
}

gptype right() {
    short T = antXD;
    clockTic();
    antXD = -1 * antYD;
    antYD = T;
    return 0;
}

// Function Implementations
gptype ifFoodAhead(FunctionNode* this) {
    assert(this->type == 'f');
// if food exists at next step execute branch 0 else branch 1
// ignore any return values
    short newX = antX + antXD;
    short newY = antY + antYD;

    if (newX < 0 || newX > MAXX || newY < 0 || newY > MAXY)
        return execNode(this->branch[1]);

    if (grid[newY][newX] == 1)
        return execNode(this->branch[0]);

    return execNode(this->branch[1]);
}

gptype progN2(FunctionNode* this) {
    assert(this->type == 'f');
// unconditionally execute branch 0 & 1
// ignore any return values
    execNode(this->branch[0]);
    return execNode(this->branch[1]);
}

gptype progN3(FunctionNode* this) {
    assert(this->type == 'f');
// unconditionally execute branch 0, 1 & 2
// ignore any return values
    execNode(this->branch[0]);
    execNode(this->branch[1]);
    return execNode(this->branch[2]);
}

Terminal terminal[] = {
    {
        't',
        TERM_MOVE,
        "move",
        move
    },
    {
        't',
        TERM_LEFT,
        "left",
        left
    },
    {
        't',
        TERM_RIGHT,
        "right",
        right
    },
    {
        't',
        -1,
        NULL,
        NULL
    },
};

Function function[] = {
    {
        FUNC_IF_FOOD_AHEAD,
        "if_food_ahead",
        2,
        ifFoodAhead
    },
    {
        FUNC_PROG_N2,
        "prog_n2",
        2,
        progN2
    },
    {
        FUNC_PROG_N3,
        "prog_n3",
        3,
        progN3
    },
    {
        -1,
        NULL,
        0,
        NULL
    }
};

//// sample trees /////
Tree* createWinnerSampleTree() {
    FunctionNode* f1 = newFunctionNode(FUNC_IF_FOOD_AHEAD);
    FunctionNode* f2 = newFunctionNode(FUNC_PROG_N3);
    FunctionNode* f3 = newFunctionNode(FUNC_PROG_N2);
    FunctionNode* f4 = newFunctionNode(FUNC_PROG_N2);
    FunctionNode* f5 = newFunctionNode(FUNC_IF_FOOD_AHEAD);
    FunctionNode* f6 = newFunctionNode(FUNC_PROG_N2);
    FunctionNode* f7 = newFunctionNode(FUNC_IF_FOOD_AHEAD);
    FunctionNode* f8 = newFunctionNode(FUNC_PROG_N2);

    setArg(f1, 1, (Node*)&terminal[TERM_MOVE]);
    setArg(f1, 2, (Node*)f2);

    setArg(f2, 1, (Node*)&terminal[TERM_LEFT]);
    setArg(f2, 2, (Node*)f3);
    setArg(f2, 3, (Node*)f4);

    setArg(f3, 1, (Node*)f5);
    setArg(f3, 2, (Node*)f6);

    setArg(f4, 1, (Node*)f7);
    setArg(f4, 2, (Node*)&terminal[TERM_MOVE]);

    setArg(f5, 1, (Node*)&terminal[TERM_MOVE]);
    setArg(f5, 2, (Node*)&terminal[TERM_RIGHT]);

    setArg(f6, 1, (Node*)&terminal[TERM_RIGHT]);
    setArg(f6, 2, (Node*)f8);

    setArg(f7, 1, (Node*)&terminal[TERM_MOVE]);
    setArg(f7, 2, (Node*)&terminal[TERM_LEFT]);

    setArg(f8, 1, (Node*)&terminal[TERM_LEFT]);
    setArg(f8, 2, (Node*)&terminal[TERM_RIGHT]);

    Tree* t = newTree((Node*)f1);
    t->tfid = t->tcid = 0;
    return t;
}
