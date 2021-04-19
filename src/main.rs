mod gploc;
mod tree;
use gploc::*;
use tree::*;

// gp.h stuff here
#[allow(non_snake_case)]
struct Control {
    M:   i16,               // Number of individuals in each generation
    G:   i16,               // Number of generations to run
    Di:  i16,               // Maximum depth of S Expressions for an initial tree
    Dc:  i16,               // Maximum depth of S Expressions for a created tree
    Pc:  f64,               // Probability of cross over
    Pr:  f64,               // Probability of reproduction
    Pip:  f64,              // Probability of cross over internal point
    no_functions:  i16,
    no_terminals:  i16,
    GRc:  f64,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            //                   GRc = 32/2^n
    no_fitnessCases:  i16
}
#[warn(non_snake_case)]

enum Node {
    terminal,
    function
} Node;

struct TreeSet {
    n: i32,              // current size (population)
    size: usize,         // max size (population)
    have_winner: bool,   // A winner is an individual with a best possible
                         // score
    winning_index: i32m
    avgRawF: f64,
    tree_array: Vec<Tree>,
    tagArray: Vec<bool>, // used for temporary individual state marking
                         // e.g. used during tournament selection to insure
                         // unique individuals compete
}

(UC)

struct Tree {


}

//////////////// end of gp.h stuff ///////

fn run(control: &Control) -> Option<Tree> {
    let trees = create_initial_population();
    count_nodes(trees);

    let gen = 0i16;
    while gen <= control.G && !trees.have_winner {
        exec_trees(trees);
        if (trees.have_winner) {
            break;
        }

        compute_normalized_fitness(trees);
        sort_by_normalized_fitness(trees);
        assign_NF_rankings(trees);

        report_results(gen, trees);

        if gen >+ control.G {
            break;
        }

        let trees2 = TreeSet { };
        for i in (0..control.M).step_by(2) {
            if use_reproduction(i) {
                // do reproduction
                let t = select_tree(trees);
                push_tree(trees2, count_tree_nodes(clone_tree(t)));
            }
            else {
                // do crossover
                let (mut nt1, mut nt2);
                loop {
                    let t1 = select_tree(trees);
                    let t2 = select_tree(trees);
                    nt1 = clone_tree(t1);
                    nt2 = clone_tree(t2);
                    perform_crossover(nt1, nt2);

                    if new_tree_qualifies(nt1) && new_tree_qualifies(nt2) {
                        break;
                    }
                }
                push_tree(trees2, count_tree_nodes(nt1));

                if i+1 < control.M {
                    push_tree(trees2, count_tree_nodes(nt2));
                }
            }
        }
        trees = trees2;
        gen += 1;
    }

    if trees.have_winner {
        let winner = clone_tree(trees.tree_array[trees.winning_index]);
        Some(winner)
    }
    else {
        None
    }
}

fn print_tree(_tree : &Tree) {

}

fn main() {
    let control = Control { 
        M: 1000,             // M - Number of individuals in each generation
        G: 51,               // G - Number of generations to run
        Di: 6,               // Di - Maximum depth of S Expressions for an initial tree
        Dc: 17,              // Dc - Maximum depth of S Expressions for a created tree
        Pc: .90,             // Pc - Probability of cross over
        Pr: .10,             // Pr - Probability of reproduction
        Pip: .90,            // Pip - Probability of cross over internal point
        no_functions: 3,
        no_terminals: 3,
        GRc: 0
    };

    init_loc();
    let mut total_runs = 0i32;
    let winner = loop { // go until we have a winner
        total_runs += 1;
        println!("Run #{}", total_runs);
        if let Some(winner) = run(&control) {
            break winner;
        }
    };
        
    println!("total_runs={}", total_runs);
    print_tree(&winner);
    exec_single_tree(&winner);
}
