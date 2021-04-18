mod gploc;
mod tree;

#[allow(non_snake_case)]
struct Control {
    M : i16,                // Number of individuals in each generation
    G : i16,                // Number of generations to run
    Di : i16,               // Maximum depth of S Expressions for an initial tree
    Dc : i16,               // Maximum depth of S Expressions for a created tree
    Pc : f64,               // Probability of cross over
    Pr : f64,               // Probability of reproduction
    Pip : f64,              // Probability of cross over internal point
    no_functions : i16,
    no_terminals : i16,
    GRc : f64,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            //                   GRc = 32/2^n
    no_fitnessCases : i16
}
#[warn(non_snake_case)]

use gploc::*;
use tree::*;

fn run(control: &Control) -> Option<Tree> {
    return Some(Tree { });

    let trees = create_initial_population();
    while (gen <= control.G && !trees.have_winner) {
        exec_trees(trees);
        if (trees.have_winner) {
            break;
        }

        compute_normalized_fitness(trees);
        sort_by_normalized_fitness(trees);
        assign_NF_rankings(trees);

        report_results(gen, trees);

        if (gen >+ control.G) {
            break;
        }

        let trees2 = TreeSet { };
        for i in 0..control.M {

        }


    }
}

fn print_tree(_tree : &Tree) {

}

fn main() {
    let control = Control { };
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
