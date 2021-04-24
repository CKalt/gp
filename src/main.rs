mod gploc;
mod tree;
use gploc::*;
use tree::*;

/////////////////////// STUBS /////////////////////////////

fn push_tree(mut _trees: &TreeSet, mut _tree: Tree) -> () {
}

fn create_unique_tree(_trees: &TreeSet, _d: u16) -> Tree {
    Tree::new()
}

fn count_nodes(mut _trees: &TreeSet) -> () {

}
        
fn compute_normalized_fitness(mut _trees: &TreeSet) ->() {

}

fn sort_by_normalized_fitness(mut _trees: &TreeSet) ->() {

}

fn assign_nf_rankings(mut _trees: &TreeSet) ->() {

}
        
fn report_results(_gen: u16, _trees: &TreeSet) -> () {

}

fn select_tree(trees: &TreeSet) -> &Tree {
    assert!(trees.tree_vec.len() > 0);
    &trees.tree_vec[0]
}

fn count_tree_nodes(tree: Tree) -> Tree {
    tree
}

fn clone_tree(_tree: &Tree) -> Tree {
    Tree::new()
}
                    
fn perform_crossover(mut _t1: &Tree, mut _t2: &Tree) -> () {
}

fn new_tree_qualifies(_tree: &Tree) -> bool {
    true
}

fn print_tree(_tree : &Tree) {

}

/////////////////////// NOT STUBS /////////////////////////////

fn create_initial_population() -> TreeSet {
    // Following Koza's recipe, he calls "ramped half-and-half",
    // we will evenly produce population segments starting with 2
    // up to the maxium depth (CONTROL.Di) and alternate
    // between Full Method and Grow Method for S Expressions.
    //
    let trees = TreeSet::new();
    let seg = trees.tree_vec.len() as f64 / (CONTROL.Di as f64 - 1.0f64);
    let mut bdr = 0.0f64;

    for d in 2..=CONTROL.Di {
        bdr += seg;
        while trees.tree_vec.len() < bdr as usize {
            push_tree(&trees, create_unique_tree(&trees, d));
        }
    }

    // fill out to end in case there are "left-overs" due to rounding
    while trees.tree_vec.len() < trees.tree_vec.len() {
        push_tree(&trees, create_unique_tree(&trees, CONTROL.Di));
    }
    trees
}

fn use_reproduction(index: usize ) -> bool {
    let result = (index as f64 / CONTROL.M as f64) < CONTROL.Pr;
    result
}

fn run() -> Option<Tree> {
    let mut trees = create_initial_population();
    count_nodes(&trees);

    let mut gen = 0u16;
    while gen <= CONTROL.G && trees.winning_index == None {
        exec_trees(&trees);
        if trees.winning_index != None {
            break;
        }

        compute_normalized_fitness(&trees);
        sort_by_normalized_fitness(&trees);
        assign_nf_rankings(&trees);

        report_results(gen, &trees);

        if gen >= CONTROL.G {
            break;
        }

        let trees2 = TreeSet::new();
        for i in (0..CONTROL.M).step_by(2) {
            if use_reproduction(i) {
                // do reproduction
                let t = select_tree(&trees);
                push_tree(&trees2, count_tree_nodes(clone_tree(t)));
            }
            else {
                // do crossover
                let (mut nt1, mut nt2);
                loop {
                    let t1 = select_tree(&trees);
                    let t2 = select_tree(&trees);
                    nt1 = clone_tree(t1);
                    nt2 = clone_tree(t2);
                    perform_crossover(&nt1, &nt2);

                    if new_tree_qualifies(&nt1) && new_tree_qualifies(&nt2) {
                        break;
                    }
                }
                push_tree(&trees2, count_tree_nodes(nt1));

                if i+1 < CONTROL.M {
                    push_tree(&trees2, count_tree_nodes(nt2));
                }
            }
        }
        trees = trees2;
        gen += 1;
    }

    match trees.winning_index {
        Some(i) => {
            let winner = clone_tree(&trees.tree_vec[i]);
            Some(winner)
        },
        _ => None,
    }
}


fn main() {
    init_loc();
    let mut total_runs = 0i32;
    let winner = loop { // go until we have a winner
        total_runs += 1;
        println!("Run #{}", total_runs);
        if let Some(winner) = run() {
            break winner;
        }
    };
        
    println!("total_runs={}", total_runs);
    print_tree(&winner);
    exec_single_tree(&winner);
}
