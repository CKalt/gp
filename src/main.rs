mod gprun;
mod tree;
mod control;

use gprun::*;
use tree::*;
use control::*;
use Node::*;

fn push_tree(trees: &mut TreeSet, tree: Tree) -> () {
    assert!(trees.tree_vec.len() < CONTROL.M);
    trees.tree_vec.push(tree);
}

fn gen_tree(generate_method : GenerateMethod, d: u16) -> Tree {
    match generate_method {
        GenerateMethod::Full => gen_tree_full_method(d),
        GenerateMethod::Grow => gen_tree_grow_method(d),
    }
}

fn tree_match(t1: &Tree, t2: &Tree) -> bool {
    t1.root.deep_match(&t2.root)
}

fn tree_match_exists(trees: &TreeSet, target_tree: &Tree) -> bool {
    for tree in trees.tree_vec.iter() {
        if tree_match(tree, target_tree) {
            return true;
        }
    }
    return false;
}

fn create_unique_tree(trees: &TreeSet, mut d: u16) -> Tree {
    let mut i = 1u16;
    let gen_method = 
        if (trees.tree_vec.len() % 2) == 0 {
            GenerateMethod::Full
        }
        else {
            GenerateMethod::Grow
        };

    let mut t = gen_tree(gen_method, d);

    while tree_match_exists(trees, &t) {
        i += 1;
        if i > 11 {
            d += 1;
            i = 1;
        }
        t = gen_tree(gen_method, d);
    }
    return t;
}

fn gen_tree_full_method(depth: u16) -> Tree {
    let mut root = FunctionNode::new_rnd();
    gen_tree_full_method_r(&mut root, 2, depth);
    Tree::new(root)
}

fn gen_tree_full_method_r(func_node: &mut FunctionNode, level: u16, depth: u16) {
    if level >= depth {
        for i in 0..func_node.fnc.arity {
            let rnd_tref = Terminal::get_rnd_ref(); // Always a Terminal Node
            func_node.set_arg(i, TNode(rnd_tref));
        }
    }
    else {
        let c_depth = level+1;
        for i in 0..func_node.fnc.arity {
            let rnd_fn = FunctionNode::new_rnd(); // Always a Funciton Node
            let node: &mut Node = func_node.set_arg(i, FNode(rnd_fn));

            match *node {
                FNode(ref mut fn_ref) =>
                    gen_tree_full_method_r(fn_ref, c_depth, depth),
                _ => panic!("expected FunctionNode"),
            }
        }
    }
}

fn gen_tree_grow_method(depth: u16) -> Tree {
    let mut root = FunctionNode::new_rnd();
    gen_tree_grow_method_r(&mut root, 2, depth);

    Tree::new(root)
}
    
fn gen_tree_grow_method_r(func_node: &mut FunctionNode, level: u16,
    depth: u16) {
    if level >= depth {
        for i in 0..func_node.fnc.arity {
            let rnd_tref = Terminal::get_rnd_ref(); // Always a Terminal Node
            func_node.set_arg(i, TNode(rnd_tref));
        }
    }
    else {
        let c_depth = level+1;
        for i in 0..func_node.fnc.arity {
            let rnd_ft_node = Node::new_rnd(); // Either a Function or Terminal Node
            let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
            if let FNode(ref mut fn_ref) = node {
                gen_tree_grow_method_r(fn_ref, c_depth, depth);
            }
        }
    }
}

/////////////////////// STUBS /////////////////////////////
fn report_results(_gen: u16, _trees: &TreeSet) -> () {

}

fn select_tree(trees: &TreeSet) -> &Tree {
    assert!(trees.tree_vec.len() > 0);
    &trees.tree_vec[0]
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
    let mut trees = TreeSet::new();
    let seg = trees.tree_vec.len() as f64 / (CONTROL.Di as f64 - 1.0f64);
    let mut bdr = 0.0f64;

    for d in 2..=CONTROL.Di {
        bdr += seg;
        while trees.tree_vec.len() < bdr as usize {
            let new_tree = create_unique_tree(&trees, d);
            push_tree(&mut trees, new_tree);
        }
    }

    // fill out to end in case there are "left-overs" due to rounding
    while trees.tree_vec.len() < trees.tree_vec.len() {
        let new_tree = create_unique_tree(&trees, CONTROL.Di);
        push_tree(&mut trees, new_tree);
    }
    trees
}

fn use_reproduction(index: usize ) -> bool {
    let result = (index as f64 / CONTROL.M as f64) < CONTROL.Pr;
    result
}

fn run() -> Option<Tree> {
    let mut trees = create_initial_population();
    trees.count_nodes();

    let mut gen = 0u16;
    while gen <= CONTROL.G && trees.winning_index == None {
        exec_trees(&mut trees);
        if trees.winning_index != None {
            break;
        }

        trees.compute_normalized_fitness()
             .sort_by_normalized_fitness();

        report_results(gen, &mut trees);

        if gen >= CONTROL.G {
            break;
        }

        let mut trees2 = TreeSet::new();
        for i in (0..CONTROL.M).step_by(2) {
            if use_reproduction(i) {
                // do reproduction
                let mut t = select_tree(&trees).clone();
                t.count_nodes();
                push_tree(&mut trees2, t);
            }
            else {
                // do crossover
                let (mut nt1, mut nt2);
                loop {
                    let t1 = select_tree(&trees);
                    let t2 = select_tree(&trees);
                    nt1 = t1.clone();
                    nt2 = t2.clone();
                    perform_crossover(&nt1, &nt2);

                    if new_tree_qualifies(&nt1) && new_tree_qualifies(&nt2) {
                        break;
                    }
                }
                nt1.count_nodes();
                push_tree(&mut trees2, nt1);

                if i+1 < CONTROL.M {
                    nt2.count_nodes();
                    push_tree(&mut trees2, nt2);
                }
            }
        }
        trees = trees2;
        gen += 1;
    }

    match trees.winning_index {
        Some(i) => {
            let winner = trees.tree_vec[i].clone();
            Some(winner)
        },
        _ => None,
    }
}


fn main() {
    init_run();
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
