mod gprun;
mod tree;
mod control;
mod gprng;
mod util;
#[cfg(gpopt_choice_logging="write")]
mod choice_logging;

use gprun::*;
use tree::*;
use control::*;
#[cfg(gpopt_choice_logging="write")]
use choice_logging::*;
use Node::*;

#[cfg(not(gpopt_rng="FileStream"))]
use rand::Rng;

use gprng::GpRng;
use gprng::GpRngFactory;
use std::mem;

fn push_tree(trees: &mut TreeSet, mut tree: Tree) {
    let index = trees.tree_vec.len();
    assert!(index < CONTROL.M);
    tree.tcid = index;
    trees.tree_vec.push(tree);
}

fn gen_tree(rng: &mut GpRng, generate_method : GenerateMethod,
        d: u16) -> Tree {
    match generate_method {
        GenerateMethod::Full => gen_tree_full_method(rng, d),
        GenerateMethod::Grow => gen_tree_grow_method(rng, d),
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

fn create_unique_tree(rng: &mut GpRng, trees: &TreeSet,
        mut d: u16) -> Tree {
    let mut i = 1u16;
    let gen_method = 
        if (trees.tree_vec.len() % 2) == 0 {
            GenerateMethod::Full
        }
        else {
            GenerateMethod::Grow
        };

    let mut t = gen_tree(rng, gen_method, d);

    while tree_match_exists(trees, &t) {
        i += 1;
        if i > 11 {
            d += 1;
            i = 1;
        }
        t = gen_tree(rng, gen_method, d);
    }

    return t;
}

fn gen_tree_full_method(rng: &mut GpRng, depth: u16) -> Tree {
    let mut root = FunctionNode::new_rnd(rng);
    gen_tree_full_method_r(rng, &mut root, 2, depth);
    Tree::new(root)
}

fn gen_tree_full_method_r(rng: &mut GpRng,
        func_node: &mut FunctionNode, level: u16, depth: u16) {
    if level >= depth {
        for i in 0..func_node.fnc.arity {
            let rnd_tref = Terminal::get_rnd_ref(rng); // Always a Terminal Node
            func_node.set_arg(i, TNode(rnd_tref));
        }
    }
    else {
        let c_depth = level+1;
        for i in 0..func_node.fnc.arity {
            let rnd_fn = FunctionNode::new_rnd(rng); // Always a Funciton Node
            let node: &mut Node = func_node.set_arg(i, FNode(rnd_fn));

            match *node {
                FNode(ref mut fn_ref) =>
                    gen_tree_full_method_r(rng, fn_ref, c_depth, depth),
                _ => panic!("expected FunctionNode"),
            }
        }
    }
}

fn gen_tree_grow_method(rng: &mut GpRng, depth: u16) -> Tree {
    let mut root = FunctionNode::new_rnd(rng);
    gen_tree_grow_method_r(rng, &mut root, 2, depth);

    Tree::new(root)
}
    
fn gen_tree_grow_method_r(rng: &mut GpRng, 
        func_node: &mut FunctionNode, level: u16, depth: u16) {
    if level >= depth {
        for i in 0..func_node.fnc.arity {
            let rnd_tref = Terminal::get_rnd_ref(rng); // Always a Terminal Node
            func_node.set_arg(i, TNode(rnd_tref));
        }
    }
    else {
        let c_depth = level+1;
        for i in 0..func_node.fnc.arity {
            let rnd_ft_node = Node::new_rnd(rng); // Either a Function or Terminal Node
            let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
            if let FNode(ref mut fn_ref) = node {
                gen_tree_grow_method_r(rng, fn_ref, c_depth, depth);
            }
        }
    }
}

fn report_results(rng: &mut GpRng, trees: &mut TreeSet,header_need: &mut bool,
        n_pellets: &u16) -> () {
    if CONTROL.show_all_trees {
        println!("Generation {}", trees.gen);
        trees.print();
    }
    if CONTROL.show_all_tree_results {
        println!("gen {}", trees.gen);
        tree_result_header(None, n_pellets);
        for (i,t) in trees.tree_vec.iter().enumerate() {
            report_tree_result(t, Some(i), None, -1.0);
        }
    }
    if CONTROL.show_best_tree_results {
        if *header_need {
            tree_result_header(Some(trees.gen), n_pellets);
            *header_need = false;
        }
        let i = trees.tree_vec.len()-1;
        let best_tree = &trees.tree_vec[i];
        report_tree_result(best_tree, Some(i), Some(trees.gen), trees.avg_raw_f);
    }
    if CONTROL.run_tests {
        run_tests(rng, trees);
    }
}

// rnd_internal_point - randomly decides whether to do crossover at an
// internal point (function) or terminal based on control Pip value.
#[cfg(gpopt_choice_logging="write")]
fn rnd_internal_point(rng: &mut GpRng) -> bool {
    #[cfg(gpopt_rng="FileStream")]
    let num: GpFloat = rng.gen_float();

    #[cfg(not(gpopt_rng="FileStream"))]
    let num: GpFloat = rng.gen_range(0.0..1.0);

    let result = num < CONTROL.Pip;
    choice_log(1, if result { "1" } else { "0" });

    result
}
#[cfg(gpopt_choice_logging="off")]
fn rnd_internal_point(rng: &mut GpRng) -> bool {
    #[cfg(gpopt_rng="FileStream")]
    let num: GpFloat = rng.gen_float();

    #[cfg(not(gpopt_rng="FileStream"))]
    let num: GpFloat = rng.gen_range(0.0..1.0);

    num < CONTROL.Pip // if Pip is .90 then true for all values less than .90.
}

fn perform_crossover(rng: &mut GpRng,
        t1: &mut Tree, t2: &mut Tree) {
    assert_ne!(t1.num_terminal_nodes, None);
    assert_ne!(t2.num_terminal_nodes, None);

    let swap_target1 : &mut Node =
        if t1.num_function_nodes.unwrap() > 0 && rnd_internal_point(rng) {
            t1.get_rnd_function_node_ref(rng)
        } else {
            t1.get_rnd_terminal_node_ref(rng)
        };

    let swap_target2 : &mut Node =
        if t2.num_function_nodes.unwrap() > 0 && rnd_internal_point(rng) {
            t2.get_rnd_function_node_ref(rng)
        } else {
            t2.get_rnd_terminal_node_ref(rng)
        };

    mem::swap(swap_target1, swap_target2);
}

fn create_initial_population(rng: &mut GpRng) -> TreeSet {
    // Following Koza's recipe, he calls "ramped half-and-half",
    // we will evenly produce population segments starting with 2
    // up to the maxium depth (CONTROL.Di) and alternate
    // between Full Method and Grow Method for S Expressions.
    let mut trees = TreeSet::new(0);
    let seg = CONTROL.M as GpFloat / (CONTROL.Di as GpFloat - 1.0);
    let mut bdr = 0.0;

    #[cfg(gpopt_trace="on")]
    println!("TP001:create_init_pop start");
    for d in 2..=CONTROL.Di {
        bdr += seg;
        while trees.tree_vec.len() < bdr as usize &&
              trees.tree_vec.len() < CONTROL.M {
            let mut new_tree = create_unique_tree(rng, &trees, d);
            new_tree.count_nodes();

            push_tree(&mut trees, new_tree);

            #[cfg(gpopt_trace="on")]
            trees.tree_vec[trees.tree_vec.len()-1].print();
        }
    }

    // fill out to end in case there are "left-overs" due to rounding
    while trees.tree_vec.len() < CONTROL.M {
        let mut new_tree = create_unique_tree(rng, &trees, CONTROL.Di);
        new_tree.count_nodes();
        push_tree(&mut trees, new_tree);
            
        #[cfg(gpopt_trace="on")]
        trees.tree_vec[trees.tree_vec.len()-1].print();
    }
    #[cfg(gpopt_trace="on")]
    println!("TP002:create_init_pop done");

    trees
}

//fn use_reproduction(index: usize ) -> bool {
//    let float_index_ratio = index as GpFloat / CONTROL.M as GpFloat;
//    let int_index_ratio = Fitness::float_to_int(
//        float_index_ratio
//    );
//    let int_control_ratio = Fitness::float_to_int(CONTROL.Pr);
//    let result = int_index_ratio < int_control_ratio;
//    result
//}

#[cfg(gpopt_choice_logging="write")]
fn use_reproduction(index: usize ) -> bool {
    let result = (index as GpFloat / CONTROL.M as GpFloat) < CONTROL.Pr;
    choice_log(9, if result { "1" } else { "0" });
    result
}
#[cfg(gpopt_choice_logging="off")]
fn use_reproduction(index: usize ) -> bool {
    let result = (index as GpFloat / CONTROL.M as GpFloat) < CONTROL.Pr;
    result
}

fn run(rng: &mut GpRng, run_number: i32) -> Option<Tree> {
    if CONTROL.show_controls {
        println!("M = {}, G = {}, D = {}", CONTROL.M, CONTROL.G, CONTROL.Di);
    }

    let mut trees = create_initial_population(rng);
    trees.count_nodes();
    trees.gen = 0u16;
    let mut header_need: bool = true;
    while trees.gen <= CONTROL.G && trees.winning_index == None {
        let n_pellets = exec_trees(&mut trees, run_number);
        if trees.winning_index != None {
            break;
        }

        #[cfg(gpopt_trace="on")]
        println!("TP0:trace log for gen={}", trees.gen);

        trees.compute_normalized_fitness()
             .sort_by_normalized_fitness();

//        #[cfg(gpopt_trace="on")]
//        if trees.gen == 1 {
//            panic!("pause");
//        }
//        else {
//            println!("at pause point gen={} skipping until gen=1 here.", trees.gen);
//        }

        report_results(rng, &mut trees, &mut header_need, &n_pellets);


//if trees.gen == 1 {
//    println!("pause");
//}
        if trees.gen >= CONTROL.G {
            break;
        }

        let mut trees2 = TreeSet::new(trees.gen);
        #[cfg(gpopt_trace="on")]
        println!("TP003:breed start");
        while trees2.tree_vec.len() < CONTROL.M {
            if use_reproduction(trees2.tree_vec.len()) {
                // do reproduction
                #[cfg(gpopt_trace="on")]
                println!("TP003.1:breeding reproduction branch");

                let mut t = trees.select_tree(rng).clone();
                t.clear_node_counts();
                t.count_nodes();
                push_tree(&mut trees2, t);

                #[cfg(gpopt_trace="on")]
                trees2.tree_vec[trees2.tree_vec.len()-1].print();
            }
            else {
                // do crossover
                #[cfg(gpopt_trace="on")]
                println!("TP003.1:breeding crossover branch");

                let (mut nt1, mut nt2);
                loop {
                    let t1 = trees.select_tree(rng);
                    let t2 = trees.select_tree(rng);
                    nt1 = t1.clone();
                    nt2 = t2.clone();
                    perform_crossover(rng, &mut nt1, &mut nt2);
                    if nt1.qualifies() && nt2.qualifies() {
                        break;
                    }
                }
                nt1.clear_node_counts();
                nt1.count_nodes();

                push_tree(&mut trees2, nt1);

                #[cfg(gpopt_trace="on")]
                trees2.tree_vec[trees2.tree_vec.len()-1].print();

                if trees2.tree_vec.len() < CONTROL.M {
                    nt2.clear_node_counts();
                    nt2.count_nodes();
                    push_tree(&mut trees2, nt2);

                    #[cfg(gpopt_trace="on")]
                    trees2.tree_vec[trees2.tree_vec.len()-1].print();
                }
            }
        }

        #[cfg(gpopt_trace="on")]
        println!("TP004:breed done");

        assert_eq!(trees2.tree_vec.len(), CONTROL.M);
        trees = trees2;
        trees.gen += 1;
    }

    match trees.winning_index {
        Some(i) => {
            let winner = trees.tree_vec[i].clone();
            Some(winner)
        },
        _ => None,
    }
}

fn run_tests(rng: &mut GpRng, trees: &mut TreeSet) {
    trees.count_nodes();
    println!("Testing getRndFunctionIndex and getRndTerminalIndex");
    for i in 1..=10 {
        println!("=============\n| Trial {:3} |\n=============", i);
        let t = trees.get_rnd_tree(rng);
        t.print();

        let (fi, fnode) = t.get_rnd_function_node_ref_i(rng);
        println!("\n--------\nRnd fi={} Subtree -->\n", fi);
        fnode.print();

        let tnode = t.get_rnd_terminal_node_ref(rng);
        tnode.print();
        println!("\n--------");
    }
}

fn main() {
    init_run();
    #[cfg(gpopt_choice_logging="write")]
    init_choice_log();

    let mut run_number = 0i32;
    let mut rng = GpRngFactory::new();
    let opt_winner = loop { // go until we have a winner
        run_number += 1;
        println!("Run #{}", run_number);
        if let Some(winner) = run(&mut rng, run_number) {
            break Some(winner);
        }
        else if run_number == CONTROL.R && CONTROL.R != 0 {
            break None;
        }
    };
        
    if let Some(mut winner) = opt_winner {
        println!("run_number={}", run_number);
        winner.print();
        exec_single_tree(&mut winner);
    } else {
        println!("Exceeded CONTROL.R ({}) runs without finding a winner.",
            CONTROL.R);
    }
}
