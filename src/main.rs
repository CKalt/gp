extern crate gp;
mod gprun;
mod fitness;
mod tree;
mod tree_set;
mod control;
mod gprng;
mod util;

#[cfg(test)]
mod tests;

use gprun::*;
use tree::*;
use tree_set::*;
use control::*;

use gprng::GpRng;
use gprng::GpRngFactory;
use fitness::GpHits;

use BranchType::*;

fn report_results(rng: &mut GpRng, trees: &mut TreeSet,header_need: &mut bool,
        hits: &GpHits) -> () {
    if CONTROL.show_all_trees {
        println!("Generation {}", trees.gen);
        trees.print();
    }
    if CONTROL.show_all_tree_results {
        println!("gen {}", trees.gen);
        Tree::print_result_header(None, hits);
        for t in trees.tree_vec.iter() {
            t.print_result(None, -1.0);
        }
    }
    if CONTROL.show_best_tree_results {
        if *header_need {
            Tree::print_result_header(Some(trees.gen), hits);
            *header_need = false;
        }
        let i = trees.tree_vec.len()-1;
        let best_tree = &trees.tree_vec[i];
        best_tree.print_result(Some(trees.gen), trees.avg_raw_f);
    }
    if CONTROL.run_tests {
        run_tests(rng, trees);
    }
}

fn run(rng: &mut GpRng, run_number: i32) -> Option<Winner> {
    if CONTROL.show_controls {
        println!("M = {}, G = {}, D = {}", CONTROL.M, CONTROL.G, CONTROL.Di);
    }
    let mut trees = TreeSet::create_initial_population(rng);

    trees.count_nodes();
    trees.gen = 0u16;
    let mut header_need: bool = true;
    // when gen == CONTROL.G we are done because gen starts with 0
    // ie. we do not process gen == 51 if G==51
    while trees.gen < CONTROL.G && trees.winning_index == None {
        let hits = trees.exec_all();
        if trees.winning_index != None {
            break;
        }

        #[cfg(gpopt_trace="on")]
        println!("TP0:trace log for gen={}", trees.gen);

        trees.compute_normalized_fitness();
        trees.sort_by_normalized_fitness();
        report_results(rng, &mut trees, &mut header_need, &hits);

        if trees.gen >= CONTROL.G {
            break;
        }

        let trees2 = trees.breed_new_generation(rng);

        trees = trees2;
        trees.gen += 1;
    }

    // return some optional Winner or None
    match trees.winning_index {
        Some(i) => {
            let winner = Winner{
                tree: trees.tree_vec[i].clone(),
                run: run_number,
                gen: trees.gen,
                e: CONTROL.computational_effort(run_number, trees.gen),
            };
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

        let (fi, bt, fnode) = t.get_rnd_function_node_ref_i(rng);

        println!("\n--------\nRnd fi={} Subtree of {} Branch Type -->\n", fi,
            match bt {
                Result0 => "Result",
                FunctionDef0 => "Function Def",
            });
        fnode.print();

        let (_, tnode) = t.get_rnd_terminal_node_ref(rng);
        tnode.print();
        println!("\n--------");
    }
}

use gp::parsers::parse_func_call;
use gp::parsers::parse_sexpr;

fn main() {
    init_run();

    if true {
        assert_eq!(
            format!("{:?}", parse_sexpr("( ADD  X    Y )"
                    )), 
                    r#"Ok(("", Func("ADD", [Term("X"), Term("Y")])))"#);
        assert_eq!(
            format!("{:?}", parse_func_call("  ABCD123 X1 (Y 1) Z  ")),
                    r#"Ok(("  ", Func("ABCD123", [Term("X1"), Func("Y", [Term("1")]), Term("Z")])))"#);
    }

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
        let rc = RunContext::new();
        println!("run={}, gen={}", winner.run, winner.gen);
        println!("Computational effort (e) = {}", winner.e);
        winner.print_result();
        rc.print_run_illustration(&format!("Have Winner! - Run# {} Gen# {}", run_number,
            winner.gen));
        winner.tree.print_exec_one();
    } else {
        println!("Exceeded CONTROL.R ({}) runs without finding a winner.",
            CONTROL.R);
    }
}
