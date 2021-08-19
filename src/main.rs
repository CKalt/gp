extern crate gp;
extern crate lazy_static;

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

use std::fs::OpenOptions;
use std::path::Path;
use std::io::prelude::*;

    
#[test]
fn test_tee_to_run_log() {
    use std::fs::File;
    use std::io::BufReader;
    use std::fs::remove_file;
    use sscanf::scanf;

    // remove file if it already exists
    if Path::new(CONTROL.run_log_file).exists() {
        if let Err(e) = remove_file(CONTROL.run_log_file) {
            panic!("failed remove file: {}", e);
        }
    }
    let win: i32 = 5;
    let run: i32 = 10;
    let gen: u16 = 23;
    let e: u64 = 233235235;

    tee_to_run_log(win, run, gen, e);

    // Read
    let mut file = BufReader::new(File::open(CONTROL.run_log_file).unwrap());

    // read the first line and extract the number from it
    let mut input = String::new();
    file.read_line(&mut input).unwrap();
    input.truncate(input.len() - 1);            // remove new line
    let parsed = scanf!(input, "win={}, run={}, gen={}, effort={}",
        i32, i32, u16, u64);

    assert_eq!(parsed, Some((win, run, gen, e)));
}

fn tee_to_run_log(winner: &Winner) {
    let msg = format!("win={}, run={}, gen={}, effort={}, S={}", 
        winner.win, winner.run, winner.gen, winner.e,
        winner.structural_complexity());

    let create_bool = !Path::new(CONTROL.run_log_file).exists();
    let mut file =
         OpenOptions::new()
            .write(true)
            .create(create_bool)
            .append(!create_bool)
            .open(CONTROL.run_log_file)
            .unwrap();

    println!("{}", msg);
    if let Err(_e) = writeln!(file, "{}", msg) {
        panic!("Couldn't write to file");
    }
}

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
                win: 0,
                tree: trees.tree_vec[i].clone(),
                run: run_number,
                gen: trees.gen,
                e: EvalCount::num_evals_for_winner()
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
                ResultProducing => "Result".to_string(),
                FunctionDefining(adf_num) => 
                    format!("Function Def ADF{}", adf_num),
            });
        fnode.print();

        let (_, tnode) = t.get_rnd_terminal_node_ref(rng);
        tnode.print();
        println!("\n--------");
    }
}

fn main() {
    init_run();

    #[cfg(gpopt_mult_threaded="yes")]
    println!("CURRENT NUM THREADS= {}", rayon::current_num_threads());

    let mut win_number = 0i32;
    let mut run_number = 0i32;
    let mut rng = GpRngFactory::new();
    loop { // outer loop checks off winners until CONTROL.W 
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
            win_number += 1;
            winner.win = win_number;
            run_number = 0i32;
            let rc = RunContext::new();
            tee_to_run_log(&winner);

            winner.print_result();
            rc.print_run_illustration(
                &format!("Have Winner #{} - Run# {} Gen# {}",
                    win_number, run_number, winner.gen));
            winner.tree.print_exec_one();
            assert!(CONTROL.W >= win_number);
            if CONTROL.W > 0 && CONTROL.W == win_number {
                println!("Win has reached max of {} runs ending.", CONTROL.W);
                break;
            } else {
                EvalCount::reset();
            }
        } else {
            println!("Exceeded CONTROL.R ({}) runs without finding a winner.",
                CONTROL.R);
        }
    }
}
