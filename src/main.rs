#[macro_use]
extern crate lazy_static;
extern crate mut_static;

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
    let control = CONTROL.read().unwrap();
    let log_file = control.run_log_file;

    if Path::new(log_file).exists() {
        if let Err(e) = remove_file(log_file) {
            panic!("failed remove file: {}", e);
        }
    }
    let win: i32 = 5;
    let run: i32 = 10;
    let gen: u16 = 23;
    let e: i64 = 233235235;

    tee_to_run_log(win, run, gen, e);

    // Read
    let mut file = BufReader::new(File::open(log_file).unwrap());

    // read the first line and extract the number from it
    let mut input = String::new();
    file.read_line(&mut input).unwrap();
    input.truncate(input.len() - 1);            // remove new line
    let parsed = scanf!(input, "win={}, run={}, gen={}, effort={}",
        i32, i32, u16, i64);

    assert_eq!(parsed, Some((win, run, gen, e)));
}

fn tee_to_run_log(win: i32, run: i32, gen: u16, e: i64) {
    let msg = format!("win={}, run={}, gen={}, effort={}", win, run, gen, e);

    let control = CONTROL.read().unwrap();
    let log_file = control.run_log_file;
    let create_bool = !Path::new(log_file).exists();
    let mut file =
         OpenOptions::new()
            .write(true)
            .create(create_bool)
            .append(!create_bool)
            .open(log_file)
            .unwrap();

    println!("{}", msg);
    if let Err(_e) = writeln!(file, "{}", msg) {
        panic!("Couldn't write to file");
    }
}

fn report_results(rng: &mut GpRng, trees: &mut TreeSet,header_need: &mut bool,
        hits: &GpHits) -> () {
    let control = CONTROL.read().unwrap();

    if control.show_all_trees {
        println!("Generation {}", trees.gen);
        trees.print();
    }
    if control.show_all_tree_results {
        println!("gen {}", trees.gen);
        Tree::print_result_header(None, hits);
        for t in trees.tree_vec.iter() {
            t.print_result(None, -1.0);
        }
    }
    if control.show_best_tree_results {
        if *header_need {
            Tree::print_result_header(Some(trees.gen), hits);
            *header_need = false;
        }
        let i = trees.tree_vec.len()-1;
        let best_tree = &trees.tree_vec[i];
        best_tree.print_result(Some(trees.gen), trees.avg_raw_f);
    }
    if control.run_tests {
        run_tests(rng, trees);
    }
}

fn run(rng: &mut GpRng, run_number: i32) -> Option<Winner> {
    let control = CONTROL.read().unwrap();
    if control.show_controls {
        println!("M = {}, G = {}, D = {}", control.M, control.G, control.Di);
    }
    let mut trees = TreeSet::create_initial_population(rng);

    trees.count_nodes();
    trees.gen = 0u16;
    let mut header_need: bool = true;
    // when gen == control.G we are done because gen starts with 0
    // ie. we do not process gen == 51 if G==51
    while trees.gen < control.G && trees.winning_index == None {
        let hits = trees.exec_all();
        if trees.winning_index != None {
            break;
        }

        #[cfg(gpopt_trace="on")]
        println!("TP0:trace log for gen={}", trees.gen);

        trees.compute_normalized_fitness();
        trees.sort_by_normalized_fitness();
        report_results(rng, &mut trees, &mut header_need, &hits);

        if trees.gen >= control.G {
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
                e: control.computational_effort(run_number, trees.gen),
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

pub static FUNCS_RPB: &[&[Function]] = &[&FUNCTIONS_RESULT_BRANCH];
pub static TERMS_RPB: &[&[Terminal]] = &[&TERMINALS_RESULT_BRANCH];
pub static FUNCS_FDB: &[&[Function]] = &[&FUNCTIONS_FUNC_DEF_BRANCH];
pub static TERMS_FDB: &[&[Terminal]] = &[&TERMINALS_FUNC_DEF_BRANCH];

fn main() {
    CONTROL.set(Control::new((&FUNCS_RPB, &TERMS_RPB),
                              Some((&FUNCS_FDB, &TERMS_FDB))))
                .expect("error during control/config setup.");
    let control = CONTROL.read().unwrap();

    init_run();

    let mut win_number = 0i32;
    let mut run_number = 0i32;
    let mut rng = GpRngFactory::new();
    loop { // outer loop checks off winners until control.W 
        let opt_winner = loop { // go until we have a winner
            run_number += 1;
            println!("Run #{}", run_number);
            if let Some(winner) = run(&mut rng, run_number) {
                break Some(winner);
            }
            else if run_number == control.R && control.R != 0 {
                break None;
            }
        };
            
        if let Some(mut winner) = opt_winner {
            win_number += 1;
            let rc = RunContext::new();
            tee_to_run_log(win_number, winner.run, winner.gen, winner.e);

            winner.print_result();
            rc.print_run_illustration(
                &format!("Have Winner #{} - Run# {} Gen# {}",
                    win_number, run_number, winner.gen));
            winner.tree.print_exec_one();
            assert!(control.W >= win_number);
            if control.W > 0 && control.W == win_number {
                println!("Win has reached max of {} runs ending.", control.W);
                break;
            }
        } else {
            println!("Exceeded control.R ({}) runs without finding a winner.",
                control.R);
        }
    }
}
