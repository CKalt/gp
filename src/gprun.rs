use crate::tree::Tree;
use crate::tree::TreeSet;
use crate::tree::exec_node;
use crate::tree::GpFloat;

use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::control::CONTROL;

pub enum GpType {
    Value(f32),
    Terminate,
}

// TERMINAL SPECIFICS
pub static TERMINAL: [Terminal; CONTROL.num_terminals as usize] = [
    Terminal {
        tid:  0u8,
        name: "L0",
        code: terminal_l0
    },
    Terminal {
        tid:  1u8,
        name: "W0",
        code: terminal_w0,
    },
    Terminal {
        tid:  2u8,
        name: "H0",
        code: terminal_h0,
    },
    Terminal {
        tid:  3u8,
        name: "L1",
        code: terminal_l1
    },
    Terminal {
        tid:  4u8,
        name: "W1",
        code: terminal_w1,
    },
    Terminal {
        tid:  5u8,
        name: "H1",
        code: terminal_h1,
    },
    Terminal {
        tid:  5u8,
        name: "H1",
        code: terminal_d,
    },
];

fn terminal_l0(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].l0);
}

fn terminal_w0(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].w0);
}

fn terminal_h0(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].h0);
}

fn terminal_l1(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].l1);
}

fn terminal_w1(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].w1);
}

fn terminal_h1(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].h1);
}

fn terminal_d(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].d);
}

pub static FUNCTION: [Function; CONTROL.num_functions as usize] = [
    Function {
        fid:  0u8,
        name: "+",
        arity: 2,
        code: func_add,
    },
    Function {
        fid:  1u8,
        name: "-",
        arity: 2,
        code: func_sub,
    },
    Function {
        fid:  2u8,
        name: "*",
        arity: 2,
        code: func_mult,
    },
    Function {
        fid:  2u8,
        name: "%",
        arity: 2,
        code: func_prot_div,
    },
];

fn eval_binary_op(val1: GpType, val2: GpType, op: fn(i32,i32)->i32 ) -> GpType {
    match val1 {
        Terminate => val1,
        Value(num1) =>
            match val2 {
                Terminate => val2,
                Value(num2) => Value( op(num1, num2) )
            }
    }
}

fn func_add(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = exec_node(rc, &func.branch[0]);
    let val2 = exec_node(rc, &func.branch[1]);
    eval_binary_op(val1, val2, |a,b| a+b)
}

fn func_sub(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = exec_node(rc, &func.branch[0]);
    let val2 = exec_node(rc, &func.branch[1]);
    eval_binary_op(val1, val2, |a,b| a-b)
}

fn func_mult(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = exec_node(rc, &func.branch[0]);
    let val2 = exec_node(rc, &func.branch[1]);
    eval_binary_op(val1, val2, |a,b| a*b)
}

fn func_prot_div(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = exec_node(rc, &func.branch[0]);
    let val2 = exec_node(rc, &func.branch[1]);
    eval_binary_op(val1, val2, |a,b| if b == 0.0 {1.0} else {a/b} )
}


fn prog_n2(rc: &mut RunContext, func: &FunctionNode) -> GpType {
// unconditionally execute branch 0 & 1
// ignore any return values
    exec_node(rc, &func.branch[0]);
    exec_node(rc, &func.branch[1])
}

fn prog_n3(rc: &mut RunContext, func: &FunctionNode) -> GpType {
// unconditionally execute branch 0, 1 & 2
// ignore any return values
    exec_node(rc, &func.branch[0]);
    exec_node(rc, &func.branch[1]);
    exec_node(rc, &func.branch[2])
}

#[allow(dead_code)]
fn clock_reset(rc: &mut RunContext) {
    rc.clock = 0;
}

fn clock_tic(rc: &mut RunContext) -> GpType {
    if rc.clock < RUN_CONTROL.max_clock {
        rc.clock += 1;
        GpType::Continue
    }
    else {
        GpType::Terminate
    }
}

fn terminate(rc: &mut RunContext) -> GpType {
    rc.clock = RUN_CONTROL.max_clock;
    GpType::Terminate
}

/// RunControl defines parameters controlling the running of individuals.
pub struct RunControl {
    pub max_clock: u16, // sets limit for entire population run duration
}

pub static RUN_CONTROL: RunControl = RunControl{
    max_clock: 3000,
};

////////////////////////////////////////////////////////////////////
type FoodCoord = (usize, usize); // (x, y)
#[derive(Copy, Clone)]
pub enum GridCellState {
    Clear,
    Food,
    FoodEaten,
    NoFoodFound,
}
use GridCellState::*;

struct FitnessCase {
    l0: f32,
    w0: f32,
    h0: f32,
    l1: f32,
    w1: f32,
    h1: f32,
    d:  f32,
}

pub const RUN_CONTROL_NUM_FITNESS_CASES: usize = 10;

struct RunContext {
    fitness_cases: [FitnessCase; RUN_CONTROL_NUM_FITNESS_CASES],
    cur_fc_index: usize,
    clock: u16,

}
impl RunContext {
    fn new() -> RunContext {
        let rc = RunContext {
            clock: 0,
            cur_fc_index: 0,
            fitness_cases:
                [
                    FitnessCase{
                        l0:   3.0,
                        w0:   4.0,
                        h0:   7.0,
                        l1:   2.0,
                        w1:   5.0,
                        h1:   3.0,
                        d:   54.0,
                    },
                    FitnessCase{
                        l0:   7.0,
                        w0:  10.0,
                        h0:   9.0,
                        l1:  10.0,
                        w1:   3.0,
                        h1:   1.0,
                        d:  600.0,
                    },
                    FitnessCase{
                        l0: 10.0,
                        w0: 9.0,
                        h0: 4.0,
                        l1: 8.0,
                        w1: 1.0,
                        h1: 6.0,
                        d:  312.0,
                    },
                    FitnessCase{
                        l0: 3.0,
                        w0: 9.0,
                        h0: 5.0,
                        l1: 1.0,
                        w1: 6.0,
                        h1: 4.0,
                        d:  111.0,
                    },
                    FitnessCase{
                        l0: 4.0,
                        w0: 3.0,
                        h0: 2.0,
                        l1: 7.0,
                        w1: 6.0,
                        h1: 1.0,
                        d: -18.0,
                    },
                    FitnessCase{
                        l0:   3.0,
                        w0:   3.0,
                        h0:   1.0,
                        l1:   9.0,
                        w1:   5.0,
                        h1:   4.0,
                        d: -171.0,
                    },
                    FitnessCase{
                        l0:   5.0,
                        w0:   9.0,
                        h0:   9.0,
                        l1:   1.0,
                        w1:   7.0,
                        h1:   6.0,
                        d:  363.0,
                    },
                    FitnessCase{
                        l0:   1.0,
                        w0:   2.0,
                        h0:   9.0,
                        l1:   3.0,
                        w1:   9.0,
                        h1:   2.0,
                        d:  -36.0,
                    },
                    FitnessCase{
                        l0: 2.0,
                        w0: 6.0,
                        h0: 8.0,
                        l1: 2.0,
                        w1: 6.0,
                        h1: 10.0,
                        d: -24.0,
                    },
                    FitnessCase{
                        l0:  8.0,
                        w0:  1.0,
                        h0:  0.0,
                        l1:  7.0,
                        w1:  5.0,
                        h1:  1.0,
                        d:  45.0,
                    },
                ],
            };

        for (i, fc) in rc.iter() {
            if fc.d != (fc.w*fc.h*fc.l) - (fc.w*fc.h*c.l) {
                panic!("Fitness case #{} is invalid.", i);
            }
        }

        rc
    }
    fn prepare_run(&mut self) { 
    }
}

pub fn exec_trees(mut trees: &mut TreeSet, run_number: i32) -> u16 {
    let mut rc = RunContext::new();

    for (i, tree) in trees.tree_vec.iter_mut().enumerate() {
        rc.prepare_run();
        while rc.clock < RUN_CONTROL.max_clock {
            exec_node(&mut rc, &mut tree.root);
        }

        if tree.compute_fitness(&rc) {
            report_tree_result(tree, tree.tfid, None, -1.0);
            rc.print_grid(&format!("Have Winner! - Run# {} Gen# {}", run_number,
                trees.gen));
            trees.winning_index = Some(i);
            break;
        }
    }
    rc.n_pellets
}

pub fn exec_single_tree(tree : &mut Tree) {
    let mut rc = RunContext::new();
    rc.prepare_run();
    rc.print_grid("Before Run");
    while rc.clock < RUN_CONTROL.max_clock {
        exec_node(&mut rc, &mut tree.root);
    }

    if tree.compute_fitness(&rc) {
        println!("Have Winner");
    }
    rc.print_grid("After Run");
}

pub fn report_tree_result(t: &Tree, i: Option<usize> , opt_gen: Option<u16>, avg_raw_f: GpFloat) {
    assert_eq!(i, t.tfid);
    let f = &t.fitness;
    let tfid = if let Some(num) = t.tfid { num } else { 0 };
    if let Some(gen) = opt_gen {
        println!("{:6} {:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                 gen, tfid, t.tcid, t.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
    } else {
        println!("{:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                tfid, t.tcid, t.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
    }
}

pub fn tree_result_header(opt_gen: Option<u16>, n_pellets: &u16) {
    println!("n_pellets={}", *n_pellets);
    if let Some(_) = opt_gen {
        println!("{:>6} {:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8} {:>6}", 
            "gen", "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr", "avgRawF");
        println!("{:>6} {:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
            "----", "---", "---", "-----", "---", "---", "------", "------", "------");
    } else {
        println!("{:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
            "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr");
        println!("{:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
            "---", "---", "-----", "---", "---", "------", "------", "------");
    }
}

