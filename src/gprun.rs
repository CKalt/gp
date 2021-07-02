use crate::tree::Tree;
use crate::tree::exec_node;
use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::control::CONTROL;

use crate::tree::GpFloat;
use crate::tree::GpHits;
use crate::tree::GpRaw;
    
#[cfg(gpopt_fitness_type="int")]
use crate::tree::GpInt;
#[cfg(gpopt_fitness_type="int")]
use crate::tree::DL_SHIFT;

#[cfg(gpopt_fitness_type="int")]
use crate::Fitness;

pub enum GpType {
    Init,
    Value(GpRaw),
    Terminate,
}
use GpType::*;

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
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].l0);
}

fn terminal_w0(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].w0);
}

fn terminal_h0(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].h0);
}

fn terminal_l1(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].l1);
}

fn terminal_w1(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].w1);
}

fn terminal_h1(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    return GpType::Value(rc.fitness_cases[rc.cur_fc_index].h1);
}

fn terminal_d(rc: &mut RunContext) -> GpType {
    #[cfg(gpopt_clock_termination="on")]
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
        fid:  3u8,
        name: "%",
        arity: 2,
        code: func_prot_div,
    },
];

fn eval_binary_op(val1: GpType, val2: GpType, op: fn(GpRaw,GpRaw)->GpRaw) -> GpType {
    match val1 {
        Terminate => val1,
        Value(num1) =>
            match val2 {
                Terminate => val2,
                Value(num2) => Value( op(num1, num2) ),
                Init => panic!("attempt to eval an invalid Init GpType"),
            },
        Init => panic!("attempt to eval an invalid Init GpType"),
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

#[cfg(gpopt_clock_termination="on")]
fn clock_reset(rc: &mut RunContext) {
    rc.clock = 0;
}

#[cfg(gpopt_clock_termination="on")]
fn clock_tic(rc: &mut RunContext) -> GpType {
    if rc.clock < RUN_CONTROL.max_clock {
        rc.clock += 1;
        GpType::Continue
    }
    else {
        GpType::Terminate
    }
}

#[cfg(gpopt_clock_termination="on")]
fn terminate(rc: &mut RunContext) -> GpType {
    rc.clock = RUN_CONTROL.max_clock;
    GpType::Terminate
}

/// RunControl defines parameters controlling the running of individuals.
pub struct RunControl {
    #[cfg(gpopt_clock_termination="on")]
    pub max_clock: u16, // sets limit for entire population run duration
}

pub static RUN_CONTROL: RunControl = RunControl {
#[cfg(gpopt_clock_termination="on")]
    max_clock: 3000,
};

pub const RUN_CONTROL_NUM_FITNESS_CASES: usize = 10;

// Fitness
struct FitnessCase {
    l0: GpRaw,
    w0: GpRaw,
    h0: GpRaw,
    l1: GpRaw,
    w1: GpRaw,
    h1: GpRaw,
    d:  GpRaw,
}
impl FitnessCase {
    pub fn compute_error(&self, result: GpType) -> GpRaw {
        match result {
            Value(result_value) => {
                (result_value - self.d).abs()
            },
            Init =>
                panic!("Unable to eval fitness case for Non-Value result=Init"),
            Terminate =>
                panic!("Unable to eval fitness case for Non-Value result=Terminate"),

        }
    }
}

/// RunContext provides runtime control over a running individual. Each 
/// node and terminal exec call recieves a reference to its RunContext
/// where it can then access it's fitness case data and currency values.
pub struct RunContext {
    pub fitness_cases: [FitnessCase; RUN_CONTROL_NUM_FITNESS_CASES],
    #[cfg(gpopt_clock_termination="on")]
    clock: u16,
    pub cur_fc_index: usize,

    pub hits: GpHits,
    pub error: GpRaw,
}
impl RunContext {
    pub fn new() -> RunContext {
        let rc = RunContext {
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
            #[cfg(gpopt_clock_termination="on")]
            clock: 0,
            hits: 0,
            error: 0.0,
        };

        for (i, fc) in rc.fitness_cases.iter().enumerate() {
            if fc.d != (fc.w0 * fc.h0 * fc.l0) - (fc.w1 * fc.h1 * fc.l1) {
                panic!("Fitness case #{} is invalid.", i);
            }
        }
        rc
    }
    pub fn print_run_illustration(&self, _label: &str) { }
    pub fn prepare_run(&mut self) { }
    pub fn get_hits_label() -> &'static str {
        "num cases w/error lt 0.01"
    }

    /// computes fitness returns true if winner
    #[cfg(gpopt_fitness_type="float")]
    pub fn compute_fitness(&self, tree: &mut Tree) -> bool {
        let mut f = &mut tree.fitness;
        f.r = self.error;
        tree.hits = self.hits;

        // init fitness "base" values
        f.n = -1.0;
        f.a = -1.0;
        f.nfr = -1.0;
        f.raw = f.r;

        // average over generation
        f.a = 1.0 / (1.0 + f.s as GpFloat);

        // if each fitness case was a hit then we have a winner.
        return self.hits == (self.fitness_cases.len() as GpHits);
    }
    #[cfg(gpopt_fitness_type="int")]
    pub fn compute_fitness(&self, tree: &mut Tree) -> bool {
        panic!("int fitness_type is not implemented.");
    }
}

pub fn init_run() { }
