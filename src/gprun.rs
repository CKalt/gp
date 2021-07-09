use crate::tree::Tree;
use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::tree::TreeBranch;
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

#[cfg(gpopt_exec_criteria="each_fitness_case")]
pub type GpType = GpRaw;

fn function_add(fc: &FitnessCase, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(fc, &func.branch[0]);
    let val2 = Tree::exec_node(fc, &func.branch[1]);
    val1 + val2
}

fn function_sub(fc: &FitnessCase, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(fc, &func.branch[0]);
    let val2 = Tree::exec_node(fc, &func.branch[1]);
    val1 - val2
}

fn function_mult(fc: &FitnessCase, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(fc, &func.branch[0]);
    let val2 = Tree::exec_node(fc, &func.branch[1]);
    val1 * val2
}

fn function_prot_div(fc: &FitnessCase, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(fc, &func.branch[0]);
    let val2 = Tree::exec_node(fc, &func.branch[1]);
    if val2 == 0.0 {1.0} else {val1/val2}
}

fn function_adf0(fc: &mut FitnessCase, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(fc, &func.branch[0]);
    let arg2 = Tree::exec_node(fc, &func.branch[1]);
    let arg3 = Tree::exec_node(fc, &func.branch[2]);

    fc.exec_adf0(arg1, arg2, arg3)
}

fn terminal_l0(fc: &FitnessCase) -> GpType {
    fc.l0
}

fn terminal_w0(fc: &FitnessCase) -> GpType {
    fc.w0
}

fn terminal_h0(fc: &FitnessCase) -> GpType {
    fc.h0
}

fn terminal_l1(fc: &FitnessCase) -> GpType {
    fc.l1
}

fn terminal_w1(fc: &FitnessCase) -> GpType {
    fc.w1
}

fn terminal_h1(fc: &FitnessCase) -> GpType {
    fc.h1
}

fn terminal_arg0(fc: &FitnessCase) -> GpType {
    match fc.adf0_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[0],
    }
}

fn terminal_arg1(fc: &FitnessCase) -> GpType {
    match fc.adf0_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[1],
    }
}

fn terminal_arg2(fc: &FitnessCase) -> GpType {
    match fc.adf0_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[2],
    }
}

pub static FUNCTIONS_RESULT_BRANCH: [Function; CONTROL.num_functions_result_branch as usize] = [
    Function {
        fid:  0u8,
        name: "+",
        arity: 2,
        code: function_add,
    },
    Function {
        fid:  1u8,
        name: "-",
        arity: 2,
        code: function_sub,
    },
    Function {
        fid:  2u8,
        name: "*",
        arity: 2,
        code: function_mult,
    },
    Function {
        fid:  3u8,
        name: "%",
        arity: 2,
        code: function_prot_div,
    },
    Function {
        fid:  4u8,
        name: "ADF0",
        arity: 3,
        code: function_adf0,
    },
];

pub static FUNCTIONS_FUNC_DEF_BRANCH: [Function; CONTROL.num_functions_func_def_branch as usize] = [
    Function {
        fid:  0u8,
        name: "+",
        arity: 2,
        code: function_add,
    },
    Function {
        fid:  1u8,
        name: "-",
        arity: 2,
        code: function_sub,
    },
    Function {
        fid:  2u8,
        name: "*",
        arity: 2,
        code: function_mult,
    },
    Function {
        fid:  3u8,
        name: "%",
        arity: 2,
        code: function_prot_div,
    },
];

// TERMINAL SPECIFICS - RESULT PRODUCING BRANCH - result_branch
pub static TERMINALS_RESULT_BRANCH: [Terminal; CONTROL.num_terminals_result_branch as usize] = [
    Terminal {
        tid:  0u8,
        name: "L0",
        code: terminal_l0,
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
];

// TERMINAL SPECIFICS FUNCTION DEFINING BRANCH - func_def_branch
pub static TERMINALS_FUNC_DEF_BRANCH: [Terminal; CONTROL.num_terminals_func_def_branch as usize] = [
    Terminal {
        tid:  0u8,
        name: "ARG0",
        code: terminal_arg0,
    },
    Terminal {
        tid:  1u8,
        name: "ARG1",
        code: terminal_arg1,
    },
    Terminal {
        tid:  2u8,
        name: "ARG2",
        code: terminal_arg2,
    },
];

pub const RUN_CONTROL_NUM_FITNESS_CASES: usize = 10;

// Fitness
pub struct FitnessCase<'a> {
    l0: GpRaw,
    w0: GpRaw,
    h0: GpRaw,
    l1: GpRaw,
    w1: GpRaw,
    h1: GpRaw,
    pub d:  GpRaw,
    pub func_def_branch: Option<&'a TreeBranch>, // used for calling adf0
    adf0_args: Option<Vec<GpType>>,
}
impl FitnessCase<'_> {
    #[cfg(gpopt_exec_criteria="each_fitness_case")]
    pub fn compute_error(&self, result: GpType) -> GpRaw {
        (result- self.d).abs()
    }
    pub fn exec_adf0(&mut self, arg1: GpType, arg2: GpType, arg3: GpType) -> GpType {
        let func_def_branch = self.func_def_branch.expect("branch not assigned for exec_adf0");

        match self.adf0_args {
            None    => self.adf0_args = Some(vec![ arg1, arg2, arg3 ]),
            Some(_) => panic!("exec_adf0 with non empty as args list"),
        }

        let result = Tree::exec_node(self, &func_def_branch.root);
        self.adf0_args = None;
        result
    }
}

/// RunContext provides runtime control over a running individual. Each 
/// node and terminal exec call recieves a reference to its RunContext
/// where it can then access it's fitness case data and currency values.
pub struct RunContext<'a> {
    pub fitness_cases: [FitnessCase<'a>; RUN_CONTROL_NUM_FITNESS_CASES],
    pub hits: GpHits,
    pub error: GpRaw,
}
impl RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let rc = RunContext {
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
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0:   7.0,
                        w0:  10.0,
                        h0:   9.0,
                        l1:  10.0,
                        w1:   3.0,
                        h1:   1.0,
                        d:  600.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0: 10.0,
                        w0: 9.0,
                        h0: 4.0,
                        l1: 8.0,
                        w1: 1.0,
                        h1: 6.0,
                        d:  312.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0: 3.0,
                        w0: 9.0,
                        h0: 5.0,
                        l1: 1.0,
                        w1: 6.0,
                        h1: 4.0,
                        func_def_branch: None,
                        adf0_args: None,
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
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0:   3.0,
                        w0:   3.0,
                        h0:   1.0,
                        l1:   9.0,
                        w1:   5.0,
                        h1:   4.0,
                        d: -171.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0:   5.0,
                        w0:   9.0,
                        h0:   9.0,
                        l1:   1.0,
                        w1:   7.0,
                        h1:   6.0,
                        d:  363.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0:   1.0,
                        w0:   2.0,
                        h0:   9.0,
                        l1:   3.0,
                        w1:   9.0,
                        h1:   2.0,
                        d:  -36.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0: 2.0,
                        w0: 6.0,
                        h0: 8.0,
                        l1: 2.0,
                        w1: 6.0,
                        h1: 10.0,
                        d: -24.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                    FitnessCase{
                        l0:  8.0,
                        w0:  1.0,
                        h0: 10.0,
                        l1:  7.0,
                        w1:  5.0,
                        h1:  1.0,
                        d:  45.0,
                        func_def_branch: None,
                        adf0_args: None,
                    },
                ],
            hits: 0,
            error: 0.0,
        };

        for (i, fc) in rc.fitness_cases.iter().enumerate() {
            if fc.d != (fc.w0 * fc.h0 * fc.l0) - (fc.w1 * fc.h1 * fc.l1) {
                println!("({}*{}*{}) - ({}*{}*{}) = {}, d={}",
                    fc.w0, fc.h0, fc.l0,
                    fc.w1, fc.h1, fc.l1,
                    (fc.w0 * fc.h0 * fc.l0) - (fc.w1 * fc.h1 * fc.l1),
                    fc.d);
                panic!("Fitness case #{} is invalid.", i+1);
            }
        }
        rc
    }
    pub fn print_run_illustration(&self, _label: &str) {
    }
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
        f.s = f.r;

        // average over generation
        f.a = 1.0 / (1.0 + f.s as GpFloat);

        // if each fitness case was a hit then we have a winner.
        let max_possible_hits = self.fitness_cases.len() as GpHits;
        let result = self.hits == max_possible_hits;
        return result;
    }
    #[cfg(gpopt_fitness_type="int")]
    pub fn compute_fitness(&self, tree: &mut Tree) -> bool {
        panic!("int fitness_type is not implemented.");
    }
}

pub fn init_run() { }
