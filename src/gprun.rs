use crate::fitness::Fitness;
use crate::tree::Tree;
use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::tree::TreeBranch;

use crate::fitness::GpFloat;
use crate::fitness::GpFitness;
use crate::fitness::GpHits;
use crate::fitness::GpRaw;

pub type GpType = bool;

//fn fix_nan(num: GpType) -> GpType {
//    if num.is_nan() {
//        0.0
//    } else {
//        num
//    }
//}

fn function_and(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(rc, &func.branch[0]);
    let val2 = Tree::exec_node(rc, &func.branch[1]);
    val1 & val2
}

fn function_or(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(rc, &func.branch[0]);
    let val2 = Tree::exec_node(rc, &func.branch[1]);
    val1 | val2
}

fn function_nand(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(rc, &func.branch[0]);
    let val2 = Tree::exec_node(rc, &func.branch[1]);
    !(val1 & val2)
}

fn function_nor(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let val1 = Tree::exec_node(rc, &func.branch[0]);
    let val2 = Tree::exec_node(rc, &func.branch[1]);
    !(val1 | val2)
}

#[cfg(gpopt_adf="yes")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, arg1, arg2, arg3, arg4)
}

fn terminal_d0(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[0]
}

fn terminal_d1(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[1]
}

fn terminal_d2(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[2]
}

fn terminal_d3(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[3]
}

fn terminal_d4(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[4]
}

fn terminal_d5(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[5]
}

fn terminal_d6(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[6]
}

#[cfg(any(gpopt_even_parity_k="8",
          gpopt_even_parity_k="9",
          gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
fn terminal_d7(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[7]
}

#[cfg(any(gpopt_even_parity_k="9",
          gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
fn terminal_d8(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[8]
}

#[cfg(any(gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
fn terminal_d9(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[9]
}

#[cfg(gpopt_even_parity_k="11")]
fn terminal_d10(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[10]
}

/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg0(rc: &RunContext) -> GpType {
     match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[0],
    }
}

/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg1(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[1],
    }
}

/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg2(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[2],
    }
}
            
/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg3(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[3],
    }
}

#[cfg(gpopt_adf="no")]
pub static FUNCTIONS_RESULT_BRANCH_NO_ADF: [Function; 4] = [
    Function {
        fid:  1u8,
        name: "AND",
        arity: 2,
        code: function_and,
        opt_adf_num: None,
    },
    Function {
        fid:  2u8,
        name: "OR",
        arity: 2,
        code: function_or,
        opt_adf_num: None,
    },
    Function {
        fid:  3u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
        opt_adf_num: None,
    },
    Function {
        fid:  4u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
        opt_adf_num: None,
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_RESULT_BRANCH_ADF: [Function; 6] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
        opt_adf_num: None,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
        opt_adf_num: None,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
        opt_adf_num: None,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
        opt_adf_num: None,
    },
    Function {
        fid:  4u8,
        name: "ADF0",
        arity: EVEN_PARITY_K_VALUE as u8 - 1u8,
        code: function_adf,
        opt_adf_num: Some(0),   // ideintifies ADF0
    },
    Function {
        fid:  5u8,
        name: "ADF1",
        arity: EVEN_PARITY_K_VALUE as u8 - 1u8,
        code: function_adf,
        opt_adf_num: Some(1),   // ideintifies ADF0
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_FUNC_DEF_BRANCH_ADF0: [Function; 4] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
        opt_adf_num: None,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
        opt_adf_num: None,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
        opt_adf_num: None,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
        opt_adf_num: None,
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_FUNC_DEF_BRANCH_ADF1: [Function; 5] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
        opt_adf_num: None,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
        opt_adf_num: None,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
        opt_adf_num: None,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
        opt_adf_num: None,
    },
    Function {
        fid:  4u8,
        name: "ADF0",
        arity: EVEN_PARITY_K_VALUE as u8 - 1u8,
        code: function_adf,
        opt_adf_num: Some(0),   // ideintifies ADF0
    },
];

// TERMINAL SPECIFICS - RESULT PRODUCING BRANCH - result_branch

#[cfg(gpopt_even_parity_k="7")]
pub const EVEN_PARITY_K_VALUE: usize = 7;
#[cfg(gpopt_even_parity_k="8")]
pub const EVEN_PARITY_K_VALUE: usize = 8;
#[cfg(gpopt_even_parity_k="9")]
pub const EVEN_PARITY_K_VALUE: usize = 9;
#[cfg(gpopt_even_parity_k="10")]
pub const EVEN_PARITY_K_VALUE: usize = 10;
#[cfg(gpopt_even_parity_k="11")]
pub const EVEN_PARITY_K_VALUE: usize = 11;

pub static TERMINALS_RESULT_BRANCH: [Terminal; EVEN_PARITY_K_VALUE] = [
    Terminal {
        tid:  0u8,
        name: "D0",
        code: terminal_d0,
    },
    Terminal {
        tid:  1u8,
        name: "D1",
        code: terminal_d1,
    },
    Terminal {
        tid:  2u8,
        name: "D2",
        code: terminal_d2,
    },
    Terminal {
        tid:  3u8,
        name: "D3",
        code: terminal_d3,
    },
    Terminal {
        tid:  4u8,
        name: "D4",
        code: terminal_d4,
    },
    Terminal {
        tid:  5u8,
        name: "D5",
        code: terminal_d5,
    },
    Terminal {
        tid:  6u8,
        name: "D6",
        code: terminal_d6,
    },
#[cfg(any(gpopt_even_parity_k="8",
          gpopt_even_parity_k="9",
          gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
    Terminal {
        tid:  7u8,
        name: "D7",
        code: terminal_d7,
    },
#[cfg(any(gpopt_even_parity_k="9",
          gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
    Terminal {
        tid:  8u8,
        name: "D8",
        code: terminal_d8,
    },
#[cfg(any(gpopt_even_parity_k="10",
          gpopt_even_parity_k="11"))]
    Terminal {
        tid:  9u8,
        name: "D9",
        code: terminal_d9,
    },
#[cfg(gpopt_even_parity_k="11")]
    Terminal {
        tid:  10u8,
        name: "D10",
        code: terminal_d10,
    },
];

// TERMINAL SPECIFICS FUNCTION DEFINING BRANCH - func_def_branch
// Same for ADF0 and ADF1
#[cfg(gpopt_adf="yes")]
pub static TERMINALS_FUNC_DEF_BRANCH_ADF_0_1: [Terminal; 4] = [
    Terminal {
        tid:  0u8,
        name: "ARG0",
        code: terminal_adf_arg0,
    },
    Terminal {
        tid:  1u8,
        name: "ARG1",
        code: terminal_adf_arg1,
    },
    Terminal {
        tid:  2u8,
        name: "ARG2",
        code: terminal_adf_arg2,
    },
    Terminal {
        tid:  3u8,
        name: "ARG3",
        code: terminal_adf_arg3,
    },
];

pub const RUN_CONTROL_NUM_FITNESS_CASES: usize =
    2usize.pow(EVEN_PARITY_K_VALUE as u32);

// FitnessCase
pub struct FitnessCase {
    input_bits: [bool; EVEN_PARITY_K_VALUE],
    output_bit: bool,
}
impl FitnessCase {
    fn new() -> FitnessCase {
        FitnessCase {
            input_bits: [false; EVEN_PARITY_K_VALUE],
            output_bit: false,
        }
    }
}

/// RunContext provides runtime control over a running individual. Each 
/// node and terminal exec call recieves a reference to its RunContext
/// where it can then access it's fitness case data and currency values.
pub struct RunContext<'a> {
    pub fitness_cases: Vec::<FitnessCase>,
    pub opt_func_def_branches: Option<Vec<&'a TreeBranch>>, // adf0, adf1,...
    pub opt_adf_args: Option<Vec<GpType>>,
    pub cur_fc: usize,
    pub hits: GpHits,
    pub error: GpRaw,
}
impl<'a> RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let mut rc = RunContext {
            cur_fc: 0,
            opt_func_def_branches: None,
            opt_adf_args: None,
            fitness_cases: Vec::new(),
            hits: 0,
            error: 0
        };

        for i in 0..RUN_CONTROL_NUM_FITNESS_CASES {
            let mut fc = FitnessCase::new();
            let mut sum_of_bits = 0u16;
            for bit in 0..EVEN_PARITY_K_VALUE {
                let bit_value: bool = (i & 2u8.pow(bit as u32) as usize) != 0;
                fc.input_bits[bit] = bit_value;
                if bit_value {
                    sum_of_bits += 1;
                }
            }
            fc.output_bit = ( sum_of_bits % 2 ) == 0;
            rc.fitness_cases.push(fc);
        }

        rc
    }
    pub fn get_cur_fc(&self) -> &FitnessCase {
        &self.fitness_cases[self.cur_fc]
    }
    pub fn print_run_illustration(&self, _label: &str) { }
    pub fn prepare_run(&mut self) { }
    pub fn get_hits_label() -> &'static str {
        "num cases w/error lt 0.01"
    }
    /// computes fitness returns true if winner
    pub fn compute_fitness(&self) -> (Fitness, bool) {
        let mut f = Fitness::new();

        f.r = self.error;
        f.hits = self.hits;

        f.n = -1.0;
        f.nfr = -1.0;
        f.raw = self.hits as GpFitness;
        f.s = f.r;
        f.a = 1.0 / (1.0 + f.s as GpFloat);

        // if each fitness case was a hit then we have a winner.
        let max_possible_hits = self.fitness_cases.len() as GpHits;
        let is_winner = self.hits == max_possible_hits;
        (f, is_winner)
    }

    /// exec_adf: executes and automatically defined function specified
    /// by the func_def_branch arg.
    #[cfg(gpopt_adf="yes")]
    pub fn exec_adf(&mut self, adf_num: usize, arg1: GpType,
                arg2: GpType, arg3: GpType, arg4: GpType)
            -> GpType {
        let func_def_branch = 
            self
                .opt_func_def_branches
                .as_ref()
                .expect("exec_adf with None set for branches.")[adf_num];

        match self.opt_adf_args {
            None    => {
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3, arg4 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref args) => {
                let (a,b,c,d) = (args[0], args[1], args[2], args[3]);
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3, arg4 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(vec![ a, b, c, d ]);
                result
            }
        }
    }

    // compute error for a single fitness case
    pub fn compute_error(&self, result: GpType) -> GpRaw {
        // since GpType is bool
        // only two possible values 0 or 1
        // 0 means output matched the target, 1 means it missed.
        let correct_result = self.fitness_cases[self.cur_fc].output_bit;
        if correct_result == result {
            0
        } else {
            1
        }
    }
}

pub fn init_run() { }

#[cfg(test)]
pub mod tests {
    #[test]
    #[cfg(gpopt_adf="yes")]
    pub fn test_print_exec_one() {
        use crate::tree::*;
        // first build tree:
        let (rb_str, fd_vec) = get_k_n_test();
        let mut tree = Tree::parse(rb_str, fd_vec);

        assert_eq!(tree.print_exec_one(), true);
    }

    #[cfg(gpopt_even_parity_k="7")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
"#;
        let adf0_str = r#"
"#;
        let adf1_str = r#"
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="8")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
"#;
        let adf0_str = r#"
"#;
        let adf1_str = r#"
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="9")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
"#;
        let adf0_str = r#"
"#;
        let adf1_str = r#"
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="10")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
"#;
        let adf0_str = r#"
"#;
        let adf1_str = r#"
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="11")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
"#;
        let adf0_str = r#"
"#;
        let adf1_str = r#"
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }
}
