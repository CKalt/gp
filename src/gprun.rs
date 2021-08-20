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

pub const EVEN_PARITY_K_VALUE: usize = 5;

#[cfg(gpopt_num_adf="1")]
pub const NUM_ADF: usize = 1;
#[cfg(gpopt_num_adf="2")]
pub const NUM_ADF: usize = 2;
#[cfg(gpopt_num_adf="3")]
pub const NUM_ADF: usize = 3;
#[cfg(gpopt_num_adf="4")]
pub const NUM_ADF: usize = 4;
#[cfg(gpopt_num_adf="5")]
pub const NUM_ADF: usize = 5;

#[cfg(gpopt_adf_parity="2")]
pub const ADF_PARITY: usize = 2;
#[cfg(gpopt_adf_parity="3")]
pub const ADF_PARITY: usize = 3;
#[cfg(gpopt_adf_parity="4")]
pub const ADF_PARITY: usize = 4;

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
    let args: Vec<GpTyp> = Vec::new();
    for b_i in 0..ADF_PARITY {
        args.push(Tree::exec_node(rc, &func.branch[b_i]));
    }

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, &args);
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
pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    vec![
        vec![
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
        ]
    ]
}

#[cfg(gpopt_adf="yes")]
pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    let mut funcs = vec![
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

    for a_i in 0..NUM_ADF {
        funcs.push(
            Function {
                fid:  a_i as u8 + 4u8,
                name: &format!("ADF{}", a_i),
                arity: ADF_ARITY as u8,
                code: function_adf,
                opt_adf_num: Some(a_i),   // ideintifies ADFn
            }
        );
    }
    vec![
        // FUNCTIONS FOR RESULT BRANCH 0 (only one result branch)
        funcs
    ]
}

pub fn get_terminals_for_result_branches() -> Vec<Vec<Terminal>> {
    vec![ 
        // TERMINALS FOR RESULT BRANCH 0 (only one result branch)
        vec![
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
        ]
    ]
}

#[cfg(gpopt_adf="yes")]
pub fn get_functions_for_func_def_branches() -> Vec<Vec<Function>> {
    let funcs = vec![
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
}

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
    pub opt_adf_args: Option<&'a Vec<GpType>>,
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
    pub fn exec_adf(&mut self, adf_num: usize, args: &Vec<GpType>)
            -> GpType {
        let func_def_branch = 
            self
                .opt_func_def_branches
                .as_ref()
                .expect("exec_adf with None set for branches.")[adf_num];

        match self.opt_adf_args {
            None    => {
                self.opt_adf_args = Some(&args);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref orig_args) => {
                let save_args = orig_args.clone();
                self.opt_adf_args = Some(args.clone());
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(save_args.clone());
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
