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
#[cfg(gpopt_even_parity_k="3")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, arg1, arg2)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="4")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, arg1, arg2, arg3)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="5")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, arg1, arg2, arg3, arg4)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="6")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);
    let arg5 = Tree::exec_node(rc, &func.branch[4]);

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, arg1, arg2, arg3, arg4, arg5)
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

#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_d3(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[3]
}

#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_d4(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[4]
}

#[cfg(gpopt_even_parity_k="6")]
fn terminal_d5(rc: &RunContext) -> GpType {
    rc.get_cur_fc().input_bits[5]
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
#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_adf_arg2(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[2],
    }
}
            
/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_adf_arg3(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[3],
    }
}

/// The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="6")]
fn terminal_adf_arg4(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[4],
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
        opt_adf_num: Some(1),   // ideintifies ADF1
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

#[cfg(gpopt_even_parity_k="3")]
pub const EVEN_PARITY_K_VALUE: usize = 3;
#[cfg(gpopt_even_parity_k="4")]
pub const EVEN_PARITY_K_VALUE: usize = 4;
#[cfg(gpopt_even_parity_k="5")]
pub const EVEN_PARITY_K_VALUE: usize = 5;
#[cfg(gpopt_even_parity_k="6")]
pub const EVEN_PARITY_K_VALUE: usize = 6;

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
#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  3u8,
        name: "D3",
        code: terminal_d3,
    },
#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  4u8,
        name: "D4",
        code: terminal_d4,
    },
#[cfg(gpopt_even_parity_k="6")]
    Terminal {
        tid:  5u8,
        name: "D5",
        code: terminal_d5,
    },
];

// TERMINAL SPECIFICS FUNCTION DEFINING BRANCH - func_def_branch
// Same for ADF0 and ADF1
#[cfg(gpopt_adf="yes")]
pub static TERMINALS_FUNC_DEF_BRANCH_ADF_0_1: [Terminal; EVEN_PARITY_K_VALUE-1] = [
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
#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  2u8,
        name: "ARG2",
        code: terminal_adf_arg2,
    },
#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  3u8,
        name: "ARG3",
        code: terminal_adf_arg3,
    },
#[cfg(gpopt_even_parity_k="6")]
    Terminal {
        tid:  4u8,
        name: "ARG4",
        code: terminal_adf_arg4,
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
    #[cfg(gpopt_even_parity_k="3")]
    pub fn exec_adf(&mut self, adf_num: usize, arg1: GpType,
                arg2: GpType)
            -> GpType {
        let func_def_branch = 
            self
                .opt_func_def_branches
                .as_ref()
                .expect("exec_adf with None set for branches.")[adf_num];

        match self.opt_adf_args {
            None    => {
                self.opt_adf_args = Some(vec![ arg1, arg2 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref args) => {
                let (a,b) = (args[0], args[1]);
                self.opt_adf_args = Some(vec![ arg1, arg2 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(vec![ a, b ]);
                result
            }
        }
    }
    /// exec_adf: executes and automatically defined function specified
    /// by the func_def_branch arg.
    #[cfg(gpopt_adf="yes")]
    #[cfg(gpopt_even_parity_k="4")]
    pub fn exec_adf(&mut self, adf_num: usize, arg1: GpType,
            arg2: GpType, arg3: GpType)
            -> GpType {
        let func_def_branch = 
            self
                .opt_func_def_branches
                .as_ref()
                .expect("exec_adf with None set for branches.")[adf_num];

        match self.opt_adf_args {
            None    => {
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref args) => {
                let (a,b,c) = (args[0], args[1], args[2]);
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(vec![ a, b, c ]);
                result
            }
        }
    }
    /// exec_adf: executes and automatically defined function specified
    /// by the func_def_branch arg.
    #[cfg(gpopt_adf="yes")]
    #[cfg(gpopt_even_parity_k="5")]
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
    /// exec_adf: executes and automatically defined function specified
    /// by the func_def_branch arg.
    #[cfg(gpopt_adf="yes")]
    #[cfg(gpopt_even_parity_k="6")]
    pub fn exec_adf(&mut self, adf_num: usize, arg1: GpType,
                arg2: GpType, arg3: GpType, arg4: GpType, arg5: GpType)
            -> GpType {
        let func_def_branch = 
            self
                .opt_func_def_branches
                .as_ref()
                .expect("exec_adf with None set for branches.")[adf_num];

        match self.opt_adf_args {
            None    => {
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3, arg4, arg5 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref args) => {
                let (a,b,c,d,e) = (args[0], args[1], args[2], args[3], args[4]);
                self.opt_adf_args = Some(vec![ arg1, arg2, arg3, arg4, arg5 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(vec![ a, b, c, d, e ]);
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

    #[cfg(gpopt_even_parity_k="3")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
(NAND
  (NAND
    (ADF1
      (NOR D2 D2)
      (ADF1 D0 D1))
    (ADF0
      (ADF0 D2 D0)
      (AND D1 D1)))
  (ADF0
    (OR
      (ADF0 D2 D1)
      (ADF1 D2 D0))
    (NAND
      (NAND D1 D0)
      (ADF0 D2 D2)))
)
"#;
        let adf0_str = r#"
(AND
  (NAND
    (NAND
      (OR ARG1 ARG0)
      (OR ARG0 ARG0))
    (NAND
      (NAND ARG0 ARG1)
      (NOR ARG1 ARG1)))
  (NAND
    (AND
      (OR ARG0 ARG1)
      (AND ARG1 ARG0))
    (AND
      (NAND ARG0 ARG0)
      (AND ARG1 ARG1)))
)
"#;
        let adf1_str = r#"
(ADF0
  (OR
    (AND
      (NOR ARG0 ARG1)
      (OR ARG0 ARG0))
    (NAND
      (OR ARG1 ARG0)
      (NAND ARG1 ARG0)))
  (AND
    (ADF0
      (ADF0 ARG1 ARG0)
      (AND ARG1 ARG1))
    (OR
      (NOR ARG0 ARG1)
      (ADF0 ARG0 ARG0)))
)
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="4")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
(ADF0
  (ADF0
    (ADF1
      (AND
        (OR D0 D2)
        (ADF1 D2 D1 D1))
      (ADF0
        (ADF1 D1 D3 D3)
        (OR D3 D2)
        (NOR D2 D3))
      (AND D3 D1))
    (OR
      (ADF1
        (NOR D0 D2)
        (ADF1 D0 D1 D1)
        (AND D0 D2))
      (NAND
        (ADF0 D1 D2 D3)
        (AND D2 D1)))
    (NOR
      (ADF0
        (NAND D3 D3)
        (OR D3 D1)
        (ADF1 D0 D1 D3))
      (NAND
        (ADF1 D3 D0 D0)
        (OR D0 D2))))
  (NAND
    (ADF1
      (NOR
        (ADF1 D0 D0 D1)
        (NAND D0 D0))
      (ADF0
        (ADF0 D0 D1 D3)
        (NOR D1 D0)
        (AND D0 D2))
      (OR
        (ADF1 D3 D3 D1)
        (NOR D2 D1)))
    (NAND
      (OR
        (ADF1
          (NOR D3 D2) D0 D2)
        (AND D2 D1))
      (ADF1
        (AND D1 D0)
        (ADF1 D0 D2 D1)
        (NAND D3 D0))))
  (NAND
    (NAND
      (ADF1
        (OR D3 D0)
        (AND D3 D2)
        (ADF0
          (ADF1 D1 D3 D3)
          (OR D3 D2)
          (NOR D2 D3)))
      (ADF0
        (AND D3 D2)
        (OR D3 D1)
        (ADF1 D3 D0 D0)))
    (ADF0
      (NAND
        (AND D1 D1)
        (OR D2 D1))
      (NAND
        (NAND D3 D0)
        (AND D0 D3))
      (ADF0
        (NAND D2 D0)
        (OR D0 D2)
        (NOR D1 D1))))
)
"#;
        let adf0_str = r#"
(NOR
  (NOR
    (AND
      (OR
        (AND ARG0 ARG2)
        (NAND ARG0 ARG0))
      (NAND
        (NAND ARG0 ARG2)
        (AND ARG1 ARG2)))
    (AND
      (NOR
        (NAND ARG0 ARG2)
        (OR ARG1 ARG0))
      (NAND
        (OR ARG0 ARG1)
        (OR ARG0 ARG0))))
  (AND
    (NOR
      (AND
        (AND ARG2 ARG1)
        (NAND ARG1 ARG0))
      (OR
        (AND ARG1 ARG2)
        (NAND ARG0 ARG0)))
    (OR
      (AND
        (NAND
          (OR ARG2 ARG2)
          (NAND ARG1 ARG1))
        (NOR ARG0 ARG2))
      (NAND
        (NAND ARG2 ARG0)
        (NOR ARG2 ARG0))))
)
"#;
        let adf1_str = r#"
(AND
  (OR
    (NAND
      (ADF0
        (AND ARG1 ARG0)
        (NAND ARG0 ARG2)
        (AND ARG2 ARG2))
      (NAND
        (NAND ARG1 ARG2)
        (OR ARG2 ARG1)))
    (OR
      (NOR
        (ADF0 ARG2 ARG0 ARG2)
        (OR ARG0 ARG2))
      (ADF0
        (OR ARG2 ARG2)
        (ADF0 ARG2 ARG1 ARG0)
        (NOR ARG2 ARG0))))
  (OR
    (AND
      (OR
        (ADF0 ARG0 ARG2 ARG2)
        (OR ARG0 ARG0))
      (OR
        (NOR ARG1 ARG2) ARG0))
    (ADF0
      (OR
        (NOR ARG2 ARG0)
        (NOR
          (ADF0
            (OR ARG2 ARG0)
            (NAND ARG2 ARG2)
            (ADF0 ARG0 ARG1 ARG2))
          (NAND ARG1 ARG0)))
      (ADF0
        (ADF0 ARG1 ARG1 ARG2)
        (OR ARG1 ARG1)
        (NAND ARG1 ARG2))
      (OR
        (OR ARG1 ARG2)
        (NOR ARG0 ARG0))))
)
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }


    #[cfg(gpopt_even_parity_k="5")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
(NAND
  (NAND
    (NAND
      (NAND
        (ADF1 D1
          (NOR D3 D0)
          (ADF1 D3 D4 D0 D3)
          (OR D1 D4))
        (NOR D2 D2))
      (ADF0
        (ADF0
          (OR D1 D3) D2
          (NOR D1 D2)
          (AND D0 D2))
        (AND
          (AND D0 D2)
          (ADF1 D3 D4 D0 D3))
        (AND D3 D1)
        (NAND
          (ADF0
            (ADF1 D4
              (ADF1 D0 D2
                (AND D0 D2)
                (OR D4 D0))
              (ADF1 D1
                (AND D1 D4)
                (AND D0 D4)
                (NAND D2 D4)) D0) D1
            (ADF1 D0 D2
              (AND D0 D2)
              (OR D4 D0))
            (NOR
              (NOR
                (ADF0 D1 D2 D2 D4)
                (ADF0 D3 D4 D3 D0))
              (ADF1 D0
                (AND D3 D2)
                (ADF0 D4 D4 D0 D3) D2)))
          (OR D2 D3))))
    (NOR D2 D2))
  (ADF0
    (ADF0
      (OR D1 D3) D2
      (NOR D1 D2)
      (NOR D0 D1))
    (AND
      (AND D0 D2)
      (ADF1 D3 D4 D0 D3))
    (AND D3 D1)
    (NAND
      (ADF0
        (ADF1 D4
          (ADF1 D0 D2
            (NAND D2 D4)
            (OR D4 D0))
          (ADF1 D1
            (AND D1 D4)
            (AND D0 D4)
            (NAND D2 D4)) D0) D1
        (AND D0
          (NAND
            (NAND D0 D0) D3))
        (NOR
          (ADF0
            (NOR D1 D1)
            (NOR D3 D4) D1
            (AND D0 D2))
          (ADF1 D0
            (AND D3 D2)
            (ADF0 D4 D4 D0 D3) D2)))
      (OR D2 D0)))
)
"#;
        let adf0_str = r#"
(NAND
  (OR
    (AND
      (OR ARG0 ARG3)
      (AND ARG1 ARG3))
    (OR
      (NOR ARG3 ARG1)
      (AND ARG3 ARG1)))
  (NAND
    (NOR
      (NOR ARG1 ARG2)
      (OR ARG3 ARG3))
    (NOR
      (OR ARG0 ARG2)
      (AND ARG2 ARG0)))
)
"#;
        let adf1_str = r#"
(ADF0
  (NAND
    (NOR
      (OR ARG2 ARG3)
      (NAND ARG1 ARG2))
    (NOR
      (OR ARG2 ARG3)
      (NAND ARG1 ARG2)))
  (NOR
    (NOR
      (ADF0 ARG1 ARG0 ARG2 ARG2)
      (ADF0 ARG3 ARG1 ARG2 ARG0))
    (NOR ARG2 ARG0))
  (ADF0
    (OR
      (NOR ARG3 ARG3)
      (ADF0 ARG0 ARG3 ARG3 ARG1))
    (AND
      (OR ARG3 ARG2)
      (ADF0 ARG2 ARG1 ARG0 ARG3))
    (ADF0
      (NAND ARG0 ARG2)
      (ADF0 ARG3 ARG0 ARG3 ARG3)
      (NOR ARG3 ARG1)
      (NAND ARG0 ARG1))
    (NAND
      (OR ARG0 ARG3)
      (OR ARG2 ARG1)))
  (NAND
    (NOR
      (OR ARG1 ARG2)
      (ADF0 ARG0 ARG3 ARG1 ARG1))
    (NOR
      (OR ARG2 ARG3) ARG3))
)
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }

    #[cfg(gpopt_even_parity_k="6")]
    fn get_k_n_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
(OR
  (ADF1
    (OR
      (NOR D2 D5)
      (NAND D5 D0))
    (NAND
      (OR D0 D2)
      (ADF1 D3 D5 D2 D5 D5))
    (ADF0
      (AND D2 D4)
      (ADF1 D5 D2 D4 D1 D1)
      (ADF1 D2 D1 D4 D3 D3)
      (OR D5 D5)
      (AND D5 D3))
    (OR D0 D2)
    (NAND
      (NAND D4
        (OR D2 D0))
      (ADF0 D5 D4 D3 D0 D4)))
  (ADF0
    (NOR
      (NOR D0 D0)
      (ADF1 D5 D0 D0 D0 D2))
    (ADF0
      (OR D2 D3)
      (NOR D3 D1)
      (AND D2 D3)
      (ADF0 D1 D1 D5 D1 D2)
      (NAND D5 D0))
    (OR
      (ADF1 D3 D3 D0 D0 D4)
      (NOR D0 D3))
    (NOR
      (ADF0 D3 D2 D0 D3 D5)
      (ADF0 D4 D5 D2 D1 D1))
    (ADF0
      (ADF0 D4 D4 D2 D2 D1)
      (NOR D3 D5)
      (AND D3 D0)
      (ADF1 D4 D4 D1 D1 D0)
      (ADF1 D3 D2 D1 D0 D3)))
)
"#;
        let adf0_str = r#"
(AND
  (NOR
    (NOR
      (OR ARG1 ARG3)
      (NAND ARG3 ARG3))
    (OR
      (AND ARG4 ARG3)
      (NOR ARG4 ARG3)))
  (NAND
    (NOR
      (NOR ARG2 ARG1)
      (NAND ARG3 ARG0))
    (NOR
      (AND ARG3 ARG3)
      (AND ARG4 ARG4)))
)
"#;
        let adf1_str = r#"
(OR
  (AND
    (AND
      (OR ARG0 ARG1)
      (NOR ARG4 ARG0))
    (AND
      (ADF0
        (NOR ARG3 ARG1)
        (OR ARG3 ARG3)
        (ADF0 ARG1 ARG3 ARG1 ARG0 ARG3)
        (NAND ARG4 ARG1)
        (NAND ARG3 ARG4))
      (NOR ARG4 ARG0)))
  (ADF0
    (OR
      (NAND ARG0 ARG0)
      (NAND ARG1 ARG1))
    (ADF0
      (NOR ARG3 ARG1)
      (OR ARG3 ARG3)
      (ADF0 ARG1 ARG3 ARG1 ARG0 ARG3)
      (NAND ARG4 ARG1)
      (NAND ARG3 ARG4))
    (ADF0
      (NAND ARG3 ARG3)
      (NAND ARG1 ARG0)
      (NAND ARG1 ARG3)
      (AND ARG0 ARG0)
      (NOR ARG2 ARG4))
    (OR
      (AND ARG3 ARG4) ARG3)
    (NOR
      (OR ARG1 ARG1) ARG1))
)
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }
}
