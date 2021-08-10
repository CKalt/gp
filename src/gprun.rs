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
fn function_adf0(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);

    let func_def_branch =
        self.opt_func_def_branch_adf0
            .expect("branch not assigned for exec_adf0");

    rc.exec_adf(func_def_branch, arg1, arg2)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="4")]
fn function_adf0(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);

    let func_def_branch =
        self.opt_func_def_branch_adf0
            .expect("branch not assigned for exec_adf0");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="5")]
fn function_adf0(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);

    let func_def_branch =
        self.opt_func_def_branch_adf0
            .expect("branch not assigned for exec_adf0");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3, arg4)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="6")]
fn function_adf0(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);
    let arg5 = Tree::exec_node(rc, &func.branch[4]);

    let func_def_branch =
        self.opt_func_def_branch_adf0
            .expect("branch not assigned for exec_adf0");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3, arg4, arg5)
}



///  
///  Note that the following set of functions defined for adf1
///  vary from the set defined for adf0 only by the fact that
///  self.opt_func_def_branch_adf1 is used as the func_def_branch.
///  An improved architecture being considered is to have the 
///  func_def_branch made as an optional variable with the
///  func FunctionNode struct, that way only one set will be 
///  required.
///  
///  This almost certainly will be the architecture used for point
///  typing when the number of args will also vary by individual tree.

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="3")]
fn function_adf1(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);

    let func_def_branch =
        self.opt_func_def_branch_adf1
            .expect("branch not assigned for exec_adf1");

    rc.exec_adf(func_def_branch, arg1, arg2)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="4")]
fn function_adf1(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);

    let func_def_branch =
        self.opt_func_def_branch_adf1
            .expect("branch not assigned for exec_adf1");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="5")]
fn function_adf1(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);

    let func_def_branch =
        self.opt_func_def_branch_adf1
            .expect("branch not assigned for exec_adf1");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3, arg4)
}

#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="6")]
fn function_adf1(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);
    let arg3 = Tree::exec_node(rc, &func.branch[2]);
    let arg4 = Tree::exec_node(rc, &func.branch[3]);
    let arg5 = Tree::exec_node(rc, &func.branch[4]);

    let func_def_branch =
        self.opt_func_def_branch_adf1
            .expect("branch not assigned for exec_adf1");

    rc.exec_adf(func_def_branch, arg1, arg2, arg3, arg4, arg5)
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

# The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf0_arg0(rc: &RunContext) -> GpType {
     match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[0],
    }
}

# The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
fn terminal_adf0_arg1(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[1],
    }
}

# The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_adf0_arg2(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[2],
    }
}
            
# The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
fn terminal_adf0_arg3(rc: &RunContext) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[3],
    }
}

# The number of args for adf for even k parity is k-1.
#[cfg(gpopt_adf="yes")]
#[cfg(gpopt_even_parity_k="6")]
fn terminal_adf0_arg4(rc: &RunContext) -> GpType {
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
    },
    Function {
        fid:  2u8,
        name: "OR",
        arity: 2,
        code: function_or,
    },
    Function {
        fid:  3u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
    },
    Function {
        fid:  4u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_RESULT_BRANCH_ADF: [Function; 6] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
    },
    Function {
        fid:  4u8,
        name: "ADF0",
        arity: 2,
        code: function_adf0,
    },
    Function {
        fid:  5u8,
        name: "ADF1",
        arity: 2,
        code: function_adf1,
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_FUNC_DEF_BRANCH_ADF0: [Function; 4] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
    },
];

#[cfg(gpopt_adf="yes")]
pub static FUNCTIONS_FUNC_DEF_BRANCH_ADF1: [Function; 5] = [
    Function {
        fid:  0u8,
        name: "AND",
        arity: 2,
        code: function_and,
    },
    Function {
        fid:  1u8,
        name: "OR",
        arity: 2,
        code: function_or,
    },
    Function {
        fid:  2u8,
        name: "NAND",
        arity: 2,
        code: function_nand,
    },
    Function {
        fid:  3u8,
        name: "NOR",
        arity: 2,
        code: function_nor,
    },
    Function {
        fid:  4u8,
        name: "ADF0",
        arity: 2,
        code: function_adf0,
    },
];

// TERMINAL SPECIFICS - RESULT PRODUCING BRANCH - result_branch

#[cfg(gpopt_even_parity_k="3")]
const EVEN_PARITY_K_VALUE: usize = 3;
#[cfg(gpopt_even_parity_k="4")]
const EVEN_PARITY_K_VALUE: usize = 4;
#[cfg(gpopt_even_parity_k="5")]
const EVEN_PARITY_K_VALUE: usize = 5;
#[cfg(gpopt_even_parity_k="6")]
const EVEN_PARITY_K_VALUE: usize = 6;

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
#[cfg(gpopt_adf="yes")]
pub static TERMINALS_FUNC_DEF_BRANCH: [Terminal; EVEN_PARITY_K_VALUE-1] = [
    Terminal {
        tid:  0u8,
        name: "ARG0",
        code: terminal_adf0_arg0,
    },
    Terminal {
        tid:  1u8,
        name: "ARG1",
        code: terminal_adf0_arg1,
    },
#[cfg(any(gpopt_even_parity_k="4",
          gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  2u8,
        name: "ARG2",
        code: terminal_adf0_arg2,
    },
#[cfg(any(gpopt_even_parity_k="5",
          gpopt_even_parity_k="6"))]
    Terminal {
        tid:  3u8,
        name: "ARG3",
        code: terminal_adf0_arg3,
    },
#[cfg(gpopt_even_parity_k="6")]
    Terminal {
        tid:  4u8,
        name: "ARG4",
        code: terminal_adf0_arg4,
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
    pub opt_func_def_branch_adf0: Option<&'a TreeBranch>, // used for calling adf0
    pub opt_func_def_branch_adf1: Option<&'a TreeBranch>, // used for calling adf1
    pub opt_adf_args: Option<Vec<GpType>>,
    pub cur_fc: usize,
    pub hits: GpHits,
    pub error: GpRaw,
}
impl RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let mut rc = RunContext {
            cur_fc: 0,
            opt_func_def_branch_adf0: None,
            opt_func_def_branch_adf1: None,
            opt_adf_args: None,
            fitness_cases: Vec::new(),
            hits: 0,
            error: 0
        };

        for i in 0..RUN_CONTROL_NUM_FITNESS_CASES {
            let mut fc = FitnessCase::new();
            for bit in 0..6 {
                fc.input_bits[bit] = (i & 2u8.pow(bit as u32) as usize) != 0;
            }
            fc.output_bit = Self::is_bool6_sym(i as u8);
            rc.fitness_cases.push(fc);
        }

        rc
    }
    fn is_bool6_sym(val: u8) -> bool {
        // bit @pos 2^0 must match bit @pos 2^5
        // bit @pos 2^1 must match bit @pos 2^4
        // bit @pos 2^2 must match bit @pos 2^3

        ((val & 2u8.pow(0)) != 0) == ((val & 2u8.pow(5)) != 0) &&
        ((val & 2u8.pow(1)) != 0) == ((val & 2u8.pow(4)) != 0) &&
        ((val & 2u8.pow(2)) != 0) == ((val & 2u8.pow(3)) != 0)
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
    pub fn exec_adf(&mut self, func_def_branch: &'a TreeBranch, arg1: GpType,
                arg2: GpType)
            -> GpType {
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
    pub fn exec_adf(&mut self, func_def_branch: &'a TreeBranch, arg1: GpType,
            arg2: GpType, arg3: GpType)
            -> GpType {
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
    pub fn exec_adf(&mut self, func_def_branch: &'a TreeBranch, arg1: GpType,
                arg2: GpType, arg3: GpType, arg4: GpType)
            -> GpType {
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
    pub fn exec_adf(&mut self, func_def_branch: &'a TreeBranch, arg1: GpType,
                arg2: GpType, arg3: GpType, arg4: GpType, arg5: GpType)
            -> GpType {
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
    use crate::Tree;
    use super::*;

    #[test]
    pub fn test_print_exec_one() {
        // first build this test tree:
        // func def branch:
        //     (AND ARG0 ARG1)
        // result branch:
        //     (OR (ADF0 true true) false)
        // exec_tree s/b: true
        let mut tree = Tree::parse(("(OR (ADF0 D0 D1) D2)",
                                    CONTROL.funcs_rpb[0],
                                    CONTROL.terms_rpb[0]), // result branch 0
                               Some(("(AND ARG0 ARG1)",
                                    CONTROL.funcs_fdb[0],
                                    CONTROL.terms_fdb[0]))); // func def branch 0
        assert_eq!(tree.print_exec_one(), false);

        let rb0_s = r#"
(AND
  (NOR
    (OR
      (ADF0 D3 D2)
      (ADF0 D4 D1))
    (NOR D0
      (NAND D0 D5)))
  (NOR
    (NOR
      (NOR D0 D5)
      (AND D0 D5))
    (AND
      (OR
        (ADF0 D3 D2)
        (ADF0 D4 D1))
      (OR D3 D2)))
)
        "#;

        let fd0_s = r#"
(AND
  (AND
    (AND
      (NAND ARG1 ARG0)
      (OR ARG0 ARG1))
    (OR
      (NAND ARG1 ARG0)
      (NOR ARG0 ARG0)))
  (NAND
    (OR
      (NOR ARG1 ARG0)
      (NOR
        (NOR ARG0
          (NAND
            (OR ARG1 ARG1)
            (NOR ARG0 ARG0)))
        (OR
          (NOR
            (OR
              (NOR
                (OR ARG1 ARG1)
                (NOR
                  (OR
                    (NAND ARG0 ARG0) ARG0)
                  (NAND ARG0 ARG1)))
              (NOR
                (NAND ARG0 ARG1)
                (NOR ARG0 ARG0)))
            (NAND ARG1 ARG0))
          (NAND
            (AND ARG0 ARG0)
            (NAND ARG1 ARG0)))))
    (NOR ARG1 ARG0))
)
        "#;

        let mut tree = Tree::parse((rb0_s,
                                    CONTROL.funcs_rpb[0],
                                    CONTROL.terms_rpb[0]), // result branch 0
                               Some((fd0_s,
                                    CONTROL.funcs_fdb[0],
                                    CONTROL.terms_fdb[0]))); // func def branch 0
        assert_eq!(tree.print_exec_one(), true);

    }
}
