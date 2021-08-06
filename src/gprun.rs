use crate::fitness::Fitness;
use crate::tree::Tree;
use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::tree::TreeBranch;
use crate::control::CONTROL;

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

fn function_adf0(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let arg1 = Tree::exec_node(rc, &func.branch[0]);
    let arg2 = Tree::exec_node(rc, &func.branch[1]);

    rc.exec_adf0(arg1, arg2)
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

fn terminal_arg0(rc: &RunContext) -> GpType {
     match rc.opt_adf0_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[0],
    }
}

fn terminal_arg1(rc: &RunContext) -> GpType {
    match rc.opt_adf0_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[1],
    }
}

pub static FUNCTIONS_RESULT_BRANCH: [Function; CONTROL.num_functions_result_branch as usize] = [
    Function {
        fid:  0u8,
        name: "ADF0",
        arity: 2,
        code: function_adf0,
    },
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

pub static FUNCTIONS_FUNC_DEF_BRANCH: [Function; CONTROL.num_functions_func_def_branch as usize] = [
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

// TERMINAL SPECIFICS - RESULT PRODUCING BRANCH - result_branch
pub static TERMINALS_RESULT_BRANCH: [Terminal; CONTROL.num_terminals_result_branch as usize] = [
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
];

pub const RUN_CONTROL_NUM_FITNESS_CASES: usize = 64;

// FitnessCase
pub struct FitnessCase {
    input_bits: [bool; 6],
    output_bit: bool,
}
impl FitnessCase {
    fn new() -> FitnessCase {
        FitnessCase {
            input_bits: [false; 6],
            output_bit: false,
        }
    }
}

/// RunContext provides runtime control over a running individual. Each 
/// node and terminal exec call recieves a reference to its RunContext
/// where it can then access it's fitness case data and currency values.
pub struct RunContext<'a> {
    pub fitness_cases: Vec::<FitnessCase>,
    pub opt_func_def_branch: Option<&'a TreeBranch>, // used for calling adf0
    pub opt_adf0_args: Option<Vec<GpType>>,
    pub cur_fc: usize,
    pub hits: GpHits,
    pub error: GpRaw,
}
impl RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let mut rc = RunContext {
            cur_fc: 0,
            opt_func_def_branch: None,
            opt_adf0_args: None,
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
    pub fn exec_adf0(&mut self, arg1: GpType, arg2: GpType)
            -> GpType {
        let func_def_branch =
            self.opt_func_def_branch.expect("branch not assigned for exec_adf0");

        match self.opt_adf0_args {
            None    => {
                self.opt_adf0_args = Some(vec![ arg1, arg2 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf0_args = None;
                result
            },
            Some(ref args) => {
                let (a,b) = (args[0], args[1]);
                self.opt_adf0_args = Some(vec![ arg1, arg2 ]);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf0_args = Some(vec![ a, b ]);
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
                                    &FUNCTIONS_RESULT_BRANCH,
                                    &TERMINALS_RESULT_BRANCH), // result branch 0
                               Some(("(AND ARG0 ARG1)",
                                    &FUNCTIONS_FUNC_DEF_BRANCH,
                                    &TERMINALS_FUNC_DEF_BRANCH))); // func def branch 0
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
                                    &FUNCTIONS_RESULT_BRANCH,
                                    &TERMINALS_RESULT_BRANCH), // result branch 0
                               Some((fd0_s,
                                    &FUNCTIONS_FUNC_DEF_BRANCH,
                                    &TERMINALS_FUNC_DEF_BRANCH))); // func def branch 0
        assert_eq!(tree.print_exec_one(), true);

    }
}
