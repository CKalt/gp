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

fn function_if(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    if Tree::exec_node(rc, &func.branch[0]) {
        Tree::exec_node(rc, &func.branch[1])
    } else {
        Tree::exec_node(rc, &func.branch[2])
    }
}

fn function_and(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    if Tree::exec_node(rc, &func.branch[0]) {
        Tree::exec_node(rc, &func.branch[1]) 
    } else {
        false
    }
}

fn function_or(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    if Tree::exec_node(rc, &func.branch[0]) {
        true
    } else {
        Tree::exec_node(rc, &func.branch[1]) 
    }
}

fn function_not(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    !Tree::exec_node(rc, &func.branch[0])
}

fn function_homing(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    Tree::exec_node(rc, &func.branch[0])
}

#[cfg(gpopt_adf="yes")]
fn function_adf(rc: &mut RunContext, func: &FunctionNode) -> GpType {
    let mut args: Vec<GpType> = Vec::new();
    for b_i in 0..ADF_ARITY {
        args.push(Tree::exec_node(rc, &func.branch[b_i]));
    }

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, &args)
}

fn terminal_data(rc: &RunContext, term: &Terminal) -> GpType {
    rc.get_cur_fc().input_bits[term.index]
}

fn terminal_adf_arg(rc: &RunContext, term: &Terminal) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[term.index],
    }
}

pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    let mut funcs = vec![
        Function {
            fid:  0u8,
            name: "AND".to_string(),
            arity: 2,
            code: function_and,
            opt_adf_num: None,
        },
        Function {
            fid:  1u8,
            name: "OR".to_string(),
            arity: 2,
            code: function_or,
            opt_adf_num: None,
        },
        Function {
            fid:  2u8,
            name: "NAND".to_string(),
            arity: 2,
            code: function_nand,
            opt_adf_num: None,
        },
        Function {
            fid:  3u8,
            name: "NOR".to_string(),
            arity: 2,
            code: function_nor,
            opt_adf_num: None,
        },
    ];

    for a_i in 0..NUM_ADF {
        funcs.push(
            Function {
                fid:  a_i as u8 + 4u8,
                name: format!("ADF{}", a_i),
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
    let mut terms: Vec<Terminal> = Vec::new();
    for tid in 0..EVEN_PARITY_K_VALUE {
        terms.push(
            Terminal {
                tid:  tid as u8,
                name: format!("D{}", tid),
                code: terminal_data,
                index: tid,
            },
        );
    }
    vec![terms]
}

pub fn get_functions_for_func_def_branches() -> Vec<Vec<Function>> {
    let mut branches: Vec<Vec<Function>> = Vec::new();
    for b_i in 0..NUM_ADF {
        let mut funcs = vec![
            Function {
                fid:  0u8,
                name: "AND".to_string(),
                arity: 2,
                code: function_and,
                opt_adf_num: None,
            },
            Function {
                fid:  1u8,
                name: "OR".to_string(),
                arity: 2,
                code: function_or,
                opt_adf_num: None,
            },
            Function {
                fid:  2u8,
                name: "NAND".to_string(),
                arity: 2,
                code: function_nand,
                opt_adf_num: None,
            },
            Function {
                fid:  3u8,
                name: "NOR".to_string(),
                arity: 2,
                code: function_nor,
                opt_adf_num: None,
            },
        ];
        for f_i in 0..b_i {
            funcs.push(
                Function {
                    fid:  4u8 + f_i as u8,
                    name: format!("ADF{}", f_i),
                    arity: ADF_ARITY as u8,
                    code: function_adf,
                    opt_adf_num: Some(f_i),
                },
            );
        }

        branches.push(funcs);
    }
    branches
}

pub fn get_terminals_for_func_def_branches() -> Vec<Vec<Terminal>> {
    let mut branches: Vec<Vec<Terminal>> = Vec::new();
    for _b_i in 0..NUM_ADF {
        let mut terms = Vec::new();
        for t_i in 0..ADF_ARITY {
            terms.push(
                Terminal {
                    tid:  t_i as u8,
                    name: format!("ARG{}", t_i),
                    code: terminal_adf_arg,
                    index: t_i
                },
            );
        }
        branches.push(terms);
    }
    branches
}

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
    pub cur_pos: (u8, u8), // x,y position on 4x6 character grid
    pub hits: GpHits,
    pub error: GpRaw,
}
impl<'a> RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let mut rc = RunContext {
            cur_fc: 0,
            cur_pos: (0, 0),
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
            None => {
                let passed_args = args.clone();
                self.opt_adf_args = Some(passed_args);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = None;
                result
            },
            Some(ref orig_args) => {
                let save_args = orig_args.clone();
                let passed_args = args.clone();
                self.opt_adf_args = Some(passed_args);
                let result = Tree::exec_node(self, &func_def_branch.root);
                self.opt_adf_args = Some(save_args);
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
    /// It is designed to accept a tree that has been already evolved
    /// and declared a winner. It will parse the tree and run
    /// the tree for all fitness cases to confirm that it is indeed 
    /// a valid tree that wins.
    /// 
    /// Success confirms that the parser works and has the ability
    /// to run the tree independently against the complete set of
    /// fitness cases.
    /// It also confirms the fitness measure correctly establishes
    /// tree is a Winner.
    pub fn test_eval_one_tree() {
        use crate::tree::*;
        // first build tree:
        let (rb_str, fd_vec) = get_k_5_test();
        let mut tree = Tree::parse(rb_str, fd_vec);

        assert_eq!(tree.eval_one_tree(), true);
    }

    fn get_k_5_test() -> (&'static str, Vec<&'static str>) {
        let rb_str = r#"
(ADF0
  (ADF0
    (AND D3
      (OR D0
        (OR
          (AND D4 D1)
          (OR D3 D4))))
    (ADF0
      (OR D0 D0)
      (ADF1 D2 D1)))
  (ADF1
    (ADF0
      (NAND D1 D1)
      (ADF1 D4 D1))
    (OR
      (NOR D2 D1)
      (NAND D3 D4)))
)
"#;
        let adf0_str = r#"
(AND
  (NAND
    (NAND
      (AND ARG0 ARG1)
      (AND ARG0 ARG1))
    (OR
      (AND ARG1 ARG1)
      (OR ARG0 ARG0)))
  (NAND
    (AND
      (OR ARG0 ARG0)
      (NOR ARG0 ARG1))
    (NOR
      (AND ARG0 ARG1)
      (AND ARG0 ARG0)))
)
"#;
        let adf1_str = r#"
(ADF0
  (ADF0
    (NAND
      (ADF0 ARG1 ARG0)
      (ADF0 ARG1 ARG0))
    (AND
      (NOR ARG1 ARG1)
      (ADF0 ARG0 ARG1)))
  (OR ARG1
    (OR
      (NAND ARG0 ARG0)
      (NOR ARG0 ARG0)))
)
"#;
        (rb_str, vec![adf0_str, adf1_str])
    }
}
