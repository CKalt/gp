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
pub const NUM_ADF: u8 = 5;
pub const ADF_ARITY: u8 = 0;

/// RunContext provides runtime control over a running individual. Each 
/// node and terminal exec call recieves a reference to its RunContext
/// where it can then access it's fitness case data and currency values.
/// Since there is a one to one correspondence between threads and 
/// running individual there will be one RunContext for each thread
/// running for an Individual. 
pub struct RunContext<'a> {
    pub opt_func_def_branches: Option<Vec<&'a TreeBranch>>, // adf0, adf1,...
    pub opt_adf_args: Option<Vec<GpType>>,
    pub cur_fc_index: usize,     // Current fitness case being processed.
    pub cur_pos: (i8, i8), // x,y position on 4x6 character grid
    pub hits: GpHits,
    pub error: GpRaw,
    pub opt_run_result: Option<IndividualRunResult>,   // Program is done when not None
}
impl<'a> RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let mut rc = RunContext {
            cur_fc_index: 0,
            cur_pos: (0, 0),
            opt_func_def_branches: None,
            opt_adf_args: None,
            hits: 0,
            error: 0,
            opt_run_result: None,
        };

        rc
    }
    pub fn cur_fc(&self) -> &FitnessCase {
        &FITNESS_CASES[self.cur_fc_index]
    }
    pub fn print_run_illustration(&self, _label: &str) { }
    pub fn prepare_run(&mut self) { }
    pub fn get_hits_label() -> &'static str {
        "num hits"
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
    /// by the adf_num arg
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

//fn fix_nan(num: GpType) -> GpType {
//    if num.is_nan() {
//        0.0
//    } else {
//        num
//    }
//}

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

#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg(rc: &RunContext, term: &Terminal) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[term.index],
    }
}

pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    let funcs = vec![
        Function {
            fid:  0u8,
            name: "IF".to_string(),
            arity: 3,
            code: function_if,
            opt_adf_num: None,
        },
        Function {
            fid:  1u8,
            name: "AND".to_string(),
            arity: 2,
            code: function_or,
            opt_adf_num: None,
        },
        Function {
            fid:  2u8,
            name: "OR".to_string(),
            arity: 2,
            code: function_or,
            opt_adf_num: None,
        },
        Function {
            fid:  3u8,
            name: "NOT".to_string(),
            arity: 1,
            code: function_not,
            opt_adf_num: None,
        },
        Function {
            fid:  4u8,
            name: "HOMING".to_string(),
            arity: 1,
            code: function_homing,
            opt_adf_num: None,
        },
    ];
    vec![funcs]
}

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
    let save_pos = rc.cur_pos;
    let val = Tree::exec_node(rc, &func.branch[0])
    rc.cur_pos = save_pos;
    val
}


pub fn get_terminals_for_result_branches() -> Vec<Vec<Terminal>> {
    // Terminal Set: I,L, NIL, X, N, NE, E, SE, S, SW, W, NW, (GO-N),
    //    (GO-NE), (GO-E), (GO-SE), (GO-S), (GO-SW), (GO-W), (G-W),
    //    (GO-NW)
    let terms = vec![
        Terminal::new(0, "I", terminal_i),
        Terminal::new(1, "L", terminal_l),
        Terminal::new(2, "NIL", terminal_nil),
        Terminal::new(3, "X", terminal_x),
        Terminal::new(4, "N", terminal_n), 
        Terminal::new(5, "NE", terminal_ne), 
        Terminal::new(6, "E", terminal_e),
        Terminal::new(7, "SE", terminal_se),
        Terminal::new(8, "S", terminal_s),
        Terminal::new(9, "SW", terminal_sw),
        Terminal::new(10, "W", terminal_w),
        Terminal::new(11, "NW", terminal_nw),
        Terminal::new(12, "(GO-N)", terminal_go_n),
        Terminal::new(13, "(GO-NE)", terminal_go_ne),
        Terminal::new(14, "(GO-E)", terminal_go_e),
        Terminal::new(15, "(GO-SE)", terminal_go_se),
        Terminal::new(16, "(GO-S)", terminal_go_s),
        Terminal::new(17, "(GO-SW)", terminal_go_sw),
        Terminal::new(18, "(GO-W)", terminal_go_w),
        Terminal::new(19, "(GO-NW)", terminal_go_nw),
    ];
    vec![terms]
}

/// Ends program with run_result indicating Letter I
fn terminal_i(rc: &mut RunContext, term: &Terminal) -> GpType {
    rc.opt_run_result = Some(IsLetter('I'));
    true
}

/// Ends program with run_result indicating Letter L
fn terminal_l(rc: &mut RunContext, term: &Terminal) -> GpType {
    rc.opt_run_result = Some(IsLetter('L'));
    true
}

/// Ends program with run_result indicating Not Letter
fn terminal_nil(rc: &mut RunContext, term: &Terminal) -> GpType {
    rc.run_result = NotLetter;
    true
}

/// Returns pixel at current position for current fitness case
fn terminal_x(rc: &mut RunContext, term: &Terminal) -> GpType {
    rc.cur_fc().get_absolute_pixel(rc.cur_pos)
}

fn terminal_n(rc: &mut RunContext, term: &Terminal) -> GpType {
    rc.cur_fc().get_relative_pixel(rc, (0, -1))
}

fn terminal_ne(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (1, -1))
}

fn terminal_e(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (1, 0))
}

fn terminal_se(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (1, 1))
}

fn terminal_s(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (0, 1))
}

fn terminal_sw(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (-1, 1))
}

fn terminal_w(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (-1, 0))
}

fn terminal_nw(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().get_relative_pixel(rc, (-1, -1))
}

fn terminal_go_n(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (0, -1))
}

fn terminal_go_ne(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (1, -1))
}

fn terminal_go_e(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (1, 0))
}

fn terminal_go_se(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (1, 1))
}

fn terminal_go_s(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (0, 1))
}

fn terminal_go_sw(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (-1, 1))
}

fn terminal_go_w(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (-1, 0))
}

fn terminal_go_nw(fc: &mut RunContext, term: &Terminal) -> GpType {
    rc.get_fc().move_relative_pixel(rc, (-1, -1))
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

// FitnessCase
type FitnessCase = [[bool; 4]; 6];
impl FitnessCaseTrait for FitnessCase {
    /// get pixel value for fixed point at `pos` tuple (x,y)
    fn get_absolute_pixel(&self, pos(i8, i8)) -> bool {
        self[pos.1 as usize][pos.0 as usize]
    }
    /// get pixel relatative to RunContext cur_pos
    fn get_relative_pixel(&self, rc: &RunContext, delta(i8, i8)) -> bool {
        let pos = (rc.cur_pos.0 + delta.0, rc.cur_pos.1 + delta.1);
        self.get_pixel(pos)
    }
    /// move cur_pos for RunContext to point relative to its current position
    /// get pixel at that new location.
    fn move_relative_pixel(&self, rc: &mut RunContext, delta(i8, i8)) -> bool {
        rc.cur_pos.0 += delta.0;
        rc.cur_pos.1 += delta.1;
        self.get_pixel(rc.cur_pos)
    }
}

struct FitnessCases {
    fc: [FitnessCase; 78],
}

const FITNESS_CASES: FitnessCases = FitnessCases {
    fc: [
            [ // 00
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 01
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 02
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 03
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 04
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 05
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 06
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, true, true],
                    [ false, false, false, false],
            ],
            [ // 07
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, true],
                    [ false, false, false, false],
            ],
            [ // 08
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, false, false, false],
            ],
            [ // 09
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, true, false, false],
            ],
            [ // 10
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, true, false],
            ],
            [ // 11
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 12
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 13
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 14
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 15
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, true],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 16
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 17
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 18
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ true, false, false, false],
            ],
            [ // 19
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, true],
            ],
            [ // 20
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 21
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 22
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 23
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ true, true, true, true],
            ],
            [ // 24
                    [ false, false, true, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, true],
            ],
            [ // 25
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, true, false, false],
            ],
            [ // 26
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, true],
                    [ false, true, false, false],
            ],
            [ // 27
                    [ false, true, false, true],
                    [ false, true, false, true],
                    [ false, true, false, true],
                    [ false, true, false, true],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 28
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, true],
                    [ false, false, false, false],
            ],
            [ // 29
                    [ false, true, false, true],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 30
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ true, false, false, false],
            ],
            [ // 31
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 32
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, true, false],
                    [ false, true, true, true],
                    [ true, false, false, false],
            ],
            [ // 33
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 34
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 35
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 36
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 37
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 38
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
            ],
            [ // 39
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, true, false],
            ],
            [ // 40
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 41
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 42
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 43
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 44
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, false, false, false],
            ],
            [ // 45
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 46
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 47
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 48
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 49 (Note: Koza II book shows copy of 45 here which seems like an error)
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 50
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 51
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
            ],
            [ // 52
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, true, true],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 53
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, true, false],
                    [ false, false, false, false],
            ],
            [ // 54
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, false, true, false],
            ],
            [ // 55
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, true, false, true],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 56
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ true, true, true, true],
            ],
            [ // 57
                    [ false, false, true, false],
                    [ false, true, false, false],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, false, false, true],
            ],
            [ // 58
                    [ false, true, true, true],
                    [ false, true, false, true],
                    [ false, true, false, true],
                    [ false, true, false, true],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 59
                    [ false, true, false, false],
                    [ true, true, true, false],
                    [ false, true, false, false],
                    [ true, true, true, false],
                    [ false, true, false, false],
                    [ false, false, false, false],
            ],
            [ // 60
                    [ true, false, true, false],
                    [ true, false, true, false],
                    [ true, true, true, false],
                    [ false, false, true, false],
                    [ false, false, true, false],
                    [ false, false, false, false],
            ],
            [ // 61
                    [ false, false, true, false],
                    [ false, false, true, false],
                    [ false, false, true, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, true, false],
            ],
            [ // 62
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ true, true, true, true],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 63
                    [ false, false, true, false],
                    [ false, false, false, true],
                    [ true, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, true, false, false],
            ],
            [ // 64
                    [ true, false, true, false],
                    [ false, true, false, true],
                    [ true, false, true, false],
                    [ false, true, false, true],
                    [ true, false, true, false],
                    [ false, true, false, true],
            ],
            [ // 65
                    [ false, true, false, true],
                    [ true, false, true, false],
                    [ false, true, false, true],
                    [ true, false, true, false],
                    [ false, true, false, true],
                    [ true, false, true, false],
            ],
            [ // 66
                    [ false, false, false, true],
                    [ false, false, false, true],
                    [ false, false, false, true],
                    [ false, false, false, true],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 67
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, false, false, false],
            ],
            [ // 68
                    [ false, false, true, true],
                    [ false, false, true, true],
                    [ false, false, true, true],
                    [ false, false, true, true],
                    [ false, false, true, true],
                    [ false, false, false, false],
            ],
            [ // 69
                    [ false, true, true, true],
                    [ false, true, true, true],
                    [ false, true, true, true],
                    [ false, true, true, true],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 70
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 71
                    [ false, false, false, false],
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 72
                    [ true, true, true, true],
                    [ true, true, true, true],
                    [ true, true, true, true],
                    [ true, true, true, true],
                    [ true, true, true, true],
                    [ true, true, true, true],
            ],
            [ // 73
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
                    [ false, false, false, false],
            ],
            [ // 74
                    [ false, true, true, true],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, true, false],
                    [ false, true, true, true],
                    [ false, false, false, false],
            ],
            [ // 75
                    [ false, true, true, false],
                    [ false, true, true, false],
                    [ false, true, false, false],
                    [ false, true, false, true],
                    [ false, true, true, true],
                    [ true, false, false, false],
            ],
            [ // 76
                    [ false, true, true, false],
                    [ false, true, false, true],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, true, true, true],
                    [ true, false, false, false],
            ],
            [ // 77
                    [ false, false, true, false],
                    [ false, true, false, true],
                    [ true, true, false, false],
                    [ false, true, false, false],
                    [ false, false, true, true],
                    [ true, false, false, false],
            ],
    ],
};

enum IndividualRunResult {
    IsLetter(char),
    NotLetter,
};

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
