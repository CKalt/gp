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
    
#[cfg(gpopt_adf="yes")]
pub const NUM_ADF: u8 = 5;

#[cfg(gpopt_adf="yes")]
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
    pub cur_fc_index: usize,   // Index of current fitness case being processed.
    pub cur_pos: (i8, i8), // x,y position on 4x6 character grid
    pub hits: GpHits,
    pub error: GpRaw,
    pub opt_run_result: Option<IndividualRunResult>, // Program is done when not None
}
impl<'a> RunContext<'_> {
    pub fn new() -> RunContext<'static> {
        let rc = RunContext {
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
    pub fn init_new_fitness_case(&mut self, fc_i: usize) {
        self.cur_fc_index = fc_i;
        self.opt_run_result = None;
        self.cur_pos = (0, 0);
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
        let max_possible_hits = FITNESS_CASES.fc.len() as GpHits;
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
    // a value of 0 is treated as a hit
    pub fn compute_error(&self, _result: GpType) -> GpRaw {
        if let Some(run_result) = &self.opt_run_result {
            match run_result {
                IsLetter('I') => {
                    if self.cur_fc_index == FITNESS_CASE_I_INDEX {
                        0     // true positive
                    } else {
                        1     // false positive or wrong positive
                    }
                },
                IsLetter('L') => {
                    if self.cur_fc_index == FITNESS_CASE_L_INDEX {
                        0   // true positive
                    } else {
                        1   // false positive or wrong positive
                    }
                },
                _ => {
                    if self.cur_fc_index == FITNESS_CASE_I_INDEX {
                        23  // false negative
                    } else if self.cur_fc_index == FITNESS_CASE_L_INDEX {
                        23  // false negative
                    } else {
                        0   // true negative
                    }
                }
            }
        } else {
            // treat like Nil
            if self.cur_fc_index == FITNESS_CASE_I_INDEX {
                23  // false negative
            } else if self.cur_fc_index == FITNESS_CASE_L_INDEX {
                23  // false negative
            } else {
                0   // true negative
            }
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
        args.push(Tree::exec_node(rc, &func.branch[b_i as usize]));
    }

    let adf_num =
        func.fnc.opt_adf_num.expect("branch not assigned for exec_adf0");

    rc.exec_adf(adf_num, &args)
}

#[cfg(not(gpopt_adf_arity="0"))]
#[cfg(gpopt_adf="yes")]
fn terminal_adf_arg(rc: &RunContext, term: &Terminal) -> GpType {
    match rc.opt_adf_args {
        None => panic!("terminal arg access with empty args list"),
        Some(ref args) => args[term.index],
    }
}

#[cfg(gpopt_adf="no")]
pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    let funcs = vec![
        Function {
            fid:  0u8,
            name: "IF".to_string(),
            arity: 3,
            code: function_if,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: Some((
                vec![
                    vec![1,2,3,4],   // arg0: AND,OR,NOT,HOMING
                    vec![0],         // arg1: IF
                    vec![0],         // arg2: IF
                ],
                vec![
                    vec![12,13,14,15,16,17,18,19], // arg0: GON, GONE...GONW
                    vec![0,1,2],                   // arg1: I,L,NIL
                    vec![0,1,2],                   // arg2: I,L,NIL
                ])), 
        },
        Function {
            fid:  1u8,
            name: "AND".to_string(),
            arity: 2,
            code: function_and,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  2u8,
            name: "OR".to_string(),
            arity: 2,
            code: function_or,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  3u8,
            name: "NOT".to_string(),
            arity: 1,
            code: function_not,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  4u8,
            name: "HOMING".to_string(),
            arity: 1,
            code: function_homing,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
    ];

    vec![funcs]
}

#[cfg(gpopt_adf="yes")]
pub fn get_functions_for_result_branches() -> Vec<Vec<Function>> {
    let mut funcs = vec![
        Function {
            fid:  0u8,
            name: "IF".to_string(),
            arity: 3,
            code: function_if,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: Some((
                vec![
                    vec![1,2,3,4],   // arg0: AND,OR,NOT,HOMING
                    vec![0],         // arg1: IF
                    vec![0],         // arg2: IF
                ],
                vec![
                    vec![12,13,14,15,16,17,18,19], // arg0: GON, GONE...GONW
                    vec![0,1,2],                   // arg1: I,L,NIL
                    vec![0,1,2],                   // arg2: I,L,NIL
                ])), 
        },
        Function {
            fid:  1u8,
            name: "AND".to_string(),
            arity: 2,
            code: function_and,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  2u8,
            name: "OR".to_string(),
            arity: 2,
            code: function_or,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  3u8,
            name: "NOT".to_string(),
            arity: 1,
            code: function_not,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
        Function {
            fid:  4u8,
            name: "HOMING".to_string(),
            arity: 1,
            code: function_homing,
            opt_adf_num: None,
            #[cfg(gpopt_syntactic_constraints="yes")] 
            opt_incl_constraints: None,
        },
    ];

    for adf_num in 0..NUM_ADF {
        funcs.push(
            Function {
                fid:  5u8 + adf_num,
                name: format!("ADF{}", adf_num),
                arity: 0,
                code: function_adf,
                opt_adf_num: Some(adf_num as usize),
                #[cfg(gpopt_syntactic_constraints="yes")] 
                opt_incl_constraints: None,
            }
        );
    }
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
    let val = Tree::exec_node(rc, &func.branch[0]);
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
        // TODO: improve parser code to accept hyphens embeded in symbol names
        // and then improve below names (e.g. GON => GO-N)
        Terminal::new(12, "GON", terminal_go_n),
        Terminal::new(13, "GONE", terminal_go_ne),
        Terminal::new(14, "GOE", terminal_go_e),
        Terminal::new(15, "GOSE", terminal_go_se),
        Terminal::new(16, "GOS", terminal_go_s),
        Terminal::new(17, "GOSW", terminal_go_sw),
        Terminal::new(18, "GOW", terminal_go_w),
        Terminal::new(19, "GONW", terminal_go_nw),
    ];
    vec![terms]
}

/// assigns final run_result for program
/// TODO : terminate program since result is known.
fn set_run_result(rc: &mut RunContext, rr: IndividualRunResult) {
    // set result if not already set, otherwise do nothing
    if let None = rc.opt_run_result {
        rc.opt_run_result = Some(rr);
    }
}

fn terminal_i(rc: &mut RunContext, _term: &Terminal) -> GpType {
    set_run_result(rc, IsLetter('I'));
    true
}

/// Ends program with run_result indicating Letter L
fn terminal_l(rc: &mut RunContext, _term: &Terminal) -> GpType {
    set_run_result(rc, IsLetter('L'));
    true
}

/// Ends program with run_result indicating Not Letter
fn terminal_nil(rc: &mut RunContext, _term: &Terminal) -> GpType {
    set_run_result(rc, NotLetter);
    true
}

/// Returns pixel at current position for current fitness case
fn terminal_x(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_absolute_pixel(rc.cur_pos)
}

fn terminal_n(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (0, -1))
}

fn terminal_ne(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (1, -1))
}

fn terminal_e(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (1, 0))
}

fn terminal_se(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (1, 1))
}

fn terminal_s(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (0, 1))
}

fn terminal_sw(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (-1, 1))
}

fn terminal_w(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (-1, 0))
}

fn terminal_nw(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.get_relative_pixel(rc, (-1, -1))
}

fn terminal_go_n(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (0, -1))
}

fn terminal_go_ne(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (1, -1))
}

fn terminal_go_e(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (1, 0))
}

fn terminal_go_se(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (1, 1))
}

fn terminal_go_s(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (0, 1))
}

fn terminal_go_sw(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (-1, 1))
}

fn terminal_go_w(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (-1, 0))
}

fn terminal_go_nw(rc: &mut RunContext, _term: &Terminal) -> GpType {
    let cur_fc = &FITNESS_CASES.fc[rc.cur_fc_index];
    cur_fc.move_relative_pixel(rc, (-1, -1))
}

#[cfg(gpopt_adf="no")]
pub fn get_functions_for_func_def_branches() -> Vec<Vec<Function>> {
    let branches: Vec<Vec<Function>> = Vec::new();
    branches
}

#[cfg(gpopt_adf="no")]
pub fn get_terminals_for_func_def_branches() -> Vec<Vec<Terminal>> {
    let branches: Vec<Vec<Terminal>> = Vec::new();
    branches
}

#[cfg(gpopt_adf="yes")]
pub fn get_functions_for_func_def_branches() -> Vec<Vec<Function>> {
    let mut branches: Vec<Vec<Function>> = Vec::new();
    for _b_i in 0..NUM_ADF {
        let funcs = vec![
            Function {
                fid:  0u8,
                name: "AND".to_string(),
                arity: 2,
                code: function_or,
                opt_adf_num: None,
                #[cfg(gpopt_syntactic_constraints="yes")] 
                opt_incl_constraints: None,
            },
            Function {
                fid:  1u8,
                name: "OR".to_string(),
                arity: 2,
                code: function_or,
                opt_adf_num: None,
                #[cfg(gpopt_syntactic_constraints="yes")] 
                opt_incl_constraints: None,
            },
            Function {
                fid:  2u8,
                name: "NOT".to_string(),
                arity: 1,
                code: function_not,
                opt_adf_num: None,
                #[cfg(gpopt_syntactic_constraints="yes")] 
                opt_incl_constraints: None,
            },
        ];

        branches.push(funcs);
    }
    branches
}

#[cfg(gpopt_adf="yes")]
pub fn get_terminals_for_func_def_branches() -> Vec<Vec<Terminal>> {
    let mut branches: Vec<Vec<Terminal>> = Vec::new();
    for _b_i in 0..NUM_ADF {
        // X, N, NE, E, SE, S, SW, W and NW
        let funcs = vec![
            Terminal::new(0, "X", terminal_x),
            Terminal::new(0, "N", terminal_n),
            Terminal::new(0, "NE", terminal_ne),
            Terminal::new(0, "E", terminal_e),
            Terminal::new(0, "SE", terminal_se),
            Terminal::new(0, "S", terminal_s),
            Terminal::new(0, "SW", terminal_sw),
            Terminal::new(0, "W", terminal_w),
            Terminal::new(0, "NW", terminal_nw),
        ];

        branches.push(funcs);
    }
    branches
}

// FitnessCase
const LETTER_WIDTH:  i8 = 4;
const LETTER_HEIGHT: i8 = 6;
pub type FitnessCase =
    [[bool; LETTER_WIDTH as usize]; LETTER_HEIGHT as usize];
pub trait FitnessCaseTrait {
    /// get pixel value for fixed point at `pos` tuple (x,y)
    fn get_absolute_pixel(&self, pos: (i8, i8)) -> bool;
    /// get pixel relatative to RunContext cur_pos
    fn get_relative_pixel(&self, rc: &RunContext, delta: (i8, i8)) -> bool;
    /// move cur_pos for RunContext to point relative to its current position
    /// get pixel at that new location.
    fn move_relative_pixel(&self, rc: &mut RunContext, delta: (i8, i8))
            -> bool;
}

impl FitnessCaseTrait for FitnessCase {
    /// get pixel value for fixed point at `pos` tuple (x,y)
    fn get_absolute_pixel(&self, pos: (i8, i8)) -> bool {
        self[pos.1 as usize][pos.0 as usize]
    }
    /// get pixel relatative to RunContext cur_pos
    fn get_relative_pixel(&self, rc: &RunContext, delta: (i8, i8)) -> bool {
        let pos =
            ((rc.cur_pos.0 + delta.0 + LETTER_WIDTH) % LETTER_WIDTH,
            (rc.cur_pos.1 + delta.1 + LETTER_HEIGHT) % LETTER_HEIGHT);

        self.get_absolute_pixel(pos)
    }
    /// move cur_pos for RunContext to point relative to its current position
    /// get pixel at that new location.
    fn move_relative_pixel(&self, rc: &mut RunContext, delta: (i8, i8))
            -> bool {
        rc.cur_pos.0 = (rc.cur_pos.0 + delta.0 + LETTER_WIDTH)
            % LETTER_WIDTH;
        rc.cur_pos.1 = (rc.cur_pos.1 + delta.1 + LETTER_HEIGHT) 
            % LETTER_HEIGHT;

        self.get_absolute_pixel(rc.cur_pos)
    }
}

pub struct FitnessCases {
    pub fc: [FitnessCase; 78],
}

const FITNESS_CASE_I_INDEX: usize = 0;
const FITNESS_CASE_L_INDEX: usize = 1;
pub const FITNESS_CASES: FitnessCases = FitnessCases {
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
            [ // 49 (Note: Koza II book shows copy of 45 here which seems like an error,
              //  esepecically because this obvous case used instead below is missing.)
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

/// IndividualRunResult returns results as type for RunContext::run_result
/// as a side affect not requiring for recursive return result.
pub enum IndividualRunResult {
    IsLetter(char),
    NotLetter,
}
use IndividualRunResult::*;

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
