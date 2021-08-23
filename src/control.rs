use std::sync::Mutex;
use crate::fitness::GpFloat;
use lazy_static::*;
use crate::gprun::*;
use crate::tree::*;

pub type TreeDepth = u16;

#[allow(non_snake_case)]
pub struct Control {
    pub funcs_rpb:              Box<Vec<Vec<Function>>>,
    pub terms_rpb:              Box<Vec<Vec<Terminal>>>,
    pub funcs_fdb:              Box<Vec<Vec<Function>>>,
    pub terms_fdb:              Box<Vec<Vec<Terminal>>>,
    pub M:   usize,             // Number of individuals in each generation
    pub G:   u16,               // Number of generations to run
    pub Di:  TreeDepth,         // Maximum depth of S Expressions for an initial tree
    pub Dc:  TreeDepth,         // Maximum depth of S Expressions for a created tree
    pub Pc:  GpFloat,               // Probability of cross over
    pub Pr:  GpFloat,               // Probability of reproduction
    pub Pip:  GpFloat,              // Probability of cross over internal point
    pub GRc:  GpFloat,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            // GRc = 32/2^n
    pub R:  i32,            // stop after this many runs (0=no limit) (inner loop)
    pub W:  i32,            // stop after this many wins (0=no limit) (outer loop)
    pub no_fitness_cases:  u16,
    pub show_all_trees:    bool,
    pub show_all_tree_results: bool,
    pub show_best_tree_results: bool,
    pub show_controls: bool,
    pub run_tests: bool,
    pub run_log_file:       &'static str,
}

pub struct EvalCount {
    winner_found: bool,
    num_evals_for_winner: u64,
}
impl EvalCount {
    fn new() -> Self {
        EvalCount {
            winner_found: false,
            num_evals_for_winner: 0,
        }
    }
    pub fn inc(is_winner: bool) -> bool {
        let mut ec = EVAL_COUNT.lock().unwrap();
        if !ec.winner_found {
            if is_winner {
                ec.winner_found = true;
            }
            ec.num_evals_for_winner += 1;
        }
        is_winner
    }
    pub fn num_evals_for_winner() -> u64 {
        let ec = EVAL_COUNT.lock().unwrap();
        ec.num_evals_for_winner
    }
    pub fn reset() {
        let mut ec = EVAL_COUNT.lock().unwrap();
        ec.num_evals_for_winner = 0;
        ec.winner_found = false;
    }
}

lazy_static! {
    pub static ref EVAL_COUNT: Mutex<EvalCount> = Mutex::new(EvalCount::new());
    pub static ref RUN_LOG_FNAME: String = run_log_fname();
    #[warn(non_snake_case)]
    pub static ref CONTROL: Control = {
        Control {
            funcs_rpb:          Box::new(get_functions_for_result_branches()),
            terms_rpb:          Box::new(get_terminals_for_result_branches()),
            funcs_fdb:          Box::new(get_functions_for_func_def_branches()),
            terms_fdb:          Box::new(get_terminals_for_func_def_branches()),
            M:                  4000,      // Number of individuals in each generation
            G:                  51,         // Number of generations to run
            Di:                 6,          // Maximum depth of S Expressions for an initial tree
            Dc:                 17,         // Maximum depth of S Expressions for a created tree
            Pc:                 0.90,       // Probability of cross over
            Pr:                 0.10,       // Probability of reproduction
            Pip:                0.90,       // Probability of cross over internal point
            GRc:                0.16,
            R:                  0,
            W:                  10,
            no_fitness_cases:   0,
            show_all_trees:     false,
            show_all_tree_results: false,
            show_best_tree_results: true,
            show_controls: false,
            run_tests: false,
            run_log_file:       &RUN_LOG_FNAME,
        }
    };
}

/// Called by lazy_static! macro to set RUN_LOG_FNAME
pub fn run_log_fname() -> String {
    #[cfg(gpopt_adf="yes")]
    let val = format!(
            "gp_run_even_{}_parity_{}_adf_{}_arity.log", EVEN_PARITY_K_VALUE,
            NUM_ADF, ADF_ARITY);
    #[cfg(gpopt_adf="no")]
    let val = format!(
            "gp_run_even_{}_parity_no_adf_{}_arity.log", 
            EVEN_PARITY_K_VALUE, ADF_ARITY);

    val
}
