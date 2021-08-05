use mut_static::MutStatic;

use crate::fitness::GpFloat;
use crate::tree::*;

pub type TreeDepth = u16;

#[allow(non_snake_case)]
pub struct Control<'a> {
    #[allow(dead_code)]
    ftpair_rpb: FTSetPair<'a>,
    #[allow(dead_code)]
    opt_ftpair_fdb: Option<FTSetPair<'a>>,
    pub M:   usize,             // Number of individuals in each generation
    pub G:   u16,               // Number of generations to run
    pub Di:  TreeDepth,         // Maximum depth of S Expressions for an initial tree
    pub Dc:  TreeDepth,         // Maximum depth of S Expressions for a created tree
    pub Pc:  GpFloat,               // Probability of cross over
    pub Pr:  GpFloat,               // Probability of reproduction
    pub Pip:  GpFloat,              // Probability of cross over internal point
    pub num_functions_result_branch: u8,          // num functions result producing branch 0
    pub num_terminals_result_branch: u8,          // num terminals result producing branch 0
    pub num_functions_func_def_branch: u8,          // num functions function defining branch 0
    pub num_terminals_func_def_branch: u8,          // num terminals function defining branch 0
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
impl Control<'_> {
    pub fn new<'a>(ftpair_rpb: FTSetPair<'a>, 
               opt_ftpair_fdb: Option<FTSetPair<'a>>) -> Control<'a> {
        Control {
            ftpair_rpb,
            opt_ftpair_fdb,
            M:                  16000,      // Number of individuals in each generation
            G:                  51,         // Number of generations to run
            Di:                 6,          // Maximum depth of S Expressions for an initial tree
            Dc:                 17,         // Maximum depth of S Expressions for a created tree
            Pc:                 0.90,       // Probability of cross over
            Pr:                 0.10,       // Probability of reproduction
            Pip:                0.90,       // Probability of cross over internal point
            num_functions_result_branch: 5,          // num functions result producing branch 0
            num_terminals_result_branch: 6,          // num terminals result producing branch 0
            num_functions_func_def_branch: 4,          // num functions function defining branch 0
            num_terminals_func_def_branch: 2,          // num terminals function defining branch 0
            GRc:                0.16,
            R:                  0,
            W:                  1,
            no_fitness_cases:   0,
            show_all_trees:     false,
            show_all_tree_results: false,
            show_best_tree_results: true,
            show_controls: true,
            run_tests: false,
            run_log_file:       "gp_run.log",
        }
    }
    pub fn computational_effort(&self, runs: i32, gen: u16) -> i64 {
        self.M as i64 * (
            (self.G as i64 * (runs-1) as i64) + (gen+1) as i64
        )
    }
}

lazy_static! {
    pub static ref CONTROL: MutStatic<Control<'static>> = MutStatic::new();
}
