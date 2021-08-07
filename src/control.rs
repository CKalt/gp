use crate::fitness::GpFloat;
use lazy_static::*;
use crate::gprun::*;
use crate::tree::*;

pub type TreeDepth = u16;

#[allow(non_snake_case)]
pub struct Control<'a> {
    pub funcs_rpb:              Box<[&'a[Function]]>,
    pub terms_rpb:              Box<[&'a[Terminal]]>,
    pub funcs_fdb:              Box<[&'a[Function]]>,
    pub terms_fdb:              Box<[&'a[Terminal]]>,
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
impl Control<'_> {
    pub fn computational_effort(&self, runs: i32, gen: u16) -> i64 {
        self.M as i64 * (
            (self.G as i64 * (runs-1) as i64) + (gen+1) as i64
        )
    }
}

lazy_static! {
#[warn(non_snake_case)]
    pub static ref CONTROL: Control<'static> = {
        Control {
            #[cfg(gpopt_adf="yes")]
            funcs_rpb:          Box::new([&FUNCTIONS_RESULT_BRANCH_ADF]),
            #[cfg(gpopt_adf="no")]
            funcs_rpb:          Box::new([&FUNCTIONS_RESULT_BRANCH_NO_ADF]),
            terms_rpb:          Box::new([&TERMINALS_RESULT_BRANCH]),
            #[cfg(gpopt_adf="yes")]
            funcs_fdb:          Box::new([&FUNCTIONS_FUNC_DEF_BRANCH]),
            #[cfg(gpopt_adf="yes")]
            terms_fdb:          Box::new([&TERMINALS_FUNC_DEF_BRANCH]),
            #[cfg(gpopt_adf="no")]
            funcs_fdb:          Box::new([]),
            #[cfg(gpopt_adf="no")]
            terms_fdb:          Box::new([]),
            M:                  16000,      // Number of individuals in each generation
            G:                  51,         // Number of generations to run
            Di:                 6,          // Maximum depth of S Expressions for an initial tree
            Dc:                 17,         // Maximum depth of S Expressions for a created tree
            Pc:                 0.90,       // Probability of cross over
            Pr:                 0.10,       // Probability of reproduction
            Pip:                0.90,       // Probability of cross over internal point
            GRc:                0.16,
            R:                  0,
            W:                  1,
            no_fitness_cases:   0,
            show_all_trees:     false,
            show_all_tree_results: false,
            show_best_tree_results: true,
            show_controls: false,
            run_tests: false,
            run_log_file:       "gp_run.log",
        }
    };
}
