use crate::tree::GpFloat;

#[cfg(gpopt_trail="santa_fe")]
const CONTROL_M: usize = 1000;            /// SANTE FE (M)
#[cfg(gpopt_trail="los_altos")]
const CONTROL_M: usize = 2000;            /// LOS ALTOS (M)

pub type TreeDepth = u16;

#[allow(non_snake_case)]
pub struct Control {
    pub M:   usize,             // Number of individuals in each generation
    pub G:   u16,               // Number of generations to run
    pub Di:  TreeDepth,         // Maximum depth of S Expressions for an initial tree
    pub Dc:  TreeDepth,         // Maximum depth of S Expressions for a created tree
    pub Pc:  GpFloat,               // Probability of cross over
    pub Pr:  GpFloat,               // Probability of reproduction
    pub Pip:  GpFloat,              // Probability of cross over internal point
    pub num_functions:  u8,
    pub num_terminals:  u8,
    pub GRc:  GpFloat,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            // GRc = 32/2^n
    pub R:  i32,            // stop after this many runs (0=no limit)
    pub no_fitness_cases:  u16,
    pub show_all_trees:    bool,
    pub show_all_tree_results: bool,
    pub show_best_tree_results: bool,
    pub show_controls: bool,
    pub run_tests: bool,
}

#[warn(non_snake_case)]
pub const CONTROL: Control = Control {
    M:                  CONTROL_M,  // Number of individuals in each generation
    G:                  51,         // Number of generations to run
    Di:                 6,          // Maximum depth of S Expressions for an initial tree
    Dc:                 17,         // Maximum depth of S Expressions for a created tree
    Pc:                 0.90,       // Probability of cross over
    Pr:                 0.10,       // Probability of reproduction
    Pip:                0.90,       // Probability of cross over internal point
    num_functions:      3,
    num_terminals:      3,
    GRc:                0.16,
    R:                  0,
    no_fitness_cases:   0,
    show_all_trees:     false,
    show_all_tree_results: false,
    show_best_tree_results: true,
    show_controls: true,
    run_tests: false,
};
