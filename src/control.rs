#[allow(non_snake_case)]
pub struct Control {
    pub M:   usize,             // Number of individuals in each generation
    pub G:   u16,               // Number of generations to run
    pub Di:  u16,               // Maximum depth of S Expressions for an initial tree
    pub Dc:  u16,               // Maximum depth of S Expressions for a created tree
    pub Pc:  f64,               // Probability of cross over
    pub Pr:  f64,               // Probability of reproduction
    pub Pip:  f64,              // Probability of cross over internal point
    pub num_functions:  u8,
    pub num_terminals:  u8,
    pub GRc:  f64,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            // GRc = 32/2^n
    pub no_fitness_cases:  u16,
    pub show_all_trees:    bool,
}

#[warn(non_snake_case)]
pub const CONTROL: Control = Control {
    M:                  1000,       // Number of individuals in each generation
    G:                  51,         // Number of generations to run
    Di:                 6,          // Maximum depth of S Expressions for an initial tree
    Dc:                 17,         // Maximum depth of S Expressions for a created tree
    Pc:                 0.90,       // Probability of cross over
    Pr:                 0.10,       // Probability of reproduction
    Pip:                0.90,       // Probability of cross over internal point
    num_functions:       3,
    num_terminals:       3,
    GRc:                0.0,
    no_fitness_cases:   0,
    show_all_trees:     false,
};
