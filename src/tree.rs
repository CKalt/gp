#[allow(non_snake_case)]
pub struct Control {
    pub M:   usize,             // Number of individuals in each generation
    pub G:   u16,               // Number of generations to run
    pub Di:  u16,               // Maximum depth of S Expressions for an initial tree
    pub Dc:  u16,               // Maximum depth of S Expressions for a created tree
    pub Pc:  f64,               // Probability of cross over
    pub Pr:  f64,               // Probability of reproduction
    pub Pip:  f64,              // Probability of cross over internal point
    pub no_functions:  u16,
    pub no_terminals:  u16,
    pub GRc:  f64,              // Greedy C Value (Function of M)
                            // n   M     C   GRc 
                            // 0  <1000  -    -
                            // 0   1000  1   .32
                            // 1   2000  2   .16
                            // 2   4000  4   .08
                            // 3   8000  8   .04
                            //                   GRc = 32/2^n
    pub no_fitness_cases:  u16
}

#[warn(non_snake_case)]
pub const CONTROL: Control = Control {
    M:                  1000,       // M - Number of individuals in each generation
    G:                  51,         // G - Number of generations to run
    Di:                 6,          // Di - Maximum depth of S Expressions for an initial tree
    Dc:                 17,         // Dc - Maximum depth of S Expressions for a created tree
    Pc:                 0.90,       // Pc - Probability of cross over
    Pr:                 0.10,       // Pr - Probability of reproduction
    Pip:                0.90,       // Pip - Probability of cross over internal point
    no_functions:       3,
    no_terminals:       3,
    GRc:                0.0,
    no_fitness_cases:   0,
};

#[allow(dead_code)]
enum Node {
    Terminal,
    Function
}

#[allow(dead_code)]
pub struct TreeSet {
    pub winning_index:  Option<usize>,
    avg_raw_f:          f64,
    pub tree_vec:       Vec<Tree>,
    tag_array:          [bool; CONTROL.M], // used for temporary individual state
                        // marking during tournament selection to insure
                        // unique individuals compete by avoiding duplicate
                        // usage.
}
impl TreeSet {
    pub fn new() -> TreeSet {
        TreeSet {
            winning_index:  None,
            avg_raw_f:      0.0,
            tree_vec:       Vec::new(),
            tag_array:      [false; CONTROL.M],
        }
    }
}

#[allow(dead_code)]
struct BaseFitness {
    nfr: f64,
    n:   f64,
    a:   f64,
    raw: f64,           // note that if loc Fitness might uses short type
                        // and then this just contains a casted copy
}
impl BaseFitness {
    fn new() -> BaseFitness {
        BaseFitness {
            nfr: 0.0,
            n:   0.0,
            a:   0.0,
            raw: 0.0,
        }
    }
}

#[allow(dead_code)]
pub struct Tree {
    tfid: u32,     // This is Tree's zero based index within TreeSet.tree_vec
                   // after sorting for fitness (least fit are lower valued)
    tcid: u32,     // The id of the Tree when first created and put into the array
    root: Node,
    fitness: BaseFitness,
    no_functions: u32,
    no_terminals: u32,
    hits: u32,
}
impl Tree {
    pub fn new() -> Tree {
        Tree { 
            tfid: 0,
            tcid: 0,
            root: Node::Terminal,
            fitness: BaseFitness::new(),
            no_functions: 0,
            no_terminals: 0,
            hits: 0,
        }
    }
}
