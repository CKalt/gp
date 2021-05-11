use rand::Rng;
use crate::gprun::*;
use Node::*;

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
    no_fitness_cases:   0
};

pub enum GpType {
    Continue,
    Terminate,
}

type FuncNodeCode = fn (rc: &mut RunContext, fnc: &FunctionNode) -> GpType;
type TermNodeCode = fn (rc: &mut RunContext) -> GpType;

#[allow(dead_code)]
pub enum Node {
    TNode(& 'static Terminal), // Terminal nodes are borrowed references
    FNode(FunctionNode),  // Function nodes are not references, they are owners.
}
impl Node {
    pub fn new_rnd() -> Node {
        let mut rng = rand::thread_rng();
        let num_ft = CONTROL.num_functions + CONTROL.num_terminals;
        let r = rng.gen_range(0..num_ft);
        if r < CONTROL.num_terminals {
            TNode(&TERMINAL[r as usize])
        }
        else {
            let rand_fid = r - CONTROL.num_terminals;
            FNode(FunctionNode::new(rand_fid))
        }
    }
    fn clone(&self) -> Node {
        match self {
            TNode(tn_ref) => TNode(tn_ref),
            FNode(func_node) => FNode(func_node.clone())
        }
    }
    pub fn deep_match(&self, other: &Node) -> bool {
        match (self, other) {
            (TNode(tn1), TNode(tn2)) => tn1.tid == tn2.tid,
            (FNode(fn1), FNode(fn2)) => {
                if fn1.fid != fn2.fid {
                    false
                } 
                else {
                    for (sub1, sub2) in fn1.branch.iter().zip( 
                                        fn2.branch.iter()) {
                        if !sub1.deep_match(sub2) {
                            return false;
                        }
                    }
                    true
                }
            },
            _ => false,
        }
    }
    fn count_nodes(&self, counts: &mut (u32, u32)) {
        match self {
            TNode(_) => counts.0 += 1,
            FNode(ref fn_ref) => {
                counts.1 += 1;
                for node in fn_ref.branch.iter() {
                    node.count_nodes(counts);
                }
            }
        }
    }
}

#[allow(dead_code)]
pub struct Function {
    fid:    u8,
    name:   & 'static str,
    pub arity:  u8,
    code:   FuncNodeCode,
}

static FUNCTION: [Function; CONTROL.num_functions as usize] = [
    Function {
        fid:  0u8,
        name: "IfFoodAhead",
        arity: 2,
        code: if_food_ahead,
    },
    Function {
        fid:  1u8,
        name: "ProgN2",
        arity: 2,
        code: prog_n2
    },
    Function {
        fid:  2u8,
        name: "ProgN3",
        arity: 3,
        code: prog_n3
    },
];

#[allow(dead_code)]
fn clock_reset(rc: &mut RunContext) {
    rc.clock = 0;
}

fn clock_tic(rc: &mut RunContext) -> GpType {
    if rc.clock < RUN_CONTROL.max_clock {
        rc.clock += 1;
        GpType::Continue
    }
    else {
        GpType::Terminate
    }
}

fn terminate(rc: &mut RunContext) -> GpType {
    rc.clock = RUN_CONTROL.max_clock;
    GpType::Terminate
}

fn if_food_ahead(rc: &mut RunContext, func: &FunctionNode) -> GpType {
// if food exists at next step execute branch 0 else branch 1
// ignore any return values
    let new_x = rc.ant_x + rc.ant_xd;
    let new_y = rc.ant_y + rc.ant_yd;

    if new_x < 0 || new_x > RUN_CONTROL_MAX_X ||
       new_y < 0 || new_y > RUN_CONTROL_MAX_Y {
        return exec_node(rc, &func.branch[1]);
    }

    if rc.grid[new_y as usize][new_y as usize] == 1 {
        return exec_node(rc, &func.branch[0]);
    }

    exec_node(rc, &func.branch[1])
}

fn prog_n2(rc: &mut RunContext, func: &FunctionNode) -> GpType {
// unconditionally execute branch 0 & 1
// ignore any return values
    exec_node(rc, &func.branch[0]);
    exec_node(rc, &func.branch[1])
}

fn prog_n3(rc: &mut RunContext, func: &FunctionNode) -> GpType {
// unconditionally execute branch 0, 1 & 2
// ignore any return values
    exec_node(rc, &func.branch[0]);
    exec_node(rc, &func.branch[1]);
    exec_node(rc, &func.branch[2])
}

#[allow(dead_code)]
pub struct Terminal {
    tid:    u8,
    name:   & 'static str,
    code:   TermNodeCode,
}

impl Terminal {
    pub fn get_rnd_ref() -> & 'static Terminal {
        let mut rng = rand::thread_rng();
        let t_id: u8 = rng.gen_range(0..CONTROL.num_terminals);
        &TERMINAL[t_id as usize]
    }
}

static TERMINAL: [Terminal; CONTROL.num_terminals as usize] = [
    Terminal {
        tid:  0u8,
        name: "Move",
        code: move_ant,
    },
    Terminal {
        tid:  1u8,
        name: "Left",
        code: left,
    },
    Terminal {
        tid:  2u8,
        name: "Right",
        code: right,
    },
];

fn move_ant(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    rc.ant_x += rc.ant_xd;
    rc.ant_y += rc.ant_yd;

    // check for out of bounds
    if rc.ant_x < 0 || rc.ant_x > RUN_CONTROL_MAX_X ||
      rc.ant_y < 0 || rc.ant_y > RUN_CONTROL_MAX_Y {
        return GpType::Continue;
    }
    // check for food to eat
    else if rc.grid[rc.ant_y as usize][rc.ant_y as usize] == 1 {
        rc.grid[rc.ant_y as usize][rc.ant_y as usize] = 2; // mark as eaten

        rc.eat_count += 1;
        if rc.eat_count == rc.n_pellets {
            return terminate(rc);
        }
    }
    // check for never been here before
    else if rc.grid[rc.ant_y as usize][rc.ant_y as usize] == 0 {
        rc.grid[rc.ant_y as usize][rc.ant_y as usize] = 3;    // mark as been here
    }
        
    // don't do anything but continue
    return GpType::Continue;
}

fn left(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }
    let t = rc.ant_xd;
    rc.ant_xd = -1 * rc.ant_yd;
    rc.ant_yd = t;

    return GpType::Continue;
}

fn right(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }
    let t = rc.ant_xd;
    rc.ant_xd = -1 * rc.ant_yd;
    rc.ant_yd = t;

    return GpType::Continue;
}

pub struct FunctionNode {
    fid: u8,
    pub fnc: & 'static Function,
    branch: Vec<Node>, 
}
impl FunctionNode {
    fn new(fid: u8) -> FunctionNode {
        // fid maps to index into function array
        let fref = &FUNCTION[fid as usize];
        FunctionNode {
            fid:    fid,
            fnc:    fref,
            branch: Vec::new()
        }
    }

    fn clone(&self) -> FunctionNode {
        let mut new_func_node = FunctionNode::new(self.fid);
        for i in 0..self.fnc.arity {
            new_func_node.set_arg(i, self.branch[i as usize].clone());
        }

        new_func_node
    }

    /// create a new FunctionNode choosing which one at random.
    pub fn new_rnd() -> FunctionNode {
        let mut rng = rand::thread_rng();
        let rand_fid: u8 = rng.gen_range(0..CONTROL.num_functions);
        FunctionNode::new(rand_fid)
    }

    /// index - 0 based arg index, node ownership is moved here
    /// to internal branch vector.
    /// will panic if attempt is made to assign index out of 
    /// 0..N sequence.
    pub fn set_arg(&mut self, index: u8, node: Node) -> &mut Node {
        assert!(index < self.fnc.arity);

        if self.branch.len() > index.into() {
            self.branch[index as usize] = node;
        }
        else if self.branch.len() == index.into() {
            self.branch.push(node);
        }
        else {
            panic!("attempt to FunctionNode::set_arg out of sequence.");
        }

        &mut self.branch[index as usize]
    }
}

#[derive(Copy, Clone)]
pub enum GenerateMethod {
    Full,
    Grow
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
            tree_vec:       Vec::new(), // TODO: for better performance change to array (must find good init method)
            tag_array:      [false; CONTROL.M], // This is easy to init
        }
    }

    /// for each tree compute the n_functions, n_terminals
    /// values.
    pub fn count_nodes(&mut self) {
        for tree in self.tree_vec.iter_mut() {
            tree.count_nodes();
        }
    }

    pub fn compute_normalized_fitness(&mut self) {
        let mut sum_a = 0.0f64;
        let mut sum_raw = 0.0f64;

        for t in self.tree_vec.iter() {
            sum_a += t.fitness.a;
            sum_raw += t.fitness.raw;
        }

        self.avg_raw_f = sum_raw / (self.tree_vec.len() as f64);

        for t in self.tree_vec.iter_mut() {
            t.fitness.n = t.fitness.a / sum_a;
        }
    }
}

#[allow(dead_code)]
struct Fitness {
    // Base values
    nfr: f64,
    n:   f64,
    a:   f64,
    raw: f64,           // note that if run Fitness uses integer type
                        // and then this just contains a converted copy

    // Ren Values
    r: u16,
    s: u16,
}
impl Fitness {
    fn new() -> Fitness {
        Fitness {
            nfr: 0.0,
            n:   0.0,
            a:   0.0,
            raw: 0.0,
            r: 0,
            s: 0,
        }
    }
}

#[allow(dead_code)]
pub struct Tree {
    pub tfid: usize,  // This is Tree's zero based index within TreeSet.tree_vec
                      // after sorting for fitness (least fit are lower valued)
    tcid: usize,      // The id of the Tree when first created and put into the array
    pub root: Node,
    fitness: Fitness,
    num_functions: u32,
    num_terminals: u32,
    hits: u32,
}
impl Tree {
    pub fn new(root: FunctionNode) -> Tree {
        Tree { 
            tfid: 0,
            tcid: 0,
            root: FNode(root),
            fitness: Fitness::new(),
            num_functions: 0,
            num_terminals: 0,
            hits: 0,
        }
    }
    pub fn clone(&self) -> Tree {
        let new_root = self.root.clone(); // clone Node

        Tree { 
            tfid: 0,
            tcid: 0,
            root: new_root,
            fitness: Fitness::new(),
            num_functions: self.num_functions,
            num_terminals: self.num_terminals,
            hits: 0,
        }
    }
    pub fn count_nodes(&mut self) {
        let mut counts = (0u32, 0u32); // (num_terms, num_funcs)
        self.root.count_nodes(&mut counts);
        self.num_terminals = counts.0;
        self.num_functions = counts.1;
    }

    /// computes fitness returns true if winner
    pub fn compute_fitness(&mut self, rc: &RunContext) -> bool {
        let mut f = &mut self.fitness;
        f.r = rc.eat_count;
        self.hits = f.r as u32;

        // init fitness "base" values
        f.n = -1.0;
        f.a = -1.0;
        f.nfr = -1.0;
        f.raw = f.r.into();

        // average over generation
        f.s = rc.n_pellets - rc.eat_count;
        f.a = 1.0f64 / (1.0f64 + (f.s as f64));

        return f.s == 0;
   }
}

pub fn exec_node(rc: &mut RunContext, node: &Node) -> GpType {
    match node {
        TNode(t) => {
            (t.code)(rc)
        }

        FNode(f) => {
            (f.fnc.code)(rc, &f)
        }
    }
}

