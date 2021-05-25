use rand::Rng;
use crate::gprun::*;
use crate::control::CONTROL;
use crate::control::TreeDepth;
use Node::*;

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
    pub fn print(&self) {
        self.print_r(0)
    }
    fn print_r(&self, depth: TreeDepth) {
        match self {
            TNode(tn_ref) => print!(" {}", tn_ref.name),
            FNode(func_node) => {
                println!("");
                Node::print_tab(depth * 2);
                for node in &func_node.branch {
                    node.print_r(depth+1);
                }
                if depth == 0 {
                    // if there are one or more functions here
                    // we skip a line here for readabiility
                    for node in &func_node.branch {
                        if let FNode(_) = node {
                            println!("");
                            break;
                        }
                    }
                }
                println!(")")
            }
        }
    }
    fn print_tab(count: u16) {
        for _ in 0..count {
            print!(" ");
        }
    }
    fn find_function_node_ref(tree: &mut Tree, fi: TreeNodeIndex) -> &mut Node {
        match tree.root {
            TNode(_) => panic!(
              "Invalid attempt to find a Function node with a terminal tree."),
            FNode(_) => {
                let mut cur_fi: TreeNodeIndex = 0;
                tree.root.find_function_node_ref_r(fi, &mut cur_fi).unwrap()
            }
        }
    }
    /// recursively iterate tree using depth first search to locate function
    /// node at index `fi`. At each entry `self` is a candidate a function node.
    /// `*cur_fi` is incremented for each candidate when when equal to `fi`
    /// `self` is the found result returned up the recursive call stack.
    fn find_function_node_ref_r(&mut self, fi: TreeNodeIndex,
            cur_fi: &mut TreeNodeIndex) -> Option<&mut Node> {
        if let FNode(_) = self {
            if fi == *cur_fi {
                return Some(self);
            }
            else {
                // not found yet, so recursively descend into each
                // child function node of this function node (skiping
                // terminal nodes). Before each call (descent) the `self`
                // currencly values are set according to the location in tree.
                if let FNode(ref mut fn_ref) = self {
                    for cn in fn_ref.branch.iter_mut() {
                        if let FNode(_) = cn {
                            *cur_fi += 1;
                            let result = cn.find_function_node_ref_r(fi, cur_fi);
                            if let Some(_) = result {
                                return result;
                            }
                        }
                    }
                }
                else {
                    panic!("expected FNode.")
                }
            }
            return None;
        }
        else {
            panic!("expected function node");
        }
    }
    fn find_terminal_node_ref(tree: &mut Tree, ti: TreeNodeIndex) -> &mut Node {
        match tree.root {
            TNode(_) => panic!(
              "Invalid recursive attempt find a node for a single node root."),
            FNode(_) => {
                let mut cur_ti: TreeNodeIndex = 0;
                tree.root.find_terminal_node_ref_r(ti, &mut cur_ti).unwrap()
            }
        }
    }
    /// recursively iterate tree using depth first search to locate terminal
    /// node at index `ti`. At each entry `self.node` should be a function node.
    /// candidate at `self.ni`
    fn find_terminal_node_ref_r(&mut self, ti: TreeNodeIndex,
            cur_ti: &mut TreeNodeIndex) -> Option<&mut Node> {
        match self {
            TNode(_) => 
                if ti == *cur_ti {
                    Some(self)
                } else {
                    *cur_ti += 1;
                    None // not found yet
                },
            FNode(fn_ref) => {
                    // not found yet, so recursively descend into each
                    // child function node 
                    for cn in fn_ref.branch.iter_mut() {
                        let result = cn.find_terminal_node_ref_r(ti, cur_ti);
                        if let Some(_) = result {
                            return result;
                        }
                    }
                    None
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

type TreeNodeIndex = u32;
#[allow(dead_code)]
pub struct TreeSet {
    pub winning_index:  Option<usize>,
    pub avg_raw_f:          f64,
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

    pub fn get_rnd_tree(&mut self) -> &mut Tree {
        let mut rng = rand::thread_rng();
        let rnd_index: usize = rng.gen_range(0..self.tree_vec.len());
        &mut self.tree_vec[rnd_index]
    }

    pub fn print(&self) {
        for tree in self.tree_vec.iter() {
            tree.print();
        }
    }

    /// for each tree compute the n_functions, n_terminals
    /// values.
    pub fn count_nodes(&mut self) {
        for tree in self.tree_vec.iter_mut() {
            tree.count_nodes();
        }
    }

    pub fn compute_normalized_fitness(&mut self) -> &mut TreeSet {
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

        self
    }

    pub fn sort_by_normalized_fitness(&mut self) -> &mut TreeSet {
        self.tree_vec
            .sort_by(|a, b| a.fitness.n.partial_cmp(&b.fitness.n).unwrap());
            
        self.assign_nf_rankings()
    }

    fn assign_nf_rankings(&mut self) -> &mut TreeSet {
        let mut nfr = 0.0f64;
        for (i, t) in self.tree_vec.iter_mut().enumerate() {
            nfr += t.fitness.n;
            t.fitness.nfr = nfr;
            t.tfid = Some(i);
        }
        self
    }
}

#[allow(dead_code)]
pub struct Fitness {
    // Base values
    pub nfr: f64,
    pub n:   f64,
    pub a:   f64,
    pub raw: f64,           // note that if run Fitness uses integer type
                        // and then this just contains a converted copy

    // Ren Values
    pub r: u16,
    pub s: u16,
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
    pub tfid: Option<usize>,  // None until sorted, then this is Tree's zero
                      // based index within TreeSet.tree_vec after sorting for
                      // fitness (least fit are lower valued)
    pub tcid: usize,  // The id of the Tree when first created and put into the array
    pub root: Node,
    pub fitness: Fitness,
    pub num_function_nodes: Option<TreeNodeIndex>,
    pub num_terminal_nodes: Option<TreeNodeIndex>,
    pub hits: u32,
}
impl Tree {
    pub fn new(root: FunctionNode) -> Tree {
        Tree { 
            tfid: None,
            tcid: 0,
            root: FNode(root),
            fitness: Fitness::new(),
            num_function_nodes: None,     // None until count_nodes is called
            num_terminal_nodes: None,     // None until count_nodes is called
            hits: 0,
        }
    }
    pub fn clone(&self) -> Tree {
        let new_root = self.root.clone(); // clone Node

        Tree { 
            tfid: None,
            tcid: 0,
            root: new_root,
            fitness: Fitness::new(),
            num_function_nodes: self.num_function_nodes,
            num_terminal_nodes: self.num_terminal_nodes,
            hits: 0,
        }
    }
    /// count_nodes descends through tree and computes counts.
    /// To insure code performs efficiently we assert that counts are not
    /// counted twice. If this is infact needed
    /// I.e. Expects that nodes have not already been counted as represented by
    /// num_terminal_nodes and num_function_nodes values of None.
    pub fn count_nodes(&mut self) {
        assert_eq!(self.num_terminal_nodes, None);
        assert_eq!(self.num_function_nodes, None);
        let mut counts = (0u32, 0u32); // (num_terms, num_funcs)
        self.root.count_nodes(&mut counts);
        self.num_terminal_nodes = Some(counts.0);
        self.num_function_nodes = Some(counts.1);
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
    pub fn print(&self) {
        println!("-------------tree ({:?}/{})----------------",
            self.tfid, self.tcid);
        println!("nFunctions = {:?}\nnTerminals= {:?}", self.num_function_nodes,
            self.num_terminal_nodes);
        self.root.print();
        println!("");
    }
    pub fn get_rnd_function_node_ref_i(&mut self)
        -> (TreeNodeIndex, &mut Node) {
        let fi = self.get_rnd_function_index();
        (fi, Node::find_function_node_ref(self, fi))
    }
    pub fn get_rnd_function_node_ref(&mut self) -> &mut Node {
        let fi = self.get_rnd_function_index();
        Node::find_function_node_ref(self, fi)
    }
    fn get_rnd_function_index(&self) -> TreeNodeIndex {
        let num_fnodes = self
            .num_function_nodes
            .expect("Tree does not have count of function nodes. ");

        let mut rng = rand::thread_rng();
        rng.gen_range(0..num_fnodes)
    }
    pub fn get_rnd_terminal_node_ref(&mut self) -> &mut Node {
        let ti = self.get_rnd_terminal_index();
        Node::find_terminal_node_ref(self, ti)
    }
    fn get_rnd_terminal_index(&self) -> TreeNodeIndex {
        let num_tnodes = self
            .num_terminal_nodes
            .expect("Tree does not have count of terminal nodes. ");

        let mut rng = rand::thread_rng();
        rng.gen_range(0..num_tnodes)
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

