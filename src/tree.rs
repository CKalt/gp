use crate::gprun::*;
use crate::gprun::GridCellState::*;
use crate::control::CONTROL;
use crate::control::TreeDepth;
use crate::gprng::GpRng;
    
#[cfg(gpopt_trace="on")]
use crate::gprng::TRACE_COUNT;

#[cfg(not(gpopt_rng="FileStream"))]
use rand::Rng;

#[cfg(gpopt_choice_logging="write")]
use crate::choice_logging::*;
use Node::*;

pub enum GpType {
    Continue,
    Terminate,
}

type FuncNodeCode = fn (rc: &mut RunContext, fnc: &FunctionNode) -> GpType;
type TermNodeCode = fn (rc: &mut RunContext) -> GpType;

pub enum Node {
    TNode(& 'static Terminal), // Terminal nodes are borrowed references
    FNode(FunctionNode),  // Function nodes are not references, they are owners.
}
impl Node {
#[cfg(gpopt_choice_logging="write")]
    /// (Formerly called newRndFTNode() in gp.c)
    pub fn new_rnd(rng: &mut GpRng) -> Node {
        let num_ft = CONTROL.num_functions + CONTROL.num_terminals;
        let r = rng.gen_range(0i32..(num_ft as i32)) as u8;
        choice_log(2, &r.to_string());
        if r < CONTROL.num_terminals {
            TNode(&TERMINAL[r as usize])
        }
        else {
            let rand_fid = r - CONTROL.num_terminals;
            FNode(FunctionNode::new(rand_fid))
        }
    }
#[cfg(gpopt_choice_logging="off")]
    pub fn new_rnd(rng: &mut GpRng) -> Node {
        let num_ft = CONTROL.num_functions + CONTROL.num_terminals;
        let r = rng.gen_range(0..num_ft as i32) as u8;
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
    fn count_nodes(&self, counts: &mut (i32, i32)) {
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
    /// Print node recursively decending tree.
    fn print_r(&self, depth: TreeDepth) {
        match self {
            TNode(tn_ref) => print!(" {}", tn_ref.name),
            FNode(func_node) => {
                println!("");
                Node::print_tab(depth * 2);
                print!("({}", func_node.fnc.name);
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
                print!(")")
            }
        }
    }
    fn print_tab(count: u16) {
        for _ in 0..count {
            print!(" ");
        }
    }
    /// always performed against a Tree.root node
    /// which is why a tree is used instead of a Node here.
    /// It sets up the recursive call to find_function_node_ref_r.
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
    /// recursively iterate Node tree using depth first search to locate function
    /// node at index `fi`. At each entry `self` is a candidate function node.
    /// `*cur_fi` is incremented for each candidate and when equal to `fi`
    /// `self` is the found Some(result) returned up the recursive call stack,
    /// recursive iteration continues for None return values.
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
    /// always performed against a Tree.root node
    /// which is why a tree is used instead of a Node here.
    /// It sets up the recursive call to find_terminal_node_ref_r.
    fn find_terminal_node_ref(tree: &mut Tree, ti: TreeNodeIndex) -> &mut Node {
        match tree.root {
            TNode(_) => &mut tree.root,  // tree is just a node
            FNode(_) => {
                let mut cur_ti: TreeNodeIndex = 0;
                tree.root.find_terminal_node_ref_r(ti, &mut cur_ti).unwrap()
            }
        }
    }
    /// recursively iterate tree using depth first search to locate terminal
    /// node at index `ti`. Function children are recursiviely called, and 
    /// Terminal Nodes either return Some(self) if `ti` == `*cur_ti` or None
    /// causing iteration to continue.
    fn find_terminal_node_ref_r(&mut self, ti: TreeNodeIndex,
            cur_ti: &mut TreeNodeIndex) -> Option<&mut Node> {
        match self {
            TNode(_) => 
                if ti == *cur_ti {
                    Some(self) // found terminal we're looking for.
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
    fn depth_gt(&self, d: TreeDepth, so_far: TreeDepth) -> bool {
        if so_far > d {
            return true;
        }

        match self {
            TNode(_) => false,
            FNode(func_node) => {
                for node in func_node.branch.iter() {
                    if node.depth_gt(d, so_far + 1) {
                        return true;
                    }
                }
                false
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
        name: "if_food_ahead",
        arity: 2,
        code: if_food_ahead,
    },
    Function {
        fid:  1u8,
        name: "prog_n2",
        arity: 2,
        code: prog_n2
    },
    Function {
        fid:  2u8,
        name: "prog_n3",
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

    if let Food = rc.grid[new_y as usize][new_x as usize] {
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

pub struct Terminal {
    tid:    u8,
    name:   & 'static str,
    code:   TermNodeCode,
}

impl Terminal {
#[cfg(gpopt_choice_logging="write")]
    pub fn get_rnd_ref(rng: &mut GpRng) -> & 'static Terminal {
        let t_id: u8 = rng.gen_range(0i32..CONTROL.num_terminals as i32) as u8;
        choice_log(3, &t_id.to_string());
        &TERMINAL[t_id as usize]
    }
#[cfg(gpopt_choice_logging="off")]
    /// (Formerly called getRndTNode() in gp.c)
    pub fn get_rnd_ref(rng: &mut GpRng) -> & 'static Terminal {
        let t_id: u8 = rng.gen_range(0..CONTROL.num_terminals as i32) as u8;
        &TERMINAL[t_id as usize]
    }
}

static TERMINAL: [Terminal; CONTROL.num_terminals as usize] = [
    Terminal {
        tid:  0u8,
        name: "move",
        code: move_ant,
    },
    Terminal {
        tid:  1u8,
        name: "left",
        code: left,
    },
    Terminal {
        tid:  2u8,
        name: "right",
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
    else if let Food = rc.grid[rc.ant_y as usize][rc.ant_x as usize] {
        rc.grid[rc.ant_y as usize][rc.ant_x as usize] = FoodEaten;

        rc.eat_count += 1;
        if rc.eat_count == rc.n_pellets {
            return terminate(rc);
        }
    }
    // check for never been here before and no food
    else if let Clear = rc.grid[rc.ant_y as usize][rc.ant_x as usize] {
        rc.grid[rc.ant_y as usize][rc.ant_x as usize] = NoFoodFound;
    }
        
    // don't do anything but continue
    return GpType::Continue;
}

fn left(rc: &mut RunContext) -> GpType {
    if let GpType::Terminate = clock_tic(rc) {
        return GpType::Terminate
    }

    let t = rc.ant_xd;
    rc.ant_xd = rc.ant_yd;
    rc.ant_yd = -1 * t;

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
    /// (formerly called newRndFNode() in gp.c)
#[cfg(gpopt_choice_logging="write")]
    pub fn new_rnd(rng: &mut GpRng) -> FunctionNode {
        let rand_fid: u8 = rng.gen_range(0..CONTROL.num_functions as i32) as u8;
        choice_log(4, &rand_fid.to_string());

        FunctionNode::new(rand_fid)
    }
#[cfg(gpopt_choice_logging="off")]
    pub fn new_rnd(rng: &mut GpRng) -> FunctionNode {
        let rand_fid: u8 = rng.gen_range(0..CONTROL.num_functions as i32) as u8;
        FunctionNode::new(rand_fid)
    }

    /// index - 0 based arg index, node ownership is moved here
    /// to internal branch vector.
    /// will panic if attempt is made to assign index out of 
    /// 0..N sequence.
    pub fn set_arg(&mut self, index_u8: u8, node: Node) -> &mut Node {
        let index = index_u8 as usize;
        assert!(index < self.fnc.arity.into());

        if self.branch.len() > index {
            self.branch[index] = node;
        }
        else if self.branch.len() == index {
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

type TreeNodeIndex = i32;
#[allow(dead_code)]
pub struct TreeSet {
    pub winning_index:  Option<usize>,
    pub avg_raw_f:      GpFloat,
    pub tree_vec:       Vec<Tree>,
    tag_array:          [bool; CONTROL.M], // used for temporary individual state
                        // marking during tournament selection to insure
                        // unique individuals compete by avoiding duplicate
                        // usage.
    pub gen: u16,       // current generation for treeset
}
impl TreeSet {
    pub fn new(gen: u16) -> TreeSet {
        TreeSet {
            winning_index:  None,
            avg_raw_f:      0.0,
            tree_vec:       Vec::new(), // TODO: for better performance change to array (must find good init method)
            tag_array:      [false; CONTROL.M], // This is easy to init,
            gen:            gen
        }
    }

#[cfg(gpopt_choice_logging="write")]
    pub fn get_rnd_tree(&mut self, rng: &mut GpRng)
            -> &mut Tree {
        let rnd_index: usize = rng.gen_range(0..self.tree_vec.len() as i32) as usize;
        choice_log(5, &rnd_index.to_string());
        &mut self.tree_vec[rnd_index]
    }
#[cfg(gpopt_choice_logging="off")]
    pub fn get_rnd_tree(&mut self, rng: &mut GpRng)
            -> &mut Tree {
        let rnd_index: usize = rng.gen_range(0..self.tree_vec.len() as i32) as usize;
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
        let mut sum_a: GpInt = 0;
        let mut sum_raw: GpInt = 0;

        for t in self.tree_vec.iter() {
            sum_a += t.fitness.lng_a;
            sum_raw += t.fitness.lng_raw;
        }

        let int_avg_raw_f: GpInt = sum_raw / (self.tree_vec.len() as GpInt);
        self.avg_raw_f = Fitness::int_to_float(int_avg_raw_f);
        let flt_sum_a = Fitness::int_to_float(sum_a);
            
        #[cfg(gpopt_trace="on")]
        for (i,t) in self.tree_vec.iter_mut().enumerate() {
            let dbl_n = t.fitness.a() / flt_sum_a;
            t.fitness.lng_n = Fitness::float_to_int(dbl_n);

            {
                println!("TP1:i={}, tcid={}, hits={}, t.fitness.n={:10.9} a={:10.9} sum_a={:10.9}",
                    i, t.tcid, t.hits, t.fitness.n(), t.fitness.a(), flt_sum_a);

//                if t.tcid == 2 {
//                    t.print();
//                    exec_single_tree(t);
//                }
            }
        }
        #[cfg(gpopt_trace="off")]
        for t in self.tree_vec.iter_mut() {
            let dbl_n = t.fitness.a() / flt_sum_a;
            t.fitness.lng_n = Fitness::float_to_int(dbl_n);
        }

        #[allow(unreachable_code)]
        self
    }

    pub fn sort_by_normalized_fitness(&mut self) -> &mut TreeSet {
        self.tree_vec
            .sort_by(|a, b| a.fitness.lng_n.partial_cmp(&b.fitness.lng_n).unwrap());
            
        self.assign_nfr_rankings()
    }

    fn assign_nfr_rankings(&mut self) -> &mut TreeSet {
        let mut nfr: GpInt = 0;
        for (i, t) in self.tree_vec.iter_mut().enumerate() {
            nfr += t.fitness.lng_n;
            t.fitness.lng_nfr = nfr;
            t.tfid = Some(i);
        }
        self
    }
}

#[cfg(gpopt_select_method="fpb")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpInt, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpInt) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
}

#[cfg(gpopt_select_method="fpb")]
impl SelectMethod for TreeSet {
    fn select_ind_bin_r(&self, level: GpInt, lo: usize, hi: usize) -> usize {
        let gap = ((hi - lo) / 2) as usize;
        let guess = lo + gap + 1;

        if self.tree_vec[guess-1].fitness.lng_nfr < level &&
                self.tree_vec[guess].fitness.lng_nfr >= level {
            guess
        }
        else if self.tree_vec[guess].fitness.lng_nfr < level {
            // new_lo = guess
            self.select_ind_bin_r(level, guess, hi)
        }
        else {
            // new_hi = guess
            self.select_ind_bin_r(level, lo, guess-1)
        }
    }

    fn select_ind_bin(&self, level: GpInt) -> usize {
        assert_ne!(self.tree_vec.len(),0);
        let result = 
            if self.tree_vec[0].fitness.lng_nfr >= level ||
                    self.tree_vec.len() == 1 {
                0
            } else if level >
                    self.tree_vec[self.tree_vec.len() - 2].fitness.lng_nfr {
                self.tree_vec.len() - 1
            } else {
                self.select_ind_bin_r(level, 0, self.tree_vec.len() - 1)
            };
        result
    }

#[cfg(gpopt_choice_logging="write")]
    fn select_tree(&self, rng: &mut GpRng) -> &Tree {
        let greedy_val = rnd_greedy_val(rng);
        let i = self.select_ind_bin(greedy_val);
        choice_log(6, &i.to_string());
        return &self.tree_vec[i];
    }
#[cfg(gpopt_choice_logging="off")]
    fn select_tree(&self, rng: &mut GpRng) -> &Tree {
        let greedy_val = rnd_greedy_val(rng);
        let i = self.select_ind_bin(greedy_val);
        return &self.tree_vec[i];
    }
}

fn rnd_greedy_val(rng: &mut GpRng) -> GpInt {
    #[cfg(not(gpopt_rng="FileStream"))]
    let r = rng.gen::<GpFloat>();  // Note optional choice logging for this function
                               // done in TreeSet::select_tree, which is
                               // the only function that calls here. This allows
                               // integer logging thereby removing floating point variences.

    #[cfg(gpopt_rng="FileStream")]
    let r = rng.gen_float();

//    #[cfg(gpopt_trace="on")]
//    unsafe {
//        println!("TP005.2:(rnd_greedy_dbl) rlog={},r={}", TRACE_COUNT, r);
//    }

    let dbl_val: GpFloat = 
        if CONTROL.GRc < 0.0001 {
            r
        } else if rng.gen_range(0..10) > 1 {
            // select from top set
            1.0 - (CONTROL.GRc * r)
        } else {
            // select from lower set
            (1.0 - CONTROL.GRc) * r
        };
    Fitness::float_to_int(dbl_val)
}

pub type GpInt = i64;
pub type GpFloat = f64;
const DL_SHIFT: GpFloat = 1000000000.0;
pub struct Fitness {
    // Base values - real values are stored after multiplying by DL_SHIFT
    pub lng_nfr: GpInt,
    pub lng_n:   GpInt,
    pub lng_a:   GpInt,
    pub lng_raw: GpInt,   // note that if run Fitness uses integer type
                        // and then this just contains a converted copy

    // Ren Values
    pub r: u16,
    pub s: u16,
}
impl Fitness {
    fn new() -> Fitness {
        Fitness {
            lng_nfr: 0,
            lng_n:   0,
            lng_a:   0,
            lng_raw: 0,
            r: 0,
            s: 0,
        }
    }
    #[inline(always)]
    pub fn int_to_float(val: GpInt) -> GpFloat {
        (val as GpFloat) / DL_SHIFT
    }
    #[inline(always)]
    pub fn float_to_int(fval: GpFloat) -> GpInt {
        let fval: GpFloat = DL_SHIFT * fval;
        let r_fval = if fval > 0.0 { fval + 0.5 } else { fval - 0.5 };

        r_fval as GpInt
    }
    #[inline(always)]
    pub fn nfr(&self) -> GpFloat {
        Self::int_to_float(self.lng_nfr)
    }
    #[inline(always)]
    pub fn n(&self) -> GpFloat {
        Self::int_to_float(self.lng_n)
    }
    #[inline(always)]
    pub fn a(&self) -> GpFloat {
        Self::int_to_float(self.lng_a)
    }
}

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
    pub fn clear_node_counts(&mut self) {
        self.num_terminal_nodes = None;
        self.num_function_nodes = None;
    }
    /// count_nodes descends through tree and computes counts.
    /// To insure code performs efficiently we assert that counts are not
    /// counted twice. If this is infact needed
    /// I.e. Expects that nodes have not already been counted as represented by
    /// num_terminal_nodes and num_function_nodes values of None.
    pub fn count_nodes(&mut self) {
//        assert_eq!(self.num_terminal_nodes, None);
//        assert_eq!(self.num_function_nodes, None);
        let mut counts = (0i32, 0i32); // (num_terms, num_funcs)
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
        f.lng_n = -1;
        f.lng_a = -1;
        f.lng_nfr = -1;
        f.lng_raw = f.r as GpInt * DL_SHIFT as GpInt;

        // average over generation
        f.s = rc.n_pellets - rc.eat_count;
        let a = 1.0 / (1.0 + (f.s as GpFloat));
        f.lng_a = Fitness::float_to_int(a);

        return f.s == 0;
    }
            
    pub fn print(&self) {
        #[cfg(gpopt_trace="on")]
        {
            let count = get_log_count();

            println!("-log={:<10}------tree ({:?}/{})----------------",
                count, self.tfid, self.tcid);
        }
        #[cfg(gpopt_trace="off")]
        println!("-------------tree ({:?}/{})----------------",
            self.tfid, self.tcid);

        println!("nFunctions = {:?}\nnTerminals= {:?}", self.num_function_nodes,
            self.num_terminal_nodes);
        self.root.print();
        println!("");
    }

    pub fn get_rnd_function_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex, &mut Node) {
        let fi = self.get_rnd_function_index(rng);
        (fi, Node::find_function_node_ref(self, fi))
    }
    pub fn get_rnd_function_node_ref(&mut self,
            rng: &mut GpRng) -> &mut Node {
        let fi = self.get_rnd_function_index(rng);
        Node::find_function_node_ref(self, fi)
    }
#[cfg(gpopt_choice_logging="write")]
    fn get_rnd_function_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_fnodes = self
            .num_function_nodes
            .expect("Tree does not have count of function nodes. ");

        let result = rng.gen_range(0..num_fnodes as i32) as TreeNodeIndex;
        choice_log(7, &result.to_string());
        result
    }
#[cfg(gpopt_choice_logging="off")]
    fn get_rnd_function_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_fnodes = self
            .num_function_nodes
            .expect("Tree does not have count of function nodes. ");

        rng.gen_range(0..num_fnodes)
    }
    pub fn get_rnd_terminal_node_ref(&mut self,
            rng: &mut GpRng) -> &mut Node {
        let ti = self.get_rnd_terminal_index(rng);
        Node::find_terminal_node_ref(self, ti)
    }
#[cfg(gpopt_choice_logging="write")]
    fn get_rnd_terminal_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_tnodes = self
            .num_terminal_nodes
            .expect("Tree does not have count of terminal nodes. ");

        let result = rng.gen_range(0..num_tnodes) as TreeNodeIndex;
        choice_log(8, &result.to_string());
        result
    }
#[cfg(gpopt_choice_logging="off")]
    fn get_rnd_terminal_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_tnodes = self
            .num_terminal_nodes
            .expect("Tree does not have count of terminal nodes. ");

        rng.gen_range(0..num_tnodes)
    }
    fn depth_gt(&self, d: TreeDepth) -> bool {
        self.root.depth_gt(d, 1)
    }
    pub fn qualifies(&self) -> bool {
        !self.depth_gt(CONTROL.Dc)
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

