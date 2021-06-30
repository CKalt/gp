use crate::gprun::*;
use crate::control::CONTROL;
use crate::control::TreeDepth;
use crate::gprng::GpRng;

#[cfg(not(gpopt_rng="FileStream"))]
use rand::Rng;

use Node::*;

type FuncNodeCode = fn (rc: &mut RunContext, fnc: &FunctionNode) -> GpType;
pub type TermNodeCode = fn (rc: &mut RunContext) -> GpType;

pub enum Node {
    TNode(& 'static Terminal), // Terminal nodes are borrowed references
    FNode(FunctionNode),  // Function nodes are not references, they are owners.
}
impl Node {
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

pub struct Function {
    #[allow(dead_code)]
    pub fid:    u8,
    pub name:   & 'static str,
    pub arity:  u8,
    pub code:   FuncNodeCode,
}

pub struct Terminal {
    pub tid:    u8,
    pub name:   & 'static str,
    pub code:   TermNodeCode,
}
impl Terminal {
    /// (Formerly called getRndTNode() in gp.c)
    pub fn get_rnd_ref(rng: &mut GpRng) -> & 'static Terminal {
        let t_id: u8 = rng.gen_range(0..CONTROL.num_terminals as i32) as u8;
        &TERMINAL[t_id as usize]
    }
}

pub struct FunctionNode {
    pub fid: u8,
    pub fnc: & 'static Function,
    pub branch: Vec<Node>, 
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

    #[cfg(gpopt_fitness_type="int")]
    pub fn compute_normalized_fitness(&mut self) -> &mut TreeSet {
        let mut sum_a: GpInt = 0;
        let mut sum_raw: GpInt = 0;

        for t in self.tree_vec.iter() {
            sum_a += t.fitness.a;
            sum_raw += t.fitness.raw;
        }

        let int_avg_raw_f: GpInt = sum_raw / (self.tree_vec.len() as GpInt);
        self.avg_raw_f = Fitness::int_to_float(int_avg_raw_f);
        let flt_sum_a = Fitness::int_to_float(sum_a);
            
        #[cfg(gpopt_trace="on")]
        for (i,t) in self.tree_vec.iter_mut().enumerate() {
            let dbl_n = t.fitness.a() / flt_sum_a;

            t.fitness.n = Fitness::float_to_int(dbl_n);
            println!("TP1:i={}, tcid={}, hits={}, t.fitness.n={:10.9} a={:10.9} sum_a={:10.9}",
                i, t.tcid, t.hits, t.fitness.n(), t.fitness.a(), flt_sum_a);

        }
        #[cfg(gpopt_trace="off")]
        for t in self.tree_vec.iter_mut() {
            let dbl_n = t.fitness.a() / flt_sum_a;
            t.fitness.n = Fitness::float_to_int(dbl_n);
        }

        self
    }
    #[cfg(gpopt_fitness_type="float")]
    pub fn compute_normalized_fitness(&mut self) -> &mut TreeSet {
        let mut sum_a: GpFloat = 0.0;
        let mut sum_raw: GpFloat = 0.0;

        for t in self.tree_vec.iter() {
            sum_a += t.fitness.a;
            sum_raw += t.fitness.raw;
        }

        self.avg_raw_f = sum_raw / (self.tree_vec.len() as GpFloat);

        #[cfg(gpopt_trace="on")]
        for (i,t) in self.tree_vec.iter_mut().enumerate() {
            t.fitness.n = t.fitness.a() / sum_a;
            println!("TP1:i={}, tcid={}, hits={}, t.fitness.n={:10.9} a={:10.9} sum_a={:10.9}",
                i, t.tcid, t.hits, t.fitness.n(), t.fitness.a(), sum_a);
        }
        #[cfg(gpopt_trace="off")]
        for t in self.tree_vec.iter_mut() {
            t.fitness.n = t.fitness.a() / sum_a;
        }

        self
    }

    pub fn sort_by_normalized_fitness(&mut self) -> &mut TreeSet {
        #[cfg(gpopt_fitness_type="int")]
        self.tree_vec
            .sort_by(|a, b| a.fitness.n.partial_cmp(&b.fitness.n).unwrap());
        #[cfg(gpopt_fitness_type="float")]
        self.tree_vec
            .sort_by(|a, b| a.fitness.n.partial_cmp(&b.fitness.n).unwrap());
            
        self.assign_nfr_rankings()
    }

    #[cfg(gpopt_fitness_type="int")]
    fn assign_nfr_rankings(&mut self) -> &mut TreeSet {
        let mut nfr: GpInt = 0;
        for (i, t) in self.tree_vec.iter_mut().enumerate() {
            nfr += t.fitness.n;
            t.fitness.nfr = nfr;
            t.tfid = Some(i);
        }
        self
    }
    #[cfg(gpopt_fitness_type="float")]
    fn assign_nfr_rankings(&mut self) -> &mut TreeSet {
        let mut nfr: GpFloat = 0.0;
        for (i, t) in self.tree_vec.iter_mut().enumerate() {
            nfr += t.fitness.n;
            t.fitness.nfr = nfr;
            t.tfid = Some(i);
        }
        self
    }

    pub fn exec(&mut self, run_number: i32) -> u16 {
        let mut rc = RunContext::new();

        for (i, tree) in self.tree_vec.iter_mut().enumerate() {
            rc.prepare_run();
            while rc.clock < RUN_CONTROL.max_clock {
                exec_node(&mut rc, &mut tree.root);
            }

            if rc.compute_fitness(tree) {
                tree.print_result(None, -1.0);
                rc.print_run_illustration(&format!("Have Winner! - Run# {} Gen# {}", run_number,
                    self.gen));
                self.winning_index = Some(i);
                break;
            }
        }
        rc.hits
    }
}

#[cfg(gpopt_select_method="fpb")]
#[cfg(gpopt_fitness_type="int")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpInt, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpInt) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
}
#[cfg(gpopt_select_method="fpb")]
#[cfg(gpopt_fitness_type="float")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpFloat, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpFloat) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
}


#[cfg(gpopt_select_method="fpb")]
impl SelectMethod for TreeSet {
    #[cfg(gpopt_fitness_type="int")]
    fn select_ind_bin_r(&self, level: GpInt, lo: usize, hi: usize) -> usize {
        let gap = ((hi - lo) / 2) as usize;
        let guess = lo + gap + 1;

        if self.tree_vec[guess-1].fitness.nfr < level &&
                self.tree_vec[guess].fitness.nfr >= level {
            guess
        }
        else if self.tree_vec[guess].fitness.nfr < level {
            // new_lo = guess
            self.select_ind_bin_r(level, guess, hi)
        }
        else {
            // new_hi = guess
            self.select_ind_bin_r(level, lo, guess-1)
        }
    }
    #[cfg(gpopt_fitness_type="float")]
    fn select_ind_bin_r(&self, level: GpFloat, lo: usize, hi: usize) -> usize {
        let gap = ((hi - lo) / 2) as usize;
        let guess = lo + gap + 1;

        if self.tree_vec[guess-1].fitness.nfr < level &&
                self.tree_vec[guess].fitness.nfr >= level {
            guess
        }
        else if self.tree_vec[guess].fitness.nfr < level {
            // new_lo = guess
            self.select_ind_bin_r(level, guess, hi)
        }
        else {
            // new_hi = guess
            self.select_ind_bin_r(level, lo, guess-1)
        }
    }

    #[cfg(gpopt_fitness_type="int")]
    fn select_ind_bin(&self, level: GpInt) -> usize {
        assert_ne!(self.tree_vec.len(),0);
        let result = 
            if self.tree_vec[0].fitness.nfr >= level ||
                    self.tree_vec.len() == 1 {
                0
            } else if level >
                    self.tree_vec[self.tree_vec.len() - 2].fitness.nfr {
                self.tree_vec.len() - 1
            } else {
                self.select_ind_bin_r(level, 0, self.tree_vec.len() - 1)
            };
        result
    }
    #[cfg(gpopt_fitness_type="float")]
    fn select_ind_bin(&self, level: GpFloat) -> usize {
        assert_ne!(self.tree_vec.len(),0);
        let result = 
            if self.tree_vec[0].fitness.nfr >= level ||
                    self.tree_vec.len() == 1 {
                0
            } else if level >
                    self.tree_vec[self.tree_vec.len() - 2].fitness.nfr {
                self.tree_vec.len() - 1
            } else {
                self.select_ind_bin_r(level, 0, self.tree_vec.len() - 1)
            };
        result
    }

    fn select_tree(&self, rng: &mut GpRng) -> &Tree {
        let greedy_val = rnd_greedy_val(rng);
        let i = self.select_ind_bin(greedy_val);
        return &self.tree_vec[i];
    }
}

#[cfg(gpopt_fitness_type="int")]
fn rnd_greedy_val(rng: &mut GpRng) -> GpInt {
    #[cfg(not(gpopt_rng="FileStream"))]
    let r = rng.gen::<GpFloat>();  // Note optional choice logging for this function
                               // done in TreeSet::select_tree, which is
                               // the only function that calls here. This allows
                               // integer logging thereby removing floating point variences.

    #[cfg(gpopt_rng="FileStream")]
    let r = rng.gen_float();

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
#[cfg(gpopt_fitness_type="float")]
fn rnd_greedy_val(rng: &mut GpRng) -> GpFloat {
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
    dbl_val
}


#[cfg(gpopt_fitness_type="int")]
pub type GpInt = i64;
pub type GpFloat = f64;

#[cfg(gpopt_fitness_type="int")]
pub type GpFitness = GpInt;

#[cfg(gpopt_fitness_type="float")]
pub type GpFitness = GpFloat;

#[cfg(gpopt_fitness_type="int")]
pub const DL_SHIFT: GpFloat = 1000000000.0;
        
pub struct Fitness {
    // Base values - real values are stored after multiplying by DL_SHIFT
    pub nfr: GpFitness,
    pub n:   GpFitness,
    pub a:   GpFitness,
    pub raw: GpFitness,   // note that if run Fitness uses integer type
                        // and then this just contains a converted copy

    // Ren Values
    pub r: u16,
    pub s: u16,
}

#[cfg(gpopt_fitness_type="int")]
impl Fitness {
    fn new() -> Fitness {
        Fitness {
            nfr: 0,
            n:   0,
            a:   0,
            raw: 0,

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
        Self::int_to_float(self.nfr)
    }
    #[inline(always)]
    pub fn n(&self) -> GpFloat {
        Self::int_to_float(self.n)
    }
    #[inline(always)]
    pub fn a(&self) -> GpFloat {
        Self::int_to_float(self.a)
    }
}

#[cfg(gpopt_fitness_type="float")]
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
    #[inline(always)]
    pub fn nfr(&self) -> GpFloat {
        self.nfr
    }
    #[inline(always)]
    pub fn n(&self) -> GpFloat {
        self.n
    }
    #[inline(always)]
    pub fn a(&self) -> GpFloat {
        self.a
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
    #[allow(dead_code)]
    pub fn get_rnd_function_node_ref(&mut self,
            rng: &mut GpRng) -> &mut Node {
        let fi = self.get_rnd_function_index(rng);
        Node::find_function_node_ref(self, fi)
    }
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
    #[allow(dead_code)]
    pub fn get_rnd_terminal_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex, &mut Node) {
        let ti = self.get_rnd_terminal_index(rng);
        (ti, Node::find_terminal_node_ref(self, ti))
    }
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
    pub fn exec(&mut self) {
        let mut rc = RunContext::new();
        rc.prepare_run();
        rc.print_run_illustration("Before Run");
        while rc.clock < RUN_CONTROL.max_clock {
            exec_node(&mut rc, &mut self.root);
        }

        if rc.compute_fitness(self) {
            println!("Have Winner");
        }
        rc.print_run_illustration("After Run");
    }
    pub fn print_result(&self, opt_gen: Option<u16>, avg_raw_f: GpFloat) {
        let f = &self.fitness;
        let tfid = if let Some(num) = self.tfid { num } else { 0 };
        if let Some(gen) = opt_gen {
            println!("{:6} {:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                     gen, tfid, self.tcid, self.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
        } else {
            println!("{:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                    tfid, self.tcid, self.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
        }
    }
    pub fn print_result_header(opt_gen: Option<u16>, hits: &u16) {
        println!("{}={}", RunContext::get_hits_label(), *hits);
        if let Some(_) = opt_gen {
            println!("{:>6} {:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8} {:>6}", 
                "gen", "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr", "avgRawF");
            println!("{:>6} {:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
                "----", "---", "---", "-----", "---", "---", "------", "------", "------");
        } else {
            println!("{:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
                "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr");
            println!("{:>4} {:>4} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8}", 
                "---", "---", "-----", "---", "---", "------", "------", "------");
        }
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

