use format_num::NumberFormat;

#[cfg(gpopt_select_method="tournament")]
use std::collections::HashMap;

use crate::gprun::*;
use crate::control::CONTROL;
use crate::control::TreeDepth;
use crate::gprng::GpRng;
use rand::Rng;

use Node::*;

type FuncNodeCode = fn (rc: &mut RunContext, fnc: &FunctionNode) -> GpType;
pub type TermNodeCode = fn (rc: &RunContext) -> GpType;

pub struct Winner {
    pub tree: Tree,
    pub run: i32,
    pub gen: u16,
}
impl Winner {
    pub fn print_result(&self) {
        Tree::print_result_header(None, &self.tree.fitness.hits);
        self.tree.print_result(None, -1.0);
    }
}

pub enum Node {
    TNode(& 'static Terminal), // Terminal nodes are borrowed references
    FNode(FunctionNode),  // Function nodes are not references, they are owners.
}
impl Node {
    pub fn new_rnd(rng: &mut GpRng,
            funcs: &'static [Function], terms: &'static [Terminal]) -> Node {
        let num_ft = funcs.len() + terms.len();
        let r = rng.gen_range(0..num_ft as i32) as u8;
        if r < terms.len() as u8 {
            TNode(&terms[r as usize])
        }
        else {
            let rand_fid = r - terms.len() as u8;
            FNode(FunctionNode::new(rand_fid, funcs))
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
    /// always performed against a TreeBranch.root node
    /// which is why a tree is used instead of a Node here.
    /// It sets up the recursive call to find_function_node_ref_r.
    fn find_function_node_ref(&mut self, fi: TreeNodeIndex) -> &mut Node {
        match self {
            TNode(_) => panic!(
              "Invalid attempt to find a Function node with a terminal tree."),
            FNode(_) => {
                let mut cur_fi: TreeNodeIndex = 0;
                self.find_function_node_ref_r(fi, &mut cur_fi).unwrap()
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
    /// always performed against a TreeBranch::root node
    /// which is why a tree is used instead of a Node here.
    /// It sets up the recursive call to find_terminal_node_ref_r.
    fn find_terminal_node_ref(&mut self, ti: TreeNodeIndex) -> &mut Node {
        match self {
            TNode(_) => self,  // tree is just a node
            FNode(_) => {
                let mut cur_ti: TreeNodeIndex = 0;
                self.find_terminal_node_ref_r(ti, &mut cur_ti).unwrap()
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
    fn node_depth_gt(&self, d: TreeDepth, so_far: TreeDepth) -> bool {
        if so_far > d {
            return true;
        }

        match self {
            TNode(_) => false,
            FNode(func_node) => {
                for node in func_node.branch.iter() {
                    if node.node_depth_gt(d, so_far + 1) {
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
    pub fn get_rnd_ref(rng: &mut GpRng, terms: &'static [Terminal]) -> &'static Terminal {
        let t_id: u8 = rng.gen_range(0..terms.len() as i32) as u8;
        &terms[t_id as usize]
    }
}

pub struct FunctionNode {
    pub fid: u8,
    pub fnc: & 'static Function,
    pub branch: Vec<Node>, 
}
impl FunctionNode {
    fn new(fid: u8, funcs: &'static [Function]) -> FunctionNode {
        // fid maps to index into function array
        let fref = &funcs[fid as usize];
        FunctionNode {
            fid:    fid,
            fnc:    fref,
            branch: Vec::new()
        }
    }

    fn clone(&self) -> FunctionNode {
        let mut new_func_node = FunctionNode {
                fid: self.fid,
                fnc: self.fnc,
                branch: Vec::new(),
            };
        for i in 0..self.fnc.arity {
            new_func_node.set_arg(i, self.branch[i as usize].clone());
        }

        new_func_node
    }

    /// create a new FunctionNode choosing which one at random.
    /// (formerly called newRndFNode() in gp.c)
    pub fn new_rnd(rng: &mut GpRng, funcs: &'static [Function]) -> FunctionNode {
        let rand_fid: u8 = rng.gen_range(0..funcs.len() as i32) as u8;
        FunctionNode::new(rand_fid, funcs)
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

pub type GpHits = u16;
pub type GpRaw = f32;
pub type GpStandardized = f32;

pub type GpFloat = f32;

#[cfg(gpopt_fitness_type="int")]
pub type GpInt = i64;

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
                          // and then this just contains a converted copy of r

    pub hits: GpHits,
    pub r: GpRaw,
    pub s: GpStandardized,
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
            hits: 0,
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
    pub fn new() -> Fitness {
        Fitness {
            nfr: 0.0,
            n:   0.0,
            a:   0.0,
            raw: 0.0,

            r: 0.0,
            s: 0.0,
            hits: 0,
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
    fn clone(&self) -> Fitness {
        Fitness {
            nfr: self.nfr,
            n:   self.n,
            a:   self.a,
            raw: self.raw,

            r: self.r,
            s: self.s,
            hits: self.hits,
        }
    }
}

pub enum BranchType {
    Result0,
    FunctionDef0,
}
use BranchType::*;

pub struct TreeBranch {
    pub root: Node,
    pub num_function_nodes: Option<TreeNodeIndex>,
    pub num_terminal_nodes: Option<TreeNodeIndex>,
}
impl TreeBranch {
    fn new(root: FunctionNode) -> TreeBranch {
        TreeBranch{
            root: FNode(root),
            num_function_nodes: None,     // None until count_nodes is called
            num_terminal_nodes: None,     // None until count_nodes is called
        }
    }
    fn clone(&self) -> TreeBranch {
        TreeBranch{
            root: self.root.clone(),
            num_function_nodes: self.num_function_nodes,
            num_terminal_nodes: self.num_terminal_nodes,
        }
    }
    pub fn clear_node_counts(&mut self) {
        self.num_terminal_nodes = None;
        self.num_function_nodes = None;
    }
    pub fn count_nodes(&mut self) {
        let mut counts = (0i32, 0i32); // (num_terms, num_funcs)
        self.root.count_nodes(&mut counts);
        self.num_terminal_nodes = Some(counts.0);
        self.num_function_nodes = Some(counts.1);
    }
}

pub struct Tree {
    pub tfid: Option<usize>,  // None until sorted, then this is Tree's zero
                      // based index within TreeSet.tree_vec after sorting for
                      // fitness (least fit are lower valued)
    pub tcid: usize,  // The id of the Tree when first created and put into the array
    pub fitness: Fitness,
    pub result_branch: TreeBranch,
    pub func_def_branch: TreeBranch,
}
impl Tree {
    pub fn new(result_branch_root: FunctionNode, func_def_branch_root: FunctionNode) -> Tree {
        Tree { 
            tfid: None,
            tcid: 0,
            fitness: Fitness::new(),
            result_branch: TreeBranch::new(result_branch_root),
            func_def_branch: TreeBranch::new(func_def_branch_root),
        }
    }
    pub fn clone(&self) -> Tree {
        Tree { 
            tfid: None,
            tcid: 0,
            fitness: self.fitness.clone(),
            result_branch: self.result_branch.clone(),
            func_def_branch: self.func_def_branch.clone(),
        }
    }
    pub fn get_num_function_nodes_bt(&self, b_type: &BranchType) -> Option<TreeNodeIndex> {
        match b_type {
            Result0      =>  self.result_branch.num_function_nodes,
            FunctionDef0 =>  self.func_def_branch.num_function_nodes,
        }
    }
    pub fn get_num_function_nodes(&self) -> Option<TreeNodeIndex> {
        if self.func_def_branch.num_function_nodes == None &&
            self.result_branch.num_function_nodes == None {
            None
        } else {
            let num_fnodes1 = self.func_def_branch.num_function_nodes.unwrap();
            let num_fnodes2 = self.result_branch.num_function_nodes.unwrap();
            Some(num_fnodes1 + num_fnodes2)
        }
    }
    pub fn get_num_terminal_nodes_bt(&self, b_type: &BranchType) -> Option<TreeNodeIndex> {
        match b_type {
            Result0      =>  self.result_branch.num_terminal_nodes,
            FunctionDef0 =>  self.func_def_branch.num_terminal_nodes,
        }
    }
    pub fn get_num_terminal_nodes(&self) -> Option<TreeNodeIndex> {
        let num_tnodes1 = self.func_def_branch.num_terminal_nodes.unwrap();
        let num_tnodes2 = self.result_branch.num_terminal_nodes.unwrap();
        Some(num_tnodes1 + num_tnodes2)
    }
    pub fn clear_node_counts(&mut self) {
        self.result_branch.clear_node_counts();
        self.func_def_branch.clear_node_counts();
    }
    pub fn count_nodes(&mut self) {
        self.result_branch.count_nodes();
        self.func_def_branch.count_nodes();
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

        println!("nFunctions = {:?}\nnTerminals= {:?}", self.get_num_function_nodes(),
            self.get_num_terminal_nodes());
        println!("Function Def Branch0:");
        self.func_def_branch.root.print();
        println!("\nResult Branch0:");
        self.result_branch.root.print();
        println!("");
    }
    pub fn get_rnd_function_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex, BranchType, &mut Node) {
        let fi = self.get_rnd_function_index(rng);
        let num_fd_branch_fnodes = self.func_def_branch.num_function_nodes.unwrap();

        if fi < num_fd_branch_fnodes {
            (fi,FunctionDef0, 
                self.func_def_branch.root.find_function_node_ref(fi))
        } else {
            let adj_fi = fi - num_fd_branch_fnodes;
            (adj_fi, Result0,
                self.result_branch.root.find_function_node_ref(adj_fi))
        }
    }
    pub fn get_rnd_function_node_ref(&mut self,
            rng: &mut GpRng) -> (BranchType, &mut Node) {
        let fi = self.get_rnd_function_index(rng);
        let num_fd_branch_fnodes = self.func_def_branch.num_function_nodes.unwrap();

        if fi < num_fd_branch_fnodes {
            (FunctionDef0, 
                self.func_def_branch.root.find_function_node_ref(fi))
        } else {
            (Result0,
                self.result_branch.root.find_function_node_ref(
                        fi - num_fd_branch_fnodes))
        }
    }
    pub fn get_rnd_function_node_ref_bt(&mut self,
            rng: &mut GpRng, b_type: &BranchType) -> &mut Node {
        let fi = self.get_rnd_function_index_bt(rng, b_type);
        match b_type {
            Result0 => self.result_branch.root.find_function_node_ref(fi),
            FunctionDef0 => self.func_def_branch.root.find_function_node_ref(fi),
        }
    }
    fn get_rnd_function_index_bt(&self, rng: &mut GpRng, b_type: &BranchType)
            -> TreeNodeIndex {
        let num_fnodes = self
            .get_num_function_nodes_bt(b_type)
            .expect("Tree does not have count of function nodes.");

        rng.gen_range(0..num_fnodes)
    }
    fn get_rnd_function_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_fnodes = self
            .get_num_function_nodes()
            .expect("Tree does not have count of function nodes.");

        rng.gen_range(0..num_fnodes)
    }
    pub fn get_rnd_terminal_node_ref(&mut self,
            rng: &mut GpRng) -> (BranchType, &mut Node) {
        let ti = self.get_rnd_terminal_index(rng);
        let num_fd_branch_tnodes = self.func_def_branch.num_terminal_nodes.unwrap();

        if ti < num_fd_branch_tnodes {
            (FunctionDef0, 
                self.func_def_branch.root.find_terminal_node_ref(ti))
        } else {
            (Result0,
                self.result_branch.root.find_terminal_node_ref(
                        ti - num_fd_branch_tnodes))
        }
    }
    pub fn get_rnd_terminal_node_ref_bt(&mut self,
            rng: &mut GpRng, b_type: &BranchType) -> &mut Node {
        let ti = self.get_rnd_terminal_index_bt(rng, b_type);
        match b_type {
            Result0 => self.result_branch.root.find_terminal_node_ref(ti),
            FunctionDef0 => self.func_def_branch.root.find_terminal_node_ref(ti),
        }
    }
    #[allow(dead_code)]
    pub fn get_rnd_terminal_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex,  BranchType, &mut Node) {
        let ti = self.get_rnd_terminal_index(rng);
        let num_fd_branch_tnodes = self.func_def_branch.num_terminal_nodes.unwrap();

        if ti < num_fd_branch_tnodes {
            (ti, FunctionDef0, 
                self.func_def_branch.root.find_terminal_node_ref(ti))
        } else {
            let adj_ti = ti - num_fd_branch_tnodes;
            (adj_ti, Result0,
                self.result_branch.root.find_terminal_node_ref(adj_ti))
        }
    }
    fn get_rnd_terminal_index_bt(&self, rng: &mut GpRng, b_type: &BranchType)
            -> TreeNodeIndex {
        let num_tnodes = self
            .get_num_terminal_nodes_bt(b_type)
            .expect("Tree does not have count of terminal nodes.");

        rng.gen_range(0..num_tnodes)
    }
    fn get_rnd_terminal_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_tnodes = self
            .get_num_terminal_nodes()
            .expect("Tree does not have count of terminal nodes. ");

        rng.gen_range(0..num_tnodes)
    }
    fn tree_depth_gt(&self, d: TreeDepth) -> bool {
        self.result_branch.root.node_depth_gt(d, 1) ||
            self.func_def_branch.root.node_depth_gt(d, 1)
    }
    pub fn qualifies(&self) -> bool {
        !self.tree_depth_gt(CONTROL.Dc)
    }
    pub fn print_exec_one(&mut self) {
        self.print();
        let mut rc = RunContext::new();
        rc.prepare_run();
        rc.print_run_illustration("Before Run");
            
        #[cfg(gpopt_exec_criteria="clock")]
        panic!("clock exec not used this problem.");

        #[cfg(gpopt_exec_criteria="each_fitness_case")]
        let num = NumberFormat::new();
        #[cfg(gpopt_exec_criteria="each_fitness_case")]
        {
            let mut sum_hits: GpHits = 0;
            let mut sum_error: GpRaw = 0.0;

            rc.func_def_branch = Some(&self.func_def_branch);
            for i in 0..rc.fitness_cases.len() {
                let result = Tree::exec_node(&mut rc, &self.result_branch.root);

                let error = rc.compute_error(result);
                let hit = error < 0.01;
                sum_error += error;
                if hit {
                    sum_hits += 1;
                }

                println!("i={},d={},result={},error={},hit?={},sum_hits={}",
                    num.format("2d", i as f64),
                    num.format("10.5f", rc.get_cur_fc().d),
                    num.format("10.5f", result),
                    num.format("10.5f", error),
                    num.format("1d", hit as u8),
                    num.format("2d", sum_hits as f64));
            }
            rc.func_def_branch = None;

            rc.hits = sum_hits;
            rc.error = sum_error;
        }

        let (f, is_winner) = rc.compute_fitness();
        self.fitness = f;
        if is_winner {
            println!("Have Winner");
        }
        rc.print_run_illustration("After Run");
    }
    pub fn print_result(&self, opt_gen: Option<u16>, avg_raw_f: GpFloat) {
        let f = &self.fitness;
        let tfid = if let Some(num) = self.tfid { num } else { 0 };
        if let Some(gen) = opt_gen {
            println!("{:6} {:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                     gen, tfid, self.tcid, self.fitness.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
        } else {
            println!("{:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                    tfid, self.tcid, self.fitness.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
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
}

