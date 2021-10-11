use crate::fitness::Fitness;
use crate::fitness::GpHits;
use crate::fitness::GpRaw;
use crate::fitness::GpFloat;
use crate::gprun::*;
use crate::control::CONTROL;
use crate::control::TreeDepth;
use crate::gprng::GpRng;
use rand::Rng;

#[cfg(test)]
use gp::parsers::parse_sexpr;
#[cfg(test)]
use gp::parsers::ParsedSNode;
#[cfg(test)]
use gp::parsers::ParsedSNode::*;

use Node::*;

type FuncNodeCode = fn (rc: &mut RunContext, fnc: &FunctionNode) -> GpType;
pub type TermNodeCode = fn (rc: &mut RunContext, &Terminal) -> GpType;

pub struct Winner {
    pub win: i32,   // Win number (among sequence of runs)
    pub tree: Tree,
    pub run: i32,
    pub gen: u16,
    pub e: u64,    // computational effort
}
impl Winner {
    pub fn print_result(&self) {
        Tree::print_result_header(None, &self.tree.fitness.hits);
        self.tree.print_result(None, -1.0);
    }
    pub fn structural_complexity(&self) -> TreeNodeIndex {
        self.tree.get_num_terminal_nodes().unwrap() 
            + self.tree.get_num_function_nodes().unwrap()
    }
}

pub enum Node {
    TNode(& 'static Terminal), // Terminal nodes are borrowed static references
    FNode(FunctionNode),  // FunctionNodes are not references, they are owners.
                          // They contain borrowed static references to
                          // Function structs
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
    #[cfg(gpopt_syntactic_constraints="yes")] 
    pub fn new_rnd_with_constraints(rng: &mut GpRng,
            funcs: &'static [Function], terms: &'static [Terminal],
            term_constraints: &Vec<u8>, func_constraints: &Vec<u8>) -> Node {
        let num_ft = funcs.len() + terms.len();
        loop {
            let r = rng.gen_range(0..num_ft as i32) as u8;
            if r < terms.len() as u8 {
                if term_constraints.iter().any(|&x| x == terms[r as usize].tid) {
                    return TNode(&terms[r as usize])
                }
            }
            else {
                let rand_fid = r - terms.len() as u8;
                if func_constraints.iter().any(|&x| x == rand_fid) {
                    return FNode(FunctionNode::new(rand_fid, funcs))
                }
            }
        }
    }
    // performs a deep (recusrive) build of node tree from parsed input.
    #[cfg(test)]
    pub fn deep_new_from_parse_tree(parsed_node: &ParsedSNode,
            funcs: &'static [Function], terms: &'static [Terminal]) -> Node {
        match parsed_node {
            Term(t_str) => {
                TNode(Terminal::get_ref_by_name(t_str, terms))
            },
            Func(f_str, f_args) => {
                let mut func_node = FunctionNode::new_by_name(f_str, funcs);
                for (i, arg) in f_args.iter().enumerate() {
                    if i > func_node.fnc.arity.into() {
                        panic!("More args parsed than function supports. fname={}",
                            func_node.fnc.name);
                    } else {
                        use std::convert::TryInto;
                        func_node.set_arg(i.try_into().unwrap(), 
                            Self::deep_new_from_parse_tree(arg, funcs, terms));
                    }
                }
                FNode(func_node)
            },
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
    fn count_nodes(&self, counts: &mut (usize, usize)) {
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
    pub fn print_node(&self) {
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
    #[cfg(gpopt_syntactic_constraints="yes")] 
    /// always performed against a TreeBranch.root node
    /// which is why a tree is used instead of a Node here.
    /// It sets up the recursive call to find_function_node_ref_ploc_r.
    /// This is an extended version of find_function_node_ref() that returns
    /// an optional ploc (parent locater value) that defines the found node's
    /// parent and child arg position that is required for syntactic
    /// constraints.
    fn find_function_node_ref_ploc(&mut self, fi: TreeNodeIndex) -> 
            (&mut Node, Option<(&Function, usize)>) {
        match self {
            TNode(_) => panic!(
              "Invalid attempt to find a Function node with a terminal tree."),
            FNode(_) => {
                let mut cur_fi: TreeNodeIndex = 0;
                // Since root is not a parent, passing None for parent 
                // locater (ploc) pair
                self.find_function_node_ref_ploc_r(fi, &mut cur_fi, None)
                    .unwrap()
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
    #[cfg(gpopt_syntactic_constraints="yes")] 
    /// recursively iterate Node tree using depth first search to locate function
    /// node at index `fi`. At each entry `self` is a candidate function node.
    /// `*cur_fi` is incremented for each candidate and when equal to `fi`
    /// `self` is the found Some(result) returned up the recursive call stack,
    /// recursive iteration continues for None return values.
    /// ploc: optional tuple consisting of : (parent node, zero based child
    /// arg position of node inside parent node branch vector)
    /// This is an extended version of find_function_node_ref_r() that returns
    /// an optional ploc (parent locater value) that defines the found node's
    /// parent and child arg position that is required for syntactic
    /// constraints.
    fn find_function_node_ref_ploc_r<'a>(&'a mut self, fi: TreeNodeIndex,
            cur_fi: &mut TreeNodeIndex,
            ploc: Option<(&'a Function, usize)>) -> 
                  Option<(&'a mut Node, Option<(&'a Function, usize)>)> {
        if let FNode(_) = self {
            if fi == *cur_fi {
                return Some((self, ploc));
            }
            else {
                // not found yet, so recursively descend into each
                // child function node of this function node (skiping
                // terminal nodes). Before each call (descent) the `self`
                // currencly values are set according to the location in tree.

                if let FNode(ref mut fn_ref) = self {
                    for (cn_i, cn) in fn_ref.branch.iter_mut().enumerate() {
                        if let FNode(_) = cn {
                            *cur_fi += 1;
                            let result = 
                                cn.find_function_node_ref_ploc_r(fi,
                                    cur_fi, Some((fn_ref.fnc, cn_i)));
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

// Types used in Function and Terminal Set constraint
// management.
//
// FSet and TSet define vector of vectors.
// Outside Vector indexed by branch number.
//   for Function defining branches this would be for
//   ADF0, ADF1...ADFN
//   for Result producint branches this would be for
//   RPB0, RPB1...RPBN
pub type FSet = Vec<Vec<Function>>; // [branch_i][func_i]
pub type TSet = Vec<Vec<Terminal>>; // [branch_i][term_i]
//pub type FTSet = (FSet, TSet);
#[cfg(gpopt_syntactic_constraints="yes")] 
type ArgNodeConstraints = Vec<NodeConstraints>; // constraints [argnum][ftids]
type NodeConstraints = Vec<u8>; // constraints [argnum][ftids]
#[cfg(gpopt_syntactic_constraints="yes")] 
pub type ArgNodeConsFTPair = (ArgNodeConstraints, ArgNodeConstraints);
pub type NodeConsFTPair = (NodeConstraints, NodeConstraints);

pub struct Function {
    #[allow(dead_code)]
    pub fid:    u8,
    pub name:   String,
    pub arity:  u8,
    pub code:   FuncNodeCode,
    pub opt_adf_num: Option<usize>, // if present, identifies adf number [0..n]
                                    // index into RunContext::opt_func_def_branches

    // following optional pair defines syntactic constraints that
    // when exist, define vector for each arg (i.e. 0..arity)
    // each holding a vector of the subset of fids and tids
    // allowed (included in the function or terminal set)
    // This filtering will be applied during tree_gen and crossover 
    // to insure only child args conforming to these contstrains exist.
    #[cfg(gpopt_syntactic_constraints="yes")] 
    pub opt_constraints: Option<ArgNodeConsFTPair>,
}

pub struct Terminal {
    pub tid:    u8,
    pub name:   String,
    pub code:   TermNodeCode,
    pub index:  usize,
}
impl Terminal {
    pub fn new(tid: u8, name: &str, code: TermNodeCode) -> Terminal {
        Terminal {
            tid,
            name: name.to_string(),
            code,
            index: 0,
        }
    }
    /// (Formerly called getRndTNode() in gp.c)
    pub fn get_rnd_ref(rng: &mut GpRng, terms: &'static [Terminal]) -> &'static Terminal {
        let t_id: u8 = rng.gen_range(0..terms.len() as i32) as u8;
        &terms[t_id as usize]
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    pub fn get_rnd_ref_with_constraints(rng: &mut GpRng,
            terms: &'static [Terminal], constraints: &Vec<u8>)
                -> &'static Terminal {
        let tc_id = rng.gen_range(0..constraints.len() as i32) as usize;
        let t_id = constraints[tc_id] as usize;
        &terms[t_id]
    }
    #[cfg(test)]
    pub fn get_ref_by_name(t_name: &str, terms: &'static [Terminal])
        -> &'static Terminal {
        for t in terms.iter() {
            if t.name.eq(t_name) {
                return t
            }
        }
        panic!("failed to get term ref by name={}", t_name);
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
    #[cfg(test)]
    fn new_by_name(f_name: &str, funcs: &'static [Function]) -> FunctionNode {
        for f in funcs.iter() {
            if f.name.eq(f_name) {
                return
                    FunctionNode {
                        fid:    f.fid,
                        fnc:    f,
                        branch: Vec::new()
                    }
            }
        }
        panic!("failed to get func by name={}", f_name);
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

    /// create a new FunctionNode choosing which Function at random.
    /// (formerly called newRndFNode() in gp.c)
    pub fn new_rnd(rng: &mut GpRng, funcs: &'static Vec<Function>) -> FunctionNode {
        assert_ne!(funcs.len(), 0);
        let rand_fid: u8 = rng.gen_range(0..funcs.len() as i32) as u8;
        FunctionNode::new(rand_fid, funcs)
    }

    /// create a new FunctionNode choosing which Function from a subset of 
    /// available Function's indicated by constraints vector.
    /// (formerly called newRndFNode() in gp.c)
    pub fn new_rnd_with_constraints(rng: &mut GpRng,
            funcs: &'static Vec<Function>,
            constraints: &NodeConstraints) -> FunctionNode {
        assert_ne!(funcs.len(), 0);
        let rand_cfid = rng.gen_range(0..constraints.len());
        let rand_fid = constraints[rand_cfid];
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

pub type TreeNodeIndex = usize;

pub enum BranchType {
    ResultProducing,
    FunctionDefining(TreeNodeIndex),
}
use BranchType::*;

pub struct TreeBranch {
    pub root: Node,
    pub num_function_nodes: Option<TreeNodeIndex>,
    pub num_terminal_nodes: Option<TreeNodeIndex>,
}
impl TreeBranch {
    pub fn new(root: Node) -> TreeBranch {
        TreeBranch{
            root,
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
        let mut counts = (0usize, 0usize);
        self.root.count_nodes(&mut counts);
        self.num_terminal_nodes = Some(counts.0);
        self.num_function_nodes = Some(counts.1);
    }
}

type TreeBranches = Vec<TreeBranch>;
trait CloneBranches {
    fn clone(&self) -> Self;
}

impl CloneBranches for TreeBranches {
    fn clone(&self) -> Self {
        let mut tree_branches: TreeBranches = Vec::new();
        for tree_branch in self {
            tree_branches.push(tree_branch.clone());
        }
        tree_branches
    }
}

pub struct Tree {
    pub tfid: Option<usize>,  // None until sorted, then this is Tree's zero
                      // based index within TreeSet.tree_vec after sorting for
                      // fitness (least fit are lower valued)
    pub tcid: usize,  // The id of the Tree when first created and put into the array
    pub fitness: Fitness,
    pub result_branch: TreeBranch,
    pub opt_func_def_branches: Option<Vec<TreeBranch>>, // If None no-adf, else adf.
    pub is_winner: bool,        // true if known winner
}
impl Tree {
    /// new Tree is constructed with a optional set of func def branches. 
    pub fn new(result_branch: TreeBranch,
            opt_func_def_branches: Option<Vec<TreeBranch>>) -> Tree {
        Tree { 
            tfid: None,
            tcid: 0,
            fitness: Fitness::new(),
            result_branch,
            opt_func_def_branches,
            is_winner: false,
        }
    }
    #[cfg(test)]
    /// mainly used for testing
    pub fn new_empty()  -> Tree {
        let terms = &CONTROL.terms_rpb[0];
        let node = Node::TNode(&terms[0]);
        let mut tree = Tree {
            tfid: None,
            tcid: 0,
            fitness: Fitness::new(),
            result_branch: TreeBranch::new(node),
            opt_func_def_branches: None,
            is_winner: true,
        };
        tree.count_nodes();
        tree
    }
    pub fn clone(&self) -> Tree {
        match &self.opt_func_def_branches {
            // no adf case
            None =>
                Tree { 
                    tfid: None,
                    tcid: 0,
                    fitness: self.fitness.clone(),
                    result_branch: self.result_branch.clone(),
                    opt_func_def_branches: None,
                    is_winner: self.is_winner,
                },
            // adf case
            Some(func_def_branches) =>
                Tree { 
                    tfid: None,
                    tcid: 0,
                    fitness: self.fitness.clone(),
                    result_branch: self.result_branch.clone(),
                    opt_func_def_branches: Some(func_def_branches.clone()),
                    is_winner: self.is_winner,
                },
        }
    }
    /// Get num function nodes for a specific branch type
    pub fn get_num_function_nodes_bt(&self, b_type: &BranchType)
            -> Option<TreeNodeIndex> {
        match b_type {
            ResultProducing => self.result_branch.num_function_nodes,
            FunctionDefining(adf_num) =>  
                self.opt_func_def_branches
                    .as_ref().unwrap()[*adf_num].num_function_nodes,
        }
    }
    /// Get total num function nodes across all branches
    pub fn get_num_function_nodes(&self) -> Option<TreeNodeIndex> {
        match &self.opt_func_def_branches {
            // no adf case (just pass through result branch node count)
            None => self.result_branch.num_function_nodes,
            // adf case - unwrap branches and reassemble optional (some) total
            Some(func_def_branches) =>
                if self.result_branch.num_function_nodes == None {
                    None // a pretty useless tree here.
                } else {
                    let mut num_fnodes_total: TreeNodeIndex =
                            self.result_branch.num_function_nodes.unwrap();
                    for func_def_branch in func_def_branches.iter() {
                        num_fnodes_total += func_def_branch.num_function_nodes.unwrap();
                    }
                    Some(num_fnodes_total)
                },
        }
    }
    /// Get num terminal nodes for a specific branch 
    pub fn get_num_terminal_nodes_bt(&self, b_type: &BranchType) -> Option<TreeNodeIndex> {
        match b_type {
            ResultProducing   => self.result_branch.num_terminal_nodes,
            FunctionDefining(adf_num) => self.opt_func_def_branches.as_ref()
                .unwrap()[*adf_num].num_terminal_nodes,
        }
    }
    /// Get total num terminal nodes across all branches
    pub fn get_num_terminal_nodes(&self) -> Option<TreeNodeIndex> {
        let num_rb_terminal_nodes = self.result_branch.num_terminal_nodes.unwrap();
        match &self.opt_func_def_branches {
            // no adf case - only need count of terminal nodes in result branch
            None => Some(num_rb_terminal_nodes),
            // adf case - need to add up result branch and func def branches
            Some(func_def_branches) => {
                let mut sum_tnodes = num_rb_terminal_nodes;
                for func_def_branch in func_def_branches.iter() {
                    sum_tnodes += func_def_branch.num_terminal_nodes.unwrap();
                }
                Some(sum_tnodes)
            },
        }
    }
    /// clear node counts across all branches
    pub fn clear_node_counts(&mut self) {
        self.result_branch.clear_node_counts(); // for adf and no-adf cases
        if let Some(func_def_branches) = &mut self.opt_func_def_branches {
            // adf case only
            for func_def_branch in func_def_branches.iter_mut() {
                func_def_branch.clear_node_counts();
            }
        }
    }
    /// count nodes across all branches
    pub fn count_nodes(&mut self) {
        self.result_branch.count_nodes(); // for adf and no-adf cases
        if let Some(func_def_branches) = &mut self.opt_func_def_branches {
            // adf case only
            for func_def_branch in func_def_branches.iter_mut() {
                func_def_branch.count_nodes();
            }
        }
    }
    /// print tree
    pub fn print_tree(&self) {
        #[cfg(gpopt_trace="on")]
        {
        //    let count = get_log_count();
            let count = 0;

            println!("-log={:<10}------tree ({:?}/{})----------------",
                count, self.tfid, self.tcid);
        }
        #[cfg(gpopt_trace="off")]
        println!("-------------tree ({:?}/{})----------------",
            self.tfid, self.tcid);

        println!("nFunctions = {:?}\nnTerminals= {:?}", self.get_num_function_nodes(),
            self.get_num_terminal_nodes());
        if let Some(func_def_branches) = &self.opt_func_def_branches {
            // adf case only
            for (adf_num, func_def_branch) in
                    func_def_branches.iter().enumerate() {
                println!("\n-------------------\n");
                println!("Function Def ADF{}:", adf_num);
                func_def_branch.root.print_node();
            }
        }
                
        println!("\n-------------------\n");
        println!("Result Branch0:");
        self.result_branch.root.print_node();
        println!("");
    }
    /// get a random function node over all branches. return node index, branch
    /// type and ref to node. Index sequence is in branch order major
    /// adf0..adfN, rb0, rb0..rbN and and minor within each branch tree, depth
    /// first.
    pub fn get_rnd_function_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex, BranchType, &mut Node) {
        let fi = self.get_rnd_function_index(rng);

        match &mut self.opt_func_def_branches {
            // no adf case
            None =>
                (fi, ResultProducing,
                    self.result_branch.root.find_function_node_ref(fi)),

            // adf case
            Some(func_def_branches) => {
                let mut sum_fnodes: TreeNodeIndex = 0;
                let mut fi_offset: TreeNodeIndex = 0;
                for (adf_num, func_def_branch) in
                    func_def_branches.iter_mut().enumerate() {
                    sum_fnodes += func_def_branch.num_function_nodes.unwrap();
                    if fi < sum_fnodes {
                        return (fi, FunctionDefining(adf_num), 
                            func_def_branch.root
                                .find_function_node_ref(fi - fi_offset))
                    }
                    fi_offset = sum_fnodes;
                }
                let adj_fi = fi - fi_offset;
                (adj_fi, ResultProducing, 
                    self.result_branch.root
                        .find_function_node_ref(adj_fi))
            }
        }
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    /// get a random function node over all branches. return branch type and
    /// mutable ref to node. (use get_rnd_function_node_ref if parent node 
    /// locater is also required.)
    pub fn get_rnd_function_node_ref(&mut self,
            rng: &mut GpRng) -> (BranchType, &mut Node) {
        let fi = self.get_rnd_function_index(rng);
        match &mut self.opt_func_def_branches {
            // no adf case
            None =>
                (ResultProducing,
                    self.result_branch.root.find_function_node_ref(fi)),
            // adf case
            Some(func_def_branches) => {
                let mut sum_fnodes: TreeNodeIndex = 0;
                let mut fi_offset: TreeNodeIndex = 0;
                for (adf_num, func_def_branch) in
                    func_def_branches.iter_mut().enumerate() {
                    sum_fnodes += func_def_branch.num_function_nodes.unwrap();
                    if fi < sum_fnodes {
                        return (FunctionDefining(adf_num), 
                            func_def_branch.root
                                .find_function_node_ref(fi - fi_offset))
                    }
                    fi_offset = sum_fnodes;
                }
                (ResultProducing, 
                    self.result_branch.root
                        .find_function_node_ref(fi - fi_offset))
            }
        }
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    /// get a random function node over all branches. return branch type and
    /// mutable ref to node and optional parent ref and arg_num
    /// (use get_rnd_function_node_ref if parent node locater is not needed.)
    /// This is an extended version of get_rnd_function_node_ref() that returns
    /// an optional ploc (parent locater value) that defines the found node's
    /// parent and child arg position that is required for syntactic
    /// constraints.
    pub fn get_rnd_function_node_ref_ploc(&mut self,
            rng: &mut GpRng) ->
                (BranchType, (&mut Node, Option<(&Function, usize)>)) {
        let fi = self.get_rnd_function_index(rng);
        match &mut self.opt_func_def_branches {
            // no adf case
            None =>
                (ResultProducing,
                    self.result_branch.root.find_function_node_ref_ploc(fi)),
            // adf case
            Some(func_def_branches) => {
                let mut sum_fnodes: TreeNodeIndex = 0;
                let mut fi_offset: TreeNodeIndex = 0;
                for (adf_num, func_def_branch) in
                    func_def_branches.iter_mut().enumerate() {
                    sum_fnodes += func_def_branch.num_function_nodes.unwrap();
                    if fi < sum_fnodes {
                        return (FunctionDefining(adf_num), 
                            func_def_branch.root
                                .find_function_node_ref_ploc(fi - fi_offset))
                    }
                    fi_offset = sum_fnodes;
                }
                (ResultProducing, 
                    self.result_branch.root
                        .find_function_node_ref_ploc(fi - fi_offset))
            }
        }
    }
    /// get a random function node by branch type (i.e. for a specific branch)
    /// return mutable ref to node.
    pub fn get_rnd_function_node_ref_bt(&mut self,
            rng: &mut GpRng, b_type: &BranchType) -> &mut Node {
        let fi = self.get_rnd_function_index_bt(rng, b_type);
        match b_type {
            ResultProducing => self
                .result_branch.root.find_function_node_ref(fi),
            FunctionDefining(adf_num) => self.opt_func_def_branches
                .as_mut().unwrap()[*adf_num].root.find_function_node_ref(fi),
        }
    }
    /// get a random function index for a specific branch type.
    fn get_rnd_function_index_bt(&self, rng: &mut GpRng, b_type: &BranchType)
            -> TreeNodeIndex {
        let num_fnodes = self
            .get_num_function_nodes_bt(b_type)
            .expect("Tree does not have count of function nodes.");
        assert_ne!(num_fnodes, 0);
        rng.gen_range(0..num_fnodes)
    }
    /// get a random function index across all nodes
    fn get_rnd_function_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_fnodes = self
            .get_num_function_nodes()
            .expect("Tree does not have count of function nodes.");

        rng.gen_range(0..num_fnodes)
    }
    /// get a random terminal node over all branches. return node index, branch
    /// type and ref to node. Index sequence is in branch order major
    /// adf0..adf1, rb0, rb0..rbN and and minor within each branch tree, depth
    /// first.
    pub fn get_rnd_terminal_node_ref(&mut self,
            rng: &mut GpRng) -> (BranchType, &mut Node) {
        let ti = self.get_rnd_terminal_index(rng);
        match &mut self.opt_func_def_branches {
            // no adf case
            None => (ResultProducing,
                self.result_branch.root.find_terminal_node_ref(ti)),
            // adf case
            Some(func_def_branches) => {
                let mut sum_tnodes: TreeNodeIndex = 0;
                let mut ti_offset: TreeNodeIndex = 0;
                for (adf_num, func_def_branch) in
                    func_def_branches.iter_mut().enumerate() {
                    sum_tnodes += func_def_branch.num_terminal_nodes.unwrap();

                    if ti < sum_tnodes {
                        return (FunctionDefining(adf_num), 
                            func_def_branch.root
                                .find_terminal_node_ref(ti - ti_offset))
                    }
                    ti_offset = sum_tnodes;
                }
                (ResultProducing, 
                    self.result_branch.root
                        .find_terminal_node_ref(ti - ti_offset))
            }
        }
    }
    /// get a random terminal node for specific branch. return ref to node.
    pub fn get_rnd_terminal_node_ref_bt(&mut self,
            rng: &mut GpRng, b_type: &BranchType) -> &mut Node {
        let ti = self.get_rnd_terminal_index_bt(rng, b_type);
        match b_type {
            ResultProducing => self.result_branch.root.find_terminal_node_ref(ti),
            FunctionDefining(adf_num) => self.opt_func_def_branches.
                as_mut().unwrap()[*adf_num].root.find_terminal_node_ref(ti),
        }
    }
    /// get a random terminal node over all branches. return node index, branch
    /// type and ref to node. Index sequence is in branch order major
    /// adfN..adf0, rb0, rb0..rbN and and minor within each branch tree, depth
    /// first.
    #[allow(dead_code)]
    pub fn get_rnd_terminal_node_ref_i(&mut self,
            rng: &mut GpRng)
        -> (TreeNodeIndex,  BranchType, &mut Node) {
        let ti = self.get_rnd_terminal_index(rng);
        match &mut self.opt_func_def_branches {
            // no adf case
            None => (ti, ResultProducing,
                self.result_branch.root.find_terminal_node_ref(ti)),
            // adf case
            Some(func_def_branches) => {
                let mut sum_tnodes: TreeNodeIndex = 0;
                let mut ti_offset: TreeNodeIndex = 0;
                for (adf_num, func_def_branch) in
                    func_def_branches.iter_mut().enumerate() {
                    sum_tnodes += func_def_branch.num_terminal_nodes.unwrap();

                    if ti < sum_tnodes {
                        let adj_ti = ti - ti_offset;
                        return (adj_ti, FunctionDefining(adf_num), 
                            func_def_branch.root
                                .find_terminal_node_ref(adj_ti))
                    }
                    ti_offset = sum_tnodes;
                }
                let adj_ti = ti - ti_offset;
                (adj_ti, ResultProducing, 
                    self.result_branch.root
                        .find_terminal_node_ref(adj_ti))
            }
        }
    }
    /// get a random terminal index for a specific branch type.
    fn get_rnd_terminal_index_bt(&self, rng: &mut GpRng, b_type: &BranchType)
            -> TreeNodeIndex {
        let num_tnodes = self
            .get_num_terminal_nodes_bt(b_type)
            .expect("Tree does not have count of terminal nodes.");

        rng.gen_range(0..num_tnodes)
    }
    /// get a random function index across all nodes
    fn get_rnd_terminal_index(&self, rng: &mut GpRng)
            -> TreeNodeIndex {
        let num_tnodes = self
            .get_num_terminal_nodes()
            .expect("Tree does not have count of terminal nodes. ");

        rng.gen_range(0..num_tnodes)
    }
    /// return true if depth of any branch in tree is greater than d.
    fn tree_depth_gt(&self, d: TreeDepth) -> bool {
        if self.result_branch.root.node_depth_gt(d, 1) {
            true
        } else {
            if let Some(func_def_branches) = &self.opt_func_def_branches {
                for func_def_branch in func_def_branches.iter() {
                    if func_def_branch.root.node_depth_gt(d, 1) {
                        return true;
                    }
                }
                return false;
            }
            else {
                false
            }
        }
    }
    /// return true tree qualifies in population based on local (self only)
    /// metrics
    pub fn qualifies(&self) -> bool {
        !self.tree_depth_gt(CONTROL.Dc)
    }

    /// Execute a single tree for each fitness case, compute fitness and print
    /// results. Return true when tree is a winner, otherwise false.
    pub fn eval_one_tree(&mut self) -> bool {
        self.print_tree();
        let mut rc = RunContext::new();
        rc.prepare_run();
        rc.print_run_illustration("Before Run");
            
        let mut sum_hits: GpHits = 0;
        let mut sum_error: GpRaw = 0;

        if let Some(func_def_branches) = &self.opt_func_def_branches {
            let mut func_def_branch_refs:  Vec<&TreeBranch> = Vec::new();
            for func_def_branch in func_def_branches {
                func_def_branch_refs.push(func_def_branch);
            }
            rc.opt_func_def_branches = Some(func_def_branch_refs);
        } else {
            rc.opt_func_def_branches = None;
        }

        for fc_i in 0..FITNESS_CASES.fc.len() {
            rc.init_new_fitness_case(fc_i);
            let result = Tree::exec_node(&mut rc, &self.result_branch.root);

            let error = rc.compute_error(result);
            let hit = error == 0;
            sum_error += error;
            if hit {
                sum_hits += 1;
            }
        }
        if let Some(_) = rc.opt_func_def_branches {
            rc.opt_func_def_branches = None;
        }

        rc.hits = sum_hits;
        rc.error = sum_error;

        let (f, is_winner) = rc.compute_fitness();
        self.fitness = f;
        if is_winner {
            println!("Have Winner");
        }
        rc.print_run_illustration("After Run");
        is_winner
    }
    /// print results for a single tree.
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
    /// print tree result header
    pub fn print_result_header(opt_gen: Option<u16>, hits: &GpHits) {
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
                (t.code)(rc, &t)
            }

            FNode(f) => {
                (f.fnc.code)(rc, &f)
            }
        }
    }
    #[inline(always)]
    pub fn exec_tree(&self, rc: &mut RunContext) -> GpType {
        Self::exec_node(rc, &self.result_branch.root)
    }

    #[cfg(test)]
    pub fn parse(
            rb0_s: &str,
            fdb_ses: Vec<&str>) -> Tree {
        // There always is a result branch and we parse that first.
        let rb_result = parse_sexpr(rb0_s);
        let rb_node =
            match rb_result {
                Ok(parsed_tuple) => parsed_tuple.1,
                Err(error) => panic!("cannot parse result branch: {:?}", error),
            };

        let result_branch_root = 
            Node::deep_new_from_parse_tree(&rb_node, &CONTROL.funcs_rpb[0],
                &CONTROL.terms_rpb[0]);

        // There may or may not be a func def branch0
        let opt_func_def_branches = 
            if CONTROL.terms_fdb.len() > 0 {
                let mut branches: Vec<TreeBranch> = Vec::new();
                for (b_i, fdb_s) in fdb_ses.iter().enumerate() {
                    let fd_result = parse_sexpr(fdb_s);
                    let fd_node = 
                        match fd_result {
                            Ok(parsed_tuple) => parsed_tuple.1,
                            Err(error) => panic!("cannot parse result branch: {:?}", error),
                        };
                    let tb = TreeBranch::new(
                                Node::deep_new_from_parse_tree(
                                    &fd_node,
                                    &CONTROL.funcs_fdb[b_i],
                                    &CONTROL.terms_fdb[b_i]));
                    branches.push(tb);
                }
                Some(branches)
            } else {
                None
            };

        let mut tree = Tree::new(
            TreeBranch::new(result_branch_root), opt_func_def_branches);
        tree.count_nodes();
        tree
    }
}
