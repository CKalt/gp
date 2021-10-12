use std::mem;
use rand::Rng;
    
#[cfg(gpopt_mult_threaded="yes")]
use rayon::prelude::*;

#[cfg(gpopt_select_method="tournament")]
use std::collections::HashMap;

use crate::gprun::*;
use crate::tree::*;
use crate::control::*;
use crate::gprng::GpRng;
use crate::fitness::GpFloat;
use crate::fitness::GpHits;
use crate::fitness::GpRaw;

use Node::*;

pub struct TreeSet {
    pub winning_index:  Option<usize>,
    pub avg_raw_f:      GpFloat,
    pub tree_vec:       Vec<Tree>,
    pub gen: u16,       // current generation for treeset
}

impl TreeSet {
    pub fn new(gen: u16) -> TreeSet {
        TreeSet {
            winning_index:  None,
            avg_raw_f:      0.0,
            tree_vec:       Vec::new(),
            gen:            gen
        }
    }
    fn tree_match(t1: &Tree, t2: &Tree) -> bool {
        if t1.result_branch.root.deep_match(&t2.result_branch.root) {
            match &t1.opt_func_def_branches {
                // no adf case
                None => true,
                // adf case
                Some(func_def_branches) => {
                    // if any fd branch doesn't match false
                    for (b_i, func_def_branch) in
                        func_def_branches.iter().enumerate() {
                        if !func_def_branch.root.deep_match(
                            &t2.opt_func_def_branches.as_ref().unwrap()[b_i].root
                        ) {
                            return false;
                        }
                    }
                    true
                },
            }
        } else {
            false
        }
    }
    fn tree_match_exists(&self, target_tree: &Tree) -> bool {
        for tree in self.tree_vec.iter() {
            if Self::tree_match(tree, target_tree) {
                return true;
            }
        }
        return false;
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    fn gen_tree_full_method(rng: &mut GpRng, depth: u16) -> Tree {
        assert_eq!(CONTROL.terms_fdb.len(), CONTROL.funcs_fdb.len());

        let mut result_branch_root =
            FunctionNode::new_rnd(rng, &CONTROL.funcs_rpb[0]);

        // call version without constraints arg
        Self::gen_tree_full_method_r(rng, &mut result_branch_root, 2, depth,
                &CONTROL.funcs_rpb[0], &CONTROL.terms_rpb[0]);

        if CONTROL.terms_fdb.len() > 0 {
            let mut func_def_branches: Vec<TreeBranch> = Vec::new();
            for i in 0..CONTROL.funcs_fdb.len() {
                let mut func_def_branch_root =
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_fdb[i]);

                Self::gen_tree_full_method_r(rng, &mut func_def_branch_root,
                    2, depth, &CONTROL.funcs_fdb[i], &CONTROL.terms_fdb[i]);
                func_def_branches.push(TreeBranch::new(FNode(func_def_branch_root)));
            }
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), 
                Some(func_def_branches))
        } else {
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), None)
        }
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    fn gen_tree_full_method(rng: &mut GpRng, depth: u16) -> Tree {
        assert_eq!(CONTROL.terms_fdb.len(), CONTROL.funcs_fdb.len());

        let mut result_branch_root = {
                let root_constraints: Box<Vec<NodeConsFTPair>>;
                if let Some(root_constraints) = &CONTROL.opt_rpb_root_constraints {
                    let (f_set, _) = root_constraints[0];
                    FunctionNode::new_rnd_constrained(
                        rng, &CONTROL.funcs_rpb[0], &f_set)
                } else {
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_rpb[0])
                }
            };

        // call version with constraints arg set to None
        Self::gen_tree_full_method_r(rng, &mut result_branch_root, 2, depth,
                &CONTROL.funcs_rpb[0], &CONTROL.terms_rpb[0], &None);

        if CONTROL.terms_fdb.len() > 0 {
            let mut func_def_branches: Vec<TreeBranch> = Vec::new();
            for i in 0..CONTROL.funcs_fdb.len() {
                let mut func_def_branch_root =
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_fdb[i]);

                Self::gen_tree_full_method_r(rng, &mut func_def_branch_root,
                    2, depth, &CONTROL.funcs_fdb[i], &CONTROL.terms_fdb[i], &None);
                func_def_branches.push(TreeBranch::new(FNode(func_def_branch_root)));
            }
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), 
                Some(func_def_branches))
        } else {
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), None)
        }
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    fn gen_tree_full_method_r(rng: &mut GpRng,
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static Vec<Function>, terms: &'static Vec<Terminal>) {
        if level >= depth {
            for i in 0..func_node.fnc.arity {
                // Always a Terminal Node
                let rnd_tref = Terminal::get_rnd_ref(rng, terms);
                func_node.set_arg(i, TNode(rnd_tref));
            }
        } else {
            let c_depth = level+1;
            for i in 0..func_node.fnc.arity { 
                // Always a Funciton Node
                let rnd_fn = FunctionNode::new_rnd(rng, funcs);
                let node: &mut Node = func_node.set_arg(i, FNode(rnd_fn));

                match *node {
                    FNode(ref mut fn_ref) =>
                        Self::gen_tree_full_method_r(rng, fn_ref, c_depth,
                                depth, funcs, terms),
                    _ => panic!("expected FunctionNode"),
                }
            }
        }
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    // opt_prev_constraints is used here to implement a deep application
    // of a constraint. Initial state is None, and subsequent states are
    // either a continuance of the previous state, or a change to a new one.
    fn gen_tree_full_method_r(rng: &mut GpRng,
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static Vec<Function>, terms: &'static Vec<Terminal>,
            opt_prev_constraints: &Option<ArgNodeConsFTPair>) {

        // If the current func_node has constraints, use them,
        // otherwise fallback to previously used constraints if any.
        let opt_constraints =
            if let Some(_) = func_node.fnc.opt_constraints {
                &func_node.fnc.opt_constraints
            } else if let Some(_) = *opt_prev_constraints {
                opt_prev_constraints
            }
            else {
                &None
            };
            
        if level >= depth {
            // Always a Terminal Node
            if let Some((_, ref term_constraints)) = opt_constraints {
                for i in 0..func_node.fnc.arity {
                    // We have syntactic terminal constraints for this function
                    // so loop until the randomly selected terminal qualifies.
                    let rnd_tref =
                        Terminal::get_rnd_ref_with_constraints(rng, terms,
                            &term_constraints[i as usize]);
                    func_node.set_arg(i, TNode(rnd_tref));
                }
            } else {
                // No terminal constraints so anything in terminal set (terms)
                // is allowed.
                for i in 0..func_node.fnc.arity {
                    let rnd_tref = Terminal::get_rnd_ref(rng, terms);
                    func_node.set_arg(i, TNode(rnd_tref));
                }
            }
        } else {
            let c_depth = level+1;
            // Always a Function Node
            if let Some((ref func_constraints,_)) = opt_constraints {
                for i in 0..func_node.fnc.arity { 
                    // We have syntactic function constraints for this function
                    // so loop until the randomly selected FunctionNode qualifies
                    let rnd_fn = FunctionNode::new_rnd_with_constraints(
                                rng, funcs, &func_constraints[i as usize]);

                    let node: &mut Node = func_node.set_arg(i, FNode(rnd_fn));
                    match *node {
                        FNode(ref mut fn_ref) =>
                            Self::gen_tree_full_method_r(rng, fn_ref, c_depth,
                                    depth, funcs, terms, opt_constraints),
                        _ => panic!("expected FunctionNode"),
                    }
                }
            }
            else {
                // No function constraints so anything in function set (funcs)
                // is allowed.
                for i in 0..func_node.fnc.arity { 
                    let rnd_fn = FunctionNode::new_rnd(rng, funcs);
                    let node: &mut Node = func_node.set_arg(i, FNode(rnd_fn));

                    match *node {
                        FNode(ref mut fn_ref) =>
                            Self::gen_tree_full_method_r(rng, fn_ref, c_depth,
                                    depth, funcs, terms, opt_constraints),
                        _ => panic!("expected FunctionNode"),
                    }
                }
            }
        }
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    fn gen_tree_grow_method(rng: &mut GpRng, depth: u16) -> Tree {
        assert_eq!(CONTROL.terms_fdb.len(), CONTROL.funcs_fdb.len());

        let mut result_branch_root =
            FunctionNode::new_rnd(rng, &CONTROL.funcs_rpb[0]);

        Self::gen_tree_grow_method_r(rng, &mut result_branch_root, 2, depth,
                &CONTROL.funcs_rpb[0], &CONTROL.terms_rpb[0]);

        if CONTROL.terms_fdb.len() > 0 {
            let mut func_def_branches: Vec<TreeBranch> = Vec::new();
            for i in 0..CONTROL.funcs_fdb.len() {
                let mut func_def_branch_root =
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_fdb[i]);
                Self::gen_tree_grow_method_r(rng, &mut func_def_branch_root,
                        2, depth, &CONTROL.funcs_fdb[i], &CONTROL.terms_fdb[i]);
                func_def_branches
                    .push(TreeBranch::new(FNode(func_def_branch_root)));
            }
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)),
                Some(func_def_branches))
        } else {
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), None)
        }
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    fn gen_tree_grow_method(rng: &mut GpRng, depth: u16) -> Tree {
        assert_eq!(CONTROL.terms_fdb.len(), CONTROL.funcs_fdb.len());

        let mut result_branch_root = {
                let root_constraints: Box<Vec<NodeConsFTPair>>;
                if let Some(root_constraints) = &CONTROL.opt_rpb_root_constraints {
                    let (f_set, _) = root_constraints[0];
                    FunctionNode::new_rnd_constrained(
                        rng, &CONTROL.funcs_rpb[0], &f_set)
                } else {
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_rpb[0])
                }
            };

        if CONTROL.terms_fdb.len() > 0 {
            let mut func_def_branches: Vec<TreeBranch> = Vec::new();
            for i in 0..CONTROL.funcs_fdb.len() {
                let mut func_def_branch_root =
                    FunctionNode::new_rnd(rng, &CONTROL.funcs_fdb[i]);
                Self::gen_tree_grow_method_r(rng, &mut func_def_branch_root,
                        2, depth, &CONTROL.funcs_fdb[i], &CONTROL.terms_fdb[i]);
                func_def_branches
                    .push(TreeBranch::new(FNode(func_def_branch_root)));
            }
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)),
                Some(func_def_branches))
        } else {
            Tree::new(
                TreeBranch::new(FNode(result_branch_root)), None)
        }
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    fn gen_tree_grow_method_r(rng: &mut GpRng, 
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static [Function], terms: &'static [Terminal]) {
        if level >= depth {
            for i in 0..func_node.fnc.arity {
                // Always a Terminal Node
                let rnd_tref = Terminal::get_rnd_ref(rng, terms);
                func_node.set_arg(i, TNode(rnd_tref));
            }
        }
        else {
            let c_depth = level+1;
            for i in 0..func_node.fnc.arity {
                // Either a Function or Terminal Node
                let rnd_ft_node = Node::new_rnd(rng, funcs, terms);
                let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
                if let FNode(ref mut fn_ref) = node {
                    Self::gen_tree_grow_method_r(rng, fn_ref, c_depth, depth,
                                    funcs, terms);
                }
            }
        }
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    fn gen_tree_grow_method_r(rng: &mut GpRng, 
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static [Function], terms: &'static [Terminal]) {
        if level >= depth {
            // Always a Terminal Node
            if let Some((_, ref term_constraints)) =
                    func_node.fnc.opt_constraints {
                for i in 0..func_node.fnc.arity {
                    loop {
                        // pick random terminals until one qualifies
                        let rnd_tref = Terminal::get_rnd_ref(rng, terms);
                        if term_constraints[i as usize].iter().any(|&x| x == rnd_tref.tid) {
                            func_node.set_arg(i, TNode(rnd_tref));
                            break;
                        }
                    }
                }
            } else {
                for i in 0..func_node.fnc.arity {
                    let rnd_tref = Terminal::get_rnd_ref(rng, terms);
                    func_node.set_arg(i, TNode(rnd_tref));
                }
            }
        }
        else {
            let c_depth = level+1;

            if let Some((ref term_constraints, ref func_constraints)) =
                func_node.fnc.opt_constraints {
                for i in 0..func_node.fnc.arity {
                    // Either a Function or Terminal Node
                    let rnd_ft_node =
                        Node::new_rnd_with_constraints(rng, funcs, terms,
                            &term_constraints[i as usize],
                            &func_constraints[i as usize]);
                    let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
                    if let FNode(ref mut fn_ref) = node {
                        Self::gen_tree_grow_method_r(rng, fn_ref, c_depth, depth,
                                        funcs, terms);
                    }
                }
            } else {
                for i in 0..func_node.fnc.arity {
                    // Either a Function or Terminal Node
                    let rnd_ft_node = Node::new_rnd(rng, funcs, terms);
                    let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
                    if let FNode(ref mut fn_ref) = node {
                        Self::gen_tree_grow_method_r(rng, fn_ref, c_depth, depth,
                                        funcs, terms);
                    }
                }
            }
        }
    }
    fn gen_tree(rng: &mut GpRng, generate_method : GenerateMethod,
            d: u16) -> Tree {
        match generate_method {
            GenerateMethod::Full => Self::gen_tree_full_method(rng, d),
            GenerateMethod::Grow => Self::gen_tree_grow_method(rng, d),
        }
    }
    fn create_unique_tree(&self, rng: &mut GpRng, mut d: u16) -> Tree {
        let mut i = 1u16;
        let gen_method = 
            if (self.tree_vec.len() % 2) == 0 {
                GenerateMethod::Full
            }
            else {
                GenerateMethod::Grow
            };

        let mut t = Self::gen_tree(rng, gen_method, d);

        while self.tree_match_exists(&t) {
            i += 1;
            if i > 11 {
                d += 1;
                i = 1;
            }
            t = Self::gen_tree(rng, gen_method, d);
        }

        return t;
    }
    pub fn create_initial_population(rng: &mut GpRng) -> TreeSet {
        // Following Koza's recipe, he calls "ramped half-and-half",
        // we will evenly produce population segments starting with 2
        // up to the maxium depth (CONTROL.Di) and alternate
        // between Full Method and Grow Method for S Expressions.
        let mut trees = Self::new(0);
        let seg = CONTROL.M as GpFloat / (CONTROL.Di as GpFloat - 1.0);
        let mut bdr = 0.0;

        #[cfg(gpopt_trace="on")]
        println!("TP001:create_init_pop start");
        for d in 2..=CONTROL.Di {
            bdr += seg;
            while trees.tree_vec.len() < bdr as usize &&
                  trees.tree_vec.len() < CONTROL.M {
                let mut new_tree = trees.create_unique_tree(rng, d);
                new_tree.count_nodes();

                trees.push_tree(new_tree);

                #[cfg(gpopt_trace="on")]
                trees.tree_vec[trees.tree_vec.len()-1].print();
            }
        }

        // fill out to end in case there are "left-overs" due to rounding
        while trees.tree_vec.len() < CONTROL.M {
            let mut new_tree = trees.create_unique_tree(rng, CONTROL.Di);
            new_tree.count_nodes();
            trees.push_tree(new_tree);
                
            #[cfg(gpopt_trace="on")]
            trees.tree_vec[trees.tree_vec.len()-1].print();
        }
        #[cfg(gpopt_trace="on")]
        println!("TP002:create_init_pop done");

        trees
    }
    pub fn get_rnd_tree(&mut self, rng: &mut GpRng)
            -> &mut Tree {
        let rnd_index: usize = rng.gen_range(0..self.tree_vec.len() as i32) as usize;
        &mut self.tree_vec[rnd_index]
    }
    pub fn print(&self) {
        for tree in self.tree_vec.iter() {
            tree.print_tree();
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
                i, t.tcid, t.fitness.hits, t.fitness.n(), t.fitness.a(), sum_a);
        }
        #[cfg(gpopt_trace="off")]
        for t in self.tree_vec.iter_mut() {
            t.fitness.n = t.fitness.a() / sum_a;
        }

        self
    }
    pub fn sort_by_normalized_fitness(&mut self) -> &mut TreeSet {
        self.tree_vec
            .sort_by(|a, b| a.fitness.n.partial_cmp(&b.fitness.n).unwrap());
            
        self.assign_nfr_rankings()
    }
    fn assign_nfr_rankings(&mut self) -> &mut TreeSet {
        let mut nfr: GpFloat = 0.0;
        for (i, t) in self.tree_vec.iter_mut().enumerate() {
            nfr += t.fitness.n;
            t.fitness.nfr = nfr;
            t.tfid = Some(i);
        }
        self
    }
    pub fn exec_slice(trees: &mut [Tree]) -> GpHits {
        let mut rc = RunContext::new();
        for tree in trees.iter_mut() {
            rc.prepare_run();

            // init accumulators
            let mut sum_hits: GpHits = 0;
            let mut sum_error: GpRaw = 0;
            if let Some(func_def_branches) = &tree.opt_func_def_branches {
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
                let result = tree.exec_tree(&mut rc);

                let error = rc.compute_error(result);
                sum_error += error;
                if error == 0 {
                    sum_hits += 1;
                }
            }
            if let Some(_) = rc.opt_func_def_branches {
                rc.opt_func_def_branches = None;
            }
            rc.hits = sum_hits;
            rc.error = sum_error;

            let (f, is_winner) = rc.compute_fitness();
            tree.fitness = f;
            tree.is_winner = EvalCount::inc(is_winner);
        }
        rc.hits
    }
    #[cfg(gpopt_mult_threaded="yes")]
    fn exec_chunks(chunks: &mut [&mut [Tree]]) -> GpHits {
        chunks.par_iter_mut()
            .map(|chunk| Self::exec_slice(chunk))
            .sum()
    }
    #[cfg(gpopt_mult_threaded="yes")]
    pub fn exec_all(&mut self) -> GpHits {
        let num_chunks: usize = 50;
        let chunk_size: usize = (self.tree_vec.len() / num_chunks) + 1;

        let mut chunks: Vec<&mut [Tree]> = Vec::new();
        for chunk in self.tree_vec.chunks_mut(chunk_size) {
            chunks.push(chunk);
        }

        let total_hits = Self::exec_chunks(&mut chunks);

        self.winning_index = None;
        for (t_i, t) in self.tree_vec.iter().enumerate() {
            if t.is_winner {
                self.winning_index = Some(t_i);
            }
        }

        total_hits
    }
    #[cfg(gpopt_mult_threaded="no")]
    pub fn exec_all(&mut self) -> GpHits {
        let total_hits = Self::exec_slice(&mut self.tree_vec);

        self.winning_index = None;
        for (t_i, t) in self.tree_vec.iter().enumerate() {
            if t.is_winner {
                self.winning_index = Some(t_i);
            }
        }

        total_hits
    }
    //fn use_reproduction(index: usize ) -> bool {
    //    let float_index_ratio = index as GpFloat / CONTROL.M as GpFloat;
    //    let int_index_ratio = Fitness::float_to_int(
    //        float_index_ratio
    //    );
    //    let int_control_ratio = Fitness::float_to_int(CONTROL.Pr);
    //    let result = int_index_ratio < int_control_ratio;
    //    result
    //}
    fn use_reproduction(&self) -> bool {
        let index = self.tree_vec.len();
        let result = (index as GpFloat / CONTROL.M as GpFloat) < CONTROL.Pr;
        result
    }
    fn push_tree(&mut self, mut tree: Tree) {
        let index = self.tree_vec.len();
        assert!(index < CONTROL.M);
        tree.tcid = index;
        self.tree_vec.push(tree);
    }
    // rnd_int_pt_decide - randomly decides whether to do crossover at an
    // internal point (function) or terminal based on control Pip value.
    fn rnd_int_pt_decide(rng: &mut GpRng) -> bool {
        #[cfg(gpopt_rng="file_stream")]
        let num: GpFloat = rng.gen_float();

        #[cfg(not(gpopt_rng="file_stream"))]
        let num: GpFloat = rng.gen_range(0.0..1.0);

        num < CONTROL.Pip // if Pip is .90 then true for all values less than .90.
    }
    #[cfg(gpopt_syntactic_constraints="no")] 
    fn perform_crossover(rng: &mut GpRng, t1: &mut Tree, t2: &mut Tree) {
        assert_ne!(t1.get_num_terminal_nodes(), None);
        assert_ne!(t2.get_num_terminal_nodes(), None);

        let node:  &mut Node;
        let b_type: BranchType;

        let swap_target1 =
            if t1.get_num_function_nodes().unwrap() > 0 && Self::rnd_int_pt_decide(rng) {
                let nref_pair = t1.get_rnd_function_node_ref(rng);
                b_type = nref_pair.0;
                node = nref_pair.1;
                node
            } else {
                let nref_pair = t1.get_rnd_terminal_node_ref(rng);
                b_type = nref_pair.0;
                node = nref_pair.1;
                node
            };

        let swap_target2 =
            if t2.get_num_function_nodes_bt(&b_type).unwrap() > 0 && Self::rnd_int_pt_decide(rng) {
                t2.get_rnd_function_node_ref_bt(rng, &b_type)
            } else {
                t2.get_rnd_terminal_node_ref_bt(rng, &b_type)
            };

        mem::swap(swap_target1, swap_target2);
    }
    #[cfg(gpopt_syntactic_constraints="yes")] 
    fn perform_crossover_with_syntatic_constraints(rng: &mut GpRng,
                t1: &mut Tree, t2: &mut Tree) -> bool {
        assert_ne!(t1.get_num_terminal_nodes(), None);
        assert_ne!(t2.get_num_terminal_nodes(), None);

        let mut opt_f_set: Option<&NodeConstraints> = None;
        let mut opt_t_set: Option<&NodeConstraints> = None;

        // Swap target1 Node is randomly chosen with across all t1 nodes ( favoring internal nodes
        // per rnd_int_pt_decide control. If the selection is t1's root node, then the Syntactic
        // constraints come from Control.opt_rpb_root_constraints, otherwise the parent FunctionNode along
        // with the Node's argument position determine the Function and Terminal set (FTSet) of
        // qualifying candidates to confirm.
        let node:  &mut Node;
        let btype: BranchType;
        let swap_target1 =
            if t1.get_num_function_nodes().unwrap() > 0 && Self::rnd_int_pt_decide(rng) {
                let fnode_loc: (BranchType, (&mut Node, Option<(&Function, usize)>))
                            = t1.get_rnd_function_node_ref_ploc(rng);

                let opt_ploc: Option<(&Function, usize)>;
                {
                    // In future following unpack will only require one
                    // statement.
                    // (Currently only unstable Rust permits doing this without
                    //  temporary variables.)
                    let (t_btype, (t_node, t_opt_ploc)) = fnode_loc;
                    btype = t_btype;
                    node = t_node;
                    opt_ploc = t_opt_ploc;
                }

                if let Some((fnc, arg_num)) = opt_ploc {
                    // Note: this branch is taken unless node is a root node.
                    let constraints: &ArgNodeConsFTPair;
                    if let Some(ref constraints) = fnc.opt_constraints {
                        let (f_cons, t_cons) = *constraints;
                        opt_f_set = Some(&f_cons[arg_num]);
                        opt_t_set = Some(&t_cons[arg_num]);
                    }
                } else {
                    // Note: node is a root node
                    // only rpb can have root level constraints (for now)
                    // eventually all roots should have a possible fset in Control.
                    let rpb_root_constraints: Box<Vec<NodeConsFTPair>>;
                    if let ResultProducing = btype {
                        if let Some(rpb_root_constraints) =
                                CONTROL.opt_rpb_root_constraints {

                            let (ref f_cons, ref t_cons)
                                = rpb_root_constraints[0];
                            if f_cons.len() > 0 {
                                opt_f_set = Some(f_cons);
                            }
                            if t_cons.len() > 0 {
                                opt_t_set = Some(t_cons);
                            }
                        }
                    }
                }

                node
            } else {
                let nref_pair = t1.get_rnd_terminal_node_ref(rng);
                btype = nref_pair.0;
                node = nref_pair.1;
                node
            };

        let swap_target2 =
            if let Some(f_set) = opt_f_set {
                // apply syntactic constraint
                if t2.get_num_function_nodes_bt(&btype).unwrap() > 0
                   && Self::rnd_int_pt_decide(rng) {
                    // loop until a random function node qulaifies as swap point 2
                    // note that we currently do a shallow (avoid deep) check as it's not
                    // required for this program. Future versions may require
                    // both a deep and a shallow check at this point depending on 
                    // the makeup of the constraints.
                    // OR until attempts_left hits zero.
                    let mut attempts_left = t2.get_num_function_nodes_bt(&btype)
                            .unwrap() * 2usize;
                    loop {
                        let node = t2.get_rnd_function_node_ref_bt(rng, &btype);
                        if let FNode(func_node) = node {
                            if f_set.iter()
                                    .any(|&x| x == func_node.fnc.fid) {
                                break node;
                            }
                        }
                        else {
                            panic!("not func_node");
                        }
                        attempts_left -= 1;
                        if attempts_left == 0 {
                            return false;
                        }
                    }
                } else {
                    if let Some(t_set) = opt_t_set {
                        // loop until a random terminal node qulaifies as swap point 2
                        // OR until attempts_left hits zero.
                        let mut attempts_left = t2.get_num_terminal_nodes_bt(&btype)
                            .unwrap() * 2usize;
                        loop {
                            let node = t2.get_rnd_terminal_node_ref_bt(rng, &btype);
                            if let TNode(tn_ref) = node {
                                if t_set.iter()
                                        .any(|&x| x == tn_ref.tid) {
                                    break node;
                                }
                            }
                            else {
                                panic!("not func_node");
                            }
                            attempts_left -= 1;
                            if attempts_left == 0 {
                                return false;
                            }
                        }
                    } else {
                        panic!("opt_f_set not None, but opt_t_set None");
                    }
                }
            } else if t2.get_num_function_nodes_bt(&btype).unwrap() > 0
                && Self::rnd_int_pt_decide(rng) {
                    t2.get_rnd_function_node_ref_bt(rng, &btype)
            } else {
                t2.get_rnd_terminal_node_ref_bt(rng, &btype)
            };

        mem::swap(swap_target1, swap_target2);
        true
    }
    pub fn breed_new_generation(&mut self, rng: &mut GpRng) -> TreeSet {
        let mut new_trees = Self::new(self.gen);
        #[cfg(gpopt_trace="on")]
        println!("TP003:breed start");
        while new_trees.tree_vec.len() < CONTROL.M {
            if new_trees.use_reproduction() {
                // do reproduction
                #[cfg(gpopt_trace="on")]
                println!("TP003.1:breeding reproduction branch");

                let mut t = self.select_tree(rng).clone();
                
                t.clear_node_counts();
                t.count_nodes();
                new_trees.push_tree(t);

                #[cfg(gpopt_trace="on")]
                new_trees.tree_vec[new_trees.tree_vec.len()-1].print();
            }
            else {
                // do crossover
                #[cfg(gpopt_trace="on")]
                println!("TP003.1:breeding crossover branch");

                let (mut nt1, mut nt2);
                loop {
                    let t1 = self.select_tree(rng);
                    let t2 = self.select_tree(rng);

                    nt1 = t1.clone();
                    nt2 = t2.clone();

                    // when syntactic_constraints is off
                    // we use both trees made with crossover
                    // when syntactic_constraints is on
                    // we toss the second, because constraints 
                    // are only enforced for the first one.
                    #[cfg(gpopt_syntactic_constraints="no")] 
                    {
                        Self::perform_crossover(rng, &mut nt1, &mut nt2);
                        if nt1.qualifies() && nt2.qualifies() {
                            break;
                        }
                    }
                    #[cfg(gpopt_syntactic_constraints="yes")] 
                    {
                        Self::perform_crossover_with_syntatic_constraints(rng,
                            &mut nt1, &mut nt2);
                        if nt1.qualifies() {
                            break;
                        }
                    }
                }
                nt1.clear_node_counts();
                nt1.count_nodes();

                new_trees.push_tree(nt1);

                #[cfg(gpopt_trace="on")]
                new_trees.tree_vec[new_trees.tree_vec.len()-1].print();

                #[cfg(gpopt_syntactic_constraints="no")] 
                if new_trees.tree_vec.len() < CONTROL.M {
                    nt2.clear_node_counts();
                    nt2.count_nodes();
                    new_trees.push_tree(nt2);

                    #[cfg(gpopt_trace="on")]
                    new_trees.tree_vec[new_trees.tree_vec.len()-1].print();
                }
            }
        }
        #[cfg(gpopt_trace="on")]
        println!("TP004:breed done");

        assert_eq!(new_trees.tree_vec.len(), CONTROL.M);
        new_trees
    }
}

#[cfg(gpopt_select_method="fpb")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpFloat, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpFloat) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
    fn rnd_greedy_val(rng: &mut GpRng) -> GpFloat;
}

#[cfg(gpopt_select_method="tournament")]
pub trait SelectMethod {
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
}

#[cfg(gpopt_select_method="fpb")]
impl SelectMethod for TreeSet {
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
        let greedy_val = Self::rnd_greedy_val(rng);
        let i = self.select_ind_bin(greedy_val);
        return &self.tree_vec[i];
    }

    #[cfg(gpopt_select_method="fpb")]
    fn rnd_greedy_val(rng: &mut GpRng) -> GpFloat {
        #[cfg(not(gpopt_rng="file_stream"))]
        let r = rng.gen::<GpFloat>();  // Note optional choice logging for this function
                                   // done in TreeSet::select_tree, which is
                                   // the only function that calls here. This allows
                                   // integer logging thereby removing floating point variences.

        #[cfg(gpopt_rng="file_stream")]
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
        dbl_val
    }
}

#[cfg(gpopt_select_method="tournament")]
impl SelectMethod for TreeSet {
    fn select_tree(&self, rng: &mut GpRng) -> &Tree {
        // perform tournament selection by choosing the tree with best 
        // normalized fitness among a random set of 7 trees.
        assert!(!self.tree_vec.len() > 50);

        let last_ti = self.tree_vec.len() - 1;
        let mut r: usize = rng.gen_range(0..=last_ti);
        let mut best_tree = &self.tree_vec[r];

        let mut t_group: HashMap<usize, bool> = HashMap::new();

        t_group.insert(r, true); // tag the slot as used
        for _ in 2..=7 {
            let mut sanity: u8 = 100;

            // loop until we find an r not already in t_group
            while t_group.contains_key(&r) {
                r = rng.gen_range(0..=last_ti);

                sanity -= 1u8;
                if sanity < 1 {
                    panic!("sanity limit reached while doing tournament selection");
                }
            }
            t_group.insert(r, true); // tag the slot as used

            if best_tree.fitness.n < self.tree_vec[r].fitness.n {
                best_tree = &self.tree_vec[r];
            }
        }
        best_tree
    }
}
