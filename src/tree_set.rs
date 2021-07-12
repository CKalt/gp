use std::mem;

#[cfg(gpopt_select_method="tournament")]
use std::collections::HashMap;

use crate::gprun::*;
use crate::tree::*;
use crate::control::CONTROL;
use crate::gprng::GpRng;
use rand::Rng;

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
            tree_vec:       Vec::new(), // TODO: for better performance change to array (must find good init method)
            gen:            gen
        }
    }
    fn tree_match(t1: &Tree, t2: &Tree) -> bool {
        t1.result_branch.root.deep_match(&t2.result_branch.root) &&
            t1.func_def_branch.root.deep_match(&t2.func_def_branch.root)
    }
    fn tree_match_exists(&self, target_tree: &Tree) -> bool {
        for tree in self.tree_vec.iter() {
            if Self::tree_match(tree, target_tree) {
                return true;
            }
        }
        return false;
    }
    fn gen_tree_full_method(rng: &mut GpRng, depth: u16) -> Tree {
        let mut result_branch_root =
            FunctionNode::new_rnd(rng, &FUNCTIONS_RESULT_BRANCH);
        Self::gen_tree_full_method_r(rng, &mut result_branch_root, 2, depth,
                &FUNCTIONS_RESULT_BRANCH, &TERMINALS_RESULT_BRANCH);

        let mut func_def_branch_root =
            FunctionNode::new_rnd(rng, &FUNCTIONS_FUNC_DEF_BRANCH);
        Self::gen_tree_full_method_r(rng, &mut func_def_branch_root, 2, depth,
                &FUNCTIONS_FUNC_DEF_BRANCH, &TERMINALS_FUNC_DEF_BRANCH);

        Tree::new(result_branch_root, func_def_branch_root)
    }
    fn gen_tree_full_method_r(rng: &mut GpRng,
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static [Function], terms: &'static [Terminal]) {
        if level >= depth {
            for i in 0..func_node.fnc.arity {
                let rnd_tref = Terminal::get_rnd_ref(rng, terms); // Always a Terminal Node
                func_node.set_arg(i, TNode(rnd_tref));
            }
        }
        else {
            let c_depth = level+1;
            for i in 0..func_node.fnc.arity {
                let rnd_fn = FunctionNode::new_rnd(rng, funcs); // Always a Funciton Node
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
    fn gen_tree_grow_method(rng: &mut GpRng, depth: u16) -> Tree {
        let mut result_branch_root =
            FunctionNode::new_rnd(rng, &FUNCTIONS_RESULT_BRANCH);
        Self::gen_tree_grow_method_r(rng, &mut result_branch_root, 2, depth,
                &FUNCTIONS_RESULT_BRANCH, &TERMINALS_RESULT_BRANCH);

        let mut func_def_branch_root =
            FunctionNode::new_rnd(rng, &FUNCTIONS_FUNC_DEF_BRANCH);
        Self::gen_tree_grow_method_r(rng, &mut func_def_branch_root, 2, depth,
                &FUNCTIONS_FUNC_DEF_BRANCH, &TERMINALS_FUNC_DEF_BRANCH);

        Tree::new(result_branch_root, func_def_branch_root)
    }
    fn gen_tree_grow_method_r(rng: &mut GpRng, 
            func_node: &mut FunctionNode, level: u16, depth: u16,
            funcs: &'static [Function], terms: &'static [Terminal]) {
        if level >= depth {
            for i in 0..func_node.fnc.arity {
                let rnd_tref = Terminal::get_rnd_ref(rng, terms); // Always a Terminal Node
                func_node.set_arg(i, TNode(rnd_tref));
            }
        }
        else {
            let c_depth = level+1;
            for i in 0..func_node.fnc.arity {
                let rnd_ft_node = Node::new_rnd(rng, funcs, terms); // Either a Function or Terminal Node
                let node: &mut Node = func_node.set_arg(i, rnd_ft_node);
                if let FNode(ref mut fn_ref) = node {
                    Self::gen_tree_grow_method_r(rng, fn_ref, c_depth, depth,
                                    funcs, terms);
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
    pub fn exec_all(&mut self) -> u16 {
        let mut rc = RunContext::new();
        for (t_i, tree) in self.tree_vec.iter_mut().enumerate() {
            rc.prepare_run();

            #[cfg(gpopt_exec_criteria="clock")]
            panic!("clock exec not used this problem.");
            #[cfg(gpopt_exec_criteria="each_fitness_case")]
            {
                // init accumulators
                let mut sum_hits: GpHits = 0;
                let mut sum_error: GpRaw = 0.0;
                rc.func_def_branch = Some(&tree.func_def_branch);
                for fc_i in 0..rc.fitness_cases.len() {
                    rc.cur_fc = fc_i;
                    let result = Tree::exec_node(&mut rc, &tree.result_branch.root);

                    let error = rc.compute_error(result);
                    sum_error += error;
                    if error < 0.01 {
                        sum_hits += 1;
                    }
                }
                rc.func_def_branch = None;
                rc.hits = sum_hits;
                rc.error = sum_error;
            }

            let (f, is_winner) = rc.compute_fitness();
            tree.fitness = f;
            if is_winner {
                self.winning_index = Some(t_i);
                break;
            }
        }
        rc.hits
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
    // rnd_internal_point - randomly decides whether to do crossover at an
    // internal point (function) or terminal based on control Pip value.
    fn rnd_internal_point(rng: &mut GpRng) -> bool {
        #[cfg(gpopt_rng="file_stream")]
        let num: GpFloat = rng.gen_float();

        #[cfg(not(gpopt_rng="file_stream"))]
        let num: GpFloat = rng.gen_range(0.0..1.0);

        num < CONTROL.Pip // if Pip is .90 then true for all values less than .90.
    }

    fn perform_crossover(rng: &mut GpRng, t1: &mut Tree, t2: &mut Tree) {
        assert_ne!(t1.get_num_terminal_nodes(), None);
        assert_ne!(t2.get_num_terminal_nodes(), None);

        let node:  &mut Node;
        let b_type: BranchType;

        let swap_target1 =
            if t1.get_num_function_nodes().unwrap() > 0 && Self::rnd_internal_point(rng) {
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
            if t2.get_num_function_nodes_bt(&b_type).unwrap() > 0 && Self::rnd_internal_point(rng) {
                t2.get_rnd_function_node_ref_bt(rng, &b_type)
            } else {
                t2.get_rnd_terminal_node_ref_bt(rng, &b_type)
            };

        mem::swap(swap_target1, swap_target2);
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

                    Self::perform_crossover(rng, &mut nt1, &mut nt2);
                    if nt1.qualifies() && nt2.qualifies() {
                        break;
                    }
                }
                nt1.clear_node_counts();
                nt1.count_nodes();

                new_trees.push_tree(nt1);

                #[cfg(gpopt_trace="on")]
                new_trees.tree_vec[new_trees.tree_vec.len()-1].print();

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
#[cfg(gpopt_fitness_type="int")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpInt, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpInt) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
    fn rnd_greedy_val(rng: &mut GpRng) -> GpFloat;
}
#[cfg(gpopt_select_method="fpb")]
#[cfg(gpopt_fitness_type="float")]
pub trait SelectMethod {
    fn select_ind_bin_r(&self, r: GpFloat, lo: usize, hi: usize) -> usize;
    fn select_ind_bin(&self, r: GpFloat) -> usize;
    fn select_tree(&self, rng: &mut GpRng) -> &Tree;
    fn rnd_greedy_val(rng: &mut GpRng) -> GpFloat;
}

#[cfg(gpopt_select_method="tournament")]
#[cfg(gpopt_fitness_type="float")]
pub trait SelectMethod {
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
        let greedy_val = Self::rnd_greedy_val(rng);
        let i = self.select_ind_bin(greedy_val);
        return &self.tree_vec[i];
    }

    #[cfg(gpopt_select_method="fpb")]
    #[cfg(gpopt_fitness_type="int")]
    fn rnd_greedy_val(rng: &mut GpRng) -> GpInt {
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
        Fitness::float_to_int(dbl_val)
    }

    #[cfg(gpopt_select_method="fpb")]
    #[cfg(gpopt_fitness_type="float")]
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

