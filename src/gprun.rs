use crate::tree::Tree;
use crate::tree::TreeSet;
use crate::tree::exec_node;

/// RunControl defines parameters controlling the running of individuals.
pub struct RunControl {
    pub max_clock: u16, // sets limit for entire population run duration
}

// LOS_ALTOS_TRAIL

// MAX_X & MAX_Y are zero based max_y (may be equal but cannot exceed)
pub const RUN_CONTROL_MAX_X: i16 = 66;
pub const RUN_CONTROL_MAX_Y: i16 = 43;

type FoodCoord = (usize, usize); // (x, y)

pub static RUN_CONTROL: RunControl = RunControl{
    max_clock: 3000,
};

#[derive(Copy, Clone)]
pub enum GridCellState {
    Clear,
    Food,
    FoodEaten,
    NoFoodFound,
}
use GridCellState::*;

pub struct RunContext {
    pub grid: [[GridCellState; RUN_CONTROL_MAX_X as usize + 1usize]; 
                               RUN_CONTROL_MAX_Y as usize + 1usize],
    pub ant_x: i16,
    pub ant_y: i16,
    pub ant_xd: i16,
    pub ant_yd: i16,
    pub eat_count: u16,
    pub clock: u16,
    pub n_pellets: u16,
    food: Vec<FoodCoord>, // stores pattern for food trail transfered 
                          // to grid during init_grid.
}
impl RunContext {
    fn new() -> RunContext {
        RunContext {
            grid: [[Clear; RUN_CONTROL_MAX_X as usize + 1usize];
                                   RUN_CONTROL_MAX_Y as usize + 1usize],
            ant_x: 0,
            ant_y: 0,
            ant_xd: 1,
            ant_yd: 0,
            eat_count: 0,
            clock: 0,
            n_pellets: 0,
            // Food for Koza's Los Altos Trail
            food: vec![ // set of FoodCoord (x,y) tuples
    (1,0), (2,0), (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (4,5), (5,5), 
    (6,5), (8,5), (9,5), (10,5), (11,5), (12,5), (12,6), (12,7), (12,8),
    (12,9), (12,11), (12,12), (12,13), (12,14), (12,17), (12,18), (12,19),
    (12,20), (12,21), (12,22), (12,23), (11,24), (10,24), (9,24), (8,24),
    (7,24), (4,24), (3,24), (1,25), (1,26), (1,27), (1,28), (2,30), (3,30),
    (4,30), (5,30), (7,29), (7,28), (8,27), (9,27), (10,27), (11,27), (12,27),
    (13,27), (14,27), (16,26), (16,25), (16,24), (16,21), (16,20), (16,19),
    (16,18), (17,15), (20,14), (20,13), (20,10), (20,9), (20,8), (20,7), 
    (21,5), (22,5), (24,4), (24,3), (25,2), (26,2), (27,2), (29,3), (29,4),
    (29,6), (29,9), (29,12), (28,14), (27,14), (26,14), (23,15), (24,18), 
    (27,19), (26,22), (23,23),
            ]
        }
    }
    fn init_grid(&mut self) {
        // clear grid
        self.grid = [[Clear; RUN_CONTROL_MAX_X as usize +1];
                                     RUN_CONTROL_MAX_Y as usize +1];

        // sprinkle in the food!
        for food in self.food.iter() {
            self.grid[food.1][food.0] = Food;
        }
        self.n_pellets = self.food.len() as u16;
    }
    fn prepare_run(&mut self) {
        self.init_grid();
        self.ant_x = 0;
        self.ant_y = 0;
        self.ant_xd = 1;
        self.ant_yd = 0;
        self.eat_count = 0;
        self.clock = 0;
    }
    fn print_grid(&self, label: &str) {
        println!("{}", label);
        for y in 0..=RUN_CONTROL_MAX_Y as usize {
            for x in 0..=RUN_CONTROL_MAX_X as usize {
                print!("{}",
                    match self.grid[y][x] {
                        Clear => ".",
                        Food => "X",
                        FoodEaten => "@",
                        NoFoodFound => "O",
                    });
            }
            println!("-----------------------------------\n");
        }
    }
}

pub fn init_run() {
    println!("init_run not implemented");
}

pub fn exec_trees(mut trees: &mut TreeSet) {
    let mut rc = RunContext::new();

    for (i, tree) in trees.tree_vec.iter_mut().enumerate() {
        rc.prepare_run();
        while rc.clock < RUN_CONTROL.max_clock {
            exec_node(&mut rc, &mut tree.root);
        }

        if tree.compute_fitness(&rc) {
            report_tree_result(tree, tree.tfid, None, -1.0);
            rc.print_grid("Have Winner!");
            trees.winning_index = Some(i);
            break;
        }
    }
}

pub fn exec_single_tree(tree : &mut Tree) {
    let mut rc = RunContext::new();
    rc.prepare_run();
    rc.print_grid("Before Run");
    while rc.clock < RUN_CONTROL.max_clock {
        exec_node(&mut rc, &mut tree.root);
    }

    if tree.compute_fitness(&rc) {
        println!("Have Winner");
    }
    rc.print_grid("After Run");
}

pub fn report_tree_result(t: &Tree, i: Option<usize> , opt_gen: Option<u16>, avg_raw_f: f64) {
    assert_eq!(i, t.tfid);
    let f = &t.fitness;
    let tfid = if let Some(num) = t.tfid { num } else { 0 };
    if let Some(gen) = opt_gen {
        println!("{:6} {:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                 gen, tfid, t.tcid, t.hits, f.r, f.s, f.a, f.n, f.nfr, avg_raw_f);
    } else {
        println!("{:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                tfid, t.tcid, t.hits, f.r, f.s, f.a, f.n, f.nfr, avg_raw_f);
    }
}

pub fn tree_result_header(opt_gen: Option<u16>) {
    if let Some(_) = opt_gen {
        println!("{:6} {:4} {:4} {:6} {:6} {:6} {:8} {:8} {:8} {:6}", 
            "gen", "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr", "avgRawF");
        println!("{:6} {:4} {:4} {:6} {:6} {:6} {:8} {:8} {:8}", 
            "----", "---", "---", "-----", "---", "---", "------", "------", "------");
    } else {
        println!("{:4} {:4} {:6} {:6} {:6} {:8} {:8} {:8}", 
            "tfid", "tcid", "hits", "r", "s", "a", "n", "nfr");
        println!("{:4} {:4} {:6} {:6} {:6} {:8} {:8} {:8}", 
            "---", "---", "-----", "---", "---", "------", "------", "------");
    }
}

