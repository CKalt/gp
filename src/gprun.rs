use crate::tree::Tree;
use crate::tree::TreeSet;
use crate::tree::exec_node;
use crate::tree::GpFloat;

use crate::tree::Function;
use crate::tree::FunctionNode;
use crate::tree::Terminal;
use crate::control::CONTROL;

pub enum GpType {
    Continue,
    Terminate,
}

// TERMINAL SPECIFICS
pub static TERMINAL: [Terminal; CONTROL.num_terminals as usize] = [
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

pub static FUNCTION: [Function; CONTROL.num_functions as usize] = [
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


/// RunControl defines parameters controlling the running of individuals.
pub struct RunControl {
    pub max_clock: u16, // sets limit for entire population run duration
}

//////////////////////
// SANTA_FE TRAIL   //
//////////////////////

// MAX_X & MAX_Y are zero based max_y (may be equal but cannot exceed)
//////////////////////
// LOS_ALTOS TRAIL  //
//////////////////////

// MAX_X & MAX_Y are zero based max_y (may be equal but cannot exceed)
pub const RUN_CONTROL_MAX_X: i16 = 43;
pub const RUN_CONTROL_MAX_Y: i16 = 66;

fn trail_factory() -> Vec<FoodCoord> {
    vec![ // set of FoodCoord (x,y) tuples
    (1,0), (2,0), (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (4,5), 
    (5,5), (6,5), (8,5), (9,5), (10,5), (11,5), (12,5), (12,6), (12,7), 
    (12,8), (12,9), (12,11), (12,12), (12,13), (12,14), (12,17), (12,18), 
    (12,19), (12,20), (12,21), (12,22), (12,23), (11,24), (10,24), (9,24), 
    (8,24), (7,24), (4,24), (3,24), (1,25), (1,26), (1,27), (1,28), 
    (1,29), (1,30), (1,31), (2,33), (3,33), (4,33), (5,33), (7,32), (7,31),
    (8,30), (9,30), (10,30), (11,30), (12,30), (17,30), (18,30), (20,29),
    (20,28), (20,27), (20,26), (20,25), (20,24), (20,21), (20,20), (20,19),
    (20,18), (21,15), (24,14), (24,13), (24,10), (24,9), (24,8), (24,7),
    (25,5), (26,5), (27,5), (28,5), (29,5), (30,5), (31,5), (32,5), (33,5), 
    (34,5), (35,5), (36,5), (38,4), (38,3), (39,2), (40,2), (41,2), (43,3), 
    (43,4), (43,6), (43,9), (43,12), (42,14), (41,14), (40,14), (37,15), 
    (38,18), (41,19), (40,22), (37,23), (37,24), (37,25), (37,26), (37,27), 
    (37,28), (37,29), (37,30), (37,31), (37,32), (37,33), (37,34), (35,34), 
    (35,35), (35,36), (35,37), (35,38), (35,39), (35,40), (35,41), (35,42), 
    (34,42), (33,42), (32,42), (31,42), (30,42), (30,44), (30,45), (30,46), 
    (30,47), (30,48), (30,49), (28,50), (28,51), (28,52), (28,53), (28,54), 
    (28,55), (28,56), (27,56), (26,56), (25,56), (24,56), (23,56), (22,58), 
    (22,59), (22,60), (22,61), (22,62), (22,63), (22,64), (22,65), (22,66), 
    ]
}

pub static RUN_CONTROL: RunControl = RunControl{
    max_clock: 3000,
};

////////////////////////////////////////////////////////////////////
type FoodCoord = (usize, usize); // (x, y)
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
            food: trail_factory()
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
        println!("---------------------------------------------------------------------------------------");
        for y in 0..=RUN_CONTROL_MAX_Y as usize {
            for x in 0..=RUN_CONTROL_MAX_X as usize {
                print!("{}",
                    match self.grid[y][x] {
                        Clear => ". ",
                        Food => "X ",
                        FoodEaten => "@ ",
                        NoFoodFound => "O ",
                    });
            }
            println!();
        }
        println!("---------------------------------------------------------------------------------------");
    }
}

pub fn init_run() { }

pub fn exec_trees(mut trees: &mut TreeSet, run_number: i32) -> u16 {
    let mut rc = RunContext::new();

    for (i, tree) in trees.tree_vec.iter_mut().enumerate() {
        rc.prepare_run();
        while rc.clock < RUN_CONTROL.max_clock {
            exec_node(&mut rc, &mut tree.root);
        }

        if tree.compute_fitness(&rc) {
            report_tree_result(tree, tree.tfid, None, -1.0);
            rc.print_grid(&format!("Have Winner! - Run# {} Gen# {}", run_number,
                trees.gen));
            trees.winning_index = Some(i);
            break;
        }
    }
    rc.n_pellets
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

pub fn report_tree_result(t: &Tree, i: Option<usize> , opt_gen: Option<u16>, avg_raw_f: GpFloat) {
    assert_eq!(i, t.tfid);
    let f = &t.fitness;
    let tfid = if let Some(num) = t.tfid { num } else { 0 };
    if let Some(gen) = opt_gen {
        println!("{:6} {:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                 gen, tfid, t.tcid, t.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
    } else {
        println!("{:4} {:4} {:6} {:6} {:6} {:6.6} {:6.6} {:6.6} {:6.2}", 
                tfid, t.tcid, t.hits, f.r, f.s, f.a(), f.n(), f.nfr(), avg_raw_f);
    }
}

pub fn tree_result_header(opt_gen: Option<u16>, n_pellets: &u16) {
    println!("n_pellets={}", *n_pellets);
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

