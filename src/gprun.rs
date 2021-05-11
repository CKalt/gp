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

pub static RUN_CONTROL: RunControl = RunControl{
    max_clock: 3000,
};

pub struct RunContext {
    pub grid: [[u8; RUN_CONTROL_MAX_X as usize + 1usize]; 
                    RUN_CONTROL_MAX_Y as usize + 1usize ],
    pub ant_x: i16,
    pub ant_y: i16,
    pub ant_xd: i16,
    pub ant_yd: i16,
    pub eat_count: u16,
    pub clock: u16,
    pub n_pellets: u16,
}
impl RunContext {
    fn new() -> RunContext {
        RunContext {
            grid: [[0; RUN_CONTROL_MAX_X as usize + 1usize];
                       RUN_CONTROL_MAX_Y as usize + 1usize],
            ant_x: 0,
            ant_y: 0,
            ant_xd: 1,
            ant_yd: 0,
            eat_count: 0,
            clock: 0,
            n_pellets: 0,
        }
    }
}

pub fn init_run() {
    println!("init_run not implemented");
}

pub fn prepare_run(_rc: &mut RunContext) {
    // STUB
}

pub fn exec_trees(mut trees: &mut TreeSet) {
    let mut rc = RunContext::new();

    for (i, tree_ref) in trees.tree_vec.iter_mut().enumerate() {
        prepare_run(&mut rc);
        while rc.clock < RUN_CONTROL.max_clock {
            exec_node(&mut rc, &mut tree_ref.root);
        }

        if compute_fitness(tree_ref) {
            report_tree_result(tree_ref, tree_ref.tfid, 0, -1.0);
            print_grid(&rc, "Have Winner!");
            trees.winning_index = Some(i);
            break;
        }
    }
}

pub fn exec_single_tree(_tree : &Tree) {
    println!("exec_single_tree not implemented");
}

fn compute_fitness(_t: &mut Tree) -> bool {
    true
}


fn report_tree_result(_t: &Tree, _i: usize , _gen: u16, _avg_raw_f: f64) {
    // STUB
}

fn print_grid(_rc: &RunContext, _label: &str) {
    // STUB
}

