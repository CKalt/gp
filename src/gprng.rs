#[cfg(gpopt_rng="file_stream")]
use crate::util::read_i32_pair_from_fbuf_rdr;

#[cfg(gpopt_rng="file_stream")]
use std::sync::atomic::{AtomicI32, Ordering};

#[cfg(gpopt_rng="file_stream")]
use crate::util::open_fbuf_rdr;

#[cfg(gpopt_rng="seedable")]
use rand::{SeedableRng};

#[cfg(gpopt_rng="thread")]
pub type GpRng = rand::rngs::ThreadRng;

#[cfg(gpopt_rng="seedable")]
pub type GpRng = rand::rngs::StdRng;

#[cfg(gpopt_rng="file_stream")]
pub type GpRng = FileStreamRng ;

#[cfg(gpopt_rng="file_stream")]
use std::io::BufReader;

#[cfg(gpopt_rng="file_stream")]
use std::ops::Range;

#[cfg(gpopt_rng="file_stream")]
use std::fs::File;

#[cfg(gpopt_rng="file_stream")]
pub struct FileStreamRng {
    fbuf_rdr: BufReader<File>,
}

#[cfg(gpopt_rng="file_stream")]
impl FileStreamRng {
    pub fn new() -> FileStreamRng {
        FileStreamRng {
            fbuf_rdr: open_fbuf_rdr("rng.log"),
        }
    }
    pub fn gen_range(&mut self, r: Range<i32>) -> i32 {
        assert_ne!(r.start, r.end);
        let result = rnd(r.end as i16 - 1, &mut self.fbuf_rdr) as i32;
        if (result < r.start) || (result >= r.end) {
            panic!("value ({}) from input stream out of range [{}..{}]",
                result, r.start, r.end);
        }
        result
    }
    pub fn gen_float(&mut self) -> f64 {
        rnd_dbl(&mut self.fbuf_rdr)
    }
}

pub struct GpRngFactory {
}
impl GpRngFactory {
    #[cfg(gpopt_rng="thread")]
    pub fn new() -> rand::rngs::ThreadRng {
        rand::thread_rng()
    }

    #[cfg(gpopt_rng="seedable")]
    pub fn new() -> rand::rngs::StdRng {
        #[cfg(gpopt_fitness_type="int")]
        let seed: [u8; 32] = [151; 32];     // Ant Winner during Run# 6 Gen# 24

        #[cfg(gpopt_fitness_type="float")]
        let seed: [u8; 32] = [73; 32];      // Ant Winner during Run# 3 Gen# 48
                                            // two-boxes-adf winner during run=23, gen=13

        SeedableRng::from_seed(seed)
    }

    #[cfg(gpopt_rng="file_stream")]
    pub fn new() -> FileStreamRng {
        FileStreamRng::new()
    }
}

#[cfg(gpopt_rng="file_stream")]
const RAND_MAX: i32 = 2147483647;    // defined in gcc <stdlib.h>

#[cfg(gpopt_rng="file_stream")]
fn rnd(max: i16, rng_rdr: &mut BufReader<File>) -> i16 {
    #[cfg(gpopt_trace="on")]
    let (c,r) = read_i32_pair_from_fbuf_rdr(rng_rdr);
    #[cfg(gpopt_trace="off")]
    let (c,r) = read_i32_pair_from_fbuf_rdr(rng_rdr);

    #[cfg(gpopt_trace="on")]
    println!("TP006:c={},r={}", c, r);

    TRACE_COUNT.store(c, Ordering::SeqCst);

    let d = (r as f64) / (RAND_MAX as f64);
    let f= (d * (max as f64)) + 0.5f64;
    let result = f.floor() as i16;
    result
}
    
#[cfg(gpopt_rng="file_stream")]
pub static TRACE_COUNT: AtomicI32 = AtomicI32::new(0);

// rnd - return random double value between 0.0 and 1.0
#[cfg(gpopt_rng="file_stream")]
fn rnd_dbl(rng_rdr: &mut BufReader<File>) -> f64 {
    #[cfg(gpopt_trace="on")]
    let (c,r) = read_i32_pair_from_fbuf_rdr(rng_rdr);

    #[cfg(gpopt_trace="off")]
    let (c,r) = read_i32_pair_from_fbuf_rdr(rng_rdr);

    TRACE_COUNT.store(c, Ordering::SeqCst);

    #[cfg(gpopt_trace="on")]
    println!("TP006:c={},r={}", c, r);

    let r_float = (r as f64)/(RAND_MAX as f64);
    r_float
}

#[allow(dead_code)]
#[cfg(gpopt_rng="file_stream")]
pub fn get_trace_count() -> i32 {
    TRACE_COUNT.load(Ordering::SeqCst)
}
