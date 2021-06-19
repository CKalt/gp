#[cfg(gpopt_rng="FileStream")]
use util::read_i32_from_file_buf_rdr;
#[cfg(gpopt_rng="FileStream")]
use util::open_file_buf_rdr;

#[cfg(gpopt_rng="Seedable")]
use rand::{SeedableRng};

#[cfg(gpopt_rng="Thread")]
pub type GpRng = rand::rngs::ThreadRng;

#[cfg(gpopt_rng="Seedable")]
pub type GpRng = rand::rngs::StdRng;

#[cfg(gpopt_rng="FileStream")]
pub type GpRng = FileStreamRng ;

#[cfg(gpopt_rng="FileStream")]
pub struct FileStreamRng {
    fbuf_rdr: BufReader<File>,
}

#[cfg(gpopt_rng="FileStream")]
impl FileStreamRng {
    fn new() -> FileStreamRng {
        FileStreamRng {
            fbuf_rdr: open_fbuf_rdr("rng.log"),
        }
    }
    fn gen_range(&mut self, r: Range<i32>) -> i32 {
        let result = rnd(r.end as i16, self.fbuf_rdr);
        if (result < r.start) || (result >= r.end) {
            panic!("value ({}) from input stream out of range [{}..{}]",
                result, r.start, r.end);
        }
        result
    }
    fn gen<f32>(&mut self) -> f32 {
        rnd_dbl(self.fbuf_rdr)
    }
}

pub struct GpRngFactory {
}
impl GpRngFactory {
    #[cfg(gpopt_rng="Thread")]
    pub fn new() -> rand::rngs::ThreadRng {
        rand::thread_rng()
    }

    #[cfg(gpopt_rng="Seedable")]
    pub fn new() -> rand::rngs::StdRng {
        let seed: [u8; 32] = [151; 32];
        SeedableRng::from_seed(seed)
    }

    #[cfg(gpopt_rng="FileStream")]
    pub fn new() -> BufReader<File> {
        FileStreamRng::new()
    }
}

#[cfg(gpopt_rng="FileStream")]
use math::round;

#[cfg(gpopt_rng="FileStream")]
const RAND_MAX: i32 = 2147483647;    // defined in gcc <stdlib.h>

#[cfg(gpopt_rng="FileStream")]
fn rnd(max: i16, rng_rdr: &mut BufReader<File>) -> i16 {
    let r = read_i32_from_fbuf_rdr(rng_rdr);
    let d = (r as f64) / (RAND_MAX as f64);
    floor((d * (max as f64)) + 0.5f64, 0) as i16
}

// rnd - return random double value between 0.0 and 1.0
#[cfg(gpopt_rng="FileStream")]
fn rnd_dbl(rng_rdr: &mut BufReader<File>) -> f64 {
    (read_i32_from_fbuf_rdr(rng_rdr) as f64)/(RAND_MAX as f64)
}

#[cfg(gpopt_rng="FileStream")]
fn rnd_greedy_dbl_as_ll() -> GpInt {
    let r = rnd_dbl();
    let result: f64;

    if CONTROL.GRc < 0.0001 {
        result = r;
    } else if rnd(9) > 1 {
        // select from top set
        result = 1.0f64 - (CONTROL.GRc * r);
    } else {
        // select from lower set
        result = (1.0f64 - CONTROL.GRc)*r;
    }
    ll_result = float_to_int(result);

    return ll_result;
}

