#[cfg(gpopt_rng="Seedable")]
use rand::{SeedableRng};

#[cfg(gpopt_rng="Thread")]
pub type GpRng = rand::rngs::ThreadRng;

#[cfg(gpopt_rng="Seedable")]
pub type GpRng = rand::rngs::StdRng;

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
}
