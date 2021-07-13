pub type GpHits = u8;
pub type GpRaw = u8;
pub type GpStandardized = u8;
pub type GpFloat = f32;
pub type GpFitness = GpFloat;

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

impl Fitness {
    pub fn new() -> Fitness {
        Fitness {
            nfr: 0.0,
            n:   0.0,
            a:   0.0,
            raw: 0.0,

            r: 0,
            s: 0,
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
    pub fn clone(&self) -> Fitness {
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

