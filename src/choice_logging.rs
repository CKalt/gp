#![cfg(gpopt_choice_logging="write")]
use std::convert::TryInto;
use std::sync::atomic::{AtomicUsize, Ordering};
use crate::util::*;
    
static LOG_COUNT: AtomicUsize = AtomicUsize::new(0);
const CHOICE_LOG_FNAME: &str = "choice.log";

pub fn init_choice_log() {
    remove_file_if_exists(CHOICE_LOG_FNAME);
}

pub fn inc_log_counter() -> i32 {
    let old_count = LOG_COUNT.fetch_add(1, Ordering::SeqCst);

    // convert from usize, return if able, panic otherwise.
    (old_count+1).try_into().unwrap()
}

#[cfg(gpopt_trace="on")]
pub fn get_log_count() ->i32 {
    LOG_COUNT.load(Ordering::Relaxed) as i32
}

pub fn choice_log(tp: u8, choice_value: &str) -> i32 {
    let counter = inc_log_counter();
    let msg = format!("{},{},{}", counter, tp, choice_value);
    append_line_to_file(CHOICE_LOG_FNAME, &msg);
    counter
}
