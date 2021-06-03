#![allow(dead_code)]
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::fs;
use std::path::Path;

pub fn append_line_to_file(fname: &str, text: &str) {
    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(fname)
        .unwrap();

    if let Err(e) = writeln!(file, "{}", text) {
        eprintln!("Couldn't write to file ({}): {}", fname, e);
    }
}

pub fn remove_file_if_exists(fname: &str) {
    if Path::new(fname).exists() {
        fs::remove_file(fname).expect("remove_file failed");
    }
}
