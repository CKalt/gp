#![allow(dead_code)]
use std::fs::OpenOptions;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufReader;
use std::io::BufRead;
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

pub fn write_i32_to_file(fname: &str, num: i32) {
    let mut file = File::create(fname).unwrap();
    if let Err(e) = writeln!(file, "{}", num) {
        eprintln!("Couldn't write to file ({}): {}", fname, e);
    }
}

pub fn read_i32_from_file(fname: &str) -> i32 {
    let f = File::open(fname).unwrap();
    let mut file = BufReader::new(&f);
    let mut line = String::new();
    let len = file.read_line(&mut line).unwrap();
    if len < 1 {
        panic!("empty string from file.");
    }
    line.truncate(line.len() - 1);
    let count = line.parse::<i32>().unwrap();
    count
}

pub fn inc_file_counter(fname: &str) -> i32 {
    let count = read_i32_from_file(fname);
    let new_count = count+1;
    write_i32_to_file(fname, new_count);
    new_count
}

pub fn reset_file_counter(fname: &str) {
    remove_file_if_exists(fname);
    write_i32_to_file(fname, 0);
}

