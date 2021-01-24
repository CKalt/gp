struct Tree {

}

fn run() -> Result <Tree, String> {


}

fn main() {
    let winner = loop {
        let result = run();

        match result {
            Ok(winner) => break winner,
            Err(e) => println!("no winner: {}", e),
        }
    };
}
