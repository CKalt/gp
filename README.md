# gp
Configuration
    The file .cargo/config.toml is used to specify conditional compile options that 
    are passed to the rustc --cfg command
    
    gpopt_trail="santa_fe" | "los_altos"
        Default: "santa_fe"

Example contents for .cargo/config.toml

[build]
#rustflags = "--cfg gpopt_trail=\"santa_fe\""
rustflags = "--cfg gpopt_trail=\"los_altos\""


Building
    $ cargo build [--verbose]

Running
    $ cargo run



# Testing for ant run against gp.c

