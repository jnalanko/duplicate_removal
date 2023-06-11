use clap::{Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    
    Command::new("horizontal_coverage") // Todo: rename to just coverage
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .arg(
            Arg::new("reads")
                .short('r')
                .long("reads")
                .help("Input reads in fasta or fastq format, possibly gzipped")
                .global(true),
        )
        .arg(
            Arg::new("seed-length")
                .short('k')
                .long("seed-length")
                .help("The length of the alignment seeds")
                .global(true),
        )
}
