use clap::{ArgGroup, Parser};
use opt_bbwt::{
    lyndon::{lyndon_suffix_table, lyndon_table},
    optimal_bbwt::optimal_bbwt,
    string::rotate_string,
    suffix_array::{circular_suffix_array, sa_to_rank, suffix_array},
    transform::{bbwt::bbwt, bwt::bwt},
};

#[derive(Debug, Parser)]
#[command(
    about = "Run string algorithms on a Pizza&Chili corpus or a local file",
    group = ArgGroup::new("input").args(["corpus", "file"])
)]
struct Args {
    /// Operation to run: bwt, bbwt, optimal_bbwt, lyndon, lyndon_suffix,
    /// suffix_array, or circular_suffix_array.
    #[arg(short, long)]
    action: String,

    /// Pizza&Chili corpus name.
    #[arg(short, long)]
    corpus: Option<String>,

    /// Local input file. Cannot be combined with --corpus.
    #[arg(short, long)]
    file: Option<std::path::PathBuf>,

    /// Maximum number of bytes to read from the corpus.
    #[arg(short, long, default_value_t = 50_000)]
    length: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let text = match (&args.corpus, &args.file) {
        (Some(corpus), None) => opt_bbwt::file::load_corpus(corpus, args.length)?,
        (None, Some(file)) => opt_bbwt::file::load_file_prefix(
            file.to_str().ok_or("input file path must be valid UTF-8")?,
            args.length,
        )?,
        (None, None) => opt_bbwt::file::load_corpus("dna", args.length)?,
        (Some(_), Some(_)) => unreachable!("clap prevents combining corpus and file"),
    };
    if text.is_empty() {
        return Err("the selected corpus is empty".into());
    }

    match args.action.as_str() {
        "bwt" => {
            let transformed = bwt(&text);
            println!("BWT: {:?}", transformed);
        }
        "bbwt" => {
            let transformed = bbwt(&text);
            println!("BBWT: {:?}", transformed);
        }
        "optimal_bbwt" => {
            let rotation_index = optimal_bbwt(&text);
            let rotated = rotate_string(&text, rotation_index);
            println!("Optimal BBWT rotation index: {}", rotation_index);
            println!("Rotated string: {:?}", rotated);
        }
        "lyndon" => {
            let suffix_array = suffix_array(&text);
            let rank = sa_to_rank(&suffix_array);
            let lyndon_table = lyndon_table(&text, &rank);
            println!("Lyndon table: {:?}", lyndon_table);
        }
        "lyndon_suffix" => {
            let lyndon_suffix_table = lyndon_suffix_table(&text);
            println!("Lyndon suffix table: {:?}", lyndon_suffix_table);
        }
        "suffix_array" => {
            let suffix_array = suffix_array(&text);
            println!("Suffix array: {:?}", suffix_array);
        }
        "circular_suffix_array" => {
            let csa = circular_suffix_array(&text);
            println!("Circular suffix array: {:?}", csa);
        }
        _ => {
            return Err(format!("Unknown action: {}", args.action).into());
        }
    }

    Ok(())
}
