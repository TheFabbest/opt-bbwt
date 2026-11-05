use clap::Parser;
use opt_bbwt::{
    bbwt::bbwt,
    bwt::bwt,
    lyndon::{lyndon_suffix_table, lyndon_table},
    optimalbbwt::optimal_bbwt_chain,
    rank::sa_to_rank,
    saca::{circular_suffix_array, sais},
    string::rotate_string,
};

#[derive(Debug, Parser)]
#[command(about = "Run BBWT and heuristic measurements on a corpus")]
struct Args {
    /// What to do
    #[arg(short, long)]
    action: String,

    /// Name of the file under ./corpora.
    #[arg(short, long, default_value = "dna")]
    corpus: String,

    /// Maximum number of bytes to read from the corpus.
    #[arg(short, long, default_value_t = 50_000)]
    length: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let text = opt_bbwt::file::load_corpus(&args.corpus, args.length)?;
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
            let rotation_index = optimal_bbwt_chain(&text, true);
            let rotated = rotate_string(&text, rotation_index);
            println!("Optimal BBWT rotation index: {}", rotation_index);
            println!("Rotated string: {:?}", rotated);
        }
        "lyndon" => {
            let suffix_array = sais(&text);
            let rank = sa_to_rank(&suffix_array);
            let lyndon_table = lyndon_table(&text, &rank);
            println!("Lyndon table: {:?}", lyndon_table);
        }
        "lyndon_suffix" => {
            let lyndon_suffix_table = lyndon_suffix_table(&text);
            println!("Lyndon suffix table: {:?}", lyndon_suffix_table);
        }
        "suffix_array" => {
            let suffix_array = sais(&text);
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
