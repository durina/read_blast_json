use clap::{Parser};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to Alignment file stored in fasta format
    #[arg(short='i', long="infile", required = true, action=clap::ArgAction::Append)]
    pub input_json: PathBuf,
    /// Path to store output
    #[arg(short='o', long="out-path", required = false, default_value_t=String::from("."))]
    pub out_path: String,
    /// Path to store output
    #[arg(short='p', long="prefix", required = false, default_value_t=String::from("out"))]
    pub prefix: String
}