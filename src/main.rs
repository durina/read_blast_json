/*
    Read BLAST json output from BLAST

    For Sensor project
    go through the json file and
        get subject sequence
        get taxid
        get coverage in the qseq
        get assiciated e-value
        get sseqid and sstart and send

    summarise the variations from the blast result
        sseq and the number of instances
            taxids associated with it
            corrected sstart and send

    Planned modules
        clap input
            get json file
            *get evalue cutoff
            *get coverage cutoff

        process json
        ouptut
 */
mod utils;
use utils::get_input::Cli;
use utils::process_json::receive_json;
use clap::Parser;
use crate::utils::summarise_output::summarise_output;

fn main() {
    let cli = Cli::parse();
    let parsed_output = receive_json(cli.input_json);
    summarise_output(parsed_output, cli.out_path, cli.prefix);
}