use serde_json::{Value};
use std::path::PathBuf;
use std::io::BufReader;
use std::fs::File;
use crate::utils::hit_struct::BlastQuery;

pub fn receive_json(json_file: PathBuf) -> Vec<BlastQuery> {
    let blast_file =  File::open(json_file).expect("Unable to open input json file");
    let blast_file_buffer = BufReader::new(blast_file);
    process_json(blast_file_buffer)
}

fn process_json(blast_file_buffer: BufReader<File>) -> Vec<BlastQuery> {
    // extract the items from the file
    let blast_data: Value = serde_json::from_reader(blast_file_buffer).expect("Not a valid json file");
    let collection_data = blast_data["BlastOutput2"].as_array().unwrap();
    let mut hit_collection: Vec<BlastQuery> = Vec::new();
    collection_data
        .iter()
        .for_each(| subject| {
            let hit_array = &subject["report"]["results"];
            hit_collection.push(BlastQuery::new(hit_array))
            }
        );
    hit_collection
}

