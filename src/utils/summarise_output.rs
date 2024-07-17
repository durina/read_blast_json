/*

    list as the hits as tables
    hit_match   taxid   title   evalue  coverage    query_start - query_match - query_end   query_strand hit_start - hit_match - hit_end    hit_strand  repeat_times

    hit_match
                title   accessionid     start   end     revised start   revised end
*/

// link: https://www.ncbi.nlm.nih.gov/nuccore/OV838944?report=fasta&log$=seqview&format=text&from=1243360&to=1243325&strand=true

use std::fs::File;
use std::io::{BufWriter, Write};
use crate::utils::hit_struct::BlastQuery;

pub fn summarise_output(parsed_output: Vec<BlastQuery>, outpath: String, prefix: String) {
    // write summary table
    for (x, query) in parsed_output.iter().enumerate() {
        let title = &query.query_title;
        let query_len = &query.query_len;
        let output_table_name = format!("{}/{}_output_table_{}_{}.tsv", outpath.strip_suffix('/').unwrap(), prefix, x+1, title);
        let mut blast_table_writer = BufWriter::new(File::create(output_table_name).unwrap());
        let headers = "\
        hit_seq\t\
        coverage\t\
        %identity\t\
        times_hit_seq_repeats\t\
        taxid\t\
        sequence_title\t\
        evalue\t\
        query_match_start\t\
        query_match_seq\t\
        query_match_end\t\
        query_match_strand\t\
        ".to_string();
        writeln!(blast_table_writer, "List of hit sequences that matched {}", title).unwrap();
        writeln!(blast_table_writer, "{}", headers).unwrap();
        let output_list_name = format!("{}/{}_seq_summary_{}_{}.tsv", outpath.strip_suffix('/').unwrap(), prefix, x+1, title);
        let mut blast_seq_list = BufWriter::new(File::create(output_list_name).unwrap());
        let headers = "\
        hit_seq\t\
        hit_accession\t\
        hit_actual_start\t\
        hit_actual_end\t\
        hit_NCBI_link\t\
        hit_strand\t\
        hit_revised_start\t\
        hit_revised_end\t\
        revised_hit_NCBI_link\t\
        hit_taxid
        ".to_string();
        writeln!(blast_seq_list, "Sequences that are similar to {}, grouped by taxonomy.", title).unwrap();
        writeln!(blast_seq_list, "{}", headers).unwrap();
        query.blast_results.iter().for_each( |hit_seq |
            {
                let (query_match_seq, query_match_start, query_match_end, query_match_strand, iden) = &hit_seq.query_extent;
                let query_strand = return_strand_char(query_match_strand);
                // write seq list
                writeln!(blast_table_writer,"{hitseq}\t{cov}\t{iden}\t{repeat}\t{taxid}\t{title}\t{evalue}\t{query_match_start}\t{query_match_seq}\t{query_match_end}\t{query_strand}",
                         hitseq=hit_seq.hitseq, cov=hit_seq.coverage, repeat=hit_seq.count,taxid=hit_seq.taxid.unwrap_or(1),
                         title=hit_seq.title, evalue=hit_seq.eval).unwrap();
                writeln!(blast_seq_list, "{}", hit_seq.hitseq).unwrap();
                hit_seq.target_extent
                    .iter()
                    .for_each(|(accession, (hit_start, hit_end, hit_strand))|
                        {
                            let hit_strand = return_strand_char(hit_strand);
                            let (revised_hit_start, revised_hit_end) = get_revised_coordinates(hit_start, hit_end, &hit_strand, query_match_start, query_match_end, query_len);
                            let hit_link = view_link_generator(accession, hit_start, hit_end, &hit_strand);
                            let revised_hit_link = view_link_generator(accession, &revised_hit_start, &revised_hit_end, &hit_strand);
                            writeln!(blast_seq_list, "\t{accession}\t{hit_start}\t{hit_end}\t{hit_link}\t{hit_strand}\t{revised_hit_start}\t{revised_hit_end}\t{revised_hit_link}\t{taxid}",
                                     taxid = hit_seq.taxid.unwrap_or(1)).unwrap();
                        }
                    );


            }
        );
    }

}

fn get_revised_coordinates(hit_start: &u64, hit_end: &u64, hit_strand: &char, query_start: &u64, query_end: &u64, query_len: &u64) -> (u64, u64) {
    let sequence_start = 1u64;
    let five_prime_offset = query_start - sequence_start;
    let three_prime_offset = query_len - query_end;
    let revised_hit_start: u64;
    let revised_hit_end: u64;
    if hit_strand == &'+' {
        revised_hit_start = if &five_prime_offset >= hit_start {
            println!("fp_offset {five_prime_offset} hit_start{hit_start}");
            1
        } else {
            hit_start - five_prime_offset
        };
        revised_hit_end = hit_end + three_prime_offset;
    } else if hit_strand == &'-' {
        revised_hit_start = hit_start + five_prime_offset;
        revised_hit_end = if &three_prime_offset >= hit_end {
            println!("3p_offset {three_prime_offset} {hit_end}");
            1
        } else {
            hit_end - three_prime_offset
        }
    } else {
        panic!("Strand information incomplete")
    }
    (revised_hit_start, revised_hit_end)
}

fn view_link_generator(accession: &str, from: &u64, to:&u64, strand: &char ) -> String {
    let k = if strand == &'+' {
        "off"
    } else {
        "on"
    };
    format!("https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&\
        id={accession}&from={from}&to={to}&strand={k}&extrafeat=null&conwithfeat=on&hide-cdd=on&ncbi_phid=CE8A5A904B610751000000000436035A")
}


fn return_strand_char(strand: &bool) -> char {
    if *strand {
        '+'
    } else {
        '-'
    }
}