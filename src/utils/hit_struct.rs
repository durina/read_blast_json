// go through the json file and
// get subject sequence
// get taxid
// get coverage in the qseq
// get assiciated e-value
// get sseqid and sstart and send strand

use std::collections::HashMap;
use serde_json::Value;

pub struct BlastQuery {
    pub query_title: String,
    pub query_len: u64,
    pub blast_results: Vec<HitSeq>
}

pub struct HitSeq {
    pub taxid: Option<u64>,
    pub hitseq: String,
    pub coverage: f64,
    pub eval: f64,
    pub query_extent: (String, u64, u64, bool, f64),
    pub target_extent: HashMap<String, (u64, u64, bool)>,
    pub gaps: u64,
    pub midline: String,
    pub count: u64,
    pub title: String
}

impl BlastQuery {
    pub fn new(hit_array: &Value) -> BlastQuery {
        // hit_array is json file ["BlastOutput2"][n]["report"]["results"]
        let mut results: Vec<HitSeq> = Vec::new();
        let query_len = hit_array["search"]["query_len"].as_f64().unwrap();
        hit_array["search"]["hits"].as_array().unwrap()
            .iter()
            .for_each(|hit| {
                extract_hit_data(hit, &mut results, &query_len)
        });
        BlastQuery {
            query_title: hit_array["search"]["query_title"].as_str().unwrap().to_string(),
            query_len: hit_array["search"]["query_len"].as_u64().unwrap(),
            blast_results: results
        }
    }
}

fn extract_hit_data(hit: &Value, hit_seq_vec: &mut Vec<HitSeq>, query_len: &f64) {
    let hit_seq_meta = &hit["description"][0];
    let hit_seq_info = &hit["hsps"][0];
    let this_taxid = hit_seq_meta["taxid"].as_u64();
    let this_hitseq = hit_seq_info["hseq"].as_str().unwrap().to_string();
    let hit = hit_seq_vec.iter()
        .position(|hitseq| hitseq.taxid == this_taxid && hitseq.hitseq == this_hitseq);
    match hit {
        Some(x) => {
            hit_seq_vec[x].update(hit_seq_info, hit_seq_meta);
        }
        None => {
            let new_hit_seq = HitSeq::new_from(hit_seq_info, hit_seq_meta, query_len);
            hit_seq_vec.push(new_hit_seq);
        }
    }
}


impl HitSeq {
    pub fn new_from(hit_info: &Value, hit_meta: &Value, qlen: &f64) -> HitSeq {
        let query_start = hit_info["query_from"].as_u64().unwrap();
        let query_end = hit_info["query_to"].as_u64().unwrap();
        let query_match_seq = hit_info["qseq"].as_str().unwrap().to_string();
        let query_strand = get_strand_info(hit_info["query_strand"].as_str().unwrap());
        let identity = hit_info["identity"].as_f64().unwrap();
        let align_len = hit_info["align_len"].as_f64().unwrap();

        let coverage = hit_info["align_len"].as_f64().unwrap()/qlen*100.0;
        let evalue = hit_info["evalue"].as_f64().unwrap();

        let hit_start = hit_info["hit_from"].as_u64().unwrap();
        let hit_end = hit_info["hit_to"].as_u64().unwrap();
        let hit_strand = get_strand_info(hit_info["hit_strand"].as_str().unwrap());
        let hit_seqid = hit_meta["accession"].as_str().unwrap().to_string();
        let hit_title = hit_meta["title"].as_str().unwrap().to_string();

        let taxid = hit_meta["taxid"].as_u64();
        let hitseq = hit_info["hseq"].as_str().unwrap().to_string();

        let mut target_map: HashMap<String, (u64, u64, bool)> = HashMap::new();
        target_map.insert(hit_seqid, (hit_start, hit_end, hit_strand));

        let gaps = hit_info["gaps"].as_u64().unwrap();
        let midline = hit_info["midline"].as_str().unwrap();

        let percent_identity = identity/align_len*100.0;

        HitSeq {
            taxid,
            hitseq,
            coverage,
            eval: evalue,
            query_extent: (query_match_seq, query_start, query_end, query_strand, percent_identity),
            target_extent: target_map,
            gaps,
            midline: midline.to_string(),
            count: 1,
            title: hit_title
        }
    }

    fn update(&mut self, hit_info: &Value, hit_meta: &Value) {
        let hit_start = hit_info["hit_from"].as_u64().unwrap();
        let hit_end = hit_info["hit_to"].as_u64().unwrap();
        let hit_strand = get_strand_info(hit_info["hit_strand"].as_str().unwrap());
        let hit_seqid = hit_meta["accession"].as_str().unwrap().to_string();
        self.count += 1;
        self.target_extent.entry(hit_seqid).or_insert((hit_start, hit_end, hit_strand));
    }

}

fn get_strand_info(strand: &str) -> bool {
    if strand == "Plus" {
        true
    } else if strand == "Minus"{
        false
    } else {
        panic!()
    }
}