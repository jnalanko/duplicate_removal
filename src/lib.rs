/*use std::fs::File;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::str;
use my_seqio::reader::{DynamicFastXReader, FastXReader};
use my_seqio::writer::DynamicFastXWriter;

pub fn parse_fasta(filename: &String) -> BTreeMap::<String, Vec<u64>>{
    let mut input = DynamicFastXReader::new_from_file(&filename);

    // A hash table with string keys and Vec<u8> values
    let mut seq_to_coverage = BTreeMap::<String, Vec<u64>>::new();
    while let Some(mut rec) = input.read_next(){
        let name = str::from_utf8(rec.head).unwrap();

        // Split name by whitespace and take first token
        let name = name.split_whitespace().next().unwrap();
        let len = rec.seq.len();

        // Create a vector of u8 of length len
        let coverage = vec![0u64; len];
        seq_to_coverage.insert(name.to_string(), coverage);
    }
    return seq_to_coverage;
}
*/