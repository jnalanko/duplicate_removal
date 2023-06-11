use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::io::{BufWriter};
use std::str;
use my_seqio::reader::{DynamicFastXReader, FastXReader};
use my_seqio::writer::DynamicFastXWriter;
use core::cmp::{min,max};

mod cli;
mod partial_suffix_sort;

fn main() {

    let matches = cli::build_cli().get_matches();

    let reads_in = matches.get_one::<String>("reads").unwrap();
    let k: usize = matches.get_one::<String>("seed-length").unwrap().parse::<usize>().unwrap();

    let mut input = DynamicFastXReader::new_from_file(&reads_in);
    let mut seq_concat = Vec::<u8>::new();
    let mut seq_lengths = Vec::<usize>::new();
    let mut headers_concat = Vec::<u8>::new();
    let mut header_lengths = Vec::<usize>::new();    
    while let Some(rec) = input.read_next(){
        let head = rec.head;
        let seq = rec.seq;

        // Append name to seq_concat
        seq_concat.extend_from_slice(&seq);
        seq_lengths.push(seq.len());

        // Append header to header_concat
        headers_concat.extend_from_slice(&head);
        header_lengths.push(head.len());
    }

    eprintln!("Running partial suffix sort");
    let partial_SA = partial_suffix_sort::partial_suffix_sort(&seq_concat, k);
    let n = partial_SA.len();

    eprintln!("Inverting partial SA");
    let mut inverse_partial_SA = vec![0usize; n];

    for i in 0..n{
        inverse_partial_SA[partial_SA[i]] = i;
    }

    eprintln!("Finding k-mer runs");
    for _ in 0..k{
        seq_concat.push(b'$'); // Padding so that k-mers can go over the end
    }
    let mut kmer_run_start = 0 as usize;
    let mut n_runs: usize = 0;
    let mut max_run: usize = 0;
    for i in 0..n{
        let text_pos = partial_SA[i];
        if i > 0 {
            // Do we have the start of a new run?
            let prev_text_pos = partial_SA[i-1];
            if seq_concat[text_pos..text_pos+k] != seq_concat[prev_text_pos..prev_text_pos+k]{
                n_runs += 1;
                max_run = max(max_run, i - kmer_run_start);
                kmer_run_start = i;
            }
        }
    }
    // Last run
    n_runs += 1;
    max_run = max(max_run, n - kmer_run_start);

    eprintln!("Found {} k-mer runs", n_runs);
    eprintln!("Average run lengths: {}", n as f64 / n_runs as f64);
    eprintln!("Longest run is {} k-mers", max_run);

}

