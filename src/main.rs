use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::io::{BufWriter};
use std::str;
use my_seqio::reader::{DynamicFastXReader, FastXReader};
use my_seqio::writer::DynamicFastXWriter;
use edit_distance;
use core::cmp::{min,max};
use std::collections::HashSet;

mod cli;
mod partial_suffix_sort;

fn extract_sequence<'a>(id: usize, seqs_concat: &'a Vec<u8>, cumul_seq_lengths: &'a Vec<usize>) -> &'a [u8]{
    let start = match id {
        0 => 0,
        _ => cumul_seq_lengths[id-1]
    };
    let end = cumul_seq_lengths[id];
    return &seqs_concat[start..end];
}

// End is one past the end
fn check_for_duplicates(doc_array: &Vec<usize>, cumul_seq_lengths: &Vec<usize>, seqs_concat: &Vec<u8>, already_tested: &mut HashSet<(usize,usize)>, SA_start: usize, SA_end: usize){
    // A hash set of integers
    if SA_end - SA_start > 100 {
        return; // Too many occurrences to check all pairs
    }

    for i in SA_start..SA_end{
        for j in i+1..SA_end{
            let id1 = doc_array[i];
            let id2 = doc_array[j];
            if id1 == id2{
                continue; // Same read
            }

            let s1 = extract_sequence(id1, seqs_concat, cumul_seq_lengths);
            let s2 = extract_sequence(id2, seqs_concat, cumul_seq_lengths);
            if !already_tested.contains(&(id1,id2)) {
                let result = edit_distance::edit_distance(str::from_utf8(s1).unwrap(), str::from_utf8(s2).unwrap());
                already_tested.insert((id1,id2));
                let d = result as f64; 
                if d < max(s1.len(), s2.len()) as f64 * 0.05 { // Less than 5% edit distance compared to the length of the longer sequence 
                    eprintln!("read {} (length {}) vs read {} (length {}): {}", id1, s1.len(), id2, s2.len(), result);
                }
            }
        }
    }

}

fn get_doc_array(partial_SA: &Vec<usize>, cumul_seq_lengths: &Vec<usize>) -> Vec<usize>{
    eprintln!("Inverting partial SA");
    let n = partial_SA.len();
    let mut inverse_partial_SA = vec![0usize; n];

    for i in 0..n{
        inverse_partial_SA[partial_SA[i]] = i;
    }

    eprintln!("Building doc array");

    let mut doc_array: Vec<usize> = vec![0; n];
    let mut seq_id = 0 as usize;
    for i in 0..n{
        while i == cumul_seq_lengths[seq_id]{
            seq_id += 1;
        }
        doc_array[inverse_partial_SA[i]] = seq_id; 
    }

    return doc_array;
}

fn main() {

    let matches = cli::build_cli().get_matches();

    let reads_in = matches.get_one::<String>("reads").unwrap();
    let k: usize = matches.get_one::<String>("seed-length").unwrap().parse::<usize>().unwrap();

    let mut input = DynamicFastXReader::new_from_file(&reads_in);
    let mut seqs_concat = Vec::<u8>::new();
    let mut cumul_seq_lengths = Vec::<usize>::new();
    let mut headers_concat = Vec::<u8>::new();
    let mut cumul_header_lengths = Vec::<usize>::new();
    while let Some(rec) = input.read_next(){
        let head = rec.head;
        let seq = rec.seq;

        // Append name to seq_concat
        seqs_concat.extend_from_slice(&seq);
        if cumul_seq_lengths.is_empty(){
            cumul_seq_lengths.push(seq.len());
        } else {
            cumul_seq_lengths.push(seq.len() + cumul_seq_lengths.last().unwrap());
        };

        // Append header to header_concat
        headers_concat.extend_from_slice(&head);
        if cumul_header_lengths.is_empty(){
            cumul_header_lengths.push(head.len());
        } else {
            cumul_header_lengths.push(head.len() + cumul_header_lengths.last().unwrap());
        };
    }
    let n = seqs_concat.len();

    eprintln!("Running partial suffix sort");
    let partial_SA = partial_suffix_sort::partial_suffix_sort(&seqs_concat, k);
    let doc_array = get_doc_array(&partial_SA, &cumul_seq_lengths);

    eprintln!("Finding k-mer runs");
    for _ in 0..k{
        seqs_concat.push(b'$'); // Padding so that k-mers can go over the end
    }
    let mut kmer_run_start = 0 as usize;
    let mut n_runs: usize = 0;
    let mut max_run: usize = 0;
    let mut already_tested = HashSet::<(usize,usize)>::new();
    for i in 0..n{
        let text_pos = partial_SA[i];
        if i > 0 {
            // Do we have the start of a new run?
            let prev_text_pos = partial_SA[i-1];
            if seqs_concat[text_pos..text_pos+k] != seqs_concat[prev_text_pos..prev_text_pos+k]{
                check_for_duplicates(&doc_array, &cumul_seq_lengths, &seqs_concat, &mut already_tested, kmer_run_start, i);
                n_runs += 1;
                max_run = max(max_run, i - kmer_run_start);
                kmer_run_start = i;
            }
        }
    }
    // Last run
    check_for_duplicates(&doc_array, &cumul_seq_lengths, &seqs_concat, &mut already_tested, kmer_run_start, n);
    n_runs += 1;
    max_run = max(max_run, n - kmer_run_start);

    eprintln!("Found {} k-mer runs", n_runs);
    eprintln!("Average run lengths: {}", n as f64 / n_runs as f64);
    eprintln!("Longest run is {} k-mers", max_run);

}

