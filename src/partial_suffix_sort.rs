fn get_cumulative(v: &Vec<u64>) -> Vec<u64>{
    let mut C = vec![0u64; v.len()];
    for i in 1..v.len(){
        C[i] = v[i-1] + C[i-1];
    }
    return C;
}

#[allow(non_snake_case)]
pub fn partial_suffix_sort(S: &Vec<u8>, mut k: usize) -> Vec<usize>{

    if k > S.len(){
        k = S.len();
    }

    // Create a vector with elements 0,1,2...
    let mut order: Vec<usize> = (0..S.len()).collect();

    // Find out the counts of charactes considered on the first round
    let mut counts: Vec<u64> = vec![0u64; 256];
    for i in k-1..S.len(){
        counts[S[i] as usize] += 1;
    }

    counts[0] += k as u64 - 1; // Past the end padding charaters

    // k rounds of LSB radix sort
    for round in 0..k{
        let mut C = get_cumulative(&counts);
        // Sort the order vector by the k-mer starting at position i
        let mut new_order: Vec<usize> = vec![0; S.len()];
        for i in order{
            let mut c = 0 as u8;
            if i + k - 1 - round < S.len(){
                c = S[i+k-1-round];
            }
            new_order[C[c as usize] as usize] = i;
            C[c as usize] += 1;
        }
        order = new_order;

        if round < k-1{
            // Remove one padding characters from the counts and add one real character
            counts[0] -= 1;
            counts[S[k-1-round-1] as usize] += 1;
        }
    }

    return order;
 
}

// Unit test
#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;

    #[test]
    fn test_partial_suffix_sort_banana(){
        let S: Vec<u8> = b"banana$".to_vec();	
        let k = 10;
        let order = partial_suffix_sort(&S, k);
        assert_eq!(order, vec![6, 5, 3, 1, 0, 4, 2]);
    }

    #[test]
    fn test_partial_suffix_sort_vs_rust_std_sort(){
        let S: Vec<u8> = b"ababbbabababbbababbbaababaaabaababbbaababbabababbbabbaabbabababbababab".to_vec();
        let k = 3;
        let n = S.len();
        let mut suffixes = Vec::<Vec<u8>>::new();
        for i in 0..S.len(){
            suffixes.push(S[i..].to_vec());
        }

        suffixes.sort();
        let partial_SA = partial_suffix_sort(&S, k);

        for i in 0..suffixes.len(){
            let A = &suffixes[i][0..min(k, suffixes[i].len())];
            let B = &S[partial_SA[i]..min(partial_SA[i]+k, n)];
            assert_eq!(A,B);
        }
    }
}