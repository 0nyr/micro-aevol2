//[
// Created by arrouan on 01/10/18.
//
#include "Dna.h"

#include <iostream>
#include <cassert>

Dna::Dna(
    Kokkos::View<
        char*, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > & DNA_seqs, 
    int length, 
    int seq_start,
    Threefry::Gen &&rng
) : DNA_seqs(DNA_seqs), seq_start(seq_start), 
    max_seq_length(length), seq_length(length) {
    // Generate a random genome
    for (int32_t i = seq_start; i < seq_start + length; i++) {
        DNA_seqs(i) = '0' + rng.random(NB_BASE);
    }
}

Dna::Dna(const Dna &clone, const size_t new_seq_start):
    DNA_seqs(clone.DNA_seqs), 
    max_seq_length(clone.max_seq_length),
    seq_length(clone.seq_length)
{
    seq_start = new_seq_start;
}

int Dna::length() const {
    return seq_length;
}

void Dna::save(gzFile backup_file) {
    // store seq_length, seq_start, max_seq_length
    gzwrite(backup_file, &seq_length, sizeof(seq_length));
    gzwrite(backup_file, &seq_start, sizeof(seq_start));
    gzwrite(backup_file, &max_seq_length, sizeof(max_seq_length));

    // copy the sequence from DNA_seqs to temporary seq
    char tmp_seq[seq_length];
    for (size_t i = seq_start; i < seq_start + seq_length; i++) {
        tmp_seq[i - seq_start] = DNA_seqs(i);
    }

    gzwrite(backup_file, tmp_seq, seq_length * sizeof(tmp_seq[0]));
}

void Dna::load(gzFile backup_file) {
    // read seq_length, seq_start, max_seq_length
    gzread(backup_file, &seq_length, sizeof(seq_length));
    gzread(backup_file, &seq_start, sizeof(seq_start));
    gzread(backup_file, &max_seq_length, sizeof(max_seq_length));

    char tmp_seq[seq_length];
    gzread(backup_file, tmp_seq, seq_length * sizeof(tmp_seq[0]));

    // copy the sequence from temporary seq_ to DNA_seqs
    for (int32_t i = 0; i < seq_length; i++) {
        DNA_seqs(seq_start + i) = tmp_seq[i];
    }
}

void Dna::set(int pos, char c) {
    DNA_seqs(seq_start + pos) = c;
}

/**
 * Remove the DNA inbetween pos_1 included and pos_2 included
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
    assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_length);
    // TODO: Kokkos
    for(size_t i = 0; pos_2+i+1 < seq_length; ++i) {
        DNA_seqs(pos_1+i) = DNA_seqs(pos_2+i+1);
    }
    // update real useful length
    seq_length = seq_length-(pos_2-pos_1+1);
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, char ** insertedSeq, size_t insertedSeqLength) {
    // Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_length);
    assert(pos + insertedSeqLength < seq_length);

    for (size_t i = 0; i < insertedSeqLength; i++) {
        DNA_seqs(seq_start + pos + i) = (*insertedSeq)[i];
    }
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna *dna) {
    assert(pos >= 0 && pos < seq_length);
    assert(pos + dna->seq_length < seq_length);

    for (size_t i = 0; i < dna->seq_length; i++) {
        DNA_seqs(seq_start + pos + i) = DNA_seqs(dna->seq_start + i);
    }
}

void Dna::do_switch(int pos) {
    if (DNA_seqs(seq_start + pos) == '0') DNA_seqs(seq_start + pos) = '1';
    else DNA_seqs(seq_start + pos) = '0';
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
    size_t seq_dupl_size = pos_2 - pos_1; // pos_2 is excluded
    char * seq_dupl[seq_dupl_size];

    if (pos_1 < pos_2) {
        //
        //       pos_1         pos_2                   -> 0-
        //         |             |                   -       -
        // 0--------------------------------->      -         -
        //         ===============                  -         - pos_1
        //           tmp (copy)                      -       -
        //                                             -----      |
        //                                             pos_2    <-'
        //
        // TODO: Kokkos
        for(size_t i = pos_1; i < pos_2; ++i) { // pos_2 is excluded
            *seq_dupl[i-pos_1] = DNA_seqs(seq_start + i);
        }
        insert(pos_3, seq_dupl, seq_dupl_size);
    } else { // if (pos_1 >= pos_2)
        // The segment to duplicate includes the origin of replication.
        // The copying process will be done in two steps.
        //
        //                                            ,->
        //    pos_2                 pos_1            |      -> 0-
        //      |                     |                   -       - pos_2
        // 0--------------------------------->     pos_1 -         -
        // ======                     =======            -         -
        //  tmp2                        tmp1              -       -
        //                                                  -----
        //
        //
        // TODO: Kokkos
        for(size_t i = pos_1; i < seq_length; ++i) { 
            *seq_dupl[i-pos_1] = DNA_seqs(seq_start + i);
        }
        size_t startShift = seq_length - pos_1;
        for(size_t i = 0; i < pos_2; ++i) { // pos_2 is excluded
            *seq_dupl[i+startShift] = DNA_seqs(seq_start + i);
        }
        insert(pos_3, seq_dupl, seq_dupl_size);
    }
}

int Dna::promoter_at(int pos) {
    //std::cout << "begin" << "Dna::promoter_at" << std::endl;
    int dist_lead = 0;

    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        int search_pos = pos + motif_id;
        if (search_pos >= seq_length) {
            search_pos -= seq_length;
        }
            
        // Searching for the promoter
        // do sum of PROM_SIZE first elements
        dist_lead += PROM_SEQ[motif_id] == DNA_seqs(seq_start + search_pos) ? 0 : 1;
    }

    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
    int dist_term_lead = 0;
    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;
        int left = pos + (TERM_SIZE - 1) - motif_id;

        // loop back the dna inf needed
        if (right >= length()) right -= length();
        if (left >= length()) left -= length();

        // Search for the terminators
        dist_term_lead += DNA_seqs(seq_start + right) != DNA_seqs(seq_start + left) ? 1 : 0;
    }
    return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
    bool start = false;
    int t_pos, k_t;

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        t_pos = pos + k_t;
        if (t_pos >= seq_length) {
            t_pos -= seq_length;
        }
        if (DNA_seqs(seq_start + t_pos) == SHINE_DAL_SEQ[k_t]) {
            start = true;
        } else {
            start = false;
            break;
        }
    }

    return start;
}

bool Dna::protein_stop(int pos) {
    bool is_protein;
    int t_k;

    for (int k = 0; k < CODON_SIZE; k++) {
        t_k = pos + k;
        if (t_k >= seq_length)
            t_k -= seq_length;

        if (DNA_seqs(seq_start + t_k) == PROTEIN_END[k]) {
            is_protein = true;
        } else {
            is_protein = false;
            break;
        }
    }

    return is_protein;
}

int Dna::codon_at(int pos) {
    int value = 0;

    int t_pos;

    for (int i = 0; i < CODON_SIZE; i++) {
        t_pos = pos + i;
        if (t_pos >= seq_length)
            t_pos -= seq_length;
        if (DNA_seqs(seq_start + t_pos) == '1')
            value += 1 << (CODON_SIZE - i - 1);
    }

    return value;
}