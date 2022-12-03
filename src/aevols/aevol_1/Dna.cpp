//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <iostream>
#include <cassert>

Dna::Dna(int length, Threefry::Gen &&rng) {
    // Generate a random genome
    seq_ = boost::dynamic_bitset<>(length);
    // TODO: Kokkos
    //std::cout << "before for" << "Dna::Dna" << std::endl;
    for (int32_t i = 0; i < length; i++) {
        seq_[i] = rng.random(NB_BASE);
    }
    //std::cout << "after for" << "Dna::Dna" << std::endl;
    // print seq_
    // for (size_t i = 0; i < length; i++) {
    //     std::cout << seq_[i];
    // }
    // std::cout << std::endl;
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    //std::cout << "begin" << "Dna::save" << std::endl;
    int dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));

    // need to convert bitset to array of char   boost::dynamic_bitset::to_string(seq_).data()
    std::string dna_string;
    boost::to_string(seq_, dna_string);
    gzwrite(backup_file, dna_string.data(), dna_length * sizeof(seq_[0]));
    //std::cout << "end" << "Dna::save" << std::endl;
}

void Dna::load(gzFile backup_file) {
    //std::cout << "begin" << "Dna::load" << std::endl;
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_ = boost::dynamic_bitset<>(tmp_seq, tmp_seq + dna_length);
    //std::cout << "end" << "Dna::load" << std::endl;
}

void Dna::set(int pos, bool c) {
    seq_[pos] = c;
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
    //std::cout << "begin" << "Dna::remove" << std::endl;
    assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());
    boost::dynamic_bitset<> newseq_(seq_.size()-(pos_2-pos_1+1));
    // TODO: Kokkos
    for(size_t i = 0; i < seq_.size(); i++) {
        if(i < pos_1 || i > pos_2) {
            newseq_[i] = seq_[i];
        }
    }
    seq_ = newseq_;
    //std::cout << "end" << "Dna::remove" << std::endl;
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, boost::dynamic_bitset<>& insertedSeq) {
    //std::cout << "begin" << "Dna::insert" << std::endl;
    // Insert sequence 'seq' at position 'pos'
    assert(pos >= 0 && pos < seq_.size());

    boost::dynamic_bitset<> newseq_(seq_.size()+insertedSeq.size());
    // TODO: Kokkos
    for (size_t i = 0; i < seq_.size()+insertedSeq.size(); i++)
    {
        if(i < pos) {
            newseq_[i] = seq_[i];
        } else if(i >= pos && i < pos+insertedSeq.size()) {
            newseq_[i] = insertedSeq[i-pos];
        } else {
            newseq_[i] = seq_[i-insertedSeq.size()];
        }
    }
    seq_ = newseq_;
    //std::cout << "end" << "Dna::insert" << std::endl;
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna *dna) {
    Dna::insert(pos, dna->seq_);
}

void Dna::do_switch(int pos) {
    if (seq_[pos] == 0) seq_[pos] = 1;
    else seq_[pos] = 0;
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
    //std::cout << "begin" << "Dna::do_duplication" << std::endl;
    // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
    char *duplicate_segment = NULL;

    int32_t seg_length;

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
        boost::dynamic_bitset<> seq_dupl(pos_2 - pos_1 + 1);
        // TODO: Kokkos
        for(size_t i=pos_1; i < pos_2; ++i) {
            seq_dupl[i-pos_1] = seq_[i+pos_1];
        }
        insert(pos_3, seq_dupl);
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
        boost::dynamic_bitset<> seq_dupl(pos_2 - pos_1 + 1);
        // TODO: Kokkos
        for(size_t i=pos_1; i < seq_.size(); ++i) {
            seq_dupl[i-pos_1] = seq_[i];
        }
        size_t startShift = seq_.size()-pos_1;
        for(size_t i=0; i < pos_2; ++i) {
            seq_dupl[i+startShift] = seq_[i];
        }
        insert(pos_3, seq_dupl);
    }
    //std::cout << "end" << "Dna::do_duplication" << std::endl;
}

int Dna::promoter_at(int pos) {
    //std::cout << "begin" << "Dna::promoter_at" << std::endl;
    int dist_lead = 0;

    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        int search_pos = pos + motif_id;
        if (search_pos >= seq_.size()) {
            search_pos -= seq_.size();
        }
            
        // Searching for the promoter
        // do sum of PROM_SIZE first elements
        bool PROM_SEQvalue = PROM_SEQ[motif_id] == '0' ? false : true;
        dist_lead += (int)PROM_SEQ[motif_id] == (int)seq_[search_pos] ? 0 : 1;
    }

    //std::cout << "end" << "Dna::promoter_at" << std::endl;
    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) {
    //std::cout << "begin" << "Dna::terminator_at" << std::endl;
    int term_dist[TERM_STEM_SIZE];
    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;
        int left = pos + (TERM_SIZE - 1) - motif_id;

        // loop back the dna inf needed
        if (right >= length()) right -= length();
        if (left >= length()) left -= length();

        // Search for the terminators
        term_dist[motif_id] = seq_[right] != seq_[left] ? 1 : 0;
    }
    int dist_term_lead = term_dist[0] +
                         term_dist[1] +
                         term_dist[2] +
                         term_dist[3];

    //std::cout << "end" << "Dna::terminator_at" << std::endl;
    return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
    //std::cout << "begin" << "Dna::shine_dal_start" << std::endl;
    bool start = false;
    int t_pos, k_t;

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        t_pos = pos + k_t;
        if (t_pos >= seq_.size()) {
            t_pos -= seq_.size();
        }
        if ((int)seq_[t_pos] == (int)SHINE_DAL_SEQ[k_t]) {
            start = true;
        } else {
            start = false;
            break;
        }
    }

    //std::cout << "end" << "Dna::shine_dal_start" << std::endl;
    return start;
}

bool Dna::protein_stop(int pos) {
    //std::cout << "begin" << "Dna::protein_stop" << std::endl;
    bool is_protein;
    int t_k;

    for (int k = 0; k < CODON_SIZE; k++) {
        t_k = pos + k;
        if (t_k >= seq_.size())
            t_k -= seq_.size();

        if ((int)seq_[t_k] == (int)PROTEIN_END[k]) {
            is_protein = true;
        } else {
            is_protein = false;
            break;
        }
    }

    //std::cout << "end" << "Dna::protein_stop" << std::endl;
    return is_protein;
}

int Dna::codon_at(int pos) {
    //std::cout << "begin" << "Dna::codon_at" << std::endl;
    int value = 0;

    int t_pos;

    for (int i = 0; i < CODON_SIZE; i++) {
        t_pos = pos + i;
        if (t_pos >= seq_.size())
            t_pos -= seq_.size();
        if (seq_[t_pos] == 1)
            value += 1 << (CODON_SIZE - i - 1);
    }

    //std::cout << "end" << "Dna::codon_at" << std::endl;
    return value;
}