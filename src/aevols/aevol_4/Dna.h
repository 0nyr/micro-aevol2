//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <bitset>
#include <zlib.h>
#include <boost/dynamic_bitset.hpp>

#include "Threefry.h"
#include "aevol_constants.h"

class Dna {

public:
    Dna(
        boost::dynamic_bitset<> & DNA_seqs
    ): DNA_seqs(DNA_seqs) {};

    Dna(const Dna &clone, const size_t new_seq_start);

    Dna(
        boost::dynamic_bitset<> & DNA_seqs, 
        int length, 
        int seq_start, 
        Threefry::Gen &&rng
    );

    ~Dna() = default;

    int length() const;

    void save(gzFile backup_file);

    void load(gzFile backup_file);

    void set(int pos, bool c);

    /// Remove the DNA inbetween pos_1 and pos_2
    void remove(int pos_1, int pos_2);

    /// Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, boost::dynamic_bitset<>& insertedSeq);

    /// Insert a sequence of a given length at a given position into the DNA of the Organism
    void insert(int pos, Dna *seq);

    void do_switch(int pos);

    void do_duplication(int pos_1, int pos_2, int pos_3);

    int promoter_at(int pos);

    int terminator_at(int pos);

    bool shine_dal_start(int pos);

    bool protein_stop(int pos);

    int codon_at(int pos);

    boost::dynamic_bitset<> & DNA_seqs;
    size_t seq_start;
    size_t max_seq_length;
    size_t seq_length;

};
