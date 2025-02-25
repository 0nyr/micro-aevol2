// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************

#pragma once

#include <memory>

#include "Abstract_ExpManager.h"
#include "Threefry.h"
#include "DnaMutator.h"
#include "Organism.h"
#include "Stats.h"

/**
 * Main class of the simulator.
 * ExpManager is in charge of running the simulation and maintaining all the data.
 * It is also that class that implements checkpointing and restore mechanisms.
 */
class ExpManager : public Abstract_ExpManager {
#ifdef USE_CUDA
    friend class cuExpManager;
#endif
public:
    ExpManager(int grid_height, int grid_width, int seed, double mutation_rate, int init_length_dna,
               int backup_step, int nb_host_threads);

    explicit ExpManager(int time, int nb_host_threads);

    ~ExpManager() override;

    void save(int t) const final;

    void load(int t) final;

    void run_evolution(int nb_gen) override;

private:
    __device__ void run_a_step(
        size_t indiv_id, 
        Kokkos::View<
            double*, 
            Kokkos::DefaultExecutionSpace::memory_space, 
            Kokkos::MemoryTraits<Kokkos::Atomic>
        > bestFitnessYet
    );

    __device__ void prepare_mutation(int indiv_id) const;

    __device__ void selection(int indiv_id) const;

    // kokkos CPU & GPU
    Kokkos::View<
        Organism, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > internal_organisms_;
    Kokkos::View<
        Organism, 
        Kokkos::DefaultExecutionSpace::memory_space
    > internal_organisms_gpu;

    Kokkos::View<
        Organism, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > prev_internal_organisms_;
    Kokkos::View<
        Organism, 
        Kokkos::DefaultExecutionSpace::memory_space
    > prev_internal_organisms_gpu;

    Kokkos::View<
        Organism, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > best_indiv;
    Kokkos::View<
        Organism, 
        Kokkos::DefaultExecutionSpace::memory_space
    > best_indiv_gpu;

    Kokkos::View<
        char*, 
        Kokkos::DefaultExecutionSpace::memory_space
    > DNA_seqs_gpu;
    Kokkos::View<
        char*, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > DNA_seqs;

    int *next_generation_reproducer_;

    Kokkos::View<
        DnaMutator, 
        Kokkos::DefaultExecutionSpace::memory_space
    > dna_mutator_array_gpu;
    Kokkos::View<
        DnaMutator, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > dna_mutator_array_;

    int nb_indivs_;

    int seed_;
    std::unique_ptr<Threefry> rng_;

    double *target;

    Stats *stats_best = nullptr;
    Stats *stats_mean = nullptr;

    int grid_height_;
    int grid_width_;

    double mutation_rate_;

    int backup_step_;
    int nb_host_threads_;
};
