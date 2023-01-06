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


#include <iostream>
#include <zlib.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>

using namespace std;

#include "ExpManager.h"
#include "AeTime.h"
#include "Gaussian.h"

// For time tracing
#include "Timetracer.h"

/**
 * Constructor for initializing a new simulation
 *
 * @param grid_height : Height of the grid containing the organisms
 * @param grid_width : Width of the grid containing the organisms
 * @param seed : Global seed for all the PRNG of the simulation
 * @param mutation_rate : Mutation rates for all the organism during the simulation
 * @param init_length_dna : Size of the randomly generated DNA at the initialization of the simulation
 * @param selection_pressure : Selection pressure used during the selection process
 * @param backup_step : How much often checkpoint must be done
 */
ExpManager::ExpManager(int grid_height, int grid_width, int seed, double mutation_rate, int init_length_dna,
                       int backup_step, int nb_host_threads)
        : seed_(seed), rng_(new Threefry(grid_width, grid_height, seed)), nb_host_threads_(nb_host_threads) {
    // Initializing the data structure
    grid_height_ = grid_height;
    grid_width_ = grid_width;

    backup_step_ = backup_step;

    nb_indivs_ = grid_height * grid_width;

    internal_organisms_gpu = Kokkos::View<
            Organism, 
            Kokkos::DefaultExecutionSpace::memory_space
        > ("internal_organisms_gpu", nb_indivs_);
    Kokkos::View<
        Organism, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > internal_organisms_ = Kokkos::create_mirror_view(internal_organisms_gpu);

    prev_internal_organisms_gpu = Kokkos::View<
            Organism, 
            Kokkos::DefaultExecutionSpace::memory_space
        > ("prev_internal_organisms_gpu", nb_indivs_);
    Kokkos::View<
        Organism, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > prev_internal_organisms_ = Kokkos::create_mirror_view(prev_internal_organisms_gpu);


    next_generation_reproducer_ = new int[nb_indivs_]();
    dna_mutator_array_ = new DnaMutator *[nb_indivs_];

    mutation_rate_ = mutation_rate;

    DNA_seqs_gpu = Kokkos::View<
            char*, 
            Kokkos::DefaultExecutionSpace::memory_space
        > ("DNA_seqs", nb_indivs_ * init_length_dna);
    Kokkos::View<
        char*, 
        Kokkos::DefaultHostExecutionSpace::memory_space
    > DNA_seqs = Kokkos::create_mirror_view(DNA_seqs_gpu);

    // Building the target environment
    // OPTI: Remove useless malloc
    Gaussian g1 (1.2, 0.52, 0.12);
    Gaussian g2 (-1.4, 0.5, 0.07);
    Gaussian g3 (0.3, 0.8, 0.03);

    target = new double[FUZZY_SAMPLING];
    double geometric_area = 0.0;
    for (int i = 0; i < FUZZY_SAMPLING; i++) {
        double pt_i = ((double) i) / (double) FUZZY_SAMPLING;

        double tmp = g1.compute_y(pt_i);
        tmp += g2.compute_y(pt_i);
        tmp += g3.compute_y(pt_i);

        tmp = tmp > Y_MAX ? Y_MAX : tmp;
        tmp = tmp < Y_MIN ? Y_MIN : tmp;

        target[i] = tmp;
        geometric_area += tmp / (double)FUZZY_SAMPLING;
    }

    printf("Initialized environmental target %f\n", geometric_area);


    // Initializing the PRNGs
    for (int indiv_id = 0; indiv_id < nb_indivs_; ++indiv_id) {
        dna_mutator_array_[indiv_id] = nullptr;
    }

    // Generate a random organism that is better than nothing
    Kokkos::View<
        bool, 
        Kokkos::DefaultHostExecutionSpace::memory_space, 
        Kokkos::MemoryTraits<Kokkos::Atomic>
    > found_organism("found_organism"); // kokkos atomic bool
    found_organism() = false;

    // kokkos parallel for while found_organism is false
    Kokkos::parallel_for(
        "ExpManager::ExpManager find organism",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nb_host_threads_),
        [=](const size_t i) 
    {
        while(!found_organism()) {
            Organism random_organism(
                DNA_seqs, // pass by reference the pointed object
                0,
                init_length_dna,
                rng_->gen(0, Threefry::MUTATION)
            );
            random_organism.locate_promoters();
            random_organism.evaluate(target);
            double r_compare = round((random_organism.metaerror - geometric_area) * 1E10) / 1E10;

            if (r_compare < 0.0 && found_organism() == false) {
                // good random organism found
                    internal_organisms_(0) = random_organism;
                    found_organism() = true;
            }
        }
    });

    printf("Populating the environment\n");

    // Create a population of clones based on the randomly generated organism
    // kokkos parallel for: create nb_indivs_/NB_THREADS clones
    Kokkos::parallel_for(
        "ExpManager::ExpManager populate",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
            0, nb_host_threads_
        ),
        [=](const int thread_id) 
    {
        size_t start = thread_id * nb_indivs_ / nb_host_threads_;
        size_t end = (thread_id + 1) * nb_indivs_ / nb_host_threads_;
        for (int indiv_id = start; indiv_id < end; ++indiv_id) {
            prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id] =
                std::make_shared<Organism>(internal_organisms_[0], indiv_id);
        }
    });

    // Create backup and stats directory
    create_directory();

    // copy the data from the host to the device
    Kokkos::deep_copy(*DNA_seqs_gpu, *DNA_seqs);
}

/**
 * Constructor to resume/restore a simulation from a given backup/checkpoint file
 *
 * @param time : resume from this generation
 */
ExpManager::ExpManager(int time, int nb_host_threads) : nb_host_threads_(nb_host_threads) {
    target = new double[FUZZY_SAMPLING];

    load(time);

    double geometric_area = 0;
    for (int i = 0; i < FUZZY_SAMPLING - 1; i++) {
        // Computing a trapezoid area
        geometric_area += ((fabs(target[i]) + fabs(target[i + 1])) / (2 * (double) FUZZY_SAMPLING));
    }

    printf("Initialized environmental target %f\n", geometric_area);

    dna_mutator_array_ = new DnaMutator *[nb_indivs_];
    for (int indiv_id = 0; indiv_id < nb_indivs_; ++indiv_id) {
        dna_mutator_array_[indiv_id] = nullptr;
    }

}

/**
 * Checkpointing/Backup of the population of organisms
 *
 * @param t : simulated time of the checkpoint
 */
void ExpManager::save(int t) const {

    // get data back from the device
    Kokkos::deep_copy(*DNA_seqs, *DNA_seqs_gpu);

    char exp_backup_file_name[255];

    sprintf(exp_backup_file_name, "backup/backup_%d.zae", t);

    // -------------------------------------------------------------------------
    // Open backup files
    // -------------------------------------------------------------------------
    gzFile exp_backup_file = gzopen(exp_backup_file_name, "w");


    // -------------------------------------------------------------------------
    // Check that files were correctly opened
    // -------------------------------------------------------------------------
    if (exp_backup_file == Z_NULL) {
        printf("Error: could not open backup file %s\n",
               exp_backup_file_name);
        exit(EXIT_FAILURE);
    }


    // -------------------------------------------------------------------------
    // Write the backup file
    // -------------------------------------------------------------------------
    gzwrite(exp_backup_file, &t, sizeof(t));

    gzwrite(exp_backup_file, &grid_height_, sizeof(grid_height_));
    gzwrite(exp_backup_file, &grid_width_, sizeof(grid_width_));

    gzwrite(exp_backup_file, &backup_step_, sizeof(backup_step_));

    gzwrite(exp_backup_file, &mutation_rate_, sizeof(mutation_rate_));

    for (int i = 0; i < FUZZY_SAMPLING; i++) {
        double tmp = target[i];
        gzwrite(exp_backup_file, &tmp, sizeof(tmp));
    }

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        prev_internal_organisms_[indiv_id]->save(exp_backup_file);
    }

    rng_->save(exp_backup_file);

    if (gzclose(exp_backup_file) != Z_OK) {
        cerr << "Error while closing backup file" << endl;
    }
}

/**
 * Loading a simulation from a checkpoint/backup file
 *
 * @param t : resuming the simulation at this generation
 */
void ExpManager::load(int t) {

    char exp_backup_file_name[255];

    sprintf(exp_backup_file_name, "backup/backup_%d.zae", t);

    // -------------------------------------------------------------------------
    // Open backup files
    // -------------------------------------------------------------------------
    gzFile exp_backup_file = gzopen(exp_backup_file_name, "r");


    // -------------------------------------------------------------------------
    // Check that files were correctly opened
    // -------------------------------------------------------------------------
    if (exp_backup_file == Z_NULL) {
        printf("Error: could not open backup file %s\n",
               exp_backup_file_name);
        exit(EXIT_FAILURE);
    }


    // -------------------------------------------------------------------------
    // Write the backup file
    // -------------------------------------------------------------------------
    int time;
    gzread(exp_backup_file, &time, sizeof(time));
    AeTime::set_time(time);

    gzread(exp_backup_file, &grid_height_, sizeof(grid_height_));

    gzread(exp_backup_file, &grid_width_, sizeof(grid_width_));

    nb_indivs_ = grid_height_ * grid_width_;

    internal_organisms_ = new std::shared_ptr<Organism>[nb_indivs_];
    prev_internal_organisms_ = new std::shared_ptr<Organism>[nb_indivs_];

    // No need to save/load this field from the backup because it will be set at selection()
    next_generation_reproducer_ = new int[nb_indivs_]();

    gzread(exp_backup_file, &backup_step_, sizeof(backup_step_));

    gzread(exp_backup_file, &mutation_rate_, sizeof(mutation_rate_));

    for (int i = 0; i < FUZZY_SAMPLING; i++) {
        double tmp;
        gzread(exp_backup_file, &tmp, sizeof(tmp));
        target[i] = tmp;
    }

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id] =
                std::make_shared<Organism>(exp_backup_file, *DNA_seqs);
    }

    rng_ = std::move(std::make_unique<Threefry>(grid_width_, grid_height_, exp_backup_file));
    seed_ = rng_->get_seed();

    if (gzclose(exp_backup_file) != Z_OK) {
        cerr << "Error while closing backup file" << endl;
    }
}

/**
 * Destructor of the ExpManager class
 */
ExpManager::~ExpManager() {
    delete stats_best;
    delete stats_mean;

    delete[] dna_mutator_array_;

    delete[] internal_organisms_;
    delete[] prev_internal_organisms_;
    delete[] next_generation_reproducer_;
    delete[] target;
}

/**
 * Selection process: for a given cell in the grid of the population, compute which organism win the computation
 * 
  * @param indiv_id : Unique identification number of the cell
  * 
  * TODO: Optimize function call (parallel calls)
 */
void ExpManager::selection(int indiv_id) const {
    double local_fit_array[NEIGHBORHOOD_SIZE];
    double probs[NEIGHBORHOOD_SIZE];
    int count = 0;
    double sum_local_fit = 0.0;

    int32_t x = indiv_id / grid_height_;
    int32_t y = indiv_id % grid_height_;

    int cur_x, cur_y;

    for (int8_t i = -1; i < NEIGHBORHOOD_WIDTH - 1; i++) {
        for (int8_t j = -1; j < NEIGHBORHOOD_HEIGHT - 1; j++) {
            cur_x = (x + i + grid_width_) % grid_width_;
            cur_y = (y + j + grid_height_) % grid_height_;

            local_fit_array[count] = prev_internal_organisms_[cur_x * grid_width_ + cur_y]->fitness;
            sum_local_fit += local_fit_array[count];

            count++;
        }
    }

    for (int8_t i = 0; i < NEIGHBORHOOD_SIZE; i++) {
        probs[i] = local_fit_array[i] / sum_local_fit;
    }

    auto rng = std::move(rng_->gen(indiv_id, Threefry::REPROD));
    int found_org = rng.roulette_random(probs, NEIGHBORHOOD_SIZE);

    int x_offset = (found_org / NEIGHBORHOOD_WIDTH) - 1;
    int y_offset = (found_org % NEIGHBORHOOD_HEIGHT) - 1;

    next_generation_reproducer_[indiv_id] = ((x + x_offset + grid_width_) % grid_width_) * grid_height_ +
                                            ((y + y_offset + grid_height_) % grid_height_);
}

/**
 * Prepare the mutation generation of an organism
 *
 * @param indiv_id : Organism unique id
 */
void ExpManager::prepare_mutation(int indiv_id) const {
    auto *rng = new Threefry::Gen(std::move(rng_->gen(indiv_id, Threefry::MUTATION)));
    const shared_ptr<Organism> &parent = prev_internal_organisms_[next_generation_reproducer_[indiv_id]];
    dna_mutator_array_[indiv_id] = new DnaMutator(
            rng,
            parent->length(),
            mutation_rate_);
    dna_mutator_array_[indiv_id]->generate_mutations();

    if (dna_mutator_array_[indiv_id]->hasMutate()) {
        internal_organisms_[indiv_id] = std::make_shared<Organism>(parent, indiv_id);
    } else {
        int parent_id = next_generation_reproducer_[indiv_id];

        internal_organisms_[indiv_id] = prev_internal_organisms_[parent_id];
        internal_organisms_[indiv_id]->reset_mutation_stats();
    }
}

/**
 * Execute a generation of the simulation for all the Organisms
 *
 */
__device__ void ExpManager::run_a_step(
    size_t indiv_id, 
    Kokkos::View<
        double*, 
        Kokkos::DefaultExecutionSpace::memory_space, 
        Kokkos::MemoryTraits<Kokkos::Atomic>
    > bestFitnessYet
    ) {

    // Running the simulation process for each organism
    selection(indiv_id);
    prepare_mutation(indiv_id);

    if (dna_mutator_array_[indiv_id]->hasMutate()) {
        auto &mutant = internal_organisms_[indiv_id];
        mutant->apply_mutations(dna_mutator_array_[indiv_id]->mutation_list_);
        mutant->evaluate(target);
    }

    Kokkos::fence("Fence 1");

    // Swap Population
    prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id];
    internal_organisms_[indiv_id] = nullptr;

    if(indiv_id == 0) {
        bestFitnessYet(0) = prev_internal_organisms_[0]->fitness;
        // Stats
        stats_best->reinit(AeTime::time());
        stats_mean->reinit(AeTime::time());
    }

    Kokkos::fence("Fence 2");

    // Search for the best
    if (prev_internal_organisms_[indiv_id]->fitness > bestFitnessYet(0)) {
        bestFitnessYet(0) = prev_internal_organisms_[indiv_id]->fitness;
    }

    Kokkos::fence("Fence 3");
    if (dna_mutator_array_[indiv_id]->hasMutate())
        prev_internal_organisms_[indiv_id]->compute_protein_stats();

    Kokkos::fence("Fence 4");
    if (indiv_id == 0) {
        stats_best->write_best(best_indiv);
        stats_mean->write_average(prev_internal_organisms_, nb_indivs_);
    }
}


/**
 * Run the evolution for a given number of generation
 *
 * @param nb_gen : Number of generations to simulate
 */
void ExpManager::run_evolution(int nb_gen) {
    INIT_TRACER("trace.csv", {"FirstEvaluation", "STEP"});
    Kokkos::parallel_for(
            "ExpManager::run_a_step 3",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, nb_indivs_),
            [=] (const size_t indiv_id) 
        {

            internal_organisms_[indiv_id]->locate_promoters();
            prev_internal_organisms_[indiv_id]->evaluate(target);
            prev_internal_organisms_[indiv_id]->compute_protein_stats();

            // Stats
            if (indiv_id == 0) {
                stats_best = new Stats(AeTime::time(), true);
                stats_mean = new Stats(AeTime::time(), false);

                printf("Running evolution from %d to %d\n", AeTime::time(), AeTime::time() + nb_gen);
            }
        }
    );

    FLUSH_TRACES(0)

    std::cout << "just before bestFitnessYet() " << std::endl;
    Kokkos::View<
        double*, 
        Kokkos::DefaultExecutionSpace::memory_space, 
        Kokkos::MemoryTraits<Kokkos::Atomic>
    > bestFitnessYet(
        "best fitness", 1
    );
    std::cout << "declared bestFitnessYet() " << std::endl;

    // Stats
    stats_best = new Stats(AeTime::time(), true);
    stats_mean = new Stats(AeTime::time(), false);

    for (int gen = 0; gen < nb_gen; gen++) {
        // each generation is dependent on the previous one
        AeTime::plusplus();

        Kokkos::parallel_for(
            "ExpManager::run_a_step 3",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, nb_indivs_),
            [=] (const size_t indiv_id) 
            {
                run_a_step(indiv_id, bestFitnessYet);

                if(indiv_id == 0) {
                    printf("Generation %d : Best individual fitness %e\n", AeTime::time(), best_indiv->fitness);
                }

                delete dna_mutator_array_[indiv_id];
                dna_mutator_array_[indiv_id] = nullptr;

                if(indiv_id == 0) {
                    if (AeTime::time() % backup_step_ == 0) {
                        save(AeTime::time());
                        cout << "Backup for generation " << AeTime::time() << " done !" << endl;
                    }
                }
            }
        );
    }
    STOP_TRACER
}
