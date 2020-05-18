"""
Created: 2018
Updated: may 2019
Commented: may 2020
Author: Matteo Tomasini
Paper: When does gene flow facilitate evolutionary rescue?
doi: https://doi.org/10.1101/622142

In this script, we simulate evolutionary rescue in a habitat with 2 demes.
To work, the script needs a file called "parameters.txt" that contains the parameters (see example in the folder).

Modify the variables in the beginning of the script to change the parameters of the simulation.

Python version: 3.X (last tested: 3.6)
"""

import numpy as np
import math
import time

BURNING_PHASE = 0                       # no burning phase as we normally start with initial frequency f_0
BURNING_OFF = 500                       # number of generations simulated after whole deme is deteriorated
STANDING_GENETIC_VARIATION = True       # do we simulate standing genetic variation? (If false sets initial_frequency = 0)
IS_BEVERTON_HOLT = False                # if True, uses Beverton-Hold dynamics in undeteriorated deme
DENSITY_REGULATION = True               # if True, regulates density to carrying capacity
THRESHOLD = 0.5                         # percentage of carrying capacity necessary to declare rescue


# FUNCTIONS


def zero_maker(length):
    """
    Custom function to make list of zeros (for readability)
    """
    list_of_zeros = [0] * length
    return list_of_zeros


def positive(number):
    """
    Custom function to set to zero negative numbers, does nothing with numbers >= 0
    """
    if number < 0:
        return 0
    else:
        return number


def deteriorate(what_time_is_it):
    """
    Sets is_deteriorated of a deme to True according to the time in the simulation
    """
    for it in iteration_list:
        if what_time_is_it >= (BURNING_PHASE + it * epoch):
            deme_deterioration_state[it] = True


def expected_offspring(is_deteriorated, deme_population, deme_capacity, is_mutant):
    """
    Calculates the expected mutant or wildtype offspring in a deme according to whether that deme is deteriorated, its population.
    """

    if is_deteriorated and is_mutant:
        fitness = 1. + selection_for
    elif is_deteriorated and not is_mutant:
        fitness = 1. - stress
    elif not is_deteriorated and is_mutant:
        fitness = 1. - selection_against
    else:
        fitness = 1.

    if not IS_BEVERTON_HOLT:
        average_offspring = fitness
    else:
        if not is_deteriorated:
            average_offspring = growth_rate * fitness / (1 + (growth_rate - 1) * deme_population / deme_capacity)
        else:
            average_offspring = fitness

    return average_offspring


def reproduce(mutants_in_world, wildtypes_in_world):
    """
    Extract poisson distributed numbers for offsprings according to the expected offspring
    """
    for i in iteration_list:
        total_population = mutants_in_world[i] + wildtypes_in_world[i]
        expected_mutants = expected_offspring(deme_deterioration_state[i], total_population,
                                              carrying_capacities[i], True)
        expected_wildtypes = expected_offspring(deme_deterioration_state[i], total_population,
                                                carrying_capacities[i], False)
        mutants_in_world[i] = np.random.poisson(mutants_in_world[i] * expected_mutants)
        wildtypes_in_world[i] = np.random.poisson(wildtypes_in_world[i] * expected_wildtypes)


def generate_new_mutants(wildtypes):
    """
    If we are within the time limits set by the variable WINDOW_OF_MUTATION, generates new mutants from wildtype with binomial sampling
    """
    new_mutants = zero_maker(number_of_demes)
    if WINDOW_OF_MUTATION[0] < generation < WINDOW_OF_MUTATION[1]:
        for i in iteration_list:
            if wildtypes[i] > 0:
                new_mutants[i] = np.random.binomial(wildtypes[i], mutation_rate)

    return new_mutants


def migrate(mutants_in_world, wildtypes_in_world):
    """
    Models migration as binomial sampling: first selects wildtypes going right from left deme, then viceversa, then moves them.
    Distinction between wildtypes and mutants is done according to current mutant_frequencies
    """
    total_population = [0, 0]
    mutant_frequencies = [0., 0.]
    for i in iteration_list:
        total_population[i] = mutants_in_world[i] + wildtypes_in_world[i]
        if total_population[i] > 0:
            mutant_frequencies[i] = float(mutants_in_world[i]) / total_population[i]

    wildtypes_going_right = np.random.binomial(wildtypes_in_world[0], migration_rate_1to2)
    mutants_going_right = np.random.binomial(mutants_in_world[0], migration_rate_1to2)
    wildtypes_going_left = np.random.binomial(wildtypes_in_world[1], migration_rate_2to1)
    mutants_going_left = np.random.binomial(mutants_in_world[1], migration_rate_2to1)

    mutants_in_world[0] = mutants_in_world[0] - mutants_going_right + mutants_going_left
    mutants_in_world[1] = mutants_in_world[1] - mutants_going_left + mutants_going_right
    wildtypes_in_world[0] = wildtypes_in_world[0] - wildtypes_going_right + wildtypes_going_left
    wildtypes_in_world[1] = wildtypes_in_world[1] - wildtypes_going_left + wildtypes_going_right


def down_regulate_density(mutants_in_world, wildtypes_in_world):
    """
    If population in non-deteriorated deme is too high, down-regulates to carrying capacity, according to frequencies (with poisson sampling)
    """
    total_population = [0, 0]
    mutant_frequencies = [0., 0.]
    for i in iteration_list:
        total_population[i] = mutants_in_world[i] + wildtypes_in_world[i]
        if total_population[i] > 0:
            mutant_frequencies[i] = float(mutants_in_world[i]) / total_population[i]

        if total_population[i] > carrying_capacities[i]:
            number_of_dead = total_population[i] - carrying_capacities[i]
            dead_mutants = np.random.binomial(number_of_dead, mutant_frequencies[i])
            dead_wildtypes = number_of_dead - dead_mutants
            mutants_in_world[i] = max(mutants_in_world[i] - dead_mutants, 0)
            wildtypes_in_world[i] = max(wildtypes_in_world[i] - dead_wildtypes, 0)


def up_regulate_density(mutants_in_world, wildtypes_in_world):
    """
    If population in non-deteriorated deme is too low, down-regulates to carrying capacity, according to frequencies (with poisson sampling).
    """
    total_population = [0, 0]
    mutant_frequencies = [0., 0.]
    for i in iteration_list:
        total_population[i] = mutants_in_world[i] + wildtypes_in_world[i]
        if total_population[i] > 0:
            mutant_frequencies[i] = float(mutants_in_world[i]) / total_population[i]

        if not deme_deterioration_state[i]:           # Up regulation only occurs in non-deteriorated state (when False)
            if total_population[i] < carrying_capacities[i]:
                new_total_population = np.random.poisson(carrying_capacities[i])
                mutants_in_world[i] = np.random.binomial(new_total_population, mutant_frequencies[i])
                wildtypes_in_world[i] = max(new_total_population - mutants_in_world[i], 0)


def is_rescued(percentage_threshold):
    """
    Changes flag is_rescued to True at the end of generation, so that simulation can stop when rescue is achieved.
    Condition for rescue: if mutants have attained a certain threshold (in terms of percentage carrying capacity in the whole habitat (cf Uecker et al. 2014))
    """
    if generation > BURNING_PHASE:   # rescue can only occur after the borning phase
        if sum(list_of_mutants) >= percentage_threshold * total_carrying_capacity:
            return True
        else:
            return False
    else:
        return False


# START SIMULATION

start_time = time.time()

# download the parameters
file_of_parameters = open('Parameters.txt', 'r')
parameters_to_read = np.loadtxt(file_of_parameters, skiprows=1, delimiter='\t')
file_of_parameters.close()

# write header into output file
header = [['OneEpoch', 'Demes', 'TotalKappa', 'migration', 'mutation', 'frequency', 'growth', 's', 'r',
           'z', 'RatioCapacities', 'RatioMigration', 'replicates', 'Rescue', 'Error']]
output_file = 'Simulation_output'
output = open('%s.txt' % output_file, 'w')
np.savetxt(output, header, fmt='%s')
output.close()

for line in parameters_to_read:
    [number_of_generations, number_of_demes, total_carrying_capacity, migration_rate, mutation_rate,
     initial_frequency, growth_rate, selection_against, stress, selection_for, ratio_between_capacities,
     ratio_between_migrations, replicates, probability_of_rescue, confidence_interval] = line

    number_of_demes = int(number_of_demes)

    if number_of_demes != 2:
        raise ValueError('There are more than 2 demes!')

    epoch = number_of_generations

    iterate_over_replicates = np.linspace(1., replicates, replicates)

    total_replicate_time = BURNING_PHASE + epoch + BURNING_OFF
    iterate_over_generations = np.linspace(1., total_replicate_time, total_replicate_time)

    # This is the window of time within which a mutation can occur. It is set to a very large value, but can be reduced.
    # (useful if one wants to look only at standing genetic variation)
    WINDOW_OF_MUTATION = [0, 100000000]

    deme_carrying_capacity_1 = math.ceil(total_carrying_capacity * ratio_between_capacities)
    deme_carrying_capacity_2 = total_carrying_capacity - deme_carrying_capacity_1
    carrying_capacities = [deme_carrying_capacity_1, deme_carrying_capacity_2]
    iteration_list = list(range(number_of_demes))

    migration_rate_1to2 = migration_rate * ratio_between_migrations
    migration_rate_2to1 = migration_rate * (1 - ratio_between_migrations)

    if not STANDING_GENETIC_VARIATION:
        initial_frequency = 0.0

    result_of_replicate = []  # this keeps track of rescue events.

    for replicate in iterate_over_replicates:
        has_survived_yet = False

        list_of_populations = [deme_carrying_capacity_1, deme_carrying_capacity_2]

        list_of_mutants = []
        for i in iteration_list:
            mutants_in_deme = math.ceil(list_of_populations[i] * initial_frequency)
            list_of_mutants.append(mutants_in_deme)

        list_of_wildtypes = []
        for i in iteration_list:
            wildtypes_in_deme = list_of_populations[i] - list_of_mutants[i]
            list_of_wildtypes.append(wildtypes_in_deme)

        deme_deterioration_state = [False] * number_of_demes

        for generation in iterate_over_generations:
            if has_survived_yet:
                break

            if DENSITY_REGULATION:
                down_regulate_density(list_of_mutants, list_of_wildtypes)
                up_regulate_density(list_of_mutants, list_of_wildtypes)

            deteriorate(generation)

            reproduce(list_of_mutants, list_of_wildtypes)

            list_of_new_mutants = generate_new_mutants(list_of_wildtypes)

            total_after_reproduction = []

            for i in iteration_list:
                if list_of_new_mutants[i] > list_of_wildtypes[i]:
                    list_of_new_mutants[i] = list_of_wildtypes[i]

                list_of_mutants[i] += int(list_of_new_mutants[i])
                list_of_wildtypes[i] -= int(list_of_new_mutants[i])
                total_after_reproduction.append(list_of_mutants[i] + list_of_wildtypes[i])

            migrate(list_of_mutants, list_of_wildtypes)

            for i in iteration_list:
                list_of_populations[i] = list_of_mutants[i] + list_of_wildtypes[i]

            if sum(list_of_populations) == 0:
                break

            has_survived_yet = is_rescued(THRESHOLD)

            if generation == total_replicate_time:
                has_survived_yet = True

        if has_survived_yet:
            result_of_replicate.append(True)

    final_probability_of_rescue = sum(result_of_replicate) / replicates
    final_confidence_interval = math.fabs(
        -2.57 * math.sqrt(final_probability_of_rescue * (1 - final_probability_of_rescue) / replicates))

    parameters_to_append = [[number_of_generations, number_of_demes, total_carrying_capacity, migration_rate,
                             mutation_rate, initial_frequency, growth_rate, selection_against,
                             stress, selection_for, ratio_between_capacities, ratio_between_migrations, replicates,
                             final_probability_of_rescue, final_confidence_interval]]

    output = open('%s.txt' % output_file, 'a')
    np.savetxt(output, parameters_to_append, delimiter='\t', fmt='%.5f')
    output.close()

print('time of execution: {}'.format(time.time() - start_time))
