"""
This file generates a file of parameters to be used in "main_simulation.py".
"""

import numpy as np

header = [['Theta', 'Demes', 'TotalKappa', 'migration', 'isGlobal', 'mutation', 'frequency', 'growth', 's2', 'r', 's1',
           'RatioCapacities', 'RatioMigration', 'replicates', 'Rescue', 'Error']]
output = open('Parameters.txt', 'w')
np.savetxt(output, header, fmt='%s')
output.close()

replicates = 2000
probability_of_rescue = 0.0
confidence_interval = 0.0

# structure of the habitat #
list_of_demes = [2]
# epoch = 50.
list_of_epoch = [500]
list_of_carrying_capacities = [20000]
list_of_migration_rates = np.logspace(-4, 0, 25)
list_of_mutation_rates = [1. / list_of_carrying_capacities[0]]
initial_frequencies = [0.0000]   # These are redefined later, in the loop.
list_of_growth_rates = [1.5]  # Growth rates must be > 1
list_of_capacity_ratios = [0.5]
list_of_migration_ratios = [0.5]

# selection coefficients
list_of_selection_for = [0.02]    # Selection in favor of mutants in deteriorated environment, > 0
list_of_selection_against = [0.1, 0.9]  # Selection against mutants in non-deteriorated environment, > 0
list_of_stress = [0.3]  # Selection against wildtype in deteriorated environment , > 0

for total_number_of_generations in list_of_epoch:
    for total_carrying_capacity in list_of_carrying_capacities:
        for demes in list_of_demes:
            for migration_rate in list_of_migration_rates:
                for mutation_rate in list_of_mutation_rates:
                    for mutant_frequency in initial_frequencies:
                        for growth_rate in list_of_growth_rates:
                            for selection_against in list_of_selection_against:
                                mutant_frequency = mutation_rate / selection_against
                                for stress in list_of_stress:
                                    for selection_for in list_of_selection_for:
                                        for ratio_between_capacities in list_of_capacity_ratios:
                                            for ratio_between_migrations in list_of_migration_ratios:
                                                parameters = [[total_number_of_generations, demes,
                                                               total_carrying_capacity, migration_rate,
                                                               mutation_rate, mutant_frequency, growth_rate,
                                                               selection_against, stress, selection_for,
                                                               ratio_between_capacities, ratio_between_migrations,
                                                               replicates, probability_of_rescue, confidence_interval]]
                                                output = open('Parameters.txt', 'a')
                                                np.savetxt(output, parameters, delimiter='\t', fmt='%.5f')
                                                output.close()
