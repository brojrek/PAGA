from Chromosome import *

#TODO Validate somewhere -> number of parents (got from selections - will be added in the future) 
# plus number of parents times a number of mutations must be LESS than population_size.
# i.e. if number of parents = 10
# and number of  mutations = 7
# population_size must be greater than 10 + 10*7 = 80 (in that case, 100 would be a good minimum population_size)

class Population:
    '''
    Class representing a population of chromosomes. <br>
    Objects contain a list of objects of Chromosome class.
    '''

    def __init__(self, primogenitor_chromosome, population_size):

        self.primogenitor_chromosome = primogenitor_chromosome
        self.population = [primogenitor_chromosome]
        self.population_size = population_size
        self.parent_number = 5

        self.fill_population()

    @property
    def population(self):
        return self._population
    
    @property
    def primogenitor_chromosome(self):
        return self._primogenitor_chromosome
    
    @property
    def parent_number(self):
        return self._parent_number
    
    @population.setter
    def population(self, population):
        self._population = population
    
    @primogenitor_chromosome.setter
    def primogenitor_chromosome(self, primogenitor_chromosome):
        if not type(primogenitor_chromosome) is Chromosome: raise Exception("primogenitor_chromosome must be an instance of Chromosome")
        self._primogenitor_chromosome = primogenitor_chromosome

    @parent_number.setter
    def parent_number(self, parent_number):
        self._parent_number = parent_number

    def fill_population(self):
        while len(self.population) < self.population_size:
            self.population.append(self.primogenitor_chromosome.offspring_mutate_add_gaps(False))
            self.population = sorted(self.population, reverse = True)
            sorted(self.population, reverse = True)

    def new_generation(self, parents): 
        self.population = parents[:]

        # Prolongation
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_prolongation())

        # Add one gap
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_add_gaps(True))

        # Add multiple gaps
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_add_gaps(False))

        # Shuffle gaps
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_shuffle_gaps())

        # Move gap
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_move_gap())

        # Move a section of gaps
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_move_section())

        # Remove gaps
        for chromosome in parents:
            self._population.append(chromosome.offspring_mutate_remove_gaps())

        if len(self.population) > self.population_size: raise Exception("new_generation: generated population size (" + len(self.population) + ") exceeds population_size value (" + self.population_size + "). Please, adjust input parameters to avoid this.")

        self.fill_population()

    def selection_simple(self):
        return self.population[:self.parent_number]