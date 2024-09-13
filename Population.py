from Chromosome import *
import numpy as np

class Population:
    """
    Class representing a population of chromosomes. <br>
    Objects contain a list of objects of Chromosome class.

    Attributes
    ----------
    primogenitor_chromosome : Chromosome
        a chromosome that will be a parent of whole population
    population_size : int
        total number of chromosomes in population

    prolongation_active : boolean
        default True
    add_one_gap_active : boolean
        default True
    add_multiple_gaps_active : boolean
        default True
    gap_shuffle_gaps_active : boolean
        default True
    move_gap_active : boolean
        default True
    move_section_active : boolean
        default True
    gap_remove_active : boolean
        default True
    Methods
    -------
    new_generation():
        Does this and that and then returns sth

    fill_population():
        Does this and that and then returns sth

    selection_truncation():
        Does this and that and then returns sth
 
    selection_roulette():
        Does this and that and then returns sth
 
    selection_tournament():
        Does this and that and then returns sth
 
    """

    def __init__(self, primogenitor_chromosome, population_size):
        """
        Constructs all the necessary attributes for the population object.

        Parameters
        ----------
        primogenitor_chromosome : Chromosome
            

        population_size : list
            list which first element is a scoring matrix (NDarray) and second is an alphabet (dict)
        penalty : int
            penalty for gaps in the sequence, used in calculationg score
        gap_prolong_odds : int, optional
            odds 1:gap_prolong_odds for each gap in the sequence to be prolonged in offspring_mutate_prolongation mutation
        gap_shuffle_odds : int, optional
            odds 1:gap_shuffle_odds for each already existing gap to be randomly shuffled in sequence.
        gap_remove_odds : int, optional
            odds 1:gap_remove_odds for each gap in the sequence to be removed in offspring_mutaute_remove_gaps mutation
        new_gaps_limit_factor : int, optional 
            number of gaps to be added randomly is a number between (1) and (len(sequence)/(new_gaps_limit_factor))
        """
        self.primogenitor_chromosome = primogenitor_chromosome
        self.population = [primogenitor_chromosome]
        self.population_size = population_size
        self.parent_number = 5

        self.prolongation_active = True
        self.add_one_gap_active = True
        self.add_multiple_gaps_active = True
        self.gap_shuffle_gaps_active = True
        self.move_gap_active = True
        self.move_section_active = True
        self.gap_remove_active = True

        self.fill_population()

    def __len__(self):
        return len(self.population)
    
    @property
    def population(self):
        """Get or set the population"""
        return self._population
    
    @property
    def primogenitor_chromosome(self):
        """Get or set the primogenitor_chromosome"""
        return self._primogenitor_chromosome

    @property
    def parent_number(self):
        """Get or set the parent_number"""
        return self._parent_number
    
    @property
    def prolongation_active(self):
        """Get or set the prolongation_active"""
        return self._prolongation_active

    @property
    def add_one_gap_active(self):
        """Get or set the add_one_gap_active"""
        return self._add_one_gap_active
    
    @property
    def add_multiple_gaps_active(self):
        """Get or set the add_multiple_gaps_active"""
        return self._add_multiple_gaps_active
    
    @property
    def gap_shuffle_gaps_active(self):
        """Get or set the gap_shuffle_gaps_active"""
        return self._gap_shuffle_gaps_active
    
    @property
    def move_gap_active(self):
        """Get or set the move_gap_active"""
        return self._move_gap_active

    @property
    def move_section_active(self):
        """Get or set the move_section_active"""
        return self._move_section_active
        
    @property
    def gap_remove_active(self):
        """Get or set the gap_remove_active"""
        return self._gap_remove_active
    
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

    @prolongation_active.setter
    def prolongation_active(self, prolongation_active):
        self._prolongation_active = prolongation_active

    @add_one_gap_active.setter
    def add_one_gap_active(self, add_one_gap_active):
        self._add_one_gap_active = add_one_gap_active
    
    @gap_shuffle_gaps_active.setter
    def add_multiple_gaps_active(self, add_multiple_gaps_active):
        self._add_multiple_gaps_active = add_multiple_gaps_active
    
    @gap_shuffle_gaps_active.setter
    def gap_shuffle_gaps_active(self, gap_shuffle_gaps_active):
        self._gap_shuffle_gaps_active = gap_shuffle_gaps_active
    
    @move_gap_active.setter
    def move_gap_active(self, move_gap_active):
        self._move_gap_active = move_gap_active

    @move_section_active.setter
    def move_section_active(self, move_section_active):
        self._move_section_active = move_section_active
        
    @gap_remove_active.setter
    def gap_remove_active(self, gap_remove_active):
        self._gap_remove_active = gap_remove_active

    def fill_population(self):
        while len(self.population) < self.population_size:
            self.population.append(self.primogenitor_chromosome.offspring_mutate_add_gaps(False))
            self.population = sorted(self.population, reverse = True)
            sorted(self.population, reverse = True)

    def new_generation(self, parents): 
        """
        Method replaces current population with new one based on selected parents and their offspring (generated by enabled mutations).
        After appending offspring to parent list, fill the population with random childreen of progenitor_chromosome.

        Parameters
        ----------
        parents : list[Chromosome]
            A list of chromosomes, selected to be a parents to new generation.
        -------
        """
        self.population = parents[:]

        # Prolongation
        if self.prolongation_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_prolongation())
        # Add one gap
        if self.add_one_gap_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_add_gaps(True))

        # Add multiple gaps
        if self.add_multiple_gaps_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_add_gaps(False))

        # Shuffle gaps
        if self.gap_shuffle_gaps_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_shuffle_gaps())

        # Move gap
        if self.move_gap_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_move_gap())

        # Move a section of gaps
        if self.move_section_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_move_section())

        # Remove gaps
        if self.gap_remove_active:
            for chromosome in parents:
                self.population.append(chromosome.offspring_mutate_remove_gaps())

        if len(self.population) > self.population_size: raise Exception("new_generation: generated population size (" + len(self.population) + ") exceeds population_size value (" + self.population_size + "). Please, adjust input parameters to avoid this.")

        self.fill_population()

    def selection_truncation(self):
        """
        Function performs a truncation selection on self.population.

        Truncation method - select self.parent_number number of best chromosomes from population.
        
        Returns
        ----------
        list
            A list of self.parent_number chromosomes with a highest score
        """
        return self.population[:self.parent_number]
    
    def selection_roulette(self):
        """
        Function performs a roulette selection, always keeping the single highscore Chromosome.

        Roulette method - based on it's score, every chromosome has proportional chance to be selected. 
        Perform random, weighted selection of self.parent_number chromosomes from population.

        Returns
        ----------
        list
            A list of size self.parent_number containing chromosome with a highest score and chromosomes selected with a roulette method
        """
        current_population = sorted(self.population[:], reverse=True)
        best = [current_population.pop(0)]
        min_score = current_population[-1].score

        if min_score <= 0:
            weights = sorted([(x.score + abs(min_score) + 1) for x in current_population], reverse=True)
        elif min_score == 1:
            weights = sorted([x.score for x in current_population], reverse=True)
        else: #min_score > 1
            weights = sorted([x.score-min_score for x in current_population], reverse=True)
        sum_score = sum(weights)
        weights = [x/sum_score for x in weights]

        return best + list(np.random.choice(a = current_population, size = self.parent_number-1, p = weights, replace = False))

    def selection_tournament(self):
        """
        Function performs a tournament selection.
        
        Tournament method - split population to self.parent_number equal groups. Select best chromosome from each group.
        Returns
        ----------
        list
            A list of size self.parent_number containing chromosomes selected with a tournament method
        """
        current_population = self.population[:]
        tournament_size = len(current_population)//self.parent_number
        new_population = []
        
        while len(current_population) > tournament_size:
            tournament_group = []
            while len(tournament_group) < tournament_size:
                tournament_group.append(current_population.pop(random.randint(0, len(current_population)-1)))
            new_population.append(max(tournament_group))
            
        if len(current_population) > 0:
            new_population.append(max(current_population))
        
        return new_population