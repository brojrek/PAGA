import random
import numpy

class Chromosome:
    '''
    Class representing a chromosome - a single unit in the population.
    Objects contain two sequences, a pair (to be) aligned.
    '''

    gap_prolong_odds = 6 # Odds 1:gap_prolong_odds for each gap in the sequence to be prolonged in offspring_mutate_prolongation mutation
    gap_shuffle_odds = 3 # Odds 1:gap_shuffle_odds for each already existing gap to be randomly shuffled in sequence.
    gap_remove_odds = 4  # Odds 1:gap_remove_odds
    new_gaps_limit_factor = 10 # offspring_mutate_add_gaps -> number of gaps to be added randomly is a number between (1) and (len(sequence)/(new_gaps_limit_factor))


    def __init__(self, sequence_A, sequence_B, scoring, penalty):
        self.sequence_A = sequence_A
        self.sequence_B = sequence_B

        self.matrix = scoring[0]
        self.alphabet = scoring[1]
        self.penalty = penalty

        self.equalize()  
        self.calculate_score()

    def __eq__(self, other):
        return self.score == other.score and self.sequence_A == other.sequence_A and self.sequence_B == self.sequence_B

    def __lt__(self, other):
        return self.score < other.score
    
    def __gt__(self, other):
        return self.score > other.score

    def calculate_score(self):
        score = 0
        for idx in range(0, self.get_length()):
        
            character_A = self.sequence_A[idx]
            character_B = self.sequence_B[idx]

            if character_A == '-' or character_B == '-':
                score -=  self.penalty
            
            else:
                row = self.alphabet[character_A]
                col = self.alphabet[character_B]
                score = score + self.matrix[row, col]
        self.score = int(score)

    @property
    def sequence_A(self):
        return self._sequence_A
    
    @property
    def sequence_B(self):
        return self._sequence_B
    
    @property
    def penalty(self):
        return self._penalty
    
    @sequence_A.setter
    def sequence_A(self, sequence_A):
        if not sequence_A or sequence_A == "": raise Exception("Sequence_A cannot be empty")
        self._sequence_A = sequence_A

    @sequence_B.setter
    def sequence_B(self, sequence_B):
        if not sequence_B or sequence_B == "": raise Exception("Sequence_B cannot be empty")
        self._sequence_B = sequence_B
    
    @penalty.setter
    def penalty(self, penalty):
        if penalty <= 0: raise Exception("Penalty must be greater than 0")
        self._penalty = penalty

    def get_length(self):
        """
        Asserts that both sequences are even length. Gets sequences length.
        returns: int representing a length of sequences on a chromosone 
        """

        assert(len(self.sequence_A) == len(self.sequence_B))
        return len(self.sequence_A)
    
    def get_score(self):
        """
        returns: int representing a score of sequences alignment (for given scoring matrix, alphabet and gap penalty)
        """

        return self.score
    
    def set_sequences(self, sequence_A, sequence_B):
        """
        Sets a value of a pair of sequences on a chromosome.
        """

        self.sequence_A = sequence_A
        self.sequence_B = sequence_B
        self.equalize()

    def equalize(self):
        '''
        Adds gaps on the end of shorter sequence to match the length of the second sequence.
        '''

        len_A = len(self.sequence_A)
        len_B = len(self.sequence_B)

        if len_A > len_B:
            self.sequence_B = self.sequence_B + ("-" * (len_A-len_B))
            
        if len_A < len_B:
            self.sequence_A = self.sequence_A + ("-" * (len_B-len_A))

        self.gap_reduction()

        assert len(self.sequence_A) == len(self.sequence_B)

    def gap_reduction(self):
        '''
        Function that removes unnecessary gaps (when both sequences have gap on the same index).
        '''
        
        for idx in range(self.get_length()-1, 0, -1):
            if self.sequence_A[idx] == '-' and self.sequence_B[idx] == '-':
                self.sequence_A = self.sequence_A[:idx] + self.sequence_A[idx+1:]
                self.sequence_B = self.sequence_B[:idx] + self.sequence_B[idx+1:]

    def gap_indexes(self, seq):
        indexes = []
        for idx in range(len(seq)-1, 0, -1):
            if seq[idx] == '-':
                indexes.append(idx)
        return indexes
    
    def gap_sections_indexes(self, seq):
        gap_sections = []

        recording = False
        for idx in range(len(seq)-1, 0, -1):
            if seq[idx] == '-':
                if not recording:
                    recording = True
                    idx_end = idx
            else:
                if recording:
                    recording = False
                    gap_sections.append((idx+1, idx_end))
        return gap_sections
    '''
    def randomly_swap_index(idx, seq):
    # if idx is the first index, move it right (only option)
        if idx == 0:
            seq = seq[1] + "-" + seq[2:]
        # if idx is the last index, move it left (only option)
        elif idx == len(seq) - 1:
            seq = seq[:-2] + "-" + seq[-2]
        # if idx is not at beginning nor at end of a sequence, randomly decide whether to move it left or right
        else:
            coin_flip = random.randint(0, 1)
            # move right
            if coin_flip == 0: 
                seq = seq[:idx] + seq[idx+1] + "-" + seq[idx+2:] 
            # move left
            else: 
                seq = seq[:idx-1] + "-" + seq[idx-1] + seq[idx+1:] 
        return seq
    '''
    def randomly_swap_indexes(self, seq, idx_start, idx_end = None):
        '''
        Asserting attempt to move whole string (idx_start = 0 and idx_end = len(seq))
        This scenario is not expected in our class.
        '''
        if idx_end == None:
            idx_end = idx_start

        assert not (idx_start == 0 and idx_end == len(seq))
        sliced = seq[idx_start:idx_end+1]

        # if idx_start is the first index, move it right (only option)
        if idx_start == 0:
            seq = seq[idx_end+1] + sliced + seq[idx_end+2:]

        # if idx_end is the last index, move it left (only option)
        elif idx_end == len(seq) - 1:
            seq = seq[:idx_start-1] + sliced + seq[idx_start-1]

        # in other cases, randomly decide whether to move it left or right
        else:
            coin_flip = random.randint(0, 1)
            # move right
            if coin_flip == 0: 
                seq = seq[:idx_start] + seq[idx_end+1] + sliced + seq[idx_end+2:]
            # move left
            else: 
                seq = seq[:idx_start-1] + sliced + seq[idx_start-1] + seq[idx_end+1:]
        return seq

    # Mutation: Prolong existing gaps. [Look for local optimum.]
    def offspring_mutate_prolongation(self):
        '''
        Function performes a prolongation mutation (by prolonging gaps) on the chromosome and returns the result as new chromosome (an offspring).
        Every gap existing in the sequence has a fixed probability to be doubled.
    
        Description: 
            1. randomly choose one of two sequences in a chromosome (chosen_sequence)
            2. iterate through chosen sequence looking for gaps
                With fixed odds of 1:(self.prolongation_odds) add current gap index to list (indexes_to_insert)
            3. For each index remembered on the list, add a new gap to chosen_sequence.
            4. Create and equalize offspring - a new Chromosome of mutated sequences
            5. Return an offspring (Chromosome)

        returns: new Chromosome 
        '''

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        coin_flip = random.randint(0, 1) # Randomly select sequence_A or sequence_B to mutate
        if coin_flip == 0:
            chosen_sequence = sequence_A
        else:
            chosen_sequence = sequence_B
        
        indexes_to_insert = []
        for k in range(len(chosen_sequence)):                                                  
            if chosen_sequence[k] == '-' and random.randint(0, self.gap_prolong_odds) == 1:  
                indexes_to_insert.append(k)                                                     
        
        while len(indexes_to_insert) > 0:
            idx = indexes_to_insert.pop(-1)
            chosen_sequence = chosen_sequence[:idx] + "-" + chosen_sequence[idx:]

        if coin_flip == 0:
            return type(self)(chosen_sequence, sequence_B, (self.matrix, self.alphabet), self.penalty)
        else:
            return type(self)(sequence_A, chosen_sequence, (self.matrix, self.alphabet), self.penalty)
    
    # Mutation: Add new gaps. [Look for global optimum]
    def offspring_mutate_add_gaps(self, random_gaps = False):
        '''
        Function performes a mutation (by adding gaps) on the chromosome and returns the result as new chromosome (an offspring).
          
        Args:
                random_gaps (boolean): if <i>False</i> (default), add a single new gap. If <i>True</i>, add a random number of
                                      new gaps but not more than <br>
                                      <i>(Length of sequence) // self.new_gaps_limit_factor+1 </i><br>
                                           where <i>self.new_gaps_limit_factor</i> equals 10 by default <br>
                                           (for sequence of length 100, the upper limit of new gaps is 100/10+1=11)
        Returns:
                offspring (Chromosome) : a new Chromosome object
        '''

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        assert len(sequence_A) == len(sequence_B)
        
        coin_flip = random.randint(0, 1) 
        if coin_flip == 0:
            chosen_sequence = sequence_A
        else:
            chosen_sequence = sequence_B
            
        if random_gaps == True:
            ngaps = random.randint(1, (len(chosen_sequence)//self.new_gaps_limit_factor)+1) # Randomly determine a number of gaps to be added.
        else:
            ngaps = 1
            
        for _ in range(ngaps): # Loop to randomly add gaps to the chosen sequence
            idx = random.randint(0, len(chosen_sequence))
            chosen_sequence = chosen_sequence[:idx] + "-" + chosen_sequence[idx:]

        if coin_flip == 0:
            return type(self)(chosen_sequence, sequence_B, (self.matrix, self.alphabet), self.penalty)
        else:
            return type(self)(sequence_A, chosen_sequence, (self.matrix, self.alphabet), self.penalty)

    
    # Mutation: Shuffle gaps. [Look for global optimum]
    def offspring_mutate_shuffle_gaps(self):
        '''
        Function performes a mutation of shuffling random (existing) gaps on the chromosome and returns the result as a new chromosome (an offspring).

        description: 1. Randomly choose one of two sequences in a chromosome (chosen_sequence)
                     2. Each gap in the chosen sequence may be removed with odds 1:gap_shuffle_odds. Count removed gaps.
                     3. For each removed gap, add new one at the random index.
                     4. Create and equalize offspring - a new Chromosome of mutated sequences
                     5. Return an offspring (Chromosome)

        returns: new Chromosome
        '''

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        length = self.get_length()
        gaps_removed = 0
        
        coin_flip = random.randint(0, 1) 
        if coin_flip == 0:
            chosen_sequence = sequence_A
        else:
            chosen_sequence = sequence_B

        for idx in range(length-1, 0, -1):
            if chosen_sequence[idx] == '-' and random.randint(0, self.gap_shuffle_odds) == 1:
                chosen_sequence = chosen_sequence[:idx] + chosen_sequence[idx+1:]
                gaps_removed += 1
                
        while gaps_removed > 0:
            idx = random.randint(0, len(chosen_sequence))
            chosen_sequence = chosen_sequence[:idx] + "-" + chosen_sequence[idx:]
            gaps_removed -= 1

        if coin_flip == 0:
            offspring = type(self)(chosen_sequence, sequence_B, (self.matrix, self.alphabet), self.penalty)
        else:
            offspring = type(self)(sequence_A, chosen_sequence, (self.matrix, self.alphabet), self.penalty)

        return offspring

    # Mutation: Remove gaps. [Look for global optimum]
    def offspring_mutate_remove_gaps(self):
        '''
        Function performes a mutation of removing random (existing) gaps on the chromosome and returns the result as a new chromosome (an offspring).
        
        description: 1. Both sequences are mutated
                     2. Each gap in each sequence may be removed with odds 1:gap_remove_odds.
                     3. Create and equalize offspring - a new Chromosome of mutated sequences
                     4. Return an offspring (Chromosome)

        returns: new Chromosome
        '''

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        length = self.get_length()
        
        for idx in range(length-1, 0, -1):
            if sequence_A[idx] == '-' and random.randint(0, self.gap_remove_odds) == 1:
                sequence_A = sequence_A[:idx] + sequence_A[idx+1:]
            
        for idx in range(length-1, 0, -1):
            if sequence_B[idx] == '-' and random.randint(0, self.gap_remove_odds) == 1:
                sequence_B = sequence_B[:idx] + sequence_B[idx+1:]
            
        offspring = type(self)(sequence_A, sequence_B, (self.matrix, self.alphabet), self.penalty)

        return offspring
    
    # Mutation: Move a single gap by one index to the left or right.
    def offspring_mutate_move_gap(self):
        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        gaps_A = self.gap_indexes(sequence_A)
        gaps_B = self.gap_indexes(sequence_B)

        if len(gaps_A) > 0:
            # randomly select one gap to be moved
            idx = random.choice(gaps_A)
            sequence_A = self.randomly_swap_indexes(sequence_A, idx)

        if len(gaps_B) > 0:
            # randomly select one gap to be moved
            idx = random.choice(gaps_B)
            sequence_B = self.randomly_swap_indexes(sequence_B, idx)

        offspring = type(self)(sequence_A, sequence_B, (self.matrix, self.alphabet), self.penalty)
        
        return offspring
    
        # Mutation: Move a whole section of gaps by one index to the left or right.
    def offspring_mutate_move_section(self):
        '''
        Function performes a mutation of moving whole randomly selected gaps section on the chromosome and returns the result as a new chromosome (an offspring).
        
        description: 1. List all gaps sections indexes for both sequences.
                     2. Randomly select one sequence to be prioritised for mutation.
                     3. If prioritised sequence has no gaps, try with a second sequence. If both have no gaps, return #TODO result of add gaps mutation.
                     3. Create and equalize offspring - a new Chromosome of mutated sequences
                     4. Return an offspring (Chromosome)

        returns: new Chromosome
        '''
        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        gap_sections_A = self.gap_sections_indexes(sequence_A)
        gap_sections_B = self.gap_sections_indexes(sequence_B)

        # Flip a coin to determine which sequence will be prioritised for mutation.
        coin_flip = random.randint(0, 1)

        if len(gap_sections_A) > 0 and (coin_flip == 0 or (coin_flip == 1 and len(gap_sections_B) == 0)):
            # randomly select one gap to be moved
            idx = random.choice(gap_sections_A)
            sequence_A = self.randomly_swap_indexes(sequence_A, idx[0], idx[1])

        elif len(gap_sections_B) > 0 and (coin_flip == 1 or (coin_flip == 0 and len(gap_sections_A) == 0)):
            # randomly select one gap to be moved
            idx = random.choice(gap_sections_B)
            sequence_B = self.randomly_swap_indexes(sequence_B, idx[0], idx[1])
        
        offspring = type(self)(sequence_A, sequence_B, (self.matrix, self.alphabet), self.penalty)
        return offspring

