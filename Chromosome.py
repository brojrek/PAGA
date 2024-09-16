import random

class Chromosome:
    """
    A class to represent a chromosome - a single unit in the population.

    Attributes
    ----------
    sequence_A : str
        first sequence of the chromosome
    sequence_B : str
        second sequence of the chromosome
    penalty : int
        penalty for gaps in the sequence, used in calculationg score
    matrix : NDarray
        scoring matrix used to calculate alignment score
    alphabet : dict
        alphabet which keys are a characters and values are a corresponding column/row number in matrix
    gap_prolong_odds : int
        odds 1:gap_prolong_odds for each gap in the sequence to be prolonged in offspring_mutate_prolongation mutation
    gap_shuffle_odds : int
        odds 1:gap_shuffle_odds for each already existing gap to be randomly shuffled in sequence.
    gap_remove_odds : int
        odds 1:gap_remove_odds for each gap in the sequence to be removed in offspring_mutaute_remove_gaps mutation
    new_gaps_limit_factor : int
        number of gaps to be added randomly is a number between (1) and (len(sequence)/(new_gaps_limit_factor))

    Methods
    -------
    offspring_mutate_add_gaps(multiple_gaps=False):
        Returns a chromosome with one gap randomly added.
        If multiple_gaps is True, randomly adds multiple gaps instead

    offspring_mutate_prolongation():
        Returns a chromosome with some gaps prolongated

    offspring_mutate_shuffle_gaps():
        Returns a chromosome with some gaps shuffled

    offspring_mutate_move_gap():
        Returns a chromosome with random gap moved

    offspring_mutate_move_section():
        Returns a chromosome with random section of gaps moved

    offspring_mutate_remove_gaps():
        Returns a chromosome with random gaps removed   
    """


    def __init__(self, sequence_A, sequence_B, scoring, penalty, gap_prolong_odds=6, gap_shuffle_odds=3, gap_remove_odds=4, new_gaps_limit_factor=10):
        """
        Constructs all the necessary attributes for the chromosome object.

        Parameters
        ----------
        sequence_A : str
            first sequence of the chromosome
        sequence_B : str
            second sequence of the chromosome
        scoring : list[NDarray, dict]
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


        self.sequence_A = sequence_A
        self.sequence_B = sequence_B

        self.matrix = scoring[0]
        self.alphabet = scoring[1]
        self.penalty = penalty

        self.gap_prolong_odds = gap_prolong_odds #gap_prolong_odds # Odds 1:gap_prolong_odds for each gap in the sequence to be prolonged in offspring_mutate_prolongation mutation
        self.gap_shuffle_odds = gap_shuffle_odds #gap_shuffle_odds # Odds 1:gap_shuffle_odds for each already existing gap to be randomly shuffled in sequence.
        self.gap_remove_odds = gap_remove_odds #gap_remove_odds  # Odds 1:gap_remove_odds
        self.new_gaps_limit_factor = new_gaps_limit_factor #new_gaps_limit_factor # offspring_mutate_add_gaps -> number of gaps to be added randomly is a number between (1) and (len(sequence)/(new_gaps_limit_factor))

        self.equalize()  
        self.calculate_score()

    def __eq__(self, other):
        return self.score == other.score and self.sequence_A == other.sequence_A and self.sequence_B == self.sequence_B

    def __lt__(self, other):
        return self.score < other.score
    
    def __gt__(self, other):
        return self.score > other.score

    #@property
    #def sequence_A(self):
    #    """Get or set the sequence_A. After setting, sequence_A and sequence_B will be automatically equalized."""

    #    return self._sequence_A
    
    #@property
    #def sequence_B(self):
    #    """Get or set the sequence_B. After setting, sequence_A and sequence_B will be automatically equalized."""

    #    return self._sequence_B
    
    @property
    def penalty(self):
        """Get or set the penalty"""
        return self._penalty
    
    @property
    def score(self):
        """Get or set the score"""
        return self._score

    #@sequence_A.setter
    #def sequence_A(self, sequence_A):
    #    if sequence_A == "": raise Exception("Sequence_A cannot be empty")
    #    self._sequence_A = sequence_A
    #    self.equalize()

    #@sequence_B.setter
    #def sequence_B(self, sequence_B):
    #    if sequence_B == "": raise Exception("Sequence_B cannot be empty")
    #    self._sequence_B = sequence_B
    #    self.equalize()

    @penalty.setter
    def penalty(self, penalty):
        if penalty <= 0: raise Exception("Penalty must be greater than 0")
        self._penalty = penalty

    
    def calculate_score(self):
        """
        Calculates and sets score property of the chromosome.
        """

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
        self._score = int(score)

    def get_length(self):
        """
        Gets sequences length. Both sequences in the chromosome are the same length.

        Returns
        -------
        len(self.sequence_A) : int
            Length of sequences on a chromosone 
        """

        return len(self.sequence_A)

    def equalize(self):
        """
        Adds gaps on the end of shorter sequence to match the length of the second sequence.
        """

        len_A = len(self.sequence_A)
        len_B = len(self.sequence_B)
        coin_flip = random.randint(0, 1)
        if len_A > len_B:
            if coin_flip == 0:
                self.sequence_B = self.sequence_B + ("-" * (len_A-len_B))
            else:
                self.sequence_B = ("-" * (len_A-len_B)) + self.sequence_B 
        if len_A < len_B:
            if coin_flip == 0:
                self.sequence_A = self.sequence_A + ("-" * (len_B-len_A))
            else:
                self.sequence_A = ("-" * (len_B-len_A)) + self.sequence_A

        self.gap_reduction()

        assert len(self.sequence_A) == len(self.sequence_B)

    def gap_reduction(self):
        """
        Removes unnecessary gaps (when both sequences have gap on the same index).
        """
        
        for idx in range(self.get_length()-1, -1, -1):
            if self.sequence_A[idx] == '-' and self.sequence_B[idx] == '-':
                self.sequence_A = self.sequence_A[:idx] + self.sequence_A[idx+1:]
                self.sequence_B = self.sequence_B[:idx] + self.sequence_B[idx+1:]

    def gap_indexes(self, seq):
        """
        Gets indexes of the gaps in the given sequence.

        Returns
        -------
        indexes : list[int]
            List of indexes of gaps in the seq.
        """

        indexes = []
        for idx in range(len(seq)-1, 0, -1):
            if seq[idx] == '-':
                indexes.append(idx)
        return indexes
    
    def gap_sections_indexes(self, seq):
        """
        Gets indexes of the gaps in the given sequence.

        Returns
        -------
        indexes : list[touple(int, int)]
            List of indexes of gaps in the seq
        """

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

    def randomly_swap_indexes(self, seq, idx_start, idx_end = None):
        """
        Moves a character on the given index by 1 randomly to left/right. If idx_end is provided,
        a whole section between idx_start and idx_end will be moved.

        Parameters
        ----------
        seq : str
            Sequence in which a character/section will be moved
        idx_start : int
            Start index of the section to be moved
        idx_end : int, optional
            End index of the section to be moved (default is None)

        Returns
        -------
        seq: str
            A sequence with moved gaps.
        """

        if idx_end == None:
            idx_end = idx_start

        if (idx_start == 0 and idx_end == len(seq)):
            pass

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
        """
        Function performes a prolongation mutation (by prolonging gaps) on the chromosome and returns the result as new chromosome (an offspring).
        Every gap existing in the sequence has a fixed odds to be doubled.
        Odds are 1:(self.prolongation_odds)

        Returns
        -------
        Chromosome
            A new chromosome with prolongated gaps.
        """

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
    def offspring_mutate_add_gaps(self, multiple_gaps = False):
        """
        Function performes a mutation (by adding gaps) on the chromosome and returns the result as new chromosome (an offspring).
        
        Parameters
        ----------
        multiple_gaps : boolean, optional
            On False, a signle gap will be added. On True, random number of gaps, up to len(sequence)/new_gaps_limit_factor will be added

        Returns
        -------
        Chromosome
            A new chromosome with added gaps
        """

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        assert len(sequence_A) == len(sequence_B)
        
        coin_flip = random.randint(0, 1) 
        if coin_flip == 0:
            chosen_sequence = sequence_A
        else:
            chosen_sequence = sequence_B
            
        if multiple_gaps == True:
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
        """
        Function performes a mutation of shuffling random (existing) gaps on the chromosome and returns the result as a new chromosome (an offspring).
        Each gap in the randomly chosen sequence may be shuffled with odds 1:gap_shuffle_odds.

        Returns
        -------
        Chromosome
            A new chromosome with shuffled gaps
        """

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
        """
        Function performes a mutation of removing random (existing) gaps on the chromosome and returns the result as a new chromosome (an offspring).
        Each gap in each sequence may be removed with odds 1:gap_remove_odds.

        Returns
        -------
        Chromosome
            A new chromosome with removed gaps
        """

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
        """
        Function performes a mutation of moving whole randomly selected gaps section on the chromosome and returns the result as a new chromosome (an offspring).
        
        Returns
        -------
        Chromosome
            A new chromosome with moved section of gaps
        """

        sequence_A = self.sequence_A
        sequence_B = self.sequence_B

        gap_sections_A = self.gap_sections_indexes(sequence_A)
        gap_sections_B = self.gap_sections_indexes(sequence_B)

        # Flip a coin to determine which sequence will be prioritised for mutation.
        coin_flip = random.randint(0, 1)

        if len(gap_sections_A) > 0 and (coin_flip == 0 or (coin_flip == 1 and len(gap_sections_B) == 0)):
            # randomly select one gap section to be moved
            idx = random.choice(gap_sections_A)
            sequence_A = self.randomly_swap_indexes(sequence_A, idx[0], idx[1])

        elif len(gap_sections_B) > 0 and (coin_flip == 1 or (coin_flip == 0 and len(gap_sections_A) == 0)):
            # randomly select one gap section to be moved
            idx = random.choice(gap_sections_B)
            sequence_B = self.randomly_swap_indexes(sequence_B, idx[0], idx[1])
        
        offspring = type(self)(sequence_A, sequence_B, (self.matrix, self.alphabet), self.penalty)
        return offspring