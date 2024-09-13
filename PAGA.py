import numpy
from Chromosome import *
from Population import *

import sys
import time

def read_fasta(filename):
    '''
    Reads a sequence from .fasta file containing exactly one sequence. First line is a header.

    filename: string representing a path or name of the .fasta file

    returns: string representing sequence read
    '''

    seq_file = open(filename, "r")
    f=(seq_file.read()).split('\n')
    seq=''
    for i in range(1, len(f)):
        seq=seq+f[i]
    return seq

def load_score_matrix(filename):
    '''
    #TODO documentation
    '''
    f = open(filename,'r')
    content = f.read()
    content = content.splitlines()
    n = len(content)-1
    matrix=numpy.zeros((n,n),'i') #d - float i - int  
    f.close()
    alphabet={}
    row=content[0].split()
    for j in range(0,len(row)):
        alphabet[row[j]]=j

    for i in range(1,len(content)):
        row=content[i].split()
        for j in range(1,len(row)):
            if row[j]!='':
                matrix[i-1,j-1] = float(row[j])
                matrix[j-1,i-1] = float(matrix[i-1,j-1])

    return matrix, alphabet

seq1 = read_fasta("ACE2_[Homo_sapiens].fasta")
seq2 = read_fasta("ACE2_[Rhinolophus_affinis].fasta")
print(len(seq1))
print(len(seq2))
ch = Chromosome(seq1, seq2, load_score_matrix("BLOSUM62.txt"), 5)
population = Population(ch, 100)


start_time = time.time()
for n in range(50):
    print('Iteration: ' + str(n) + ', ' + str(round(time.time() - start_time, 2)) + ' seconds passed. Highscore: ' + str(population.population[0].score))
    sys.stdout.flush()
    time.sleep(1)
    population.new_generation(population.selection_tournament())
    print(len(population))
end_time = time.time()
print(population.population[0].score, population.population[0].sequence_A, population.population[0].sequence_B)
#print([x.score for x in population.population])
print("Total execution time: " + str(round(end_time - start_time, 2)) + " seconds.")