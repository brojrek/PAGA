def read_fasta(filename):
    seq_file = open(filename, "r")
    f=(seq_file.read()).split('\n')
    seq=''
    for i in range(1, len(f)):
        seq=seq+f[i]
    return seq

seq1 = read_fasta('ACE2_[Homo sapiens].fasta')
seq2 = read_fasta('ACE2_[Rhinolophus affinis].fasta')

print(seq1)
print(seq2)