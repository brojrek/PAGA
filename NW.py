import numpy as np

def nw_alignment(seq_A, seq_B, scoring, gap_penalty):
    nx = len(seq_A)
    ny = len(seq_B)
    matrix, alphabet = scoring[0], scoring[1]
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx, nx + 1)
    F[0,:] = np.linspace(0, -ny, ny + 1)

    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4

    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if not (seq_A[i] == '-' or seq_B[j] == '-'):
                t[0] = F[i,j] + matrix[alphabet[seq_A[i]], alphabet[seq_B[j]]] #matrix[alphabet[seq_A[i]], alphabet[seq_B[j]]]

            t[1] = F[i,j+1] - gap_penalty
            t[2] = F[i+1,j] - gap_penalty
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rseq_A = []
    rseq_B = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rseq_A.append(seq_A[i-1])
            rseq_B.append(seq_B[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rseq_A.append(seq_A[i-1])
            rseq_B.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rseq_A.append('-')
            rseq_B.append(seq_B[j-1])
            j -= 1
    # Reverse the strings.
    rseq_A = ''.join(rseq_A)[::-1]
    rseq_B = ''.join(rseq_B)[::-1]

    score = 0
    assert len(rseq_A) == len(rseq_B)
    for pos in range(0, len(rseq_A)):
        if rseq_A[pos] == '-' or rseq_B[pos] == '-':
            score -= gap_penalty
        else:
            score += matrix[alphabet[rseq_A[pos]], alphabet[rseq_B[pos]]]
    return [[rseq_A, rseq_B], score]
