import numpy as np

# Algoritmo de alineamiento global
def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)
    ptr = np.zeros((n+1, m+1), dtype=int)  # 0: diag, 1: left, 2: up

    # inicialización
    for i in range(1, n+1):
        score[i][0] = score[i-1, 0] + gap
        ptr[i, 0] = 2
    for j in range(1, m+1):
        score[0][j] =  score[0, j-1] + gap
        ptr[0, j] = 3

    # llenado
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1,j] + gap
            best = score[i,j-1] + gap
            score[i][j] = max(diag, up, best)
            if best == diag:
                ptr[i][j] = 1
            elif best == up:
                ptr[i][j] = 2
            else:
                ptr[i][j] = 3
               
    i, j = n, m
    a1, a2 = [], []
    while i > 0 or j > 0:
        if ptr[i, j] == 1:
            a1.append(seq1[i-1]); a2.append(seq2[j-1])
            i -= 1; j -= 1
        elif ptr[i, j] == 2:
            a1.append(seq1[i-1]); a2.append('-')
            i -= 1
        else:
            a1.append('-'); a2.append(seq2[j-1])
            j -= 1

    a1, a2 = ''.join(reversed(a1)), ''.join(reversed(a2))
    return score, ptr, a1, a2, score[n, m]

# Algoritmo de alineamiento local
def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)
    ptr = np.zeros((n+1, m+1), dtype=int) # 1: diag, 2: left, 3: up, 0: stop

    max_score = 0
    max_i = 0
    max_j = 0

    # llenado
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1][j] + gap
            left = score[i, j-1] + gap
            best = max(0, diag, up, left)
            score[i][j] = best
            if best == 0:
                ptr[i, j] = 0
            elif best == diag:
                ptr[i, j] = 1
            elif best == up:
                ptr[i, j] = 2
            else:
                ptr[i, j] = 3
            if best > max_score:
                max_score = best
                max_i, max_j = i, j
    
    i, j = max_i, max_j
    a1, a2 = [], []
    while i > 0 and j > 0 and ptr[i, j] != 0:
        if ptr[i, j] == 1:
            a1.append(seq1[i-1]); a2.append(seq2[j-1])
            i -= 1; j -= 1
        elif ptr[i, j] == 2:
            a1.append(seq1[i-1]); a2.append('-')
            i -= 1
        else:
            a1.append('-'); a2.append(seq2[j-1])
            j -= 1

    a1, a2 = ''.join(reversed(a1)), ''.join(reversed(a2))
    return score, ptr, a1, a2, max_score


def markers(a1, a2):
    marks = ''.join('|' if c1 == c2 else ' ' for c1, c2 in zip(a1, a2))
    identity = sum(c == '|' for c in marks) / max(1, len(marks))
    return marks, identity