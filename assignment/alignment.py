import numpy as np

import numpy as np

# Algoritmo de alineamiento global
def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)

    # inicialización
    for i in range(n+1):
        score[i][0] = i * gap
    for j in range(m+1):
        score[0][j] = j * gap

    # llenado
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(diag, delete, insert)

    return score

# Algoritmo de alineamiento local
def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)

    max_score = 0
    max_pos = None

    # llenado
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(0, diag, delete, insert)
            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    return score, max_score, max_pos
