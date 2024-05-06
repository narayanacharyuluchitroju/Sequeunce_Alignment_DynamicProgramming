import numpy as np
def backtrack_alignment(X, Y, B, match=-1, mismatch=3, gap=2):
    m, n = len(X), len(Y)
    i, j = m, n
    X_aligned = ''
    Y_aligned = ''

    while i > 0 and j > 0:
        if B[i, 1] == B[i-1, 0] + (match if X[i-1] == Y[j-1] else mismatch):
            X_aligned = X[i-1] + X_aligned
            Y_aligned = Y[j-1] + Y_aligned
            i -= 1
            j -= 1
        elif B[i, 1] == B[i, 0] + gap:
            X_aligned = X[i-1] + X_aligned
            Y_aligned = '-' + Y_aligned
            i -= 1
        else:
            X_aligned = '-' + X_aligned
            Y_aligned = Y[j-1] + Y_aligned
            j -= 1

    while i > 0:
        X_aligned = X[i-1] + X_aligned
        Y_aligned = '-' + Y_aligned
        i -= 1
    while j > 0:
        X_aligned = '-' + X_aligned
        Y_aligned = Y[j-1] + Y_aligned
        j -= 1

    return X_aligned, Y_aligned

# Given data
X = "STAY"
Y = "SOTY"
B_last_two_columns = np.array([[6., 8.], [3., 5.], [0., 2.], [2., 3.], [4., 1.]])

# Backtrack alignment
X_aligned, Y_aligned = backtrack_alignment(X, Y, B_last_two_columns)

# Output aligned sequences
print("Sequence 1 after alignment:", X_aligned)
print("Sequence 2 after alignment:", Y_aligned)
