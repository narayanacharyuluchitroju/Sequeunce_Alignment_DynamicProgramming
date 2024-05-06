def hirschberg(X, Y, match_score, mismatch_score, gap_penalty):
    m, n = len(X), len(Y)

    # Initialize the DP table
    dp = [[0] * (n + 1) for _ in range(2)]

    # Perform the alignment
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if X[i - 1] == Y[j - 1]:
                dp[i % 2][j] = dp[(i - 1) % 2][j - 1] + match_score
            else:
                dp[i % 2][j] = max(dp[(i - 1) % 2][j] + gap_penalty,
                                   dp[i % 2][j - 1] + gap_penalty,
                                   dp[(i - 1) % 2][j - 1] + mismatch_score)

    print(dp)
    # Traceback to find the alignment
    align_X, align_Y = [], []
    i, j = m, n
    while i > 0 and j > 0:
        if X[i - 1] == Y[j - 1]:
            align_X.insert(0, X[i - 1])
            align_Y.insert(0, Y[j - 1])
            i -= 1
            j -= 1
        elif dp[(i - 1) % 2][j] + gap_penalty == dp[i % 2][j]:
            align_X.insert(0, X[i - 1])
            align_Y.insert(0, '-')
            i -= 1
        elif dp[i % 2][j - 1] + gap_penalty == dp[i % 2][j]:
            align_X.insert(0, '-')
            align_Y.insert(0, Y[j - 1])
            j -= 1
        else:
            align_X.insert(0, X[i - 1])
            align_Y.insert(0, Y[j - 1])
            i -= 1
            j -= 1

    while i > 0:
        align_X.insert(0, X[i - 1])
        align_Y.insert(0, '-')
        i -= 1

    while j > 0:
        align_X.insert(0, '-')
        align_Y.insert(0, Y[j - 1])
        j -= 1

    return align_X, align_Y

# Example usage:
X = "STAY"
Y = "SOTY"
match_score = 2
mismatch_score = -1
gap_penalty = 3

align_X, align_Y = hirschberg(X, Y, match_score, mismatch_score, gap_penalty)
print("".join(align_X))
print("".join(align_Y))
