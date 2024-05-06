def hirschberg(X, Y, space_cost, match_cost, mismatch_cost):
    def needleman_wunsch(X, Y):
        # Needleman-Wunsch algorithm for alignment
        m = len(X)
        n = len(Y)
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        for i in range(1, m + 1):
            dp[i][0] = dp[i - 1][0] + space_cost

        for j in range(1, n + 1):
            dp[0][j] = dp[0][j - 1] + space_cost

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if X[i - 1] == Y[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1] + match_cost
                else:
                    dp[i][j] = max(dp[i - 1][j - 1] + mismatch_cost,
                                   dp[i - 1][j] + space_cost,
                                   dp[i][j - 1] + space_cost)

        for i in dp:
            print(i)
        print()
        return dp

    def backtrack(dp, X, Y):
        # Backtracking to find the alignment
        m = len(X)
        n = len(Y)
        i = m
        j = n
        align_X = ""
        align_Y = ""

        while i > 0 or j > 0:
            if i > 0 and j > 0 and X[i - 1] == Y[j - 1]:
                align_X = X[i - 1] + align_X
                align_Y = Y[j - 1] + align_Y
                i -= 1
                j -= 1
            elif i > 0 and dp[i][j] == dp[i - 1][j] + space_cost:
                align_X = X[i - 1] + align_X
                align_Y = "-" + align_Y
                i -= 1
            else:
                align_X = "-" + align_X
                align_Y = Y[j - 1] + align_Y
                j -= 1

        return align_X, align_Y

    # Main function
    m = len(X)
    n = len(Y)

    print(X, m, Y)

    if m == 0:
        return n * space_cost, "-" * n, Y
    if n == 0:
        return m * space_cost, X, "-" * m
    if m == 1:
        return (n - 1) * space_cost + (match_cost if X[0] == Y[0] else mismatch_cost), X, Y if X == Y else Y.replace(X[0], '_')
    if n == 1:
        return (m - 1) * space_cost + (match_cost if X[0] == Y[0] else mismatch_cost), X if X == Y else X.replace(Y[0], '_'), Y

    mid = m // 2
    score_left = needleman_wunsch(X[:mid], Y)
    score_right = needleman_wunsch(X[mid:][::-1], Y[::-1])

    max_score = -float('inf')
    split_point = 0

    for j in range(n + 1):
        score = score_left[mid][j] + score_right[m - mid][n - j]
        if score > max_score:
            max_score = score
            split_point = j

    left_score, align_X_left, align_Y_left = hirschberg(X[:mid], Y[:split_point], space_cost, match_cost, mismatch_cost)
    right_score, align_X_right, align_Y_right = hirschberg(X[mid:], Y[split_point:], space_cost, match_cost, mismatch_cost)

    return left_score + right_score, align_X_left + align_X_right, align_Y_left + align_Y_right


def read_input(file_path):
    with open(file_path, 'r') as file:
        X = file.readline().strip()
        Y = file.readline().strip()
        space_cost = int(file.readline().strip())
        match_cost = int(file.readline().strip())
        mismatch_cost = int(file.readline().strip())
    return X, Y, space_cost, match_cost, mismatch_cost


def write_output(score, align_X, align_Y, file_path):
    with open(file_path, 'w') as file:
        file.write(str(score) + '\n')
        file.write(align_X + '\n')
        file.write(align_Y + '\n')


if __name__ == "__main__":
    input_file = "input.txt"
    output_file = "output.txt"

    X, Y, space_cost, match_cost, mismatch_cost = read_input(input_file)
    score, align_X, align_Y = hirschberg(X, Y, space_cost, match_cost, mismatch_cost)

    write_output(score, align_X, align_Y, output_file)
