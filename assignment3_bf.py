############
# Project : String Matching
# Course: CS7200
# Team:
# CHITROJU KODANDA SAIAYYAPPA RAGHAVENDRA SESHA NARAYANACHARYULU
# VARUN GRANDHI
# SUDHEER DANIEL MEGHAVARAM
############


# Read input from file
with open('input.txt', 'r') as file:
    seq1 = file.readline().strip()
    seq2 = file.readline().strip()
    gap_cost = int(file.readline().strip())
    match_cost = int(file.readline().strip())
    mismatch_cost = int(file.readline().strip())


def alignment_cost(seq1, seq2):
    alignments = []

    def backtrack(align1='', align2='', i=0, j=0):
        if i == len(seq1) and j == len(seq2):
            alignments.append((align1, align2))
            return
        if i < len(seq1):
            backtrack(align1 + seq1[i], align2 + '_', i + 1, j)
        if j < len(seq2):
            backtrack(align1 + '_', align2 + seq2[j], i, j + 1)
        if i < len(seq1) and j < len(seq2):
            backtrack(align1 + seq1[i], align2 + seq2[j], i + 1, j + 1)

    backtrack()
    return alignments


def calculate_score(align1, align2, gap_cost, match_cost, mismatch_cost):
    score = 0
    for char1, char2 in zip(align1, align2):
        if char1 == '_' or char2 == '_':
            score += gap_cost
        elif char1 == char2:
            score += match_cost
        else:
            score += mismatch_cost
    return score


def brute_force_alignment(seq1, seq2, gap_cost, match_cost, mismatch_cost):
    alignments = alignment_cost(seq1, seq2)
    min_score = float('inf')
    optimal_alignment = None

    for align1, align2 in alignments:
        score = calculate_score(align1, align2, gap_cost, match_cost, mismatch_cost)
        if score < min_score:
            min_score = score
            optimal_alignment = (align1, align2)

    return optimal_alignment, min_score


# Perform brute force sequence alignment
optimal_alignment, alignment_score = brute_force_alignment(seq1, seq2, gap_cost, match_cost, mismatch_cost)

# Output the optimal alignment and its score
print(alignment_score)
print(optimal_alignment[0])
print(optimal_alignment[1])
