class Alignment(object):
    def __init__(self, gap_penalty=2, match_score=-1, mismatch_score=3):
        self.seq_a = None
        self.seq_b = None
        self.len_a = None
        self.len_b = None
        self.gap_penalty = gap_penalty
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.separator = '_'

    def match(self, a, b):
        return self.match_score if a == b else self.mismatch_score

    def delete(self, a):
        return self.gap_penalty

    def insert(self, a):
        return self.gap_penalty

    def score(self, aligned_seq_a, aligned_seq_b):
        score = 0
        for a, b in zip(aligned_seq_a, aligned_seq_b):
            if a == b:
                score += self.match_score
            else:
                if a == self.separator or b == self.separator:
                    score += self.gap_penalty
                else:
                    score += self.mismatch_score
        return score

    def map_alignment(self, aligned_seq_a, aligned_seq_b):
        map_b2a = []
        idx = 0
        for x, y in zip(aligned_seq_a, aligned_seq_b):
            if x == y:
                map_b2a.append(idx)
                idx += 1
            elif x == self.separator:
                map_b2a.append(idx)
            elif y == self.separator:
                idx += 1
        return map_b2a


class Needleman(Alignment):
    def __init__(self, *args):
        super(Needleman, self).__init__()
        self.semi_global = False
        self.matrix = None

    def init_matrix(self):
        rows = self.len_a + 1
        cols = self.len_b + 1
        self.matrix = [[0] * cols for i in range(rows)]

    def compute_matrix(self):
        seq_a = self.seq_a
        seq_b = self.seq_b
        len_a = self.len_a
        len_b = self.len_b

        if not self.semi_global:
            for i in range(1, len_a + 1):
                self.matrix[i][0] = self.delete(seq_a[i - 1]) + self.matrix[i - 1][0]
            for i in range(1, len_b + 1):
                self.matrix[0][i] = self.insert(seq_b[i - 1]) + self.matrix[0][i - 1]

        for i in range(1, len_a + 1):
            for j in range(1, len_b + 1):
                score_sub = self.matrix[i - 1][j - 1] + self.match(seq_a[i - 1], seq_b[j - 1])
                score_del = self.matrix[i - 1][j] + self.delete(seq_a[i - 1])
                score_ins = self.matrix[i][j - 1] + self.insert(seq_b[j - 1])
                self.matrix[i][j] = max(score_sub, score_del, score_ins)

    def backtrack(self):
        aligned_seq_a, aligned_seq_b = [], []
        seq_a, seq_b = self.seq_a, self.seq_b

        if self.semi_global:
            last_col_max, val = max(enumerate([row[-1] for row in self.matrix]), key=lambda a: a[1])
            last_row_max, val = max(enumerate([col for col in self.matrix[-1]]), key=lambda a: a[1])

            if self.len_a < self.len_b:
                i, j = self.len_a, last_row_max
                aligned_seq_a = [self.separator] * (self.len_b - last_row_max)
                aligned_seq_b = seq_b[last_row_max:]
            else:
                i, j = last_col_max, self.len_b
                aligned_seq_a = seq_a[last_col_max:]
                aligned_seq_b = [self.separator] * (self.len_a - last_col_max)
        else:
            i, j = self.len_a, self.len_b

        mat = self.matrix

        while i > 0 or j > 0:
            if self.semi_global and (i == 0 or j == 0):
                if i == 0 and j > 0:
                    aligned_seq_a = [self.separator] * j + aligned_seq_a
                    aligned_seq_b = seq_b[:j] + aligned_seq_b
                elif i > 0 and j == 0:
                    aligned_seq_a = seq_a[:i] + aligned_seq_a
                    aligned_seq_b = [self.separator] * i + aligned_seq_b
                break

            if j > 0 and mat[i][j] == mat[i][j - 1] + self.insert(seq_b[j - 1]):
                aligned_seq_a.insert(0, self.separator * len(seq_b[j - 1]))
                aligned_seq_b.insert(0, seq_b[j - 1])
                j -= 1

            elif i > 0 and mat[i][j] == mat[i - 1][j] + self.delete(seq_a[i - 1]):
                aligned_seq_a.insert(0, seq_a[i - 1])
                aligned_seq_b.insert(0, self.separator * len(seq_a[i - 1]))
                i -= 1

            elif i > 0 and j > 0 and mat[i][j] == mat[i - 1][j - 1] + self.match(seq_a[i - 1], seq_b[j - 1]):
                aligned_seq_a.insert(0, seq_a[i - 1])
                aligned_seq_b.insert(0, seq_b[j - 1])
                i -= 1
                j -= 1

            else:
                raise Exception('Backtrack error', i, j, seq_a[i - 2:i + 1], seq_b[j - 2:j + 1])

        return aligned_seq_a, aligned_seq_b

    def align(self, seq_a, seq_b, semi_global=True):
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.len_a = len(self.seq_a)
        self.len_b = len(self.seq_b)

        self.semi_global = semi_global

        self.init_matrix()
        self.compute_matrix()
        return self.backtrack()


class Hirschberg(Alignment):
    def __init__(self):
        super(Hirschberg, self).__init__()
        self.needleman = Needleman()

    def last_row(self, seqa, seqb):
        lena = len(seqa)
        lenb = len(seqb)
        pre_row = [0] * (lenb + 1)
        cur_row = [0] * (lenb + 1)

        for j in range(1, lenb + 1):
            pre_row[j] = pre_row[j - 1] + self.insert(seqb[j - 1])

        for i in range(1, lena + 1):
            cur_row[0] = self.delete(seqa[i - 1]) + pre_row[0]
            for j in range(1, lenb + 1):
                score_sub = pre_row[j - 1] + self.match(seqa[i - 1], seqb[j - 1])
                score_del = pre_row[j] + self.delete(seqa[i - 1])
                score_ins = cur_row[j - 1] + self.insert(seqb[j - 1])
                cur_row[j] = max(score_sub, score_del, score_ins)

            pre_row = cur_row
            cur_row = [0] * (lenb + 1)

        return pre_row

    def align_rec(self, seq_a, seq_b):
        aligned_a, aligned_b = [], []
        len_a, len_b = len(seq_a), len(seq_b)

        if len_a == 0:
            for i in range(len_b):
                aligned_a.append(self.separator * len(seq_b[i]))
                aligned_b.append(seq_b[i])
        elif len_b == 0:
            for i in range(len_a):
                aligned_a.append(seq_a[i])
                aligned_b.append(self.separator * len(seq_a[i]))

        elif len(seq_a) == 1:
            aligned_a, aligned_b = self.needleman.align(seq_a, seq_b)

        else:
            mid_a = int(len_a / 2)

            rowleft = self.last_row(seq_a[:mid_a], seq_b)
            rowright = self.last_row(seq_a[mid_a:][::-1], seq_b[::-1])

            rowright.reverse()

            row = [l + r for l, r in zip(rowleft, rowright)]
            maxidx, maxval = max(enumerate(row), key=lambda a: a[1])

            mid_b = maxidx

            aligned_a_left, aligned_b_left = self.align_rec(seq_a[:mid_a], seq_b[:mid_b])
            aligned_a_right, aligned_b_right = self.align_rec(seq_a[mid_a:], seq_b[mid_b:])
            aligned_a = aligned_a_left + aligned_a_right
            aligned_b = aligned_b_left + aligned_b_right

        return aligned_a, aligned_b

    def align(self, seq_a, seq_b):
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.len_a = len(self.seq_a)
        self.len_b = len(self.seq_b)
        return self.align_rec(self.seq_a, self.seq_b)


def compare(str1, str2):
    seq1 = list(str1)
    seq2 = list(str2)
    for algorithm in [Needleman(), Hirschberg()]:
        a, b = algorithm.align(seq1, seq2)
        print(algorithm.score(a, b))
        print("".join(a))
        print("".join(b))
        print()


if __name__ == "__main__":
    compare("STAY", "SOTY")
