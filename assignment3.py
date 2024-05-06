import sys


class SequenceAlignment:
    def __init__(self, input_file_path):
        self.input_file_path = input_file_path
        self.str1, self.str2, self.insertion_cost, self.match_cost, self.mismatch_cost = self.read_input()
        self.dp = None
        self.aligned_str1 = None
        self.aligned_str2 = None

    def read_input(self):
        with open(self.input_file_path, 'r') as file:
            str1 = file.readline().strip()
            str2 = file.readline().strip()
            insertion_cost = int(file.readline().strip())
            match_cost = int(file.readline().strip())
            mismatch_cost = int(file.readline().strip())
        return str1, str2, insertion_cost, match_cost, mismatch_cost

    def write_output(self, cost):
        output_file_path = self.input_file_path.replace('input', 'output')
        with open(output_file_path, 'w') as file:
            file.write(str(cost) + '\n')
            file.write(self.aligned_str1.lower() + '\n')
            file.write(self.aligned_str2.lower() + '\n')

    def initialize_dp(self):
        rows = len(self.str1) + 1
        cols = len(self.str2) + 1
        self.dp = [[0] * cols for _ in range(rows)]
        for i in range(rows):
            self.dp[i][0] = i * self.insertion_cost
        for j in range(cols):
            self.dp[0][j] = j * self.insertion_cost

    def SpaceEfficientAlignment(self, s1, s2):
        m, n = len(s1), len(s2)
        array2d = [[0] * 2 for _ in range(m + 5)]

        for i in range(m + 1):
            array2d[i][0] = i * self.insertion_cost

        for j in range(1, n + 1):
            array2d[0][1] = j * self.insertion_cost
            for i in range(1, m + 1):
                array2d[i][1] = min(self.match_cost if s1[i - 1] == s2[j - 1] else self.mismatch_cost + array2d[i - 1][0],
                                    self.insertion_cost + array2d[i - 1][1], self.insertion_cost + array2d[i][0])

            for i in range(m + 1):
                array2d[i][0] = array2d[i][1]

        vec = [array2d[i][1] for i in range(m + 1)]
        return vec

    def fill_dp_table(self):
        for i in range(1, len(self.str1) + 1):
            for j in range(1, len(self.str2) + 1):
                if self.str1[i - 1] == self.str2[j - 1]:
                    match_cost = self.dp[i - 1][j - 1] + self.match_cost
                else:
                    match_cost = self.dp[i - 1][j - 1] + self.mismatch_cost
                insert_cost = self.dp[i][j - 1] + self.insertion_cost
                delete_cost = self.dp[i - 1][j] + self.insertion_cost
                self.dp[i][j] = min(insert_cost, delete_cost, match_cost)

    def traceback(self):
        aligned_str1 = ''
        aligned_str2 = ''
        i, j = len(self.str1), len(self.str2)
        while i > 0 or j > 0:
            if i > 0 and j > 0 and self.dp[i][j] == self.dp[i - 1][j - 1] + (
                    self.match_cost if self.str1[i - 1] == self.str2[j - 1] else self.mismatch_cost):
                aligned_str1 = self.str1[i - 1] + aligned_str1
                aligned_str2 = self.str2[j - 1] + aligned_str2
                i -= 1
                j -= 1
            elif i > 0 and self.dp[i][j] == self.dp[i - 1][j] + self.insertion_cost:  # Deletion
                aligned_str1 = self.str1[i - 1] + aligned_str1
                aligned_str2 = '_' + aligned_str2
                i -= 1
            elif j > 0 and self.dp[i][j] == self.dp[i][j - 1] + self.insertion_cost:  # Insertion
                aligned_str1 = '_' + aligned_str1
                aligned_str2 = self.str2[j - 1] + aligned_str2
                j -= 1
            else:
                aligned_str1 = self.str1[i - 1] + aligned_str1
                aligned_str2 = self.str2[j - 1] + aligned_str2
                i -= 1
                j -= 1

        while i > 0:
            aligned_str1 = self.str1[i - 1] + aligned_str1
            aligned_str2 = '_' + aligned_str2
            i -= 1
        while j > 0:
            aligned_str1 = '_' + aligned_str1
            aligned_str2 = self.str2[j - 1] + aligned_str2
            j -= 1

        self.aligned_str1 = aligned_str1
        self.aligned_str2 = aligned_str2

    def calculate_alignment(self):
        self.initialize_dp()
        self.fill_dp_table()
        self.traceback()
        cost = self.dp[len(self.str1)][len(self.str2)]
        return cost, self.aligned_str1, self.aligned_str2


if __name__ == "__main__":
    # Check if input file is provided
    if len(sys.argv) != 2:
        print("Usage: python3 script.py input_file_name.txt")
        sys.exit(1)

    # Read the input file name from the command line
    input_file_path = sys.argv[1]

    alignment = SequenceAlignment(input_file_path)
    cost, aligned_str1, aligned_str2 = alignment.calculate_alignment()
    alignment.write_output(cost)

    print("Least cost:", cost)
    print("Aligned string 1:", aligned_str1)
    print("Aligned string 2:", aligned_str2)
