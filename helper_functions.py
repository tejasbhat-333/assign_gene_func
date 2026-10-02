
GAP_PENALTY = -1


def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman-Wunsch algorithm."""
    n, m = len(seq1), len(seq2)

    scores = [[0.0] * (m + 1) for _ in range(n + 1)]
    traceback = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialise first row and column
    for i in range(1, n + 1):
        scores[i][0] = i * GAP_PENALTY
        traceback[i][0] = "up"

    for j in range(1, m + 1):
        scores[0][j] = j * GAP_PENALTY
        traceback[0][j] = "left"

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = (
                scores[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )
            up = scores[i - 1][j] + GAP_PENALTY
            left = scores[i][j - 1] + GAP_PENALTY

            best = max(diagonal, up, left)
            scores[i][j] = best

            # Tie-breaking: diagonal, then up, then left
            if best == diagonal:
                traceback[i][j] = "diagonal"
            elif best == up:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

    # Trace back from the bottom-right cell
    aligned1, aligned2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        direction = traceback[i][j]

        if direction == "diagonal":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif direction == "left":
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), scores[n][m]


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm."""
    n, m = len(seq1), len(seq2)

    scores = [[0.0] * (m + 1) for _ in range(n + 1)]
    traceback = [[None] * (m + 1) for _ in range(n + 1)]

    best_score = 0.0
    best_position = (0, 0)

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = (
                scores[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )
            up = scores[i - 1][j] + GAP_PENALTY
            left = scores[i][j - 1] + GAP_PENALTY

            best = max(0.0, diagonal, up, left)
            scores[i][j] = best

            # Zero means the local alignment starts afresh here
            if best == 0:
                traceback[i][j] = None
            elif best == diagonal:
                traceback[i][j] = "diagonal"
            elif best == up:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

            if best > best_score:
                best_score = best
                best_position = (i, j)

    # Trace back from the highest-scoring cell until score reaches zero
    aligned1, aligned2 = [], []
    i, j = best_position

    while i > 0 and j > 0 and traceback[i][j] is not None:
        direction = traceback[i][j]

        if direction == "diagonal":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif direction == "left":
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), best_score


def scoring_function_simple(aa_i, aa_j):
    """Simple scoring: +1 for a match and -1 for a mismatch."""
    return 1 if aa_i == aa_j else -1