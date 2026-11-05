"""1426 (JB): shared definitions (tapes, moves, strip, windows), written from Paper B Section 2 and SPEC Section 2."""
import secrets

K = 12
MOVES = ["L", "R", "S1", "S2", "S3", "M"]
MIDX = {m: i for i, m in enumerate(MOVES)}
DRIFT = {"L": (-1, 1), "R": (1, -1), "S1": (1, 1), "S2": (1, 1), "S3": (1, 1), "M": (-1, -1)}
SCALE = 88  # fixed-point scale of densities: integer X stands for X / 2^88


def heights(level, word):
    hs = [level]
    for ch in word:
        hs.append(hs[-1] + (1 if ch == "+" else -1))
    return hs


def in_strip(level, word):
    return all(0 <= h < K for h in heights(level, word))


def apply_move(level, word, h, mv):
    """Apply a move to the tape (level, word, head). Returns (word, head) or None if not applicable/legal.
    The head must be interior. Strip legality of the result is checked."""
    n = len(word)
    if not (0 < h < n - 1):
        return None
    c = word[h]
    if mv == "L":
        nw, nh = word, h - 1
    elif mv == "R":
        nw, nh = word, h + 1
    elif mv == "S1":
        if c != "+": return None
        nw, nh = word[:h] + "+-+" + word[h + 1:], h + 1
    elif mv == "S2":
        if c != "-": return None
        nw, nh = word[:h] + "+--" + word[h + 1:], h + 1
    elif mv == "S3":
        if c != "-": return None
        nw, nh = word[:h] + "--+" + word[h + 1:], h + 1
    elif mv == "M":
        if c != "+": return None
        a, b = word[h - 1], word[h + 1]
        if a == "+" and b == "+": return None
        s = (1 if a == "+" else -1) + (1 if b == "+" else -1) + 1
        x = "+" if s == 1 else "-"
        assert s in (1, -1)
        nw, nh = word[:h - 1] + x + word[h + 2:], h - 1
    else:
        raise ValueError(mv)
    if not in_strip(level, nw):
        return None
    return nw, nh


def legal_moves_window(level, word):
    """Legal moves at a window (head at the middle letter), in the order L,R,S1,S2,S3,M (SPEC 2.3)."""
    r = (len(word) - 1) // 2
    out = []
    for mv in MOVES:
        if apply_move(level, word, r, mv) is not None:
            out.append(mv)
    return out


def windows(n):
    """All (level, word) of n letters in the strip, levels increasing, words lexicographic with + before -."""
    res = []
    def rec(level, h, pre):
        if len(pre) == n:
            res.append((level, "".join(pre))); return
        for ch, s in (("+", 1), ("-", -1)):
            if 0 <= h + s < K:
                pre.append(ch); rec(level, h + s, pre); pre.pop()
    for lev in range(K):
        rec(lev, lev, [])
    return res


def code(word):
    v = 0
    for ch in word:
        v = (v << 1) | (1 if ch == "+" else 0)
    return v


def decode(c, n):
    return "".join("+" if (c >> (n - 1 - j)) & 1 else "-" for j in range(n))


def is_probable_prime(n):
    if n < 2: return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in small:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in small:  # deterministic for n < 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True


def random_prime(bits=62):
    while True:
        c = secrets.randbits(bits) | (1 << (bits - 1)) | 1
        if is_probable_prime(c):
            return c
