from math import isqrt


def is_prime(n):
    """Return True if n is prime, using trial division."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2

    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def ceil_div(a, b):
    """Return ceil(a / b) for integers b > 0."""
    return -((-a) // b)


def linear_bounds(inequalities, lower_bound=3):
    """Solve finitely many inequalities u*alpha + v >= 0.

    Returns (lower, upper), where upper=None means there is no upper bound.
    Returns None if the inequalities are inconsistent.
    """
    lower = lower_bound
    upper = None

    for u, v in inequalities:
        if u > 0:
            lower = max(lower, ceil_div(-v, u))
        elif u < 0:
            new_upper = v // (-u)
            if upper is None:
                upper = new_upper
            else:
                upper = min(upper, new_upper)
        elif v < 0:
            return None

    if upper is not None and lower > upper:
        return None

    return lower, upper


def pair_has_solution(eta, j):
    """Check whether some alpha >= 3 satisfies the preliminary conditions."""
    A = 4 * j - eta
    B = 12 - j

    inequalities = [
        (A, B - 1),                    # Equation (9)
        (eta - j, j + eta - 12),       # Equation (10)
        (eta, j - 12),                 # S >= 3
        (8 - eta, -j - 4),             # S <= 2*alpha - 1
    ]

    bounds = linear_bounds(inequalities)
    if bounds is None:
        return None

    lower, upper = bounds

    # Integrality of S and t depends only on alpha modulo 4.
    for residue in range(4):
        alpha = lower + ((residue - lower) % 4)

        if upper is not None and alpha > upper:
            continue

        S_numerator = eta * alpha + j
        t_numerator = A * alpha + B

        if S_numerator % 4 != 0 or t_numerator % 4 != 0:
            continue

        S = S_numerator // 4
        t = t_numerator // 4

        if (
            3 <= S <= 2 * alpha - 1
            and t >= 1
            and (eta - j) * alpha + j + eta - 12 >= 0
        ):
            return alpha

    return None


def P(alpha, x):
    """The polynomial P(alpha,x) from Proposition 5.1."""
    return (
        3 * (x * x - 3 * x - 2) * alpha * alpha
        + 4 * (x * x + 2) * alpha
        - (x * x - 3 * x + 2)
    )


def main():
    # First finite check: determine the possible pairs (eta,j).
    possible_pairs = []

    for eta in range(1, 8):
        for j in range(1, 9):
            witness_alpha = pair_has_solution(eta, j)
            if witness_alpha is not None:
                possible_pairs.append((eta, j, witness_alpha))

    print("Possible pairs (eta,j), with one witness alpha:")
    for eta, j, alpha in possible_pairs:
        print(f"eta={eta}, j={j}, witness alpha={alpha}")

    # Second finite check: impose equation (8), integrality of M,
    # and primality of p=M*alpha^2+1.
    surviving_cases = []

    for eta, j, _ in possible_pairs:
        A = 4 * j - eta
        B = 12 - j
        C = 3 * eta
        E = 3 * j - 36 + 4 * eta

        # Equation (8) is
        # M*(A*alpha+B) = C*alpha+E.
        remainder = A * E - C * B

        if remainder == 0:
            # This occurs only for (eta,j)=(6,6).
            # Equation (8) gives M=1, while integrality forces alpha odd,
            # so p=alpha^2+1 is even and greater than 2.
            assert (eta, j) == (6, 6)
            continue

        if A < 0:
            alpha_max = (B - 1) // (-A)
        elif A > 0:
            alpha_max = (abs(remainder) - B) // A
        else:
            continue

        for alpha in range(3, alpha_max + 1):
            S_numerator = eta * alpha + j
            t_numerator = A * alpha + B

            if S_numerator % 4 != 0 or t_numerator % 4 != 0:
                continue

            S = S_numerator // 4
            t = t_numerator // 4

            if not (3 <= S <= 2 * alpha - 1):
                continue
            if t < 1:
                continue
            if (eta - j) * alpha + j + eta - 12 < 0:
                continue

            numerator = C * alpha + E
            denominator = A * alpha + B

            if numerator % denominator != 0:
                continue

            M = numerator // denominator
            if M < 1:
                continue

            p = M * alpha * alpha + 1
            if is_prime(p):
                surviving_cases.append((alpha, M, p, S, t, eta, j))

    surviving_cases.sort()

    print("\nCases surviving equation (8) and primality:")
    print("alpha  M    p    S    t   eta   j")
    for alpha, M, p, S, t, eta, j in surviving_cases:
        print(
            f"{alpha:5d} {M:2d} {p:4d} {S:4d} "
            f"{t:4d} {eta:5d} {j:3d}"
        )

    print("\nFinal root check:")
    for alpha, M, p, S, t, eta, j in surviving_cases:
        print(f"\nalpha={alpha}, p={p}, S={S}")
        found_zero = False

        for r in range(1, alpha + 1):
            s = S - r
            if not (r < s <= alpha):
                continue

            residue_r = P(alpha, r) % p
            residue_s = P(alpha, s) % p

            print(
                f"(r,s)=({r},{s}): "
                f"P(alpha,r)={residue_r}, "
                f"P(alpha,s)={residue_s} mod p"
            )

            if residue_r == 0 and residue_s == 0:
                found_zero = True

        assert not found_zero


if __name__ == "__main__":
    main()
