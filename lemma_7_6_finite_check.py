from fractions import Fraction
from math import isqrt


def positive_roots(A, B, C):
    """Return the positive rational roots of A*x^2+B*x+C."""
    if A == 0:
        if B == 0:
            return []
        root = Fraction(-C, B)
        return [root] if root > 0 else []

    discriminant = B * B - 4 * A * C
    if discriminant < 0:
        return []

    square_root = isqrt(discriminant)
    if square_root * square_root != discriminant:
        return []

    roots = {
        Fraction(-B + square_root, 2 * A),
        Fraction(-B - square_root, 2 * A),
    }
    return sorted(root for root in roots if root > 0)


def main():
    # Finite check for 5 <= ell <= 15.
    square_discriminant_cases = []

    for ell in range(5, 16):
        lam = 1

        # This is exactly lambda < 3ell/(16-ell).
        while (16 - ell) * lam < 3 * ell:
            A = 2 * lam * (16 - ell) - 6 * ell
            B = 3 * ell * (lam + 3)
            C = ell * (3 - ell - lam)

            roots = positive_roots(A, B, C)
            if roots:
                square_discriminant_cases.append((ell, lam, roots))

            lam += 1

    print("Square-discriminant cases for 5 <= ell <= 15:")
    for ell, lam, roots in square_discriminant_cases:
        roots_text = ", ".join(str(root) for root in roots)
        print(f"ell={ell}, lambda={lam}, positive roots K: {roots_text}")

    print("\nAdmissibility check ell | 4K:")
    for ell, lam, roots in square_discriminant_cases:
        for K in roots:
            admissible = (
                K.denominator == 1
                and (4 * K.numerator) % ell == 0
            )
            print(
                f"ell={ell}, lambda={lam}, K={K}: "
                f"admissible={admissible}"
            )
            assert not admissible

    # The ell=16 case.
    ell_16_cases = []
    for K in range(4, 12, 4):
        if 32 % (3 * K - 1) == 0:
            ell_16_cases.append(K)

    print("\nell=16 cases satisfying 4 | K and 3K-1 | 32:")
    print(ell_16_cases)
    assert ell_16_cases == []

    # Final finite check for 17 <= ell <= 20.
    final_cases = []

    for ell in range(17, 21):
        # Positivity of
        # ((32-2ell)K^2 + 3ell*K - ell)/ell
        # implies K < 3ell/(2ell-32).
        K_max = (3 * ell - 1) // (2 * ell - 32)

        for K in range(2, K_max + 1):
            if (4 * K) % ell != 0:
                continue

            numerator = (
                (32 - 2 * ell) * K * K
                + 3 * ell * K
                - ell
            )

            if numerator <= 0:
                continue

            assert numerator % ell == 0
            mu_over_2 = numerator // ell
            right_side = 6 * K * K - 9 * K - 3 + ell
            remainder = right_side % mu_over_2

            final_cases.append((ell, K, mu_over_2, remainder))

    print("\nCases for 17 <= ell <= 20:")
    print("ell   K   mu/2   remainder")
    for ell, K, mu_over_2, remainder in final_cases:
        print(f"{ell:3d} {K:3d} {mu_over_2:6d} {remainder:11d}")
        assert remainder != 0


if __name__ == "__main__":
    main()
