from typing import Literal, Tuple

import matplotlib.pyplot as plt
import numpy as np
from tqdm.auto import tqdm


# Safe exponential function
def safe_exp(x):
    # return np.exp(np.minimum(x, 700))
    return np.exp(x)


# Precompute P_s(x, y) and Q_s(x, y) for all s in [1, 30]
def precompute_PQ(
    width: int, height: int, X: np.ndarray, Y: np.ndarray, n_range: int
) -> Tuple[np.ndarray, np.ndarray]:
    P_table = np.zeros((n_range, height, width))
    Q_table = np.zeros((n_range, height, width))

    for s in range(1, n_range + 1):
        P_table[s - 1] = np.arctan(
            np.tan(
                2 * np.sin(5 * s) * X
                - 2 * np.cos(5 * s) * Y
                + 3 * np.cos(5 * s)
                + (3 * np.cos(14 * X - 19 * Y + 5 * s)) / 200
            )
        )
        Q_table[s - 1] = np.arctan(
            np.tan(
                2 * (np.cos(5 * s) * X + np.sin(5 * s) * Y + 2 * np.cos(4 * s))
                + (3 * np.cos(18 * X + 15 * Y + 4 * s)) / 200
            )
        )

    return P_table, Q_table


# Precompute W(x, y) for the entire image
def precompute_W(
    width: int, height: int, X: np.ndarray, Y: np.ndarray, n_range: int
) -> np.ndarray:
    W_table = np.zeros((n_range, height, width))

    for s in range(1, n_range + 1):
        left = np.cos(
            28**s * 25**-s * (np.cos(2 * s) * X + np.sin(2 * s) * Y) + 2 * np.sin(5 * s)
        )
        right = np.cos(
            28**s * 25**-s * (np.cos(2 * s) * Y - np.sin(2 * s) * X) + 2 * np.sin(6 * s)
        )
        exp_exponent = -3 * (left**2 * right**2 + (-97 / 100))
        W_table[s - 1] = safe_exp(-safe_exp(exp_exponent))

    return W_table


# Function F (Vectorized)
def F(x: np.ndarray) -> np.ndarray:
    return (
        255
        * safe_exp(-safe_exp(-1000 * x))
        * np.abs(x) ** (safe_exp(-safe_exp(1000 * (x - 1))))
    )


def A(v: int, s: int, P_table: np.ndarray, Q_table: np.ndarray) -> np.ndarray:
    P_s = P_table[s - 1]
    Q_s = Q_table[s - 1]

    s_const = -safe_exp(-1000 * (s - 1 / 2))
    weighted = -safe_exp(
        v
        * (
            (5 / 4) * (1 - P_s) * (Q_s**2)
            + (P_s**2)
            - 11 / 20
            + (np.arctan(100 * (v - 100)) / (10 * np.pi))
        )
    )
    penalty = -safe_exp(v * (Q_s**2 + P_s**2 - 1))
    return safe_exp(s_const + weighted + penalty)


def U(
    s: int, X: np.ndarray, Y: np.ndarray, P_table: np.ndarray, Q_table: np.ndarray
) -> np.ndarray:
    P_s = P_table[s - 1]
    Q_s = Q_table[s - 1]

    def _inner_fn(mode: Literal["M", "N"]):
        R_v = R(0 if mode == "M" else 1, s, P_table, Q_table)

        if mode == "M":
            left_exponent = -safe_exp(
                -100
                * (
                    P_s
                    - 57 / 100
                    - (3 / 20 + np.cos(7 * Q_s + 2 * s) / 10)
                    * np.cos(
                        (10 + 3 * np.cos(14 * s)) * np.arccos(R_v)
                        + 3 / 10 * np.cos(45 * X + 47 * Y + np.cos(17 * X))
                        + 2 * np.cos(5 * s)
                    )
                )
            )
            right_exponent = -safe_exp(1000 * (P_s - 18 / 25 + (3 / 2 * Q_s) ** 8))
        else:
            left_exponent = -safe_exp(
                100
                * (
                    P_s
                    - 37 / 50
                    - (3 / 20 + np.cos(8 * Q_s + 5 * s) / 10)
                    * np.cos(
                        (10 + 3 * np.cos(16 * s)) * np.arccos(R_v)
                        + 3 / 10 * np.cos(38 * X - 47 * Y + np.cos(19 * X))
                        + 2 * np.cos(4 * s)
                    )
                )
            )
            right_exponent = -safe_exp(-1000 * (P_s - 71 / 100 - (3 / 2 * Q_s) ** 8))
        return safe_exp(left_exponent + right_exponent)

    M_s = _inner_fn(mode="M")
    N_s = _inner_fn(mode="N")

    return 1 - (1 - M_s) * (1 - N_s)


def R(t: int, s: int, P_table: np.ndarray, Q_table: np.ndarray) -> np.ndarray:

    def E(t, s, P_table, Q_table):
        P_s = P_table[s - 1]
        Q_s = Q_table[s - 1]

        left_term = 1000 / np.sqrt(20) * Q_s
        center_term = np.sqrt(5 * np.abs(20 - 20 * (1 - 2 * t) * P_s - 27 * t))
        right_term = (
            1 + 50 * np.sqrt(np.abs(4 * (200 - (20 * (1 - 2 * t) * P_s + 27 * t) ** 2)))
        ) ** -1

        return left_term + center_term + right_term

    E_t_s = E(t, s, P_table, Q_table)
    return E_t_s * safe_exp(-safe_exp(1000 * (np.abs(E_t_s) - 1)))


def L(
    v: int, s: int, P_table: np.ndarray, Q_table: np.ndarray, W_xy: np.ndarray
) -> np.ndarray:
    R_0_s = R(t=0, s=s, P_table=P_table, Q_table=Q_table)
    P_s = P_table[s - 1]
    C_20_s = C(v=20, s=s, R_s=R_0_s, P_table=P_table, Q_table=Q_table, W_xy=W_xy)
    C_10_s = C(v=10, s=s, R_s=R_0_s, P_table=P_table, Q_table=Q_table, W_xy=W_xy)
    A_4_s = A(v=4, s=s, P_table=P_table, Q_table=Q_table)
    A_1000_s = A(v=1000, s=s, P_table=P_table, Q_table=Q_table)
    B_s = B(s=s, R_s=R_0_s, P_table=P_table)

    left_term_1 = 1 / 10 - (1 / 40) * np.cos(20 * np.arccos(R_0_s)) * np.cos(25 * P_s)
    left_term_2 = (
        4 * v**2
        - 13 * v
        + 11
        + np.cos(7 * s + v * s)
        + 20 * safe_exp(-safe_exp(-70 * (C_20_s - 1 / 2)))
        + 20 * safe_exp(-safe_exp(-10 * (C_10_s - 1 / 2)))
    )

    return left_term_1 * left_term_2 * A_4_s * A_1000_s + B_s


def C(
    v: int,
    s: int,
    R_s: np.ndarray,
    P_table: np.ndarray,
    Q_table: np.ndarray,
    W_xy: np.ndarray,
) -> np.ndarray:
    P_s = P_table[s - 1]
    Q_s = Q_table[s - 1]

    cos_term = np.cos(10 * np.arccos(R_s)) * np.cos(25 / 2 * P_s)
    sin_term = np.sin(10 * np.arccos(R_s)) * np.sin(25 / 2 * P_s)

    left_exponent = -(
        safe_exp(v * (cos_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (cos_term + 7 / 10 + W_xy / 5))
        + safe_exp(v * (sin_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (sin_term + 7 / 10 + W_xy / 5))
    )
    right_exponent = -safe_exp(
        3 / 2 * ((Q_s) ** 2 + (P_s - 1 / 4) ** 2 - 21 / 50 + W_xy / 5)
    )

    return safe_exp(left_exponent + right_exponent)


def B(s: int, R_s: np.ndarray, P_table: np.ndarray) -> np.ndarray:
    P_s = P_table[s - 1]
    return safe_exp(
        -safe_exp(-70 * (np.cos(20 * np.arccos(R_s)) * np.cos(25 * P_s) - 47 / 50))
    )


# H function (Vectorized version)
def H(
    v: int,
    X: np.ndarray,
    Y: np.ndarray,
    P_table: np.ndarray,
    Q_table: np.ndarray,
    W_xy: np.ndarray,
    n_range: int = 30,
) -> np.ndarray:
    sum_term = np.zeros_like(X)

    prod_cache = {}
    for s in tqdm(range(1, n_range + 1)):
        prod_term = np.ones_like(X)

        r = s - 1
        if prod_cache.get(r, None) is not None:
            prod_term *= prod_cache[r]

        prod_term *= 1 - A(v=1000, s=r, P_table=P_table, Q_table=Q_table)
        U_r = U(s=r, X=X, Y=Y, P_table=P_table, Q_table=Q_table)
        prod_term *= 1 - safe_exp(-safe_exp(-1000 * (r - 1 / 2))) * U_r
        prod_term *= 1 - 5 / 4 * A(v=4, s=r, P_table=P_table, Q_table=Q_table)

        U_s = U(s=s, X=X, Y=Y, P_table=P_table, Q_table=Q_table)
        L_v_s = L(v=v, s=s, P_table=P_table, Q_table=Q_table, W_xy=W_xy)
        A_1000_s = A(v=1000, s=s, P_table=P_table, Q_table=Q_table)

        right_term = (5 - 3 * (v - 1) ** 2 + W_xy) / 10
        right_term *= safe_exp(-safe_exp(-np.abs(71 / 10 - 10 * P_table[s - 1]))) * U_s
        right_term += A_1000_s * (1 - U_s) * L_v_s

        sum_term += prod_term + right_term

    return sum_term


# Vectorized RGB function
def compute_RGB(
    X: np.ndarray,
    Y: np.ndarray,
    P_table: np.ndarray,
    Q_table: np.ndarray,
    W_xy: np.ndarray,
    n_range: int = 30,
) -> np.ndarray:
    R = F(H(0, X, Y, P_table, Q_table, W_xy, n_range))
    G = F(H(1, X, Y, P_table, Q_table, W_xy, n_range))
    B = F(H(2, X, Y, P_table, Q_table, W_xy, n_range))

    img = np.stack([R, G, B], axis=-1).astype(np.uint8)
    return img


# Generate image with vectorized operations
def generate_image(
    width: int = 2000, height: int = 1200, denominator: int | None = None
):
    x_denominator = denominator or width // 2
    x_vals = (np.arange(width) + 1 - width // 2) / x_denominator
    y_denominator = denominator or height // 2
    y_vals = (height // 2 - np.arange(height)) / y_denominator
    X, Y = np.meshgrid(x_vals, y_vals)

    print("Precomputing P_s, Q_s tables...")
    P_table, Q_table = precompute_PQ(width, height, X, Y, n_range=30)

    print("Precomputing W(x, y) table...")
    W_table = precompute_W(width, height, X, Y, n_range=40)
    W_xy = W_table.sum(axis=0)

    print("Generating image using vectorized RGB computation...")
    image = compute_RGB(
        X=X, Y=Y, P_table=P_table, Q_table=Q_table, W_xy=W_xy, n_range=30
    )

    return image


# Run the optimized function
if __name__ == "__main__":
    width, height = 2000, 1200
    image = generate_image(width=width, height=height, denominator=height // 2)
    plt.imsave("strawberries_optimized_250317.png", image)
    plt.imshow(image)
    plt.axis("off")
    plt.show()
