from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from tqdm.auto import tqdm


# Safe exponential function
def safe_exp(x) -> np.ndarray:
    return np.exp(np.minimum(x, 700))
    # return np.exp(x)


# Function F (Vectorized)
def F(x: np.ndarray) -> np.ndarray:
    return (
        255
        * safe_exp(-safe_exp(-1000 * x))
        * np.abs(x) ** (safe_exp(-safe_exp(1000 * (x - 1))))
    )


def L(v: int, Table: Dict[str, np.ndarray]) -> np.ndarray:
    S = Table["S"]
    R0 = Table["R0"]
    P = Table["P"]
    C20 = Table["C20"]
    C10 = Table["C10"]
    A4 = Table["A4"]
    A1000 = Table["A1000"]
    B = Table["B"]

    left_term_1 = 1 / 10 - (1 / 40) * np.cos(20 * np.arccos(R0)) * np.cos(25 * P)
    left_term_2 = (
        4 * v**2
        - 13 * v
        + 11
        + np.cos(7 * S + v * S)
        + 20 * safe_exp(-safe_exp(-70 * (C20 - 1 / 2)))
        + 20 * safe_exp(-safe_exp(-10 * (C10 - 1 / 2)))
    )

    return left_term_1 * left_term_2 * A4 * A1000 + B


def H(v: int, Table: Dict[str, np.ndarray]) -> np.ndarray:

    def f(r: int) -> np.ndarray:
        return (
            (1 - Table["A1000"][r])
            * (1 - safe_exp(-safe_exp(-1000 * (r - 0.5))) * Table["U"][r])
            * (1 - 1.25 * Table["A4"][r])
        )

    def g(v: int, s: int) -> np.ndarray:
        return (
            (5 - 3 * (v - 1) ** 2 + Table["W_xy"])
            / 10
            * safe_exp(-safe_exp(-np.abs(71 / 10 - 10 * Table["P"][s])))
            * Table["U"][s]
        ) + Table["A1000"][s] * (1 - Table["U"][s]) * L_val[s]

    L_val = L(v, Table)
    f_cumprod = np.ones_like(Table["W_xy"])
    H_val = np.zeros_like(Table["W_xy"])
    for s in tqdm(range(1, 31)):
        f_cumprod *= f(s - 1)
        H_val += f_cumprod * g(v, s)

    return H_val


# Vectorized RGB function
def compute_RGB(Table: Dict[str, np.ndarray]) -> np.ndarray:
    R = F(H(0, Table))
    G = F(H(1, Table))
    B = F(H(2, Table))

    img = np.stack([R, G, B], axis=-1).astype(np.uint8)
    return img


# Generate image with vectorized operations
def generate_image(
    width: int = 2000, height: int = 1200, denominator: int | None = None
) -> np.ndarray:
    X, Y = _create_meshgrid(width, height, denominator)

    Table = _prepare_tables(width, height, X, Y)

    print("Generating image using vectorized RGB computation...")
    image = compute_RGB(Table)
    return image


def _create_meshgrid(
    width: int, height: int, denominator: int
) -> Tuple[np.ndarray, np.ndarray]:
    x_denominator = denominator or width // 2
    x_vals = (np.arange(width) + 1 - width // 2) / x_denominator
    y_denominator = denominator or height // 2
    y_vals = (height // 2 - np.arange(height)) / y_denominator
    X, Y = np.meshgrid(x_vals, y_vals)
    return X, Y


def _prepare_tables(
    width: int, height: int, X: np.ndarray, Y: np.ndarray
) -> Dict[str, np.ndarray]:
    print("Prepare 'S' Table")
    S = _prepare_s_range(width, height, n_range=31)
    print("Prepare 'P' Table")
    P = _prepare_p_table(X, Y, S)
    print("Prepare 'Q' Table")
    Q = _prepare_q_table(X, Y, S)
    print("Prepare 'W' Table")
    W = _prepare_w_table(X, Y, width, height, n_range=40)
    W_xy = W.sum(axis=0)
    print("Prepare 'A' Table")
    A4 = _prepare_a_table(4, S, P, Q)
    A1000 = _prepare_a_table(1000, S, P, Q)
    print("Prepare 'R' Table")
    R0 = _prepare_r_table(0, P, Q)
    R1 = _prepare_r_table(1, P, Q)
    print("Prepare 'U' Table")
    U = _prepare_u_table(X, Y, S, P, Q, R0, R1)
    print("Prepare 'C' Table")
    C10 = _prepare_c_table(10, R0, P, Q, W_xy)
    C20 = _prepare_c_table(20, R0, P, Q, W_xy)
    print("Prepare 'B' Table")
    B = _prepare_b_table(R0, P)

    return {
        "S": S,
        "P": P,
        "Q": Q,
        "W_xy": W_xy,
        "A4": A4,
        "A1000": A1000,
        "R0": R0,
        "R1": R1,
        "U": U,
        "C10": C10,
        "C20": C20,
        "B": B,
    }


def _prepare_s_range(width: int, height: int, n_range: int) -> np.ndarray:
    return np.tile(np.arange(n_range).reshape(n_range, 1, 1), (1, height, width))


def _prepare_p_table(X: np.ndarray, Y: np.ndarray, S: np.ndarray) -> np.ndarray:
    return np.arctan(
        np.tan(
            2 * np.sin(5 * S) * X
            - 2 * np.cos(5 * S) * Y
            + 3 * np.cos(5 * S)
            + (3 * np.cos(14 * X - 19 * Y + 5 * S)) / 200
        )
    )


def _prepare_q_table(X: np.ndarray, Y: np.ndarray, S: np.ndarray) -> np.ndarray:
    return np.arctan(
        np.tan(
            2 * (np.cos(5 * S) * X + np.sin(5 * S) * Y + 2 * np.cos(4 * S))
            + (3 * np.cos(18 * X + 15 * Y + 4 * S)) / 200
        )
    )


def _prepare_w_table(
    X: np.ndarray, Y: np.ndarray, width: int, height: int, n_range: int
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


def _prepare_a_table(v: int, S: np.ndarray, P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    s_const = -safe_exp(-1000 * (S - 1 / 2))
    weighted = -safe_exp(
        v
        * (
            (5 / 4) * (1 - P) * (Q**2)
            + (P) ** 2
            - 11 / 20
            + (np.arctan(100 * (v - 100)) / (10 * np.pi))
        )
    )
    penalty = -safe_exp(v * (Q**2 + P**2 - 1))
    return safe_exp(s_const + weighted + penalty)


def _prepare_u_table(
    X: np.ndarray,
    Y: np.ndarray,
    S: np.ndarray,
    P: np.ndarray,
    Q: np.ndarray,
    R0: np.ndarray,
    R1: np.ndarray,
) -> np.ndarray:
    left = -safe_exp(
        -100
        * (
            P
            - 57 / 100
            - (3 / 20 + np.cos(7 * Q + 2 * S) / 10)
            * np.cos(
                (10 + 3 * np.cos(14 * S)) * np.arccos(R0)
                + 3 / 10 * np.cos(45 * X + 47 * Y + np.cos(17 * X))
                + 2 * np.cos(5 * S)
            )
        )
    )
    right = -safe_exp(1000 * (P - 18 / 25 + (3 / 2 * Q) ** 8))
    M = safe_exp(left + right)

    left = -safe_exp(
        100
        * (
            P
            - 37 / 50
            - (3 / 20 + np.cos(8 * Q + 5 * S) / 10)
            * np.cos(
                (10 + 3 * np.cos(16 * S)) * np.arccos(R1)
                + 3 / 10 * np.cos(38 * X - 47 * Y + np.cos(19 * X))
                + 2 * np.cos(4 * S)
            )
        )
    )
    right = -safe_exp(-1000 * (P - 71 / 100 - (3 / 2 * Q) ** 8))
    N = safe_exp(left + right)

    return 1 - (1 - M) * (1 - N)


def _prepare_r_table(t: int, P: np.ndarray, Q: np.ndarray) -> np.ndarray:

    def E(t, P, Q):
        left_term = 1000 / np.sqrt(20) * Q
        center_term = np.sqrt(5 * np.abs(20 - 20 * (1 - 2 * t) * P - 27 * t))
        right_term = (
            1 + 50 * np.sqrt(np.abs(4 * (200 - (20 * (1 - 2 * t) * P + 27 * t) ** 2)))
        ) ** -1

        return left_term + center_term + right_term

    E_t_s = E(t, P, Q)
    return E_t_s * safe_exp(-safe_exp(1000 * (np.abs(E_t_s) - 1)))


def _prepare_c_table(
    v: int,
    R: np.ndarray,
    P: np.ndarray,
    Q: np.ndarray,
    W_xy: np.ndarray,
) -> np.ndarray:
    cos_term = np.cos(10 * np.arccos(R)) * np.cos(25 / 2 * P)
    sin_term = np.sin(10 * np.arccos(R)) * np.sin(25 / 2 * P)

    left = -(
        safe_exp(v * (cos_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (cos_term + 7 / 10 + W_xy / 5))
        + safe_exp(v * (sin_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (sin_term + 7 / 10 + W_xy / 5))
    )
    right = -safe_exp(3 / 2 * ((Q) ** 2 + (P - 1 / 4) ** 2 - 21 / 50 + W_xy / 5))

    return safe_exp(left + right)


def _prepare_b_table(R: np.ndarray, P: np.ndarray) -> np.ndarray:
    return safe_exp(
        -safe_exp(-70 * (np.cos(20 * np.arccos(R)) * np.cos(25 * P) - 47 / 50))
    )


# Run the optimized function
if __name__ == "__main__":
    width, height = 2000, 1200
    image = generate_image(width=width, height=height, denominator=height // 2)
    plt.imsave("strawberries_optimized_250317_1.png", image)
    plt.imshow(image)
    plt.axis("off")
    plt.show()
