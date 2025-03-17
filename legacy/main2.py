"""
Strawberries by Hamid Naderi Yeganeh

2000 X 1200 Image.
For m = 1, 2, 3, ..., 2000 and n = 1, 2, 3, ..., 1200,
the color of the pixel or the row n and the column m is ::

RGB(
    F(
        H0(
            (m - 1000) / 600, (601 - n) / 600
        )
    ),
    F(
        H1(
            (m - 1000) / 600, (601 - n) / 600
        )
    ),
    F(
        H2(
            (m - 1000) / 600, (601 - n) / 600
        )
    ),
)
where

F(x) = 255e^{-e^{-1000x}} \cdot |x|^{e^{-e^{1000(x-1)}}}

H_v(x, y) =
\sum_{s=1}^{30} \left[
    \left(
        \prod_{r=0}^{s-1}
        \left( 1 - A_{1000,r}(x,y) \right)
        \left( 1 - e^{-e^{-1000\left(r - \frac{1}{2} \right)}} U_r(x,y) \right)
        \left( 1 - \frac{5}{4} A_{4,r}(x, y) \right)
    \right)
\right.
\]
\[
\left.
    \times \left(
        \frac{5 - 3(v - 1)^2 + W(x,y)}{10}
        e^{-e^{-\left|\frac{71}{10} - 10P_s(x,y)\right|}} U_s(x,y)
        + A_{1000,s}(x,y) (1 - U_s(x,y)) L_{v,s}(x,y)
    \right)
\right].

L_{v,s}(x,y) =
\left( \frac{1}{10} - \frac{1}{40} \cos \left( 20 \arccos(R_{0,s}(x,y)) \right) \cos(25 P_s(x,y)) \right)
\]
\[
\times \left( 4v^2 - 13v + 11 + \cos(7s + vs) + 20 e^{-e^{-70 (c_{20,s}(x,y) - \frac{1}{2})}} + 20 e^{-e^{-10 (c_{10,s}(x,y) - \frac{1}{2})}} \right)
\]
\[
\times A_{4,s}(x,y) A_{1000,s}(x,y) + B_s(x,y).

C_{v,s}(x,y) =
e^{-e^{v \left( \cos \left( 10 \arccos(R_{0,s}(x,y)) \right) \cos \left( \frac{25}{2} P_s(x,y) \right) - \frac{7}{10} - \frac{W(x,y)}{5} \right)}}
\]
\[
- e^{-e^{v \left( \cos \left( 10 \arccos(R_{0,s}(x,y)) \right) \cos \left( \frac{25}{2} P_s(x,y) \right) + \frac{7}{10} + \frac{W(x,y)}{5} \right)}}
\]
\[
- e^{-e^{v \left( \sin \left( 10 \arccos(R_{0,s}(x,y)) \right) \sin \left( \frac{25}{2} P_s(x,y) \right) - \frac{7}{10} - \frac{W(x,y)}{5} \right)}}
\]
\[
- e^{-e^{v \left( \sin \left( 10 \arccos(R_{0,s}(x,y)) \right) \sin \left( \frac{25}{2} P_s(x,y) \right) + \frac{7}{10} + \frac{W(x,y)}{5} \right)}}
\]
\[
\times e^{-e^{-\frac{3}{2} \left( (Q_s(x,y))^2 + \left( P_s(x,y) - \frac{1}{4} \right)^2 - \frac{21}{50} + \frac{W(x,y)}{5} \right)}}.


B_s(x,y) = e^{-e^{-70 \left( \cos \left( 20 \arccos(R_{0,s}(x,y)) \right) \cos \left( 25 P_s(x,y) \right) - \frac{47}{50} \right)}}

A_{v,s}(x,y) =
e^{-e^{-1000(s - \frac{1}{2})} - e^{v \left( \frac{5}{4} (1 - P_s(x,y)) (Q_s(x,y))^2 + (P_s(x,y))^2 - \frac{11}{20} + \frac{\arctan(100(v - 100))}{10\pi} \right)} - e^{v \left( (Q_s(x,y))^2 + (P_s(x,y))^2 - 1 \right)}}.

U_s(x,y) = 1 - (1 - M_s(x,y))(1 - N_s(x,y)).

M_s(x,y) =
e^{-e^{-100 \left( P_s(x,y) - \frac{57}{100} - \left( \frac{3}{20} + \frac{\cos(7Q_s(x,y) + 2s)}{10} \right) \cos \left( (10 + 3\cos(14s)) \arccos(R_{0,s}(x,y)) + \frac{3}{10} \cos(45x + 47y + \cos(17x)) + 2\cos(5s) \right) \right) }}
\]
\[
- e^{-e^{1000 \left( P_s(x,y) - \frac{18}{25} + \left(\frac{3}{2} Q_s(x,y) \right)^8 \right)}}.

N_s(x,y) =
e^{-e^{100 \left( P_s(x,y) - \frac{37}{50} - \left(\frac{3}{20} + \frac{\cos(8Q_s(x,y) + 5s)}{10} \right) \cos \left( (10 + 3\cos(16s)) \arccos(R_{1,s}(x,y)) + \frac{3}{10} \cos(38x - 47y + \cos(19x)) + 2\cos(4s) \right) \right)} - e^{-1000 \left( P_s(x,y) - \frac{71}{100} - \left(\frac{3}{2} Q_s(x,y) \right)^8 \right)}}.

R_{t,s}(x,y) = E_{t,s}(x,y) e^{-e^{1000(|E_{t,s}(x,y)| - 1)}}

E_{t,s}(x,y) =
\frac{1000}{\sqrt{20}} Q_s(x,y) \sqrt{5|20 - 20(1 - 2t) P_s(x,y) - 27t|}
\left( 1 + 50\sqrt{|4(200 - (20(1 - 2t) P_s(x,y) + 27t)^2)|} \right)^{-1}

P_s(x,y) =
\arctan \left( \tan \left( 2\sin(5s) x - 2\cos(5s) y + 3\cos(5s) + \frac{3\cos(14x - 19y + 5s)}{200} \right) \right).

Q_s(x,y) =
\arctan \left( \tan \left( 2(\cos(5s) x + \sin(5s) y + 2\cos(4s)) + \frac{3\cos(18x + 15y + 4s)}{200} \right) \right).

W(x,y) =
\sum_{s=1}^{40} e^{-e^{-3 \cos^2 \left( 28^s 25^{-s} (\cos(2s)x + \sin(2s)y) + 2 \sin(5s) \right) \cos^2 \left( 28^s 25^{-s} (\cos(2s)y - \sin(2s)x) + 2 \sin(6s) \right) - \frac{97}{100} }}.

"""

from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
from tqdm.auto import tqdm


def safe_exp(x):
    return np.exp(np.minimum(x, 700))


def F(x):
    return (
        255
        * safe_exp(-safe_exp(-1000 * x))
        * np.abs(x) ** (safe_exp(-safe_exp(1000 * (x - 1))))
    )


def L(v, s, x, y):
    R_0 = R(0, s, x, y)
    P_s = P(s, x, y)
    C_20 = C(20, s, x, y)
    C_10 = C(10, s, x, y)
    A_4 = A(4, s, x, y)
    A_1000 = A(1000, s, x, y)
    B_s = B(s, x, y)

    left_term_1 = 1 / 10 - (1 / 40) * np.cos(20 * np.arccos(R_0)) * np.cos(25 * P_s)
    left_term_2 = (
        4 * v**2
        - 13 * v
        + 11
        + np.cos(7 * s + v * s)
        + 20 * safe_exp(-safe_exp(-70 * (C_20 - 1 / 2)))
        + 20 * safe_exp(-safe_exp(-10 * (C_10 - 1 / 2)))
    )

    return left_term_1 * left_term_2 * A_4 * A_1000 + B_s


def C(v, s, x, y):
    R_0 = R(0, s, x, y)
    P_s = P(s, x, y)
    Q_s = Q(s, x, y)
    W_xy = W(x, y)
    cos_term = np.cos(10 * np.arccos(R_0)) * np.cos(25 / 2 * P_s)
    sin_term = np.sin(10 * np.arccos(R_0)) * np.sin(25 / 2 * P_s)

    left_exponent = -(
        safe_exp(v * (cos_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (cos_term + 7 / 10 + W_xy / 5))
        + safe_exp(v * (sin_term - 7 / 10 - W_xy / 5))
        + safe_exp(-v * (sin_term + 7 / 10 + W_xy / 5))
    )
    right_exponent = -safe_exp(
        -3 / 2 * ((Q_s) ** 2 + (P_s - 1 / 4) ** 2 - 21 / 50 + W_xy / 5)
    )

    return safe_exp(left_exponent + right_exponent)


@lru_cache
def B(s, x, y):
    R_0 = R(0, s, x, y)
    P_s = P(s, x, y)
    exponent = -safe_exp(
        -70 * (np.cos(20 * np.arccos(R_0)) * np.cos(25 * P_s) - 47 / 50)
    )
    return safe_exp(exponent)


@lru_cache
def A(v, s, x, y):
    P_s = P(s, x, y)
    Q_s = Q(s, x, y)
    left_exponent = -safe_exp(-1000 * (s - 1 / 2))
    center_exponent = -safe_exp(
        v
        * (
            (5 / 4) * (1 - P_s) * (Q_s) ** 2
            + (P_s) ** 2
            - 11 / 20
            + np.arctan(100 * (v - 100)) / (10 * np.pi)
        )
    )
    right_exponent = -safe_exp(v * (Q_s**2 + P_s**2 - 1))

    return safe_exp(left_exponent + center_exponent + right_exponent)


def U(s, x, y):
    return 1 - (1 - M(s, x, y)) * (1 - N(s, x, y))


def M(s, x, y):
    P_s = P(s, x, y)
    Q_s = Q(s, x, y)
    R_0 = R(0, s, x, y)

    left_exponent = -safe_exp(
        -100
        * (
            P_s
            - 57 / 100
            - (3 / 20 + np.cos(7 * Q_s + 2 * s) / 10)
            * np.cos(
                (10 + 3 * np.cos(14 * s)) * np.arccos(R_0)
                + 3 / 10 * np.cos(45 * x + 47 * y + np.cos(17 * x))
                + 2 * np.cos(5 * s)
            )
        )
    )
    right_exponent = -safe_exp(1000 * (P_s - 18 / 25 + (3 / 2 * Q_s) ** 8))
    return safe_exp(left_exponent + right_exponent)


def N(s, x, y):
    P_s = P(s, x, y)
    Q_s = Q(s, x, y)
    R_1 = R(1, s, x, y)

    left_exponent = -safe_exp(
        -100
        * (
            P_s
            - 37 / 50
            - (3 / 20 + np.cos(8 * Q_s + 5 * s) / 10)
            * np.cos(
                (10 + 3 * np.cos(16 * s)) * np.arccos(R_1)
                + 3 / 10 * np.cos(38 * x - 47 * y + np.cos(19 * x))
                + 2 * np.cos(4 * s)
            )
        )
    )
    right_exponent = -safe_exp(1000 * (P_s - 71 / 100 + (3 / 2 * Q_s) ** 8))
    return safe_exp(left_exponent + right_exponent)


@lru_cache
def R(t, s, x, y):
    E_t_s = E(t, s, x, y)
    return E_t_s * safe_exp(-safe_exp(1000 * (np.abs(E_t_s) - 1)))


def E(t, s, x, y):
    Q_s = Q(s, x, y)
    P_s = P(s, x, y)
    left_term = 1000 / np.sqrt(20) * Q_s
    center_term = np.sqrt(5 * np.abs(20 - 20 * (1 - 2 * t) * P_s - 27 * t))
    right_term = (
        1 + 50 * np.sqrt(np.abs(4 * (200 - (20 * (1 - 2 * t) * P_s + 27 * t) ** 2)))
    ) ** -1
    return left_term * center_term * right_term


@lru_cache
def P(s, x, y):
    return np.arctan(
        np.tan(
            2 * np.sin(5 * s) * x
            - 2 * np.cos(5 * s) * y
            + 3 * np.cos(5 * s)
            + 3 * np.cos(14 * x - 19 * y + 5 * s) / 200
        )
    )


@lru_cache
def Q(s, x, y):
    return np.arctan(
        np.tan(
            2 * (np.cos(5 * s) * x + np.sin(5 * s) * y + 2 * np.cos(4 * s))
            + 3 * np.cos(18 * x + 15 * y + 4 * s) / 200
        )
    )


def W(x, y):

    def W_s(s):

        left = np.cos(
            28**s * 25**-s * (np.cos(2 * s) * x + np.sin(2 * s) * y) + 2 * np.sin(5 * s)
        )
        right = np.cos(
            28**s * 25**-s * (np.cos(2 * s) * y - np.sin(2 * s) * x) + 2 * np.sin(6 * s)
        )
        exp_exponent = -3 * ((left**2) * (right**2) + (-97 / 100))
        exponent = -safe_exp(exp_exponent)
        return safe_exp(exponent)

    return np.sum([W_s(s) for s in range(1, 41)])


def H(v, x, y):
    W_xy = W(x, y)
    sum_term = 0
    for s in range(1, 31):
        prod_term = 1
        for r in range(s):
            prod_term *= (
                (1 - A(1000, r, x, y))
                * (1 - safe_exp(-safe_exp(-1000 * (r - 1 / 2))) * U(r, x, y))
                * (1 - 5 / 4 * A(4, r, x, y))
            )

        sum_term += prod_term * ((5 - 3 * (v - 1) ** 2 + W_xy) / 10) * safe_exp(
            -safe_exp(-abs(71 / 10 - 10 * P(s, x, y)))
        ) * U(s, x, y) + A(1000, s, x, y) * (1 - U(s, x, y)) * L(v, s, x, y)
    return sum_term


def RGB(m, n):
    x = (m - 1000) / 1000
    y = (601 - n) / 600
    return [F(H(0, x, y)), F(H(1, x, y)), F(H(2, x, y))]


def generate_image(width=2000, height=1200):
    img = np.zeros((height, width, 3), dtype=np.uint8)
    f = open("temp2.txt", mode="w")
    for m in tqdm(range(width)):
        temp = ""
        for n in tqdm(range(height), leave=False):
            img[n, m] = RGB(m, n)
            temp += f"{img[n, m]}\t"
        f.write(temp + "\n")
    f.close()
    return img


if __name__ == "__main__":
    image = generate_image()
    plt.imsave("strawberries_generated.png", image)
    plt.imshow(image)
    plt.axis("off")
    plt.show()
