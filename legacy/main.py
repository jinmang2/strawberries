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
e^{-e^{100 \left( P_s(x,y) - \frac{37}{50} - \left(\frac{3}{20} + \frac{\cos(8Q_s(x,y) + 5s)}{10} \right) \cos \left( (10 + 3\cos(16s)) \arccos(R_{1,s}(x,y)) + \frac{3}{10} \cos(38x - 47y + \cos(19x)) + 2\cos(45) \right) \right)} - e^{-1000 \left( P_s(x,y) - \frac{71}{100} - \left(\frac{3}{2} Q_s(x,y) \right)^8 \right)}}.

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

import matplotlib.pyplot as plt
import numpy as np


def safe_exp(x):
    return safe_exp(np.minimum(x, 700))


def F(x):
    return (
        255
        * safe_exp(-safe_exp(-1000 * x))
        * np.abs(x) ** (safe_exp(-safe_exp(1000 * (x - 1))))
    )


def L(v, s, x, y):
    return (
        1 / 10
        - (1 / 40) * np.cos(20 * np.arccos(R(0, s, x, y))) * np.cos(25 * P(s, x, y))
    ) * (
        4 * v**2
        - 13 * v
        + 11
        + np.cos(7 * s + v * s)
        + 20 * safe_exp(-safe_exp(-70 * (C(20, s, x, y) - 1 / 2)))
        + 20 * safe_exp(-safe_exp(-10 * (C(10, s, x, y) - 1 / 2)))
    ) * A(
        4, s, x, y
    ) * A(
        1000, s, x, y
    ) + B(
        s, x, y
    )


def C(v, s, x, y):
    return (
        safe_exp(
            -safe_exp(
                v
                * (
                    np.cos(10 * np.arccos(R(0, s, x, y))) * np.cos(25 / 2 * P(s, x, y))
                    - 7 / 10
                    - W(x, y) / 5
                )
            )
        )
        - safe_exp(
            -safe_exp(
                v
                * (
                    np.cos(10 * np.arccos(R(0, s, x, y))) * np.cos(25 / 2 * P(s, x, y))
                    + 7 / 10
                    + W(x, y) / 5
                )
            )
        )
        - safe_exp(
            -safe_exp(
                v
                * (
                    np.sin(10 * np.arccos(R(0, s, x, y))) * np.sin(25 / 2 * P(s, x, y))
                    - 7 / 10
                    - W(x, y) / 5
                )
            )
        )
        - safe_exp(
            -safe_exp(
                v
                * (
                    np.sin(10 * np.arccos(R(0, s, x, y))) * np.sin(25 / 2 * P(s, x, y))
                    + 7 / 10
                    + W(x, y) / 5
                )
            )
        )
        * safe_exp(
            -safe_exp(
                -3
                / 2
                * (
                    (Q(s, x, y)) ** 2
                    + (P(s, x, y) - 1 / 4) ** 2
                    - 21 / 50
                    + W(x, y) / 5
                )
            )
        )
    )


def B(s, x, y):
    return safe_exp(
        -safe_exp(
            -70
            * (
                np.cos(20 * np.arccos(R(0, s, x, y))) * np.cos(25 * P(s, x, y))
                - 47 / 50
            )
        )
    )


def A(v, s, x, y):
    result = (
        safe_exp(-safe_exp(-1000 * (s - 1 / 2)))
        - safe_exp(
            v
            * (
                (5 / 4) * (1 - P(s, x, y)) * (Q(s, x, y)) ** 2
                + (P(s, x, y)) ** 2
                - 11 / 20
                + np.arctan(100 * (v - 100)) / (10 * np.pi)
            )
        )
        - safe_exp(v * ((Q(s, x, y)) ** 2 + (P(s, x, y)) ** 2 - 1))
    )

    if np.isnan(result) or np.isinf(result):
        print(f"A({v}, {s}, {x}, {y}) resulted in NaN or Inf")

    return result


def U(s, x, y):
    return 1 - (1 - M(s, x, y)) * (1 - N(s, x, y))


def M(s, x, y):
    return safe_exp(
        -safe_exp(
            -100
            * (
                P(s, x, y)
                - 57 / 100
                - (3 / 20 + np.cos(7 * Q(s, x, y) + 2 * s) / 10)
                * np.cos(
                    (10 + 3 * np.cos(14 * s)) * np.arccos(R(0, s, x, y))
                    + 3 / 10 * np.cos(45 * x + 47 * y + np.cos(17 * x))
                    + 2 * np.cos(5 * s)
                )
            )
        )
    )


def N(s, x, y):
    return safe_exp(
        -safe_exp(
            100
            * (
                P(s, x, y)
                - 37 / 50
                - (3 / 20 + np.cos(8 * Q(s, x, y) + 5 * s) / 10)
                * np.cos(
                    (10 + 3 * np.cos(16 * s)) * np.arccos(R(1, s, x, y))
                    + 3 / 10 * np.cos(38 * x - 47 * y + np.cos(19 * x))
                    + 2 * np.cos(45)
                )
            )
        )
    )


def R(t, s, x, y):
    return E(t, s, x, y) * safe_exp(-safe_exp(1000 * (np.abs(E(t, s, x, y)) - 1)))


def E(t, s, x, y):
    return (
        (1000 / np.sqrt(20))
        * Q(s, x, y)
        * np.sqrt(5 * np.abs(20 - 20 * (1 - 2 * t) * P(s, x, y) - 27 * t))
        * (
            1
            + 50
            * np.sqrt(np.abs(4 * (200 - (20 * (1 - 2 * t) * P(s, x, y) + 27 * t) ** 2)))
        )
        ** -1
    )


def P(s, x, y):
    result = np.arctan(
        np.tan(
            2 * np.sin(5 * s) * x
            - 2 * np.cos(5 * s) * y
            + 3 * np.cos(5 * s)
            + 3 * np.cos(14 * x - 19 * y + 5 * s) / 200
        )
    )

    if np.isnan(result) or np.isinf(result):
        print(f"P({s}, {x}, {y}) resulted in NaN or Inf")

    return result


def Q(s, x, y):
    result = np.arctan(
        np.tan(
            2 * (np.cos(5 * s) * x + np.sin(5 * s) * y + 2 * np.cos(4 * s))
            + 3 * np.cos(18 * x + 15 * y + 4 * s) / 200
        )
    )

    if np.isnan(result) or np.isinf(result):
        print(f"Q({s}, {x}, {y}) resulted in NaN or Inf")

    return result


def W(x, y):
    return np.sum(
        [
            safe_exp(
                -safe_exp(
                    -3
                    * np.cos(
                        28**s * 25**-s * (np.cos(2 * s) * x + np.sin(2 * s) * y)
                        + 2 * np.sin(5 * s)
                    )
                    ** 2
                    * np.cos(
                        28**s * 25**-s * (np.cos(2 * s) * y - np.sin(2 * s) * x)
                        + 2 * np.sin(6 * s)
                    )
                    ** 2
                    - 97 / 100
                )
            )
            for s in range(1, 41)
        ]
    )


def H(v, x, y):
    W_val = W(x, y)
    sum_term = 0
    for s in range(1, 31):
        prod_term = 1
        for r in range(s):
            prod_term *= (
                (1 - A(1000, r, x, y))
                * (1 - safe_exp(-safe_exp(-1000 * (r - 1 / 2))) * U(r, x, y))
                * (1 - 5 / 4 * A(4, r, x, y))
            )

        sum_term += prod_term * ((5 - 3 * (v - 1) ** 2 + W_val) / 10) * safe_exp(
            -safe_exp(-abs(71 / 10 - 10 * P(s, x, y)))
        ) * U(s, x, y) + A(1000, s, x, y) * (1 - U(s, x, y)) * L(v, s, x, y)
    return sum_term


def RGB(m, n):
    x = (m - 1000) / 600
    y = (601 - n) / 600
    return [F(H(0, x, y)), F(H(1, x, y)), F(H(2, x, y))]


def generate_image(width=2000, height=1200):
    img = np.zeros((height, width, 3), dtype=np.uint8)
    for m in range(width):
        for n in range(height):
            img[n, m] = RGB(m, n)
    return img


if __name__ == "__main__":
    # image = generate_image()
    # plt.imshow(image)
    # plt.axis("off")
    # plt.show()

    print(RGB(0, 1))
