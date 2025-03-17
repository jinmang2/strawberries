# Strawberries

Strawberries by Hamid Naderi Yeganeh

```shell
$ python generate.py --filename strawberries.png
```

2000 X 1200 Image.
For $m = 1, 2, 3, ..., 2000$ and $n = 1, 2, 3, ..., 1200$,
the color of the pixel or the row n and the column m is ::

$$
RGB\Bigg(
F\bigg(
H0\Big(
\cfrac{m - 1000}{600}, \cfrac{601 - n}{600}
\Big)
\bigg),
F\bigg(
H1\Big(
\cfrac{m - 1000}{600}, \frac{601 - n}{600}
\Big)
\bigg),
F\bigg(
H2\Big(
\cfrac{m - 1000}{600}, \cfrac{601 - n}{600}
\Big)
\bigg)
\Bigg)
$$

where

$$F(x) = 255e^{-e^{-1000x}} \cdot |x|^{e^{-e^{1000(x-1)}}}$$

$$
H_v(x, y) =
\sum_{s=1}^{30} \left[ \left(
\prod_{r=0}^{s-1}
\left( 1 - A_{1000,r}(x,y) \right)
\left( 1 - e^{-e^{-1000\left(r - \frac{1}{2} \right)}} U_r(x,y) \right)
\left( 1 - \frac{5}{4} A_{4,r}(x, y) \right)
\right) \right.
$$

$$
\left. \times \left(
\frac{5 - 3(v - 1)^2 + W(x,y)}{10}
 e^{-e^{-\left|\frac{71}{10} - 10P_s(x,y)\right|}} U_s(x,y) +
A_{1000,s}(x,y) \Big(1 - U_s(x,y)\Big) L_{v,s}(x,y)
\right) \right].
$$

$$
L_{v,s}(x,y) =
\left( \frac{1}{10} - \frac{1}{40} \cos \left( 20 \arccos(R_{0,s}(x,y)) \right) \cos(25 P_s(x,y)) \right)
$$

$$
\times \left( 4v^2 - 13v + 11 + \cos(7s + vs) + 20 e^{-e^{-70 (c_{20,s}(x,y) - \frac{1}{2})}} + 20 e^{-e^{-10 (c_{10,s}(x,y) - \frac{1}{2})}} \right)
$$

$$
\times A_{4,s}(x,y) A_{1000,s}(x,y) + B_s(x,y).
$$

$$
C_{v,s}(x,y) =
 e^{-e^{v \left( \cos \left( 10 \arccos(R_{0,s}(x,y)) \right) \cos \left( \frac{25}{2} P_s(x,y) \right) - \frac{7}{10} - \frac{W(x,y)}{5} \right)}
 -e^{v \left( \cos \left( 10 \arccos(R_{0,s}(x,y)) \right) \cos \left( \frac{25}{2} P_s(x,y) \right) + \frac{7}{10} + \frac{W(x,y)}{5} \right)}
 -e^{v \left( \sin \left( 10 \arccos(R_{0,s}(x,y)) \right) \sin \left( \frac{25}{2} P_s(x,y) \right) - \frac{7}{10} - \frac{W(x,y)}{5} \right)}
 -e^{v \left( \sin \left( 10 \arccos(R_{0,s}(x,y)) \right) \sin \left( \frac{25}{2} P_s(x,y) \right) + \frac{7}{10} + \frac{W(x,y)}{5} \right)}}
\times e^{-e^{-\frac{3}{2} \left( (Q_s(x,y))^2 + \left( P_s(x,y) - \frac{1}{4} \right)^2 - \frac{21}{50} + \frac{W(x,y)}{5} \right)}}.
$$

$$

B_s(x,y) = e^{-e^{-70 \left( \cos \left( 20 \arccos(R*{0,s}(x,y)) \right) \cos \left( 25 P_s(x,y) \right) - \frac{47}{50} \right)}}


$$

$$

A_{v,s}(x,y) =
e^{-e^{-1000(s - \frac{1}{2})} - e^{v \left( \frac{5}{4} (1 - P_s(x,y)) (Q_s(x,y))^2 + (P_s(x,y))^2 - \frac{11}{20} + \frac{\arctan(100(v - 100))}{10\pi} \right)} - e^{v \left( (Q_s(x,y))^2 + (P_s(x,y))^2 - 1 \right)}}.


$$

$$

U_s(x,y) = 1 - (1 - M_s(x,y))(1 - N_s(x,y)).


$$

$$

M_s(x,y) =
e^{-e^{-100 \left( P_s(x,y) - \frac{57}{100} - \left( \frac{3}{20} + \frac{\cos(7Q_s(x,y) + 2s)}{10} \right) \cos \left( (10 + 3\cos(14s)) \arccos(R*{0,s}(x,y)) + \frac{3}{10} \cos(45x + 47y + \cos(17x)) + 2\cos(5s) \right) \right) }}


$$

$$

-   e^{-e^{1000 \left( P_s(x,y) - \frac{18}{25} + \left(\frac{3}{2} Q_s(x,y) \right)^8 \right)}}.
$$

$$
N_s(x,y) =
 e^{-e^{100 \left( P_s(x,y) - \frac{37}{50} - \left(\frac{3}{20} + \frac{\cos(8Q_s(x,y) + 5s)}{10} \right) \cos \left( (10 + 3\cos(16s)) \arccos(R_{1,s}(x,y)) + \frac{3}{10} \cos(38x - 47y + \cos(19x)) + 2\cos(4s) \right) \right)}}.
$$

$$

R_{t,s}(x,y) = E_{t,s}(x,y) e^{-e^{1000(|E_{t,s}(x,y)| - 1)}}


$$

$$

E_{t,s}(x,y) =
\frac{1000}{\sqrt{20}} Q_s(x,y) \sqrt{5|20 - 20(1 - 2t) P_s(x,y) - 27t|}
\left( 1 + 50\sqrt{|4(200 - (20(1 - 2t) P_s(x,y) + 27t)^2)|} \right)^{-1}


$$

$$

P_s(x,y) =
\arctan \left( \tan \left( 2\sin(5s) x - 2\cos(5s) y + 3\cos(5s) + \frac{3\cos(14x - 19y + 5s)}{200} \right) \right).


$$

$$

Q_s(x,y) =
\arctan \left( \tan \left( 2(\cos(5s) x + \sin(5s) y + 2\cos(4s)) + \frac{3\cos(18x + 15y + 4s)}{200} \right) \right).


$$

$$

W(x,y) =
\sum_{s=1}^{40} e^{-e^{-3 \cos^2 \left( 28^s 25^{-s} (\cos(2s)x + \sin(2s)y) + 2 \sin(5s) \right) \cos^2 \left( 28^s 25^{-s} (\cos(2s)y - \sin(2s)x) + 2 \sin(6s) \right) - \frac{97}{100} }}.


$$
