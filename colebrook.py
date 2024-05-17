# By John Nguyen <knightofthealtar64@gmail.com>
# Python 3.10.11
# MIT License

from math import sqrt, log

# Program for calculating Darcy-Weisbach friction factor using Colebrook-White equation
# in other words, friction inside a pipe


k = 3.28e-6 * 12 # inch, from https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html
D = 0.5 # inch, from fluids system design

Re = 723250 # dimensionless, from ideal gas calculation

max_err = 0.000000000001

def error(f):
    if f > 0:
        return 1/sqrt(f) - (-2 * log(k / D / 3.7 + 2.51 / Re / sqrt(f)))
    else:
        return 10000


guess = 1
a = 1
b = 0
err = error(guess)

if error(a) > 0 or error(b) < 0:
    print("variables disordered")
    exit(1)

i = 0

while i < 100 :
    mid = (a + b) / 2
    fmid = error(mid)
    err_lim = (-b + a) / 2

    if fmid > 0:
        b = mid
    elif fmid < 0:
        a = mid

    if err_lim < max_err:
        print(f"Converges at {mid} after {i} tries")
        break

    i += 1
