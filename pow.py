#!/usr/bin/python2

"""
Python implementation of

"Randomizing the Montgomery Powering Ladder"
Duc-Phong Le, Chik How Tan and Michael Tunstall
http://eprint.iacr.org/2015/657

(C) 2015 Jack Lloyd

"""

import random
import sys

p = 2**521 - 1
p = 2**127 - 1
#p = 503

test_rounds = 500

def words(k, word_bits):
    orig = k
    v = []
    b = (1 << word_bits)
    while k > 0:
        w = k % b
        v.insert(0, w)
        k /= b

    cmp = 0
    for z in v:
        cmp = (cmp << word_bits) + z
    #print cmp, orig
    assert cmp == orig
    return v

def bits(k):
    return words(k, 1)

def bits0(k):
    kbits = bits(k)
    return (kbits[0:-1], kbits[-1])

def group_op(x, y):
    return (x * y) % p

def random_bit():
    return random.randint(0, 1)

def group_inv(x):
    xinv = pow(x, p - 2, p)
    if group_op(x, xinv) != 1:
        print x, xinv, group_op(x, xinv)
        sys.exit(1)
    return xinv

def algo2(x, k):
    R = [1, x]
    for ki in bits(k):
        R[0] = group_op(R[0], R[ki])
        R[1] = group_op(R[0], x)
    return R[0]

def algo3(x, k):
    R = [1, 1]
    for ki in bits(k):
        R[1] = group_op(R[0], R[ki ^ 1])
        R[0] = group_op(R[1], x)
    return R[0]

def algo4(x, k):
    b = random_bit()
    R = [1, x if b == 0 else 1]

    for ki in bits(k):
        if b ^ ki == 1:
            b = random_bit()
        R[b] = group_op(R[0], R[b ^ ki])
        R[b ^ 1] = group_op(R[b], x)
    return R[0]

def algo5(x, k):
    R = [1,1]
    for ki in bits(k):
        R[ki] = group_op(R[0], R[0])
        R[ki ^ 1] = group_op(R[ki], x)
    return R[0]

def algo6(x, k):
    R = [1,x]
    U = [group_inv(x), x]
    for ki in bits(k):
        R[0] = group_op(R[0], R[1])
        R[1] = group_op(R[0], U[ki])
        l = ki ^ 1
    return R[l]

def precompute_7(x, h):
    U = [0 for i in range(6*h+3)]
    assert len(U) == 6*h+3

    U[3*h+1] = 1
    xinv = group_inv(x)

    for i in range(1, 3*h + 2):
        U[3*h+1+i] = group_op(U[3*h+i], x)
        U[3*h+1-i] = group_op(U[3*h+2-i], xinv)

    #for i in range(len(U) - 1):
    #    assert group_op(U[i], x) == U[i+1]

    return U

def algo7(x, k, h):
    U = precompute_7(x, h)

    kn, k0 = bits0(k)

    #R = [1, 1]
    R = [1, x]
    alpha = 0
    for ki in kn:
        gamma = random.randint(-h, h)
        l = gamma - 2*alpha + ki - (ki ^ 1)
        alpha = gamma
        R[0] = group_op(R[0], R[1])
        R[1] = group_op(R[0], U[3*h + 1 + l])

    R[0] = group_op(R[0], R[1])
    R[0] = group_op(R[0], U[3*h + 1 - alpha - (k0 ^ 1)])
    return R[0]

def algo7a(x, k, h):
    U = precompute_7(x, h)

    kn, k0 = bits0(k)

    #R = [1, 1]
    R = [x, x]
    alpha = 0
    for ki in kn:
        gamma = random.randint(-h, h)
        l = gamma - 2*alpha + ki - (ki ^ 1)
        alpha = gamma
        R[1] = group_op(R[0], U[3*h + 1 + l])
        R[0] = group_op(R[0], R[1])

    R[0] = group_op(R[0], U[3*h + 1 - alpha - (k0 ^ 1)])
    return R[0]


def algo7b(x, k, h):
    U = precompute_7(x, h)

    kn, k0 = bits0(k)

    R = x
    alpha = 0
    for ki in kn:
        gamma = random.randint(-h, h)
        l = gamma - 2*alpha + ki - (ki ^ 1)
        alpha = gamma
        R = group_op(R, R)
        R = group_op(R, U[3*h + 1 + l])

    R = group_op(R, U[3*h + 1 - alpha - (k0 ^ 1)])
    return R

def algo8(x, k, h):
    return x
    U = [0 for i in range(2*h+2)]

    U[h] = 1
    xinv = group_inv(x)

    for i in range(1, h+1):
        U[h+i] = group_op(U[h-1+i], x)
        U[h-i] = group_op(U[h+1-i], xinv)
    U[2*h+1] = group_op(U[2*h], x)

    for i in range(len(U) - 1):
        assert group_op(U[i], x) == U[i+1]

    print h, U

    R = [1, 1]
    alpha = 0
    kbits = bits(k)
    k0 = kbits[len(kbits)-1]
    kbits = kbits[0:len(kbits)-1]
    for ki in kbits:
        R[0] = group_op(R[0], R[1])

        if alpha >= 0:
            lb = 2*alpha - h
            ub = h + 1 - ki
        else:
            lb = -h - ki
            ub = 2*alpha + h - 1
        if ub < lb:
            print "swap"
            (lb,ub) = (ub,lb)
        gamma = random.randint(lb, ub)
        print "rnd", lb, ub, gamma

        R[1] = group_op(R[0], U[h + gamma + ki])
        alpha = 2*alpha - gamma

    R[0] = group_op(R[0], R[1])
    R[0] = group_op(R[0], U[h + alpha + k0])
    return R[0]

def precompute_9(x, h, m):
    T = [0 for i in range(2*(m+1)*h + m - 1)]

    half = len(T)/2
    T[half] = 1
    xinv = group_inv(x)

    for i in range(half + 1, len(T)):
        #print "+1", i, i-1
        T[i] = group_op(T[i-1], x)

    for i in reversed(range(0, half)):
        T[i] = group_op(T[i+1], xinv)

    print "T", T
    for i in range(len(T) - 1):
        assert T[i] > 0
        assert group_op(T[i], x) == T[i+1]

    return T

def algo9(x, k, h, word_bits):
    m = 2**word_bits
    print "bs", word_bits, m

    T = precompute_9(x, h, m)

    R = 1
    alpha = 0

    kw = words(k, word_bits)

    print "k", k
    print "kw", kw

    for w in kw[0:-1]:
        print "w", w
        for i in range(word_bits):
            R = group_op(R, R)
        gamma = random.randrange(-h, h)
        idx = m * alpha - gamma + w
        R = group_op(R, T[idx])
        alpha = gamma

    for i in range(word_bits):
        R = group_op(R, R)
    idx = m * alpha + kw[0]
    R = group_op(R, T[idx])
    return R

def precompute_10(x, h, m):
    T = [0 for i in range((m+1)*h - 1)]

    T[0] = x

    for i in range(1, len(T)):
        T[i] = group_op(T[i-1], x)

    print "T", T
    for i in range(len(T) - 1):
        assert T[i] > 0
        assert group_op(T[i], x) == T[i+1]

    return T

def algo10(x, k, h):
    word_bits = 2
    m = 2**word_bits
    print "bs", word_bits, m

    T = precompute_10(x, h, m)

    R = 1
    alpha = 0

    kw = words(k, word_bits)

    print "k", k
    print "kw", kw

    for w in kw[0:-1]:
        print "w", w
        for i in range(word_bits):
            print "dbl"
            R = group_op(R, R)
        bnd = m * alpha + w

        if min(h, bnd) == 0:
            gamma = 0
        else:
            gamma = random.randrange(0, min(h, bnd))
        idx = bnd - gamma
        R = group_op(R, T[idx])
        alpha = gamma

    for i in range(word_bits):
        R = group_op(R, R)
    idx = m * alpha + kw[0]
    R = group_op(R, T[idx])
    return R

def check(x, k, ref, who, val):
    if val != ref:
        print "Algo", who, "failed:"
        print "x", x
        print "k", k
        print "ref", ref
        print "val", val
        sys.exit(1)

def main():

    for i in range(test_rounds):
        x = random.randint(2, p-1)
        k = random.randint(2, p-1)

        h = random.randint(1, 5)

        bits = random.choice([1,2,4,6,8])

        v0 = pow(x, k, p)

        check(x, k, v0, "2", algo2(x, k))
        check(x, k, v0, "3", algo3(x, k))
        check(x, k, v0, "4", algo4(x, k))
        check(x, k, v0, "5", algo5(x, k))
        check(x, k, v0, "6", algo6(x, k))
        check(x, k, v0, "7", algo7(x, k, h))
        check(x, k, v0, "7a", algo7a(x, k, h))
        check(x, k, v0, "7b", algo7b(x, k, h))

        # Not working:
        #check(x, k, v0, "8", algo8(x, k, h))
        #check(x, k, v0, "9", algo8(x, k, h, bits))
        #check(x, k, v0, "10", algo8(x, k, h, bits))

if __name__ == '__main__':
    main()
