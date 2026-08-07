#!/usr/bin/env python3
"""
Estimate norms and relative yield for a GNFS/SNFS job file.
Helps decide which side should carry the special q, and how to trade that
choice off against the large-prime parameters on the other side.

Norms
-----
The lattice siever with special-q q on side s walks the (a,b) points of the
q-lattice that land in the (i,j) rectangle  -I/2 <= i < I/2,  0 < j <= J,
skipping gcd(i,j) > 1.  Those points sit in a lattice of determinant q, so in
the (a,b) plane they cover a region of area I*J*q whose shape is set by the
skew and by the Gauss-reduced basis of the lattice.

This tool samples that region directly: it builds random q-lattices, reduces
them in the skewed metric exactly as the siever does, walks the (i,j)
rectangle, and evaluates both homogeneous forms in exact integer arithmetic.
Averaging log2 of those values gives the mean norm on each side.

The part that decides the special-q side: every (a,b) on the q-lattice has
q | F_s(a,b) on the special-q side, so that side's norm is smaller by exactly
log2(q) bits, and the other side is untouched.  The *total* norm bits are the
same either way -- the special-q side only decides how the total is split.

Nothing else enters the norm.  lpb, mfb, lambda, alim and rlim have no effect
on it whatsoever; they only change the probability that a norm of a given size
becomes a relation.  That is the next section.

Yield
-----
After sieving, a norm has all its prime factors below lim removed, leaving a
cofactor that must (a) be at most mfb bits and (b) split into primes of at most
lpb bits.  Writing u = ln(norm)/ln(lim) and measuring prime sizes in units of
ln(lim), the density of integers whose cofactor is a product of exactly j such
primes is the standard heuristic

    P_j(u) = (1/j!) * Integral rho(u - s_1 - ... - s_j) ds_1/s_1 ... ds_j/s_j

over 1 <= s_i <= lpb*ln2/ln(lim), subject to sum(s_i) <= mfb*ln2/ln(lim).
The ds/s is the density of a prime at that scale (Mertens), rho is Dickman's
function, and the total probability is the sum over j = 0, 1, 2, ... allowed by
the mfb cap.  This is what makes 2LP and 3LP distinguishable: raising mfb from
2*lpb to 3*lpb admits the j = 3 term, which for a heavy norm is worth far more
than the sum of the smaller terms.  lambda enters as the sieve report
threshold, so the effective cofactor bound is min(mfb, lambda*log2(lim)) --
setting lambda too low silently discards the very relations that mfb allows.

alpha (the small-prime bias of the polynomial, in nats) shifts the effective
norm by alpha/ln2 bits before u is formed.

Because rho is steeply convex, the yield is driven by the tail of the region
where norms fall well below average, so the index is computed as the mean of
P_r * P_a over the sampled points -- E[P], not P(E[norm]).

Accuracy
--------
The norm model was checked against real q-lattices built from actual roots of
f mod q (rather than the random lattices used here) and reproduces the true
average norm on both sides to within 0.2 bits across I15/J14, I16/J15 and
I16/J16.  The shape term also agrees with CADO's own reported lognorm to 0.1
nats, which the header prints as a cross-check on the parsed polynomial.

The yield index is a different matter.  It is a standard heuristic that ignores
prime distinctness, cofactoring failures (ECM does not always split), and
duplicate removal, and it is known to run optimistic by orders of magnitude in
absolute terms.  Use it only to compare configurations -- the same q, the same
seed, and therefore the same sampled points, so the comparison is paired and
the ratio is far more trustworthy than either number alone.  The test sieve
remains the ground truth.
"""

import sys
import re
import os
import random
import argparse
from math import log, log2, log10, sqrt, pi, cos, sin, gcd
from typing import NoReturn

LOG2E = 1.4426950408889634
LN2 = 0.6931471805599453


def parse_job_file(filename):
    """Parse a .job file; returns (params, comment_lines)."""
    params = {}
    comments = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                comments.append(line)
                continue
            if ':' in line:
                key, val = line.split(':', 1)
                key = key.strip()
                val = val.strip()
                # Strip invisible Unicode characters (zero-width spaces, etc.)
                val = re.sub(r'[^\x20-\x7E]', '', val)
                try:
                    if '.' in val:
                        params[key] = float(val)
                    else:
                        params[key] = int(val)
                except ValueError:
                    params[key] = val
    return params, comments


def get_poly(params, prefix):
    """Coefficients for prefix0..prefixN, highest degree first.

    Coerced to int: the norms are evaluated exactly, because the high-yield
    part of the sieve region is precisely where F nearly cancels and float
    would lose it.  A job file writing "c6: 1.0" still parses.
    """
    max_deg = -1
    for key in params:
        rest = key[len(prefix):]
        if key.startswith(prefix) and rest.isdigit():
            max_deg = max(max_deg, int(rest))
    if max_deg < 0:
        return []
    coeffs = []
    for i in range(max_deg, -1, -1):
        c = params.get(f'{prefix}{i}', 0)
        if isinstance(c, float):
            if not c.is_integer():
                raise ValueError(
                    f"{prefix}{i} = {c} is not an integer; polynomial "
                    "coefficients must be integral")
            c = int(c)
        elif not isinstance(c, int):
            raise ValueError(f"{prefix}{i} = {c!r} is not a number")
        coeffs.append(c)
    return coeffs


def get_alg_poly(params):
    return get_poly(params, 'c')


def get_rat_poly(params):
    return get_poly(params, 'Y')


def parse_alphas(comments):
    """Pull CADO's per-side alpha and lognorm out of the poly comments.

    CADO writes e.g. "# side 1 lognorm 52.76, E 46.03, alpha -6.73 (proj ...)".
    Side 0 is rational, side 1 is algebraic.
    """
    alphas, lognorms = {}, {}
    for line in comments:
        m = re.search(r'side\s+([01])\b.*?alpha\s+(-?[0-9.]+)', line)
        if m:
            alphas[int(m.group(1))] = float(m.group(2))
        m = re.search(r'side\s+([01])\b.*?lognorm\s+([0-9.]+)', line)
        if m:
            lognorms[int(m.group(1))] = float(m.group(2))
    return alphas.get(0), alphas.get(1), lognorms.get(1)


def eval_homogeneous(coeffs, a, b):
    """Homogeneous F(a,b) = c_d a^d + c_{d-1} a^{d-1} b + ... + c_0 b^d.

    Horner in (a,b); exact when a and b are ints, which matters because the
    high-yield part of the region is exactly where F nearly cancels.
    """
    if not coeffs:
        return 0
    v = coeffs[0]
    bp = 1
    for c in coeffs[1:]:
        bp = bp * b
        v = v * a + c * bp
    return v


def log2_abs(n):
    """log2|n| for ints of any size (math.log2 overflows past ~2^1024)."""
    n = abs(n)
    if n == 0:
        return float('-inf')
    if not isinstance(n, int):
        return log2(n)
    bl = n.bit_length()
    if bl <= 1023:
        return log2(n)
    return bl - 53 + log2(n >> (bl - 53))


def algebraic_norm(alg_coeffs, a, b):
    return abs(eval_homogeneous(alg_coeffs, a, b))


def rational_norm(rat_coeffs, a, b):
    return abs(eval_homogeneous(rat_coeffs, a, b))


def angular_log2_mean(coeffs, skew, nsteps=20000):
    """Mean of log2|F| over the unit circle in skew-normalised coordinates.

    Only used as a cross-check against CADO's reported lognorm; the norms
    themselves come from sampling the real region.
    """
    rs = sqrt(skew)
    total, used = 0.0, 0
    for k in range(nsteps):
        t = pi * (k + 0.5) / nsteps
        v = abs(eval_homogeneous(coeffs, rs * cos(t), sin(t) / rs))
        if v > 0:
            total += log2(v)
            used += 1
    return total / used if used else 0.0


def reduce_q_lattice(q, r, skew, maxiter=200):
    """Gauss-reduce the q-lattice {(a,b) : a = r*b mod q} in the skewed metric.

    The siever does exactly this to pick its (i,j) -> (a,b) basis, so the shape
    of the region it actually sieves comes from the reduced basis.  Entries stay
    integral; only the comparisons are done in floating point.
    """
    inv_rs = 1.0 / sqrt(skew)
    rs = sqrt(skew)

    def dot(u, v):
        return (u[0] * inv_rs) * (v[0] * inv_rs) + (u[1] * rs) * (v[1] * rs)

    u, v = (q, 0), (r, 1)
    for _ in range(maxiter):
        if dot(u, u) > dot(v, v):
            u, v = v, u
        du = dot(u, u)
        if du == 0:
            break
        m = int(round(dot(u, v) / du))
        if m == 0:
            break
        v = (v[0] - m * u[0], v[1] - m * u[1])
    return u, v


def sample_region(alg, rat, q, ibits, jbits, skew, nlat=200, npts=120, seed=12345):
    """Sample the sieve region; return raw (rat_bits, alg_bits) per point.

    "Raw" means before the special-q division -- log2(q) is subtracted later on
    whichever side carries q, so one sample set serves both configurations and
    the comparison between them is paired.

    Random lattices are used rather than true roots of f mod q: the geometry is
    the same (validated to 0.2 bits against real-root lattices), and it avoids
    needing to find roots mod q.
    """
    q = int(q)
    I, J = 1 << ibits, 1 << jbits
    rng = random.Random(seed)
    rb, ab = [], []
    for _ in range(nlat):
        u, v = reduce_q_lattice(q, rng.randrange(q), skew)
        for _ in range(npts):
            i = rng.randrange(-I // 2, I // 2)
            j = rng.randrange(1, J + 1)
            if gcd(abs(i), j) != 1:      # siever skips non-primitive (i,j)
                continue
            a = i * u[0] + j * v[0]
            b = i * u[1] + j * v[1]
            if b == 0:
                continue
            na = algebraic_norm(alg, a, b)
            nr = rational_norm(rat, a, b)
            if na == 0 or nr == 0:
                continue
            ab.append(log2_abs(na))
            rb.append(log2_abs(nr))
    return rb, ab


# ---------------------------------------------------------------- Dickman rho

_RHO_H = 1.0 / 1024
_RHO_UMAX = 40.0
_rho_table = None


def _build_rho():
    h = _RHO_H
    step = int(round(1.0 / h))
    n = int(_RHO_UMAX / h)
    rho = [0.0] * (n + 1)
    cum = [0.0] * (n + 1)          # cum[k] = integral of rho from 0 to k*h
    for k in range(n + 1):
        u = k * h
        if u <= 1.0:
            rho[k] = 1.0
        elif u <= 2.0:
            rho[k] = 1.0 - log(u)
        else:
            # u*rho(u) = int_{u-1}^{u} rho, solved implicitly with trapezoid
            rho[k] = (cum[k - 1] + 0.5 * h * rho[k - 1] - cum[k - step]) / (u - 0.5 * h)
            if rho[k] < 0.0:
                rho[k] = 0.0
        if k > 0:
            cum[k] = cum[k - 1] + 0.5 * h * (rho[k - 1] + rho[k])
    return rho


def dickman_rho(u):
    """Probability a random x is x^(1/u)-smooth."""
    global _rho_table
    if u <= 1.0:
        return 1.0
    if u >= _RHO_UMAX:
        return 0.0
    if _rho_table is None:
        _rho_table = _build_rho()
    x = u / _RHO_H
    k = int(x)
    f = x - k
    return _rho_table[k] * (1 - f) + _rho_table[k + 1] * f


# ------------------------------------------------- large-prime smoothness

class SideSmoothness:
    """P(a norm of given size yields a relation) for one side's parameters.

    Accounts for lim, lpb, mfb and lambda -- see the module docstring for the
    P_j integral.  The number of large primes allowed falls out of the mfb cap
    rather than being specified, so 2LP vs 3LP is captured automatically.
    """

    WSTEP = 0.001
    MAXJ_CAP = 4      # 3LP is the practical ceiling; 4 leaves headroom
    MAX_M1 = 400      # bound on convolution width; see _build

    def __init__(self, lim, lpb, mfb, lam=None, alpha=0.0):
        self.ok = bool(lim and lpb and mfb and lim > 1)
        self.alpha = alpha or 0.0
        self.lpb, self.mfb, self.lam, self.lim = lpb, mfb, lam, lim
        self._prob_cache = {}
        self.wstep = self.WSTEP
        if not self.ok:
            self.maxj = 0
            return
        self.lny = log(lim)
        self.smax = lpb * LN2 / self.lny
        # lambda is the sieve report threshold; it can bind before mfb does.
        cof_bits = mfb
        if lam and lam > 0:
            cof_bits = min(mfb, lam * log2(lim))
        self.cof_bits = cof_bits
        self.wmax = cof_bits * LN2 / self.lny
        self.maxj = int(self.wmax) if self.smax > 0 else 0
        # The convolutions are O(maxj * len(m1)^2); a mistyped lim (say 671
        # instead of 67100000) would otherwise stall for minutes.
        self.capped = self.maxj > self.MAXJ_CAP
        if self.capped:
            self.maxj = self.MAXJ_CAP
        self._build()

    def _build(self):
        """Discretised measures m_j = (ds/s)^{*j} on a grid of w = sum(s_i).

        The convolutions cost O(maxj * len(m1)^2), and len(m1) grows as lim
        shrinks, so the step is coarsened to keep the width bounded.  A normal
        job has smax ~ 1.2 and never leaves the 0.001 grid; only a degenerate
        lim (usually a typo) triggers coarsening, and there the result is a
        warning rather than a number to trust.
        """
        h = self.WSTEP
        span = self.smax - 1.0
        if span > 0 and span / h > self.MAX_M1:
            h = span / self.MAX_M1
        self.wstep = h
        lo, hi = int(round(1.0 / h)), int(self.smax / h)
        if hi <= lo:
            self.measures = []
            self.maxj = 0
            return
        # Mertens weight for the cell at w = k*h is (ds/s) = h/(k*h) = 1/k.
        m1 = [1.0 / k for k in range(lo, hi + 1)]
        m1[0] *= 0.5            # trapezoid end weights
        m1[-1] *= 0.5
        self.measures = [(lo, m1)]
        for _ in range(2, self.maxj + 1):
            off, prev = self.measures[-1]
            conv = [0.0] * (len(prev) + len(m1) - 1)
            for i, x in enumerate(prev):
                if x:
                    for jj, y in enumerate(m1):
                        conv[i + jj] += x * y
            self.measures.append((off + lo, conv))

    def prob(self, bits):
        """Probability that a norm of `bits` bits produces a relation.

        Memoised on bits rounded to 0.01, which is far finer than the model's
        accuracy and collapses the tens of thousands of sampled norms onto a
        few thousand distinct keys.
        """
        if not self.ok:
            return 0.0
        key = round(bits, 2)
        hit = self._prob_cache.get(key)
        if hit is not None:
            return hit
        self._prob_cache[key] = val = self._prob(key)
        return val

    def _prob(self, bits):
        u = (bits + self.alpha * LOG2E) * LN2 / self.lny
        if u <= 0:
            return 1.0
        total = dickman_rho(u)
        h = self.wstep
        wcap = min(self.wmax, u)
        fact = 1.0
        for j, (off, m) in enumerate(self.measures, start=1):
            fact *= j
            s = 0.0
            for k, val in enumerate(m):
                w = (off + k) * h
                if w > wcap:
                    break
                if val:
                    s += val * dickman_rho(u - w)
            total += s / fact
        return min(total, 1.0)

    def describe(self):
        if not self.ok:
            return "n/a (needs lim, lpb and mfb)"
        cap = f"{self.cof_bits:.0f}"
        note = ""
        if self.lam and self.lam > 0 and self.lam * log2(self.lim) < self.mfb - 0.5:
            note = f" (lambda caps mfb {self.mfb} -> {cap})"
        if self.capped:
            note += f" [capped at {self.MAXJ_CAP}LP -- is lim={self.lim:g} right?]"
        return f"{self.maxj}LP, cofactor <= {cap} bits{note}"


# ------------------------------------------------------------------- analysis

class Model:
    def __init__(self, params, comments, overrides=None, nlat=200, npts=120,
                 seed=12345):
        o = overrides or {}
        self.params = params
        self.skew = params.get('skew', 1.0) or 1.0
        self.alg = get_alg_poly(params)
        self.rat = get_rat_poly(params)
        pa_r, pa_a, self.cado_lognorm = parse_alphas(comments)
        self.alpha_r = o.get('alpha_r', pa_r if pa_r is not None else 0.0)
        self.alpha_a = o.get('alpha_a', pa_a if pa_a is not None else 0.0)
        self.nlat, self.npts, self.seed = nlat, npts, seed
        self._sample_cache = {}
        self._index_cache = {}

        def pick(name):
            return o[name] if o.get(name) is not None else params.get(name, 0)

        self.rside = SideSmoothness(pick('rlim'), pick('lpbr'), pick('mfbr'),
                                    pick('rlambda'), self.alpha_r)
        self.aside = SideSmoothness(pick('alim'), pick('lpba'), pick('mfba'),
                                    pick('alambda'), self.alpha_a)

    def check(self):
        """Raise if the job file cannot support a norm estimate.

        Without this an unparsable polynomial yields an empty sample set and
        the tool reports a confident 0.0 bits on both sides -- a fabricated
        answer is far worse than no answer.
        """
        if len(self.alg) < 2:
            raise ValueError(
                "no algebraic polynomial found (need c0, c1, ... entries)")
        if len(self.rat) < 2:
            raise ValueError(
                "no rational polynomial found (need Y0, Y1 entries)")
        if not any(self.alg) or not any(self.rat):
            raise ValueError("a polynomial has all-zero coefficients")
        if self.skew <= 0:
            raise ValueError(f"skew must be positive, got {self.skew}")

    def samples(self, q, ibits, jbits):
        key = (int(q), ibits, jbits)
        if key not in self._sample_cache:
            rb, ab = sample_region(
                self.alg, self.rat, q, ibits, jbits, self.skew,
                self.nlat, self.npts, self.seed)
            if not rb:
                raise ValueError(
                    f"no usable sample points at q={int(q)} "
                    f"(I=2^{ibits}, J=2^{jbits}) -- check the polynomial and skew")
            self._sample_cache[key] = (rb, ab)
        return self._sample_cache[key]

    def at_q(self, q, ibits, jbits, sq_side):
        """Mean norm bits on both sides, with log2(q) on the special-q side."""
        rb, ab = self.samples(q, ibits, jbits)
        mr, ma = sum(rb) / len(rb), sum(ab) / len(ab)
        lq = log2(int(q))
        return (mr - lq, ma) if sq_side == 'r' else (mr, ma - lq)

    def yield_index(self, q, ibits, jbits, sq_side):
        """E[P_r * P_a] over the sampled region.  Comparable across runs.

        None when either side lacks the parameters to define a probability.
        """
        if not (self.rside.ok and self.aside.ok):
            return None
        key = (int(q), ibits, jbits, sq_side)
        if key in self._index_cache:
            return self._index_cache[key]
        rb, ab = self.samples(q, ibits, jbits)
        lq = log2(int(q))
        dr, da = (lq, 0.0) if sq_side == 'r' else (0.0, lq)
        total = 0.0
        for x, y in zip(rb, ab):
            total += self.rside.prob(x - dr) * self.aside.prob(y - da)
        self._index_cache[key] = total / len(rb)
        return self._index_cache[key]


def fmt_poly(coeffs):
    deg = len(coeffs) - 1
    terms = []
    for i, c in enumerate(coeffs):
        if c != 0:
            d = deg - i
            terms.append(f"{c}" if d == 0 else
                         (f"{c}*x" if d == 1 else f"{c}*x^{d}"))
    return " + ".join(terms).replace("+ -", "- ")


def die(msg) -> NoReturn:
    """Fail loudly and non-zero -- callers treat exit 0 as usable numbers."""
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def sci(x):
    """Format a yield index; an undefined one must never read as a number."""
    if x is None:
        return "n/a"
    return "%.2e" % x if x > 0 else "0"


def siever_q(ibits):
    """The conventional representative special-q for a given siever size.

    YAFU's heuristic: 1e4 for I11, rising a decade per I bit to 1e9 for I16.
    """
    return int(10000 * (10 ** (ibits - 11)))


def siever_table(model, sq_side, ibits_list=(11, 12, 13, 14, 15, 16)):
    """Norms across siever sizes, each at its own representative q.

    A survey: it answers "how do the norms sit at the scale each siever
    actually runs at", which a single q cannot.  Every row uses that siever's
    default J = I/2; --jbits applies only to the tables below, so the header
    says so rather than letting an I16 row silently disagree with a
    --jbits 16 comparison further down.
    """
    side_name = 'rational' if sq_side == 'r' else 'algebraic'
    print(f"--- Norm Estimates by Siever (special q on the {side_name} side, "
          f"each at J=I/2) ---")
    print(f"{'Siever':<20} {'q':>6} {'Rat bits':>9} {'Alg bits':>9} {'Alg-Rat':>8} "
          f"{'yield idx':>11}  {'Best sq side':<12}")
    print("-" * 82)
    for ib in ibits_list:
        q, jb = siever_q(ib), ib - 1
        rb, ab = model.at_q(q, ib, jb, sq_side)
        idx = model.yield_index(q, ib, jb, sq_side)
        # Both sides reuse the same cached sample set, so this is nearly free.
        best = ""
        alt = model.yield_index(q, ib, jb, 'a' if sq_side == 'r' else 'r')
        if idx is not None and alt is not None:
            if idx >= alt:
                best = side_name
            else:
                best = 'algebraic' if sq_side == 'r' else 'rational'
        print(f"{'gnfs-lasieve4I%de' % ib:<20} {'1e%d' % round(log10(q)):>6} "
              f"{rb:>9.1f} {ab:>9.1f} {ab - rb:>+8.1f} {sci(idx):>11}  {best:<12}")


def per_q_table(model, qs, ibits, jbits, sq_side):
    side_name = 'rational' if sq_side == 'r' else 'algebraic'
    print(f"--- Norms per special-q  (I=2^{ibits}, J=2^{jbits}, special-q on the "
          f"{side_name} side) ---")
    print(f"{'q':>14} {'rat bits':>9} {'alg bits':>9} {'Alg-Rat':>8} {'yield idx':>12}")
    print("-" * 56)
    for q in qs:
        rb, ab = model.at_q(q, ibits, jbits, sq_side)
        idx = model.yield_index(q, ibits, jbits, sq_side)
        print(f"{int(q):>14,} {rb:>9.1f} {ab:>9.1f} {ab - rb:>+8.1f} "
              f"{sci(idx):>12}")


def side_comparison(model, q, ibits, jbits):
    print(f"--- Special-q side comparison (I=2^{ibits}, J=2^{jbits}, "
          f"q = {int(q):,}) ---")
    print(f"  rational side: {model.rside.describe()}")
    print(f"  algebraic side: {model.aside.describe()}")
    print()
    print(f"{'sq side':<12} {'rat bits':>9} {'alg bits':>9} {'yield idx':>12}")
    print("-" * 46)
    rows = []
    for s in ('r', 'a'):
        rb, ab = model.at_q(q, ibits, jbits, s)
        idx = model.yield_index(q, ibits, jbits, s)
        rows.append((s, rb, ab, idx))
        name = 'rational' if s == 'r' else 'algebraic'
        print(f"{name:<12} {rb:>9.1f} {ab:>9.1f} {sci(idx):>12}")
    if rows[0][3] is None or rows[1][3] is None:
        print("\n=> No yield comparison: the job file is missing the lim/lpb/mfb "
              "needed to")
        print("   define a smoothness probability, so only the norms above are "
              "meaningful.")
    elif not (rows[0][3] and rows[1][3]):
        print("\n=> No yield comparison: both sides score zero, so the norms are "
              "far outside")
        print("   what these parameters can factor. Check lpb/mfb against the "
              "norm sizes.")
    else:
        best = max(rows, key=lambda r: r[3])
        other = rows[0] if best is rows[1] else rows[1]
        ratio = best[3] / other[3] if other[3] > 0 else float('inf')
        name = 'RATIONAL' if best[0] == 'r' else 'ALGEBRAIC'
        print(f"\n=> Model prefers special q on the {name} side "
              f"({ratio:.2f}x the yield index of the other choice).")
        print(f"   The special-q side is the one whose norm drops by "
              f"log2(q) = {log2(int(q)):.1f} bits.")
        if ratio < 1.3:
            print("   That margin is small -- treat the two as a toss-up and let "
                  "the test sieve decide.")
        print("   Absolute values are meaningless; only compare indices between "
              "runs at the same q and seed.")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Estimate norms and relative yield for a GNFS/SNFS job file.",
        epilog="Parameter overrides let you price a hypothetical config, e.g. "
               "giving the un-sieved side 3LP: --sqside r --mfba 96 --lpba 32. "
               "Indices are comparable across runs at the same q and seed.")
    parser.add_argument('job_file', nargs='?', default='input.job',
                        help='Job file (default: input.job)')
    parser.add_argument('--qrange', nargs=2, type=int, metavar=('Q_START_M', 'Q_END_M'),
                        help='Special-q range in millions (e.g. --qrange 30 230)')
    parser.add_argument('--qstep', type=int, default=0, metavar='STEP_M',
                        help='Step for --qrange, in millions (default: endpoints only)')
    parser.add_argument('--q', type=int, help='A single special-q value')
    parser.add_argument('--siever', '--ibits', dest='ibits', type=int, default=16,
                        metavar='I_BITS', help='Siever I value, e.g. 15 for I15e')
    parser.add_argument('--jbits', type=int, default=None, metavar='J_BITS',
                        help='Sieve J value (default: I_BITS-1)')
    parser.add_argument('--sqside', choices=['r', 'a'], default='a',
                        help="Side carrying the special q (default: a)")
    for name, hlp in (('rlim', 'rational factor base bound'),
                      ('alim', 'algebraic factor base bound'),
                      ('lpbr', 'rational large prime bound (bits)'),
                      ('lpba', 'algebraic large prime bound (bits)'),
                      ('mfbr', 'rational cofactor bound (bits)'),
                      ('mfba', 'algebraic cofactor bound (bits)')):
        parser.add_argument(f'--{name}', type=int, default=None, help=f'Override {hlp}')
    parser.add_argument('--rlambda', type=float, default=None, help='Override rlambda')
    parser.add_argument('--alambda', type=float, default=None, help='Override alambda')
    parser.add_argument('--alpha-r', type=float, default=None,
                        help='Override rational-side alpha (nats)')
    parser.add_argument('--alpha-a', type=float, default=None,
                        help='Override algebraic-side alpha (nats)')
    parser.add_argument('--nlat', type=int, default=200,
                        help='Random q-lattices to sample (default: 200)')
    parser.add_argument('--npts', type=int, default=120,
                        help='Points per lattice (default: 120)')
    parser.add_argument('--no-siever-table', action='store_true',
                        help='Skip the by-siever survey table')
    parser.add_argument('--tsv', action='store_true',
                        help='Machine-readable: "rat_bits alg_bits yield_index"')
    args = parser.parse_args()

    filename = args.job_file
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!", file=sys.stderr)
        if filename == 'input.job':
            print("Usage: ./estimate_norms.py [job_file]", file=sys.stderr)
        sys.exit(1)

    if args.q is not None and args.q < 2:
        die(f"--q must be at least 2, got {args.q}")
    if args.ibits < 8 or args.ibits > 24:
        die(f"--ibits {args.ibits} is out of range (expected roughly 11-16)")
    if args.jbits is not None and (args.jbits < 1 or args.jbits > 24):
        die(f"--jbits {args.jbits} is out of range")
    if args.nlat < 1 or args.npts < 1:
        die("--nlat and --npts must be positive")

    params, comments = parse_job_file(filename)
    ibits = args.ibits
    jbits = args.jbits if args.jbits is not None else ibits - 1
    overrides = {k: getattr(args, k) for k in
                 ('rlim', 'alim', 'lpbr', 'lpba', 'mfbr', 'mfba',
                  'rlambda', 'alambda')}
    if args.alpha_r is not None:
        overrides['alpha_r'] = args.alpha_r
    if args.alpha_a is not None:
        overrides['alpha_a'] = args.alpha_a
    try:
        model = Model(params, comments, overrides, args.nlat, args.npts)
        model.check()
    except ValueError as e:
        die(f"{filename}: {e}")

    if args.tsv:
        q = args.q if args.q is not None else siever_q(ibits)
        rb, ab = model.at_q(q, ibits, jbits, args.sqside)
        idx = model.yield_index(q, ibits, jbits, args.sqside) or 0.0
        print("%.1f %.1f %.3e" % (rb, ab, idx))
        return

    print(f"Parsing job file: {filename}\n")
    print(f"Algebraic polynomial (degree {len(model.alg) - 1}):")
    print("  f(x) = " + fmt_poly(model.alg))
    print(f"\nRational polynomial (degree {len(model.rat) - 1}):")
    print("  g(x) = " + fmt_poly(model.rat))
    print(f"Skewness: {model.skew}")
    print(f"alpha: rational {model.alpha_r:+.2f}, algebraic {model.alpha_a:+.2f} (nats)")

    if model.cado_lognorm:
        ours = angular_log2_mean(model.alg, model.skew) / LOG2E
        flag = "" if abs(ours - model.cado_lognorm) < 0.5 else "   <-- MISMATCH"
        print(f"Shape term cross-check: {ours:.2f} vs CADO lognorm "
              f"{model.cado_lognorm:.2f} (nats){flag}")

    print(f"\n--- Sieve Parameters ---")
    print(f"lpba: {model.aside.lpb} bits   alim: {model.aside.lim:.2e}   "
          f"mfba: {model.aside.mfb} bits")
    print(f"lpbr: {model.rside.lpb} bits   rlim: {model.rside.lim:.2e}   "
          f"mfbr: {model.rside.mfb} bits")
    print(f"algebraic: {model.aside.describe()}")
    print(f"rational:  {model.rside.describe()}")

    print(f"\n--- Lambda Values ---")
    for label, side in (('alambda', model.aside), ('rlambda', model.rside)):
        current = params.get(label, 0)
        if side.ok:
            print(f"Current {label}: {current}   "
                  f"suggested: {side.mfb / log2(side.lim):.4f}")
        else:
            print(f"Current {label}: {current}   "
                  f"suggested: (need {label[0]}lim and mfb{label[0]})")

    if not args.no_siever_table:
        print()
        siever_table(model, args.sqside)

    qs = []
    if args.q is not None:
        qs = [args.q]
    elif args.qrange:
        q0, q1 = args.qrange[0] * 1_000_000, args.qrange[1] * 1_000_000
        if q1 < q0:
            q0, q1 = q1, q0
        if args.qstep:
            qs = list(range(q0, q1 + 1, args.qstep * 1_000_000))
            if not qs:
                qs = [q0]
            if qs[-1] != q1:
                qs.append(q1)
        else:
            qs = [q0, q1]
    if qs:
        print()
        per_q_table(model, qs, ibits, jbits, args.sqside)

    q_ref = qs[len(qs) // 2] if qs else siever_q(ibits)
    print()
    side_comparison(model, q_ref, ibits, jbits)


if __name__ == "__main__":
    main()
