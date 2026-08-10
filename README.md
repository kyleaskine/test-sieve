# GNFS Test Sieve Tools

Tools for testing and optimizing GNFS (General Number Field Sieve) sieving parameters.

## Setup

### Prerequisites

1. **GGNFS binaries**: You need compiled versions of the GGNFS lattice sievers
   - Required versions: `gnfs-lasieve4I14e`, `gnfs-lasieve4I15e`, `gnfs-lasieve4I16e`
   - Build from source: [lasieve in YAFU on GitHub](https://github.com/bbuhrow/yafu)
   - Copy the compiled binaries into this directory

2. **Job file**: Create an `input.job` file containing your polynomial and sieve parameters
   - Must include polynomial coefficients (c0-c6, Y0, Y1)
   - Must include skewness and sieve parameters (alim, rlim, lpba, lpbr, mfba, mfbr, alambda, rlambda)
   - See example format in existing input.job

### File Structure

After setup, your directory should look like:
```
test-sieve/
├── gnfs-lasieve4I14e     # Binary executable
├── gnfs-lasieve4I15e     # Binary executable
├── gnfs-lasieve4I16e     # Binary executable
├── input.job             # Your polynomial and parameters
├── estimate_norms.py     # Norm estimation tool
└── test_sieve.sh         # Test sieve runner
```

## Usage

### 1. Estimate Norms (estimate_norms.py)

Estimates the average norm on each side over the region the siever actually
covers for a given special-q, and prices the trade-off between which side
carries the special q and how the large-prime parameters are set.

```bash
./estimate_norms.py [job_file] [--ibits 15] [--jbits 14] [--sqside a|r]
                    [--q Q | --qrange START_M END_M [--qstep STEP_M]]
                    [--lpbr N --lpba N --mfbr N --mfba N --rlambda X ...]
```

- If no file is specified, defaults to `input.job`
- `--ibits` is the siever I value (15 for `gnfs-lasieve4I15e`); `--jbits`
  defaults to `ibits-1`, so pass `--jbits 16` when sieving I16 with `-J 16`
- `Raw rat`/`Raw alg` are the side-independent norms before choosing which
  side carries the special q. `Sq rat`/`Sq alg` are the effective norms after
  subtracting `log2(q)` from the side selected by `--sqside`.

**Norms vs. yield -- two separate layers.** The norm depends only on the
polynomials, skew, `I`/`J`, `q`, and which side carries `q`. `lpb`, `mfb`,
`lambda`, `alim`/`rlim` do not enter it at all; they only change the
probability that a norm of a given size becomes a relation. That probability
is the `yield idx`, and it is where 2LP vs 3LP shows up.

**Example output:**

With no `--q`/`--qrange`, a survey table shows each siever at the
representative q for its size (1e4 for I11 rising a decade per bit to 1e9 for
I16), so you are not reading a single arbitrary q. Skip it with
`--no-siever-table`.

```
--- Norm Estimates by Siever (special q on the algebraic side) ---
Siever                    q  Raw rat  Raw alg   Sq rat   Sq alg  Alg-Rat   yield idx  Best sq side
-----------------------------------------------------------------------------------------------------
gnfs-lasieve4I14e       1e7    160.2    144.2    160.2    120.9    -39.3    1.06e-06  algebraic
gnfs-lasieve4I15e       1e8    162.9    160.3    162.9    133.7    -29.2    2.71e-07  algebraic
gnfs-lasieve4I16e       1e9    165.6    176.1    165.6    146.2    -19.4    6.44e-08  algebraic

--- Special-q side comparison (I=2^16, J=2^15, q = 1,000,000,000) ---
  rational side: 3LP, cofactor <= 89 bits
  algebraic side: 2LP, cofactor <= 59 bits
  raw norms (before special q): rational 165.6 bits, algebraic 176.1 bits

sq side       rat bits  alg bits    yield idx
----------------------------------------------
rational         135.7     176.1     4.63e-08
algebraic        165.6     146.2     6.44e-08

=> Model prefers special q on the ALGEBRAIC side (1.39x the yield index).
```

`Best sq side` can change with q, so the survey is worth a glance before
committing to a side for a long range.

**Pricing a compensation trade.** When the norms are nearly even (typical for
degree-6 SNFS) you may want to sieve one side and give the *other* side 3LP and
a higher `lpb` to make up for it. The overrides let you cost that out; indices
are comparable across runs at the same `q` and seed:

```bash
# sieve algebraic, rational stays as-is
./estimate_norms.py --ibits 15 --q 130000000 --sqside a --tsv
# sieve rational instead, compensating that side with 3LP and lpb+1
./estimate_norms.py --ibits 15 --q 130000000 --sqside r --lpbr 32 --mfbr 93 --tsv
```

Note that `mfb` only buys an extra large prime once it clears `k*log2(lim)` --
with `alim` at 1.34e8, three large primes need at least 81 bits, so `mfba` 64
and 80 are worth exactly the same. `lambda` can also bind before `mfb` does
(the effective cap is `min(mfb, lambda*log2(lim))`), silently discarding the
relations `mfb` was meant to allow; the tool reports the effective cap.

`yield idx` is a heuristic: validated to track *gains* from adding a large
prime or raising `lpb` to within 1-4%, but optimistic by orders of magnitude in
absolute terms. Compare indices, never read one on its own. See the docstring
in `estimate_norms.py` for the model and its validation.

### 2. Test Sieve (test_sieve.sh)

Runs test sieves across specified q0 ranges to measure yields and optimize parameters.

```bash
./test_sieve.sh
```

The script will prompt for:
- Side to sieve (algebraic or rational)
- Siever version (14, 15, or 16)
- Start value in millions (e.g., 100 for 100M)
- End value in millions
- Interval in millions

**Outputs:**
- Detailed sieve statistics for each q0 value
- Expected relation counts (exp_rel)
- Side-independent `raw-rat`/`raw-alg` norms per q0, via
  `estimate_norms.py` (skipped if `python3` is unavailable). Run the estimator
  directly for the full special-q-adjusted comparison.
- Suggested lambda values
- Creates `result.job` with optimal parameters

## Tips

- Run `estimate_norms.py` first to determine which side to sieve; it is a
  cheap average-case model, so let the measured yields overrule it
- Use smaller intervals (1M-10M) for initial testing
- Compare yields across different siever versions (14 vs 15 vs 16)
- Adjust lambda values based on suggested output
- Monitor sec/rel to balance speed vs yield

## Example Workflow

```bash
# 1. Create your input.job file with polynomial and parameters

# 2. Estimate which side to sieve
./estimate_norms.py

# 3. Run test sieves
./test_sieve.sh
# Choose rational/algebraic based on step 2
# Select siever version (usually 15 is good balance)
# Test range, e.g., 100M to 150M with 10M intervals

# 4. Review output and adjust parameters in input.job if needed
```
