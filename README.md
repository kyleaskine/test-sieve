# GNFS Test Sieve Tools

Tools for testing and optimizing GNFS (General Number Field Sieve) sieving parameters.

## Setup

### Prerequisites

1. **GGNFS binaries**: You need compiled versions of the GGNFS lattice sievers
   - Required versions: `gnfs-lasieve4I14e`, `gnfs-lasieve4I15e`, `gnfs-lasieve4I16e`
   - Build from source: [GGNFS on SourceForge](https://sourceforge.net/projects/ggnfs/)
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

Analyzes your polynomial to determine which side (algebraic or rational) to sieve.

```bash
./estimate_norms.py [job_file]
```

- If no file is specified, defaults to `input.job`
- Shows norm estimates for each siever version (14, 15, 16)
- Recommends which side to sieve based on norm sizes

**Example output:**
```
Siever               Alg bits   Rat bits     Diff Sieve side
--------------------------------------------------------------
gnfs-lasieve4I14e        92.3       95.1     +2.8 rational
gnfs-lasieve4I15e        94.5       97.3     +2.8 rational
gnfs-lasieve4I16e        96.7       99.5     +2.8 rational

=> Algebraic side has SMALLER norms
   Sieve the RATIONAL side
```

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
- Suggested lambda values
- Creates `result.job` with optimal parameters

## Tips

- Run `estimate_norms.py` first to determine which side to sieve
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
