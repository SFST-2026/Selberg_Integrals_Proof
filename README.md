# Selberg A₂ Moment Verification Scripts

Scripts verifying the theorems in "Exact Trigonometric Moments of the A₂ Selberg Integral."

## Scripts

| Script | Verifies | Runtime |
|--------|----------|---------|
| `selberg_moments_verify.py` | Theorems 2–3 (7 formulas, k=1..30); Table 1 | ~5 min |
| `dixon_selberg_verify.py` | Theorem 1 (Dixon, k=1..15) | ~30 s |
| `selberg_30digit_verify.py` | Conjecture 6 (30-digit central equation) | ~30 s |
| `selberg_verify_20digits.py` | Full pipeline (20 digits, detailed output) | ~5 s |
| `selberg_precision_results.py` | Table 2 (collected 30-digit values) | reference |
| `selberg_heat_kernel_analysis.py` | Heat kernel analysis (needs GPU data) | varies |

## Requirements

- Python 3.8+
- `mpmath` >= 1.3 (`pip install mpmath`)

## Numerical Parameters

- Gauss–Legendre: N = 60–120 points per dimension
- Internal precision: `mp.dps` = 40–90 decimal digits
- Lattice cutoff: N_max = 10 (307 shells on Z⁴)
- Fourier modes: M = 12–15

## Usage

```bash
python3 selberg_moments_verify.py     # Verify 7 formulas
python3 dixon_selberg_verify.py       # Verify Dixon identity
python3 selberg_30digit_verify.py     # 30-digit central equation
```

Each script is self-contained and prints explicit pass/fail output.

## Author

M.W. Le Borgne, Independent Researcher, Schmitten, Germany

Computational assistance by Claude (Anthropic). All results independently verified.
