[README.md](https://github.com/user-attachments/files/28675413/README.md)
# Q-quadratic convergence of the centralized circumcentered-reflection method

Companion code and figures for the paper

> Y. Bello-Cruz,
> *Q-quadratic convergence of the centralized circumcentered-reflection
> method under a relative interior condition.*

The experiments are implemented in **Julia**; a **Python** port reproduces the
same numbers and produces the figures. All scripts regenerate every numerical
claim and figure in the paper.

## What the code shows

For two equality-constrained feasibility problems — overlapping ellipsoids in a
plane, and a fixed-trace spectral feasibility problem in `S^8` — the classical
full-dimensional hypothesis of [Behling–Bello-Cruz–Iusem–Santos, *Math.
Program.* 2024] fails, because the constraints confine the sets to a proper
affine subspace. Only the relative-interior analysis of the paper applies.
The centralized circumcentered-reflection method (cCRM) converges
**Q-quadratically**: the ratio `dist(z^{k+1}) / dist(z^k)^2` stabilizes at a
finite value below the predicted bound `2 max(kappa_X, kappa_Y) / omega`, and
the method reaches machine precision in a handful of steps where the method of
alternating projections (MAP) and Douglas–Rachford (DR) take many.

## Repository layout

```
.
├── julia/                  # primary implementation (Julia, stdlib only)
│   ├── verify_numerics.jl  # high-precision checks of the lemmas and examples
│   └── run_experiments.jl  # cCRM vs MAP vs DR on both test problems
├── python/                 # Python port (NumPy/SciPy/Matplotlib)
│   ├── verify_numerics.py
│   ├── experiments_ccrm.py # projections, circumcenter, cCRM/MAP/DR operators
│   └── run_experiments.py  # runs the experiments and writes the figures
├── figures/                # figures used in the paper
│   ├── convergence_ellipsoid.png
│   ├── convergence_spectral.png
│   └── ccrm_convergence.png
└── paper/
    └── cCRM_revised.pdf
```

## Running

### Julia (primary)

```bash
cd julia
julia verify_numerics.jl     # lemma / example checks (BigFloat precision)
julia run_experiments.jl     # cCRM vs MAP vs DR; prints distance histories
```

Requires Julia ≥ 1.6. Uses only the standard library (`LinearAlgebra`,
`Printf`); see `julia/Project.toml`.

### Python (port, also produces the figures)

```bash
cd python
pip install -r requirements.txt
python verify_numerics.py    # high-precision checks (needs mpmath)
python run_experiments.py    # writes convergence_*.png into the working dir
```

Requires Python ≥ 3.9 with `numpy`, `scipy`, `matplotlib`, `mpmath`.

## Reproducibility

The spectral experiment uses a fixed random seed (`20260606` in Python) so the
reported numbers reproduce exactly. The Julia port uses Julia's default RNG;
set a seed with `Random.seed!` if bit-identical output across languages is
needed. All distances are `dist(z^k, X ∩ Y)`, computed by projecting each
iterate onto the intersection with Dykstra's algorithm.

## Citation

If you use this code, please cite the paper and Julia:

- Y. Bello-Cruz, *Q-quadratic convergence of the centralized
  circumcentered-reflection method under a relative interior condition.*
- J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah,
  *Julia: a fresh approach to numerical computing*, SIAM Rev. 59 (2017) 65–98.

## License

MIT. See [LICENSE](LICENSE).
