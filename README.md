# Difference-of-Log-Normals distribution

MATLAB functions for working with the Difference-of-Log-Normals (DLN) distribution: the difference between two jointly log-normal random variables. The distribution supports positive and negative outcomes and is used in my research on firm fundamentals, growth and returns.

## Functions

| Function | Purpose |
|---|---|
| `dlnpdf`, `dlncdf`, `dlninv` | Density, cumulative distribution and inverse CDF |
| `dlnrnd` | Random draws |
| `dlnmom` | Moments |
| `dlnfit` | Parameter fitting |
| `dlnpar` | Parameter transformations |
| `mvdlnpdf`, `mvdlncdf`, `mvdlnrnd` | Multivariate DLN functions |

Add the repository directory to your MATLAB path. Each function's source header documents its arguments and parameter conventions.

The univariate parameter order is `[muP, sigP, muN, sigN, rPN]`: the means and standard deviations of the underlying normals, followed by their correlation.

Illustrative calls, checked against the function signatures but not executed in this documentation update:

```matlab
p = [0, 1, 0, 1, 0];
w = dlnrnd(p(1), p(2), p(3), p(4), p(5), [1, 1000]);
f = dlnpdf([-1, 0, 1], p);
```

Dependencies vary by function. Distribution and simulation routines use MATLAB statistical functions. `dlnfit` uses Optimization Toolbox and Global Optimization Toolbox functions and enables parallel multistart. It requires at least 100 positive and 100 negative input observations. Consult the source before choosing the execution environment; this update does not establish a supported MATLAB release or validate numerical behavior.

## Research

- [The Difference-of-Log-Normals Distribution: Properties, Estimation, and Growth](https://arxiv.org/abs/2302.02486) — 2023 technical note on the distribution and its estimators.
- [Growth and Differences of Log-Normals](https://kn.owled.ge/research/#growth) — related research; selected public manuscript and abstract.
- [Research software overview](https://kn.owled.ge/software/#dln).

The code snapshot and linked papers have different dates; this is not a claim that the repository reproduces every result in the later manuscripts. Existing license files and notices are retained.
