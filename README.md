# QuantumIQT — Breidbart Attack on IQT1

Simulation and notes for the IQT1 quantum token protocol (Kent et al. 2022)
under a Breidbart intercept-and-resend attack.

---

## What this is

The IQT1 protocol lets a user (Alice) prove ownership of a quantum token at
one of two spacelike-separated points, without storing any quantum state.
Bob (the bank) sends `N` BB84 single photons; Alice measures them, commits
classically, and presents at her chosen point.

This project studies what happens when Alice is **dishonest** (= Eve) and
measures photons in the **Breidbart (π/8) basis** — the optimal intercept-and-resend
strategy for BB84. We analyse both the **full Breidbart attack** (p=1, every pulse)
and the **partial Breidbart attack** (p=3/4, QKD boundary).

---

## Files

| File | Description |
|---|---|
| `iqt_implementation.py` | IQT1 protocol simulation — QBER distributions, Panels A & B |
| `final_analysis.py` | Full analysis — P_forge vs delta, ε_unf vs N, I_BE vs N |
| `breidbart_notes.pdf` | Self-contained lecture notes (6 pages) |
| `breidbart_notes.tex` | LaTeX source for the notes |

---

## Key quantities

| Symbol | Value | Nature |
|---|---|---|
| `p_BE = 1/2 + 1/(2√2)` | 0.8536 | Born rule — Breidbart success prob per qubit |
| `e_BE = sin²(π/8)` | 14.64% | QBER Bob observes in IQT (Eve reports directly) |
| `Q_tol` | 11.00% | QKD threshold — solution of `h(q) = 1/2` |
| margin = `e_BE - Q_tol` | 3.64% | Fundamental gap |
| `delta = e_BE - gamma_err` | (0, e_BE) | Bank's strictness — chosen by designer |
| `gamma_err` | [0, Q_tol] | Bank's acceptance threshold — chosen (protocol) |
| `n_min` | ~3,400 at delta=3.64% | Min photons for P_forge negligible |
| `ε_unf = p_BE^N` | ~10⁻⁶⁹ at N=1000 | Paper eq. (D1) — full quantum bound |

---

## Security chain

```
Q_tol (11%)  <  e_BE (14.64%)  <  25% (individual attack)
```

- margin = 3.64% → gap between Eve's error rate and the detection boundary
- `ε_unf = p_BE^N` → forging probability exponentially suppressed in N
- At N=1000: Eve would need 10⁶⁹ attempts — more than atoms on Earth

---

## Three-panel analysis (`final_analysis.py`)

**Panel A** — `P_forge(delta, n)` vs `delta`: how strict must the bank be?

**Panel B** — `ε_unf = p_BE^N` vs N (eq. D1): the paper's quantum bound.

**Panel C** — `I_BE` (Eve's total MI) vs N, with `ε_unf` scale on right axis:
Eve accumulates ~200 bits at N=1000 yet `ε_unf = 10⁻⁶⁹`.
Information does not equal forgery.

---

## QBER: IQT vs QKD

| Scenario | QBER | Formula |
|---|---|---|
| QKD (Eve re-sends photon) | 25% | `1 - (p_BE² + e_BE²)` |
| IQT (Eve reports directly) | 14.64% | `e_BE = sin²(π/8)` |

---

## Mutual information

| Attack | I_BE per qubit | Total at N=10,000 |
|---|---|---|
| Full Breidbart (p=1) | 0.399 bits/qubit | ~2,000 bits |
| Partial Breidbart (p=3/4) | 0.299 bits/qubit | ~1,497 bits |

---

## Minimum token size

```
n_min = [Φ⁻¹(P_forge_target)]² * e_BE*(1-e_BE) / delta²
```

Scales as `1/delta²` — halving `delta` quadruples `n_min`.

| delta | gamma_err | n_min (P_forge=10⁻⁹) |
|---|---|---|
| 3.64% (margin) | 11% = Q_tol | ~3,400 |
| 7.32% = e_BE/2 | 7.32% | ~800 |
| 14.64% = e_BE | 0% (ideal) | ~210 |

---

## Quick start

```bash
pip install numpy scipy matplotlib
python final_analysis.py      # main analysis + terminal table
python iqt_implementation.py  # protocol simulation
```

---

## Reference

A. Kent, D. Lowndes, D. Pitalúa-García, J. Rarity —
*Practical quantum tokens without quantum memories and experimental tests*,
arXiv:2104.11717v4 (2022).

Equation (D1): `ε_unf = p_BE^N` — the Breidbart attack is optimal and its
per-qubit success rate directly gives the unforgeability bound.
