# QuantumIQT — Breidbart Attack on IQT1

Simulation and notes for the IQT1 quantum token protocol (Kent et al. 2022)
under a Breidbart intercept-and-resend attack.

---

## What this is

The IQT1 protocol lets a user (Alice) prove ownership of a quantum token at
one of two spacelike-separated points, without storing any quantum state.
Bob (the bank) sends `N` BB84 single photons; Alice measures them, commits
classically, and presents at her chosen point.

This project studies what happens when Alice is **dishonest** and measures
every photon in the **Breidbart (π/8) basis** — the optimal intercept-and-resend
strategy for BB84 — instead of a random BB84 basis.

---

## Files

| File | Description |
|---|---|
| `iqt_implementation.py` | Full IQT1 protocol simulation + Panels A & B |
| `iqt_panels_ABC.py` | Extended simulation including eq. (D1) unforgeability bound (Panel C) |
| `breidbart_notes.pdf` | Self-contained lecture notes (7 pages) |
| `breidbart_notes.tex` | LaTeX source for the notes |

---

## Key quantities

| Symbol | Value | How obtained |
|---|---|---|
| `p_BE = 1/2 + 1/(2√2)` | 0.8536 | Born rule — Breidbart success prob per qubit |
| `e_BE = 1/2 - 1/(2√2)` | 0.1464 | QBER Bob observes in IQT (Eve reports directly) |
| `gamma_err` | ~5% | **Chosen** by the bank as acceptance threshold |
| `Q_tol` | 17.05% | **Calculated** from `0.5*(1-h(q)) = q` |
| `ε_unf = p_BE^N` | ~10⁻⁶⁹ at N=1000 | Paper eq. (D1), Lemma 1 |

The chain of inequalities that makes the protocol secure:

```
gamma_err (5%)  <  e_BE (14.64%)  <  Q_tol (17.05%)
```

- `e_BE > gamma_err` → Breidbart attack **always detected** by the bank's QBER check  
- `e_BE < Q_tol` → Eve is in the information-advantage zone, but the `gamma_err` check catches her first  
- `ε_unf = p_BE^N` → forging probability **exponentially suppressed** in N

---

## Protocol flow (IQT1, Kent et al. Table I)

```
Stage I
  Step 1.  Bob generates N BB84 photons |ψ_k⟩ = |φ_{t_k u_k}⟩
  Step 2.  Alice measures each photon
             honest  → random BB84 basis y_k,  outcome x_k
             Breidbart Eve → π/8 basis,  outcome x_k (correct w.p. p_BE)
  Steps 3-6. Classical commitment: Alice sends d = y ⊕ z to Bob

Stage II
  Steps 7-9. Alice reveals b (chosen point) and c = b ⊕ z
  Step 10.   Alice's agent presents token x to Bob's agent at Q_b
  Step 11.   Bob validates: checks x[Δ_b] == t[Δ_b]
             where Δ_b = {k : ỹ_{b,k} = u_k}  (≈ N/2 positions)
             rejects if QBER > gamma_err
```

Key algebraic fact: after cancellation of z and c, `d̃_b = y` always.
So `Δ_b = {k : y_k = u_k}` — positions where Alice claimed Bob's basis.

---

## QBER: IQT vs QKD

The Breidbart QBER depends on whether Eve re-sends a photon (QKD) or
reports outcomes directly (IQT):

| Scenario | QBER | Formula |
|---|---|---|
| QKD (Eve re-sends photon to Bob) | 25% | `1 - (p_BE² + e_BE²)` |
| IQT (Eve reports outcomes directly) | 14.64% | `e_BE = 1/2 - 1/(2√2)` |

In IQT, Bob checks Alice's classical outcomes `x[Δ_b]` against `t[Δ_b]`
directly — there is no second photon transmission — so the relevant QBER
is `e_BE`, not 0.25.

---

## Simulation results

Run with `N = 2000`, `n_trials = 1000`, `gamma_err = 5%`:

```
Honest Alice     QBER = 0.0000  → always accepted
Breidbart Eve    QBER = 0.1462 ± 0.0117  → always rejected
Theory e_BE             = 0.1464  ✓
```

Forging probability (both sub-tokens pass simultaneously):

```
N =   50 :  P_forge ≈ 8%      (finite-size fluctuations)
N =  200 :  P_forge ≈ 0.2%
N =  300 :  P_forge = 0       (consistent with 10⁻²¹)
N = 1000 :  P_forge = 0       (consistent with 10⁻⁶⁹)
```

---

## Quick start

```bash
pip install numpy scipy matplotlib
python iqt_implementation.py   # produces iqt_panels_AB.png
python iqt_panels_ABC.py       # produces iqt_panels_ABC.png (includes D1 bound)
```

---

## Reference

A. Kent, D. Lowndes, D. Pitalúa-García, J. Rarity —
*Practical quantum tokens without quantum memories and experimental tests*,
arXiv:2104.11717v4 (2022).

Equation (D1): `ε_unf = (1/2 + 1/(2√2))^N`  — proved in Appendix D via
reduction to a state-discrimination task; the optimal strategy is exactly
the Breidbart measurement.
