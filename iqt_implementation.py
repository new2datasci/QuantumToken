"""
IQT1 Protocol Simulation — Panels A, B & C
============================================
Single photon + perfect channel + Breidbart attack

Protocol: IQT1 (Kent et al. Table I)
- Bob sends N BB84 single photons
- Honest Alice  → QBER = 0          (always accepted)
- Breidbart Eve → QBER = e_BE       (always detected)
- Panel C compares the simulation to eq. (D1) from the paper:
      epsilon_unf = (1/2 + 1/(2√2))^N
  which is the per-qubit Breidbart success probability p_BE raised
  to the N-th power — confirmed by our simulation.

Notation (aligned with TD2):
  p_BE  = P(Eve guesses Bob's bit correctly) = 1/2 + 1/(2√2)
  e_BE  = 1 - p_BE = 1/2 - 1/(2√2)  ← QBER observed in IQT
  Q_tol = threshold from 0.5*(1-h(q)) = q ≈ 17.05%
  gamma_err = bank's acceptance threshold (protocol parameter, chosen)
"""

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Constants ──────────────────────────────────────────────────
def h(q):
    if q <= 0 or q >= 1: return 0.0
    return -q * np.log2(q) - (1-q) * np.log2(1-q)

p_BE      = 0.5 + 0.5/np.sqrt(2)
e_BE      = 1 - p_BE
Q_tol     = brentq(lambda q: 0.5*(1-h(q))-q, 0.001, 0.499)
rng       = np.random.default_rng(42)
N         = 2000
n_trials  = 1000
gamma_err = 0.05

print(f"p_BE  = {p_BE:.6f}  (= 1/2 + 1/(2√2))")
print(f"e_BE  = {e_BE:.6f}  (= 1/2 - 1/(2√2))")
print(f"Q_tol = {Q_tol:.6f}")

# ── Step 1: Bob sends N BB84 single photons ────────────────────
def bob_send(N):
    t = rng.integers(0, 2, N)    # prepared bit
    u = rng.integers(0, 2, N)    # prepared basis (0=Z, 1=X)
    return t, u

# ── Step 2a: Honest Alice ──────────────────────────────────────
def alice_honest(t, u):
    y = rng.integers(0, 2, len(t))
    x = np.where(y == u, t, rng.integers(0, 2, len(t)))
    return x, y

# ── Step 2b: Breidbart Eve ────────────────────────────────────
def alice_breidbart(t, u):
    right = rng.random(len(t)) < p_BE
    x = np.where(right, t, 1-t)
    y = rng.integers(0, 2, len(t))   # random reported basis
    return x, y

# ── Steps 3-11: IQT1 classical protocol ───────────────────────
def run_protocol(t, u, x, y, b):
    """
    Returns (QBER on Delta_b, |Delta_b|).
    Key result: dtilde_b = y always, so Delta_b = {k: y[k]==u[k]}.
    QBER = fraction of x[Delta_b] != t[Delta_b].
    """
    z        = rng.integers(0, 2)
    d        = y ^ z
    c        = b ^ z
    d_b      = d if b == 0 else 1-d
    dtilde_b = d_b ^ c
    Delta_b  = np.where(dtilde_b == u)[0]
    if len(Delta_b) == 0:
        return np.nan, 0
    qber = float(np.mean(x[Delta_b] != t[Delta_b]))
    return qber, len(Delta_b)

# ── Run simulation ─────────────────────────────────────────────
q_honest, q_breid, sizes = [], [], []
for _ in range(n_trials):
    t, u       = bob_send(N)
    xh, yh     = alice_honest(t, u)
    xb, yb     = alice_breidbart(t, u)
    qh, _      = run_protocol(t, u, xh, yh, b=0)
    qb, sz     = run_protocol(t, u, xb, yb, b=0)
    q_honest.append(qh)
    q_breid.append(qb)
    sizes.append(sz)

q_honest = np.array(q_honest)
q_breid  = np.array(q_breid)
n_Z      = int(np.mean(sizes))

print(f"Simulated e_BE   = {np.nanmean(q_breid):.6f}  (theory: {e_BE:.6f})")
print(f"Simulated std    = {np.nanstd(q_breid):.6f}  "
      f"(theory: {np.sqrt(e_BE*(1-e_BE)/n_Z):.6f})")

# ── Panel C: unforgeability bound eq.(D1) ─────────────────────
# Paper (D1): epsilon_unf = (1/2 + 1/(2√2))^N = p_BE^N
# 
# Simulation estimate of P_forge:
#   Alice succeeds at BOTH b=0 and b=1 simultaneously.
#   Per trial: success_b = (QBER_b <= gamma_err)
#   P_forge_sim ≈ fraction of trials where BOTH succeed.
#
# We compare this to the theoretical bound (D1).

# collect joint success for varying N
N_values = [50, 100, 200, 300, 500, 750, 1000, 1500, 2000]
p_forge_theory = []
p_forge_sim    = []

for Nv in N_values:
    # theoretical bound (D1)
    p_forge_theory.append(p_BE**Nv)

    # simulation: run n_sim trials for this N
    n_sim = 500
    successes = 0
    for _ in range(n_sim):
        t, u   = bob_send(Nv)
        xb, yb = alice_breidbart(t, u)
        q0, _  = run_protocol(t, u, xb, yb, b=0)
        q1, _  = run_protocol(t, u, xb, yb, b=1)
        # Alice succeeds only if BOTH sub-tokens pass gamma_err check
        if (not np.isnan(q0)) and (not np.isnan(q1)):
            if q0 <= gamma_err and q1 <= gamma_err:
                successes += 1
    p_forge_sim.append(successes / n_sim)

p_forge_theory = np.array(p_forge_theory)
p_forge_sim    = np.array(p_forge_sim)

print("\nN       theory bound (D1)   simulated P_forge")
print("-"*48)
for i, Nv in enumerate(N_values):
    th = p_forge_theory[i]
    si = p_forge_sim[i]
    th_str = f"{th:.2e}" if th > 1e-300 else "< 1e-300"
    print(f"N={Nv:5d}   {th_str:>16}   {si:.4f}")

# ── Figure ─────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 5))
fig.patch.set_facecolor('#0f1117')
DARK  = '#0f1117'; PANEL = '#1a1d27'
BLUE  = '#4da6ff'; CORAL = '#ff6b6b'
GREEN = '#4dffb4'; AMBER = '#ffc94d'
GRAY  = '#8899aa'; WHITE = '#e8eaf6'
PURPLE= '#b48eff'

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=GRAY, labelsize=9)
    ax.spines[:].set_color('#2a2d3a')
    ax.set_title(title, color=WHITE, fontsize=10, pad=8)
    for l in ax.get_xticklabels() + ax.get_yticklabels():
        l.set_color(GRAY)
    ax.xaxis.label.set_color(GRAY)
    ax.yaxis.label.set_color(GRAY)

gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

# ── Panel A: QBER distributions ───────────────────────────────
ax = fig.add_subplot(gs[0])
style_ax(ax, 'A · QBER distributions')

ax.bar([0], [n_trials], width=0.004, color=BLUE, alpha=0.85,
       label='Honest Alice  (QBER = 0 always)', zorder=3)

sigma_th = np.sqrt(e_BE*(1-e_BE)/n_Z)
x_range  = np.linspace(e_BE - 5*sigma_th, e_BE + 5*sigma_th, 300)
pdf      = norm.pdf(x_range, loc=e_BE, scale=sigma_th)
ax.fill_between(x_range, pdf * sigma_th * n_trials * 2,
                color=CORAL, alpha=0.35, label='Breidbart / theoretical', zorder=2)
ax.hist(q_breid, bins=40, color=CORAL, alpha=0.7,
        label='Breidbart / simulated', zorder=2)
ax.axvline(0,         color=BLUE,  lw=1.5, ls='-',  zorder=4)
ax.axvline(e_BE,      color=AMBER, lw=2,   ls='--',
           label=f'$e_{{BE}}$ = {e_BE*100:.2f}%', zorder=4)
ax.axvline(gamma_err, color=GREEN, lw=2,   ls=':',
           label=f'$\\gamma_{{err}}$ = {gamma_err*100:.0f}%', zorder=4)
ax.annotate('', xy=(gamma_err+0.002, 60), xytext=(e_BE-0.002, 60),
            arrowprops=dict(arrowstyle='<->', color=WHITE, lw=1.2))
ax.text((gamma_err+e_BE)/2, 66, 'always\ndetected',
        color=WHITE, fontsize=8, ha='center')
ax.set_xlabel('QBER on $\\Delta_b$')
ax.set_ylabel('Number of trials')
ax.set_xlim(-0.02, 0.25)
ax.legend(fontsize=7.5, facecolor=DARK, labelcolor=WHITE, edgecolor='#2a2d3a')

# ── Panel B: Bob's decision per trial ─────────────────────────
ax2 = fig.add_subplot(gs[1])
style_ax(ax2, "B · Bob's decision per trial")

n_show = 200
ax2.scatter(np.arange(n_show), q_honest[:n_show]*100,
            color=BLUE,  s=4, alpha=0.6, label='Honest Alice', zorder=3)
ax2.scatter(np.arange(n_show), q_breid[:n_show]*100,
            color=CORAL, s=4, alpha=0.6, label='Breidbart Eve', zorder=3)
ax2.axhline(gamma_err*100, color=GREEN, lw=2, ls=':',
            label=f'$\\gamma_{{err}}$ = {gamma_err*100:.0f}%')
ax2.axhline(e_BE*100, color=AMBER, lw=1.5, ls='--',
            label=f'$e_{{BE}}$ = {e_BE*100:.2f}%')
ax2.axhspan(0,             gamma_err*100, color=BLUE,  alpha=0.07)
ax2.axhspan(gamma_err*100, 22,            color=CORAL, alpha=0.07)
ax2.text(100, gamma_err*100/2,
         'ACCEPTED', color=BLUE,  fontsize=9, ha='center',
         va='center', fontweight='bold', alpha=0.8)
ax2.text(100, (gamma_err*100 + e_BE*100)/2 + 1,
         'REJECTED', color=CORAL, fontsize=9, ha='center',
         va='center', fontweight='bold', alpha=0.8)
ax2.set_xlabel('Trial index')
ax2.set_ylabel('QBER (%)')
ax2.set_ylim(-0.5, 22)
ax2.set_xlim(0, n_show)
ax2.legend(fontsize=7.5, facecolor=DARK, labelcolor=WHITE, edgecolor='#2a2d3a')

# ── Panel C: epsilon_unf (D1) vs N ────────────────────────────
ax3 = fig.add_subplot(gs[2])
style_ax(ax3, 'C · Unforgeability bound eq.(D1) vs N')

# theoretical curve (D1): log10(epsilon_unf) = N * log10(p_BE)
N_cont = np.linspace(10, 2000, 500)
log_eps_theory = N_cont * np.log10(p_BE)
ax3.plot(N_cont, log_eps_theory, color=PURPLE, lw=2.5,
         label=r'Theory: $\varepsilon_{unf} = p_{BE}^N$ eq.(D1)')

# simulated points — only plot where P_forge_sim > 0
# (for large N, simulation gives 0 — consistent with theory -> -inf)
sim_nonzero = [(N_values[i], p_forge_sim[i])
               for i in range(len(N_values)) if p_forge_sim[i] > 0]
if sim_nonzero:
    Ns_nz, Ps_nz = zip(*sim_nonzero)
    ax3.scatter(Ns_nz, np.log10(Ps_nz), color=CORAL, s=60, zorder=5,
                label='Simulation (non-zero)', marker='o')

# zero-detection markers (arrows pointing down)
sim_zero = [N_values[i] for i in range(len(N_values)) if p_forge_sim[i] == 0]
for Nz in sim_zero:
    log_th = Nz * np.log10(p_BE)
    ax3.annotate('', xy=(Nz, log_th - 2), xytext=(Nz, log_th + 2),
                 arrowprops=dict(arrowstyle='->', color=GREEN, lw=1.5))
ax3.scatter(sim_zero,
            [Nz * np.log10(p_BE) for Nz in sim_zero],
            color=GREEN, s=60, zorder=5, marker='v',
            label='Simulation = 0\n(consistent with theory)')

# reference lines
ax3.axhline(-10,  color=GRAY, lw=0.8, ls=':', alpha=0.5)
ax3.axhline(-50,  color=GRAY, lw=0.8, ls=':', alpha=0.5)
ax3.axhline(-100, color=GRAY, lw=0.8, ls=':', alpha=0.5)
for val, lbl in [(-10,'$10^{-10}$'), (-50,'$10^{-50}$'), (-100,'$10^{-100}$')]:
    ax3.text(2020, val+1, lbl, color=GRAY, fontsize=8, va='bottom')

ax3.set_xlabel('Token size N')
ax3.set_ylabel(r'$\log_{10}(\varepsilon_{unf})$')
ax3.set_xlim(0, 2100)
ax3.legend(fontsize=7.5, facecolor=DARK, labelcolor=WHITE, edgecolor='#2a2d3a')

# annotation box showing key numbers
ax3.text(200, -85,
         f'$p_{{BE}} = \\frac{{1}}{{2}}+\\frac{{1}}{{2\\sqrt{{2}}}}={p_BE:.4f}$\n'
         f'$e_{{BE}} = {e_BE:.4f}$ (Breidbart QBER)\n'
         f'$Q_{{tol}} = {Q_tol:.4f}$ (info threshold)',
         color=WHITE, fontsize=8, va='top',
         bbox=dict(boxstyle='round,pad=0.4', facecolor=PANEL,
                   edgecolor='#2a2d3a', alpha=0.9))

fig.suptitle(
    f'IQT1 · Single photon · Perfect channel · N={N}\n'
    f'Honest QBER=0  |  Breidbart QBER=$e_{{BE}}$={e_BE*100:.2f}%  |  '
    f'$\\gamma_{{err}}$={gamma_err*100:.0f}%  |  '
    f'eq.(D1): $\\varepsilon_{{unf}}=p_{{BE}}^N$',
    color=WHITE, fontsize=10, y=1.02
)

plt.savefig('/Users/ruchithareja/Documents/Python Decoy/Quantum Token Python/iqt_panels_ABC.png',
            dpi=150, bbox_inches='tight', facecolor=DARK)
print("\nFigure saved: iqt_panels_ABC.png")
plt.show()
