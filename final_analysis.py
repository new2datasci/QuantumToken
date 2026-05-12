"""
QuantumIQT — Final Analysis
============================
Single photon, no loss, full channel.

Three panels:
  A — P_forge(delta, n) vs delta  [bank varies gamma_err = e_BE - delta]
  B — eps_unf = p_BE^N vs N       [paper's bound, eq. D1]
  C — Mutual information I_BE     [full vs partial Breidbart, vs N]
"""

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ══════════════════════════════════════════════════════════════
# Protocol parameters
# ══════════════════════════════════════════════════════════════
N_PAPER   = 742_491          # Alice detections (Kent et al.)
N_Z_PAPER = N_PAPER // 2    # sub-token size

# ══════════════════════════════════════════════════════════════
# Derived constants
# ══════════════════════════════════════════════════════════════
p_BE = 0.5 + 0.5/np.sqrt(2)
e_BE = 1 - p_BE              # = sin^2(pi/8) ~ 0.1464

def h(q):
    if q<=0 or q>=1: return 0.0
    return -q*np.log2(q) - (1-q)*np.log2(1-q)

Q_tol  = brentq(lambda q: h(q)-0.5, 0.001, 0.499)
margin = e_BE - Q_tol        # ~ 3.64%

# Mutual information per qubit
I_BE_full    = 1 - h(e_BE)              # full Breidbart p=1
I_BE_partial = 0.75 * (1 - h(e_BE))    # partial Breidbart p=3/4

print("="*60)
print("  Final Analysis")
print("="*60)
print(f"  e_BE          = {e_BE:.6f}  ({e_BE*100:.2f}%)")
print(f"  Q_tol         = {Q_tol:.6f}  ({Q_tol*100:.2f}%)")
print(f"  margin        = {margin*100:.4f}%")
print(f"  I_BE full     = {I_BE_full:.4f} bits/qubit  (p=1)")
print(f"  I_BE partial  = {I_BE_partial:.4f} bits/qubit  (p=3/4)")
print(f"  eps_unf(N=1000) = p_BE^1000 = 10^({1000*np.log10(p_BE):.1f})")
# ══════════════════════════════════════════════════════════════
# Delta sweep table (printed to terminal)
# ══════════════════════════════════════════════════════════════
from scipy.stats import norm as _norm

def _n_min_approx(eps, delta):
    z = _norm.ppf(eps)
    return int(np.ceil(z**2 * e_BE*(1-e_BE) / delta**2))

def _n_min_exact(eps, delta):
    from scipy.stats import binom as _binom
    for n in range(100, 5_000_001, 100):
        k = int((e_BE-delta)*n)
        if k < 0: return n
        if _binom.cdf(k, n, e_BE) <= eps:
            return n
    return np.nan

named_deltas = [
    ("margin = e_BE - Q_tol",    margin),
    ("1.5 x margin",             margin*1.5),
    ("2 x margin",               margin*2.0),
    ("delta = 5%",               0.05),
    ("delta = e_BE/2",           e_BE/2),
    ("delta = 0.75*e_BE",        e_BE*0.75),
    ("delta = 0.9*e_BE",         e_BE*0.9),
    ("delta = e_BE  (ideal IQT)", e_BE*0.999),
]
eps_list = [1e-6, 1e-9, 1e-20]

print()
print("="*90)
print("  Delta sweep: delta = e_BE - gamma_err")
print("="*90)
print(f"  {'delta':>8}  {'gamma_err':>10}  {'label':28}", end="")
for eps in eps_list:
    print(f"  n_min(1e{int(np.log10(eps))})", end="")
print()
print("  "+"-"*90)
for label, delta in named_deltas:
    ge = e_BE - delta
    print(f"  {delta*100:>7.3f}%  {ge*100:>9.3f}%  {label:28}", end="")
    for eps in eps_list:
        ne = _n_min_exact(eps, delta)
        na = _n_min_approx(eps, delta)
        if np.isnan(ne):
            print(f"  {'N/A':>14}", end="")
        else:
            print(f"  {ne:>6,}/{na:>6,}", end="")
    print()
print("  Format: exact/approx")
print()



# ══════════════════════════════════════════════════════════════
# Figure
# ══════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 6))
fig.patch.set_facecolor('#0f1117')
DARK='#0f1117'; PANEL='#1a1d27'
BLUE='#4da6ff'; CORAL='#ff6b6b'; GREEN='#4dffb4'
AMBER='#ffc94d'; GRAY='#8899aa'; WHITE='#e8eaf6'; PURPLE='#b48eff'
TEAL='#4dffee'

def sa(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=GRAY, labelsize=9)
    ax.spines[:].set_color('#2a2d3a')
    ax.set_title(title, color=WHITE, fontsize=10, pad=8)
    for l in ax.get_xticklabels()+ax.get_yticklabels():
        l.set_color(GRAY)
    ax.xaxis.label.set_color(GRAY)
    ax.yaxis.label.set_color(GRAY)

gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.40)

# ── Panel A: P_forge vs delta ──────────────────────────────────
ax = fig.add_subplot(gs[0])
sa(ax, r'A · $P_{forge}(\delta, n)$ vs $\delta$'
       r'   [$\gamma_{err} = e_{BE} - \delta$]')

delta_arr = np.linspace(1e-4, e_BE*0.999, 400)
n_values  = [1000, 2000, 5000, 10000, 50000, N_Z_PAPER]
colors_n  = [CORAL, AMBER, GREEN, TEAL, BLUE, PURPLE]

for n_val, col in zip(n_values, colors_n):
    pf = np.array([binom.cdf(int((e_BE-d)*n_val), n_val, e_BE)
                   for d in delta_arr])
    log_pf = np.log10(np.maximum(pf, 1e-310))
    lbl = f'$n={n_val//1000}$k' if n_val < N_Z_PAPER else f'$n_Z$={n_val:,} (paper)'
    ax.plot(delta_arr*100, log_pf, color=col, lw=2, label=lbl)

# key vertical lines
ax.axvline(margin*100, color=GREEN, lw=2, ls='--',
           label=f'margin={margin*100:.2f}%\n'
                 r'($\gamma_{err}=Q_{tol}$)')
ax.axvline(e_BE*100,   color=WHITE, lw=1.5, ls=':',
           label=r'$\delta=e_{BE}$  ($\gamma_{err}=0$, ideal)')

# security references
for val, lbl in [(-9,'$10^{-9}$'), (-20,'$10^{-20}$'), (-50,'$10^{-50}$')]:
    ax.axhline(val, color=GRAY, lw=0.7, ls=':', alpha=0.5)
    ax.text(e_BE*100, val+0.5, lbl, color=GRAY, fontsize=8)

ax.set_xlabel(r'$\delta = e_{BE} - \gamma_{err}$ (%)')
ax.set_ylabel(r'$\log_{10}(P_{forge})$')
ax.set_xlim(0, e_BE*100)
ax.legend(fontsize=7.5, facecolor=DARK, labelcolor=WHITE,
          edgecolor='#2a2d3a')

# ── Panel B: eps_unf = p_BE^N vs N ────────────────────────────
ax2 = fig.add_subplot(gs[1])
sa(ax2, r'B · $\varepsilon_{unf} = p_{BE}^N$ vs $N$   [eq. D1, paper bound]')
 
N_cont   = np.logspace(1, 5, 400)
log_eu   = N_cont * np.log10(p_BE)
ax2.plot(N_cont, log_eu, color=PURPLE, lw=3,
         label=r'$\varepsilon_{unf}=p_{BE}^N$')
 
# shade region
ax2.fill_between(N_cont, log_eu, -400,
                 alpha=0.08, color=PURPLE)
 
# mark key N values
for N_val, col, lbl in [
    (100,       CORAL,  '$N=100$'),
    (1000,      AMBER,  '$N=1000$'),
    (10000,     GREEN,  '$N=10^4$'),
    (N_Z_PAPER, BLUE,   f'Paper $n_Z$={N_Z_PAPER:,}'),
]:
    log_val = N_val * np.log10(p_BE)
    ax2.scatter([N_val], [log_val], color=col, s=80, zorder=6)
    ax2.annotate(f'{lbl}\n$10^{{{log_val:.0f}}}$',
             xy=(N_val, log_val),
             xytext=(N_val*1.4, log_val+15),
             color=col, fontsize=7.5,
             arrowprops=dict(arrowstyle='->', color=col, lw=1),
             bbox=dict(boxstyle='round,pad=0.4', facecolor=DARK,
                       edgecolor='none', alpha=0.9))
 
ax2.set_xscale('log')
ax2.set_xlabel('Token size $N$ (log scale)')
ax2.set_ylabel(r'$\log_{10}(\varepsilon_{unf})$')
ax2.set_xlim(N_cont[0], N_cont[-1])
ax2.legend(fontsize=8, facecolor=DARK, labelcolor=WHITE,
           edgecolor='#2a2d3a')

# ── Panel C: I_BE + eps_unf vs N ─────────────────────────────
ax3 = fig.add_subplot(gs[2])
sa(ax3, r'C · $I_{BE}$ (Eve knows more) vs $\varepsilon_{unf}$ (cannot forge)')

N_mi = np.linspace(100, 20000, 500)
I_full_total    = I_BE_full    * N_mi / 2
I_partial_total = I_BE_partial * N_mi / 2
log_eu          = N_mi * np.log10(p_BE)

# Left axis: total I_BE
ax3.plot(N_mi, I_full_total,    color=CORAL, lw=2.5,
         label=f'$I_{{BE}}$ full ($p=1$): {I_BE_full:.3f} bits/qubit')
ax3.plot(N_mi, I_partial_total, color=BLUE,  lw=2.5, ls='--',
         label=f'$I_{{BE}}$ partial ($p=3/4$): {I_BE_partial:.3f} bits/qubit')

ax3.fill_between(N_mi, I_partial_total, I_full_total,
                 alpha=0.15, color=AMBER, label='extra info: full vs partial')
ax3.fill_between(N_mi, 0, I_partial_total, alpha=0.06, color=BLUE)
ax3.set_ylabel('Total $I_{BE}$ (bits) — Eve\'s information $\\uparrow$',
               color=CORAL)
ax3.tick_params(axis='y', labelcolor=CORAL)

# Right axis: eps_unf as tick labels only — no line
ax3r = ax3.twinx()
ax3r.set_facecolor(PANEL)
ax3r.spines[:].set_color('#2a2d3a')
ax3r.set_ylim(0, abs(log_eu[-1]))
ax3r.set_yticks([0, 200, 400, 600, 800, 1000, 1200, 1400])
ax3r.set_yticklabels(
    ['$10^{0}$', '$10^{-200}$', '$10^{-400}$',
     '$10^{-600}$', '$10^{-800}$', '$10^{-1000}$',
     '$10^{-1200}$', '$10^{-1400}$'],
    color=PURPLE, fontsize=8)
    
ax3r.tick_params(axis='y', labelcolor=PURPLE)
ax3r.set_ylabel(r'$\varepsilon_{unf}=p_{BE}^N$ — forging bound $\downarrow$',
                color=PURPLE)

# annotate the contrast at N=5000
N_ann = 5000
I_f   = I_BE_full * N_ann/2
eu_a  = N_ann * np.log10(p_BE)
ax3.annotate(
    f'$N={N_ann:,}$\n'
    f'$I_{{BE}}^{{full}}={I_f:.0f}$ bits\n'
    f'$\\varepsilon_{{unf}}=10^{{{eu_a:.0f}}}$',
    xy=(N_ann, I_f), xytext=(N_ann+2500, I_f-600),
    color=WHITE, fontsize=8,
    arrowprops=dict(arrowstyle='->', color=WHITE, lw=1.2),
    bbox=dict(boxstyle='round,pad=0.3', facecolor=PANEL,
              edgecolor=GRAY, alpha=0.9))

ax3.set_xlim(0, 20000)
ax3.set_xticks([0, 5000, 10000, 15000, 20000])
ax3.set_xticklabels(['0', '5k', '10k', '15k', '20k'])
ax3.set_xlabel('Token size $N$ = pulses sent = detections')
ax3.set_ylim(0)

l1, n1 = ax3.get_legend_handles_labels()
l2, n2 = ax3r.get_legend_handles_labels()
ax3.legend(l1+l2, n1+n2, fontsize=7.5, facecolor=DARK,
           labelcolor=WHITE, edgecolor='#2a2d3a', loc='upper left')

fig.suptitle(
    r'QuantumIQT — Final Analysis  '
    r'[single photon, no loss, full Breidbart $p=1$ / partial $p=3/4$]'
    '\n'
    r'$e_{BE}=\sin^2(\pi/8)\approx14.64\%$  $\cdot$  '
    r'$Q_{tol}\approx11\%$  $\cdot$  '
    r'margin$=3.64\%$  $\cdot$  '
    r'$\varepsilon_{unf}=p_{BE}^N$ (eq. D1)',
    color=WHITE, fontsize=10, y=1.02
)
 
plt.savefig('/Users/ruchithareja/Documents/Python Decoy/Quantum Token Python/IQT/final_analysis.png',
            dpi=150, bbox_inches='tight', facecolor=DARK)
print("\nFigure saved.")
plt.show()
plt.close()