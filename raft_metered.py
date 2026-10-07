"""
raft_metered.py — RAFT polymerization with metered CTA addition (method of moments)

Kinetic model of RAFT polymerization with metered chain-transfer-agent (CTA)
feeding. Predicts Mn and dispersity (D) evolution under varied CTA feed
profiles, to guide experimental design of metered-addition polymerization.

Model: ideal degenerative-transfer limit (fast RAFT equilibrium).
State: M, I, lam0 (total radicals), mu0 (dormant chains), CTA,
       T1/T2 = total 1st/2nd moments (active + dormant), D0/D1/D2 = dead chains.

Key idea: chain transfer conserves the TOTAL moments (active + dormant). In the
fast-exchange limit every chain grows at the same average rate kp*M*fa, where
fa = lam0/(lam0+mu0) is the fraction of time a chain spends active. That is
what makes RAFT "living": dormant chains keep growing via reactivation.

Assumptions:
  * RAFT exchange is fast vs propagation -> active and dormant chains share the
    same chain-length distribution at all times (ideal living limit).
  * Termination lumped into a single kt
    (d[rad]/dt = 2*f*kd*[I] - kt*[rad]^2); dead-chain moments use the matching
    lumped form. Teaching-grade model: trust trends, not absolute D.
  * Rate constants are illustrative (styrene-like, 70 C, toluene); replace with
    literature values for a specific system.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- parameters
# Illustrative values for a VBPX-like styrenic monomer at 70 C.
PARAMS = dict(
    kd=3.2e-5,      # s^-1, AIBN decomposition at 70 C
    f=0.6,          # initiator efficiency
    kp=480.0,       # L/mol/s, propagation (styrene-like)
    kt=1.0e8,       # L/mol/s, termination (lumped; see docstring)
    ktr=5.0e5,      # L/mol/s, radical addition to CTA
    M0=1.0,         # mol/L, initial monomer
    I0=15e-3,       # mol/L, AIBN ([M]:[I] = 200:3)
    CTA_total=15e-3,  # mol/L, total CTA ([M]:[CTA] = 200:3)
    M_monomer=299.0,  # g/mol, VBPX ~ C21H17NO (approx)
    t_end=72 * 3600,  # s, 72 h
)

# ---------------------------------------------------------------- feed profile
def feed_rate(t, mode, CTA_total, t_feed):
    """CTA feed rate (mol/L/s): 'batch' = all at t=0; 'metered' = constant
    feed of CTA_total spread over t_feed seconds."""
    if mode == "batch":
        return 0.0  # all CTA present at t=0 (set in y0)
    elif mode == "metered":
        return CTA_total / t_feed if t <= t_feed else 0.0
    else:
        raise ValueError(f"unknown feed mode: {mode}")


def odes(t, y, p, mode, t_feed):
    M, I, lam0, mu0, CTA, T1, T2, D0, D1, D2 = y
    kd, f, kp, kt, ktr = p["kd"], p["f"], p["kp"], p["kt"], p["ktr"]

    T0 = lam0 + mu0
    fa = lam0 / T0 if T0 > 0 else 0.0  # fast-exchange: fraction of time active
    lam1 = T1 * fa
    lam2 = T2 * fa

    dM = -kp * M * lam0
    dI = -kd * I
    dlam0 = 2 * f * kd * I - kt * lam0**2
    dmu0 = ktr * lam0 * CTA          # pre-equilibrium capping builds dormant pool
    dCTA = -ktr * lam0 * CTA + feed_rate(t, mode, p["CTA_total"], t_feed)
    # Total moments: transfer conserves T1/T2; every chain grows at kp*M*fa.
    dT1 = kp * M * lam0 - kt * lam0 * lam1
    dT2 = kp * M * fa * (2 * T1 + T0) - kt * lam0 * lam2
    dD0 = kt * lam0**2               # dead chains (lumped: 1 per event)
    dD1 = kt * lam0 * lam1
    dD2 = kt * lam0 * lam2
    return [dM, dI, dlam0, dmu0, dCTA, dT1, dT2, dD0, dD1, dD2]


def run(mode, p, t_feed=None):
    """Run one scenario; return dict of time series."""
    t_feed = t_feed or p["t_end"]
    CTA_0 = p["CTA_total"] if mode == "batch" else 0.0
    y0 = [p["M0"], p["I0"], 0, 0, CTA_0, 0, 0, 0, 0, 0]
    sol = solve_ivp(odes, (0, p["t_end"]), y0, args=(p, mode, t_feed),
                    method="BDF", rtol=1e-8, atol=1e-14,
                    t_eval=np.linspace(0, p["t_end"], 400))
    M, I, lam0, mu0, CTA, T1, T2, D0, D1, D2 = sol.y
    Mm = p["M_monomer"]
    conv = 1 - M / p["M0"]
    N_tot, W1_tot, W2_tot = (lam0 + mu0) + D0, T1 + D1, T2 + D2
    # np.where evaluates both branches; silence the harmless 0/0 at t=0.
    with np.errstate(divide="ignore", invalid="ignore"):
        Mn = np.where(N_tot > 0, Mm * W1_tot / N_tot, np.nan)
        Mw = np.where(W1_tot > 0, Mm * W2_tot / W1_tot, np.nan)
    return dict(t=sol.t, conv=conv, Mn=Mn, Mw=Mw, D=Mw / Mn, CTA=CTA)


def main():
    p = PARAMS

    # ---- validation: FRP limit and [CTA]/[I] series (batch) ----
    print("== validation ==")
    frp = run("batch", {**p, "ktr": 0.0, "CTA_total": 0.0})
    print("FRP (ktr=0, no CTA):      final D = %.2f (expect ~1.5-2.0)" % frp["D"][-1])
    for ratio in [1, 3, 10]:
        s = run("batch", {**p, "CTA_total": p["I0"] * ratio})
        print("[CTA]/[I]=%-2d (batch):      final D = %.2f, Mn = %.0f g/mol"
              % (ratio, s["D"][-1], s["Mn"][-1]))

    # ---- feed-profile comparison at fixed total CTA ----
    print("== feed-profile comparison (total CTA fixed) ==")
    scenarios = {
        "batch": run("batch", p),
        "metered_24h": run("metered", p, t_feed=24 * 3600),
        "metered_6h": run("metered", p, t_feed=6 * 3600),
    }
    for name, s in scenarios.items():
        print("%-12s: final D = %.2f, Mn = %.0f g/mol"
              % (name, s["D"][-1], s["Mn"][-1]))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for name, s in scenarios.items():
        ax[0].plot(s["conv"], s["Mn"] / 1000, label=name)
        ax[1].plot(s["conv"], s["D"], label=name)
    ax[0].set(xlabel="conversion", ylabel="Mn (kg/mol)")
    ax[1].set(xlabel="conversion", ylabel="dispersity D")
    for a in ax:
        a.legend()
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("mn_d_vs_conversion.png", dpi=150)
    print("saved mn_d_vs_conversion.png")


if __name__ == "__main__":
    main()
