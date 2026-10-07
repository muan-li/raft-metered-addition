# RAFT Polymerization with Metered CTA Addition — Kinetic Simulator

A method-of-moments ODE model of RAFT polymerization with metered
chain-transfer-agent (CTA) feeding. Predicts number-average molecular weight
(Mn) and dispersity (D) evolution under varied CTA feed profiles, to guide
experimental design of metered-addition polymerization.

## Model

Ideal degenerative-transfer limit (fast RAFT equilibrium). State variables:
monomer, initiator, total radicals, dormant chains, CTA, total 1st/2nd moments
(active + dormant), dead-chain moments. Chain transfer conserves the total
moments; in the fast-exchange limit every chain (active or dormant) grows at
the same average rate — that is what makes RAFT "living".

Assumptions:
- RAFT exchange is fast vs propagation: active and dormant chains share the
  same chain-length distribution at all times (ideal living limit).
- Termination lumped into a single kt (convention:
  d[rad]/dt = 2*f*kd*[I] - kt*[rad]^2).
- Rate constants are illustrative (styrene-like, 70 C, toluene); replace with
  literature values for a specific system.

## Run

```bash
pip install -r requirements.txt
python raft_metered.py
```

Prints validation numbers, then compares three feed profiles
(`batch`, `metered_24h`, `metered_6h`) and saves `mn_d_vs_conversion.png`:
Mn and dispersity vs conversion. Feed profiles are implemented in
`feed_rate()`; scenarios are configured in `main()`.

## Validation

- FRP limit (ktr = 0, no CTA): D -> 2.14 at high conversion (uncontrolled) ✓
- Batch RAFT, [CTA]/[I] = 1 / 3 / 10: D = 1.55 / 1.29 / 1.26,
  Mn = 8098 / 4242 / 1590 g/mol ✓ (textbook RAFT behavior: more CTA ->
  more chains, lower Mn, narrower distribution)

## Results

**Feed-profile comparison** (total CTA fixed at 15 mM, [M]:[CTA]:[I] = 200:3:3):

| feed profile | final D | final Mn (g/mol) |
|---|---|---|
| batch (all CTA at t=0) | 1.55 | 8098 |
| metered over 6 h | 1.82 | 8098 |
| metered over 24 h | 2.50 | 8098 |

- Mn is identical in all cases: set by total CTA, not by feed profile.
- Slower feeding broadens the distribution: chains initiated early grow
  uncontrolled before CTA arrives. The transient shows the metered-addition
  fingerprint — Mn initially *decreases* (newly created chains dilute the
  average) and D goes through a hump (long uncontrolled chains mixed with
  newly born short chains); the hump's turning point marks the end of feeding.

**ktr sensitivity** (batch, [CTA]/[I] = 1):

| ktr (L/mol/s) | minimum D | final D |
|---|---|---|
| 1e4 | 1.34 | 1.55 |
| 1e5 | 1.22 | 1.55 |
| 5e5 | 1.22 | 1.55 |
| 1e6 | 1.22 | 1.55 |

- Threshold behavior: 1e4 -> 1e5 deepens control (Dmin 1.34 -> 1.22); above
  ~1e5 the system is in the fast-exchange limit and faster transfer changes
  nothing.
- Final D is independent of ktr (1.55 in all cases): the endpoint is set by the
  dead-chain fraction, i.e. by [CTA]/[I], not by transfer rate.

**Takeaway for experimental design** — three knobs, three jobs:
- `[CTA]/[I]` sets the endpoint dispersity;
- `ktr` sets how fast/deep control is established (saturates above ~1e5);
- the feed profile shapes the distribution (metered addition trades some
  narrowness for control over shape and low instantaneous [CTA]).

## Reference

Method of moments for RAFT polymerization: standard textbook treatment of
living-polymerization moments plus degenerative chain transfer.
