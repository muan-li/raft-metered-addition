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
python raft_metered_skeleton.py
```

Outputs `mn_d_vs_conversion.png`: Mn and dispersity vs conversion for each
feed profile. Feed profiles are implemented in `feed_rate()`
(`"batch"`, `"metered"`, ...); scenarios are configured in `main()`.

## Validation

- FRP limit (ktr = 0, no CTA): D -> ~2.14 at high conversion (uncontrolled) ✓
- Batch RAFT, [CTA]/[I] = 1 / 3 / 10: D = 1.55 / 1.30 / 1.26,
  Mn = 8098 / 4242 / 1590 g/mol ✓ (textbook RAFT behavior: more CTA ->
  more chains, lower Mn, narrower distribution)

## Results

Work in progress — feed-profile optimization (batch vs metered addition at
fixed total CTA) coming.

## Reference

Method of moments for RAFT polymerization: standard textbook treatment of
living-polymerization moments plus degenerative chain transfer.
