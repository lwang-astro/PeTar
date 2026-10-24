# T9 micro-TDE reference data (author-provided)

- `micro_TDE_events.csv`, `micro_TDE_snapshots.csv`: reference micro-TDE sample
  from S. Rastello (private communication, MODEST26 follow-up, 2026-10;
  PeTar_TDE_SR runs, Z = 0.02). Ten BH-star events in five configuration
  classes (single-single, wide binary, binary-single, triple, 2+2).
- `micro_TDE_snapshots.csv` supplies the last pre-event snapshot state
  (positions/velocities in pc and pc/Myr, masses in Msun, BSE type) used by
  the T9 event-replay runs; `micro_TDE_events.csv` supplies the
  author-measured trigger states (a, e, r_peri, beta, t_to_peri) used as
  comparison values.
- Column documentation is in the original `README_micro_TDE.txt`
  (summarised in the T9 section of `test/validation/README.md`).
