# How every number is computed

Each column binderqc reports, with its definition, its source, and where it
departs from that source. Everything is biotite + numpy on coordinates — no MD,
no PAE, no trained model, no external binary.

Throughout: **SASA** is solvent-accessible surface area from
`biotite.structure.sasa`; **relative SASA (relSASA)** is residue SASA divided by
the Tien et al. 2013 theoretical maximum for that residue type; **Cβ** is the
side-chain centroid proxy (Cα for glycine).

## The two exposure conventions

This is the one convention worth understanding before reading any number, because
mixing the two silently is the failure mode.

| convention | SASA computed on | correct for |
|---|---|---|
| `complex` | binder **with the target present** | anything that happens to the *bound* molecule |
| `free` | binder **alone** | anything that happens to the binder in solution, before or between binding events |

A hydrophobic residue buried against the target reads as safe under `complex` and
as a liability under `free`. Both are right, for different questions.

**Aggregation and patch columns (`sap_*`, `a3d_*`, `charge_patch_*`) are always
computed on the free binder**, because that is how their source definitions are
written (Chennamsetty 2009, Kuriata 2019, Raybould 2019 are all surface
properties of the molecule itself). The oxidation and deamidation checks take
`--degradation-context` (`complex` by default).

## Interface and pose

**`binder_bsa`** — area the binder buries on binding: its SASA alone minus its
SASA in the complex. A tiny interface is the clearest tell of a junk binder.
Validated against PDBePISA on two-chain complexes (`tests/pisa_correctness.py`).

**`interface_packing`** — interface heavy-atom contacts within 4 Å per 100 Å² of
buried area. A packing-density *proxy* for Rosetta's contact molecular surface;
denser contacts per buried area means better packing. Not the validated Rosetta
metric.

**`n_hbonds`** — interface N/O···N/O pairs within 3.5 Å. A geometric proxy: with
no hydrogens we cannot check donor/acceptor identity or angle, so this
over-counts.

**`n_salt_bridges`** — distinct cross-chain residue pairs with a cationic
side-chain N (Arg/Lys/His) within 4 Å of an anionic side-chain O (Asp/Glu).

**`buried_unsat_polar`** — binder polar heavy atoms that are (a) in the interface
zone, (b) buried in the complex, yet (c) have no polar heavy atom in a *different*
residue close enough to act as an H-bond partner. Burying a polar group costs the
H-bond it made to water, so an unsatisfied buried polar is destabilising. A known
separator of successful from failed interface designs (Stranges & Kuhlman, *Protein
Sci* 2013) and penalised in de novo binder design (Cao et al., *Nature* 2022).
Same no-hydrogens caveat as `n_hbonds`.

**`approach_angle`** — angle (0–90°) between the binder's long axis and the
paratope→epitope axis. Near 0° the binder points end-on into the target; near 90°
it lies across the surface. STCRpy measures a docking angle for TCRs against a
fixed MHC frame; this is the same idea without needing a reference frame, so it
works for any binder.

## Grippability

**`epitope_planarity`** — RMSD of epitope points to their best-fit plane. Small
means a flat patch with little to grip; larger means concave or knobby. NaN below
three points.

**`glyco_occluded_epitope`** — N-glycosylation sequons (N-X-[S/T], X≠P) on the
*target* that are both solvent-exposed and at or near the epitope, where an
installed glycan would sterically mask the patch. Two things make this
SASA-aware rather than a sequence grep: exposure is measured on the **target
alone** (glycans are installed before binding, so complex SASA would wrongly hide
the very sites of interest), and sequons are read off consecutive modelled
residue ids, so an unresolved loop cannot fabricate a phantom sequon.

**`grippability_consensus()`** — combines the physical read above with an
*optional* learned score (e.g. IARA hotspot probability), passed in by the caller
rather than invoked here, so binderqc stays dependency-free and never vendors an
externally-licensed model.

## Did it bind where you asked?

**`epitope_contact_fraction`, `epitope_precision`** — every other epitope column is
descriptive: it reports whatever patch the complex happens to show, so a design
that drifted to a second site passes all of them, because the site it landed on
*becomes* the epitope they describe. These two compare contacted residues against
an epitope supplied from outside the structure. Residue-level over target residues
in contact, not contact-weighted — a count of atom pairs would let one deeply
buried residue outvote five shallow ones.

## Tag-site selection

Per terminus: **relative SASA** (buried termini make bad tag sites),
**distance** from the terminal Cα to the nearest paratope Cα, and
**orientation** — the cosine between the direction the chain extends past the
terminus and the direction to the paratope. A tag grows outward along
(terminal Cα − adjacent Cα); if that points at the paratope the tag heads
straight for the target, which is bad even when the terminus is far away.

**`sg_sasa`** — for a cysteine terminus, the SASA of the SG atom itself.
≈0 Å² means buried and unreactive; >~15 Å² means exposed and conjugatable.

## Aggregation and surface patches

All on the **isolated binder**.

**`sap_score` / `sap_total`** — Spatial Aggregation Propensity (Chennamsetty et
al. 2009). Per residue, sum over residues whose Cβ lies within 10 Å of
(relSASA × Black-Mould hydrophobicity, shifted so Gly = 0). `sap_score` is the
peak (hottest patch); `sap_total` is the sum of positive windowed values
(whole-surface load). *Simplifications vs the original:* single-structure SASA
rather than an MD average, and whole-residue relSASA (Tien max-ASA) rather than
side-chain SAA against an Ala-X-Ala reference.

**`a3d_score` / `a3d_total_positive`** — Aggrescan3D (Zambrano et al. 2015),
ported from Aggrescan3D 1.0.2. Per residue, with Cα as centre:

```
agg_sup = 0                                             if RSA% < 10  (buried)
        = a3v[aa] · 0.0599 · exp(0.0521 · min(RSA, 55))  otherwise
agg_fin = agg_sup + Σ over residues within 10 Å (Cα–Cα, excluding self)
                     of  agg_sup_neighbour · 1.2915 · exp(−2.56/10 · dist)
```

Buried central residues keep `agg_fin = 0`. `a3d_score` is the peak;
`a3d_total_positive` is the sum of positive values, which is what the Aggrescan3D
server reports as its "A3D score". The `a3v` scale and every constant are taken
verbatim from Aggrescan3D 1.0.2 and are regression-tested
(`tests/test_constant_tables.py`). The scale was checked against upstream's own
`aggrescan/data/matrices/aggrescan.mat` on 2026-11-20 — **zero mismatches across
all twenty residues**. That file is not vendored here, because Aggrescan3D is
distributed "free for non-commercial users" and binderqc is MIT.

RSA for this column uses the **naccess** reference max-ASA that Aggrescan3D's own
freeSASA call uses (`_A3D_REF_MAX_ASA`), not the Tien 2013 table the rest of the
tool uses. That distinction is load-bearing: Tien's values are ~25% larger, which
moves every partially-exposed residue down a band and made the port run about 30%
hot. *Remaining caveat:* absolute SASA still comes from biotite rather than
freeSASA with naccess radii, so values are close but not identical.

Measured end to end against real Aggrescan3D 1.0.2 on 40 single-chain structures
(`tests/a3d_correctness.py`): Pearson r = 0.9996 for `a3d_score`, 0.9998 for
`a3d_total_positive`, median relative error 0.9%, mean ratio 0.99.

**`charge_patch_pos` / `charge_patch_neg`** — TAP-style surface charge patches
(Raybould et al. 2019). Per residue, the SASA-weighted formal charge summed over
residues whose Cβ is within 10 Å, keeping the positive and negative sides
separate so a basic cluster is not cancelled by a nearby acidic one. Charges are
R/K (+1) and D/E (−1); His is left neutral (~10% protonated at pH 7.4).
Clustered like-charge drives viscosity and self-association, which net charge, pI
and salt-bridge counts all miss. **Thresholds are deliberately not baked in** —
TAP's amber/red cutoffs are calibrated on antibody Fv surfaces and do not
transfer to de novo minibinders, so this is a reported number, not a gate.

**`paratope_hydrophobicity` / `paratope_charge`** — the patch idea applied where
TAP actually applies it: to the binding region rather than the whole surface.
Over binder residues within 10 Å (Cβ) of an interface residue, the SASA-weighted
Black-Mould hydrophobicity (Gly-shifted) and the SASA-weighted net formal charge.
Hydrophobicity is a positive load; charge is signed (>0 basic, <0 acidic).

## Sequence

**`sequence_liabilities`** — motif red flags, with thresholds from Adaptyv's
protein-qc skill (MIT). Things to look at, not automatic rejections.

Deamidation (N-[G/S/T]) and Asp isomerisation (D-[G/S/T/D/H]) are deliberately
**not** sequence-only here: Sydow et al. 2014 found the real risk is set by local
backbone flexibility and exposure as well as the C+1 flank, so a motif hit
over-predicts. binderqc has the structure, so both are **exposure-gated** on the
N/D residue's own relSASA, alongside the oxidation check.

**`gravy`, `net_charge`, `pi`, `lys_arg_pct`, `instability_index`** — standard
ProtParam-style quantities. `net_charge` and `lys_arg_pct` are **whole-binder**,
unlike `paratope_charge` which is interface-scoped. That distinction is load-bearing:
in a large de novo campaign (211 designs, two CROs' SPR/BLI plus BVP-ELISA and an
ovalbumin polyreactivity counter-screen), whole-molecule net charge (ρ = −0.57) and
Lys+Arg content (ρ = 0.50) were the two strongest correlates of nonspecific binding
out of 25 candidate measures.

## Self-association

**`self_association()`** — takes a predicted binder homodimer plus a `score_structure`
row for the same binder on its target, and reports how much of the paratope the
self-interface covers. Interface PAE, ipTM and ipSAE all answer *whether* a
self-interface is predicted, never *where*. The distinction decides what to do:
a binder that self-associates away from the paratope is a formulation problem;
one that self-associates *through* the paratope competes with its own target
binding. Residue numbering must match between the two structures;
`residues_matched` is false when it does not.

## What binderqc deliberately does not compute

**ipSAE and other confidence metrics.** Those need PAE, not coordinates. binderqc
is model-agnostic physical QC: it takes a finished complex from any source and
asks physical questions about it.

## References

- Tien et al. (2013) *PLoS ONE* — maximum accessible surface areas
- Chennamsetty et al. (2009) *PNAS* — Spatial Aggregation Propensity
- Zambrano et al. (2015) *Nucleic Acids Res* — Aggrescan3D; Kuriata et al. (2019) A3D 2.0
- Raybould et al. (2019) *PNAS* — Therapeutic Antibody Profiler patches
- Stranges & Kuhlman (2013) *Protein Sci* — buried unsatisfied polars
- Cao et al. (2022) *Nature* — de novo binder design
- Sydow et al. (2014) *PLoS ONE* — structure-dependent deamidation/isomerisation risk
- Black & Mould (1991) — residue hydrophobicity scale
- Guruprasad et al. (1990) — instability index
