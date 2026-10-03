# Developability metrics: provenance and evidence

This documents the literature basis for every developability / interface-quality
metric `binderqc` reports, and — honestly — which ones are well-supported versus
which are weak heuristics you should report but not gate on. References were
checked against primary sources (not from memory); where a citation is marked
*(unverified)* the metric is standard but a specific primary reference still needs
to be pinned before it is quoted as fact.

## Empirical anchor

The de facto ground truth for "what makes a biologic developable" is the
clinical-stage landscape:

- **Jain et al. 2017, PNAS** — *Biophysical properties of the clinical-stage
  antibody landscape*. 12 biophysical assays measured on 137 clinical-stage mAbs;
  the empirical basis for which in-silico flags track real liabilities and where
  the acceptable limits sit. [doi:10.1073/pnas.1616408114](https://doi.org/10.1073/pnas.1616408114)
- **Raybould et al. 2019, PNAS** — *Five computational developability guidelines
  for therapeutic antibody profiling* (the Therapeutic Antibody Profiler, TAP):
  CDR length, patches of surface hydrophobicity (PSH), patches of positive/negative
  charge (PPC/PNC), structural Fv charge symmetry. Source for the charge/hydrophobic
  patch metrics. [doi:10.1073/pnas.1810576116](https://doi.org/10.1073/pnas.1810576116)

> Caveat that applies to all of these: TAP/Jain thresholds were calibrated on
> antibody Fv regions. `binderqc` targets small de novo binders (SAP is noted as
> calibrated on 52–65-residue binders), so treat absolute thresholds as
> provisional and prefer *relative ranking within a campaign* over hard cutoffs.

## Metric-by-metric

| Metric (column) | Measures | Primary reference | Evidence | Use |
|---|---|---|---|---|
| `sap_score`, `sap_total`, `sap_per_res` | Spatial aggregation propensity — dynamically exposed, spatially-clustered hydrophobics | Chennamsetty et al., PNAS 2009 [doi:10.1073/pnas.0904191106](https://doi.org/10.1073/pnas.0904191106) | **Strong** — validated, widely used aggregation predictor | Keep; primary aggregation flag |
| `a3d_score`, `a3d_total_positive` | Structure-aware aggregation propensity (Aggrescan3D) | Kuriata et al., NAR 2019 (A3D 2.0) [doi:10.1093/nar/gkz321](https://doi.org/10.1093/nar/gkz321); orig. Zambrano et al. 2015 | **Strong–Moderate** — established; complements SAP | Keep; cross-check with SAP |
| `charge_patch_pos`, `charge_patch_neg` | Surface charge patches (TAP PPC/PNC style) | Raybould et al., PNAS 2019 | **Strong** — TAP validated on clinical mAbs | Keep; thresholds provisional for minibinders |
| `paratope_hydrophobicity`, `paratope_charge` | Interface patch composition (TAP PSH-style, paratope-scoped) | Raybould et al., PNAS 2019 | **Moderate** — right concept, interface-scoped variant not independently validated here | Keep as ranking signal |
| `buns_interface`, `interface_polar_satisfied_frac` | Buried unsatisfied polar atoms at the interface | Stranges & Kuhlman, Protein Sci 2013 [doi:10.1002/pro.2187](https://doi.org/10.1002/pro.2187); penalized in Cao et al., Nature 2022 [doi:10.1038/s41586-022-04654-9](https://doi.org/10.1038/s41586-022-04654-9) | **Strong** for interface quality — buried-unsat separates successful from failed interface designs | Keep; the key polar-interface-quality signal (new) |
| `n_hbonds` | Interface N/O–N/O contacts ≤3.5 Å (H-bond proxy, no H/angle) | geometric proxy; H-bond design context Boyken et al., Science 2016 [doi:10.1126/science.aad8865](https://doi.org/10.1126/science.aad8865) | **Moderate** — over-counts (no donor/acceptor or angle check) | Report; pair with `buns_interface` |
| `n_salt_bridges` | Cross-interface cation–anion pairs ≤4 Å | geometric | **Moderate** — geometric proxy | Report |
| `gravy` | Mean Kyte–Doolittle hydrophobicity | Kyte & Doolittle, J Mol Biol 1982 [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0) | **Moderate** — crude 1-D hydrophobicity; 3-D SAP/A3D are better for aggregation/solubility | Report; defer to SAP/A3D |
| `pi`, `mw`, `ext_coeff_280` | ProtParam isoelectric point, mass, extinction | Gasteiger et al. 2005, ProtParam (ExPASy) | **Moderate** (pI: solubility dips near pI) / descriptive (MW, ε) | Report |
| `instability_index` | Guruprasad dipeptide-composition instability | Guruprasad et al., Protein Eng 1990 [doi:10.1093/protein/4.2.155](https://doi.org/10.1093/protein/4.2.155) | **Weak** — a 1990 dipeptide heuristic; poorly predictive of real stability (e.g. applicability study, Int. J. 2019) | **Report only — never gate** |
| `sequence_liabilities`: deamidation (NG/NS/NT), Asp isomerization (D-[G/S/T/D/H]) | Chemical degradation motifs | Sydow et al., PLoS ONE 2014 [doi:10.1371/journal.pone.0100736](https://doi.org/10.1371/journal.pone.0100736) | **Strong motif, but sequence-only over-predicts** — actual risk needs local flexibility + C-flank size | Keep; **refinement: gate on exposure/flexibility** (binderqc has structure) |
| `sequence_liabilities`: unpaired Cys | Free-thiol (scrambling/aggregation/conjugation) risk | developability liability literature (reviewed in Jain 2017 context) *(unverified specific ref)* | **Moderate** — odd-count heuristic; even counts can still expose a free thiol | Keep; refine with SG SASA |
| `sequence_liabilities`: exposed Met/Trp oxidation | Solvent-exposed oxidation-prone residues (SASA-gated) | oxidation liability literature *(unverified specific ref)* | **Moderate–Strong** — exposure-gating is the correct refinement | Keep; elevate severity for radioligand (RLT) candidates |
| `epitope_glyco_occluded`, `epitope_glyco_sites` | Target N-X-[S/T] sequons (X≠P), SASA-aware, near interface | canonical sequon rule *(specific ref unverified)* | **Strong motif / Moderate for actual occupancy** — structure-aware check is the right refinement | Keep |
| `interface_packing` | Heavy-atom contacts per 100 Å² BSA (CMS proxy) | proxy for Rosetta contact-molecular-surface | **Weak** — unvalidated proxy (already labelled in code) | Report, don't gate |
| `approach_angle`, `epitope_planarity`, grippability | Pose / designability heuristics (IARA grippability validated in-house: EDB 9/9, TNC 0/8) | design heuristics | **Moderate** — designability, not developability | Report |

## Verdict — help vs. report-only

**Help (well-grounded, keep as quality signals / gates):** `sap_score`, `a3d_score`,
`charge_patch_*`, **`buns_interface` / `interface_polar_satisfied_frac`** (new),
deamidation & Asp-isomerization motifs *(with structural gating)*, unpaired-Cys,
exposed Met/Trp oxidation, structure-aware glyco occlusion.

**Report, do not gate (weak or crude):** `instability_index` (weak 1990 heuristic —
the clearest "don't gate" case), `gravy` (defer to SAP/A3D for aggregation),
`interface_packing` (unvalidated proxy), `pi` (context-dependent).

**Highest-value refinement:** the deamidation/isomerization flags are currently
sequence-motif-only, which over-predicts. Sydow 2014 shows the real risk is set by
local backbone flexibility and the C-terminal flanking residue — and `binderqc`
already has the structure, so gating these flags on exposure/flexibility is the
single biggest accuracy gain available.

## Not currently in binderqc (candidates, require model dependencies)

These are well-supported but each adds a model dependency beyond the
biotite/numpy/pandas core, so they belong as optional extras rather than core:
NetSolP solubility (Thumuluri et al. 2022), WALTZ amyloid stretches
(Maurer-Stroh et al. 2010), MHC-II immunogenicity (IEDB), and an inverse-folder
interface NLL (Frame2Seq / ProteinMPNN) as an orthogonal "does the sequence fit
the backbone" check.
