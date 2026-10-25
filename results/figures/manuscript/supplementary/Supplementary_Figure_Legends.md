# Additional file 1. Supplementary Figures S1-S6

Single-nucleus transcriptomic reanalysis prioritizes NFE2L2 and THRB regulatory programs in
focal cortical dysplasia type II astrocytes.

Each figure is supplied as a TIFF file at submission resolution. Panel letters are lower case
and all figures use the same palette as the main figures. The seven prioritized candidate
transcription factors are NFE2L2, THRB, SATB2, SOX2, SMAD3, SOX10 and CUX1.

**Figure_S1_discovery_object_construction_qc.tiff**

Figure S1. Discovery object construction and quality control. (a) Nuclei retained per source in
the discovery object (523, 426, 911 and 462; blue, internal control; orange, lesion). (b)
Source-level filtering and subsampling trace on a logarithmic scale, from input through quality
control and subsampling to the final object; colour is the condition and line style is the donor
(solid, G120; dashed, G133), and each trace is labelled with its sample at the right. (c) Number
of genes detected per nucleus in each source and (d) mitochondrial read fraction per nucleus in
each source; the dashed lines mark the filtering thresholds (200 and 8,000 genes, 20%
mitochondrial reads) and each violin shows the median as a horizontal line. Filtering additionally
capped each source at 3,000 nuclei. The merged object contained 2,322 nuclei, 36,601 genes and
105 motif-supported regulons, and the AUCell activity matrix matched the object by barcode with
no missing nucleus.

**Figure_S2_full_differential_regulon_landscape.tiff**

Figure S2. Full differential regulon landscape. (a) Row z-scored mean regulon activity of the top
40 differential regulons across the four samples (RdBu, range -1.7 to 1.7, matching Figure 2);
the seven prioritized candidates are labelled in teal. (b) The two ends of the differential
landscape: the 14 most internal-control-associated and the 14 most lesion-associated regulons
ordered by the lesion minus internal-control mean AUCell difference, with the group boundary
marked by a dashed line and the two groups labelled in the panel. (c) Co-activity of the
candidate regulons, the complete symmetric Spearman correlation matrix over all nuclei, with
blue indicating negative and orange positive correlation; CEBPD, JUNB and TCF4 are shown as
non-prioritized context. (d) Prioritisation funnel from 105 motif-supported regulons to 76 with
regulon FDR < 0.05, 30 meeting every numerical criterion, 9 literature-qualified factors and the
7 prioritized candidates.

**Figure_S3_CellOracle_THRB_vector_field.tiff**

Figure S3. CellOracle supplementary perturbation view. (a) Observed mean displacement of each
perturbed network against the randomised control (2 standard deviations); bars are the observed
network, open circles the randomised control and z the paired permutation statistic. A
permutation reference distribution was available for four transcription factors, and BHLHE40
(marked with a dagger) is not one of the prioritized candidates. (b) THRB field drawn as the
transcription-factor-specific residual, that is the local mean displacement of the THRB network
minus the mean displacement of the seven perturbed candidate networks, which is the component not
shared across candidates; it is drawn on its own scale and carries its own scale bar. (c, d) NFE2L2
network-specific field, the local mean observed displacement minus the paired randomised-control
displacement, computed separately in donors G120 (n = 949) and G133 (n = 1,373). Panels c and d
share the vector scale of the NFE2L2 field in Figure 3D. All three fields use the same reference
embedding and the same bounding box, and each carries a scale bar in the lower left corner.

**Figure_S4_supportive_external_evidence.tiff**

Figure S4. Supportive external evidence for all seven prioritized candidates, drawn as bubble
plots in which area carries the value and colour its interpretation. (a) GSE190452, two-group
cohort. Rows are the candidates and columns are the group log2 fold change (lesion minus control),
the cluster direction support, the significant cluster support and the detection fraction in the
expected-higher group. Bubble area is the magnitude of the value normalised within each column and
colour indicates whether the external direction matches the candidate's expected direction (teal,
as expected; orange, not reproduced); the sign of the fold change is therefore not drawn, and the
per-column area normalisation means areas are comparable within a column only. (b) GSE140393,
lesion-only cohort. Detection fraction per nucleus in the five lesion samples; area is the
detection fraction and colour indicates whether the 5% evaluability rule is met (teal, detected;
grey, below 5%). In both panels a light dashed guide line is drawn at every column position, and a
dashed line at the same height separates the two leading axes from the five supporting candidates.
Each panel carries a size key and a colour key. The archived support score of GSE140393 is not
shown because in a single-group cohort it rewards low detection. SATB2 stayed below the 5% rule in
every sample of GSE140393 and is therefore not evaluable in that cohort. The archived support
scores that accompany these two cohorts, which cover the three candidates with archived external
expression support (NFE2L2, THRB and SOX2), are given in Additional file 2, Table ST07; the
remaining candidates in panels a and b are shown from the re-run support tables.

**Figure_S5_functional_enrichment_convergence.tiff**

Figure S5. Functional enrichment and program convergence. (a) Regulon-target set size against the
number of direction-consistent intersection genes with the differential expression program. The
x axis is logarithmic because four candidates have target sets of only 7 to 18 genes, and each
label is joined to its point by a leader line. (b) The same convergence expressed as the
percentage of regulon target genes that are differentially expressed, ordered by value (SOX2
80.0%, THRB 76.5%, SATB2 59.5%, SOX10 38.9%, SMAD3 31.6%, NFE2L2 29.5%, CUX1 28.6%). (c)
Regulon-target set size on a logarithmic axis against the number of significant GO and KEGG terms
after full-family Benjamini-Hochberg correction. (d) GO terms shared by at least two candidates,
with colour giving the number of sharing candidates (2 or 3). Row labels are shortened terms and
the complete GO term names are given in Additional file 2, Table ST08.

**Figure_S6_drug_signature_enrichment.tiff**

Figure S6. Drug-signature enrichment of the candidate programs. (a) Significant DSigDB terms per
query set: 893 for the NFE2L2 regulon targets and 492 for the NFE2L2 core program. Both THRB
query sets returned no significant term, so drug-signature enrichment is driven by the NFE2L2
program alone. (b) Top ten NFE2L2 core-program terms by adjusted p value. (c) mTOR-axis primary
leads and (d) supportive clinical leads; bar length is the best adjusted p value across the
candidate query sets on a log10 scale and bar colour is the curated clinical-support tier, where
tier A is the highest and tier E the lowest level of clinical support. A compound that reached no
significant term is marked by a short tick at zero instead of a bar, which is the case for
temsirolimus and everolimus in panel c. (e) Blood-brain-barrier annotation of the twelve nominated
compounds grouped into high, moderate and low penetration, with bar length the best adjusted p
value and colour the blood-brain-barrier class. The annotation separates compounds that could reach
a central nervous system target from those that would need delivery optimisation.
