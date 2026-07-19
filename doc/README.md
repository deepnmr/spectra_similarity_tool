# Source papers

The Shift-Tolerant Correlation Coefficient method (this work) is documented in full, with
equations and proofs, in the Angewandte manuscript [`LCC_angewandte.md`](LCC_angewandte.md)
and its Supporting Information [`LCC_angewandte_SI.md`](LCC_angewandte_SI.md).

The documents also report the experimental Local-Contrast STCC option added in
`hsqc_lcc.py`: a square-root intensity transform followed by fixed-$3\sigma$
local background subtraction. STCC (`--method lcc`) remains the default pending
an improvement in held-out threshold error or rejection false-positive rate.

The methods in this repository come from the following papers. Only the
open-access one is stored here; the two Elsevier papers are kept locally but
excluded from git (see `.gitignore`) to respect their copyright — use the DOIs
to obtain them.

| file | paper | access |
| --- | --- | --- |
| `Pierens2012_JCheminformatics_HSQC_NN_DGA.pdf` | Pierens et al. 2012 | Open Access (CC-BY), included |
| `Bodis2009_Talanta_2D_HSQC_bins.pdf` | Bodis et al. 2009 | Elsevier, local only |
| `Castillo2013_Chemometrics_tree.pdf` | Castillo et al. 2013 | Elsevier, local only |

## Citations

1. L. Bodis, A. Ross, E. Pretsch, *A novel spectra similarity measure*,
   Chemometrics and Intelligent Laboratory Systems 85 (2007) 1-8.
   doi:10.1016/j.chemolab.2006.03.006 — the 1D bin method.
2. L. Bodis, A. Ross, J. Bodis, E. Pretsch, *Automatic compatibility tests of
   HSQC NMR spectra with proposed structures of chemical compounds*, Talanta 79
   (2009) 1379-1386. doi:10.1016/j.talanta.2009.06.017 — the 2D bin extension.
3. A.M. Castillo, L. Uribe, L. Patiny, J. Wist, *Fast and shift-insensitive
   similarity comparisons of NMR using a tree-representation of spectra*,
   Chemometrics and Intelligent Laboratory Systems 127 (2013) 1-6.
   doi:10.1016/j.chemolab.2013.05.009 — the quad-tree method.
4. G.K. Pierens, S. Brossi, Z. Yang, D.C. Reutens, V. Vegh, *HSQC spectral based
   similarity matching of compounds using nearest neighbours and a fast discrete
   genetic algorithm*, Journal of Cheminformatics 4 (2012) 25.
   doi:10.1186/1758-2946-4-25 — the nearest-neighbour peak method.
