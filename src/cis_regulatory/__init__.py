"""Q4.3 cis-regulatory analysis (Lynch/Wagner transposable-element hypothesis).

The Q4.3 question is whether the human decidual *cis*-regulatory landscape
is disproportionately built from ancient transposable elements (TEs) —
specifically the MER20 and MER41 families that Lynch et al. (2011, 2015)
implicated in wiring the progesterone/cAMP decidual program. The substrate
is GSE61793 (Lynch lab), which ships already-processed hg19 MACS peak
calls for three assays:

- ``h3k27ac`` — active-enhancer peaks (the primary TE-rewiring substrate),
- ``h3k4me3`` — promoter peaks,
- ``dnasei`` — DNaseI open-chromatin union.

Modules:

- :mod:`src.cis_regulatory.peaks` — load BED peak calls into tidy,
  bioframe-compatible DataFrames and summarise them.
- :mod:`src.cis_regulatory.te_overlap` — overlap peaks with the UCSC
  RepeatMasker hg19 annotation and quantify the TE-derived fraction,
  flagging MER20 / MER41.

All heavy IO (downloading RepeatMasker, reading BEDs) is kept in thin
loader functions so the analytic functions stay pure and testable.
"""
