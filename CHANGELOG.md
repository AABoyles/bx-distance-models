# Changelog

All notable changes to this project are recorded here. HIV-TRACE's `tn93`
([veg/tn93](https://github.com/veg/tn93), `src/tn93_shared.cc`) is the reference for what the
ambiguity modes mean.

## Unreleased

### Added

- `SequenceSet` ambiguity mode `'gapmm'` (`tn93 -a gapmm`): as `'average'`, and a gap facing
  anything other than a gap is compared as an `N`, so a gap facing a base adds a quarter to each of
  the four (base, X) cells. Gap facing gap is skipped. Every metric supports it, not only TN93, and
  `links()` keeps its early exit. Sites compared as `N` are counted with popcounts over a separate
  bit-plane rather than one by one (see `docs/performance-analysis.md`).
- `?` is read as `N` (any base), as in tn93. Outside `'gapmm'` both are still skipped.
- Reference tests against tn93 v1.1.0 for TN93 in every ambiguity mode
  (`tests/reference/tn93_distances.tsv`, from `tests/reference/tn93/generate.py`), including
  MicrobeTrace's golden alignments.
