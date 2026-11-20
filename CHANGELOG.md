# Changelog

All notable changes to this project are recorded here. HIV-TRACE's `tn93`
([veg/tn93](https://github.com/veg/tn93), `src/tn93_shared.cc`) is the reference for what the
ambiguity modes mean.

## Unreleased

### Changed

- **`N` and `?` are now fourfold ambiguity codes in `'resolve'` and `'average'` modes**, as in
  tn93: `N` facing `A` counts as a match (A-A += 1) when resolving, and adds a quarter to each of
  A-A, A-C, A-G and A-T when averaging. They were skipped before, so distances in these modes
  change for any sequence that contains `N` or `?`: resolving usually lowers them (more matching
  sites), averaging usually raises them. They now agree with HIV-TRACE's `tn93 -a resolve|average`.
  `'skip'` mode and the single-pair functions (`tn93_distance`, ...) are unchanged and still skip
  `N`. Pairs whose sequences hold `N` also take an extra popcount pass in these modes (~40 ns/pair
  in `bench:api`, where every sequence has ~3 N).

### Added

- `SequenceSet` option `maxAmbiguityFraction` for `'resolve'` mode: HIV-TRACE's
  `tn93 -a resolve -g FRACTION`, natively. A pair is resolved unless the sites resolution touched
  are at least that fraction of its compared sites; then it is averaged. As in tn93, a resolved
  site is a base facing a code that contains it, or two ambiguity codes (N included) whether or
  not they share a base. One scan per pair, and `links()` keeps its exact early exit:
  `new SequenceSet(seqs, 'resolve', { maxAmbiguityFraction: 0.015 })`. Unknown options throw.
- TypeScript types `AmbiguityMode` and `SequenceSetOptions`.
- `SequenceSet` ambiguity mode `'gapmm'` (`tn93 -a gapmm`): as `'average'`, and a gap facing
  anything other than a gap is compared as an `N`, so a gap facing a base adds a quarter to each of
  the four (base, X) cells. Gap facing gap is skipped. Every metric supports it, not only TN93, and
  `links()` keeps its early exit. Sites compared as `N` are counted with popcounts over a separate
  bit-plane rather than one by one (see `docs/performance-analysis.md`).
- `?` is read as `N` (any base), as in tn93.
- Reference tests against tn93 v1.1.0 for TN93 in every ambiguity mode
  (`tests/reference/tn93_distances.tsv`, from `tests/reference/tn93/generate.py`), including
  MicrobeTrace's golden alignments.
