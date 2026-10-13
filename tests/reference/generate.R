# Generate reference distances with ape::dist.dna for the Rust test `ape_reference`.
#
#   Rscript tests/reference/generate.R      (from the repository root)
#
# Writes tests/reference/ape_distances.tsv. The output is deterministic for a given ape version;
# regenerate it only when adding cases or models, and commit the script and the table together.
#
# Each row is one aligned pair, scored on its own (a two-sequence alignment), with pairwise
# deletion: a site counts only if both sequences hold A/C/G/T there. That is the library's rule
# (any case, U = T). ape has no U and uses global base frequencies for F81, T92 and TN93, so:
#   - U is read as T before the pair reaches ape;
#   - base.freq is the composition of the compared sites only, which is what the library uses.
# ape's T92 sets h = 2 * GC * (1 - GC) from the pooled G+C content, which is Tamura's (1992)
# h = GC1 + GC2 - 2 * GC1 * GC2 only when both sequences have the same G+C content at the compared
# sites. T92 is therefore NA for every other pair.
#
# When a base is missing from the compared sites, ape's F81 and TN93 are NaN (0/0 terms) while the
# library drops the zero-weight terms; likewise T92 when the G+C content is 0 or 1. Those values are
# NA here: src/tests.rs covers the degenerate compositions.
#
# ape also computes log arguments from rounded proportions, so a pair exactly on the JC69 or K80
# saturation boundary can come out finite (K80 with 1 transition and 1 transversion in 3 sites:
# 1 - 2/3 - 1/3 = 5.6e-17, a distance of 19). The library treats such arguments as zero; the
# table records Inf, the exact answer, for those pairs (found from the integer counts).
#
# Saturated models give NaN in ape (the log of a non-positive number) and Infinity in the library;
# no comparable sites gives NaN in both.
suppressPackageStartupMessages(library(ape))
set.seed(20261013)

BASES <- c("A", "C", "G", "T")
TRANSITION <- c(A = "G", G = "A", C = "T", T = "C")
TRANSVERSIONS <- list(A = c("C", "T"), G = c("C", "T"), C = c("A", "G"), T = c("A", "G"))
NOISE <- c("-", "-", "-", "N", "N", "?", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V")
COMPOSITIONS <- list(
  uniform = c(0.25, 0.25, 0.25, 0.25),
  gc_rich = c(0.15, 0.35, 0.35, 0.15),
  at_rich = c(0.35, 0.15, 0.15, 0.35),
  skewed  = c(0.40, 0.10, 0.30, 0.20)
)
LENGTHS <- c(5, 31, 63, 64, 65, 127, 128, 129, 200, 500, 1000)

random_seq <- function(len, comp) sample(BASES, len, replace = TRUE, prob = comp)

# Substitute transitions at rate `ts` and transversions at rate `tv`.
mutate <- function(s, ts, tv) {
  u <- runif(length(s))
  for (k in seq_along(s)) {
    if (u[k] < ts) s[k] <- TRANSITION[[s[k]]]
    else if (u[k] < ts + tv) s[k] <- sample(TRANSVERSIONS[[s[k]]], 1)
  }
  s
}

# Swap `k` disjoint pairs of sites: the composition is unchanged, so T92 is comparable with ape.
shuffle_sites <- function(s, k) {
  idx <- sample(length(s), 2 * min(k, length(s) %/% 2))
  half <- length(idx) / 2
  a <- idx[seq_len(half)]
  b <- idx[half + seq_len(half)]
  tmp <- s[a]; s[a] <- s[b]; s[b] <- tmp
  s
}

# Replace a fraction of sites with gaps, N and ambiguity codes.
add_noise <- function(s, rate) {
  hit <- runif(length(s)) < rate
  s[hit] <- sample(NOISE, sum(hit), replace = TRUE)
  s
}

# Spell some A/C/G/T in lower case and some T as U; neither changes the distance.
respell <- function(s, rate) {
  low <- runif(length(s)) < rate
  s[low] <- tolower(s[low])
  u <- s %in% c("T", "t") & runif(length(s)) < rate
  s[u] <- ifelse(s[u] == "T", "U", "u")
  s
}

cases <- list()
add_case <- function(kind, s1, s2) {
  cases[[length(cases) + 1]] <<- list(kind = kind, s1 = paste(s1, collapse = ""), s2 = paste(s2, collapse = ""))
}

# Random pairs over the range of divergences, Ts/Tv ratios, compositions and lengths
for (comp in names(COMPOSITIONS)) {
  for (len in LENGTHS) {
    for (div in c(0.01, 0.1, 0.3)) {
      kappa <- sample(c(0.5, 2, 8), 1)
      s1 <- random_seq(len, COMPOSITIONS[[comp]])
      add_case(paste("mutated", comp, sep = "/"), s1, mutate(s1, div * kappa / (kappa + 1), div / (kappa + 1)))
    }
  }
}

# The same with gaps, N and ambiguity codes in both sequences, and mixed spelling
for (comp in names(COMPOSITIONS)) {
  for (len in c(65, 200, 1000)) {
    for (rate in c(0.02, 0.2)) {
      s1 <- random_seq(len, COMPOSITIONS[[comp]])
      s2 <- mutate(s1, 0.08, 0.04)
      add_case(paste("noisy", comp, sep = "/"), respell(add_noise(s1, rate), 0.1), respell(add_noise(s2, rate), 0.1))
    }
  }
}

# Equal composition at the compared sites, so ape's T92 is Tamura's; gaps at shared sites only
for (comp in names(COMPOSITIONS)) {
  for (len in c(64, 129, 500)) {
    for (k in c(1, 10, 60)) {
      s1 <- random_seq(len, COMPOSITIONS[[comp]])
      s2 <- shuffle_sites(s1, k)
      gaps <- runif(len) < 0.05
      s1[gaps] <- "-"
      s2[gaps] <- sample(c("-", "N"), sum(gaps), replace = TRUE)
      add_case(paste("shuffled", comp, sep = "/"), s1, s2)
    }
  }
}

# Unrelated sequences: most models saturate
for (len in c(20, 64, 300)) {
  for (comp in names(COMPOSITIONS)) {
    add_case(paste("unrelated", comp, sep = "/"), random_seq(len, COMPOSITIONS[[comp]]), random_seq(len, COMPOSITIONS[[comp]]))
  }
}

# Identical sequences, and pairs without a single comparable site
s <- random_seq(150, COMPOSITIONS$uniform)
add_case("identical", s, s)
add_case("identical", add_noise(s, 0.3), s)
add_case("no_sites", strsplit("----NNNN", "")[[1]], strsplit("ACGTACGT", "")[[1]])
add_case("no_sites", strsplit("RYKMSWRY", "")[[1]], strsplit("ACGTACGT", "")[[1]])

MODELS <- c("N", "raw", "TS", "TV", "JC69", "K80", "F81", "T92", "TN93", "logdet")

score <- function(s1, s2) {
  a <- toupper(strsplit(s1, "")[[1]])
  b <- toupper(strsplit(s2, "")[[1]])
  a[a == "U"] <- "T"
  b[b == "U"] <- "T"
  keep <- a %in% BASES & b %in% BASES
  x <- as.DNAbin(rbind(a, b))
  bf <- if (any(keep)) base.freq(as.DNAbin(rbind(a[keep], b[keep]))) else rep(0.25, 4)
  gc1 <- sum(a[keep] %in% c("C", "G"))
  gc2 <- sum(b[keep] %in% c("C", "G"))
  gc_usable <- gc1 == gc2 && gc1 > 0 && gc1 < sum(keep)
  s <- sum(keep)
  ts <- sum(keep & a != b & (a %in% c("A", "G")) == (b %in% c("A", "G")))
  tv <- sum(keep & a != b) - ts
  on_boundary <- c(JC69 = s > 0 && 3 * s == 4 * (ts + tv),
                   K80 = s > 0 && (s == 2 * ts + tv || s == 2 * tv))
  d <- sapply(MODELS, function(m) {
    if (m %in% names(on_boundary) && on_boundary[[m]]) return(Inf)
    if (m == "T92" && !gc_usable) return(NA)
    if (m %in% c("F81", "TN93") && any(bf == 0)) return(NA)
    suppressWarnings(as.numeric(dist.dna(x, model = m, pairwise.deletion = TRUE, base.freq = bf)))
  })
  c(sites = sum(keep), d)
}

fmt <- function(v) ifelse(is.na(v) & !is.nan(v), "NA", sprintf("%.17g", v))

rows <- lapply(seq_along(cases), function(i) {
  cs <- cases[[i]]
  v <- score(cs$s1, cs$s2)
  c(id = i, kind = cs$kind, fmt(v), seq1 = cs$s1, seq2 = cs$s2)
})
table <- do.call(rbind, rows)
colnames(table) <- c("id", "kind", "sites", MODELS, "seq1", "seq2")

out <- file.path("tests", "reference", "ape_distances.tsv")
con <- file(out, "w")
writeLines(c(
  sprintf("# Generated by tests/reference/generate.R with ape %s (%s). Do not edit.",
          packageVersion("ape"), R.version.string),
  paste(colnames(table), collapse = "\t"),
  apply(table, 1, paste, collapse = "\t")
), con)
close(con)
cat(sprintf("Wrote %d pairs to %s\n", nrow(table), out))
