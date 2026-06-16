# Shared snapshot scrub used by the print/format snapshot tests.
#
# Masks every decimal number, together with the sign and any padding in front of
# it (so a sign flip of a near-zero value cannot change the masked line), and
# collapses the value-driven inter-column spacing of factor-label matrix headers
# (lines of the form "F1    F2    F3"). A matrix header carries no decimal for the
# number mask to absorb, so its column alignment would otherwise be pinned verbatim
# and shift whenever a body value changes width. Print snapshots recorded under
# local_reproducible_output() therefore pin down layout, section order, headers
# (labels and their order), and wording, but not computed decimals or the column
# widths they drive, both of which can drift across BLAS implementations, platforms,
# and rotation-package versions. Integer output (factor counts, df, N, residual
# counts) is kept verbatim.
scrub_num <- function(lines) {
  lines <- gsub("\\s*-?(\\d+)?\\.\\d+", " <num>", lines)
  hdr <- grepl("^\\s*F[0-9]+( +F[0-9]+)+ *$", lines)
  lines[hdr] <- gsub("(F[0-9]+) +", "\\1 ", lines[hdr])
  lines
}

# Like scrub_num, but also masks integer percentages. The EFA_AVERAGE print reports
# error/convergence/Heywood/admissibility rates as `round(mean(...) * 100)` over many inner
# solutions; those percentages can flip across BLAS implementations, so the snapshot pins the
# wording but not the rate. (The deterministic "N EFAs" grid size is left verbatim.)
scrub_num_pct <- function(lines) scrub_num(gsub("[0-9]+%", "<pct>", lines))
