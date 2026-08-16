# One emitted report item: the line matching `start` plus the continuation lines the
# console-width packer indented under it. Used by the width tests of every printer that packs
# a long line, so they all agree on what counts as one item (the indent is .efa_wrap_chunks()'s
# `exdent`, and a blank line or any unindented line ends the item).
wrapped_item <- function(lines, start) {
  i <- grep(start, lines)[1L]
  j <- i
  while (j < length(lines) && startsWith(lines[j + 1L], "  ")) j <- j + 1L
  lines[i:j]
}

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

# Like scrub_num, but also masks integer percentages and the effective number of solutions
# each fit index is averaged over. The efa_average() print reports
# error/convergence/Heywood/admissibility rates as `round(mean(...) * 100)` over many inner
# solutions, and the Model Fit block counts the solutions that reached a usable fit index;
# both are driven by the same per-solution convergence and Heywood outcomes, which can flip
# across BLAS implementations. The snapshot therefore pins the wording but not those counts.
# (The deterministic grid size - the "N EFAs" sentence and the "of N" denominators - is left
# verbatim.)
scrub_num_pct <- function(lines) {
  lines <- gsub("[0-9]+%", "<pct>", lines)
  lines <- gsub("averaged over [0-9]+ of", "averaged over <n> of", lines)
  scrub_num(lines)
}
