# Shared snapshot scrub used by the print/format snapshot tests.
#
# Masks every decimal number, together with the sign and any padding in front of
# it (so a sign flip of a near-zero value cannot change the masked line). Print
# snapshots recorded under local_reproducible_output() therefore pin down layout,
# section order, headers, and wording, but not computed decimals, whose digits can
# drift across BLAS implementations and platforms. Integer output (factor counts,
# df, N, residual counts) is kept verbatim.
scrub_num <- function(lines) gsub("\\s*-?(\\d+)?\\.\\d+", " <num>", lines)

# Like scrub_num, but also masks integer percentages. The EFA_AVERAGE print reports
# error/convergence/Heywood/admissibility rates as `round(mean(...) * 100)` over many inner
# solutions; those percentages can flip across BLAS implementations, so the snapshot pins the
# wording but not the rate. (The deterministic "N EFAs" grid size is left verbatim.)
scrub_num_pct <- function(lines) scrub_num(gsub("[0-9]+%", "<pct>", lines))
