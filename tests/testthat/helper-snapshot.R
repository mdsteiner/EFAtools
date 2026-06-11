# Shared snapshot scrub used by the print/format snapshot tests.
#
# Masks every decimal number, together with the sign and any padding in front of
# it (so a sign flip of a near-zero value cannot change the masked line). Print
# snapshots recorded under local_reproducible_output() therefore pin down layout,
# section order, headers, and wording, but not computed decimals, whose digits can
# drift across BLAS implementations and platforms. Integer output (factor counts,
# df, N, residual counts) is kept verbatim.
scrub_num <- function(lines) gsub("\\s*-?(\\d+)?\\.\\d+", " <num>", lines)
