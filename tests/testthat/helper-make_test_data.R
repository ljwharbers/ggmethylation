# Helpers shared by the unit tests. testthat sources helper-*.R before the
# test files.

# Minimal methylation_data constructor (used by test-build_read_panel.R,
# test-binary_call_mode.R, test-smooth_ci_ribbon.R and test-delta_track.R).
make_test_data = function(reads_df, sites_df,
                          region_start = 1000L, region_end = 2000L,
                          cigar_features = data.frame(
                            read_name  = character(0),
                            type       = character(0),
                            ref_start  = integer(0),
                            ref_end    = integer(0),
                            length     = integer(0),
                            stringsAsFactors = FALSE
                          )) {
  gr = GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = region_start, end = region_end)
  )
  structure(
    list(
      reads          = reads_df,
      sites          = sites_df,
      region         = gr,
      mod_code       = "m",
      group_tag      = NULL,
      cigar_features = cigar_features
    ),
    class = "methylation_data"
  )
}

# Paths to the synthetic fixture modBAM and VCF (data-raw/make_test_fixture.R).
fixture_bam = function() testthat::test_path("fixtures", "synthetic.bam")
fixture_vcf = function() testthat::test_path("fixtures", "synthetic.vcf.gz")
