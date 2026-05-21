refseq2fasta <- function(
  phyloseq_obj,
  out_dir,
  phyloseq_list = FALSE,
  extra_id = NULL,
  out_ext = ".fa"
) {
  if (!phyloseq_list) {
    # # Names for files
    nms <- deparse(substitute(phyloseq_obj))
    nms <- gsub("_physeq", "", nms)

    seqs <- phyloseq::refseq(phyloseq_obj)

    # build a safe output path
    fname <- paste0(nms, extra_id, out_ext)
    fname <- gsub("[[:space:]]+", "_", fname)
    fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

    outfile <- file.path(out_dir, fname)

    # Getting Biostring to talk to phylotools
    new_seqs <- Biostrings::as.data.frame(seqs) %>%
      tibble::rownames_to_column(., var = "seq.name") %>%
      rename(seq.text = x)

    tryCatch(
      phylotools::dat2fasta(new_seqs, outfile = outfile),
      error = function(e) {
        message("Failed to write ", outfile, ": ", conditionMessage(e))
      }
    )
  } else if (phyloseq_list && !is.list(phyloseq_list)) {
    message(
      "Not a list of phyloseq objects: class('",
      class(phyloseq_obj),
      "')"
    )
  } else if (is.list(phyloseq_list)) {
    # Names for files
    nms <- names(phyloseq_list)
    nms <- gsub("_physeq", "", nms)

    # Refseqs to a FASTA named by the list name
    purrr::walk2(
      phyloseq_list,
      nms,
      ~ {
        seqs <- phyloseq::refseq(.x)

        # build a safe output path
        fname <- paste0(.y, extra_id, out_ext)
        fname <- gsub("[[:space:]]+", "_", fname)
        fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

        outfile <- file.path(out_dir, fname)

        # Getting Biostring to talk to phylotools
        new_seqs <- Biostrings::as.data.frame(seqs) %>%
          tibble::rownames_to_column(., var = "seq.name") %>%
          rename(seq.text = x)

        tryCatch(
          phylotools::dat2fasta(new_seqs, outfile = outfile),
          error = function(e) {
            message("Failed to write ", outfile, ": ", conditionMessage(e))
          }
        )
      }
    )
  }
}
