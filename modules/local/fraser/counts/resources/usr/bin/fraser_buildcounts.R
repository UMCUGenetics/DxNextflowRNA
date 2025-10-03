#!/usr/bin/env Rscript

library(data.table)
library(FRASER)
library(argparse)

get_arguments <- function() {
  parser <- ArgumentParser()

  parser$add_argument("--refset_junctions", nargs = "+")
  parser$add_argument("--sample_junctions", nargs = "+")
  parser$add_argument("--prefix")

  return(parser$parse_args())
}


merge_tsvs <- function(files, sample_names = NULL) {
  if (is.null(sample_names)) {
    # standaard sample namen uit bestandsnaam
    sample_names <- sub("\\.tsv\\.gz$", "", basename(files))
  }
  stopifnot(length(files) == length(sample_names))

  dt_list <- lapply(seq_along(files), function(i) {
    f <- files[i]
    s <- sample_names[i]

    dt <- fread(f)

    # alle kolommen behalve "sample" zijn de keys
    key_cols <- setdiff(names(dt), "sample")

    setnames(dt, old = "sample", new = s)
    setkeyv(dt, key_cols)
    dt
  })

  # bepaal merge keys uit de eerste tabel
  key_cols <- setdiff(names(dt_list[[1]]), sample_names[1])

  merged <- Reduce(function(x, y) merge(x, y, by = key_cols, all = TRUE), dt_list)

  # NA → 0
  for (j in (length(key_cols) + 1):ncol(merged)) {
    set(merged, which(is.na(merged[[j]])), j, 0)
  }

  return(merged)
}

# Voeg startID en endID toe
add_junction_ids <- function(junc_df){
  dt <- as.data.table(junc_df)
  setorder(dt, seqnames, start, end, strand)

  # StartIDs
  donors <- unique(dt[, .(seqnames, start, strand)])
  donors[, startID := .I]
  dt <- merge(dt, donors, by=c("seqnames","start","strand"), all.x=TRUE, sort=FALSE)

  # EndIDs
  acceptors <- unique(dt[, .(seqnames, end, strand)])
  acceptors[, endID := .I + nrow(donors)]
  dt <- merge(dt, acceptors, by=c("seqnames","end","strand"), all.x=TRUE, sort=FALSE)

  setcolorder(dt, c("seqnames","start","end","width","strand",
                    setdiff(names(dt), c("seqnames","start","end","width","strand"))))
  return(as.data.frame(dt))
}




# Leid splice sites af
make_splice_sites <- function(junc_df){
  dt <- as.data.table(junc_df)
  sample_cols <- setdiff(names(dt), c("seqnames","start","end","width","strand","startID","endID"))

  # Donors
  donor_sites <- unique(dt[, .(seqnames, pos = start, strand)])
  donor_sites[, `:=`(
    start = pos - 1L,
    end   = pos,
    width = 2L,
    type  = "Donor"
  )]
  donor_sites[, spliceSiteID := .I]
  donor_sites[, pos := NULL]

  donor_counts <- dt[, lapply(.SD, sum), by=.(seqnames, start, strand), .SDcols=sample_cols]
  donors <- merge(donor_sites, donor_counts,
                  by.x=c("seqnames","end","strand"),
                  by.y=c("seqnames","start","strand"), all.x=TRUE)

  # Acceptors
  acceptor_sites <- unique(dt[, .(seqnames, pos = end, strand)])
  acceptor_sites[, `:=`(
    start = pos,
    end   = pos + 1L,
    width = 2L,
    type  = "Acceptor"
  )]
  acceptor_sites[, spliceSiteID := .I + nrow(donors)]
  acceptor_sites[, pos := NULL]

  acceptor_counts <- dt[, lapply(.SD, sum), by=.(seqnames, end, strand), .SDcols=sample_cols]
  acceptors <- merge(acceptor_sites, acceptor_counts,
                     by.x=c("seqnames","start","strand"),
                     by.y=c("seqnames","end","strand"), all.x=TRUE)

  splice_dt <- rbind(donors, acceptors, fill=TRUE)
  setcolorder(splice_dt, c("seqnames","start","end","width","strand",
                           "spliceSiteID","type", sample_cols))
  return(as.data.frame(splice_dt))
}

# Sample table
make_sample_table <- function(junc_df){
  sample_cols <- setdiff(names(junc_df),
    c("seqnames","start","end","width","strand","startID","endID","spliceSiteID","type"))
  data.frame(sampleID = sample_cols)
}


args <- get_arguments()

junction_files <- c(args$refset_junctions, args$sample_junctions)
junction_dt <- merge_tsvs(junction_files)
junction_dt <- add_junction_ids(junction_dt)
splice_dt <- make_splice_sites(junction_dt)
sampleTable <- make_sample_table(junction_dt)


### Resultaten bekijken
cat("Junctions:\n")
print(junction_dt)

cat("\nSplice sites:\n")
print(splice_dt)

cat("\nSample table:\n")
print(sampleTable)


fwrite(junction_dt, paste0(args$prefix, "_junction_counts.tsv.gz"), sep="\t")
fwrite(splice_dt,   paste0(args$prefix, "_spliceSite_counts.tsv.gz"), sep="\t")
fwrite(sampleTable, paste0(args$prefix, "_sampleTable_countTable.tsv"), sep="\t")



# Onderstaande moet naar de algemene fraser module
## library(FRASER)
## library(data.table)

## # Lees tabellen
## sampleTable <- fread("sampleTable_countTable.tsv")
## junctions   <- fread("junction_counts.tsv.gz")
## spliceSites <- fread("spliceSite_counts.tsv.gz")

## # Maak FraserDataSet
## fds <- FraserDataSet(
##   colData    = sampleTable,
##   junctions  = junctions,
##   spliceSites= spliceSites,
##   workingDir = "FRASER_output"
## )
