###############################################################################
##### wrapper to handle FASTA (possibly containing >1 sequences) as input #####
##############################################################################

is.enriched.wrapper = function(input.fasta, target.kmer, output.csv,
                               significance = 0.05, multiple.testing = NULL){
  # input:
  # - input.fasta: fasta file containing (possibly multiple) sequences; case-insensitive
  # - target.kmer: a string of k-mer of interest; case-insensitive
  # - output.csv: filename of output csv file
  # - significance: significance level; between 0 and 1
  # - multiple.testing: if NULL, no correction for multiple testing; otherwise, one of 
  #   the possible arguments for 'method' in p.adjust()
  # output:
  # a csv file containing the following columns for each sequence
  # - kmer.hits: number of times the k-mer of interest appears in sequence
  # - kmer.p: empirical p-value for enrichment of k-mer of interest
  # - kmer.p.adj: if applicable, adjusted p-value for enrichment of k-mer of interest
  # - count.q1/mean/median/q3/max: summary statistics of counts of all k-mers present
  # - significnt: 0/1; 1 if kmer.p(.adj) < significance
  
  ##### read in fasta file
  fasta = read.delim(input.fasta, header=F, as.is=T)
  
  ##### extract sequences and store in a list
  # fasta sequences are preceded by a single-line description starting with >
  # get row numbers of these descriptions 
  fasta.desc.idx = which(grepl(pattern=">", fasta[, ])==T)
  # actual sequences start at the next line after these descriptions
  seq.start.idx = fasta.desc.idx + 1
  # actual sequences end at the line immediately before these descriptions 
  # take care of edge case (first and last sequence)
  seq.end.idx = c(fasta.desc.idx[-1]-1, nrow(fasta))
  
  n.seq = length(fasta.desc.idx)
  seqs = vector('list', n.seq)
  for (i in 1:n.seq){
    seqs[[i]] = paste(fasta[seq.start.idx[i]:seq.end.idx[i], ], collapse="")
  }
  
  ##### assess k-mer enrichment for each sequence
  enrich = lapply(seqs, is.enriched, kmer = target.kmer)
  enrich = do.call(rbind, enrich)
  
  ##### if indicated, correct for multiple testing
  if (!is.null(multiple.testing)){
    enrich = cbind(enrich, kmer.p.adj = p.adjust(p = enrich[, 'kmer.p'], 
                                                 method = multiple.testing))
  }
  
  ##### significant?
  if (!is.null(multiple.testing)){
    # if correcting for multiple testing, base significance on adjusted p-values
    enrich = cbind(enrich, significant = enrich[, 'kmer.p.adj'] < significance)
  } else {
    # if not, base significance on (unadjusted) p-values
    enrich = cbind(enrich, significant = enrich[, 'kmer.p'] < significance)
  }
  
  ##### export
  write.table(enrich, file=output.csv, quote=F, sep=",", row.names=F)
}


######################################################################
##### compute empirical p-value as a measure of k-mer enrichment #####
######################################################################

is.enriched = function(input.seq, kmer){
  # input:
  # - input.seq: a string of input sequence; case-insentitive
  # - k-mer: a string of k-mer; case-insensitive
  # output:
  # a vector containing 
  # - kmer.hits: number of times the k-mer appears in input.seq
  # - kmer.p: empirical p-value for k-mer enrichment
  # - count.q1/mean/median/q3/max: summary statistics of counts of all k-mers present
  
  kmer = tolower(kmer)
  k = nchar(kmer)
  
  seq.leng = nchar(input.seq)
  input.seq = tolower(input.seq)

  ### initiate count vector
  count.vec = c()
  
  ### go through input sequence, count all the different k-mers present
  for (i in 1:(seq.leng-k+1)){
    cur.kmer = substr(input.seq, start=i, stop=i+k-1)
    if ( is.null(count.vec[cur.kmer]) || is.na(count.vec[cur.kmer]) ) {
      # if count.vec is null (when i=1) or cur.kmer not yet in count.vec
      count.vec = c(count.vec, 1)
      names(count.vec)[length(count.vec)] = cur.kmer
    } else {
      # if cur.kmer already in count.vec
      count.vec[cur.kmer] = count.vec[cur.kmer] + 1
    }
  }
  
  # each kmer combination must be unique
  stopifnot( length(unique(names(count.vec))) == length(count.vec) )
  # no 0 or NA allowed
  stopifnot( sum( is.na(count.vec) | count.vec==0 )==0 )
  
  ### empirical p-value
  p.val = sum(count.vec[kmer] > count.vec) / length(count.vec)
  
  return(c(kmer.hits = as.numeric(count.vec[kmer]), 
           kmer.p = p.val,
           count.q1 = as.numeric(quantile(count.vec, 0.25)),
           count.mean = mean(count.vec),
           count.median = median(count.vec),
           count.q3 = as.numeric(quantile(count.vec, 0.75)),
           count.max = max(count.vec)))
}