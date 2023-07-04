#'Plot convergence graphs and seqlogo
#'@name plot_data
#'@param new_alpha PWM model.
#'@export
plot_data = function(new_alpha) {
  convergence = new_alpha[[2]][-c(1)]
  changes = new_alpha[[3]][-c(1:2)]
  convergence = (convergence - min(convergence)) / (max(convergence) - min(convergence))
  changes = (changes - min(changes)) / (max(changes) - min(changes))
  changes = c(0, changes)
  par(mfrow = c(1,1), mar = c(5,4,4,6))
  plot(as.double(convergence), type='l', lwd = 1, xlab = "Iterations", ylab = "")
  points(as.double(changes), type='l', lty = 2, col=2)
  legend("right", y = 0.75, legend = c("IC", "C"), fill = c(1, 2))
  rownames(new_alpha[[1]]) = c("A", "C", "G", "T")
  show(universalmotif::view_motifs(new_alpha[[1]]))
}

#'@title conservation_level
#'@description computes de the conservation level from consensus string.
#'@param consensus String of the consensus multiple local alignment.
#'@return The conservation level of the consensus string.
#'@export
conservation_level = function(consensus) {
  k = nchar(consensus)
  clevel = 0
  for (i in 1:k) {
    n = Biostrings::substr(x = consensus, start = i, stop = i)
    if (n %in% c("A", "C", "G", "T")) {
      clevel = clevel + 1
    }
    
    else if (n %in% c("R", "Y", "S", "W", "K", "M")) {
      clevel = clevel + 0.5
    }
    
    else if (n %in% c("B", "D", "H", "V")) {
      clevel = clevel + 0.25
    }
  }
  
  clevel / k
}

#'@title extract_jaspar_metadata
#'@description extract JASPAR 2020 or 2022 metadata of all motifs.
#'@param version This is 2020 or 2022.
#'@return The metadata of all motifs in the JASPAR repository.
#'@export
extract_jaspar_metadata = function(version) {
  
  if (version == "2022") matrices = TFBSTools::getMatrixSet(x = JASPAR2022, opts = list())
  else if (version == "2020") matrices = TFBSTools::getMatrixSet(x = JASPAR2020, opts = list())
  else {
    show("Version can be 2022 or 2021")
    return (-1)
  }
  
  metadata = lapply(matrices, function(x) {
    motif = universalmotif::create_motif(x@profileMatrix)
    new_line = data.frame (
      ID = x@ID,
      K = dim(x@profileMatrix)[2],
      NAME = x@name,
      CLASS = x@matrixClass,
      STRAND = x@strand,
      FAMILY = paste(x@tags$family, collapse = ""),
      MEDLINE = x@tags$medline,
      TYPE = paste(x@tags$type, collapse = ""),
      SPECIES = x@tags$species,
      CONSENSUS = motif@consensus,
      IC = motif@icscore,
      N = sum(x@profileMatrix[,1]),
      conservation_level = conservation_level(motif@consensus)
    )
  })
  
  df_jaspar = data.frame(do.call(rbind, metadata), stringsAsFactors = F)
  df_jaspar
}

#'@title write_jaspar_metadata
#'@description Write file JASPAR 2020 or 2022 metadata of all motifs.
#'@param metadata Data Frame with JASPAR metadata.
#'@return Write file meta_jaspar.txt in the current directory.
#'@export
write_jaspar_metadata = function(metadata) {
  write.table(metadata, "meta_jaspar.txt", sep = "\t")
}

#'@title read_jaspar_metadata
#'@description Read file JASPAR 2020 or 2022 metadata of all motifs.
#'@param file Pat with JASPAR metadata file.
#'@return Data Frame with JASPAR metadata.
#'@export
read_jaspar_metadata = function(file) {
  read.table(file)
}

#'@title generate_OOPS
#'@description Generate synthetic OOPS dataset.
#'@param clevel Conservation level.
#'@param k Size of motif to be inserted.
#'@param n Number of sequences.
#'@param t Size of each sequence.
#'@param pos Position of the motif in the sequences.
#'@return Synthetic OOPS dataset.
#'@export
generate_OOPS = function(clevel, k, n = 20000, t = 105, pos = 50) {
  consensus = sample(x = c("A", "C", "G", "T"), size = k, replace = T)
  fasta = c()
  m = t - k + 1
  kmers = c()
  for (i in 1:n) {
    bin_vector = extraDistr::rbern(n = k, prob = clevel)
    kmer = consensus
    kmer[bin_vector == 0] = sample(x = sample(x = c("A", "C", "G", "T"), size = sum(bin_vector == 0), replace = T))
    seqi = sample(x = c("A", "C", "G", "T"), size = t, replace = T)
    seqi[pos:(pos + k - 1)] = kmer
    fasta[i] = paste0(seqi, collapse = "")
    kmers = append(kmers, paste0(kmer, collapse = ""))
  }
    
  alpha = kmers2alpha(kmers)
  list("fasta" = fasta, "alpha" = alpha) 
}

#'@title generate_ZOOPS
#'@description Generate synthetic ZOOPS dataset.
#'@param clevel Conservation level.
#'@param k Size of motif to be inserted.
#'@param n Number of sequences.
#'@param t Size of each sequence.
#'@param pos Position of the motif in the sequences.
#'@param noise Level of noise.
#'@return Synthetic ZOOPS dataset.
#'@export
generate_ZOOPS = function(clevel, k, n = 20000, t = 105, pos = 50, noise = 0.05) {
  consensus = sample(x = c("A", "C", "G", "T"), size = k, replace = T)
  fasta = c()
  m = t - k + 1
  kmers = c()
  for (i in 1:n) {
    bin_vector = extraDistr::rbern(n = k, prob = clevel)
    kmer = consensus
    kmer[bin_vector == 0] = sample(x = sample(x = c("A", "C", "G", "T"), size = sum(bin_vector == 0), replace = T))
    seqi = sample(x = c("A", "C", "G", "T"), size = t, replace = T)
    if (runif(1) > noise) {
      seqi[pos:(pos + k - 1)] = kmer
      kmers = append(kmers, paste0(kmer, collapse = ""))
    } 
    fasta[i] = paste0(seqi, collapse = "")
  }
  alpha = kmers2alpha(kmers)
  list("fasta" = fasta, "alpha" = alpha)
}

#'@title generate_ANR
#'@description Generate synthetic ANR dataset.
#'@param clevel Conservation level.
#'@param k Size of motif to be inserted.
#'@param n Number of sequences.
#'@param t Size of each sequence.
#'@param w Probability of the kmer be a motif.
#'@return Synthetic ANR dataset.
#'@export
generate_ANR = function(clevel, k, n = 20000, t = 105, w = 0.1) {
  consensus = sample(x = c("A", "C", "G", "T"), size = k, replace = T)
  fasta = c()
  m = t - k + 1
  kmers = c()
  for (i in 1:n) {
    seqi = sample(x = c("A", "C", "G", "T"), size = t, replace = T)
    for (pos in 1:m) {
      if (extraDistr::rbern(1, p = w) == 1) {
        bin_vector = extraDistr::rbern(n = k, prob = clevel)
        kmer = consensus
        kmer[bin_vector == 0] = sample(x = sample(x = c("A", "C", "G", "T"), size = sum(bin_vector == 0), replace = T))
        seqi[pos:(pos + k - 1)] = kmer
        kmers = append(kmers, paste0(kmer, collapse = ""))
      }
    }
    fasta[i] = paste0(seqi, collapse = "")
  }
  
  alpha = kmers2alpha(kmers)
  list("fasta" = fasta, "alpha" = alpha)
}

#'@title data_generate
#'@description Generate synthetic dataset.
#'@param clevel Conservation level.
#'@param k Size of motif to be inserted.
#'@param n Number of sequences.
#'@param t Size of each sequence.
#'@param noise Level of noise.
#'@param w Probability of the kmer be a motif.
#'@return Synthetic ZOOPS dataset.
#'@export
data_generate = function(clevel, k, n = 1000, t = 105, noise = 0.05, w = 0.1) {
  consensus = sample(x = c("A", "C", "G", "T"), size = k, replace = T)
  fasta = c()
  m = t - k + 1
  kmers = c()
  for (i in 1:n) {
    seqi = sample(x = c("A", "C", "G", "T"), size = t, replace = T)
    
    if (runif(1) > noise) {
      for (pos in 1:m) {
        if (runif(1) < w) {
          bin_vector = extraDistr::rbern(n = k, prob = clevel)
          kmer = consensus
          kmer[bin_vector == 0] = sample(x = sample(x = c("A", "C", "G", "T"), size = sum(bin_vector == 0), replace = T))
          kmers = append(kmers, paste0(kmer, collapse = ""))
          seqi[pos:(pos + k - 1)] = kmer
        }
      }
    } 
    
    fasta[i] = paste0(seqi, collapse = "")
    
  }
  
  list("fasta" = fasta, "alpha" = kmers2alpha(kmers))
}
