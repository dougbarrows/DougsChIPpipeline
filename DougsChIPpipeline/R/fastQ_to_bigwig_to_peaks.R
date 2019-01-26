
#' ChIPseq Pipeline (fastq to bw to peaks)
#'
#' Dougs pipeline that takes a fastq, aligns, creates a bigwig, then calls peaks if desired, with a lot of fun QC along the way
#'
#' @rdname fastQ_to_bigwig_to_peaks
#' @param fastq_path path to the folder wehre the fastqs are
#' @param mapq mapping quality score (mapq) used for filtering BAMs. The default is mapq = 15, which is what ChIPQC uses as its cut off
#' @param buildIndex if TRUE this will build the index for the Rsubread alignment, only need to do this once. 
#' @param sample_sheet_path Path to the sample sheet for peak calling. If this argument is not included, no peaks will be called.
#' @param bw MACS2 argument if calling peaks, default is 300
#' @param q MACS2 argument if calling peaks, default is 0.05
#' @param genome gemone used for alignments, etc
#' @param blacklist_path path to blacklist bed file (include this as data once you know how to do that?)
#' @details NA
#' @return many things
#' @examples
#' message("coming soon!")
#' @export
#' @import ShortRead Rsubread Rsamtools BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg19 BSgenome.Mmusculus.UCSC.mm9 BSgenome.Mmusculus.UCSC.mm10 rtracklayer GenomicAlignments ggplot2 ChIPQC devtools soGGi TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Mmusculus.UCSC.mm9.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene ChIPseeker


fastQ_to_bigwig_to_peaks <- function(fastq_path, mapq = 15, buildIndex = FALSE, sample_sheet_path,  bw = 300, q = 0.05, genome, blacklist_path) {

  ################
  # must have BSgenome index files for Rsubread in the same folder as this markdown
  # if not, then put buildIndex = TRUE and it will get built in the pipeline 
  ##############
  
  ##########################
  #####if calling peaks!!!
  # must open r studio from terminal with command "open -a RStudio" for the system command to inherit the right path for macs2
  # can check you are okay by running command "system("echo $PATH)", which should return the same path as if you were to run "echo $PATH" in the terminal and should include path to macs2 (in my case it is in the anaconda folder)
  # also make sure you are in python2 environment before you open Rstudio from terminal - "source activate python2"
  
  # NOTE: the pipeline will not call peaks if you do not put a sample sheet path in the command
  ###########################
  
  if(missing(buildIndex)) stop("No 'buildIndex' argument specified, indicate whether you want an index built")
  
  if(missing(blacklist_path)) stop("path to blacklist ('blacklist_path') bed file is required")
  
  if(missing(genome)) stop("genome ('genome') specification is required")
  
  if(missing(fastq_path)) stop("path to fastq files ('fastq_path') is required")
  
  if (missing(sample_sheet_path)) {
    print("No sample sheet included, no peaks will be called")
  } 
  
  fastq <- list.files(fastq_path)
  
  ShortRead_QA_path <- paste0(fastq_path, "ShortRead_QA")
  dir.create(ShortRead_QA_path)
  
  BAM_path <- paste0(fastq_path, "BAM_files/")
  dir.create(BAM_path)
  
  ChIPQC_path <- paste0(fastq_path, "ChIP_QC")
  dir.create(ChIPQC_path)
  
  bigwig_path <- paste0(fastq_path, "bigwigs/")
  dir.create(bigwig_path)
  ###########
  
  # Now we loop through file to analyze each one
  
  bam_files <- vector(mode = "character", length = length(fastq))
  
  for (i in seq_along(fastq)) {
    
    ##########
    # do some QC with shortread package
    ##########
    
    print(paste0("ShortRead QA analysis for ", fastq[i], "..."))
    
    # calculate the percent of duplicates in the fastq file
    shortread_fastq <- readFastq(paste0(fastq_path, fastq[i]))
    dupLogical <- srduplicated(shortread_fastq)
    numberOfDups <- table(dupLogical)
    perc_dups <- numberOfDups[2]/(numberOfDups[1] + numberOfDups[2])
    names(perc_dups) <- "Percent_Duplicates"
    out <- capture.output(perc_dups)
    cat( out, file = paste0(ShortRead_QA_path, "/", fastq[i], "_duplication_perc.txt"), sep = "\n", append = TRUE)
    
    # generate a QA report (similar to FastQC)
    myQA <- qa(paste0(fastq_path, fastq[i]))
    myReport <- report(myQA, dest = paste0(ShortRead_QA_path, "/", fastq[i], "_ShortreadQA.html"))
    ######
    
    # filter reads and write out a FastQ of our filtered reads
    
    # this will pass the fastq file in a certain number of reads at a time (can be reduced if computer can't handle this number)
    fqStreamer <- FastqStreamer(paste0(fastq_path, fastq[i]),
                                n=5000000)
    filtered_fastq <- paste0( fastq_path, "filtered_",fastq[i])
    
    # this loop will use the Fastq streamer to loop through the whole file and filter out reads that have more thatn 10 "N" which means the sequencer couldnt identify the base, and those that have a quality score ('alphabetScore') below a certain value (300 is what Tom used in the course)
    TotalReads <- 0
    TotalReadsFilt <- 0
    while (length(fq <- yield(fqStreamer))>0) {
      TotalReads <- TotalReads+length(fq)
      filt1 <- fq[alphabetScore(fq) > 300 ]
      filt2 <- filt1[alphabetFrequency(sread(filt1))[,"N"] < 10]
      TotalReadsFilt <- TotalReadsFilt+length(filt2)
      writeFastq(filt2, filtered_fastq ,mode="a")
    }
    TotalReads
    TotalReadsFilt
    
    
    ########## build the index for Rsubread alignment - only need to do this once
    
    # Here we cycle through the major chromosomes and create a DNAStringSet object from the retrieved sequences.
    if (buildIndex == TRUE) {
      if (i == 1) {
        print("Building index...")
        
        if (genome == "hg19") {
          mainChromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
          mainChrSeq <-
            lapply(mainChromosomes, function(x)
              BSgenome.Hsapiens.UCSC.hg19[[x]])
          names(mainChrSeq) <- mainChromosomes
          mainChrSeqSet <- DNAStringSet(mainChrSeq)
          
          # Now we have a DNAStringSet object we can use the writeXStringSet to create our FASTA file of sequences to align to.
          
          writeXStringSet(mainChrSeqSet,
                          "BSgenome.Hsapiens.UCSC.hg19.fa")
          
          #From Toms course - The Rsubread package offers a faster aligner than the QuasR package (bowtie for R) although the Rsubread package only works on Macs. For alignment with the Rsubread package we must first build our genome index for Rsubread using the buildindex() function.The buildindex() function simply takes the parameters of our desired index name and the FASTA file to build index from.
          
          buildindex("BSgenome.Hsapiens.UCSC.hg19",
                     "BSgenome.Hsapiens.UCSC.hg19.fa")
        }
        
        if (genome == "hg38") {
          mainChromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
          mainChrSeq <-
            lapply(mainChromosomes, function(x)
              BSgenome.Hsapiens.UCSC.hg38[[x]])
          names(mainChrSeq) <- mainChromosomes
          mainChrSeqSet <- DNAStringSet(mainChrSeq)
          
          writeXStringSet(mainChrSeqSet,
                          "BSgenome.Hsapiens.UCSC.hg38.fa")
          
          buildindex("BSgenome.Hsapiens.UCSC.hg38",
                     "BSgenome.Hsapiens.UCSC.hg38.fa")
        }
        
        if (genome == "mm9") {
          mainChromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
          mainChrSeq <-
            lapply(mainChromosomes, function(x)
              BSgenome.Mmusculus.UCSC.mm9[[x]])
          names(mainChrSeq) <- mainChromosomes
          mainChrSeqSet <- DNAStringSet(mainChrSeq)
          
          writeXStringSet(mainChrSeqSet,
                          "BSgenome.Mmusculus.UCSC.mm9.fa")
          
          buildindex("BSgenome.Mmusculus.UCSC.mm9",
                     "BSgenome.Mmusculus.UCSC.mm9.fa")
        }
        
        if (genome == "mm10") {
          mainChromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
          mainChrSeq <-
            lapply(mainChromosomes, function(x)
              BSgenome.Mmusculus.UCSC.mm10[[x]])
          names(mainChrSeq) <- mainChromosomes
          mainChrSeqSet <- DNAStringSet(mainChrSeq)
          
          writeXStringSet(mainChrSeqSet,
                          "BSgenome.Mmusculus.UCSC.mm10.fa")
          
          buildindex("BSgenome.Mmusculus.UCSC.mm10",
                     "BSgenome.Mmusculus.UCSC.mm10.fa")
        }
      }
    }
    
    ###### now to the alignment using Rsubread
    
    # extract the base part of the name of each fastq file so that this can be used as the base for output files generated below
    if (grepl(".fastq", filtered_fastq)) {
      base <- strsplit(fastq[i], ".fastq", fixed = TRUE)[[1]][1]
    } else if (grepl(".fq", filtered_fastq)){
      base <- strsplit(fastq[i], ".fq", fixed = TRUE)[[1]][1]
    } 
    
    #align with Rsubread
    
    bam <- paste0(base, ".BAM")
    
    # here I am using the filtered fastq, do we need to do this given that we filter out low quality alignments? I suppose we don't, but it might make the alignments faster if we are gettign rid of bad reads befroehand
    
    print(paste0("Aligning ", fastq[i], "..."))
    
    if (genome == "hg19") {
      out <- capture.output(align(index = "BSgenome.Hsapiens.UCSC.hg19",
                                  readfile1 = filtered_fastq,
                                  output_file = paste0(BAM_path, bam),
                                  type="dna"))  # in Tom's class he had "phredOffset = 64" in this command too, but Rsubread kept throwing a warning saying it was wrong, so I'll just keep the default for now
      
      cat( out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    }
    
    if (genome == "hg38") {
      out <- capture.output(align(index = "BSgenome.Hsapiens.UCSC.hg38",
                                  readfile1 = filtered_fastq,
                                  output_file = paste0(BAM_path, bam),
                                  type="dna"))  
      
      cat( out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    }
    
    if (genome == "mm9") {
      out <- capture.output(align(index = "BSgenome.Mmusculus.UCSC.mm9",
                                  readfile1 = filtered_fastq,
                                  output_file = paste0(BAM_path, bam),
                                  type="dna")) 
      
      cat( out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    }
    
    if (genome == "mm10") {
      out <- capture.output(align(index = "BSgenome.Mmusculus.UCSC.mm10",
                                  readfile1 = filtered_fastq,
                                  output_file = paste0(BAM_path, bam),
                                  type="dna")) 
      
      cat( out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    }
    
    print(paste0("Sorting ", fastq[i], "..."))
    sorted_bam_pre <- paste0( BAM_path, base, "_sorted")
    sortBam(paste0(BAM_path, bam), sorted_bam_pre) # don't put ".BAM" at end of output file as it adds it anyway
    
    sorted_bam_post <- paste0(BAM_path, base, "_sorted.BAM")
    indexBam(sorted_bam_post)
    
    file.remove(paste0(BAM_path, bam))
    
    out <- capture.output(quickBamFlagSummary(sorted_bam_post))
    cat("BamFlagSummary - Before Filtering", out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    
    # filter out low quality reads
    # mapping quality score (mapq) can be specified, but the default is mapq = 15, which is what ChIPQC uses as its cut off, so seems reasonable
    
    print(paste0("Filtering out low quality reads and duplicates for ", fastq[i], "..."))
    param <- ScanBamParam(mapqFilter = mapq, 
                          flag = scanBamFlag(isDuplicate = FALSE))
    sorted_bam_post_mapq <- paste0(sorted_bam_pre, "nodup_q", as.character(mapq), ".BAM")
    filterBam(file = sorted_bam_post, destination = sorted_bam_post_mapq, param = param)
    
    bam_files[i] <- sorted_bam_post_mapq # this will eventualy be used to get the exact names of the BAM files generated to call peaks
    
    out <- capture.output(quickBamFlagSummary(sorted_bam_post_mapq))
    cat("BamFlagSummary - After Filtering", out, file = paste0(BAM_path, fastq[i], "_BAMlog.txt"), sep = "\n", append = TRUE)
    
    # make figure showing all mapped reads to each chromosome
    mappedReads <- idxstatsBam(sorted_bam_post_mapq)
    TotalMapped <- sum(mappedReads[,"mapped"])
    ggplot(mappedReads,aes(x=seqnames,y=mapped))+
      geom_bar(stat="identity") + 
      ggtitle(fastq[i]) + 
      coord_flip()
    ggsave(paste0(BAM_path, fastq[i], "_chromDistr.pdf"))
    
    #######
    # ChIP QC
    ########
    
    # using the pre-filtered BAM
    print(paste0("ChIP QC for unfiltered ", fastq[i], "..."))
    
    QCresult <- ChIPQCsample(reads = sorted_bam_post,
                             annotation = genome,
                             blacklist = blacklist_path)
    out <- capture.output(QCmetrics(QCresult))
    cat("ChIP QC Metrics for unfiltered", out, file = paste0(ChIPQC_path, "/", fastq[i], "_ChIPQC_output.txt"), sep = "\n", append = TRUE)
    
    myFlags <- flagtagcounts(QCresult)
    out <- capture.output(myFlags["DuplicateByChIPQC"]/myFlags["Mapped"])
    cat( out, file = paste0(ChIPQC_path, "/", fastq[i], "_ChIPQC_output.txt"), sep = "\n", append = TRUE)
    
    plotCC(QCresult) + ggtitle(fastq[i])
    ggsave(paste0(ChIPQC_path, "/", fastq[i], "_ChIPQC_CCplot.pdf")) 
    
    plotSSD(QCresult)+xlim(0,5) + ggtitle(fastq[i])
    ggsave(paste0(ChIPQC_path, "/", fastq[i], "_ChIPQC_SSDplot.pdf"))
    
    ######
    
    # make big wig
    
    print(paste0("Making BigWig for ", fastq[i], "..."))
    alignment <- readGAlignments(sorted_bam_post_mapq)
    reads_coverage <- coverage(alignment)
    
    export.bw(reads_coverage, con = paste0(bigwig_path, base, "_bigWig.bw"))
  }
  
  # this is from Tom's course and will produce a normalized bigwig where reads are scaled to the number of mapped reads multiplied by a million (reads per million)
  coverage_norm <- coverage(alignment,
                            weight = (10^6)/TotalMapped)
  export.bw(coverage_norm, paste0(bigwig_path, base, "_bigWig_norm.bw"))
  
  ##########
  # get plots over genes using soGGi
  ##########
  
  print(paste0("Making gene region signal plots from bigwigs..."))
  
  
  # get GRanges object of all genes
  
  
  if (genome == "hg19") {
    whole_gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  }
  
  if (genome == "hg38") {
    whole_gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
  
  if (genome == "mm9") {
    whole_gene <- genes(TxDb.Mmusculus.UCSC.mm9.knownGene)
  }
  
  if (genome == "mm10") {
    whole_gene <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
  }
  
  bigWigs <- list.files(paste0(fastq_path, "bigwigs"), full.names = TRUE)
  whole_gene_plots <- vector(mode = "list", length = length(bigWigs))
  
  names(whole_gene) <- NULL # needed to do this or the "percentOfRegion" setting below threw an error
  
  for (i in seq_along(bigWigs)) {
    
    whole_gene_plots[[i]] <- regionPlot(bigWigs[i],
                                        testRanges = whole_gene,
                                        style = "percentOfRegion",
                                        format = "bigwig")
    
  }
  
  # this will make a figure for each chip separately, need to work on code below in the chunk devoted to soGGi with either rbind, or the concatenation to get all in same plot. Basically it wouldn't properly concatenate in the for loop as I had it set up below. 
  for(i in seq_along(whole_gene_plots)) {
    plotRegion(whole_gene_plots[[i]])
    ggsave(paste0(strsplit(metadata(whole_gene_plots[[i]])$names, ".bw", fixed = TRUE)[[1]], ".pdf"))
  }
  
  
  ##############
  
  # Peak Calling
  
  # this will only run peaks if a sample sheet path is specified in the original call to the function
  
  ###############
  
  
  if (!missing(sample_sheet_path)) {
    # call peaks with macs2
    print("Calling Peaks...")
    
    data <- read.delim(sample_sheet_path, header = TRUE, sep = ",")
    
    if (missing(genome)) {
      macs2_genome <- "hs"
    }
    if (genome == "hg19") {
      macs2_genome <- "hs"
    }
    if (genome == "hg38") {
      macs2_genome <- "hs"
    }
    if (genome == "mm9") {
      macs2_genome <- "mm"
    }
    if (genome == "mm10") {
      macs2_genome <- "mm"
    }
    
    for (i in 1:nrow(data)) {
      if (data[i,4] == 0) {
        dir.create(paste0(fastq_path, data[i,1], "_peaks/"))
        path_output <- paste0(fastq_path, data[i,1], "_peaks/")
        input <- bam_files[grep(data[i,3], bam_files)]
        peaks <- bam_files[grep(data[i,2], bam_files)]
        
        macsCommand <- paste0("macs2 callpeak", 
                              " -t ", peaks, 
                              " -c ", input, 
                              " -g ", macs2_genome,
                              " --outdir ", path_output, 
                              " -n ",  data[i,1], 
                              " --bw ", bw, " -q ", q)  
      }
      else {
        dir.create(paste0(fastq_path, data[i,1], "_broad_peaks/"))
        path_output <- paste0(fastq_path, data[i,1], "_broad_peaks/")
        input <- bam_files[grep(data[i,3], bam_files)]
        peaks <- bam_files[grep(data[i,2], bam_files)]
        
        macsCommand <- paste0("macs2 callpeak", 
                              " -t ", peaks, 
                              " -c ", input, 
                              " -g ", macs2_genome, 
                              " --outdir ", path_output, 
                              " -n ",  data[i,1], 
                              " --broad", 
                              " --bw ", bw, " -q ", q)  
      }
      system(macsCommand)
      
      #genomic annotation using the annotatePeak function
      if (genome == "hg19") {
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      } else if (genome == "hg38") {
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      } else if (genome == "mm9") {
        txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
      } else if (genome == "mm10") {
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      }
      
      if (data[i,4] == 0) {
        peakAnnoList_marklist <- annotatePeak(peak = paste0(path_output, data[i,1], "_peaks.narrowPeak" ),
                                              TxDb = txdb,
                                              tssRegion = c(-5000, 3000),
                                              verbose = FALSE)
        plotAnnoBar(peakAnnoList_marklist)
        ggsave(paste0(path_output, data[i,1], "_annotation.pdf"), height = 4, width =  6)
      } else {
        peakAnnoList_marklist <- annotatePeak(peak = paste0(path_output, data[i,1], "_peaks.broadPeak" ),
                                              TxDb = txdb,
                                              tssRegion = c(-5000, 3000),
                                              verbose = FALSE)
        plotAnnoBar(peakAnnoList_marklist)
        ggsave(paste0(path_output, data[i,1], "_annotation.pdf"), height = 4, width =  6)
        
      }
    }
    
    
  }
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
}