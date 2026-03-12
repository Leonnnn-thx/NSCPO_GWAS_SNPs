rm(list = ls())
gc()

setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/new_input_atac")
library(gkmSVM) 
library(BSgenome.Hsapiens.UCSC.hg19.masked) 

# util for getting negative group 
genNullSeqs_HEPM <- function (inputBedFN, genomeVersion = "hg19", outputBedFN = "negSet.bed", outputPosFastaFN = "posSet.fa", outputNegFastaFN = "negSet.fa", 
                              xfold = 1, repeat_match_tol = 0.02, GC_match_tol = 0.02, 
                              length_match_tol = 0.02, batchsize = 5000, nMaxTrials = 20, 
                              genome = NULL) 
{
  if (is.null(genome)) {
    
    if (toupper(genomeVersion) == "DANRER7") {
      if (requireNamespace("BSgenome.Drerio.UCSC.danRer7.masked", 
                           quietly = TRUE)) {
        genome <- BSgenome.Drerio.UCSC.danRer7.masked::BSgenome.Drerio.UCSC.danRer7.masked
      }
    }
    else {
      if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19.masked", 
                           quietly = TRUE)) {
        genome <- BSgenome.Hsapiens.UCSC.hg19.masked::BSgenome.Hsapiens.UCSC.hg19.masked
      }
    }
  }
  seqnams = GenomeInfoDb::seqnames(genome)
  chrlens = GenomeInfoDb::seqlengths(genome)
  chrpos = cumsum(as.numeric(chrlens))
  pmax = max(chrpos)
  chrpos = c(chrpos, 1e+12)
  chrpos0 = c(0, chrpos)
  ichrA = as.character(names(chrlens))
  getichrpos = function(ipos) {
    j = order(ipos)
    ipos = sort(ipos)
    ci = 1
    res = rep(NA, length(ipos))
    for (i in 1:length(ipos)) {
      while (ipos[i] > chrpos[ci]) {
        ci = ci + 1
      }
      res[j[i]] = ci
    }
    return(res)
  }
  generateRandomGenSeqs = function(seqlens) {
    rpos = sample(pmax, length(seqlens), replace = TRUE)
    ichr1 = getichrpos(rpos)
    ichr2 = getichrpos(rpos + seqlens)
    jj = which(ichr1 != ichr2)
    while (length(jj) > 0) {
      rpos[jj] = sample(pmax, length(jj), replace = TRUE)
      ichr1 = getichrpos(rpos)
      ichr2 = getichrpos(rpos + seqlens)
      jj = which(ichr1 != ichr2)
    }
    chr = ichrA[ichr1]
    start = rpos - chrpos0[ichr1]
    names <- chr
    ranges <- IRanges::IRanges(start = start, width = seqlens)
    strand <- BiocGenerics::strand(sample(c("+", "-"), length(names), 
                                          replace = TRUE))
    gr <- GenomicRanges::GRanges(seqnames = names, ranges = ranges, 
                                 strand = strand)
  }
  inBed = rtracklayer::import.bed(inputBedFN)
  inbed = GenomicRanges::as.data.frame(inBed)
  jj = which(is.na(match(as.character(inbed$seqnames), as.character(seqnams))))
  if (length(jj) > 0) {
    cat(paste("ERROR: Chromosome name not recognized for", 
              length(jj), "sequences.\n"))
    cat(unique(as.character(inbed$seqnames[jj])))
    return(NULL)
  }
  jj = which(inbed$end > GenomeInfoDb::seqlengths(genome)[as.character(inbed$seqnames)])
  if (length(jj) > 0) {
    cat("ERROR: Region outside chromosome. (Check the genome version) \n")
    print(inbed[jj, ])
    return(NULL)
  }
  gcContent <- function(seqs) {
    alf <- Biostrings::alphabetFrequency(seqs, as.prob = TRUE)
    gc = rowSums(alf[, c("G", "C"), drop = FALSE])
  }
  matchSeqs = function(gc1, gc2, len1, len2, rpt1, rpt2, gc_th = 0.02, 
                       len_th = 0.02, rpt_th = 0.02) {
    len_th = len_th * len1
    i1 = order(gc1)
    i2 = order(gc2)
    gc1 = gc1[i1]
    gc2 = gc2[i2]
    len1 = len1[i1]
    len2 = len2[i2]
    rpt1 = rpt1[i1]
    rpt2 = rpt2[i2]
    gc2 = c(gc2, 1e+10)
    len_th = len_th[i1]
    m2 = 1
    N = length(i1)
    N2 = length(i2)
    mtc1 = rep(NA, N)
    mtc2 = rep(0, length(i2))
    for (i in 1:N) {
      gc1i = gc1[i]
      len1i = len1[i]
      rpt1i = rpt1[i]
      len_thi = len_th[i]
      while (gc1i - gc2[m2] > gc_th) {
        m2 = m2 + 1
      }
      if (m2 <= N2) {
        m2b = m2
        while (gc2[m2b] - gc1i <= gc_th) {
          if ((mtc2[m2b] == 0) & (abs(len1i - len2[m2b]) <= 
                                  len_thi) & (abs(rpt1i - rpt2[m2b]) <= rpt_th)) {
            mtc2[m2b] = i
            mtc1[i] = m2b
            if (m2b == m2) {
              m2 = m2 + 1
            }
            break
          }
          m2b = m2b + 1
        }
      }
      else {
        break
      }
    }
    mtc1 = i2[mtc1]
    res = rep(NA, N)
    res[i1] = mtc1
    return(res)
  }
  repeatRat = function(bed) {
    chrs = unique(GenomeInfoDb::seqnames(bed))
    rpts = rep(0, length(bed))
    for (ichr in as.character(chrs)) {
      seq = genome[[ichr]]
      rpt = Biostrings::masks(seq)[["TRF"]]
      jj = which(as.character(GenomeInfoDb::seqnames(bed)) == 
                   ichr)
      if (length(jj) > 0) {
        jbed = bed[jj]
        jrpts = rep(0, length(jj))
        olaps <- IRanges::findOverlaps(rpt, jbed@ranges)
        qdf = GenomicRanges::as.data.frame(rpt)[S4Vectors::queryHits(olaps), 
        ]
        isect <- IRanges::pintersect(IRanges::IRanges(start = qdf$start, 
                                                      end = qdf$end), jbed@ranges[S4Vectors::subjectHits(olaps)])
        jres = S4Vectors::subjectHits(olaps)
        olap_width = BiocGenerics::width(isect)
        for (i in 1:length(jres)) {
          jrpts[jres[i]] = jrpts[jres[i]] + olap_width[i]
        }
        rpts[jj] = jrpts
      }
    }
    rpts = rpts/BiocGenerics::width(bed)
  }
  inbed = GenomicRanges::as.data.frame(inBed)
  cat(" importing sequences for", inputBedFN, "from", GenomeInfoDb::bsgenomeName(genome), 
      "\n")
  inSeqs = Biostrings::getSeq(genome, inBed)
  seqlens = inbed$width
  inGC = gcContent(inSeqs)
  cat(" calculating repeat distributions\n")
  inRpt = repeatRat(inBed)
  nout = round(nrow(inbed) * xfold)
  outbed = matrix(ncol = ncol(inbed), nrow = nout)
  outSeq = rep(inSeqs, length = nout)
  colnames(outbed) = colnames(inbed)
  unmatched = 1:length(outSeq)
  desGC = rep(inGC, length = nout)
  desRpt = rep(inRpt, length = nout)
  desLens = rep(seqlens, length = nout)
  for (iter in 1:nMaxTrials) {
    if (length(unmatched) > 0) {
      cat(" Trial", iter, "out of", nMaxTrials, "\n")
      rndBed = generateRandomGenSeqs(rep(desLens[unmatched], 
                                         length.out = batchsize))
      rndbed = GenomicRanges::as.data.frame(rndBed)
      cat(" importing sequences\n")
      rndSeqs = Biostrings::getSeq(genome, rndBed)
      rndGC = gcContent(rndSeqs)
      cat(" calculating repeat distributions\n")
      rndRpt = repeatRat(rndBed)
      cat(" matching sequences\n")
      mtc = matchSeqs(desGC[unmatched], rndGC, desLens[unmatched], 
                      BiocGenerics::width(rndBed), desRpt[unmatched], 
                      rndRpt, gc_th = GC_match_tol, len_th = length_match_tol, 
                      rpt_th = repeat_match_tol)
      jj = which(!is.na(mtc))
      if (length(jj) > 0) {
        outbed[unmatched[jj], 1:5] = as.matrix(rndbed[mtc[jj], 
        ])
        outSeq[unmatched[jj], ] = rndSeqs[mtc[jj], ]
        unmatched = unmatched[-jj]
      }
      cat(nrow(outbed) - length(unmatched), "sequences found so far, ", 
          length(unmatched), " remaining.\n")
    }
  }
  if (length(unmatched) > 0) {
    outbed = outbed[-unmatched, ]
    outSeq = outSeq[-unmatched, ]
  }
  outbed = gsub(" ", "", outbed)
  write.table(as.matrix(outbed[, 1:3]), quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE, file = outputBedFN)
  if (requireNamespace("seqinr", quietly = TRUE)) {
    outseqnams = paste(outbed[, 1], outbed[, 2], outbed[, 
                                                        3], "neg", 1:nrow(outbed), sep = "_")
    seqinr::write.fasta(sequences = sapply(as.character(outSeq), 
                                           strsplit, ""), names = outseqnams, file.out = outputNegFastaFN)
    inseqnams = paste(as.character(inbed[, 1]), inbed[, 2], 
                      inbed[, 3], "pos", 1:nrow(inbed), sep = "_")
    seqinr::write.fasta(sequences = sapply(as.character(inSeqs), 
                                           strsplit, ""), names = inseqnams, file.out = outputPosFastaFN)
  }
  return(outputNegFastaFN)
}

# get negative control group
genNullSeqs_HEPM('HEPM_AEs_FINAL.bed',nMaxTrials=10,xfold=1,genomeVersion='hg19', outputPosFastaFN='HEPM_AEs_FINAL_pos.fa', outputBedFN='HEPM_AEs_FINAL_neg1x.bed', outputNegFastaFN='HEPM_AEs_FINAL_neg1x.fa')

gkmsvm_kernel('HEPM_AEs_FINAL_pos.fa','HEPM_AEs_FINAL_neg1x.fa', 'HEPM_AEs_kernel.out')

gkmsvm_train('HEPM_AEs_kernel.out', 'HEPM_AEs_FINAL_pos.fa', 'HEPM_AEs_FINAL_neg1x.bed', svmfnprfx= 'HEPM_AEs'); #trains SVM

gkmsvm_trainCV('HEPM_AEs_kernel.out','HEPM_AEs_FINAL_pos.fa','HEPM_AEs_FINAL_neg1x.bed',svmfnprfx='HEPM_AEs',outputCVpredfn='HEPM_AEs_cvpred.out',outputROCfn='HEPM_AEs_roc.out')

# Figure 3f

library(ROCR)
library(readr)

### Prepare the function for calculating auPRC (subfunction of gkmsvm_trainCV())
auPRC <- function(perf) {
  rec <- perf@x.values
  prec <- perf@y.values
  result <- list()
  for (i in 1:length(perf@x.values)) {
    result[i] <- list(sum((rec[[i]][2:length(rec[[i]])] - 
                             rec[[i]][2:length(rec[[i]]) - 1]) * prec[[i]][-1]))
  }
  return(result)
}
## Evaluate HEPM_training kernel ####
HEPM_training <- read_delim("HEPM_AEs_cvpred.out", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
colnames(HEPM_training) <-c("sequenceid", "preds","labs","nCV")
pred <- prediction( HEPM_training$preds, HEPM_training$labs)
perf <- performance(pred,"tpr","fpr")  
pdf(file="HEPM_AEs_ROC.pdf", width=4.5, height=7)
plot(perf, colorize=TRUE) 
dev.off()
perf <- performance(pred, "auc") # calculate auROC
perf@y.values[[1]]   # return the auROC=0.856238

## precision/recall curve (x-axis: recall, y-axis: precision)
perf1 <- performance(pred, "prec", "rec")
pdf(file="HEPM_AEs_PRC.pdf", width=4.5, height=7)
plot(perf1,colorize=TRUE)  
dev.off()

auPRC(perf1)   # return the auPRC=0.8738478


# Figure 3g

library(ggplot2)
library(RColorBrewer)
library(ggpubr)

gkmsvm_classify('VistaEnh_human.fa',svmfnprfx='HEPM_AEs', 'VistaEnh_human_weight.out') 
gkmsvm_classify('VistaNoEnh_human_unique.fa',svmfnprfx='HEPM_AEs', 'VistaNoEnh_human_weight.out')
gkmsvm_classify('Vista_Human_positive_facialmesenchymeEnh.fa',svmfnprfx='HEPM_AEs', 'Vista_facial_enh_score.csv')


plot.format=theme(plot.background=element_blank(),panel.grid=element_blank(),panel.background=element_blank(),panel.border=element_rect(color="black",linewidth=0.5,fill=NA),axis.line=element_blank(),axis.ticks=element_line(color="black",linewidth=0.5),axis.text=element_text(color="black",size=9),axis.title=element_text(color="black",size=9),plot.title=element_text(color="black",size=12),legend.background=element_blank(),legend.key=element_blank(),legend.text=element_text(color="black",size=9),legend.title=element_text(color="black",size=7))

mydata <-read.csv("Vista_score_3group.csv", header = T, sep = ",")

compaired_all <- list(
  c("VistaEnh", "VistaNoEnh"),
  c("VistaEnh", "Vista_facial_enh"),
  c("VistaNoEnh", "Vista_facial_enh")
)

num_classes <- length(unique(mydata$class))
colors <- brewer.pal(max(3, num_classes), "Set2")[1:num_classes]


P1 <- ggplot(mydata, aes(x = class, y = value)) +
  geom_violin(aes(fill = class), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_manual(values = colors) +
  stat_compare_means(comparisons = compaired_all, 
                     method = "wilcox.test",
                     label = "p.format") +
  scale_y_continuous(breaks = c(100, 125, 150, 175, 200)) +
  theme_light() + 
  theme(
    panel.grid.major.y = element_line(color = "gray70", 
                                      size = 0.5, 
                                      linetype = "dashed"),  
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  plot.format

P1
ggsave("Fig3g.vista_facialenh_totalenh_noenh_violin.pdf", plot = P1, width = 10, height = 5)
