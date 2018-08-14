### R code from vignette source 'exomeCopy.Rnw'

###################################################
### code chunk number 1: exomeCopy.Rnw:18-19
###################################################
options(width = 70)

rm(list =ls())
###################################################
### code chunk number 2: exomeCopy.Rnw:55-77 (eval = FALSE)
###################################################
# library(exomeCopy)
# target.file <- "targets.bed"
# bam.files <- c("/path/to/file1.bam", "/path/to/file2.bam", "/path/to/file3.bam")
# sample.names <- c("sample1","sample2","sample3")
# reference.file <- "/path/to/reference_genome.fa"
# target.df <- read.delim(target.file, header = FALSE)
# target <- GRanges(seqname = target.df[,1], IRanges(start = target.df[,2] + 1, end = target.df[,3]))
# counts <- target
# for (i in 1:length(bam.files)) {
#   mcols(counts)[[sample.names[i]]] <- countBamInGRanges(bam.files[i], target)
# }
# counts$GC <- getGCcontent(target, reference.file)
# counts$GC.sq <- counts$GC^2
# counts$bg <- generateBackground(sample.names, counts, median)
# counts$log.bg <- log(counts$bg + .1)
# counts$width <- width(counts)
# fit.list <- lapply(sample.names, function(sample.name) {
#   lapply(seqlevels(target), function(seq.name) {
#     exomeCopy(counts[seqnames(counts) == seq.name], sample.name, X.names = c("log.bg", "GC", "GC.sq","width"), S = 0:4, d = 2)
#   })
# })
# compiled.segments <- compileCopyCountSegments(fit.list)


###################################################
### code chunk number 3: exomeCopy.Rnw:90-93
###################################################
library(exomeCopy)
gr <- GRanges(seqname="seq1",IRanges(start=1,end=345))
subdivideGRanges(gr)


###################################################
### code chunk number 4: exomeCopy.Rnw:99-106
###################################################
  plot(0,0,xlim=c(0,500),ylim=c(0,25),type="n",yaxt="n",ylab="",xlab="width of input GRanges object",main="Effect of subdivideGRanges")
  abline(v=1:5*100,col="grey")
  for (i in 1:24) {
    gr <- GRanges(seqname="chr1",IRanges(start=1,width=(i*20)))
    sbd.gr <- subdivideGRanges(gr)
    arrows(start(sbd.gr),rep(i,length(sbd.gr)),end(sbd.gr),rep(i,length(sbd.gr)),length=.04,angle=90,code=3)
  }


###################################################
### code chunk number 5: exomeCopy.Rnw:112-119
###################################################
target.file <- system.file("extdata","targets.bed",package="exomeCopy")
target.df <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end")) 
target <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
target
target <- reduce(target, min.gapwidth=0)
target.sub <- subdivideGRanges(target)
target.sub


###################################################
### code chunk number 6: exomeCopy.Rnw:128-137
###################################################
bam.file <- system.file("extdata","mapping.bam",package="exomeCopy")
scanBamHeader(bam.file)[[1]]$targets
seqlevels(target.sub)
toy.counts <- target.sub
sample.df <- data.frame(samples="sample1",bam.files=bam.file,stringsAsFactors=FALSE)
for (i in 1:nrow(sample.df)) {
  mcols(toy.counts)[[sample.df$samples[i]]] <- countBamInGRanges(sample.df$bam.files[i],target.sub)
}
toy.counts


###################################################
### code chunk number 7: exomeCopy.Rnw:143-146
###################################################
reference.file <- system.file("extdata","reference.fa",package="exomeCopy")
toy.counts$GC <- getGCcontent(target.sub, reference.file)
toy.counts


###################################################
### code chunk number 8: exomeCopy.Rnw:155-158
###################################################
data(exomecounts)
length(exomecounts)
head(exomecounts)


###################################################
### code chunk number 9: exomeCopy.Rnw:164-165
###################################################
plot(start(exomecounts),exomecounts$HG00551,xlim=c(0.8e6,1.8e6),xlab="genomic position",ylab="counts",main="HG00551 read counts in exonic ranges")


###################################################
### code chunk number 10: exomeCopy.Rnw:173-176
###################################################
chr1a <- exomecounts
seqlevels(chr1a) <- "chr1a"
example.counts <- c(exomecounts, chr1a)


###################################################
### code chunk number 11: exomeCopy.Rnw:183-187
###################################################
exome.samples <- grep("HG.+",colnames(mcols(example.counts)),value=TRUE)
example.counts$bg <- generateBackground(exome.samples, example.counts, median)
example.counts$log.bg <- log(example.counts$bg + .1)
example.counts$bg.var <- generateBackground(exome.samples, example.counts, var)


###################################################
### code chunk number 12: exomeCopy.Rnw:192-194
###################################################
summary(example.counts$bg)
example.counts <- example.counts[example.counts$bg > 0,]


###################################################
### code chunk number 13: exomeCopy.Rnw:199-201
###################################################
example.counts$GC.sq <- example.counts$GC^2
example.counts$width <- width(example.counts)


###################################################
### code chunk number 14: exomeCopy.Rnw:245-257
###################################################
simulateCNV <- function(x,indices,multiply,prob) {
  x[indices] <- x[indices] + multiply * rbinom(length(indices),prob=prob,size=x[indices])
  return(x)
}
set.seed(2)
cnv.probs <- rep(c(.99,.5,.5,.95),each=2)
cnv.mult <- rep(c(-1,1),each=4)
cnv.starts <- rep(c(1,301,601,901),each=2)
for (i in 1:8) {
  mcols(example.counts)[[exome.samples[i]]] <- simulateCNV(mcols(example.counts)[[exome.samples[i]]],cnv.starts[i]:(cnv.starts[i] + 99),multiply=cnv.mult[i],prob=cnv.probs[i])
  mcols(example.counts)[[exome.samples[i]]] <- simulateCNV(mcols(example.counts)[[exome.samples[i]]],1000 + cnv.starts[i]:(cnv.starts[i] + 99),multiply=cnv.mult[i],prob=cnv.probs[i])
}

example.counts

###################################################
### code chunk number 15: exomeCopy.Rnw:266-268
###################################################
  fit <- exomeCopy(example.counts[seqnames(example.counts) == "chr1",],sample.name=exome.samples[3],X.names=c("log.bg","GC","GC.sq","width"),S=0:6,d=2)
  show(fit)

  function (rdata =example.counts[seqnames(example.counts) == "chr1",] , sample.name=exome.samples[3],
            X.names=c("log.bg","GC","GC.sq","width"), Y.names, fit.var = FALSE, 
            reltol = 1e-04, S = 0:4, d = 2, goto.cnv = 1e-04, goto.normal = 1/20, 
            init.phi = "norm") 
  {
    if (!sample.name %in% colnames(rdata)) {
      stop("sample.name is not a column in rdata.")
    }
    O <- rdata[[sample.name]]
    if (any(O != round(O) | O < 0)) {
      stop("Sample counts must be non-negative integers")
    }
    if (mean(O == 0) > 0.9) {
      warning("More than 90% of sample counts are zero, exomeCopy will return NULL")
      return(NULL)
    }
    if (any(S != round(S) | S < 0) | any(d != round(d) | d < 
                                         0)) {
      stop("S and d must be non-negative integers")
    }
    if (!all(d %in% S)) {
      stop("The normal state, d, must be one of the possible copy states in S")
    }
    if (!(all(X.names %in% colnames(rdata)))) {
      stop("all X.names must be columns in rdata")
    }
    if (length(rdata) != 1) {
      stop("The current rdata argument has ranges over the following chromosomes/spaces: ", 
           paste(names(rdata), collapse = ", "), ".\n The rdata argument should contain ranged data over a single space.\n You can pass exomeCopy a ranged data from a single space with this syntax: rdata['chr1'].\n See vignette for example of running on multiple chromosomes.")
    }
    if (nrow(rdata) < 100) {
      warning("exomeCopy was tested for thousands of ranges covering a targeted region on a single chromosome.  The results might not be reliable for less than a hundred ranges.")
    }
    if (is.unsorted(start(rdata))) {
      stop("Genomic ranges for ranged data must be sorted.")
    }
    normal.state = which(S == d)
    X <- as.matrix(as.data.frame(values(rdata))[, X.names])
    colnames(X) <- X.names
    if (fit.var) {
      if (!(all(Y.names %in% colnames(rdata)))) {
        stop("Y.names must be variable names in rdata")
      }
      Y <- as.matrix(as.data.frame(values(rdata))[, Y.names])
      colnames(Y) <- Y.names
    }
    controls <- list(reltol = reltol, maxit = 10000)
    X.full <- cbind(intercept = rep(1, nrow(X)), scale(X))
    lmfit <- lm(log(O + 0.1) ~ X.full + 0)
    beta.hat <- lmfit$coefficients
    names(beta.hat) <- colnames(X.full)
    if (init.phi == "norm") {
      phi.hat <- (var(O - exp(lmfit$fitted)) - mean(O))/mean(O)^2
    }
    else if (init.phi == "counts") {
      phi.hat <- (var(O) - mean(O))/mean(O)^2
    }
    phi.hat <- ifelse(phi.hat < 1e-06, 1e-06, phi.hat)
    init.par <- list(goto.cnv = goto.cnv, goto.normal = goto.normal, 
                     beta.hat = beta.hat, phi.hat = phi.hat)
    fx.par <- list(S = S, d = d, normal.state = normal.state, 
                   fit.var = fit.var)
    nstates <- length(S)
    if (!fit.var) {
      data <- list(O = O, X = X.full)
      finite.test <- is.finite(negLogLike(c(logit(goto.cnv), 
                                            logit(goto.normal), beta.hat, log(phi.hat)), fx.par, 
                                          data, nstates, stFn, trFn, emFn))
      if (!finite.test) {
        beta.hat <- rep(0, length(beta.hat))
        beta.hat[1] <- log(mean(data$O) + 0.1)
      }
      nm.fit <- optim(c(logit(goto.cnv), logit(goto.normal), 
                        beta.hat, log(phi.hat)), function(par) negLogLike(par, 
                                                                          fx.par, data, nstates, stFn, trFn, emFn), method = "Nelder-Mead", 
                      control = controls)
    }
    else {
      Y.full <- cbind(intercept = rep(1, nrow(Y)), scale(Y))
      data <- list(O = O, X = X.full, Y = Y.full)
      gamma.hat <- c(log(phi.hat), rep(0, ncol(Y)))
      nm.fit <- optim(c(logit(goto.cnv), logit(goto.normal), 
                        beta.hat, gamma.hat), function(par) negLogLike(par, 
                                                                       fx.par, data, nstates, stFn, trFn, emFn), method = "Nelder-Mead", 
                      control = controls)
    }
    goto.cnv.hat <- logistic(nm.fit$par[1])
    goto.normal.hat <- logistic(nm.fit$par[2])
    A <- trFn(nm.fit$par[1:2], fx.par, data, nstates)
    beta.hat <- nm.fit$par[3:(2 + ncol(X.full))]
    names(beta.hat) <- colnames(X.full)
    mu.hat <- exp(as.numeric(X.full %*% beta.hat))
    if (!fit.var) {
      phi.hat <- exp(nm.fit$par[(3 + ncol(X.full))])
      gamma.hat <- NULL
      type = "exomeCopy"
    }
    else {
      gamma.hat <- nm.fit$par[(3 + ncol(X.full)):length(nm.fit$par)]
      names(gamma.hat) <- colnames(Y.full)
      phi.hat <- exp(as.numeric(Y.full %*% gamma.hat))
      phi.hat[phi.hat < 1e-06] <- 1e-06
      type = "exomeCopyVar"
    }
    path <- viterbiPath(nm.fit$par, fx.par, data, nstates, stFn, 
                        trFn, emFn)
    emit.probs <- emFn(nm.fit$par, fx.par, data, nstates)
    log.odds <- log(emit.probs[cbind(path, seq(path))] + 1e-06) - 
      log(emit.probs[normal.state, ] + 1e-06)
    final.par <- list(goto.cnv = goto.cnv.hat, goto.normal = goto.normal.hat, 
                      beta = beta.hat, gamma = gamma.hat, phi = phi.hat)
    fit <- new("ExomeCopy", sample.name = sample.name, type = type, 
               path = Rle(path), ranges = ranges(rdata), O.norm = as.numeric(O/mu.hat), 
               log.odds = log.odds, fx.par = fx.par, init.par = init.par, 
               final.par = final.par, counts = nm.fit$counts, convergence = nm.fit$convergence, 
               nll = nm.fit$value)
    return(fit)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
###################################################
### code chunk number 16: exomeCopy.Rnw:273-274
###################################################
  copyCountSegments(fit)


###################################################
### code chunk number 17: exomeCopy.Rnw:280-282
###################################################
  cnv.cols <- c("red","orange","black","deepskyblue","blue","blue2","blue4")
  plot(fit,col=cnv.cols)


###################################################
### code chunk number 18: exomeCopy.Rnw:290-293
###################################################
runExomeCopy <- function(sample.name,seqs) {
  lapply(seqs,function(seq.name) exomeCopy(example.counts[seqnames(example.counts) == seq.name],sample.name,X.names=c("log.bg","GC","GC.sq","width"),S=0:4,d=2))
}


###################################################
### code chunk number 19: exomeCopy.Rnw:298-308
###################################################
seqs <- c("chr1","chr1a")
names(seqs) <- seqs
samples <- exome.samples[1:8]
names(samples) <- samples
t0 <- as.numeric(proc.time()[3])
fit.list <- lapply(samples,runExomeCopy,seqs)
t1 <- as.numeric(proc.time()[3])
time.elapsed <- as.numeric(t1 - t0)
paste(round(time.elapsed),"seconds for",length(samples),"samples,",round(sum(width(example.counts))/1e3),"kb of target")
paste("~",round(time.elapsed / 60 / (8 * sum(width(example.counts))) * 32e6,1)," minutes for 1 sample, 32 Mb of target",sep="")


###################################################
### code chunk number 20: exomeCopy.Rnw:315-317
###################################################
compiled.segments <- compileCopyCountSegments(fit.list)
CNV.segments <- compiled.segments[compiled.segments$copy.count != 2]


###################################################
### code chunk number 21: exomeCopy.Rnw:325-328
###################################################
CNV.segments[1:6]
table(CNV.segments$nranges)
CNV.segments <- CNV.segments[CNV.segments$nranges > 5]


###################################################
### code chunk number 22: exomeCopy.Rnw:334-336
###################################################
CNV.overlaps.matrix <- as.matrix(findOverlaps(CNV.segments,drop.self=TRUE))
head(CNV.overlaps.matrix)


###################################################
### code chunk number 23: exomeCopy.Rnw:342-346
###################################################
par(mfrow=c(2,1),mar=c(4,3,2,1))
cnv.cols <- c("red","orange","black","deepskyblue","blue")
plotCompiledCNV(CNV.segments=CNV.segments, seq.name="chr1", col=cnv.cols)
plotCompiledCNV(CNV.segments=CNV.segments, seq.name="chr1a", col=cnv.cols)


###################################################
### code chunk number 24: exomeCopy.Rnw:354-356
###################################################
fit.list[[1]][[1]]@init.par$beta.hat
fit.list[[1]][[1]]@final.par$beta


###################################################
### code chunk number 25: exomeCopy.Rnw:361-362
###################################################
sessionInfo()


