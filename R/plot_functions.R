# Function to plot FDP-TPR curves
PlotFDPTPR <- function(res.list, contrast.vec, colors, sort.list, TPcol = "UPS", plotSVG = FALSE){
  
  cex.main <- 1.2
  cex.lab <- 1
  cex <- 1
  lwd <- 1
  
  if(plotSVG){
    cex.main <- 4.3
    cex.lab <- 3.3
    cex <- 3
    lwd <- 6
  }
  
  if(length(TPcol) == 1) {TPcol <- rep(TPcol, length(sort.list))}
  
  silent <- map(contrast.vec,
                
                ~{
                  contr <- .x
                  
                  if(plotSVG){
                    grDevices::svg(paste0(contr, ".svg"), width=10, height=10)
                  }
                  
                  main <- names(contrast.vec[contrast.vec == contr])
                  if(is.null(main)){main <- contr}
                  
                  plot(1, type="n", main = main, cex.main=cex.main, xlim=c(0,0.15), ylim=c(0,1), xlab="False discovery proportion", ylab="True positive rate", cex.lab=cex.lab, cex.axis=cex, cex=cex, las=1, frame.plot=FALSE)
                  
                  pmap(list(res.list, colors, sort.list, TPcol), function(x,y,z,TPcol)  {
                    
                    TP <- rlang::parse_quo(TPcol, env = rlang::caller_env())
                    
                    plot.obj <- subset(x, x$contrast == contr)
                    
                    if(!is.na(z[3])){
                      plot.obj[,z[3]] <- plot.obj %>% pull(z[3]) %>% abs %>% `::` ("dplyr","desc")(.)
                    }
                    
                    if(!is.na(z[4])){
                      plot.obj[,z[4]] <- plot.obj %>% pull(z[4]) %>% abs %>% `::` ("dplyr","desc")(.)
                    }
                    
                    z <- na.omit(z)
                    
                    max_ups <- plot.obj %>% pull(!!TP) %>% sum #/length(unique(plot.obj$contrast)) # 40 UPS proteins
                    plot.obj <- filter(plot.obj, !is.na(!!rlang::sym(z[1]))) %>% #group_by(contrast) %>% 
                      arrange(!!! rlang::syms(z)) %>%
                      mutate(FDP = cummean(!(!!TP)),
                             FPR = cumsum(!(!!TP))/sum(!(!!TP)),
                             TPR  = cumsum((!!TP))/max_ups,
                             UPS_hits = cumsum((!!TP)))
                    
                    lines(x = c(0, plot.obj$FDP), y = c(0, plot.obj$TPR), col = y, cex.lab = cex.lab, cex.axis = cex, cex = lwd, lwd = lwd, pch = 1)
                    
                    ### 1% FDR ###
                    FDP1 <- plot.obj$FDP[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.01])]
                    TPR1 <- plot.obj$TPR[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.01])]
                    if(length(FDP1) == 0){FDP1 <- 0}
                    if(length(TPR1) == 0){TPR1 <- 0}
                    if(FDP1 < 0.01){
                      pch = 17
                    } else{pch = 2}
                    points(FDP1, TPR1, col = y, pch = pch, cex = lwd, lwd = lwd)
                    ###
                    
                    ### 5% FDR ###
                    FDP5 <- plot.obj$FDP[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.05])]
                    TPR5 <- plot.obj$TPR[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.05])]
                    if(length(FDP5) == 0){FDP5 <- 0}
                    if(length(TPR5) == 0){TPR5 <- 0}
                    if(FDP5 < 0.05){
                      pch = 15
                    } else{pch = 0}
                    points(FDP5, TPR5, col = y, pch = pch, cex = lwd, lwd = lwd)
                    ###
                    
                    ### 10% FDR ###
                    FDP10 <- plot.obj$FDP[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.1])]
                    TPR10 <- plot.obj$TPR[which.max(pull(plot.obj, !!z[1])[pull(plot.obj, !!z[1]) < 0.1])]
                    if(length(FDP10) == 0){FDP10 <- 0}
                    if(length(TPR10) == 0){TPR10 <- 0}
                    if(FDP10 < 0.1){
                      pch = 19
                    } else{pch = 1}
                    points(FDP10, TPR10, col = y, pch = pch, cex = lwd, lwd = lwd)
                    ###
                    cat(paste0(contr), "\n")
                    cat(paste0(c(100*FDP1, 100*FDP5, 100*FDP10), "% \n"))
                  })
                  
                  abline(v = 0.01, lty = "dotted", lwd = lwd)
                  abline(v = 0.05, lty = "dotted", lwd = lwd)
                  abline(v = 0.1, lty = "dotted", lwd = lwd)
                  
                
                  if(plotSVG){
                    dev.off()
                  }
                  }
                
  )
  
  
}

# Function to plot scatter plots with histogram on the axis
scatterHist <- function(x, y, x2, y2, main = "", col1 = "#E15E9E", col2 = "#1B2944", col3 = "green", pch1 = 1, pch2 = 4, xlab = "Log2 odds ratio", ylab = "Log2 fold change", trueFC = NULL, square = FALSE, plotSVG = FALSE, filename = "scatterhist.svg"){
  
  cex.main <- 1.2
  cex.lab <- 1
  cex <- 1
  bottom <- 3
  left <- 3
  
  if(plotSVG){
    cex.main <- 4
    cex.lab <- 3
    cex <- 2.5
    bottom <- 5
    left <- 5
    
    grDevices::svg(filename, width=10, height=10)
  }
  
  def.par <- par(no.readonly = TRUE)
  
  zones = matrix(c(2,0,1,3), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(4/5,1/5), heights = c(1/5,4/5))
  
  xlim = c(min(density(c(x, x2), na.rm = TRUE)$x), max(density(c(x, x2), na.rm = TRUE)$x))
  ylim = c(min(density(c(y, y2), na.rm = TRUE)$x), max(density(c(y, y2), na.rm = TRUE)$x))
  
  xNA <- x[is.na(y)]
  x2NA <- x2[is.na(y2)]
  x <- x[!is.na(y)]
  x2 <- x2[!is.na(y2)]
  y <- na.omit(y)
  y2 <- na.omit(y2)
  
  xhist = hist(xNA, plot=FALSE, breaks = seq(xlim[1], xlim[2], by = (xlim[2]-xlim[1])/50))
  xhist2 = hist(x2NA, plot=FALSE, breaks = seq(xlim[1], xlim[2], by = (xlim[2]-xlim[1])/50))

  top <- max(c(xhist$density + xhist2$density), na.rm = TRUE)
  
  if(square){
    xlim[1] <- min(c(xlim[1], ylim[1]))
    xlim[2] <- max(c(xlim[2], ylim[2]))
    ylim[1] <- min(c(xlim[1], ylim[1]))
    ylim[2] <- max(c(xlim[2], ylim[2]))
  }
  
  # Main plot
  par(mar=c(bottom,left,1,1))
  plot(c(x,x2),c(y,y2), las=1, frame.plot=FALSE, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab
       , cex.lab=cex.lab, cex.axis=cex, cex.sub=cex, lwd=cex)
  
  if(!is.null(trueFC)){
    abline(h = trueFC, col = col3, lwd = cex)
    abline(v = trueFC, col = col3, lwd = cex)
  }
  
  points(x, y, las=1, col = col1, pch = pch1, cex.lab=cex.lab, cex.axis=cex, cex.sub=cex, lwd=cex)
  points(x2, y2, las=1, col = col2, pch = pch2, cex.lab=cex.lab, cex.axis=cex, cex.sub=cex, lwd=cex)
  
  # x-axis histogram
  par(mar=c(0,left,1,1))
  
  if(is.finite(top)){
    densTot <- na.omit(xhist$density)+na.omit(xhist2$density)
    densFP <- na.omit(xhist$density)
    if(length(densTot) != 0){barplot(densTot, axes=FALSE, ylim=c(0, top), space=0, col = col2, border = NA, main = main, cex.main=cex.main)}
    if(length(densFP) != 0){barplot(densFP, axes=FALSE, ylim=c(0, top), space=0, col = col1, border = NA, add = TRUE)}
  abline(h=0)
  }
  
  par(oma=c(bottom,left,0,0))
  
  par(def.par)
  if(plotSVG){
    dev.off()
  }
}

# Function to plot the underlying data for one or more proteins, with an option to save all plots in a PDF file
plot_proteins <- function(gene.names, title, msnset, pdf = TRUE){
  
  cex.main <- 1.2
  cex.lab <- 1
  cex <- 1
  lwd <- 2
  
  if(pdf){
    pdf(title, width=10, height=10)
  }
  
  colors_all <- rep(c("#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0",
                      "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#800000",
                      "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#000000"), 10000)
  pchs_all <- rep(c(0:25,33,35:38,47,48:57,60:90,97:122,126), 10000)
  
  for(i in 1:length(gene.names)){
    
    # Select expression data for protein
    data <- Biobase::exprs(msnset)[fData(msnset)$gene.name == gene.names[i],]
    data <- data[,c("LA3", "LA4", "LA8", "RA3", "RA4", "RA8", "SepA3", "SepA4", "SepA8", "LV3", "LV4", "LV8", "RV3", "RV4", "RV8", "SepV3", "SepV4", "SepV8")]
    data <- data.frame(data, sequence = row.names(data))
    
    # Convert data to long format
    long <- data %>% gather(sample, expression, LA3:SepV8, factor_key=TRUE) %>% na.omit
    
    colors <- colors_all[as.numeric(long$sequence)]
    pchs <- pchs_all[as.numeric(long$sequence)]
    
    if(any(is.na(colors))){stop("Colors NA!")}
    if(any(is.na(pchs))){stop("Pchs NA!")}
    
    
    plot(jitter((as.numeric(long$sample)), factor=0.5), long$expression, col = colors, pch=pchs, xaxt = "n", xlim = c(0, 19), ylab="Log2 peptide intensity", xlab="sample", main = gene.names[i], las=2, frame.plot=FALSE, cex = cex, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex, lwd = cex)
    axis(side = 1, at = 1:18, labels = levels(long$sample), cex = cex, cex.lab = cex.lab, cex.axis = cex, las = 2)
    
    abline(v = 9.5, lwd = lwd)
    
  }
  
  if(pdf){
    dev.off()
  }
  
}

# Function to plot the underlying data for four proteins and to save them in an SVG file

plot_proteinsSVG <- function(gene.names, msnset, title){
  
  grDevices::svg(title, width=40, height=40)
  cex.main <- 8.3
  cex.lab <- 6.3
  cex <- 5
  lwd <- 5
  par(mfrow=c(2,2))
  par(mar=c(32.1,17.1,10.1,2.1))
  
  
  colors_all <- rep(c("#e6194b", "#3cb44b", "#0082c8", "#f58231", "#911eb4", "#46f0f0",
                      "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#800000",
                      "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#000000"), 10000)
  pchs_all <- rep(c(0:25,33,35:38,47,48:57,60:90,97:122,126), 10000)
  
  for(i in 1:length(gene.names)){
    
    # Select expression data for protein
    data <- Biobase::exprs(msnset)[fData(msnset)$gene.name == gene.names[i],]
    data <- data[,c("LA3", "LA4", "LA8", "RA3", "RA4", "RA8", "SepA3", "SepA4", "SepA8", "LV3", "LV4", "LV8", "RV3", "RV4", "RV8", "SepV3", "SepV4", "SepV8")]
    data <- data.frame(data, sequence = row.names(data))
    
    # Convert data to long format
    long <- data %>% gather(sample, expression, LA3:SepV8, factor_key=TRUE) %>% na.omit
    
    colors <- colors_all[as.numeric(long$sequence)]
    pchs <- pchs_all[as.numeric(long$sequence)]
    
    if(any(is.na(colors))){stop("Colors NA!")}
    if(any(is.na(pchs))){stop("Pchs NA!")}
    
    
    plot(jitter((as.numeric(long$sample)), factor=0.5), long$expression, col = colors, pch=pchs, xaxt = "n", xlim = c(0, 19), ylab="Log2 peptide intensity", xlab="sample", main = gene.names[i], las=2, frame.plot=FALSE, cex = cex, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex, lwd = cex, ann = FALSE)
    axis(side = 1, at = 1:18, labels = levels(long$sample), cex = cex, cex.lab = cex.lab, cex.axis = cex, las = 2)
    
    mtext(side = 3, text = gene.names[i], line = 3, cex = cex.main, font=2)
    mtext(side = 1, text = "sample", line = 3, cex = cex.lab, padj = 4)
    mtext(side = 2, text = "Log2 peptide intensity", line = 3, cex = cex.lab, padj = -1.5)
    
    abline(v = 9.5, lwd = lwd)
  }
  dev.off()
}