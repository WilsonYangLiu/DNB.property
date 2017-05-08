# This script is used to create a Boxplot with some facinating property:
# 1. Instead of showing only average of SD, PCC_in, PCC_out respectively, 
#    we recalculate all SD, PCC_in, PCC_out for known DNB, compare with control. 
# 2. We use these values to create the Boxplot, and map values of all SD, PCC_in, 
#    PCC_out to this plot.
# The `maps` is how the controls should be chosen for a sepecific treatment, say, 
# w3 to w3_c, w5 to w5_c (w* is the treatment sample, w*_c is the coresponding control sample)

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/')

#1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
Gene.count <- read.csv(file = './0.Backup/[!]fpkm/FPKM_after.txt', sep = '\t', row.names = 1)
wx.count <- list() # For disease 
wx.count.c <- list() # For control
wx.count$w3 <- Gene.count[, grep(pattern = 'X3', x = colnames(Gene.count), value = TRUE)]
wx.count$w5 <- Gene.count[, grep(pattern = 'X5', x = colnames(Gene.count), value = TRUE)]
wx.count$w9 <- Gene.count[, grep(pattern = 'X9', x = colnames(Gene.count), value = TRUE)]
wx.count$w14 <- Gene.count[, grep(pattern = 'X14', x = colnames(Gene.count), value = TRUE)]
wx.count$w17 <- Gene.count[, grep(pattern = 'X17', x = colnames(Gene.count), value = TRUE)]
wx.count.c$wx <- Gene.count[, grep(pattern = 'c', x = colnames(Gene.count), value = TRUE)]
#1-----#

# Gene.count <- read.csv(file = './2+.Count2FPKM/fpkm.csv', row.names = 1)
# Gene.count <- read.csv(file = './2.count/htseq/count.txt', sep = '\t', row.names = 1)
# wx.count <- list() # For disease 
# wx.count.c <- list() # For control
# wx.count$w3 <- Gene.count[, grep(pattern = 'X3WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w5 <- Gene.count[, grep(pattern = 'X5WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w9 <- Gene.count[, grep(pattern = 'X9WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w14 <- Gene.count[, grep(pattern = 'X14WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w17 <- Gene.count[, grep(pattern = 'X17WT', x = colnames(Gene.count), value = TRUE)]
# wx.count.c$wx <- Gene.count[, grep(pattern = 'N|C', x = colnames(Gene.count), value = TRUE)]

## mapping between disease and control
maps <- rep('wx', 5); names(maps) <- names(wx.count)

DNB <- read.table(file = './3.2.dnb/module_gene.txt', stringsAsFactors = FALSE)$V1
#-- for week x (3, 5, 9, 14, 17)
SD <- list()
cor.in <- list()
cor.out <- list()
for (w in names(wx.count)) {
  wx.DNB <- wx.count[[w]][DNB, ]
  wx.others <- wx.count[[w]][!(rownames(wx.count[[w]]) %in% DNB), ]
  
  SD.tmp <- as.numeric(apply(wx.DNB, MARGIN = 1, FUN = function(x) sd(x)))
  names(SD.tmp) <- rownames(wx.DNB)
  SD[[w]] <- SD.tmp
  
  cor.in.tmp <- abs(cor(t(wx.DNB)))
  cor.in.tmp[!(lower.tri(cor.in.tmp))] <- 0
  cor.in[[w]] <- cor.in.tmp
  
  cor.out.tmp <- abs(cor(x = t(wx.DNB), y = t(wx.others), method = "spearman"))
  cor.out[[w]] <- cor.out.tmp
}

SD.c <- list()
cor.in.c <- list()
cor.out.c <- list()
for (w in names(wx.count.c)) {
  wx.DNB <- wx.count.c[[w]][DNB, ]
  wx.others <- wx.count.c[[w]][!(rownames(wx.count.c[[w]]) %in% DNB), ]
  
  SD.tmp <- as.numeric(apply(wx.DNB, MARGIN = 1, FUN = function(x) sd(x)))
  names(SD.tmp) <- rownames(wx.DNB)
  SD.c[[w]] <- SD.tmp
  
  cor.in.tmp <- abs(cor(t(wx.DNB)))
  cor.in.tmp[!(lower.tri(cor.in.tmp))] <- 0
  cor.in.c[[w]] <- cor.in.tmp
  
  cor.out.tmp <- abs(cor(x = t(wx.DNB), y = t(wx.others), method = "spearman"))
  cor.out.c[[w]] <- cor.out.tmp
}

SD.adj <- list()
cor.in.adj <- list()
cor.out.adj <- list()
for (m in 1:length(maps)) {
  SD.adj[[names(maps[m])]] <- SD[[names(maps[m])]] / (SD.c[[maps[m]]] + 0.01)
  
  cor.in.tmp <- cor.in[[names(maps[m])]] / (cor.in.c[[maps[m]]] + 0.1)
  cor.in.tmp <- cor.in.tmp[lower.tri(cor.in.tmp)]
  cor.in.tmp <- cor.in.tmp[!(is.infinite(cor.in.tmp)) & !(is.nan(cor.in.tmp))]
  cor.in.adj[[names(maps[m])]] <- cor.in.tmp
  
  cor.out.tmp <- c(cor.out[[names(maps[m])]] / (cor.out.c[[maps[m]]]) + 0.1)
  cor.out.tmp <- cor.out.tmp[!(is.infinite(cor.out.tmp)) & !(is.nan(cor.out.tmp))]
  cor.out.adj[[names(maps[m])]] <- cor.out.tmp
}


### For SD
critical.time <- 3
nTime <- length(SD.adj)
boxplot(SD.adj, lwd = 2,
        names = c(names(wx.count)),
        xlab = 'Time', ylab = '', 
        main = 'SD_adj',
        outline = FALSE, add = FALSE)
# stripchart(SD.adj, vertical = TRUE, method = "jitter", add = TRUE, 
#            pch = 20, col = c('grey', 'blue')[(1:nTime == critical.time) + 1])
lines(1:nTime, 
      sapply(names(wx.count), FUN = function(x) mean(SD.adj[[x]])), 
      type = 'b', 
      lwd = 3, col = 'red', pch = 1)

### For cor.in ?
boxplot(cor.in.adj, lwd = 2, 
        names = c(names(wx.count)),
        xlab = 'Time', ylab = '', 
        main = 'PCC_in_adj',
        outline = FALSE, add = FALSE)
# stripchart(cor.in.adj, vertical = TRUE, method = "jitter", add = TRUE, 
#            pch = 20, col = c('grey', 'blue')[(1:nTime == critical.time) + 1])
lines(1:nTime, 
      sapply(names(wx.count), FUN = function(x) mean(cor.in.adj[[x]])), 
      type = 'b', 
      lwd = 3, col = 'red', pch = 1)

### For cor.out ?
boxplot(cor.out.adj, lwd = 2, 
        names = c(names(wx.count)),
        xlab = 'Time', ylab = '', 
        main = 'PCC_out_adj',
        outline = FALSE, add = FALSE)
# stripchart(cor.out.adj, vertical = TRUE, method = "jitter", add = TRUE, 
#            pch = 20, col = c('grey', 'blue')[(1:nTime == critical.time) + 1])
lines(1:nTime, 
      sapply(names(wx.count), FUN = function(x) mean(cor.out.adj[[x]])), 
      type = 'b', 
      lwd = 3, col = 'red', pch = 1)


### Optional
mean.SD <- sapply(names(wx.count), FUN = function(x) mean(SD.adj[[x]]))
mean.PCC.in <- sapply(names(wx.count), FUN = function(x) mean(cor.in.adj[[x]]))
mean.PCC.out <- sapply(names(wx.count), FUN = function(x) mean(cor.out.adj[[x]]))
DNB.score <- mean.SD * mean.PCC.in / mean.PCC.out

plot(DNB.score,
     xlab = '', ylab = '', main = 'CI',
     type = 'b', 
     lwd = 3, col = 'red', 
     xaxt = "n", axes = TRUE)
axis(1, at=1:5, labels=names(wx.count))
