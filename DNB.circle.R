# 

rm(list = ls())
setwd(dir = '/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/')

# Load Expresion Value
#1-- This scirpt is used for './0.Backup/[!]fpkm/FPKM_after.txt' ---#
Gene.count <- read.csv(file = './0.Backup/[!]fpkm/FPKM_after.txt', sep = '\t', row.names = 1)
rownames(Gene.count) <- sapply(rownames(Gene.count), FUN = function(x) stringr::str_to_upper(x))
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
# rownames(Gene.count) <- sapply(rownames(Gene.count), FUN = function(x) stringr::str_to_upper(x))
# wx.count <- list() # For disease
# wx.count.c <- list() # For control
# wx.count$w3 <- Gene.count[, grep(pattern = 'X3WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w5 <- Gene.count[, grep(pattern = 'X5WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w9 <- Gene.count[, grep(pattern = 'X9WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w14 <- Gene.count[, grep(pattern = 'X14WT', x = colnames(Gene.count), value = TRUE)]
# wx.count$w17 <- Gene.count[, grep(pattern = 'X17WT', x = colnames(Gene.count), value = TRUE)]
# wx.count.c$wx <- Gene.count[, grep(pattern = 'N|C', x = colnames(Gene.count), value = TRUE)]

# Load DNB
DNB <- read.table(file = './3.2.dnb/module_gene.txt', stringsAsFactors = FALSE)$V1
DNB <- sapply(DNB, FUN = function(x) stringr::str_to_upper(x))

# Construct Network
Nodes <- data.frame(Node = rownames(Gene.count), 
                    isDNB = '-', 
                    stringsAsFactors = FALSE)
Nodes$isDNB[Nodes$Node %in% DNB] <- '+'

Edges <- data.frame()
for (dnb in DNB) {
  Edges.tmp <- data.frame(node1 = dnb, 
                          node2 = rownames(Gene.count)[rownames(Gene.count) != dnb], 
                          stringsAsFactors = FALSE)
  Edges <- rbind(Edges, Edges.tmp)
}

#save(Edges, Nodes, wx.count, wx.count.c, file = './3.2.dnb/var_used.as.sd.pcc.RData')
#---- Start Here ----
load(file = './3.2.dnb/var_used.as.sd.pcc.RData')

# mapping between disease and control
maps <- rep('wx', 5); names(maps) <- names(wx.count)

#-- for week x (3, 5, 9, 14, 17)
for (w in names(wx.count)) {
  wx <- wx.count[[w]]
  Nodes[w] <- sapply(Nodes$Node, FUN = function(x) {
    sd(wx[x, ])
  })
  
  wx <- t(wx)
  tmp <- apply(Edges[, c('node1', 'node2')], MARGIN = 1, FUN = function(x) {
    cor(x = wx[, x[1]], y = wx[, x[2]], method = "spearman")
  })
  Edges[w] <- tmp
}

for (w in names(wx.count.c)) {
  wx <- wx.count.c[[w]]
  Nodes[w] <- sapply(Nodes$Node, FUN = function(x) {
    sd(wx[x, ])
  })
  
  wx <- t(wx)
  tmp <- apply(Edges[, c('node1', 'node2')], MARGIN = 1, FUN = function(x) {
    cor(x = wx[, x[1]], y = wx[, x[2]], method = "spearman")
  })
  Edges[w] <- tmp
}

for (m in 1:length(maps)) {
  Nodes[paste(names(maps[m]), 'adj', sep = '_')] <- Nodes[names(maps[m])] / (Nodes[maps[m]] + 0.001)
  
  Edges[paste(names(maps[m]), 'adj', sep = '_')] <- abs(Edges[names(maps[m])]) / (abs(Edges[maps[m]]) + 0.001)
}

#save(Edges, Nodes, file = './3.2.dnb/Final.net.RData')

#---- Start Here ----
load(file = './3.2.dnb/Final.net.RData')
write.table(Nodes, file = './3.2.dnb/Circle_node.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(Edges, file = './3.2.dnb/Circle_edge.tsv', quote = FALSE, sep = '\t', row.names = FALSE)


