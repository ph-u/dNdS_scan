#!/bin/env Rscript
# author: ph-u
# script: string_network_Coloring.r
# desc: recolor string network
# in: Rscript string_network_Coloring.r
# out: res/string_network_Coloring(V).*
# arg: 0
# date: 20250323, 20250716

source("p_src.r")
library("WeightedTreemaps")
cBp0 = cBpOri[-c(1,9)]
cons.df = read.csv(paste0(pT[1],"p_conservationTables.csv"), header = T)
cons.df$gNam = GTF$gNam[match(cons.df$gene, GTF$Locus.Tag)]
cons.df$stringAnn = string.gNam$annotation[match(cons.df$gene, string.gNam$queryItem)]
cons.df$annotation = paste(cons.df$pRod,cons.df$stringAnn)
cons.df$pRod = cons.df$stringAnn = NULL
string.int = read.table(paste0(pT[1],"string_interactions_20250323.tsv"), comment.char = "", header = T, quote = "", sep = "\t")
string.map = read.table(paste0(pT[1],"string_mapping_20250323.tsv"), comment.char = "", header = T, quote = "", sep = "\t")
string.net = read.table(paste0(pT[1],"string_network_coordinates_20250323.tsv"), comment.char = "", header = T, quote = "", sep = "\t")
string.net$PAnum = read.table(text=sub("[.]PA","@PA",string.net$identifier), sep = "@")[,2]

##### k-means clusters #####
#f.kmeans = list.files(paste0(pT[1],"string_diff"),"kmeans_clusters", full.names = T)
#for(i in 1:length(f.kmeans)){
#  i0 = read.table(f.kmeans[i], comment.char = "", header = T, quote = "", sep = "\t")
#  colnames(i0)[2] = paste0("cluster.number.",i)
#  if(i>1){f.kmeansDF = merge(f.kmeansDF, i0[,c(6,2)], by = "protein.identifier")}else{f.kmeansDF = i0[,c(6,2)]}
#};rm(i, i0)
f.kmeansDF = read.csv(paste0(pT[2], "string_diff.csv"), header = T)
string.net$clusColor = f.kmeansDF$cID[match(string.net$PAnum, f.kmeansDF$PAnum)]

##### interaction line segment coordinates #####
string.int$x0 = string.net$x_position[match(string.int$node1_string_id, string.net$identifier)]
string.int$y0 = string.net$y_position[match(string.int$node1_string_id, string.net$identifier)]
string.int$x1 = string.net$x_position[match(string.int$node2_string_id, string.net$identifier)]
string.int$y1 = string.net$y_position[match(string.int$node2_string_id, string.net$identifier)]

##### plot map #####
pdf(paste0(pT[2],"string_network_Coloring.pdf"), width = 20, height = 20)
plot(x = string.net$x_position, y = string.net$y_position, pch = 20, cex = 5, col = "#00000000")
segments(x0 = string.int$x0, y0 = string.int$y0, x1=string.int$x1, y1=string.int$y1, lwd = string.int$combined_score*10, col = "#00000011")
points(x = string.net$x_position, y = string.net$y_position, pch = 20, cex = 5, col = cBp0[string.net$clusColor])
text(x = string.net$x_position, y = string.net$y_position, labels = cons.df$gNam[match(string.net$PAnum, cons.df$gene)], offset = 1)
legend("topleft", legend = unique(f.kmeansDF$cID), fill = cBp0[1:max(f.kmeansDF$cID)], ncol = 5)
invisible(dev.off())

##### Voronoi map #####
jpeg(paste0(pT[2],"string_network_ColoringV.jpeg"), width = 1500, height = 1500, res = 300)
par(mar = c(0,0,0,0)+.1)
drawTreemap(voronoiTreemap(
  data = as.data.frame(table(string.net$clusColor)),
  levels = "Var1",
  cell_size = "Freq",
  shape = "circle",
  positioning = "regular",
  seed = 20250716
), label_size = 1, label_color = "#ffffffff", border_color = "#00000000", color_palette = cBp0[1:max(string.net$clusColor, na.rm = T)])
invisible(dev.off())

##### Word distribution PCA analysis #####
#string.map$annotation = gsub("\\s+", " ",gsub("[]]"," ",gsub("[[]"," ",gsub("[']"," ",tolower(gsub("[!@#$%^&*()_+=,./<>?|`~;:]","",gsub("--","-",gsub(" ","-",cons.df$annotation[match(string.map$queryItem,cons.df$gene)]))))))))

## Add in extra data from PseudoCap Function/Pathways/GO page
#string.map$annotation = gsub("pppa-",gsub(" ","-",paste0("pppa-",tolower("a key regulatory component of T6SS in Pseudomonas aeruginosa protein dephosphorylation protein serine threonine phosphatase activity Biofilm formation Bacterial secretion system PROTEIN PHOSPHATASE 2C")),string.map$annotation)
#string.map$annotation = gsub("pfpi-",gsub(" ","-",paste0("pfpi-",tolower("single-species biofilm formation on inanimate substrate cellular response to antibiotic single-species biofilm formation bacterial-type flagellum-dependent swarming motility cellular protein metabolic process Class I glutamine amidotransferase-like Glutamine amidotransferase JCVI: DJ-1/PfpI/YhbO family deglycase protease")),string.map$annotation)
#string.map$annotation = gsub("pfpi-",gsub(" ","-",paste0("pfpi-",tolower("")),string.map$annotation)
#keyBio = c() # !!!!!

#string.annoDB = txtConcat(paste(unique(string.map$annotation), collapse = "-"))

#i0 = c();for(i in 1:length(string.annoDB)){if(length(grep(paste0(string.annoDB[i],"s"),string.annoDB))>0){i0 = c(i0,grep(paste0(string.annoDB[i],"s"),string.annoDB))}};rm(i)
#string.annoDB = string.annoDB[-i0];rm(i0)
