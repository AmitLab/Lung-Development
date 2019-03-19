# Functions for further calculations
library(zoo)
library(pcaMethods)
library(reshape2)
library(Rgraphviz)
library(plyr)
library(metacell)


## Creating a list of the most significant genes for each metacell
#' @param mc mc object 
#' @param mat mat object 
#' @param good_mcs names of relevant metacells 
#' @param nms_per_mc number of top expressed genes we would like to choose from each metacell 
#' @param nms_thresh include only genes with at least one footprint score above nms_thresh
#' @param max_num maximal length for the final gene list
#' @param bad_genes genes to exclude
#' @param must_haves genes to include (no matter what)
#' @param ord ordering method for the gene list (none - for no order, hc - for hierarchial clustering or max.col for ordering by footprint score)

choose_genes_from_mc = function(mc,mat,good_mcs = colnames(mc@mc_fp),nms_per_mc = 5, nms_thresh = 5, max_num = Inf, bad_genes = c(), must_haves = c(),ord = "none") {
  lfp = log2(mc@mc_fp[,good_mcs])
  nms = unique(as.vector(apply(lfp,2,function(x){ head(names(sort(x, decreasing = T)),nms_per_mc)}))) # choosing top nms_per_mc genes for each mc
  nms = setdiff(nms, c(bad_genes, names(which(apply(lfp[nms,],1,max) < nms_thresh)))) # filter out bad genes and genes with maximum footprint below nms thresh
  nms = union(must_haves, head(names(sort(apply(lfp[nms, ], 1, max),T)), max_num - length(must_haves))) # add must_have genes
  if (ord == "hc") {
    nms = nms[ hclust(dist(cor(t(lfp[nms,]))), "ward.D2")$order]
  } else if (ord == "max.col") {
    nms = nms[ order(max.col(lfp[nms,]), rowSums(as.matrix(mat@mat[nms,])))]
  }
  nms
}


## Creating an heatmap with time point ordering
#' @param mc_id 
#' @param gset_id
#' @param mat_id
#' @param fig_fn file name for figure
#' @param mc_ord order for metacells 
#' @param tps factor object with ordinal metadata for each cell (in this case - time point) 
#' @param tp_cols color key for each level of tps 

heatmap_tp <-function(mc_id, gset_id, mat_id = mc_id,fig_fn = NULL,mc_ord = NULL,tps = NULL,tp_cols = NULL){
  mcp_heatmap_height = tgconfig::get_param("mcp_heatmap_height",package = "metacell")
  mcp_heatmap_width = tgconfig::get_param("mcp_heatmap_width",package = "metacell")
  mcp_heatmap_ideal_umi = tgconfig::get_param("mcp_heatmap_ideal_umi",package = "metacell")
  mcp_heatmap_fp_shades = colorRampPalette(tgconfig::get_param("mcp_heatmap_seq_shades",package = "metacell"))(1000)
  mcp_heatmap_text_cex = tgconfig::get_param("mcp_heatmap_text_cex",package = "metacell")
  mcp_heatmap_alt_side_text = tgconfig::get_param("mcp_heatmap_alt_side_text",package = "metacell")
  mcp_heatmap_latgene_color = tgconfig::get_param("mcp_heatmap_latgene_color",package = "metacell")
  mcp_heatmap_lwd = tgconfig::get_param("mcp_heatmap_lwd",package = "metacell")
  mcp_2d_legend_cex = tgconfig::get_param("mcell_mc2d_legend_cex",package = "metacell")
  mcp_heatmap_text_cex_2 = mcp_heatmap_text_cex + 2
  
  mc = scdb_mc(mc_id)
  if(is.null(mc)) {
    stop("undefined meta cell object " , mc_id)
  }
  n_mc = ncol(mc@mc_fp)
  gset = scdb_gset(gset_id)
  if(is.null(gset)) {
    stop("undefined gset object when trying to plot markers, id " , gset_id)
  }
  scmat = scdb_mat(mat_id)
  if(is.null(scmat)) {
    stop("undefined mat object when trying to plot markers, id " , mat_id)
  }
  if(length(intersect(names(mc@mc), colnames(scmat@mat))) != length(names(mc@mc))) {
    stop("cells in meta cell are missing from provided matrix in mc_plot_marks")
  }
  if(is.null(tps)){
    stop("undefined tps")
  }  
  
  if(is.null(tp_cols)){
    tp_cols = colorRampPalette(tgconfig::get_param("mcp_heatmap_fp_shades",package = "metacell"))(length(levels(tps)))
  }
  
  if(is.null(mcp_heatmap_ideal_umi)) {
    mcp_heatmap_ideal_umi = quantile(colSums(as.matrix(scmat@mat)), 0.25)
  }
  if(is.null(fig_fn)) {
    fig_fn = scfigs_fn(mc_id, "cells_tp_heat_marks")
  }
  if (is.null(mc_ord)) {
    mc_ord = 1:ncol(mc@mc_fp)
  }
  cell_ord = names(mc@mc)[order(order(mc_ord)[mc@mc])]
  good_marks = intersect(names(gset@gene_set), rownames(mc@mc_fp))
  if(length(good_marks) < 2) {
    stop("Could not get >=2 markers to plot marker matrix")
  }
  gene_folds = mc@mc_fp
  if(!is.null(fig_fn)) {
    if(is.null(mcp_heatmap_height)){
      mcp_heatmap_height = 30*length(good_marks) + 250
    }
    if(is.null(mcp_heatmap_width)){
      mcp_heatmap_width= max(min(2500,length(cell_ord)+200),800)
    }
    png(fig_fn, w=mcp_heatmap_width + 220,h=mcp_heatmap_height + 100)
  }
  
  layout(matrix(c(rep(1,45),2,3,rep(4,22),rep(5,23),6,6),ncol=2,byrow = FALSE),width = c(mcp_heatmap_width,220))
  top_marg=c(0,20,5,5)
  par(mar=top_marg)
  mat = as.matrix(as.matrix(scmat@mat)[good_marks, cell_ord])
  totu = colSums(as.matrix(scmat@mat)[, cell_ord])
  mat = t(t(mat)/totu)*mcp_heatmap_ideal_umi
  lus_1 = log2(1+7*mat)
  lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
  if (length(cell_ord) < mcp_heatmap_width) {
    smooth_n = 1
  }
  smooth_n = max(2,ceiling(2*length(cell_ord)/mcp_heatmap_width))
  lus_smoo = t(apply(lus, 1, function(x) rollmean(x,smooth_n, fill=0)))
  image(t(lus_smoo), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n')
  vals_range = range(lus_smoo)
  cell_x = rep(NA, time=length(cell_ord))
  names(cell_x) = names(mc@mc)
  cell_x[cell_ord] = 1:length(cell_ord)
  cl_x = tapply(cell_x, mc@mc, mean)/length(cell_ord)
  mtext(1:n_mc, side = 3, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
  cl_x_b = tapply(cell_x, mc@mc, max)/length(cell_ord)
  abline(v=cl_x_b, lwd=mcp_heatmap_lwd)
  all_marks = good_marks
  g_n = length(all_marks)
  gene_cols = rep("black", g_n)
  if(mcp_heatmap_alt_side_text) {
    odd = seq(1,g_n,2)
    even = seq(2,g_n,2)
    mtext(substr(all_marks[odd],1,8),
          at=seq(0,1,length.out=g_n)[odd],
          side=2, las=2, cex=mcp_heatmap_text_cex_2, col=gene_cols[odd])
    mtext(substr(all_marks[even],1,8),
          at=seq(0,1,length.out=g_n)[even],
          side=4, las=2, cex=mcp_heatmap_text_cex_2, col=gene_cols[even])
  } else {
    mtext(substr(all_marks,1,8),
          at=seq(0,1,length.out=g_n),
          side=2, las=2, cex=mcp_heatmap_text_cex_2, col=gene_cols)
  }
  lower_marg=c(0,20,0,5)
  par(mar=lower_marg)
  image(as.matrix(1:length(cell_ord),nrow=1), col=mc@colors[mc@mc[cell_ord]], xaxt='n', yaxt='n')
  par(mar=lower_marg)
  image(matrix(as.numeric(tps[cell_ord])), col = tp_cols, axes = F)
  right_top_marg = c(30,1,20,15)
  par(mar=right_top_marg)
  image(t(matrix(seq(vals_range[1], vals_range[2], length=1000))),col = mcp_heatmap_fp_shades, axes =F)
  mtext(as.integer(vals_range[1]),side=1,at = 0.25,las = 1,line=1, cex=mcp_heatmap_text_cex_2+1)
  mtext(as.integer(vals_range[2]),side=3,at = 0.25,las = 1,line=1, cex=mcp_heatmap_text_cex_2+1)
  mtext("log normalized UMI",side = 4, at = 0.5, cex=mcp_heatmap_text_cex_2 + 1, line = 3)
  right_buttom_marg = c(1,0,20,1)
  par(mar=right_buttom_marg)
  plot.new()
  legend("bottomleft", legend = levels(tps), pch = 22, cex = mcp_2d_legend_cex + 3,pt.cex = mcp_2d_legend_cex + 4, pt.bg = tp_cols,col="black", bty = "n",y.intersp = 2)
  par(mar=c(0,0,0,0))
  plot.new()
  dev.off()
}


## Summing and normalizing umi counts for each group of cells
#' @param mat_id id of matrix object
#' @param comb factor with categorial metadata (annotations) for each cell
#' @param bad_genes genes to exclude
#' @param cells well id's (cell names)
#' @param min_comb minimal size of group (minimal number of cells in a group) otherwise we will assign NA to the group
#' @param choose_genes boolean - if T we choose only differential genes by chi square test, else - we then choose all genes
#' @param normalize boolean - if T we normalize by cell size and multiply by 1000, else - we use the original matrix
#' @param g1 first group of cells (we use it only if comb is NULL)
#' @param g2 second group of cells (we use it only if comb is NULL)

sc_to_bulk = function(mat_id, comb=NULL, bad_genes = c(), cells = names(comb), min_comb = 0, choose_genes = T, normalize = T,
                      g1 = NULL, g2 = NULL) {
  
  if (is.null(comb)) { cells = union(g1,g2)}
  mat = scdb_mat(mat_id)
  umis = as.matrix(mat@mat)[,cells]
  umis_n = sweep(umis,2,colSums(umis),"/") * 1000
  if (!normalize) { umis_n = umis}
  if (is.null(comb)) {
    comb = 1 + (cells %in% g1)
    names(comb) = cells
  }
  MAP = as.numeric(factor(comb[cells])); names(MAP) = cells
  if (choose_genes) {
    genes = setdiff(scr_chi_square_diff_genes(umis[,cells], MAP = MAP[cells], fdr = T, pval = 1e-3), bad_genes)
  } else {
    genes = rownames(umis)
  }
  m = t(apply(umis_n[genes, cells],1,tapply,comb[cells],sum))
  sizes = table(comb[cells])
  good_clusts = names(which(sizes >= min_comb))
  m = m[,good_clusts]; sizes = sizes[good_clusts]
  m = sweep(m,2,as.vector(sizes),"/") * min(sizes)
  m
}


## Chi square test for differential gene expression
#' @param umis umi count table (normalized or not)
#' @param MAP a factor with cells mapping to groups
#' @param g1 first group of cells (we use it only if MAP is NULL)
#' @param g2 second group of cells (we use it only if MAP is NULL)
#' @param pval p.value 
#' @param should we use FDR correction?

scr_chi_square_diff_genes = function(umis, MAP = NULL, g1 = NULL, g2 = NULL, pval, fdr = F) {
  if (is.null(MAP)) {
    MAP = c(rep(1,length(g1)), rep(2, length(g2)))
    names(MAP) = c(g1, g2)
  }
  cells = names(MAP)
  umis = umis[,cells]
  uniform_a = rowSums(umis)/sum(umis)
  exp_count = matrix(uniform_a, ncol = 1) %*% matrix(colSums(umis),1) # exp_counts per cell
  dimnames(exp_count)  = dimnames(umis)
  ex = t(daply(.data= data.frame(cbind(V1 = MAP, t(exp_count)), check.names = F), .(V1), colSums))[-1,]
  obs = t(daply(.data= data.frame(cbind(V1 = MAP, t(umis)), check.names = F), .(V1), colSums))[-1,]
  x2 = rowSums(((obs-ex)^2 )/ex ) # chi^2 with df = ncluster-1
  if (!fdr) {
    sig_genes = x2 > qchisq(1-pval,df= length(unique(MAP)) - 1)
  } else {
    pvals = p.adjust(1 - pchisq(x2, df = length(unique(MAP)) - 1), "fdr")
    sig_genes = pvals < pval
  }
  sig_genes[ is.na(sig_genes)] = F
  return (names(sig_genes)[sig_genes])
  
}


