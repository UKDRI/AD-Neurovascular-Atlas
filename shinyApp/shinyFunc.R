library(Matrix) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(ggpubr) 
library(grid) 
library(gridExtra) 
library(tools) 



##### Functions for aesthetic control #####
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}
ggsav <- function(filename, ...){  
  if(file_ext(filename) == "pdf"){  
    ggsave(filename = filename, useDingbats = FALSE, ...)  
  } else {  
    ggsave(filename = filename, ...)  
  }
}

# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 



##### Common plotting functions #####
# Plot cell information / gene expression on dimred 
sc2Ddimr <- function(inpConf, inpMeta, inpDimr, inpdr, inp1, 
                     inpH5, inpGene, inpDtyp, inpsub1, inpsub2, inpmin, inpmax, 
                     inpsiz, inpord, inpcol, inpfsz, inpasp, inptxt, inplab){ 
  # Tidy inputs 
  inpDtyp = gsub("^Assay: ", "", inpDtyp) 
  inpH5 = paste0(inpH5, inpDtyp, ".h5") 
  inpGene = inpGene[[inpDtyp]] 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = data.table(inpDimr[[inpdr]]); colnames(ggData) = c("X","Y")
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  ggData$sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE]
  if(inpDtyp == "Cell Information"){
    ggData$val = inpMeta[, inpConf[UI == inp1]$ID, with = FALSE]
  } else {
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    # ggData[val < 0]$val = 0 
    h5file$close_all() 
  }
  # Split ggData into foreground / background
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  # Reorder the cells by value
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  # Do factoring if required 
  ggCont = TRUE    # Mark is plotted value is continuous
  if(inpDtyp == "Cell Information"){
    if(!is.na(inpConf[UI == inp1]$fCL)){ 
      ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
      names(ggCol) = levels(ggData$val) 
      ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
      ggData$val = factor(ggData$val, levels = ggLvl) 
      ggCol = ggCol[ggLvl] 
      ggCont = FALSE
    } 
  }
  # Apply min/max.cutoff if required 
  if(ggCont){ 
    min.cutoff = quantile(ggData$val, inpmin/100) 
    max.cutoff = quantile(ggData$val, inpmax/100) 
    ggData[val <= min.cutoff]$val = min.cutoff 
    ggData[val >= max.cutoff]$val = max.cutoff 
  }
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + 
    xlab(paste0(inpdr,"1")) + ylab(paste0(inpdr,"2")) + 
    sctheme(base_size = inpfsz, XYval = inptxt) 
  if(ggCont){ 
    ggOut = ggOut + scale_color_gradientn(paste0(inp1, "  "), colours = inpcol) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    legfsz = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    legfsz = 0.75 * (inpfsz - (1.5 * floor(legfsz/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = legfsz)) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = (inpfsz/4), seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

# Function to get cell numbers 
sc2Dnum <- function(inpConf, inpMeta, inpDimr, inpdr, inpdrX, inpdrY, inp1, 
                    inpH5, inpGene, inpDtyp, inpsub1, inpsub2, inpsplt){ 
  # Tidy inputs 
  inpDtyp = gsub("^Assay: ", "", inpDtyp) 
  inpH5 = paste0(inpH5, inpDtyp, ".h5") 
  inpGene = inpGene[[inpDtyp]] 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = data.table(inpDimr[[inpdr]]); colnames(ggData) = c("X","Y") 
  ggData$sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE] 
  if(inpDtyp == "Cell Information"){ 
    ggData$group = inpMeta[, inpConf[UI == inp1]$ID, with = FALSE] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$group = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    h5file$close_all() 
  } 
  # Split ggData into foreground / background 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  # Split inp1 if necessary 
  if(length(inpConf[UI == inp1]$fCL) == 0){ 
    ggData$group = cut(ggData$group, breaks = round(inpsplt, 0)) 
  } else if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggData$group = cut(ggData$group, breaks = round(inpsplt, 0)) 
  } 
  # Actual data.table 
  ggData1 = ggData
  if(!is.null(inpdrX[1])){
    ggData1 = ggData1[X > inpdrX[1]][X < inpdrX[2]]
    ggData1 = ggData1[Y > inpdrY[1]][Y < inpdrY[2]]
  } 
  ggData1 = ggData1[, .(nZoom = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nZoom"), with = FALSE] 
  ggData[is.na(nZoom)]$nZoom = 0 
  ggData$pctZoom = 100 * ggData$nZoom / ggData$nCells 
  ggData = ggData[order(group)] 
  return(ggData) 
} 

# Plot X/Y relationship (either confusion matrix or scatter plot) 
sc2Dcomp <- function(inpConf, inpMeta, inp1, inp2, inpH5, inpGene,  
                     inpDtyp1, inpDtyp2, inpsub1, inpsub2,  
                     inpmin1, inpmax1, inpmin2, inpmax2, 
                     inpsiz, inpcol, inpfsz){  
  # Tidy inputs  
  inpDtyp1 = gsub("^Assay: ", "", inpDtyp1)  
  inpDtyp2 = gsub("^Assay: ", "", inpDtyp2)  
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}  
  # Prepare ggData  
  ggData = data.table(sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE]) 
  if(inpDtyp1 == "Cell Information"){ 
    ggData$val1 = inpMeta[, inpConf[UI == inp1]$ID, with = FALSE] 
  } else { 
    h5file <- H5File$new(paste0(inpH5, inpDtyp1, ".h5") , mode = "r")  
    h5data <- h5file[["grp"]][["data"]]  
    ggData$val1 = h5data$read(args = list(inpGene[[inpDtyp1]][inp1], quote(expr=)))  
    h5file$close_all()  
  } 
  if(inpDtyp2 == "Cell Information"){ 
    ggData$val2 = inpMeta[, inpConf[UI == inp2]$ID, with = FALSE] 
  } else { 
    h5file <- H5File$new(paste0(inpH5, inpDtyp2, ".h5") , mode = "r")  
    h5data <- h5file[["grp"]][["data"]]  
    ggData$val2 = h5data$read(args = list(inpGene[[inpDtyp2]][inp2], quote(expr=)))  
    h5file$close_all()  
  } 
  # Split ggData into foreground 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){  
    ggData = ggData[sub %in% inpsub2]  
  }  
  # Identify continuous covariates and Apply min/max.cutoff if required 
  ggCont1 = TRUE    # Mark is plotted value is continuous 
  if(inpDtyp1 == "Cell Information"){ 
    if(!is.na(inpConf[UI == inp1]$fCL)){  
      ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]  
      names(ggCol) = levels(ggData$val1)  
      ggLvl = levels(ggData$val1)[levels(ggData$val1) %in% unique(ggData$val1)]  
      ggData$val1 = factor(ggData$val1, levels = ggLvl)  
      # ggCol = ggCol[ggLvl]  
      ggCont1 = FALSE 
    }  
  } 
  if(ggCont1){  
    min.cutoff1 = quantile(ggData$val1, inpmin1/100) 
    max.cutoff1 = quantile(ggData$val1, inpmax1/100) 
    ggData[val1 <= min.cutoff1]$val1 = min.cutoff1 
    ggData[val1 >= max.cutoff1]$val1 = max.cutoff1 
  } 
  # Repeat for inp2 
  ggCont2 = TRUE    # Mark is plotted value is continuous 
  if(inpDtyp2 == "Cell Information"){ 
    if(!is.na(inpConf[UI == inp2]$fCL)){  
      ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]  
      names(ggCol) = levels(ggData$val2)  
      ggLvl = levels(ggData$val2)[levels(ggData$val2) %in% unique(ggData$val2)]  
      ggData$val2 = factor(ggData$val2, levels = ggLvl)  
      # ggCol = ggCol[ggLvl]  
      ggCont2 = FALSE 
    }  
  } 
  if(ggCont2){  
    min.cutoff2 = quantile(ggData$val2, inpmin2/100) 
    max.cutoff2 = quantile(ggData$val2, inpmax2/100) 
    ggData[val2 <= min.cutoff2]$val2 = min.cutoff2 
    ggData[val2 >= max.cutoff2]$val2 = max.cutoff2 
  } 
   
  # Start ggplot 
  if(ggCont1 == TRUE & ggCont2 == TRUE){ 
    # Both cont: Scatter plot 
    ggOut <- ggplot(ggData, aes(val1, val2)) +  
      geom_point() + sctheme(base_size = inpfsz) + xlab(inp1) + ylab(inp2) 
  } else if(ggCont1 == FALSE & ggCont2 == FALSE){ 
    # Both cate: Confusion matrix 
    allCombi <- CJ(val1 = unique(ggData$val1), val2 = unique(ggData$val2)) 
    ggData <- ggData[, .N, by = .(val1, val2)] 
    ggData <- merge(allCombi, ggData, by = c("val1", "val2"), all.x = TRUE) 
    ggData[is.na(N), N := 0]  # Replace NA with 0  
    ggOut <- ggplot(ggData, aes(val1, val2, fill = N, label = N)) + 
      geom_tile(color = "white") + geom_text(size = inpfsz/5) +  
      scale_fill_gradient(low = "white", high = "steelblue") + 
      theme_minimal(base_size = inpfsz) + xlab(inp1) + ylab(inp2) +  
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
      ggtitle("Confusion matrix") 
  } else { 
    # Cate-cont: Ask user to use VlnBoxP directly! 
    ggOut <- ggplot() + sctheme(base_size = inpfsz) +  
      labs(title = paste0("Please use violin / boxplot to interrogate the data! \n ", 
                          "Or refer to `Cell numbers / statistics` tab on the right...")) 
  } 
  return(ggOut) 
} 

# Plot X/Y relationship (either confusion matrix or scatter plot) 
sc2Dcnum <- function(inpConf, inpMeta, inp1, inp2, inpH5, inpGene,  
                     inpDtyp1, inpDtyp2, inpsub1, inpsub2,  
                     inpmin1, inpmax1, inpmin2, inpmax2, inpcut = 0){  
  # Tidy inputs  
  inpDtyp1 = gsub("^Assay: ", "", inpDtyp1)  
  inpDtyp2 = gsub("^Assay: ", "", inpDtyp2)  
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}  
  # Prepare ggData  
  ggData = data.table(sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE]) 
  if(inpDtyp1 == "Cell Information"){ 
    ggData$val1 = inpMeta[, inpConf[UI == inp1]$ID, with = FALSE] 
  } else { 
    h5file <- H5File$new(paste0(inpH5, inpDtyp1, ".h5") , mode = "r")  
    h5data <- h5file[["grp"]][["data"]]  
    ggData$val1 = h5data$read(args = list(inpGene[[inpDtyp1]][inp1], quote(expr=)))  
    h5file$close_all()  
  } 
  if(inpDtyp2 == "Cell Information"){ 
    ggData$val2 = inpMeta[, inpConf[UI == inp2]$ID, with = FALSE] 
  } else { 
    h5file <- H5File$new(paste0(inpH5, inpDtyp2, ".h5") , mode = "r")  
    h5data <- h5file[["grp"]][["data"]]  
    ggData$val2 = h5data$read(args = list(inpGene[[inpDtyp2]][inp2], quote(expr=)))  
    h5file$close_all()  
  } 
  # Split ggData into foreground 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){  
    ggData = ggData[sub %in% inpsub2]  
  }  
  # Identify continuous covariates and Apply min/max.cutoff if required 
  ggCont1 = TRUE    # Mark is plotted value is continuous 
  if(inpDtyp1 == "Cell Information"){ 
    if(!is.na(inpConf[UI == inp1]$fCL)){  
      ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]  
      names(ggCol) = levels(ggData$val1)  
      ggLvl = levels(ggData$val1)[levels(ggData$val1) %in% unique(ggData$val1)]  
      ggData$val1 = factor(ggData$val1, levels = ggLvl)  
      # ggCol = ggCol[ggLvl]  
      ggCont1 = FALSE 
    }  
  } 
  if(ggCont1){  
    min.cutoff1 = quantile(ggData$val1, inpmin1/100) 
    max.cutoff1 = quantile(ggData$val1, inpmax1/100) 
    ggData[val1 <= min.cutoff1]$val1 = min.cutoff1 
    ggData[val1 >= max.cutoff1]$val1 = max.cutoff1 
  } 
  # Repeat for inp2 
  ggCont2 = TRUE    # Mark is plotted value is continuous 
  if(inpDtyp2 == "Cell Information"){ 
    if(!is.na(inpConf[UI == inp2]$fCL)){  
      ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]  
      names(ggCol) = levels(ggData$val2)  
      ggLvl = levels(ggData$val2)[levels(ggData$val2) %in% unique(ggData$val2)]  
      ggData$val2 = factor(ggData$val2, levels = ggLvl)  
      # ggCol = ggCol[ggLvl]  
      ggCont2 = FALSE 
    }  
  } 
  if(ggCont2){  
    min.cutoff2 = quantile(ggData$val2, inpmin2/100) 
    max.cutoff2 = quantile(ggData$val2, inpmax2/100) 
    ggData[val2 <= min.cutoff2]$val2 = min.cutoff2 
    ggData[val2 >= max.cutoff2]$val2 = max.cutoff2 
  } 
  
  # Actual data.table 
  if(ggCont1 == TRUE & ggCont2 == TRUE){ 
    # Both cont: Pearson / Spearman corr 
    ggP <- round(cor(ggData$val1, ggData$val2, method = "pearson"), 6) 
    ggS <- round(cor(ggData$val1, ggData$val2, method = "spearman"), 6) 
    ggK <- round(cor(ggData$val1, ggData$val2, method = "kendall"), 6) 
    ggModel <- lm(val2 ~ val1, data = ggData)
    ggCoeff <- signif(coef(ggModel), 6)
    ggData = data.table(measures = c("Pearson Corr.",
                                     "Spearman Corr.",
                                     "Kendall Corr.",
                                     "LM (Y/RHS~X/LHS) Gradient",
                                     "LM (Y/RHS~X/LHS) Intercept"),
                        value = c(ggP, ggS, ggK, ggCoeff[2], ggCoeff[1]))
    
  } else if(ggCont1 == FALSE & ggCont2 == FALSE){ 
    # Both cate: Confusion matrix 
    ggData <- ggData[, .N, by = .(val1, val2)][order(val1, val2)]
    colnames(ggData) <- c("value_left", "value_right", "nCells")
    
  } else { 
    # Cate-cont: Ask user to use VlnBoxP directly! 
    if(ggCont1 == TRUE){
      colnames(ggData) <- c("sub", "val2", "val1")
    } 
    ggData$express = FALSE 
    ggData[val2 > inpcut]$express = TRUE 
    ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "val1"] 
    ggData = ggData[, .(nCells = .N), by = "val1"] 
    ggData = ggData1[ggData, on = "val1"] 
    ggData = ggData[, c("val1", "nCells", "nExpress"), with = FALSE] 
    ggData[is.na(nExpress)]$nExpress = 0 
    ggData$pctExpress = round(100 * ggData$nExpress / ggData$nCells, 2)
    colnames(ggData)[1] = "group"
    colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
    ggData = ggData[order(group)] 
  } 
  return(ggData) 
} 

# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpDimr, inpdr, inp1, inp2, 
                     inpH5, inpGene, inpsub1, inpsub2, 
                     inpmin1, inpmax1, inpmin2, inpmax2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = data.table(inpDimr[[inpdr]]); colnames(ggData) = c("X","Y") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  ggData$sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE] 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  # Apply min/max.cutoff if required 
  min.cutoff1 = quantile(ggData$val1, inpmin1/100) 
  max.cutoff1 = quantile(ggData$val1, inpmax1/100) 
  ggData[val1 <= min.cutoff1]$val1 = min.cutoff1 
  ggData[val1 >= max.cutoff1]$val1 = max.cutoff1 
  min.cutoff2 = quantile(ggData$val2, inpmin2/100) 
  max.cutoff2 = quantile(ggData$val2, inpmax2/100) 
  ggData[val2 <= min.cutoff1]$val2 = min.cutoff2 
  ggData[val2 >= max.cutoff1]$val2 = max.cutoff2 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(paste0(inpdr,"1")) + ylab(paste0(inpdr,"2")) + 
    sctheme(base_size = inpfsz, XYval = inptxt) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = inpfsz, XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpH5, inpGene, inpsub1, inpsub2){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 

# Plot violin / boxplot 
filter_list <- function(item, fullList, reqLen = 2) { 
  # Filter the item based on the master list 
  filtered_item <- item[item %in% fullList] 
  # Return the item only if its length is 2 
  if (length(filtered_item) == reqLen) { 
    return(filtered_item) 
  } else { 
   return(NULL) # Return NULL if the condition is not met 
  } 
} 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, inpH5, inpGene, inpDtyp, 
                     inpsub1, inpsub2, inptyp, inppts, inpstg, inpstp1, inpstp2, 
                     inpsiz, inpfsz, inpnoi){ 
  # Tidy inputs 
  inpDtyp = gsub("^Assay: ", "", inpDtyp) 
  inpH5 = paste0(inpH5, inpDtyp, ".h5") 
  inpGene = inpGene[[inpDtyp]] 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inpDtyp == "Cell Information"){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    if(inpnoi){ 
      set.seed(42) 
      tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
      ggData$val = ggData$val + tmpNoise 
    } 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  # Add stats test if specified (stg: stats-global, stp: stats-pairwise) 
  if(inpstg != "none"){ 
    ggOut = ggOut + stat_compare_means(size = inpfsz/3, label.x.npc = "centre", method = inpstg) 
  } 
  if(inpstp1 != "none"){ 
    compList = strsplit(gsub("\"|'| ", "", inpstp2), "\n")[[1]] 
    compList = lapply(compList, function(x) { 
      strsplit(x, ",|;")[[1]] 
    }) 
    compList = lapply(compList, filter_list, fullList = ggData$X) 
    compList = compList[!sapply(compList, is.null)] 
    if(length(compList) > 0){ 
      ggOut = ggOut + stat_compare_means(size = inpfsz/3, method = inpstp1, comparisons = compList) 
    } 
  } 
  # Continue ggplot 
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = inpfsz, Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut) 
} 

# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inpord1, inpord2, inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Reorder X-axis if required 
  if(inpord1 != "Original order"){ 
    if(inptyp == "Proportion"){ 
      ggOrder = ggData[grp == inpord1][order(pctCells)] 
      if(inpord2 == "Decreasing"){ggOrder = ggOrder[order(-pctCells)]} 
    } else { 
      ggOrder = ggData[grp == inpord1][order(nCells)] 
      if(inpord2 == "Decreasing"){ggOrder = ggOrder[order(-nCells)]} 
    }  
    if(inpord2 == "Decreasing"){ 
      finalOrder = c(as.character(ggOrder$X), 
                     setdiff(levels(ggOrder$X), as.character(ggOrder$X))) 
    } else { 
      finalOrder = c(setdiff(levels(ggOrder$X), as.character(ggOrder$X)), 
                     as.character(ggOrder$X)) 
    } 
    ggData$X = factor(ggData$X, levels = finalOrder) 
    ggLvl = c(setdiff(ggLvl, inpord1), inpord1) 
    ggData$grp = factor(ggData$grp, levels = ggLvl) 
  } 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = inpfsz, Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 

# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ',|;|"|\n')[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 

# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpH5, inpGene, inpsub1, inpsub2, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, inpExp, inpMax, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
   
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("cellID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("cellID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
   
  # Aggregate 
  if(inpExp){ggData$val = expm1(ggData$val)} 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(cellID)), 
                  by = c("geneName", "grpBy")] 
  if(inpExp){ggData$val = log1p(ggData$val)} 
   
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    ggData[val >  inpMax]$val =  inpMax 
    ggData[val < -inpMax]$val = -inpMax 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
   
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = inpfsz) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = inpfsz, Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
   
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = inpfsz, Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = inpcols) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = inpfsz, Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = inpcols) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
     
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 

# Plot cell information / gene expression on spatial 
sc2Dspat <- function(inpConf, inpMeta, inpImg, inpAlp, inp1, 
                     inpH5, inpGene, inpDtyp, inpsub1, inpsub2, inpmin, inpmax, 
                     inpsiz, inpcol, inpfsz, inplab){ 
  # Tidy inputs 
  inpDtyp = gsub("^Assay: ", "", inpDtyp) 
  inpH5 = paste0(inpH5, inpDtyp, ".h5") 
  inpGene = inpGene[[inpDtyp]] 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Setup image
  orgXY = c(0, dim(inpImg[["bg_image"]])[1], 0, dim(inpImg[["bg_image"]])[2]) 
  orgXY = orgXY / inpImg$lowres 
  inpImgGrob <- rasterGrob(inpImg[["bg_image"]], interpolate = FALSE, 
                           width=unit(1,"npc"), height=unit(1,"npc")) 
  # Prepare ggData 
  ggData = data.table(inpImg$coord[, c("imagecol", "imagerow")])
  colnames(ggData) = c("X","Y")
  ggData$Y = orgXY[4] - ggData$Y
  ggData$sub = inpMeta[, inpConf[UI == inpsub1]$ID, with = FALSE]
  ggData$siz = 1
  if(inpDtyp == "Cell Information"){
    ggData$val = inpMeta[, inpConf[UI == inp1]$ID, with = FALSE]
  } else {
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    ggData[val < 0]$val = 0 
    h5file$close_all() 
  }
  # Split ggData into foreground / background
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  # Do factoring if required 
  ggCont = TRUE    # Mark is plotted value is continuous
  if(inpDtyp == "Cell Information"){
    if(!is.na(inpConf[UI == inp1]$fCL)){ 
      ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
      names(ggCol) = levels(ggData$val) 
      ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
      ggData$val = factor(ggData$val, levels = ggLvl) 
      ggCol = ggCol[ggLvl] 
      ggCont = FALSE
    } 
  }
  # Apply min/max.cutoff if required 
  if(ggCont){ 
    min.cutoff = quantile(ggData$val, inpmin/100) 
    max.cutoff = quantile(ggData$val, inpmax/100) 
    ggData[val <= min.cutoff]$val = min.cutoff 
    ggData[val >= max.cutoff]$val = max.cutoff 
  }
  
  # Actual ggplot 
  if(ggCont){
    ggOut = ggplot(ggData, aes(X, Y, color = val, size = siz, alpha = val)) +
      annotation_custom(inpImgGrob, xmin = orgXY[1], xmax = orgXY[2], 
                                    ymin = orgXY[3], ymax = orgXY[4])
  } else {
    ggOut = ggplot(ggData, aes(X, Y, color = val, size = siz)) +
      annotation_custom(inpImgGrob, xmin = orgXY[1], xmax = orgXY[2], 
                                    ymin = orgXY[3], ymax = orgXY[4])
  }
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16, alpha = inpAlp[1]) 
  } 
  if(ggCont){ 
    ggOut = ggOut + 
      geom_point(shape = 16) + 
      scale_color_gradientn(paste0(inp1, "  "), colours = inpcol) + 
      scale_size_continuous(range = c(0, inpsiz), limits = c(0,1)) + 
      guides(color = guide_colorbar(barwidth = 15),
             size = "none", alpha = "none") + scale_alpha(range = inpAlp)
  } else { 
    legfsz = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    legfsz = 0.75 * (inpfsz - (1.5 * floor(legfsz/50))) 
    ggOut = ggOut + 
      geom_point(shape = 16, alpha = inpAlp[2]) + 
      scale_color_manual("", values = ggCol) + 
      scale_size_continuous(range = c(0, inpsiz), limits = c(0,1)) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow),
             size = "none") + 
      theme(legend.text = element_text(size = legfsz)) 
    if(inplab){
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = (inpfsz/4), seed = 42) 
    } 
  } 
  ggOut = ggOut + coord_fixed() + 
    theme_void(base_size = inpfsz) + theme(
      axis.text = element_blank(), axis.line = element_blank(), 
      legend.position = "bottom", legend.key = element_rect(colour = NA, fill = NA))
  return(ggOut) 
} 


