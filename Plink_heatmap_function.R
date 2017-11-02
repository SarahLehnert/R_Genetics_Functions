#Function to go from genepop file to heatmap of LD (r2) values
#Function converts genepop to ped/map and runs plink to get R2 matrix (LD matrix)
#Matrix is then used to create a heatmap of LD values

Plink_heatmap <- function(genepop, where.plink, where.PGDspider) {
  
  
  #Using the genepop files  - run through plink
  #Create spid file
  spidTop <- "# spid-file generated: Thu Mar 10 09:40:22 AST 2016
  # GENEPOP Parser questions
  PARSER_FORMAT=GENEPOP
  # Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
  GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
  # Select the type of the data:
  GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
  # How are Microsat alleles coded?
  GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS
  # PED Writer questions
  WRITER_FORMAT=PED
  # Save MAP file"
  map.loc <- paste0("PED_WRITER_MAP_FILE_QUESTION= ", "PGDtest")
  spidBottom <- "# Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
  PED_WRITER_ZERO_CHAR_QUESTION=
  # Specify the locus/locus combination you want to write to the PED file:
  PED_WRITER_LOCUS_COMBINATION_QUESTION=
  # Do you want to save an additional MAP file with loci information?
  PED_WRITER_MAP_QUESTION=true"
  spid.file <- c(spidTop, map.loc, spidBottom)
  
  
  where.PLINK<- where.plink
  where.PGDspider<- where.PGDspider
  
  ###write spid file genepop conversion 
  write(x = spid.file, file = paste0(where.PGDspider, "/", "hyb.spid"))
  
  file.copy(from = genepop, to = paste0(where.PGDspider, "/", "genepop_test.txt"), overwrite = T)
  
  ###convert Genepop to PED using PGD spider
  where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)
  
  
  #### OS X and LINUX CALL
  
  if(Sys.info()["sysname"] != "Windows"){
    
    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_test.txt"
    execute.SPIDER <- paste0("java -Xmx1024", "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    
    goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
    ### run PGDspider through system
    system(run.PGDspider)
    
  } # End MAC LINUX IF
  
  #### Windows call
  
  if(Sys.info()["sysname"] == "Windows"){
    
    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_test.txt"
    execute.SPIDER <- paste0("java -Xmx1024", "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
    ### run PGDspider through system
    shell(run.PGDspider)
    
  } # End WINDOWS IF
  
  
  remember.PEDpath.PGD <- paste0(where.PGDspider, "PGDtest.ped")
  remember.MAPpath.PGD <- paste0(where.PGDspider, "PGDtest.map")
  
  
  ###move to plink folder
  ped.path <- paste0(where.PGDspider, "/", "PGDtest.ped")
  map.path <- paste0(where.PGDspider, "/", "PGDtest.map")
  file.copy(from = ped.path, to = where.PLINK, overwrite = TRUE)
  
  file.copy(from = map.path, to = where.PLINK, overwrite = TRUE)
  
  
  plink_ped_path <- paste0(where.PLINK, "/", "PGDtest.ped")
  plink_map_path <- paste0(where.PLINK, "/", "PGDtest.map")
  
  ##Run plink for LD
  ####Plink creates file in plink.ld with LD values. Rename this file and move it to folder to have results later (will continue to over write plink.ld for each time plink is run)
  
  writeLines("Running plink for LD ")
  writeLines("
             ")
  
  where.PLINK.go <- gsub(x = where.PLINK, pattern = " ", replacement = "\\ ", fixed = TRUE)
  go.to.PLINK <- paste0("cd ", where.PLINK.go)    
  
  
  ### OSX LINUX PLINK call
  if(Sys.info()["sysname"] != "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file PGDtest --r2 square")
    ### run PLINK through system
    system(execute.PLINK)
  }
  
  ### Windows PLINK CALL
  if(Sys.info()["sysname"] == "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file PGDtest --r2 square")
    ### run PLINK through system
    shell(execute.PLINK)
  }
  
  writeLines("Reading R2 matrix")
  writeLines("
             ")
  
  
  #Read plink results Pop2 (LD matrix)
  LD_matrix=data.table::fread(paste(where.plink, "plink.ld", sep="/"), header=F)
  LD_matrix=as.matrix(LD_matrix)
  
  
  writeLines("Creating Heatmap from R2 matrix")
  writeLines("
             ")
  
  
  #colour palette
  mypalette4<-colorRampPalette(c("gray95","dodgerblue", "blue","midnightblue", "black"))
  
  
  
  pdf(file="LD_matrix.pdf", width = 13, height=12, bg = "white")
  gplots::heatmap.2(LD_matrix,
                    Rowv=FALSE, #rows should be reordered as required
                    Colv = "Rowv", #columns should be treated as rows
                    dendrogram="none", #no trees
                    scale="none",
                    breaks=100, #number of break points used to bin into colours
                    col=mypalette4,
                    trace="none", #whether lines should be drawn between cols, rows,
                    margins=c(5,5),#margins for column names and row names
                    labRow= " ",
                    labCol= " ",
                    key=TRUE,
                    keysize = 1,
                    density.info="none"
                    #lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
                    #lhei = c(20,5), #row height for plot elements
                    #lwid = c(8,30)  #column width for plot elements)
  )
  dev.off()
}
