
##Function to generate LD-Fst relationship plot and LD heatmap for top 500 high Fst loci for two groups/populations/morphs
##Require R packages: diveRsity, gplots, genepopedit, data.table
#Input genepop file with only two groups/populations and include path for plink and PGDspider


LD_Fst_Heatmap <- function(genepop, where.plink, where.PGDspider) {
  
  #Console message
  writeLines("Require R packages: diveRsity, gplots, genepopedit, data.table\n\n")
  
  
  #Console message
  writeLines("Calculating pairwise Fst\n\n")
  
  #generate locus specific Fst values between the two populations
  diveRsity::diffCalc(genepop, "fst", fst=T)
  
  #Read in Fst file and remove last row of Fst file (which is Global Fst)
  fst=read.table("fst-[diffCalc]/std_stats.txt", header=T)
  rows=length(fst$loci)-1
  Fst_1=fst[1:rows, ]
  
  #Sort Fst from lowest to highest value
  sort_fst <- Fst_1[order(as.numeric(as.character(Fst_1$Fst))), ]
  sort_fst2 <- sort_fst[!(is.na(sort_fst$Fst)),]
  sort_fst2$loci <- factor(sort_fst2$loci) #remove "Global" level
  
  #order of Loci to keep for genepop file (in order from lowest to highest)
  loci_keep=as.character(sort_fst2$loci)
  
  #Get population group names
  pop_groups= genepopedit::genepop_detective(genepop, "Pops")
  
  writeLines("Ordering loci from lowest to highest FST and subsetting genepop files")
  writeLines("
             ")
  
  #Subset genepop files for loci in order for all indivdiuals and within each population group
  genepopedit::subset_genepop(genepop, keep = TRUE, subs = loci_keep, path ="genepop_all_ordered_fst.txt")
  genepopedit::subset_genepop(genepop, keep = TRUE, spop = pop_groups[1], subs = loci_keep, path ="genepop_Pop1_ordered_fst.txt")
  genepopedit::subset_genepop(genepop, keep = TRUE, spop = pop_groups[2], subs = loci_keep, path ="genepop_Pop2_ordered_fst.txt")
  
  
  
  #Using the new genepop files generated - run through plink
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
  
  #Indicate where plink and PGDspider are
  where.PLINK<- where.plink
  where.PGDspider<- where.PGDspider
  
  ###write spid file genepop conversion 
  write(x = spid.file, file = paste0(where.PGDspider, "/", "hyb.spid"))
  file.copy(from = "genepop_all_ordered_fst.txt", to = paste0(where.PGDspider, "/", "genepop_all_ordered_fst.txt"), overwrite = T) 
  
  
  ###convert Genepop to PED using PGD spider
  where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)
  
  
  #### OS X and LINUX CALL
  
  if(Sys.info()["sysname"] != "Windows"){
    
    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_all_ordered_fst.txt"
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
    input.file.call <- "-inputfile genepop_all_ordered_fst.txt"
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
  
  ###move files to plink folder to run plink
  ped.path <- paste0(where.PGDspider, "/", "PGDtest.ped")
  map.path <- paste0(where.PGDspider, "/", "PGDtest.map")
  file.copy(from = ped.path, to = where.PLINK, overwrite = TRUE)
  
  file.copy(from = map.path, to = where.PLINK, overwrite = TRUE)
  
  plink_ped_path <- paste0(where.PLINK, "/", "PGDtest.ped")
  plink_map_path <- paste0(where.PLINK, "/", "PGDtest.map")
  
  ##Run plink for LD (for all groups together - LD overall individuals)
  
  writeLines("Running plink for LD")
  writeLines("
             ")
  
  ####Plink creates file in plink.ld with LD values. Read in this file (will continue to over write plink.ld for each time plink is run)
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
  
  
  
  #Console message
  
  writeLines("\n\n  Reading LD matrix \n\n ")
  
  
  ##Since probably don't want to save LD matrix file as it is very large just read it in and use it
  LD_all=data.table::fread(paste0(where.plink, "/", "plink.ld"), header=F)
  
  #Take row means for LD matrix (=mean LD across all loci)
  LD_all$LD_means=rowMeans(LD_all, na.rm = T)
  
  #Add loci names to table
  loci_names_all=genepopedit::genepop_detective("genepop_all_ordered_fst.txt", "Loci")
  LD_all$loci=loci_names_all
  
  #Create new dataframe with LD means and Loci names
  Mean_LD=cbind(LD_all$loci, LD_all$LD_means)
  colnames(Mean_LD)<-c("loci", "mean_LD")
  
  #merge mean LD values with Fst values based on loci names to plot LD vs Fst relationship
  merged=merge(fst, Mean_LD, by=1) 
  merged2=merged[,c(1,7,9)] #remove extra statistics generated by diveRsity that aren't needed now
  
  #ensure LD and Fst are read as numeric 
  merged2$mean_LD=as.numeric(as.character(merged2$mean_LD))
  merged2$Fst=as.numeric(as.character(merged2$Fst))
  
  #save PDF plot of LD-Fst relationship for the species
  pdf(file="LD_FST_relationship.pdf", width = 8, height=7, bg = "white")
  plot(merged2$Fst, merged2$mean_LD, xlab="Locus specific FST", ylab="Mean LD (overall)", las=1, pch=19)
  dev.off()
  
  #save text file with FST and LD values for later
  write.table(merged2, "LD_FST_results.txt", row.names = F, quote = F)
  
  #Now run analyses for within each population group and generate Heat Map
  
  ###Population/Group 1
  
  writeLines("Converting files from genepop to PED format")
  writeLines("
             ")
  
  file.remove(paste0(where.PLINK, "/", "PGDtest.ped"))
  file.remove(paste0(where.PLINK, "/", "PGDtest.map"))
  
  #Using the new genepop files generated - run through plink
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
  file.copy(from = "genepop_Pop1_ordered_fst.txt",to = paste0(where.PGDspider, "/", "genepop_Pop1_ordered_fst.txt"), overwrite = T)
  
  ###convert Genepop to PED using PGD spider
  where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)
  
  #### OS X and LINUX CALL
  
  if(Sys.info()["sysname"] != "Windows"){
    
    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_Pop1_ordered_fst.txt"
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
    input.file.call <- "-inputfile genepop_Pop1_ordered_fst.txt"
    execute.SPIDER <- paste0("java -Xmx1024", "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    #  ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
    # ### run PGDspider through system
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
  ####Plink creates file in plink.ld with LD values. 
  where.PLINK.go <- gsub(x = where.PLINK, pattern = " ", replacement = "\\ ", fixed = TRUE)
  go.to.PLINK <- paste0("cd ", where.PLINK.go)    
  
  writeLines("Running plink for LD for pop 1")
  writeLines("
             ") 
  
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
  
  
  
  #Console message
  writeLines("
             ")
  writeLines("Reading LD matrix for Pop 1")
  writeLines("
             ")
  
  #Read plink results. Use this LD matrix for plot below (both matrix results together for both pops)
  LD_pop1=data.table::fread(paste(where.plink, "plink.ld", sep="/"), header=F)
  LD_pop1=as.matrix(LD_pop1)
  
  
  LD_pop1_a=data.table::fread(paste(where.plink, "plink.ld", sep="/"), header=F)
  #Add loci names to table
  
  loci_names_all=genepopedit::genepop_detective("genepop_Pop1_ordered_fst.txt", "Loci")
  LD_pop1_a$LD_means=rowMeans(LD_pop1_a, na.rm = T)
  LD_pop1_a$loci=loci_names_all
  
  #Create new dataframe with LD means and Loci names
  Mean_LD=cbind(LD_pop1_a$loci, LD_pop1_a$LD_means)
  colnames(Mean_LD)<-c("loci", "mean_LD")
  
  #merge mean LD values with Fst values based on loci names to plot LD vs Fst relationship
  merged=merge(fst, Mean_LD, by=1) 
  merged2=merged[,c(1,7,9)] #remove extra statistics generated by diveRsity that aren't needed now
  
  #ensure LD and Fst are read as numeric 
  merged2$mean_LD=as.numeric(as.character(merged2$mean_LD))
  merged2$Fst=as.numeric(as.character(merged2$Fst))
  
  #save PDF plot of LD-Fst relationship for the species
  pdf(file="LD_FST_relationship_pop1.pdf", width = 8, height=7, bg = "white")
  plot(merged2$Fst, merged2$mean_LD, xlab="Locus specific FST", ylab="Mean LD (within pop 1)", las=1, pch=19)
  dev.off()
  
  #Now run analyses for within each population group and generate Heat Map
  
  ##Population/group 2
  
  file.remove(paste0(where.PLINK, "/", "PGDtest.ped"))
  file.remove(paste0(where.PLINK, "/", "PGDtest.map"))
  
  writeLines("Converting files from genepop to PED format")
  writeLines("
             ")
  
  
  #Using the new genepop files generated - run through plink
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
  
  file.copy(from = "genepop_Pop2_ordered_fst.txt",to = paste0(where.PGDspider, "/", "genepop_Pop2_ordered_fst.txt"), overwrite = T)
  
  ###convert Genepop to PED using PGD spider
  where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)
  
  
  #### OS X and LINUX CALL
  
  if(Sys.info()["sysname"] != "Windows"){
    
    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_Pop2_ordered_fst.txt"
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
    input.file.call <- "-inputfile genepop_Pop2_ordered_fst.txt"
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
  
  writeLines("Running plink for LD for pop 2")
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
  
  
  #Console message
  writeLines("
             ")
  writeLines("Reading LD matrix for Pop 2")
  writeLines("
             ")
  
  #Read plink results Pop2 (LD matrix)
  LD_pop2=data.table::fread(paste(where.plink, "plink.ld", sep="/"), header=F)
  LD_pop2=as.matrix(LD_pop2)
  
  
  
  #Add loci names to table and calculate row means
  LD_pop2_a=data.table::fread(paste(where.plink, "plink.ld", sep="/"), header=F)
  #Add loci names to table
  
  loci_names_all=genepopedit::genepop_detective("genepop_Pop2_ordered_fst.txt", "Loci")
  LD_pop2_a$LD_means=rowMeans(LD_pop2_a, na.rm = T)
  
  LD_pop2_a$loci=loci_names_all
  
  #Create new dataframe with LD means and Loci names
  Mean_LD=cbind(LD_pop2_a$loci, LD_pop2_a$LD_means)
  colnames(Mean_LD)<-c("loci", "mean_LD")
  
  #merge mean LD values with Fst values based on loci names to plot LD vs Fst relationship
  merged=merge(fst, Mean_LD, by=1) 
  merged2=merged[,c(1,7,9)] #remove extra statistics generated by diveRsity that aren't needed now
  
  #ensure LD and Fst are read as numeric 
  merged2$mean_LD=as.numeric(as.character(merged2$mean_LD))
  merged2$Fst=as.numeric(as.character(merged2$Fst))
  
  #save PDF plot of LD-Fst relationship for the species
  pdf(file="LD_FST_relationship_pop2.pdf", width = 8, height=7, bg = "white")
  plot(merged2$Fst, merged2$mean_LD, xlab="Locus specific FST", ylab="Mean LD (within pop 2)", las=1, pch=19)
  dev.off()
  
  #Create matrix with upper tri as Pop1 and lower tri as Pop2
  new<-LD_pop2
  diag(new) <- 1
  new[upper.tri(new)] <- LD_pop1[upper.tri(LD_pop1)]
  
  #for top 500 loci
  end=ncol(new)
  start=end-500
  
  
  #colour palette
  mypalette4<-colorRampPalette(c("gray95","dodgerblue", "blue","midnightblue", "black"))
  
  #Console message
  writeLines("Generating LD heatmap\n\n")
  
  
  pdf(file="LDmatrix_TwoPops_FstOrdered.pdf", width = 13, height=12, bg = "white")
  gplots::heatmap.2(new[start:end, start:end],
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
                    xlab=pop_groups[2],
                    ylab=pop_groups[1],
                    keysize = 1,
                    density.info="none"
                    #lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
                    #lhei = c(20,5), #row height for plot elements
                    #lwid = c(8,30)  #column width for plot elements)
  )
  dev.off()
  
  
  
  #remove extra files created in working directory
  file.remove("fst-[diffCalc]/std_stats.txt")
  file.remove("genepop_all_ordered_fst.txt")
  file.remove("genepop_Pop1_ordered_fst.txt")
  file.remove("genepop_Pop2_ordered_fst.txt")
  file.remove("fst-[diffCalc]/")
  file.remove(paste0(where.PLINK, "/", "PGDtest.ped"))
  file.remove(paste0(where.PLINK, "/", "PGDtest.map"))
  
}

