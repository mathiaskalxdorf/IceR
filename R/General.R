.onAttach <- function(libname, pkgname) {

  library(utils)
  if (!("cleaver" %in% rownames(installed.packages()))) {
    packageStartupMessage(
      paste0(
        "Please install `cleaver` by",
        " `BiocManager::install('cleaver')`"
      )
    )
  }
  if (!("PECA" %in% rownames(installed.packages()))) {
    packageStartupMessage(
      paste0(
        "Please install `PECA` by",
        " `BiocManager::install('PECA')`"
      )
    )
  }

}



'%!in%' <- function(x,y)!('%in%'(x,y))
'%not in%' <- function(x,y)!('%in%'(x,y))


runUI <- function() {
  appDir <- system.file("shiny", "DDAiceR_UI", package = "DDAiceR")
  if (appDir == "") {
    stop("Could not find UI directory. Try re-installing `DDAiceR`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


###Save data.frames or list of data.frames as xlsx file (with separate tabs per list item)
###Data = Data.frame or list of data.frames to be saved
###File = File path where xlsx file should be saved
SaveExcel = function(Data,File)
{
  wb <- createWorkbook("Data")
  if(is.data.frame(Data))
  {
    addWorksheet(wb, "Data")
    writeData(wb, "Data",Data)
  }else
  {
    for(t in 1:length(Data))
    {
      addWorksheet(wb, names(Data)[t])
      writeData(wb, names(Data)[t],Data[[t]])
    }
  }

  saveWorkbook(wb, File, overwrite = TRUE)
}

####################################################################
########Handling MaxQ and Requant data##############################
####################################################################
####Load MaxQ data from summary folder
####Returns prefiltered Protein and Peptide data as separate lists
####if path is directly specified it has to end with \\
load_MaxQ_data <- function(path=NA,remove_contaminants=T,remove_reverse=T,intensity_used=c("LFQ intensity","iBAQ","Intensity"))
{
  options(warn=-1)
  if(is.na(path))
  {
    print("Select a file in the MaxQ output folder")
    path_to_MaxQ <- file.choose()
    temp <- unlist(gregexpr("\\\\",path_to_MaxQ))
    path_to_MaxQ <- substr(path_to_MaxQ,1,temp[length(temp)])
    data_protein <- read.table(paste(path_to_MaxQ,"proteinGroups.txt",sep=""),sep="\t",header=T)
    data_peptide <- read.table(paste(path_to_MaxQ,"peptides.txt",sep=""),sep="\t",header=T)
    print(paste("Selected path to MaxQuant output:",path_to_MaxQ))
  }else
  {
    path_to_MaxQ <- path
    data_protein <- read.table(paste(path_to_MaxQ,"proteinGroups.txt",sep=""),sep="\t",header=T)
    data_peptide <- read.table(paste(path_to_MaxQ,"peptides.txt",sep=""),sep="\t",header=T)
  }


  ###Remove reverse hits and contaminants?
  data_protein <- data_protein[which(data_protein$Only.identified.by.site != "+"),]

  if(remove_contaminants == T)
  {
    data_protein <- data_protein[which(data_protein$Potential.contaminant != "+"),]
    data_peptide <- data_peptide[which(data_peptide$Potential.contaminant != "+"),]
  }
  if(remove_reverse == T)
  {
    data_protein <- data_protein[which(data_protein$Reverse != "+"),]
    data_peptide <- data_peptide[which(data_peptide$Reverse != "+"),]
  }
  ###Extract quant columns
  intensity_used <- gsub(" ","\\\\.",intensity_used[1])
  data_protein_quant <- data_protein[,which(grepl(paste(intensity_used,"\\.",sep=""),colnames(data_protein)))]
  data_peptide_quant <- data_peptide[,which(grepl("Intensity\\.",colnames(data_peptide)))]

  ###Extract number of quantified peptides per protein and sample
  data_protein_peptide_count <- data_protein[,which(grepl("^Peptides\\.",colnames(data_protein)))]

  ###missing values are marked as NA
  data_protein_quant[data_protein_quant == 0] <- NA
  data_peptide_quant[data_peptide_quant == 0] <- NA

  ###log2 transform quant data
  data_protein_quant <- log2(data_protein_quant)
  data_peptide_quant <- log2(data_peptide_quant)

  ###prepare additional information per row on protein level
  data_protein_info <- data.frame(Gene_name=str_split(data_protein$Gene.names,";",simplify = T)[,1],
                                  ID=str_split(data_protein$Majority.protein.IDs,";",simplify = T)[,1],
                                  Organism=substr(data_protein$Fasta.headers,regexpr("OS=",data_protein$Fasta.headers)+3,regexpr("GN=",data_protein$Fasta.headers)-2),
                                  num_peptides=data_protein$Peptides,
                                  num_unique_peptides=data_protein$Unique.peptides,
                                  Gene_names_all=data_protein$Gene.names,
                                  IDs_major=data_protein$Majority.protein.IDs)
  data_protein_info$Gene_name <- as.character(data_protein_info$Gene_name)
  data_protein_info$ID <- as.character(data_protein_info$ID)
  data_protein_info$Organism <- as.character(data_protein_info$Organism)
  data_protein_info$Gene_names_all <- as.character(data_protein_info$Gene_names_all)
  data_protein_info$IDs_major <- as.character(data_protein_info$IDs_major)
  if(any(data_protein_info$Gene_name == ""))
  {
    data_protein_info$Gene_name <- as.character(data_protein_info$Gene_name)
    data_protein_info$Gene_name[which(data_protein_info$Gene_name == "")] <- "Unassigned"
  }
  data_protein_info$Gene_name <- make.unique(data_protein_info$Gene_name)

  ###prepare additional information per row on peptide level
  data_peptide_info <- data.frame(Sequence=data_peptide$Sequence,
                                  Gene_name=str_split(data_peptide$Gene.names,";",simplify = T)[,1],
                                  ID=str_split(data_peptide$Leading.razor.protein,";",simplify = T)[,1],
                                  Organism=data_protein_info$Organism[match(str_split(data_peptide$Leading.razor.protein,";",simplify = T)[,1],data_protein_info$ID)],
                                  Gene_names_all=data_peptide$Gene.names,
                                  IDs_major=data_peptide$Proteins,
                                  Start.position=data_peptide$Start.position,
                                  End.position=data_peptide$End.position,
                                  Score=data_peptide$Score,
                                  PEP=data_peptide$PEP)

  if(any(data_peptide_info$Gene_name == ""))
  {
    data_peptide_info$Gene_name <- as.character(data_peptide_info$Gene_name)
    data_peptide_info$Gene_name[which(data_peptide_info$Gene_name == "")] <- "Unassigned"
  }

  if(any(is.na(data_peptide_info$Organism)))
  {
    data_peptide_info$Organism[is.na(data_peptide_info$Organism)] <- ""
  }

  ###make rownames of peptide and protein quant more readable
  rownames(data_protein_quant) <- data_protein_info$Gene_name
  rownames(data_protein_peptide_count) <- data_protein_info$Gene_name
  rownames(data_peptide_quant) <- data_peptide_info$Sequence

  return(list(Protein_level=list(Quant_data=data_protein_quant,Meta_data=data_protein_info,Num_peptides_per_quant=data_protein_peptide_count),Peptide_level=list(Quant_data=data_peptide_quant,Meta_data=data_peptide_info)))

  options(warn=0)
}

####Load Requant data from output folder
####Returns prefiltered Protein and feature data as separate lists
####if paths are directly specified it has to end with \\
load_Requant_data <- function(path_to_requant_folder=NA,file_name_extension=NA,path_MaxQ=NA,quant_value=c("LFQ","Total","Top3"),imputed=T)
{
  options(warn=-1)
  library(openxlsx)
  quant_value <- quant_value[1]
  if(is.na(path_to_requant_folder) | is.na(file_name_extension))
  {
    print("Select Features_quantification.xlsx file")
    path_to_requant <- file.choose()
    temp <- unlist(gregexpr("\\\\",path_to_requant))
    path_to_requant_folder <- substr(path_to_requant,1,temp[length(temp)])
    file_name_extension <- substr(path_to_requant,temp[length(temp)]+1,200)
    file_name_extension <- gsub("Features_quantification|\\.xlsx","",file_name_extension)
    print(paste("Selected path to MaxRequanteR output:",path_to_requant_folder))
    print(paste("Selected MaxRequanteR output name:",file_name_extension))
  }else
  {
    file_name_extension <- paste("_",file_name_extension,sep="")
    print(paste("Selected path to MaxRequanteR output:",path_to_requant_folder))
    print(paste("Selected MaxRequanteR output name:",file_name_extension))
  }

  if(is.na(path_MaxQ))
  {
    print("Select a file in the corresponding MaxQuant output folder")
    path_to_MaxQ <- file.choose()
    temp <- unlist(gregexpr("\\\\",path_to_MaxQ))
    path_to_MaxQ <- substr(path_to_MaxQ,1,temp[length(temp)])
    MaxQ_data <- read.table(paste(path_to_MaxQ,"proteinGroups.txt",sep=""),sep="\t",header=T)
    MaxQ_data$Organism <- substr(MaxQ_data$Fasta.headers,regexpr("OS=",MaxQ_data$Fasta.headers)+3,regexpr("GN=",MaxQ_data$Fasta.headers)-2)
    MaxQ_data$ID <- str_split(MaxQ_data$Majority.protein.IDs,";",simplify = T)[,1]
    print(paste("Selected path to MaxQuant output:",path_to_MaxQ))
  }else
  {
    path_to_MaxQ <- path_MaxQ
    MaxQ_data <- read.table(paste(path_to_MaxQ,"proteinGroups.txt",sep=""),sep="\t",header=T)
    MaxQ_data$Organism <- substr(MaxQ_data$Fasta.headers,regexpr("OS=",MaxQ_data$Fasta.headers)+3,regexpr("GN=",MaxQ_data$Fasta.headers)-2)
    MaxQ_data$ID <- str_split(MaxQ_data$Majority.protein.IDs,";",simplify = T)[,1]
  }

  if(quant_value == "Total" & imputed == F)protein_quant_tab <- "Total"
  if(quant_value == "Total" & imputed == T)protein_quant_tab <- "Total_imputed"
  if(quant_value == "Top3" & imputed == F)protein_quant_tab <- "Top3"
  if(quant_value == "Top3" & imputed == T)protein_quant_tab <- "Top3_imputed"
  if(quant_value == "LFQ" & imputed == F)protein_quant_tab <- "LFQ"
  if(quant_value == "LFQ" & imputed == T)protein_quant_tab <- "LFQ_imputed"



  if(imputed == F)feature_quant_tab <- "Quantification"
  if(imputed == T)feature_quant_tab <- "Quantification_imputed"


  sheets_protein <- getSheetNames(paste(path_to_requant_folder,"Proteins_quantification",file_name_extension,".xlsx",sep=""))
  sheets_features <- getSheetNames(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""))

  ###protein level data
  if(imputed == T & protein_quant_tab %not in% sheets_protein)protein_quant_tab <- gsub("_imputed","",protein_quant_tab)

  data_protein <- read.xlsx(paste(path_to_requant_folder,"Proteins_quantification",file_name_extension,".xlsx",sep=""),protein_quant_tab)
  colnames(data_protein) <- gsub("^X","",colnames(data_protein))

  ###feature level data
  if(imputed == T & feature_quant_tab %not in% sheets_features)feature_quant_tab <- gsub("_imputed","",feature_quant_tab)
  features <- read.xlsx(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""),"features",skipEmptyRows = F)
  features_sample_matrix <- read.xlsx(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""),feature_quant_tab,skipEmptyRows = F)
  features_sample_matrix_counts <- read.xlsx(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""),"Ion_count",skipEmptyRows = F)
  features_sample_matrix_pvals <- read.xlsx(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""),"Quantification_pvals",skipEmptyRows = F)
  features_sample_matrix_S2B <- read.xlsx(paste(path_to_requant_folder,"Features_quantification",file_name_extension,".xlsx",sep=""),"Signal_to_background",skipEmptyRows = F)

  data_protein_quant <- data_protein[,which(!grepl("median_alignment_score|median_quant_pvals|median_S2B",colnames(data_protein)))[-c(1:3)]]
  data_protein_pvals <- data_protein[,which(grepl("median_quant_pvals",colnames(data_protein)))]
  data_protein_S2B <- data_protein[,which(grepl("median_S2B",colnames(data_protein)))]

  ##prepare additional information per row on protein level
  data_protein_info <- data.frame(Gene_name=str_split(data_protein$Gene_Name,";",simplify = T)[,1],
                                  ID=str_split(data_protein$UniProt_Identifier,";",simplify = T)[,1],
                                  Organism=MaxQ_data$Organism[match(str_split(data_protein$UniProt_Identifier,";",simplify = T)[,1],MaxQ_data$ID)],
                                  num_quant_features=data_protein$num_quant_features,
                                  Gene_names_all=data_protein$Gene_Name)
  data_protein_info$Gene_name <- as.character(data_protein_info$Gene_name)
  data_protein_info$ID <- as.character(data_protein_info$ID)
  data_protein_info$Organism <- as.character(data_protein_info$Organism)
  data_protein_info$Gene_names_all <- as.character(data_protein_info$Gene_names_all)

  if(any(is.na(data_protein_info$Gene_name)))
  {
    data_protein_info$Gene_name <- as.character(data_protein_info$Gene_name)
    data_protein_info$Gene_name[which(is.na(data_protein_info$Gene_name))] <- "Unassigned"
  }

  ##keep IDs which correspond to Major protein id from MaxQ
  keep <- which(data_protein_info$ID %in% features$Protein)###quantified by at least 1 unique feature
  keep <- unique(append(keep,which(data_protein_info$ID %in% str_split(MaxQ_data$Protein.IDs[which(grepl(";",MaxQ_data$Protein.IDs))],";",simplify = T)[,1]))) ###also keep leading protein id if only quantified by overlapping peptides
  data_protein <- data_protein[keep,]
  data_protein_quant <- data_protein_quant[keep,]
  data_protein_pvals <- data_protein_pvals[keep,]
  data_protein_S2B <- data_protein_S2B[keep,]
  data_protein_info <- data_protein_info[keep,]

  ##make gene names unique if some are existing duplicated
  data_protein_info$Gene_name <- make.unique(data_protein_info$Gene_name)

  ##filter feature data for either unknown feature or specific peptide sequence (no overlapping features)
  features_sample_matrix <- features_sample_matrix[which(!grepl(";|\\|",features$Sequence)),]
  features_sample_matrix_counts <- features_sample_matrix_counts[which(!grepl(";|\\|",features$Sequence)),]
  features_sample_matrix_pvals <- features_sample_matrix_pvals[which(!grepl(";|\\|",features$Sequence)),]
  features_sample_matrix_S2B <- features_sample_matrix_S2B[which(!grepl(";|\\|",features$Sequence)),]

  features <- features[which(!grepl(";|\\|",features$Sequence)),]

  ###prepare additional information per row on feature level
  ##which protein ID should be kept in case if a peptide belongs to two or more protein IDs
  multi_IDs <- str_split(features$Protein,";|\\|",simplify = T)
  IDs_used <- vector(mode="character",length(nrow(features)))
  for(i in 1:nrow(multi_IDs))
  {
    IDs_used[i] <- multi_IDs[i,which(multi_IDs[i,] %in% MaxQ_data$ID)[1]]
  }
  IDs_used[is.na(IDs_used)] <- ""

  data_peptide_info <- data.frame(Sequence=features$Sequence,
                                  Gene_name=data_protein_info$Gene_name[match(str_split(features$Protein,";|\\|",simplify = T)[,1],data_protein_info$ID)],
                                  ID=IDs_used,
                                  Feature_name=features$Feature_name,
                                  Organism=data_protein_info$Organism[match(str_split(features$Protein,";|\\|",simplify = T)[,1],data_protein_info$ID)],
                                  Charge=features$Charge,
                                  IDs_major=features$Protein,
                                  Score=features$mean_Scores)
  data_peptide_info$Sequence <- as.character(data_peptide_info$Sequence)
  data_peptide_info$Gene_name <- as.character(data_peptide_info$Gene_name)
  data_peptide_info$ID <- as.character(data_peptide_info$ID)
  data_peptide_info$Feature_name <- as.character(data_peptide_info$Feature_name)
  data_peptide_info$Organism <- as.character(data_peptide_info$Organism)
  data_peptide_info$IDs_major <- as.character(data_peptide_info$IDs_major)

  if(any(is.na(data_peptide_info$Gene_name)))
  {
    data_peptide_info$Gene_name <- as.character(data_peptide_info$Gene_name)
    data_peptide_info$Gene_name[which(is.na(data_peptide_info$Gene_name))] <- "Unassigned"
  }

  if(any(is.na(data_peptide_info$Organism)))
  {
    data_peptide_info$Organism[is.na(data_peptide_info$Organism)] <- ""
  }
  ##keep IDs which correspond to Major protein id from MaxQ
  features <- features[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]

  features_sample_matrix <- features_sample_matrix[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]
  features_sample_matrix_counts <- features_sample_matrix_counts[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]
  features_sample_matrix_pvals <- features_sample_matrix_pvals[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]
  features_sample_matrix_S2B <- features_sample_matrix_S2B[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]

  data_peptide_info <- data_peptide_info[which(data_peptide_info$ID %in% MaxQ_data$ID | is.na(data_peptide_info$ID) | data_peptide_info$ID == ""),]

  ###make rownames of peptide and protein quant more readable
  rownames(data_protein_quant) <- data_protein_info$Gene_name
  rownames(data_protein_pvals) <- data_protein_info$Gene_name
  rownames(data_protein_S2B) <- data_protein_info$Gene_name

  rownames(features_sample_matrix) <- data_peptide_info$Feature_name
  rownames(features_sample_matrix_counts) <- data_peptide_info$Feature_name
  rownames(features_sample_matrix_pvals) <- data_peptide_info$Feature_name
  rownames(features_sample_matrix_S2B) <- data_peptide_info$Feature_name

  return(list(Protein_level=list(Quant_data=data_protein_quant,
                                 Meta_data=data_protein_info,
                                 Quant_pVal=data_protein_pvals,
                                 Quant_S2B=data_protein_S2B),
              Peptide_level=list(Quant_data=features_sample_matrix,
                                 Ion_counts=features_sample_matrix_counts,
                                 Quant_pVal=features_sample_matrix_pvals,
                                 S2B=features_sample_matrix_S2B,
                                 Meta_data=data_peptide_info,
                                 Meta_data_full=features)))
  options(warn=0)
}

###Adds sample annotation information to Protein and Peptide quant data
###Annotation information in rows with nrow = ncol of quant data in data_list containing protein and peptide level data
add_annotations <- function(data_list,Annotations)
{
  data_list$Annotations <- Annotations
  return(data_list)
}

###Change column names (sample names)
set_sample_names <- function(data_list,sample_names)
{
  if("Quant_data" %in% names(data_list$Protein_level))colnames(data_list$Protein_level$Quant_data) <- sample_names
  if("Quant_data_norm" %in% names(data_list$Protein_level))colnames(data_list$Protein_level$Quant_data_norm) <- sample_names

  if("Quant_data" %in% names(data_list$Peptide_level))colnames(data_list$Peptide_level$Quant_data) <- sample_names
  if("Quant_data_norm" %in% names(data_list$Peptide_level))colnames(data_list$Peptide_level$Quant_data_norm) <- sample_names
  return(data_list)
}

###Determine general numbers of the data set like number of proteins, number of peptides and number of missing values
determine_general_numbers <- function(data_list)
{
  if("Protein_level" %in% names(data_list))
  {
    ###number of proteins
    temp <- plyr::count(data_list$Protein_level$Meta_data$Organism)
    colnames(temp) <- c("Organism","Count")
    data_list$Protein_level$num_prots <- temp

    ###number of proteins per samples
    temp <- as.data.frame(t(as.matrix((colSums(!is.na(data_list$Protein_level$Quant_data))))))
    data_list$Protein_level$num_prots_per_sample <- temp

    ###missing values on protein level - absolute and relative
    temp <- c(length(which(is.na(data_list$Protein_level$Quant_data))),
              length(which(is.na(data_list$Protein_level$Quant_data)))/(nrow(data_list$Protein_level$Quant_data)*ncol(data_list$Protein_level$Quant_data))*100)
    names(temp) <- c("absolute","relative")
    data_list$Protein_level$missing_values <- temp

  }

  if("Peptide_level" %in% names(data_list))
  {
    #check if we are currently looking at requant data
    if("Meta_data_full" %in% names(data_list$Peptide_level))#requant data
    {
      ###missing values on peptide level - absolute and relative
      temp <- c(length(which(is.na(data_list$Peptide_level$Quant_data[which(!grepl("_i|_d|_pmp",data_list$Peptide_level$Meta_data$Feature_name)),]))),
                length(which(is.na(data_list$Peptide_level$Quant_data[which(!grepl("_i|_d|_pmp",data_list$Peptide_level$Meta_data$Feature_name)),])))/(nrow(data_list$Peptide_level$Quant_data[which(!grepl("_i|_d|_pmp",data_list$Peptide_level$Meta_data$Feature_name)),])*ncol(data_list$Peptide_level$Quant_data))*100)
    }else #any other data
    {
      ###missing values on peptide level - absolute and relative
      temp <- c(length(which(is.na(data_list$Peptide_level$Quant_data))),
                length(which(is.na(data_list$Peptide_level$Quant_data)))/(nrow(data_list$Peptide_level$Quant_data)*ncol(data_list$Peptide_level$Quant_data))*100)
    }

    names(temp) <- c("absolute","relative")
    data_list$Peptide_level$missing_values <- temp

    ###number of peptides
    temp <- length(unique(data_list$Peptide_level$Meta_data$Sequence))
    names(temp) <- "num_peptides"
    data_list$Peptide_level$num_peptides <- temp

    ###number of peptides per samples
    temp <- as.data.frame(t(as.matrix((colSums(!is.na(data_list$Peptide_level$Quant_data))))))
    data_list$Peptide_level$num_peptides_per_sample <- temp
  }

  return(data_list)
}

###Compare general numbers in plots
###list_of_data_lists = named ist of data_lists e.g. list(MaxQ=MaxQ_data,Requant=Requant_data)
###colors = colors which should be used per sample
compare_general_numbers <- function(list_of_data_lists,colors=c("darkgrey","chocolate3"),Legendpos = "top",margins=c(12,4,4,9),inset=c(-0.435,0))
{
  ###protein numbers
  dat_prot <- NULL
  names_available=NULL
  for(i in names(list_of_data_lists))
  {
    if("Protein_level" %in% names(list_of_data_lists[[i]]))
    {
      if(is.null(dat_prot))
      {
        dat_prot <- list_of_data_lists[[i]]$Protein_level$num_prots
      }else
      {
        dat_prot <- full_join(dat_prot,list_of_data_lists[[i]]$Protein_level$num_prots,by="Organism")
      }
      names_available <- append(names_available,i)
    }

  }
  colnames(dat_prot)[2:ncol(dat_prot)] <- names_available
  orgs <- dat_prot$Organism
  dat_prot <- dat_prot[,-1]
  dat_prot <- as.data.frame(t(dat_prot))
  rownames(dat_prot) <- names_available
  colnames(dat_prot) <- orgs

  if(length(colors) > length(names_available))colors <- colors[1:length(names_available)]

  ###peptide numbers
  dat_peps <- NULL
  names_available=NULL
  for(i in names(list_of_data_lists))
  {
    if("Peptide_level" %in% names(list_of_data_lists[[i]]))
    {
      if(is.null(dat_peps))
      {
        dat_peps <- list_of_data_lists[[i]]$Peptide_level$num_peptides
      }else
      {
        dat_peps <- append(dat_peps,list_of_data_lists[[i]]$Peptide_level$num_peptides)
      }
      names_available <- append(names_available,i)
    }else
    {
      dat_peps <- append(dat_peps,NA)
      names_available <- append(names_available,i)
    }
  }
  names(dat_peps) <- names_available

  ###missing values
  missing_values_prot <- NULL
  for(i in names(list_of_data_lists)) ###protein level
  {
    if("Protein_level" %in% names(list_of_data_lists[[i]]))
    {
      missing_values_prot <- append(missing_values_prot,list_of_data_lists[[i]]$Protein_level$missing_values[2])
    }else
    {
      missing_values_prot <- append(missing_values_prot,NA)
    }
  }
  missing_values_peps <- NULL
  for(i in names(list_of_data_lists)) ###peptide level
  {
    if("Peptide_level" %in% names(list_of_data_lists[[i]]))
    {
      missing_values_peps <- append(missing_values_peps,list_of_data_lists[[i]]$Peptide_level$missing_values[2])
    }else
    {
      missing_values_peps <- append(missing_values_peps,NA)
    }
  }
  missing_values <- cbind(missing_values_prot,missing_values_peps)
  rownames(missing_values) <- names(list_of_data_lists)
  colnames(missing_values) <- c("Protein","Peptide")


  ###plot protein count
  p <- BarplotsSBS(dat_prot,main="Number of identified proteins",ylab="Count",col=colors,AvgLine = F,shownumbers = T,Legendtitle = "Data",Legends = rownames(dat_prot),Legendpos = Legendpos,margins=margins,inset=inset)

  ###plot peptide count
  p <- Barplots(dat_peps,main="Number of identified peptides",ylab="Count",col=colors,AvgLine = F,shownumbers = T,Legendtitle = "Data",Legends = names(dat_peps),Legendpos = Legendpos,Legendscol = colors,Name = names(dat_peps),xlab="",margins=margins,inset=inset)

  ###plot missing value rates
  p <- BarplotsSBS(missing_values,AvgLine = F,Name = colnames(missing_values),main="Missing value rate",ylab="Fraction missing values [%]",col=colors,Legendtitle = "Data",Legends = names(list_of_data_lists),Legendpos = Legendpos,shownumbers = T,margins=margins,inset=inset)

  ###plot number of quantified proteins per sample and data_list

  for(i in names(list_of_data_lists))
  {
    if("Protein_level" %in% names(list_of_data_lists[[i]]))
    {
      Barplots(list_of_data_lists[[i]]$Protein_level$num_prots_per_sample[,order(as.numeric(as.matrix(list_of_data_lists[[i]]$Protein_level$num_prots_per_sample)))],
               shownumbers = F,
               col=colors[which(names(list_of_data_lists) == i)],
               xlab = "",
               ylab="Count",
               main=paste(i,"- Number of quantified proteins"))
    }
    if("Peptide_level" %in% names(list_of_data_lists[[i]]))
    {
      Barplots(list_of_data_lists[[i]]$Peptide_level$num_peptides_per_sample[,order(as.numeric(as.matrix(list_of_data_lists[[i]]$Peptide_level$num_peptides_per_sample)))],
               shownumbers = F,
               col=colors[which(names(list_of_data_lists) == i)],
               xlab = "",
               ylab="Count",
               main=paste(i,"- Number of quantified peptides"))
    }
  }
}

###Plot accuracy of quantification on protein and peptide level
###list_of_data_lists = named ist of data_lists e.g. list(MaxQ=MaxQ_data,Requant=Requant_data)
plot_accuracy <- function(list_of_data_lists,colors=c("darkgrey","chocolate3"),Legendpos = "topleft",margins=c(2,4,4,10),inset=c(-0.435,0),Annotation_column=1,plot_for=c("both","protein","peptide"),allow_missing_values=T,representation_of=c("total","per group"),show_numbers=T,numbers_size=1,round_to_k=T)
{
  plot_for <- plot_for[1]
  representation_of <- representation_of[1]

  SDs_protein_list <- list()
  SDs_peptide_list <- list()
  for(n in names(list_of_data_lists))
  {
    if(plot_for %in% c("both","protein"))
    {
      ##protein
      temp_dat <- list_of_data_lists[[n]]$Protein_level$Quant_data_norm
      temp_anno <- list_of_data_lists[[n]]$Annotations[,Annotation_column]
      SDs <- as.data.frame(matrix(ncol=4,nrow=nrow(temp_dat)*length(unique(temp_anno))))
      colnames(SDs) <- c("SD","Dilution","Method","Intensity")
      SDs$Dilution <- sort(rep(unique(temp_anno),nrow(temp_dat)))
      SDs$Method <- paste("Protein_",n,sep="")
      for(d in unique(temp_anno))
      {
        temp_dat_2 <- 2^as.matrix(temp_dat[,which(temp_anno==d)])
        temp_sds <- rowSds(temp_dat_2,na.rm = allow_missing_values)
        temp_mean <- rowMeans(temp_dat_2,na.rm = allow_missing_values)
        temp_sds <- (temp_sds/temp_mean)*100
        missing_val <- which(rowSums(is.na(temp_dat_2))>0)
        temp_sds[missing_val] <- NA
        SDs[which(SDs$Dilution == d),1] <- temp_sds
        SDs[which(SDs$Dilution == d),4] <- rowMedians(temp_dat_2,na.rm=T)
      }
      SDs_protein_list[[n]] <- SDs
    }
    ##peptide
    if(plot_for %in% c("both","peptide"))
    {
      temp_dat <- list_of_data_lists[[n]]$Peptide_level$Quant_data_norm
      temp_anno <- list_of_data_lists[[n]]$Annotations[,Annotation_column]
      SDs <- as.data.frame(matrix(ncol=4,nrow=nrow(temp_dat)*length(unique(temp_anno))))
      colnames(SDs) <- c("SD","Dilution","Method","Intensity")
      SDs$Dilution <- sort(rep(unique(temp_anno),nrow(temp_dat)))
      SDs$Method <- paste("Peptide_",n,sep="")
      for(d in unique(temp_anno))
      {
        temp_dat_2 <- 2^as.matrix(temp_dat[,which(temp_anno==d)])
        temp_sds <- rowSds(temp_dat_2,na.rm = allow_missing_values)
        temp_mean <- rowMeans(temp_dat_2,na.rm = allow_missing_values)
        temp_sds <- (temp_sds/temp_mean)*100
        missing_val <- which(rowSums(is.na(temp_dat_2))>0)
        temp_sds[missing_val] <- NA
        SDs[which(SDs$Dilution == d),1] <- temp_sds
        SDs[which(SDs$Dilution == d),4] <- rowMedians(temp_dat_2,na.rm=T)
      }
      SDs_peptide_list[[n]] <- SDs
    }
  }

  ###Combine and plot
  SDs_combined_protein <- NULL
  SDs_combined_peptide <- NULL
  ##protein level
  for(n in unique(c(names(SDs_protein_list),names(SDs_peptide_list))))
  {
    if(plot_for %in% c("both","protein"))
    {
      if(is.null(SDs_combined_protein))
      {
        SDs_combined_protein <- SDs_protein_list[[n]]
      }else
      {
        SDs_combined_protein <- rbind(SDs_combined_protein,SDs_protein_list[[n]])
      }
    }

    if(plot_for %in% c("both","peptide"))
    {
      if(is.null(SDs_combined_peptide))
      {
        SDs_combined_peptide <- SDs_peptide_list[[n]]
      }else
      {
        SDs_combined_peptide <- rbind(SDs_combined_peptide,SDs_peptide_list[[n]])
      }
    }

  }

  if(representation_of == "total")
  {
    unit_count_proteins <- ""
    unit_count_peptides <- ""
    ##determine number of quantifications on protein and peptide level
    if(plot_for %in% c("both","protein"))
    {
      count_protein <- plyr::count(SDs_combined_protein$Method[which(!is.na(SDs_combined_protein$SD))])
      if(any(count_protein$freq>2000) & round_to_k == T)
      {
        count_protein$freq <- count_protein$freq/1000
        unit_count_proteins <- "K"
      }
    }
    if(plot_for %in% c("both","peptide"))
    {
      count_peptide <- plyr::count(SDs_combined_peptide$Method[which(!is.na(SDs_combined_peptide$SD))])
      if(any(count_peptide$freq>2000) & round_to_k == T)
      {
        count_peptide$freq <- count_peptide$freq/1000
        unit_count_peptides <- "K"
      }
    }

    par_save <- par()
    par(mar=margins)
    ###plot
    if(plot_for %in% c("both","protein"))
    {
      p <- boxplot(SDs_combined_protein$SD~SDs_combined_protein$Method,outline=F,las=2,col=colors,xaxt = "n",xlab="",ylab="CV [%]",main="Protein - Precision of quantification")
      par(xpd=T)
      legend(Legendpos,inset=inset, legend = names(list_of_data_lists), fill = colors,cex=0.8,text.font=1,horiz=F,border = NA,bg="transparent",box.col = NA,ncol=1)
      par(xpd=F)
      range_y <- par("usr")[4]-par("usr")[3]
      if(show_numbers == T)
      {
        for(n in 1:length(names(SDs_protein_list)))
        {
          y <- p$stats[3,n]+(range_y*0.05)
          text(n,y,paste(round(count_protein$freq[which(count_protein$x==paste("Protein_",names(SDs_protein_list)[n],sep=""))],digits = 0),unit_count_proteins,sep=""),cex=numbers_size)
        }
      }

    }

    if(plot_for %in% c("both","peptide"))
    {
      p <- boxplot(SDs_combined_peptide$SD~SDs_combined_peptide$Method,outline=F,las=2,col=colors,xaxt = "n",xlab="",ylab="CV [%]",main="Peptide - Precision of quantification")
      par(xpd=T)
      legend(Legendpos,inset=inset, legend = names(list_of_data_lists), fill = colors,cex=0.8,text.font=1,horiz=F,border = NA,bg="transparent",box.col = NA,ncol=1)
      par(xpd=F)
      range_y <- par("usr")[4]-par("usr")[3]
      if(show_numbers == T)
      {
        for(n in 1:length(names(SDs_peptide_list)))
        {
          y <- p$stats[3,n]+(range_y*0.05)
          text(n,y,paste(round(count_peptide$freq[which(count_peptide$x==paste("Peptide_",names(SDs_peptide_list)[n],sep=""))],digits = 0),unit_count_peptides,sep=""),cex=numbers_size)
        }
      }
    }
    par <- par_save
  }

  if(representation_of == "per group")
  {
    unit_count_proteins <- ""
    unit_count_peptides <- ""

    if(plot_for %in% c("both","protein"))col_name <- colnames(SDs_combined_protein)[2]
    if(plot_for %in% c("both","peptide"))col_name <- colnames(SDs_combined_peptide)[2]
    ##determine number of quantifications on protein and peptide level
    if(plot_for %in% c("both","protein"))
    {
      count_protein <- plyr::count(SDs_combined_protein[which(!is.na(SDs_combined_protein$SD)),],vars = c("Method",col_name))
      if(any(count_protein$freq>2000) & round_to_k == T)
      {
        count_protein$freq <- count_protein$freq/1000
        unit_count_proteins <- "K"
      }
      count_protein <- count_protein[order(count_protein[,2]),]
    }
    if(plot_for %in% c("both","peptide"))
    {
      count_peptide <- plyr::count(SDs_combined_peptide[which(!is.na(SDs_combined_peptide$SD)),],vars = c("Method",col_name))
      if(any(count_peptide$freq>2000) & round_to_k == T)
      {
        count_peptide$freq <- count_peptide$freq/1000
        unit_count_peptides <- "K"
      }
      count_peptide <- count_peptide[order(count_peptide[,2]),]
    }
    par_save <- par()
    par(mar=margins)
    if(plot_for %in% c("both","protein"))names <- SDs_combined_protein[!duplicated(SDs_combined_protein[,2:3]),]
    if(plot_for %in% c("both","peptide"))names <- SDs_combined_peptide[!duplicated(SDs_combined_peptide[,2:3]),]
    names <- names[order(names[,2]),2]
    ###plot
    if(plot_for %in% c("both","protein"))
    {
      p <- boxplot(SDs_combined_protein$SD~SDs_combined_protein$Method+SDs_combined_protein[,2],outline=F,las=2,col=colors,xlab="",names=names,ylab="CV [%]",main="Protein - Precision of quantification")
      par(xpd=T)
      legend(Legendpos,inset=inset, legend = names(list_of_data_lists), fill = colors,cex=0.8,text.font=1,horiz=F,border = NA,bg="transparent",box.col = NA,ncol=1)
      par(xpd=F)
      range_y <- par("usr")[4]-par("usr")[3]
      for(n in 1:nrow(count_protein))
      {
        y <- p$stats[3,n]+(range_y*0.02)
        text(n,y,paste(round(count_protein$freq[n],digits = 0),unit_count_proteins,sep=""),cex=numbers_size/2)
      }
    }

    if(plot_for %in% c("both","peptide"))
    {
      p <- boxplot(SDs_combined_peptide$SD~SDs_combined_peptide$Method+SDs_combined_peptide[,2],outline=F,las=2,col=colors,names=names,xlab="",ylab="CV [%]",main="Peptide - Precision of quantification")
      par(xpd=T)
      legend(Legendpos,inset=inset, legend = names(list_of_data_lists), fill = colors,cex=0.8,text.font=1,horiz=F,border = NA,bg="transparent",box.col = NA,ncol=1)
      par(xpd=F)
      range_y <- par("usr")[4]-par("usr")[3]
      for(n in 1:nrow(count_peptide))
      {
        y <- p$stats[3,n]+(range_y*0.02)
        text(n,y,paste(round(count_peptide$freq[n],digits = 0),unit_count_peptides,sep=""),cex=numbers_size/2)
      }
    }
    par <- par_save

  }



  #n <- length(names(list_of_data_lists))
  #axis(1, at = 0:1*n + (n-(0.5*(n-1))), labels = c("Peptide-level","Protein-level"), tick = TRUE)

  #legend(Legendpos,inset=c(-0.12,0), legend = names(list_of_data_lists), fill = colors,cex=0.8,ncol=1,text.font=1,horiz=T,border = NA,bg="transparent",box.col = NA)
}

###Visualize data set on protein or peptide level in a heatmap
###data_list = output from load_MaxQ_data or load_Requant_data
###annotation = vector indicating how samples should be grouped or data frame with different annotations per cols and sampels in rows
###colors_annotation = named vector of colors with names corresponding to annotation group names
###subset_indices = integer vector indicating which rows of protein quant or peptide quant table should be plotted
###used_quant should be norm or raw indicating of normalized or raw intensities should be plotted
###quant_level should be protein or peptide. indicates if quant data on protein or peptide level should be plotted
###breaks should be a numerical vector indicating min and max value between which the colors of the plotting are ranging
###colors_plot should be a vector of 2 or 3 colors which should be used to color data in the range of breaks
###row_order integer vector indicating in which order rows should be plotted
visualize_dataset_heatmap <- function(data_list,annotation=NULL,colors_annotation=NULL,subset_rows=NULL,subset_columns=NULL,used_quant=c("norm","raw"),quant_level=c("protein","peptide"),breaks=NULL,colors_plot=c("deepskyblue","firebrick"),main="",annotation_name="Sample",row_order = NULL,show_rownames=F,rownames=NULL,remove_NA_rows=T,hclust_col=F,hclust_row=F)
{
  library(grDevices)
  used_quant <- used_quant[1]
  quant_level <- quant_level[1]
  ##Only plot a subset of features?
  if(is.null(subset_rows) & quant_level == "protein")subset_rows <- 1:nrow(data_list$Protein_level$Quant_data)
  if(is.null(subset_rows) & quant_level == "peptide")subset_rows <- 1:nrow(data_list$Peptide_level$Quant_data)

  if(is.null(subset_columns) & quant_level == "protein")subset_columns <- 1:ncol(data_list$Protein_level$Quant_data)
  if(is.null(subset_columns) & quant_level == "peptide")subset_columns <- 1:ncol(data_list$Peptide_level$Quant_data)

  if(!is.data.frame(annotation))
  {
    annotation <- annotation[subset_columns]
    annotation <- data.frame(anno=annotation)

  }else
  {
    annotation <- annotation[subset_columns,]
  }

  ##which quantification data should be plotted, raw or normalized. Standard is normalized
  if(used_quant == "norm" & quant_level == "protein")temp_dat <- data_list$Protein_level$Quant_data_norm[subset_rows,subset_columns]
  if(used_quant == "raw" & quant_level == "protein")temp_dat <- data_list$Protein_level$Quant_data[subset_rows,subset_columns]

  if(used_quant == "norm" & quant_level == "peptide")temp_dat <- data_list$Peptide_level$Quant_data_norm[subset_rows,subset_columns]
  if(used_quant == "raw" & quant_level == "peptide")temp_dat <- data_list$Peptide_level$Quant_data[subset_rows,subset_columns]

  if(used_quant == "norm" & quant_level == "protein" & main == "")main="Protein level - Normalized"
  if(used_quant == "norm" & quant_level == "peptide" & main == "")main="Peptide level - Normalized"
  if(used_quant == "raw" & quant_level == "protein" & main == "")main="Protein level - Raw"
  if(used_quant == "raw" & quant_level == "peptide" & main == "")main="Peptide level - Raw"
  ##order samples according to annotation?
  if(is.null(breaks))breaks <- c(min(temp_dat,na.rm=T),max(temp_dat,na.rm=T))

  if(!is.null(annotation))
  {
    if(ncol(annotation) == 1)
    {
      temp_dat <- temp_dat[,do.call(order, as.data.frame(annotation[, 1:ncol(annotation)]))]
      annotation <- as.data.frame(annotation[do.call(order, as.data.frame(annotation[, 1:ncol(annotation)])),])
      colnames(annotation) <- "anno"
    }else
    {
      temp_dat <- temp_dat[,do.call(order, annotation[, 1:ncol(annotation)])]
      annotation <- annotation[do.call(order, annotation[, 1:ncol(annotation)]),]
    }
  }else
  {
    annotation$V1 <- colnames(temp_dat)
  }

  if(is.null(colors_annotation))
  {
    if(is.null(annotation))
    {
      colors_annotation <- hcl.colors(ncol(temp_dat), alpha = 1, rev = FALSE)
      names(colors_annotation) <- colnames(temp_dat)
    }else
    {
      colors_annotation <- list()

      for(c in 1:ncol(annotation))
      {
        colors_annotation[[c]] <- hcl.colors(length(unique(annotation[,c])), alpha = 1, rev = FALSE)
        names(colors_annotation[[c]]) <- unique(annotation[,c])
      }
      names(colors_annotation) <- colnames(annotation)
    }
  }
  if(!is.list(colors_annotation))
  {
    colors_annotation <- list(anno=colors_annotation)
  }

  ###order rows by sum of intensities if no ordering is specified
  if(is.null(row_order))row_order <- order(rowSums(temp_dat,na.rm=T))
  temp_dat <- temp_dat[row_order,]

  ###rownames
  if(show_rownames == T)
  {
    if(!is.null(rownames))
    {
      rownames <- as.character(rownames)[row_order]
    }else
    {
      rownames <- rownames(temp_dat)
    }
  }

  ###remove rows only containing NAs?
  if(remove_NA_rows==T)
  {
    if(!is.null(rownames))rownames <- rownames[which(rowSums(is.na(temp_dat)) < ncol(temp_dat))]
    temp_dat <- temp_dat[which(rowSums(is.na(temp_dat)) < ncol(temp_dat)),]
  }

  ###if names occur several times repeated then only keep the rowname in the middle of all rownames in a block
  if(show_rownames == T)
  {
    last_changed = 0
    for(i in 1:(length(rownames)-1))
    {
      if(rownames[i] != rownames[i+1] | i == (length(rownames)-1))
      {
        if(last_changed == i - 1)
        {
          rownames[i] <- rownames[i]
          last_changed <- i
        }else
        {
          temp_name <- rownames[i]
          if(i != (length(rownames)-1))
          {
            rownames[(last_changed+1):i] <- ""
            delta <- i - (last_changed)
          }else
          {
            rownames[(last_changed+1):(i+1)] <- ""
            delta <- i+1 - (last_changed)
          }

          rownames[last_changed+ceiling(delta/2)] <- temp_name
          last_changed <- i
        }

      }
    }
  }

  ##plot
  Heatmap(data = temp_dat,annotation = annotation,colors_annotation = colors_annotation,hclust_col = hclust_col,hclust_row = hclust_row,show_rownames = show_rownames,colors_plot = colors_plot,color_breaks = breaks,main=main,annotation_name = annotation_name,row_labels=rownames)

}

###Normalize data based on median or density maxima
###data=data frame with features in rows and samples in cols
###method=how the data should be normalized. either "density" --> use density maxium or "median" --> use median for normalization
###norm_to: can be set if values should be normalized to a specific value
###norm_on_subset: numeric vector indicating which rows should be only used for normlization
normalize_data <- function(data,method="density",norm_to=NULL,norm_on_subset=NULL,main="Data")
{
  if(is.null(norm_on_subset))norm_on_subset <- 1:nrow(data)
  densityplots(data[norm_on_subset,],main=paste(main,"- Raw"),xlab="Intensity, log2",col = "black")
  maxima <- NULL
  if(method == "density")
  {
    for(c in 1:ncol(data))
    {
      maxima <- append(maxima,maxDensity(data[norm_on_subset,c]))
    }
  }
  if(method == "median")
  {
    maxima <- colMedians(as.matrix(data[norm_on_subset,]),na.rm=T)
  }
  if(is.null(norm_to))norm_to <- mean(maxima,na.rm=T)
  norm_factors <- norm_to/maxima
  data_norm <- data
  for(c in 1:ncol(data))
  {
    data_norm[,c] <- data_norm[,c]*norm_factors[c]
  }
  densityplots(data_norm,main=paste(main,"- Normalized"),xlab="Intensity, log2",col = "black")
  return(data_norm)
}

###Barplot stacked
#order_groups --> order table according to group names and plot with trellis
#group_names --> character vector of length ncol(Data) indicating which col belongs to which group
Barplotsstacked = function(Data,Name="",ylab="Y-axis",logy=F,main="Titel",col="lightblue",AvgLine=T,Legends=NA,Legendtitle="Legend",Legendpos = "topright",shownumbers=T,shownumbers_total=T,order_groups=F,group_names="",ylim=NULL,margins=c(12,4,4,9),inset=c(-0.3,0))
{
  orig_par <- par()
  curdata <- as.matrix(Data)
  curdata[which(is.na(curdata))] <- 0
  if(any(is.na(Name))){
    Name[is.na(Name)] <- ""
  }

  if(any(Name != ""))
  {
    names <- Name
  }else
  {
    names <- colnames(Data)
  }

  ###order samples by groups
  if(order_groups == T)
  {
    Data <- Data[,order(group_names)]
  }
  if(order_groups == T)
  {
    par(mar=margins,xpd=F)
  }else
  {
    par(mar=margins,xpd=F)
  }
  if(is.null(ylim))
  {
    ymin <- ifelse(min(colSums(curdata,na.rm=T),na.rm=T) < 0,min(colSums(curdata,na.rm=T),na.rm=T) - (0.1*min(colSums(curdata,na.rm=T),na.rm=T)),0)
    ymax <- ifelse(max(colSums(curdata,na.rm=T),na.rm=T) > 0,max(colSums(curdata,na.rm=T),na.rm=T) + (0.1*max(colSums(curdata,na.rm=T),na.rm=T)),0)
  }else
  {
    ymin <- ifelse(ylim[1] < 0,ylim[1] - (0.1*ylim[1]),0)
    ymax <- ifelse(ylim[2] > 0,ylim[2] + (0.1*ylim[2]),0)
  }

  if(logy == T & ymin == 0 & ymax > 1)
  {
    ymin = 1
  }else
  {
    logy = F
  }

  if(logy == F)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",ylim=c(ymin,ymax))
  if(logy == T)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",ylim=c(ymin,ymax),log="y")


  axis(1, at=mids, labels=names, las=3,cex.axis=0.6)
  if(AvgLine == T)
  {
    abline(h=mean(curdata,na.rm=T),lty=2)
    par(xpd=T)
    text(par("usr")[2]+1,mean(curdata,na.rm=T),round(mean(curdata,na.rm=T),digits=0))
  }

  if(!is.na(Legends))
  {
    par(xpd=T)
    defaultpar <- par()
    par(font=2)
    legend(Legendpos,inset=inset, legend = Legends, fill = col,cex=0.8, title=Legendtitle,title,ncol=1,text.font=1)
    par(defaultpar)
    par(xpd=F)
  }

  if(shownumbers == T)
  {
    for(i in 1:ncol(curdata))
    {
      x <- mids[i]
      for(j in 1:nrow(curdata))
      {
        if(abs(curdata[j,i]) >= 0.1*(ymax+ymin)) ###only add label if value >= 10% of yplot area otherwise bar too small
        {
          if(j == 1)
          {
            y <- 0 ###start y value
          }else
          {
            y <- sum(curdata[1:j-1,i],na.rm=T)
          }
          y <- y + (curdata[j,i]/2)
          text(x,y,labels = round(curdata[j,i],digits=1),cex=0.7)
        }
      }
    }
  }

  if(shownumbers_total == T)
  {
    for(i in 1:ncol(curdata))
    {
      x <- mids[i]
      ###plot label complete bar
      y <- sum(curdata[1:nrow(curdata),i],na.rm=T)
      text(x,y,labels = round(y,digits=1),cex=0.7,pos=3)
    }
  }





  ###order samples by groups trellis
  if(order_groups == T & logy == F)
  {
    abline(h=par("usr")[4])
    for(g in sort(unique(group_names)))
    {
      indices <- which(sort(group_names) == g)

      if(min(indices) == 1)
      {
        abline(v=mids[min(indices)]-0.625)
      }
      abline(v=mids[max(indices)]+0.625)
      par(xpd=T)
      text(mean(c(mids[min(indices)],mids[max(indices)])),par("usr")[4]+(par("usr")[4]-par("usr")[3])*0.05,g)
      par(xpd=F)

    }
  }





  par(orig_par)
  return(mids)
}
BarplotsSBS = function(Data,ErrbarData=NA,Name="",ylab="Y-axis",main="Titel",col="lightblue",AvgLine=T,Legends=NA,Legendtitle="Legend",Legendpos = "topright",ylim=NA,logy=F,shownumbers=F,shownumbers_digits=1,separation=T,horiz_line=NULL,margins=c(8,4,4,4),inset=c(-0.1,0))
{
  orig_par <- par()
  curdata <- as.matrix(Data)

  if(any(is.na(Name))){
    Name[is.na(Name)] <- ""
  }

  if(any(Name != ""))
  {
    names <- Name
  }else
  {
    names <- colnames(Data)
  }
  par(mar=margins,xpd=F)
  if(is.na(ylim))
  {
    if(logy == F)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",beside=T)
    if(logy == T)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",beside=T,log="y")
  }else
  {
    if(logy == F)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",beside=T,ylim = ylim)
    if(logy == T)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",beside=T,ylim = ylim,log="y")
  }

  if(separation == T)
  {
    for(i in 1:ncol(mids)-1)
    {
      abline(v=mids[nrow(mids),i]+1,lty=2)
    }

  }

  if(!is.na(ErrbarData))
  {
    for(i in 1:ncol(curdata))
    {
      errbar(mids[,i],curdata[,i], curdata[,i]+ErrbarData[,i], curdata[,i]-ErrbarData[,i], add=T, pch=26, cap=.01)
    }

  }

  axis(1, at=colMeans(mids), labels=names, las=3,cex.axis=0.9)
  if(AvgLine == T)
  {
    abline(h=mean(curdata,na.rm=T),lty=2)
    par(xpd=T)
    text(par("usr")[2]+1,mean(curdata,na.rm=T),round(mean(curdata,na.rm=T),digits=0))
  }

  if(!is.null(horiz_line))
  {
    par(xpd=F)
    abline(h=horiz_line,lty=2)
  }

  if(!is.na(Legends))
  {
    par(xpd=T)
    defaultpar <- par()
    par(font=2)
    legend(Legendpos,inset=inset, legend = Legends, fill = col,cex=0.8, title=Legendtitle,title,ncol=1,text.font=1)
    par(defaultpar)
  }
  par(xpd=T)
  if(shownumbers == T)
  {
    for(i in 1:length(curdata))
    {
      x <- mids[i]
      ###plot label complete bar
      y <- curdata[i]
      text(x,y,labels = round(y,digits = shownumbers_digits),cex=0.7,pos=3,srt=90,offset = 1)
    }
  }
  par(xpd=F)
  par(orig_par)
  return(mids)
}
Barplots = function(Data,ErrbarData=NA,Name="",xlab="X-axis",ylab="Y-axis",main="Titel",col="lightblue",AvgLine=T,digits_average=0,Legends=NA,Legendscol=NA,Legendtitle="Legend",Legendpos = "topright",shownumbers=T,shownumbers_digits=1,ylim=NA,logy=F,margins=c(10.1,4.1,4.1,4.1),inset=c(-0.1,0))
{
  #orig_par <- par()
  curdata <- as.numeric(Data)

  if(any(is.na(Name))){
    Name[is.na(Name)] <- ""
  }

  if(any(Name != ""))
  {
    names <- Name
  }else
  {
    names <- colnames(Data)
  }
  par(mar=margins,xpd=F)

  if(is.na(ylim))
  {
    if(logy == F)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n")
    if(logy == T)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",log="y")
  }else
  {
    if(logy == F)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",ylim = ylim)
    if(logy == T)mids <- barplot(curdata,ylab=ylab,col=col,las=2,main=main,xaxt = "n",ylim = ylim,log="y")
  }



  title(xlab = xlab)
  if(!is.na(ErrbarData))
  {
    errbar(mids,curdata, as.numeric(curdata+ErrbarData), as.numeric(curdata-ErrbarData), add=T, pch=26, cap=.01)
  }
  axis(1, at=mids, labels=names, las=3,cex.axis=0.9)
  if(AvgLine == T)
  {
    abline(h=mean(curdata,na.rm=T),lty=2)
    par(xpd=T)
    text(par("usr")[2]+1,mean(curdata,na.rm=T),round(mean(curdata,na.rm=T),digits=digits_average))
  }
  if(!is.na(Legends))
  {
    par(xpd=T)
    defaultpar <- par()
    par(font=2)
    legend(Legendpos,inset=inset, legend = Legends, fill = Legendscol,cex=0.8, title=Legendtitle,title,ncol=1,text.font=1)
    par(defaultpar)
  }
  par(xpd=T)
  if(shownumbers == T)
  {
    for(i in 1:length(curdata))
    {
      x <- mids[i]
      ###plot label complete bar
      y <- curdata[i]
      text(x,y,labels = round(y,digits = shownumbers_digits),cex=0.7,pos=3)
    }
  }
  par(xpd=F)
  #par(orig_par)
  return(mids)
}

###plot density plots for data.frame with columns = samples and rows = observations
densityplots = function(datafcs,names=NA,col=NA,main="",xlab="",xlim=NA,lwd=2)
{
  #default_par <- par() #save par
  par(mar=c(5.1, 4.1, 4.1, 6.1))
  #####datafcs should contain only cols containg data for which density should be displayed
  densitys <- list()
  datafcs[sapply(datafcs, is.infinite)] <- NA
  densmax <- 0
  for(i in 1:ncol(datafcs))
  {
    if(!is.na(sum(unique(datafcs[,i]),na.rm=T)) & length(unique(datafcs[,i])) > 1 | length(unique(datafcs[,i])) == 1 & !is.na(unique(datafcs[,i])[1]))
    {
      densitys[[i]] <- density(datafcs[,i],na.rm=T)
      if(max(densitys[[i]]$y) > densmax & length(unique(datafcs[,i])) > 1) ### only if its not the reference
      {
        densmax <- max(densitys[[i]]$y)
      }
    }else
    {
      densitys[[i]] <- NA
    }

  }
  if(is.na(col)){col <- randomColor(ncol(datafcs))}
  if(length(col) == 1){col <- rep(col,ncol(datafcs))}
  counter = 0
  for(i in 1:ncol(datafcs))
  {
    if(!is.na(densitys[[i]][3]))
    {
      counter <- counter + 1
      if(counter == 1)
      {
        if(!is.na(xlim))
        {
          plot(densitys[[1]],xlim = xlim,ylim=c(0,densmax),main=main,xlab=xlab,col=col[i],lwd=lwd)
        }else
        {
          plot(densitys[[1]],ylim=c(0,densmax),main=main,xlab=xlab,col=col[i],lwd=lwd)
        }

      }else
      {
        lines(densitys[[i]],col=col[i],lwd=lwd)
      }
      x <- maxDensity(datafcs[,i])
      segments(x0 = x,x1 = x, y0 = 0, y1 = max(densitys[[i]]$y),lty=2,col="grey")
    }
  }

  if(!is.na(names))
  {
    par(xpd=T)
    par(font=2)
    legend("topright",inset=c(-0.25,0), legend = names, fill = col,cex=0.8, title="Samples",title,ncol=1,text.font=1)
    par(xpd=F)
  }
  #par(default_par)
  #par(xpd=F)
}

####Heatmap visualization of data
####data = data.frame, rows = observations, cols = samples
####annotation = vector of length ncol(data), classifying each sample to a class or data frame with different annotations in cols
####annotation_name = a name for the classification e.g. Tumor or a vector of names if more than one annotation should be plotted
####colors_annotation = a named vector with colors if only one annotation should be used or a list of named color vectors
####colors = named vector with names = unique(annotation)
####colors_plot = color vector of length 2 or 3
Heatmap <- function(data,annotation,annotation_name="Anno",main="Heatmap",colors_annotation,colors_plot=NA,color_breaks=NA,hclust_col=T,hclust_row=T,show_rownames=T,show_colnames=F,fontsize_rows=10,interactive=F,row_labels=NULL)
{
  library(gplots)
  library(pheatmap)
  library(heatmaply)

  if(any(is.na(colors_plot)))
  {
    colors_plot <- bluered(100)
  }else
  {
    if(length(colors_plot) == 2)
    {
      colors_plot <- colorpanel(100,colors_plot[1],colors_plot[2])
    }
    if(length(colors_plot) == 2)
    {
      colors_plot <- colorpanel(100,colors_plot[1],colors_plot[2],colors_plot[3])
    }
  }

  if(is.null(row_labels)) row_labels <- rownames(data)

  if(interactive == F)
  {
    if(!is.data.frame(annotation))
    {
      annotation_col = data.frame(class=annotation)
      rownames(annotation_col) <- colnames(data)
      ann_colors = list(class = colors_annotation)
      colnames(annotation_col) <- annotation_name
      names(ann_colors) <- annotation_name
    }else
    {
      annotation_col = annotation
      rownames(annotation_col) <- colnames(data)
      ann_colors = colors_annotation
      colnames(annotation_col) <- colnames(annotation_col)
      if(length(annotation_name) == ncol(annotation_col))
      {
        colnames(annotation_col) <- annotation_name
        names(ann_colors) <- annotation_name
      }

    }

    if(any(!is.na(color_breaks)))
    {
      color_breaks <- seq(from=color_breaks[1],to = color_breaks[2],length.out = length(colors_plot))
      pheatmap(mat = data,color = colors_plot,breaks = color_breaks,annotation_col = annotation_col,annotation_colors=ann_colors,main=main,show_colnames = show_colnames,show_rownames = show_rownames,fontsize_row = fontsize_rows,cluster_rows = hclust_row,cluster_cols = hclust_col,border_color=NA,labels_row=row_labels)
    }else
    {
      max <- max(abs(data),na.rm=T)
      color_breaks = seq(from=-max,to = max,length.out = length(colors_plot))

      pheatmap(mat = data,color = colors_plot,annotation_col = annotation_col,breaks = color_breaks,annotation_colors=ann_colors,main=main,show_colnames = show_colnames,show_rownames = show_rownames,fontsize_row = fontsize_rows,cluster_rows = hclust_row,cluster_cols = hclust_col,border_color=NA,labels_row=row_labels)
    }
  }
  if(interactive == T)
  {
    library(plotly)

    if(any(!is.na(color_breaks)))
    {
      data[data < color_breaks[1]] <- color_breaks[1]
      data[data > color_breaks[2]] <- color_breaks[2]

      heatmaply(as.matrix(data),col = colors_plot,col_side_colors = annotation,col_side_palette = colors_annotation,main=main,showticklabels = c(show_colnames,show_rownames),Rowv = hclust_row,Colv = hclust_col,labRow = row_labels)
    }else
    {
      max <- max(abs(data),na.rm=T)
      color_breaks <- c(-max,max)

      heatmaply(as.matrix(data),col = colors_plot,limits = color_breaks,col_side_colors = annotation,col_side_palette = colors_annotation,main=main,showticklabels = c(show_colnames,show_rownames),Rowv = hclust_row,Colv = hclust_col,labRow = row_labels)
    }
  }

}

####Perform DE analysis using probe-level expression change averaging, e.g. performing DE on protein level by integrating peptide level data
####peptide_data = samples in cols, observations in rows, intensities in log2
####peptide_data_quant_significance = Data.frame with ncol and nrow equal to peptide_data representing significance for every quantification
####ids = e.g. protein names to which respective feature corresponds to, aggregation per protein
####anno = annotation of grouping per sample, 2 groups
####goup1_name = name of group1 in anno
####group2_name = name of group2 in anno
####TopN_ratio = indicates based on how many top-abundant peptides the actual ratio between both groups should be calculated. If NA, ratio is calculated based on all peptides
####pvalue_cutoff = indicates which quantification values should be used for DE. If NA, all quantifications are used
PECA_analysis <- function(peptide_data,peptide_data_quant_significance,ids,anno=NULL,group1_name,group2_name,TopN_ratio=5,pvalue_cutoff=NA)
{
  library(PECA)
  library(matrixStats)
  colnames(peptide_data) <- gsub("^X","",colnames(peptide_data))
  #get names of samples in groups to be compared
  group_1 <- colnames(peptide_data)[which(anno == group1_name)]
  group_2 <- colnames(peptide_data)[which(anno == group2_name)]
  #max_pvalue per feature under investigation
  if(!is.na(pvalue_cutoff))
  {
    max_pval <- rowMaxs(as.matrix(peptide_data_quant_significance[,append(group_1,group_2)]),na.rm=T)
    peptide_data <- peptide_data[which(max_pval < pvalue_cutoff),]
    ids <- ids[which(max_pval < pvalue_cutoff)]
  }

  peptide_data <- data.frame(id=ids,2^peptide_data)
  colnames(peptide_data) <- gsub("^X","",colnames(peptide_data))

  DE_results <- PECA_df(df = peptide_data,id=id,samplenames1 = group_1,samplenames2 = group_2,test = "t",progress = F)
  colnames(DE_results)[c(1,5,6)] <- c("logFC","P.Value","adj.P.Val")

  #calculate ratio based on TopN abundant peptides per protein
  if(!is.na(TopN_ratio))
  {
    peptide_data <- log2(peptide_data[,-1])
    for(i in 1:nrow(DE_results))
    {
      temp_quants <- peptide_data[which(ids == rownames(DE_results)[i]),]
      if(nrow(temp_quants)>TopN_ratio)
      {
        temp_quants <- temp_quants[order(rowSums(temp_quants,na.rm = T),decreasing = T)[1:TopN_ratio],]
        ratio <- median(rowMeans(temp_quants[,which(anno == group1_name)],na.rm=T)-rowMeans(temp_quants[,which(anno == group2_name)],na.rm=T),na.rm=T)
        set(DE_results,as.integer(i),1L,ratio)
      }
    }
  }


  return(DE_results)
}
