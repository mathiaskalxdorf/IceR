#' Function to determine at which value a density maximum is reached
#' @param data Numeric vector
#' @details Uses kernel density estimation function from R-package stats removing missing values
#' @return Numeric indicating at which value a density maximum is reached
#' @export
maxDensity <- function(data)
{
  dens <- stats::density(data,na.rm=T)
  return(dens$x[which(dens$y == max(dens$y))])
}

#' Convert thermo raw files to mzXML with centroided ms1 scans using the ProteoWizard tool msConvert
#' @param path_to_raw Path to folder containing raw files which should be converted
#' @details Requires installation of ProteoWizard (http://proteowizard.sourceforge.net/download.html). Pay attention to installation requirements.
#' @return Resulting mzXML files are stored in a sub-directory within specified raw file folder
#' @export
run_msconvert_raw_mzXML <- function(path_to_raw=NULL)
{
  #suppressWarnings(suppressMessages(library(ff,quietly = T)))
  #suppressWarnings(suppressMessages(library(rChoiceDialogs,quietly = T)))

  if(is.null(path_to_raw))path_to_raw <- rChoiceDialogs::rchoose.dir(caption = "Select folder containing raw files")

  raw_files <- list.files(path_to_raw)
  raw_files <- raw_files[which(grepl("\\.raw",raw_files))]

  ###check which raw files still have to be converted
  mzXMLs_available <- list.files(base::paste(path_to_raw,"/mzXML",sep=""))

  files_to_be_converted <- raw_files[which(base::gsub("\\.raw","",raw_files) %not in% base::gsub("\\.mzXML","",mzXMLs_available))]

  if(length(files_to_be_converted)>0)
  {
    ###get home directory
    home_folder <- Sys.getenv("HOME")
    home_folder <- base::gsub("Documents","AppData/Local/Apps/",home_folder)

    ###find MSConvert folder
    folders <- list.dirs(path = home_folder, full.names = TRUE, recursive = F)
    folders <- folders[which(grepl("ProteoWizard",folders))]
    folders <- folders[length(folders)]

    if(file.exists(base::paste(folders,"\\msconvert.exe",sep="")))
    {
      path_to_msconvert <- base::paste(folders,"\\msconvert.exe",sep="")
    }else
    {
      print("Could not find msconvert.exe. Can be usually found in C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version")
      pb <- tcltk::tkProgressBar("Warning!",min = 0,max = 3,initial = 0,label = "Could not find msConvert.exe.",width = 500)

      counter <- 1
      label <- c("On Windows typically located in C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version",
                 "Please specify location of msConvert.exe")
      while(T)
      {
        Sys.sleep(3)
        tcltk::setTkProgressBar(pb,value = counter,label = label[counter])
        counter <- counter + 1
        if(counter == 4)break
      }
      close(pb)
      path_to_msconvert <- file.choose()
    }

    ###get user folder
    win_user_folder <- path.expand('~')

    ###create temporary folder
    dir.create(base::paste(win_user_folder,"\\temp_msconvert",sep=""),showWarnings = F)
    dir.create(base::paste(win_user_folder,"\\temp_msconvert\\mzXML",sep=""),showWarnings = F)

    ###create temporary config.txt and files.txt
    temp_path <- base::paste(win_user_folder,"\\temp_msconvert",sep="")
    setwd(temp_path)
    ####config
    fileConn<-file("config.txt")
    writeLines(c("mzXML=true",
                 "64=true",
                 "noindex=false",
                 "zlib=true",
                 "filter=\"peakPicking vendor msLevel=1\"",
                 "filter=\"msLevel 1\""), fileConn)
    close(fileConn)
    ####files
    fileConn<-file("files.txt")
    writeLines(base::paste(path_to_raw,"\\",files_to_be_converted,sep=""), fileConn)
    close(fileConn)

    ###prepare arguments for msconvert
    arg <- base::paste("-f ",temp_path,"\\files.txt",
                 " -o ",temp_path,"\\mzXML",
                 " -c ",temp_path,"\\config.txt",sep="")

    ###run msconvert
    system2(path_to_msconvert, args = arg)

    ###move mzXMLs to original raw folder
    from <- temp_path
    to   <- path_to_raw
    path1 <- base::paste0(from,"\\mzXML")
    path2 <- base::paste0(to,"\\mzXML")

    dir.create(path2,showWarnings = F)

    for(f in base::gsub("\\.raw",".mzXML",files_to_be_converted))
    {
      ff::file.move(base::paste(path1,"\\",f,sep=""),path2)
    }

    ###remove temp msconvert folder
    r <- file.remove(c("files.txt","config.txt"))
  }

}

#' Prepare mzXML files for IceR workflow
#' @param path_to_mzXML Path to folder containing mzXML files
#' @param n_cores Numbers of CPU cores which should be used to perform conversion
#' @details Converts MS1 spectra in mzXML files into tables containing m/z, RT and intensity information per ion
#' @return Resulting ion tables are stored in a sub-directory (all_ion_lists) of the mzXML folder as .RData files
#' @export
mzxml_to_list <- function(path_to_mzXML,n_cores=2)
{
  #suppressWarnings(suppressMessages(library(doParallel,quietly = T)))
  i <- 0
  convert <- function(mzXMLfile,path_to_mzXML)
  {
    #suppressWarnings(suppressMessages(library(readMzXmlData,quietly = T)))

    data <- base::paste(path_to_mzXML,"/",mzXMLfile,sep="")
    pb <- tcltk::tkProgressBar(title = "Read mzXML",label=base::paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    ms <- readMzXmlData::readMzXmlFile(data)
    close(pb)

    #get total rowcount
    rowcount <- 0
    for(i in 1:length(ms))#
    {
      if(ms[[i]]$metaData$msLevel == 1)
      {
        rowcount <- rowcount + length(ms[[i]]$spectrum$mass)
      }
    }

    ###now extract data
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    dat <- data.table::as.data.table(matrix(ncol=3,nrow=rowcount))
    colnames(dat) <- c("m.z","RT","Intensity")
    sample <- mzXMLfile
    sample <- base::substr(sample,1,regexpr(".mzXML",sample)-1)
    dat$m.z <- as.numeric(dat$m.z)
    dat$RT <- as.numeric(dat$RT)
    dat$Intensity <- as.numeric(dat$Intensity)
    ind <- 1
    max <- length(ms)
    pb <- tcltk::tkProgressBar(title = "Extract data",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)

    for(i in 1:length(ms))
    {
      if(ms[[i]]$metaData$msLevel == 1)
      {
        start <- ind
        stop <- ind + length(ms[[i]]$spectrum$mass) - 1
        data.table::set(x = dat,i = start:stop,j = (1L),value = ms[[i]]$spectrum$mass)
        data.table::set(x = dat,i = start:stop,j = (2L),value = ms[[i]]$metaData$retentionTime/60)
        data.table::set(x = dat,i = start:stop,j = (3L),value = ms[[i]]$spectrum$intensity)
        ind <- ind + length(ms[[i]]$spectrum$mass)
      }
      tcltk::setTkProgressBar(pb, i, label = base::paste( round(i/max*100, 0),"% done (",i,"/",max,")",sep=""))
    }
    close(pb)
    save(dat,file=base::paste(path_to_mzXML,"/all_ion_lists/",sample,"_all_ions.RData",sep=""))
  }

  ##Step - Extract all ions per ms1 spectra

  mzXMLfiles <- list.files(path_to_mzXML)
  mzXMLfiles <- mzXMLfiles[which(grepl(".mzXML",mzXMLfiles))]

  dir.create(base::paste(path_to_mzXML,"/all_ion_lists",sep=""),showWarnings = F)

  setwd(base::paste(path_to_mzXML,"/all_ion_lists",sep=""))

  ###check if all ion list are already available and only generate those which are required
  available_all_ion_lists <- list.files()[which(grepl("\\.RData",list.files()))]

  missing_all_ion_lists <- mzXMLfiles[which(base::gsub("\\.mzXML","",mzXMLfiles) %not in% base::gsub("_all_ions\\.RData","",available_all_ion_lists))]
  if(length(missing_all_ion_lists)>0)
  {
    mzXMLfiles <- missing_all_ion_lists

    cl <- parallel::makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=mzXMLfiles) %dopar%
      {
        convert(i,path_to_mzXML)
      }
    parallel::stopCluster(cl)
  }

}

#' Extract MS1-spectra from raw TIMS-ToF Pro data
#' @param path_to_raw Path to folder containing d-files
#' @details Convert MS1 spectra with TIMS information in raw files of Bruker TIMS-ToF Pro Mass Spectrometers into tables containing m/z, RT, inverse ion mobility and intensity information per ion. This process is currently very slow and can take several days even for a small data set. It also requires enough space on the home drive of the PC as well as at least 50 Gb of memory.
#' @return Resulting ion tables are stored in a sub-directory (all_ion_lists) of the raw folder as .RData files
#' @export
convert_rawTIMS <- function(path_to_raw=NULL)
{
  ###extract spectra from raw text
  extract_spectra <- function(path,filename)
  {
    #library(data.table)
    con <- file(base::paste(path,filename,sep=""), "r")
    #read until spectrumList length is given
    while(T)
    {
      line = readLines(con, 1)
      if(grepl("    spectrumList \\(",line))
      {
        break
      }
    }
    #get number of total spectra
    n_spectra <- as.numeric(base::gsub("    spectrumList \\(| spectra):","",line[1]))

    #prepare table summarizing relevant information
    spectra <- base::as.data.frame(matrix(ncol=5,nrow=n_spectra))
    colnames(spectra) <- c("RT","1/KO","mz","int","num_ions")
    spectra$RT <- as.numeric(spectra$RT)
    spectra$`1/KO` <- as.numeric(spectra$`1/KO`)
    spectra$mz <- as.character(spectra$mz)
    spectra$int <- as.character(spectra$int)
    spectra$num_ions <- as.numeric(spectra$num_ions)
    count_spectra <- 0

    #jump to section where individual spectra are starting
    while(T)
    {
      line = readLines(con, 1)
      if(grepl("      spectrum:",line))
      {
        break
      }
    }
    #now read individual spectra
    print(base::paste(base::gsub(".txt","",filename),": Prepare TIMS-ToF data",sep=""))
    pb <- utils::txtProgressBar(min = 0,max = n_spectra,style = 3)
    for(i in 1:n_spectra)
    {
      #check if ion data is available
      if(i == 1)
      {
        line = readLines(con, 3)
        ion_count <- as.numeric(base::gsub("        defaultArrayLength: ","",line[3]))
      }else
      {
        line = readLines(con, 4)
        ion_count <- as.numeric(base::gsub("        defaultArrayLength: ","",line[4]))
      }


      if(ion_count > 0)
      {
        #read data for this spectrum
        line = readLines(con, 43)

        #increment number of available spectra with ions
        count_spectra <- count_spectra + 1

        #get RT and 1/K0
        RT <- as.numeric(base::gsub("            cvParam: scan start time, |, second","",line[32]))
        ims <- as.numeric(base::gsub("            cvParam: inverse reduced ion mobility, |, volt-second per square centimeter","",line[33]))

        #get m/z and intensities
        mz <- base::gsub(base::paste("          binary: \\[",ion_count,"\\] ",sep=""),"",line[40])
        int <- base::gsub(base::paste("          binary: \\[",ion_count,"\\] ",sep=""),"",line[43])

        #store data in table
        data.table::set(spectra,as.integer(count_spectra),as.integer(1:5),list(RT,ims,mz,int,ion_count))

      }else
      {
        #jump to the end of this spectrum
        line = readLines(con, 41)
      }
      utils::setTxtProgressBar(pb, i)

    }
    close(pb)

    spectra <- spectra[1:count_spectra,]

    close(con)

    ##prepare final table
    table_store <- base::as.data.frame(matrix(ncol=4,nrow=sum(spectra$num_ions)))
    colnames(table_store) <- c("RT","1/KO","mz","int")
    table_store$RT <- as.numeric(table_store$RT)
    table_store$`1/KO` <- as.numeric(table_store$`1/KO`)
    table_store$mz <- as.numeric(table_store$mz)
    table_store$int <- as.numeric(table_store$int)
    counter <- 1
    print(base::paste(base::gsub(".txt","",filename),": Finalize conversion of TIMS-ToF data",sep=""))
    pb <- utils::txtProgressBar(min = 0,max = nrow(spectra),style = 3)
    for(i in 1:nrow(spectra))
    {
      start <- counter
      end <- counter + spectra$num_ions[i] - 1
      mzs <- as.numeric(unlist(strsplit(spectra$mz[i]," ")))
      ints <- as.numeric(unlist(strsplit(spectra$int[i]," ")))
      data.table::set(table_store,as.integer(start:end),as.integer(1:4),list(rep(spectra$RT[i],spectra$num_ions[i]),
                                                                 rep(spectra$`1/KO`[i],spectra$num_ions[i]),
                                                                 mzs,
                                                                 ints))
      counter <- counter + spectra$num_ions[i]

      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    dir.create(base::paste(path,"all_ion_lists",sep=""),showWarnings = F)
    print(base::paste(base::gsub(".txt","",filename),": Store extracted spectra data",sep=""))
    save(table_store,file = base::paste(path,"all_ion_lists/",base::gsub(".txt","_all_ions.RData",filename),sep=""))

    rm(spectra,table_store)
    gc()
  }
  '%!in%' <- function(x,y)!('%in%'(x,y))
  '%not in%' <- function(x,y)!('%in%'(x,y))

  #suppressWarnings(suppressMessages(library(ff,quietly = T)))
  #suppressWarnings(suppressMessages(library(rChoiceDialogs,quietly = T)))

  if(is.null(path_to_raw))path_to_raw <- rChoiceDialogs::rchoose.dir(caption = "Select folder containing raw files")

  ###get home directory
  home_folder <- Sys.getenv("HOME")
  home_folder <- base::gsub("Documents","AppData/Local/Apps/",home_folder)

  ###find MSConvert folder
  folders <- list.dirs(path = home_folder, full.names = TRUE, recursive = F)
  folders <- folders[which(grepl("ProteoWizard",folders))]
  folders <- folders[length(folders)]

  raw_files <- list.files(path_to_raw)
  raw_files <- raw_files[which(grepl("\\.d",raw_files))]

  ###check which raw files still have to be converted
  text_available <- list.files(base::paste(path_to_raw,"/all_ion_lists",sep=""))

  files_to_be_converted <- raw_files[which(base::gsub("\\.d","",raw_files) %not in% base::gsub("_all_ions.RData","",text_available))]
  if(length(files_to_be_converted)>0)
  {
    ##find location of msconvert.exe
    if(file.exists(base::paste(folders,"\\msconvert.exe",sep="")))
    {
      path_to_msconvert <- base::paste(folders,"\\msconvert.exe",sep="")
    }else
    {
      print("Could not find msconvert.exe. Can be usually found in C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version")
      pb <- tcltk::tkProgressBar("Warning!",min = 0,max = 3,initial = 0,label = "Could not find msConvert.exe.",width = 500)

      counter <- 1
      label <- c("Typically located in C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version",
                 "Please specify location of msConvert.exe")
      while(T)
      {
        Sys.sleep(3)
        tcltk::setTkProgressBar(pb,value = counter,label = label[counter])
        counter <- counter + 1
        if(counter == 4)break
      }
      close(pb)
      path_to_msconvert <- file.choose()
    }

    ##prepare conversion for each file into raw text
    ###get user folder
    win_user_folder <- path.expand('~')

    ###create temporary folder
    dir.create(base::paste(win_user_folder,"\\temp_msconvert",sep=""),showWarnings = F)
    dir.create(base::paste(win_user_folder,"\\temp_msconvert\\temp",sep=""),showWarnings = F)

    ###create temporary config.txt and files.txt
    temp_path <- base::paste(win_user_folder,"\\temp_msconvert",sep="")
    setwd(temp_path)

    ####config
    fileConn<-file("config.txt")
    writeLines(c("text=true",
                 "64=true",
                 "zlib=true",
                 "filter=\"peakPicking vendor msLevel=1\"",
                 "filter=\"msLevel 1\""), fileConn)
    close(fileConn)

    ###now perform conversion using msconvert and subsequently by R for each file step by step

    for(i in 1:length(files_to_be_converted))
    {
      ####specify files.txt for msconvert containing the path to the current file
      fileConn<-file("files.txt")
      writeLines(base::paste(path_to_raw,"\\",files_to_be_converted[i],sep=""), fileConn)
      close(fileConn)

      ###prepare arguments for msconvert
      arg <- base::paste("-f ",temp_path,"\\files.txt",
                   " -o ",temp_path,"\\temp",
                   " -c ",temp_path,"\\config.txt",sep="")

      ###run msconvert
      if(file.exists(base::paste(temp_path,"\\temp\\",base::gsub(".d",".txt",files_to_be_converted[i]),sep="")) == F)#check if file is already converted to raw text
      {#if not use msconvert
        system2(path_to_msconvert, args = arg)
      }else
      {
        print("Converted txt file already available")
      }

      ###when finished, trigger conversion into spectra RData file
      if(file.exists(base::paste(temp_path,"\\temp\\all_ion_lists\\",base::gsub(".d","_all_ions.RData",files_to_be_converted[i]),sep="")) == F)#check if file is already converted to raw text
      {#if not use msconvert
        extract_spectra(path = base::paste(temp_path,"\\temp\\",sep=""),filename = base::gsub(".d",".txt",files_to_be_converted[i]))
      }else
      {
        print("Spectra were already extracted")
      }

      ###copy final spectra .RData file to destination and remove raw txt file
      from <- temp_path
      to   <- path_to_raw
      path1 <- base::paste0(from,"\\temp\\all_ion_lists")
      path2 <- base::paste0(to,"\\all_ion_lists")

      dir.create(path2,showWarnings = F)

      f <- base::gsub("\\.d","_all_ions.RData",files_to_be_converted[i])
      r <- ff::file.move(base::paste(path1,"\\",f,sep=""),path2)
      #remove raw txt file
      r <- file.remove(base::paste(temp_path,"\\temp\\",base::gsub(".d",".txt",files_to_be_converted[i]),sep=""))

    }
  }else
  {
    print("All files are already converted.")
  }
}

#' Perform alignment of pre-determined MS1-features by MaxQuant over proteomics samples either analyzed on an Orbitrap machine (e.g. Q-Exactive or Orbitrap Fusion) or a TIMS-ToF Pro
#' @param path_to_MaxQ_output Path to folder containing MaxQuant outputs (txt folder containing at least allpeptides.txt, evidence.txt, peptides.txt and proteinGroups.txt)
#' @param path_to_output Path to folder where IceR results should be stored
#' @param output_file_names_add IceR result name tag. By default IceR_analysis
#' @param mz_window Numeric value indicating maximal m/z deviation around a center. By default set to NA which indicates that the function determines automatically this parameter based on sd of m/z of identified peptides between samples.
#' @param min_mz_window Numberic value indicating how large the automatically determined m/z-window should be. By default set to 0.001 Da. Only required if m/z alignment window should be automatically determined. If set to NA, no minimal window size will be required.
#' @param RT_window Numeric value indicating maximal RT deviation around a center. By default set to NA which indicates that the function determines automatically this parameter based on sd of RT of identified peptides between samples.
#' @param min_RT_window Numberic value indicating how large the automatically determined RT-window should be. By default set to 1 min. Only required if RT alignment window should be automatically determined. If set to NA, no minimal window size will be required.
#' @param feature_mass_deviation_collapse Numeric value indicating which minimal mass deviation is required to distinguish IceR features. By default set to 0.002 Da. IceR features with overlapping RT-windows and mass differences smaller than the specified value are merged.
#' @param only_unmodified_peptides Boolean value indicating if only unmodified peptide sequences are used for alignment. By default set to F.
#' @param remove_contaminants Boolean value indicating if peptide features labeled as contaminants should be removed. By default set to T.
#' @param sample_list Character vector (raw file names) listing which samples should be aligned. By default all samples occuring in MaxQuant outputs are aligned.
#' @param align_unknown Boolean value indicating if only peptide features or also unsequenced features should be aligend over samples. By default set to F.
#' @param min_num_ions_collapse Numeric value indicating how many unsequenced MaxQuant features have to be at least detected over all samples to result in an IceR feature. Only required if align_unknown is set to T. By default set to 10.
#' @param MassSpec_mode String being either "Orbitrap" or "TIMSToF" specifying by which type of Mass Spectrometer the data was generated. By default it expects Thermo Orbitrap data.
#' @param IM_window Numeric value indicating maximal ion mobility (inverse K0) deviation around a center. By default set to NA which indicates that the function determines automatically this parameter based on sd of ion mobility of identified peptides between samples.
#' @param min_IM_window Numberic value indicating how large the automatically determined ion mobility-window (inverse K0) should be. By default set to 0.002. Only required if ion mobility alignment window should be automatically determined. If set to NA, no minimal window size will be required.
#' @param multiplicity Numeric value between 1 - 3 indicating if multiple channels are mixed on MS1-level. Label-free corresponds to 1 and SILAC to 2 (light and heavy) or 3 (light, medium, heavy)
#' @param SILAC_settings List of 3 string vectors named light, medium and heavy. Each string vector should indicate which heavy isotope of Lys and/or Arg was used for respective channel. String vectors can contain Lys2, Lys4, Lys6, Lys8, Arg6 and Arg10.
#' @details Performs the first steps of the IceR workflow: 1) Alignment window determination if not specified. 2) Alignment of MaxQuant features into IceR features. 3) Transfer of sequence information between MaxQuant features aligned into IceR features. 4) Extraction, modelling and prediction of RT- and m/z-correction factor per IceR feature and sample.
#' @return Outputs are stored in the sub-directory Temporary_files within specified output folder. MaxQuant allpeptides.txt and evidence.txt are converted to RData files. QC plots of estimated alignment windows as well as of random forest modesl and generalized additive models are stored in a QC_plots.pdf. Relevant QC data is stored in Feature_alignment_QC_data.RData. Aligned IceR features are stored in Features_aligned_merged.txt
#' @export
align_features <- function(path_to_MaxQ_output,path_to_output,align_unknown=F,output_file_names_add="IceR_analysis",mz_window=NA,min_mz_window = 0.001,RT_window=NA,min_RT_window=1,min_num_ions_collapse=10,feature_mass_deviation_collapse=0.002,only_unmodified_peptides=F,sample_list=NA,remove_contaminants=T,MassSpec_mode=c("Orbitrap","TIMSToF"),IM_window=NA,min_IM_window=0.002,multiplicity=c(1,2,3),SILAC_settings=list(light=c(""),medium=c("Lys4","Arg6"),heavy=c("Lys8","Arg10")))
{
  # path_to_MaxQ_output <- "D:\\Publication\\IceR\\test IceR example\\MaxQ"
  # path_to_output <- "D:\\Publication\\IceR\\test IceR example\\IceR_V1001"
  # align_unknown=F
  # output_file_names_add="IceR_analysis_V1001"
  # mz_window=NA
  # min_mz_window = 0.001
  # RT_window=NA
  # min_RT_window=1
  # min_num_ions_collapse=10
  # feature_mass_deviation_collapse=0.002
  # only_unmodified_peptides=F
  # sample_list=c("20200110_QE1_DDA_1H25_T5_E3_R1","20200110_QE1_DDA_1H25_T5_E3_R2","20200110_QE1_DDA_1H25_T5_E9_R1","20200110_QE1_DDA_1H25_T5_E9_R2")
  # remove_contaminants=T
  # MassSpec_mode="Orbitrap"
  # IM_window=NA
  # min_IM_window=0.002
  # multiplicity=1
  # SILAC_settings=list(light=c(""),medium=c("Lys4","Arg6"),heavy=c("Lys8","Arg10"))

  options(warn=-1)
  # suppressWarnings(suppressMessages(library(data.table,quietly = T)))
  # suppressWarnings(suppressMessages(library(stringr,quietly = T)))
  # suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
  # suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
  # suppressWarnings(suppressMessages(library(ff,quietly = T)))
  # suppressWarnings(suppressMessages(library(matrixStats,quietly = T)))
  use_mz_at_max_int_for_correction=F

  #check if multiplicity is 2 or 3 (e.g. SILAC samples)
  multiplicity <- multiplicity[1] #take first definition if a vector
  #if multiplicity is two, specification of light and heavy will be used
  #if multiplicity is three, specifications of light, medium and heavy will be used
  #SILAC_settings is a list of light, medium and heavy and each is a character vector
  #this character vector should contain one or multiple of the following names: Lys4, Lys8, Arg6, Arg10 indicating which heavy isotope was used per channel
  #can be also combinations like Lys4 and Arg6 for medium channel --> medium=c("Lys4","Arg6")
  #convert isotope AA labels into expected mass shifts
  if(multiplicity > 1)
  {
    if(multiplicity == 2)SILAC_settings$medium <- NULL #exclude medium channel in case of only 2 channels
    require_lys <- F
    require_arg <- F
    for(i in 1:length(SILAC_settings)) ##convert label names into expected mass shifts
    {
      current_mass_shift_Lys <- 0
      current_mass_shift_Arg <- 0

      if(any(SILAC_settings[[i]] == "Lys2"))current_mass_shift_Lys <- current_mass_shift_Lys + 1.994070
      if(any(SILAC_settings[[i]] == "Lys4"))current_mass_shift_Lys <- current_mass_shift_Lys + 4.025108
      if(any(SILAC_settings[[i]] == "Lys6"))current_mass_shift_Lys <- current_mass_shift_Lys + 6.020130
      if(any(SILAC_settings[[i]] == "Lys8"))current_mass_shift_Lys <- current_mass_shift_Lys + 8.014200
      if(any(SILAC_settings[[i]] == "Arg6"))current_mass_shift_Arg <- current_mass_shift_Arg + 6.020130
      if(any(SILAC_settings[[i]] == "Arg10"))current_mass_shift_Arg <- current_mass_shift_Arg + 10.00827
      SILAC_settings[[i]] <- list(labels = SILAC_settings[[i]],
                                  mass_shift_Lys = current_mass_shift_Lys,
                                  mass_shift_Arg = current_mass_shift_Arg)
      if(current_mass_shift_Lys != 0)require_lys <- T
      if(current_mass_shift_Arg != 0)require_arg <- T

    }
  }

  #Mass spec Mode can only be Orbitrap or TIMSToF
  MassSpec_mode <- MassSpec_mode[1]

  if(is.na(sample_list))sample_list<-NULL

  if(output_file_names_add != "")output_file_names_add <- base::paste("_",output_file_names_add,sep="")

  setwd(path_to_output)
  dir.create("Temporary_files")

  if(file.exists(base::paste("Temporary_files/Features_aligned_merged",output_file_names_add,".txt",sep="")))
  {
    print("Alignment already done")
  }else
  {
    if(file.exists("Temporary_files/allPeptides.RData"))
    {
      load("Temporary_files/allPeptides.RData")
      load("Temporary_files/evidence.RData")

      if(is.null(sample_list))
      {
        sample_list <- sort(unique(allpeptides$Raw.file))
      }

    }else
    {
      setwd(path_to_MaxQ_output)
      options(fftempdir = path_to_MaxQ_output)
      print("Read MaxQ results")
      temp <- utils::read.csv(file = "allPeptides.txt",sep='\t',nrows = 2044,header=T)
      tempclasses = sapply(temp, class)
      if(MassSpec_mode == "Orbitrap")
      {
        if(multiplicity == 1)
        {
          tempclasses["Intensity"] = "numeric"
          tempclasses[which(tempclasses == "logical" | tempclasses == "character")] <- "factor"
        }else
        {
          tempclasses["Intensity"] = "numeric"
          tempclasses["Intensity.L"] = "numeric"
          if(multiplicity == 3)tempclasses["Intensity.M"] = "numeric"
          tempclasses["Intensity.H"] = "numeric"
          tempclasses["Score"] = "numeric"
          if(any(names(tempclasses) == "PEP"))tempclasses["PEP"] = "factor"
          if(any(colnames(temp) == "Retention.Length..FWHM.")) tempclasses["Retention.Length..FWHM."] = "numeric"
          if(any(colnames(temp) == "Retention.length..FWHM.")) tempclasses["Retention.length..FWHM."] = "numeric"
          tempclasses[which(tempclasses == "logical" | tempclasses == "character")] <- "factor"
        }

      }else
      {
        if(multiplicity == 1)
        {
          tempclasses["m.z"] = "numeric"
          tempclasses["Mass"] = "numeric"
          tempclasses["Retention.time" ] = "numeric"
          tempclasses["Retention.length"] = "numeric"
          tempclasses["Ion.mobility.index"] = "numeric"
          tempclasses["Ion.mobility.index.length"] = "numeric"
          tempclasses[which(tempclasses == "logical" | tempclasses == "character")] <- "factor"
        }
      }
      allpeptides_save <- ff::read.csv.ffdf(file = "allPeptides.txt",sep='\t',VERBOSE = F,colClasses=tempclasses,next.rows = 100000)##read in data in chunks of 100000 rows

      allpeptides <- base::as.data.frame(allpeptides_save)
      allpeptides <- allpeptides[order(allpeptides$Mass),]

      if(multiplicity > 1)allpeptides <- allpeptides[which(allpeptides$Type != "ISO"),] ###remove ISO features as they are not adding anything at the moment

      #allpeptidessave <- allpeptides
      ###free some memory
      rm(allpeptides_save)
      gc()

      ##unify column names
      #RT
      if(length(which(colnames(allpeptides) == "Retention.Length"))==0) ###possibly upper and lower case problem
      {
        if(length(which(colnames(allpeptides) == "Retention.length"))>0) ###possibly upper and lower case problem
        {
          colnames(allpeptides)[which(colnames(allpeptides) == "Retention.length")] <- "Retention.Length"
        }else
        {
          allpeptides$Retention.Length = 0
        }
      }

      if(MassSpec_mode == "Orbitrap")allpeptides$Retention.Length <- allpeptides$Retention.Length/60 ##in Orbitrap allpeptides.txt the retention length is given in sec instead of min like in e.g. evidence.txt

      #MSMS scan numbers
      if(length(which(colnames(allpeptides) == "MSMS.Scan.Numbers"))==0)
      {
        if(length(which(colnames(allpeptides) == "MS.MS.scan.number"))>0)
        {
          colnames(allpeptides)[which(colnames(allpeptides) == "MS.MS.scan.number")] <- "MSMS.Scan.Numbers"
        }else
        {
          allpeptides$MSMS.Scan.Numbers = 0
        }
      }

      ###Some peptides are not available in allpeptides.txt thus we get this additional information from the evidence.txt file
      evidence <- utils::read.csv(file = "evidence.txt",sep='\t',header=T)
      if(any(colnames(evidence) == "MS.MS.Scan.Number"))colnames(evidence)[which(colnames(evidence) == "MS.MS.Scan.Number")] <- "MSMS.Scan.Numbers"
      if(any(colnames(evidence) == "MS.MS.scan.number"))colnames(evidence)[which(colnames(evidence) == "MS.MS.scan.number")] <- "MSMS.Scan.Numbers"
      if(any(colnames(evidence) == "MSMS.Scan.Numbers"))colnames(evidence)[which(colnames(evidence) == "MSMS.Scan.Numbers")] <- "MSMS.Scan.Numbers"
      if(any(colnames(evidence) == "Retention.length"))colnames(evidence)[which(colnames(evidence) == "Retention.length")] <- "Retention.Length"
      if(any(colnames(evidence) == "Ion.mobility.length"))colnames(evidence)[which(colnames(evidence) == "Ion.mobility.length")] <- "Ion.mobility.index.length"
      if(any(colnames(evidence) == "Number.of.scans"))colnames(evidence)[which(colnames(evidence) == "Number.of.scans")] <- "Number.of.frames"

      #load msms info in case of SILAC data
      if(multiplicity > 1)
      {
        temp <- utils::read.csv(file = "msms.txt",sep='\t',nrows = 2044,header=T)
        tempclasses = sapply(temp, class)
        tempclasses[which(tempclasses != "factor")] <- "factor"
        msms_save <- ff::read.csv.ffdf(file = "msms.txt",sep='\t',VERBOSE = F,colClasses=tempclasses,next.rows = 100000)##read in data in chunks of 100000 rows
        msms <- base::as.data.frame(msms_save)

        if(length(which(colnames(msms) == "Labeling.state"))==0) ###possibly upper and lower case problem
        {
          colnames(msms)[which(colnames(msms) == "Labeling.State")] <- "Labeling.state"
        }

        #match msms to evidence table
        match_order <- match(evidence$Best.MS.MS,msms$id)

        #add labeling state info from msms to evidence table
        evidence$Labeling.State <- as.numeric(as.character(msms$Labeling.state))[match_order]
      }


      print("Read MaxQ results finished")

      if(MassSpec_mode == "Orbitrap")
      {
        #in case of SILAC only use true identification per SILAC state as given (intensities could be requantified by MaxQ if setting enabled)
        if(multiplicity > 1)
        {
          #remove unclear labeling states
          sel <- which(is.na(evidence$Labeling.State) | evidence$Labeling.State == -1)
          if(length(sel) > 0)evidence <- evidence[-sel,]

          #remove intensities of silac channels which were not directly identified
          evidence$Intensity.H[which(evidence$Labeling.State == 0)] <- NA
          evidence$Intensity.M[which(evidence$Labeling.State == 0)] <- NA
          evidence$Intensity.H[which(evidence$Labeling.State == 1)] <- NA
          evidence$Intensity.L[which(evidence$Labeling.State == 1)] <- NA
          evidence$Intensity.L[which(evidence$Labeling.State == 2)] <- NA
          evidence$Intensity.M[which(evidence$Labeling.State == 2)] <- NA
        }

        add_data <- base::as.data.frame(matrix(ncol=ncol(allpeptides),nrow=nrow(evidence),NA))
        colnames(add_data) <- colnames(allpeptides)
        match_col_names <- match(colnames(allpeptides),colnames(evidence))
        data.table::set(add_data,j=as.integer(which(!is.na(match_col_names))),value = evidence[,match_col_names[which(!is.na(match_col_names))]])

        ####now find rows (peptides + modification + sample) which are not yet present in allpeptides.txt
        temp1 <- base::paste(add_data$Sequence,add_data$Modifications,add_data$Raw.file,sep="_")
        temp2 <- base::paste(allpeptides$Sequence,allpeptides$Modifications,allpeptides$Raw.file,sep="_")
        add_data <- add_data[which(temp1 %not in% temp2),]
        allpeptides <- rbind(allpeptides,add_data)
        rm(temp1,temp2,add_data)
        gc()
      }

      if(MassSpec_mode == "TIMSToF")
      {
        add_data <- base::as.data.frame(matrix(ncol=ncol(allpeptides),nrow=nrow(evidence),NA))
        colnames(add_data) <- colnames(allpeptides)
        match_col_names <- match(colnames(allpeptides),colnames(evidence))
        data.table::set(add_data,j=as.integer(which(!is.na(match_col_names))),value = evidence[,match_col_names[which(!is.na(match_col_names))]])

        ###in TIMS allpeptides.txt column sequence and modification is missing thus take evidence table completeley
        add_data$Sequence <- evidence$Sequence
        add_data$Modifications <- evidence$Modifications
        add_data$Proteins <- evidence$Proteins
        add_data$Score <- evidence$Score

        allpeptides$Sequence <- ""
        allpeptides$Modifications <- ""
        allpeptides$Proteins <- ""
        allpeptides$Score <- NA
        allpeptides$Score <- as.numeric(allpeptides$Score)

        allpeptides <- rbind(allpeptides,add_data)
        rm(add_data)
        gc()

        ###convert Ion.mobility.index into 1/K0
        fit <- stats::lm(evidence$X1.K0~evidence$Ion.mobility.index)
        allpeptides$inverse_K0 <- (allpeptides$Ion.mobility.index*fit$coefficients[2])+fit$coefficients[1]
        allpeptides$inverse_K0_length <- allpeptides$Ion.mobility.index.length/1000

      }

      #if in SILAC mode, treat each channel as individual sample
      if(multiplicity > 1)
      {
        temp_allpeptides <- base::as.data.frame(matrix(ncol=ncol(allpeptides)+4,nrow=0))
        colnames(temp_allpeptides) <- c(colnames(allpeptides),"m.z_SILAC","m.z_SILAC_uncalibrated","Lys_count","Arg_count")

        #change intensities set to 0 to NA
        allpeptides$Intensity.L[allpeptides$Intensity.L == 0] <- NA
        allpeptides$Intensity.M[allpeptides$Intensity.M == 0] <- NA
        allpeptides$Intensity.H[allpeptides$Intensity.H == 0] <- NA

        #Now expand every evidence row to an individual row per SILAC state where intensity != NA
        for(s in unique(allpeptides$Raw.file))
        {
          for(label in names(SILAC_settings))
          {
            if(label == "light")sel <- which(allpeptides$Raw.file == s & !is.na(allpeptides$Intensity.L))
            if(label == "medium")sel <- which(allpeptides$Raw.file == s & !is.na(allpeptides$Intensity.M))
            if(label == "heavy")sel <- which(allpeptides$Raw.file == s & !is.na(allpeptides$Intensity.H))

            #set intensity column to observed intensity in respective channel
            temp <- allpeptides[sel,]
            if(label == "light")temp$Intensity <- temp$Intensity.L
            if(label == "medium")temp$Intensity <- temp$Intensity.M
            if(label == "heavy")temp$Intensity <- temp$Intensity.H

            #set observed m/z at max intensity
            if(label == "light")temp$Max.intensity.m.z.0 <- temp$Max.intensity.m.z.0
            if(label == "medium")temp$Max.intensity.m.z.0 <- temp$Max.intensity.m.z.1
            if(label == "heavy")temp$Max.intensity.m.z.0 <- temp$Max.intensity.m.z.2

            #get number of Lys and Arg per feature
            count_AA_K <- ifelse(temp$Sequence != " ",stringr::str_count(temp$Sequence, "K"),NA)
            count_AA_R <- ifelse(temp$Sequence != " ",stringr::str_count(temp$Sequence, "R"),NA)
            if(any(colnames(temp) == "Lys.Count"))
            {
              temp$Lys_count <- ifelse(!is.na(temp$Lys.Count),temp$Lys.Count,count_AA_K)
            }else
            {
              temp$Lys_count <- count_AA_K
            }
            if(any(colnames(temp) == "Arg.Count"))
            {
              temp$Arg_count <- ifelse(!is.na(temp$Arg.Count),temp$Arg.Count,count_AA_R)
            }else
            {
              temp$Arg_count <- count_AA_R
            }

            #check if m/z shifts come from only Lys, only Arg or both
            #if only coming from Lys then only keep features for which Lys count is known
            #if only coming from Arg then only keep features for which Arg count is known
            #in case of both, require Lys and Arg count to be known
            #further filter for features with at least one of each required AA otherwise L, M and H channels can not be distinguished
            if(require_lys == T & require_arg == T)
            {
              temp <- temp[which(!is.na(temp$Lys_count) & temp$Lys_count > 0 | !is.na(temp$Arg_count) & temp$Arg_count > 0),]
            }else if(require_lys == T)
            {
              temp <- temp[which(!is.na(temp$Lys_count) & temp$Lys_count > 0),]
            }else if(require_arg == T)
            {
              temp <- temp[which(!is.na(temp$Arg_count) & temp$Arg_count > 0),]
            }
            temp$Lys_count[is.na(temp$Lys_count)] <- 0
            temp$Arg_count[is.na(temp$Arg_count)] <- 0

            #adjust m/z to m/z including expected mass shift
            temp$m.z_SILAC <- temp$m.z + ((SILAC_settings[[label]]$mass_shift_Lys*temp$Lys_count)/temp$Charge) + ((SILAC_settings[[label]]$mass_shift_Arg*temp$Arg_count)/temp$Charge)
            temp$m.z_SILAC_uncalibrated <- temp$Uncalibrated.m.z + + ((SILAC_settings[[label]]$mass_shift_Lys*temp$Lys_count)/temp$Charge) + ((SILAC_settings[[label]]$mass_shift_Arg*temp$Arg_count)/temp$Charge)

            #adjust raw file name
            temp$Raw.file <- base::paste(s,"_Channel_",label,sep="")

            temp_allpeptides <- rbind(temp_allpeptides,temp)
          }
        }

        allpeptides <- temp_allpeptides
      }

      ###if no specific sample list is defined which should be used, use all raw files
      if(is.null(sample_list))
      {
        sample_list <- sort(unique(allpeptides$Raw.file))
      }

      ###clean allpeptides table to reduce required memory space
      if(MassSpec_mode == "Orbitrap")
      {
        if(multiplicity == 1)
        {
          if(any(colnames(allpeptides) == "Resolution"))
          {
            allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","Mass","Uncalibrated.m.z","Resolution","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity")]
          }else
          {
            allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","Mass","Uncalibrated.m.z","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity")]
            allpeptides$Resolution <- 0
            print("MS resolution seems to be missing in MaxQ outputs !!!")
          }
        }else if(multiplicity > 1)
        {
          if(any(colnames(allpeptides) == "Resolution"))
          {
            allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","m.z_SILAC","Mass","Uncalibrated.m.z","m.z_SILAC_uncalibrated","Resolution","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity","Lys_count","Arg_count")]
          }else
          {
            allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","m.z_SILAC","Mass","Uncalibrated.m.z","m.z_SILAC_uncalibrated","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity","Lys_count","Arg_count")]
            allpeptides$Resolution <- 0
            print("MS resolution seems to be missing in MaxQ outputs !!!")
          }

        }

      }
      if(MassSpec_mode == "TIMSToF")
      {
        allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","Mass","Retention.time","Retention.Length","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity","inverse_K0","inverse_K0_length")]
        allpeptides$Resolution <- 0 ###no resolution information available
      }

      if(remove_contaminants == T)exclude <- which(grepl("CON",allpeptides$Proteins) | grepl("REV",allpeptides$Proteins))
      if(remove_contaminants == F)exclude <- which(grepl("REV",allpeptides$Proteins))
      if(length(exclude) > 0)allpeptides <- allpeptides[-exclude,] ###remove potential reverse and contaminant peptides


      ###add calibrated RT to all peptides
      if(multiplicity == 1)
      {
        temp_evidence <- evidence[,c("Sequence","Raw.file","Charge","Modifications","Calibrated.retention.time")]
        temp_evidence <- stats::aggregate(temp_evidence$Calibrated.retention.time,list(Sequence=evidence$Sequence,Raw.file=evidence$Raw.file,Charge=evidence$Charge,Modifications=evidence$Modifications),FUN=mean,na.rm=T)
        colnames(temp_evidence)[5] <- "Calibrated.retention.time"
        temp <- dplyr::left_join(allpeptides,temp_evidence,by=c("Sequence"="Sequence","Raw.file"="Raw.file","Charge"="Charge","Modifications"="Modifications"))

        missing_RT_calibration_indices <- which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")
        if(length(missing_RT_calibration_indices)>0)
        {
          for(ind in missing_RT_calibration_indices)
          {
            sub <- temp[which(temp$Sequence == temp$Sequence[ind] & temp$Modifications == temp$Modifications[ind]),]
            temp$Calibrated.retention.time[ind] <- stats::median(sub$Calibrated.retention.time,na.rm=T)
          }
        }

        ###if any peptide feature still doesn?t have a valid calibrated RT --> just use observed RT
        temp$Calibrated.retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")] <- temp$Retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")]
        allpeptides <- temp
      }else if(multiplicity > 1)
      {
        temp_evidence <- evidence[,c("Sequence","Raw.file","Charge","Modifications","Calibrated.retention.time")]
        temp_evidence <- stats::aggregate(temp_evidence$Calibrated.retention.time,list(Sequence=evidence$Sequence,Raw.file=evidence$Raw.file,Charge=evidence$Charge,Modifications=evidence$Modifications),FUN=mean,na.rm=T)
        colnames(temp_evidence)[5] <- "Calibrated.retention.time"
        allpeptides$Raw.file_orig <- base::gsub("_Channel_light|_Channel_medium|_Channel_heavy","",allpeptides$Raw.file)
        temp <- dplyr::left_join(allpeptides,temp_evidence,by=c("Sequence"="Sequence","Raw.file_orig"="Raw.file","Charge"="Charge","Modifications"="Modifications"))

        missing_RT_calibration_indices <- which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")
        if(length(missing_RT_calibration_indices)>0)
        {
          for(ind in missing_RT_calibration_indices)
          {
            sub <- temp[which(temp$Sequence == temp$Sequence[ind] & temp$Modifications == temp$Modifications[ind]),]
            temp$Calibrated.retention.time[ind] <- stats::median(sub$Calibrated.retention.time,na.rm=T)
          }
        }

        ###if any peptide feature still doesn?t have a valid calibrated RT --> just use observed RT
        temp$Calibrated.retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")] <- temp$Retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")]
        allpeptides <- temp
      }

      #determine isotope corrected m/z
      if(MassSpec_mode == "Orbitrap" & multiplicity == 1)
      {
        ###calculate deviations between samples
        add <- base::as.data.frame(matrix(nrow=nrow(allpeptides),ncol=3))
        colnames(add) <- c("isotope_at_max_int","isotope_corrected_mz_at_max_int","delta_mz_to_mz_at_max_int")
        add[,1] <- as.numeric(add[,1])
        add[,2] <- as.numeric(add[,2])
        add[,3] <- as.numeric(add[,3])

        isotopes <- matrix(nrow=nrow(allpeptides),ncol=4)
        isotopes[,1] <- allpeptides$m.z
        isotopes[,2] <- ((allpeptides$m.z*allpeptides$Charge)+(1*1.002054))/allpeptides$Charge
        isotopes[,3] <- ((allpeptides$m.z*allpeptides$Charge)+(2*1.002054))/allpeptides$Charge
        isotopes[,4] <- ((allpeptides$m.z*allpeptides$Charge)+(3*1.002054))/allpeptides$Charge

        max <- nrow(allpeptides)
        pb <- tcltk::tkProgressBar(title = "Determine deviations of true m/z from observed m/z",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        for(i in 1:nrow(allpeptides))
        {
          if(!is.na(allpeptides$Max.intensity.m.z.0[i]))
          {
            mz_max_int <- allpeptides$Max.intensity.m.z.0[i]
            deltas <- abs(isotopes[i,] - mz_max_int)
            closest_iso <- which(deltas == min(deltas,na.rm=T))-1
            if(closest_iso > 0) ###+1 or +2 or +3 isotope shows highest intensity
            { ##correct back to isotope +0 m/z
              iso_mult <- closest_iso*1.002054
              mz_max_int <- ((mz_max_int*allpeptides$Charge[i])-(iso_mult))/allpeptides$Charge[i]
            }
            delta_mz_to_mz_at_max_int <- mz_max_int - allpeptides$m.z[i]
            data.table::set(add,as.integer(i),as.integer(1:3),value=as.list(c(closest_iso,mz_max_int,delta_mz_to_mz_at_max_int)))
          }

          updatecounter <- updatecounter + 1
          if(updatecounter >= 100)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)

        allpeptides <- cbind(allpeptides,add)

        ###evidence table

        add <- base::as.data.frame(matrix(nrow=nrow(evidence),ncol=3))
        colnames(add) <- c("isotope_at_max_int","isotope_corrected_mz_at_max_int","delta_mz_to_mz_at_max_int")
        add[,1] <- as.numeric(add[,1])
        add[,2] <- as.numeric(add[,2])
        add[,3] <- as.numeric(add[,3])

        isotopes <- matrix(nrow=nrow(evidence),ncol=4)
        isotopes[,1] <- evidence$m.z
        isotopes[,2] <- ((evidence$m.z*evidence$Charge)+(1*1.002054))/evidence$Charge
        isotopes[,3] <- ((evidence$m.z*evidence$Charge)+(2*1.002054))/evidence$Charge
        isotopes[,4] <- ((evidence$m.z*evidence$Charge)+(3*1.002054))/evidence$Charge

        max <- nrow(evidence)
        pb <- tcltk::tkProgressBar(title = "Determine deviations of true m/z from observed m/z",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0

        for(i in 1:nrow(evidence))
        {
          if(!is.na(evidence$Max.intensity.m.z.0[i]))
          {
            mz_max_int <- evidence$Max.intensity.m.z.0[i]
            deltas <- abs(isotopes[i,] - mz_max_int)
            closest_iso <- which(deltas == min(deltas,na.rm=T))-1
            if(closest_iso > 0) ###+1 or +2 or +3 isotope shows highest intensity
            { ##correct back to isotope +0 m/z
              iso_mult <- closest_iso*1.002054
              mz_max_int <- ((mz_max_int*evidence$Charge[i])-(iso_mult))/evidence$Charge[i]
            }
            delta_mz_to_mz_at_max_int <- mz_max_int - evidence$m.z[i]
            data.table::set(add,as.integer(i),as.integer(1:3),value=as.list(c(closest_iso,mz_max_int,delta_mz_to_mz_at_max_int)))
          }

          updatecounter <- updatecounter + 1
          if(updatecounter >= 100)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)

        evidence <- cbind(evidence,add)
      }

      ###finally save prepared data
      setwd(path_to_output)

      save(allpeptides,file = "Temporary_files/allPeptides.RData")
      save(evidence,file = "Temporary_files/evidence.RData")
    }

    QC_data <- list() ##here relevant qc data is stored and finally saved as RData which can be used for re-generating plots

    grDevices::pdf(base::paste("Temporary_files/QC_plots",output_file_names_add,".pdf",sep=""))

    RT_calibration <- T
    mz_calibration <- T

    if(MassSpec_mode == "TIMSToF")
    {
      IM_calibration <- T
    }else
    {
      IM_calibration <- F
    }

    if(MassSpec_mode == "Orbitrap")QC_data[["MaxQ_calibrations"]] <- evidence[,c("Retention.time.calibration","Uncalibrated...Calibrated.m.z..Da.")]

    allpeptides <- allpeptides[order(allpeptides$m.z),]
    rownames(allpeptides) <- c(1:nrow(allpeptides))

    if(only_unmodified_peptides == T)
    {
      identified_ions <- subset(allpeptides,Sequence != " " & Sequence != "" & Modifications == "Unmodified")
    }else
    {
      identified_ions <- subset(allpeptides,Sequence != " " & Sequence != "")
    }

    identified_ions <- identified_ions[order(identified_ions$Sequence,identified_ions$Intensity,decreasing = T),]

    identified_ions$Sequence <- as.character(identified_ions$Sequence)
    unique_peptides <- unique(identified_ions$Sequence)

    if(MassSpec_mode == "Orbitrap")
    {
      windows <- base::as.data.frame(matrix(ncol=11,nrow=nrow(identified_ions)))
      colnames(windows) <- c("Peptide","Mod","Charge","mean_m.z","sd_m.z","mean_RT","sd_RT","num_ions","mean_RT_calibrated","sd_RT_calibrated","sd_m.z_uncalibrated")

      windows$Peptide <- as.character(windows$Peptide)
      windows$Mod <- as.character(windows$Mod)
      windows$Charge <- as.numeric(windows$Charge)
      windows$mean_m.z <- as.numeric(windows$mean_m.z)
      windows$sd_m.z <- as.numeric(windows$sd_m.z)
      windows$mean_RT <- as.numeric(windows$mean_RT)
      windows$sd_RT <- as.numeric(windows$sd_RT)
      windows$num_ions <- as.numeric(windows$num_ions)
      windows$mean_RT_calibrated <- as.numeric(windows$mean_RT_calibrated)
      windows$sd_RT_calibrated <- as.numeric(windows$sd_RT_calibrated)
      windows$sd_m.z_uncalibrated <- as.numeric(windows$sd_m.z_uncalibrated)
    }

    if(MassSpec_mode == "TIMSToF")
    {
      windows <- base::as.data.frame(matrix(ncol=13,nrow=nrow(identified_ions)))
      colnames(windows) <- c("Peptide","Mod","Charge","mean_m.z","sd_m.z","mean_RT","sd_RT","num_ions","mean_RT_calibrated","sd_RT_calibrated","sd_m.z_uncalibrated","mean_IM","sd_IM")

      windows$Peptide <- as.character(windows$Peptide)
      windows$Mod <- as.character(windows$Mod)
      windows$Charge <- as.numeric(windows$Charge)
      windows$mean_m.z <- as.numeric(windows$mean_m.z)
      windows$sd_m.z <- as.numeric(windows$sd_m.z)
      windows$mean_RT <- as.numeric(windows$mean_RT)
      windows$sd_RT <- as.numeric(windows$sd_RT)
      windows$num_ions <- as.numeric(windows$num_ions)
      windows$mean_RT_calibrated <- as.numeric(windows$mean_RT_calibrated)
      windows$sd_RT_calibrated <- as.numeric(windows$sd_RT_calibrated)
      windows$sd_m.z_uncalibrated <- as.numeric(windows$sd_m.z_uncalibrated)
      windows$mean_IM <- as.numeric(windows$mean_IM)
      windows$sd_IM <- as.numeric(windows$sd_IM)
    }

    #if multiplicity > 1 then only use light channel for the alignment parameter determination as M and H will have exact same m/z and RT
    if(multiplicity > 1)
    {
      identified_ions <- identified_ions[!duplicated(identified_ions[,c("Sequence","Raw.file_orig","Charge","Modifications")]),]
    }

    ###remove outlier samples where RT of a peptide is very different from all other samples
    outlier_RT_deviation <- (max(allpeptides$Retention.time,na.rm=T)/100)*5 ###deviation should not be larger than 5 % of the total chromatographic retention length

    max <- length(unique_peptides)
    pb <- tcltk::tkProgressBar(title = "Determine matching parameter windows",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    start <- 1
    max_ind <- nrow(identified_ions)
    count_features <- 0
    for(i in 1:length(unique_peptides))
    {
      ind <- start
      while(T)
      {
        ind = ind + 100
        if(ind < max_ind)
        {
          if(identified_ions[ind,"Sequence"] != unique_peptides[i])
          {
            break
          }
        }else
        {
          ind <- max_ind
          break
        }

      }
      sub <- identified_ions[start:ind,]

      ####reset for next start index
      bins <- (ind-start)/100
      start <- ifelse(bins==1,start,start+((bins-1)*100))

      sub <- sub[which(sub$Sequence == unique_peptides[i]),]
      if(length(which(!is.na(sub$Intensity)))>0)sub<-sub[which(!is.na(sub$Intensity)),]
      sub <- sub[!duplicated(sub$Raw.file),]

      for(c in unique(sub$Charge))
      {
        sub1 <- sub[which(sub$Charge == c),]
        for(m in unique(sub1$Modifications))
        {
          sub1 <- sub[which(sub$Charge == c & sub$Modifications == m),]
          count_features <- count_features + 1

          # if(multiplicity > 1)
          # {
          #   sub1$m.z <- sub1$Max.intensity.m.z.0 - (sub1$m.z_SILAC-sub1$m.z)
          # }

          median_RT <- stats::median(sub1$Retention.time,na.rm=T)
          remove <- which(sub1$Retention.time > median_RT + outlier_RT_deviation)
          if(length(remove)>0)sub1 <- sub1[-remove,]
          ###determine mean m/z at peak maximum and mean RT
          if(MassSpec_mode == "Orbitrap")
          {
            data.table::set(windows,as.integer(count_features),as.integer(c(3:11)),value=as.list(c(c,
                                                                                       mean(sub1$m.z,na.rm=T),
                                                                                       stats::sd(sub1$m.z,na.rm=T),
                                                                                       mean(sub1$Retention.time,na.rm=T),
                                                                                       stats::sd(sub1$Retention.time,na.rm=T),
                                                                                       nrow(sub1),mean(sub1$Calibrated.retention.time,na.rm=T),
                                                                                       stats::sd(sub1$Calibrated.retention.time,na.rm=T),
                                                                                       stats::sd(sub1$isotope_corrected_mz_at_max_int,na.rm=T))))
          }
          if(MassSpec_mode == "TIMSToF")
          {
            data.table::set(windows,as.integer(count_features),as.integer(c(3:13)),value=as.list(c(c,
                                                                                       mean(sub1$m.z,na.rm=T),
                                                                                       stats::sd(sub1$m.z,na.rm=T),
                                                                                       mean(sub1$Retention.time,na.rm=T),
                                                                                       stats::sd(sub1$Retention.time,na.rm=T),
                                                                                       nrow(sub1),mean(sub1$Calibrated.retention.time,na.rm=T),
                                                                                       stats::sd(sub1$Calibrated.retention.time,na.rm=T),
                                                                                       stats::sd(sub1$isotope_corrected_mz_at_max_int,na.rm=T),
                                                                                       mean(sub1$inverse_K0,na.rm=T),
                                                                                       stats::sd(sub1$inverse_K0,na.rm=T))))
          }
          data.table::set(windows,as.integer(count_features),as.integer(1:2),value=as.list(c(unique_peptides[i],m)))
        }

      }

      updatecounter <- updatecounter + 1
      if(updatecounter >= 10)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- lubridate::seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    windows <- windows[1:count_features,]

    #generate some QC plots indicating variation of RT and m/z (and IM) for identified features across samples
    if(multiplicity == 1)graphics::boxplot(windows$sd_m.z,main="Standard deviation of peptide features m/z",outline=F)
    graphics::boxplot(windows$sd_RT,main="Standard deviation of peptide features RT",outline=F)
    if(MassSpec_mode == "TIMSToF")graphics::boxplot(windows$sd_IM,main="Standard deviation of peptide features inverse K0",outline=F)
    if(is.na(mz_window))
    {
      ###based on boxplots
      if(!is.na(min_mz_window))
      {
        if(grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5]>min_mz_window)
        {
          borders_m.z <- c(-grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5],grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5]) ###upper whisker
        }else
        {
          borders_m.z <- c(-min_mz_window,min_mz_window)
        }
      }else
      {
        borders_m.z <- c(-grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5],grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5]) ###upper whisker
      }


    }else
    {
      ###used defined parameter
      borders_m.z <- c(-mz_window,mz_window)
    }
    if(is.na(RT_window))
    {
      ###based on boxplots
      if(!is.na(min_RT_window))
      {
        if(grDevices::boxplot.stats(windows$sd_RT)$`stats`[5]>min_RT_window)
        {
          borders_RT <- c(-grDevices::boxplot.stats(windows$sd_RT)$`stats`[5],grDevices::boxplot.stats(windows$sd_RT)$`stats`[5]) ##75% quantil based on boxplots
        }else
        {
          borders_RT <- c(-min_RT_window,min_RT_window)
        }
      }else
      {
        borders_RT <- c(-grDevices::boxplot.stats(windows$sd_RT)$`stats`[5],grDevices::boxplot.stats(windows$sd_RT)$`stats`[5]) ##75% quantil based on boxplots
      }

    }else
    {
      ###use defined parameter
      borders_RT <- c(-RT_window,RT_window)
    }

    if(MassSpec_mode == "TIMSToF")
    {
      if(is.na(IM_window))
      {
        ###based on boxplots
        if(!is.na(min_IM_window))
        {
          if(grDevices::boxplot.stats(windows$sd_IM)$`stats`[5]>min_IM_window)
          {
            borders_IM <- c(-grDevices::boxplot.stats(windows$sd_IM)$`stats`[5],grDevices::boxplot.stats(windows$sd_IM)$`stats`[5]) ##75% quantil based on boxplots
          }else
          {
            borders_IM <- c(-min_IM_window,min_IM_window)
          }
        }else
        {
          borders_IM <- c(-grDevices::boxplot.stats(windows$sd_IM)$`stats`[5],grDevices::boxplot.stats(windows$sd_IM)$`stats`[5]) ##75% quantil based on boxplots
        }

      }else
      {
        ###use defined parameter
        borders_IM <- c(-IM_window,IM_window)
      }
    }

    borders_RT_use_save = borders_RT

    if(MassSpec_mode == "Orbitrap")QC_data[["Feature_alignment_windows"]] <- list(RT_window=borders_RT,mz_window=borders_m.z)
    if(MassSpec_mode == "TIMSToF")QC_data[["Feature_alignment_windows"]] <- list(RT_window=borders_RT,mz_window=borders_m.z,IM_window=borders_IM)


    ###determine for how many windows we expect an overlap of ions for different peptide sequences based on the chosen parameters
    windows <- windows[order(windows$mean_m.z),]
    max_i <- nrow(windows)
    windows$overlap <- 0

    max <- nrow(windows)
    pb <- tcltk::tkProgressBar(title = "Determine overlapping peptide features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    if(MassSpec_mode == "Orbitrap")
    {
      for(i in 1:nrow(windows))
      {
        start <- ifelse((i - 5) < 1,1,(i - 5))
        end <- ifelse((i + 5) > max_i,max_i,(i + 5))
        sub <- windows[start:end,]
        sub <- sub[which(sub$mean_m.z >= windows$mean_m.z[i] + borders_m.z[1] & sub$mean_m.z <= windows$mean_m.z[i] + borders_m.z[2] & sub$mean_RT >= windows$mean_RT[i] + borders_RT[1] & sub$mean_RT <= windows$mean_RT[i] + borders_RT[2] & sub$Charge == windows$Charge[i]),]

        data.table::set(windows,as.integer(i),12L,value=(nrow(sub)-1))

        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }
    if(MassSpec_mode == "TIMSToF")
    {
      for(i in 1:nrow(windows))
      {
        start <- ifelse((i - 5) < 1,1,(i - 5))
        end <- ifelse((i + 5) > max_i,max_i,(i + 5))
        sub <- windows[start:end,]
        sub <- sub[which(sub$mean_m.z >= windows$mean_m.z[i] + borders_m.z[1] & sub$mean_m.z <= windows$mean_m.z[i] + borders_m.z[2] &
                           sub$mean_RT >= windows$mean_RT[i] + borders_RT[1] & sub$mean_RT <= windows$mean_RT[i] + borders_RT[2] &
                           sub$mean_IM >= windows$mean_IM[i] + borders_IM[1] & sub$mean_IM <= windows$mean_IM[i] + borders_IM[2] &
                           sub$Charge == windows$Charge[i]),]

        data.table::set(windows,as.integer(i),12L,value=(nrow(sub)-1))

        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }
    close(pb)

    p <- Barplots(plyr::count(windows$overlap)[,2],main="Number of expected feature overlaps",ylab="Count",Name = plyr::count(windows$overlap)[,1],xlab = "Overlaps",AvgLine = F)
    QC_data[["Alignment_deviations_overlap"]] <- windows

    #####Perform alignment of peptide features

    allpeptides <- allpeptides[order(allpeptides$m.z),]
    rownames(allpeptides) <- c(1:nrow(allpeptides))
    allpeptides$Sequence <- as.character(allpeptides$Sequence)
    allpeptides$Modifications <- as.character(allpeptides$Modifications)
    allpeptides$Proteins <- base::gsub(";","|",allpeptides$Proteins)
    allpeptides$Raw.file <- as.character(allpeptides$Raw.file)
    ###presubset all ions based on charge, here use only rows of allions which are coming from samples which should be also used
    ####also create a library of ions per charge state per mz_window of 0.5 Da
    print("Prepare for peptide feature matching - Indexing all ions")
    rownames(allpeptides) <- c(1:nrow(allpeptides))
    allpeptides_frag <- list()
    allpeptides_frag_indices_per_mz_window <- list()

    for(i in 1:max(allpeptides$Charge))
    {
      if(length(which(allpeptides$Charge == i))>0)
      {
        allpeptides_frag[[i]] <- allpeptides[which(allpeptides$Charge == i),]
        min_max_mz <- c(floor(min(allpeptides_frag[[i]]$m.z)),ceiling(max(allpeptides_frag[[i]]$m.z)))
        indices <- matrix(ncol=4,nrow=min_max_mz[2]-min_max_mz[1])
        indices[,1] <- min_max_mz[1]:(min_max_mz[2]-1)
        indices[,2] <- (min_max_mz[1]+1):(min_max_mz[2])
        for(j in 1:nrow(indices))
        {
          indx <- which(allpeptides_frag[[i]]$m.z >= indices[j,1] & allpeptides_frag[[i]]$m.z <= indices[j,2])
          if(length(indx) > 0)
          {
            indices[j,3:4] <- c(min(indx),max(indx))
          }else
          {
            indices[j,3:4] <- c(NA,NA)
          }

        }
        allpeptides_frag_indices_per_mz_window[[i]] <- indices
      }else
      {
        allpeptides_frag[[i]] <- NA
        allpeptides_frag_indices_per_mz_window[[i]] <- NA
      }

    }
    rm(indices)
    gc()

    temp_data <- base::as.data.frame(matrix(ncol=1,nrow=nrow(allpeptides))) ###here already used ions will be marked
    temp_data[,1] <- as.numeric(temp_data[,1])

    temp <- subset(allpeptides,Sequence != " " & Sequence != "")

    if(MassSpec_mode == "Orbitrap")
    {
      if(multiplicity == 1)
      {
        features <- base::as.data.frame(matrix(ncol=21,nrow=length(unique(base::paste(temp$Sequence,temp$Charge,temp$Modifications)))))
        colnames(features) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")

        features <- data.table::as.data.table(features)
        features$Feature_name <- as.character(features$Feature_name)
        features$Mass <- as.numeric(features$Mass)
        features$RT <- as.numeric(features$RT)
        features$Sequence <- as.character(features$Sequence)
        features$Protein <- as.character(features$Protein)
        features$MSMS.Scan.Numbers <- as.character(features$MSMS.Scan.Numbers)
        features$m.z_range_min <- as.numeric(features$m.z_range_min)
        features$m.z_range_max <- as.numeric(features$m.z_range_max)
        features$RT_range_min <- as.numeric(features$RT_range_min)
        features$RT_range_max <- as.numeric(features$RT_range_max)
        features$ion_indices <- as.character(features$ion_indices)
        features$count_ion_indices <- as.numeric(features$count_ion_indices)
        features$m.z <- as.numeric(features$m.z)
        features$Charge <- as.numeric(features$Charge)
        features$mean_Scores <- as.character(features$mean_Scores)
        features$num_matches <- as.character(features$num_matches)
        features$Modifications <- as.character(features$Modifications)
        features$RT_length <- as.numeric(features$RT_length)
        features$Observed_RT <- as.character(features$Observed_RT)
        features$Observed_mz <- as.character(features$Observed_mz)
        features$Observed_score <- as.character(features$Observed_score)
      }else if(multiplicity > 1)
      {
        features <- base::as.data.frame(matrix(ncol=24,nrow=length(unique(base::paste(temp$Sequence,temp$Charge,temp$Modifications)))))
        colnames(features) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score","m.z._shift_light","m.z._shift_medium","m.z._shift_heavy")

        features <- data.table::as.data.table(features)
        features$Feature_name <- as.character(features$Feature_name)
        features$Mass <- as.numeric(features$Mass)
        features$RT <- as.numeric(features$RT)
        features$Sequence <- as.character(features$Sequence)
        features$Protein <- as.character(features$Protein)
        features$MSMS.Scan.Numbers <- as.character(features$MSMS.Scan.Numbers)
        features$m.z_range_min <- as.numeric(features$m.z_range_min)
        features$m.z_range_max <- as.numeric(features$m.z_range_max)
        features$RT_range_min <- as.numeric(features$RT_range_min)
        features$RT_range_max <- as.numeric(features$RT_range_max)
        features$ion_indices <- as.character(features$ion_indices)
        features$count_ion_indices <- as.numeric(features$count_ion_indices)
        features$m.z <- as.numeric(features$m.z)
        features$Charge <- as.numeric(features$Charge)
        features$mean_Scores <- as.character(features$mean_Scores)
        features$num_matches <- as.character(features$num_matches)
        features$Modifications <- as.character(features$Modifications)
        features$RT_length <- as.numeric(features$RT_length)
        features$Observed_RT <- as.character(features$Observed_RT)
        features$Observed_mz <- as.character(features$Observed_mz)
        features$Observed_score <- as.character(features$Observed_score)
        features$m.z._shift_light <- as.numeric(features$m.z._shift_light)
        features$m.z._shift_medium <- as.numeric(features$m.z._shift_medium)
        features$m.z._shift_heavy <- as.numeric(features$m.z._shift_heavy)
      }

    }
    if(MassSpec_mode == "TIMSToF")
    {
      features <- base::as.data.frame(matrix(ncol=26,nrow=length(unique(base::paste(temp$Sequence,temp$Charge,temp$Modifications)))))
      colnames(features) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score","Inv_K0","Inv_K0_length","Inv_K0_range_min","Inv_K0_range_max","Observed_IM")

      features <- data.table::as.data.table(features)
      features$Feature_name <- as.character(features$Feature_name)
      features$Mass <- as.numeric(features$Mass)
      features$RT <- as.numeric(features$RT)
      features$Sequence <- as.character(features$Sequence)
      features$Protein <- as.character(features$Protein)
      features$MSMS.Scan.Numbers <- as.character(features$MSMS.Scan.Numbers)
      features$m.z_range_min <- as.numeric(features$m.z_range_min)
      features$m.z_range_max <- as.numeric(features$m.z_range_max)
      features$RT_range_min <- as.numeric(features$RT_range_min)
      features$RT_range_max <- as.numeric(features$RT_range_max)
      features$ion_indices <- as.character(features$ion_indices)
      features$count_ion_indices <- as.numeric(features$count_ion_indices)
      features$m.z <- as.numeric(features$m.z)
      features$Charge <- as.numeric(features$Charge)
      features$mean_Scores <- as.character(features$mean_Scores)
      features$num_matches <- as.character(features$num_matches)
      features$Modifications <- as.character(features$Modifications)
      features$RT_length <- as.numeric(features$RT_length)
      features$Observed_RT <- as.character(features$Observed_RT)
      features$Observed_mz <- as.character(features$Observed_mz)
      features$Observed_score <- as.character(features$Observed_score)
      features$Inv_K0 <- as.numeric(features$Inv_K0)
      features$Inv_K0_length <- as.numeric(features$Inv_K0_length)
      features$Inv_K0_range_min <- as.numeric(features$Inv_K0_range_min)
      features$Inv_K0_range_max <- as.numeric(features$Inv_K0_range_max)
      features$Observed_IM <- as.character(features$Observed_IM)
    }

    ###get ordering of allpeptides based on sequence
    order_by_Sequence <- order(allpeptides$Sequence,decreasing = T)
    unique_peptides <- sort(unique_peptides,decreasing = T)
    number_of_rows_per_peptide_sequence <- plyr::count(allpeptides$Sequence[which(allpeptides$Sequence != " " & allpeptides$Sequence != "")])
    number_of_rows_per_peptide_sequence <- number_of_rows_per_peptide_sequence[order(number_of_rows_per_peptide_sequence$x,decreasing = T),]

    max <- length(unique_peptides)
    pb <- tcltk::tkProgressBar(title = "Matching known features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    last_index <- 1 ###required for subsetting allpeptides per peptide sequence
    count_features <- 0

    if(MassSpec_mode == "Orbitrap")
    {
      for(i in 1:length(unique_peptides))#
      {
        sub_peptide <- allpeptides[order_by_Sequence[last_index:(last_index-1+number_of_rows_per_peptide_sequence[i,2])],]
        last_index <- last_index + number_of_rows_per_peptide_sequence[i,2]
        sub_peptide <- sub_peptide[which(sub_peptide$Sequence == unique_peptides[i]),]

        for(c in unique(sub_peptide$Charge))
        {
          for(m in unique(sub_peptide$Modifications))
          {
            sub <- sub_peptide[which(sub_peptide$Charge == c & sub_peptide$Modifications == m),]
            if(nrow(sub) > 0)
            {
              sub <- sub[order(sub$Intensity,decreasing = T),]

              ###only use features close to median feature over all samples for this peptide
              median_m.z <- stats::median(sub$m.z,na.rm=T)##
              median_RT <- stats::median(sub$Retention.time,na.rm=T)
              selection <- which(abs(sub$Retention.time-median_RT) < borders_RT[2] & abs(sub$m.z-median_m.z) < borders_m.z[2])
              sub <- sub[selection,]

              if(nrow(sub) > 0)
              {
                ###check if peptide was sequenced several times. If this is the case, use feature closest to median
                if(any(duplicated(sub$Raw.file)))
                {
                  sub <- sub[!duplicated(sub$Raw.file),]
                }

                median_m.z <- stats::median(sub$m.z,na.rm=T)##
                median_RT <- stats::median(sub$Retention.time,na.rm=T)

                RT_max_min <- max(sub$Retention.time,na.rm=T) - min(sub$Retention.time,na.rm=T)
                too_large_RT_variation <- F

                ###get subset of allpeptides based on charge and mz window indexing
                start_indices <- which(allpeptides_frag_indices_per_mz_window[[c]][,1] > median_m.z + borders_m.z[1])[1] - 1
                end_indices <- which(allpeptides_frag_indices_per_mz_window[[c]][,2] > median_m.z + borders_m.z[2])[1] + 1

                if(is.na(end_indices) | end_indices > nrow(allpeptides_frag_indices_per_mz_window[[c]]))end_indices <- nrow(allpeptides_frag_indices_per_mz_window[[c]])

                if(is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & !is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices+1,3]))
                {
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices+1,3]

                }else if(!is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & is.na(allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]))
                {
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,4]
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]
                }else if(is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & !is.na(allpeptides_frag_indices_per_mz_window[[c]][end_indices,3]))
                {
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,3]
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                }else
                {
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                }

                temp_sub <- allpeptides_frag[[c]][start_indices:end_indices,]
                sub <- temp_sub[which(temp_sub$m.z >= median_m.z + borders_m.z[1] & temp_sub$m.z <= median_m.z + borders_m.z[2] & temp_sub$Retention.time >= median_RT + borders_RT[1] & temp_sub$Retention.time <= median_RT + borders_RT[2] | temp_sub$Sequence == unique_peptides[i] & temp_sub$Charge == c & temp_sub$Modifications == m),]

                sub <- sub[order(sub$Intensity,sub$Sequence,decreasing = T),]

                ###check if selected MaxQ features in expected RT and mz window contains outliers (too high deviation in RT or mz compared to all other features)
                box_stats_rt <- grDevices::boxplot.stats(sub$Retention.time)$stats
                box_stats_mz <- grDevices::boxplot.stats(sub$m.z)$stats
                outlier <- which(sub$Retention.time < box_stats_rt[2] | sub$Retention.time > box_stats_rt[4] | sub$m.z < box_stats_mz[2] | sub$m.z > box_stats_mz[4])
                #features with sequence identification will be never regarded as outlier
                if(any(sub$Sequence[outlier] == unique_peptides[i]))
                {
                  outlier <- outlier[-which(sub$Sequence[outlier] == unique_peptides[i])]
                }
                if(length(outlier)>0) ###any outlier detected? exclude them
                {
                  sub <- sub[-outlier,]
                }

                ###check if the fraction of features with not currently searched sequence is > 25 %
                ###if yes, keep all sequences and indicate this with multiple sequences for this feature otherwise remove other sequences
                count_sequences <- plyr::count(sub$Sequence[which(sub$Sequence != " " & sub$Sequence != "")])
                if(count_sequences$freq[which(count_sequences$x == unique_peptides[i])] >= 0.8*sum(count_sequences$freq))
                {
                  sub <- sub[which(sub$Sequence == unique_peptides[i] | sub$Sequence == " " | sub$Sequence == ""),]
                }

                if(any(duplicated(sub$Raw.file)))
                {
                  sub <- sub[!duplicated(sub$Raw.file),]
                }

                if(nrow(sub) == 0) ###variation in RT for same peptide sequence between samples is too high
                {
                  borders_RT_save <- borders_RT
                  too_large_RT_variation <- T
                  borders_RT <- c(-RT_max_min/2,RT_max_min/2)
                  sub <- temp_sub[which(temp_sub$m.z >= median_m.z + borders_m.z[1] & temp_sub$m.z <= median_m.z + borders_m.z[2] & temp_sub$Retention.time >= median_RT + borders_RT[1] & temp_sub$Retention.time <= median_RT + borders_RT[2]),]
                  if(nrow(sub) == 0)###still no match
                  {
                    too_large_RT_variation <- F
                    borders_RT <- borders_RT_save
                  }
                }

                if(nrow(sub) > 0)
                {
                  count_features <- count_features + 1

                  ###new feature - give a new name
                  name <- base::paste("Feature",count_features,sep="_")

                  sequences_relevant <- which(nchar(as.character(sub$Sequence)) > 1)
                  if(length(sequences_relevant) > 0)
                  {
                    sequence <- base::paste(unique(c(sub_peptide$Sequence[1],sub[sequences_relevant,"Sequence"])),collapse=";")
                  }else
                  {
                    sequence <- sub_peptide$Sequence[1]
                  }

                  proteins_relevant <- which(nchar(as.character(sub$Proteins)) > 1)
                  if(length(proteins_relevant) > 0)
                  {
                    protein <- base::paste(unique(c(sub_peptide$Proteins[1],sub[proteins_relevant,"Proteins"])),collapse=";")
                  }else
                  {
                    protein <- sub_peptide$Proteins[1]
                  }

                  msmsscansrelevant_relevant <- which(nchar(as.character(sub$MSMS.Scan.Numbers)) > 1)
                  if(length(msmsscansrelevant_relevant) > 0)
                  {
                    msmsscan <- base::paste(unique(base::paste(sub[msmsscansrelevant_relevant,"Raw.file"],"_msms_",sub[msmsscansrelevant_relevant,"MSMS.Scan.Numbers"],sep="")),collapse=";")
                  }else
                  {
                    msmsscan <- ""
                  }

                  median_mass <- stats::median(sub$Mass,na.rm=T)
                  Charge <- c

                  scores <- NULL
                  num_matches <- NULL
                  for(s in unique(unlist(stringr::str_split(sequence,";"))))
                  {
                    scores <- append(scores,mean(sub$Score[which(sub$Sequence == s)],na.rm=T))
                    num_matches <- append(num_matches,nrow(sub[which(sub$Sequence == s),]))
                  }

                  Modifications <- ifelse(m != "Unmodified",m,"")

                  RT_length <- max(sub$Retention.Length,na.rm = T)

                  temp_RT <- stats::aggregate(sub$Retention.time,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                  Observed_RT <- base::paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

                  if(multiplicity == 1) #use m.z corrected for isotope pattern
                  {
                    temp_mz <- stats::aggregate(sub$isotope_corrected_mz_at_max_int,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                    Observed_mz <- base::paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list
                  }else if(multiplicity > 1)
                  {
                    temp_mz <- stats::aggregate(sub$m.z_SILAC,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                    Observed_mz <- base::paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

                    #calculate expected heavy isotope induced shifts
                    lys_count <- stats::median(sub$Lys_count,na.rm=T)
                    arg_count <- stats::median(sub$Arg_count,na.rm=T)
                    isotope_shifts <- c(((SILAC_settings$light$mass_shift_Lys*lys_count) + (SILAC_settings$light$mass_shift_Arg*arg_count))/c,
                                        ((SILAC_settings$medium$mass_shift_Lys*lys_count) + (SILAC_settings$medium$mass_shift_Arg*arg_count))/c,
                                        ((SILAC_settings$heavy$mass_shift_Lys*lys_count) + (SILAC_settings$heavy$mass_shift_Arg*arg_count))/c)
                  }

                  temp_score <- stats::aggregate(sub$Score,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                  Observed_score <- base::paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed scores according to sample_list


                  data.table::set(x = features,i = as.integer(count_features),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
                                                                                                                                    sequence,
                                                                                                                                    protein,
                                                                                                                                    msmsscan,
                                                                                                                                    base::paste(rownames(sub),collapse = ","),
                                                                                                                                    base::paste(scores,collapse=";"),
                                                                                                                                    base::paste(num_matches,collapse=";"),
                                                                                                                                    Modifications,
                                                                                                                                    Observed_RT,
                                                                                                                                    Observed_mz,
                                                                                                                                    Observed_score
                  )))

                  if(multiplicity == 1)
                  {
                    data.table::set(x = features,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
                      as.numeric(median_m.z),
                      as.numeric(median_RT),
                      as.numeric(median_m.z+borders_m.z[1]),
                      as.numeric(median_m.z+borders_m.z[2]),
                      as.numeric(median_RT+borders_RT[1]),
                      as.numeric(median_RT+borders_RT[2]),
                      as.numeric(nrow(sub)),
                      as.numeric(median_mass),
                      as.numeric(Charge),
                      RT_length
                    )))
                  }else if(multiplicity > 1)
                  {
                    data.table::set(x = features,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18,22:24)),value = as.list(c(
                      as.numeric(median_m.z),
                      as.numeric(median_RT),
                      as.numeric(median_m.z+borders_m.z[1]),
                      as.numeric(median_m.z+borders_m.z[2]),
                      as.numeric(median_RT+borders_RT[1]),
                      as.numeric(median_RT+borders_RT[2]),
                      as.numeric(nrow(sub)),
                      as.numeric(median_mass),
                      as.numeric(Charge),
                      RT_length,
                      isotope_shifts[1],
                      isotope_shifts[2],
                      isotope_shifts[3]
                    )))
                  }


                  ###now mark matched features
                  data.table::set(temp_data,i = as.integer(rownames(sub)),j=1L,value=1)
                  #temp_data[as.numeric(rownames(sub)),1] <- 1

                  ###if RT_window had to be increased, reset borders_RT
                  if(too_large_RT_variation == T)
                  {
                    too_large_RT_variation <- F
                    borders_RT <- borders_RT_save
                  }
                }
              }
            }
          }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 50)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }
    if(MassSpec_mode == "TIMSToF")
    {
      for(i in 1:length(unique_peptides))
      {
        sub_peptide <- allpeptides[order_by_Sequence[last_index:(last_index-1+number_of_rows_per_peptide_sequence[i,2])],]
        last_index <- last_index + number_of_rows_per_peptide_sequence[i,2]
        sub_peptide <- sub_peptide[which(sub_peptide$Sequence == unique_peptides[i]),]

        for(c in unique(sub_peptide$Charge))
        {
          for(m in unique(sub_peptide$Modifications))
          {
            sub <- sub_peptide[which(sub_peptide$Charge == c & sub_peptide$Modifications == m),]
            if(nrow(sub) > 0)
            {
              sub <- sub[order(sub$Intensity,decreasing = T),]

              ###only use features close to median feature over all samples for this peptide
              median_m.z <- stats::median(sub$m.z,na.rm=T)##
              median_RT <- stats::median(sub$Retention.time,na.rm=T)
              median_IM <- stats::median(sub$inverse_K0,na.rm=T)
              selection <- which(abs(sub$inverse_K0-median_IM) < borders_IM[2] & abs(sub$Retention.time-median_RT) < borders_RT[2] & abs(sub$m.z-median_m.z) < borders_m.z[2])
              sub <- sub[selection,]

              if(nrow(sub) > 0)
              {
                ###check if peptide was sequenced several times. If this is the case, use feature closest to median
                if(any(duplicated(sub$Raw.file)))
                {
                  sub <- sub[!duplicated(sub$Raw.file),]
                }

                median_m.z <- stats::median(sub$m.z,na.rm=T)##
                median_RT <- stats::median(sub$Retention.time,na.rm=T)
                median_IM <- stats::median(sub$inverse_K0,na.rm=T)

                RT_max_min <- max(sub$Retention.time,na.rm=T) - min(sub$Retention.time,na.rm=T)
                too_large_RT_variation <- F

                ###get subset of allpeptides based on charge and mz window indexing
                start_indices <- which(allpeptides_frag_indices_per_mz_window[[c]][,1] > median_m.z + borders_m.z[1])[1] - 1
                end_indices <- which(allpeptides_frag_indices_per_mz_window[[c]][,2] > median_m.z + borders_m.z[2])[1] + 1

                if(is.na(end_indices) | end_indices > nrow(allpeptides_frag_indices_per_mz_window[[c]]))end_indices <- nrow(allpeptides_frag_indices_per_mz_window[[c]])

                if(is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & !is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices+1,3]))
                {
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices+1,3]

                }else if(!is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & is.na(allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]))
                {
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,4]
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]
                }else if(is.na(allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]) & !is.na(allpeptides_frag_indices_per_mz_window[[c]][end_indices,3]))
                {
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,3]
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                }else
                {
                  start_indices <- allpeptides_frag_indices_per_mz_window[[c]][start_indices,3]
                  end_indices <- allpeptides_frag_indices_per_mz_window[[c]][end_indices,4]
                }

                temp_sub <- allpeptides_frag[[c]][start_indices:end_indices,]
                sub <- temp_sub[which(temp_sub$m.z >= median_m.z + borders_m.z[1] & temp_sub$m.z <= median_m.z + borders_m.z[2] &
                                        temp_sub$Retention.time >= median_RT + borders_RT[1] & temp_sub$Retention.time <= median_RT + borders_RT[2] &
                                        temp_sub$inverse_K0 >= median_IM + borders_IM[1] & temp_sub$inverse_K0 <= median_IM + borders_IM[2] |
                                        temp_sub$Sequence == unique_peptides[i] & temp_sub$Charge == c & temp_sub$Modifications == m),]

                sub <- sub[order(sub$Intensity,sub$Sequence,decreasing = T),]

                ###check if selected MaxQ features in expected RT, IM and mz window contains outliers (too high deviation in RT, IM or mz compared to all other features)
                box_stats_rt <- grDevices::boxplot.stats(sub$Retention.time)$stats
                box_stats_mz <- grDevices::boxplot.stats(sub$m.z)$stats
                box_stats_im <- grDevices::boxplot.stats(sub$inverse_K0)$stats
                outlier <- which(sub$Retention.time < box_stats_rt[2] | sub$Retention.time > box_stats_rt[4] |
                                   sub$m.z < box_stats_mz[2] | sub$m.z > box_stats_mz[4] |
                                   sub$inverse_K0 < box_stats_im[2] | sub$inverse_K0 > box_stats_im[4])
                #features with sequence identification will be never regarded as outlier
                if(any(sub$Sequence[outlier] == unique_peptides[i]))
                {
                  outlier <- outlier[-which(sub$Sequence[outlier] == unique_peptides[i])]
                }
                if(length(outlier)>0) ###any outlier detected? exclude them
                {
                  sub <- sub[-outlier,]
                }

                ###check if the fraction of features with not currently searched sequence is > 25 %
                ###if yes, keep all sequences and indicate this with multiple sequences for this feature otherwise remove other sequences
                count_sequences <- plyr::count(sub$Sequence[which(sub$Sequence != " " & sub$Sequence != "")])
                if(count_sequences$freq[which(count_sequences$x == unique_peptides[i])] >= 0.8*sum(count_sequences$freq))
                {
                  sub <- sub[which(sub$Sequence == unique_peptides[i] | sub$Sequence == " " | sub$Sequence == ""),]
                }

                if(any(duplicated(sub$Raw.file)))
                {
                  sub <- sub[!duplicated(sub$Raw.file),]
                }

                if(nrow(sub) == 0) ###variation in RT for same peptide sequence between samples is too high
                {
                  borders_RT_save <- borders_RT
                  too_large_RT_variation <- T
                  borders_RT <- c(-RT_max_min/2,RT_max_min/2)
                  sub <- temp_sub[which(temp_sub$m.z >= median_m.z + borders_m.z[1] & temp_sub$m.z <= median_m.z + borders_m.z[2] &
                                          temp_sub$Retention.time >= median_RT + borders_RT[1] & temp_sub$Retention.time <= median_RT + borders_RT[2] &
                                          temp_sub$inverse_K0 >= median_IM + borders_IM[1] & temp_sub$inverse_K0 <= median_IM + borders_IM[2]),]
                  if(nrow(sub) == 0)###still no match
                  {
                    too_large_RT_variation <- F
                    borders_RT <- borders_RT_save
                  }
                }

                if(nrow(sub) > 0)
                {
                  count_features <- count_features + 1

                  ###new feature - give a new name
                  name <- base::paste("Feature",count_features,sep="_")

                  sequences_relevant <- which(nchar(as.character(sub$Sequence)) > 1)
                  if(length(sequences_relevant) > 0)
                  {
                    sequence <- base::paste(unique(c(sub_peptide$Sequence[1],sub[sequences_relevant,"Sequence"])),collapse=";")
                  }else
                  {
                    sequence <- sub_peptide$Sequence[1]
                  }

                  proteins_relevant <- which(nchar(as.character(sub$Proteins)) > 1)
                  if(length(proteins_relevant) > 0)
                  {
                    protein <- base::paste(unique(c(sub_peptide$Proteins[1],sub[proteins_relevant,"Proteins"])),collapse=";")
                  }else
                  {
                    protein <- sub_peptide$Proteins[1]
                  }

                  msmsscansrelevant_relevant <- which(nchar(as.character(sub$MSMS.Scan.Numbers)) > 1)
                  if(length(msmsscansrelevant_relevant) > 0)
                  {
                    msmsscan <- base::paste(unique(base::paste(sub[msmsscansrelevant_relevant,"Raw.file"],"_msms_",sub[msmsscansrelevant_relevant,"MSMS.Scan.Numbers"],sep="")),collapse=";")
                  }else
                  {
                    msmsscan <- ""
                  }

                  median_mass <- stats::median(sub$Mass,na.rm=T)
                  Charge <- c

                  scores <- NULL
                  num_matches <- NULL
                  for(s in unique(unlist(stringr::str_split(sequence,";"))))
                  {
                    scores <- append(scores,mean(sub$Score[which(sub$Sequence == s)],na.rm=T))
                    num_matches <- append(num_matches,nrow(sub[which(sub$Sequence == s),]))
                  }

                  Modifications <- ifelse(m != "Unmodified",m,"")

                  RT_length <- max(sub$Retention.Length,na.rm = T)
                  temp_RT <- stats::aggregate(sub$Retention.time,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
                  Observed_RT <- base::paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

                  temp_mz <- stats::aggregate(sub$m.z,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
                  Observed_mz <- base::paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

                  temp_score <- stats::aggregate(sub$Score,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
                  Observed_score <- base::paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed scores according to sample_list

                  IM_length <- max(sub$inverse_K0_length,na.rm = T)
                  temp_IM <- stats::aggregate(sub$inverse_K0,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
                  Observed_IM <- base::paste(temp_IM[match(sample_list,temp_IM$Raw.file),2],collapse=";") ###order all observed IM according to sample_list



                  data.table::set(x = features,i = as.integer(count_features),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21,26)),value = as.list(c(name,
                                                                                                                                       sequence,
                                                                                                                                       protein,
                                                                                                                                       msmsscan,
                                                                                                                                       base::paste(rownames(sub),collapse = ","),
                                                                                                                                       base::paste(scores,collapse=";"),
                                                                                                                                       base::paste(num_matches,collapse=";"),
                                                                                                                                       Modifications,
                                                                                                                                       Observed_RT,
                                                                                                                                       Observed_mz,
                                                                                                                                       Observed_score,
                                                                                                                                       Observed_IM
                  )))
                  data.table::set(x = features,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18,22,23,24,25)),value = as.list(c(
                    as.numeric(median_m.z),
                    as.numeric(median_RT),
                    as.numeric(median_m.z+borders_m.z[1]),
                    as.numeric(median_m.z+borders_m.z[2]),
                    as.numeric(median_RT+borders_RT[1]),
                    as.numeric(median_RT+borders_RT[2]),
                    as.numeric(nrow(sub)),
                    as.numeric(median_mass),
                    as.numeric(Charge),
                    RT_length,
                    median_IM,
                    IM_length,
                    as.numeric(median_IM+borders_IM[1]),
                    as.numeric(median_IM+borders_IM[2])
                  )))

                  ###now mark matched features
                  data.table::set(temp_data,i = as.integer(rownames(sub)),j=1L,value=1)
                  #temp_data[as.numeric(rownames(sub)),1] <- 1

                  ###if RT_window had to be increased, reset borders_RT
                  if(too_large_RT_variation == T)
                  {
                    too_large_RT_variation <- F
                    borders_RT <- borders_RT_save
                  }
                }
              }
            }
          }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 50)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }

    close(pb)

    features <- features[1:count_features,]
    #write.table(features,base::paste("Temporary_files/Features_aligned_peptides",output_file_names_add,".txt",sep=""),row.names = F)

    borders_RT <- borders_RT_use_save

    if(MassSpec_mode == "Orbitrap")save(temp_data,features,allpeptides,borders_m.z,borders_RT,windows,file=base::paste("Temporary_files/features_aligned_step_1",output_file_names_add,".RData",sep=""))
    if(MassSpec_mode == "TIMSToF")save(temp_data,features,allpeptides,borders_m.z,borders_RT,borders_IM,windows,file=base::paste("Temporary_files/features_aligned_step_1",output_file_names_add,".RData",sep=""))
    #####Perform alignment of undefined features

    # if(align_unknown == T)
    # {
    #   ###2.) step - align remaining features
    #
    #   #load(base::paste("features_aligned_step_1",output_file_names_add,".RData",sep=""))
    #
    #   count_features <- nrow(features)
    #
    #   ###start by trying to find unique features over all samples defined by the pair Mass and RT
    #   border_factor <- 4
    #
    #   ####use this matrix for faster subsetting
    #   temp_data_m <- as.matrix(subset(allpeptides,select=c("m.z","Retention.time","Charge")))
    #   temp_data_m <- cbind(temp_data_m,temp_data)
    #   colnames(temp_data_m)[4] <- "Matched"
    #   temp_data_m[which(is.na(temp_data_m[,4])),4] <- 0
    #   #temp_data_m <- cbind(temp_data_m,as.matrix(subset(allpeptides,select=c("Retention.Length"))))
    #
    #   features_unknown <- base::as.data.frame(matrix(ncol=21,nrow=(nrow(temp_data)/min_num_ions_collapse)))
    #   colnames(features_unknown) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")
    #
    #
    #   features_unknown <- data.table::as.data.table(features_unknown)
    #   features_unknown$Feature_name <- as.character(features_unknown$Feature_name)
    #   features_unknown$Mass <- as.numeric(features_unknown$Mass)
    #   features_unknown$RT <- as.numeric(features_unknown$RT)
    #   features_unknown$Sequence <- as.character(features_unknown$Sequence)
    #   features_unknown$Protein <- as.character(features_unknown$Protein)
    #   features_unknown$MSMS.Scan.Numbers <- as.character(features_unknown$MSMS.Scan.Numbers)
    #   features_unknown$m.z_range_min <- as.numeric(features_unknown$m.z_range_min)
    #   features_unknown$m.z_range_max <- as.numeric(features_unknown$m.z_range_max)
    #   features_unknown$RT_range_min <- as.numeric(features_unknown$RT_range_min)
    #   features_unknown$RT_range_max <- as.numeric(features_unknown$RT_range_max)
    #   features_unknown$ion_indices <- as.character(features_unknown$ion_indices)
    #   features_unknown$count_ion_indices <- as.numeric(features_unknown$count_ion_indices)
    #   features_unknown$m.z <- as.numeric(features_unknown$m.z)
    #   features_unknown$Charge <- as.numeric(features_unknown$Charge)
    #   features_unknown$mean_Scores <- as.character(features_unknown$mean_Scores)
    #   features_unknown$num_matches <- as.character(features_unknown$num_matches)
    #   features_unknown$Modifications <- as.character(features_unknown$Modifications)
    #   features_unknown$RT_length <- as.numeric(features_unknown$RT_length)
    #   features_unknown$Observed_RT <- as.character(features_unknown$Observed_RT)
    #   features_unknown$Observed_mz <- as.character(features_unknown$Observed_mz)
    #   features_unknown$Observed_score <- as.character(features_unknown$Observed_score)
    #
    #   max <- nrow(temp_data_m)
    #   pb <- tcltk::tkProgressBar(title = "Aligning unknown features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    #   start_time <- Sys.time()
    #   updatecounter <- 0
    #   time_require <- 0
    #   start <- 1
    #   end <- max-1
    #
    #   count_features_unknown <- 0
    #
    #   for(i in start:end) ###go through all features and try to match with others in the data set
    #   {
    #     if(temp_data_m[i,4] == 0)###check if already matched
    #     {
    #       ind <- i
    #       while(T)
    #       {
    #         ind = ind + 1000
    #         if(ind > nrow(temp_data_m))
    #         {
    #           ind <- nrow(temp_data_m)
    #           break
    #         }
    #         if(temp_data_m[ind,1] > (temp_data_m[i,1]+(border_factor*borders_m.z[2])))
    #         {
    #           break
    #         }
    #       }
    #
    #       sub <- temp_data_m[i:ind,]
    #       sub <- sub[which(sub[,4] == 0),]
    #
    #       ###determine which other features might fit to the current unmatched feature
    #       ###here use a wide mass and RT range to get a first preselection for Mass and RT distribution estimation                                  )
    #
    #       matches_preselection <- which(sub[,1] >= (temp_data_m[i,1]+(border_factor*borders_m.z[1])) & sub[,1] <= (temp_data_m[i,1]+(border_factor*borders_m.z[2])) & sub[,2] >= (temp_data_m[i,2]+(border_factor*borders_RT[1])) & sub[,2] <= (temp_data_m[i,2]+(border_factor*borders_RT[2])) & as.numeric(sub[,3]) == temp_data_m[i,3])
    #
    #       #Sys.time()-start
    #       if(length(matches_preselection) >= min_num_ions_collapse) ###at least n other features
    #       {
    #         ###now get estimate of true Mass and RT by using the median
    #         median_m.z <- stats::median(sub[matches_preselection,1],na.rm=T)
    #         median_RT <- stats::median(sub[matches_preselection,2],na.rm=T)
    #
    #         ###determine which features will be finally matched with the normal m/z and RT range
    #         matches_temp <- which(as.numeric(sub[,1]) >= (median_m.z+borders_m.z[1]) & as.numeric(sub[,1]) <= (median_m.z+borders_m.z[2]) & as.numeric(sub[,2]) >= (median_RT+borders_RT[1]) & as.numeric(sub[,2]) <= (median_RT+borders_RT[2]) & as.numeric(sub[,3]) == temp_data_m[i,3])
    #
    #         matches <- as.numeric(rownames(sub)[matches_temp]) ###convert back to original data
    #
    #         if(length(matches) >= min_num_ions_collapse)
    #         {
    #           count_features <- count_features + 1
    #           count_features_unknown <- count_features_unknown + 1
    #           ###new feature - give a new name
    #           name <- base::paste("Feature",count_features,sep="_")
    #
    #           sequences_relevant <- which(nchar(as.character(allpeptides[matches,"Sequence"])) > 1)
    #           if(length(sequences_relevant) > 0)
    #           {
    #             sequence <- base::paste(unique(allpeptides[matches[sequences_relevant],"Sequence"]),collapse=";")
    #           }else
    #           {
    #             sequence <- ""
    #           }
    #
    #           proteins_relevant <- which(nchar(as.character(allpeptides[matches,"Proteins"])) > 1)
    #           if(length(proteins_relevant) > 0)
    #           {
    #             protein <- base::paste(unique(allpeptides[matches[proteins_relevant],"Proteins"]),collapse=";")
    #           }else
    #           {
    #             protein <- ""
    #           }
    #
    #           msmsscansrelevant_relevant <- which(nchar(as.character(allpeptides[matches,"MSMS.Scan.Numbers"])) > 1)
    #           if(length(msmsscansrelevant_relevant) > 0)
    #           {
    #             msmsscan <- base::paste(unique(base::paste(allpeptides[matches[msmsscansrelevant_relevant],"Raw.file"],"_msms_",allpeptides[matches[msmsscansrelevant_relevant],"MSMS.Scan.Numbers"],sep="")),collapse=";")
    #           }else
    #           {
    #             msmsscan <- ""
    #           }
    #
    #           median_mass <- stats::median(allpeptides$Mass[matches],na.rm=T)
    #           Charge <- temp_data_m[i,3]
    #
    #
    #           if(sequence != "")
    #           {
    #             scores <- NULL
    #             num_matches <- NULL
    #             for(s in unique(unlist(stringr::str_split(sequence,";"))))
    #             {
    #               scores <- append(scores,mean(allpeptides[matches,"Score"][which(allpeptides[matches,"Sequence"] == s)],na.rm=T))
    #               num_matches <- append(num_matches,length(which(allpeptides[matches,"Sequence"] == s)))
    #             }
    #
    #             Modifications_relevant <- which(nchar(as.character(allpeptides[matches,"Modifications"])) > 1)
    #             if(length(Modifications_relevant) > 0)
    #             {
    #               Modifications <- base::paste(unique(allpeptides[matches[Modifications_relevant],"Modifications"]),collapse=";")
    #               Modifications <- base::gsub("Unmodified","",Modifications)
    #
    #             }else
    #             {
    #               Modifications <- ""
    #             }
    #           }else
    #           {
    #             scores <- NA
    #             num_matches <- NA
    #             Modifications <- ""
    #           }
    #
    #           RT_length <- stats::median(allpeptides$Retention.Length[matches],na.rm=T)
    #
    #           temp_RT <- stats::aggregate(allpeptides$Retention.time[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
    #           Observed_RT <- base::paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list
    #
    #           temp_mz <- stats::aggregate(allpeptides$isotope_corrected_mz_at_max_int[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
    #           Observed_mz <- base::paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list
    #
    #           temp_score <- stats::aggregate(allpeptides$Score[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
    #           Observed_score <- base::paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed scores according to sample_list
    #
    #
    #           data.table::set(x = features_unknown,i = as.integer(count_features_unknown),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
    #                                                                                                                                             sequence,
    #                                                                                                                                             protein,
    #                                                                                                                                             msmsscan,
    #                                                                                                                                             base::paste(matches,collapse = ","),
    #                                                                                                                                             base::paste(scores,collapse=";"),
    #                                                                                                                                             base::paste(num_matches,collapse=";"),
    #                                                                                                                                             Modifications,
    #                                                                                                                                             Observed_RT,
    #                                                                                                                                             Observed_mz,
    #                                                                                                                                             Observed_score
    #           )))
    #           data.table::set(x = features_unknown,i = as.integer(count_features_unknown),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
    #             as.numeric(median_m.z),
    #             as.numeric(median_RT),
    #             as.numeric(median_m.z+borders_m.z[1]),
    #             as.numeric(median_m.z+borders_m.z[2]),
    #             as.numeric(median_RT+borders_RT[1]),
    #             as.numeric(median_RT+borders_RT[2]),
    #             as.numeric(length(matches)),
    #             as.numeric(median_mass),
    #             as.numeric(Charge),
    #             RT_length
    #           )))
    #
    #
    #           ###now mark matched features
    #           # starttime <- Sys.time()
    #           # temp_data_m[matches,4] <- 1
    #           # print(Sys.time()-starttime)
    #
    #           data.table::set(temp_data_m,i=as.integer(matches),j=4L,value=1)
    #
    #         }
    #
    #
    #       }
    #     }
    #
    #     updatecounter <- updatecounter + 1
    #     if(updatecounter >= 5000)
    #     {
    #       time_elapsed <- difftime(Sys.time(),start_time,units="secs")
    #       time_require <- (time_elapsed/(i/max))*(1-(i/max))
    #       td <- lubridate::seconds_to_period(time_require)
    #       time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))
    #
    #       updatecounter <- 0
    #       tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    #     }
    #
    #   }
    #   close(pb)
    #
    #   features_unknown <- features_unknown[1:count_features_unknown,]
    #
    #   #write.table(features_unknown,base::paste("Temporary_files/Features_aligned_unknown",output_file_names_add,".txt",sep=""),row.names = F)
    #
    #   ####Combine known and unknown features
    #   features <- rbind(features,features_unknown)
    #   rownames(features) <- c()
    #
    #   ###number of matched to number of unmatched ions
    #   value <- length(which(temp_data_m$Matched == 1))/nrow(temp_data_m) * 100
    #   barplot(value,main="Relative fraction of features which are matched",ylab="Fraction [%]")
    #
    #   QC_data[["Numbers_matched_features"]] <- value
    #
    #   ###Store all results
    #
    #   save(temp_data_m,features,count_features_unknown,file=base::paste("features_aligned_step_2",output_file_names_add,".RData",sep=""))
    # }

    #write.table(features,base::paste("Temporary_files/Features_aligned",output_file_names_add,".txt",sep=""),row.names = F)

    ####Merge features which were currently regarded as separated although same RT and m/z
    #features <- utils::read.table("Features_alligned.txt",header = T)

    #####prepare columns to indicate number of charges and with which other features the respective feature was merged
    features$num_diff_charges <- NA
    features$num_diff_charges <- as.numeric(features$num_diff_charges)
    features$merged_with <- ""
    ncolumns <- ncol(features)

    ###collapse features with same m/z and RT parameters and count number of different charge states for same feature
    borders_mass <- c(-feature_mass_deviation_collapse,feature_mass_deviation_collapse)

    max <- nrow(features)
    pb <- tcltk::tkProgressBar(title = "Merge features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    if(MassSpec_mode == "Orbitrap")
    {
      for(i in 1:nrow(features))
      {
        if(is.na(features[i,"num_diff_charges"]))
        {
          selection <- which(features$Mass >= features$Mass[i]+borders_mass[1] & features$Mass <= features$Mass[i]+borders_mass[2] &
                               features$RT >= features$RT[i]+borders_RT[1] & features$RT <= features$RT[i]+borders_RT[2])
          if(length(selection) > 1)
          {
            ###merge same charges
            for(c in unique(features[selection,]$Charge))
            {
              if(length(which(features[selection,]$Charge == c)) > 1) ###more than one feature with same charge -> merge
              {
                m.z_range_min <- min(features[selection[which(features[selection,]$Charge == c)],]$m.z_range_min,na.rm=T)
                m.z_range_max <- max(features[selection[which(features[selection,]$Charge == c)],]$m.z_range_max,na.rm=T)
                RT_range_min <- min(features[selection[which(features[selection,]$Charge == c)],]$RT_range_min,na.rm=T)
                RT_range_max <- max(features[selection[which(features[selection,]$Charge == c)],]$RT_range_max,na.rm=T)
                ion_indices <- base::paste(unlist(features[selection[which(features[selection,]$Charge == c)],]$ion_indices),collapse=",")
                count_ion_indices <- sum(features[selection[which(features[selection,]$Charge == c)],]$count_ion_indices,na.rm = T)
                mass <- mean(unlist(features[selection[which(features[selection,]$Charge == c)],]$Mass),na.rm=T)
                m.z <- mean(c(m.z_range_min,m.z_range_max))
                RT <- mean(c(RT_range_min,RT_range_max))

                Sequence <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Sequence,";")))

                scores <- NULL
                num_matches <- NULL
                Modifications <- NULL
                for(s in unique(Sequence))
                {
                  scores_temp <- as.numeric(unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$mean_Scores,";")))
                  scores <- append(scores,mean(scores_temp,na.rm=T))

                  num_matches_temp <- as.numeric(unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$num_matches,";")))
                  num_matches <- append(num_matches,sum(num_matches_temp,na.rm=T))

                  Modifications_temp <- unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$Modifications,";"))
                  Modifications <- append(Modifications,base::paste(Modifications_temp,collapse=","))
                }
                scores <- base::paste(scores,collapse=",")
                num_matches <- base::paste(num_matches,collapse=",")
                Modifications <- base::paste(Modifications,collapse=",")

                if(any(Sequence == ""))Sequence <- Sequence[-which(Sequence == "")]
                Sequence <- base::paste(Sequence,collapse=";")

                Protein <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Protein,";")))
                if(any(Protein == ""))Protein <- Protein[-which(Protein == "")]
                Protein <- base::paste(Protein,collapse=";")

                MSMS.Scan.Numbers <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$MSMS.Scan.Numbers,";")))
                if(any(MSMS.Scan.Numbers == ""))MSMS.Scan.Numbers <- MSMS.Scan.Numbers[-which(MSMS.Scan.Numbers == "")]
                MSMS.Scan.Numbers <- base::paste(MSMS.Scan.Numbers,collapse=";")

                Observed_RTs <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_RT,pattern = ";",simplify = T)
                Observed_RTs <- apply(Observed_RTs, 2, as.numeric)
                Observed_RTs <- colMeans(Observed_RTs,na.rm=T)
                Observed_RTs[is.na(Observed_RTs)] <- NA
                Observed_RTs <- base::paste(Observed_RTs,collapse=";")

                Observed_mz <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_mz,pattern = ";",simplify = T)
                Observed_mz <- apply(Observed_mz, 2, as.numeric)
                Observed_mz <- colMeans(Observed_mz,na.rm=T)
                Observed_mz[is.na(Observed_mz)] <- NA
                Observed_mz <- base::paste(Observed_mz,collapse=";")

                Observed_score <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_score,pattern = ";",simplify = T)
                Observed_score <- apply(Observed_score, 2, as.numeric)
                Observed_score <- colMeans(Observed_score,na.rm=T)
                Observed_score[is.na(Observed_score)] <- NA
                Observed_score <- base::paste(Observed_score,collapse=";")

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(2,3,7,8,9,10,12,13,14)),value=as.list(c(
                  m.z,
                  RT,
                  m.z_range_min,
                  m.z_range_max,
                  RT_range_min,
                  RT_range_max,
                  count_ion_indices,
                  mass,
                  c)))

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(4,5,6,11,15,16,17,19,20,21)),value=as.list(c(
                  Sequence,
                  Protein,
                  MSMS.Scan.Numbers,
                  ion_indices,
                  scores,
                  num_matches,
                  Modifications,
                  Observed_RTs,
                  Observed_mz,
                  Observed_score)))

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[-1]],j=as.integer(ncolumns),value=features$Feature_name[selection[which(features[selection,]$Charge == c)[1]]])
              }
            }
            ###add number of different charges count
            data.table::set(x=features,i=selection,j=as.integer(ncolumns-1),value=length(unique(features[selection,]$Charge)))
          }
        }
        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }
    if(MassSpec_mode == "TIMSToF")
    {
      for(i in 1:nrow(features))
      {
        if(is.na(features[i,"num_diff_charges"]))
        {
          selection <- which(features$Mass >= features$Mass[i]+borders_mass[1] & features$Mass <= features$Mass[i]+borders_mass[2] &
                               features$RT >= features$RT[i]+borders_RT[1] & features$RT <= features$RT[i]+borders_RT[2] &
                               features$Inv_K0 >= features$Inv_K0[i]+borders_IM[1] & features$Inv_K0 <= features$Inv_K0[i]+borders_IM[2])
          if(length(selection) > 1)
          {
            ###merge same charges
            for(c in unique(features[selection,]$Charge))
            {
              if(length(which(features[selection,]$Charge == c)) > 1) ###more than one feature with same charge -> merge
              {
                m.z_range_min <- min(features[selection[which(features[selection,]$Charge == c)],]$m.z_range_min,na.rm=T)
                m.z_range_max <- max(features[selection[which(features[selection,]$Charge == c)],]$m.z_range_max,na.rm=T)
                RT_range_min <- min(features[selection[which(features[selection,]$Charge == c)],]$RT_range_min,na.rm=T)
                RT_range_max <- max(features[selection[which(features[selection,]$Charge == c)],]$RT_range_max,na.rm=T)
                IM_range_min <- min(features[selection[which(features[selection,]$Charge == c)],]$Inv_K0_range_min,na.rm=T)
                IM_range_max <- max(features[selection[which(features[selection,]$Charge == c)],]$Inv_K0_range_max,na.rm=T)
                ion_indices <- base::paste(unlist(features[selection[which(features[selection,]$Charge == c)],]$ion_indices),collapse=",")
                count_ion_indices <- sum(features[selection[which(features[selection,]$Charge == c)],]$count_ion_indices,na.rm = T)
                mass <- mean(unlist(features[selection[which(features[selection,]$Charge == c)],]$Mass),na.rm=T)
                m.z <- mean(c(m.z_range_min,m.z_range_max))
                RT <- mean(c(RT_range_min,RT_range_max))
                IM <- mean(c(IM_range_min,IM_range_max))

                Sequence <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Sequence,";")))

                scores <- NULL
                num_matches <- NULL
                Modifications <- NULL
                for(s in unique(Sequence))
                {
                  scores_temp <- as.numeric(unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$mean_Scores,";")))
                  scores <- append(scores,mean(scores_temp,na.rm=T))

                  num_matches_temp <- as.numeric(unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$num_matches,";")))
                  num_matches <- append(num_matches,sum(num_matches_temp,na.rm=T))

                  Modifications_temp <- unlist(stringr::str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$Modifications,";"))
                  Modifications <- append(Modifications,base::paste(Modifications_temp,collapse=","))
                }
                scores <- base::paste(scores,collapse=",")
                num_matches <- base::paste(num_matches,collapse=",")
                Modifications <- base::paste(Modifications,collapse=",")

                if(any(Sequence == ""))Sequence <- Sequence[-which(Sequence == "")]
                Sequence <- base::paste(Sequence,collapse=";")

                Protein <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Protein,";")))
                if(any(Protein == ""))Protein <- Protein[-which(Protein == "")]
                Protein <- base::paste(Protein,collapse=";")

                MSMS.Scan.Numbers <- unique(unlist(stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$MSMS.Scan.Numbers,";")))
                if(any(MSMS.Scan.Numbers == ""))MSMS.Scan.Numbers <- MSMS.Scan.Numbers[-which(MSMS.Scan.Numbers == "")]
                MSMS.Scan.Numbers <- base::paste(MSMS.Scan.Numbers,collapse=";")

                Observed_RTs <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_RT,pattern = ";",simplify = T)
                Observed_RTs <- apply(Observed_RTs, 2, as.numeric)
                Observed_RTs <- colMeans(Observed_RTs,na.rm=T)
                Observed_RTs[is.na(Observed_RTs)] <- NA
                Observed_RTs <- base::paste(Observed_RTs,collapse=";")

                Observed_mz <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_mz,pattern = ";",simplify = T)
                Observed_mz <- apply(Observed_mz, 2, as.numeric)
                Observed_mz <- colMeans(Observed_mz,na.rm=T)
                Observed_mz[is.na(Observed_mz)] <- NA
                Observed_mz <- base::paste(Observed_mz,collapse=";")

                Observed_score <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_score,pattern = ";",simplify = T)
                Observed_score <- apply(Observed_score, 2, as.numeric)
                Observed_score <- colMeans(Observed_score,na.rm=T)
                Observed_score[is.na(Observed_score)] <- NA
                Observed_score <- base::paste(Observed_score,collapse=";")

                Observed_IMs <- stringr::str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_IM,pattern = ";",simplify = T)
                Observed_IMs <- apply(Observed_IMs, 2, as.numeric)
                Observed_IMs <- colMeans(Observed_IMs,na.rm=T)
                Observed_IMs[is.na(Observed_IMs)] <- NA
                Observed_IMs <- base::paste(Observed_IMs,collapse=";")

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(2,3,7,8,9,10,12,13,14,22,24,25)),value=as.list(c(
                  m.z,
                  RT,
                  m.z_range_min,
                  m.z_range_max,
                  RT_range_min,
                  RT_range_max,
                  count_ion_indices,
                  mass,
                  c,
                  IM,
                  IM_range_min,
                  IM_range_max)))

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(4,5,6,11,15,16,17,19,20,21,26)),value=as.list(c(
                  Sequence,
                  Protein,
                  MSMS.Scan.Numbers,
                  ion_indices,
                  scores,
                  num_matches,
                  Modifications,
                  Observed_RTs,
                  Observed_mz,
                  Observed_score,
                  Observed_IMs)))

                data.table::set(x=features,i=selection[which(features[selection,]$Charge == c)[-1]],j=as.integer(ncolumns),value=features$Feature_name[selection[which(features[selection,]$Charge == c)[1]]])
              }
            }
            ###add number of different charges count
            data.table::set(x=features,i=selection,j=as.integer(ncolumns-1),value=length(unique(features[selection,]$Charge)))
          }
        }
        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
    }

    close(pb)

    features$num_diff_charges[is.na(features$num_diff_charges)] <- 1

    features <- features[which(features$merged_with == ""),]
    features <- base::as.data.frame(features)[,-ncolumns]


    ####now add mz and RT calibration information per feature

    RT_calibration_vals <- base::as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features),NA))
    colnames(RT_calibration_vals) <- base::paste("RT_calibration",sample_list,sep=" ")
    for(c in 1:ncol(RT_calibration_vals))
    {
      RT_calibration_vals[,c] <- as.numeric(RT_calibration_vals[,c])
    }

    mz_calibration_vals <- base::as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features),NA))
    colnames(mz_calibration_vals) <- base::paste("mz_calibration",sample_list,sep=" ")
    for(c in 1:ncol(mz_calibration_vals))
    {
      mz_calibration_vals[,c] <- as.numeric(mz_calibration_vals[,c])
    }

    if(MassSpec_mode == "TIMSToF")
    {
      IM_calibration_vals <- base::as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features),NA))
      colnames(IM_calibration_vals) <- base::paste("IM_calibration",sample_list,sep=" ")
      for(c in 1:ncol(IM_calibration_vals))
      {
        IM_calibration_vals[,c] <- as.numeric(IM_calibration_vals[,c])
      }
    }

    if(any(c(RT_calibration,mz_calibration) == T)) ###any of both calibrations should be done?
    {
      ###first: start by taking already available calibration information

      peptide_features <- which(features$Sequence != "")

      pep_seq <- stringr::str_split(features$Sequence,";",simplify = T)

      if(multiplicity == 1)evidence <- evidence[which(evidence$Raw.file %in% sample_list),]
      if(multiplicity > 1)
      {
        evidence <- evidence[which(evidence$Raw.file %in% base::gsub("_Channel_heavy|_Channel_light|_Channel_medium","",sample_list)),]
      }

      ###evidence_peptides
      evidence_peptides <- unique(evidence$Sequence)
      evidence_peptides_indices <- vector(mode = "list", length = length(evidence_peptides))
      names(evidence_peptides_indices) <- evidence_peptides

      max <- length(evidence_peptides)
      pb <- tcltk::tkProgressBar(title = "Indexing peptides",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      start <- 1
      for(p in evidence_peptides)
      {
        i <- which(evidence_peptides == p)
        evidence_peptides_indices[[p]] <- which(evidence$Sequence == p)

        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)

      ###grap already available information for RT and mz calibration from MaxQ output
      max <- length(peptide_features)
      pb <- tcltk::tkProgressBar(title = "Preparing calibrations per feature",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      counter <- 0
      for(i in peptide_features)
      {
        cur_mods <- stringr::str_split(features$Modifications[i],";|,",simplify = T)
        cur_mods[cur_mods == ""] <- "Unmodified"
        cur_mods <- unique(cur_mods)

        indices <- NULL
        for(p in pep_seq[i,which(pep_seq[i,] != "")])
        {
          indices <- append(indices,evidence_peptides_indices[[p]])
        }
        evidence_sub <- evidence[indices,]
        indices <- which(evidence_sub$Charge == features$Charge[i] & evidence_sub$Modifications %in% cur_mods)

        if(length(indices) > 0)
        {
          if(RT_calibration == T)
          {
            calc_RT <- as.numeric(stringr::str_split(features$Observed_RT[i],pattern = ";",simplify = T)) - features$RT[i]
            data.table::set(RT_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_RT)))
          }
          if(IM_calibration == T)
          {
            calc_IM <- as.numeric(stringr::str_split(features$Observed_IM[i],pattern = ";",simplify = T)) - features$Inv_K0[i]
            data.table::set(IM_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_IM)))
          }
          if(mz_calibration == T)
          {
            if(use_mz_at_max_int_for_correction == T) ###determine mz correction based on observed mz peaks
            {
              calc_mz <- as.numeric(stringr::str_split(features$Observed_mz[i],pattern = ";",simplify = T)) - features$m.z[i]
              data.table::set(mz_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_mz)))
            }else ###use MaxQ correction
            {
              if(multiplicity == 1)
              {
                calc_mz <- stats::aggregate(evidence_sub[indices,"Uncalibrated...Calibrated.m.z..Da."],by=list(Sample=evidence_sub$Raw.file[indices]),FUN="mean",na.rm=T)
                data.table::set(mz_calibration_vals,as.integer(i),as.integer(match(calc_mz$Sample,sample_list)),value = as.list(c(calc_mz$x)))
              }else if(multiplicity > 1)
              {
                evidence_sub[indices,"Uncalibrated...Calibrated.m.z..Da."][is.na(evidence_sub[indices,"Uncalibrated...Calibrated.m.z..Da."])] <- 0
                calc_mz <- stats::aggregate(evidence_sub[indices,"Uncalibrated...Calibrated.m.z..Da."],by=list(Sample=evidence_sub$Raw.file[indices]),FUN="mean",na.rm=T)

                if(multiplicity == 2)
                {
                  names <- c(base::paste(calc_mz$Sample,rep("_Channel_light",length(unique(evidence_sub$Raw.file[indices]))),sep=""),
                             base::paste(calc_mz$Sample,rep("_Channel_heavy",length(unique(evidence_sub$Raw.file[indices]))),sep=""))
                  calc_mz <- calc_mz[rep(1:nrow(calc_mz),multiplicity),]
                  calc_mz$Sample <- names
                }else if(multiplicity == 3)
                {
                  names <- c(base::paste(calc_mz$Sample,rep("_Channel_light",length(unique(evidence_sub$Raw.file[indices]))),sep=""),
                             base::paste(calc_mz$Sample,rep("_Channel_medium",length(unique(evidence_sub$Raw.file[indices]))),sep=""),
                             base::paste(calc_mz$Sample,rep("_Channel_heavy",length(unique(evidence_sub$Raw.file[indices]))),sep=""))
                  calc_mz <- calc_mz[rep(1:nrow(calc_mz),multiplicity),]
                  calc_mz$Sample <- names
                }

                #order by sample_list
                correction_mz <- calc_mz$x[match(calc_mz$Sample,sample_list)]

                #get SILAC-dependent m/z differences
                observed <- as.numeric(stringr::str_split(features$Observed_mz[i],pattern = ";",simplify = T)) - features$m.z[i]

                #final correction
                correction <- observed + correction_mz

                data.table::set(mz_calibration_vals,as.integer(i),as.integer(1:ncol(mz_calibration_vals)),value = as.list(c(correction)))
              }
            }
          }
        }
        counter <- counter + 1
        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(counter/max))*(1-(counter/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, counter, label=base::paste( round(counter/max*100, 0)," % done (",counter,"/",max,", Time require: ",time_require,")",sep = ""))
        }

      }
      close(pb)

      RT_calibration_vals_save <- RT_calibration_vals

      if(mz_calibration == T & use_mz_at_max_int_for_correction==F)
      {
        ####train random forest model for predicting mz calibration

        #suppressWarnings(suppressMessages(library(randomForest,quietly = T)))

        ###train a model per sample
        models <- list()
        trainsets <- list()
        evalsets <- list()
        if(any(colnames(evidence) == "Resolution"))
        {
          temp_data <- evidence[which(rowSums(is.na(evidence[,c("Raw.file","Retention.time","m.z","Charge","Resolution","Uncalibrated...Calibrated.m.z..Da.")])) == 0),c("Raw.file","Retention.time","m.z","Charge","Resolution","Uncalibrated...Calibrated.m.z..Da.")]
          temp_data$Resolution <- 0
        }else
        {
          temp_data <- evidence[which(rowSums(is.na(evidence[,c("Raw.file","Retention.time","m.z","Charge","Uncalibrated...Calibrated.m.z..Da.")])) == 0),c("Raw.file","Retention.time","m.z","Charge","Uncalibrated...Calibrated.m.z..Da.")]
          temp_data$Resolution <- 0
        }

        #replicate rows in case of multiplicity > 1
        if(multiplicity > 1)
        {
          temp_data_temp <- NULL
          for(labels in names(SILAC_settings))
          {
            if(labels == "light")
            {
              temp <- temp_data
              temp$Raw.file <- base::paste(temp$Raw.file,"_Channel_light",sep="")
              temp_data_temp <- temp
            }
            if(labels == "medium")
            {
              temp <- temp_data
              temp$Raw.file <- base::paste(temp$Raw.file,"_Channel_medium",sep="")
              temp_data_temp <- rbind(temp_data_temp,temp)
            }
            if(labels == "heavy")
            {
              temp <- temp_data
              temp$Raw.file <- base::paste(temp$Raw.file,"_Channel_heavy",sep="")
              temp_data_temp <- rbind(temp_data_temp,temp)
            }
          }
          temp_data <- temp_data_temp
        }

        print("Train RF models for mz-corrections")
        for(s in sample_list)
        {
          selec_rows <- which(temp_data$Raw.file == s)
          temp_data2 <- temp_data[selec_rows,]
          random_selection <- sample(1:nrow(temp_data2),size = 0.8*nrow(temp_data2),replace = F)
          trainset <- temp_data2[random_selection,]
          evalset <- temp_data2[-random_selection,]

          if(length(unique(trainset[,"Uncalibrated...Calibrated.m.z..Da."]))>5)
          {
            model <- randomForest::randomForest(trainset[,c("Retention.time","m.z","Charge","Resolution")],trainset[,"Uncalibrated...Calibrated.m.z..Da."] , importance = TRUE,ntree = 500, mtry = 4,do.trace=F,nodesize = 100)
            evalset$predicted <- stats::predict(model, evalset, type = "response")
            trainset$predicted <- stats::predict(model, trainset, type = "response")
            models[[s]] <- model
            trainsets[[s]] <- trainset
            evalsets[[s]] <- evalset

            graphics::plot(trainset$predicted,trainset$Uncalibrated...Calibrated.m.z..Da.,main=base::paste(s,"- Train set (80% of data)"),xlab="Predicted m/z calibration",ylab="MaxQ determined m/z calibration")
            fit <- stats::lm(trainset$Uncalibrated...Calibrated.m.z..Da.~trainset$predicted)
            graphics::abline(fit)
            posx <- graphics::par("usr")[1] + (graphics::par("usr")[2]-graphics::par("usr")[1])*0.8
            posy <- graphics::par("usr")[3] + (graphics::par("usr")[4]-graphics::par("usr")[3])*0.2
            graphics::text(posx,posy,base::paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))
            fit_train <- fit

            graphics::plot(evalset$predicted,evalset$Uncalibrated...Calibrated.m.z..Da.,main=base::paste(s,"- Validation set (20% of data)"),xlab="Predicted m/z calibration",ylab="MaxQ determined m/z calibration")
            fit <- stats::lm(evalset$Uncalibrated...Calibrated.m.z..Da.~evalset$predicted)
            graphics::abline(fit)
            posx <- graphics::par("usr")[1] + (graphics::par("usr")[2]-graphics::par("usr")[1])*0.8
            posy <- graphics::par("usr")[3] + (graphics::par("usr")[4]-graphics::par("usr")[3])*0.2
            graphics::text(posx,posy,base::paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))
            fit_eval <- fit

            print(base::paste(s,": R^2 train-set=",round(as.numeric(summary(fit_train)[8]),digits=2),", R^2 eval-set=",round(as.numeric(summary(fit_eval)[8]),digits=2),sep=""))
          }
        }

        QC_data[["mz_calibration_models"]] <- models
        QC_data[["mz_calibration_train_sets"]] <- trainsets
        QC_data[["mz_calibration_evaluation_sets"]] <- evalsets

        ####now predict mz for missing mzcalibrations of features

        ####get ion indices per feature
        all_ion_indices <- base::as.data.frame(matrix(ncol=6,nrow=nrow(allpeptides)))
        colnames(all_ion_indices) <- c("Feature","Raw.file","Retention.time","m.z","Charge","Resolution")
        all_ion_indices$Feature <- as.numeric(all_ion_indices$Feature)
        all_ion_indices$Raw.file <- as.character(all_ion_indices$Raw.file)
        all_ion_indices$Retention.time <- as.numeric(all_ion_indices$Retention.time)
        all_ion_indices$m.z <- as.numeric(all_ion_indices$m.z)
        all_ion_indices$Charge <- as.integer(all_ion_indices$Charge)
        all_ion_indices$Resolution <- as.numeric(all_ion_indices$Resolution)
        cur_index = 1

        max <- nrow(features)
        pb <- tcltk::tkProgressBar(title = "Prepare features for m/z-calibration prediction by RF",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        for(i in 1:nrow(features))
        {
          ###start with extracting ion_indices per feature
          cur_ion_indices <- as.numeric(stringr::str_split(features$ion_indices[i],",",simplify = T))
          data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),1L,i)
          data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),2L,allpeptides[cur_ion_indices,"Raw.file"])
          data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),as.integer(3:6),value = as.list(allpeptides[cur_ion_indices,c("Retention.time","m.z","Charge","Resolution")]))

          cur_index <- cur_index+length(cur_ion_indices)

          updatecounter <- updatecounter + 1
          if(updatecounter >= 100)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)

        all_ion_indices <- all_ion_indices[1:(cur_index-1),]

        all_ion_indices$Resolution <- 0

        raw_files <- sample_list

        ####get mean ion propperties per feature to predict m/z calibration per feature for each sample
        median_feature_properties <- all_ion_indices[,-2]#[which(all_ion_indices$Raw.file == raw_files[c]),-2]

        ###determine mean values per feature
        median_feature_properties <- stats::aggregate(median_feature_properties[,-1],by=list(Feature=median_feature_properties$Feature),FUN="median",na.rm=T)

        QC_data[["mz_calibration_median_feature_properties"]] <- median_feature_properties

        max <- ncol(mz_calibration_vals)
        pb <- tcltk::tkProgressBar(title = "Predict mz calibrations using RandomForest models",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0

        for(c in 1:ncol(mz_calibration_vals))
        {
          select_model <- which(names(models) == raw_files[c])

          if(length(select_model)>0)
          {
            features_select <- which(is.na(mz_calibration_vals[,c]))

            temp_data <- median_feature_properties[which(median_feature_properties$Feature %in% features_select),]

            temp_data <- temp_data[which(rowSums(is.na(temp_data[,c("Retention.time","m.z","Charge","Resolution")])) == 0),]

            if(nrow(temp_data)> 0)
            {
              prediction <- stats::predict(models[[select_model]], temp_data[,-1], type = "response")
              #add expected SILAC isotope shifts
              if(multiplicity > 1)
              {
                if(grepl("light",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features$m.z._shift_light[temp_data$Feature]
                if(grepl("medium",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features$m.z._shift_medium[temp_data$Feature]
                if(grepl("heavy",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features$m.z._shift_heavy[temp_data$Feature]
                prediction <- prediction + SILAC_label_mz_shift
              }

              data.table::set(mz_calibration_vals,as.integer(temp_data$Feature),as.integer(c),value = prediction)
            }
          }
          updatecounter <- updatecounter + 1
          if(updatecounter >= 1)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(c/max))*(1-(c/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, c, label=base::paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)
      }

      # }else if(mz_calibration == T & use_mz_at_max_int_for_correction==T) ###use observed mz at max intensity calibration
      # {
      #   ###get rowMedians of mz_calibrations
      #   feature_specific_correction_factor <- matrixStats::rowMedians(as.matrix(mz_calibration_vals),na.rm=T)
      #   ###examine deviation of expected feature_specific shift from observed shift
      #   expected_from_observed <- mz_calibration_vals/feature_specific_correction_factor
      #   ###sample specific shift
      #   sample_specific_correction_factor <- matrixStats::colMedians(log2(as.matrix(expected_from_observed)),na.rm=T)
      #   ###examine deviation of expected feature_specific shift from observed shift
      #   expected_from_observed <- mz_calibration_vals/((2^sample_specific_correction_factor)*feature_specific_correction_factor)
      #   ###plot per sample how well expected mz correction correlates with observed correction
      #   for(c in 1:length(sample_list))
      #   {
      #     sel <- which(abs(mz_calibration_vals[,c]) < 0.01)
      #     graphics::smoothScatter(mz_calibration_vals[sel,c],(2^sample_specific_correction_factor[c])*feature_specific_correction_factor[sel],xlim=c(-0.005,0.005),ylim=c(-0.005,0.005),main=sample_list[c],xlab="Observed mz - feature mz",ylab="Expected mz - feature mz")
      #     graphics::abline(a=0,b=1)
      #     fit <- stats::lm((2^sample_specific_correction_factor[c])*feature_specific_correction_factor[sel]~mz_calibration_vals[sel,c])
      #     graphics::abline(fit,lty=2)
      #     graphics::text(0,0.004,base::paste("Rsq:",round(as.numeric(summary(fit)[8]),digits=2)))
      #   }
      #   ###store calibration factor
      #   QC_data[["mz_calibration_feature_specific_cor_factor"]] <- feature_specific_correction_factor
      #   QC_data[["mz_calibration_sample_specific_cor_factor"]] <- sample_specific_correction_factor
      #
      #   ###replace missing values with expected correction factors
      #   for(c in 1:ncol(mz_calibration_vals))
      #   {
      #     selection <- which(is.na(mz_calibration_vals[,c]))
      #     mz_calibration_vals[,c][selection] <- ((2^sample_specific_correction_factor[c])*feature_specific_correction_factor[selection])
      #   }
      # }

      #if there are still missing values just fill with 0 or expected SILAC shifts
      if(any(is.na(mz_calibration_vals)))
      {
        if(multiplicity == 1)
        {
          mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
        }else
        {
          for(c in 1:ncol(mz_calibration_vals))
          {
            if(any(is.na(mz_calibration_vals[,c])))
            {
              if(grepl("light",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_light[is.na(mz_calibration_vals[,c])]
              if(grepl("medium",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_medium[is.na(mz_calibration_vals[,c])]
              if(grepl("heavy",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_heavy[is.na(mz_calibration_vals[,c])]
            }
          }
        }
      }

      if(RT_calibration == T)
      {
        graphics::plot(1, type="n", axes=F, xlab="", ylab="")
        graphics::text(1,1,"RT GAM fitting - First round of estimation")

        RT_alignment_GAM_models <- list()
        ###fit GAM and predict per sample.
        for(c in 1:ncol(RT_calibration_vals))
        {
          x <- features$RT[which(RT_calibration_vals[,c] != 0)]
          y <- RT_calibration_vals[,c][which(RT_calibration_vals[,c] != 0)]
          if(length(x) > 10) ###require at least 500 data points to perform fitting
          {
            ##try to fit an average generalised additive model to determine a RT dependent calibration curve
            gam <- mgcv::gam(y ~ s(x), method = "REML")
            RT_alignment_GAM_models[[c]] <- gam
            x_pred <- seq(min(features$RT,na.rm=T), max(features$RT,na.rm=T), length.out = nrow(features))
            y_pred <- stats::predict(gam, base::data.frame(x = x_pred))
            lim_stats <- grDevices::boxplot.stats(y_pred)
            delta_y <- lim_stats$stats[4] - lim_stats$stats[2]
            ylim <- c(-max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))-delta_y,max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))+delta_y)
            if(ylim[2] < 2)ylim <- c(-2,2)
            graphics::smoothScatter(x,y,ylab="Observed RT calibration",main=sample_list[c],xlab="RT [min]",ylim=ylim)
            graphics::lines(x_pred,y_pred,col="red")
            graphics::legend("topright",legend="GAM",lty=1,col="red")

            ##predict RT correction for all features in current sample
            y_pred <- stats::predict(gam, base::data.frame(x = features$RT))

            ###use GAM to predict RT correction for missing features in current sample
            selection <- which(is.na(RT_calibration_vals[,c]))
            RT_calibration_vals[,c][selection] <- y_pred[selection]

          }else ###not enough observations ... skip
          {
            print(base::paste(sample_list[c],"- Not enough peptide observations for RT-GAM fitting..."))
          }
        }
        QC_data[["RT_calibration_GAM_models"]] <- RT_alignment_GAM_models

        #perform a second round of RT calibration
        #currently it is possible that several features are observed only in a few samples
        #hence observed median RT is driven by these few samples
        #based on GAM models adjust median feature RT and rerun fitting
        expected_raw_RTs <- features$RT + RT_calibration_vals
        expected_raw_RTs[is.na(expected_raw_RTs)] <- 0

        delta_rt_temp <- features$RT_range_max-features$RT
        features$RT <- matrixStats::rowMedians(as.matrix(expected_raw_RTs),na.rm=T)
        features$RT_range_min <- features$RT - delta_rt_temp
        features$RT_range_max <- features$RT + delta_rt_temp

        #restart RT calibration
        #grap known RT once again
        RT_calibration_vals <- RT_calibration_vals_save

        graphics::plot(1, type="n", axes=F, xlab="", ylab="")
        graphics::text(1,1,"RT GAM fitting - Second round of estimation")

        RT_alignment_GAM_models <- list()
        ###fit GAM and predict per sample.
        for(c in 1:ncol(RT_calibration_vals))
        {
          x <- features$RT[which(RT_calibration_vals[,c] != 0)]
          y <- RT_calibration_vals[,c][which(RT_calibration_vals[,c] != 0)]
          if(length(x) > 10) ###require at least 500 data points to perform fitting
          {
            ##try to fit an average generalised additive model to determine a RT dependent calibration curve
            gam <- mgcv::gam(y ~ s(x), method = "REML")
            RT_alignment_GAM_models[[c]] <- gam
            x_pred <- seq(min(features$RT,na.rm=T), max(features$RT,na.rm=T), length.out = nrow(features))
            y_pred <- stats::predict(gam, base::data.frame(x = x_pred))
            lim_stats <- grDevices::boxplot.stats(y_pred)
            delta_y <- lim_stats$stats[4] - lim_stats$stats[2]
            ylim <- c(-max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))-delta_y,max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))+delta_y)
            if(ylim[2] < 2)ylim <- c(-2,2)
            graphics::smoothScatter(x,y,ylab="Observed RT calibration",main=sample_list[c],xlab="RT [min]",ylim=ylim)
            graphics::lines(x_pred,y_pred,col="red")
            graphics::legend("topright",legend="GAM",lty=1,col="red")

            ##predict RT correction for all features in current sample
            y_pred <- stats::predict(gam, base::data.frame(x = features$RT))

            ###use GAM to predict RT correction for missing features in current sample
            selection <- which(is.na(RT_calibration_vals[,c]))
            RT_calibration_vals[,c][selection] <- y_pred[selection]
          }else ###not enough observations ... skip
          {
            print(base::paste(sample_list[c],"- Not enough peptide observations for RT-GAM fitting..."))
          }
        }
        QC_data[["RT_calibration_GAM_models"]] <- RT_alignment_GAM_models
      }
      RT_calibration_vals[is.na(RT_calibration_vals)] <- 0

      if(MassSpec_mode == "TIMSToF" & IM_calibration == T)
      {
        ###use median of IM_calibration per sample for NAs
        for(c in 1:ncol(IM_calibration_vals))
        {
          IM_calibration_vals[,c][is.na(IM_calibration_vals[,c])] <- stats::median(IM_calibration_vals[,c],na.rm=T)
        }
      }

    }else ###no calibration should be done
    {
      RT_calibration_vals[is.na(RT_calibration_vals)] <- 0

      if(multiplicity == 1)
      {
        mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
      }else if(multiplicity > 1)
      {
        for(c in 1:ncol(mz_calibration_vals))
        {
          if(any(is.na(mz_calibration_vals[,c])))
          {
            if(grepl("light",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_light[is.na(mz_calibration_vals[,c])]
            if(grepl("medium",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_medium[is.na(mz_calibration_vals[,c])]
            if(grepl("heavy",colnames(mz_calibration_vals)[c]))mz_calibration_vals[is.na(mz_calibration_vals[,c]),c] <- features$m.z._shift_heavy[is.na(mz_calibration_vals[,c])]
          }
        }
      }
      if(MassSpec_mode == "TIMSToF")IM_calibration_vals[is.na(IM_calibration_vals)] <- 0
    }

    features <- cbind(features,RT_calibration_vals,mz_calibration_vals)
    if(MassSpec_mode == "TIMSToF")features <- cbind(features,IM_calibration_vals)

    utils::write.table(features,base::paste("Temporary_files/Features_aligned_merged",output_file_names_add,".txt",sep=""),row.names = F)

    save(QC_data,file = base::paste("Temporary_files/Feature_alignment_QC_data.RData",sep=""))

    grDevices::dev.off()
    options(warn=0)
    crap <- gc(F)
  }
}

#' Extend IceR features by expected +1-isotopic features
#' @param path_to_features Path to folder where results of align_features() are stored.
#' @param feature_table_file_name File name which contains align_features() results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param min_observations Specifying how many MaxQuant features had to be previously detected for monoisotopic IceR features to be extended by +1-isotopic features. By default set to 0 indicating that all for all monoisotopic IceR features a +1-isotopic feature is added.
#' @details Optional step of the IceR workflow appending the list of IceR features by respective expected +1-isotopic features. Isotopic features are expected to show an m/z shift by +1.002 Da/z relative to the monoisotopic IceR feature. Only IceR features with at least min_observations and with mean Andromeda score > 25 % quantile of all mean Andromeda scores are considered.
#' @return Extended list ist stored in the specified feature_table_file_name.
#' @export
add_isotope_features <- function(path_to_features,feature_table_file_name="Features_aligned_merged_IceR_analysis.txt",min_observations=0)
{
  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  features <- utils::read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  if(!any(grepl("_i",features$Feature_name)))
  {
    ###isotope + 1 features
    isotope_features <- features
    isotope_features <- isotope_features[which(!grepl(";",isotope_features$Sequence)),]
    isotope_features <- isotope_features[which(as.numeric(as.character(isotope_features$num_matches)) >= min_observations),]

    if(nrow(isotope_features)>0)
    {
      quantile25_score <- grDevices::boxplot.stats(as.numeric(as.character(isotope_features$mean_Scores)))$stats[2]
      isotope_features <- isotope_features[which(as.numeric(as.character(isotope_features$mean_Scores)) >= quantile25_score),]
      isotope_features$m.z <- ((isotope_features$m.z*isotope_features$Charge)+1.002054)/isotope_features$Charge
      delta_mz <- isotope_features$m.z_range_max-isotope_features$m.z_range_min
      isotope_features$m.z_range_max <- isotope_features$m.z+(delta_mz/2)
      isotope_features$m.z_range_min <- isotope_features$m.z-(delta_mz/2)
      isotope_features$Feature_name <- base::paste(isotope_features$Feature_name,"_i",sep="")
      features <- rbind(features,isotope_features)
      utils::write.table(features,base::paste("Temporary_files/",feature_table_file_name,sep=""),row.names = F)
      print(base::paste("Added",nrow(isotope_features),"isotope features."))
    }else
    {
      print(base::paste("None of the peptides was observed in at least",min_observations,"samples. No isotope features were added."))
    }
  }else
  {
    print("Isotope features were already added")
  }
}


#' Extend IceR features by expected monoisotopic features which correspond to totally missed tryptic peptides of identified proteins.
#' @param path_to_features Path to folder where results of align_features() are stored.
#' @param feature_table_file_name File name which contains align_features() results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param path_to_fasta Path to fasta file used during MaxQuant data preprocessing.
#' @param min_AA_length Minimal length of a theoretical peptide which is considered to be added. By default set to 7.
#' @param max_AA_length Maximal length of a theoretical peptide which is considered to be added. By default set to 25.
#' @param max_add_per_protein Maximal number of potentially missed peptide features which should be added. By default set to 10.
#' @param min_observations Minimal number of MaxQuant features which have to fall into the expected RT- and m/z-windows of potentially missed peptide features.
#' @details Optional step of the IceR workflow appending the list of IceR features by potentially missed tryptic peptide IceR features. Requires installation of the Bioconductor package "cleaver". Expected m/z windows are determined based on masses of respective tryptic peptides. Expected RT windows are estimated based on random forest modelling with KyteDoolittle hydrophobicity, peptide length, peptide charges and amino acid composition as predictors and observed RTs as response factors. Potentially missed peptide (PMP) features are expected to show +2 or +3 charge.
#' @return Extended list ist stored in the specified feature_table_file_name.
add_missed_peptides <- function(path_to_features,feature_table_file_name="Features_aligned_merged.txt",path_to_fasta,min_AA_length=7,max_AA_length=25,max_add_per_protein=10,min_observations=5)
{
  ###to do
  #add expected mz shifts for SILAC (m.z._shift_light,m.z._shift_medium,m.z._shift_heavy) to features_unknown if multiplicity > 1

  options(warn = -1)
  #suppressWarnings(suppressMessages(library(Peptides,quietly = T)))
  #suppressWarnings(suppressMessages(library(randomForest,quietly = T)))
  #suppressWarnings(suppressMessages(library(seqinr,quietly = T)))
  #suppressWarnings(suppressMessages(library(cleaver,quietly = T)))
  RT_calibration=T
  mz_calibration=T

  grDevices::pdf("Temporary_files/Potentially_missed_peptides_QC_plots.pdf")

  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  features <- utils::read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  features <- features[which(!grepl("_i",features$Feature_name)),]

  all_peptides <- unique(as.character(stringr::str_split(features$Sequence,";",simplify = T)))
  all_proteins <- unique(features$Protein)
  all_proteins <- all_proteins[which(all_proteins != "" & !grepl(";|\\|",all_proteins))]

  features <- features[which(!grepl(";|\\|",features$Sequence)),]

  ###filter for high qualtiy identifications
  features <- features[which(as.numeric(as.character(features$mean_Scores)) >= 50),]

  ###summarize relevant parameters
  print("Prepare data for modelling peptide RT")
  aaComp <- aaComp(features$Sequence)
  aaComp <- base::data.frame(matrix(unlist(aaComp), nrow=length(aaComp), byrow=T))[,1:9]
  colnames(aaComp) <- c("Tiny","Small","Aliphatic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")

  peptide_parameters <- base::data.frame(hydrophobicities=Peptides::hydrophobicity(features$Sequence),
                                   length=nchar(as.character(features$Sequence)),
                                   charge=as.numeric(as.character(features$Charge)),
                                   aaComp,
                                   RT=as.numeric(as.character(features$RT)))
  ###Train random forest model on the parameters
  random_selection <- sample(1:nrow(peptide_parameters),size = 0.8*nrow(peptide_parameters),replace = F)
  trainset <- peptide_parameters[random_selection,]
  evalset <- peptide_parameters[-random_selection,]
  print("Perform RF modelling for peptide RT")
  model <- randomForest::randomForest(trainset[,1:12],trainset[,13], importance = TRUE,ntree = 200, mtry = 4,do.trace=F,nodesize = floor(nrow(trainset)/1000))

  evalset$predicted <- stats::predict(model, evalset[,1:12], type = "response")
  trainset$predicted <- stats::predict(model, trainset[,1:12], type = "response")

  graphics::smoothScatter(trainset$predicted,trainset$RT,main="Train set (80% of data)",xlab="Predicted RT",ylab="Observed RT")
  fit <- stats::lm(trainset$RT~trainset$predicted)
  graphics::abline(fit)
  posx <- graphics::par("usr")[1] + (graphics::par("usr")[2]-graphics::par("usr")[1])*0.8
  posy <- graphics::par("usr")[3] + (graphics::par("usr")[4]-graphics::par("usr")[3])*0.2
  graphics::text(posx,posy,base::paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))

  graphics::smoothScatter(evalset$predicted,evalset$RT,main="Eval set (80% of data)",xlab="Predicted RT",ylab="Observed RT")
  fit <- stats::lm(evalset$RT~evalset$predicted)
  graphics::abline(fit)
  posx <- graphics::par("usr")[1] + (graphics::par("usr")[2]-graphics::par("usr")[1])*0.8
  posy <- graphics::par("usr")[3] + (graphics::par("usr")[4]-graphics::par("usr")[3])*0.2
  graphics::text(posx,posy,base::paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))

  deviation_predicted_true_RT <- matrixStats::rowSds(cbind(evalset$RT,evalset$predicted),na.rm=T)
  RT_deviation_cutoff <- grDevices::boxplot.stats(deviation_predicted_true_RT)$stats[4]

  ###Perform in-silico digestion of all identified proteins
  fasta <- seqinr::read.fasta(path_to_fasta,seqtype = "AA",as.string = T)
  Sequences <- base::data.frame(matrix(unlist(fasta), nrow=length(fasta), byrow=T))
  colnames(Sequences) <- "Sequence"
  Sequences$Sequence <- as.character(Sequences$Sequence)
  fasta_headers_boundaries_uid <- gregexpr("\\|",names(fasta))
  fasta_headers_boundaries_uid <- base::data.frame(matrix(unlist(fasta_headers_boundaries_uid), nrow=length(fasta_headers_boundaries_uid), byrow=T))
  UIDs <- base::substr(names(fasta),start = fasta_headers_boundaries_uid[,1]+1,stop = fasta_headers_boundaries_uid[,2]-1)
  rownames(Sequences) <- UIDs
  Sequences <- Sequences[all_proteins,,drop=F]

  insilicodigest <- cleaver::cleave(Sequences$Sequence,"trypsin")
  names(insilicodigest) <- rownames(Sequences)

  ###Summarize potential peptides with a length of at least n and at max m amino acids
  potential_new_feature_peptides <- vector(mode="character",length = nrow(features))
  potential_new_feature_peptides_ID <- vector(mode="character",length = nrow(features))
  count <- 1
  max <- length(insilicodigest)
  pb <- tcltk::tkProgressBar(title = "Determine potential missed peptide features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0
  for(i in 1:length(insilicodigest))
  {
    temp <- insilicodigest[[i]]
    temp <- temp[-1] ###remove the first peptide as it contains signal peptide sequence which is cleaved off and thus will never be detectable
    temp <- temp[which(temp %not in% all_peptides & nchar(temp) >= min_AA_length & nchar(temp) <= max_AA_length)]
    if(length(temp)>0)
    {
      if(length(temp)>max_add_per_protein)temp <- temp[1:max_add_per_protein]
      potential_new_feature_peptides[count:(count+length(temp)-1)] <- temp
      potential_new_feature_peptides_ID[count:(count+length(temp)-1)] <- names(insilicodigest)[i]
      count <- count+length(temp)
    }

    updatecounter <- updatecounter + 1
    if(updatecounter >= 100)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)
  Seq_to_ID <- base::data.frame(Sequence=potential_new_feature_peptides,ID=potential_new_feature_peptides_ID)
  potential_new_feature_peptides <- unique(potential_new_feature_peptides)
  wrong_seq <- which(regexpr("U",potential_new_feature_peptides) != -1)
  if(length(wrong_seq)>0)potential_new_feature_peptides <- potential_new_feature_peptides[-wrong_seq]
  ###Determine Mass and RT per potential new peptide feature
  aaComp <- aaComp(potential_new_feature_peptides)
  aaComp <- base::data.frame(matrix(unlist(aaComp), nrow=length(aaComp), byrow=T))[,1:9]
  colnames(aaComp) <- c("Tiny","Small","Aliphatic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")
  new_peptide_parameters <- base::data.frame(hydrophobicities=Peptides::hydrophobicity(potential_new_feature_peptides),
                                       length=nchar(potential_new_feature_peptides),
                                       charge=2,
                                       aaComp)

  new_peptide_parameters$predicted <- stats::predict(model, new_peptide_parameters[,1:12], type = "response")

  potential_new_feature_peptides <- base::data.frame(Sequence=potential_new_feature_peptides,
                                               ID=Seq_to_ID$ID[match(potential_new_feature_peptides,Seq_to_ID$Sequence)],
                                               Mass=Peptides::mw(potential_new_feature_peptides, monoisotopic = T),
                                               RT=new_peptide_parameters$predicted)

  ###Now try to see if corresponding features with +2 or +3 charge wera already detected by MaxQ
  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  load("allpeptides.RData")
  setwd(path_to_features)

  select_dataframe_rows = function(ds, sel) {
    cnames = colnames(ds)
    rnames = rownames(ds)
    ds = base::data.frame(ds[sel,])
    colnames(ds) = cnames
    rownames(ds) = rnames[sel]
    return (ds)
  }

  allpeptides <- select_dataframe_rows(allpeptides,which(allpeptides$Sequence == " " & allpeptides$Charge %in% c(2,3) & allpeptides$Modifications == " " | allpeptides$Sequence == "" & allpeptides$Charge %in% c(2,3) & allpeptides$Modifications == " "))
  allpeptides <- allpeptides[order(allpeptides$Retention.time),]
  allpeptides <- select_dataframe_rows(allpeptides,which(allpeptides$Retention.time >= min(potential_new_feature_peptides$RT) & allpeptides$Retention.time <= max(potential_new_feature_peptides$RT)))
  allpeptides$allpeptides_rowname <- rownames(allpeptides)

  ###Indexing dat by RT windows to improve speed for subsetting
  indexing_RT_window <- 0.5

  num_windows <- ceiling(ceiling(max(allpeptides$Retention.time,na.rm=T)-min(allpeptides$Retention.time,na.rm=T))*(1/indexing_RT_window))

  Indices <- base::as.data.frame(matrix(ncol=4,nrow=num_windows))
  colnames(Indices) <- c("RT_start","RT_end","Row_start","Row_end")
  Indices$RT_start <- as.numeric(Indices$RT_start)
  Indices$RT_end <- as.numeric(Indices$RT_end)
  Indices$Row_start <- as.numeric(Indices$Row_start)
  Indices$Row_end <- as.numeric(Indices$Row_end)

  print("Indexing MaxQ features")
  max <- nrow(Indices)
  pb <- tcltk::tkProgressBar(title = "Indexing MaxQ features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0
  start <- 1
  end <- max

  cur_index <- 1
  start_RT <- floor(min(allpeptides$Retention.time,na.rm=T))

  ###presubset allpeptides into RT windows
  incr <- 1000
  for(i in 1:nrow(Indices))
  {
    start <- (i-1)*indexing_RT_window+start_RT
    end <- (i*indexing_RT_window)+start_RT

    indx <- cur_index
    indx_prev <- cur_index
    while(T)
    {
      indx <- indx + incr

      if(allpeptides$Retention.time[indx] >= end | indx > nrow(allpeptides))
      {
        break
      }else
      {
        indx_prev <- indx
      }
    }
    if(indx > nrow(allpeptides))indx <- nrow(allpeptides)

    #inds <- which(allpeptides$Retention.time > start & allpeptides$Retention.time < end)
    if(indx > indx_prev)
    {
      temp <- allpeptides[indx_prev:indx,]

      if(length(which(temp$Retention.time < end)) > 0)
      {
        inds <- cur_index:(indx_prev-1+max(which(temp$Retention.time < end)))

        incr <- max(inds) + 1 - cur_index

        cur_index <- max(inds) + 1
        if(length(inds)>0)
        {
          data.table::set(Indices,i = as.integer(i),j=as.integer(1:4),value=as.list(c(start,end,min(inds),max(inds))))
        }
      }
    }


    updatecounter <- updatecounter + 1
    if(updatecounter >= 1)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    }

  }
  close(pb)

  if(any(rowSums(is.na(Indices)) > 0))
  {
    Indices <- Indices[-which(rowSums(is.na(Indices)) > 0),]
  }

  data_frags <- list()
  for(i in 1:nrow(Indices))
  {
    data_frags[[i]] <- allpeptides[Indices$Row_start[i]:Indices$Row_end[i],]
  }

  ###now search for potential fitting MaxQ features, if MaxQ feature mass fits to potential peptide and more than one species (e.g. two and three charges)
  ###where detected, take this charge for which more MaxQ features were found
  QC_data <- NULL
  load("Temporary_files/Feature_alignment_QC_data.RData")
  borders_RT <- QC_data$Feature_alignment_windows$RT_window
  borders_m.z <- QC_data$Feature_alignment_windows$mz_window

  features_unknown <- base::as.data.frame(matrix(ncol=21,nrow=nrow(potential_new_feature_peptides)))
  colnames(features_unknown) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")
  features_unknown <- data.table::as.data.table(features_unknown)
  features_unknown$Feature_name <- as.character(features_unknown$Feature_name)
  features_unknown$Mass <- as.numeric(features_unknown$Mass)
  features_unknown$RT <- as.numeric(features_unknown$RT)
  features_unknown$Sequence <- as.character(features_unknown$Sequence)
  features_unknown$Protein <- as.character(features_unknown$Protein)
  features_unknown$MSMS.Scan.Numbers <- as.character(features_unknown$MSMS.Scan.Numbers)
  features_unknown$m.z_range_min <- as.numeric(features_unknown$m.z_range_min)
  features_unknown$m.z_range_max <- as.numeric(features_unknown$m.z_range_max)
  features_unknown$RT_range_min <- as.numeric(features_unknown$RT_range_min)
  features_unknown$RT_range_max <- as.numeric(features_unknown$RT_range_max)
  features_unknown$ion_indices <- as.character(features_unknown$ion_indices)
  features_unknown$count_ion_indices <- as.numeric(features_unknown$count_ion_indices)
  features_unknown$m.z <- as.numeric(features_unknown$m.z)
  features_unknown$Charge <- as.numeric(features_unknown$Charge)
  features_unknown$mean_Scores <- as.character(features_unknown$mean_Scores)
  features_unknown$num_matches <- as.character(features_unknown$num_matches)
  features_unknown$Modifications <- as.character(features_unknown$Modifications)
  features_unknown$RT_length <- as.numeric(features_unknown$RT_length)
  features_unknown$Observed_RT <- as.character(features_unknown$Observed_RT)
  features_unknown$Observed_mz <- as.character(features_unknown$Observed_mz)
  features_unknown$Observed_score <- as.character(features_unknown$Observed_score)

  count_features <- 0
  total_count_features <- as.numeric(base::gsub("Feature_","",features$Feature_name[nrow(features)]))

  sample_list <- base::gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration",colnames(features)))])

  print("Search for potential missed peptide features")
  max <- nrow(potential_new_feature_peptides)
  pb <- tcltk::tkProgressBar(title = "Search for missed peptide features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 350)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0
  for(i in 1:nrow(potential_new_feature_peptides))
  {
    RT_window <- c(potential_new_feature_peptides$RT[i]-(RT_deviation_cutoff/2),potential_new_feature_peptides$RT[i]+(RT_deviation_cutoff/2))

    if(RT_window[1] <= max(Indices$RT_start) & RT_window[2] <= max(Indices$RT_end))
    {
      search_range <- (which(Indices$RT_start>=RT_window[1])[1]-1):(which(Indices$RT_end>=RT_window[2])[1])
    }else
    {
      search_range <- nrow(Indices)
    }
    if(any(search_range < 1))
    {
      search_range <- search_range[-which(search_range < 1)]
    }

    if(length(search_range)>0)
    {
      sub <- data.table::rbindlist(data_frags[search_range])

      sub <- sub[which(sub$Mass >= (potential_new_feature_peptides$Mass[i] - 0.001) & sub$Mass <= (potential_new_feature_peptides$Mass[i] + 0.001)),]

      if(length(unique(sub$Raw.file))>=min_observations)
      {
        median_RT <- stats::median(sub$Retention.time,na.rm=T)
        median_mz <- stats::median(sub$m.z,na.rm=T)
        count_charges <- plyr::count(sub$Charge)
        max_count_charge <- count_charges[which(count_charges$freq == max(count_charges$freq))[1],1]

        sub <- sub[which(sub$m.z >= median_mz+borders_m.z[1] & sub$m.z <= median_mz+borders_m.z[2] & sub$Retention.time >= median_RT+borders_RT[1] & sub$Retention.time <= median_RT+borders_RT[2] & sub$Charge == max_count_charge),]
        if(length(unique(sub$Raw.file))>=min_observations)
        {
          count_features <- count_features + 1
          total_count_features <- total_count_features + 1
          name <- base::paste("Feature",total_count_features,"pmp",sep="_")

          sequence <- as.character(potential_new_feature_peptides$Sequence[i])
          protein <- as.character(potential_new_feature_peptides$ID[i])
          msmsscan <- ""

          median_mass <- stats::median(sub$Mass,na.rm=T)
          Charge <- max_count_charge

          scores <- 0
          num_matches <- length(unique(sub$Raw.file))
          Modifications <- ""

          RT_length <- stats::median(sub$Retention.Length,na.rm=T)

          temp_RT <- stats::aggregate(sub$Retention.time,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_RT <- base::paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

          temp_mz <- stats::aggregate(sub$isotope_corrected_mz_at_max_int,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_mz <- base::paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

          temp_score <- stats::aggregate(sub$Score,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_score <- base::paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed score according to sample_list

          matches <- sub$allpeptides_rowname

          data.table::set(x = features_unknown,i = as.integer(count_features),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
                                                                                                                                    sequence,
                                                                                                                                    protein,
                                                                                                                                    msmsscan,
                                                                                                                                    base::paste(matches,collapse = ","),
                                                                                                                                    base::paste(scores,collapse=";"),
                                                                                                                                    base::paste(num_matches,collapse=";"),
                                                                                                                                    Modifications,
                                                                                                                                    Observed_RT,
                                                                                                                                    Observed_mz,
                                                                                                                                    Observed_score
          )))
          data.table::set(x = features_unknown,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
            as.numeric(median_mz),
            as.numeric(median_RT),
            as.numeric(median_mz+borders_m.z[1]),
            as.numeric(median_mz+borders_m.z[2]),
            as.numeric(median_RT+borders_RT[1]),
            as.numeric(median_RT+borders_RT[2]),
            as.numeric(nrow(sub)),
            as.numeric(median_mass),
            as.numeric(Charge),
            RT_length
          )))

        }

      }
    }

    updatecounter <- updatecounter + 1
    if(updatecounter >= 10)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,") Found: ",count_features,sep = ""))
    }

  }
  close(pb)

  features_unknown <- features_unknown[1:count_features,]
  print(base::paste("Detected",count_features,"potentially missed peptides (PMPs)"))
  features_unknown$num_diff_charges <- 1

  ####now add mz and RT calibration information per feature

  RT_calibration_vals <- base::as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features_unknown),NA))
  colnames(RT_calibration_vals) <- base::paste("RT_calibration",sample_list,sep=".")
  for(c in 1:ncol(RT_calibration_vals))
  {
    RT_calibration_vals[,c] <- as.numeric(RT_calibration_vals[,c])
  }

  mz_calibration_vals <- base::as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features_unknown),NA))
  colnames(mz_calibration_vals) <- base::paste("mz_calibration",sample_list,sep=".")
  for(c in 1:ncol(mz_calibration_vals))
  {
    mz_calibration_vals[,c] <- as.numeric(mz_calibration_vals[,c])
  }


  #determine multiplicity
  if(any(grepl("_Channel_light|_Channel_medium|_Channel_heavy",sample_list)))
  {
    if(any(grepl("_Channel_medium",sample_list)))
    {
      multiplicity <- 3
    }else
    {
      multiplicity <- 2
    }
  }else
  {
    multiplicity <- 1
  }

  ###use observed RTs and fill missing values with median RT calibrations per sample
  if(RT_calibration == T)
  {
    for(i in 1:nrow(features_unknown))
    {
      calc_RT <- as.numeric(stringr::str_split(features_unknown$Observed_RT[i],pattern = ";",simplify = T)) - features_unknown$RT[i]
      data.table::set(RT_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_RT)))
    }

    ###fill missing values with median RT_calibration per sample
    for(c in 1:ncol(RT_calibration_vals))
    {
      RT_calibration_vals[,c][is.na(RT_calibration_vals[,c])] <- stats::median(RT_calibration_vals[,c],na.rm=T)
    }

    ###if anything is left with NA fill with 0
    RT_calibration_vals[is.na(RT_calibration_vals)] <- 0
  }else ###no calibration
  {
    RT_calibration_vals[is.na(RT_calibration_vals)] <- 0
  }
  ### use observed m/z calibrated to uncalibrated values per sample
  ### fill missing values with trained m/z calibrations models per sample
  if(mz_calibration == T)
  {
    ###fill with known deviations
    max <- nrow(features_unknown)
    pb <- tcltk::tkProgressBar(title = "Preparing m/z calibrations per missed peptide feature",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(features_unknown))
    {
      observed_RT <- as.numeric(stringr::str_split(features_unknown$Observed_RT[i],";",simplify = T))
      RT_window <- c(min(observed_RT,na.rm=T),max(observed_RT,na.rm=T))
      if(RT_window[1] <= max(Indices$RT_start) & RT_window[2] <= max(Indices$RT_end))
      {
        search_range <- (which(Indices$RT_start>=RT_window[1])[1]-1):(which(Indices$RT_end>=RT_window[2])[1])
      }else
      {
        search_range <- nrow(Indices)
      }
      if(any(search_range < 1))
      {
        search_range <- search_range[-which(search_range < 1)]
      }

      if(length(search_range)>0)
      {
        sub <- data.table::rbindlist(data_frags[search_range])
        temp <- sub[which(sub$allpeptides_rowname %in% as.character(stringr::str_split(features_unknown$ion_indices[i],",",simplify = T))),]
        temp$Uncalibrated...Calibrated.m.z..Da. <- temp$Uncalibrated.m.z-temp$m.z
        calc_mz <- stats::aggregate(temp$Uncalibrated...Calibrated.m.z..Da.,by=list(Sample=temp$Raw.file),FUN="mean",na.rm=T)
        data.table::set(mz_calibration_vals,as.integer(i),as.integer(match(calc_mz$Sample,sample_list)),value = as.list(c(calc_mz$x)))

      }

      updatecounter <- updatecounter + 1
      if(updatecounter >= 10)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- lubridate::seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }

    }
    close(pb)

    ####now predict mz for missing mzcalibrations of features
    setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
    load("allpeptides.RData")
    setwd(path_to_features)

    ###all following steps are only required to get the median resolution per feature
    ###we require median RT, mz, resolution and charge as a rough estimation of the corresponding parameters in the samples for which we don?t have the respective feature identified
    ####get ion indices per feature
    all_ion_indices <- base::as.data.frame(matrix(ncol=6,nrow=nrow(allpeptides)))
    colnames(all_ion_indices) <- c("Feature","Raw.file","Retention.time","m.z","Charge","Resolution")
    all_ion_indices$Feature <- as.numeric(all_ion_indices$Feature)
    all_ion_indices$Raw.file <- as.character(all_ion_indices$Raw.file)
    all_ion_indices$Retention.time <- as.numeric(all_ion_indices$Retention.time)
    all_ion_indices$m.z <- as.numeric(all_ion_indices$m.z)
    all_ion_indices$Charge <- as.integer(all_ion_indices$Charge)
    all_ion_indices$Resolution <- as.numeric(all_ion_indices$Resolution)
    cur_index = 1

    max <- nrow(features_unknown)
    pb <- tcltk::tkProgressBar(title = "Prepare missed peptide features for m/z-calibration prediction by RF",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(features_unknown))
    {
      ###start with extracting ion_indices per feature
      cur_ion_indices <- as.numeric(stringr::str_split(features_unknown$ion_indices[i],",",simplify = T))
      data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),1L,i)
      data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),2L,allpeptides[cur_ion_indices,"Raw.file"])
      data.table::set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),as.integer(3:6),value = as.list(allpeptides[cur_ion_indices,c("Retention.time","m.z","Charge","Resolution")]))

      cur_index <- cur_index+length(cur_ion_indices)

      updatecounter <- updatecounter + 1
      if(updatecounter >= 100)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- lubridate::seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    all_ion_indices <- all_ion_indices[1:(cur_index-1),]

    raw_files <- sample_list

    ####get mean ion propperties per feature to predict m/z calibration per feature for each sample
    median_feature_properties <- all_ion_indices[,-2]#[which(all_ion_indices$Raw.file == raw_files[c]),-2]

    ###determine mean values per feature
    median_feature_properties <- stats::aggregate(median_feature_properties[,-1],by=list(Feature=median_feature_properties$Feature),FUN="median",na.rm=T)

    ###perform prediction

    models <- QC_data$mz_calibration_models

    max <- ncol(mz_calibration_vals)
    pb <- tcltk::tkProgressBar(title = "Predict mz calibrations using RandomForest models",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(c in 1:ncol(mz_calibration_vals))
    {
      select_model <- which(names(models) == raw_files[c])

      if(length(select_model)>0)
      {
        features_select <- which(is.na(mz_calibration_vals[,c]))

        temp_data <- median_feature_properties[which(median_feature_properties$Feature %in% features_select),]
        temp_data$Resolution <- 0
        temp_data <- temp_data[which(rowSums(is.na(temp_data[,c("Retention.time","m.z","Charge","Resolution")])) == 0),]

        if(nrow(temp_data)> 0)
        {
          prediction <- stats::predict(models[[select_model]], temp_data[,-1], type = "response")

          if(any(is.na(prediction)))prediction[which(is.na(prediction))] <- 0
          #add expected SILAC isotope shifts
          if(multiplicity > 1)
          {
            if(grepl("light",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features_unknown$m.z._shift_light[temp_data$Feature]
            if(grepl("medium",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features_unknown$m.z._shift_medium[temp_data$Feature]
            if(grepl("heavy",colnames(mz_calibration_vals)[c]))SILAC_label_mz_shift <- features_unknown$m.z._shift_heavy[temp_data$Feature]
            prediction <- prediction + SILAC_label_mz_shift
          }
          data.table::set(mz_calibration_vals,as.integer(temp_data$Feature),as.integer(c),value = prediction)
        }
      }
      updatecounter <- updatecounter + 1
      if(updatecounter >= 1)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(c/max))*(1-(c/max))
        td <- lubridate::seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        tcltk::setTkProgressBar(pb, c, label=base::paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    ###fill all remaining NAs
    mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
  }else ###no calibration
  {
    mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
  }
  features_unknown <- cbind(features_unknown,RT_calibration_vals,mz_calibration_vals)

  ###get full features list again
  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  features <- utils::read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  ##add new PMPs
  features <- rbind(features,features_unknown)

  ###Save features list
  utils::write.table(features,base::paste("Temporary_files/",feature_table_file_name,sep=""),row.names = F)

  options(warn=0)
  grDevices::dev.off()
  crap <- gc(F)
}

#' Perform quantification of IceR features
#' @param path_to_features Path to folder where results of align_features() are stored
#' @param path_to_mzXML Path to folder containing mzXML files of samples in case of Orbitrap data.
#' @param path_to_MaxQ_output Path to folder containing MaxQuant outputs (txt folder containing at least allpeptides.txt, evidence.txt, peptides.txt and proteinGroups.txt)
#' @param feature_table_file_name File name which contains align_features() results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param output_file_names_add IceR result name tag. By default IceR_analysis
#' @param RT_calibration Boolean value indicating if corrected RT should be used during peak detection, selection and DICE, By default set to T.
#' @param mz_calibration Boolean value indicating if corrected m/z should be used during peak detection, selection and DICE, By default set to T.
#' @param abundance_estimation_correction Boolean value indicating if resulting peptide abundances should be corrected using MaxQuant results as a reference. By default set to T.
#' @param Quant_pVal_cut Numeric value used as diagnostic cutoff border for visualization of significances of ion accumulation per IceR feature quantification. Furthermore, used as cutoff to filter +1-isotopic IceR features with significant accumulation of ions. By default set to 0.05.
#' @param n_cores Numeric value specifying on how many CPU cores tasks should be parallelized. By default set to 2.
#' @param kde_resolution Numeric value specifying number of grid points per dimension. By default set to 50.
#' @param num_peaks_store Numeric value specifying number of 2D peaks to be stored during peak detection. By default set to 5.
#' @param plot_peak_detection Boolean value indicating if for every feature quantification the determined kernel density estimations and detected peaks should be visualized and stored. By default set to F.
#' @param alignment_variability_score_cutoff Numeric value specifying significance cutoff to distinguish which features show high RT- or m/z-variability of selected peaks between samples. By default set to 0.05. All features showing significant general variability (variability score < alignment_variability_score_cutoff) are excluded.
#' @param alignment_scores_cutoff Numeric value specifying significance cutoff to distinguish which samples show high RT- or m/z-variability of selected peaks for respective IceR feature. By default set to 0.05. All samples showing significant peak variability for respective IceR feature (variability score < alignment_scores_cutoff) are excluded (quantification set to NA).
#' @param mono_iso_alignment_cutoff Numeric value specifying significance cutoff to distinguish which samples show high RT- or m/z-variability of selected peaks for +1-isotopic from corresponding monoisotopic IceR feature. By default set to 0.05. All samples showing significant peak variability between selected +1-isotopic and monoisotopic IceR features (variability score < mono_iso_alignment_cutoff) are excluded (quantification of +1-isotopic feature set to NA).
#' @param calc_peptide_LFQ Boolean value specifying if multiply peptide quantification data for same peptide sequence (multiply charge states, isotope-states) should be aggregated using the MaxLFQ algorithm. By default set to F.
#' @param calc_protein_LFQ Boolean value specifying if protein quantification should be additionally performed by peptide quantification aggregation using the MaxLFQ algorithm. By default set to T.
#' @param MassSpec_mode String being either "Orbitrap" or "TIMSToF" specifying by which type of Mass Spectrometer the data was generated. By default it expects Thermo Orbitrap data.
#' @param use_IM_data Boolean value indicating if ion mobility information should be used during feature quantification in case of TIMS-ToF data. By default set to T.
#' @param path_to_extracted_spectra Path to folder containing extracted spectra files of samples in case of TIMS-ToF data.
#' @import foreach
#' @export
#' @details Performs final steps of the IceR workflow:
#' 1) Estimation of background noise per IceR feature quantification.
#' 2) 2D Kernel density estimation-based peak detection, selection and DICE-based quantification of IceR features.
#' 3) Determination of significances of ion accumulations per IceR feature quantification.
#' 4) Quality control of peak selections.
#' 5) IceR peak selection accuracy estimations per sample.
#' 6) Optional: Imputation of missing IceR feature quantifications using estimated backgroudn noise models per sample.
#' 7) Protein quantification by aggregating available peptide quantifications.
#' @return Outputs are stored in the specified output folder and intermediate results are stored in the sub-directory Temporary_files.
#' Quantification results are stored in tab-delimited text files (.tab):
#' Feature information - Features_DDAiceR_Analysis.tab
#' Peak alignment scores - Features_quantification_alignment_score_DDAiceR_Analysis.tab
#' Feature quantifications - Features_quantification_DDAiceR_Analysis.tab
#' Feature quantifications after imputation - Features_quantification_imputed_DDAiceR_Analysis.tab
#' Numbers of observed ions per feature quantification - Features_quantification_ioncount_DDAiceR_Analysis.tab
#' Alignment scores between monoisotopic and corresponding +1isotopic IceR features - Features_quantification_mono_iso_alignment_score_DDAiceR_Analysis.tab
#' Significance of ion accumulations per feature quantification - Features_quantification_pvals_DDAiceR_Analysis.tab
#' Signal to background intensity ratios per feature quantification - Features_quantification_S2B_DDAiceR_Analysis.tab
#' General variability score of peak selections - Features_quantification_variability_score_DDAiceR_Analysis.tab
#' Protein quantification - Proteins_quantification_LFQ_DDAiceR_Analysis.tab
#' Protein quantification after imputation - Proteins_quantification_LFQ_imputed_DDAiceR_Analysis.tab
#' QC results of background noise estimations are visualized in "Decoy feature quantification parameters.pdf".
#' QC results of peak selections are visualized in "Alignment and quantification scores.pdf".
#' The performance of RT- and m/z-alignments over samples is visualized in "Performance of feature alignment.pdf"
#' Estimation of required peptide abundance correction factors is visualized in "Correct feature abundance estimations Signal_Background_intensity.pdf".
#' Intermediate results of the function are stored in RData files.
requantify_features <- function(path_to_features,path_to_mzXML=NA,path_to_MaxQ_output,feature_table_file_name="Features_aligned_merged_IceR_analysis.txt",output_file_names_add="IceR_analysis",RT_calibration=T,mz_calibration=T,abundance_estimation_correction = T,Quant_pVal_cut=0.05,n_cores=2,kde_resolution = 50,num_peaks_store = 5,plot_peak_detection=F,alignment_variability_score_cutoff=0.05,alignment_scores_cutoff=0.05,mono_iso_alignment_cutoff=0.05,calc_peptide_LFQ=F,calc_protein_LFQ=T,MassSpec_mode=c("Orbitrap","TIMSToF"),use_IM_data = T,path_to_extracted_spectra=NA)
{
  # path_to_features <- "D:\\Publication\\IceR\\test IceR example\\IceR_V1001"
  # path_to_mzXML <- "D:\\Publication\\IceR\\test IceR example\\raw\\mzXML"
  # path_to_MaxQ_output <- "D:\\Publication\\IceR\\test IceR example\\MaxQ"
  # feature_table_file_name="Features_aligned_merged_IceR_analysis_V1001.txt"
  # output_file_names_add="IceR_analysis_V1001"
  # RT_calibration=T
  # mz_calibration=T
  # abundance_estimation_correction = T
  # Quant_pVal_cut=0.05
  # n_cores=4
  # kde_resolution = 50
  # num_peaks_store = 5
  # plot_peak_detection=F
  # alignment_variability_score_cutoff=0.05
  # alignment_scores_cutoff=0.05
  # mono_iso_alignment_cutoff=0.05
  # calc_peptide_LFQ=F
  # calc_protein_LFQ=T
  # MassSpec_mode="Orbitrap"
  # use_IM_data = T
  # path_to_extracted_spectra=NA

  options(warn=-1)
  #suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
  multiply_intensity_count=F
  peak_detection=T

  QC_data <- list() ##here relevant qc data is stored and finally saved as RData which can be used for re-generating plots

  mean_background_intensity <- NA
  sd_background_intensity <- NA
  n_background_intensity <- NA

  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  features <- utils::read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  ###Check peak widths and increase very small peak widths while reduce extreme long outlier peak widths
  stats <- grDevices::boxplot.stats(features$RT_length)
  features$RT_length[which(features$RT_length < stats$stats[2])] <- stats$stats[2] ##lower than 25 % quantile
  features$RT_length[which(features$RT_length > stats$stats[5])] <- stats$stats[5] ##if RT length > upper whisker --> shorten to upper whisker
  features$RT_length[is.na(features$RT_length)] <- stats$stats[3]

  if(MassSpec_mode == "TIMSToF")
  {
    stats <- grDevices::boxplot.stats(features$Inv_K0_length)
    features$Inv_K0_length[which(features$Inv_K0_length < stats$stats[2])] <- stats$stats[2] ##lower than 25 % quantile
    features$Inv_K0_length[which(features$Inv_K0_length > stats$stats[5])] <- stats$stats[5] ##if IM length > upper whisker --> shorten to upper whisker
    features$Inv_K0_length[is.na(features$Inv_K0_length)] <- stats$stats[3]
  }

  ###Add decoy features
  features_decoy <- features[which(!grepl("_pmp|_i",features$Feature_name)),]
  features_decoy$Feature_name <- base::paste(features_decoy$Feature_name,"_d",sep="")

  RT_ranges <- features_decoy$RT_range_max-features_decoy$RT_range_min
  mz_ranges <- features_decoy$m.z_range_max-features_decoy$m.z_range_min
  median_RT_window <- stats::median(features_decoy$RT_range_max-features_decoy$RT_range_min,na.rm=T)
  median_mz_window <- stats::median(features_decoy$m.z_range_max-features_decoy$m.z_range_min,na.rm=T)

  features_decoy$RT <- features_decoy$RT + 5*median_RT_window
  features_decoy$m.z <- features_decoy$m.z + 5*median_mz_window

  features_decoy$m.z_range_min <- features_decoy$m.z-(mz_ranges/2)
  features_decoy$m.z_range_max <- features_decoy$m.z+(mz_ranges/2)

  features_decoy$RT_range_min <- features_decoy$RT-(RT_ranges/2)
  features_decoy$RT_range_max <- features_decoy$RT+(RT_ranges/2)

  ###adjust observed RT and mz accordingly
  Observed_RT <- as.matrix(stringr::str_split(features_decoy$Observed_RT,";",simplify = T))
  class(Observed_RT) <- "numeric"
  Observed_RT <- Observed_RT + 5*median_RT_window
  features_decoy$Observed_RT <- apply(Observed_RT, 1, base::paste, collapse=";")

  Observed_mz <- as.matrix(stringr::str_split(features_decoy$Observed_mz,";",simplify = T))
  class(Observed_mz) <- "numeric"
  Observed_mz <- Observed_mz + 5*median_mz_window
  features_decoy$Observed_mz <- apply(Observed_mz, 1, base::paste, collapse=";")


  if(MassSpec_mode == "TIMSToF")
  {
    IM_ranges <- features_decoy$Inv_K0_range_max-features_decoy$Inv_K0_range_min
    median_IM_window <- stats::median(IM_ranges,na.rm=T)
    features_decoy$Inv_K0 <- features_decoy$Inv_K0 + 5*median_IM_window

    Observed_IM <- as.matrix(stringr::str_split(features_decoy$Observed_IM,";",simplify = T))
    class(Observed_IM) <- "numeric"
    Observed_IM <- Observed_IM + 5*median_IM_window
    features_decoy$Observed_IM <- apply(Observed_IM, 1, base::paste, collapse=";")
  }

  features <- rbind(features,features_decoy)

  if(output_file_names_add != "")output_file_names_add <- base::paste("_",output_file_names_add,sep="")

  ##Step 1 - Summarize ion intensities per aligned MaxQuant feature
  ###Function to perform extraction on multiple threads. Depending on the available ram, the extraction can be run on several threads. However, the task is using much memory so that it is recommended to just run 2 threads in parallel if only 16 gb of ram are available.
  extract_intensities_worker <- function(Sample_IDs,features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection,n_cores,ion_intensity_cutoff=F,mean_background_ion_intensity_model = NA,sd_background_ion_intensity = NA,peak_min_ion_count=NA,kde_resolution=25,num_peaks_store=5,plots=F,MassSpec_mode="Orbitrap",use_IM_data=T)
  {
    #suppressWarnings(suppressMessages(library(doParallel,quietly = T)))

    ####Function to extract intensities in respective extracted .RData for a list of selected features (or all features).
    ####indexing_RT_window defines how large each RT indexing window is to speed up subsetting. 0.5 min was observed to be good
    get_intensities <- function(Sample_ID,path,features_select,indexing_RT_window=0.1,RT_calibration,mz_calibration,peak_detection,ion_intensity_cutoff,mean_background_ion_intensity_model,sd_background_ion_intensity,peak_min_ion_count,kde_resolution,num_peaks_store,plots,MassSpec_mode,use_IM_data)
    {
      #suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
      #suppressWarnings(suppressMessages(library(stringr,quietly = T)))

      if(any(!is.na(mean_background_ion_intensity_model)))
      {
        colnames(mean_background_ion_intensity_model) <- base::gsub("-",".",colnames(mean_background_ion_intensity_model))
        names(sd_background_ion_intensity) <- base::gsub("-",".",names(sd_background_ion_intensity))
      }

      ###if no peak selection should be performed, extract all ions within the expected feature windows
      feature_no_peak_selection <- function(all_ion_data,selected_features,cur_sample,include_isotope_patterns = F,num_peaks_store=0,MassSpec_mode="Orbitrap",use_IM_data=T)
      {
        peak_ion_data_list <- list()

        ###define RT and mz window
        RT_expected <- selected_features$RT + selected_features[,base::paste("RT_calibration.",cur_sample,sep="")]
        mz_expected <- selected_features$m.z + selected_features[,base::paste("mz_calibration.",cur_sample,sep="")]

        RT_window <- c(RT_expected - (selected_features$RT_length/2),
                       RT_expected + (selected_features$RT_length/2))

        mz_window <- c(selected_features$m.z_range_min + selected_features[,base::paste("mz_calibration.",cur_sample,sep="")],
                       selected_features$m.z_range_max + selected_features[,base::paste("mz_calibration.",cur_sample,sep="")])

        if(MassSpec_mode == "TIMSToF")
        {
          im_expected <- selected_features$Inv_K0 + selected_features[,base::paste("IM_calibration.",cur_sample,sep="")]

          IM_window <- c(im_expected - (selected_features$Inv_K0_length/2),
                         im_expected + (selected_features$Inv_K0_length/2))
        }

        if(MassSpec_mode == "Orbitrap")
        {
          ###select relevant ions
          ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window[1] &
                                           all_ion_data$m.z <= mz_window[2] &
                                           all_ion_data$RT >= RT_window[1] &
                                           all_ion_data$RT <= RT_window[2]),]
          ion_data <- stats::na.omit(cbind(ion_data,data.frame(isotope=0)))

          #ion_data$isotope=0

          ##add isotope ions if wanted
          if(include_isotope_patterns == T)
          {
            iso_mult <- (1:3)*1.002054
            cur_charge <- selected_features$Charge
            cur_mz_min = mz_window[1]
            cur_mz_max = mz_window[2]
            for(isos in c(1,2,3))###isotope +1,+2,+3
            {
              mz_range <- c(((cur_mz_min*cur_charge)+(iso_mult[isos]))/cur_charge,((cur_mz_max*cur_charge)+(iso_mult[isos]))/cur_charge)

              selection <- which(all_ion_data$m.z >= mz_range[1] &
                                   all_ion_data$m.z <= mz_range[2] &
                                   all_ion_data$RT >= RT_window[1] &
                                   all_ion_data$RT <= RT_window[2])
              if(length(selection)>0)
              {
                temp <- all_ion_data[selection,]
                temp$isotope=isos
                ion_data <- rbind(ion_data,temp)
              }else
              {
                break
              }
            }

            ion_data <- ion_data[order(ion_data$RT),]
            ###only keep isotope ions which are present in spectra where also the previous isotope ions was detected
            for(iso in c(1,2,3))
            {
              if(length(which(ion_data$isotope == iso))>0 & length(which(ion_data$isotope == iso-1))>0)
              {
                iso_current <- which(ion_data$isotope == iso)

                iso_minus1 <- ion_data[which(ion_data$isotope == iso-1),]

                remove <- which(ion_data$RT[iso_current] %not in% iso_minus1$RT)
                if(length(remove)>0)iso_current <- iso_current[remove]
                if(length(iso_current)>0)ion_data <- ion_data[-iso_current,]
              }else
              {
                remove <- which(ion_data$isotope == iso)
                if(length(remove)>0)ion_data <- ion_data[-remove,]
              }
            }

            ###correct m/z of +1,+2 and +3 isotope ions to its theoretical +0 isotope m/z
            ion_data$m.z <- ((ion_data$m.z*cur_charge)-(ion_data$isotope*1.002054))/cur_charge
          }

          peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data,
                                                                        Peak_info=base::data.frame(RT=RT_expected,
                                                                                             RT_window_lower=RT_window[1],
                                                                                             RT_window_upper=RT_window[2],
                                                                                             mz=mz_expected,
                                                                                             mz_window_lower=mz_window[1],
                                                                                             mz_window_upper=mz_window[2],
                                                                                             density=0,
                                                                                             known_peak=0,
                                                                                             Peak=num_peaks_store+1,
                                                                                             ion_count=nrow(ion_data)))
        }
        if(MassSpec_mode == "TIMSToF")
        {
          ###select relevant ions
          if(use_IM_data == T)
          {
            ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window[1] &
                                             all_ion_data$m.z <= mz_window[2] &
                                             all_ion_data$RT >= RT_window[1] &
                                             all_ion_data$RT <= RT_window[2] &
                                             all_ion_data$`1/K0` >= IM_window[1] &
                                             all_ion_data$`1/K0` <= IM_window[2]),]
          }else
          {
            ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window[1] &
                                             all_ion_data$m.z <= mz_window[2] &
                                             all_ion_data$RT >= RT_window[1] &
                                             all_ion_data$RT <= RT_window[2]),]
          }
          ion_data <- stats::na.omit(cbind(ion_data,data.frame(isotope=0)))

          #ion_data$isotope=0

          peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data,
                                                                        Peak_info=base::data.frame(RT=RT_expected,
                                                                                             RT_window_lower=RT_window[1],
                                                                                             RT_window_upper=RT_window[2],
                                                                                             mz=mz_expected,
                                                                                             mz_window_lower=mz_window[1],
                                                                                             mz_window_upper=mz_window[2],
                                                                                             density=0,
                                                                                             known_peak=0,
                                                                                             Peak=num_peaks_store+1,
                                                                                             ion_count=nrow(ion_data)))
          #IM=im_expected,
          #IM_window_lower=IM_window[1],
          #IM_window_upper=IM_window[2]))
        }

        names(peak_ion_data_list) <- "Standard"
        return(peak_ion_data_list)

      }

      ###if peak selection should be performed
      feature_2D_peak_selection <- function(all_ion_data,selected_features,cur_sample,known_RT,delta_mz,delta_rt,RT_window_expand_factor=5,IM_window_expand_factor=5,mz_window_expand_factor=4,include_isotope_patterns = F,n_raster_dens_matrix=50,local_maxima_k=3,max_delta_RT=2,max_delta_mz=0.005,peak_min_ion_count=5,RT_bw=0.5,mz_bw=0.002,num_peaks_store=5,plot=F,auto_adjust_kde_resolution=T,MassSpec_mode="Orbitrap",delta_im=0.001,close_peak_merging=F,use_IM_data=T)
      {
        #suppressWarnings(suppressMessages(library(MASS,quietly = T)))
        #suppressWarnings(suppressMessages(library(raster,quietly = T)))
        #suppressWarnings(suppressMessages(library(ggtern,quietly = T)))
        peak_ion_data_list <- list()

        graph <- NA ##here we store graphical output if wanted

        ###define RT and mz window
        RT_correction <- as.numeric(selected_features[base::paste("RT_calibration.",cur_sample,sep="")])
        mz_correction <- as.numeric(selected_features[base::paste("mz_calibration.",cur_sample,sep="")])

        RT_expected <- as.numeric(selected_features$RT+RT_correction)
        mz_expected <- as.numeric(selected_features$m.z+mz_correction)


        if(delta_rt > selected_features$RT_length/2) ###use standard window as long as the peak width is < standard RT window
        {
          RT_window <- c(selected_features$RT_range_min + RT_correction,
                         selected_features$RT_range_max + RT_correction)
        }else
        {
          RT_window <- c(selected_features$RT + RT_correction - (selected_features$RT_length/2),
                         selected_features$RT + RT_correction + (selected_features$RT_length/2))
        }

        mz_window <- c(selected_features$m.z_range_min + mz_correction,
                       selected_features$m.z_range_max + mz_correction)

        RT_window_expanded <- c(RT_window[1]-((delta_rt)*RT_window_expand_factor),
                                RT_window[2]+((delta_rt)*RT_window_expand_factor))

        if(RT_window_expanded[2]-RT_window_expanded[1] > 10) ###maximum 10 min deviation
        {
          RT_window_expanded[1] <- RT_window[1]-5
          RT_window_expanded[2] <- RT_window[2]+5
        }

        mz_window_expanded <- c(mz_window[1]-(delta_mz*mz_window_expand_factor),
                                mz_window[2]+(delta_mz*mz_window_expand_factor))

        if(MassSpec_mode == "TIMSToF")
        {
          IM_correction <- as.numeric(selected_features[base::paste("IM_calibration.",cur_sample,sep="")])
          IM_expected <- as.numeric(selected_features$Inv_K0+IM_correction)

          if(delta_im > selected_features$Inv_K0_length/2) ###use standard window as long as the peak width is < standard IM window
          {
            IM_window <- c(selected_features$Inv_K0_range_min + IM_correction,
                           selected_features$Inv_K0_range_max + IM_correction)
          }else
          {
            IM_window <- c(selected_features$Inv_K0 + IM_correction - (selected_features$Inv_K0_length/2),
                           selected_features$Inv_K0 + IM_correction + (selected_features$Inv_K0_length/2))
          }
          IM_window_expanded <- c(IM_window[1]-(delta_im*IM_window_expand_factor),
                                  IM_window[2]+(delta_im*IM_window_expand_factor))
        }


        ###filter for ions in the expanded window
        if(MassSpec_mode == "Orbitrap")
        {
          ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window_expanded[1] &
                                           all_ion_data$m.z <= mz_window_expanded[2] &
                                           all_ion_data$RT >= RT_window_expanded[1] &
                                           all_ion_data$RT <= RT_window_expanded[2]),]
        }else
        {
          if(use_IM_data == T)
          {
            ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window_expanded[1] &
                                             all_ion_data$m.z <= mz_window_expanded[2] &
                                             all_ion_data$RT >= RT_window_expanded[1] &
                                             all_ion_data$RT <= RT_window_expanded[2] &
                                             all_ion_data$`1/K0` >= IM_window[1] &
                                             all_ion_data$`1/K0` <= IM_window[2]),]
          }else
          {
            ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window_expanded[1] &
                                             all_ion_data$m.z <= mz_window_expanded[2] &
                                             all_ion_data$RT >= RT_window_expanded[1] &
                                             all_ion_data$RT <= RT_window_expanded[2]),]
          }
        }

        ion_data <- stats::na.omit(cbind(ion_data,data.frame(isotope=0)))
        #ion_data$isotope=0

        ##add isotope ions if wanted
        if(include_isotope_patterns == T & nrow(ion_data)>0)
        {
          iso_mult <- (1:3)*1.002054

          cur_charge <- selected_features$Charge
          cur_mz_min = mz_window_expanded[1]
          cur_mz_max = mz_window_expanded[2]
          for(isos in c(1,2,3))###isotope +1,+2,+3
          {
            mz_range <- c(((cur_mz_min*cur_charge)+(iso_mult[isos]))/cur_charge,((cur_mz_max*cur_charge)+(iso_mult[isos]))/cur_charge)

            selection <- which(all_ion_data$m.z >= mz_range[1] &
                                 all_ion_data$m.z <= mz_range[2] &
                                 all_ion_data$RT >= RT_window_expanded[1] &
                                 all_ion_data$RT <= RT_window_expanded[2])
            if(length(selection)>0)
            {
              temp <- all_ion_data[selection,]
              temp$isotope=isos
              ion_data <- rbind(ion_data,temp)
            }else
            {
              break
            }
          }

          ion_data <- ion_data[order(ion_data$RT),]
          ###only keep isotope ions which are present in spectra where also the previous isotope ions was detected
          for(iso in c(1,2,3))
          {
            if(length(which(ion_data$isotope == iso))>0 & length(which(ion_data$isotope == iso-1))>0)
            {
              iso_current <- which(ion_data$isotope == iso)

              iso_minus1 <- ion_data[which(ion_data$isotope == iso-1),]

              remove <- which(ion_data$RT[iso_current] %not in% iso_minus1$RT)
              if(length(remove)>0)iso_current <- iso_current[remove]
              if(length(iso_current)>0)ion_data <- ion_data[-iso_current,]
            }else
            {
              remove <- which(ion_data$isotope == iso)
              if(length(remove)>0)ion_data <- ion_data[-remove,]
            }
          }

          ###correct m/z of +1,+2 and +3 isotope ions to its theoretical +0 isotope m/z
          ion_data$m.z <- ((ion_data$m.z*cur_charge)-(ion_data$isotope*1.002054))/cur_charge
        }

        if(nrow(ion_data) >= peak_min_ion_count)
        {
          ###determine 2D density matrix
          #consider different resolutions required for peaks with expected small RT lengths
          if(auto_adjust_kde_resolution==T)
          {
            required_res <- ceiling(((RT_window_expanded[2]-RT_window_expanded[1])/((selected_features$RT_length/1.9)*2))+1)
            if(required_res > n_raster_dens_matrix)n_raster_dens_matrix <- required_res
          }

          f2 <- MASS::kde2d(ion_data$RT, ion_data$m.z, n = n_raster_dens_matrix,
                      h = c(RT_bw, mz_bw),lims = c(RT_window_expanded[1],RT_window_expanded[2],mz_window_expanded[1],mz_window_expanded[2]))

          #area_per_cell <- ((max(f2$x)-min(f2$x)))*((max(f2$y)-min(f2$y)))/2

          ###determine which density (=z) will be selected to be at least exceeded --> upper whisker
          outlier_densities <- grDevices::boxplot.stats(f2$z)$stats[5]

          ## Convert it to a raster object
          r <- raster::raster(f2$z)
          raster::extent(r) <- raster::extent(c(0, length(f2$x), 0, length(f2$y)) + 0.5)

          ## Find the maximum value within the k-cell neighborhood of each cell
          f <- function(X) max(X, na.rm=TRUE)
          ww <- matrix(1, nrow=local_maxima_k, ncol=local_maxima_k) ## Weight matrix for cells in moving window
          localmax <- raster::focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)

          ## Get x-y coordinates of those cells that are local maxima
          r2 <- r==localmax
          maxXY <- raster::xyFromCell(r2, raster::Which(r2==1, cells=TRUE))

          ## Remove maxima which are below the minimal outlier density
          maxima <- base::data.frame(RT=f2$x[n_raster_dens_matrix-maxXY[,2]+1],mz=f2$y[maxXY[,1]],density=f2$z[((maxXY[,1]-1)*n_raster_dens_matrix)+(n_raster_dens_matrix-maxXY[,2]+1)])
          maxima <- maxima[which(maxima$density > outlier_densities),]
          ###next filter for maxima with at least peak_min_ion_count ions
          maxima$estimated_count <- (maxima$density/sum(as.matrix(maxima$density)))*nrow(ion_data)

          RT_cut <- selected_features$RT_length/2#ifelse(!is.na(known_RT),selected_features$RT_length/2,delta_rt)
          selection <- which(maxima$estimated_count >= peak_min_ion_count | abs(maxima$RT-RT_expected)<=RT_cut & abs(maxima$mz-mz_expected)<=delta_mz)
          ###if no maxima are left after filtering still take topN closest peak further
          if(length(selection)==0 & nrow(maxima) > 0)
          {
            ###no maxima would be left see keep up to N peaks
            distances <- base::data.frame(dist_rt=maxima$RT-RT_expected,
                                    dist_mz=maxima$mz-mz_expected,
                                    dist_total_scaled=sqrt(((maxima$RT-RT_expected))^2 + ((maxima$mz-mz_expected)*500)^2))


            ###at maximum the top N closest maxima
            ordering <- order(abs(distances$dist_total_scaled))
            maxima <- maxima[ordering,]
            if(nrow(maxima)>num_peaks_store) maxima <- maxima[1:num_peaks_store,]
          }else if(nrow(maxima) > 0)
          {
            maxima <- maxima[selection,]
          }

          ###next merge peaks which are very close together
          if(nrow(maxima)>1 & close_peak_merging == T)
          {
            close_info <- base::as.data.frame(matrix(ncol=2,nrow=nrow(maxima),0))

            for(cl in 1:nrow(maxima))
            {
              temp <- sqrt((maxima$RT[cl]-maxima$RT)^2+((maxima$mz[cl]-maxima$mz)*500)^2)
              sel <- which(temp == min(temp[-cl]))[1]
              data.table::set(close_info,as.integer(cl),as.integer(1:2),value = list(temp[sel],sel))
            }

            dist_cut <- sqrt((delta_mz*500)^2+(selected_features$RT_length/4)^2) ###merge peaks with d_mz < 0.001 and d_RT < peak_width/2
            if(any(close_info$V1 < dist_cut))
            {
              sel <- which(close_info$V1 < dist_cut)

              maxima_add <- base::as.data.frame(matrix(nrow=length(sel),ncol=4,0))
              colnames(maxima_add) <- colnames(maxima)
              ###find mean RT and mz of maxima which should be merged wheigted by density
              for(cl in 1:length(sel))
              {
                RT_add <- stats::weighted.mean(c(maxima$RT[sel[cl]],maxima$RT[close_info[sel[cl],2]]),
                                        c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                mz_add <- stats::weighted.mean(c(maxima$mz[sel[cl]],maxima$mz[close_info[sel[cl],2]]),
                                        c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                density_add <- stats::weighted.mean(c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]),
                                             c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                count_add <- stats::weighted.mean(c(maxima$estimated_count[sel[cl]],maxima$estimated_count[close_info[sel[cl],2]]),
                                           c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                data.table::set(maxima_add,as.integer(cl),as.integer(1:4),value=list(RT_add,mz_add,density_add,count_add))
              }
              maxima_add <- unique(maxima_add)
              ###remove maxima which are merged
              remove <- unique(append(sel,close_info[sel,2]))
              maxima <- maxima[-remove,]
              ###add new merged peaks
              maxima <- rbind(maxima,maxima_add)
            }
          }

          ###now check which density maxima is closest to the expected feature
          if(nrow(maxima)>0)
          {
            distances <- base::data.frame(dist_rt=maxima$RT-RT_expected,
                                    dist_mz=maxima$mz-mz_expected,
                                    dist_total_scaled=sqrt(((maxima$RT-RT_expected))^2 + ((maxima$mz-mz_expected)*500)^2))
            ###at maximum the top N closest maxima
            ordering <- order(abs(distances$dist_total_scaled))
            maxima_select <- maxima[ordering,]#which(abs(distances$dist_rt[ordering]) <= max_delta_RT & abs(distances$dist_mz[ordering]) <= max_delta_mz),]

            if(nrow(maxima_select)>num_peaks_store) maxima_select <- maxima_select[1:num_peaks_store,]

            if(nrow(maxima_select)>0)
            {
              ###label which peak was closest to expected feature RT and mz if known
              maxima_select$known_peak <- 0
              if(!is.na(known_RT))
              {
                ###RT is known
                ###check if any peak was detected close to known RT (max deviation by peak width/2)
                peak_close_rt <- which(abs(maxima_select$RT - RT_expected) <= selected_features$RT_length/1.9 & abs(maxima_select$mz - mz_expected) <= delta_mz/1.9)
                ###and if the closest peak actually also shows a good number of ions
                ###otherwise no peak will be selected here but selection will be done later
                if(length(peak_close_rt)>0 & maxima_select$estimated_count[1] >= 2* peak_min_ion_count)
                {
                  maxima_select$known_peak[peak_close_rt[1]] <- 1
                  ###reorder maxima selected such that known peak maximum is at index 1
                  maxima_select <- maxima_select[c(peak_close_rt[1],c(1:nrow(maxima_select))[-peak_close_rt[1]]),]
                }else ###peak RT is known but none of the detected peaks shows RT close to the known RT
                {
                  ###in this case none of the maxima are set to be the known peak and later we will try to decide which one is the correct one
                }
              }
              ###if any maximum was detected, return table of ions as a list containing up to 3 peak area ion lists from closest to further away peaks.
              RT_window_width <- selected_features$RT_length
              mz_window_width <- delta_mz

              peak_counter=0

              for(i in 1:nrow(maxima_select))
              {
                cur_mz_window <- c(maxima_select$mz[i]-(mz_window_width),
                                   maxima_select$mz[i]+(mz_window_width))
                cur_RT_window <- c(maxima_select$RT[i]-(RT_window_width/2),
                                   maxima_select$RT[i]+(RT_window_width/2))

                selection <- which(ion_data$m.z >= cur_mz_window[1] &
                                     ion_data$m.z <= cur_mz_window[2] &
                                     ion_data$RT >= cur_RT_window[1] &
                                     ion_data$RT <= cur_RT_window[2])

                if(length(selection)<peak_min_ion_count & maxima_select$known_peak[i] == 1)maxima_select$known_peak[i] <- 0 ###if this peak should be the correct peak but number of ions are very low then set to unknown correct peak

                peak_counter <- peak_counter + 1

                peak_ion_data_list[[peak_counter]] <- list(ion_data=ion_data[selection,],
                                                           Peak_info=base::data.frame(RT=maxima_select$RT[i],
                                                                                RT_window_lower=cur_RT_window[1],
                                                                                RT_window_upper=cur_RT_window[2],
                                                                                mz=maxima_select$mz[i],
                                                                                mz_window_lower=cur_mz_window[1],
                                                                                mz_window_upper=cur_mz_window[2],
                                                                                density=maxima_select$density[i],
                                                                                known_peak=maxima_select$known_peak[i],
                                                                                Peak=peak_counter,
                                                                                ion_count=length(selection)))

              }
              if(peak_counter>0)names(peak_ion_data_list) <- as.character(1:peak_counter)
              cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))
              ###also add quantification if standard windows are used
              selection <- which(ion_data$m.z >= mz_window[1] &
                                   ion_data$m.z <= mz_window[2] &
                                   ion_data$RT >= cur_RT_window[1] &
                                   ion_data$RT <= cur_RT_window[2])
              peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data[selection,],
                                                                            Peak_info=base::data.frame(RT=RT_expected,
                                                                                                 RT_window_lower=cur_RT_window[1],
                                                                                                 RT_window_upper=cur_RT_window[2],
                                                                                                 mz=mz_expected,
                                                                                                 mz_window_lower=mz_window[1],
                                                                                                 mz_window_upper=mz_window[2],
                                                                                                 density=0,
                                                                                                 known_peak=0,
                                                                                                 Peak=num_peaks_store+1,
                                                                                                 ion_count=length(selection)))
              ###plot should be directly generated
              if(plot==T)
              {
                graphics::image(f2,main=base::paste(cur_sample,"-",selected_features$Feature_name),xlab="RT",ylab="m/z",xlim=RT_window_expanded,ylim=mz_window_expanded)
                #graphics::points(maxima_select$RT,maxima_select$mz)
                graphics::text(maxima_select$RT,maxima_select$mz,1:nrow(maxima_select),col="black")
                ###expected window
                graphics::rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2)
                graphics::points(RT_expected,mz_expected,pch=4)
                ###adjusted window

                graphics::rect(maxima_select$RT-RT_window_width/2,maxima_select$mz-mz_window_width,maxima_select$RT+RT_window_width/2,maxima_select$mz+mz_window_width,lty=2,border=ifelse(maxima_select$known_peak==1,"green","darkgrey"))

                ###add label indicating intensity within the window
                for(i in 1:length(peak_ion_data_list))
                {
                  sum_intensity <- sum(sum(peak_ion_data_list[[i]]$ion_data$Intensity))
                  if(sum_intensity > 0)sum_intensity <- base::log2(sum_intensity)
                  if(peak_ion_data_list[[i]]$Peak_info$Peak != num_peaks_store+1)
                  {
                    graphics::text(peak_ion_data_list[[i]]$Peak_info$RT,peak_ion_data_list[[i]]$Peak_info$mz-(1.1*delta_mz),round(sum_intensity,2),col=ifelse(peak_ion_data_list[[i]]$Peak_info$known_peak==1,"green","darkgrey"))
                  }else
                  {
                    graphics::text(peak_ion_data_list[[i]]$Peak_info$RT,peak_ion_data_list[[i]]$Peak_info$mz+(1.1*delta_mz),round(sum_intensity,2),col="black")
                  }

                }
              }else ###store all relevant data for performing plotting later
              {
                peak_intensities <- base::as.data.frame(matrix(ncol=1,nrow=length(peak_ion_data_list)))
                peak_intensities[,1] <- as.numeric(peak_intensities[,1])

                for(i in 1:length(peak_ion_data_list))
                {
                  sum_intensity <- sum(sum(peak_ion_data_list[[i]]$ion_data$Intensity))
                  if(sum_intensity > 0)sum_intensity <- base::log2(sum_intensity)
                  data.table::set(peak_intensities,as.integer(i),as.integer(1),sum_intensity)
                }
                graph <- list(kdemap=f2,
                              maxima_select=maxima_select,
                              peak_intensities=peak_intensities)
              }
            }else ###no maxima with a minmal density found
            {
              ###in this case simply extract ions which are within the expected window
              RT_window_width <- selected_features$RT_length
              cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))

              selection <- which(ion_data$m.z >= mz_window[1] &
                                   ion_data$m.z <= mz_window[2] &
                                   ion_data$RT >= cur_RT_window[1] &
                                   ion_data$RT <= cur_RT_window[2])
              peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data[selection,],
                                                                            Peak_info=base::data.frame(RT=RT_expected,
                                                                                                 RT_window_lower=cur_RT_window[1],
                                                                                                 RT_window_upper=cur_RT_window[2],
                                                                                                 mz=mz_expected,
                                                                                                 mz_window_lower=mz_window[1],
                                                                                                 mz_window_upper=mz_window[2],
                                                                                                 density=0,
                                                                                                 known_peak=0,
                                                                                                 Peak=4,
                                                                                                 ion_count=length(selection)))
            }

          }else ###no maxima found
          {
            ###in this case simply extract ions which are within the expected window
            RT_window_width <- selected_features$RT_length
            cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))

            selection <- which(ion_data$m.z >= mz_window[1] &
                                 ion_data$m.z <= mz_window[2] &
                                 ion_data$RT >= cur_RT_window[1] &
                                 ion_data$RT <= cur_RT_window[2])
            peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data[selection,],
                                                                          Peak_info=base::data.frame(RT=RT_expected,
                                                                                               RT_window_lower=cur_RT_window[1],
                                                                                               RT_window_upper=cur_RT_window[2],
                                                                                               mz=mz_expected,
                                                                                               mz_window_lower=mz_window[1],
                                                                                               mz_window_upper=mz_window[2],
                                                                                               density=0,
                                                                                               known_peak=0,
                                                                                               Peak=4,
                                                                                               ion_count=length(selection)))
          }
        }else
        {
          ###not enough ions even in expanded window
          ###in this case simply extract ions which are within the expected window
          RT_window_width <- selected_features$RT_length
          cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))

          selection <- which(ion_data$m.z >= mz_window[1] &
                               ion_data$m.z <= mz_window[2] &
                               ion_data$RT >= cur_RT_window[1] &
                               ion_data$RT <= cur_RT_window[2])
          peak_ion_data_list[[as.character(num_peaks_store+1)]] <- list(ion_data=ion_data[selection,],
                                                                        Peak_info=base::data.frame(RT=RT_expected,
                                                                                             RT_window_lower=cur_RT_window[1],
                                                                                             RT_window_upper=cur_RT_window[2],
                                                                                             mz=mz_expected,
                                                                                             mz_window_lower=mz_window[1],
                                                                                             mz_window_upper=mz_window[2],
                                                                                             density=0,
                                                                                             known_peak=0,
                                                                                             Peak=4,
                                                                                             ion_count=length(selection)))
        }

        names(peak_ion_data_list)[which(names(peak_ion_data_list) != as.character(num_peaks_store+1))] <- base::paste("Peak_",names(peak_ion_data_list)[which(names(peak_ion_data_list) != as.character(num_peaks_store+1))],sep="")
        names(peak_ion_data_list)[which(names(peak_ion_data_list) == as.character(num_peaks_store+1))] <- "Standard"
        peak_ion_data_list$graph <- graph
        return(peak_ion_data_list)
      }

      `%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

      t.test2 <- function(...) ###T-test
      {
        obj<-try(stats::t.test(...), silent=TRUE)
        if (methods::is(obj, "try-error")) return(NA) else return(obj$p.value)
      }

      setwd(path)


      max <- 1
      print("Load spectra data")
      pb <- tcltk::tkProgressBar(title = base::paste("Load spectra data:",Sample_ID),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      table_store <- NULL
      if(grepl("_Channel_light|_Channel_medium|_Channel_heavy",Sample_ID))#SILAC mode --- all channels in one raw file
      {
        load(base::paste(base::gsub("_Channel_light|_Channel_medium|_Channel_heavy","",Sample_ID),"_all_ions.RData",sep="")) ###ions with RT and intensity in variable dat
      }else
      {
        load(base::paste(Sample_ID,"_all_ions.RData",sep="")) ###ions with RT and intensity in variable dat
      }

      tcltk::setTkProgressBar(pb, 1, label=base::paste( round(1/max*100, 0)," % done (",1,"/",max,")",sep = ""))
      close(pb)

      if(MassSpec_mode == "TIMSToF")
      {
        dat <- table_store
        rm(table_store)
        gc()

        dat$RT <- dat$RT/60
        dat <- dat[,c(3,1,4,2)]
        colnames(dat) <- c("m.z","RT","Intensity","1/K0")
      }

      ###Indexing dat by RT windows to improve speed for subsetting
      num_windows <- ceiling(ceiling(max(dat$RT))*(1/indexing_RT_window))

      ###if in TIMSTF mode, further index by IM to further improve speed
      if(MassSpec_mode == "TIMSToF")
      {
        indexing_IM_window <- 0.025
        num_IM_windows <- ceiling(ceiling(max(dat$`1/K0`)-min(dat$`1/K0`))*(1/indexing_IM_window))
        num_windows_RT_only <- num_windows
        min_im <- min(dat$`1/K0`)
        num_windows <- num_windows * num_IM_windows
      }

      if(MassSpec_mode == "Orbitrap")
      {
        Indices <- base::as.data.frame(matrix(ncol=4,nrow=num_windows))
        colnames(Indices) <- c("RT_start","RT_end","Row_start","Row_end")
        Indices$RT_start <- as.numeric(Indices$RT_start)
        Indices$RT_end <- as.numeric(Indices$RT_end)
        Indices$Row_start <- as.numeric(Indices$Row_start)
        Indices$Row_end <- as.numeric(Indices$Row_end)

        print("Indexing intensities")
        max <- nrow(Indices)
        pb <- tcltk::tkProgressBar(title = base::paste("Indexing intensities:",Sample_ID),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        start <- 1
        end <- max

        cur_index <- 1

        ###increment per indexing
        incr <- 10000

        for(i in 1:nrow(Indices))#
        {
          start <- (i-1)*indexing_RT_window
          end <- i*indexing_RT_window

          indx <- cur_index
          indx_prev <- cur_index
          while(T)
          {
            indx <- indx + incr

            if(dat$RT[indx] >= end | indx > nrow(dat))
            {
              break
            }else
            {
              indx_prev <- indx
            }
          }
          if(indx > nrow(dat))indx <- nrow(dat)

          #inds <- which(dat$RT > start & dat$RT < end)
          if(indx > indx_prev)
          {
            temp <- dat[indx_prev:indx,]

            if(length(which(temp$RT < end)) > 0)
            {
              inds <- cur_index:(indx_prev-1+max(which(temp$RT < end)))

              incr <- max(inds) + 1 - cur_index

              cur_index <- max(inds) + 1
              if(length(inds)>0)
              {
                data.table::set(Indices,i = as.integer(i),j=as.integer(1:4),value=as.list(c(start,end,min(inds),max(inds))))
              }
            }
          }


          updatecounter <- updatecounter + 1
          if(updatecounter >= 1)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }

        }
        close(pb)

        if(any(rowSums(is.na(Indices)) > 0))
        {
          Indices <- Indices[-which(rowSums(is.na(Indices)) > 0),]
        }

        if(min(features_select$RT_range_min) < Indices$RT_start[1])Indices$RT_start[1] <- min(features_select$RT_range_min)
        if(max(features_select$RT_range_max) > Indices$RT_end[nrow(Indices)])Indices$RT_end[nrow(Indices)] <- max(features_select$RT_range_max)

        ###Fragment data into indexed RT windows to further improve subsetting speed
        data_frags <- list()

        for(i in 1:nrow(Indices))
        {
          data_frags[[i]] <- dat[Indices$Row_start[i]:Indices$Row_end[i],]
        }

        rm(dat)
        crap <- gc(F)

      }else
      {
        Indices <- base::as.data.frame(matrix(ncol=4,nrow=num_windows))
        colnames(Indices) <- c("RT_start","RT_end","IM_start","IM_end")
        Indices$RT_start <- sort(rep(((0:(num_windows_RT_only-1))*indexing_RT_window),num_IM_windows))
        Indices$RT_end <- sort(rep(((1:num_windows_RT_only)*indexing_RT_window),num_IM_windows))
        Indices$IM_start <- min_im + ((0:(num_IM_windows-1))*indexing_IM_window)
        Indices$IM_end <- min_im + ((1:num_IM_windows)*indexing_IM_window)

        Indices_list <- list()
        Indices_list_count <- 0

        print("Indexing intensities")
        max <- num_windows_RT_only
        pb <- tcltk::tkProgressBar(title = base::paste("Indexing intensities:",Sample_ID),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        start <- 1
        end <- max

        cur_index <- 1

        ###increment per indexing
        incr <- 10000

        for(i in 1:num_windows_RT_only)
        {
          start <- (i-1)*indexing_RT_window
          end <- i*indexing_RT_window

          indx <- cur_index
          indx_prev <- cur_index
          while(T)
          {
            indx <- indx + incr

            if(dat$RT[indx] >= end | indx > nrow(dat))
            {
              break
            }else
            {
              indx_prev <- indx
            }
          }
          if(indx > nrow(dat))indx <- nrow(dat)

          if(indx > indx_prev)
          {
            temp <- dat[indx_prev:indx,]

            if(length(which(temp$RT < end)) > 0)
            {
              inds <- cur_index:(indx_prev-1+max(which(temp$RT < end)))
              incr <- max(inds) + 1 - cur_index
              temp <- dat[inds,]

              ##now further subset based on IM
              for(j in 1:num_IM_windows)
              {
                start_im <- min_im + ((j-1)*indexing_IM_window)
                end_im <- min_im + (j*indexing_IM_window)

                indx_temp <- inds[which(temp$`1/K0` >= start_im & temp$`1/K0` <= end_im)]
                Indices_list_count <- Indices_list_count + 1
                if(length(indx_temp)>0)
                {
                  Indices_list[[Indices_list_count]] <- indx_temp
                }else
                {
                  Indices_list[[Indices_list_count]] <- NA
                }
              }

              cur_index <- max(inds) + 1
            }
          }

          if(Indices_list_count != (i*num_IM_windows))
          {
            missing <- (i*num_IM_windows) - Indices_list_count

            for(j in 1:missing)
            {
              Indices_list_count <- Indices_list_count + 1
              Indices_list[[Indices_list_count]] <- NA
            }
          }

          updatecounter <- updatecounter + 1
          if(updatecounter >= 1)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }

        }
        close(pb)

        if(any(is.na(Indices_list)))
        {
          remove <- which(is.na(Indices_list))
          Indices <- Indices[-remove,]
          na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
          Indices_list <- na.omit.list(Indices_list)
        }

        ###Fragment data into indexed RT windows to further improve subsetting speed
        data_frags <- list()

        print("Fragment data")
        max <- nrow(Indices)
        pb <- tcltk::tkProgressBar(title = base::paste("Fragment data:",Sample_ID),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        for(i in 1:nrow(Indices))
        {
          data_frags[[i]] <- dat[Indices_list[[i]],]
          updatecounter <- updatecounter + 1
          if(updatecounter >= 100)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }

        }
        close(pb)
      }

      rm(dat)
      crap <- gc(F)


      ###Now extract intensities,
      ###row1 = Intensity summed,
      ###row2 = num ions,
      ###row3 = mean intensity,
      ###row4 = sd intensity,
      ###row5 = detected optimum mz window minimum
      ###row6 = detected optimum mz window maximum
      ###row7 = detected optimum RT window minimum
      ###row8 = detected optimum RT window maximum
      ###row9 = number of detected peaks

      if(peak_detection == F)
      {
        graph_peaks <- list()
        peaks_quant <- list()

        #if(MassSpec_mode == "Orbitrap")
        {
          Intensities <- base::as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
          colnames(Intensities) <- features_select$Feature_name

          ###if cut off pvalue is defined then store background quantifications separately
          if(ion_intensity_cutoff == T)
          {
            Intensities_signal_background <- base::as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
            colnames(Intensities_signal_background ) <- features_select$Feature_name
          }
        }
        # else
        # {
        #   Intensities <- base::as.data.frame(matrix(nrow=12,ncol=nrow(features_select),0))
        #   colnames(Intensities) <- features_select$Feature_name
        #
        #   ###if cut off pvalue is defined then store background quantifications separately
        #   if(ion_intensity_cutoff == T)
        #   {
        #     Intensities_signal_background <- base::as.data.frame(matrix(nrow=12,ncol=nrow(features_select),0))
        #     colnames(Intensities_signal_background ) <- features_select$Feature_name
        #   }
        # }

        ###store data for calculating scoring per sample and feature
        delta_T1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_T1) <- features_select$Feature_name
        delta_T2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_T2) <- features_select$Feature_name
        delta_M1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_M1) <- features_select$Feature_name
        delta_M2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_M2) <- features_select$Feature_name
      }else
      {
        graph_peaks <- list()
        peaks_quant <- list()
        #if(MassSpec_mode == "Orbitrap")
        {
          for(p in 1:(num_peaks_store+1)) ###1-num_peaks_store stores the quantification for top closest peaks while last stores total quantification of the window
          {
            Intensities <- base::as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
            colnames(Intensities) <- features_select$Feature_name

            ###if cut off pvalue is defined then store background quantifications separately
            if(ion_intensity_cutoff == T)
            {
              Intensities_signal_background <- base::as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
              colnames(Intensities_signal_background ) <- features_select$Feature_name
            }
            ###store data for calculating scoring per sample and feature
            delta_T1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
            colnames(delta_T1) <- features_select$Feature_name
            delta_T2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
            colnames(delta_T2) <- features_select$Feature_name
            delta_M1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
            colnames(delta_M1) <- features_select$Feature_name
            delta_M2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
            colnames(delta_M2) <- features_select$Feature_name

            peaks_quant[[p]] <- list(Intensities=Intensities,
                                     Intensities_signal_background=Intensities_signal_background,
                                     delta_T1=delta_T1,
                                     delta_T2=delta_T2,
                                     delta_M1=delta_M1,
                                     delta_M2=delta_M2)
          }
          names(peaks_quant)[1:num_peaks_store] <- base::paste("Peak_",1:num_peaks_store,sep="")
          names(peaks_quant)[num_peaks_store+1] <- "Standard"
        }
        # if(MassSpec_mode == "TIMSToF")
        # {
        #
        #   for(p in 1:(num_peaks_store+1)) ###1-num_peaks_store stores the quantification for top closest peaks while last stores total quantification of the window
        #   {
        #     Intensities <- base::as.data.frame(matrix(nrow=12,ncol=nrow(features_select),0))
        #     colnames(Intensities) <- features_select$Feature_name
        #
        #     ###if cut off pvalue is defined then store background quantifications separately
        #     if(ion_intensity_cutoff == T)
        #     {
        #       Intensities_signal_background <- base::as.data.frame(matrix(nrow=12,ncol=nrow(features_select),0))
        #       colnames(Intensities_signal_background ) <- features_select$Feature_name
        #     }
        #     ###store data for calculating scoring per sample and feature
        #     delta_T1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        #     colnames(delta_T1) <- features_select$Feature_name
        #     delta_T2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        #     colnames(delta_T2) <- features_select$Feature_name
        #     delta_M1 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        #     colnames(delta_M1) <- features_select$Feature_name
        #     delta_M2 <- base::as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        #     colnames(delta_M2) <- features_select$Feature_name
        #
        #     peaks_quant[[p]] <- list(Intensities=Intensities,
        #                              Intensities_signal_background=Intensities_signal_background,
        #                              delta_T1=delta_T1,
        #                              delta_T2=delta_T2,
        #                              delta_M1=delta_M1,
        #                              delta_M2=delta_M2)
        #   }
        #   names(peaks_quant)[1:num_peaks_store] <- base::paste("Peak_",1:num_peaks_store,sep="")
        #   names(peaks_quant)[num_peaks_store+1] <- "Standard"
        # }
      }

      print("Extracting intensities")
      max <- nrow(features_select)
      pb <- tcltk::tkProgressBar(title = base::paste("Extracting intensities:",Sample_ID),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      end <- max

      Sample_ID_save <- Sample_ID
      Sample_ID <- base::gsub("-",".",Sample_ID)

      iso_mult <- (1:3)*1.002054

      peak_width_stats <- grDevices::boxplot.stats(features_select$RT_length)$stats
      delta_mz <- stats::median(features_select$m.z - features_select$m.z_range_min,na.rm=T)
      delta_rt <- stats::median(features_select$RT - features_select$RT_range_min,na.rm=T)/2

      known_RTs <- base::as.data.frame(stringr::str_split(features_select$Observed_RT,";",simplify = T))
      colnames(known_RTs) <- base::substr(colnames(features_select)[which(grepl("RT_calibration",colnames(features_select)))],16,1000)
      known_RTs[] <- lapply(known_RTs, function(x) as.numeric(as.character(x)))
      known_RTs_all <- known_RTs
      known_RTs <- known_RTs[,Sample_ID]

      RT_window_expand_factor <- 5 ###expand RT window around expected window for 2Dpeak detection

      if(MassSpec_mode == "TIMSToF")
      {
        IM_width_stats <- grDevices::boxplot.stats(features_select$Inv_K0_length)$stats
        delta_im <- stats::median(features_select$Inv_K0 - features_select$Inv_K0_range_min,na.rm=T)/2

        known_IMs <- base::as.data.frame(stringr::str_split(features_select$Observed_IM,";",simplify = T))
        colnames(known_IMs) <- base::substr(colnames(features_select)[which(grepl("IM_calibration",colnames(features_select)))],16,1000)
        known_IMs[] <- lapply(known_IMs, function(x) as.numeric(as.character(x)))
        known_IMs_all <- known_IMs
        known_IMs <- known_IMs[,Sample_ID]

        IM_window_expand_factor <- 5
        delta_im <- stats::median(features_select$Inv_K0 - features_select$Inv_K0_range_min,na.rm=T)/2

      }

      if(plots == T)
      {
        dir.create(base::paste(path,"/2Dpeakselection",sep=""))
        dir.create(base::paste(path,"/2Dpeakselection/",Sample_ID,sep=""))
      }

      for(i in 1:nrow(features_select))
      {
        if(MassSpec_mode == "Orbitrap")
        {
          ###extract all relevant ions in RT window
          RT_correction <- ifelse(RT_calibration == T,features_select[i,base::paste("RT_calibration.",Sample_ID,sep="")],0)

          if(peak_detection == T)# & !grepl("_d",features_select$Feature_name[i]))
          {
            RT_window <- c(features_select$RT_range_min[i] + RT_correction,
                           features_select$RT_range_max[i] + RT_correction)
            RT_window <- c(RT_window[1]-(delta_rt*RT_window_expand_factor),
                           RT_window[2]+(delta_rt*RT_window_expand_factor))
          }else
          {
            RT_window <- c(features_select$RT[i]+RT_correction-(features_select$RT_length[i]),features_select$RT[i]+RT_correction+(features_select$RT_length[i]))
          }

          ###get search range in spectra
          if(RT_window[1] <= max(Indices$RT_start) & RT_window[2] <= max(Indices$RT_end))
          {
            search_range <- (which(Indices$RT_start>=RT_window[1])[1]-1):(which(Indices$RT_end>=RT_window[2])[1])
          }else
          {
            search_range <- nrow(Indices)
          }
          if(any(search_range < 1))
          {
            search_range <- search_range[-which(search_range < 1)]
          }
        }else
        {
          ###extract all relevant ions in RT window
          RT_correction <- ifelse(RT_calibration == T,features_select[i,base::paste("RT_calibration.",Sample_ID,sep="")],0)
          IM_correction <- features_select[i,base::paste("IM_calibration.",Sample_ID,sep="")]

          if(peak_detection == T)# & !grepl("_d",features_select$Feature_name[i]))
          {
            RT_window <- c(features_select$RT_range_min[i] + RT_correction,
                           features_select$RT_range_max[i] + RT_correction)
            RT_window <- c(RT_window[1]-(delta_rt*RT_window_expand_factor),
                           RT_window[2]+(delta_rt*RT_window_expand_factor))

            IM_window <- c(features_select$Inv_K0_range_min[i] + IM_correction,
                           features_select$Inv_K0_range_max[i] + IM_correction)
            IM_window <- c(IM_window[1]-(delta_im*IM_window_expand_factor),
                           IM_window[2]+(delta_im*IM_window_expand_factor))
          }else
          {
            RT_window <- c(features_select$RT[i]+RT_correction-(features_select$RT_length[i]),features_select$RT[i]+RT_correction+(features_select$RT_length[i]))
            IM_window <- c(features_select$Inv_K0[i]+IM_correction-(features_select$Inv_K0_length[i]),features_select$Inv_K0[i]+IM_correction+(features_select$Inv_K0_length[i]))
          }

          ###get search range in spectra
          if(RT_window[1] <= max(Indices$RT_start) & RT_window[2] <= max(Indices$RT_end))
          {
            search_range <- which(Indices$RT_start <= RT_window[1] & Indices$RT_end >= RT_window[1] | Indices$RT_start >= RT_window[1] & Indices$RT_end <= RT_window[2] | Indices$RT_start >= RT_window[1] & Indices$RT_start <= RT_window[2])
            search_range_IM <- which(Indices$IM_start <= IM_window[1] & Indices$IM_end >= IM_window[1] | Indices$IM_start >= IM_window[1] & Indices$IM_end <= IM_window[2] | Indices$IM_start >= IM_window[1] & Indices$IM_start <= IM_window[2])
            search_range <- search_range[which(search_range %in% search_range_IM)]

          }else
          {
            search_range <- nrow(Indices)
          }
          if(any(search_range < 1))
          {
            search_range <- search_range[-which(search_range < 1)]
          }
        }

        if(length(search_range)>0)
        {
          sub <- data.table::rbindlist(data_frags[search_range])

          res <- NULL
          if(peak_detection == T)#& !grepl("_d",features_select$Feature_name[i])) ###use 2D peak detection
          {
            if(plots == T)
            {
              grDevices::pdf(base::paste(path,"/2Dpeakselection/",Sample_ID,"/",features_select$Feature_name[i],".pdf",sep=""))
            }

            if(MassSpec_mode == "Orbitrap")
            {
              res <- feature_2D_peak_selection(mz_bw=0.002,all_ion_data = sub,selected_features = features_select[i,],cur_sample = Sample_ID,known_RT = known_RTs[i],delta_mz = delta_mz,delta_rt=delta_rt,max_delta_RT = 2*delta_rt,max_delta_mz = 3*delta_mz,peak_min_ion_count = peak_min_ion_count,n_raster_dens_matrix = kde_resolution,num_peaks_store=num_peaks_store,plot = plots,MassSpec_mode = MassSpec_mode)
            }else
            {
              res <- feature_2D_peak_selection(mz_bw=0.01,all_ion_data = sub,selected_features = features_select[i,],cur_sample = Sample_ID,known_RT = known_RTs[i],delta_mz = delta_mz,delta_rt=delta_rt,max_delta_RT = 2*delta_rt,max_delta_mz = 3*delta_mz,peak_min_ion_count = peak_min_ion_count,n_raster_dens_matrix = kde_resolution,num_peaks_store=num_peaks_store,plot = plots,MassSpec_mode = MassSpec_mode,delta_im = delta_im,IM_window_expand_factor=IM_window_expand_factor,use_IM_data=use_IM_data)
            }

            if(plots == T)
            {
              grDevices::dev.off()
            }
            if(plots==F) ###store graphs in list
            {
              #graph_peaks[[as.character(features_select$Feature_name[i])]] <- res$graph
            }
            res$graph <- NULL
          }else ###no peak detection
          {
            res <- feature_no_peak_selection(all_ion_data = sub,selected_features = features_select[i,],cur_sample = Sample_ID,num_peaks_store=num_peaks_store,MassSpec_mode = MassSpec_mode,use_IM_data=use_IM_data)
          }

          ###check if any ions are available in any window
          temp <- unlist(res)

          if(any(temp[grepl("ion_count",names(temp))]>0))
          {
            if(peak_detection == F) ##no peak detection - simply sum all ion intensities within expected window
            {
              ####all ions
              m.z_window_final <- c(res$Standard$Peak_info$mz_window_lower,res$Standard$Peak_info$mz_window_upper)
              RT_window_final <- c(res$Standard$Peak_info$RT_window_lower,res$Standard$Peak_info$RT_window_upper)
              res_int_total <- res$Standard$ion_data$Intensity
              res_mz_total <- res$Standard$ion_data$m.z
              res_rt_total <- res$Standard$ion_data$RT
              res_isotope_total <- res$Standard$ion_data$isotope
              length_data_total <- length(res_int_total)
              if(MassSpec_mode == "TIMSToF")
              {
                IM_window_final <- c(res$Standard$Peak_info$IM_window_lower,res$Standard$Peak_info$IM_window_upper)
                res_im_total <- res$Standard$ion_data$`1/K0`
              }

              ###determine which ions are showing an intensity above the background signal intensity
              if(ion_intensity_cutoff == T)
              {
                res$Standard$ion_data$signal_background <- ifelse(base::log2(res$Standard$ion_data$Intensity) > mean_background_ion_intensity_model[i,Sample_ID]+2*sd_background_ion_intensity[1,Sample_ID],1,0)

                res_int_signal <- res_int_total[which(res$Standard$ion_data$signal_background == 1)]
                res_mz_signal <- res_mz_total[which(res$Standard$ion_data$signal_background == 1)]
                res_rt_signal <- res_rt_total[which(res$Standard$ion_data$signal_background == 1)]
                res_isotope_signal <- res_isotope_total[which(res$Standard$ion_data$signal_background == 1)]
                length_data_signal <- length(res_int_signal)

                ###save signal Intensities
                if(length(res_int_signal)>0)
                {
                  #if(MassSpec_mode == "Orbitrap")
                  {
                    data.table::set(Intensities_signal_background,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_signal)),
                                                                                                     length(res_int_signal),
                                                                                                     mean(log10(res_int_signal),na.rm=T),
                                                                                                     stats::sd(log10(res_int_signal),na.rm=T),
                                                                                                     m.z_window_final[1],
                                                                                                     m.z_window_final[2],
                                                                                                     RT_window_final[1],
                                                                                                     RT_window_final[2])))
                  }
                  # else
                  # {
                  #   data.table::set(Intensities_signal_background,i=as.integer(c(1:8,11:12)),j=as.integer(i),value=list(c(log10(sum(res_int_signal)),
                  #                                                                                             length(res_int_signal),
                  #                                                                                             mean(log10(res_int_signal),na.rm=T),
                  #                                                                                             stats::sd(log10(res_int_signal),na.rm=T),
                  #                                                                                             m.z_window_final[1],
                  #                                                                                             m.z_window_final[2],
                  #                                                                                             RT_window_final[1],
                  #                                                                                             RT_window_final[2],
                  #                                                                                             IM_window_final[1],
                  #                                                                                             IM_window_final[2])))
                  # }
                }

                ###save signal+background Intensities
                if(length(res_int_total)>0)
                {
                  #if(MassSpec_mode == "Orbitrap")
                  {
                    data.table::set(Intensities_signal_background,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                                     length(res_int_total),
                                                                                                     mean(log10(res_int_total),na.rm=T),
                                                                                                     stats::sd(log10(res_int_total),na.rm=T),
                                                                                                     m.z_window_final[1],
                                                                                                     m.z_window_final[2],
                                                                                                     RT_window_final[1],
                                                                                                     RT_window_final[2])))
                  }
                  # else
                  # {
                  #   data.table::set(Intensities_signal_background,i=as.integer(c(1:8,11:12)),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                  #                                                                                    length(res_int_total),
                  #                                                                                    mean(log10(res_int_total),na.rm=T),
                  #                                                                                    stats::sd(log10(res_int_total),na.rm=T),
                  #                                                                                    m.z_window_final[1],
                  #                                                                                    m.z_window_final[2],
                  #                                                                                    RT_window_final[1],
                  #                                                                                    RT_window_final[2],
                  #                                                                                    IM_window_final[1],
                  #                                                                                    IM_window_final[2])))
                  # }
                }
              }else ###no discrimination between signal and background ions ## save all ions
              {
                ###save signal+background Intensities
                if(length(res_int_total)>0)
                {
                  #if(MassSpec_mode == "Orbitrap")
                  {
                    data.table::set(Intensities,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                   length(res_int_total),
                                                                                   mean(log10(res_int_total),na.rm=T),
                                                                                   stats::sd(log10(res_int_total),na.rm=T),
                                                                                   m.z_window_final[1],
                                                                                   m.z_window_final[2],
                                                                                   RT_window_final[1],
                                                                                   RT_window_final[2])))
                  }
                  # else
                  # {
                  #   data.table::set(Intensities,i=as.integer(c(1:8,11:12)),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                  #                                                                  length(res_int_total),
                  #                                                                  mean(log10(res_int_total),na.rm=T),
                  #                                                                  stats::sd(log10(res_int_total),na.rm=T),
                  #                                                                  m.z_window_final[1],
                  #                                                                  m.z_window_final[2],
                  #                                                                  RT_window_final[1],
                  #                                                                  RT_window_final[2],
                  #                                                                  IM_window_final[1],
                  #                                                                  IM_window_final[2])))
                  # }
                }
              }

              ###feature alignment Scoreing based on algorithm from DeMix-Q (Zhang, 2016) - use all ions
              if(length_data_total > 0)
              {
                df <- base::data.frame(mz=res_mz_total,rt=res_rt_total,iso=res_isotope_total,int=res_int_total)
                df_mono_isotope <- df[which(df$iso==0),]
                df_isotope_1 <- df[which(df$iso==1),]
                if(nrow(df_mono_isotope)>0)
                {
                  ###deviation from consensus feature
                  max_mono <- which(df_mono_isotope$int == max(df_mono_isotope$int,na.rm=T))
                  data.table::set(delta_T1,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - (features_select$RT[i]))
                  data.table::set(delta_M1,i=1L,j=as.integer(i),mean(df_mono_isotope$mz[max_mono],na.rm=T) - (features_select$m.z[i]))
                  if(nrow(df_isotope_1)>0)
                  {
                    ###deviation from monoisotopic ion to M+1 isotope ion
                    max_iso_1 <- which(df_isotope_1$int == max(df_isotope_1$int,na.rm=T))
                    data.table::set(delta_T2,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - mean(df_isotope_1$rt[max_iso_1],na.rm=T))
                    data.table::set(delta_M2,i=1L,j=as.integer(i),(mean(df_mono_isotope$mz[max_mono],na.rm=T)+(iso_mult[1]/features_select$Charge[i])) - mean(df_isotope_1$mz[max_iso_1],na.rm=T))
                  }
                }
              }
            }else ###with peak detection
            {
              ###store up to num_peaks_store potential peak windows + the original expected window
              num_peaks <- length(which(names(res) != "Standard" & names(res) != "graph"))
              for(p in names(res)[which(names(res) != "graph")])
              {
                m.z_window_final <- c(res[[p]]$Peak_info$mz_window_lower,res[[p]]$Peak_info$mz_window_upper)
                RT_window_final <- c(res[[p]]$Peak_info$RT_window_lower,res[[p]]$Peak_info$RT_window_upper)

                res_int_total <- res[[p]]$ion_data$Intensity
                res_mz_total <- res[[p]]$ion_data$m.z
                res_rt_total <- res[[p]]$ion_data$RT
                res_isotope_total <- res[[p]]$ion_data$isotope
                length_data_total <- length(res_int_total)
                if(MassSpec_mode == "TIMSToF")
                {
                  IM_window_final <- c(res$Standard$Peak_info$IM_window_lower,res$Standard$Peak_info$IM_window_upper)
                  res_im_total <- res$Standard$ion_data$`1/K0`
                }

                ###determine which ions are showing an intensity above the background signal intensity
                if(ion_intensity_cutoff == T)
                {
                  res[[p]]$ion_data$signal_background <- ifelse(base::log2(res[[p]]$ion_data$Intensity) > mean_background_ion_intensity_model[i,Sample_ID]+2*sd_background_ion_intensity[1,Sample_ID],1,0)

                  res_int_signal <- res_int_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_mz_signal <- res_mz_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_rt_signal <- res_rt_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_isotope_signal <- res_isotope_total[which(res[[p]]$ion_data$signal_background == 1)]
                  length_data_signal <- length(res_int_signal)

                  ###save signal Intensities
                  if(length(res_int_signal)>0)
                  {
                    data.table::set(peaks_quant[[p]]$Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_signal)),
                                                                                                     length(res_int_signal),
                                                                                                     mean(log10(res_int_signal),na.rm=T),
                                                                                                     stats::sd(log10(res_int_signal),na.rm=T),
                                                                                                     m.z_window_final[1],
                                                                                                     m.z_window_final[2],
                                                                                                     RT_window_final[1],
                                                                                                     RT_window_final[2],
                                                                                                     num_peaks,
                                                                                                     res[[p]]$Peak_info$known_peak)))
                  }

                  ###save signal+background Intensities
                  if(length(res_int_total)>0)
                  {
                    data.table::set(peaks_quant[[p]]$Intensities_signal_background,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                                                       length(res_int_total),
                                                                                                                       mean(log10(res_int_total),na.rm=T),
                                                                                                                       stats::sd(log10(res_int_total),na.rm=T),
                                                                                                                       m.z_window_final[1],
                                                                                                                       m.z_window_final[2],
                                                                                                                       RT_window_final[1],
                                                                                                                       RT_window_final[2],
                                                                                                                       num_peaks,
                                                                                                                       res[[p]]$Peak_info$known_peak)))
                  }else ###no intensity at all then at least save standard information
                  {
                    data.table::set(peaks_quant[[p]]$Intensities_signal_background,i=as.integer(1:10),j=as.integer(i),value=list(c(NA,
                                                                                                                       length(res_int_total),
                                                                                                                       NA,
                                                                                                                       NA,
                                                                                                                       m.z_window_final[1],
                                                                                                                       m.z_window_final[2],
                                                                                                                       RT_window_final[1],
                                                                                                                       RT_window_final[2],
                                                                                                                       num_peaks,
                                                                                                                       res[[p]]$Peak_info$known_peak)))
                  }
                }else ###no discrimination between signal and background ions ## save all ions
                {
                  ###save signal+background Intensities
                  if(length(res_int_total)>0)
                  {
                    data.table::set(Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                    length(res_int_total),
                                                                                    mean(log10(res_int_total),na.rm=T),
                                                                                    stats::sd(log10(res_int_total),na.rm=T),
                                                                                    m.z_window_final[1],
                                                                                    m.z_window_final[2],
                                                                                    RT_window_final[1],
                                                                                    RT_window_final[2],
                                                                                    num_peaks,
                                                                                    res[[p]]$Peak_info$known_peak)))
                  }else ###no intensity at all then at least save standard information
                  {
                    data.table::set(Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(NA,
                                                                                    length(res_int_total),
                                                                                    NA,
                                                                                    NA,
                                                                                    m.z_window_final[1],
                                                                                    m.z_window_final[2],
                                                                                    RT_window_final[1],
                                                                                    RT_window_final[2],
                                                                                    num_peaks,
                                                                                    res[[p]]$Peak_info$known_peak)))
                  }
                }

                ###feature alignment Scoreing based on algorithm from DeMix-Q (Zhang, 2016) - use all ions
                if(length_data_total > 0)
                {
                  df <- base::data.frame(mz=res_mz_total,rt=res_rt_total,iso=res_isotope_total,int=res_int_total)
                  df_mono_isotope <- df[which(df$iso==0),]
                  df_isotope_1 <- df[which(df$iso==1),]
                  if(nrow(df_mono_isotope)>0)
                  {
                    ###deviation from consensus feature
                    max_mono <- which(df_mono_isotope$int == max(df_mono_isotope$int,na.rm=T))
                    data.table::set(peaks_quant[[p]]$delta_T1,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - (features_select$RT[i]))
                    data.table::set(peaks_quant[[p]]$delta_M1,i=1L,j=as.integer(i),mean(df_mono_isotope$mz[max_mono],na.rm=T) - (features_select$m.z[i]))
                    if(nrow(df_isotope_1)>0)
                    {
                      ###deviation from monoisotopic ion to M+1 isotope ion
                      max_iso_1 <- which(df_isotope_1$int == max(df_isotope_1$int,na.rm=T))
                      data.table::set(peaks_quant[[p]]$delta_T2,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - mean(df_isotope_1$rt[max_iso_1],na.rm=T))
                      data.table::set(peaks_quant[[p]]$delta_M2,i=1L,j=as.integer(i),(mean(df_mono_isotope$mz[max_mono],na.rm=T)+(iso_mult[1]/features_select$Charge[i])) - mean(df_isotope_1$mz[max_iso_1],na.rm=T))
                    }
                  }
                }

              }
            }
          }

        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 5)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d:%02d',  td@day, td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }

      }
      close(pb)

      if(peak_detection == F)
      {
        Intensities <- t(Intensities)
        Intensities[Intensities == 0] <- NA

        if(ion_intensity_cutoff == T)
        {
          Intensities_signal_background <- t(Intensities_signal_background)
          Intensities_signal_background[Intensities_signal_background == 0] <- NA
        }

        delta_T1 <- t(delta_T1)
        delta_T1[delta_T1 == 0] <- NA
        delta_M1 <- t(delta_M1)
        delta_M1[delta_M1 == 0] <- NA
        delta_T2 <- t(delta_T2)
        delta_T2[delta_T2 == 0] <- NA
        delta_M2 <- t(delta_M2)
        delta_M2[delta_M2 == 0] <- NA

        peaks_quant <- list()
        if(ion_intensity_cutoff == T)
        {
          peaks_quant[[1]] <- list(Intensities=Intensities,
                                   Intensities_signal_background=Intensities_signal_background,
                                   delta_T1=delta_T1,
                                   delta_T2=delta_T2,
                                   delta_M1=delta_M1,
                                   delta_M2=delta_M2)
        }else
        {
          peaks_quant[[1]] <- list(Intensities=Intensities,
                                   delta_T1=delta_T1,
                                   delta_T2=delta_T2,
                                   delta_M1=delta_M1,
                                   delta_M2=delta_M2)
        }

        return(list(peaks_quant=peaks_quant,graph_peaks=graph_peaks))

      }else
      {
        for(p in 1:(num_peaks_store+1))
        {
          peaks_quant[[p]]$Intensities <- t(peaks_quant[[p]]$Intensities)
          peaks_quant[[p]]$Intensities[peaks_quant[[p]]$Intensities == 0] <- NA

          peaks_quant[[p]]$Intensities_signal_background <- t(peaks_quant[[p]]$Intensities_signal_background)
          peaks_quant[[p]]$Intensities_signal_background[peaks_quant[[p]]$Intensities_signal_background == 0] <- NA

          peaks_quant[[p]]$delta_T1 <- t(peaks_quant[[p]]$delta_T1)
          peaks_quant[[p]]$delta_T1[peaks_quant[[p]]$delta_T1 == 0] <- NA
          peaks_quant[[p]]$delta_M1 <- t(peaks_quant[[p]]$delta_M1)
          peaks_quant[[p]]$delta_M1[peaks_quant[[p]]$delta_M1 == 0] <- NA
          peaks_quant[[p]]$delta_T2 <- t(peaks_quant[[p]]$delta_T2)
          peaks_quant[[p]]$delta_T2[peaks_quant[[p]]$delta_T2 == 0] <- NA
          peaks_quant[[p]]$delta_M2 <- t(peaks_quant[[p]]$delta_M2)
          peaks_quant[[p]]$delta_M2[peaks_quant[[p]]$delta_M2 == 0] <- NA
        }

        return(list(peaks_quant=peaks_quant,graph_peaks=graph_peaks))

      }

    }

    ####Function to call get_intensities and finally save resulting table as .tab data table
    extract_intensities <- function(Sample_ID,features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection,ion_intensity_cutoff,mean_background_ion_intensity_model,sd_background_ion_intensity,peak_min_ion_count,kde_resolution,num_peaks_store,plots,MassSpec_mode,use_IM_data)
    {
      #suppressWarnings(suppressMessages(library(data.table,quietly = T)))

      res <- get_intensities(Sample_ID,path = path_to_raw,features_select=features_select,RT_calibration=RT_calibration,mz_calibration=mz_calibration,peak_detection=peak_detection,ion_intensity_cutoff = ion_intensity_cutoff,mean_background_ion_intensity_model=mean_background_ion_intensity_model,sd_background_ion_intensity=sd_background_ion_intensity,peak_min_ion_count=peak_min_ion_count,kde_resolution=kde_resolution,num_peaks_store = num_peaks_store,plots=plots,MassSpec_mode=MassSpec_mode,use_IM_data=use_IM_data)

      peaks_quant <- res$peaks_quant
      #peaks_graph <- res$graph_peaks

      save(peaks_quant,file = base::paste(path_to_output_folder,"/",Sample_ID,"_feature_quant.RData",sep=""))
      #if(length(peaks_graph)>0)save(peaks_graph,file = base::paste(path_to_output_folder,"/",Sample_ID,"_feature_graphs.RData",sep=""))

    }

    ####Prepare threads to run extraction. The task is using much memory so that more than 2 threads in parallel on a pc with 16 gb of ram
    ####results in slower performance than for just 2 threads
    cl <- parallel::makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=Sample_IDs) %dopar%
      {
        extract_intensities(Sample_ID = i,features_select = features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection=peak_detection,ion_intensity_cutoff = ion_intensity_cutoff,mean_background_ion_intensity_model=mean_background_ion_intensity_model,sd_background_ion_intensity=sd_background_ion_intensity,peak_min_ion_count=peak_min_ion_count,kde_resolution=kde_resolution,num_peaks_store=num_peaks_store,plots=plots,MassSpec_mode=MassSpec_mode,use_IM_data=use_IM_data)
      }
    parallel::stopCluster(cl)

    # for(i in Sample_IDs)
    # {
    #   extract_intensities(Sample_ID = i,features_select = features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection=peak_detection,ion_intensity_cutoff = ion_intensity_cutoff,mean_background_ion_intensity_model=mean_background_ion_intensity_model,sd_background_ion_intensity=sd_background_ion_intensity,peak_min_ion_count=peak_min_ion_count,kde_resolution=kde_resolution,num_peaks_store=num_peaks_store,plots=plots,MassSpec_mode=MassSpec_mode,use_IM_data=use_IM_data)
    # }

  }

  if(MassSpec_mode == "Orbitrap")
  {
    ###check for which samples mzXML files are available
    mzXMLfiles <- list.files(base::paste(path_to_mzXML,"/all_ion_lists",sep=""))
    mzXMLfiles <- mzXMLfiles[which(grepl(".RData",mzXMLfiles))]
    samples <- mzXMLfiles
    samples <- base::substr(samples,1,regexpr("_all_ions.RData",samples)-1)
  }
  if(MassSpec_mode == "TIMSToF")
  {
    ###check for which samples extracted spectra files are available
    mzXMLfiles <- list.files(path_to_extracted_spectra)
    mzXMLfiles <- mzXMLfiles[which(grepl(".RData",mzXMLfiles))]
    samples <- mzXMLfiles
    samples <- base::substr(samples,1,regexpr("_all_ions.RData",samples)-1)
    path_to_mzXML <- base::gsub("\\\\all_ion_lists|/all_ion_lists","",path_to_extracted_spectra)
  }

  ###keep samples which should be actually requantified
  Requant_samples <- base::gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration\\.",colnames(features)))])
  Requant_samples <- base::gsub("\\.","-",Requant_samples)

  #check in which quant mode ... LFQ or SILAC
  if(any(grepl("_Channel_light|_Channel_medium|_Channel_heavy",Requant_samples)))
  {
    if(any(grepl("_Channel_medium",Requant_samples)))
    {
      multiplicity <- 3
      labels <- c("_Channel_light","_Channel_medium","_Channel_heavy")
    }else
    {
      multiplicity <- 2
      labels <- c("_Channel_light","Channel_heavy")
    }
    samples <- base::paste(sort(rep(samples,multiplicity)),labels,sep="")
  }else
  {
    multiplicity <- 1
  }

  samples <- samples[which(samples %in% Requant_samples)]

  ###Perform quantification of decoy features

  dir.create(base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,sep=""),showWarnings = F)
  available <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,sep=""))
  available <- base::gsub("_feature_quant.RData","",available)
  if(length(which(samples %not in% available))>0)
  {
    selected_decoys <- which(grepl("_d",features$Feature_name))

    samples <- samples[which(samples %not in% available)]
    #single core in case of TIMSToF data as it fills memory too much
    extract_intensities_worker(Sample_IDs = as.character(samples),
                               features_select = features[selected_decoys,],
                               path_to_raw = base::paste(path_to_mzXML,"/all_ion_lists",sep=""),
                               path_to_output_folder = base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,sep=""),
                               RT_calibration=RT_calibration,
                               mz_calibration=mz_calibration,
                               peak_detection=F,
                               n_cores=ifelse(MassSpec_mode == "Orbitrap",n_cores,ifelse(n_cores >= 3,3,n_cores)),
                               MassSpec_mode = MassSpec_mode,
                               use_IM_data=use_IM_data)
  }

  ###determine distribution of background ion intensities
  samples <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,sep=""))
  samples <- samples[which(grepl("_feature_quant.RData",samples))]
  samples <- base::substr(samples,1,regexpr("_feature_quant.RData",samples)-1)

  files <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,sep=""))
  if(length(which(!grepl("_feature_quant.RData",files)))>0)
  {
    files <- files[-which(!grepl("_feature_quant.RData",files))]
  }

  features_select = features[which(grepl("_d",features$Feature_name)),]

  decoy_intensities <- base::as.data.frame(matrix(ncol=length(samples),nrow=nrow(features_select)))
  colnames(decoy_intensities) <- samples
  decoy_intensities <- sapply( decoy_intensities, as.numeric )
  decoy_intensities <- data.table::as.data.table(decoy_intensities)
  rownames(decoy_intensities) <- features_select$Feature_name

  decoy_ioncount <- decoy_intensities
  decoy_mean_intensity <- decoy_intensities
  decoy_sd_intensity <- decoy_intensities

  for(c in 1:ncol(decoy_intensities))
  {
    #load stored data into variable peaks_quant
    peaks_quant <- NULL
    load(base::paste(path_to_mzXML,"/all_ion_lists/Extracted decoy intensities",output_file_names_add,"/",colnames(decoy_intensities)[c],"_feature_quant.RData",sep=""))

    signal=peaks_quant[[1]]$Intensities

    decoy_intensities[,c] <- base::log2(10^signal[,1])
    decoy_ioncount[,c] <- signal[,2]
    decoy_mean_intensity[,c] <- base::log2(10^signal[,3])
    decoy_sd_intensity[,c] <- base::log2(10^signal[,4])
  }

  ###plot general numbers of quantifications of decoy features
  setwd(path_to_features)
  grDevices::pdf("Temporary_files/Decoy feature quantification parameters.pdf")

  ###mean decoy intensity (intensity of a single decoy ion)
  RT_all <- rep(as.numeric(features_select$RT),ncol(decoy_mean_intensity))
  x_all <- as.numeric(as.matrix(decoy_mean_intensity))
  graphics::smoothScatter(RT_all,x_all,ylab="Intensity, log2",main="All samples - Decoy feature mean intensity",xlab="RT [min]")

  ##try to fit an average generalised additive model to determine a RT dependent mean intensity and sd of intensity
  fit_gam_mean <- mgcv::gam(x_all ~ s(RT_all), method = "REML")
  x_pred <- seq(min(features_select$RT,na.rm=T), max(features_select$RT,na.rm=T), length.out = nrow(features_select))
  y_pred <- stats::predict(fit_gam_mean, base::data.frame(RT_all = x_pred))
  graphics::lines(x_pred,y_pred,col="red")
  graphics::legend("topright",legend="GAM",lty=1,col="red")

  ##now fit gam models per sample
  fit_gam_per_sample <- list()
  for(c in 1:ncol(decoy_mean_intensity))
  {
    RT <- as.numeric(features_select$RT)
    x <- as.numeric(as.matrix(decoy_mean_intensity)[,c])
    graphics::smoothScatter(RT,x,ylab="Mean intensity, log2",main=base::paste(colnames(decoy_mean_intensity)[c],"Decoy feature intensity"),xlab="RT [min]")

    ##try to fit an average generalised additive model to determine a RT dependent mean intensity and sd of intensity
    ##if not enough data points available, use the average gam
    gam <- tryCatch({mgcv::gam(x ~ s(RT), method = "REML")}, error = function(error_condition) {
      RT <- rep(as.numeric(features_select$RT),ncol(decoy_mean_intensity))
      x <- as.numeric(as.matrix(decoy_mean_intensity))
      mgcv::gam(x ~ s(RT), method = "REML")
    })
    x_pred <- seq(min(features_select$RT,na.rm=T), max(features_select$RT,na.rm=T), length.out = nrow(features_select))
    y_pred <- stats::predict(gam, base::data.frame(RT = x_pred))
    graphics::lines(x_pred,y_pred,col="red")
    graphics::legend("topright",legend="GAM",lty=1,col="red")
    fit_gam_per_sample[[colnames(decoy_mean_intensity)[c]]] <- gam
  }

  graphics::par(mfrow=c(2,2))
  graphics::boxplot(as.numeric(as.matrix(decoy_intensities)),outline=F,ylab="Summed intensity, log2",main="Summed intensity of decoy ions")
  graphics::boxplot(as.numeric(as.matrix(decoy_mean_intensity)),outline=F,ylab="Mean intensity, log2",main="Mean intensity of decoy ions")
  graphics::boxplot(as.numeric(as.matrix(decoy_sd_intensity)),outline=F,ylab="SD of intensity, log2",main="SD of intensity of decoy ions")
  graphics::boxplot(as.numeric(as.matrix(decoy_ioncount)),outline=F,ylab="Count",main="Number of ions in decoys with quantification")
  graphics::par(mfrow=c(1,1))
  grDevices::dev.off()

  ###Now predict background intensities per sample and feature by using individually fitted GAMs, multiply with median decoy ion count and add some noise using observed decoy intensity sds per sample
  background_intensity_GAM_table_per_feature <- NULL ###used for determining which ions are background ions and which are signal ions
  background_intensity_GAM_table_per_feature_sum <- NULL ###used to later impute missing values

  median_decoy_sd_per_sample <- base::as.data.frame(t(matrixStats::colMedians(as.matrix(decoy_sd_intensity),na.rm=T)))
  colnames(median_decoy_sd_per_sample) <- colnames(decoy_sd_intensity)

  median_decoy_ion_count_per_sample <- base::as.data.frame(t(matrixStats::colMedians(as.matrix(decoy_ioncount),na.rm=T)))
  colnames(median_decoy_ion_count_per_sample) <- colnames(decoy_ioncount)

  max <- ncol(decoy_mean_intensity)
  pb <- tcltk::tkProgressBar(title = "Model background intensities per feature",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0

  for(c in 1:ncol(decoy_mean_intensity))
  {
    cur_sample <- colnames(decoy_mean_intensity)[c]
    if(c != 1)
    {
      background_intensity_GAM_table_per_feature <- cbind(background_intensity_GAM_table_per_feature,
                                                          base::data.frame(mean_intensity=stats::predict(fit_gam_per_sample[[cur_sample]], base::data.frame(RT = features$RT))))
      background_intensity_GAM_table_per_feature_sum <- cbind(background_intensity_GAM_table_per_feature_sum,
                                                              base::data.frame(mean_intensity=stats::predict(fit_gam_per_sample[[cur_sample]], base::data.frame(RT = features$RT))))
    }else
    {
      background_intensity_GAM_table_per_feature <- base::data.frame(mean_intensity=stats::predict(fit_gam_per_sample[[cur_sample]], base::data.frame(RT = features$RT)))
      background_intensity_GAM_table_per_feature_sum <- base::data.frame(mean_intensity=stats::predict(fit_gam_per_sample[[cur_sample]], base::data.frame(RT = features$RT)))
    }
    colnames(background_intensity_GAM_table_per_feature)[ncol(background_intensity_GAM_table_per_feature)] <- cur_sample
    colnames(background_intensity_GAM_table_per_feature_sum)[ncol(background_intensity_GAM_table_per_feature_sum)] <- cur_sample

    background_intensity_GAM_table_per_feature_sum[,c] <- stats::rnorm(n=nrow(background_intensity_GAM_table_per_feature_sum),mean = background_intensity_GAM_table_per_feature_sum[,c],sd = median_decoy_sd_per_sample[1,c])

    updatecounter <- updatecounter + 1
    if(updatecounter >= 1)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(c/max))*(1-(c/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, c, label=base::paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)

  sd_background_intensity <- median_decoy_sd_per_sample

  rownames(background_intensity_GAM_table_per_feature) <- features$Feature_name
  rownames(background_intensity_GAM_table_per_feature_sum) <- features$Feature_name

  QC_data[["Decoy_feature_parameters"]] <- list(Decoy_mean_intensity_per_RT = base::data.frame(RT_all=RT_all,
                                                                                         mean_intensity=x_all),
                                                GAM_model_average=fit_gam_mean,
                                                GAM_model_per_sample=fit_gam_per_sample,
                                                background_intensity_GAM_table_per_feature=background_intensity_GAM_table_per_feature,
                                                background_intensity_GAM_table_per_feature_sum=background_intensity_GAM_table_per_feature_sum,
                                                sd_background_intensity=sd_background_intensity,
                                                median_decoy_ion_count_per_sample=median_decoy_ion_count_per_sample,
                                                decoy_intensities=decoy_intensities,
                                                decoy_mean_intensity=decoy_mean_intensity,
                                                decoy_sd_intensity=decoy_sd_intensity,
                                                decoy_ioncount=decoy_ioncount)

  ###now select 1000 decoy features randomly for which we will perform peak selection and quantification
  decoy_features <- which(grepl("_d",features$Feature_name))
  if(length(decoy_features)>1000)
  {
    set.seed(1)
    select_decoys <- sample(decoy_features,1000)
    decoy_features_keep <- features[select_decoys,]
    ###now add isotope peaks of these decoy features
    isotope_features <- decoy_features_keep
    isotope_features$m.z <- ((isotope_features$m.z*isotope_features$Charge)+1.002054)/isotope_features$Charge
    delta_mz <- isotope_features$m.z_range_max-isotope_features$m.z_range_min
    isotope_features$m.z_range_max <- isotope_features$m.z+(delta_mz/2)
    isotope_features$m.z_range_min <- isotope_features$m.z-(delta_mz/2)
    isotope_features$Feature_name <- base::paste(isotope_features$Feature_name,"_i",sep="")
    decoy_features_keep <- rbind(decoy_features_keep,isotope_features)

  }else ##less than 1000 decoys so just continue
  {
    decoy_features_keep <- features[decoy_features,]
    ###now add isotope peaks of these decoy features
    isotope_features <- decoy_features_keep
    isotope_features$m.z <- ((isotope_features$m.z*isotope_features$Charge)+1.002054)/isotope_features$Charge
    delta_mz <- isotope_features$m.z_range_max-isotope_features$m.z_range_min
    isotope_features$m.z_range_max <- isotope_features$m.z+(delta_mz/2)
    isotope_features$m.z_range_min <- isotope_features$m.z-(delta_mz/2)
    isotope_features$Feature_name <- base::paste(isotope_features$Feature_name,"_i",sep="")
    decoy_features_keep <- rbind(decoy_features_keep,isotope_features)
  }
  ###remove all decoy features
  features <- features[-decoy_features,]
  ###add selected decoy features plus isotope decoy features
  features <- rbind(features,decoy_features_keep)

  ###see if already samples were converted
  samples <- mzXMLfiles
  samples <- base::substr(samples,1,regexpr("_all_ions.RData",samples)-1)
  Requant_samples <- base::gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration\\.",colnames(features)))])
  Requant_samples <- base::gsub("\\.","-",Requant_samples)
  if(multiplicity > 1)
  {
    samples <- base::paste(sort(rep(samples,multiplicity)),labels,sep="")
  }
  samples <- samples[which(samples %in% Requant_samples)]

  dir.create(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,sep=""),showWarnings = F)
  available <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,sep=""))
  available <- available[which(grepl("feature_quant.RData",available))]
  available <- base::gsub("_feature_quant.RData","",available)
  peak_min_ion_count <- grDevices::boxplot.stats(as.numeric(as.matrix(decoy_ioncount)))$stats[4]
  if(length(which(samples %not in% available))>0)
  {
    samples <- samples[which(samples %not in% available)]
    ##use decoy defined cut of to distinguish background intensity from signal intensity per feature
    extract_intensities_worker(Sample_IDs = as.character(samples),
                               features_select = features,
                               path_to_raw = base::paste(path_to_mzXML,"/all_ion_lists",sep=""),
                               path_to_output_folder = base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,sep=""),
                               RT_calibration=RT_calibration,
                               mz_calibration=mz_calibration,
                               peak_detection=peak_detection,
                               n_cores=ifelse(MassSpec_mode == "Orbitrap",n_cores,ifelse(n_cores >= 3,3,n_cores)),
                               ion_intensity_cutoff = T,
                               mean_background_ion_intensity_model = background_intensity_GAM_table_per_feature,
                               sd_background_ion_intensity = sd_background_intensity,
                               peak_min_ion_count=peak_min_ion_count,
                               kde_resolution = kde_resolution,
                               num_peaks_store = num_peaks_store,
                               plots = plot_peak_detection,
                               MassSpec_mode = MassSpec_mode,
                               use_IM_data=use_IM_data) ###define background peaks during 2DKDE as peaks with <= 75% quantile
  }

  #####Summarize data for all features and samples
  samples <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,sep=""))
  samples <- samples[which(grepl("feature_quant.RData",samples))]
  samples <- base::substr(samples,1,regexpr("_feature_quant.RData",samples)-1)

  ###Available sample data
  files <- list.files(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,sep=""))
  if(length(which(!grepl("feature_quant.RData",files)))>0)
  {
    files <- files[-which(!grepl("feature_quant.RData",files))]
  }

  features_select <- features

  features_intensity <- base::as.data.frame(matrix(ncol=length(samples),nrow=nrow(features_select)))#
  colnames(features_intensity) <- samples
  features_intensity <- sapply( features_intensity, as.numeric )
  features_intensity <- data.table::as.data.table(features_intensity)
  rownames(features_intensity) <- features_select$Feature_name
  Ioncount_feature_sample_matrix <- features_intensity

  feature_with_background_intensity <- features_intensity
  Ioncount_feature_with_background_intensity <- features_intensity

  ###prepare matrices to store data
  if(peak_detection == T)
  {
    peak_quant <- list()
    for(p in 1:(num_peaks_store+1))
    {
      features_intensity_temp <- features_intensity
      feature_with_background_intensity_temp <- feature_with_background_intensity
      Ioncount_feature_sample_matrix_temp <- Ioncount_feature_sample_matrix
      Ioncount_feature_with_background_intensity_temp <- Ioncount_feature_with_background_intensity
      #dT1_temp <- dT1
      #dM1_temp <- dM1
      Peak_rt_temp <- features_intensity
      Peak_mz_temp <- features_intensity
      num_peaks_temp <- features_intensity
      correct_peak_temp <- features_intensity
      Peak_rt_with_background_temp <- features_intensity
      Peak_mz_with_background_temp <- features_intensity
      num_peaks_with_background_temp <- features_intensity
      correct_peak_with_background_temp <- features_intensity
      peak_quant[[p]] <- list(features_intensity=features_intensity_temp,
                              feature_with_background_intensity=feature_with_background_intensity_temp,
                              Ioncount_feature_sample_matrix=Ioncount_feature_sample_matrix_temp,
                              Ioncount_feature_with_background_intensity=Ioncount_feature_with_background_intensity_temp,
                              #dT1=dT1_temp,
                              #dM1=dM1_temp,
                              Peak_rt=Peak_rt_temp,
                              Peak_mz=Peak_mz_temp,
                              num_peaks=num_peaks_temp,
                              correct_peak=correct_peak_temp,
                              Peak_rt_with_background=Peak_rt_with_background_temp,
                              Peak_mz_with_background=Peak_mz_with_background_temp,
                              num_peaks_with_background=num_peaks_with_background_temp,
                              correct_peak_with_background=correct_peak_with_background_temp)
    }
    names(peak_quant) <- c(base::paste("Peak_",1:(num_peaks_store),sep=""),"Standard")
  }

  max <- ncol(features_intensity)
  pb <- tcltk::tkProgressBar(title = "Merge quantification results",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0

  ###read the data
  for(c in 1:ncol(features_intensity))
  {
    if(peak_detection == F)
    {
      ###Now extract intensities, row1 = Intensity summed, row2 = num ions, row3 = mean intensity, row4 = sd intensity,
      ###row5 = detected optimum mz window minimum
      ###row6 = detected optimum mz window maximum
      ###row7 = detected optimum RT window minimum
      ###row8 = detected optimum RT window maximum
      ###row9 = number of detected peaks
      ###row10 = index which of the peaks should be the correct one (based on known observed RT for the respective feature)

      #load stored data into variable peaks_quant
      peaks_quant <- NULL
      load(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,"/",colnames(features_intensity)[c],"_feature_quant.RData",sep=""))

      signal=peaks_quant[[1]]$Intensities
      signal_background=peaks_quant[[1]]$Intensities_signal_background

      features_intensity[,c] <- signal[,1]
      Ioncount_feature_sample_matrix[,c] <- signal[,2]

      feature_with_background_intensity[,c] <- signal_background[,1]
      Ioncount_feature_with_background_intensity[,c] <- signal_background[,2]

    }else
    {
      #load stored data into variable peaks_quant
      load(base::paste(path_to_mzXML,"/all_ion_lists/Extracted feature intensities",output_file_names_add,"/",colnames(features_intensity)[c],"_feature_quant.RData",sep=""))

      for(p in 1:(num_peaks_store+1))
      {
        ###Now extract intensities, row1 = Intensity summed, row2 = num ions, row3 = mean intensity, row4 = sd intensity,
        ###row5 = detected optimum mz window minimum
        ###row6 = detected optimum mz window maximum
        ###row7 = detected optimum RT window minimum
        ###row8 = detected optimum RT window maximum
        ###row9 = number of detected peaks
        ###row10 = index which of the peaks should be the correct one (based on known observed RT for the respective feature)

        signal=peaks_quant[[p]]$Intensities
        signal_background=peaks_quant[[p]]$Intensities_signal_background

        peak_quant[[p]]$features_intensity[,c] <- signal[,1]
        peak_quant[[p]]$Ioncount_feature_sample_matrix[,c] <- signal[,2]
        peak_quant[[p]]$Peak_rt[,c] <- (signal[,7]+signal[,8])/2
        peak_quant[[p]]$Peak_mz[,c] <- (signal[,5]+signal[,6])/2
        peak_quant[[p]]$num_peaks[,c] <- signal[,9]
        peak_quant[[p]]$correct_peak[,c] <- signal[,10]

        peak_quant[[p]]$feature_with_background_intensity[,c] <- signal_background[,1]
        peak_quant[[p]]$Ioncount_feature_with_background_intensity[,c] <- signal_background[,2]
        peak_quant[[p]]$Peak_rt_with_background[,c] <- (signal_background[,7]+signal_background[,8])/2
        peak_quant[[p]]$Peak_mz_with_background[,c] <- (signal_background[,5]+signal_background[,6])/2
        peak_quant[[p]]$num_peaks_with_background[,c] <- signal_background[,9]
        peak_quant[[p]]$correct_peak_with_background[,c] <- signal_background[,10]

      }
    }
    updatecounter <- updatecounter + 1
    if(updatecounter >= 1)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(c/max))*(1-(c/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, c, label=base::paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)

  #if peak detection was performed, then check if peak selection was already done. If yes, load stored data, if not, perform peak selection
  if(file.exists(base::paste(path_to_features,"/Temporary_files/Quantification_raw_results.RData",sep="")))
  {
    load(base::paste(path_to_features,"/Temporary_files/Quantification_raw_results.RData",sep=""))
    crap <- gc(F)
  }else
  {
    ###if no peak detection was done, directly continue
    ###otherwise select per feature best peak
    if(peak_detection == T)
    {
      features_intensity <- as.matrix(features_intensity)
      Ioncount_feature_sample_matrix <- as.matrix(Ioncount_feature_sample_matrix)
      feature_with_background_intensity <- as.matrix(feature_with_background_intensity)
      Ioncount_feature_with_background_intensity <- as.matrix(Ioncount_feature_with_background_intensity)

      ###prepare matrix in which we store decisions which peaks where used for quantification
      peak_selected <- matrix(ncol=ncol(feature_with_background_intensity),nrow=nrow(feature_with_background_intensity))
      colnames(peak_selected) <- colnames(feature_with_background_intensity)

      ###transfer all quantification data from Standard for quantifications where no peak was detected
      selection <- which(is.na(as.matrix(peak_quant$Standard$num_peaks_with_background)) & is.na(as.matrix(peak_quant$Peak_1$num_peaks_with_background)))
      features_intensity[selection] <- as.matrix(peak_quant$Standard$features_intensity)[selection]
      Ioncount_feature_sample_matrix[selection] <- as.matrix(peak_quant$Standard$Ioncount_feature_sample_matrix)[selection]
      feature_with_background_intensity[selection] <- as.matrix(peak_quant$Standard$feature_with_background_intensity)[selection]
      Ioncount_feature_with_background_intensity[selection] <- as.matrix(peak_quant$Standard$Ioncount_feature_with_background_intensity)[selection]
      peak_selected[selection] <- num_peaks_store+1 ###standard window

      ###next, transfer all quantification data from Peak1 for quantifications where correct peak is known
      selection <- which(as.matrix(peak_quant$Peak_1$correct_peak_with_background == 1))
      features_intensity[selection] <- as.matrix(peak_quant$Peak_1$features_intensity)[selection]
      Ioncount_feature_sample_matrix[selection] <- as.matrix(peak_quant$Peak_1$Ioncount_feature_sample_matrix)[selection]
      feature_with_background_intensity[selection] <- as.matrix(peak_quant$Peak_1$feature_with_background_intensity)[selection]
      Ioncount_feature_with_background_intensity[selection] <- as.matrix(peak_quant$Peak_1$Ioncount_feature_with_background_intensity)[selection]
      peak_selected[selection] <- 1 ###closest peak = peak 1

      ###if no peak is know, change from NA to 0
      peak_quant$Peak_1$correct_peak[is.na(peak_quant$Peak_1$correct_peak_with_background)] <- 0
      peak_quant$Peak_1$correct_peak_with_background[is.na(peak_quant$Peak_1$correct_peak_with_background)] <- 0

      features_intensity <- base::as.data.frame(features_intensity)
      Ioncount_feature_sample_matrix <- base::as.data.frame(Ioncount_feature_sample_matrix)
      feature_with_background_intensity <- base::as.data.frame(feature_with_background_intensity)
      Ioncount_feature_with_background_intensity <- base::as.data.frame(Ioncount_feature_with_background_intensity)
      peak_selected <- base::as.data.frame(peak_selected)

      ###For the next step of peak quantifications where true peak is not known we will need RT and mz correction factors per feature
      ###Extract RT and mz correction factors per sample and feature
      indices_RT_correction <- which(grepl("RT_calibration",colnames(features_select)))
      indices_mz_correction <- which(grepl("mz_calibration",colnames(features_select)))
      ordering_indices <- match(base::gsub("-",".",samples),base::gsub("RT_calibration\\.","",colnames(features_select)[indices_RT_correction]))
      indices_RT_correction <- indices_RT_correction[ordering_indices]
      indices_mz_correction <- indices_mz_correction[ordering_indices]

      RT_correction_factors <- features_select[,indices_RT_correction]
      mz_correction_factors <- features_select[,indices_mz_correction]
      colnames(RT_correction_factors) <- samples
      rownames(RT_correction_factors) <- features_select$Feature_name
      colnames(mz_correction_factors) <- samples
      rownames(mz_correction_factors) <- features_select$Feature_name

      ###general RT and mz windows
      delta_mz <- stats::median(features_select$m.z - features_select$m.z_range_min,na.rm=T)
      delta_rt <- stats::median(features_select$RT - features_select$RT_range_min,na.rm=T)/2

      ###split up quant results into chunks of 50000 features
      ###perform individually peak selection per chunk
      ###combine results afterwards again
      ###purpose: keep peak quant file size low otherwise memory can be quickly overloaded in case of a large experiment
      chunk_size <- 50000
      num_chunks <- ceiling(nrow(features)/chunk_size)
      for(chunk in 1:num_chunks)
      {
        start <- (chunk-1)*chunk_size+1
        end <- (chunk)*chunk_size
        if(end > nrow(features))end <- nrow(features)

        features_temp <- features_select[start:end,]

        peak_quant_temp <- peak_quant
        for(p in 1:(num_peaks_store+1))
        {
          for(p2 in 1:length(peak_quant_temp[[p]]))
          {
            peak_quant_temp[[p]][[p2]] <- peak_quant_temp[[p]][[p2]][start:end,]
          }
        }

        ######perform peak decision
        #suppressWarnings(suppressMessages(library(doParallel,quietly = T)))

        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)

        res <- foreach::foreach(s=1:length(samples)) %dopar%
          {
            peak_decision(features_select = features_temp,
                          peak_quant = peak_quant_temp,
                          samples,
                          s = s,
                          RT_correction_factors = RT_correction_factors[start:end,],
                          mz_correction_factors = mz_correction_factors[start:end,],
                          features_intensity_sample = features_intensity[start:end,s,drop=F],
                          Ioncount_sample = Ioncount_feature_sample_matrix[start:end,s,drop=F],
                          feature_with_background_intensity_sample = feature_with_background_intensity[start:end,s,drop=F],
                          Ioncount_with_background_sample = Ioncount_feature_with_background_intensity[start:end,s,drop=F],
                          peak_selected_sample = peak_selected[start:end,s,drop=F],
                          delta_mz = delta_mz,
                          delta_rt = delta_rt,
                          peak_min_ion_count=peak_min_ion_count,
                          chunk=chunk,
                          num_chunks=num_chunks)
          }
        parallel::stopCluster(cl)

        ###merge results from all threads
        for(c in 1:length(samples))
        {
          data.table::set(features_intensity,as.integer(start:end),c,res[[c]]$features_intensity_sample)
          data.table::set(Ioncount_feature_sample_matrix,as.integer(start:end),c,res[[c]]$Ioncount_sample)
          data.table::set(feature_with_background_intensity,as.integer(start:end),c,res[[c]]$feature_with_background_intensity_sample)
          data.table::set(Ioncount_feature_with_background_intensity,as.integer(start:end),c,res[[c]]$Ioncount_with_background_sample)
          data.table::set(peak_selected,as.integer(start:end),c,res[[c]]$peak_selected_sample)
        }

      }

      rownames(peak_selected) <- features$Feature_name
    }else
    {
      peak_selected <- NULL
    }

    rownames(features_intensity) <- features$Feature_name
    rownames(Ioncount_feature_sample_matrix) <- features$Feature_name
    rownames(feature_with_background_intensity) <- features$Feature_name
    rownames(Ioncount_feature_with_background_intensity) <- features$Feature_name

    setwd(path_to_features)
    dir.create("Temporary_files")
    setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
    print("Save peak detection and selection results")
    save(features,
         features_intensity,
         Ioncount_feature_sample_matrix,
         feature_with_background_intensity,
         Ioncount_feature_with_background_intensity,
         peak_selected,
         peak_quant,
         QC_data,
         path_to_features,
         path_to_MaxQ_output,
         samples,
         delta_rt,
         delta_mz,
         peak_min_ion_count,
         num_peaks_store,
         RT_calibration,
         mz_calibration,
         abundance_estimation_correction,
         Quant_pVal_cut,
         alignment_variability_score_cutoff,
         alignment_scores_cutoff,
         mono_iso_alignment_cutoff,
         calc_peptide_LFQ,
         calc_protein_LFQ,
         MassSpec_mode,
         multiplicity,
         file = "Quantification_raw_results.RData")

    crap <- gc(F)
  }

  ###read previously generated outputs
  #load("Quantification_raw_results.RData")

  setwd(path_to_features)

  grDevices::pdf("Temporary_files/Alignment and quantification scores.pdf")

  features$target_decoy <- ifelse(grepl("_d",features$Feature_name),"decoy","target")
  ###Tag decoy features which overlap with real features
  remove_decoy_outlier <- vector(mode="logical",length=length(which(features$target_decoy == "decoy")))
  features_target <- features[which(features$target_decoy == "target"),]
  decoy_indices <- which(features$target_decoy == "decoy" & !grepl("_d_i",features$Feature_name))

  max <- length(decoy_indices)
  pb <- tcltk::tkProgressBar(title = "Detect target-decoy overlapping features",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0
  for(i in 1:length(decoy_indices))
  {
    index <- decoy_indices[i]
    overlap <- which(features_target$m.z_range_min > features$m.z_range_min[index] & features_target$m.z_range_min < features$m.z_range_max[index] |
                       features_target$m.z_range_max > features$m.z_range_min[index] & features_target$m.z_range_max < features$m.z_range_max[index])

    overlap <- which(features_target$RT_range_min[overlap] > features$RT_range_min[index] & features_target$RT_range_min[overlap] < features$RT_range_max[index] |
                       features_target$RT_range_max[overlap] > features$RT_range_min[index] & features_target$RT_range_max[overlap] < features$RT_range_max[index])
    if(length(overlap)>0)
    {
      remove_decoy_outlier[i] <- T
    }

    updatecounter <- updatecounter + 1
    if(updatecounter >= 10)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)
  features$overlapping_decoys <- F
  features$overlapping_decoys[decoy_indices] <- remove_decoy_outlier

  ###t-test based determination of background quantifications
  ###determine if ion count significantly deviates from decoy feature ion counts
  ###perform this for signal intensities and background+signal intensities
  #Signal+Background quantification

  decoy_ion_count <- as.numeric(as.matrix(Ioncount_feature_with_background_intensity[which(features$target_decoy == "decoy" & !grepl("_d_i",features$Feature_name) & features$overlapping_decoys == F),]))
  decoy_ion_count[is.na(decoy_ion_count)] <- 0

  ###add a count of 1 for log scale
  decoy_ion_count <- decoy_ion_count + 1

  ###determine mean and sd of log2 transformed decoy ion counts for z scoring
  mean_decoy_count <- mean(base::log2(decoy_ion_count))
  sd_decoy_count <- stats::sd(base::log2(decoy_ion_count))

  ###Prepare for testing if observed ion counts of target features is significantly higher than for decoy ions (just random peak selection)
  target_ion_counts <- Ioncount_feature_with_background_intensity
  rows <- nrow(target_ion_counts)
  rownames(target_ion_counts) <- rownames(Ioncount_feature_with_background_intensity)
  target_ion_counts <- as.matrix(target_ion_counts)
  ###add 1 for log scale
  target_ion_counts <- target_ion_counts + 1
  zscores <- base::as.data.frame((base::log2(target_ion_counts)-mean_decoy_count)/sd_decoy_count)

  zscore_to_pval <- function(z, alternative="greater")
  {

    if(alternative == "greater")
    {
      pval = stats::pnorm(-z)
    }else if(alternative == "two.sided")
    {
      pval = stats::pnorm(-abs(z))
      pval = 2 * pval
    }else if(alternative=="less")
    {
      pval = stats::pnorm(z)
    }
    return(pval)
  }

  ###perform significance test
  pval_signal_with_background_quant <- base::as.data.frame(apply(zscores,1:2,zscore_to_pval))
  temp_pval_quant <- pval_signal_with_background_quant
  temp_pval_quant[is.na(temp_pval_quant)] <- 1

  QC_data[["Decoy_ion_counts"]] <- decoy_ion_count
  QC_data[["Target_ion_counts"]] <- target_ion_counts
  QC_data[["Quant_pval"]] <- pval_signal_with_background_quant

  ###plot pvalue of quantification per sample for target features
  if(ncol(temp_pval_quant) <= 15) ###only 15 samples can be plotted in one plot
  {
    graphics::par(mar=c(10,4,4,4))
    graphics::boxplot(-log10(temp_pval_quant[which(!grepl("_d",rownames(temp_pval_quant))),]),outline=F,main="Quantification pVals per target feature",ylab="pValue, -log10",names=colnames(temp_pval_quant),las=2)
    graphics::abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
  }else
  {
    graphics::par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(temp_pval_quant)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(temp_pval_quant))]
      graphics::boxplot(-log10(temp_pval_quant[which(!grepl("_d",rownames(temp_pval_quant))),columns]),outline=F,main="Quantification pVals per target feature",ylab="pValue, -log10",names=colnames(temp_pval_quant)[columns],las=2)
      graphics::abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
    }
  }

  ###plot pvalue of quantification per sample for decoy features
  if(ncol(temp_pval_quant) <= 15) ###only 15 samples can be plotted in one plot
  {
    graphics::par(mar=c(10,4,4,4))
    graphics::boxplot(-log10(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),]),outline=F,main="Quantification pVals per decoy feature",ylab="pValue, -log10",names=colnames(temp_pval_quant),las=2)
    graphics::abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
  }else
  {
    graphics::par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(temp_pval_quant)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(temp_pval_quant))]
      graphics::boxplot(-log10(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),columns]),outline=F,main="Quantification pVals per decoy feature",ylab="pValue, -log10",names=colnames(temp_pval_quant)[columns],las=2)
      graphics::abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
    }
  }

  ###plot per sample how many target feature quantifications were significant and not significant
  plot_data <- rbind(colSums(temp_pval_quant[which(!grepl("_d|_i",rownames(temp_pval_quant))),] < Quant_pVal_cut,na.rm=T),
                     colSums(temp_pval_quant[which(!grepl("_d|_i",rownames(temp_pval_quant))),] > Quant_pVal_cut,na.rm=T))
  rownames(plot_data) <- c("certain","uncertain")

  if(ncol(plot_data) <= 15) ###only 15 samples can be plotted in one plot
  {
    p <- Barplotsstacked(plot_data,AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(plot_data),main="Certainty of ion accumulation for target features",ylab="Count",ylim=c(0,max(colSums(plot_data,na.rm=T))),Legendtitle = "Accumulation")
  }else ###if more samples, split samples up
  {
    pages <- ceiling(ncol(plot_data)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(plot_data))]
      p <- Barplotsstacked(plot_data[,columns],AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(plot_data),main="Certainty of ion accumulation for target features",ylab="Count",ylim=c(0,max(colSums(plot_data,na.rm=T))),Legendtitle = "Accumulation")
    }
  }

  ###plot per sample how many decoy feature quantifications were significant and not significant
  plot_data <- rbind(colSums(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),] < Quant_pVal_cut,na.rm=T),
                     colSums(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),] > Quant_pVal_cut,na.rm=T))
  rownames(plot_data) <- c("certain","uncertain")

  if(ncol(plot_data) <= 15) ###only 15 samples can be plotted in one plot
  {
    p <- Barplotsstacked(plot_data,AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(plot_data),main="Certainty of ion accumulation for decoy features",ylab="Count",ylim=c(0,max(colSums(plot_data,na.rm=T))),Legendtitle = "Accumulation")
  }else ###if more samples, split samples up
  {
    pages <- ceiling(ncol(plot_data)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(plot_data))]
      p <- Barplotsstacked(plot_data[,columns],AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(plot_data),main="Certainty of ion accumulation for decoy features",ylab="Count",ylim=c(0,max(colSums(plot_data,na.rm=T))),Legendtitle = "Accumulation")
    }
  }


  ###Boxplot per sample ratio of feature signal/background
  features_background_intensity <- log10(10^feature_with_background_intensity-10^features_intensity)
  S2B <- base::log2(10^features_intensity/10^features_background_intensity)
  S2B[is.infinite(as.matrix(S2B))] <- NA
  S2B <- base::as.data.frame(S2B)
  rownames(S2B) <- features$Feature_name

  QC_data[["Signal_to_background_target_decoy"]] <- list(Target_Decoy=ifelse(grepl("_d",rownames(S2B)),"decoy","target"),
                                                         S2B=S2B)

  if(ncol(S2B) <= 15) ###only 15 samples can be plotted in one plot
  {
    graphics::par(mar=c(10,4,4,4))
    graphics::boxplot(S2B[which(!grepl("_d",rownames(S2B))),],outline=F,main="Signal to background ratio per feature quantification",ylab="Signal/Background, log2",names=colnames(S2B),las=2)
    graphics::abline(h=0,lty=2,col="red")
    for(co in 1:ncol(S2B))
    {
      med <- stats::median(S2B[which(!grepl("_d",rownames(S2B))),co],na.rm=T)
      graphics::text(co,med+((graphics::par("usr")[4]-graphics::par("usr")[3])*0.05),round(med,digits=1))
    }

  }else
  {
    graphics::par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(S2B)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(S2B))]
      graphics::boxplot(S2B[which(!grepl("_d",rownames(S2B))),columns],outline=F,main="Signal to background ratio per feature quantification",ylab="Signal/Background, log2",names=colnames(S2B)[columns],las=2)
      graphics::abline(h=0,lty=2,col="red")
      for(co in 1:length(columns))
      {
        med <- stats::median(S2B[which(!grepl("_d",rownames(S2B))),columns[co]],na.rm=T)
        graphics::text(co,med+((graphics::par("usr")[4]-graphics::par("usr")[3])*0.05),round(med,digits=1))
      }
    }
  }

  pval_signal_with_background_quant <- base::as.data.frame(pval_signal_with_background_quant)
  features_intensity <- as.matrix(features_intensity)
  feature_with_background_intensity <- as.matrix(feature_with_background_intensity)
  features_intensity[is.infinite(features_intensity)] <- NA
  feature_with_background_intensity[is.infinite(feature_with_background_intensity)] <- NA
  features_intensity <- base::as.data.frame(features_intensity)
  feature_with_background_intensity <- base::as.data.frame(feature_with_background_intensity)

  ###calculate some scores which give insight in confidence of quantifications (selection of right peaks)

  #arrange selected peaks results in a long table

  peaks <- base::data.frame(RT=unlist(sapply(peak_quant, function(x) x$Peak_rt_with_background)),
                      mz=unlist(sapply(peak_quant, function(x) x$Peak_mz_with_background)),
                      known=unlist(sapply(peak_quant, function(x) x$correct_peak_with_background)),
                      ion_count=unlist(sapply(peak_quant, function(x) x$Ioncount_feature_with_background_intensity)),
                      intensity=unlist(sapply(peak_quant, function(x) x$feature_with_background_intensity)))
  peaks$sample <- rep(sort(rep(samples,nrow(features))),length(peak_quant))
  peaks$feature <- rep(features$Feature_name,length(peak_quant)*length(samples))
  peaks$peak <- sort(rep(1:length(peak_quant),nrow(features)*length(samples)))
  peaks$known[is.na(peaks$known)] <- 0

  #add information which peak was finally used
  peaks$selected <- NA
  for(s in samples)
  {
    selection_sample <- which(peaks$sample == s)
    matching_indices <- match(peaks$feature[selection_sample],rownames(peak_selected))
    peaks$selected[selection_sample] <- peak_selected[matching_indices,s]
  }

  #add corrected RT and mz per entry (eliminate expected deviations between samples)
  peaks$RT_correct <- NA
  peaks$mz_correct <- NA
  peaks$RT_correct <- as.numeric(peaks$RT_correct)
  peaks$mz_correct <- as.numeric(peaks$mz_correct)
  for(s in 1:length(samples))
  {
    corrections <- features[,c(1,which(grepl(base::paste("_calibration\\.",base::gsub("-",".",samples[s]),"$",sep=""),colnames(features))))]
    selection_peaks <- which(peaks$sample == samples[s])
    RT_corrected <- peaks$RT[selection_peaks]-corrections[match(peaks$feature[selection_peaks],corrections$Feature_name),2]
    mz_corrected <- peaks$mz[selection_peaks]-corrections[match(peaks$feature[selection_peaks],corrections$Feature_name),3]
    data.table::set(peaks,as.integer(selection_peaks),c(10L,11L),value=list(RT_corrected,mz_corrected))
  }


  #calculate general peak variability score using delta_rt and delta_mz as sd for zscoring
  ###feature with significant score here should be potentially excluded
  Variability_alignment_scoring <- function(features,peaks,delta_rt,delta_mz,alignment_variability_score_cutoff=0.05,plot=T)
  {
    #check in which quant mode ... LFQ or SILAC
    if(any(grepl("_Channel_light|_Channel_medium|_Channel_heavy",peaks$sample)))
    {
      if(any(grepl("_Channel_medium",peaks$sample)))
      {
        multiplicity <- 3
      }else
      {
        multiplicity <- 2
      }
    }else
    {
      multiplicity <- 1
    }

    ###Convert z scores (vector or matrix) to pvalues
    zscore_to_pval <- function(z, alternative="greater")
    {

      if(alternative == "greater")
      {
        pval = stats::pnorm(-z)
      }else if(alternative == "two.sided")
      {
        pval = stats::pnorm(-abs(z))
        pval = 2 * pval
      }else if(alternative=="less")
      {
        pval = stats::pnorm(z)
      }
      return(pval)
    }
    #select peaks which were selected by peak selection
    temp <- peaks[which(peaks$peak == peaks$selected),]
    #calculate interquartile range per feature over samples for RT and mz
    temp_RT <- stats::aggregate(temp$RT,by=list(Feature=temp$feature),FUN=matrixStats::iqr,na.rm=T)
    if(multiplicity == 1)temp_mz <- stats::aggregate(temp$mz,by=list(Feature=temp$feature),FUN=matrixStats::iqr,na.rm=T)
    if(multiplicity > 1)
    {
      #correct observed mz by expected shift
      mz_temp_corrfactors <- list("light" = features$m.z._shift_light[match(temp$feature,features$Feature_name)],
                                  "medium" = features$m.z._shift_medium[match(temp$feature,features$Feature_name)],
                                  "heavy" = features$m.z._shift_heavy[match(temp$feature,features$Feature_name)])
      mz_temp_corrfactor <- ifelse(grepl("_Channel_light",temp$sample),mz_temp_corrfactors$light,
                                   ifelse(grepl("_Channel_medium",temp$sample),mz_temp_corrfactors$medium,mz_temp_corrfactors$heavy))
      temp$mz <- temp$mz - mz_temp_corrfactor
      temp_mz <- stats::aggregate(temp$mz,by=list(Feature=temp$feature),FUN=matrixStats::iqr,na.rm=T)
    }

    #z score using delta RT and delta mz as sd
    temp_RT$x <- (temp_RT$x-mean(temp_RT$x,na.rm=T))/delta_rt
    temp_mz$x <- (temp_mz$x-mean(temp_mz$x,na.rm=T))/delta_mz
    ##features with a negative z-score represent features with much lower deviation between samples than expected
    ##features with a positive z-score represent features with much higher deviation between samples than expected --> have to be inspected and eventually removed
    #graphics::boxplot(temp_mz$x,main="General mz-deviation of selected peaks between samples",ylab="Z-score - mz deviation",ylim=c(-4,4))
    #graphics::boxplot(temp_RT$x,main="General RT-deviation of selected peaks between samples",ylab="Z-score - RT deviation",ylim=c(-4,4))

    alignment_variability_score <- base::data.frame(RT_variability_pval=(zscore_to_pval(temp_RT$x)),
                                              mz_variability_pval=(zscore_to_pval(temp_mz$x)))

    alignment_variability_score <- alignment_variability_score[match(features$Feature_name,temp_RT$Feature),]
    rownames(alignment_variability_score) <- features$Feature_name
    alignment_variability_score <- base::as.data.frame(alignment_variability_score)
    alignment_variability_score <- cbind(alignment_variability_score,
                                         matrixStats::rowMins(as.matrix(alignment_variability_score),na.rm=T))
    colnames(alignment_variability_score)[3] <- "combined_variability_pval"

    #add some plots showing number of features with general high variability between samples
    if(plot == T)
    {
      temp <- matrixStats::rowMins(as.matrix(alignment_variability_score[which(!grepl("_d|_i",rownames(alignment_variability_score))),]),na.rm=T)
      Variability_scores <- cbind(length(which(temp > alignment_variability_score_cutoff)),
                                  length(which(temp < alignment_variability_score_cutoff)))
      colnames(Variability_scores) <- c("normal","high")

      p <- BarplotsSBS(Variability_scores,AvgLine = F,col=c("chocolate","grey"),shownumbers = T,main="Variability in peak selection between samples",ylab="Count")

    }

    return(alignment_variability_score)
  }

  alignment_variability_score <- Variability_alignment_scoring(features,peaks,delta_rt,delta_mz,alignment_variability_score_cutoff)

  #calculate standardized (z score) deviation of peak RT/mz from mean over all samples with used sd defined by delta_RT and delta_mz
  #z scores are then translated into two sided pvalues and per feature/sample the minimum of pval of RT or mz is reported
  Alignment_scoring <- function(peaks,features_select,samples,sd_RT=0.5,sd_mz=0.001,corrected_alignment=T,alignment_scores_cutoff=0.05,plot=T)
  {
    ####determine for every sample if detected peak RT or mz is significantly deviating from all other determined peak RT or mz
    #extract only finally used peaks
    temp_peaks <- peaks[which(peaks$peak == peaks$selected),]
    #order according to feature name
    temp_peaks <- temp_peaks[order(temp_peaks$feature,temp_peaks$sample),]

    #check if we are looking at raw mz and if multiplicity != 1
    if(corrected_alignment == F)
    {
      #check in which quant mode ... LFQ or SILAC
      if(any(grepl("_Channel_light|_Channel_medium|_Channel_heavy",peaks$sample)))
      {
        if(any(grepl("_Channel_medium",peaks$sample)))
        {
          multiplicity <- 3
        }else
        {
          multiplicity <- 2
        }
      }else
      {
        multiplicity <- 1
      }

      #if multiplicity > 1 then raw mz have to be corrected for expected mz shifts
      if(multiplicity > 1)
      {
        #correct observed mz by expected shift
        mz_temp_corrfactors <- list("light" = features$m.z._shift_light[match(temp_peaks$feature,features$Feature_name)],
                                    "medium" = features$m.z._shift_medium[match(temp_peaks$feature,features$Feature_name)],
                                    "heavy" = features$m.z._shift_heavy[match(temp_peaks$feature,features$Feature_name)])
        mz_temp_corrfactor <- ifelse(grepl("_Channel_light",temp_peaks$sample),mz_temp_corrfactors$light,
                                     ifelse(grepl("_Channel_medium",temp_peaks$sample),mz_temp_corrfactors$medium,mz_temp_corrfactors$heavy))
        temp_peaks$mz <- temp_peaks$mz - mz_temp_corrfactor
      }
    }

    #prepare matrices which should store RT and mz of selected peaks per sample and feature
    if(corrected_alignment == T)selected_peaks_RT <- stats::reshape(temp_peaks[,c("feature","sample","RT_correct")], idvar = "feature", timevar = "sample", direction = "wide")
    if(corrected_alignment == F)selected_peaks_RT <- stats::reshape(temp_peaks[,c("feature","sample","RT")], idvar = "feature", timevar = "sample", direction = "wide")
    rownames(selected_peaks_RT) <- selected_peaks_RT$feature
    selected_peaks_RT <- selected_peaks_RT[,-1]
    if(corrected_alignment == T)selected_peaks_mz <- stats::reshape(temp_peaks[,c("feature","sample","mz_correct")], idvar = "feature", timevar = "sample", direction = "wide")
    if(corrected_alignment == F)selected_peaks_mz <- stats::reshape(temp_peaks[,c("feature","sample","mz")], idvar = "feature", timevar = "sample", direction = "wide")
    rownames(selected_peaks_mz) <- selected_peaks_mz$feature
    selected_peaks_mz <- selected_peaks_mz[,-1]

    #calculate mean and sd of RT and mz per feature
    selected_peaks_distributions_RT <- base::data.frame(mean=rowMeans(as.matrix(selected_peaks_RT),na.rm=T),
                                                  sd=matrixStats::rowSds(as.matrix(selected_peaks_RT,na.rm=T)))
    selected_peaks_distributions_mz <- base::data.frame(mean=rowMeans(as.matrix(selected_peaks_mz),na.rm=T),
                                                  sd=matrixStats::rowSds(as.matrix(selected_peaks_mz,na.rm=T)))
    #convert observed RTs and mzs per feature and sample into z-scores
    selected_peaks_RT_zscore <- (selected_peaks_RT-selected_peaks_distributions_RT$mean)/sd_RT
    selected_peaks_mz_zscore <- (selected_peaks_mz-selected_peaks_distributions_mz$mean)/sd_mz
    #convert zscores into pvalues
    zscore_to_pval <- function(z)
    {
      pval = stats::pnorm(-abs(z))
      pval = 2 * pval
      return(pval)
    }
    selected_peaks_RT_pval <- base::as.data.frame(apply(selected_peaks_RT_zscore,1:2,zscore_to_pval))
    selected_peaks_mz_pval <- base::as.data.frame(apply(selected_peaks_mz_zscore,1:2,zscore_to_pval))
    if(corrected_alignment == T)colnames(selected_peaks_RT_pval) <- base::gsub("RT_correct.","",colnames(selected_peaks_RT_pval))
    if(corrected_alignment == T)colnames(selected_peaks_mz_pval) <- base::gsub("mz_correct.","",colnames(selected_peaks_mz_pval))
    if(corrected_alignment == F)colnames(selected_peaks_RT_pval) <- base::gsub("RT.","",colnames(selected_peaks_RT_pval))
    if(corrected_alignment == F)colnames(selected_peaks_mz_pval) <- base::gsub("mz.","",colnames(selected_peaks_mz_pval))
    #reorder rows to match to ordering in features_select
    selected_peaks_RT_pval <- selected_peaks_RT_pval[match(features_select$Feature_name,rownames(selected_peaks_RT_pval)),]
    selected_peaks_mz_pval <- selected_peaks_mz_pval[match(features_select$Feature_name,rownames(selected_peaks_mz_pval)),]
    #summarize alignment score per feature and sample to single minimal pval (either RT or mz)
    alignment_score_peaks <- pmin(selected_peaks_RT_pval,selected_peaks_mz_pval)

    ###add some plots about for how many features we are certain/uncertain about the peak selection
    if(plot == T)
    {
      temp_score <- alignment_score_peaks[which(!grepl("_d|_d_i",rownames(alignment_score_peaks))),]
      temp_score <- rbind(colSums(temp_score > alignment_scores_cutoff,na.rm=T),
                          colSums(temp_score < alignment_scores_cutoff,na.rm=T))
      rownames(temp_score) <- c("certain","uncertain")

      main <- ifelse(corrected_alignment == T,"Certainty of peak selection (corrected)","Certainty of peak selection (raw)")

      if(ncol(temp_score) <= 15) ###only 15 samples can be plotted in one plot
      {
        p <- Barplotsstacked(temp_score,AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(temp_score),main=main,ylab="Count",ylim=c(0,max(colSums(temp_score,na.rm=T))),Legendtitle = "Peak selection")
      }else ###if more samples, split samples up
      {
        pages <- ceiling(ncol(temp_score)/15)
        for(p in 1:pages)
        {
          columns <- (((p-1)*15)+1):(p*15)
          columns <- columns[which(columns <= ncol(temp_score))]
          p <- Barplotsstacked(temp_score[,columns],AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(temp_score),main=main,ylab="Count",ylim=c(0,max(colSums(temp_score,na.rm=T))),Legendtitle = "Peak selection")
        }
      }
    }

    return(alignment_score_peaks)
  }
  #perform alignment scoring based on corrected mz and RT per feature/sample
  alignment_scores_peaks_correct <- Alignment_scoring(peaks,features,samples,sd_RT = delta_rt,sd_mz = delta_mz,corrected_alignment = T,alignment_scores_cutoff=alignment_scores_cutoff)
  #perform alignment scoring based on raw mz and RT per feature/sample to detect unexpected complete outliers
  alignment_scores_peaks_raw <- Alignment_scoring(peaks,features,samples,sd_RT = delta_rt,sd_mz = delta_mz,corrected_alignment = F,alignment_scores_cutoff=alignment_scores_cutoff)

  ###Calculate score based on cooccurence of isotope peaks
  if(any(grepl("_i",features$Feature_name) & !grepl("_d_i",features$Feature_name)))
  {
    Mono_iso_alignment_scoring <- function(features,samples,peaks,pval_signal_with_background_quant,delta_rt,delta_mz,mono_iso_alignment_cutoff,plot=T)
    {
      #prepare dataframe into which all results are saved
      mono_iso_alignment_summary <- base::as.data.frame(matrix(nrow=nrow(features),ncol=length(samples)))
      rownames(mono_iso_alignment_summary) <- features$Feature_name
      for(c in 1:ncol(mono_iso_alignment_summary))
      {
        mono_iso_alignment_summary[,c] <- as.numeric(mono_iso_alignment_summary[,c])
      }

      if(length(which(peaks$known == 1)) >= 10)
      {
        zscore_to_pval <- function(z, alternative="greater")
        {

          if(alternative == "greater")
          {
            pval = stats::pnorm(-z)
          }else if(alternative == "two.sided")
          {
            pval = stats::pnorm(-abs(z))
            pval = 2 * pval
          }else if(alternative=="less")
          {
            pval = stats::pnorm(z)
          }
          return(pval)
        }

        ###estimate expected distributions of deviations in RT and mz between M and M+1 peaks
        mean_RT_distribution <- vector("numeric",length(samples))
        sd_RT_distribution <- vector("numeric",length(samples))
        mean_mz_distribution <- vector("numeric",length(samples))
        sd_mz_distribution <- vector("numeric",length(samples))
        for(s in 1:length(samples))
        {
          temp <- peaks[which(peaks$sample == samples[s]),]
          #add quantification scores
          temp$quant_score <- pval_signal_with_background_quant[match(temp$feature,rownames(pval_signal_with_background_quant)),s]
          #add charge state
          temp$charge <- features$Charge[match(features$Feature_name,temp$feature)]
          #selected peaks RT and mz
          temp_selected <- temp[which(temp$peak == temp$selected),]
          #split features into mono and +1 isotopes
          temp_selected_mono <- temp_selected[which(!grepl("_i",temp_selected$feature)),]
          temp_selected_iso <- temp_selected[which(grepl("_i",temp_selected$feature)),]
          temp_selected_iso$feature <- base::gsub("_i","",temp_selected_iso$feature)
          #combine both and only keep features for which we also have +1 isotope features
          combined <- dplyr::inner_join(temp_selected_mono,temp_selected_iso,by="feature")
          #determine expected variation between mono and +1 isotopes based on true identifications and significantly quantified
          #mono and +1 isotope features
          true <- combined[which(combined$known.x == 1 & combined$quant_score.x < 0.05 & combined$quant_score.y < 0.05),]
          true$delta_rt_iso_mono <- true$RT.y-true$RT.x
          true$delta_mz_iso_mono <- (((true$mz.y*true$charge.y)-1.002054)/true$charge.y)-true$mz.x

          mean_RT_distribution[s] <- mean(true$delta_rt_iso_mono,na.rm=T)
          sd_RT_distribution[s] <- stats::sd(true$delta_rt_iso_mono,na.rm=T)
          mean_mz_distribution[s] <- mean(true$delta_mz_iso_mono,na.rm=T)
          sd_mz_distribution[s] <- stats::sd(true$delta_mz_iso_mono,na.rm=T)
        }
        #estimate normal distribution of deviations in RT and mz between mono and +1 isotope peak
        RT_deviation_norm <- base::data.frame(x=seq(-delta_rt*3, delta_rt*3, length=1000),
                                        y=stats::dnorm(seq(-delta_rt*3, delta_rt*3, length=1000), mean=mean(mean_RT_distribution,na.rm=T), sd=mean(sd_RT_distribution,na.rm=T)))
        mz_deviation_norm <- base::data.frame(x=seq(-delta_mz*3,delta_mz*3, length=1000),
                                        y=stats::dnorm(seq(-delta_mz*3,delta_mz*3, length=1000), mean=mean(mean_mz_distribution,na.rm=T), sd=mean(sd_mz_distribution,na.rm=T)))
        graphics::plot(RT_deviation_norm$x,RT_deviation_norm$y,type="l",xlab="RT(M+1) - RT(M)",main="RT-deviation M vs M+1 peaks - true quantifications",ylab="Density")
        graphics::abline(v=mean(true$delta_rt_iso_mono),lty=2)
        graphics::plot(mz_deviation_norm$x,mz_deviation_norm$y,type="l",xlab="mz(M+1) - mz(M)",main="mz-deviation M vs M+1 peaks - true quantifications",ylab="Density")
        graphics::abline(v=mean(true$delta_mz_iso_mono),lty=2)
        #now calculate for every available feature for which an isotope feature is available how well both peaks are aligned
        mean_RT_distribution_mean <- mean(mean_RT_distribution,na.rm=T)
        sd_RT_distribution_mean <- mean(sd_RT_distribution,na.rm=T)
        mean_mz_distribution_mean <- mean(mean_mz_distribution,na.rm=T)
        sd_mz_distribution_mean <- mean(sd_mz_distribution,na.rm=T)
        for(s in 1:length(samples))
        {
          temp <- peaks[which(peaks$sample == samples[s]),]
          #add quantification scores
          temp$quant_score <- pval_signal_with_background_quant[match(temp$feature,rownames(pval_signal_with_background_quant)),s]
          #add charge state
          temp$charge <- features$Charge[match(features$Feature_name,temp$feature)]
          #selected peaks RT and mz
          temp_selected <- temp[which(temp$peak == temp$selected),]
          #split features into mono and +1 isotopes
          temp_selected_mono <- temp_selected[which(!grepl("_i",temp_selected$feature)),]
          temp_selected_iso <- temp_selected[which(grepl("_i",temp_selected$feature)),]
          temp_selected_iso$feature <- base::gsub("_i","",temp_selected_iso$feature)
          #combine both and only keep features for which we also have +1 isotope features
          combined <- dplyr::inner_join(temp_selected_mono,temp_selected_iso,by="feature")
          combined$delta_rt_iso_mono <- combined$RT.y-combined$RT.x
          combined$delta_mz_iso_mono <- (((combined$mz.y*combined$charge.y)-1.002054)/combined$charge.y)-combined$mz.x
          #convert to z-score
          combined$z_delta_rt_iso_mono <- (combined$delta_rt_iso_mono-mean_RT_distribution_mean)/sd_RT_distribution_mean
          combined$z_delta_mz_iso_mono <- (combined$delta_mz_iso_mono-mean_mz_distribution_mean)/sd_mz_distribution_mean
          #convert to pvalues
          combined$RT_deviation_pval <- zscore_to_pval(combined$z_delta_rt_iso_mono,"two.sided")
          combined$mz_deviation_pval <- zscore_to_pval(combined$z_delta_mz_iso_mono,"two.sided")
          #store results
          combined$score <- matrixStats::rowMins(as.matrix(combined[,c("RT_deviation_pval","mz_deviation_pval")]),na.rm=T)

          mono_iso_alignment_summary[,s] <- combined$score[match(base::gsub("_i","",rownames(mono_iso_alignment_summary)),combined$feature)]
          colnames(mono_iso_alignment_summary)[s] <- samples[s]
        }

        mono_iso_alignment_summary[is.infinite(as.matrix(mono_iso_alignment_summary))] <- NA

        ###add some plots about for how many features we are certain about peak selection based on detected deviations to isotope peak selection
        if(plot == T)
        {
          temp_alignmentscores <- mono_iso_alignment_summary
          temp_alignmentscores <- temp_alignmentscores[which(!grepl("_d|_i|_d_i",rownames(temp_alignmentscores))),]
          RT_alignment_scores <- rbind(colSums(temp_alignmentscores > mono_iso_alignment_cutoff,na.rm=T),
                                       colSums(temp_alignmentscores < mono_iso_alignment_cutoff,na.rm=T))
          rownames(RT_alignment_scores) <- c("certain","uncertain")

          if(ncol(RT_alignment_scores) <= 15) ###only 15 samples can be plotted in one plot
          {
            p <- Barplotsstacked(RT_alignment_scores,AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(RT_alignment_scores),main="Certainty of peak selection based on mono/+1 isotope peak selection",ylab="Count",ylim=c(0,max(colSums(RT_alignment_scores,na.rm=T))),Legendtitle = "Peak selection")
          }else ###if more samples, split samples up
          {
            pages <- ceiling(ncol(RT_alignment_scores)/15)
            for(p in 1:pages)
            {
              columns <- (((p-1)*15)+1):(p*15)
              columns <- columns[which(columns <= ncol(RT_alignment_scores))]
              p <- Barplotsstacked(RT_alignment_scores[,columns],AvgLine = F,col=c("chocolate","grey"),shownumbers_total = F,shownumbers = T,Legends = rownames(RT_alignment_scores),main="Certainty of peak selection based on mono/+1 isotope peak selection",ylab="Count",ylim=c(0,max(colSums(RT_alignment_scores,na.rm=T))),Legendtitle = "Peak selection")
            }
          }
        }
      }else
      {
        print("Warning: The true peak for too few features is known. Skipping +1-isotope alignment scoring.")
      }

      return(mono_iso_alignment_summary)
    }

    mono_iso_alignment_summary <- Mono_iso_alignment_scoring(features,samples,peaks,pval_signal_with_background_quant,delta_rt,delta_mz,mono_iso_alignment_cutoff)
  }else
  {
    mono_iso_alignment_summary <- NA
  }

  FDR_peak_selection <- Peak_selection_FDR(num_features=500,features,samples,peaks,path_to_features,peak_quant,feature_with_background_intensity,peak_selected,delta_mz,delta_rt,peak_min_ion_count,num_peaks_store,alignment_scores_cutoff,n_cores,peak_decision,Alignment_scoring,plot=T)

  QC_data$FDR_peak_selection <- FDR_peak_selection

  if(any(!is.na(FDR_peak_selection$Large_Intensity_delta_FDR)))
  {
    for(i in 1:length(FDR_peak_selection$Total_FDR))
    {
      print(base::paste(samples[i]," - FDR of peak selection: ",round(FDR_peak_selection$Large_Intensity_delta_FDR[i],digits=1)," %",sep=""))
    }
  }

  grDevices::dev.off()

  ###function to plot peak selections
  # Plot_feature_quantification <- function(path_to_graph_data,features,selected_feature,samples,peaks,num_peaks_store)
  # {
  #   e <- new.env()
  #
  #   max <- length(samples)
  #   pb <- tcltk::tkProgressBar(title = "Plot peak selections per sample",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  #   start_time <- Sys.time()
  #   updatecounter <- 0
  #   time_require <- 0
  #
  #   for(s in 1:length(samples))
  #   {
  #     load(base::paste(path_to_graph_data,"\\",samples[s],"_feature_graphs.RData",sep=""),envir=e)
  #
  #     max2 <- length(selected_feature)
  #     pb2 <- tcltk::tkProgressBar(title = "Plot peak selections",label=base::paste( round(0/max2*100, 0),"% done"), min = 0,max = max2, width = 300)
  #     start_time2 <- Sys.time()
  #     updatecounter2 <- 0
  #     time_require2 <- 0
  #     for(f in selected_feature)
  #     {
  #       i <- which(selected_feature == f)
  #       temp <- e$peaks_graph[[as.character(f)]]
  #       if(!is.na(temp))
  #       {
  #         ###plot KDE
  #         image(temp$kdemap,main=base::paste(f,"-",samples[s]),xlab="RT",ylab="m/z",xlim=c(min(temp$kdemap$x),max(temp$kdemap$x)),ylim=c(min(temp$kdemap$y),max(temp$kdemap$y)))
  #         ###indicate where peaks were detected
  #         graphics::text(temp$maxima_select$RT,temp$maxima_select$mz,1:nrow(temp$maxima_select),col=ifelse(temp$maxima_select$known_peak == 1,"green","black"))
  #         ###expected window
  #         sel <- which(features$Feature_name == f)
  #         RT_expected <- features$RT[sel] + features[sel,which(grepl("RT_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))]
  #         RT_window_width <- features$RT_length[sel]
  #         cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))
  #         mz_expected <- features$m.z[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))]
  #         mz_window <- c(features$m.z_range_min[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))],
  #                        features$m.z_range_max[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))])
  #         delta_mz <- (mz_window[2]-mz_window[1])/2
  #
  #         graphics::rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2)
  #         graphics::points(RT_expected,mz_expected,pch=4)
  #         ###indicate which peak was selected
  #         temp_selected <- peaks[which(peaks$sample == samples[s] & peaks$feature == f & peaks$peak == peaks$selected),]
  #
  #         if(temp_selected$peak == (num_peaks_store+1)) #selected standard window
  #         {
  #           graphics::rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2,border="blue")
  #           intensity_temp <- temp$peak_intensities[nrow(temp$peak_intensities),1]
  #           graphics::text(RT_expected,mz_expected-(1.1*delta_mz),round(intensity_temp,2),col="blue")
  #         }else
  #         {
  #           rect_pos_x_1 <- temp$maxima_select$RT[temp_selected$peak]-RT_window_width/2
  #           rect_pos_y_1 <- temp$maxima_select$mz[temp_selected$peak]-delta_mz
  #           rect_pos_x_2 <- temp$maxima_select$RT[temp_selected$peak]+RT_window_width/2
  #           rect_pos_y_2 <- temp$maxima_select$mz[temp_selected$peak]+delta_mz
  #
  #
  #           graphics::rect(rect_pos_x_1,rect_pos_y_1,rect_pos_x_2,rect_pos_y_2,lty=2,border="blue")
  #           intensity_temp <- temp$peak_intensities[temp_selected$peak,1]
  #           graphics::text(temp_selected$RT,temp_selected$mz-(1.1*delta_mz),round(intensity_temp,2),col="blue")
  #         }
  #       }
  #
  #       updatecounter2 <- updatecounter2 + 1
  #       if(updatecounter2 >= 1)
  #       {
  #         time_elapsed <- difftime(Sys.time(),start_time2,units="secs")
  #         time_require2 <- (time_elapsed/(i/max2))*(1-(i/max2))
  #         td <- lubridate::seconds_to_period(time_require2)
  #         time_require2 <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))
  #
  #         updatecounter <- 0
  #         tcltk::setTkProgressBar(pb2, i, label=base::paste( round(i/max2*100, 0)," % done (",i,"/",max2,", Time require: ",time_require2,")",sep = ""))
  #       }
  #     }
  #     close(pb2)
  #     updatecounter <- updatecounter + 1
  #     if(updatecounter >= 1)
  #     {
  #       time_elapsed <- difftime(Sys.time(),start_time,units="secs")
  #       time_require <- (time_elapsed/(s/max))*(1-(s/max))
  #       td <- lubridate::seconds_to_period(time_require)
  #       time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))
  #
  #       updatecounter <- 0
  #       tcltk::setTkProgressBar(pb, s, label=base::paste( round(s/max*100, 0)," % done (",s,"/",max,", Time require: ",time_require,")",sep = ""))
  #     }
  #   }
  #   close(pb)
  # }
  #
  # f <- c("Feature_16094","Feature_27778")
  # View(peaks[which(peaks$feature%in%f & peaks$peak == peaks$selected),])
  # features_select[which(features_select$Feature_name == f),]
  # Plot_feature_quantification(path_to_graph_data=base::paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""),
  #                             features=features,
  #                             selected_feature=f,
  #                             samples=samples[c(13:18)],
  #                             peaks=peaks,
  #                             num_peaks_store=num_peaks_store)


  ###evaluate alignment performance
  #after peak alignment
  peaks_selected <- peaks[which(peaks$peak == peaks$selected),]
  #raw peaks
  peaks_raw <- peaks[which(peaks$peak == 6),]
  select_only_known <- which(base::paste(peaks_raw$sample,peaks_raw$feature) %in% base::paste(peaks_selected$sample,peaks_selected$feature) )
  peaks_raw <- peaks_raw[select_only_known,]

  if(multiplicity > 1)
  {
    #correct observed mz by expected shift
    mz_temp_corrfactors <- list("light" = features$m.z._shift_light[match(peaks_raw$feature,features$Feature_name)],
                                "medium" = features$m.z._shift_medium[match(peaks_raw$feature,features$Feature_name)],
                                "heavy" = features$m.z._shift_heavy[match(peaks_raw$feature,features$Feature_name)])
    mz_temp_corrfactor <- ifelse(grepl("_Channel_light",peaks_raw$sample),mz_temp_corrfactors$light,
                                 ifelse(grepl("_Channel_medium",peaks_raw$sample),mz_temp_corrfactors$medium,mz_temp_corrfactors$heavy))
    peaks_raw$mz <- peaks_raw$mz - mz_temp_corrfactor


    mz_temp_corrfactors <- list("light" = features$m.z._shift_light[match(peaks_selected$feature,features$Feature_name)],
                                "medium" = features$m.z._shift_medium[match(peaks_selected$feature,features$Feature_name)],
                                "heavy" = features$m.z._shift_heavy[match(peaks_selected$feature,features$Feature_name)])
    mz_temp_corrfactor <- ifelse(grepl("_Channel_light",peaks_selected$sample),mz_temp_corrfactors$light,
                                 ifelse(grepl("_Channel_medium",peaks_selected$sample),mz_temp_corrfactors$medium,mz_temp_corrfactors$heavy))
    peaks_selected$mz <- peaks_selected$mz - mz_temp_corrfactor
  }

  #deviation in raw data
  temp_mean_raw <- stats::aggregate(peaks_raw[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_raw$feature),FUN=mean,na.rm=T)
  temp_sd_raw <- stats::aggregate(peaks_raw[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_raw$feature),FUN=stats::sd,na.rm=T)
  #deviation after alignment
  temp_mean_aligned <- stats::aggregate(peaks_selected[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_selected$feature),FUN=mean,na.rm=T)
  temp_sd_aligned <- stats::aggregate(peaks_selected[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_selected$feature),FUN=stats::sd,na.rm=T)
  #exclude features where sd in uncorrected RT > 1 min
  sel <- which(temp_sd_raw$RT > 1)
  temp_mean_raw <- temp_mean_raw[-sel,]
  temp_sd_raw <- temp_sd_raw[-sel,]
  temp_mean_aligned <- temp_mean_aligned[-sel,]
  temp_sd_aligned <- temp_sd_aligned[-sel,]
  grDevices::pdf("Performance of feature alignment.pdf",useDingbats = F)
  graphics::boxplot(temp_sd_raw$RT,temp_sd_aligned$RT_correct,outline=F,main="Variability of feature RT",names = c("Raw","Aligned"),ylab="SD of RT [min]")
  graphics::boxplot(temp_sd_raw$mz,temp_sd_aligned$mz_correct,outline=F,main="Variability of feature mz",names = c("Raw","Aligned"),ylab="SD of mz [Da]")
  grDevices::dev.off()


  ###Save temporary results
  print("Save quantification results after alignment scoring")

  temp_results <- list(features=features,
                       features_intensity=features_intensity,
                       Ioncount_feature_sample_matrix=Ioncount_feature_sample_matrix,
                       feature_with_background_intensity=feature_with_background_intensity,
                       Ioncount_feature_with_background_intensity=Ioncount_feature_with_background_intensity,
                       peak_selected=peak_selected,
                       S2B=S2B,
                       pval_signal_with_background_quant=pval_signal_with_background_quant,
                       alignment_variability_score=alignment_variability_score,
                       alignment_scores_peaks_correct=alignment_scores_peaks_correct,
                       alignment_scores_peaks_raw=alignment_scores_peaks_raw,
                       mono_iso_alignment_summary=mono_iso_alignment_summary,
                       peaks=peaks,
                       QC_data)
  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  save(temp_results,file = "Quantification_raw_results_with_scores.RData")

  crap <- gc(F)
  ###load temporary results
  # setwd(base::paste(path_to_features,"\\Temporary_files",sep=""))
  # load("Quantification_raw_results_with_scores.RData")
  # features <- temp_results$features
  # features_intensity <- temp_results$features_intensity
  # Ioncount_feature_sample_matrix <- temp_results$Ioncount_feature_sample_matrix
  # feature_with_background_intensity <- temp_results$feature_with_background_intensity
  # feature_with_background_intensity_imputed <- temp_results$feature_with_background_intensity_imputed
  # Ioncount_feature_with_background_intensity <- temp_results$Ioncount_feature_with_background_intensity
  # peak_selected <- temp_results$peak_selected
  # S2B <- temp_results$S2B
  # pval_signal_with_background_quant <- temp_results$pval_signal_with_background_quant
  # alignment_variability_score <- temp_results$alignment_variability_score
  # alignment_scores_peaks_correct <- temp_results$alignment_scores_peaks_correct
  # alignment_scores_peaks_raw <- temp_results$alignment_scores_peaks_raw
  # mono_iso_alignment_summary <- temp_results$mono_iso_alignment_summary
  # peaks <- temp_results$peaks
  ###remove decoy ions from results
  selection <- which(!grepl("_d",features$Feature_name))
  if(length(selection) > 0)
  {
    features<-features[selection,]
    Ioncount_feature_sample_matrix<-Ioncount_feature_sample_matrix[selection,]
    feature_with_background_intensity<-feature_with_background_intensity[selection,]
    Ioncount_feature_with_background_intensity<-Ioncount_feature_with_background_intensity[selection,]
    peak_selected<-peak_selected[selection,]
    S2B<-S2B[selection,]
    pval_signal_with_background_quant<-pval_signal_with_background_quant[selection,]
    alignment_variability_score<-alignment_variability_score[selection,]
    alignment_scores_peaks_correct<-alignment_scores_peaks_correct[selection,]
    alignment_scores_peaks_raw<-alignment_scores_peaks_raw[selection,]
    if(any(!is.na(mono_iso_alignment_summary)))mono_iso_alignment_summary<-mono_iso_alignment_summary[selection,]
  }

  ###remove isotope features which never show significant quantification
  selection <- which(grepl("_i",features$Feature_name) & matrixStats::rowMins(as.matrix(pval_signal_with_background_quant),na.rm=T) > Quant_pVal_cut)

  total_length <- length(which(grepl("_i",features$Feature_name)))
  if(length(selection) > 0)
  {
    features<-features[-selection,]
    Ioncount_feature_sample_matrix<-Ioncount_feature_sample_matrix[-selection,]
    feature_with_background_intensity<-feature_with_background_intensity[-selection,]
    Ioncount_feature_with_background_intensity<-Ioncount_feature_with_background_intensity[-selection,]
    peak_selected<-peak_selected[-selection,]
    S2B<-S2B[-selection,]
    pval_signal_with_background_quant<-pval_signal_with_background_quant[-selection,]
    alignment_variability_score<-alignment_variability_score[-selection,]
    alignment_scores_peaks_correct<-alignment_scores_peaks_correct[-selection,]
    alignment_scores_peaks_raw<-alignment_scores_peaks_raw[-selection,]
    if(any(!is.na(mono_iso_alignment_summary)))mono_iso_alignment_summary<-mono_iso_alignment_summary[-selection,]
    print(base::paste("Removed ",length(selection)," isotope features (",round(length(selection)/total_length*100,digits=1)," %) as they dont show significant ion accumulation in any sample.",sep=""))
  }

  ####impute missing values based on generalized additive model
  impute_feature_level_quant <- function(data,features,background_intensity_GAM_table_per_feature_sum)
  {
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    quant_data <- data.table::copy(data)
    background_intensity_GAM_table_per_feature <- background_intensity_GAM_table_per_feature_sum

    for(c in 1:ncol(quant_data))
    {
      missing_vals_rows <- which(is.na(quant_data[,c]))

      impute_vals <- background_intensity_GAM_table_per_feature[match(rownames(quant_data)[missing_vals_rows],rownames(background_intensity_GAM_table_per_feature)),c]

      impute_vals <- log10(2^impute_vals)

      data.table::set(quant_data,i=as.integer(missing_vals_rows),j=as.integer(c),value=impute_vals)
    }

    return(quant_data)
  }

  feature_with_background_intensity_imputed <- impute_feature_level_quant(data=feature_with_background_intensity,
                                                                          features=features,
                                                                          background_intensity_GAM_table_per_feature_sum=QC_data$Decoy_feature_parameters$background_intensity_GAM_table_per_feature_sum)

  ##now perform filtering of quantifications based on scores
  #remove all features with too high variability
  selection <- rownames(alignment_variability_score)[which(matrixStats::rowMins(as.matrix(alignment_variability_score),na.rm=T) < alignment_variability_score_cutoff)]
  if(length(selection)>0)
  {
    print(base::paste("Removed ",length(which(!grepl("_i",selection)))," (",length(which(grepl("_i",selection)))," isotope)"," features from quantification results due to too high variability in alignment between samples.",sep=""))

    selection <- match(selection,features$Feature_name)

    features<-features[-selection,]
    Ioncount_feature_sample_matrix<-Ioncount_feature_sample_matrix[-selection,]
    feature_with_background_intensity<-feature_with_background_intensity[-selection,]
    feature_with_background_intensity_imputed<-feature_with_background_intensity_imputed[-selection,]
    Ioncount_feature_with_background_intensity<-Ioncount_feature_with_background_intensity[-selection,]
    peak_selected<-peak_selected[-selection,]
    S2B<-S2B[-selection,]
    pval_signal_with_background_quant<-pval_signal_with_background_quant[-selection,]
    alignment_variability_score<-alignment_variability_score[-selection,]
    alignment_scores_peaks_correct<-alignment_scores_peaks_correct[-selection,]
    alignment_scores_peaks_raw<-alignment_scores_peaks_raw[-selection,]
    if(any(!is.na(mono_iso_alignment_summary)))mono_iso_alignment_summary<-mono_iso_alignment_summary[-selection,]
  }

  #remove all quantifications per sample with uncertain peak selection

  for(s in 1:length(samples))
  {
    temp_scores <- alignment_scores_peaks_correct[,which(grepl(samples[s],colnames(alignment_scores_peaks_correct))),drop=F]
    temp_scores_raw <- alignment_scores_peaks_raw[,which(grepl(samples[s],colnames(alignment_scores_peaks_correct))),drop=F]
    selection <- which(temp_scores < alignment_scores_cutoff)# | temp_scores_raw < alignment_scores_cutoff)
    names(selection) <- rownames(temp_scores)[selection]

    feature_with_background_intensity[selection,s] <- NA
    feature_with_background_intensity_imputed[selection,s] <- NA

    print(base::paste(samples[s],": Removed ",length(which(!grepl("_i",names(selection))))," (",length(which(grepl("_i",names(selection))))," isotope)"," quantifications (",round(length(selection)/nrow(alignment_scores_peaks_correct)*100,digits=1)," %) due to uncertain peak selection",sep=""))
  }

  #remove isotope quantifications for uncertain mono-+1-isos but keep mono quantification if otherwise certain of peak selection
  if(any(!is.na(mono_iso_alignment_summary)))
  {
    for(s in 1:length(samples))
    {
      if(length(which(peaks$known == 1)) >= 10)
      {
        temp_scores <- mono_iso_alignment_summary[,which(grepl(samples[s],colnames(mono_iso_alignment_summary)))]
        selection <- which(temp_scores < mono_iso_alignment_cutoff & grepl("_i",rownames(mono_iso_alignment_summary)))
        names(selection) <- rownames(mono_iso_alignment_summary)[selection]

        #how many of these still have a quantification (not already removed by previous filtering step)
        temp_quant <- feature_with_background_intensity[selection,s]

        feature_with_background_intensity[selection,s] <- NA
        feature_with_background_intensity_imputed[selection,s] <- NA

        print(base::paste(samples[s],": Removed ",length(which(!is.na(temp_quant)))," isotope quantifications (",round(length(which(!is.na(temp_quant)))/nrow(alignment_scores_peaks_correct)*100,digits=1)," %) due to discrepancy in peak selection between mono-/+1-isotope features",sep=""))
      }else ### not enough known peaks ... remove all isotopes
      {
        selection <- which(grepl("_i",rownames(feature_with_background_intensity)))
        feature_with_background_intensity[selection,s] <- NA
        feature_with_background_intensity_imputed[selection,s] <- NA

        print(base::paste(samples[s],": Removed all isotope quantifications due to too few known true peaks.",sep=""))
      }

    }
  }


  #calculate fraction of missing values before and after imputation
  total <- nrow(feature_with_background_intensity)*ncol(feature_with_background_intensity)
  missing <- length(which(is.na(as.numeric(as.matrix(feature_with_background_intensity)))))
  print(base::paste("Quantification without imputation: ",round(missing/total*100,digits=1)," % missing values",sep=""))
  missing <- length(which(is.na(as.numeric(as.matrix(feature_with_background_intensity_imputed)))))
  print(base::paste("Quantification with imputation: ",round(missing/total*100,digits=1)," % missing values",sep=""))


  ####If features contain PMPs (potentially missed peptides) then evaluate now which of them are correlating (spearman, rank) in abundance to the other peptides of the same proteins
  # if(any(grepl("_pmp",features$Feature_name)))
  # {
  #   cor_p_cut <- 0.05
  #   selection <- which(grepl("_pmp",features$Feature_name))
  #   PMPs <- list(features=features[selection,],
  #                features_intensity=features_intensity[selection,],
  #                Ioncount_feature_sample_matrix=Ioncount_feature_sample_matrix[selection,],
  #                feature_with_background_intensity=feature_with_background_intensity[selection,],
  #                Ioncount_feature_with_background_intensity=Ioncount_feature_with_background_intensity[selection,],
  #                Alignment_scores=Alignment_scores[selection,],
  #                Quant_pvals_signal_with_background=Quant_pvals_signal_with_background[selection,])
  #
  #   ###go through all proteins with PMPs and evaluate if PMPs are correlating with known peptides of the respective protein
  #   remove_list <- vector(mode="numeric",length(PMPs$features$Feature_name))
  #   remove_list[remove_list==0]<-NA
  #   count_remove <- 0
  #
  #   unique_PMP_proteins <- unique(PMPs$features$Protein)
  #
  #   max <- length(unique_PMP_proteins)
  #   pb <- tcltk::tkProgressBar(title = "Detect significantly correlating PMPs",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  #   start_time <- Sys.time()
  #   updatecounter <- 0
  #   time_require <- 0
  #   pbi <- 0
  #   for(p in unique_PMP_proteins)
  #   {
  #     selection <- which(features$Protein == p & !grepl("_i",features$Feature_name))
  #     temp <- features[selection,]
  #     temp_quant <- feature_with_background_intensity_imputed[selection,]
  #
  #     pmp_features <- which(grepl("_pmp",temp$Feature_name))
  #
  #     correlation <- rcorr(t(temp_quant),type="pearson")
  #
  #     correlation$P[correlation$P==0]<-NA
  #
  #     remove <- vector(mode="logical",length(pmp_features))
  #     for(i in pmp_features)
  #     {
  #       if(!any(correlation$P[i,-pmp_features] < cor_p_cut,na.rm=T))
  #       {
  #         remove[which(pmp_features==i)] <- T
  #       }
  #     }
  #
  #     if(any(remove)) ###should any of the pmp features be removed
  #     {
  #       remove <- selection[pmp_features[which(remove==T)]]
  #       remove_list[(count_remove+1):(count_remove+length(remove))] <- remove
  #       count_remove <- count_remove + length(remove)
  #     }
  #
  #     pbi <- pbi + 1
  #
  #     updatecounter <- updatecounter + 1
  #     if(updatecounter >= 10)
  #     {
  #       time_elapsed <- difftime(Sys.time(),start_time,units="secs")
  #       time_require <- (time_elapsed/(pbi/max))*(1-(pbi/max))
  #       td <- lubridate::seconds_to_period(time_require)
  #       time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))
  #
  #       updatecounter <- 0
  #       tcltk::setTkProgressBar(pb, pbi, label=base::paste( round(pbi/max*100, 0)," % done (",pbi,"/",max,", Time require: ",time_require,")",sep = ""))
  #     }
  #   }
  #   close(pb)
  #
  #   remove_list <- remove_list[1:count_remove]
  #
  #   print(base::paste(nrow(PMPs$features)-length(remove_list),"of PMPs show significant correlation with peptides of respective proteins."))
  #
  #   ###now remove PMPs which are not at all correlating with any peptide of the same protein
  #   if(length(remove_list)>0)
  #   {
  #     ###store excluded features and quant
  #     excluded_PMPs <- list(features=features[remove_list,],
  #                           features_intensity=features_intensity[remove_list,],
  #                           Ioncount_feature_sample_matrix=Ioncount_feature_sample_matrix[remove_list,],
  #                           feature_with_background_intensity=feature_with_background_intensity[remove_list,],
  #                           feature_with_background_intensity_imputed=feature_with_background_intensity_imputed[remove_list,],
  #                           Ioncount_feature_with_background_intensity=Ioncount_feature_with_background_intensity[remove_list,],
  #                           Alignment_scores=Alignment_scores[remove_list,],
  #                           Quant_pvals_signal_with_background=Quant_pvals_signal_with_background[remove_list,])
  #     save(excluded_PMPs,file = "excluded_PMPs.RData")
  #
  #     ###keep all the other features and continue
  #     features<-features[-remove_list,]
  #     features_intensity<-features_intensity[-remove_list,]
  #     Ioncount_feature_sample_matrix<-Ioncount_feature_sample_matrix[-remove_list,]
  #     feature_with_background_intensity<-feature_with_background_intensity[-remove_list,]
  #     feature_with_background_intensity_imputed<-feature_with_background_intensity_imputed[-remove_list,]
  #     Ioncount_feature_with_background_intensity<-Ioncount_feature_with_background_intensity[-remove_list,]
  #     Alignment_scores<-Alignment_scores[-remove_list,]
  #     Quant_pvals_signal_with_background<-Quant_pvals_signal_with_background[-remove_list,]
  #   }
  # }
  #
  # ###Filter for isotope features which at least in one sample show significant quantification
  # if(any(grepl("_i",features$Feature_name)))
  # {
  #   selection <- which(grepl("_i",features$Feature_name))
  #   quant_pvals_temp <- Quant_pvals_signal_with_background[selection,]
  #   remove_list <- which(matrixStats::rowMins(as.matrix(quant_pvals_temp),na.rm=T)>0.01)
  #   remove_list <- selection[remove_list]
  #
  #   if(length(remove_list)>0)
  #   {
  #     features<-features[-remove_list,]
  #     features_intensity<-features_intensity[-remove_list,]
  #     Ioncount_feature_sample_matrix<-Ioncount_feature_sample_matrix[-remove_list,]
  #     feature_with_background_intensity<-feature_with_background_intensity[-remove_list,]
  #     feature_with_background_intensity_imputed<-feature_with_background_intensity_imputed[-remove_list,]
  #     Ioncount_feature_with_background_intensity<-Ioncount_feature_with_background_intensity[-remove_list,]
  #     Alignment_scores<-Alignment_scores[-remove_list,]
  #     Quant_pvals_signal_with_background<-Quant_pvals_signal_with_background[-remove_list,]
  #   }
  # }

  ###store the updated results
  # setwd(base::paste(path_to_features,"\\Temporary_files",sep=""))
  # utils::write.table(features_intensity,base::paste("Features_quantification_signal_only",output_file_names_add,".txt",sep=""))
  # utils::write.table(Ioncount_feature_sample_matrix,base::paste("Features_quantification_ion_counts_signal_only",output_file_names_add,".txt",sep=""))
  # utils::write.table(feature_with_background_intensity,base::paste("Features_quantification_signal_background",output_file_names_add,".txt",sep=""))
  # utils::write.table(feature_with_background_intensity,base::paste("Features_quantification_signal_background_imputed",output_file_names_add,".txt",sep=""))
  # utils::write.table(Ioncount_feature_with_background_intensity,base::paste("Features_quantification_ion_count_signal_background",output_file_names_add,".txt",sep=""))
  # utils::write.table(Alignment_scores,base::paste("Alignment_scores",output_file_names_add,".txt",sep=""))
  # utils::write.table(Quant_pvals_signal_with_background,base::paste("Features_quantification_pvals_signal_background",output_file_names_add,".txt",sep=""))
  #

  #reorder samples if multiplicity > 1 as it should be ordered L, M, H like in MaxQ output
  if(multiplicity > 1)
  {
    #get raw sample names
    sample <- sort(unique(base::gsub("_Channel_light|_Channel_medium|_Channel_heavy","",colnames(features_intensity))))

    #define correct order
    if(multiplicity == 2)
    {
      labels <- c("_Channel_light","_Channel_heavy")
      sample <- base::paste(sort(rep(sample,2)),labels,sep="")
    }
    if(multiplicity == 3)
    {
      labels <- c("_Channel_light","_Channel_medium","_Channel_heavy")
      sample <- base::paste(sort(rep(sample,3)),labels,sep="")
    }

    #get numerical order for tables
    order <- match(sample,colnames(features_intensity))

    #change respective table column order
    features_intensity <- features_intensity[,order]
    Ioncount_feature_sample_matrix <- Ioncount_feature_sample_matrix[,order]
    feature_with_background_intensity <- feature_with_background_intensity[,order]
    feature_with_background_intensity_imputed <- feature_with_background_intensity_imputed[,order]
    Ioncount_feature_with_background_intensity <- Ioncount_feature_with_background_intensity[,order]
    peak_selected <- peak_selected[,order]
    S2B <- S2B[,order]
    pval_signal_with_background_quant <- pval_signal_with_background_quant[,order]
    alignment_scores_peaks_correct <- alignment_scores_peaks_correct[,order]
    alignment_scores_peaks_raw <- alignment_scores_peaks_raw[,order]
    mono_iso_alignment_summary <- mono_iso_alignment_summary[,order]
  }

  ###store results after filtering
  print("Save quantification results after alignment scoring and filter")

  temp_results <- list(features=features,
                       features_intensity=features_intensity,
                       Ioncount_feature_sample_matrix=Ioncount_feature_sample_matrix,
                       feature_with_background_intensity=feature_with_background_intensity,
                       feature_with_background_intensity_imputed=feature_with_background_intensity_imputed,
                       Ioncount_feature_with_background_intensity=Ioncount_feature_with_background_intensity,
                       peak_selected=peak_selected,
                       S2B=S2B,
                       pval_signal_with_background_quant=pval_signal_with_background_quant,
                       alignment_variability_score=alignment_variability_score,
                       alignment_scores_peaks_correct=alignment_scores_peaks_correct,
                       alignment_scores_peaks_raw=alignment_scores_peaks_raw,
                       mono_iso_alignment_summary=mono_iso_alignment_summary,
                       QC_data)
  setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
  save(temp_results,file = "Quantification_raw_results_with_scores_filtered.RData")

  # path_to_features <- "F:\\9_Spike_in_data_sets\\4_spike-in human Shen\\Requant\\19 - DDAicer reprocessed"
  # path_to_MaxQ_output <- "F:\\9_Spike_in_data_sets\\4_spike-in human Shen\\MaxQuant"
  # output_file_names_add <- "_DDAiceR_analysis"
  # n_cores <- 3
  # abundance_estimation_correction <- T
  # calc_protein_LFQ <- T
  # calc_peptide_LFQ <- F
  #
  # setwd(base::paste(path_to_features,"\\Temporary_files",sep=""))
  # load("Quantification_raw_results_with_scores_filtered.RData")
  # features <- temp_results$features
  # features_intensity <- temp_results$features_intensity
  # Ioncount_feature_sample_matrix <- temp_results$Ioncount_feature_sample_matrix
  # feature_with_background_intensity <- temp_results$feature_with_background_intensity
  # feature_with_background_intensity_imputed <- temp_results$feature_with_background_intensity_imputed
  # Ioncount_feature_with_background_intensity <- temp_results$Ioncount_feature_with_background_intensity
  # peak_selected <- temp_results$peak_selected
  # S2B <- temp_results$S2B
  # pval_signal_with_background_quant <- temp_results$pval_signal_with_background_quant
  # alignment_variability_score <- temp_results$alignment_variability_score
  # alignment_scores_peaks_correct <- temp_results$alignment_scores_peaks_correct
  # mono_iso_alignment_summary <- temp_results$mono_iso_alignment_summary
  # QC_data <- temp_results$QC_data
  # library(stringr)
  # library(doParallel)
  ###Perform protein level quantification

  correct_intensities <- function(features,feature_sample_matrix_requantified,pval_quant,MaxQ_peptides_quant,main="",corr_factor=NA)
  {
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    #suppressWarnings(suppressMessages(library(dplyr,quietly = T)))

    if(is.na(corr_factor)) ###no correction factor specified so correction will be determined and then applied
    {
      ###select monoisotopic features with significant quantification for determination of general trends in difference between MaxQ and Requant quantification
      select <- which(!grepl("_i|_d",features$Feature_name))
      feature_quant <- feature_sample_matrix_requantified[select,]
      feature_quant_pval <- pval_quant[select,]
      feature_quant[is.na(feature_quant_pval) | feature_quant_pval > 0.1] <- NA
      features_select <- features[select,]
      ###determine deviations between Requant and MaxQ in abundance estimations
      Requant_peptides_quant_seq <- stats::aggregate(10^feature_quant,by=list(Sequence=features_select$Sequence),FUN=sum)
      Requant_peptides_quant_seq[,-1] <- base::log2(Requant_peptides_quant_seq[,-1])
      ###bring MaxQ results table into same order as Requant output
      if(ncol(MaxQ_peptides_quant) != ncol(Requant_peptides_quant_seq))
      {
        colnames(MaxQ_peptides_quant) <- base::gsub("Intensity.","",colnames(MaxQ_peptides_quant))
        ordering <- vector("numeric",ncol(MaxQ_peptides_quant))
        for(i in 1:ncol(MaxQ_peptides_quant))
        {
          overlap <- grepl(base::paste(colnames(MaxQ_peptides_quant)[i],"$",sep=""),colnames(Requant_peptides_quant_seq))
          ordering[i] <- ifelse(any(overlap),which(overlap == T),NA)
        }
        MaxQ_peptides_quant <- MaxQ_peptides_quant[,which(!is.na(ordering))]
        ordering <- ordering[!is.na(ordering)]
        MaxQ_peptides_quant <- MaxQ_peptides_quant[,ordering]
      }
      comb <- dplyr::full_join(MaxQ_peptides_quant,Requant_peptides_quant_seq,by="Sequence")

      ncolumns <- ncol(feature_sample_matrix_requantified)

      dat <- base::as.data.frame(matrix(ncol=2,nrow=nrow(comb)*ncolumns))
      dat[,1] <- as.numeric(dat[,1])
      dat[,2] <- as.numeric(dat[,2])
      for(i in 1:ncolumns)
      {
        start_index <- (i-1)*nrow(comb)+1
        end_index <- i*nrow(comb)
        temp <- base::as.data.frame(cbind(comb[,i+1],comb[,i+1+ncolumns]))
        data.table::set(dat,as.integer(start_index:end_index),j=as.integer(1:2),value=temp)
      }
      colnames(dat) <- c("MaxQ","Requant")
      dat <- stats::na.omit(dat)
      if(main != "")main=base::paste(" ",main,sep="")
      grDevices::pdf(base::paste("Correct feature abundance estimations",main,".pdf",sep=""))
      graphics::smoothScatter(dat[,1],dat[,2],ylab="Requant, log2",xlab="MaxQ, log2",main="Correlation of MaxQ vs Requant on peptide level")
      graphics::abline(a=0,b=1)
      temp <- dat[c(order(dat[,2],decreasing = T)[1:(0.1*nrow(dat))],order(dat[,2],decreasing = F)[1:(0.1*nrow(dat))]),]
      y <- temp[,2]
      x <- temp[,1]
      fit <- stats::lm(y~x)
      graphics::abline(fit,lty=2,col="red")
      m <- summary(fit)$coefficients[2,1]
      intercept <- summary(fit)$coefficients[1,1]
      Rsq <- summary(fit)$r.squared
      posx <- graphics::par("usr")[1]+(graphics::par("usr")[2]-graphics::par("usr")[1])*0.15
      posy <- graphics::par("usr")[4]-(graphics::par("usr")[4]-graphics::par("usr")[3])*0.1
      graphics::text(posx,posy,labels = base::paste("R? =",round(Rsq,digits=2),"\nslope =",round(m,digits=2)))
      grDevices::dev.off()

      ###correction
      features_sample_matrix_corrected <- base::as.data.frame(base::log2(10^feature_sample_matrix_requantified))

      ###correct slope
      correction <- 1/m
      for(c in 1:ncol(features_sample_matrix_corrected))
      {
        features_sample_matrix_corrected[,c] <- (features_sample_matrix_corrected[,c]*correction)
      }
      return(list(features_sample_matrix_corrected,
                  correction_data=dat,
                  correction_fit=fit,
                  correction_factor=correction))
    }else ###correction factor is supplied so correct data accordingly
    {
      ###correction
      features_sample_matrix_corrected <- base::as.data.frame(base::log2(10^feature_sample_matrix_requantified))

      for(c in 1:ncol(features_sample_matrix_corrected))
      {
        features_sample_matrix_corrected[,c] <- features_sample_matrix_corrected[,c]*corr_factor
      }
      return(list(features_sample_matrix_corrected,
                  correction_factor=corr_factor))
    }

  }

  Top3_Protein_Quant <- function(features,feature_sample_matrix_requantified,Alignment_scores=NULL,Quant_pvals=NULL,S2B=NULL,use_overlapping=T,min_peps=2,quant_pvalue_cutoff=0.1,use_isotope_pmps=F)
  {
    pb <- tcltk::tkProgressBar(title = "Prepare Top3 quantification",label=base::paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    close(pb)
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    #suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
    #suppressWarnings(suppressMessages(library(matrixStats,quietly = T)))
    #suppressWarnings(suppressMessages(library(stringr,quietly = T)))

    feature_sample_matrix_requantified <- base::as.data.frame(feature_sample_matrix_requantified)
    if(!is.null(Quant_pvals))Quant_pvals <- base::as.data.frame(Quant_pvals)
    if(!is.null(Alignment_scores))Alignment_scores <- base::as.data.frame(Alignment_scores)
    if(!is.null(S2B))S2B <- base::as.data.frame(S2B)
    ##top3 method
    ###input matrix with samples in cols and rows correspond to unique peptides (log2 summed intensity over charge state and modification) of a respective protein
    Top3_quant <- function(pep_matrix,features_temp,Alignment_scores_temp=NULL,Quant_pvals_temp=NULL,S2B_temp=NULL)
    {
      sequence <- features_temp$Sequence

      if(length(sequence)>0)
      {
        pep_matrix <- stats::aggregate(2^pep_matrix, by=list(Sequence=sequence), FUN=sum,na.rm=T)
        pep_matrix <- pep_matrix[,-1]
        pep_matrix[pep_matrix==0] <- NA
        pep_matrix <- base::log2(pep_matrix)

        if(!is.null(Alignment_scores_temp))
        {
          score_matrix <- stats::aggregate(Alignment_scores_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          score_matrix <- score_matrix[,-1]
          score_matrix[score_matrix==0] <- NA
        }else
        {
          score_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(score_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(Quant_pvals_temp))
        {
          pval_matrix <- stats::aggregate(Quant_pvals_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          pval_matrix <- pval_matrix[,-1]
          pval_matrix[pval_matrix==0] <- NA
        }else
        {
          pval_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(pval_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(S2B_temp))
        {
          S2B_matrix <- stats::aggregate(S2B_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          S2B_matrix <- S2B_matrix[,-1]
          S2B_matrix[S2B_matrix==0] <- NA
        }else
        {
          S2B_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(S2B_matrix) <- colnames(pep_matrix)
        }

        top3_res <- NULL
        top3_score_res <- NULL
        top3_pval_res <- NULL
        top3_S2B_res <- NULL
        for(c in 1:ncol(pep_matrix))
        {
          top3 <- 2^pep_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
          sums <- sum(top3,na.rm=T)

          if(!is.null(Alignment_scores_temp))
          {
            top3_score <- score_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
            median_score <- stats::median(top3_score,na.rm=T)
          }else
          {
            median_score <- NA
          }

          if(!is.null(Quant_pvals_temp))
          {
            top3_pval <- pval_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
            median_pval <- stats::median(top3_pval,na.rm=T)
          }else
          {
            median_pval <- NA
          }

          if(!is.null(S2B_temp))
          {
            top3_S2B <- S2B_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
            median_S2B <- stats::median(top3_S2B,na.rm=T)
          }else
          {
            median_S2B <- NA
          }

          if(sums == 0)
          {
            sums <- NA
          }else
          {
            if(length(which(!is.na(top3)))>=min_peps) ### quantification only if at least 2 peptide quantifications are available
            {
              sums <- sums/length(which(!is.na(top3)))
            }else
            {
              sums <- NA
              median_score <- NA
              median_pval <- NA
              median_S2B <- NA
            }

          }
          top3_res <- append(top3_res,base::log2(sums))
          top3_score_res <- append(top3_score_res,median_score)
          top3_pval_res <- append(top3_pval_res,median_pval)
          top3_S2B_res <- append(top3_S2B_res,median_S2B)
        }
      }else
      {
        top3_res <- rep(NA,ncol(pep_matrix))
        top3_score_res <- rep(NA,ncol(score_matrix))
        top3_pval_res <- rep(NA,ncol(pval_matrix))
        top3_S2B_res <- rep(NA,ncol(pval_matrix))
      }

      names(top3_res) <- colnames(pep_matrix)
      names(top3_score_res) <- base::paste(colnames(score_matrix),"_median_score",sep="")
      names(top3_pval_res) <- base::paste(colnames(pval_matrix),"_median_pvals",sep="")
      names(top3_S2B_res) <- base::paste(colnames(pval_matrix),"_median_S2B",sep="")

      return(append(top3_res,append(top3_score_res,append(top3_pval_res,top3_S2B_res))))
    }

    if(use_isotope_pmps == F)
    {
      features_temp <- features[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      Quant_pvals_temp <- Quant_pvals[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      Alignment_scores_temp <- Alignment_scores[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      S2B_temp <- S2B[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
    }else
    {
      features_temp <- features[which(features$Protein != "" & !grepl(";",features$Sequence)),]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      Quant_pvals_temp <- Quant_pvals[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      Alignment_scores_temp <- Alignment_scores[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      S2B_temp <- S2B[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
    }

    if(!is.na(quant_pvalue_cutoff))
    {
      selection <- which(matrixStats::rowMins(as.matrix(Quant_pvals_temp),na.rm=T) < quant_pvalue_cutoff)

      features_temp <- features[selection,]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[selection,]
      Quant_pvals_temp <- Quant_pvals[selection,]
      Alignment_scores_temp <- Alignment_scores[selection,]
      S2B_temp <- S2B[selection,]
    }

    if(nrow(features_temp) > 0)
    {
      if(use_overlapping==T) ###use also features where a peptide is shared between 2 or more proteins
      {
        unique_proteins <- sort(unique(as.character(stringr::str_split(features_temp$Protein,"\\||;",simplify = T))))
      }else
      {
        unique_proteins <- sort(unique(features_temp$Protein[which(!grepl("\\||;",features_temp$Protein))]))
      }

      if(any(unique_proteins == ""))unique_proteins <- unique_proteins[-which(unique_proteins == "")]

      protein_TOP3 <- base::as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=length(unique_proteins),0))
      rownames(protein_TOP3) <- unique_proteins
      colnames(protein_TOP3) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),base::paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))

      max <- nrow(protein_TOP3)
      pb <- tcltk::tkProgressBar(title = "Perform Top3 quantification",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      for(i in 1:length(unique_proteins))
      {
        if(use_overlapping==T)
        {
          ind <- which(grepl(unique_proteins[i],features_temp$Protein))
        }else
        {
          ind <- which(features_temp$Protein == unique_proteins[i])
        }
        if(length(ind)>0)
        {
          data.table::set(protein_TOP3,as.integer(i),as.integer(1:ncol(protein_TOP3)),value=as.list(c(length(ind),as.numeric(Top3_quant(pep_matrix = feature_sample_matrix_requantified_temp[ind,],features_temp = features_temp[ind,],Alignment_scores_temp = Alignment_scores[ind,],Quant_pvals_temp=Quant_pvals[ind,],S2B_temp=S2B_temp[ind,])))))
        }else
        {
          data.table::set(protein_TOP3,as.integer(i),as.integer(1),value=0)
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)
    }else
    {
      protein_TOP3 <- base::as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=0,0))
      colnames(protein_TOP3) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),base::paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))
    }

    return(protein_TOP3)
  }

  Total_Protein_Quant <- function(features,feature_sample_matrix_requantified,Alignment_scores=NULL,Quant_pvals=NULL,S2B=NULL,use_overlapping=T,min_peps=2,quant_pvalue_cutoff=0.1,use_isotope_pmps=F)
  {
    pb <- tcltk::tkProgressBar(title = "Prepare Total quantification",label=base::paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    close(pb)
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    #suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
    #suppressWarnings(suppressMessages(library(matrixStats,quietly = T)))
    #suppressWarnings(suppressMessages(library(stringr,quietly = T)))

    feature_sample_matrix_requantified <- base::as.data.frame(feature_sample_matrix_requantified)
    if(!is.null(Quant_pvals))Quant_pvals <- base::as.data.frame(Quant_pvals)
    if(!is.null(Alignment_scores))Alignment_scores <- base::as.data.frame(Alignment_scores)
    if(!is.null(S2B))S2B <- base::as.data.frame(S2B)
    ##total method
    ###input matrix with samples in cols and rows correspond to unique peptides (log2 summed intensity over charge state and modification) of a respective protein
    Total_quant <- function(pep_matrix,features_temp,Alignment_scores_temp=NULL,Quant_pvals_temp=NULL,S2B_temp=NULL,Quant_cutoff=4)
    {
      sequence <- features_temp$Sequence

      if(length(sequence)>0)
      {
        pep_matrix <- stats::aggregate(2^pep_matrix, by=list(Sequence=sequence), FUN=sum,na.rm=T)
        pep_matrix <- pep_matrix[,-1]
        pep_matrix[pep_matrix==0] <- NA
        pep_matrix <- base::log2(pep_matrix)

        if(!is.null(Alignment_scores_temp))
        {
          score_matrix <- stats::aggregate(Alignment_scores_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          score_matrix <- score_matrix[,-1]
          score_matrix[score_matrix==0] <- NA
        }else
        {
          score_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(score_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(Quant_pvals_temp))
        {
          pval_matrix <- stats::aggregate(Quant_pvals_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          pval_matrix <- pval_matrix[,-1]
          pval_matrix[pval_matrix==0] <- NA
        }else
        {
          pval_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(pval_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(S2B_temp))
        {
          S2B_matrix <- stats::aggregate(S2B_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          S2B_matrix <- S2B_matrix[,-1]
          S2B_matrix[S2B_matrix==0] <- NA
        }else
        {
          S2B_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(S2B_matrix) <- colnames(pep_matrix)
        }

        total_res <- NULL
        total_score_res <- NULL
        total_pval_res <- NULL
        total_S2B_res <- NULL
        for(c in 1:ncol(pep_matrix))
        {
          total <- 2^pep_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
          sums <- sum(total,na.rm=T)

          if(!is.null(Alignment_scores_temp))
          {
            total_score <- score_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
            median_score <- stats::weighted.mean(total_score,total,na.rm=T)#stats::median(total_score,na.rm=T)
          }else
          {
            median_score <- NA
          }

          if(!is.null(Quant_pvals_temp))
          {
            total_pval <- pval_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
            median_pval <- stats::weighted.mean(total_pval,total,na.rm=T)#stats::median(total_pval,na.rm=T)
          }else
          {
            median_pval <- NA
          }

          if(!is.null(S2B_temp))
          {
            total_S2B <- S2B_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
            median_S2B <- stats::weighted.mean(total_S2B,total,na.rm=T)#stats::median(total_S2B,na.rm=T)
          }else
          {
            median_S2B <- NA
          }

          if(sums == 0)
          {
            sums <- NA
          }else
          {
            if(length(which(!is.na(total)))>=min_peps) ### quantification only if at least n peptide quantifications are available
            {
              sums <- sums/length(which(!is.na(total)))
            }else
            {
              sums <- NA
              median_score <- NA
              median_pval <- NA
              median_S2B <- NA
            }

          }

          total_res <- append(total_res,base::log2(sums))
          total_score_res <- append(total_score_res,median_score)
          total_pval_res <- append(total_pval_res,median_pval)
          total_S2B_res <- append(total_S2B_res,median_S2B)
        }
      }else
      {
        total_res <- rep(NA,ncol(pep_matrix))
        total_score_res <- rep(NA,ncol(score_matrix))
        total_pval_res <- rep(NA,ncol(pval_matrix))
        total_S2B_res <- rep(NA,ncol(pval_matrix))
      }

      names(total_res) <- colnames(pep_matrix)
      names(total_score_res) <- base::paste(colnames(score_matrix),"_median_score",sep="")
      names(total_pval_res) <- base::paste(colnames(pval_matrix),"_median_pvals",sep="")
      names(total_S2B_res) <- base::paste(colnames(pval_matrix),"_median_S2B",sep="")

      return(append(total_res,append(total_score_res,append(total_pval_res,total_S2B_res))))
    }

    if(use_isotope_pmps == F)
    {
      features_temp <- features[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      Quant_pvals_temp <- Quant_pvals[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      Alignment_scores_temp <- Alignment_scores[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
      S2B_temp <- S2B[which(features$Protein != "" & !grepl(";",features$Sequence) & !grepl("_i|_pmp",features$Feature_name)),]
    }else
    {
      features_temp <- features[which(features$Protein != "" & !grepl(";",features$Sequence)),]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      Quant_pvals_temp <- Quant_pvals[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      Alignment_scores_temp <- Alignment_scores[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
      S2B_temp <- S2B[which(features$Protein != "" & !grepl(";",features$Sequence) ),]
    }

    if(!is.na(quant_pvalue_cutoff))
    {
      selection <- which(matrixStats::rowMins(as.matrix(Quant_pvals_temp),na.rm=T) < quant_pvalue_cutoff)

      features_temp <- features[selection,]
      feature_sample_matrix_requantified_temp <- feature_sample_matrix_requantified[selection,]
      Quant_pvals_temp <- Quant_pvals[selection,]
      Alignment_scores_temp <- Alignment_scores[selection,]
      S2B_temp <- S2B[selection,]
    }

    if(nrow(features_temp) > 0)
    {
      if(use_overlapping==T) ###use also features where a peptide is shared between 2 or more proteins
      {
        unique_proteins <- sort(unique(as.character(stringr::str_split(features_temp$Protein,"\\||;",simplify = T))))
      }else
      {
        unique_proteins <- sort(unique(features_temp$Protein[which(!grepl("\\||;",features_temp$Protein))]))
      }

      if(any(unique_proteins == ""))unique_proteins <- unique_proteins[-which(unique_proteins == "")]

      protein_total <- base::as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=length(unique_proteins),0))
      rownames(protein_total) <- unique_proteins
      colnames(protein_total) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),base::paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))

      max <- nrow(protein_total)
      pb <- tcltk::tkProgressBar(title = "Perform total quantification",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      for(i in 1:length(unique_proteins))
      {
        if(use_overlapping==T)
        {
          ind <- which(grepl(unique_proteins[i],features_temp$Protein))
        }else
        {
          ind <- which(features_temp$Protein == unique_proteins[i])
        }
        if(length(ind)>0)
        {
          data.table::set(protein_total,as.integer(i),as.integer(1:ncol(protein_total)),value=as.list(c(length(ind),as.numeric(Total_quant(pep_matrix = feature_sample_matrix_requantified_temp[ind,],features_temp = features_temp[ind,],Alignment_scores_temp = Alignment_scores[ind,],Quant_pvals_temp=Quant_pvals[ind,],S2B_temp=S2B_temp[ind,])))))
        }else
        {
          data.table::set(protein_total,as.integer(i),as.integer(1),value=0)
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- lubridate::seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)
    }else
    {
      protein_total <- base::as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=0,0))
      colnames(protein_total) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),base::paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),base::paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))
    }

    return(protein_total)
  }

  crap <- gc(F)
  ###perform protein level aggregation
  if(!is.null(path_to_MaxQ_output))
  {
    ###load MaxQ peptide results
    MaxQ_peptides <- utils::read.table(base::paste(path_to_MaxQ_output,"/peptides.txt",sep=""),sep="\t",header=T)
    MaxQ_peptides <- MaxQ_peptides[-which(MaxQ_peptides$Potential.contaminant == "+" | MaxQ_peptides$Reverse == "+"),]
    if(multiplicity == 1)MaxQ_peptides_quant <- MaxQ_peptides[,which(grepl("Intensity\\.",colnames(MaxQ_peptides)))]
    if(multiplicity > 1)MaxQ_peptides_quant <- MaxQ_peptides[,which(grepl("Intensity\\.[LMH]\\.",colnames(MaxQ_peptides)))]
    MaxQ_peptides_quant[MaxQ_peptides_quant==0] <- NA
    MaxQ_peptides_quant <- base::log2(MaxQ_peptides_quant)
    MaxQ_peptides_leading_razor <- base::data.frame(Sequence=MaxQ_peptides$Sequence,Leading_razor=MaxQ_peptides$Leading.razor.protein)
    MaxQ_peptides_leading_razor$Leading_razor <- as.character(MaxQ_peptides_leading_razor$Leading_razor)
    MaxQ_peptides <- base::data.frame(Sequence=MaxQ_peptides$Sequence,MaxQ_peptides_quant)

    MaxQ_protein_groups <- utils::read.table(base::paste(path_to_MaxQ_output,"/proteinGroups.txt",sep=""),sep="\t",header=T)

    ###check if IDs were correctly parsed, if not, try to parse with SwissProt or Trembl
    if(!any(colnames(MaxQ_protein_groups) == "Gene.names") & any(colnames(MaxQ_protein_groups) == "Protein.IDs")) ##not correctly parsed but contains trembl or swissprot fasta headers
    {
      print("Fasta file was not correctly parsed during search. Try to base::paste fasta headers ...")
      if(any(grepl(">sp|>tr",MaxQ_protein_groups$Fasta.headers)))
      {
        print("Detected Swiss-Prot and/or TrEMBL fasta headers.")
        ##protein level
        ###parsing was not performed correctly so we have to try to do this here expecting swissprot or trembl fasta headers
        MaxQ_protein_groups$Gene.names <- ""
        MaxQ_protein_groups$Organism <- ""
        MaxQ_protein_groups$Protein.IDs <- as.character(MaxQ_protein_groups$Protein.IDs)

        fasta_headers <- stringr::str_split(base::gsub(">","",MaxQ_protein_groups$Fasta.headers),"\\|",simplify = T)

        ##extract gene name and species information
        GN <- vector("character",nrow(MaxQ_protein_groups))
        #OS <- vector("character",nrow(MaxQ_protein_groups))
        ID <- vector("character",nrow(MaxQ_protein_groups))
        GN_start <- unlist(lapply(base::paste(gregexpr("GN=",MaxQ_protein_groups$Fasta.headers),sep=","), `[[`, 1))
        GN_start <- base::gsub("c\\(|\\)","",GN_start)
        #OS_start <- unlist(lapply(base::paste(gregexpr("OS=",MaxQ_protein_groups$Fasta.headers),sep=","), `[[`, 1))
        #OS_start <- base::gsub("c\\(|\\)","",OS_start)
        for(i in 1:nrow(MaxQ_protein_groups))
        {
          if(GN_start[i] != "-1")
          {
            indices <- as.numeric(unlist(stringr::str_split(GN_start[i],","))) + 3
            stop <- regexpr(" ",substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+10))
            GN_temp <- substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+stop-2)
            GN[i] <- base::paste(GN_temp,collapse = ";")
          }
          # if(OS_start[i] != "-1")
          # {
          #   indices <- as.numeric(unlist(stringr::str_split(OS_start[i],","))) + 3
          #   stop <- regexpr("GN=",substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+100))
          #   OS_temp <- substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+stop-3)
          #   OS[i] <- base::paste(OS_temp,collapse = ";")
          # }
          header <- fasta_headers[i,c(2,4,6)][which(fasta_headers[i,c(2,4,6)] != "")]
          if(length(header)>0)ID[i] <- base::paste(header,collapse=";")
        }
        CON_REV <- !grepl("^CON_|^REV_",MaxQ_protein_groups$Protein.IDs)
        MaxQ_protein_groups$Protein.IDs <- ifelse(CON_REV,ID, MaxQ_protein_groups$Protein.IDs)
        MaxQ_protein_groups$Gene.names <- ifelse(CON_REV,GN,"")
        #MaxQ_protein_groups$Organism <- ifelse(CON_REV,OS,"")
        MaxQ_protein_groups$Majority.protein.IDs <- MaxQ_protein_groups$Protein.IDs

        ##peptide level
        if(any(grepl("\\|",MaxQ_peptides_leading_razor$Leading_razor)))
        {
          temp <- stringr::str_split(MaxQ_peptides_leading_razor$Leading_razor,"\\|",simplify = T)
          MaxQ_peptides_leading_razor$Leading_razor <- ifelse(temp[,2] != "",temp[,2],MaxQ_peptides_leading_razor$Leading_razor)
        }

        ##requantification features
        temp <- stringr::str_split(features$Protein,"\\||;",simplify=T)
        ID <- vector("character",nrow(temp))
        for(i in 1:nrow(temp))
        {
          sel <- which(temp[i,] == "sp")+1
          if(length(sel)>0)
          {
            ID[i] <- base::paste(temp[i,sel],collapse=";")
          }else
          {
            sel <- which(temp[i,] != "")
            ID[i] <- base::paste(temp[i,sel],collapse=";")
          }
        }
        features$Protein <- ID
      }else
      {
        print("No supported fasta headers were detected. Next steps might be not fully working.")
      }
    }

    #####Correct for overestimation of high intense features
    if(abundance_estimation_correction == F) ###no correction
    {
      feature_with_background_intensity <- base::log2(10^feature_with_background_intensity)
      feature_with_background_intensity_imputed <- base::log2(10^feature_with_background_intensity_imputed)

    }else ##correct abundance estimations based on MaxQ peptide abundance estimations
    {
      ###determine abundance correction factors based on MaxQ peptide intensities and data without imputation
      cor_res <- correct_intensities(features,feature_sample_matrix_requantified = feature_with_background_intensity,pval_quant = pval_signal_with_background_quant,MaxQ_peptides_quant = MaxQ_peptides,main="Signal_Background_intensity")

      feature_with_background_intensity <- cor_res[[1]]
      QC_data[["Abundance_correction"]] <- list(correction_data=cor_res[[2]],
                                                correction_fit=cor_res[[3]],
                                                correction_factor=cor_res[[4]])

      cor_res <- correct_intensities(features,feature_sample_matrix_requantified = feature_with_background_intensity_imputed,pval_quant = pval_signal_with_background_quant,MaxQ_peptides_quant = MaxQ_peptides,main="Signal_Background_intensity_imputed",corr_factor = cor_res[[4]])

      feature_with_background_intensity_imputed <- cor_res[[1]]
    }

    ###get leading razor ID per peptide sequence
    features$Protein <- as.character(features$Protein)
    features$all_matching_Proteins <- features$Protein
    features$Protein <- as.character(MaxQ_peptides_leading_razor$Leading_razor[match(features$Sequence,MaxQ_peptides_leading_razor$Sequence)])
    features$Protein[is.na(features$Protein)] <- features$all_matching_Proteins[is.na(features$Protein)]

    ###perform peptide level LFQ
    if(calc_peptide_LFQ == T)
    {
      LFQ_peptide_quant_process <- function(features,features_quant,n_cores,label="Perform peptide-LFQ-quantification",num_ratio_samples=NA,TopN=5,seed=1)
      {
        set.seed(seed)
        ##prepare for MaxLFQ algorithm
        #suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))
        #suppressWarnings(suppressMessages(library(data.table,quietly = T)))

        calculate_LFQ <- function(peptide_quant_data,min_num_ratios=2,num_ratio_samples=NA)
        {
          #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
          #suppressWarnings(suppressMessages(library(robustbase,quietly = T)))
          #error function which is used to optimize ratios
          least_square_error <- function(par,ratio_mat)
          {
            sum <- 0
            vals <- NULL
            if(ncol(ratio_mat)>1)
            {
              for(c in 1:(ncol(ratio_mat)-1))
              {
                for(r in (c+1):nrow(ratio_mat))
                {
                  val <- (base::log2(ratio_mat[r,c])-base::log2(par[r])+base::log2(par[c]))^2
                  vals <- append(vals,val)
                  if(!is.na(val)) sum <- sum + val
                }
              }
            }

            if(sum == 0){sum = NA}
            return(sum)
          }
          #unlog intensities
          peptide_quant_data <- 2^peptide_quant_data
          #calculate summed intensities per sample
          totalsum_per_sample <- colSums(peptide_quant_data,na.rm=T)
          #determine median ratio matrix between all samples
          ratio_mat <- base::as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data)))
          if(ncol(ratio_mat)>1)
          {
            for(c in 1:(ncol(ratio_mat)-1))
            {
              for(r in (c+1):nrow(ratio_mat))
              {
                ratios <- peptide_quant_data[,r]/peptide_quant_data[,c]
                if(length(which(!is.na(ratios))) >= min_num_ratios){ratio_mat[r,c] <- stats::median(ratios,na.rm=T)}
              }
            }
          }
          if(is.na(num_ratio_samples)) ###calculate ratios over all samples
          {
            #define start parameter
            start_par <- c(rep(1,ncol(ratio_mat)))
            #now find optimum in ratios to best recover true observed ratios between samples
            res_ratio <- NULL
            try(res_ratio <- stats::optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)
            if(!is.null(res_ratio))
            {
              #normalize ratios to sample with highest intensity
              ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample == max(totalsum_per_sample))]
              #finally calculate log2 lfq protein intensities per sample
              lfq <- base::log2(ratio_norm*totalsum_per_sample[which(totalsum_per_sample == max(totalsum_per_sample))])
              #remove quant values for samples were no ratios were available
              if(any(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
              {
                sel <- as.numeric(which(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
                lfq[sel] <- NA
              }
            }else
            {
              lfq <- rep(NA,ncol(peptide_quant_data))
            }

          }else ###calculate ratios only for a subset of samples
          {
            lfq <- base::as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data),0))
            for(s in 1:ncol(peptide_quant_data))
            {
              ###randomly select up to 6 other samples from list of samples which also contain quantifications
              samples_with_quant <- as.numeric(which(colSums(peptide_quant_data,na.rm=T) > 0))
              samples_with_quant <- samples_with_quant[which(samples_with_quant != s)]
              if(length(samples_with_quant)>0)
              {
                samples_for_comparison_per_sample <- sample(samples_with_quant,ifelse(length(samples_with_quant)>num_ratio_samples,num_ratio_samples,length(samples_with_quant)))
              }else
              {
                samples_for_comparison_per_sample <- sample(c(1:ncol(peptide_quant_data))[-i],size = num_ratio_samples,replace = F)
              }

              ratio_mat_temp <- ratio_mat[c(s,sort(samples_for_comparison_per_sample)),c(s,sort(samples_for_comparison_per_sample))]
              if(length(which(!is.na(ratio_mat_temp))) < 3 & any(!is.na(peptide_quant_data[,s])))##to few of selected random samples show an observed intensity ratio but protein is quantified in current sample
              {
                ###randomly select other samples

              }

              totalsum_per_sample_temp <- totalsum_per_sample[c(s,sort(samples_for_comparison_per_sample))]
              #define start parameter
              start_par <- c(rep(1,ncol(ratio_mat_temp)))
              #now find optimum in ratios to best recover true observed ratios between samples
              res_ratio <- NULL
              try(res_ratio <- stats::optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat_temp,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)

              if(!is.null(res_ratio))
              {
                #normalize ratios to sample with highest intensity
                ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp,na.rm=T))]
                #finally calculate log2 lfq protein intensities per sample
                temp_lfq <- base::log2(ratio_norm*totalsum_per_sample_temp[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp))])

                data.table::set(lfq,as.integer(s),as.integer(c(s,sort(samples_for_comparison_per_sample))),as.list(temp_lfq))
              }

            }
            lfq[lfq==0] <- NA
            lfq <- matrixStats::colMedians(as.matrix(lfq),na.rm=T)

            #remove quant values for samples were no ratios were available
            if(any(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
            {
              sel <- as.numeric(which(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
              lfq[sel] <- NA
            }


          }
          return(lfq)
        }

        ###txtProgressBar from package pbarETA (Francesco Napolitano) - License: LGPL-3
        txtProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA,
                                    title, label, style = 3, file = "")
        {
          formatTime <- function(seconds) {
            if (seconds == Inf || is.nan(seconds) || is.na(seconds))
              return("NA")
            seconds <- round(seconds)
            sXmin <- 60
            sXhr <- sXmin * 60
            sXday <- sXhr * 24
            sXweek <- sXday * 7
            sXmonth <- sXweek * 4.22
            sXyear <- sXmonth * 12
            years <- floor(seconds/sXyear)
            seconds <- seconds - years * sXyear
            months <- floor(seconds/sXmonth)
            seconds <- seconds - months * sXmonth
            weeks <- floor(seconds/sXweek)
            seconds <- seconds - weeks * sXweek
            days <- floor(seconds/sXday)
            seconds <- seconds - days * sXday
            hours <- floor(seconds/sXhr)
            seconds <- seconds - hours * sXhr
            minutes <- floor(seconds/sXmin)
            seconds <- seconds - minutes * sXmin
            ETA <- c(years, months, days, hours, minutes, seconds)
            startst <- which(ETA > 0)[1]
            if (is.na(startst))
              startst <- 6
            starts <- min(startst, 4)
            fmtstr <- rep("%02d", length(ETA))[startst:length(ETA)]
            fmtstr <- base::paste(fmtstr, collapse = ":")
            return(do.call(sprintf, as.list(c(as.list(fmtstr), ETA[startst:length(ETA)]))))
          }
          if (!identical(file, "") && !(inherits(file, "connection") &&
                                        isOpen(file)))
            stop("'file' must be \"\" or an open connection object")
          if (!style %in% 1L:3L)
            style <- 1
          .val <- initial
          .killed <- FALSE
          .nb <- 0L
          .pc <- -1L
          .time0 <- NA
          .timenow <- NA
          .firstUpdate <- T
          nw <- nchar(char, "w")
          if (is.na(width)) {
            width <- getOption("width")
            if (style == 3L)
              width <- width - 10L
            width <- trunc(width/nw)
          }
          if (max <= min)
            stop("must have 'max' > 'min'")
          up1 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb < nb) {
              cat(base::paste(rep.int(char, nb - .nb), collapse = ""),
                  file = file)
              utils::flush.console()
            }
            else if (.nb > nb) {
              cat("\r", base::paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            .nb <<- nb
          }
          up2 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb <= nb) {
              cat("\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            else {
              cat("\r", base::paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            .nb <<- nb
          }
          up3 <- function(value, calledOnCreation = F) {
            timenow <- proc.time()[["elapsed"]]
            if (!calledOnCreation && .firstUpdate) {
              .time0 <<- timenow
              .timenow <<- timenow
              .firstUpdate <<- F
            }
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            pc <- round(100 * (value - min)/(max - min))
            if (nb == .nb && pc == .pc && timenow - .timenow < 1)
              return()
            .timenow <<- timenow
            span <- timenow - .time0
            timeXiter <- span/(.val - min)
            ETA <- (max - .val) * timeXiter
            ETAstr <- formatTime(ETA)
            cat(base::paste(c("\r  |", rep.int(" ", nw * width + 6)),
                      collapse = ""), file = file)
            cat(base::paste(c("\r  |", rep.int(char, nb), rep.int(" ",
                                                            nw * (width - nb)), sprintf("| %3d%%", pc), ", ETA ",
                        ETAstr), collapse = ""), file = file)
            utils::flush.console()
            .nb <<- nb
            .pc <<- pc
          }
          getVal <- function() .val
          kill <- function() if (!.killed) {
            cat("\n", file = file)
            utils::flush.console()
            .killed <<- TRUE
          }
          up <- switch(style, up1, up2, up3)
          up(initial, T)
          structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
        }

        sel <- which(!grepl(";|,",features$Sequence))
        unique_peptides <- sort(unique(base::paste(features$Sequence[sel],features$Modifications[sel],sep="_")))
        unique_seq <- base::substr(unique_peptides,1,regexpr("_",unique_peptides)-1)
        unique_mod <- base::substr(unique_peptides,regexpr("_",unique_peptides)+1,nchar(unique_peptides))

        LFQ_peptide_quant <- base::as.data.frame(matrix(ncol=ncol(features_quant),nrow=length(unique_peptides),0))
        colnames(LFQ_peptide_quant) <- colnames(features_quant)

        ###Perform LFQ quantification
        cl <- snow::makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        doSNOW::registerDoSNOW(cl)
        iterations <- nrow(LFQ_peptide_quant)
        progress <- function(n) utils::setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        start <- Sys.time()
        print(base::paste(label," (",Sys.time(),")",sep=""))
        pb <- txtProgressBar(max = iterations, style = 3)
        res_LFQ <- foreach::foreach(i=1:nrow(LFQ_peptide_quant),.options.snow = opts) %dopar%
          {
            sel <- which(features$Sequence == unique_seq[i] & features$Modifications == unique_mod[i])
            testdat = features_quant[sel,]
            if(nrow(testdat) > TopN) ###if more than 5 peptides are available select TopN peptides according to intensity over all samples (select highest abundant peptides as quantifications are more accurate)
            {
              total_abundance <- rowSums(testdat,na.rm=T)
              testdat <- testdat[order(total_abundance,decreasing = T)[1:TopN],]
            }
            if(nrow(testdat)>1)
            {
              res <- calculate_LFQ(peptide_quant_data = testdat,min_num_ratios = 1,num_ratio_samples=num_ratio_samples)
            }else
            {
              res <- as.numeric(testdat)
              names(res) <- base::paste("V",1:length(res),sep="")
            }
            return(res)
          }
        snow::stopCluster(cl)
        close(pb)

        #Combine results
        for(i in 1:length(res_LFQ))
        {
          if(length(res_LFQ[[i]])>0)
          {
            data.table::set(LFQ_peptide_quant,as.integer(i),as.integer(1:ncol(LFQ_peptide_quant)),as.list(as.numeric(res_LFQ[[i]])))
          }
        }

        #add information about number of quant features per protein
        count_quant_features <- plyr::count(sort(base::paste(features$Sequence[sel],features$Modifications[sel],sep="_")))

        LFQ_peptide_quant <- base::data.frame(Sequence=unique_seq,
                                        Modifications=unique_mod,
                                        Protein=features$Protein[match(unique_seq,features$Sequence)],
                                        num_quant_features=count_quant_features$freq[match(unique_peptides,count_quant_features$x)],
                                        LFQ_peptide_quant)

        colnames(LFQ_peptide_quant) <- base::gsub("^X","",colnames(LFQ_peptide_quant))

        end <- Sys.time()
        print(base::paste("Finished peptide LFQ-quantification (",Sys.time(),")",sep=""))
        print(end - start)
        #replace 0 by NA
        LFQ_peptide_quant[LFQ_peptide_quant==0] <- NA

        return(LFQ_peptide_quant)
      }

      if(ncol(feature_with_background_intensity) <= 10)num_ratio_samples <- NA
      if(ncol(feature_with_background_intensity) > 10)num_ratio_samples <- 6

      LFQ_peptide_quant_with_background <- LFQ_peptide_quant_process(features,feature_with_background_intensity,n_cores,label="Perform peptide LFQ-quantification for non-imputed data",num_ratio_samples = num_ratio_samples)
      LFQ_peptide_quant_with_background_imputed <- LFQ_peptide_quant_process(features,feature_with_background_intensity_imputed,n_cores,label="Perform peptide LFQ-quantification for imputed data",num_ratio_samples = num_ratio_samples)
      save(LFQ_peptide_quant_with_background,LFQ_peptide_quant_with_background_imputed,file = "Peptide_LFQ_temp.RData")
    }

    ###perform protein level quantification
    if(calc_protein_LFQ == T)
    {
      ##implementation of the MaxLFQ algorithm. num_ratio_samples indicates between how many samples the ratio matrices should be determined.
      #if num_ratio_samples is set to NA it will perform least-square analysis between all samples
      #if num_ratio_samples is set to a number < number of samples least square analysis is performed for randomly picked n (=num_ratio_samples) samples to reduced computation time
      LFQ_protein_quant_process <- function(features,features_quant,n_cores,label="Perform LFQ-quantification",num_ratio_samples=NA,TopN=5,seed=1)
      {
        set.seed(seed)
        ##prepare for MaxLFQ algorithm
        #suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))
        #suppressWarnings(suppressMessages(library(data.table,quietly = T)))

        calculate_LFQ <- function(peptide_quant_data,min_num_ratios=2,num_ratio_samples=NA)
        {
          #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
          #suppressWarnings(suppressMessages(library(robustbase,quietly = T)))
          #error function which is used to optimize ratios
          least_square_error <- function(par,ratio_mat)
          {
            sum <- 0
            vals <- NULL
            if(ncol(ratio_mat)>1)
            {
              for(c in 1:(ncol(ratio_mat)-1))
              {
                for(r in (c+1):nrow(ratio_mat))
                {
                  val <- (base::log2(ratio_mat[r,c])-base::log2(par[r])+base::log2(par[c]))^2
                  vals <- append(vals,val)
                  if(!is.na(val)) sum <- sum + val
                }
              }
            }

            if(sum == 0){sum = NA}
            return(sum)
          }
          #unlog intensities
          peptide_quant_data <- 2^peptide_quant_data
          #calculate summed intensities per sample
          totalsum_per_sample <- colSums(peptide_quant_data,na.rm=T)
          #determine median ratio matrix between all samples
          ratio_mat <- base::as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data)))
          if(ncol(ratio_mat)>1)
          {
            for(c in 1:(ncol(ratio_mat)-1))
            {
              for(r in (c+1):nrow(ratio_mat))
              {
                ratios <- peptide_quant_data[,r]/peptide_quant_data[,c]
                if(length(which(!is.na(ratios))) >= min_num_ratios){ratio_mat[r,c] <- stats::median(ratios,na.rm=T)}
              }
            }
          }
          if(is.na(num_ratio_samples)) ###calculate ratios over all samples
          {
            #define start parameter
            start_par <- c(rep(1,ncol(ratio_mat)))
            #now find optimum in ratios to best recover true observed ratios between samples
            res_ratio <- NULL
            try(res_ratio <- stats::optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)
            if(!is.null(res_ratio))
            {
              #normalize ratios to sample with highest intensity
              ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample == max(totalsum_per_sample))]
              #finally calculate log2 lfq protein intensities per sample
              lfq <- base::log2(ratio_norm*totalsum_per_sample[which(totalsum_per_sample == max(totalsum_per_sample))])
              #remove quant values for samples were no ratios were available
              if(any(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
              {
                sel <- as.numeric(which(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
                lfq[sel] <- NA
              }
            }else
            {
              lfq <- rep(NA,ncol(peptide_quant_data))
            }

          }else ###calculate ratios only for a subset of samples
          {
            lfq <- base::as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data),0))
            for(s in 1:ncol(peptide_quant_data))
            {
              ###randomly select up to 6 other samples from list of samples which also contain quantifications
              samples_with_quant <- as.numeric(which(colSums(peptide_quant_data,na.rm=T) > 0))
              samples_with_quant <- samples_with_quant[which(samples_with_quant != s)]
              if(length(samples_with_quant)>0)
              {
                samples_for_comparison_per_sample <- sample(samples_with_quant,ifelse(length(samples_with_quant)>num_ratio_samples,num_ratio_samples,length(samples_with_quant)))
              }else
              {
                samples_for_comparison_per_sample <- sample(c(1:ncol(peptide_quant_data))[-i],size = num_ratio_samples,replace = F)
              }

              ratio_mat_temp <- ratio_mat[c(s,sort(samples_for_comparison_per_sample)),c(s,sort(samples_for_comparison_per_sample))]
              if(length(which(!is.na(ratio_mat_temp))) < 3 & any(!is.na(peptide_quant_data[,s])))##to few of selected random samples show an observed intensity ratio but protein is quantified in current sample
              {
                ###randomly select other samples

              }

              totalsum_per_sample_temp <- totalsum_per_sample[c(s,sort(samples_for_comparison_per_sample))]
              #define start parameter
              start_par <- c(rep(1,ncol(ratio_mat_temp)))
              #now find optimum in ratios to best recover true observed ratios between samples
              res_ratio <- NULL
              try(res_ratio <- stats::optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat_temp,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)

              if(!is.null(res_ratio))
              {
                #normalize ratios to sample with highest intensity
                ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp,na.rm=T))]
                #finally calculate log2 lfq protein intensities per sample
                temp_lfq <- base::log2(ratio_norm*totalsum_per_sample_temp[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp))])

                data.table::set(lfq,as.integer(s),as.integer(c(s,sort(samples_for_comparison_per_sample))),as.list(temp_lfq))
              }

            }
            lfq[lfq==0] <- NA
            lfq <- matrixStats::colMedians(as.matrix(lfq),na.rm=T)

            #remove quant values for samples were no ratios were available
            if(any(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
            {
              sel <- as.numeric(which(colSums(!is.na(ratio_mat)) == 0 & rowSums(!is.na(ratio_mat)) == 0))
              lfq[sel] <- NA
            }


          }
          return(lfq)
        }

        ###txtProgressBar from package pbarETA (Francesco Napolitano) - License: LGPL-3
        txtProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA,
                                    title, label, style = 3, file = "")
        {
          formatTime <- function(seconds) {
            if (seconds == Inf || is.nan(seconds) || is.na(seconds))
              return("NA")
            seconds <- round(seconds)
            sXmin <- 60
            sXhr <- sXmin * 60
            sXday <- sXhr * 24
            sXweek <- sXday * 7
            sXmonth <- sXweek * 4.22
            sXyear <- sXmonth * 12
            years <- floor(seconds/sXyear)
            seconds <- seconds - years * sXyear
            months <- floor(seconds/sXmonth)
            seconds <- seconds - months * sXmonth
            weeks <- floor(seconds/sXweek)
            seconds <- seconds - weeks * sXweek
            days <- floor(seconds/sXday)
            seconds <- seconds - days * sXday
            hours <- floor(seconds/sXhr)
            seconds <- seconds - hours * sXhr
            minutes <- floor(seconds/sXmin)
            seconds <- seconds - minutes * sXmin
            ETA <- c(years, months, days, hours, minutes, seconds)
            startst <- which(ETA > 0)[1]
            if (is.na(startst))
              startst <- 6
            starts <- min(startst, 4)
            fmtstr <- rep("%02d", length(ETA))[startst:length(ETA)]
            fmtstr <- base::paste(fmtstr, collapse = ":")
            return(do.call(sprintf, as.list(c(as.list(fmtstr), ETA[startst:length(ETA)]))))
          }
          if (!identical(file, "") && !(inherits(file, "connection") &&
                                        isOpen(file)))
            stop("'file' must be \"\" or an open connection object")
          if (!style %in% 1L:3L)
            style <- 1
          .val <- initial
          .killed <- FALSE
          .nb <- 0L
          .pc <- -1L
          .time0 <- NA
          .timenow <- NA
          .firstUpdate <- T
          nw <- nchar(char, "w")
          if (is.na(width)) {
            width <- getOption("width")
            if (style == 3L)
              width <- width - 10L
            width <- trunc(width/nw)
          }
          if (max <= min)
            stop("must have 'max' > 'min'")
          up1 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb < nb) {
              cat(base::paste(rep.int(char, nb - .nb), collapse = ""),
                  file = file)
              utils::flush.console()
            }
            else if (.nb > nb) {
              cat("\r", base::paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            .nb <<- nb
          }
          up2 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb <= nb) {
              cat("\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            else {
              cat("\r", base::paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", base::paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              utils::flush.console()
            }
            .nb <<- nb
          }
          up3 <- function(value, calledOnCreation = F) {
            timenow <- proc.time()[["elapsed"]]
            if (!calledOnCreation && .firstUpdate) {
              .time0 <<- timenow
              .timenow <<- timenow
              .firstUpdate <<- F
            }
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            pc <- round(100 * (value - min)/(max - min))
            if (nb == .nb && pc == .pc && timenow - .timenow < 1)
              return()
            .timenow <<- timenow
            span <- timenow - .time0
            timeXiter <- span/(.val - min)
            ETA <- (max - .val) * timeXiter
            ETAstr <- formatTime(ETA)
            cat(base::paste(c("\r  |", rep.int(" ", nw * width + 6)),
                      collapse = ""), file = file)
            cat(base::paste(c("\r  |", rep.int(char, nb), rep.int(" ",
                                                            nw * (width - nb)), sprintf("| %3d%%", pc), ", ETA ",
                        ETAstr), collapse = ""), file = file)
            utils::flush.console()
            .nb <<- nb
            .pc <<- pc
          }
          getVal <- function() .val
          kill <- function() if (!.killed) {
            cat("\n", file = file)
            utils::flush.console()
            .killed <<- TRUE
          }
          up <- switch(style, up1, up2, up3)
          up(initial, T)
          structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
        }

        unique_proteins <- sort(unique(features$Protein[which(!grepl(";|,",features$Protein))]))
        unique_proteins <- unique_proteins[which(unique_proteins != "")]

        LFQ_protein_quant <- base::as.data.frame(matrix(ncol=ncol(features_quant),nrow=length(unique_proteins),0))
        colnames(LFQ_protein_quant) <- colnames(features_quant)
        rownames(LFQ_protein_quant) <- unique_proteins

        ###Perform LFQ quantification
        cl <- snow::makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        doSNOW::registerDoSNOW(cl)
        iterations <- nrow(LFQ_protein_quant)
        progress <- function(n) utils::setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        start <- Sys.time()
        print(base::paste(label," (",Sys.time(),")",sep=""))
        pb <- txtProgressBar(max = iterations, style = 3)
        res_LFQ <- foreach::foreach(i=1:nrow(LFQ_protein_quant),.options.snow = opts) %dopar%
          {
            sel <- which(features$Protein == unique_proteins[i])
            testdat = features_quant[sel,]
            if(nrow(testdat) > TopN) ###if more than 5 peptides are available select TopN peptides according to intensity over all samples (select highest abundant peptides as quantifications are more accurate)
            {
              total_abundance <- rowSums(testdat,na.rm=T)
              testdat <- testdat[order(total_abundance,decreasing = T)[1:TopN],]
            }
            if(nrow(testdat) >= 2)
            {
              res <- calculate_LFQ(peptide_quant_data = testdat,min_num_ratios = 2,num_ratio_samples=num_ratio_samples)
            }else
            {
              res <- as.numeric(rep(NA,ncol(features_quant)))
              names(res) <- base::paste("V",1:length(res),sep="")
            }

            return(res)
          }
        snow::stopCluster(cl)
        close(pb)

        #Combine results
        for(i in 1:length(res_LFQ))
        {
          if(length(res_LFQ[[i]])>0)
          {
            data.table::set(LFQ_protein_quant,as.integer(i),as.integer(1:ncol(LFQ_protein_quant)),as.list(as.numeric(res_LFQ[[i]])))
          }
        }
        #add information about number of quant features per protein
        count_quant_features <- plyr::count(features$Protein[which(!grepl(";|\\||,",features$Protein))])

        LFQ_protein_quant <- base::data.frame(num_quant_features=count_quant_features$freq[match(rownames(LFQ_protein_quant),count_quant_features$x)],LFQ_protein_quant)

        end <- Sys.time()
        print(base::paste("Finished LFQ-quantification (",Sys.time(),")",sep=""))
        print(end - start)
        #replace 0 by NA
        LFQ_protein_quant[LFQ_protein_quant==0] <- NA

        return(LFQ_protein_quant)
      }

      if(ncol(feature_with_background_intensity) <= 10)num_ratio_samples <- NA
      if(ncol(feature_with_background_intensity) > 10)num_ratio_samples <- 6

      if(calc_peptide_LFQ == F)
      {
        LFQ_quant_with_background <- LFQ_protein_quant_process(features,feature_with_background_intensity,n_cores,label="Perform LFQ-quantification for non-imputed data",num_ratio_samples = num_ratio_samples)
        LFQ_quant_with_background_imputed <- LFQ_protein_quant_process(features,feature_with_background_intensity_imputed,n_cores,label="Perform LFQ-quantification for imputed data",num_ratio_samples = num_ratio_samples)
      }else ###Use peptide LFQ for calcualting protein LFQ
      {
        LFQ_quant_with_background <- LFQ_protein_quant_process(LFQ_peptide_quant_with_background[,c(1:4)],LFQ_peptide_quant_with_background[,c(5:ncol(LFQ_peptide_quant_with_background))],n_cores,label="Perform LFQ-quantification for non-imputed data",num_ratio_samples = num_ratio_samples)
        LFQ_quant_with_background_imputed <- LFQ_protein_quant_process(LFQ_peptide_quant_with_background_imputed[,c(1:4)],LFQ_peptide_quant_with_background_imputed[,c(5:ncol(LFQ_peptide_quant_with_background_imputed))],n_cores,label="Perform LFQ-quantification for imputed data",num_ratio_samples = num_ratio_samples)
      }
    }

    #Perform Top3 and Total quantification
    cl <- parallel::makeCluster(ifelse(n_cores < 4,n_cores,4))#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=1:4) %dopar%
      {
        if(i == 1)
        {
          ###Perform Top3 protein quantification
          res <- Top3_Protein_Quant(features = features,feature_sample_matrix_requantified = feature_with_background_intensity,Quant_pvals=pval_signal_with_background_quant,S2B=S2B,Alignment_scores = alignment_scores_peaks_correct,use_overlapping = T,min_peps = 1)
        }
        if(i == 2)
        {
          ###Perform Total protein quantification
          res <- Total_Protein_Quant(features = features,feature_sample_matrix_requantified = feature_with_background_intensity,Quant_pvals=pval_signal_with_background_quant,S2B=S2B,Alignment_scores = alignment_scores_peaks_correct,use_overlapping = T,min_peps = 1)
        }
        if(i == 3)
        {
          ###Perform Top3 protein quantification for imputed data
          res <- Top3_Protein_Quant(features = features,feature_sample_matrix_requantified = feature_with_background_intensity_imputed,Quant_pvals=pval_signal_with_background_quant,Alignment_scores = alignment_scores_peaks_correct,S2B=S2B,use_overlapping = T,min_peps = 1)
        }
        if(i == 4)
        {
          ###Perform Total protein quantification for imputed data
          res <- Total_Protein_Quant(features = features,feature_sample_matrix_requantified = feature_with_background_intensity_imputed,Quant_pvals=pval_signal_with_background_quant,Alignment_scores = alignment_scores_peaks_correct,S2B=S2B,use_overlapping = T,min_peps = 1)
        }
        return(res)
      }
    parallel::stopCluster(cl)

    Top3_quant_with_background <- res[[1]]
    Total_quant_with_background <- res[[2]]
    Top3_quant_with_background_imputed <- res[[3]]
    Total_quant_with_background_imputed <- res[[4]]

    ###Match Gene names to Uniprot Identifier
    if(any(colnames(MaxQ_protein_groups) == "Gene.names") & any(colnames(MaxQ_protein_groups) == "Protein.IDs"))
    {
      temp <- MaxQ_protein_groups[,c("Gene.names","Protein.IDs")]
      temp <- temp[which(temp$Gene.names != ""),]

      UniProt_to_GeneName <- base::data.frame(UniProt_ID=unique(as.character(stringr::str_split(temp$Protein.IDs,";",simplify = T))),Gene_Name="")
      UniProt_to_GeneName$Gene_Name <- as.character(UniProt_to_GeneName$Gene_Name)

      for(i in 1:nrow(UniProt_to_GeneName))
      {
        UniProt_to_GeneName$Gene_Name[i] <- as.character(temp$Gene.names[which(grepl(UniProt_to_GeneName$UniProt_ID[i],temp$Protein.IDs))])
      }

      Top3_quant_with_background <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Top3_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Top3_quant_with_background),Top3_quant_with_background)
      rownames(Top3_quant_with_background) <- c()

      Total_quant_with_background <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Total_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Total_quant_with_background),Total_quant_with_background)
      rownames(Total_quant_with_background) <- c()

      Top3_quant_with_background_imputed <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Top3_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Top3_quant_with_background_imputed),Top3_quant_with_background_imputed)
      rownames(Top3_quant_with_background_imputed) <- c()

      Total_quant_with_background_imputed <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Total_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Total_quant_with_background_imputed),Total_quant_with_background_imputed)
      rownames(Total_quant_with_background_imputed) <- c()

      if(calc_peptide_LFQ == T)
      {
        LFQ_peptide_quant_with_background <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(LFQ_peptide_quant_with_background$Protein,UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=LFQ_peptide_quant_with_background$Protein,LFQ_peptide_quant_with_background[,c(1,2,4,5:ncol(LFQ_peptide_quant_with_background))])
        rownames(LFQ_peptide_quant_with_background) <- c()
        colnames(LFQ_peptide_quant_with_background) <- base::gsub("^X","",colnames(LFQ_peptide_quant_with_background))

        LFQ_peptide_quant_with_background_imputed <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(LFQ_peptide_quant_with_background_imputed$Protein,UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=LFQ_peptide_quant_with_background_imputed$Protein,LFQ_peptide_quant_with_background_imputed[,c(1,2,4,5:ncol(LFQ_peptide_quant_with_background_imputed))])
        rownames(LFQ_peptide_quant_with_background_imputed) <- c()
        colnames(LFQ_peptide_quant_with_background_imputed) <- base::gsub("^X","",colnames(LFQ_peptide_quant_with_background_imputed))
      }

      if(calc_protein_LFQ == T)
      {
        LFQ_quant_with_background <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(LFQ_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(LFQ_quant_with_background),LFQ_quant_with_background)
        rownames(LFQ_quant_with_background) <- c()
        colnames(LFQ_quant_with_background) <- base::gsub("^X","",colnames(LFQ_quant_with_background))

        LFQ_quant_with_background_imputed <- base::data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(LFQ_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(LFQ_quant_with_background_imputed),LFQ_quant_with_background_imputed)
        rownames(LFQ_quant_with_background_imputed) <- c()
        colnames(LFQ_quant_with_background_imputed) <- base::gsub("^X","",colnames(LFQ_quant_with_background_imputed))
      }

    }

    if(calc_protein_LFQ == T)
    {
      utils::write.table(x = LFQ_quant_with_background,file = base::paste(path_to_features,"/Proteins_quantification_LFQ",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
      utils::write.table(x = LFQ_quant_with_background_imputed,file = base::paste(path_to_features,"/Proteins_quantification_LFQ_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    }
    utils::write.table(x = Top3_quant_with_background,file = base::paste(path_to_features,"/Proteins_quantification_Top3",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Top3_quant_with_background_imputed,file = base::paste(path_to_features,"/Proteins_quantification_Top3_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Total_quant_with_background,file = base::paste(path_to_features,"/Proteins_quantification_Total",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Total_quant_with_background_imputed,file = base::paste(path_to_features,"/Proteins_quantification_Total_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")

    if(calc_peptide_LFQ == T)
    {
      utils::write.table(x = LFQ_peptide_quant_with_background,file = base::paste(path_to_features,"/Peptides_quantification_LFQ",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
      utils::write.table(x = LFQ_peptide_quant_with_background_imputed,file = base::paste(path_to_features,"/Peptides_quantification_LFQ_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    }
  }

  setwd(path_to_features)

  ###finally save feature level quantification
  utils::write.table(x = features,file = base::paste(path_to_features,"/Features",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
  utils::write.table(x = feature_with_background_intensity,file = base::paste(path_to_features,"/Features_quantification",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = feature_with_background_intensity_imputed,file = base::paste(path_to_features,"/Features_quantification_imputed",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = pval_signal_with_background_quant,file = base::paste(path_to_features,"/Features_quantification_pvals",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = Ioncount_feature_with_background_intensity,file = base::paste(path_to_features,"/Features_quantification_ioncount",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = S2B,file = base::paste(path_to_features,"/Features_quantification_S2B",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = alignment_variability_score,file = base::paste(path_to_features,"/Features_quantification_variability_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = alignment_scores_peaks_correct,file = base::paste(path_to_features,"/Features_quantification_alignment_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = mono_iso_alignment_summary,file = base::paste(path_to_features,"/Features_quantification_mono_iso_alignment_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")

  save(QC_data,file = base::paste("Temporary_files/Feature_quantification_QC_data.RData",sep=""))
  options(warn=0)
}


#' Peak decision algorithm only designed for internal use
#' @param features_select Features on which peak-selection should be performed
#' @param peak_quant List of determined peaks per feature
#' @param samples String vector of all samples
#' @param s Integer indicating index of current sample
#' @param RT_correction_factors RT corrections
#' @param mz_correction_factors mz corrections
#' @param features_intensity_sample features_intensity_sample
#' @param Ioncount_sample Ioncount_sample
#' @param feature_with_background_intensity_sample feature_with_background_intensity_sample
#' @param Ioncount_with_background_sample Ioncount_with_background_sample
#' @param peak_selected_sample peak_selected_sample
#' @param delta_mz delta_mz
#' @param delta_rt delta_rt
#' @param peak_min_ion_count peak_min_ion_count
#' @param chunk chunk
#' @param num_chunks num_chunks
#' @param progress Progressbar
#' @import data.table
#' @export
#' @details Peak decision algorithm
peak_decision <- function(features_select,peak_quant,samples,s,RT_correction_factors,mz_correction_factors,features_intensity_sample,Ioncount_sample,feature_with_background_intensity_sample,Ioncount_with_background_sample,peak_selected_sample,delta_mz,delta_rt,peak_min_ion_count,chunk=NULL,num_chunks=NULL,progress=T)
{
  #prevent issues during R CMD check
  ..s <- s
  rm(..s)
  ..known_peaks_indices <- 1
  rm(..known_peaks_indices)
  ..o <- 1
  rm(..o)

  #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
  #suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
  if(progress == T)
  {
    pb <- tcltk::tkProgressBar(title = base::paste("Prepare for peak selection -",samples[s]),label=base::paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    close(pb)
  }

  ###next, decide for all other feature quantifications for which at least 1 peak was available which peak quantification should be used
  max <- nrow(features_select)

  if(progress == T)
  {
    if(!is.null(chunk) & !is.null(num_chunks))
    {
      pb <- tcltk::tkProgressBar(title = base::paste("Select peaks (",chunk,"/",num_chunks,") - ",samples[s],sep=""),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    }else
    {
      pb <- tcltk::tkProgressBar(title = base::paste("Select peaks -",samples[s]),label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    }
  }

  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0

  for(i in 1:nrow(features_select))
  {
    num_detected_peaks <- ifelse(!is.na(peak_quant$Standard$num_peaks_with_background[i,..s]),as.numeric(peak_quant$Standard$num_peaks_with_background[i,..s]),
                                 ifelse(!is.na(peak_quant$Peak_1$num_peaks_with_background[i,..s]),as.numeric(peak_quant$Peak_1$num_peaks_with_background[i,..s]),NA))
    if(!is.na(num_detected_peaks)) ##peaks detected?
    {
      ###check if more than 1 peak was detected and true peak is unknown
      if(num_detected_peaks > 1)
      {
        # if(!is.na(peak_quant$Peak_1$correct_peak_with_background[i,..s]) & peak_quant$Peak_1$correct_peak_with_background[i,..s] == 1) ###correct peak is known so take the closest peak
        # {
        #   ###already done in batch before
        #   ###so nothing to be done here
        #
        #   # data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$features_intensity[i,..s]))
        #   # data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$Ioncount_feature_sample_matrix[i,..s]))
        #   # data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$feature_with_background_intensity[i,..s]))
        #   # data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$Ioncount_feature_with_background_intensity[i,..s]))
        #   #
        #   # data.table::set(dT1_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$dT1[i,..s]))
        #   # data.table::set(dM1_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$dM1[i,..s]))
        #
        # }
        if(peak_quant$Peak_1$correct_peak_with_background[i,..s] == 0) ###more than 1 peaks but true peak not known
        {
          ###check if for any sample the correct peak is known
          num_known_peaks <- length(which(peak_quant$Peak_1$correct_peak_with_background[i,] != 0))
          if(num_known_peaks > 0)
          {
            ###now check where correct known peak usually is located in samples were true location is known
            known_peaks_indices <- which(peak_quant$Peak_1$correct_peak_with_background[i,] == 1)

            RT_in_known <- as.numeric(peak_quant$Peak_1$Peak_rt_with_background[i,..known_peaks_indices])
            mz_in_known <- as.numeric(peak_quant$Peak_1$Peak_mz_with_background[i,..known_peaks_indices])

            ###get RT and mz for detected peaks in current sample
            detected_peaks_rt <- vector("numeric",as.numeric(num_detected_peaks))
            detected_peaks_mz <- vector("numeric",as.numeric(num_detected_peaks))

            names(detected_peaks_rt) <- 1:as.numeric(num_detected_peaks)
            names(detected_peaks_mz) <- 1:as.numeric(num_detected_peaks)

            for(p in 1:length(detected_peaks_rt))
            {
              detected_peaks_rt[p] <- as.numeric(peak_quant[[p]]$Peak_rt_with_background[i,..s])
              detected_peaks_mz[p] <- as.numeric(peak_quant[[p]]$Peak_mz_with_background[i,..s])
            }

            if(length(which(!is.na(detected_peaks_rt) & !is.na(detected_peaks_mz))) > 0)
            {
              detected_peaks_rt <- detected_peaks_rt[which(!is.na(detected_peaks_rt))]
              detected_peaks_mz <- detected_peaks_mz[which(!is.na(detected_peaks_mz))]

              ###now compare observed deviations in RT and mz for known peaks with the deviation for the potential peaks in current sample
              ###for this comparison we have to correct detected peaks RT and m/z with the RT and mz correction factors per corresponding sample
              ###otherwise it would be possible that current sample shows e.g. a global RT shift which could result in a major deviation from all other known samples
              ###thus a peak wouldn?t be selected although it e.g. perfectly lies close to the expected window
              RT_in_known_corrected <- as.numeric(RT_in_known - RT_correction_factors[i,known_peaks_indices])
              mz_in_known_corrected <- as.numeric(mz_in_known - mz_correction_factors[i,known_peaks_indices])

              detected_peaks_rt_corrected <- as.numeric(detected_peaks_rt - RT_correction_factors[i,s])
              detected_peaks_mz_corrected <- as.numeric(detected_peaks_mz - mz_correction_factors[i,s])
              names(detected_peaks_rt_corrected) <- names(detected_peaks_rt)
              names(detected_peaks_mz_corrected) <- names(detected_peaks_mz)
              ###check if the detected RT and mz of peaks are likely to be belonging to the same population (compared to known peaks)
              ###99 % of observed data points lie within 2*sd range --> if delta RT or delta mZ > 2*sd then most likely this is the wrong peak
              ###so we assume all peaks within this RT and mz deviation to be valid peak candidates
              ###if these deviations are smaller then the expected RT and mz window, we expand these acceptance criteria accordingly
              delta_RT_cut <- 3*stats::sd(RT_in_known_corrected,na.rm=T)
              if(is.na(delta_RT_cut) | delta_RT_cut < 2*delta_rt)delta_RT_cut <- 2*delta_rt
              delta_mz_cut <- 3*stats::sd(mz_in_known_corrected,na.rm=T)
              if(is.na(delta_mz_cut) | delta_mz_cut < 3*delta_mz)delta_mz_cut <- 3*delta_mz

              within_range <- ifelse(abs(detected_peaks_rt_corrected-mean(RT_in_known_corrected,na.rm=T)) <= delta_RT_cut &
                                       abs(detected_peaks_mz_corrected-mean(mz_in_known_corrected,na.rm=T)) <= delta_mz_cut,T,F)

              if(any(within_range))
              {
                detected_peaks_rt <- detected_peaks_rt[within_range]
                detected_peaks_mz <- detected_peaks_mz[within_range]
                detected_peaks_rt_corrected <- detected_peaks_rt_corrected[within_range]
                detected_peaks_mz_corrected <- detected_peaks_mz_corrected[within_range]

                ###now check that at the position of the potential peak no other peak is present in samples with known peak
                count_other_peaks <- sum(peak_quant$Peak_1$num_peaks_with_background[i,..known_peaks_indices]-1)
                if(count_other_peaks > 0)
                {
                  other_peaks <- base::as.data.frame(matrix(ncol=6,nrow=count_other_peaks,0))
                  counter <- 0

                  for(o in known_peaks_indices) ###collect peak data for all additional peaks in samples were peak was known
                  {
                    cur_peak_count <- as.numeric(peak_quant$Peak_1$num_peaks_with_background[i,..o]-1)
                    if(cur_peak_count > 0)
                    {
                      for(p in 1:cur_peak_count)
                      {
                        counter <- counter + 1
                        temp <- c(as.numeric(peak_quant[[p+1]]$Peak_mz_with_background[i,..o] - mz_correction_factors[i,o]),
                                  as.numeric(peak_quant[[p+1]]$Peak_rt_with_background[i,..o] - RT_correction_factors[i,o]),
                                  as.numeric(peak_quant[[p+1]]$feature_with_background_intensity[i,..o]),
                                  as.numeric(peak_quant[[p+1]]$Ioncount_feature_with_background_intensity[i,..o]))
                        ##Add distance in RT and mz to known peak
                        temp[5:6] <- c(abs(temp[1] - as.numeric(peak_quant[[1]]$Peak_mz_with_background[i,..o] - mz_correction_factors[i,o])),
                                       abs(temp[2] - as.numeric(peak_quant[[1]]$Peak_rt_with_background[i,..o] - RT_correction_factors[i,o])))

                        data.table::set(other_peaks,as.integer(counter),as.integer(1:6),value=list(temp[1],temp[2],temp[3],temp[4],temp[5],temp[6]))
                      }
                    }
                  }
                  #disregard other peaks below significance threshold
                  other_peaks <- other_peaks[which(other_peaks$V4 > peak_min_ion_count),]
                  #disregard other peaks which are too close to the expected peak
                  #cutoff RT: > 1.1* half peak width -> half RT extraction window
                  #cutoff mz: > 1.1* delta_mz -> half mz extraction window
                  other_peaks <- other_peaks[which(other_peaks$V5 > delta_mz*1.1 | other_peaks$V6 > ((features_select$RT_length[i]/2)*1.1)),]
                  ###now check that selected peaks are not close to other peaks in sampels with known peak
                  overlap <- vector("logical",length(detected_peaks_rt))
                  for(o in 1:length(detected_peaks_rt))
                  {
                    overlap[o] <- ifelse(any(abs(detected_peaks_rt_corrected[o]-other_peaks$V2) <= (delta_rt/2) & abs(detected_peaks_mz_corrected[o]-other_peaks$V1) <= delta_mz/2,na.rm=T),T,F)
                  }
                }else
                {
                  overlap <- vector("logical",length(detected_peaks_rt))
                }

                if(any(overlap == F)) ###not overlapping with any other peak
                {
                  detected_peaks_rt <- detected_peaks_rt_corrected[which(overlap == F)]
                  detected_peaks_mz <- detected_peaks_mz_corrected[which(overlap == F)]
                  RT_in_known_other_peaks <- as.numeric(peak_quant$Peak_1$Peak_rt_with_background[i,..known_peaks_indices])
                  mz_in_known_other_peaks <- as.numeric(peak_quant$Peak_1$Peak_mz_with_background[i,..known_peaks_indices])

                  ###select the peak which is closest to all other known peaks

                  sum_delta_to_known_peaks_rt <- vector("numeric",length(detected_peaks_rt))
                  sum_delta_to_known_peaks_mz <- vector("numeric",length(detected_peaks_rt))
                  sum_delta_to_known_peaks <- vector("numeric",length(detected_peaks_rt))
                  for(p in 1:length(detected_peaks_rt))
                  {
                    sum_delta_to_known_peaks_rt[p] <- sum(detected_peaks_rt[p]-RT_in_known_corrected,na.rm=T)/length(RT_in_known_corrected)
                    sum_delta_to_known_peaks_mz[p] <- sum((detected_peaks_mz[p]-mz_in_known_corrected)*500,na.rm=T)/length(RT_in_known_corrected)
                    sum_delta_to_known_peaks[p] <- abs(sum_delta_to_known_peaks_rt[p])+abs(sum_delta_to_known_peaks_mz[p])
                  }

                  selected_peak <- which(!is.na(sum_delta_to_known_peaks) & sum_delta_to_known_peaks == min(sum_delta_to_known_peaks[which(!is.na(sum_delta_to_known_peaks))],na.rm=T))[1]
                  selected_peak <- as.numeric(names(detected_peaks_rt)[selected_peak])
                  data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$features_intensity[i,..s]))
                  data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_sample_matrix[i,..s]))

                  data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$feature_with_background_intensity[i,..s]))
                  data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_with_background_intensity[i,..s]))

                  data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(selected_peak))

                }else###none of the peaks is not overlapping with an other peak detected in other samples where correct peak was known
                {
                  data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                  data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                  data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                  data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                  data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                }

              }else ###none of the peaks is within the accepted range thus use standard window
              {
                data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
              }
            }else ###none of the peaks is valid
            {
              data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
              data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

              data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
              data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

              data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
            }


          }else ###in no other sample the correct peak is known thus take total window
          {
            data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
            data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

            data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
            data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

            data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
          }
        }
      }else if(num_detected_peaks == 1) ###one peak was detected
      {
        ### check if correct peak is unknown
        if(peak_quant$Peak_1$correct_peak_with_background[i,..s] == 0)
        {
          ###if this is the case check if detected peak is comparable to correct peaks in other samples where correct peak was known
          ###check if for any sample the correct peak is known
          num_known_peaks <- length(which(peak_quant$Peak_1$correct_peak_with_background[i,] != 0))
          if(num_known_peaks > 0)
          {
            ###now check where correct known peak usually is located in samples were true location is known
            known_peaks_indices <- which(peak_quant$Peak_1$correct_peak_with_background[i,] == 1)

            RT_in_known <- as.numeric(peak_quant$Peak_1$Peak_rt_with_background[i,..known_peaks_indices])
            mz_in_known <- as.numeric(peak_quant$Peak_1$Peak_mz_with_background[i,..known_peaks_indices])

            detected_peaks_rt <- as.numeric(peak_quant[[1]]$Peak_rt_with_background[i,..s])
            detected_peaks_mz <- as.numeric(peak_quant[[1]]$Peak_mz_with_background[i,..s])
            if(length(which(!is.na(detected_peaks_rt) & !is.na(detected_peaks_mz))) > 0)
            {
              detected_peaks_rt <- detected_peaks_rt[which(!is.na(detected_peaks_rt))]
              detected_peaks_mz <- detected_peaks_mz[which(!is.na(detected_peaks_mz))]

              ###now compare observed deviations in RT and mz for known peaks with the deviation for the potential peaks in current sample
              ###for this comparison we have to correct detected peaks RT and m/z with the RT and mz correction factors per corresponding sample
              ###otherwise it would be possible that current sample shows e.g. a global RT shift which could result in a major deviation from all other known samples
              ###thus a peak wouldn?t be selected although it e.g. perfectly lies close to the expected window
              RT_in_known_corrected <- as.numeric(RT_in_known - RT_correction_factors[i,known_peaks_indices])
              mz_in_known_corrected <- as.numeric(mz_in_known - mz_correction_factors[i,known_peaks_indices])

              detected_peaks_rt_corrected <- as.numeric(detected_peaks_rt - RT_correction_factors[i,s])
              detected_peaks_mz_corrected <- as.numeric(detected_peaks_mz - mz_correction_factors[i,s])

              ###check if the detected RT and mz of peaks are likely to be belonging to the same population (compared to known peaks)
              ###99 % of observed data points lie within 3*sd range --> if delta RT or delta mZ > 2*sd then most likely this is the wrong peak
              ###so we assume all peaks within this RT and mz deviation to be valid peak candidates
              ###if these deviations are smaller then the expected RT and mz window, we expand these acceptance criteria accordingly
              delta_RT_cut <- 3*stats::sd(RT_in_known,na.rm=T)
              if(is.na(delta_RT_cut) | delta_RT_cut < 2*delta_rt)delta_RT_cut <- 2*delta_rt
              delta_mz_cut <- 3*stats::sd(mz_in_known,na.rm=T)
              if(is.na(delta_mz_cut) | delta_mz_cut < 3*delta_mz)delta_mz_cut <- 3*delta_mz

              within_range <- ifelse(abs(detected_peaks_rt_corrected-mean(RT_in_known_corrected,na.rm=T)) <= delta_RT_cut &
                                       abs(detected_peaks_mz_corrected-mean(mz_in_known_corrected,na.rm=T)) <= delta_mz_cut,T,F)

              if(any(within_range))
              {
                detected_peaks_rt <- detected_peaks_rt[within_range]
                detected_peaks_mz <- detected_peaks_mz[within_range]
                detected_peaks_rt_corrected <- detected_peaks_rt_corrected[within_range]
                detected_peaks_mz_corrected <- detected_peaks_mz_corrected[within_range]

                ###now check that at the position of the potential peak no other peak is present in samples with known peak
                count_other_peaks <- sum(peak_quant$Peak_1$num_peaks_with_background[i,..known_peaks_indices]-1)
                if(count_other_peaks > 0)
                {
                  other_peaks <- base::as.data.frame(matrix(ncol=6,nrow=count_other_peaks,0))
                  counter <- 0

                  for(o in known_peaks_indices) ###collect peak data for all additional peaks in samples were peak was known
                  {
                    cur_peak_count <- as.numeric(peak_quant$Peak_1$num_peaks_with_background[i,..o]-1)
                    if(cur_peak_count > 0)
                    {
                      for(p in 1:cur_peak_count)
                      {
                        counter <- counter + 1
                        temp <- c(as.numeric(peak_quant[[p+1]]$Peak_mz_with_background[i,..o] - mz_correction_factors[i,o]),
                                  as.numeric(peak_quant[[p+1]]$Peak_rt_with_background[i,..o] - RT_correction_factors[i,o]),
                                  as.numeric(peak_quant[[p+1]]$feature_with_background_intensity[i,..o]),
                                  as.numeric(peak_quant[[p+1]]$Ioncount_feature_with_background_intensity[i,..o]))
                        ##Add distance in RT and mz to known peak
                        temp[5:6] <- c(abs(temp[1] - as.numeric(peak_quant[[1]]$Peak_mz_with_background[i,..o] - mz_correction_factors[i,o])),
                                       abs(temp[2] - as.numeric(peak_quant[[1]]$Peak_rt_with_background[i,..o] - RT_correction_factors[i,o])))
                      }
                    }
                  }
                  #disregard other peaks below significance threshold
                  other_peaks <- other_peaks[which(other_peaks$V4 > peak_min_ion_count),]
                  #disregard other peaks which are too close to the expected peak
                  #cutoff RT: > 1.1* half peak width -> half RT extraction window
                  #cutoff mz: > 1.1* delta_mz -> half mz extraction window
                  other_peaks <- other_peaks[which(other_peaks$V5 > delta_mz*1.1 | other_peaks$V6 > ((features_select$RT_length[i]/2)*1.1)),]
                  ###now check that selected peaks are not close to other peaks in sampels with known peak

                  overlap <- ifelse(any(abs(detected_peaks_rt_corrected-other_peaks$V2) <= delta_rt/2 & abs(detected_peaks_mz_corrected-other_peaks$V1) <= delta_mz/2,na.rm=T),T,F)
                }else
                {
                  overlap <- F
                }

                if(any(overlap == F)) ###not overlapping with any other peak
                {
                  selected_peak <- 1

                  data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$features_intensity[i,..s]))
                  data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_sample_matrix[i,..s]))

                  data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$feature_with_background_intensity[i,..s]))
                  data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_with_background_intensity[i,..s]))

                  data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(selected_peak))
                }else###the peak is overlapping with an other peak detected in other samples where correct peak was known
                {
                  data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                  data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                  data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                  data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                  data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                }
              }else ###none of the peaks is within the accepted range thus use standard window
              {
                data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
              }
            }else ###no peak was valid
            {
              data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
              data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

              data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
              data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

              data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
            }
          }else ###in no other sample the correct peak was known thus use the standard window
          {
            data.table::set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
            data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

            data.table::set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
            data.table::set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

            data.table::set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
          }


        }


      }
    }

    updatecounter <- updatecounter + 1
    if(updatecounter >= 10 & progress == T)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- lubridate::seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      tcltk::setTkProgressBar(pb, i, label=base::paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    }

  }
  if(progress == T)
  {
    close(pb)
  }

  return(list(features_intensity_sample=features_intensity_sample,
              Ioncount_sample=Ioncount_sample,
              feature_with_background_intensity_sample=feature_with_background_intensity_sample,
              Ioncount_with_background_sample=Ioncount_with_background_sample,
              peak_selected_sample=peak_selected_sample))
}

#' Peak selection FDR algorithm for internal use
#' @param num_features num_features
#' @param features features
#' @param samples samples
#' @param peaks peaks
#' @param path_to_features path_to_features
#' @param peak_quant peak_quant
#' @param feature_with_background_intensity feature_with_background_intensity
#' @param peak_selected peak_selected
#' @param delta_mz delta_mz
#' @param delta_rt delta_rt
#' @param peak_min_ion_count peak_min_ion_count
#' @param num_peaks_store num_peaks_store
#' @param alignment_scores_cutoff alignment_scores_cutoff
#' @param n_cores n_cores
#' @param peak_decision peak_decision
#' @param Alignment_scoring Alignment_scoring
#' @param seed seed
#' @param plot plot
#' @import data.table
#' @import foreach
#' @import randomForest
#' @export
#' @details Peak selection FDR algorithm for internal use
Peak_selection_FDR <- function(num_features=500,features,samples,peaks,path_to_features,peak_quant,feature_with_background_intensity,peak_selected,delta_mz,delta_rt,peak_min_ion_count,num_peaks_store,alignment_scores_cutoff,n_cores,peak_decision,Alignment_scoring,seed=110519,plot=T)
{
  if(length(which(peaks$known == 1)) >= 10)
  {
    #suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    #suppressWarnings(suppressMessages(library(randomForest,quietly = T)))
    #suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
    #select features for FDR estimation per sample
    if(!is.na(seed))set.seed(seed)
    temp_peaks <- peaks[which(peaks$peak == peaks$selected & peaks$known == 1 & !grepl("_i|_d",peaks$feature)),]
    selection_feature <- sample(unique(temp_peaks$feature),length(unique(temp_peaks$feature)),replace = F)
    selection_feature <- match(selection_feature,features$Feature_name)
    selection_feature_per_sample <- list()

    for(s in 1:length(samples))
    {
      temp_peaks_2 <- temp_peaks[which(temp_peaks$sample == samples[s]),]
      temp_peaks_2 <- temp_peaks_2[match(features$Feature_name[selection_feature],temp_peaks_2$feature),]
      sel <- which(!is.na(temp_peaks_2$known))[1:num_features]
      if(any(is.na(sel)))sel <- sel[which(!is.na(sel))]
      selection_feature_per_sample[[s]] <- selection_feature[sel]
    }

    #remove known peak location for selected features and samples and replace with predicted corrections based on models
    features_select_FDR <- list()
    features$Observed_mz <- as.character(features$Observed_mz)
    features$Observed_RT <- as.character(features$Observed_RT)

    setwd(base::paste(path_to_features,"/Temporary_files",sep=""))
    e <- new.env()
    load("Feature_alignment_QC_data.RData",envir = e)
    median_feature_properties <- e$QC_data$mz_calibration_median_feature_properties
    mz_correction_models <- e$QC_data$mz_calibration_models
    RT_alignment_GAM_models <- e$QC_data$RT_calibration_GAM_models

    #determine multiplicity
    if(any(grepl("_Channel_light|_Channel_medium|_Channel_heavy",samples)))
    {
      if(any(grepl("_Channel_medium",samples)))
      {
        multiplicity <- 3
      }else
      {
        multiplicity <- 2
      }
    }else
    {
      multiplicity <- 1
    }

    for(s in 1:length(samples))
    {
      temp <- data.table::copy(features)
      temp <- temp[selection_feature_per_sample[[s]],]
      if(nrow(temp)>10)
      {
        #remove observed RT and mz
        obs_rt <- (stringr::str_split(temp$Observed_RT,";",simplify = T))
        obs_rt[,s] <- "NA"
        obs_mz <- (stringr::str_split(temp$Observed_mz,";",simplify = T))
        obs_mz[,s] <- "NA"
        temp$Observed_mz <- apply(obs_mz,1,base::paste,collapse = ";")
        #exchange mz calibration with prediction from model
        select_model <- which(names(mz_correction_models) == samples[s])

        features_select <- selection_feature_per_sample[[s]]

        temp_data <- median_feature_properties[match(features_select,median_feature_properties$Feature),]

        temp_data$Resolution <- 0

        temp_data <- temp_data[which(rowSums(is.na(temp_data[,c("Retention.time","m.z","Charge","Resolution")])) == 0),]

        if(length(select_model)>0)
        {
          if(nrow(temp_data) > 0)
          {
            prediction <- stats::predict(mz_correction_models[[select_model]], temp_data[,-1], type = "response")
            if(any(is.na(prediction)))prediction[which(is.na(prediction))] <- 0
            #add expected SILAC isotope shifts
            if(multiplicity > 1)
            {
              if(grepl("light",samples[s]))SILAC_label_mz_shift <- features$m.z._shift_light[temp_data$Feature]
              if(grepl("medium",samples[s]))SILAC_label_mz_shift <- features$m.z._shift_medium[temp_data$Feature]
              if(grepl("heavy",samples[s]))SILAC_label_mz_shift <- features$m.z._shift_heavy[temp_data$Feature]
              prediction <- prediction + SILAC_label_mz_shift
            }
            data.table::set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(colnames(temp) == base::paste("mz_calibration.",samples[s],sep=""))),prediction)
          }else
          {
            data.table::set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(colnames(temp) == base::paste("mz_calibration.",samples[s],sep=""))),0)
          }
        }else
        {
          data.table::set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(grepl("mz_calibration.",colnames(temp)))[s]),0)
        }
        #exchange RT calibration with prediction from model
        select_model <- s
        if(!is.null(RT_alignment_GAM_models[[select_model]]))
        {
          prediction <- stats::predict(RT_alignment_GAM_models[[select_model]], base::data.frame(x = temp$RT))
          if(any(is.na(prediction)))prediction[which(is.na(prediction))] <- 0
          data.table::set(temp,as.integer(1:nrow(temp)),as.integer(which(colnames(temp) == base::paste("RT_calibration.",samples[s],sep=""))),prediction)

        }else
        {
          data.table::set(temp,as.integer(1:nrow(temp)),as.integer(which(colnames(temp) == base::paste("RT_calibration.",samples[s],sep=""))),0)
        }

        features_select_FDR[[s]] <- temp
      }else
      {
        features_select_FDR[[s]] <- temp
      }
    }

    #prepare columns containing correction information
    indices_RT_correction <- which(grepl("RT_calibration",colnames(features)))
    indices_mz_correction <- which(grepl("mz_calibration",colnames(features)))
    ordering_indices <- match(base::gsub("-",".",samples),base::gsub("RT_calibration\\.","",colnames(features)[indices_RT_correction]))
    indices_RT_correction <- indices_RT_correction[ordering_indices]
    indices_mz_correction <- indices_mz_correction[ordering_indices]

    #reduce peak_quant information to only cover features selected for FDR calculation

    #prepare only required peak quant data
    peak_quant_reduced <- data.table::copy(peak_quant)
    all_selected_features <- unique(unlist(selection_feature_per_sample))

    for(p in 1:(num_peaks_store+1))
    {
      for(p2 in 1:length(peak_quant_reduced[[p]]))
      {
        peak_quant_reduced[[p]][[p2]] <- peak_quant_reduced[[p]][[p2]][all_selected_features,]
      }
    }

    #suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))

    cl <- snow::makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    doSNOW::registerDoSNOW(cl)

    start <- Sys.time()

    res <- foreach::foreach(s=1:length(samples)) %dopar%
      {
        #suppressWarnings(suppressMessages(library(data.table,quietly = T)))

        selection_feature <- selection_feature_per_sample[[s]]
        if(length(selection_feature)>0)
        {
          max <- 1
          pb <- tcltk::tkProgressBar(title = "Evaluate peak selection FDR",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)

          features_select_FDR_cur <- features_select_FDR[[s]]
          #now perform peak decision for FDR_features and indicate if same peak was selected as before
          peak_decision_same_without_peak_known <- vector("logical",length(selection_feature))

          RT_correction_factors <- features_select_FDR_cur[,indices_RT_correction]
          mz_correction_factors <- features_select_FDR_cur[,indices_mz_correction]
          colnames(RT_correction_factors) <- samples
          rownames(RT_correction_factors) <- features_select_FDR_cur$Feature_name
          colnames(mz_correction_factors) <- samples
          rownames(mz_correction_factors) <- features_select_FDR_cur$Feature_name

          feature_with_background_intensity_select_FDR <- data.table::copy(feature_with_background_intensity)
          peak_selected_select_FDR <- data.table::copy(peak_selected)

          temp_dummy_df <-  data.table::data.table(Name=rep(0L,length(selection_feature)))
          colnames(temp_dummy_df) <- samples[s]
          close(pb)

          peak_quant_temp <- data.table::copy(peak_quant_reduced)

          for(p in 1:(num_peaks_store+1))
          {
            for(p2 in 1:length(peak_quant_temp[[p]]))
            {
              peak_quant_temp[[p]][[p2]] <- peak_quant_temp[[p]][[p2]][match(selection_feature,all_selected_features),]
            }
          }

          peak_quant_temp$Peak_1$correct_peak_with_background[,s] <- 0L

          res_temp <- peak_decision(features_select = features_select_FDR_cur,
                                    peak_quant = peak_quant_temp,
                                    samples = samples,
                                    s = s,
                                    RT_correction_factors = RT_correction_factors,
                                    mz_correction_factors = mz_correction_factors,
                                    features_intensity_sample = temp_dummy_df,
                                    Ioncount_sample = temp_dummy_df,
                                    feature_with_background_intensity_sample = feature_with_background_intensity_select_FDR[selection_feature,s,drop=F],
                                    Ioncount_with_background_sample = temp_dummy_df,
                                    peak_selected_sample = peak_selected_select_FDR[selection_feature,s,drop=F],
                                    delta_mz = delta_mz,
                                    delta_rt = delta_rt,
                                    peak_min_ion_count = peak_min_ion_count,
                                    progress = T)

          peak_decision_same_without_peak_known <- res_temp$peak_selected_sample[,1] == peak_selected[selection_feature,s]

          data.table::set(feature_with_background_intensity_select_FDR,as.integer(selection_feature),as.integer(s),res_temp$feature_with_background_intensity_sample)
          data.table::set(peak_selected_select_FDR,as.integer(selection_feature),as.integer(s),res_temp$peak_selected_sample)

          names(peak_decision_same_without_peak_known) <- base::paste(features_select_FDR_cur$Feature_name,"_Sample_",s,sep="")

          return(list(peak_decision_same_without_peak_known=peak_decision_same_without_peak_known,
                      feature_with_background_intensity_select_FDR=feature_with_background_intensity_select_FDR,
                      peak_selected_select_FDR=peak_selected_select_FDR))
        }else
        {
          return(NULL)
        }

      }
    snow::stopCluster(cl)


    #determine how many false selections show large intensity difference and would not be removed by alignment scoring
    results_peak_selection_FDR_all <- list()
    max <- length(samples)
    pb <- tcltk::tkProgressBar(title="Evaluate peak selection FDR",label=base::paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(s in 1:length(samples))
    {
      if(!is.null(res[[s]]))
      {
        peak_decision_same_without_peak_known <- res[[s]]$peak_decision_same_without_peak_known
        peak_selected_select_FDR <- res[[s]]$peak_selected_select_FDR
        feature_with_background_intensity_select_FDR <- res[[s]]$feature_with_background_intensity_select_FDR
        selection_feature <- selection_feature_per_sample[[s]]
        if(length(which(peak_decision_same_without_peak_known == F))>0L)
        {
          selections_all <- 1:length(peak_decision_same_without_peak_known)

          #compare quantification differences for all tests
          compare_quant_results_for_wrong_decisions <- base::as.data.frame(matrix(ncol=2,nrow=length(selections_all)))
          colnames(compare_quant_results_for_wrong_decisions) <- c("correct","wrong")
          rownames(compare_quant_results_for_wrong_decisions) <- features_select_FDR[[s]]$Feature_name[selections_all]
          compare_quant_results_for_wrong_decisions$correct <- as.numeric(compare_quant_results_for_wrong_decisions$correct)
          compare_quant_results_for_wrong_decisions$wrong <- as.numeric(compare_quant_results_for_wrong_decisions$wrong)
          for(i in selections_all)
          {
            ind <- i
            data.table::set(compare_quant_results_for_wrong_decisions,as.integer(ind),as.integer(1:2),list(base::log2(10^feature_with_background_intensity_select_FDR[selection_feature[i],s]),
                                                                                                           base::log2(10^feature_with_background_intensity[selection_feature[i],s])))
          }
          ##compare intensities of unbiased peak selection and correct peak
          #determine distribution of quantification deviations between correct and unbiased selected peak
          if(plot==T)
          {
            # plot(density(compare_quant_results_for_wrong_decisions$wrong-compare_quant_results_for_wrong_decisions$correct,na.rm=T),xlab="Intensity-masked / Intensity-true, log2",main=base::paste(samples[s],"- Deviation in selected peak quantification"))
            # graphics::abline(v=0)
            # graphics::abline(v=-1,lty=2,col="red")
            # graphics::abline(v=1,lty=2,col="red")

            graphics::plot(compare_quant_results_for_wrong_decisions$correct,compare_quant_results_for_wrong_decisions$wrong,xlab="True peak quantification, log2",ylab="Masked peak quantification, log2",main=base::paste(samples[s],"- Error in quantification"))
            graphics::abline(a=0,b=1)
            graphics::abline(a=1,b=1,lty=2,col="red")
            graphics::abline(a=-1,b=1,lty=2,col="red")
          }

          #how many show deviation > 2 fold
          wrong_peak_intensity_outlier <- abs(compare_quant_results_for_wrong_decisions$wrong-compare_quant_results_for_wrong_decisions$correct) > 1

          #finally visualize results
          freq <- plyr::count(peak_decision_same_without_peak_known)
          if(length(which(freq$x == F)) == 0)
          {
            freq <- rbind(freq,base::data.frame(x=F,freq=0))
          }
          if(length(which(freq$x == T)) == 0)
          {
            freq <- rbind(freq,base::data.frame(x=T,freq=0))
          }
          freq$rel <- freq$freq/sum(freq$freq)*100

          plot_data <- c(freq$rel[which(freq$x==F)],
                         length(which(wrong_peak_intensity_outlier == T))/sum(freq$freq)*100)
          names(plot_data) <- c("Total",">2-fold intensity difference")

          results_peak_selection_FDR <- list(peak_decision_same_without_peak_known=peak_decision_same_without_peak_known,
                                             compare_quant_results_for_wrong_decisions=compare_quant_results_for_wrong_decisions,
                                             wrong_peak_intensity_outlier=wrong_peak_intensity_outlier,
                                             plot_data=plot_data)

        }else
        {
          freq <- plyr::count(peak_decision_same_without_peak_known)
          if(length(which(freq$x == F)) == 0)
          {
            freq <- rbind(freq,base::data.frame(x=F,freq=0))
          }
          if(length(which(freq$x == T)) == 0)
          {
            freq <- rbind(freq,base::data.frame(x=T,freq=0))
          }
          freq$rel <- freq$freq/sum(freq$freq)*100

          plot_data <- c(freq$rel[which(freq$x==F)],
                         0/sum(freq$freq)*100)
          names(plot_data) <- c("Total",">2-fold intensity difference")


          results_peak_selection_FDR <- list(peak_decision_same_without_peak_known=peak_decision_same_without_peak_known,
                                             compare_quant_results_for_wrong_decisions=NULL,
                                             wrong_peak_intensity_outlier=NULL,
                                             plot_data=plot_data)
        }
        results_peak_selection_FDR_all[[s]] <- results_peak_selection_FDR
      }else
      {
        results_peak_selection_FDR_all[[s]] <- NULL
      }


      updatecounter <- updatecounter + 1
      if(updatecounter >= 1)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(s/max))*(1-(s/max))
        td <- lubridate::seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        tcltk::setTkProgressBar(pb, s, label=base::paste( round(s/max*100, 0)," % done (",s,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    ##plot results over all samples
    Total_FDR <- NULL
    Large_Intensity_delta_FDR <- NULL
    for(s in 1:length(samples))
    {
      if(length(results_peak_selection_FDR_all[[s]]$peak_decision_same_without_peak_known) >= 100) ###require at least 100 identified peptides in a sample
      {
        Total_FDR <- append(Total_FDR,results_peak_selection_FDR_all[[s]]$plot_data[1])
        Large_Intensity_delta_FDR <- append(Large_Intensity_delta_FDR,results_peak_selection_FDR_all[[s]]$plot_data[2])
      }else
      {
        print(base::paste(samples[s],": Not enough known peaks for peak selection FDR estimation",sep=""))
        Total_FDR <- append(Total_FDR,NA)
        Large_Intensity_delta_FDR <- append(Large_Intensity_delta_FDR,NA)
      }

    }
    names(Total_FDR) <- samples
    names(Large_Intensity_delta_FDR) <- samples

    if(plot == T)
    {
      # ylim <- c(0,ifelse(max(Total_FDR,na.rm=T)<5,5,max(Total_FDR,na.rm=T)))
      # p <- Barplots(Total_FDR,AvgLine = T,digits_average = 1,Name = samples,xlab = "",ylab="FDR [%]",main = "Peak selection FDR - Total",shownumbers = F,ylim=ylim)
      # graphics::abline(h=5,lty=2,col="red")

      ylim <- c(0,ifelse(max(Large_Intensity_delta_FDR,na.rm=T)<5,5,max(Large_Intensity_delta_FDR,na.rm=T)))
      p <- Barplots(Large_Intensity_delta_FDR,AvgLine = T,digits_average = 1,Name = samples,xlab = "",ylab="FDR [%]",main = base::paste("Peak selection FDR - > 2-fold intensity difference\nBased on n=",num_features," random draws",sep=""),shownumbers = F,ylim=ylim)
      graphics::abline(h=5,lty=2,col="red")

    }

  }else
  {
    print("Warning: The true peak for too few features is known. Skipping peak selection FDR estimation.")
    results_peak_selection_FDR_all <- NA
    Total_FDR <- NA
    Large_Intensity_delta_FDR <- NA
  }


  return(list(results_peak_selection_FDR_all=results_peak_selection_FDR_all,
              Total_FDR=Total_FDR,
              Large_Intensity_delta_FDR=Large_Intensity_delta_FDR))

}


