#' Function to determine at which value a density maximum is reached
#' @param data Numeric vector
#' @details Uses kernel density estimation function from R-package stats removing missing values
#' @return Numeric indicating at which value a density maximum is reached
#' @examples
#' data <- c(1,1,1,2,2,4,4,5,10)
#' maxDensity(data)
maxDensity <- function(data)
{
  dens <- stats::density(data,na.rm=T)
  return(dens$x[which(dens$y == max(dens$y))])
}

#' Convert thermo raw files to mzXML with centroided ms1 scans using the ProteoWizard tool msConvert
#' @param path_to_raw Path to folder containing raw files which should be converted
#' @details Requires installation of ProteoWizard (http://proteowizard.sourceforge.net/download.html). Pay attention to installation requirements.
#' @return Resulting mzXML files are stored in a sub-directory within specified raw file folder
run_msconvert_raw_mzXML <- function(path_to_raw=NULL)
{
  suppressWarnings(suppressMessages(library(ff,quietly = T)))

  if(is.null(path_to_raw))path_to_raw <- choose.dir(caption = "Select folder containing raw files")

  ###get home directory
  home_folder <- Sys.getenv("HOME")
  home_folder <- gsub("Documents","AppData/Local/Apps/",home_folder)

  ###find MSConvert folder
  folders <- list.dirs(path = home_folder, full.names = TRUE, recursive = F)
  folders <- folders[which(grepl("ProteoWizard",folders))]
  folders <- folders[length(folders)]

  if(file.exists(paste(folders,"\\msconvert.exe",sep="")))
  {
    path_to_msconvert <- paste(folders,"\\msconvert.exe",sep="")
  }else
  {
    print("Select msconvert.exe. Can be usually found in C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version")
    path_to_msconvert <- file.choose()
  }

  raw_files <- list.files(path_to_raw)
  raw_files <- raw_files[which(grepl("\\.raw",raw_files))]

  ###check which raw files still have to be converted
  mzXMLs_available <- list.files(paste(path_to_raw,"\\mzXML",sep=""))

  files_to_be_converted <- raw_files[which(gsub("\\.raw","",raw_files) %not in% gsub("\\.mzXML","",mzXMLs_available))]

  if(length(files_to_be_converted)>0)
  {
    ###get user folder
    win_user_folder <- path.expand('~')

    ###create temporary folder
    dir.create(paste(win_user_folder,"\\temp_msconvert",sep=""),showWarnings = F)
    dir.create(paste(win_user_folder,"\\temp_msconvert\\mzXML",sep=""),showWarnings = F)

    ###create temporary config.txt and files.txt
    temp_path <- paste(win_user_folder,"\\temp_msconvert",sep="")
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
    writeLines(paste(path_to_raw,"\\",files_to_be_converted,sep=""), fileConn)
    close(fileConn)

    ###prepare arguments for msconvert
    arg <- paste("-f ",temp_path,"\\files.txt",
                 " -o ",temp_path,"\\mzXML",
                 " -c ",temp_path,"\\config.txt",sep="")

    ###run msconvert
    system2(path_to_msconvert, args = arg)

    ###move mzXMLs to original raw folder
    from <- temp_path
    to   <- path_to_raw
    path1 <- paste0(from,"\\mzXML")
    path2 <- paste0(to,"\\mzXML")

    dir.create(path2,showWarnings = F)

    for(f in gsub("\\.raw",".mzXML",files_to_be_converted))
    {
      file.move(paste(path1,"\\",f,sep=""),path2)
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
mzxml_to_list <- function(path_to_mzXML,n_cores=2)
{
  suppressWarnings(suppressMessages(library(doParallel,quietly = T)))

  convert <- function(mzXMLfile,path_to_mzXML)
  {
    suppressWarnings(suppressMessages(library("readMzXmlData",quietly = T)))

    data <- paste(path_to_mzXML,"\\",mzXMLfile,sep="")
    pb <- winProgressBar(title = "Read mzXML",label=paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    ms <- readMzXmlFile(data)
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
    suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    dat <- as.data.table(matrix(ncol=3,nrow=rowcount))
    colnames(dat) <- c("m.z","RT","Intensity")
    sample <- mzXMLfile
    sample <- substr(sample,1,regexpr(".mzXML",sample)-1)
    dat$m.z <- as.numeric(dat$m.z)
    dat$RT <- as.numeric(dat$RT)
    dat$Intensity <- as.numeric(dat$Intensity)
    ind <- 1
    max <- length(ms)
    pb <- winProgressBar(title = "Extract data",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)

    for(i in 1:length(ms))
    {
      if(ms[[i]]$metaData$msLevel == 1)
      {
        start <- ind
        stop <- ind + length(ms[[i]]$spectrum$mass) - 1
        set(x = dat,i = start:stop,j = (1L),value = ms[[i]]$spectrum$mass)
        set(x = dat,i = start:stop,j = (2L),value = ms[[i]]$metaData$retentionTime/60)
        set(x = dat,i = start:stop,j = (3L),value = ms[[i]]$spectrum$intensity)
        ind <- ind + length(ms[[i]]$spectrum$mass)
      }
      setWinProgressBar(pb, i, label = paste( round(i/max*100, 0),"% done (",i,"/",max,")",sep=""))
    }
    close(pb)
    save(dat,file=paste(path_to_mzXML,"\\all_ion_lists\\",sample,"_all_ions.RData",sep=""))
  }

  ##Step - Extract all ions per ms1 spectra

  mzXMLfiles <- list.files(path_to_mzXML)
  mzXMLfiles <- mzXMLfiles[which(grepl(".mzXML",mzXMLfiles))]

  dir.create(paste(path_to_mzXML,"\\all_ion_lists",sep=""),showWarnings = F)

  setwd(paste(path_to_mzXML,"\\all_ion_lists",sep=""))

  ###check if all ion list are already available and only generate those which are required
  available_all_ion_lists <- list.files()[which(grepl("\\.RData",list.files()))]

  missing_all_ion_lists <- mzXMLfiles[which(gsub("\\.mzXML","",mzXMLfiles) %not in% gsub("_all_ions\\.RData","",available_all_ion_lists))]
  if(length(missing_all_ion_lists)>0)
  {
    mzXMLfiles <- missing_all_ion_lists

    cl <- makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    registerDoParallel(cl)
    res <- foreach(i=mzXMLfiles) %dopar%
      {
        convert(i,path_to_mzXML)
      }
    stopCluster(cl)
  }

}


#' Perform alignment of pre-determined MS1-features by MaxQuant over proteomics samples
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
#' @details Performs the first steps of the IceR workflow: 1) Alignment window determination if not specified. 2) Alignment of MaxQuant features into IceR features. 3) Transfer of sequence information between MaxQuant features aligned into IceR features. 4) Extraction, modelling and prediction of RT- and m/z-correction factor per IceR feature and sample.
#' @return Outputs are stored in the sub-directory Temporary_files within specified output folder. MaxQuant allpeptides.txt and evidence.txt are converted to RData files. QC plots of estimated alignment windows as well as of random forest modesl and generalized additive models are stored in a QC_plots.pdf. Relevant QC data is stored in Feature_alignment_QC_data.RData. Aligned IceR features are stored in Features_aligned_merged.txt
align_features <- function(path_to_MaxQ_output,path_to_output,align_unknown=F,output_file_names_add="IceR_analysis",mz_window=NA,min_mz_window = 0.001,RT_window=NA,min_RT_window=1,min_num_ions_collapse=10,feature_mass_deviation_collapse=0.002,only_unmodified_peptides=F,sample_list=NULL,remove_contaminants=T)
{
  options(warn=-1)
  suppressWarnings(suppressMessages(library(data.table,quietly = T)))
  suppressWarnings(suppressMessages(library(stringr,quietly = T)))
  suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
  suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
  suppressWarnings(suppressMessages(library(ff,quietly = T)))
  use_mz_at_max_int_for_correction=F

  if(output_file_names_add != "")output_file_names_add <- paste("_",output_file_names_add,sep="")

  setwd(path_to_output)
  dir.create("Temporary_files")

  if(file.exists(paste("Temporary_files\\Features_aligned_merged",output_file_names_add,".txt",sep="")))
  {
    print("Alignment already done")
  }else
  {
    if(file.exists("Temporary_files\\allPeptides.RData"))
    {
      load("Temporary_files\\allPeptides.RData")
      load("Temporary_files\\evidence.RData")

      if(is.null(sample_list))
      {
        sample_list <- sort(unique(evidence$Raw.file))
      }

    }else
    {
      setwd(path_to_MaxQ_output)
      options(fftempdir = path_to_MaxQ_output)
      print("Read MaxQ results")
      temp <- read.csv(file = "allPeptides.txt",sep='\t',nrows = 2044,header=T)
      tempclasses = sapply(temp, class)
      tempclasses[29] = "numeric"
      tempclasses[which(tempclasses == "logical")] <- "factor"
      tempclasses[32] = "factor"
      tempclasses[33] = "factor"
      allpeptides_save <- read.csv.ffdf(file = "allPeptides.txt",sep='\t',VERBOSE = F,colClasses=tempclasses,next.rows = 100000)##read in data in chunks of 100000 rows
      allpeptides <- as.data.frame(allpeptides_save)
      allpeptides <- allpeptides[order(allpeptides$Mass),]

      #allpeptidessave <- allpeptides
      ###free some memory
      rm(allpeptides_save)
      gc()


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

      allpeptides$Retention.Length <- allpeptides$Retention.Length/60

      ###Some peptides are not available in allpeptides.txt thus we get this additional information from the evidence.txt file
      evidence <- read.csv(file = "evidence.txt",sep='\t',header=T)
      if(any(colnames(evidence) == "MS.MS.Scan.Number"))colnames(evidence)[which(colnames(evidence) == "MS.MS.Scan.Number")] <- "MSMS.Scan.Numbers"
      if(any(colnames(evidence) == "MS.MS.scan.number"))colnames(evidence)[which(colnames(evidence) == "MS.MS.scan.number")] <- "MSMS.Scan.Numbers"
      if(any(colnames(evidence) == "Retention.length"))colnames(evidence)[which(colnames(evidence) == "Retention.length")] <- "Retention.Length"
      print("Read MaxQ results finished")

      add_data <- as.data.frame(matrix(ncol=ncol(allpeptides),nrow=nrow(evidence),NA))
      colnames(add_data) <- colnames(allpeptides)
      match_col_names <- match(colnames(allpeptides),colnames(evidence))
      set(add_data,j=as.integer(which(!is.na(match_col_names))),value = evidence[,match_col_names[which(!is.na(match_col_names))]])

      ####now find rows (peptides + modification + sample) which are not yet present in allpeptides.txt
      temp1 <- paste(add_data$Sequence,add_data$Modifications,add_data$Raw.file,sep="_")
      temp2 <- paste(allpeptides$Sequence,allpeptides$Modifications,allpeptides$Raw.file,sep="_")
      add_data <- add_data[which(temp1 %not in% temp2),]
      allpeptides <- rbind(allpeptides,add_data)

      rm(temp1,temp2,add_data)
      gc()

      ###if no specific sample list is defined which should be used, use all raw files
      if(is.null(sample_list))
      {
        sample_list <- sort(unique(evidence$Raw.file))
      }

      ###clean allpeptides table to reduce required memory space
      if(any(colnames(allpeptides) == "Resolution"))
      {
        allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","Mass","Uncalibrated.m.z","Resolution","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity")]
      }else
      {
        allpeptides <- allpeptides[,c("Raw.file","Charge","m.z","Mass","Uncalibrated.m.z","Max.intensity.m.z.0","Retention.time","Retention.Length","MS.MS.IDs","Sequence","Modifications","Proteins","Score","MSMS.Scan.Numbers","Intensity")]
        allpeptides$Resolution <- 0
        print("MS resolution seems to be missing in MaxQ outputs !!!")
      }

      if(remove_contaminants == T)exclude <- which(grepl("CON",allpeptides$Proteins) | grepl("REV",allpeptides$Proteins))
      if(remove_contaminants == F)exclude <- which(grepl("REV",allpeptides$Proteins))
      if(length(exclude) > 0)allpeptides <- allpeptides[-exclude,] ###remove potential reverse and contaminant peptides


      ###add calibrated RT to all peptides
      temp_evidence <- evidence[,c("Sequence","Raw.file","Charge","Modifications","Calibrated.retention.time")]
      temp_evidence <- aggregate(temp_evidence$Calibrated.retention.time,list(Sequence=evidence$Sequence,Raw.file=evidence$Raw.file,Charge=evidence$Charge,Modifications=evidence$Modifications),FUN=mean,na.rm=T)
      colnames(temp_evidence)[5] <- "Calibrated.retention.time"
      temp <- left_join(allpeptides,temp_evidence,by=c("Sequence"="Sequence","Raw.file"="Raw.file","Charge"="Charge","Modifications"="Modifications"))

      missing_RT_calibration_indices <- which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")
      if(length(missing_RT_calibration_indices)>0)
      {
        for(ind in missing_RT_calibration_indices)
        {
          sub <- temp[which(temp$Sequence == temp$Sequence[ind] & temp$Modifications == temp$Modifications[ind]),]
          temp$Calibrated.retention.time[ind] <- median(sub$Calibrated.retention.time,na.rm=T)
        }
      }

      ###if any peptide feature still doesn?t have a valid calibrated RT --> just use observed RT
      temp$Calibrated.retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")] <- temp$Retention.time[which(is.na(temp$Calibrated.retention.time) & temp$Sequence != " " & temp$Sequence != "")]
      allpeptides <- temp

      ###calculate deviation of mz between mz at max intensity and true mz
      ##allpeptides
      add <- as.data.frame(matrix(nrow=nrow(allpeptides),ncol=3))
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
      pb <- winProgressBar(title = "Determine deviations of true m/z from observed m/z",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          set(add,as.integer(i),as.integer(1:3),value=as.list(c(closest_iso,mz_max_int,delta_mz_to_mz_at_max_int)))
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)

      allpeptides <- cbind(allpeptides,add)

      ###evidence table

      add <- as.data.frame(matrix(nrow=nrow(evidence),ncol=3))
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
      pb <- winProgressBar(title = "Determine deviations of true m/z from observed m/z",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          set(add,as.integer(i),as.integer(1:3),value=as.list(c(closest_iso,mz_max_int,delta_mz_to_mz_at_max_int)))
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)

      evidence <- cbind(evidence,add)
      ####
      ###finally save prepared data
      setwd(path_to_output)

      save(allpeptides,file = "Temporary_files\\allPeptides.RData")
      save(evidence,file = "Temporary_files\\evidence.RData")
    }

    QC_data <- list() ##here relevant qc data is stored and finally saved as RData which can be used for re-generating plots

    pdf(paste("Temporary_files\\QC_plots",output_file_names_add,".pdf",sep=""))

    RT_calibration <- T
    mz_calibration <- T

    QC_data[["MaxQ_calibrations"]] <- evidence[,c("Retention.time.calibration","Uncalibrated...Calibrated.m.z..Da.","delta_mz_to_mz_at_max_int")]

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


    windows <- as.data.frame(matrix(ncol=11,nrow=nrow(identified_ions)))
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

    ###remove outlier samples where RT of a peptide is very different from all other samples
    outlier_RT_deviation <- (max(allpeptides$Retention.time,na.rm=T)/100)*5 ###deviation should not be larger than 5 % of the total chromatographic retention length

    max <- length(unique_peptides)
    pb <- winProgressBar(title = "Determine matching parameter windows",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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

          median_RT <- median(sub1$Retention.time,na.rm=T)
          remove <- which(sub1$Retention.time > median_RT + outlier_RT_deviation)
          if(length(remove)>0)sub1 <- sub1[-remove,]
          ###determine mean m/z at peak maximum and mean RT
          set(windows,as.integer(count_features),as.integer(c(3:11)),value=as.list(c(c,
                                                                                     mean(sub1$m.z,na.rm=T),
                                                                                     sd(sub1$m.z,na.rm=T),
                                                                                     mean(sub1$Retention.time,na.rm=T),
                                                                                     sd(sub1$Retention.time,na.rm=T),
                                                                                     nrow(sub1),mean(sub1$Calibrated.retention.time,na.rm=T),
                                                                                     sd(sub1$Calibrated.retention.time,na.rm=T),
                                                                                     sd(sub1$isotope_corrected_mz_at_max_int,na.rm=T))))
          set(windows,as.integer(count_features),as.integer(1:2),value=as.list(c(unique_peptides[i],m)))
        }

      }

      updatecounter <- updatecounter + 1
      if(updatecounter >= 10)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    windows <- windows[1:count_features,]

    boxplot(windows$sd_m.z,main="Standard deviation of peptide features m/z",outline=F)
    boxplot(windows$sd_RT,main="Standard deviation of peptide features RT",outline=F)

    if(is.na(mz_window))
    {
      ###based on boxplots
      if(!is.na(min_mz_window))
      {
        if(boxplot.stats(windows$sd_m.z)$`stats`[5]>min_mz_window)
        {
          borders_m.z <- c(-boxplot.stats(windows$sd_m.z)$`stats`[5],boxplot.stats(windows$sd_m.z)$`stats`[5]) ###upper whisker
        }else
        {
          borders_m.z <- c(-min_mz_window,min_mz_window)
        }
      }else
      {
        borders_m.z <- c(-boxplot.stats(windows$sd_m.z)$`stats`[5],boxplot.stats(windows$sd_m.z)$`stats`[5]) ###upper whisker
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
        if(boxplot.stats(windows$sd_RT)$`stats`[5]>min_RT_window)
        {
          borders_RT <- c(-boxplot.stats(windows$sd_RT)$`stats`[5],boxplot.stats(windows$sd_RT)$`stats`[5]) ##75% quantil based on boxplots
        }else
        {
          borders_RT <- c(-min_RT_window,min_RT_window)
        }
      }else
      {
        borders_RT <- c(-boxplot.stats(windows$sd_RT)$`stats`[5],boxplot.stats(windows$sd_RT)$`stats`[5]) ##75% quantil based on boxplots
      }

    }else
    {
      ###use defined parameter
      borders_RT <- c(-RT_window,RT_window)
    }

    borders_RT_use_save = borders_RT

    QC_data[["Feature_alignment_windows"]] <- list(RT_window=borders_RT,mz_window=borders_m.z)

    ###determine for how many windows we expect an overlap of ions for different peptide sequences based on the chosen parameters
    windows <- windows[order(windows$mean_m.z),]
    max_i <- nrow(windows)
    windows$overlap <- 0

    max <- nrow(windows)
    pb <- winProgressBar(title = "Determine overlapping peptide features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(windows))
    {
      start <- ifelse((i - 5) < 1,1,(i - 5))
      end <- ifelse((i + 5) > max_i,max_i,(i + 5))
      sub <- windows[start:end,]
      sub <- sub[which(sub$mean_m.z >= windows$mean_m.z[i] + borders_m.z[1] & sub$mean_m.z <= windows$mean_m.z[i] + borders_m.z[2] & sub$mean_RT >= windows$mean_RT[i] + borders_RT[1] & sub$mean_RT <= windows$mean_RT[i] + borders_RT[2] & sub$Charge == windows$Charge[i]),]

      set(windows,as.integer(i),12L,value=(nrow(sub)-1))

      updatecounter <- updatecounter + 1
      if(updatecounter >= 100)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
    allpeptides$Proteins <- gsub(";","|",allpeptides$Proteins)
    allpeptides$Raw.file <- as.character(allpeptides$Raw.file)
    ###presubset all ions based on charge, here use only rows of allions which are coming from samples which should be also used
    ####also create a library of ions per charge state per mz_window of 0.5 Da
    print("Prepare for peptide feature matching - Indexing all ions")
    rownames(allpeptides) <- c(1:nrow(allpeptides))
    allpeptides_frag <- list()
    allpeptides_frag_indices_per_mz_window <- list()

    for(i in 1:max(allpeptides$Charge))
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
    }
    rm(indices)
    gc()

    temp_data <- as.data.frame(matrix(ncol=1,nrow=nrow(allpeptides))) ###here already used ions will be marked
    temp_data[,1] <- as.numeric(temp_data[,1])

    temp <- subset(allpeptides,Sequence != " " & Sequence != "")

    features <- as.data.frame(matrix(ncol=21,nrow=length(unique(paste(temp$Sequence,temp$Charge,temp$Modifications)))))
    colnames(features) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")

    features <- as.data.table(features)
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
    features$mean_Scores <- as.character(features$ion_indices)
    features$num_matches <- as.character(features$num_matches)
    features$Modifications <- as.character(features$Modifications)
    features$RT_length <- as.numeric(features$RT_length)
    features$Observed_RT <- as.character(features$Observed_RT)
    features$Observed_mz <- as.character(features$Observed_mz)
    features$Observed_score <- as.character(features$Observed_score)
    ###get ordering of allpeptides based on sequence
    order_by_Sequence <- order(allpeptides$Sequence,decreasing = T)
    unique_peptides <- sort(unique_peptides,decreasing = T)
    number_of_rows_per_peptide_sequence <- plyr::count(allpeptides$Sequence[which(allpeptides$Sequence != " " & allpeptides$Sequence != "")])
    number_of_rows_per_peptide_sequence <- number_of_rows_per_peptide_sequence[order(number_of_rows_per_peptide_sequence$x,decreasing = T),]

    max <- length(unique_peptides)
    pb <- winProgressBar(title = "Matching known features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    last_index <- 1 ###required for subsetting allpeptides per peptide sequence
    count_features <- 0

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
            median_m.z <- median(sub$m.z,na.rm=T)##
            median_RT <- median(sub$Retention.time,na.rm=T)
            selection <- which(abs(sub$Retention.time-median_RT) < borders_RT[2] & abs(sub$m.z-median_m.z) < borders_m.z[2])
            sub <- sub[selection,]

            if(nrow(sub) > 0)
            {
              ###check if peptide was sequenced several times. If this is the case, use feature closest to median
              if(any(duplicated(sub$Raw.file)))
              {
                sub <- sub[!duplicated(sub$Raw.file),]
              }

              median_m.z <- median(sub$m.z,na.rm=T)##
              median_RT <- median(sub$Retention.time,na.rm=T)

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
              box_stats_rt <- boxplot.stats(sub$Retention.time)$stats
              box_stats_mz <- boxplot.stats(sub$m.z)$stats
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
                name <- paste("Feature",count_features,sep="_")

                sequences_relevant <- which(nchar(as.character(sub$Sequence)) > 1)
                if(length(sequences_relevant) > 0)
                {
                  sequence <- paste(unique(c(sub_peptide$Sequence[1],sub[sequences_relevant,"Sequence"])),collapse=";")
                }else
                {
                  sequence <- sub_peptide$Sequence[1]
                }

                proteins_relevant <- which(nchar(as.character(sub$Proteins)) > 1)
                if(length(proteins_relevant) > 0)
                {
                  protein <- paste(unique(c(sub_peptide$Proteins[1],sub[proteins_relevant,"Proteins"])),collapse=";")
                }else
                {
                  protein <- sub_peptide$Proteins[1]
                }

                msmsscansrelevant_relevant <- which(nchar(as.character(sub$MSMS.Scan.Numbers)) > 1)
                if(length(msmsscansrelevant_relevant) > 0)
                {
                  msmsscan <- paste(unique(paste(sub[msmsscansrelevant_relevant,"Raw.file"],"_msms_",sub[msmsscansrelevant_relevant,"MSMS.Scan.Numbers"],sep="")),collapse=";")
                }else
                {
                  msmsscan <- ""
                }

                median_mass <- median(sub$Mass,na.rm=T)
                Charge <- c

                scores <- NULL
                num_matches <- NULL
                for(s in unique(unlist(str_split(sequence,";"))))
                {
                  scores <- append(scores,mean(sub$Score[which(sub$Sequence == s)],na.rm=T))
                  num_matches <- append(num_matches,nrow(sub[which(sub$Sequence == s),]))
                }

                Modifications <- ifelse(m != "Unmodified",m,"")

                RT_length <- max(sub$Retention.Length,na.rm = T)

                temp_RT <- aggregate(sub$Retention.time,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                Observed_RT <- paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

                temp_mz <- aggregate(sub$isotope_corrected_mz_at_max_int,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                Observed_mz <- paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

                temp_score <- aggregate(sub$Score,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)

                Observed_score <- paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed scores according to sample_list


                set(x = features,i = as.integer(count_features),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
                                                                                                                                  sequence,
                                                                                                                                  protein,
                                                                                                                                  msmsscan,
                                                                                                                                  paste(rownames(sub),collapse = ","),
                                                                                                                                  paste(scores,collapse=";"),
                                                                                                                                  paste(num_matches,collapse=";"),
                                                                                                                                  Modifications,
                                                                                                                                  Observed_RT,
                                                                                                                                  Observed_mz,
                                                                                                                                  Observed_score
                )))
                set(x = features,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
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

                ###now mark matched features
                set(temp_data,i = as.integer(rownames(sub)),j=1L,value=1)
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
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    features <- features[1:count_features,]
    #write.table(features,paste("Temporary_files\\Features_aligned_peptides",output_file_names_add,".txt",sep=""),row.names = F)

    borders_RT <- borders_RT_use_save

    save(temp_data,features,allpeptides,borders_m.z,borders_RT,windows,file=paste("Temporary_files\\features_aligned_step_1",output_file_names_add,".RData",sep=""))

    #####Perform alignment of undefined features

    if(align_unknown == T)
    {
      ###2.) step - align remaining features

      #load(paste("features_aligned_step_1",output_file_names_add,".RData",sep=""))

      count_features <- nrow(features)

      ###start by trying to find unique features over all samples defined by the pair Mass and RT
      border_factor <- 4

      ####use this matrix for faster subsetting
      temp_data_m <- as.matrix(subset(allpeptides,select=c("m.z","Retention.time","Charge")))
      temp_data_m <- cbind(temp_data_m,temp_data)
      colnames(temp_data_m)[4] <- "Matched"
      temp_data_m[which(is.na(temp_data_m[,4])),4] <- 0
      #temp_data_m <- cbind(temp_data_m,as.matrix(subset(allpeptides,select=c("Retention.Length"))))

      features_unknown <- as.data.frame(matrix(ncol=21,nrow=(nrow(temp_data)/min_num_ions_collapse)))
      colnames(features_unknown) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")


      features_unknown <- as.data.table(features_unknown)
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
      features_unknown$mean_Scores <- as.character(features_unknown$ion_indices)
      features_unknown$num_matches <- as.character(features_unknown$num_matches)
      features_unknown$Modifications <- as.character(features_unknown$Modifications)
      features_unknown$RT_length <- as.numeric(features_unknown$RT_length)
      features_unknown$Observed_RT <- as.character(features_unknown$Observed_RT)
      features_unknown$Observed_mz <- as.character(features_unknown$Observed_mz)
      features_unknown$Observed_score <- as.character(features_unknown$Observed_score)

      max <- nrow(temp_data_m)
      pb <- winProgressBar(title = "Aligning unknown features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      start <- 1
      end <- max-1

      count_features_unknown <- 0

      for(i in start:end) ###go through all features and try to match with others in the data set
      {
        if(temp_data_m[i,4] == 0)###check if already matched
        {
          ind <- i
          while(T)
          {
            ind = ind + 1000
            if(ind > nrow(temp_data_m))
            {
              ind <- nrow(temp_data_m)
              break
            }
            if(temp_data_m[ind,1] > (temp_data_m[i,1]+(border_factor*borders_m.z[2])))
            {
              break
            }
          }

          sub <- temp_data_m[i:ind,]
          sub <- sub[which(sub[,4] == 0),]

          ###determine which other features might fit to the current unmatched feature
          ###here use a wide mass and RT range to get a first preselection for Mass and RT distribution estimation                                  )

          matches_preselection <- which(sub[,1] >= (temp_data_m[i,1]+(border_factor*borders_m.z[1])) & sub[,1] <= (temp_data_m[i,1]+(border_factor*borders_m.z[2])) & sub[,2] >= (temp_data_m[i,2]+(border_factor*borders_RT[1])) & sub[,2] <= (temp_data_m[i,2]+(border_factor*borders_RT[2])) & as.numeric(sub[,3]) == temp_data_m[i,3])

          #Sys.time()-start
          if(length(matches_preselection) >= min_num_ions_collapse) ###at least n other features
          {
            ###now get estimate of true Mass and RT by using the median
            median_m.z <- median(sub[matches_preselection,1],na.rm=T)
            median_RT <- median(sub[matches_preselection,2],na.rm=T)

            ###determine which features will be finally matched with the normal m/z and RT range
            matches_temp <- which(as.numeric(sub[,1]) >= (median_m.z+borders_m.z[1]) & as.numeric(sub[,1]) <= (median_m.z+borders_m.z[2]) & as.numeric(sub[,2]) >= (median_RT+borders_RT[1]) & as.numeric(sub[,2]) <= (median_RT+borders_RT[2]) & as.numeric(sub[,3]) == temp_data_m[i,3])

            matches <- as.numeric(rownames(sub)[matches_temp]) ###convert back to original data

            if(length(matches) >= min_num_ions_collapse)
            {
              count_features <- count_features + 1
              count_features_unknown <- count_features_unknown + 1
              ###new feature - give a new name
              name <- paste("Feature",count_features,sep="_")

              sequences_relevant <- which(nchar(as.character(allpeptides[matches,"Sequence"])) > 1)
              if(length(sequences_relevant) > 0)
              {
                sequence <- paste(unique(allpeptides[matches[sequences_relevant],"Sequence"]),collapse=";")
              }else
              {
                sequence <- ""
              }

              proteins_relevant <- which(nchar(as.character(allpeptides[matches,"Proteins"])) > 1)
              if(length(proteins_relevant) > 0)
              {
                protein <- paste(unique(allpeptides[matches[proteins_relevant],"Proteins"]),collapse=";")
              }else
              {
                protein <- ""
              }

              msmsscansrelevant_relevant <- which(nchar(as.character(allpeptides[matches,"MSMS.Scan.Numbers"])) > 1)
              if(length(msmsscansrelevant_relevant) > 0)
              {
                msmsscan <- paste(unique(paste(allpeptides[matches[msmsscansrelevant_relevant],"Raw.file"],"_msms_",allpeptides[matches[msmsscansrelevant_relevant],"MSMS.Scan.Numbers"],sep="")),collapse=";")
              }else
              {
                msmsscan <- ""
              }

              median_mass <- median(allpeptides$Mass[matches],na.rm=T)
              Charge <- temp_data_m[i,3]


              if(sequence != "")
              {
                scores <- NULL
                num_matches <- NULL
                for(s in unique(unlist(str_split(sequence,";"))))
                {
                  scores <- append(scores,mean(allpeptides[matches,"Score"][which(allpeptides[matches,"Sequence"] == s)],na.rm=T))
                  num_matches <- append(num_matches,length(which(allpeptides[matches,"Sequence"] == s)))
                }

                Modifications_relevant <- which(nchar(as.character(allpeptides[matches,"Modifications"])) > 1)
                if(length(Modifications_relevant) > 0)
                {
                  Modifications <- paste(unique(allpeptides[matches[Modifications_relevant],"Modifications"]),collapse=";")
                  Modifications <- gsub("Unmodified","",Modifications)

                }else
                {
                  Modifications <- ""
                }
              }else
              {
                scores <- NA
                num_matches <- NA
                Modifications <- ""
              }

              RT_length <- median(allpeptides$Retention.Length[matches],na.rm=T)

              temp_RT <- aggregate(allpeptides$Retention.time[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
              Observed_RT <- paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

              temp_mz <- aggregate(allpeptides$isotope_corrected_mz_at_max_int[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
              Observed_mz <- paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

              temp_score <- aggregate(allpeptides$Score[matches],by=list(Raw.file=allpeptides$Raw.file[matches]),FUN=mean,na.rm=T)
              Observed_score <- paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed scores according to sample_list


              set(x = features_unknown,i = as.integer(count_features_unknown),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
                                                                                                                                                sequence,
                                                                                                                                                protein,
                                                                                                                                                msmsscan,
                                                                                                                                                paste(matches,collapse = ","),
                                                                                                                                                paste(scores,collapse=";"),
                                                                                                                                                paste(num_matches,collapse=";"),
                                                                                                                                                Modifications,
                                                                                                                                                Observed_RT,
                                                                                                                                                Observed_mz,
                                                                                                                                                Observed_score
              )))
              set(x = features_unknown,i = as.integer(count_features_unknown),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
                as.numeric(median_m.z),
                as.numeric(median_RT),
                as.numeric(median_m.z+borders_m.z[1]),
                as.numeric(median_m.z+borders_m.z[2]),
                as.numeric(median_RT+borders_RT[1]),
                as.numeric(median_RT+borders_RT[2]),
                as.numeric(length(matches)),
                as.numeric(median_mass),
                as.numeric(Charge),
                RT_length
              )))


              ###now mark matched features
              # starttime <- Sys.time()
              # temp_data_m[matches,4] <- 1
              # print(Sys.time()-starttime)

              set(temp_data_m,i=as.integer(matches),j=4L,value=1)

            }


          }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 5000)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }

      }
      close(pb)

      features_unknown <- features_unknown[1:count_features_unknown,]

      #write.table(features_unknown,paste("Temporary_files\\Features_aligned_unknown",output_file_names_add,".txt",sep=""),row.names = F)

      ####Combine known and unknown features
      features <- rbind(features,features_unknown)
      rownames(features) <- c()

      ###number of matched to number of unmatched ions
      value <- length(which(temp_data_m$Matched == 1))/nrow(temp_data_m) * 100
      barplot(value,main="Relative fraction of features which are matched",ylab="Fraction [%]")

      QC_data[["Numbers_matched_features"]] <- value

      ###Store all results

      save(temp_data_m,features,count_features_unknown,file=paste("features_aligned_step_2",output_file_names_add,".RData",sep=""))
    }

    #write.table(features,paste("Temporary_files\\Features_aligned",output_file_names_add,".txt",sep=""),row.names = F)

    ####Merge features which were currently regarded as separated although same RT and m/z
    #features <- read.table("Features_alligned.txt",header = T)

    #####prepare columns to indicate number of charges and with which other features the respective feature was merged
    features$num_diff_charges <- NA
    features$num_diff_charges <- as.numeric(features$num_diff_charges)
    features$merged_with <- ""
    ncolumns <- ncol(features)

    ###collapse features with same m/z and RT parameters and count number of different charge states for same feature
    borders_mass <- c(-feature_mass_deviation_collapse,feature_mass_deviation_collapse)

    max <- nrow(features)
    pb <- winProgressBar(title = "Merge features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(features))
    {
      if(is.na(features[i,"num_diff_charges"]))
      {
        selection <- which(features$Mass >= features$Mass[i]+borders_mass[1] & features$Mass <= features$Mass[i]+borders_mass[2] & features$RT >= features$RT[i]+borders_RT[1] & features$RT <= features$RT[i]+borders_RT[2])
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
              ion_indices <- paste(unlist(features[selection[which(features[selection,]$Charge == c)],]$ion_indices),collapse=",")
              count_ion_indices <- sum(features[selection[which(features[selection,]$Charge == c)],]$count_ion_indices,na.rm = T)
              mass <- mean(unlist(features[selection[which(features[selection,]$Charge == c)],]$Mass),na.rm=T)
              m.z <- mean(c(m.z_range_min,m.z_range_max))
              RT <- mean(c(RT_range_min,RT_range_max))

              Sequence <- unique(unlist(str_split(features[selection[which(features[selection,]$Charge == c)],]$Sequence,";")))

              scores <- NULL
              num_matches <- NULL
              Modifications <- NULL
              for(s in unique(Sequence))
              {
                scores_temp <- as.numeric(unlist(str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$mean_Scores,";")))
                scores <- append(scores,mean(scores_temp,na.rm=T))

                num_matches_temp <- as.numeric(unlist(str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$num_matches,";")))
                num_matches <- append(num_matches,sum(num_matches_temp,na.rm=T))

                Modifications_temp <- unlist(str_split(features[selection[which(features[selection,"Charge"] == c & grepl(s,features[selection,"Sequence"]))],]$Modifications,";"))
                Modifications <- append(Modifications,paste(Modifications_temp,collapse=","))
              }
              scores <- paste(scores,collapse=",")
              num_matches <- paste(num_matches,collapse=",")
              Modifications <- paste(Modifications,collapse=",")

              if(any(Sequence == ""))Sequence <- Sequence[-which(Sequence == "")]
              Sequence <- paste(Sequence,collapse=";")

              Protein <- unique(unlist(str_split(features[selection[which(features[selection,]$Charge == c)],]$Protein,";")))
              if(any(Protein == ""))Protein <- Protein[-which(Protein == "")]
              Protein <- paste(Protein,collapse=";")

              MSMS.Scan.Numbers <- unique(unlist(str_split(features[selection[which(features[selection,]$Charge == c)],]$MSMS.Scan.Numbers,";")))
              if(any(MSMS.Scan.Numbers == ""))MSMS.Scan.Numbers <- MSMS.Scan.Numbers[-which(MSMS.Scan.Numbers == "")]
              MSMS.Scan.Numbers <- paste(MSMS.Scan.Numbers,collapse=";")

              Observed_RTs <- str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_RT,pattern = ";",simplify = T)
              Observed_RTs <- apply(Observed_RTs, 2, as.numeric)
              Observed_RTs <- colMeans(Observed_RTs,na.rm=T)
              Observed_RTs[is.na(Observed_RTs)] <- NA
              Observed_RTs <- paste(Observed_RTs,collapse=";")

              Observed_mz <- str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_mz,pattern = ";",simplify = T)
              Observed_mz <- apply(Observed_mz, 2, as.numeric)
              Observed_mz <- colMeans(Observed_mz,na.rm=T)
              Observed_mz[is.na(Observed_mz)] <- NA
              Observed_mz <- paste(Observed_mz,collapse=";")

              Observed_score <- str_split(features[selection[which(features[selection,]$Charge == c)],]$Observed_score,pattern = ";",simplify = T)
              Observed_score <- apply(Observed_score, 2, as.numeric)
              Observed_score <- colMeans(Observed_score,na.rm=T)
              Observed_score[is.na(Observed_score)] <- NA
              Observed_score <- paste(Observed_score,collapse=";")

              set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(2,3,7,8,9,10,12,13,14)),value=as.list(c(
                m.z,
                RT,
                m.z_range_min,
                m.z_range_max,
                RT_range_min,
                RT_range_max,
                count_ion_indices,
                mass,
                c)))

              set(x=features,i=selection[which(features[selection,]$Charge == c)[1]],j=as.integer(c(4,5,6,11,15,16,17,19,20,21)),value=as.list(c(
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

              set(x=features,i=selection[which(features[selection,]$Charge == c)[-1]],j=as.integer(ncolumns),value=features$Feature_name[selection[which(features[selection,]$Charge == c)[1]]])
            }
          }
          ###add number of different charges count
          set(x=features,i=selection,j=as.integer(ncolumns-1),value=length(unique(features[selection,]$Charge)))
        }
      }
      updatecounter <- updatecounter + 1
      if(updatecounter >= 100)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    features$num_diff_charges[is.na(features$num_diff_charges)] <- 1

    features <- features[which(features$merged_with == ""),]
    features <- as.data.frame(features)[,-ncolumns]

    ####now add mz and RT calibration information per feature

    RT_calibration_vals <- as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features),NA))
    colnames(RT_calibration_vals) <- paste("RT_calibration",sample_list,sep=" ")
    for(c in 1:ncol(RT_calibration_vals))
    {
      RT_calibration_vals[,c] <- as.numeric(RT_calibration_vals[,c])
    }

    mz_calibration_vals <- as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features),NA))
    colnames(mz_calibration_vals) <- paste("mz_calibration",sample_list,sep=" ")
    for(c in 1:ncol(mz_calibration_vals))
    {
      mz_calibration_vals[,c] <- as.numeric(mz_calibration_vals[,c])
    }

    if(any(c(RT_calibration,mz_calibration) == T)) ###any of both calibrations should be done?
    {
      ###first: start by taking already available calibration information

      peptide_features <- which(features$Sequence != "")

      pep_seq <- str_split(features$Sequence,";",simplify = T)

      evidence <- evidence[which(evidence$Raw.file %in% sample_list),]

      ###evidence_peptides
      evidence_peptides <- unique(evidence$Sequence)
      evidence_peptides_indices <- vector(mode = "list", length = length(evidence_peptides))
      names(evidence_peptides_indices) <- evidence_peptides

      max <- length(evidence_peptides)
      pb <- winProgressBar(title = "Indexing peptides",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)

      ###grap already available information for RT and mz calibration from MaxQ output
      max <- length(peptide_features)
      pb <- winProgressBar(title = "Preparing RT and m/z calibrations per feature",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      counter <- 0
      for(i in peptide_features)
      {
        cur_mods <- str_split(features$Modifications[i],";|,",simplify = T)
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
            calc_RT <- as.numeric(str_split(features$Observed_RT[i],pattern = ";",simplify = T)) - features$RT[i]
            set(RT_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_RT)))
          }
          if(mz_calibration == T)
          {
            if(use_mz_at_max_int_for_correction == T) ###determine mz correction based on observed mz peaks
            {
              calc_mz <- as.numeric(str_split(features$Observed_mz[i],pattern = ";",simplify = T)) - features$m.z[i]
              set(mz_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_mz)))
            }else ###use MaxQ correction
            {
              calc_mz <- aggregate(evidence_sub[indices,"Uncalibrated...Calibrated.m.z..Da."],by=list(Sample=evidence_sub$Raw.file[indices]),FUN="mean",na.rm=T)
              set(mz_calibration_vals,as.integer(i),as.integer(match(calc_mz$Sample,sample_list)),value = as.list(c(calc_mz$x)))
            }
          }
        }
        counter <- counter + 1
        updatecounter <- updatecounter + 1
        if(updatecounter >= 100)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(counter/max))*(1-(counter/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, counter, label=paste( round(counter/max*100, 0)," % done (",counter,"/",max,", Time require: ",time_require,")",sep = ""))
        }

      }
      close(pb)

      if(mz_calibration == T & use_mz_at_max_int_for_correction==F)
      {
        ####train random forest model for predicting mz calibration

        suppressWarnings(suppressMessages(library(randomForest,quietly = T)))

        ###train a model per sample
        models <- list()
        trainsets <- list()
        evalsets <- list()
        if(any(colnames(evidence) == "Resolution"))
        {
          temp_data <- evidence[which(rowSums(is.na(evidence[,c("Raw.file","Retention.time","m.z","Charge","Resolution","Uncalibrated...Calibrated.m.z..Da.")])) == 0),c("Raw.file","Retention.time","m.z","Charge","Resolution","Uncalibrated...Calibrated.m.z..Da.")]
        }else
        {
          temp_data <- evidence[which(rowSums(is.na(evidence[,c("Raw.file","Retention.time","m.z","Charge","Uncalibrated...Calibrated.m.z..Da.")])) == 0),c("Raw.file","Retention.time","m.z","Charge","Uncalibrated...Calibrated.m.z..Da.")]
          temp_data$Resolution <- 0
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
            model <- randomForest(trainset[,c("Retention.time","m.z","Charge","Resolution")],trainset[,"Uncalibrated...Calibrated.m.z..Da."] , importance = TRUE,ntree = 100, mtry = 4,do.trace=F,nodesize = 100)
            evalset$predicted <- stats::predict(model, evalset, type = "response")
            trainset$predicted <- stats::predict(model, trainset, type = "response")
            models[[s]] <- model
            trainsets[[s]] <- trainset
            evalsets[[s]] <- evalset

            plot(trainset$predicted,trainset$Uncalibrated...Calibrated.m.z..Da.,main=paste(s,"- Train set (80% of data)"),xlab="Predicted m/z calibration",ylab="MaxQ determined m/z calibration")
            fit <- lm(trainset$Uncalibrated...Calibrated.m.z..Da.~trainset$predicted)
            abline(fit)
            posx <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.8
            posy <- par("usr")[3] + (par("usr")[4]-par("usr")[3])*0.2
            text(posx,posy,paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))
            fit_train <- fit

            plot(evalset$predicted,evalset$Uncalibrated...Calibrated.m.z..Da.,main=paste(s,"- Validation set (20% of data)"),xlab="Predicted m/z calibration",ylab="MaxQ determined m/z calibration")
            fit <- lm(evalset$Uncalibrated...Calibrated.m.z..Da.~evalset$predicted)
            abline(fit)
            posx <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.8
            posy <- par("usr")[3] + (par("usr")[4]-par("usr")[3])*0.2
            text(posx,posy,paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))
            fit_eval <- fit

            print(paste(s,": R^2 train-set=",round(as.numeric(summary(fit_train)[8]),digits=2),", R^2 eval-set=",round(as.numeric(summary(fit_eval)[8]),digits=2),sep=""))
          }
        }

        QC_data[["mz_calibration_models"]] <- models
        QC_data[["mz_calibration_train_sets"]] <- trainsets
        QC_data[["mz_calibration_evaluation_sets"]] <- evalsets

        ####now predict mz for missing mzcalibrations of features

        ####get ion indices per feature
        all_ion_indices <- as.data.frame(matrix(ncol=6,nrow=nrow(allpeptides)))
        colnames(all_ion_indices) <- c("Feature","Raw.file","Retention.time","m.z","Charge","Resolution")
        all_ion_indices$Feature <- as.numeric(all_ion_indices$Feature)
        all_ion_indices$Raw.file <- as.character(all_ion_indices$Raw.file)
        all_ion_indices$Retention.time <- as.numeric(all_ion_indices$Retention.time)
        all_ion_indices$m.z <- as.numeric(all_ion_indices$m.z)
        all_ion_indices$Charge <- as.integer(all_ion_indices$Charge)
        all_ion_indices$Resolution <- as.numeric(all_ion_indices$Resolution)
        cur_index = 1

        max <- nrow(features)
        pb <- winProgressBar(title = "Prepare features for m/z-calibration prediction by RF",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        for(i in 1:nrow(features))
        {
          ###start with extracting ion_indices per feature
          cur_ion_indices <- as.numeric(str_split(features$ion_indices[i],",",simplify = T))
          set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),1L,i)
          set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),2L,allpeptides[cur_ion_indices,"Raw.file"])
          set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),as.integer(3:6),value = as.list(allpeptides[cur_ion_indices,c("Retention.time","m.z","Charge","Resolution")]))

          cur_index <- cur_index+length(cur_ion_indices)

          updatecounter <- updatecounter + 1
          if(updatecounter >= 100)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(i/max))*(1-(i/max))
            td <- seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)

        all_ion_indices <- all_ion_indices[1:(cur_index-1),]

        raw_files <- sample_list

        ####get mean ion propperties per feature to predict m/z calibration per feature for each sample
        median_feature_properties <- all_ion_indices[,-2]#[which(all_ion_indices$Raw.file == raw_files[c]),-2]

        ###determine mean values per feature
        median_feature_properties <- aggregate(median_feature_properties[,-1],by=list(Feature=median_feature_properties$Feature),FUN="median",na.rm=T)

        QC_data[["mz_calibration_median_feature_properties"]] <- median_feature_properties

        max <- ncol(mz_calibration_vals)
        pb <- winProgressBar(title = "Predict mz calibrations using RandomForest models",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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

              set(mz_calibration_vals,as.integer(temp_data$Feature),as.integer(c),value = prediction)
            }


          }
          updatecounter <- updatecounter + 1
          if(updatecounter >= 1)
          {
            time_elapsed <- difftime(Sys.time(),start_time,units="secs")
            time_require <- (time_elapsed/(c/max))*(1-(c/max))
            td <- seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

            updatecounter <- 0
            setWinProgressBar(pb, c, label=paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
          }
        }
        close(pb)
      }else if(mz_calibration == T & use_mz_at_max_int_for_correction==T) ###use observed mz at max intensity calibration
      {
        ###get rowMedians of mz_calibrations
        feature_specific_correction_factor <- rowMedians(as.matrix(mz_calibration_vals),na.rm=T)
        ###examine deviation of expected feature_specific shift from observed shift
        expected_from_observed <- mz_calibration_vals/feature_specific_correction_factor
        ###sample specific shift
        sample_specific_correction_factor <- colMedians(log2(as.matrix(expected_from_observed)),na.rm=T)
        ###examine deviation of expected feature_specific shift from observed shift
        expected_from_observed <- mz_calibration_vals/((2^sample_specific_correction_factor)*feature_specific_correction_factor)
        ###plot per sample how well expected mz correction correlates with observed correction
        for(c in 1:length(sample_list))
        {
          sel <- which(abs(mz_calibration_vals[,c]) < 0.01)
          smoothScatter(mz_calibration_vals[sel,c],(2^sample_specific_correction_factor[c])*feature_specific_correction_factor[sel],xlim=c(-0.005,0.005),ylim=c(-0.005,0.005),main=sample_list[c],xlab="Observed mz - feature mz",ylab="Expected mz - feature mz")
          abline(a=0,b=1)
          fit <- lm((2^sample_specific_correction_factor[c])*feature_specific_correction_factor[sel]~mz_calibration_vals[sel,c])
          abline(fit,lty=2)
          text(0,0.004,paste("Rsq:",round(as.numeric(summary(fit)[8]),digits=2)))
        }
        ###store calibration factor
        QC_data[["mz_calibration_feature_specific_cor_factor"]] <- feature_specific_correction_factor
        QC_data[["mz_calibration_sample_specific_cor_factor"]] <- sample_specific_correction_factor

        ###replace missing values with expected correction factors
        for(c in 1:ncol(mz_calibration_vals))
        {
          selection <- which(is.na(mz_calibration_vals[,c]))
          mz_calibration_vals[,c][selection] <- ((2^sample_specific_correction_factor[c])*feature_specific_correction_factor[selection])
        }
      }
      mz_calibration_vals[is.na(mz_calibration_vals)] <- 0

      if(RT_calibration == T)
      {
        RT_alignment_GAM_models <- list()
        ###fit GAM and predict per sample.
        for(c in 1:ncol(RT_calibration_vals))
        {
          x <- features$RT[which(RT_calibration_vals[,c] != 0)]
          y <- RT_calibration_vals[,c][which(RT_calibration_vals[,c] != 0)]
          if(length(x) > 10) ###require at least 500 data points to perform fitting
          {
            ##try to fit an average generalised additive model to determine a RT dependent calibration curve
            gam <- gam(y ~ s(x), method = "REML")
            RT_alignment_GAM_models[[c]] <- gam
            x_pred <- seq(min(features$RT,na.rm=T), max(features$RT,na.rm=T), length.out = nrow(features))
            y_pred <- predict(gam, data.frame(x = x_pred))
            lim_stats <- boxplot.stats(y_pred)
            delta_y <- lim_stats$stats[4] - lim_stats$stats[2]
            ylim <- c(-max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))-delta_y,max(abs(c(lim_stats$stats[2],lim_stats$stats[4])))+delta_y)
            if(ylim[2] < 2)ylim <- c(-2,2)
            smoothScatter(x,y,ylab="Observed RT calibration",main=sample_list[c],xlab="RT [min]",ylim=ylim)
            lines(x_pred,y_pred,col="red")
            legend("topright",legend="GAM",lty=1,col="red")

            ##predict RT correction for all features in current sample
            y_pred <- predict(gam, data.frame(x = features$RT))

            ###use GAM to predict RT correction for missing features in current sample
            selection <- which(is.na(RT_calibration_vals[,c]))
            RT_calibration_vals[,c][selection] <- y_pred[selection]
          }else ###not enough observations ... skiü
          {
            print(paste(sample_list[c],"- Not enough peptide observations for RT-GAM fitting..."))
          }

        }
        # ###use median of RT_calibration per sample for NAs
        # for(c in 1:ncol(RT_calibration_vals))
        # {
        #   RT_calibration_vals[,c][is.na(RT_calibration_vals[,c])] <- median(RT_calibration_vals[,c],na.rm=T)
        # }
        QC_data[["RT_calibration_GAM_models"]] <- RT_alignment_GAM_models
      }
      RT_calibration_vals[is.na(RT_calibration_vals)] <- 0

    }else ###no calibration should be done
    {
      RT_calibration_vals[is.na(RT_calibration_vals)] <- 0
      mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
    }

    features <- cbind(features,RT_calibration_vals,mz_calibration_vals)

    write.table(features,paste("Temporary_files\\Features_aligned_merged",output_file_names_add,".txt",sep=""),row.names = F)

    save(QC_data,file = paste("Temporary_files\\Feature_alignment_QC_data.RData",sep=""))

    dev.off()
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
add_isotope_features <- function(path_to_features,feature_table_file_name="Features_aligned_merged_IceR_analysis.txt",min_observations=0)
{
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  features <- read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  if(!any(grepl("_i",features$Feature_name)))
  {
    ###isotope + 1 features
    isotope_features <- features
    isotope_features <- isotope_features[which(!grepl(";",isotope_features$Sequence)),]
    isotope_features <- isotope_features[which(as.numeric(as.character(isotope_features$num_matches)) >= min_observations),]

    if(nrow(isotope_features)>0)
    {
      quantile25_score <- boxplot.stats(as.numeric(as.character(isotope_features$mean_Scores)))$stats[2]
      isotope_features <- isotope_features[which(as.numeric(as.character(isotope_features$mean_Scores)) >= quantile25_score),]
      isotope_features$m.z <- ((isotope_features$m.z*isotope_features$Charge)+1.002054)/isotope_features$Charge
      delta_mz <- isotope_features$m.z_range_max-isotope_features$m.z_range_min
      isotope_features$m.z_range_max <- isotope_features$m.z+(delta_mz/2)
      isotope_features$m.z_range_min <- isotope_features$m.z-(delta_mz/2)
      isotope_features$Feature_name <- paste(isotope_features$Feature_name,"_i",sep="")
      features <- rbind(features,isotope_features)
      write.table(features,paste("Temporary_files\\",feature_table_file_name,sep=""),row.names = F)
      print(paste("Added",nrow(isotope_features),"isotope features."))
    }else
    {
      print(paste("None of the peptides was observed in at least",min_observations,"samples. No isotope features were added."))
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
  options(warn = -1)
  suppressWarnings(suppressMessages(library(Peptides,quietly = T)))
  suppressWarnings(suppressMessages(library(randomForest,quietly = T)))
  suppressWarnings(suppressMessages(library(seqinr,quietly = T)))
  suppressWarnings(suppressMessages(library(cleaver,quietly = T)))
  RT_calibration=T
  mz_calibration=T

  pdf("Temporary_files\\Potentially_missed_peptides_QC_plots.pdf")

  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  features <- read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  features <- features[which(!grepl("_i",features$Feature_name)),]

  all_peptides <- unique(as.character(str_split(features$Sequence,";",simplify = T)))
  all_proteins <- unique(features$Protein)
  all_proteins <- all_proteins[which(all_proteins != "" & !grepl(";|\\|",all_proteins))]

  features <- features[which(!grepl(";|\\|",features$Sequence)),]

  ###filter for high qualtiy identifications
  features <- features[which(as.numeric(as.character(features$mean_Scores)) >= 50),]

  ###summarize relevant parameters
  print("Prepare data for modelling peptide RT")
  aaComp <- aaComp(features$Sequence)
  aaComp <- data.frame(matrix(unlist(aaComp), nrow=length(aaComp), byrow=T))[,1:9]
  colnames(aaComp) <- c("Tiny","Small","Aliphatic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")

  peptide_parameters <- data.frame(hydrophobicities=hydrophobicity(features$Sequence),
                                   length=nchar(as.character(features$Sequence)),
                                   charge=as.numeric(as.character(features$Charge)),
                                   aaComp,
                                   RT=as.numeric(as.character(features$RT)))
  ###Train random forest model on the parameters
  random_selection <- sample(1:nrow(peptide_parameters),size = 0.8*nrow(peptide_parameters),replace = F)
  trainset <- peptide_parameters[random_selection,]
  evalset <- peptide_parameters[-random_selection,]
  print("Perform RF modelling for peptide RT")
  model <- randomForest(trainset[,1:12],trainset[,13], importance = TRUE,ntree = 200, mtry = 4,do.trace=F,nodesize = floor(nrow(trainset)/1000))

  evalset$predicted <- stats::predict(model, evalset[,1:12], type = "response")
  trainset$predicted <- stats::predict(model, trainset[,1:12], type = "response")

  smoothScatter(trainset$predicted,trainset$RT,main="Train set (80% of data)",xlab="Predicted RT",ylab="Observed RT")
  fit <- lm(trainset$RT~trainset$predicted)
  abline(fit)
  posx <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.8
  posy <- par("usr")[3] + (par("usr")[4]-par("usr")[3])*0.2
  text(posx,posy,paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))

  smoothScatter(evalset$predicted,evalset$RT,main="Eval set (80% of data)",xlab="Predicted RT",ylab="Observed RT")
  fit <- lm(evalset$RT~evalset$predicted)
  abline(fit)
  posx <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.8
  posy <- par("usr")[3] + (par("usr")[4]-par("usr")[3])*0.2
  text(posx,posy,paste("R^2:",round(as.numeric(summary(fit)[8]),digits=2)))

  deviation_predicted_true_RT <- rowSds(cbind(evalset$RT,evalset$predicted),na.rm=T)
  RT_deviation_cutoff <- boxplot.stats(deviation_predicted_true_RT)$stats[4]

  ###Perform in-silico digestion of all identified proteins
  fasta <- read.fasta(path_to_fasta,seqtype = "AA",as.string = T)
  Sequences <- data.frame(matrix(unlist(fasta), nrow=length(fasta), byrow=T))
  colnames(Sequences) <- "Sequence"
  Sequences$Sequence <- as.character(Sequences$Sequence)
  fasta_headers_boundaries_uid <- gregexpr("\\|",names(fasta))
  fasta_headers_boundaries_uid <- data.frame(matrix(unlist(fasta_headers_boundaries_uid), nrow=length(fasta_headers_boundaries_uid), byrow=T))
  UIDs <- substr(names(fasta),start = fasta_headers_boundaries_uid[,1]+1,stop = fasta_headers_boundaries_uid[,2]-1)
  rownames(Sequences) <- UIDs
  Sequences <- Sequences[all_proteins,,drop=F]

  insilicodigest <- cleave(Sequences$Sequence,"trypsin")
  names(insilicodigest) <- rownames(Sequences)

  ###Summarize potential peptides with a length of at least n and at max m amino acids
  potential_new_feature_peptides <- vector(mode="character",length = nrow(features))
  potential_new_feature_peptides_ID <- vector(mode="character",length = nrow(features))
  count <- 1
  max <- length(insilicodigest)
  pb <- winProgressBar(title = "Determine potential missed peptide features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)
  Seq_to_ID <- data.frame(Sequence=potential_new_feature_peptides,ID=potential_new_feature_peptides_ID)
  potential_new_feature_peptides <- unique(potential_new_feature_peptides)
  wrong_seq <- which(regexpr("U",potential_new_feature_peptides) != -1)
  if(length(wrong_seq)>0)potential_new_feature_peptides <- potential_new_feature_peptides[-wrong_seq]
  ###Determine Mass and RT per potential new peptide feature
  aaComp <- aaComp(potential_new_feature_peptides)
  aaComp <- data.frame(matrix(unlist(aaComp), nrow=length(aaComp), byrow=T))[,1:9]
  colnames(aaComp) <- c("Tiny","Small","Aliphatic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")
  new_peptide_parameters <- data.frame(hydrophobicities=hydrophobicity(potential_new_feature_peptides),
                                       length=nchar(potential_new_feature_peptides),
                                       charge=2,
                                       aaComp)

  new_peptide_parameters$predicted <- stats::predict(model, new_peptide_parameters[,1:12], type = "response")

  potential_new_feature_peptides <- data.frame(Sequence=potential_new_feature_peptides,
                                               ID=Seq_to_ID$ID[match(potential_new_feature_peptides,Seq_to_ID$Sequence)],
                                               Mass=mw(potential_new_feature_peptides, monoisotopic = T),
                                               RT=new_peptide_parameters$predicted)

  ###Now try to see if corresponding features with +2 or +3 charge wera already detected by MaxQ
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  load("allpeptides.RData")
  setwd(path_to_features)

  select_dataframe_rows = function(ds, sel) {
    cnames = colnames(ds)
    rnames = rownames(ds)
    ds = data.frame(ds[sel,])
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

  Indices <- as.data.frame(matrix(ncol=4,nrow=num_windows))
  colnames(Indices) <- c("RT_start","RT_end","Row_start","Row_end")
  Indices$RT_start <- as.numeric(Indices$RT_start)
  Indices$RT_end <- as.numeric(Indices$RT_end)
  Indices$Row_start <- as.numeric(Indices$Row_start)
  Indices$Row_end <- as.numeric(Indices$Row_end)

  print("Indexing MaxQ features")
  max <- nrow(Indices)
  pb <- winProgressBar(title = "Indexing MaxQ features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          set(Indices,i = as.integer(i),j=as.integer(1:4),value=as.list(c(start,end,min(inds),max(inds))))
        }
      }
    }


    updatecounter <- updatecounter + 1
    if(updatecounter >= 1)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(i/max))*(1-(i/max))
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
  ###were detected, take this charge for which more MaxQ features were found
  load("Temporary_files\\Feature_alignment_QC_data.RData")
  borders_RT <- QC_data$Feature_alignment_windows$RT_window
  borders_m.z <- QC_data$Feature_alignment_windows$mz_window

  features_unknown <- as.data.frame(matrix(ncol=21,nrow=nrow(potential_new_feature_peptides)))
  colnames(features_unknown) <- c("Feature_name","m.z","RT","Sequence","Protein","MSMS.Scan.Numbers","m.z_range_min","m.z_range_max","RT_range_min","RT_range_max","ion_indices","count_ion_indices","Mass","Charge","mean_Scores","num_matches","Modifications","RT_length","Observed_RT","Observed_mz","Observed_score")
  features_unknown <- as.data.table(features_unknown)
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
  features_unknown$mean_Scores <- as.character(features_unknown$ion_indices)
  features_unknown$num_matches <- as.character(features_unknown$num_matches)
  features_unknown$Modifications <- as.character(features_unknown$Modifications)
  features_unknown$RT_length <- as.numeric(features_unknown$RT_length)
  features_unknown$Observed_RT <- as.character(features_unknown$Observed_RT)
  features_unknown$Observed_mz <- as.character(features_unknown$Observed_mz)
  features_unknown$Observed_score <- as.character(features_unknown$Observed_score)

  count_features <- 0
  total_count_features <- as.numeric(gsub("Feature_","",features$Feature_name[nrow(features)]))

  sample_list <- gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration",colnames(features)))])

  print("Search for potential missed peptide features")
  max <- nrow(potential_new_feature_peptides)
  pb <- winProgressBar(title = "Search for missed peptide features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 350)
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
      sub <- rbindlist(data_frags[search_range])

      sub <- sub[which(sub$Mass >= (potential_new_feature_peptides$Mass[i] - 0.001) & sub$Mass <= (potential_new_feature_peptides$Mass[i] + 0.001)),]

      if(length(unique(sub$Raw.file))>=min_observations)
      {
        median_RT <- median(sub$Retention.time,na.rm=T)
        median_mz <- median(sub$m.z,na.rm=T)
        count_charges <- plyr::count(sub$Charge)
        max_count_charge <- count_charges[which(count_charges$freq == max(count_charges$freq))[1],1]

        sub <- sub[which(sub$m.z >= median_mz+borders_m.z[1] & sub$m.z <= median_mz+borders_m.z[2] & sub$Retention.time >= median_RT+borders_RT[1] & sub$Retention.time <= median_RT+borders_RT[2] & sub$Charge == max_count_charge),]
        if(length(unique(sub$Raw.file))>=min_observations)
        {
          count_features <- count_features + 1
          total_count_features <- total_count_features + 1
          name <- paste("Feature",total_count_features,"pmp",sep="_")

          sequence <- as.character(potential_new_feature_peptides$Sequence[i])
          protein <- as.character(potential_new_feature_peptides$ID[i])
          msmsscan <- ""

          median_mass <- median(sub$Mass,na.rm=T)
          Charge <- max_count_charge

          scores <- 0
          num_matches <- length(unique(sub$Raw.file))
          Modifications <- ""

          RT_length <- median(sub$Retention.Length,na.rm=T)

          temp_RT <- aggregate(sub$Retention.time,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_RT <- paste(temp_RT[match(sample_list,temp_RT$Raw.file),2],collapse=";") ###order all observed RT according to sample_list

          temp_mz <- aggregate(sub$isotope_corrected_mz_at_max_int,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_mz <- paste(temp_mz[match(sample_list,temp_mz$Raw.file),2],collapse=";") ###order all observed mz according to sample_list

          temp_score <- aggregate(sub$Score,by=list(Raw.file=sub$Raw.file),FUN=mean,na.rm=T)
          Observed_score <- paste(temp_score[match(sample_list,temp_score$Raw.file),2],collapse=";") ###order all observed score according to sample_list

          matches <- sub$allpeptides_rowname

          set(x = features_unknown,i = as.integer(count_features),j = as.integer(c(1,4,5,6,11,15,16,17,19,20,21)),value = as.list(c(name,
                                                                                                                                    sequence,
                                                                                                                                    protein,
                                                                                                                                    msmsscan,
                                                                                                                                    paste(matches,collapse = ","),
                                                                                                                                    paste(scores,collapse=";"),
                                                                                                                                    paste(num_matches,collapse=";"),
                                                                                                                                    Modifications,
                                                                                                                                    Observed_RT,
                                                                                                                                    Observed_mz,
                                                                                                                                    Observed_score
          )))
          set(x = features_unknown,i = as.integer(count_features),j = as.integer(c(2,3,7,8,9,10,12,13,14,18)),value = as.list(c(
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
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,") Found: ",count_features,sep = ""))
    }

  }
  close(pb)

  features_unknown <- features_unknown[1:count_features,]
  print(paste("Detected",count_features,"potentially missed peptides (PMPs)"))
  features_unknown$num_diff_charges <- 1

  ####now add mz and RT calibration information per feature

  RT_calibration_vals <- as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features_unknown),NA))
  colnames(RT_calibration_vals) <- paste("RT_calibration",sample_list,sep=".")
  for(c in 1:ncol(RT_calibration_vals))
  {
    RT_calibration_vals[,c] <- as.numeric(RT_calibration_vals[,c])
  }

  mz_calibration_vals <- as.data.frame(matrix(ncol=length(sample_list),nrow=nrow(features_unknown),NA))
  colnames(mz_calibration_vals) <- paste("mz_calibration",sample_list,sep=".")
  for(c in 1:ncol(mz_calibration_vals))
  {
    mz_calibration_vals[,c] <- as.numeric(mz_calibration_vals[,c])
  }


  ###use observed RTs and fill missing values with median RT calibrations per sample
  if(RT_calibration == T)
  {
    for(i in 1:nrow(features_unknown))
    {
      calc_RT <- as.numeric(str_split(features_unknown$Observed_RT[i],pattern = ";",simplify = T)) - features_unknown$RT[i]
      set(RT_calibration_vals,as.integer(i),as.integer(1:length(sample_list)),value = as.list(c(calc_RT)))
    }

    ###fill missing values with median RT_calibration per sample
    for(c in 1:ncol(RT_calibration_vals))
    {
      RT_calibration_vals[,c][is.na(RT_calibration_vals[,c])] <- median(RT_calibration_vals[,c],na.rm=T)
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
    pb <- winProgressBar(title = "Preparing m/z calibrations per missed peptide feature",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(features_unknown))
    {
      observed_RT <- as.numeric(str_split(features_unknown$Observed_RT[i],";",simplify = T))
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
        sub <- rbindlist(data_frags[search_range])
        temp <- sub[which(sub$allpeptides_rowname %in% as.character(str_split(features_unknown$ion_indices[i],",",simplify = T))),]
        temp$Uncalibrated...Calibrated.m.z..Da. <- temp$Uncalibrated.m.z-temp$m.z
        calc_mz <- aggregate(temp$Uncalibrated...Calibrated.m.z..Da.,by=list(Sample=temp$Raw.file),FUN="mean",na.rm=T)
        set(mz_calibration_vals,as.integer(i),as.integer(match(calc_mz$Sample,sample_list)),value = as.list(c(calc_mz$x)))

      }

      updatecounter <- updatecounter + 1
      if(updatecounter >= 10)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }

    }
    close(pb)

    ####now predict mz for missing mzcalibrations of features
    setwd(paste(path_to_features,"\\Temporary_files",sep=""))
    load("allpeptides.RData")
    setwd(path_to_features)

    ###all following steps are only required to get the median resolution per feature
    ###we require median RT, mz, resolution and charge as a rough estimation of the corresponding parameters in the samples for which we don?t have the respective feature identified
    ####get ion indices per feature
    all_ion_indices <- as.data.frame(matrix(ncol=6,nrow=nrow(allpeptides)))
    colnames(all_ion_indices) <- c("Feature","Raw.file","Retention.time","m.z","Charge","Resolution")
    all_ion_indices$Feature <- as.numeric(all_ion_indices$Feature)
    all_ion_indices$Raw.file <- as.character(all_ion_indices$Raw.file)
    all_ion_indices$Retention.time <- as.numeric(all_ion_indices$Retention.time)
    all_ion_indices$m.z <- as.numeric(all_ion_indices$m.z)
    all_ion_indices$Charge <- as.integer(all_ion_indices$Charge)
    all_ion_indices$Resolution <- as.numeric(all_ion_indices$Resolution)
    cur_index = 1

    max <- nrow(features_unknown)
    pb <- winProgressBar(title = "Prepare missed peptide features for m/z-calibration prediction by RF",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0
    for(i in 1:nrow(features_unknown))
    {
      ###start with extracting ion_indices per feature
      cur_ion_indices <- as.numeric(str_split(features_unknown$ion_indices[i],",",simplify = T))
      set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),1L,i)
      set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),2L,allpeptides[cur_ion_indices,"Raw.file"])
      set(all_ion_indices,i = as.integer(cur_index:(cur_index-1+length(cur_ion_indices))),as.integer(3:6),value = as.list(allpeptides[cur_ion_indices,c("Retention.time","m.z","Charge","Resolution")]))

      cur_index <- cur_index+length(cur_ion_indices)

      updatecounter <- updatecounter + 1
      if(updatecounter >= 100)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(i/max))*(1-(i/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)

    all_ion_indices <- all_ion_indices[1:(cur_index-1),]

    raw_files <- sample_list

    ####get mean ion propperties per feature to predict m/z calibration per feature for each sample
    median_feature_properties <- all_ion_indices[,-2]#[which(all_ion_indices$Raw.file == raw_files[c]),-2]

    ###determine mean values per feature
    median_feature_properties <- aggregate(median_feature_properties[,-1],by=list(Feature=median_feature_properties$Feature),FUN="median",na.rm=T)

    ###perform prediction

    models <- QC_data$mz_calibration_models

    max <- ncol(mz_calibration_vals)
    pb <- winProgressBar(title = "Predict mz calibrations using RandomForest models",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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

          set(mz_calibration_vals,as.integer(temp_data$Feature),as.integer(c),value = prediction)
        }
      }
      updatecounter <- updatecounter + 1
      if(updatecounter >= 1)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(c/max))*(1-(c/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, c, label=paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
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
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  features <- read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  ##add new PMPs
  features <- rbind(features,features_unknown)

  ###Save features list
  write.table(features,paste("Temporary_files\\",feature_table_file_name,sep=""),row.names = F)

  options(warn=0)
  dev.off()
  crap <- gc(F)
}

#' Perform quantification of IceR features
#' @param path_to_features Path to folder where results of align_features() are stored
#' @param path_to_mzXML Path to folder containing mzXML files of samples
#' @param path_to_MaxQ_output Path to folder containing MaxQuant outputs (txt folder containing at least allpeptides.txt, evidence.txt, peptides.txt and proteinGroups.txt)
#' @param path_to_output Path to folder where IceR results should be stored
#' @param feature_table_file_name File name which contains align_features() results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param output_file_names_add IceR result name tag. By default IceR_analysis
#' @param RT_calibration Boolean value indicating if corrected RT should be used during peak detection, selection and DICE, By default set to T.
#' @param mz_calibration Boolean value indicating if corrected m/z should be used during peak detection, selection and DICE, By default set to T.
#' @param abundance_estimation_correction Boolean value indicating if resulting peptide abundances should be corrected using MaxQuant results as a reference. By default set to T.
#' @param Quant_pVal_cut Numeric value used as diagnostic cutoff border for visualization of significances of ion accumulation per IceR feature quantification. Furthermore, used as cutoff to filter +1-isotopic IceR features with significant accumulation of ions. By default set to 0.05.
#' @param n_cores Numeric value specifying on how many CPU cores tasks should be parallelized. By default set to 2.
#' @param kde_resolution Numeric value specifying number of grid points per dimension. By default set to 50.
#' @param num_peaks_store Numeric value specifying number of 2D peaks to be stored during peak detection. By default set to 5.
#' @param plot_2D_peak_detection Boolean value indicating if for every feature quantification the determined kernel density estimations and detected peaks should be visualized and stored. By default set to F.
#' @param alignment_variability_score_cutoff Numeric value specifying significance cutoff to distinguish which features show high RT- or m/z-variability of selected peaks between samples. By default set to 0.05. All features showing significant general variability (variability score < alignment_variability_score_cutoff) are excluded.
#' @param alignment_scores_cutoff Numeric value specifying significance cutoff to distinguish which samples show high RT- or m/z-variability of selected peaks for respective IceR feature. By default set to 0.05. All samples showing significant peak variability for respective IceR feature (variability score < alignment_scores_cutoff) are excluded (quantification set to NA).
#' @param mono_iso_alignment_cutoff Numeric value specifying significance cutoff to distinguish which samples show high RT- or m/z-variability of selected peaks for +1-isotopic from corresponding monoisotopic IceR feature. By default set to 0.05. All samples showing significant peak variability between selected +1-isotopic and monoisotopic IceR features (variability score < mono_iso_alignment_cutoff) are excluded (quantification of +1-isotopic feature set to NA).
#' @param calc_peptide_LFQ Boolean value specifying if multiply peptide quantification data for same peptide sequence (multiply charge states, isotope-states) should be aggregated using the MaxLFQ algorithm. By default set to F.
#' @param calc_protein_LFQ Boolean value specifying if protein quantification should be additionally performed by peptide quantification aggregation using the MaxLFQ algorithm. By default set to T.
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
requantify_features <- function(path_to_features,path_to_mzXML,path_to_MaxQ_output,feature_table_file_name="Features_aligned_merged_IceR_analysis.txt",output_file_names_add="IceR_analysis",RT_calibration=T,mz_calibration=T,abundance_estimation_correction = T,Quant_pVal_cut=0.05,n_cores=2,kde_resolution = 50,num_peaks_store = 5,plot_2D_peak_detection=F,alignment_variability_score_cutoff=0.05,alignment_scores_cutoff=0.05,mono_iso_alignment_cutoff=0.05,calc_peptide_LFQ=F,calc_protein_LFQ=T)
{
  options(warn=-1)
  suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
  multiply_intensity_count=F
  peak_detection=T

  QC_data <- list() ##here relevant qc data is stored and finally saved as RData which can be used for re-generating plots

  mean_background_intensity <- NA
  sd_background_intensity <- NA
  n_background_intensity <- NA

  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  features <- read.table(feature_table_file_name,header = T)
  setwd(path_to_features)

  ###Check peak widths and increase very small peak widths while reduce extreme long outlier peak widths
  stats <- boxplot.stats(features$RT_length)
  features$RT_length[which(features$RT_length < stats$stats[2])] <- stats$stats[2] ##lower than 25 % quantile
  features$RT_length[which(features$RT_length > stats$stats[5])] <- stats$stats[5] ##if RT length > upper whisker --> shorten to upper whisker
  features$RT_length[is.na(features$RT_length)] <- stats$stats[3]

  ###Add decoy features
  features_decoy <- features[which(!grepl("_pmp|_i",features$Feature_name)),]
  features_decoy$Feature_name <- paste(features_decoy$Feature_name,"_d",sep="")

  RT_ranges <- features_decoy$RT_range_max-features_decoy$RT_range_min
  mz_ranges <- features_decoy$m.z_range_max-features_decoy$m.z_range_min
  median_RT_window <- median(features_decoy$RT_range_max-features_decoy$RT_range_min,na.rm=T)
  median_mz_window <- median(features_decoy$m.z_range_max-features_decoy$m.z_range_min,na.rm=T)

  features_decoy$RT <- features_decoy$RT + 5*median_RT_window
  features_decoy$m.z <- features_decoy$m.z + 5*median_mz_window

  features_decoy$m.z_range_min <- features_decoy$m.z-(mz_ranges/2)
  features_decoy$m.z_range_max <- features_decoy$m.z+(mz_ranges/2)

  features_decoy$RT_range_min <- features_decoy$RT-(RT_ranges/2)
  features_decoy$RT_range_max <- features_decoy$RT+(RT_ranges/2)

  ###adjust observed RT and mz accordingly
  Observed_RT <- as.matrix(str_split(features_decoy$Observed_RT,";",simplify = T))
  class(Observed_RT) <- "numeric"
  Observed_RT <- Observed_RT + 5*median_RT_window
  features_decoy$Observed_RT <- apply(Observed_RT, 1, paste, collapse=";")

  Observed_mz <- as.matrix(str_split(features_decoy$Observed_mz,";",simplify = T))
  class(Observed_mz) <- "numeric"
  Observed_mz <- Observed_mz + 5*median_mz_window
  features_decoy$Observed_mz <- apply(Observed_mz, 1, paste, collapse=";")

  features <- rbind(features,features_decoy)

  if(output_file_names_add != "")output_file_names_add <- paste("_",output_file_names_add,sep="")

  ##Step 1 - Summarize ion intensities per aligned MaxQuant feature
  ###Function to perform extraction on multiple threads. Depending on the available ram, the extraction can be run on several threads. However, the task is using much memory so that it is recommended to just run 2 threads in parallel if only 16 gb of ram are available.
  extract_intensities_worker <- function(Sample_IDs,features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection,n_cores,ion_intensity_cutoff=F,mean_background_ion_intensity_model = NA,sd_background_ion_intensity = NA,peak_min_ion_count=NA,kde_resolution=25,num_peaks_store=5,plots=F)
  {
    suppressWarnings(suppressMessages(library(doParallel,quietly = T)))

    ####Function to extract intensities in respective extracted .RData for a list of selected features (or all features).
    ####indexing_RT_window defines how large each RT indexing window is to speed up subsetting. 0.5 min was observed to be good
    get_intensities <- function(Sample_ID,path,features_select,indexing_RT_window=0.1,RT_calibration,mz_calibration,peak_detection,ion_intensity_cutoff,mean_background_ion_intensity_model,sd_background_ion_intensity,peak_min_ion_count,kde_resolution,num_peaks_store,plots)
    {
      suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
      suppressWarnings(suppressMessages(library(stringr,quietly = T)))

      ###if no peak selection should be performed, extract all ions within the expected feature windows
      feature_no_peak_selection <- function(all_ion_data,selected_features,cur_sample,include_isotope_patterns = F,num_peaks_store=0)
      {
        peak_ion_data_list <- list()

        ###define RT and mz window
        RT_expected <- selected_features$RT + selected_features[,paste("RT_calibration.",cur_sample,sep="")]
        mz_expected <- selected_features$m.z + selected_features[,paste("mz_calibration.",cur_sample,sep="")]

        RT_window <- c(RT_expected - (selected_features$RT_length/2),
                       mz_expected + (selected_features$RT_length/2))

        mz_window <- c(selected_features$m.z_range_min + selected_features[,paste("mz_calibration.",cur_sample,sep="")],
                       selected_features$m.z_range_max + selected_features[,paste("mz_calibration.",cur_sample,sep="")])

        ###select relevant ions
        ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window[1] &
                                         all_ion_data$m.z <= mz_window[2] &
                                         all_ion_data$RT >= RT_window[1] &
                                         all_ion_data$RT <= RT_window[2]),]
        ion_data$isotope=0

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
                                                                      Peak_info=data.frame(RT=RT_expected,
                                                                                           RT_window_lower=RT_window[1],
                                                                                           RT_window_upper=RT_window[2],
                                                                                           mz=mz_expected,
                                                                                           mz_window_lower=mz_window[1],
                                                                                           mz_window_upper=mz_window[2],
                                                                                           density=0,
                                                                                           known_peak=0,
                                                                                           Peak=num_peaks_store+1,
                                                                                           ion_count=nrow(ion_data)))
        names(peak_ion_data_list) <- "Standard"
        return(peak_ion_data_list)

      }

      ###if peak selection should be performed
      feature_2D_peak_selection <- function(all_ion_data,selected_features,cur_sample,known_RT,delta_mz,delta_rt,RT_window_expand_factor=5,mz_window_expand_factor=4,include_isotope_patterns = F,n_raster_dens_matrix=50,local_maxima_k=3,max_delta_RT=2,max_delta_mz=0.005,peak_min_ion_count=5,RT_bw=0.5,mz_bw=0.002,num_peaks_store=5,plot=F,auto_adjust_kde_resolution=T)
      {

        suppressWarnings(suppressMessages(library(MASS,quietly = T)))
        suppressWarnings(suppressMessages(library(raster,quietly = T)))
        suppressWarnings(suppressMessages(library(ggtern,quietly = T)))
        peak_ion_data_list <- list()

        graph <- NA ##here we store graphical output if wanted

        ###define RT and mz window
        RT_correction <- as.numeric(selected_features[paste("RT_calibration.",cur_sample,sep="")])
        mz_correction <- as.numeric(selected_features[paste("mz_calibration.",cur_sample,sep="")])

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

        ###filter for ions in the expanded window
        ion_data <- all_ion_data[which(all_ion_data$m.z >= mz_window_expanded[1] &
                                         all_ion_data$m.z <= mz_window_expanded[2] &
                                         all_ion_data$RT >= RT_window_expanded[1] &
                                         all_ion_data$RT <= RT_window_expanded[2]),]

        ion_data$isotope=0

        ##add isotope ions if wanted
        if(include_isotope_patterns == T)
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

          f2 <- kde2d(ion_data$RT, ion_data$m.z, n = n_raster_dens_matrix,
                      h = c(RT_bw, mz_bw),lims = c(RT_window_expanded[1],RT_window_expanded[2],mz_window_expanded[1],mz_window_expanded[2]))

          #area_per_cell <- ((max(f2$x)-min(f2$x)))*((max(f2$y)-min(f2$y)))/2

          ###determine which density (=z) will be selected to be at least exceeded --> upper whisker
          outlier_densities <- boxplot.stats(f2$z)$stats[5]

          ## Convert it to a raster object
          r <- raster(f2$z)
          extent(r) <- extent(c(0, length(f2$x), 0, length(f2$y)) + 0.5)

          ## Find the maximum value within the k-cell neighborhood of each cell
          f <- function(X) max(X, na.rm=TRUE)
          ww <- matrix(1, nrow=local_maxima_k, ncol=local_maxima_k) ## Weight matrix for cells in moving window
          localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)

          ## Get x-y coordinates of those cells that are local maxima
          r2 <- r==localmax
          maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))

          ## Remove maxima which are below the minimal outlier density
          maxima <- data.frame(RT=f2$x[n_raster_dens_matrix-maxXY[,2]+1],mz=f2$y[maxXY[,1]],density=f2$z[((maxXY[,1]-1)*n_raster_dens_matrix)+(n_raster_dens_matrix-maxXY[,2]+1)])
          maxima <- maxima[which(maxima$density > outlier_densities),]
          ###next filter for maxima with at least peak_min_ion_count ions
          maxima$estimated_count <- (maxima$density/sum(as.matrix(maxima$density)))*nrow(ion_data)

          RT_cut <- selected_features$RT_length/2#ifelse(!is.na(known_RT),selected_features$RT_length/2,delta_rt)
          selection <- which(maxima$estimated_count >= peak_min_ion_count | abs(maxima$RT-RT_expected)<=RT_cut & abs(maxima$mz-mz_expected)<=delta_mz)
          ###if no maxima are left after filtering still take topN closest peak further
          if(length(selection)==0 & nrow(maxima) > 0)
          {
            ###no maxima would be left see keep up to N peaks
            distances <- data.frame(dist_rt=maxima$RT-RT_expected,
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
          if(nrow(maxima)>1)
          {
            close_info <- as.data.frame(matrix(ncol=2,nrow=nrow(maxima),0))

            for(cl in 1:nrow(maxima))
            {
              temp <- sqrt((maxima$RT[cl]-maxima$RT)^2+((maxima$mz[cl]-maxima$mz)*500)^2)
              sel <- which(temp == min(temp[-cl]))[1]
              set(close_info,as.integer(cl),as.integer(1:2),value = list(temp[sel],sel))
            }

            dist_cut <- sqrt((delta_mz*500)^2+(selected_features$RT_length/4)^2) ###merge peaks with d_mz < 0.001 and d_RT < peak_width/2
            if(any(close_info$V1 < dist_cut))
            {
              sel <- which(close_info$V1 < dist_cut)

              maxima_add <- as.data.frame(matrix(nrow=length(sel),ncol=4,0))
              colnames(maxima_add) <- colnames(maxima)
              ###find mean RT and mz of maxima which should be merged wheigted by density
              for(cl in 1:length(sel))
              {
                RT_add <- weighted.mean(c(maxima$RT[sel[cl]],maxima$RT[close_info[sel[cl],2]]),
                                        c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                mz_add <- weighted.mean(c(maxima$mz[sel[cl]],maxima$mz[close_info[sel[cl],2]]),
                                        c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                density_add <- weighted.mean(c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]),
                                             c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                count_add <- weighted.mean(c(maxima$estimated_count[sel[cl]],maxima$estimated_count[close_info[sel[cl],2]]),
                                           c(maxima$density[sel[cl]],maxima$density[close_info[sel[cl],2]]))
                set(maxima_add,as.integer(cl),as.integer(1:4),value=list(RT_add,mz_add,density_add,count_add))
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
            distances <- data.frame(dist_rt=maxima$RT-RT_expected,
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
                                                           Peak_info=data.frame(RT=maxima_select$RT[i],
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
                                                                            Peak_info=data.frame(RT=RT_expected,
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
                image(f2,main=paste(cur_sample,"-",selected_features$Feature_name),xlab="RT",ylab="m/z",xlim=RT_window_expanded,ylim=mz_window_expanded)
                #points(maxima_select$RT,maxima_select$mz)
                text(maxima_select$RT,maxima_select$mz,1:nrow(maxima_select),col="black")
                ###expected window
                rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2)
                points(RT_expected,mz_expected,pch=4)
                ###adjusted window

                rect(maxima_select$RT-RT_window_width/2,maxima_select$mz-mz_window_width,maxima_select$RT+RT_window_width/2,maxima_select$mz+mz_window_width,lty=2,border=ifelse(maxima_select$known_peak==1,"green","darkgrey"))

                ###add label indicating intensity within the window
                for(i in 1:length(peak_ion_data_list))
                {
                  sum_intensity <- sum(sum(peak_ion_data_list[[i]]$ion_data$Intensity))
                  if(sum_intensity > 0)sum_intensity <- log2(sum_intensity)
                  if(peak_ion_data_list[[i]]$Peak_info$Peak != num_peaks_store+1)
                  {
                    text(peak_ion_data_list[[i]]$Peak_info$RT,peak_ion_data_list[[i]]$Peak_info$mz-(1.1*delta_mz),round(sum_intensity,2),col=ifelse(peak_ion_data_list[[i]]$Peak_info$known_peak==1,"green","darkgrey"))
                  }else
                  {
                    text(peak_ion_data_list[[i]]$Peak_info$RT,peak_ion_data_list[[i]]$Peak_info$mz+(1.1*delta_mz),round(sum_intensity,2),col="black")
                  }

                }
              }else ###store all relevant data for performing plotting later
              {
                peak_intensities <- as.data.frame(matrix(ncol=1,nrow=length(peak_ion_data_list)))
                peak_intensities[,1] <- as.numeric(peak_intensities[,1])

                for(i in 1:length(peak_ion_data_list))
                {
                  sum_intensity <- sum(sum(peak_ion_data_list[[i]]$ion_data$Intensity))
                  if(sum_intensity > 0)sum_intensity <- log2(sum_intensity)
                  set(peak_intensities,as.integer(i),as.integer(1),sum_intensity)
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
                                                                            Peak_info=data.frame(RT=RT_expected,
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
                                                                          Peak_info=data.frame(RT=RT_expected,
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
                                                                        Peak_info=data.frame(RT=RT_expected,
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

        names(peak_ion_data_list)[which(names(peak_ion_data_list) != as.character(num_peaks_store+1))] <- paste("Peak_",names(peak_ion_data_list)[which(names(peak_ion_data_list) != as.character(num_peaks_store+1))],sep="")
        names(peak_ion_data_list)[which(names(peak_ion_data_list) == as.character(num_peaks_store+1))] <- "Standard"
        peak_ion_data_list$graph <- graph
        return(peak_ion_data_list)
      }


      `%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

      t.test2 <- function(...) ###T-test
      {
        obj<-try(t.test(...), silent=TRUE)
        if (is(obj, "try-error")) return(NA) else return(obj$p.value)
      }

      setwd(path)
      max <- 1
      print("Load spectra data")
      pb <- winProgressBar(title = paste("Load spectra data:",Sample_ID),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      load(paste(Sample_ID,"_all_ions.RData",sep="")) ###ions with RT and intensity in variable dat
      setWinProgressBar(pb, 1, label=paste( round(1/max*100, 0)," % done (",1,"/",max,")",sep = ""))
      close(pb)

      ###Indexing dat by RT windows to improve speed for subsetting
      num_windows <- ceiling(ceiling(max(dat$RT))*(1/indexing_RT_window))

      Indices <- as.data.frame(matrix(ncol=4,nrow=num_windows))
      colnames(Indices) <- c("RT_start","RT_end","Row_start","Row_end")
      Indices$RT_start <- as.numeric(Indices$RT_start)
      Indices$RT_end <- as.numeric(Indices$RT_end)
      Indices$Row_start <- as.numeric(Indices$Row_start)
      Indices$Row_end <- as.numeric(Indices$Row_end)

      print("Indexing intensities")
      max <- nrow(Indices)
      pb <- winProgressBar(title = paste("Indexing intensities:",Sample_ID),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
              set(Indices,i = as.integer(i),j=as.integer(1:4),value=as.list(c(start,end,min(inds),max(inds))))
            }
          }
        }


        updatecounter <- updatecounter + 1
        if(updatecounter >= 1)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
      ###calculate TIC per MS1 spectrum
      TIC_per_MS_spectrum <- dat[,list(sum=sum(Intensity)),by="RT"]

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

        Intensities <- as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
        colnames(Intensities) <- features_select$Feature_name

        ###if cut off pvalue is defined then store background quantifications separately
        if(ion_intensity_cutoff == T)
        {
          Intensities_signal_background <- as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
          colnames(Intensities_signal_background ) <- features_select$Feature_name
        }

        ###store data for calculating scoring per sample and feature
        delta_T1 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_T1) <- features_select$Feature_name
        delta_T2 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_T2) <- features_select$Feature_name
        delta_M1 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_M1) <- features_select$Feature_name
        delta_M2 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
        colnames(delta_M2) <- features_select$Feature_name
      }else
      {
        graph_peaks <- list()

        peaks_quant <- list()

        for(p in 1:(num_peaks_store+1)) ###1-num_peaks_store stores the quantification for top closest peaks while last stores total quantification of the window
        {
          Intensities <- as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
          colnames(Intensities) <- features_select$Feature_name

          ###if cut off pvalue is defined then store background quantifications separately
          if(ion_intensity_cutoff == T)
          {
            Intensities_signal_background <- as.data.frame(matrix(nrow=10,ncol=nrow(features_select),0))
            colnames(Intensities_signal_background ) <- features_select$Feature_name
          }
          ###store data for calculating scoring per sample and feature
          delta_T1 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
          colnames(delta_T1) <- features_select$Feature_name
          delta_T2 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
          colnames(delta_T2) <- features_select$Feature_name
          delta_M1 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
          colnames(delta_M1) <- features_select$Feature_name
          delta_M2 <- as.data.frame(matrix(nrow=1,ncol=nrow(features_select),0))
          colnames(delta_M2) <- features_select$Feature_name

          peaks_quant[[p]] <- list(Intensities=Intensities,
                                   Intensities_signal_background=Intensities_signal_background,
                                   delta_T1=delta_T1,
                                   delta_T2=delta_T2,
                                   delta_M1=delta_M1,
                                   delta_M2=delta_M2)
        }
        names(peaks_quant)[1:num_peaks_store] <- paste("Peak_",1:num_peaks_store,sep="")
        names(peaks_quant)[num_peaks_store+1] <- "Standard"
      }

      print("Extracting intensities")
      max <- nrow(features_select)
      pb <- winProgressBar(title = paste("Extracting intensities:",Sample_ID),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0
      end <- max

      iso_mult <- (1:3)*1.002054

      peak_width_stats <- boxplot.stats(features_select$RT_length)$stats
      delta_mz <- median(features_select$m.z - features_select$m.z_range_min,na.rm=T)
      delta_rt <- median(features_select$RT - features_select$RT_range_min,na.rm=T)/2

      known_RTs <- as.data.frame(str_split(features_select$Observed_RT,";",simplify = T))
      colnames(known_RTs) <- substr(colnames(features_select)[which(grepl("RT_calibration",colnames(features_select)))],16,1000)
      known_RTs[] <- lapply(known_RTs, function(x) as.numeric(as.character(x)))
      known_RTs_all <- known_RTs
      known_RTs <- known_RTs[,Sample_ID]

      RT_window_expand_factor <- 5 ###expand RT window around expected window for 2Dpeak detection

      if(plots == T)
      {
        dir.create(paste(path,"\\2Dpeakselection",sep=""))
        dir.create(paste(path,"\\2Dpeakselection\\",Sample_ID,sep=""))
      }

      for(i in 1:nrow(features_select))#nrow(features_select) #which(features_select$Protein == "P00579"))
      {
        ###extract all relevant ions in RT window
        RT_correction <- ifelse(RT_calibration == T,features_select[i,paste("RT_calibration.",Sample_ID,sep="")],0)

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

        if(length(search_range)>0)
        {
          sub <- rbindlist(data_frags[search_range])

          res <- NULL
          if(peak_detection == T )#& !grepl("_d",features_select$Feature_name[i])) ###use 2D peak detection
          {
            if(plots == T)
            {
              pdf(paste(path,"\\2Dpeakselection\\",Sample_ID,"\\",features_select$Feature_name[i],".pdf",sep=""))
            }
            res <- feature_2D_peak_selection(all_ion_data = sub,selected_features = features_select[i,],cur_sample = Sample_ID,known_RT = known_RTs[i],delta_mz = delta_mz,delta_rt=delta_rt,max_delta_RT = 2*delta_rt,max_delta_mz = 3*delta_mz,peak_min_ion_count = peak_min_ion_count,n_raster_dens_matrix = kde_resolution,num_peaks_store=num_peaks_store,plot = plots)
            if(plots == T)
            {
              dev.off()
            }
            if(plots==F) ###store graphs in list
            {
              #graph_peaks[[as.character(features_select$Feature_name[i])]] <- res$graph
            }
            res$graph <- NULL
          }else ###no peak detection
          {
            res <- feature_no_peak_selection(all_ion_data = sub,selected_features = features_select[i,],cur_sample = Sample_ID,num_peaks_store=num_peaks_store)
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

              ###determine which ions are showing an intensity above the background signal intensity
              if(ion_intensity_cutoff == T)
              {
                res$Standard$ion_data$signal_background <- ifelse(log2(res$Standard$ion_data$Intensity) > mean_background_ion_intensity_model[i,Sample_ID]+2*sd_background_ion_intensity[1,Sample_ID],1,0)

                res_int_signal <- res_int_total[which(res$Standard$ion_data$signal_background == 1)]
                res_mz_signal <- res_mz_total[which(res$Standard$ion_data$signal_background == 1)]
                res_rt_signal <- res_rt_total[which(res$Standard$ion_data$signal_background == 1)]
                res_isotope_signal <- res_isotope_total[which(res$Standard$ion_data$signal_background == 1)]
                length_data_signal <- length(res_int_signal)

                ###save signal Intensities
                if(length(res_int_signal)>0)
                {
                  set(Intensities,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_signal)),
                                                                                 length(res_int_signal),
                                                                                 mean(log10(res_int_signal),na.rm=T),
                                                                                 sd(log10(res_int_signal),na.rm=T),
                                                                                 m.z_window_final[1],
                                                                                 m.z_window_final[2],
                                                                                 RT_window_final[1],
                                                                                 RT_window_final[2])))
                }

                ###save signal+background Intensities
                if(length(res_int_total)>0)
                {
                  set(Intensities_signal_background,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                                   length(res_int_total),
                                                                                                   mean(log10(res_int_total),na.rm=T),
                                                                                                   sd(log10(res_int_total),na.rm=T),
                                                                                                   m.z_window_final[1],
                                                                                                   m.z_window_final[2],
                                                                                                   RT_window_final[1],
                                                                                                   RT_window_final[2])))
                }
              }else ###no discrimination between signal and background ions ## save all ions
              {
                ###save signal+background Intensities
                if(length(res_int_total)>0)
                {
                  set(Intensities,i=as.integer(1:8),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                 length(res_int_total),
                                                                                 mean(log10(res_int_total),na.rm=T),
                                                                                 sd(log10(res_int_total),na.rm=T),
                                                                                 m.z_window_final[1],
                                                                                 m.z_window_final[2],
                                                                                 RT_window_final[1],
                                                                                 RT_window_final[2])))
                }
              }

              ###feature alignment Scoreing based on algorithm from DeMix-Q (Zhang, 2016) - use all ions
              if(length_data_total > 0)
              {
                df <- data.frame(mz=res_mz_total,rt=res_rt_total,iso=res_isotope_total,int=res_int_total)
                df_mono_isotope <- df[which(df$iso==0),]
                df_isotope_1 <- df[which(df$iso==1),]
                if(nrow(df_mono_isotope)>0)
                {
                  ###deviation from consensus feature
                  max_mono <- which(df_mono_isotope$int == max(df_mono_isotope$int,na.rm=T))
                  set(delta_T1,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - (features_select$RT[i]))
                  set(delta_M1,i=1L,j=as.integer(i),mean(df_mono_isotope$mz[max_mono],na.rm=T) - (features_select$m.z[i]))
                  if(nrow(df_isotope_1)>0)
                  {
                    ###deviation from monoisotopic ion to M+1 isotope ion
                    max_iso_1 <- which(df_isotope_1$int == max(df_isotope_1$int,na.rm=T))
                    set(delta_T2,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - mean(df_isotope_1$rt[max_iso_1],na.rm=T))
                    set(delta_M2,i=1L,j=as.integer(i),(mean(df_mono_isotope$mz[max_mono],na.rm=T)+(iso_mult[1]/features_select$Charge[i])) - mean(df_isotope_1$mz[max_iso_1],na.rm=T))
                  }
                }
              }
            }else ###with peak detection
            {
              ###store up to num_peaks_store potential peak windows + the original expected window
              num_peaks <- length(which(names(res) != "Standard"))
              for(p in names(res))
              {
                m.z_window_final <- c(res[[p]]$Peak_info$mz_window_lower,res[[p]]$Peak_info$mz_window_upper)
                RT_window_final <- c(res[[p]]$Peak_info$RT_window_lower,res[[p]]$Peak_info$RT_window_upper)

                res_int_total <- res[[p]]$ion_data$Intensity
                res_mz_total <- res[[p]]$ion_data$m.z
                res_rt_total <- res[[p]]$ion_data$RT
                res_isotope_total <- res[[p]]$ion_data$isotope
                length_data_total <- length(res_int_total)

                ###determine which ions are showing an intensity above the background signal intensity
                if(ion_intensity_cutoff == T)
                {
                  res[[p]]$ion_data$signal_background <- ifelse(log2(res[[p]]$ion_data$Intensity) > mean_background_ion_intensity_model[i,Sample_ID]+2*sd_background_ion_intensity[1,Sample_ID],1,0)

                  res_int_signal <- res_int_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_mz_signal <- res_mz_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_rt_signal <- res_rt_total[which(res[[p]]$ion_data$signal_background == 1)]
                  res_isotope_signal <- res_isotope_total[which(res[[p]]$ion_data$signal_background == 1)]
                  length_data_signal <- length(res_int_signal)

                  ###save signal Intensities
                  if(length(res_int_signal)>0)
                  {
                    set(peaks_quant[[p]]$Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_signal)),
                                                                                                     length(res_int_signal),
                                                                                                     mean(log10(res_int_signal),na.rm=T),
                                                                                                     sd(log10(res_int_signal),na.rm=T),
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
                    set(peaks_quant[[p]]$Intensities_signal_background,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                                                       length(res_int_total),
                                                                                                                       mean(log10(res_int_total),na.rm=T),
                                                                                                                       sd(log10(res_int_total),na.rm=T),
                                                                                                                       m.z_window_final[1],
                                                                                                                       m.z_window_final[2],
                                                                                                                       RT_window_final[1],
                                                                                                                       RT_window_final[2],
                                                                                                                       num_peaks,
                                                                                                                       res[[p]]$Peak_info$known_peak)))
                  }else ###no intensity at all then at least save standard information
                  {
                    set(peaks_quant[[p]]$Intensities_signal_background,i=as.integer(1:10),j=as.integer(i),value=list(c(NA,
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
                    set(Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(log10(sum(res_int_total)),
                                                                                    length(res_int_total),
                                                                                    mean(log10(res_int_total),na.rm=T),
                                                                                    sd(log10(res_int_total),na.rm=T),
                                                                                    m.z_window_final[1],
                                                                                    m.z_window_final[2],
                                                                                    RT_window_final[1],
                                                                                    RT_window_final[2],
                                                                                    num_peaks,
                                                                                    res[[p]]$Peak_info$known_peak)))
                  }else ###no intensity at all then at least save standard information
                  {
                    set(Intensities,i=as.integer(1:10),j=as.integer(i),value=list(c(NA,
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
                  df <- data.frame(mz=res_mz_total,rt=res_rt_total,iso=res_isotope_total,int=res_int_total)
                  df_mono_isotope <- df[which(df$iso==0),]
                  df_isotope_1 <- df[which(df$iso==1),]
                  if(nrow(df_mono_isotope)>0)
                  {
                    ###deviation from consensus feature
                    max_mono <- which(df_mono_isotope$int == max(df_mono_isotope$int,na.rm=T))
                    set(peaks_quant[[p]]$delta_T1,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - (features_select$RT[i]))
                    set(peaks_quant[[p]]$delta_M1,i=1L,j=as.integer(i),mean(df_mono_isotope$mz[max_mono],na.rm=T) - (features_select$m.z[i]))
                    if(nrow(df_isotope_1)>0)
                    {
                      ###deviation from monoisotopic ion to M+1 isotope ion
                      max_iso_1 <- which(df_isotope_1$int == max(df_isotope_1$int,na.rm=T))
                      set(peaks_quant[[p]]$delta_T2,i=1L,j=as.integer(i),mean(df_mono_isotope$rt[max_mono],na.rm=T) - mean(df_isotope_1$rt[max_iso_1],na.rm=T))
                      set(peaks_quant[[p]]$delta_M2,i=1L,j=as.integer(i),(mean(df_mono_isotope$mz[max_mono],na.rm=T)+(iso_mult[1]/features_select$Charge[i])) - mean(df_isotope_1$mz[max_iso_1],na.rm=T))
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
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
    extract_intensities <- function(Sample_ID,features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection,ion_intensity_cutoff,mean_background_ion_intensity_model,sd_background_ion_intensity,peak_min_ion_count,kde_resolution,num_peaks_store,plots)
    {
      suppressWarnings(suppressMessages(library(data.table,quietly = T)))

      res <- get_intensities(Sample_ID,path = path_to_raw,features_select=features_select,RT_calibration=RT_calibration,mz_calibration=mz_calibration,peak_detection=peak_detection,ion_intensity_cutoff = ion_intensity_cutoff,mean_background_ion_intensity_model=mean_background_ion_intensity_model,sd_background_ion_intensity=sd_background_ion_intensity,peak_min_ion_count=peak_min_ion_count,kde_resolution=kde_resolution,num_peaks_store = num_peaks_store,plots=plots)

      peaks_quant <- res$peaks_quant
      #peaks_graph <- res$graph_peaks

      save(peaks_quant,file = paste(path_to_output_folder,"\\",Sample_ID,"_feature_quant.RData",sep=""))
      #if(length(peaks_graph)>0)save(peaks_graph,file = paste(path_to_output_folder,"\\",Sample_ID,"_feature_graphs.RData",sep=""))

    }

    ####Prepare threads to run extraction. The task is using much memory so that more than 2 threads in parallel on a pc with 16 gb of ram
    ####results in slower performance than for just 2 threads
    cl <- makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    registerDoParallel(cl)
    res <- foreach(i=Sample_IDs) %dopar%
      {
        extract_intensities(Sample_ID = i,features_select = features_select,path_to_raw,path_to_output_folder,RT_calibration,mz_calibration,peak_detection=peak_detection,ion_intensity_cutoff = ion_intensity_cutoff,mean_background_ion_intensity_model=mean_background_ion_intensity_model,sd_background_ion_intensity=sd_background_ion_intensity,peak_min_ion_count=peak_min_ion_count,kde_resolution=kde_resolution,num_peaks_store=num_peaks_store,plots=plots)
      }
    stopCluster(cl)

  }

  ###check for which samples mzXML files are available
  mzXMLfiles <- list.files(path_to_mzXML)
  mzXMLfiles <- mzXMLfiles[which(grepl(".mzXML",mzXMLfiles))]
  samples <- mzXMLfiles
  samples <- substr(samples,1,regexpr(".mzXML",samples)-1)

  ###keep samples which should be actually requantified
  samples <- samples[which(samples %in% gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration\\.",colnames(features)))]))]

  ###Perform quantification of decoy features

  dir.create(paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,sep=""),showWarnings = F)
  available <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,sep=""))
  available <- gsub("_feature_quant.RData","",available)
  if(length(which(samples %not in% available))>0)
  {
    samples <- samples[which(samples %not in% available)]
    extract_intensities_worker(Sample_IDs = as.character(samples),
                               features_select = features[which(grepl("_d",features$Feature_name)),],
                               path_to_raw = paste(path_to_mzXML,"\\all_ion_lists",sep=""),
                               path_to_output_folder = paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,sep=""),
                               RT_calibration=RT_calibration,
                               mz_calibration=mz_calibration,
                               peak_detection=F,
                               n_cores=n_cores)
  }

  ###determine distribution of background ion intensities
  samples <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,sep=""))
  samples <- samples[which(grepl("_feature_quant.RData",samples))]
  samples <- substr(samples,1,regexpr("_feature_quant.RData",samples)-1)

  files <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,sep=""))
  if(length(which(!grepl("_feature_quant.RData",files)))>0)
  {
    files <- files[-which(!grepl("_feature_quant.RData",files))]
  }

  features_select = features[which(grepl("_d",features$Feature_name)),]

  decoy_intensities <- as.data.frame(matrix(ncol=length(samples),nrow=nrow(features_select)))
  colnames(decoy_intensities) <- samples
  decoy_intensities <- sapply( decoy_intensities, as.numeric )
  decoy_intensities <- as.data.table(decoy_intensities)
  rownames(decoy_intensities) <- features_select$Feature_name

  decoy_ioncount <- decoy_intensities
  decoy_mean_intensity <- decoy_intensities
  decoy_sd_intensity <- decoy_intensities

  for(c in 1:ncol(decoy_intensities))
  {
    #load stored data into variable peaks_quant
    load(paste(path_to_mzXML,"\\all_ion_lists\\Extracted decoy intensities",output_file_names_add,"\\",colnames(decoy_intensities)[c],"_feature_quant.RData",sep=""))

    signal=peaks_quant[[1]]$Intensities

    decoy_intensities[,c] <- log2(10^signal[,1])
    decoy_ioncount[,c] <- signal[,2]
    decoy_mean_intensity[,c] <- log2(10^signal[,3])
    decoy_sd_intensity[,c] <- log2(10^signal[,4])
  }

  ###plot general numbers of quantifications of decoy features
  setwd(path_to_features)
  pdf("Temporary_files\\Decoy feature quantification parameters.pdf")

  ###mean decoy intensity (intensity of a single decoy ion)
  RT_all <- rep(as.numeric(features_select$RT),ncol(decoy_mean_intensity))
  x_all <- as.numeric(as.matrix(decoy_mean_intensity))
  smoothScatter(RT_all,x_all,ylab="Intensity, log2",main="All samples - Decoy feature mean intensity",xlab="RT [min]")

  ##try to fit an average generalised additive model to determine a RT dependent mean intensity and sd of intensity
  fit_gam_mean <- gam(x_all ~ s(RT_all), method = "REML")
  x_pred <- seq(min(features_select$RT,na.rm=T), max(features_select$RT,na.rm=T), length.out = nrow(features_select))
  y_pred <- predict(fit_gam_mean, data.frame(RT_all = x_pred))
  lines(x_pred,y_pred,col="red")
  legend("topright",legend="GAM",lty=1,col="red")

  ##now fit gam models per sample
  fit_gam_per_sample <- list()
  for(c in 1:ncol(decoy_mean_intensity))
  {
    RT <- as.numeric(features_select$RT)
    x <- as.numeric(as.matrix(decoy_mean_intensity)[,c])
    smoothScatter(RT,x,ylab="Mean intensity, log2",main=paste(colnames(decoy_mean_intensity)[c],"Decoy feature intensity"),xlab="RT [min]")

    ##try to fit an average generalised additive model to determine a RT dependent mean intensity and sd of intensity
    gam <- gam(x ~ s(RT), method = "REML")
    x_pred <- seq(min(features_select$RT,na.rm=T), max(features_select$RT,na.rm=T), length.out = nrow(features_select))
    y_pred <- predict(gam, data.frame(RT = x_pred))
    lines(x_pred,y_pred,col="red")
    legend("topright",legend="GAM",lty=1,col="red")
    fit_gam_per_sample[[colnames(decoy_mean_intensity)[c]]] <- gam
  }

  par(mfrow=c(2,2))
  boxplot(as.numeric(as.matrix(decoy_intensities)),outline=F,ylab="Summed intensity, log2",main="Summed intensity of decoy ions")
  boxplot(as.numeric(as.matrix(decoy_mean_intensity)),outline=F,ylab="Mean intensity, log2",main="Mean intensity of decoy ions")
  boxplot(as.numeric(as.matrix(decoy_sd_intensity)),outline=F,ylab="SD of intensity, log2",main="SD of intensity of decoy ions")
  boxplot(as.numeric(as.matrix(decoy_ioncount)),outline=F,ylab="Count",main="Number of ions in decoys with quantification")
  par(mfrow=c(1,1))
  dev.off()

  ###Now predict background intensities per sample and feature by using individually fitted GAMs, multiply with median decoy ion count and add some noise using observed decoy intensity sds per sample
  background_intensity_GAM_table_per_feature <- NULL ###used for determining which ions are background ions and which are signal ions
  background_intensity_GAM_table_per_feature_sum <- NULL ###used to later impute missing values

  median_decoy_sd_per_sample <- as.data.frame(t(colMedians(as.matrix(decoy_sd_intensity),na.rm=T)))
  colnames(median_decoy_sd_per_sample) <- colnames(decoy_sd_intensity)

  median_decoy_ion_count_per_sample <- as.data.frame(t(colMedians(as.matrix(decoy_ioncount),na.rm=T)))
  colnames(median_decoy_ion_count_per_sample) <- colnames(decoy_ioncount)

  max <- ncol(decoy_mean_intensity)
  pb <- winProgressBar(title = "Model background intensities per feature",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
  start_time <- Sys.time()
  updatecounter <- 0
  time_require <- 0

  for(c in 1:ncol(decoy_mean_intensity))
  {
    cur_sample <- colnames(decoy_mean_intensity)[c]
    if(c != 1)
    {
      background_intensity_GAM_table_per_feature <- cbind(background_intensity_GAM_table_per_feature,
                                                          data.frame(mean_intensity=predict(fit_gam_per_sample[[cur_sample]], data.frame(RT = features$RT))))
      background_intensity_GAM_table_per_feature_sum <- cbind(background_intensity_GAM_table_per_feature_sum,
                                                              data.frame(mean_intensity=predict(fit_gam_per_sample[[cur_sample]], data.frame(RT = features$RT))))
    }else
    {
      background_intensity_GAM_table_per_feature <- data.frame(mean_intensity=predict(fit_gam_per_sample[[cur_sample]], data.frame(RT = features$RT)))
      background_intensity_GAM_table_per_feature_sum <- data.frame(mean_intensity=predict(fit_gam_per_sample[[cur_sample]], data.frame(RT = features$RT)))
    }
    colnames(background_intensity_GAM_table_per_feature)[ncol(background_intensity_GAM_table_per_feature)] <- cur_sample
    colnames(background_intensity_GAM_table_per_feature_sum)[ncol(background_intensity_GAM_table_per_feature_sum)] <- cur_sample

    background_intensity_GAM_table_per_feature_sum[,c] <- rnorm(n=nrow(background_intensity_GAM_table_per_feature_sum),mean = background_intensity_GAM_table_per_feature_sum[,c],sd = median_decoy_sd_per_sample[1,c])

    updatecounter <- updatecounter + 1
    if(updatecounter >= 1)
    {
      time_elapsed <- difftime(Sys.time(),start_time,units="secs")
      time_require <- (time_elapsed/(c/max))*(1-(c/max))
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, c, label=paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)

  sd_background_intensity <- median_decoy_sd_per_sample

  rownames(background_intensity_GAM_table_per_feature) <- features$Feature_name
  rownames(background_intensity_GAM_table_per_feature_sum) <- features$Feature_name

  QC_data[["Decoy_feature_parameters"]] <- list(Decoy_mean_intensity_per_RT = data.frame(RT_all=RT_all,
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
    isotope_features$Feature_name <- paste(isotope_features$Feature_name,"_i",sep="")
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
    isotope_features$Feature_name <- paste(isotope_features$Feature_name,"_i",sep="")
    decoy_features_keep <- rbind(decoy_features_keep,isotope_features)
  }
  ###remove all decoy features
  features <- features[-decoy_features,]
  ###add selected decoy features plus isotope decoy features
  features <- rbind(features,decoy_features_keep)

  ###see if already samples were converted
  samples <- mzXMLfiles
  samples <- substr(samples,1,regexpr(".mzXML",samples)-1)
  samples <- samples[which(samples %in% gsub("RT_calibration\\.","",colnames(features)[which(grepl("RT_calibration\\.",colnames(features)))]))]

  dir.create(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""),showWarnings = F)
  available <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""))
  available <- available[which(grepl("feature_quant.RData",available))]
  available <- gsub("_feature_quant.RData","",available)
  peak_min_ion_count <- boxplot.stats(as.numeric(as.matrix(decoy_ioncount)))$stats[4]
  if(length(which(samples %not in% available))>0)
  {
    samples <- samples[which(samples %not in% available)]
    ##use decoy defined cut of to distinguish background intensity from signal intensity per feature
    extract_intensities_worker(Sample_IDs = as.character(samples),
                               features_select = features,
                               path_to_raw = paste(path_to_mzXML,"\\all_ion_lists",sep=""),
                               path_to_output_folder = paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""),
                               RT_calibration=RT_calibration,
                               mz_calibration=mz_calibration,
                               peak_detection=peak_detection,
                               n_cores=n_cores,
                               ion_intensity_cutoff = T,
                               mean_background_ion_intensity_model = background_intensity_GAM_table_per_feature,
                               sd_background_ion_intensity = sd_background_intensity,
                               peak_min_ion_count=peak_min_ion_count,
                               kde_resolution = kde_resolution,
                               num_peaks_store = num_peaks_store,
                               plots = plot_2D_peak_detection) ###define background peaks during 2DKDE as peaks with <= 75% quantile
  }

  #####Summarize data for all features and samples
  samples <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""))
  samples <- samples[which(grepl("feature_quant.RData",samples))]
  samples <- substr(samples,1,regexpr("_feature_quant.RData",samples)-1)

  ###Available sample data
  files <- list.files(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""))
  if(length(which(!grepl("feature_quant.RData",files)))>0)
  {
    files <- files[-which(!grepl("feature_quant.RData",files))]
  }

  features_select <- features

  features_intensity <- as.data.frame(matrix(ncol=length(samples),nrow=nrow(features_select)))#
  colnames(features_intensity) <- samples
  features_intensity <- sapply( features_intensity, as.numeric )
  features_intensity <- as.data.table(features_intensity)
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
    names(peak_quant) <- c(paste("Peak_",1:(num_peaks_store),sep=""),"Standard")
  }

  max <- ncol(features_intensity)
  pb <- winProgressBar(title = "Merge quantification results",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
      load(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,"\\",colnames(features_intensity)[c],"_feature_quant.RData",sep=""))

      signal=peaks_quant[[1]]$Intensities
      signal_background=peaks_quant[[1]]$Intensities_signal_background

      features_intensity[,c] <- signal[,1]
      Ioncount_feature_sample_matrix[,c] <- signal[,2]

      feature_with_background_intensity[,c] <- signal_background[,1]
      Ioncount_feature_with_background_intensity[,c] <- signal_background[,2]

    }else
    {
      #load stored data into variable peaks_quant
      load(paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,"\\",colnames(features_intensity)[c],"_feature_quant.RData",sep=""))

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
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, c, label=paste( round(c/max*100, 0)," % done (",c,"/",max,", Time require: ",time_require,")",sep = ""))
    }
  }
  close(pb)

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

    features_intensity <- as.data.frame(features_intensity)
    Ioncount_feature_sample_matrix <- as.data.frame(Ioncount_feature_sample_matrix)
    feature_with_background_intensity <- as.data.frame(feature_with_background_intensity)
    Ioncount_feature_with_background_intensity <- as.data.frame(Ioncount_feature_with_background_intensity)
    peak_selected <- as.data.frame(peak_selected)

    ###For the next step of peak quantifications where true peak is not known we will need RT and mz correction factors per feature
    ###Extract RT and mz correction factors per sample and feature
    indices_RT_correction <- which(grepl("RT_calibration",colnames(features_select)))
    indices_mz_correction <- which(grepl("mz_calibration",colnames(features_select)))
    ordering_indices <- match(samples,gsub("RT_calibration\\.","",colnames(features_select)[indices_RT_correction]))
    indices_RT_correction <- indices_RT_correction[ordering_indices]
    indices_mz_correction <- indices_mz_correction[ordering_indices]

    RT_correction_factors <- features_select[,indices_RT_correction]
    mz_correction_factors <- features_select[,indices_mz_correction]
    colnames(RT_correction_factors) <- samples
    rownames(RT_correction_factors) <- features_select$Feature_name
    colnames(mz_correction_factors) <- samples
    rownames(mz_correction_factors) <- features_select$Feature_name

    ###general RT and mz windows
    delta_mz <- median(features_select$m.z - features_select$m.z_range_min,na.rm=T)
    delta_rt <- median(features_select$RT - features_select$RT_range_min,na.rm=T)/2

    ###add step to store information which peak was finally picked
    peak_decision <- function(features_select,peak_quant,samples,s,RT_correction_factors,mz_correction_factors,features_intensity_sample,Ioncount_sample,feature_with_background_intensity_sample,Ioncount_with_background_sample,peak_selected_sample,delta_mz,delta_rt,peak_min_ion_count,chunk=NULL,num_chunks=NULL,progress=T)
    {
      suppressWarnings(suppressMessages(library(data.table,quietly = T)))
      suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
      if(progress == T)
      {
        pb <- winProgressBar(title = paste("Prepare for peak selection -",samples[s]),label=paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
        close(pb)
      }

      ###next, decide for all other feature quantifications for which at least 1 peak was available which peak quantification should be used
      max <- nrow(features_select)

      if(progress == T)
      {
        if(!is.null(chunk) & !is.null(num_chunks))
        {
          pb <- winProgressBar(title = paste("Select peaks (",chunk,"/",num_chunks,") - ",samples[s],sep=""),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
        }else
        {
          pb <- winProgressBar(title = paste("Select peaks -",samples[s]),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
            #   # set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$features_intensity[i,..s]))
            #   # set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$Ioncount_feature_sample_matrix[i,..s]))
            #   # set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$feature_with_background_intensity[i,..s]))
            #   # set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$Ioncount_feature_with_background_intensity[i,..s]))
            #   #
            #   # set(dT1_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$dT1[i,..s]))
            #   # set(dM1_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Peak_1$dM1[i,..s]))
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
                  delta_RT_cut <- 3*sd(RT_in_known_corrected,na.rm=T)
                  if(is.na(delta_RT_cut) | delta_RT_cut < 2*delta_rt)delta_RT_cut <- 2*delta_rt
                  delta_mz_cut <- 3*sd(mz_in_known_corrected,na.rm=T)
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
                      other_peaks <- as.data.frame(matrix(ncol=6,nrow=count_other_peaks,0))
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

                            set(other_peaks,as.integer(counter),as.integer(1:6),value=list(temp[1],temp[2],temp[3],temp[4],temp[5],temp[6]))
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
                      set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$features_intensity[i,..s]))
                      set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_sample_matrix[i,..s]))

                      set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$feature_with_background_intensity[i,..s]))
                      set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_with_background_intensity[i,..s]))

                      set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(selected_peak))

                    }else###none of the peaks is not overlapping with an other peak detected in other samples where correct peak was known
                    {
                      set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                      set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                      set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                      set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                      set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                    }

                  }else ###none of the peaks is within the accepted range thus use standard window
                  {
                    set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                    set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                    set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                    set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                    set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                  }
                }else ###none of the peaks is valid
                {
                  set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                  set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                  set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                  set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                  set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                }


              }else ###in no other sample the correct peak is known thus take total window
              {
                set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
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
                  delta_RT_cut <- 3*sd(RT_in_known,na.rm=T)
                  if(is.na(delta_RT_cut) | delta_RT_cut < 2*delta_rt)delta_RT_cut <- 2*delta_rt
                  delta_mz_cut <- 3*sd(mz_in_known,na.rm=T)
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
                      other_peaks <- as.data.frame(matrix(ncol=6,nrow=count_other_peaks,0))
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

                      set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$features_intensity[i,..s]))
                      set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_sample_matrix[i,..s]))

                      set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$feature_with_background_intensity[i,..s]))
                      set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant[[selected_peak]]$Ioncount_feature_with_background_intensity[i,..s]))

                      set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(selected_peak))
                    }else###the peak is overlapping with an other peak detected in other samples where correct peak was known
                    {
                      set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                      set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                      set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                      set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                      set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                    }
                  }else ###none of the peaks is within the accepted range thus use standard window
                  {
                    set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                    set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                    set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                    set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                    set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                  }
                }else ###no peak was valid
                {
                  set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                  set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                  set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                  set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                  set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
                }
              }else ###in no other sample the correct peak was known thus use the standard window
              {
                set(features_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$features_intensity[i,..s]))
                set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_sample_matrix[i,..s]))

                set(feature_with_background_intensity_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$feature_with_background_intensity[i,..s]))
                set(Ioncount_with_background_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Ioncount_feature_with_background_intensity[i,..s]))

                set(peak_selected_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(length(peak_quant)))
              }


            }


          }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10 & progress == T)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
      suppressWarnings(suppressMessages(library(doParallel,quietly = T)))

      cl <- makeCluster(n_cores)
      registerDoParallel(cl)

      res <- foreach(s=1:length(samples)) %dopar%
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
      stopCluster(cl)

      ###merge results from all threads
      for(c in 1:length(samples))
      {
        set(features_intensity,as.integer(start:end),c,res[[c]]$features_intensity_sample)
        set(Ioncount_feature_sample_matrix,as.integer(start:end),c,res[[c]]$Ioncount_sample)
        set(feature_with_background_intensity,as.integer(start:end),c,res[[c]]$feature_with_background_intensity_sample)
        set(Ioncount_feature_with_background_intensity,as.integer(start:end),c,res[[c]]$Ioncount_with_background_sample)
        set(peak_selected,as.integer(start:end),c,res[[c]]$peak_selected_sample)
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
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  print("Save peak detection and selection results")
  save(features,
       features_intensity,
       Ioncount_feature_sample_matrix,
       feature_with_background_intensity,
       Ioncount_feature_with_background_intensity,
       peak_selected,
       peak_quant,
       QC_data,
       file = "Quantification_raw_results.RData")

  crap <- gc(F)

  ###read previously generated outputs
  #load("Quantification_raw_results.RData")

  setwd(path_to_features)

  pdf("Temporary_files\\Alignment and quantification scores.pdf")

  features$target_decoy <- ifelse(grepl("_d",features$Feature_name),"decoy","target")
  ###Tag decoy features which overlap with real features
  remove_decoy_outlier <- vector(mode="logical",length=length(which(features$target_decoy == "decoy")))
  features_target <- features[which(features$target_decoy == "target"),]
  decoy_indices <- which(features$target_decoy == "decoy" & !grepl("_d_i",features$Feature_name))

  max <- length(decoy_indices)
  pb <- winProgressBar(title = "Detect target-decoy overlapping features",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
      td <- seconds_to_period(time_require)
      time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

      updatecounter <- 0
      setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
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
  mean_decoy_count <- mean(log2(decoy_ion_count))
  sd_decoy_count <- sd(log2(decoy_ion_count))

  ###Prepare for testing if observed ion counts of target features is significantly higher than for decoy ions (just random peak selection)
  target_ion_counts <- Ioncount_feature_with_background_intensity
  rows <- nrow(target_ion_counts)
  rownames(target_ion_counts) <- rownames(Ioncount_feature_with_background_intensity)
  target_ion_counts <- as.matrix(target_ion_counts)
  ###add 1 for log scale
  target_ion_counts <- target_ion_counts + 1
  zscores <- as.data.frame((log2(target_ion_counts)-mean_decoy_count)/sd_decoy_count)

  zscore_to_pval <- function(z, alternative="greater")
  {

    if(alternative == "greater")
    {
      pval = pnorm(-z)
    }else if(alternative == "two.sided")
    {
      pval = pnorm(-abs(z))
      pval = 2 * pval
    }else if(alternative=="less")
    {
      pval = pnorm(z)
    }
    return(pval)
  }

  ###perform significance test
  pval_signal_with_background_quant <- as.data.frame(apply(zscores,1:2,zscore_to_pval))
  temp_pval_quant <- pval_signal_with_background_quant
  temp_pval_quant[is.na(temp_pval_quant)] <- 1

  QC_data[["Decoy_ion_counts"]] <- decoy_ion_count
  QC_data[["Target_ion_counts"]] <- target_ion_counts
  QC_data[["Quant_pval"]] <- pval_signal_with_background_quant

  ###plot pvalue of quantification per sample for target features
  if(ncol(temp_pval_quant) <= 15) ###only 15 samples can be plotted in one plot
  {
    par(mar=c(10,4,4,4))
    boxplot(-log10(temp_pval_quant[which(!grepl("_d",rownames(temp_pval_quant))),]),outline=F,main="Quantification pVals per target feature",ylab="pValue, -log10",names=colnames(temp_pval_quant),las=2)
    abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
  }else
  {
    par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(temp_pval_quant)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(temp_pval_quant))]
      boxplot(-log10(temp_pval_quant[which(!grepl("_d",rownames(temp_pval_quant))),columns]),outline=F,main="Quantification pVals per target feature",ylab="pValue, -log10",names=colnames(temp_pval_quant)[columns],las=2)
      abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
    }
  }

  ###plot pvalue of quantification per sample for decoy features
  if(ncol(temp_pval_quant) <= 15) ###only 15 samples can be plotted in one plot
  {
    par(mar=c(10,4,4,4))
    boxplot(-log10(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),]),outline=F,main="Quantification pVals per decoy feature",ylab="pValue, -log10",names=colnames(temp_pval_quant),las=2)
    abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
  }else
  {
    par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(temp_pval_quant)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(temp_pval_quant))]
      boxplot(-log10(temp_pval_quant[which(grepl("_d$",rownames(temp_pval_quant))),columns]),outline=F,main="Quantification pVals per decoy feature",ylab="pValue, -log10",names=colnames(temp_pval_quant)[columns],las=2)
      abline(h=-log10(Quant_pVal_cut),lty=2,col="red")
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
  S2B <- log2(10^features_intensity/10^features_background_intensity)
  S2B[is.infinite(as.matrix(S2B))] <- NA
  S2B <- as.data.frame(S2B)
  rownames(S2B) <- features$Feature_name

  QC_data[["Signal_to_background_target_decoy"]] <- list(Target_Decoy=ifelse(grepl("_d",rownames(S2B)),"decoy","target"),
                                                         S2B=S2B)

  if(ncol(S2B) <= 15) ###only 15 samples can be plotted in one plot
  {
    par(mar=c(10,4,4,4))
    boxplot(S2B[which(!grepl("_d",rownames(S2B))),],outline=F,main="Signal to background ratio per feature quantification",ylab="Signal/Background, log2",names=colnames(S2B),las=2)
    abline(h=0,lty=2,col="red")
    for(co in 1:ncol(S2B))
    {
      med <- median(S2B[which(!grepl("_d",rownames(S2B))),co],na.rm=T)
      text(co,med+((par("usr")[4]-par("usr")[3])*0.05),round(med,digits=1))
    }

  }else
  {
    par(mar=c(10,4,4,4))
    pages <- ceiling(ncol(S2B)/15)
    for(p in 1:pages)
    {
      columns <- (((p-1)*15)+1):(p*15)
      columns <- columns[which(columns <= ncol(S2B))]
      boxplot(S2B[which(!grepl("_d",rownames(S2B))),columns],outline=F,main="Signal to background ratio per feature quantification",ylab="Signal/Background, log2",names=colnames(S2B)[columns],las=2)
      abline(h=0,lty=2,col="red")
      for(co in 1:length(columns))
      {
        med <- median(S2B[which(!grepl("_d",rownames(S2B))),columns[co]],na.rm=T)
        text(co,med+((par("usr")[4]-par("usr")[3])*0.05),round(med,digits=1))
      }
    }
  }

  pval_signal_with_background_quant <- as.data.frame(pval_signal_with_background_quant)
  features_intensity <- as.matrix(features_intensity)
  feature_with_background_intensity <- as.matrix(feature_with_background_intensity)
  features_intensity[is.infinite(features_intensity)] <- NA
  feature_with_background_intensity[is.infinite(feature_with_background_intensity)] <- NA
  features_intensity <- as.data.frame(features_intensity)
  feature_with_background_intensity <- as.data.frame(feature_with_background_intensity)

  ###calculate some scores which give insight in confidence of quantifications (selection of right peaks)

  #arrange selected peaks results in a long table

  peaks <- data.frame(RT=unlist(sapply(peak_quant, function(x) x$Peak_rt_with_background)),
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
    corrections <- features[,c(1,which(grepl(paste("_calibration\\.",samples[s],"$",sep=""),colnames(features))))]
    selection_peaks <- which(peaks$sample == samples[s])
    RT_corrected <- peaks$RT[selection_peaks]-corrections[match(peaks$feature[selection_peaks],corrections$Feature_name),2]
    mz_corrected <- peaks$mz[selection_peaks]-corrections[match(peaks$feature[selection_peaks],corrections$Feature_name),3]
    set(peaks,as.integer(selection_peaks),c(10L,11L),value=list(RT_corrected,mz_corrected))
  }


  #calculate general peak variability score using delta_rt and delta_mz as sd for zscoring
  ###feature with significant score here should be potentially excluded
  Variability_alignment_scoring <- function(features,peaks,delta_rt,delta_mz,alignment_variability_score_cutoff=0.05,plot=T)
  {
    ###Convert z scores (vector or matrix) to pvalues
    zscore_to_pval <- function(z, alternative="greater")
    {

      if(alternative == "greater")
      {
        pval = pnorm(-z)
      }else if(alternative == "two.sided")
      {
        pval = pnorm(-abs(z))
        pval = 2 * pval
      }else if(alternative=="less")
      {
        pval = pnorm(z)
      }
      return(pval)
    }
    #select peaks which were selected by peak selection
    temp <- peaks[which(peaks$peak == peaks$selected),]
    #calculate interquartile range per feature over samples for RT and mz
    temp_RT <- aggregate(temp$RT,by=list(Feature=temp$feature),FUN=iqr,na.rm=T)
    temp_mz <- aggregate(temp$mz,by=list(Feature=temp$feature),FUN=iqr,na.rm=T)

    #z score using delta RT and delta mz as sd
    temp_RT$x <- (temp_RT$x-mean(temp_RT$x,na.rm=T))/delta_rt
    temp_mz$x <- (temp_mz$x-mean(temp_mz$x,na.rm=T))/delta_mz
    ##features with a negative z-score represent features with much lower deviation between samples than expected
    ##features with a positive z-score represent features with much higher deviation between samples than expected --> have to be inspected and eventually removed
    #boxplot(temp_mz$x,main="General mz-deviation of selected peaks between samples",ylab="Z-score - mz deviation",ylim=c(-4,4))
    #boxplot(temp_RT$x,main="General RT-deviation of selected peaks between samples",ylab="Z-score - RT deviation",ylim=c(-4,4))

    alignment_variability_score <- data.frame(RT_variability_pval=(zscore_to_pval(temp_RT$x)),
                                              mz_variability_pval=(zscore_to_pval(temp_mz$x)))

    alignment_variability_score <- alignment_variability_score[match(features$Feature_name,temp_RT$Feature),]
    rownames(alignment_variability_score) <- features$Feature_name
    alignment_variability_score <- as.data.frame(alignment_variability_score)
    alignment_variability_score <- cbind(alignment_variability_score,
                                         rowMins(as.matrix(alignment_variability_score),na.rm=T))
    colnames(alignment_variability_score)[3] <- "combined_variability_pval"

    #add some plots showing number of features with general high variability between samples
    if(plot == T)
    {
      temp <- rowMins(as.matrix(alignment_variability_score[which(!grepl("_d|_i",rownames(alignment_variability_score))),]),na.rm=T)
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
    #prepare matrices which should store RT and mz of selected peaks per sample and feature
    if(corrected_alignment == T)selected_peaks_RT <- reshape(temp_peaks[,c("feature","sample","RT_correct")], idvar = "feature", timevar = "sample", direction = "wide")
    if(corrected_alignment == F)selected_peaks_RT <- reshape(temp_peaks[,c("feature","sample","RT")], idvar = "feature", timevar = "sample", direction = "wide")
    rownames(selected_peaks_RT) <- selected_peaks_RT$feature
    selected_peaks_RT <- selected_peaks_RT[,-1]
    if(corrected_alignment == T)selected_peaks_mz <- reshape(temp_peaks[,c("feature","sample","mz_correct")], idvar = "feature", timevar = "sample", direction = "wide")
    if(corrected_alignment == F)selected_peaks_mz <- reshape(temp_peaks[,c("feature","sample","mz")], idvar = "feature", timevar = "sample", direction = "wide")
    rownames(selected_peaks_mz) <- selected_peaks_mz$feature
    selected_peaks_mz <- selected_peaks_mz[,-1]
    #calculate mean and sd of RT and mz per feature
    selected_peaks_distributions_RT <- data.frame(mean=rowMeans(as.matrix(selected_peaks_RT),na.rm=T),
                                                  sd=rowSds(as.matrix(selected_peaks_RT,na.rm=T)))
    selected_peaks_distributions_mz <- data.frame(mean=rowMeans(as.matrix(selected_peaks_mz),na.rm=T),
                                                  sd=rowSds(as.matrix(selected_peaks_mz,na.rm=T)))
    #convert observed RTs and mzs per feature and sample into z-scores
    selected_peaks_RT_zscore <- (selected_peaks_RT-selected_peaks_distributions_RT$mean)/sd_RT
    selected_peaks_mz_zscore <- (selected_peaks_mz-selected_peaks_distributions_mz$mean)/sd_mz
    #convert zscores into pvalues
    zscore_to_pval <- function(z)
    {
      pval = pnorm(-abs(z))
      pval = 2 * pval
      return(pval)
    }
    selected_peaks_RT_pval <- as.data.frame(apply(selected_peaks_RT_zscore,1:2,zscore_to_pval))
    selected_peaks_mz_pval <- as.data.frame(apply(selected_peaks_mz_zscore,1:2,zscore_to_pval))
    if(corrected_alignment == T)colnames(selected_peaks_RT_pval) <- gsub("RT_correct.","",colnames(selected_peaks_RT_pval))
    if(corrected_alignment == T)colnames(selected_peaks_mz_pval) <- gsub("mz_correct.","",colnames(selected_peaks_mz_pval))
    if(corrected_alignment == F)colnames(selected_peaks_RT_pval) <- gsub("RT.","",colnames(selected_peaks_RT_pval))
    if(corrected_alignment == F)colnames(selected_peaks_mz_pval) <- gsub("mz.","",colnames(selected_peaks_mz_pval))
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
  Mono_iso_alignment_scoring <- function(features,samples,peaks,pval_signal_with_background_quant,delta_rt,delta_mz,mono_iso_alignment_cutoff,plot=T)
  {
    #prepare dataframe into which all results are saved
    mono_iso_alignment_summary <- as.data.frame(matrix(nrow=nrow(features),ncol=length(samples)))
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
          pval = pnorm(-z)
        }else if(alternative == "two.sided")
        {
          pval = pnorm(-abs(z))
          pval = 2 * pval
        }else if(alternative=="less")
        {
          pval = pnorm(z)
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
        temp_selected_iso$feature <- gsub("_i","",temp_selected_iso$feature)
        #combine both and only keep features for which we also have +1 isotope features
        combined <- inner_join(temp_selected_mono,temp_selected_iso,by="feature")
        #determine expected variation between mono and +1 isotopes based on true identifications and significantly quantified
        #mono and +1 isotope features
        true <- combined[which(combined$known.x == 1 & combined$quant_score.x < 0.05 & combined$quant_score.y < 0.05),]
        true$delta_rt_iso_mono <- true$RT.y-true$RT.x
        true$delta_mz_iso_mono <- (((true$mz.y*true$charge.y)-1.002054)/true$charge.y)-true$mz.x

        mean_RT_distribution[s] <- mean(true$delta_rt_iso_mono,na.rm=T)
        sd_RT_distribution[s] <- sd(true$delta_rt_iso_mono,na.rm=T)
        mean_mz_distribution[s] <- mean(true$delta_mz_iso_mono,na.rm=T)
        sd_mz_distribution[s] <- sd(true$delta_mz_iso_mono,na.rm=T)
      }
      #estimate normal distribution of deviations in RT and mz between mono and +1 isotope peak
      RT_deviation_norm <- data.frame(x=seq(-delta_rt*3, delta_rt*3, length=1000),
                                      y=dnorm(seq(-delta_rt*3, delta_rt*3, length=1000), mean=mean(mean_RT_distribution,na.rm=T), sd=mean(sd_RT_distribution,na.rm=T)))
      mz_deviation_norm <- data.frame(x=seq(-delta_mz*3,delta_mz*3, length=1000),
                                      y=dnorm(seq(-delta_mz*3,delta_mz*3, length=1000), mean=mean(mean_mz_distribution,na.rm=T), sd=mean(sd_mz_distribution,na.rm=T)))
      plot(RT_deviation_norm$x,RT_deviation_norm$y,type="l",xlab="RT(M+1) - RT(M)",main="RT-deviation M vs M+1 peaks - true quantifications",ylab="Density")
      abline(v=mean(true$delta_rt_iso_mono),lty=2)
      plot(mz_deviation_norm$x,mz_deviation_norm$y,type="l",xlab="mz(M+1) - mz(M)",main="mz-deviation M vs M+1 peaks - true quantifications",ylab="Density")
      abline(v=mean(true$delta_mz_iso_mono),lty=2)
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
        temp_selected_iso$feature <- gsub("_i","",temp_selected_iso$feature)
        #combine both and only keep features for which we also have +1 isotope features
        combined <- inner_join(temp_selected_mono,temp_selected_iso,by="feature")
        combined$delta_rt_iso_mono <- combined$RT.y-combined$RT.x
        combined$delta_mz_iso_mono <- (((combined$mz.y*combined$charge.y)-1.002054)/combined$charge.y)-combined$mz.x
        #convert to z-score
        combined$z_delta_rt_iso_mono <- (combined$delta_rt_iso_mono-mean_RT_distribution_mean)/sd_RT_distribution_mean
        combined$z_delta_mz_iso_mono <- (combined$delta_mz_iso_mono-mean_mz_distribution_mean)/sd_mz_distribution_mean
        #convert to pvalues
        combined$RT_deviation_pval <- zscore_to_pval(combined$z_delta_rt_iso_mono,"two.sided")
        combined$mz_deviation_pval <- zscore_to_pval(combined$z_delta_mz_iso_mono,"two.sided")
        #store results
        combined$score <- rowMins(as.matrix(combined[,c("RT_deviation_pval","mz_deviation_pval")]),na.rm=T)

        mono_iso_alignment_summary[,s] <- combined$score[match(gsub("_i","",rownames(mono_iso_alignment_summary)),combined$feature)]
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

  ###calculate FDR for peak selections
  ###for that purpose we select randomly 1000 quantifications in random samples where true peak was known
  ###we check if the same peak will be selected if information about true peak location is removed
  Peak_selection_FDR <- function(num_features=500,features,samples,peaks,path_to_features,peak_quant,feature_with_background_intensity,peak_selected,delta_mz,delta_rt,peak_min_ion_count,num_peaks_store,alignment_scores_cutoff,n_cores,peak_decision,Alignment_scoring,seed=1,plot=T)
  {
    if(length(which(peaks$known == 1)) >= 10)
    {
      suppressWarnings(suppressMessages(library(data.table,quietly = T)))
      suppressWarnings(suppressMessages(library(randomForest,quietly = T)))
      suppressWarnings(suppressMessages(library(mgcv,quietly = T)))
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

      setwd(paste(path_to_features,"\\Temporary_files",sep=""))
      e <- new.env()
      load("Feature_alignment_QC_data.RData",envir = e)
      median_feature_properties <- e$QC_data$mz_calibration_median_feature_properties
      mz_correction_models <- e$QC_data$mz_calibration_models
      RT_alignment_GAM_models <- e$QC_data$RT_calibration_GAM_models

      for(s in 1:length(samples))
      {
        temp <- data.table::copy(features)
        temp <- temp[selection_feature_per_sample[[s]],]

        #remove observed RT and mz
        obs_rt <- (str_split(temp$Observed_RT,";",simplify = T))
        obs_rt[,s] <- "NA"
        obs_mz <- (str_split(temp$Observed_mz,";",simplify = T))
        obs_mz[,s] <- "NA"
        temp$Observed_mz <- apply(obs_mz,1,paste,collapse = ";")
        #exchange mz calibration with prediction from model
        select_model <- which(names(mz_correction_models) == samples[s])

        features_select <- selection_feature_per_sample[[s]]

        temp_data <- median_feature_properties[match(features_select,median_feature_properties$Feature),]

        temp_data <- temp_data[which(rowSums(is.na(temp_data[,c("Retention.time","m.z","Charge","Resolution")])) == 0),]

        if(length(select_model)>0)
        {
          if(nrow(temp_data) > 0)
          {
            prediction <- stats::predict(mz_correction_models[[select_model]], temp_data[,-1], type = "response")
            if(any(is.na(prediction)))prediction[which(is.na(prediction))] <- 0

            set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(grepl("mz_calibration.",colnames(temp)))[s]),prediction)
          }else
          {
            set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(grepl("mz_calibration.",colnames(temp)))[s]),0)
          }
        }else
        {
          set(temp,as.integer(match(rownames(temp_data),rownames(temp))),as.integer(which(grepl("mz_calibration.",colnames(temp)))[s]),0)
        }
        #exchange RT calibration with prediction from model
        select_model <- s
        if(!is.null(RT_alignment_GAM_models[[select_model]]))
        {
          prediction <- stats::predict(RT_alignment_GAM_models[[select_model]], data.frame(x = temp$RT))
          if(any(is.na(prediction)))prediction[which(is.na(prediction))] <- 0
          set(temp,as.integer(1:nrow(temp)),as.integer(which(grepl("RT_calibration.",colnames(temp)))[s]),prediction)

        }else
        {
          set(temp,as.integer(1:nrow(temp)),as.integer(which(grepl("RT_calibration.",colnames(temp)))[s]),0)
        }

        features_select_FDR[[s]] <- temp
      }

      #prepare columns containing correction information
      indices_RT_correction <- which(grepl("RT_calibration",colnames(features)))
      indices_mz_correction <- which(grepl("mz_calibration",colnames(features)))
      ordering_indices <- match(samples,gsub("RT_calibration\\.","",colnames(features)[indices_RT_correction]))
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

      suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))

      cl <- makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
      registerDoSNOW(cl)

      start <- Sys.time()

      res <- foreach(s=1:length(samples)) %dopar%
        {
          suppressWarnings(suppressMessages(library(data.table,quietly = T)))
          max <- 1
          pb <- winProgressBar(title = "Evaluate peak selection FDR",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)

          selection_feature <- selection_feature_per_sample[[s]]
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

          temp_dummy_df <-  data.table(Name=rep(0L,length(selection_feature)))
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

          set(feature_with_background_intensity_select_FDR,as.integer(selection_feature),as.integer(s),res_temp$feature_with_background_intensity_sample)
          set(peak_selected_select_FDR,as.integer(selection_feature),as.integer(s),res_temp$peak_selected_sample)

          names(peak_decision_same_without_peak_known) <- paste(features_select_FDR_cur$Feature_name,"_Sample_",s,sep="")

          return(list(peak_decision_same_without_peak_known=peak_decision_same_without_peak_known,
                      feature_with_background_intensity_select_FDR=feature_with_background_intensity_select_FDR,
                      peak_selected_select_FDR=peak_selected_select_FDR))
        }
      stopCluster(cl)

      #determine how many false selections show large intensity difference and would not be removed by alignment scoring
      results_peak_selection_FDR_all <- list()
      max <- length(samples)
      pb <- winProgressBar(title="Evaluate peak selection FDR",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
      start_time <- Sys.time()
      updatecounter <- 0
      time_require <- 0

      for(s in 1:length(samples))
      {
        peak_decision_same_without_peak_known <- res[[s]]$peak_decision_same_without_peak_known
        peak_selected_select_FDR <- res[[s]]$peak_selected_select_FDR
        feature_with_background_intensity_select_FDR <- res[[s]]$feature_with_background_intensity_select_FDR
        selection_feature <- selection_feature_per_sample[[s]]
        if(length(which(peak_decision_same_without_peak_known == F))>0L)
        {
          #prepare updated peaks table
          # peaks_select_FDR <- data.table::copy(peaks)
          # peaks_select_FDR$feature <- as.character(peaks_select_FDR$feature)
          # selection_sample <- which(peaks_select_FDR$sample == samples[s] & peaks_select_FDR$feature %in% features$Feature_name[selection_feature])
          # for(i in selection_feature)
          # {
          #   f <- features$Feature_name[i]
          #   sel <- which(peaks_select_FDR$feature[selection_sample] == f)
          #   set(peaks_select_FDR,as.integer(selection_sample[sel]),9L,peak_selected_select_FDR[i,s])
          # }

          wrong_selections <- which(peak_decision_same_without_peak_known == F)
          #compare quantification differences for situation where wrong peak was selected
          compare_quant_results_for_wrong_decisions <- as.data.frame(matrix(ncol=2,nrow=length(wrong_selections)))
          colnames(compare_quant_results_for_wrong_decisions) <- c("correct","wrong")
          rownames(compare_quant_results_for_wrong_decisions) <- features_select_FDR[[s]]$Feature_name[which(peak_decision_same_without_peak_known == F)]
          compare_quant_results_for_wrong_decisions$correct <- as.numeric(compare_quant_results_for_wrong_decisions$correct)
          compare_quant_results_for_wrong_decisions$wrong <- as.numeric(compare_quant_results_for_wrong_decisions$wrong)
          for(i in wrong_selections)
          {
            ind <- which(wrong_selections == i)
            set(compare_quant_results_for_wrong_decisions,as.integer(ind),as.integer(1:2),list(log2(10^feature_with_background_intensity_select_FDR[selection_feature[i],s]),
                                                                                               log2(10^feature_with_background_intensity[selection_feature[i],s])))
          }
          ##compare intensities of correct and wrong peak selections
          #determine distribution of quantification deviations between correct and wrong
          if(plot==T & length(wrong_selections)>=3L)
          {
            plot(density(compare_quant_results_for_wrong_decisions$wrong-compare_quant_results_for_wrong_decisions$correct,na.rm=T),xlab="Intensity-wrong / Intensity-true, log2",main=paste(samples[s],"- Deviation in selected peak quantification"))
            abline(v=0)
            abline(v=-1,lty=2,col="red")
            abline(v=1,lty=2,col="red")

            plot(compare_quant_results_for_wrong_decisions$correct,compare_quant_results_for_wrong_decisions$wrong,xlab="True peak quantification, log2",ylab="Wrong peak quantification, log2",main=paste(samples[s],"- Error in quantification"))
            abline(a=0,b=1)
            abline(a=1,b=1,lty=2,col="red")
            abline(a=-1,b=1,lty=2,col="red")
          }

          #how many show deviation > 2 fold
          wrong_peak_intensity_outlier <- abs(compare_quant_results_for_wrong_decisions$wrong-compare_quant_results_for_wrong_decisions$correct) > 1

          #finally visualize results
          freq <- plyr::count(peak_decision_same_without_peak_known)
          if(length(which(freq$x == F)) == 0)
          {
            freq <- rbind(freq,data.frame(x=F,freq=0))
          }
          if(length(which(freq$x == T)) == 0)
          {
            freq <- rbind(freq,data.frame(x=T,freq=0))
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
            freq <- rbind(freq,data.frame(x=F,freq=0))
          }
          if(length(which(freq$x == T)) == 0)
          {
            freq <- rbind(freq,data.frame(x=T,freq=0))
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

        updatecounter <- updatecounter + 1
        if(updatecounter >= 1)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(s/max))*(1-(s/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, s, label=paste( round(s/max*100, 0)," % done (",s,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)

      ##plot results over all samples
      Total_FDR <- NULL
      Large_Intensity_delta_FDR <- NULL
      for(s in 1:length(samples))
      {
        Total_FDR <- append(Total_FDR,results_peak_selection_FDR_all[[s]]$plot_data[1])
        Large_Intensity_delta_FDR <- append(Large_Intensity_delta_FDR,results_peak_selection_FDR_all[[s]]$plot_data[2])
      }
      names(Total_FDR) <- samples
      names(Large_Intensity_delta_FDR) <- samples

      if(plot == T)
      {
        ylim <- c(0,ifelse(max(Total_FDR,na.rm=T)<5,5,max(Total_FDR,na.rm=T)))
        p <- Barplots(Total_FDR,AvgLine = T,digits_average = 1,Name = samples,xlab = "",ylab="FDR [%]",main = "Peak selection FDR - Total",shownumbers = F,ylim=ylim)
        abline(h=5,lty=2,col="red")

        ylim <- c(0,ifelse(max(Large_Intensity_delta_FDR,na.rm=T)<5,5,max(Large_Intensity_delta_FDR,na.rm=T)))
        p <- Barplots(Large_Intensity_delta_FDR,AvgLine = T,digits_average = 1,Name = samples,xlab = "",ylab="FDR [%]",main = "Peak selection FDR - > 2-fold intensity difference",shownumbers = F,ylim=ylim)
        abline(h=5,lty=2,col="red")

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

  FDR_peak_selection <- Peak_selection_FDR(num_features=500,features,samples,peaks,path_to_features,peak_quant,feature_with_background_intensity,peak_selected,delta_mz,delta_rt,peak_min_ion_count,num_peaks_store,alignment_scores_cutoff,n_cores,peak_decision,Alignment_scoring,seed=1,plot=T)

  QC_data$FDR_peak_selection <- FDR_peak_selection

  if(any(!is.na(FDR_peak_selection$Large_Intensity_delta_FDR)))
  {
    for(i in 1:length(FDR_peak_selection$Total_FDR))
    {
      print(paste(samples[i]," - FDR of peak selection: ",round(FDR_peak_selection$Total_FDR[i],digits=1)," % (with > 2-fold wrong abundance after filtering: ",round(FDR_peak_selection$Large_Intensity_delta_FDR[i],digits=1)," %)",sep=""))
    }
  }

  dev.off()

  ###function to plot peak selections
  Plot_feature_quantification <- function(path_to_graph_data,features,selected_feature,samples,peaks,num_peaks_store)
  {
    e <- new.env()

    max <- length(samples)
    pb <- winProgressBar(title = "Plot peak selections per sample",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(s in 1:length(samples))
    {
      load(paste(path_to_graph_data,"\\",samples[s],"_feature_graphs.RData",sep=""),envir=e)

      max2 <- length(selected_feature)
      pb2 <- winProgressBar(title = "Plot peak selections",label=paste( round(0/max2*100, 0),"% done"), min = 0,max = max2, width = 300)
      start_time2 <- Sys.time()
      updatecounter2 <- 0
      time_require2 <- 0
      for(f in selected_feature)
      {
        i <- which(selected_feature == f)
        temp <- e$peaks_graph[[as.character(f)]]
        if(!is.na(temp))
        {
          ###plot KDE
          image(temp$kdemap,main=paste(f,"-",samples[s]),xlab="RT",ylab="m/z",xlim=c(min(temp$kdemap$x),max(temp$kdemap$x)),ylim=c(min(temp$kdemap$y),max(temp$kdemap$y)))
          ###indicate where peaks were detected
          text(temp$maxima_select$RT,temp$maxima_select$mz,1:nrow(temp$maxima_select),col=ifelse(temp$maxima_select$known_peak == 1,"green","black"))
          ###expected window
          sel <- which(features$Feature_name == f)
          RT_expected <- features$RT[sel] + features[sel,which(grepl("RT_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))]
          RT_window_width <- features$RT_length[sel]
          cur_RT_window <- c(RT_expected-(RT_window_width/2),RT_expected+(RT_window_width/2))
          mz_expected <- features$m.z[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))]
          mz_window <- c(features$m.z_range_min[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))],
                         features$m.z_range_max[sel] + features[sel,which(grepl("mz_calibration.",colnames(features)) & grepl(samples[s],colnames(features)))])
          delta_mz <- (mz_window[2]-mz_window[1])/2

          rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2)
          points(RT_expected,mz_expected,pch=4)
          ###indicate which peak was selected
          temp_selected <- peaks[which(peaks$sample == samples[s] & peaks$feature == f & peaks$peak == peaks$selected),]

          if(temp_selected$peak == (num_peaks_store+1)) #selected standard window
          {
            rect(cur_RT_window[1],mz_window[1],cur_RT_window[2],mz_window[2],lty=2,border="blue")
            intensity_temp <- temp$peak_intensities[nrow(temp$peak_intensities),1]
            text(RT_expected,mz_expected-(1.1*delta_mz),round(intensity_temp,2),col="blue")
          }else
          {
            rect_pos_x_1 <- temp$maxima_select$RT[temp_selected$peak]-RT_window_width/2
            rect_pos_y_1 <- temp$maxima_select$mz[temp_selected$peak]-delta_mz
            rect_pos_x_2 <- temp$maxima_select$RT[temp_selected$peak]+RT_window_width/2
            rect_pos_y_2 <- temp$maxima_select$mz[temp_selected$peak]+delta_mz


            rect(rect_pos_x_1,rect_pos_y_1,rect_pos_x_2,rect_pos_y_2,lty=2,border="blue")
            intensity_temp <- temp$peak_intensities[temp_selected$peak,1]
            text(temp_selected$RT,temp_selected$mz-(1.1*delta_mz),round(intensity_temp,2),col="blue")
          }
        }

        updatecounter2 <- updatecounter2 + 1
        if(updatecounter2 >= 1)
        {
          time_elapsed <- difftime(Sys.time(),start_time2,units="secs")
          time_require2 <- (time_elapsed/(i/max2))*(1-(i/max2))
          td <- seconds_to_period(time_require2)
          time_require2 <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb2, i, label=paste( round(i/max2*100, 0)," % done (",i,"/",max2,", Time require: ",time_require2,")",sep = ""))
        }
      }
      close(pb2)
      updatecounter <- updatecounter + 1
      if(updatecounter >= 1)
      {
        time_elapsed <- difftime(Sys.time(),start_time,units="secs")
        time_require <- (time_elapsed/(s/max))*(1-(s/max))
        td <- seconds_to_period(time_require)
        time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

        updatecounter <- 0
        setWinProgressBar(pb, s, label=paste( round(s/max*100, 0)," % done (",s,"/",max,", Time require: ",time_require,")",sep = ""))
      }
    }
    close(pb)
  }
  #
  # f <- c("Feature_16094","Feature_27778")
  # View(peaks[which(peaks$feature%in%f & peaks$peak == peaks$selected),])
  # features_select[which(features_select$Feature_name == f),]
  # Plot_feature_quantification(path_to_graph_data=paste(path_to_mzXML,"\\all_ion_lists\\Extracted feature intensities",output_file_names_add,sep=""),
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
  select_only_known <- which(paste(peaks_raw$sample,peaks_raw$feature) %in% paste(peaks_selected$sample,peaks_selected$feature) )
  peaks_raw <- peaks_raw[select_only_known,]
  #deviation in raw data
  temp_mean_raw <- aggregate(peaks_raw[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_raw$feature),FUN=mean,na.rm=T)
  temp_sd_raw <- aggregate(peaks_raw[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_raw$feature),FUN=sd,na.rm=T)
  #deviation after alignment
  temp_mean_aligned <- aggregate(peaks_selected[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_selected$feature),FUN=mean,na.rm=T)
  temp_sd_aligned <- aggregate(peaks_selected[,c("RT","mz","RT_correct","mz_correct")],by=list(peaks_selected$feature),FUN=sd,na.rm=T)
  #exclude features where sd in uncorrected RT > 1 min
  sel <- which(temp_sd_raw$RT > 1)
  temp_mean_raw <- temp_mean_raw[-sel,]
  temp_sd_raw <- temp_sd_raw[-sel,]
  temp_mean_aligned <- temp_mean_aligned[-sel,]
  temp_sd_aligned <- temp_sd_aligned[-sel,]
  pdf("Performance of feature alignment.pdf",useDingbats = F)
  boxplot(temp_sd_raw$RT,temp_sd_aligned$RT_correct,outline=F,main="Variability of feature RT",names = c("Raw","Aligned"),ylab="SD of RT [min]")
  boxplot(temp_sd_raw$mz,temp_sd_aligned$mz_correct,outline=F,main="Variability of feature mz",names = c("Raw","Aligned"),ylab="SD of mz [Da]")
  dev.off()


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
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  save(temp_results,file = "Quantification_raw_results_with_scores.RData")

  crap <- gc(F)
  ###load temporary results
  # setwd(paste(path_to_features,"\\Temporary_files",sep=""))
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
    mono_iso_alignment_summary<-mono_iso_alignment_summary[selection,]
  }

  ###remove isotope features which never show significant quantification
  selection <- which(grepl("_i",features$Feature_name) & rowMins(as.matrix(pval_signal_with_background_quant),na.rm=T) > Quant_pVal_cut)

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
    mono_iso_alignment_summary<-mono_iso_alignment_summary[-selection,]
  }
  print(paste("Removed ",length(selection)," isotope features (",round(length(selection)/total_length*100,digits=1)," %) as they don´t show significant ion accumulation in any sample.",sep=""))

  ####impute missing values based on generalized additive model
  impute_feature_level_quant <- function(data,features,background_intensity_GAM_table_per_feature_sum)
  {
    suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    quant_data <- data.table::copy(data)
    background_intensity_GAM_table_per_feature <- background_intensity_GAM_table_per_feature_sum

    for(c in 1:ncol(quant_data))
    {
      missing_vals_rows <- which(is.na(quant_data[,c]))

      impute_vals <- background_intensity_GAM_table_per_feature[match(rownames(quant_data)[missing_vals_rows],rownames(background_intensity_GAM_table_per_feature)),c]

      impute_vals <- log10(2^impute_vals)

      set(quant_data,i=as.integer(missing_vals_rows),j=as.integer(c),value=impute_vals)
    }

    return(quant_data)
  }

  feature_with_background_intensity_imputed <- impute_feature_level_quant(data=feature_with_background_intensity,
                                                                          features=features,
                                                                          background_intensity_GAM_table_per_feature_sum=QC_data$Decoy_feature_parameters$background_intensity_GAM_table_per_feature_sum)

  ##now perform filtering of quantifications based on scores
  #remove all features with too high variability
  selection <- rownames(alignment_variability_score)[which(rowMins(as.matrix(alignment_variability_score),na.rm=T) < alignment_variability_score_cutoff)]
  if(length(selection)>0)
  {
    print(paste("Removed ",length(which(!grepl("_i",selection)))," (",length(which(grepl("_i",selection)))," isotope)"," features from quantification results due to too high variability in alignment between samples.",sep=""))

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
    mono_iso_alignment_summary<-mono_iso_alignment_summary[-selection,]
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

    print(paste(samples[s],": Removed ",length(which(!grepl("_i",names(selection))))," (",length(which(grepl("_i",names(selection))))," isotope)"," quantifications (",round(length(selection)/nrow(alignment_scores_peaks_correct)*100,digits=1)," %) due to uncertain peak selection",sep=""))
  }

  #remove isotope quantifications for uncertain mono-+1-isos but keep mono quantification if otherwise certain of peak selection
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

      print(paste(samples[s],": Removed ",length(which(!is.na(temp_quant)))," isotope quantifications (",round(length(which(!is.na(temp_quant)))/nrow(alignment_scores_peaks_correct)*100,digits=1)," %) due to discrepancy in peak selection between mono-/+1-isotope features",sep=""))
    }else ### not enough known peaks ... remove all isotopes
    {
      selection <- which(grepl("_i",rownames(feature_with_background_intensity)))
      feature_with_background_intensity[selection,s] <- NA
      feature_with_background_intensity_imputed[selection,s] <- NA

      print(paste(samples[s],": Removed all isotope quantifications due to too few known true peaks.",sep=""))
    }

  }

  #calculate fraction of missing values before and after imputation
  total <- nrow(feature_with_background_intensity)*ncol(feature_with_background_intensity)
  missing <- length(which(is.na(as.numeric(as.matrix(feature_with_background_intensity)))))
  print(paste("Quantification without imputation: ",round(missing/total*100,digits=1)," % missing values",sep=""))
  missing <- length(which(is.na(as.numeric(as.matrix(feature_with_background_intensity_imputed)))))
  print(paste("Quantification with imputation: ",round(missing/total*100,digits=1)," % missing values",sep=""))


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
  #   pb <- winProgressBar(title = "Detect significantly correlating PMPs",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
  #       td <- seconds_to_period(time_require)
  #       time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))
  #
  #       updatecounter <- 0
  #       setWinProgressBar(pb, pbi, label=paste( round(pbi/max*100, 0)," % done (",pbi,"/",max,", Time require: ",time_require,")",sep = ""))
  #     }
  #   }
  #   close(pb)
  #
  #   remove_list <- remove_list[1:count_remove]
  #
  #   print(paste(nrow(PMPs$features)-length(remove_list),"of PMPs show significant correlation with peptides of respective proteins."))
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
  #   remove_list <- which(rowMins(as.matrix(quant_pvals_temp),na.rm=T)>0.01)
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
  # setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  # write.table(features_intensity,paste("Features_quantification_signal_only",output_file_names_add,".txt",sep=""))
  # write.table(Ioncount_feature_sample_matrix,paste("Features_quantification_ion_counts_signal_only",output_file_names_add,".txt",sep=""))
  # write.table(feature_with_background_intensity,paste("Features_quantification_signal_background",output_file_names_add,".txt",sep=""))
  # write.table(feature_with_background_intensity,paste("Features_quantification_signal_background_imputed",output_file_names_add,".txt",sep=""))
  # write.table(Ioncount_feature_with_background_intensity,paste("Features_quantification_ion_count_signal_background",output_file_names_add,".txt",sep=""))
  # write.table(Alignment_scores,paste("Alignment_scores",output_file_names_add,".txt",sep=""))
  # write.table(Quant_pvals_signal_with_background,paste("Features_quantification_pvals_signal_background",output_file_names_add,".txt",sep=""))
  #

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
  setwd(paste(path_to_features,"\\Temporary_files",sep=""))
  save(temp_results,file = "Quantification_raw_results_with_scores_filtered.RData")

  # path_to_features <- "F:\\9_Spike_in_data_sets\\4_spike-in human Shen\\Requant\\19 - DDAicer reprocessed"
  # path_to_MaxQ_output <- "F:\\9_Spike_in_data_sets\\4_spike-in human Shen\\MaxQuant"
  # output_file_names_add <- "_DDAiceR_analysis"
  # n_cores <- 3
  # abundance_estimation_correction <- T
  # calc_protein_LFQ <- T
  # calc_peptide_LFQ <- F
  #
  # setwd(paste(path_to_features,"\\Temporary_files",sep=""))
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
    suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    suppressWarnings(suppressMessages(library(dplyr,quietly = T)))

    if(is.na(corr_factor)) ###no correction factor specified so correction will be determined and then applied
    {
      ###select monoisotopic features with significant quantification for determination of general trends in difference between MaxQ and Requant quantification
      select <- which(!grepl("_i|_d",features$Feature_name))
      feature_quant <- feature_sample_matrix_requantified[select,]
      feature_quant_pval <- pval_quant[select,]
      feature_quant[is.na(feature_quant_pval) | feature_quant_pval > 0.1] <- NA
      features_select <- features[select,]
      ###determine deviations between Requant and MaxQ in abundance estimations
      Requant_peptides_quant_seq <- aggregate(10^feature_quant,by=list(Sequence=features_select$Sequence),FUN=sum)
      Requant_peptides_quant_seq[,-1] <- log2(Requant_peptides_quant_seq[,-1])
      ###bring MaxQ results table into same order as Requant output
      if(ncol(MaxQ_peptides_quant) != ncol(Requant_peptides_quant_seq))
      {
        colnames(MaxQ_peptides_quant) <- gsub("Intensity.","",colnames(MaxQ_peptides_quant))
        ordering <- vector("numeric",ncol(MaxQ_peptides_quant))
        for(i in 1:ncol(MaxQ_peptides_quant))
        {
          overlap <- grepl(paste(colnames(MaxQ_peptides_quant)[i],"$",sep=""),colnames(Requant_peptides_quant_seq))
          ordering[i] <- ifelse(any(overlap),which(overlap == T),NA)
        }
        MaxQ_peptides_quant <- MaxQ_peptides_quant[,which(!is.na(ordering))]
        ordering <- ordering[!is.na(ordering)]
        MaxQ_peptides_quant <- MaxQ_peptides_quant[,ordering]
      }
      comb <- full_join(MaxQ_peptides_quant,Requant_peptides_quant_seq,by="Sequence")

      ncolumns <- ncol(feature_sample_matrix_requantified)

      dat <- as.data.frame(matrix(ncol=2,nrow=nrow(comb)*ncolumns))
      dat[,1] <- as.numeric(dat[,1])
      dat[,2] <- as.numeric(dat[,2])
      for(i in 1:ncolumns)
      {
        start_index <- (i-1)*nrow(comb)+1
        end_index <- i*nrow(comb)
        temp <- as.data.frame(cbind(comb[,i+1],comb[,i+1+ncolumns]))
        set(dat,as.integer(start_index:end_index),j=as.integer(1:2),value=temp)
      }
      colnames(dat) <- c("MaxQ","Requant")
      dat <- na.omit(dat)
      if(main != "")main=paste(" ",main,sep="")
      pdf(paste("Correct feature abundance estimations",main,".pdf",sep=""))
      smoothScatter(dat[,1],dat[,2],ylab="Requant, log2",xlab="MaxQ, log2",main="Correlation of MaxQ vs Requant on peptide level")
      abline(a=0,b=1)
      temp <- dat[c(order(dat[,2],decreasing = T)[1:(0.1*nrow(dat))],order(dat[,2],decreasing = F)[1:(0.1*nrow(dat))]),]
      y <- temp[,2]
      x <- temp[,1]
      fit <- lm(y~x)
      abline(fit,lty=2,col="red")
      m <- summary(fit)$coefficients[2,1]
      intercept <- summary(fit)$coefficients[1,1]
      Rsq <- summary(fit)$r.squared
      posx <- par("usr")[1]+(par("usr")[2]-par("usr")[1])*0.15
      posy <- par("usr")[4]-(par("usr")[4]-par("usr")[3])*0.1
      text(posx,posy,labels = paste("R? =",round(Rsq,digits=2),"\nslope =",round(m,digits=2)))
      dev.off()

      ###correction
      features_sample_matrix_corrected <- as.data.frame(log2(10^feature_sample_matrix_requantified))

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
      features_sample_matrix_corrected <- as.data.frame(log2(10^feature_sample_matrix_requantified))

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
    pb <- winProgressBar(title = "Prepare Top3 quantification",label=paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    close(pb)
    suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
    suppressWarnings(suppressMessages(library(matrixStats,quietly = T)))
    suppressWarnings(suppressMessages(library(stringr,quietly = T)))

    feature_sample_matrix_requantified <- as.data.frame(feature_sample_matrix_requantified)
    if(!is.null(Quant_pvals))Quant_pvals <- as.data.frame(Quant_pvals)
    if(!is.null(Alignment_scores))Alignment_scores <- as.data.frame(Alignment_scores)
    if(!is.null(S2B))S2B <- as.data.frame(S2B)
    ##top3 method
    ###input matrix with samples in cols and rows correspond to unique peptides (log2 summed intensity over charge state and modification) of a respective protein
    Top3_quant <- function(pep_matrix,features_temp,Alignment_scores_temp=NULL,Quant_pvals_temp=NULL,S2B_temp=NULL)
    {
      sequence <- features_temp$Sequence

      if(length(sequence)>0)
      {
        pep_matrix <- aggregate(2^pep_matrix, by=list(Sequence=sequence), FUN=sum,na.rm=T)
        pep_matrix <- pep_matrix[,-1]
        pep_matrix[pep_matrix==0] <- NA
        pep_matrix <- log2(pep_matrix)

        if(!is.null(Alignment_scores_temp))
        {
          score_matrix <- aggregate(Alignment_scores_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          score_matrix <- score_matrix[,-1]
          score_matrix[score_matrix==0] <- NA
        }else
        {
          score_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(score_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(Quant_pvals_temp))
        {
          pval_matrix <- aggregate(Quant_pvals_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          pval_matrix <- pval_matrix[,-1]
          pval_matrix[pval_matrix==0] <- NA
        }else
        {
          pval_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(pval_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(S2B_temp))
        {
          S2B_matrix <- aggregate(S2B_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
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
            median_score <- median(top3_score,na.rm=T)
          }else
          {
            median_score <- NA
          }

          if(!is.null(Quant_pvals_temp))
          {
            top3_pval <- pval_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
            median_pval <- median(top3_pval,na.rm=T)
          }else
          {
            median_pval <- NA
          }

          if(!is.null(S2B_temp))
          {
            top3_S2B <- S2B_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c][1:3]
            median_S2B <- median(top3_S2B,na.rm=T)
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
          top3_res <- append(top3_res,log2(sums))
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
      names(top3_score_res) <- paste(colnames(score_matrix),"_median_score",sep="")
      names(top3_pval_res) <- paste(colnames(pval_matrix),"_median_pvals",sep="")
      names(top3_S2B_res) <- paste(colnames(pval_matrix),"_median_S2B",sep="")

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
      selection <- which(rowMins(as.matrix(Quant_pvals_temp),na.rm=T) < quant_pvalue_cutoff)

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
        unique_proteins <- sort(unique(as.character(str_split(features_temp$Protein,"\\||;",simplify = T))))
      }else
      {
        unique_proteins <- sort(unique(features_temp$Protein[which(!grepl("\\||;",features_temp$Protein))]))
      }

      if(any(unique_proteins == ""))unique_proteins <- unique_proteins[-which(unique_proteins == "")]

      protein_TOP3 <- as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=length(unique_proteins),0))
      rownames(protein_TOP3) <- unique_proteins
      colnames(protein_TOP3) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))

      max <- nrow(protein_TOP3)
      pb <- winProgressBar(title = "Perform Top3 quantification",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          set(protein_TOP3,as.integer(i),as.integer(1:ncol(protein_TOP3)),value=as.list(c(length(ind),as.numeric(Top3_quant(pep_matrix = feature_sample_matrix_requantified_temp[ind,],features_temp = features_temp[ind,],Alignment_scores_temp = Alignment_scores[ind,],Quant_pvals_temp=Quant_pvals[ind,],S2B_temp=S2B_temp[ind,])))))
        }else
        {
          set(protein_TOP3,as.integer(i),as.integer(1),value=0)
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)
    }else
    {
      protein_TOP3 <- as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=0,0))
      colnames(protein_TOP3) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))

    }

    return(protein_TOP3)
  }

  Total_Protein_Quant <- function(features,feature_sample_matrix_requantified,Alignment_scores=NULL,Quant_pvals=NULL,S2B=NULL,use_overlapping=T,min_peps=2,quant_pvalue_cutoff=0.1,use_isotope_pmps=F)
  {
    pb <- winProgressBar(title = "Prepare Total quantification",label=paste( round(0/1*100, 0),"% done"), min = 0,max = 1, width = 300)
    close(pb)
    suppressWarnings(suppressMessages(library(data.table,quietly = T)))
    suppressWarnings(suppressMessages(library(lubridate,quietly = T)))
    suppressWarnings(suppressMessages(library(matrixStats,quietly = T)))
    suppressWarnings(suppressMessages(library(stringr,quietly = T)))

    feature_sample_matrix_requantified <- as.data.frame(feature_sample_matrix_requantified)
    if(!is.null(Quant_pvals))Quant_pvals <- as.data.frame(Quant_pvals)
    if(!is.null(Alignment_scores))Alignment_scores <- as.data.frame(Alignment_scores)
    if(!is.null(S2B))S2B <- as.data.frame(S2B)
    ##total method
    ###input matrix with samples in cols and rows correspond to unique peptides (log2 summed intensity over charge state and modification) of a respective protein
    Total_quant <- function(pep_matrix,features_temp,Alignment_scores_temp=NULL,Quant_pvals_temp=NULL,S2B_temp=NULL,Quant_cutoff=4)
    {
      sequence <- features_temp$Sequence

      if(length(sequence)>0)
      {
        pep_matrix <- aggregate(2^pep_matrix, by=list(Sequence=sequence), FUN=sum,na.rm=T)
        pep_matrix <- pep_matrix[,-1]
        pep_matrix[pep_matrix==0] <- NA
        pep_matrix <- log2(pep_matrix)

        if(!is.null(Alignment_scores_temp))
        {
          score_matrix <- aggregate(Alignment_scores_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          score_matrix <- score_matrix[,-1]
          score_matrix[score_matrix==0] <- NA
        }else
        {
          score_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(score_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(Quant_pvals_temp))
        {
          pval_matrix <- aggregate(Quant_pvals_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
          pval_matrix <- pval_matrix[,-1]
          pval_matrix[pval_matrix==0] <- NA
        }else
        {
          pval_matrix <- matrix(nrow=nrow(pep_matrix),ncol=ncol(pep_matrix),NA)
          colnames(pval_matrix) <- colnames(pep_matrix)
        }
        if(!is.null(S2B_temp))
        {
          S2B_matrix <- aggregate(S2B_temp, by=list(Sequence=sequence), FUN=mean,na.rm=T)
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
            median_score <- weighted.mean(total_score,total,na.rm=T)#median(total_score,na.rm=T)
          }else
          {
            median_score <- NA
          }

          if(!is.null(Quant_pvals_temp))
          {
            total_pval <- pval_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
            median_pval <- weighted.mean(total_pval,total,na.rm=T)#median(total_pval,na.rm=T)
          }else
          {
            median_pval <- NA
          }

          if(!is.null(S2B_temp))
          {
            total_S2B <- S2B_matrix[order(pep_matrix[,c],na.last = T,decreasing = T),c]
            median_S2B <- weighted.mean(total_S2B,total,na.rm=T)#median(total_S2B,na.rm=T)
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

          total_res <- append(total_res,log2(sums))
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
      names(total_score_res) <- paste(colnames(score_matrix),"_median_score",sep="")
      names(total_pval_res) <- paste(colnames(pval_matrix),"_median_pvals",sep="")
      names(total_S2B_res) <- paste(colnames(pval_matrix),"_median_S2B",sep="")

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
      selection <- which(rowMins(as.matrix(Quant_pvals_temp),na.rm=T) < quant_pvalue_cutoff)

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
        unique_proteins <- sort(unique(as.character(str_split(features_temp$Protein,"\\||;",simplify = T))))
      }else
      {
        unique_proteins <- sort(unique(features_temp$Protein[which(!grepl("\\||;",features_temp$Protein))]))
      }

      if(any(unique_proteins == ""))unique_proteins <- unique_proteins[-which(unique_proteins == "")]

      protein_total <- as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=length(unique_proteins),0))
      rownames(protein_total) <- unique_proteins
      colnames(protein_total) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))

      max <- nrow(protein_total)
      pb <- winProgressBar(title = "Perform total quantification",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)
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
          set(protein_total,as.integer(i),as.integer(1:ncol(protein_total)),value=as.list(c(length(ind),as.numeric(Total_quant(pep_matrix = feature_sample_matrix_requantified_temp[ind,],features_temp = features_temp[ind,],Alignment_scores_temp = Alignment_scores[ind,],Quant_pvals_temp=Quant_pvals[ind,],S2B_temp=S2B_temp[ind,])))))
        }else
        {
          set(protein_total,as.integer(i),as.integer(1),value=0)
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10)
        {
          time_elapsed <- difftime(Sys.time(),start_time,units="secs")
          time_require <- (time_elapsed/(i/max))*(1-(i/max))
          td <- seconds_to_period(time_require)
          time_require <- sprintf('%02d:%02d:%02d', td@hour, lubridate::minute(td), round(lubridate::second(td),digits=0))

          updatecounter <- 0
          setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,", Time require: ",time_require,")",sep = ""))
        }
      }
      close(pb)
    }else
    {
      protein_total <- as.data.frame(matrix(ncol=1+(4*ncol(feature_sample_matrix_requantified)),nrow=0,0))
      colnames(protein_total) <- c("num_quant_features",colnames(feature_sample_matrix_requantified),paste("median_alignment_score_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_quant_pvals_",colnames(feature_sample_matrix_requantified),sep=""),paste("median_S2B_",colnames(feature_sample_matrix_requantified),sep=""))
    }

    return(protein_total)
  }

  crap <- gc(F)
  ###perform protein level aggregation
  if(!is.null(path_to_MaxQ_output))
  {
    ###load MaxQ peptide results
    MaxQ_peptides <- read.table(paste(path_to_MaxQ_output,"\\peptides.txt",sep=""),sep="\t",header=T)
    MaxQ_peptides <- MaxQ_peptides[-which(MaxQ_peptides$Potential.contaminant == "+" | MaxQ_peptides$Reverse == "+"),]
    MaxQ_peptides_quant <- MaxQ_peptides[,which(grepl("Intensity\\.",colnames(MaxQ_peptides)))]
    MaxQ_peptides_quant[MaxQ_peptides_quant==0] <- NA
    MaxQ_peptides_quant <- log2(MaxQ_peptides_quant)
    MaxQ_peptides_leading_razor <- data.frame(Sequence=MaxQ_peptides$Sequence,Leading_razor=MaxQ_peptides$Leading.razor.protein)
    MaxQ_peptides_leading_razor$Leading_razor <- as.character(MaxQ_peptides_leading_razor$Leading_razor)
    MaxQ_peptides <- data.frame(Sequence=MaxQ_peptides$Sequence,MaxQ_peptides_quant)

    MaxQ_protein_groups <- read.table(paste(path_to_MaxQ_output,"\\proteinGroups.txt",sep=""),sep="\t",header=T)

    ###check if IDs were correctly parsed, if not, try to parse with SwissProt or Trembl
    if(!any(colnames(MaxQ_protein_groups) == "Gene.names") & any(colnames(MaxQ_protein_groups) == "Protein.IDs")) ##not correctly parsed but contains trembl or swissprot fasta headers
    {
      print("Fasta file was not correctly parsed during search. Try to paste fasta headers ...")
      if(any(grepl(">sp|>tr",MaxQ_protein_groups$Fasta.headers)))
      {
        print("Detected Swiss-Prot and/or TrEMBL fasta headers.")
        ##protein level
        ###parsing was not performed correctly so we have to try to do this here expecting swissprot or trembl fasta headers
        MaxQ_protein_groups$Gene.names <- ""
        MaxQ_protein_groups$Organism <- ""
        MaxQ_protein_groups$Protein.IDs <- as.character(MaxQ_protein_groups$Protein.IDs)

        fasta_headers <- str_split(gsub(">","",MaxQ_protein_groups$Fasta.headers),"\\|",simplify = T)

        ##extract gene name and species information
        GN <- vector("character",nrow(MaxQ_protein_groups))
        #OS <- vector("character",nrow(MaxQ_protein_groups))
        ID <- vector("character",nrow(MaxQ_protein_groups))
        GN_start <- unlist(lapply(paste(gregexpr("GN=",MaxQ_protein_groups$Fasta.headers),sep=","), `[[`, 1))
        GN_start <- gsub("c\\(|\\)","",GN_start)
        #OS_start <- unlist(lapply(paste(gregexpr("OS=",MaxQ_protein_groups$Fasta.headers),sep=","), `[[`, 1))
        #OS_start <- gsub("c\\(|\\)","",OS_start)
        for(i in 1:nrow(MaxQ_protein_groups))
        {
          if(GN_start[i] != "-1")
          {
            indices <- as.numeric(unlist(str_split(GN_start[i],","))) + 3
            stop <- regexpr(" ",substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+10))
            GN_temp <- substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+stop-2)
            GN[i] <- paste(GN_temp,collapse = ";")
          }
          # if(OS_start[i] != "-1")
          # {
          #   indices <- as.numeric(unlist(str_split(OS_start[i],","))) + 3
          #   stop <- regexpr("GN=",substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+100))
          #   OS_temp <- substring(MaxQ_protein_groups$Fasta.headers[i],indices,indices+stop-3)
          #   OS[i] <- paste(OS_temp,collapse = ";")
          # }
          header <- fasta_headers[i,c(2,4,6)][which(fasta_headers[i,c(2,4,6)] != "")]
          if(length(header)>0)ID[i] <- paste(header,collapse=";")
        }
        CON_REV <- !grepl("^CON_|^REV_",MaxQ_protein_groups$Protein.IDs)
        MaxQ_protein_groups$Protein.IDs <- ifelse(CON_REV,ID, MaxQ_protein_groups$Protein.IDs)
        MaxQ_protein_groups$Gene.names <- ifelse(CON_REV,GN,"")
        #MaxQ_protein_groups$Organism <- ifelse(CON_REV,OS,"")
        MaxQ_protein_groups$Majority.protein.IDs <- MaxQ_protein_groups$Protein.IDs

        ##peptide level
        if(any(grepl("\\|",MaxQ_peptides_leading_razor$Leading_razor)))
        {
          temp <- str_split(MaxQ_peptides_leading_razor$Leading_razor,"\\|",simplify = T)
          MaxQ_peptides_leading_razor$Leading_razor <- ifelse(temp[,2] != "",temp[,2],MaxQ_peptides_leading_razor$Leading_razor)
        }

        ##requantification features
        temp <- str_split(features$Protein,"\\||;",simplify=T)
        ID <- vector("character",nrow(temp))
        for(i in 1:nrow(temp))
        {
          sel <- which(temp[i,] == "sp")+1
          if(length(sel)>0)
          {
            ID[i] <- paste(temp[i,sel],collapse=";")
          }else
          {
            sel <- which(temp[i,] != "")
            ID[i] <- paste(temp[i,sel],collapse=";")
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
      feature_with_background_intensity <- log2(10^feature_with_background_intensity)
      feature_with_background_intensity_imputed <- log2(10^feature_with_background_intensity_imputed)

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
        suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))
        suppressWarnings(suppressMessages(library(data.table,quietly = T)))

        calculate_LFQ <- function(peptide_quant_data,min_num_ratios=2,num_ratio_samples=NA)
        {
          suppressWarnings(suppressMessages(library(data.table,quietly = T)))
          suppressWarnings(suppressMessages(library(robustbase,quietly = T)))
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
                  val <- (log2(ratio_mat[r,c])-log2(par[r])+log2(par[c]))^2
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
          ratio_mat <- as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data)))
          if(ncol(ratio_mat)>1)
          {
            for(c in 1:(ncol(ratio_mat)-1))
            {
              for(r in (c+1):nrow(ratio_mat))
              {
                ratios <- peptide_quant_data[,r]/peptide_quant_data[,c]
                if(length(which(!is.na(ratios))) >= min_num_ratios){ratio_mat[r,c] <- median(ratios,na.rm=T)}
              }
            }
          }
          if(is.na(num_ratio_samples)) ###calculate ratios over all samples
          {
            #define start parameter
            start_par <- c(rep(1,ncol(ratio_mat)))
            #now find optimum in ratios to best recover true observed ratios between samples
            res_ratio <- NULL
            try(res_ratio <- optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)
            if(!is.null(res_ratio))
            {
              #normalize ratios to sample with highest intensity
              ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample == max(totalsum_per_sample))]
              #finally calculate log2 lfq protein intensities per sample
              lfq <- log2(ratio_norm*totalsum_per_sample[which(totalsum_per_sample == max(totalsum_per_sample))])
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
            lfq <- as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data),0))
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
              try(res_ratio <- optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat_temp,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)

              if(!is.null(res_ratio))
              {
                #normalize ratios to sample with highest intensity
                ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp,na.rm=T))]
                #finally calculate log2 lfq protein intensities per sample
                temp_lfq <- log2(ratio_norm*totalsum_per_sample_temp[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp))])

                set(lfq,as.integer(s),as.integer(c(s,sort(samples_for_comparison_per_sample))),as.list(temp_lfq))
              }

            }
            lfq[lfq==0] <- NA
            lfq <- colMedians(as.matrix(lfq),na.rm=T)

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
            fmtstr <- paste(fmtstr, collapse = ":")
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
              cat(paste(rep.int(char, nb - .nb), collapse = ""),
                  file = file)
              flush.console()
            }
            else if (.nb > nb) {
              cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
            }
            .nb <<- nb
          }
          up2 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb <= nb) {
              cat("\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
            }
            else {
              cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
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
            cat(paste(c("\r  |", rep.int(" ", nw * width + 6)),
                      collapse = ""), file = file)
            cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ",
                                                            nw * (width - nb)), sprintf("| %3d%%", pc), ", ETA ",
                        ETAstr), collapse = ""), file = file)
            flush.console()
            .nb <<- nb
            .pc <<- pc
          }
          getVal <- function() .val
          kill <- function() if (!.killed) {
            cat("\n", file = file)
            flush.console()
            .killed <<- TRUE
          }
          up <- switch(style, up1, up2, up3)
          up(initial, T)
          structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
        }

        sel <- which(!grepl(";|,",features$Sequence))
        unique_peptides <- sort(unique(paste(features$Sequence[sel],features$Modifications[sel],sep="_")))
        unique_seq <- substr(unique_peptides,1,regexpr("_",unique_peptides)-1)
        unique_mod <- substr(unique_peptides,regexpr("_",unique_peptides)+1,nchar(unique_peptides))

        LFQ_peptide_quant <- as.data.frame(matrix(ncol=ncol(features_quant),nrow=length(unique_peptides),0))
        colnames(LFQ_peptide_quant) <- colnames(features_quant)

        ###Perform LFQ quantification
        cl <- makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        registerDoSNOW(cl)
        iterations <- nrow(LFQ_peptide_quant)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        start <- Sys.time()
        print(paste(label," (",Sys.time(),")",sep=""))
        pb <- txtProgressBar(max = iterations, style = 3)
        res_LFQ <- foreach(i=1:nrow(LFQ_peptide_quant),.options.snow = opts) %dopar%
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
              names(res) <- paste("V",1:length(res),sep="")
            }
            return(res)
          }
        stopCluster(cl)
        close(pb)

        #Combine results
        for(i in 1:length(res_LFQ))
        {
          if(length(res_LFQ[[i]])>0)
          {
            set(LFQ_peptide_quant,as.integer(i),as.integer(1:ncol(LFQ_peptide_quant)),as.list(as.numeric(res_LFQ[[i]])))
          }
        }

        #add information about number of quant features per protein
        count_quant_features <- plyr::count(sort(paste(features$Sequence[sel],features$Modifications[sel],sep="_")))

        LFQ_peptide_quant <- data.frame(Sequence=unique_seq,
                                        Modifications=unique_mod,
                                        Protein=features$Protein[match(unique_seq,features$Sequence)],
                                        num_quant_features=count_quant_features$freq[match(unique_peptides,count_quant_features$x)],
                                        LFQ_peptide_quant)

        colnames(LFQ_peptide_quant) <- gsub("^X","",colnames(LFQ_peptide_quant))

        end <- Sys.time()
        print(paste("Finished peptide LFQ-quantification (",Sys.time(),")",sep=""))
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
        suppressWarnings(suppressMessages(library(doSNOW,quietly = T)))
        suppressWarnings(suppressMessages(library(data.table,quietly = T)))

        calculate_LFQ <- function(peptide_quant_data,min_num_ratios=2,num_ratio_samples=NA)
        {
          suppressWarnings(suppressMessages(library(data.table,quietly = T)))
          suppressWarnings(suppressMessages(library(robustbase,quietly = T)))
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
                  val <- (log2(ratio_mat[r,c])-log2(par[r])+log2(par[c]))^2
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
          ratio_mat <- as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data)))
          if(ncol(ratio_mat)>1)
          {
            for(c in 1:(ncol(ratio_mat)-1))
            {
              for(r in (c+1):nrow(ratio_mat))
              {
                ratios <- peptide_quant_data[,r]/peptide_quant_data[,c]
                if(length(which(!is.na(ratios))) >= min_num_ratios){ratio_mat[r,c] <- median(ratios,na.rm=T)}
              }
            }
          }
          if(is.na(num_ratio_samples)) ###calculate ratios over all samples
          {
            #define start parameter
            start_par <- c(rep(1,ncol(ratio_mat)))
            #now find optimum in ratios to best recover true observed ratios between samples
            res_ratio <- NULL
            try(res_ratio <- optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)
            if(!is.null(res_ratio))
            {
              #normalize ratios to sample with highest intensity
              ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample == max(totalsum_per_sample))]
              #finally calculate log2 lfq protein intensities per sample
              lfq <- log2(ratio_norm*totalsum_per_sample[which(totalsum_per_sample == max(totalsum_per_sample))])
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
            lfq <- as.data.frame(matrix(nrow=ncol(peptide_quant_data),ncol=ncol(peptide_quant_data),0))
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
              try(res_ratio <- optim(par=start_par, fn=least_square_error, ratio_mat=ratio_mat_temp,lower = 1,upper=100,method = "L-BFGS-B"),silent = T)

              if(!is.null(res_ratio))
              {
                #normalize ratios to sample with highest intensity
                ratio_norm <- res_ratio$par / res_ratio$par[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp,na.rm=T))]
                #finally calculate log2 lfq protein intensities per sample
                temp_lfq <- log2(ratio_norm*totalsum_per_sample_temp[which(totalsum_per_sample_temp == max(totalsum_per_sample_temp))])

                set(lfq,as.integer(s),as.integer(c(s,sort(samples_for_comparison_per_sample))),as.list(temp_lfq))
              }

            }
            lfq[lfq==0] <- NA
            lfq <- colMedians(as.matrix(lfq),na.rm=T)

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
            fmtstr <- paste(fmtstr, collapse = ":")
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
              cat(paste(rep.int(char, nb - .nb), collapse = ""),
                  file = file)
              flush.console()
            }
            else if (.nb > nb) {
              cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
            }
            .nb <<- nb
          }
          up2 <- function(value) {
            if (!is.finite(value) || value < min || value > max)
              return()
            .val <<- value
            nb <- round(width * (value - min)/(max - min))
            if (.nb <= nb) {
              cat("\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
            }
            else {
              cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""),
                  "\r", paste(rep.int(char, nb), collapse = ""),
                  sep = "", file = file)
              flush.console()
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
            cat(paste(c("\r  |", rep.int(" ", nw * width + 6)),
                      collapse = ""), file = file)
            cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ",
                                                            nw * (width - nb)), sprintf("| %3d%%", pc), ", ETA ",
                        ETAstr), collapse = ""), file = file)
            flush.console()
            .nb <<- nb
            .pc <<- pc
          }
          getVal <- function() .val
          kill <- function() if (!.killed) {
            cat("\n", file = file)
            flush.console()
            .killed <<- TRUE
          }
          up <- switch(style, up1, up2, up3)
          up(initial, T)
          structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
        }

        unique_proteins <- sort(unique(features$Protein[which(!grepl(";|,",features$Protein))]))
        unique_proteins <- unique_proteins[which(unique_proteins != "")]

        LFQ_protein_quant <- as.data.frame(matrix(ncol=ncol(features_quant),nrow=length(unique_proteins),0))
        colnames(LFQ_protein_quant) <- colnames(features_quant)
        rownames(LFQ_protein_quant) <- unique_proteins

        ###Perform LFQ quantification
        cl <- makeCluster(n_cores)#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        registerDoSNOW(cl)
        iterations <- nrow(LFQ_protein_quant)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        start <- Sys.time()
        print(paste(label," (",Sys.time(),")",sep=""))
        pb <- txtProgressBar(max = iterations, style = 3)
        res_LFQ <- foreach(i=1:nrow(LFQ_protein_quant),.options.snow = opts) %dopar%
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
              names(res) <- paste("V",1:length(res),sep="")
            }

            return(res)
          }
        stopCluster(cl)
        close(pb)

        #Combine results
        for(i in 1:length(res_LFQ))
        {
          if(length(res_LFQ[[i]])>0)
          {
            set(LFQ_protein_quant,as.integer(i),as.integer(1:ncol(LFQ_protein_quant)),as.list(as.numeric(res_LFQ[[i]])))
          }
        }
        #add information about number of quant features per protein
        count_quant_features <- plyr::count(features$Protein[which(!grepl(";|\\||,",features$Protein))])

        LFQ_protein_quant <- data.frame(num_quant_features=count_quant_features$freq[match(rownames(LFQ_protein_quant),count_quant_features$x)],LFQ_protein_quant)

        end <- Sys.time()
        print(paste("Finished LFQ-quantification (",Sys.time(),")",sep=""))
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
    cl <- makeCluster(ifelse(n_cores < 4,n_cores,4))#as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
    registerDoParallel(cl)
    res <- foreach(i=1:4) %dopar%
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
    stopCluster(cl)

    Top3_quant_with_background <- res[[1]]
    Total_quant_with_background <- res[[2]]
    Top3_quant_with_background_imputed <- res[[3]]
    Total_quant_with_background_imputed <- res[[4]]

    ###Match Gene names to Uniprot Identifier
    if(any(colnames(MaxQ_protein_groups) == "Gene.names") & any(colnames(MaxQ_protein_groups) == "Protein.IDs"))
    {
      temp <- MaxQ_protein_groups[,c("Gene.names","Protein.IDs")]
      temp <- temp[which(temp$Gene.names != ""),]

      UniProt_to_GeneName <- data.frame(UniProt_ID=unique(as.character(str_split(temp$Protein.IDs,";",simplify = T))),Gene_Name="")
      UniProt_to_GeneName$Gene_Name <- as.character(UniProt_to_GeneName$Gene_Name)

      for(i in 1:nrow(UniProt_to_GeneName))
      {
        UniProt_to_GeneName$Gene_Name[i] <- as.character(temp$Gene.names[which(grepl(UniProt_to_GeneName$UniProt_ID[i],temp$Protein.IDs))])
      }

      Top3_quant_with_background <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Top3_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Top3_quant_with_background),Top3_quant_with_background)
      rownames(Top3_quant_with_background) <- c()

      Total_quant_with_background <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Total_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Total_quant_with_background),Total_quant_with_background)
      rownames(Total_quant_with_background) <- c()

      Top3_quant_with_background_imputed <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Top3_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Top3_quant_with_background_imputed),Top3_quant_with_background_imputed)
      rownames(Top3_quant_with_background_imputed) <- c()

      Total_quant_with_background_imputed <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(Total_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(Total_quant_with_background_imputed),Total_quant_with_background_imputed)
      rownames(Total_quant_with_background_imputed) <- c()

      if(calc_peptide_LFQ == T)
      {
        LFQ_peptide_quant_with_background <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(LFQ_peptide_quant_with_background$Protein,UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=LFQ_peptide_quant_with_background$Protein,LFQ_peptide_quant_with_background[,c(1,2,4,5:ncol(LFQ_peptide_quant_with_background))])
        rownames(LFQ_peptide_quant_with_background) <- c()
        colnames(LFQ_peptide_quant_with_background) <- gsub("^X","",colnames(LFQ_peptide_quant_with_background))

        LFQ_peptide_quant_with_background_imputed <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(LFQ_peptide_quant_with_background_imputed$Protein,UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=LFQ_peptide_quant_with_background_imputed$Protein,LFQ_peptide_quant_with_background_imputed[,c(1,2,4,5:ncol(LFQ_peptide_quant_with_background_imputed))])
        rownames(LFQ_peptide_quant_with_background_imputed) <- c()
        colnames(LFQ_peptide_quant_with_background_imputed) <- gsub("^X","",colnames(LFQ_peptide_quant_with_background_imputed))
      }

      if(calc_protein_LFQ == T)
      {
        LFQ_quant_with_background <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(LFQ_quant_with_background),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(LFQ_quant_with_background),LFQ_quant_with_background)
        rownames(LFQ_quant_with_background) <- c()
        colnames(LFQ_quant_with_background) <- gsub("^X","",colnames(LFQ_quant_with_background))

        LFQ_quant_with_background_imputed <- data.frame(Gene_Name=UniProt_to_GeneName$Gene_Name[match(rownames(LFQ_quant_with_background_imputed),UniProt_to_GeneName$UniProt_ID)],UniProt_Identifier=rownames(LFQ_quant_with_background_imputed),LFQ_quant_with_background_imputed)
        rownames(LFQ_quant_with_background_imputed) <- c()
        colnames(LFQ_quant_with_background_imputed) <- gsub("^X","",colnames(LFQ_quant_with_background_imputed))
      }

    }

    if(calc_protein_LFQ == T)
    {
      utils::write.table(x = LFQ_quant_with_background,file = paste(path_to_features,"\\Proteins_quantification_LFQ",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
      utils::write.table(x = LFQ_quant_with_background_imputed,file = paste(path_to_features,"\\Proteins_quantification_LFQ_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    }
    utils::write.table(x = Top3_quant_with_background,file = paste(path_to_features,"\\Proteins_quantification_Top3",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Top3_quant_with_background_imputed,file = paste(path_to_features,"\\Proteins_quantification_Top3_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Total_quant_with_background,file = paste(path_to_features,"\\Proteins_quantification_Total",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    utils::write.table(x = Total_quant_with_background_imputed,file = paste(path_to_features,"\\Proteins_quantification_Total_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")

    if(calc_peptide_LFQ == T)
    {
      utils::write.table(x = LFQ_peptide_quant_with_background,file = paste(path_to_features,"\\Peptides_quantification_LFQ",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
      utils::write.table(x = LFQ_peptide_quant_with_background_imputed,file = paste(path_to_features,"\\Peptides_quantification_LFQ_imputed",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
    }
  }

  setwd(path_to_features)

  ###finally save feature level quantification
  utils::write.table(x = features,file = paste(path_to_features,"\\Features",output_file_names_add,".tab",sep=""),row.names = F,sep = "\t")
  utils::write.table(x = feature_with_background_intensity,file = paste(path_to_features,"\\Features_quantification",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = feature_with_background_intensity_imputed,file = paste(path_to_features,"\\Features_quantification_imputed",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = pval_signal_with_background_quant,file = paste(path_to_features,"\\Features_quantification_pvals",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = Ioncount_feature_with_background_intensity,file = paste(path_to_features,"\\Features_quantification_ioncount",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = S2B,file = paste(path_to_features,"\\Features_quantification_S2B",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = alignment_variability_score,file = paste(path_to_features,"\\Features_quantification_variability_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = alignment_scores_peaks_correct,file = paste(path_to_features,"\\Features_quantification_alignment_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")
  utils::write.table(x = mono_iso_alignment_summary,file = paste(path_to_features,"\\Features_quantification_mono_iso_alignment_score",output_file_names_add,".tab",sep=""),row.names = T,sep = "\t")

  save(QC_data,file = paste("Temporary_files\\Feature_quantification_QC_data.RData",sep=""))
  options(warn=0)
}




