#' ---------------------------------------------
#' 
#' 
#' Helper functions for AMP analysis
#' 
#' 
#' ---------------------------------------------

library(tcltk2)

Compile_Raven_selns <- function(site_id = character(), 
                                dep_id = character()){
  dir_seln <- tk_choose.dir(caption = "Select folder with selection tables (.txt)")
  seln_tables <- dir_seln |>
    # list all files with .txt extension
    list.files(pattern = ".txt") |>
    map_chr(~paste0(dir_seln,"\\",  .)) |>
    # use map to iterate read.delim() function over each file in the directory
    map(~read.delim(.)) |>
    # rename columns with funky character formatting
    map(~rename(.,"Begin_Time" = "Begin.Time..s.",
                "End_Time" = "End.Time..s.",
                "Low_Freq" = "Low.Freq..Hz.",
                "High_Freq" = "High.Freq..Hz.", 
                "Begin_Date" = "Begin.Date", 
                "Begin_Clock" = "Begin.Clock.Time",
                "End_Clock" = "End.Clock.Time", 
                "Delta_Time_s" = "Delta.Time..s.")) |>
    map(~mutate(., Begin_file_date = first(Begin_Date)))
  
  
  # Now we have selection tables as a list, but we want them all together in one dataframe
  all_selns <- do.call("rbind", seln_tables) |>
    # fix any T's that got converted to logical, then convert all to TRANSIT
    mutate(Behavior = gsub(pattern = "TRUE", replacement = "T", x = Behavior),
           Behavior = gsub(pattern = "T", replacement = "Transit", fixed=TRUE, x = Behavior),
           Behavior = gsub(pattern = "M", replacement = "Maneuver", fixed=TRUE, x = Behavior))
  
  
  return = all_selns
  
}

# Test
# Compile_Raven_selns(site_id = "SWS", dep_id = "test")



Compile_Raven_bio_selns <- function(site_id = character(), 
                                    dep_id = character()){
  dir_seln <- tk_choose.dir(caption = "Select folder with selection tables (.txt)")
  seln_tables <- dir_seln |>
    # list all files with .txt extension
    list.files(pattern = ".txt") |>
    map_chr(~paste0(dir_seln,"\\",  .)) |>
    # use map to iterate read.delim() function over each file in the directory
    map(~read.delim(.)) |>
    # rename columns with funky character formatting
    map(~rename(.,"Begin_Time" = "Begin.Time..s.",
                "End_Time" = "End.Time..s.",
                "Low_Freq" = "Low.Freq..Hz.",
                "High_Freq" = "High.Freq..Hz.", 
                "Begin_Date" = "Begin.Date", 
                "Begin_Clock" = "Begin.Clock.Time",
                "End_Clock" = "End.Clock.Time", 
                "Delta_Time_s" = "Delta.Time..s.",
                "Call_type" = "Call.Type"))
  
  
  # Now we have selection tables as a list, but we want them all together in one dataframe
  all_selns <- do.call("rbind", seln_tables) 
  
  
  return = all_selns
  
}



Compile_ves_detections <- function(site_id = character()){
  dir_detect <- choose.dir(caption = paste0("Select folder with validated detections ", site_id))
  det_tables <- dir_detect |>
    # list all files with .csv extension
    list.files(pattern = ".csv") |>
    map_chr(~paste0(dir_detect,"\\",  .)) |>
    # use map to iterate read_csv() function over each file in the directory
    map(~read_csv(.) |>
          select(Detection, ISOStartTime, ISOEndTime, StartTime, EndTime,
                 start, end, Labels, `new labels`) |>
          # mutate to add deployment ID from filename 
          # use .* operator to subset after "_MA-RI_" and before "_5s_48Hz" in filename
          mutate(start_date_ISO = as_date(ISOStartTime),
                 Dep_ID = sub(".*_MA-RI_","", .x),
                 Dep_ID = sub("_5s_48Hz.*","",Dep_ID))) 
  
  
  # Now we have selection tables as a list, but we want them all together in one dataframe
  all_dets <- do.call("rbind", det_tables) |>
    # add new column for site
    mutate(SITE = site_id)
  
  
  return = all_dets
  
}



getDeploymentInfo <- function(){
  
  ##### Deployment info ###
  site_id <- dlg_input(message = "Site ID, e.g. 'TRE'")$res
  dep_id <- dlg_input(message = "Deployment ID, YYYYMM")$res
  
  # Start and end dates of deployment
  start_dep_date <- as_date(dlg_input(message = "Start date: YYYY-MM-DD")$res)
  end_dep_date <- as_date(dlg_input(message = "End date: YYYY-MM-DD")$res)
  
  # Time zone of sound files
  tz_files <- dlg_list(title = "Original file TZ",choices = c(grep("Etc/", OlsonNames(), value=TRUE)))$res
  # Local time zone
  tz_local <- dlg_list(title = "TZ for Figures",choices = c(grep("Etc/", OlsonNames(), value=TRUE)))$res
  
  return = list(site_id = site_id, 
                dep_id = dep_id, 
                start_dep_date = start_dep_date, 
                end_dep_date = end_dep_date, 
                tz_files = tz_files, 
                tz_local = tz_local)
  
}




inside_tables_to_hp <- function(ins_table, 
                                site_id = character()){
  # ves_og <- read_csv(choose.files(caption = "Choose Ship Detection Notes for given dep")) |>
  ves_og <- ins_table |>
    # change Behavior to factor and explicitly set all behaviors as levels
    mutate(Behavior = factor(Behavior, levels = c("CPA","CPAM","TA","TAM","TB","MANEUVER","TRANSIT")))
  
  
  
  #### Create Hourly Presence Table ####
  
  # Count instances of each behavior per date-hour
  hr_tally <- ves_og |>
    mutate(Begin_Date = as_date(Begin_Date, format = "%m/%d/%Y"))|>
    group_by(Begin_Date, Behavior, Begin_hour, SiteID) |>
    count() |>
    pivot_wider(names_from = Behavior,
                values_from = n,
                names_expand = TRUE,
                names_prefix = "ins_",
                values_fill = 0) |> # expand to include a column for all possible levels
    mutate(trans_inside = rowSums(across(c(ins_TB, ins_CPA, ins_TA, ins_TRANSIT))),
           man_inside = rowSums(across(c(ins_TAM, ins_CPAM, ins_MANEUVER))),
           total_inside = rowSums(across(c(trans_inside, man_inside))))
  
  # create a new df with all hours for the whole dataset
  date_range_dep <- seq.Date(from = min(hr_tally$Begin_Date), 
                             to = max(hr_tally$Begin_Date), by = "day") |>
    crossing(seq(0,23,1))
  
  # rename columns
  names(date_range_dep) <- c("Begin_Date","Begin_hour")
  
  # join 2 data frames together to add behavior tally
  hourly_pres <- date_range_dep |>
    left_join(hr_tally, by = c("Begin_Date","Begin_hour")) |>
    replace_na(list(Ves_counts_inside = 0,
                    SiteID = site_id))
  
  return = hourly_pres
}

# Taiki's stuff down here ####
library(lubridate)
library(readxl)
library(dplyr)
library(signal)

tryPm <- try(packageVersion('PAMmisc') >= '1.12.5', silent=TRUE)
if(!isTRUE(tryPm)) {
  pak::pkg_install('TaikiSan21/PAMmisc')
}
library(PAMmisc)

# x is either a POSIXct time or a wav file name
# tz is the actual timezone of these 
# result will be recast to UTC
#
tzAdjuster <- function(x, tz) {
  # step 1 convert files as if UTC time
  if(is.character(x)) {
    x <- fileToTime(x)
  }
  # special handling for UTC+-X style timezone (these are non-standard)
  if(grepl('UTC[\\+\\-][0-9]{1,2}', tz)) {
    tz <- gsub('UTC', '', tz)
    sign <- substr(tz, start=1, stop=1)
    nums <- as.numeric(substr(tz, start=2, stop=nchar(tz)))
    tz <- 'UTC'
  } else {
    sign <- '+'
    nums <- 0
  }
  # for proper Olson TZ's, need to force it then re-convert to UTC
  if(tz %in% OlsonNames()) {
    x <- force_tz(x, tzone=tz)
    x <- with_tz(x, tzone='UTC')
  }
  if(nums != 0) {
    x <- switch(sign,
                '+' = {
                  x - 3600 * nums
                },
                '-' = {
                  x + 3600 * nums
                }
    )
  }
  x
}

# Helper function to parse various styles of wav file names to times
# x is name of wav file (or many)
# result is POSIXct time in UTC
#
fileToTime <- function(x) {
  if(length(x) > 1) {
    result <- lapply(x, fileToTime)
    result <- as.POSIXct(unlist(result), origin='1970-01-01 00:00:00', tz='UTC')
    return(result)
  }
  x <- basename(x)
  x <- tools::file_path_sans_ext(x)
  format <- c('pamguard', 'pampal', 'soundtrap', 'sm3', 'icListens1', 'icListens2', 'AMAR')
  for(f in format) {
    switch(
      f,
      'pamguard' = {
        date <- gsub('.*([0-9]{8}_[0-9]{6}_[0-9]{3})$', '\\1', x)
        posix <- as.POSIXct(substr(date, 1, 15), tz = 'UTC', format = '%Y%m%d_%H%M%S')
        if(is.na(posix)) next
        millis <- as.numeric(substr(date, 17, 19)) / 1e3
        if(!is.na(posix)) {
          break
        }
      },
      'pampal' = {
        date <- gsub('.*([0-9]{14}_[0-9]{3})$', '\\1', x)
        posix <- as.POSIXct(substr(date, 1, 14), tz = 'UTC', format = '%Y%m%d%H%M%S')
        if(is.na(posix)) next
        millis <- as.numeric(substr(date, 16, 18)) / 1e3
        if(!is.na(posix)) {
          break
        }
      },
      'soundtrap' = {
        date <- gsub('.*\\.([0-9]{12})$', '\\1', x)
        posix <- as.POSIXct(date, format = '%y%m%d%H%M%S', tz='UTC')
        millis <- 0
        if(!is.na(posix)) {
          break
        }
      },
      'sm3' = {
        date <- gsub('.*\\_([0-9]{8}_[0-9]{6})Z?$', '\\1', x)
        posix <- as.POSIXct(date, format = '%Y%m%d_%H%M%S', tz='UTC')
        millis <- 0
        if(!is.na(posix)) {
          break
        }
      },
      'icListens1' = {
        date <- gsub('.*_([0-9]{8}-[0-9]{6})$', '\\1', x)
        posix <- as.POSIXct(date, format = '%Y%m%d-%H%M%S', tz='UTC')
        millis <- 0
        if(!is.na(posix)) {
          break
        }
      },
      'icListens2' = {
        date <- gsub('.*_([0-9]{6}-[0-9]{6})$', '\\1', x)
        posix <- as.POSIXct(date, format = '%y%m%d-%H%M%S', tz='UTC')
        millis <- 0
        if(!is.na(posix)) {
          break
        }
      },
      'AMAR' = {
        # 'AMAR668.9.20210823T231318Z.wav' example
        date <- gsub('.*([0-9]{8}T[0-9]{6}Z)$', '\\1', x)
        posix <- as.POSIXct(date, format='%Y%m%dT%H%M%SZ', tz='UTC')
        millis <- 0
        if(!is.na(posix)) {
          break
        }
      }
    )
  }
  posix + millis
}

# clarify which should be chooser options
arcToRaven <- function(arc=NULL, wav=NULL, wavTz='UTC', gpsTz='UTC', freq=250, duration=5, file=NULL, maxPoints=20) {
  if(is.character(arc)) {
    arc <- read_excel(arc)
  }
  if('DateTimeS' %in% names(arc)) {
    arc$UTC <- as.POSIXct(arc$DateTimeS, format='%Y-%m-%dT%H:%M:%S', tz='UTC')
  } else if(all(c('Date', 'Time') %in% names(arc))) {
    arc$UTC <- as.POSIXct(paste0(arc$Date, '_', arc$Time), 
                          format='%y%m%d_%H:%M:%S', tz='UTC')
  } else {
    stop('Column "DateTimeS" or columns "Date" and "Time" are not present')
  }
  arc$UTC <- tzAdjuster(arc$UTC, tz=gpsTz)
  wavTime <- tzAdjuster(wav, tz=wavTz)
  freq <- as.integer(freq)
  arc$OBJECTID <- as.integer(arc$OBJECTID)
  beginTime <-  as.numeric(difftime(arc$UTC, wavTime, units='secs'))
  if(any(beginTime < 0)) {
    stop('Some start times are negative, likely that timezone is incorrect')
  }
  result <- data.frame(
    Selection = arc$OBJECTID,
    View='Spectrogram',
    Channel=1L,
    'Begin Time (s)' = beginTime,
    'End Time (s)' = beginTime + duration,
    'Low Freq (Hz)' = freq,
    'High Freq (Hz)' = freq,
    'Begin Date' = NA,
    'Begin Clock Time' = NA,
    'End Clock Time' = NA,
    'Delta Time (s)' = NA,
    'Distance (m)' = arc$NEAR_DIST,
    FID_2 = arc$OBJECTID,
    # Site = arc$Comment,
    Lat = arc$Latitude,
    Long = arc$Longitude,
    # Segment = arc$Descript,
    check.names=FALSE
  )
  if('Comment' %in% names(arc)) {
    result$Site <- arc$Comment
  }
  if('Segment' %in% names(arc)) {
    result$Segment <- arc$Segment
  } else if('Descript' %in% names(arc)) {
    result$Segment <- arc$Descript
  }
  if('Segment' %in% names(result)) {
    result <- bind_rows(lapply(split(result, result$Segment), function(x) {
      if(nrow(x) < maxPoints) {
        warning('Segment ', x$Segment[1], ' had less than ', maxPoints,
                ' data points', call.=FALSE)
        keepSeq <- seq_len(nrow(x))
      } else {
        keepSeq <- seq(from=1, to=nrow(x), length.out=maxPoints)
      }
      x[keepSeq, ]
    }))
  } else {
    warning('"Segment" column not found in Arc data')
  }
  if(!is.null(file)) {
    write.table(result, file=file, sep='\t', row.names=FALSE, quote=FALSE)
  }
  result
}

# dir is wav files
# calibration is large negative number
# freqMin/Max are allowed peak ranges
calculatePeakLev <- function(dir, cal=NULL, freqMin=NULL, freqMax=NULL, progress=TRUE, file=NULL) {
  if(length(dir) == 1 &&
     dir.exists(dir)) {
    wavFiles <- list.files(dir, full.names=TRUE, pattern='\\.wav$')
  } else if all(file.exists(dir)) {
    wavFiles <- dir
  } else {
    warning('dir must be a directory or list of wav files')
    return(NULL)
  }
  WINDOW <- numeric(0)
  if(progress) {
    pb <- txtProgressBar(min=0, max=length(wavFiles), style=3)
    ix <- 0
  }
  result <- lapply(wavFiles, function(x) {
    data <- fastReadWave(x)[1,]
    data <- data - mean(data)
    data <- data #* 10^(-cal/20)
    sr <- attr(data, 'rate')
    if(length(WINDOW) < sr) {
      WINDOW <<- signal::hanning(sr)
    }
    if(is.null(freqMin)) {
      freqMin <<- 0
    }
    if(is.null(freqMax)) {
      freqMax <<- sr/2
    }
    nfft <- sr
    noverlap <- round(sr / 2)
    rms <- 20*log10(sqrt(mean(data^2))) - cal #believe this is correct cal
    hop <- nfft - noverlap
    welch <- pwelch(data, nfft=nfft, noverlap=noverlap, sr=sr, window=WINDOW, demean='long')
    peakLev <- max(abs(welch$spec[welch$freq >= freqMin & welch$freq <= freqMax]))
    peakFreq <- welch$freq[abs(welch$spec) == peakLev]
    peakLev <- 10*log10(peakLev)- cal
    result <- list(PeakFreq=peakFreq, PeakLev=peakLev, RMS=rms, Filename=basename(x))
    splitName <- strsplit(basename(x), '_')[[1]]
    
    if(length(splitName) != 3) {
      distance <- NA
    } else {
      distPart <- splitName[3]
      distPart <- gsub('\\.wav$', '', distPart)
      distance <- suppressWarnings(as.numeric(distPart))
    }
    result$Distance <- distance
    if(progress) {
      ix <<- ix + 1
      setTxtProgressBar(pb, value=ix)
    }
    result
  })
  result <- bind_rows(result)
  result$Calibration <- cal
  result$FreqMin <- freqMin
  result$FreqMax <- freqMax
  if(!is.null(file)) {
    if(!grepl('csv$', file)) {
      file <- paste0(file, '.csv')
    }
    write.csv(result, file=file, row.names=FALSE)
  }
  result
}
