# Purpose -----------------------------------------------------------------
# This script is used to process the seabird CTD files. See `seabirdProcessing.html` for a vignette of this workflow.

# Libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(segmented)
library(lubridate)
library(zoo)

# Reading in data ---------------------------------------------------------
# Specify the month, year, and total number of casts
# A cast is defined as a deployment of the Seabird into the water to collect data. Theoretically, there is one cast per station.

# Function to read in the seabird data
read_seabird <- function(year, month,
                         path = "U:\\LTM\\Bay Study\\SeaBird\\Data") {
  # Does the path exist?
  if (!file.exists(path)) stop("Make sure the file path exists.", call. = F)
  
  if (missing(year)) {
    message("Pick a year:")
    return(print(list.files(path)))
  }
  
  if (missing(month)) {
    message("Pick a month:")
    return(print(list.files(file.path(path, year))))
  }
  
  # List of files in the designated year and month
  files <- list.files(file.path(path, year, month), pattern = "*.cnv", full.names = TRUE)
  # Reading in all the files. Will return a list
  lapply(files, function(x) {
    # Pulling the cast name from the file name, which is the values between the month, year, and .cnv strings
    castName <- gsub(paste0(".*/", year, "/", month, "/", month, year, "(\\d+)\\.cnv"), "\\1", x)
    
    oce::read.ctd.sbe(x, station = castName)
  }) %>% 
    setNames(gsub(".*/(.*)\\.cnv", "\\1", files))
}

# Actually read it in
seabirdCtd <- read_seabird(2023, "January")
# Do this across ALL the casts in that file
seabirdData <- lapply(seabirdCtd, function(x) {
  slot(x, "data") %>% 
    data.frame() %>% 
    mutate(cast = x@metadata$station,
           startTime = force_tz(x@metadata$startTime, tzone = "America/Los_Angeles"),
           month = format(startTime, format = "%m"),
           # Create the depth bins that baystudy has historically used
           depthBin = cut(depth, 
                          breaks = c(0, seq(0.7, ceiling(max(depth)) + 0.2, by = 0.5)),
                          labels = seq(0.5, ceiling(max(depth)), by = 0.5)))
})

# Isolating the downcast --------------------------------------------------

# Fit the model, extract the breakpoints and various fit metrics
breakpointsSelect <- lapply(seabirdData, 
                            function(x) {
                              
                              x <<- x
                              lmMod <- lm(pressure ~ scan, data = x)
                              
                              # Change Kmax if needed
                              fitSegmented <- tryCatch(selgmented(lmMod, type = "bic", Kmax = 5),
                                                       error = function(cond) {
                                                         warning(cond, 
                                                                 call. = F)
                                                       })
                              
                              if (!isTRUE(is.character(fitSegmented))) {
                                rsme <- sqrt(mean((x$pressure - fitSegmented$fitted.values)^2))
                                breakpoints <- fitSegmented$psi[, 2]
                                
                                p <- function(y) {
                                  plot(fitSegmented)
                                  points(x[, c("scan", "pressure")])
                                  abline(v = breakpoints, col = "#FF4D00")
                                }
                                
                                # Return this
                                list(data = x,
                                     mod = fitSegmented,
                                     plot = function(){p(fitSegmented)},
                                     rsme = rsme,
                                     breakpoints = breakpoints)
                              } else {
                                NULL
                              }
                            }) %>% 
  setNames(names(seabirdData))

# Isolating the downcast: several filters
# 1. The max depth is the end of the downcast. If there are multiple scans with the same max value, the first is used
# 2. The beginning of the downcast is a scan that is at least 10 scans prior to the scan of the max depth AND the difference between the scan at the start of the downcast and the scan that occurs 15 scans later is < -1
# These scans can be interpreted as: 1) the start of the downcast occurs before the end of the downcast and right before the seabird start its descent; and 2) the end of the downcast occurs when the seabird has reached the deepest depth that it will go

downcast <- lapply(breakpointsSelect,
                              function(x) {
                                x$data %>% 
                                  mutate(breakpoints = ifelse(row_number() %in% round(x$breakpoints), 
                                                              scan, NA),
                                         maxDepth = ifelse(depth == max(depth, na.rm = T), scan, NA),
                                         fromMax = scan - max(maxDepth, na.rm = T),
                                         slope = depth - lead(depth, 15),
                                         bookEndDowncast = ifelse(!is.na(maxDepth) | (!is.na(breakpoints) & 
                                                                                        fromMax < -10 & slope < -1) &
                                                                    scan <= first(which(depth == max(depth, na.rm = T))),
                                                                  T, F)) %>% 
                                  filter(between(scan, first(which(bookEndDowncast == T)), last(which(bookEndDowncast == T))))
                              })

# Quickly look across the downcast to see if there are casts during which this procedure did not accurately ID the downcast
downcast %>% 
  bind_rows(.id = "cast") %>% 
  ggplot(aes(scan, pressure)) +
  geom_line() +
  facet_wrap(~cast, scales = "free") 

breakpointsSelect[[1]]$plot()

# Binding in the station for each cast ------------------------------------
# This is the datasheet that is entered manually of the start time of each seabird.
stationCast <- read.csv(file.path("seabird", "SeabirdLogJan-Mar2023.csv")) %>%
  group_by(Station) %>% 
  transmute(startTimeLog = as.POSIXct(paste0(Date, " ", ifelse(nchar(Time) == 3, paste0("0", Time), Time)), 
                                      format = "%m/%d/%Y %H%M"),
            month = format(as.Date(Date, format = "%m/%d/%Y"), format = "%m"),
            cast = Cast, Station, rep = 1:n(), Depth, Comments) %>% 
  right_join(downcast %>% 
               bind_rows(.id = "fileName") %>% 
               mutate(cast = as.integer(cast)),
             by = c("month", "cast")) %>% 
  group_by(cast = factor(cast), Station = factor(Station)) %>% 
  mutate(timeDifference = as.numeric(difftime(startTime, startTimeLog, units = "mins")))

# This can then be quickly plotted to visualize if the above matches are correct:
distinct(stationCast, cast, rep, Station, timeDifference) %>% 
  ggplot(aes(factor(Station), timeDifference, color = factor(rep, levels = 1:3))) +
  # geom_segment(aes(xend = factor(Station), 
  #                  yend = ifelse(timeDifference >= 0, -Inf, Inf)), 
  #              size = 1) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, size = 1, color = "#EBEBEB") +
  expand_limits(y = 0) +
  labs(x = "Station", y = "Time Difference (min)", color = "Rep") +
  theme_bw(base_size = 24) +
  # guides(color = guide_legend(override.aes = list(linetype = c(0, 0)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Detect "outlying" points during downcast --------------------------------
# Given how tight the readings are, can try a basic sd approach

betweenOutlier <- function(x, sd = 2) {
  !between(x, mean(x, na.rm = T) - 2 * sd(x, na.rm = T),
           mean(x, na.rm = T) + 2 * sd(x, na.rm = T))
}

downcast$January2023001 %>% 
  pivot_longer(c(temperature, salinity, conductivity, depth, pressure), 
               names_to = "variable", values_to = "values") %>% 
  group_by(variable = factor(variable, levels = c("conductivity", "salinity", "temperature", "depth", "pressure"))) %>% 
  mutate(Outliers = betweenOutlier(values)) %>% 
  ggplot(aes(scan, values, color = Outliers)) +
  geom_point(size = 3) +
  facet_wrap(~variable, scale = "free") +
  theme_classic(base_size = 24)

# Another approach is to use a rolling MAD, median absolute deviation from the median. For more information:
# https://stats.stackexchange.com/a/35612 and ?mad()

MAD <- function(x, window = 2, threshold = 3) {
  # ut = upper threshold; lt = lower threshold
  ut <- function(x) {m = median(x); m + threshold * median(abs(x - m))}
  lt <- function(x) {m = median(x); m - threshold * median(abs(x - m))}
  
  xZoo = zoo(x)
  
  df <- data.frame(x, 
                   lowerThreshold = rollapply(xZoo, window, lt, align = "right", fill = x[1]),
                   upperThreshold = rollapply(xZoo, window, ut, align = "right", fill = x[1]))
  
  df$outliers <- factor(((df$x < df$lowerThreshold) | (df$x > df$upperThreshold)),
                        levels = c(T, F))
  df
}

plot_outliers <- function(downcast, ...) {
  
  cols <- scales::hue_pal()(4)
  
  downcast %>% 
    pivot_longer(c(temperature, salinity, conductivity, depth, pressure), 
                 names_to = "variable", values_to = "values") %>% 
    group_by(variable = factor(variable, levels = c("conductivity", "salinity", "temperature", "depth", "pressure"))) %>% 
    mutate(Outliers = MAD(values, ...)) %>% 
    ggplot(aes(scan, values, color = Outliers$outliers)) +
    geom_point(size = 3) +
    geom_line(aes(y = Outliers$upperThreshold, color = "upperThreshold")) +
    geom_line(aes(y = Outliers$lowerThreshold, color = "lowerThreshold")) +
    facet_wrap(~variable, scale = "free") +
    labs(x = "Scan", color = "Outlier") +
    scale_color_manual(values = c("TRUE" = cols[1],
                                  "FALSE" = cols[3],
                                  "upperThreshold" = cols[2],
                                  "lowerThreshold" = cols[4])) +
    guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA, NA),
                                                    linetype = c(NA, NA, 1, 1)))) +
    theme_classic(base_size = 24)  
}

# So far, a window of 2 works the best...
plot_outliers(downcast$January2023001)

