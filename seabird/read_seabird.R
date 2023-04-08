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

what <- read_seabird(2023, "January")

what2 <- lapply(what, function(x) {
  slot(x, "data") %>% 
    data.frame() %>% 
    mutate(cast = x@metadata$station,
           startTime = force_tz(x@metadata$startTime, tzone = "America/Los_Angeles"),
           month = format(startTime, format = "%m"),
           depthBin = cut(depth, 
                          breaks = c(0, seq(0.7, ceiling(max(depth)) + 0.2, by = 0.5)),
                          labels = seq(0.5, ceiling(max(depth)), by = 0.5)))
})

# Breakpoint analysis -----------------------------------------------------
# Used to find just the downcast
# First, identify each change point MANUALLY
manualChangePoints <- read.csv(file.path("seabird", "manualStartPoint.csv"))

# Now to try different methods
library(changepoint)

changepoints <- lapply(list(cpt.mean, cpt.meanvar, cpt.var),
       function(x) {
         lapply(what2, function(data) {
           # First pressure
           model <- x(data$pressure, method = "PELT")
           p <- function(y) plot(y)
           changePoints <- model@cpts
           
           # Then depth
           modelDepth <- x(data$depth, method = "PELT")
           changePointsDepth <- modelDepth@cpts
           
           list(model = model,
                plot = function(){p(model)},
                cpts = changePoints,
                modelDepth = modelDepth,
                plotDepth = function(){p(modelDepth)},
                cptsDepth = changePointsDepth)
       }) %>% 
           setNames(names(what2))}) %>% 
  setNames(c("mean", "meanVar", "var"))

# Appears that mean works the best
tibble(cast = names(what2), 
       cptsMean = sapply(changepoints$mean, function(x) x$cpts[1]),
       cptsMeanVar = sapply(changepoints$meanVar, function(x) x$cpts[1]),
       cptsMedian = sapply(changepoints$var, function(x) x$cpts[1])) %>% 
  bind_cols(select(manualChangePoints, startScan)) %>% 
  mutate(cptsMeanDifference = cptsMean - startScan,
         cptsMeanVarDifference = cptsMeanVar - startScan,
         cptsMedianDifference = cptsMedian - startScan) %>% 
  pivot_longer(contains("Difference"), names_to = "method", values_to = "differenceToManual") %>% 
  ggplot(aes(cast, differenceToManual, fill = method)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tibble(cast = names(what2), 
       numberCptsMean = sapply(changepoints$mean, function(x) length(x$cpts)),
       numberCptsMeanVar = sapply(changepoints$meanVar, function(x) length(x$cpts)),
       numberCptsMedian = sapply(changepoints$var, function(x) length(x$cpts))) %>% 
  pivot_longer(-cast, names_to = "method", values_to = "numberOfCpts") %>% 
  ggplot(aes(cast, numberOfCpts, fill = method)) +
  geom_col(position = "dodge")
# The number of change points may not be as informative here.

tibble(cast = names(what2), 
       cptsMean = sapply(changepoints$mean, function(x) x$cpts[1]),
       cptsMeanVar = sapply(changepoints$meanVar, function(x) x$cpts[1]),
       cptsMedian = sapply(changepoints$var, function(x) x$cpts[1]),
       # Depth
       cptsMeanDepth = sapply(changepoints$mean, function(x) x$cptsDepth[1]),
       cptsMeanVarDepth = sapply(changepoints$meanVar, function(x) x$cptsDepth[1]),
       cptsMedianDepth = sapply(changepoints$var, function(x) x$cptsDepth[1])) %>% 
  bind_cols(select(manualChangePoints, startScan)) %>% 
  mutate(cptsMeanDifference = cptsMean - startScan,
         cptsMeanVarDifference = cptsMeanVar - startScan,
         cptsMedianDifference = cptsMedian - startScan,
         # Depth
         cptsMeanDifferenceDepth = cptsMeanDepth - startScan,
         cptsMeanVarDifferenceDepth = cptsMeanVarDepth - startScan,
         cptsMedianDifferenceDepth = cptsMedianDepth - startScan) %>% 
  pivot_longer(contains("Difference"), names_to = "method", values_to = "differenceToManual") %>% 
  ggplot(aes(cast, differenceToManual, fill = method)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Answer essentially the same between depth and pressure
tibble(cast = names(what2), 
       cptsMean = sapply(changepoints$mean, function(x) x$cpts[1]),
       # Depth
       cptsMeanDepth = sapply(changepoints$mean, function(x) x$cptsDepth[1])) %>% 
  bind_cols(select(manualChangePoints, startScan)) %>% 
  mutate(cptsMeanDifference = cptsMean - startScan,
         # Depth
         cptsMeanDifferenceDepth = cptsMeanDepth - startScan) %>% 
  pivot_longer(contains("Difference"), names_to = "method", values_to = "differenceToManual") %>% 
  ggplot(aes(cast, differenceToManual, fill = method)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Breakpoint analysis -----------------------------------------------------
library(segmented)
# fitTry <- lm(pressure ~ scan, data = what2$January2023029)
# fitSegmented <- selgmented(fitTry)

# Approach will be to try a bunch of different number of max breakpoints and see which one fits the best
# Then, to find the correct one, applies several filters, specifically: in order
# 1. Find max pressure
# 2. Find closest breakpoint that occurs before max pressure is acheived
# Additional filters if needed
# 3. Filter for descentRate > 0.1
# 4. Find closest to a pressure of 2

# breakpoints <- lapply(3:10, 
#                       function(x) {
#                         lapply(what2, function(data) {
#                           
#                           cast <- unique(data$cast)
#                           cat("Working on", x, "breakpoint(s) for cast", cast, "\n")
#                           lmMod <- lm(pressure ~ scan, data = data)
#                           
#                           fitSegmented <- tryCatch(segmented(lmMod, sig.Z = ~scan, npsi = x),
#                                                    error = function(cond) {
#                                                      warning(paste0(cond, " for max breakpoint ", x, " for cast ", cast), 
#                                                              call. = F)
#                                                    })
#                           
#                           if (!isTRUE(is.character(fitSegmented))) {
#                             rsme <- sqrt(mean((data$pressure - fitSegmented$fitted.values)^2))
#                             breakpoints <- fitSegmented$psi
#                             p <- function(y) {
#                               plot(y)
#                               points(data[, c("scan", "pressure")])
#                               points.segmented(y)
#                               lines.segmented(y)
#                             }
#                             
#                             list(mod = fitSegmented,
#                                  plot = function(){p(fitSegmented)},
#                                  rsme = rsme,
#                                  npsi = x,
#                                  breakpoints = breakpoints)
#                           } else {
#                             NULL
#                           }
#                         })
#                       }) %>% 
#   setNames(c(paste0("npsi_", 3:10)))
# 
# bestBreaks <- lapply(breakpoints, 
#                      function(x) {
#                        lapply(x, function(y) {
#                          nBreak <- ifelse(is.null(y$npsi), "Failed", y$npsi)
#                          tibble("rmse_{nBreak}" := y$rsme,
#                                 "breakpoints_{nBreak}" := data.frame(y$breakpoints))
#                        }) %>% 
#                          bind_rows(.id = "cast")
#                      }) %>% 
#   reduce(full_join, by = "cast") %>% 
#   pivot_longer(-cast, names_to = "bestBreak", values_to = "scan") %>% 
#   group_by(cast) %>% 
#   slice_min(scan)

# How to extract the breakpoints:
breakpointsSelect <- lapply(what2, 
                            function(x) {
                              
                              x <<- x
                              lmMod <- lm(pressure ~ scan, data = x)
                              
                              fitSegmented <- tryCatch(selgmented(lmMod, type = "bic", Kmax = 5),
                                                       error = function(cond) {
                                                         warning(cond, 
                                                                 call. = F)
                                                       })
                              
                              rsme <- sqrt(mean((x$pressure - fitSegmented$fitted.values)^2))
                              breakpoints <- fitSegmented$psi[, 2]
                              
                              p <- function(y) {
                                plot(fitSegmented)
                                points(x[, c("scan", "pressure")])
                                abline(v = breakpoints, col = "#FF4D00")
                              }
                              
                              list(data = x,
                                   mod = fitSegmented,
                                   rsme = rsme,
                                   breakpoints = breakpoints,
                                   plot = function(){p(fitSegmented)})
                            }) %>% 
  setNames(names(what2))

# breakpointsSelect$January2023021$data %>% 
#   mutate(breakpoints = ifelse(row_number() %in% round(breakpointsSelect$January2023021$breakpoints), 
#                               scan, NA),
#          maxDepth = ifelse(depth == max(depth, na.rm = T), scan, NA)) %>% 
#   mutate(fromMax = scan - max(maxDepth, na.rm = T),
#          slope = depth - lead(depth, 15),
#          bookEndDowncast = ifelse(!is.na(maxDepth) | (!is.na(breakpoints) & fromMax < -10 & slope < -1) &
#                            scan <= first(which(depth == max(depth, na.rm = T))),
#                            T, F)) %>% 
#   filter(between(scan, first(which(bookEndDowncast == T)), last(which(bookEndDowncast == T))))
# 
# breakpointsSelect$January2023046$plot()

downcastBreakpoints <- lapply(breakpointsSelect,
       function(x) {
         x$data %>% 
           mutate(breakpoints = ifelse(row_number() %in% round(x$breakpoints), 
                                           scan, NA),
                      maxDepth = ifelse(depth == max(depth, na.rm = T), scan, NA),
                  fromMax = scan - max(maxDepth, na.rm = T),
                  slope = depth - lead(depth, 15),
                  bookEndDowncast = ifelse(!is.na(maxDepth) | (!is.na(breakpoints) & fromMax < -10 & slope < -1) &
                                             scan <= first(which(depth == max(depth, na.rm = T))),
                                           T, F)) %>% 
           filter(between(scan, first(which(bookEndDowncast == T)), last(which(bookEndDowncast == T))))
       })

downcastBreakpoints %>% 
  bind_rows(.id = "cast") %>% 
  ggplot(aes(scan, pressure)) +
  geom_line() +
  facet_wrap(~cast, scales = "free") 

# Now, to bind back to the station
what3 <- read.csv(file.path("seabird", "SeabirdLogJan-Mar2023.csv")) %>%
  group_by(Station) %>% 
  transmute(startTimeLog = as.POSIXct(paste0(Date, " ", ifelse(nchar(Time) == 3, paste0("0", Time), Time)), 
                                   format = "%m/%d/%Y %H%M"),
            month = format(as.Date(Date, format = "%m/%d/%Y"), format = "%m"),
            cast = Cast, Station, rep = 1:n(), Depth, Comments) %>% 
  # IMPORTANT: Remove casts that are not valid here
  # filter(cast != 22) %>% 
  right_join(downcastBreakpoints %>% 
               bind_rows(.id = "fileName") %>% 
               mutate(cast = as.integer(cast)),
             by = c("month", "cast")) %>% 
  group_by(cast = factor(cast), Station = factor(Station)) %>% 
  mutate(timeDifference = as.numeric(difftime(startTime, startTimeLog, units = "mins")))

distinct(what3, cast, rep, Station, timeDifference) %>% 
  ggplot(aes(factor(Station), timeDifference, color = factor(rep, levels = 1:3))) +
  geom_point(size = 5) +
  scale_color_manual(values = c("1" = "#FF4D00", "2" = "#FFFFB8")) +
  expand_limits(y = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Station", y = "Time Difference (min)", color = "Rep")
  

# Correlation between pressure and depth basically the same...?
what2 %>% 
  bind_rows() %>% 
  # ggplot(aes(depth, pressure)) +
  # geom_line()
  select(depth, pressure) %>% 
  {cor.test(.$depth, .$pressure)}


