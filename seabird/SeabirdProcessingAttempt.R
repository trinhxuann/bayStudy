library(oce)
# library(ocedata)
library(tidyverse)
library(lubridate)
library(dplyr)

# load data ---------------------------------------------------------------

## ID month and year you are working with
month <- "January"
year <- 2023
# cast per station, i.e., 52 stations ideally. So there should be 52 files
casts <- 52

## determine path
pathname <- paste("U:\\LTM\\Bay Study\\SeaBird\\Data\\", year , "\\", month, sep = "")

## get file names
files<- list.files(path = pathname, pattern = "*.cnv", full.names=TRUE)

##get all the .cnv file from the directory
  
  ctd = list()

for (i in 1:casts){
  
  ctd[[i]] = read.ctd(files[i])%>%
    ctdTrim(method = "downcast")
    #%>%    ctdDecimate(p = 0.2)
}
  ## name casts; need to know how many casts, change as needed
  names(ctd) <- paste("raw_", 1:casts, sep = "")


# get data from list ------------------------------------------------------

  
getdata <-  lapply(ctd, slot, name = "data")
profiles <- purrr::map_df(getdata, data.frame, .id = 'name')  %>% 
  filter(scan > 100) %>% 
  filter(descentRate > 0.01) %>% 
  mutate(survey = month, year = year)
 

# get start time from metadata to help match station to cast -------------

## extracting time
  a <- ctd %>% 
    map(slot, name = "metadata")
   b <- lapply(a, function (x) `[`(x, c('startTime'))) 
   df <- data.frame(matrix(unlist(b), nrow=casts, byrow=TRUE),stringsAsFactors=FALSE)
   df <- df %>% rename(x = 1)
   df <- as.POSIXct(df$x, origin ="1970-01-01 00:00:00", tz = "UTC")
   df <- as.data.frame(df)

## make dataframe
starttimes <- df %>% 
  mutate(name = paste("raw_", 1:casts, sep = "")) %>% 
  rename(startTime = df) %>% 
  select(name, startTime) %>% 
  mutate(startTime = as.character(startTime))

## use shiny to enter data to table

# Start time should match the start time of the tow per station.

## save csv to match station using seabird log, when matched re-save with same name
path2 <- paste(pathname, "\\StartTimesForMatch_", month, year, ".csv", sep = "")
write_csv(starttimes, file = path2)

## read in matched cast/station

matched <- read_csv(path2) %>% 
  select(name, startTime)

## add station to profile data
## LEFT OFF HERE TRYING TO DO ROUNDED DEPTH AND DEPTH BIN
cleanprofiles <- profiles %>% 
  left_join(matched, by = "name")
  # filter(!station == 999) 
#%>% 
  # mutate(roundedDepth = cut(depth, breaks = c(0, 0.7, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7), 
  #                           #labels = c("0.5", "1", "1.5", "2", "2.5", "3", "3.5")
  #                           labels = F)) %>% 
  # mutate(case_when())
path3 <- paste(pathname, "\\CleanProfiles_", month, year, ".csv", sep = "")
write_csv(cleanprofiles, path3)


## plot

cleanprofiles %>% 
  ggplot(aes(x=temperature, y=depth)) +
  geom_point() +
  scale_y_reverse() +
  facet_wrap(~name)





# Archive -----------------------------------------------------------------

## Another way to load profiles

files<- dir("U:\\LTM\\Bay Study\\SeaBird\\Data\\2022\\October",pattern = "*.cnv") #get all the .cnv file from the directory

casts<- list() #make a list dataframe
for (ifile in 1:length(files)) {
  casts[[ifile]]<- read.ctd.sbe(files[ifile]) 
}
casts <- read.oce(dir("U:\\LTM\\Bay Study\\SeaBird\\Data\\2022\\October",pattern = "*.cnv"))


filenames <- list.files(path = "U:\\LTM\\Bay Study\\SeaBird\\Data\\2022\\October", pattern="*.cnv", full.names=TRUE)
casts <- lapply(filenames, read.ctd.sbe)
names(casts) <- paste("raw_", 1:53, sep = "")

casts$raw_1@data%>%as.data.frame()%>%head()
summary(casts$raw_1)

plot(casts$raw_1)

casts_trim <- casts %>% 
  lapply(ctdTrim)

attributes(casts_trim)

casts_tables <- lapply(casts_trim$, as.tibble)

plot(casts_clean$raw_1)
