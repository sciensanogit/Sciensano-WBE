############################################################################### #
# Aim ----
#| Prepare SARS-CoV-2, RSV, Influenza data for visualization
# Content:
#| Load package
#| Load data
#| Clean data
#| Normalization
#| Smoothening
#| Aggregation at national level
#| Export data to .csv, .xls, .rds
############################################################################### #

# load packages ----
# specify package location
# .libPaths( "//sciensano.be/fs/1150_EPIVG_EpiInfect/15_WBE/PROJECTS/CodeLibraryR/librairies/R/4.5.2" )

# select packages
pkgs <- c("dplyr", "tidyr", "zoo", "writexl", "ggplot2")
# install packages if not yet done
install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, FUN = library, character.only = TRUE))

# load data ----

# import correction factors
correction_factor <- read.csv("data-correction_factor.csv") %>%
  filter(grepl("COV_4.|PMMV_2.|RSV_1.|INF_1.", idprotocol)) %>%
  rename(labProtocolID = idprotocol, value = result)

# Belgian data are available here https://www.geo.be/catalog/details/9eec5acf-a2df-11ed-9952-186571a04de2?l=en
#| Metadata
#| siteName is the name of the treatment plant
#| collDTStart is the date of sampling
#| labName is the name of the lab analysing the sample
#| labProtocolID is the protocol used to analyse the dample
#| flowRate is the flow rate measured at the inlet of the treatment plant during sampling
#| popServ is the population covered by the treatment plant
#| measure is the target measured
#| value is the result

# sars-cov-2 data
df_sc <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantscovid&outputFormat=csv")

# load influenza data
df_inf <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantsinfluenza&outputFormat=csv")

# load RSV data
df_rsv <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantsrsv&outputFormat=csv")

# pmmv data
df_pmmv <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantspmmv&outputFormat=csv")

# join both
df <- df_sc %>%
  rbind(df_inf, df_rsv, df_pmmv)

# clean data
df <- df %>%
  select(siteName, collDTStart, labName, labProtocolID, flowRate, popServ, measure, value, quality) %>% 
  rename(date = collDTStart) %>% 
  mutate(date = as.Date(date))

# set and subset dates
date_start <- as.Date("2024-06-10")
date_graph_end <- as.Date("2026-08-28")
date_switch <- "2025-09-01" # method transition
date_reporting <- as.Date("2026-02-01", format = "%Y-%m-%d")

df <- df %>%
  filter(date > date_start & date < date_reporting)

# rename measures
# diplay existing measure
# unique(df$measure)
df[df$measure == "SARS-CoV-2 E gene", ]$measure <- "Egen"                                  
df[df$measure == "SARS-CoV-2 nucleocapsid gene, allele 2", ]$measure <- "N2gen"             
df[df$measure == "Pepper mild mottle virus capsid protein gene region", ]$measure <- "PMMV"                                  
df[df$measure == "Influenza virus type A", ]$measure <- "FluA"                                
df[df$measure == "Influenza virus type B", ]$measure <- "FluB"                       
df[df$measure == "Human respiratory syncytial virus type A", ]$measure <- "RSVA"                     
df[df$measure == "Human respiratory syncytial virus type B", ]$measure <- "RSVB"  

# translate siteName to english 
df[df$siteName == "Bruxelles-Sud", ]$siteName <- "Brussels-South"
df[df$siteName == "Bruxelles-Nord", ]$siteName <- "Brussels-North"

## Apply correction factors ----
# create df1 to keep concentrations before the date at which the analytical method was modified
# subset based on labProtocolID used betwen date_start and date_end
df1 <- df %>%
  filter(date < date_switch & grepl("_COV_4.|_PMMV_2.|_RSV_1.|_INF_1.", labProtocolID))

# correct the concentrations before switch date with correction factors
df1 <- df1 %>%
  left_join(
    correction_factor %>%
      select(labProtocolID, measure, value) %>%
      rename(factor = value),
    by = c("labProtocolID", "measure")) %>%
  mutate(value = value * factor) %>%
  select(-factor)

# create df2 with concentrations after the switch date
# subset based on labProtocolID used betwen date_start and date_end
df2 <- df %>%
  filter(date >= date_switch & grepl("_COV_5.|_PMMV_3.|_RSV_2.|_INF_2.", labProtocolID))

# append the two;
df <- rbind(df1, df2)

# rename measures
df[df$measure == "Egen", ]$measure <- "SARS-CoV-2"
df[df$measure == "N2gen", ]$measure <- "SARS-CoV-2"
df[df$measure == "FluA", ]$measure <- "Influenza"
df[df$measure == "FluB", ]$measure <- "Influenza"
df[df$measure == "RSVA", ]$measure <- "RSV"
df[df$measure == "RSVB", ]$measure <- "RSV"

# apply LOQ provided by the lab
df[df$measure == "PMMV" & df$value < 125, ]$value <- NA

# remove outliers
df[df$quality == "Quality concerns", ]$value <- NA

# normalization ----
# compute mean of replicated analysis of each measure
df <- df %>%
  select(date, siteName, labName, flowRate, popServ, measure, value) %>%
  group_by(date, siteName, labName, flowRate, popServ, measure) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>% ungroup()

# pivot to have SARS, Influenza, RSV and PMMV as variables
df <- df %>%
  pivot_wider(names_from = measure, values_from = value)

# pivot again to have SARS, Influenza, RSV in the "measure" variable
df <- df %>%
  select(date, siteName, labName, flowRate, popServ, PMMV, 'SARS-CoV-2', Influenza, RSV) %>%
  pivot_longer(cols = 'SARS-CoV-2':RSV, names_to = "measure", values_to = "value")

# compute viral load (value_load), viral ratio (value_ratio)
df <- df %>%
  mutate(value_load = value*flowRate*24*1000000/popServ*100000,
         value_pmmv = value/PMMV)

# save
df_site_raw <- df

# smoothening ----
# compute the linear extrapolation data
df <- df_site_raw %>%
  group_by(measure, siteName) %>%
  complete(date = seq(min(date), max(date), "day")) %>%
  mutate(value_avg14d_past = na.approx(value, maxgap = 14, na.rm = FALSE),
         value_load_avg14d_past = na.approx(value_load, maxgap = 14, na.rm = FALSE),
         value_pmmv_avg14d_past = na.approx(value_pmmv, maxgap = 14, na.rm = FALSE))

# compute moving average on past 14 days
df <- df %>%
  group_by(measure, siteName) %>%
  mutate(across(value_avg14d_past:value_pmmv_avg14d_past,
                ~ rollmean(.x, k = 14, fill = NA, na.rm = TRUE, align = "right")))

# save
df_site <- df

# national level ----
## aggregation ----
# compute weighted mean with factor being the population served by each site
df <- df_site_raw %>%
  select(measure, date, popServ, value, value_load, value_pmmv) %>%
  mutate(siteName = "Belgium") %>%
  group_by(measure, siteName, date) %>%
  summarise(across(value:value_pmmv, ~ weighted.mean(.x, popServ, na.rm=TRUE))) %>% ungroup()

## smoothening ----
# linear extrapolation data
df <- df %>%
  group_by(measure, siteName) %>%
  complete(date = seq(min(date), max(date), "day")) %>%
  mutate(value_avg14d_past = na.approx(value, maxgap = 14, na.rm = FALSE),
         value_load_avg14d_past = na.approx(value_load, maxgap = 14, na.rm = FALSE),
         value_pmmv_avg14d_past = na.approx(value_pmmv, maxgap = 14, na.rm = FALSE))

# moving average on past 14 days
df <- df %>%
  group_by(measure, siteName) %>%
  mutate(across(value_avg14d_past:value_pmmv_avg14d_past,
                ~ rollmean(.x, k = 14, fill = NA, na.rm = TRUE, align = "right")))

# save
df_nation <- df

# export data ----
# create folder if not existing
dir.create("./data", showWarnings = F)

# export as rds
saveRDS(list(df_site_raw, df_site, df_nation),
        file = "./data/Belgium_export.rds")

# visual ----
# create folder if not existing
dir.create("./plot", showWarnings = F)

plot <- df_nation %>%
  filter(measure == "RSV") %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = value_pmmv), na.rm = T) +
  geom_line(aes(y = value_pmmv_avg14d_past), na.rm = T)

plot

# save
ggsave(file="./plot/Graph_be-RSV-nation-viral_ratio.png",
       plot, width = 21, height = 12, dpi = 200)
