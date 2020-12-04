library(reshape2) 
library(dplyr) 
library(raster)
library(sp) 
library(tidyverse)
library(prism)

setwd("~/github/2020-pn-tmp/")
options(prism.path = "~/github/2020-pn-tmp/prism_data")

info <- read_csv("samples.csv") %>%
  select(site, ava_id, ava, lon, lat) %>%
  filter(!is.na(lon)) %>%
  distinct()

vins_spdf <- SpatialPointsDataFrame(coords = info[ , c('lon','lat')], data = info, 
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))

# ppt ---------------------------------------------------------------------

get_prism_monthlys(type = 'ppt', years=c(2017, 2019), mon = 4:9, keepZip = TRUE)

rs_ppt <- ls_prism_data()[1:nrow(ls_prism_data()), 1][grepl(pattern = "ppt", x = ls_prism_data()[1:nrow(ls_prism_data()), 1])]
rs_ppt <- prism_stack(rs_ppt)

vins_ppt <- raster::extract(rs_ppt, vins_spdf, na.rm=TRUE, sp=TRUE) %>% 
  as_tibble %>%
  dplyr::select(-lon.1, -lat.1)

colnames(vins_ppt) <- gsub("_bil", "", colnames(vins_ppt))
colnames(vins_ppt) <- gsub("PRISM_ppt", "", colnames(vins_ppt))
colnames(vins_ppt) <- gsub("_stable_4kmM3_", "", colnames(vins_ppt))
colnames(vins_ppt)
ppt <- vins_ppt %>% 
  pivot_longer(cols = "201701":"201909", names_to = "date", values_to = "ppt") %>%
  mutate(year = gsub(".{2}$", "", date)) %>%
  mutate(month = gsub("^201[79]", "", date)) %>%
  mutate(month = gsub("^0", "", month)) %>%
  mutate(month = as.numeric(month)) %>%
  filter(month %in% c(1:9))

ggplot(ppt, aes(x = year, y = ppt)) +
  geom_boxplot() +
  facet_wrap(~month) + 
  theme_minimal()

ppt_sum <- ppt %>%
  group_by(site, year) %>%
  summarize(total_ppt = sum(ppt))


# write.csv(precip_data, "precip_data_sept.csv")

# gdd ---------------------------------------------------------------------

get_prism_dailys(type = 'tmin', minDate = "2017-04-01", maxDate = "2017-09-30", keepZip = F)
get_prism_dailys(type = 'tmin', minDate = "2019-04-01", maxDate = "2019-09-30", keepZip = F)

rs_tmin <- ls_prism_data()[1:nrow(ls_prism_data()), 1][grepl(pattern = "tmin", x = ls_prism_data()[1:nrow(ls_prism_data()), 1])]
rs_tmin <- prism_stack(rs_tmin)
vins_tmin <- raster::extract(rs_tmin, vins_spdf, na.rm=TRUE, sp=TRUE) %>% 
  as_tibble %>%
  dplyr::select(-lon.1, -lat.1)

colnames(vins_tmin) <- gsub("_bil", "", colnames(vins_tmin))
colnames(vins_tmin) <- gsub("PRISM_tmin_stable_4kmD2_", "", colnames(vins_tmin))
colnames(vins_tmin)
tmin <- vins_tmin %>% 
  pivot_longer(cols = "20170401":"20190930", names_to = "date", values_to = "tmin") %>%
  mutate(year = gsub(".{4}$", "", date)) %>%
  mutate(month = gsub("^201[79]", "", date)) %>%
  mutate(month = gsub("^0", "", month)) %>%
  mutate(month = gsub(".{2}$", "", month)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = gsub("^.{6}", "", date)) %>%
  mutate(day = gsub("^0", "", day)) %>%
  mutate(day = as.numeric(day))

# tmax
get_prism_dailys(type = 'tmax', minDate = "2017-04-01", maxDate = "2017-09-30", keepZip = F)
get_prism_dailys(type = 'tmax', minDate = "2019-04-01", maxDate = "2019-09-30", keepZip = F)

rs_tmax <- ls_prism_data()[1:nrow(ls_prism_data()), 1][grepl(pattern = "tmax", x = ls_prism_data()[1:nrow(ls_prism_data()), 1])]
rs_tmax <- prism_stack(rs_tmax)
vins_tmax <- raster::extract(rs_tmax, vins_spdf, na.rm=TRUE, sp=TRUE) %>% 
  as_tibble %>%
  dplyr::select(-lon.1, -lat.1)

colnames(vins_tmax) <- gsub("_bil", "", colnames(vins_tmax))
colnames(vins_tmax) <- gsub("PRISM_tmax_stable_4kmD2_", "", colnames(vins_tmax))
colnames(vins_tmax)
tmax <- vins_tmax %>% 
  pivot_longer(cols = "20170401":"20190930", names_to = "date", values_to = "tmax") %>%
  mutate(year = gsub(".{4}$", "", date)) %>%
  mutate(month = gsub("^201[79]", "", date)) %>%
  mutate(month = gsub("^0", "", month)) %>%
  mutate(month = gsub(".{2}$", "", month)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(day = gsub("^.{6}", "", date)) %>%
  mutate(day = gsub("^0", "", day)) %>%
  mutate(day = as.numeric(day))

# GDD equation: 
# https://www.evineyardapp.com/blog/2017/03/01/why-the-need-to-calculate-growing-degree-days-in-vineyard/
# http://oregonviticulture.net/gdd/gdd.html
# http://wine.wsu.edu/extension/weather/growing-degree-days/
gdd_c <- left_join(tmin, tmax) %>%
  mutate(gdd = ((tmin + tmax) / 2) - 10) %>%
  group_by(site, year) %>%
  summarize(total_gdd_c = sum(gdd))

gdd_f <- left_join(tmin, tmax) %>%
  mutate(tmin_f = (9/5) * tmin + 32) %>%
  mutate(tmax_f = (9/5) * tmax + 32) %>%
  mutate(gdd_f = ((tmin_f + tmax_f) / 2) - 50) %>%
  group_by(site, year) %>%
  summarize(total_gdd_f = sum(gdd_f))


gdd <- left_join(gdd_c, gdd_f)


# combine -----------------------------------------------------------------

prism <- full_join(gdd, ppt_sum) %>%
  mutate(year = as.numeric(year))
View(prism)

info_all <- read_csv("samples.csv") %>%
  full_join(prism)
write_csv(info_all, "samples.csv")
