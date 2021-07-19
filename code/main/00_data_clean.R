################################################################################
# Description: Join 2016-18 diagnoses with village locations
################################################################################
################################################################################

diag_raw <- read.csv(here::here("KAMIS/KAMIS geotagged","village level data-Bi+Jh.csv"),header = T)


################################################################################
# Output cleaned datasets
################################################################################

rawdir <- "./KAMIS/Raw"
outdir <- "./KAMIS/Clean/"

version <- "current"
start <- ymd("2013-01-01")
end <- ymd("2020-12-31")

agg <- aggregate_diag(start, end, state = "BH", log = here::here(outdir,paste0("data_aggregation_log_", start, "-" ,end,".txt")))

saveRDS(agg$VL, file = here::here(outdir,paste0("vl_",start,"_",end,".rds")))
saveRDS(agg$PKDL, file = here::here(outdir,paste0("pkdl_",start,"_",end,".rds")))

png(filename = here::here(outdir,paste0("KAMIS case count ", start, "-" ,end,".png")), height = 600, width = 800, res = 150)
print(agg$plot.totals)
dev.off()

agg$VL %>%
  group_by(agg_date) %>%
  summarise(n = sum(n, na.rm = TRUE), pop = sum(pop)) %>%
  ggplot(aes(agg_date, n*1e4/pop)) +
  geom_col(fill = "grey") + 
  xlim(c(ymd("2018-12-01"), end)) +
  geom_vline(xintercept = ymd("2020-03-01"), col = "red", lty = "dashed") +
  labs(x = "Month",y = "Incidence/10,000", title = "VL diagnoses Jan 2019-Jan 2021",
       subtitle = "Declaration of COVID-19 pandemic highlighted in red") +
  theme_minimal()




# Village incidence dataset with GPS for affected villages
# [original]
village_raw <- read.csv(here::here("KAMIS/KAMIS geotagged","village level data-Bi+Jh.csv"),header = T)