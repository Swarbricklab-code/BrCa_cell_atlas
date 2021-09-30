# =========================
# Daniel Roden
# d.roden@garvan.org.au
# Garvan Institute
# -------------------------
# VERSION: 1.0
# =========================
# DESCRIPTION:
# ============
# Utility/Wrapper functions for survival analysis
#
# ------
# NOTES:
# ------
#
# --------
# TODO:
# --------
#
# ======================================================================
# ===============
# LIBRARIES
# ===============
library(survival)
library(survminer)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ===============
# FUNCTIONS
# ===============

run_survival_analysis_logrank <- function(clinical_data, risk.table = F, color = NULL, conf.int = FALSE, end_time = Inf, palette, title = "", exclude_samples=FALSE) 
{
  if(end_time != Inf)
  {
    if(exclude_samples) {
      message(paste0("filtering OS to end at time = ", end_time))
      message(paste0("BEFORE: ", nrow(clinical_data)))
      clinical_data <- clinical_data %>% dplyr::filter(OS <= end_time)
      message(paste0("AFTER: ", nrow(clinical_data)))
    } else {
      
    }
  } else {
    message(paste0("Using all OS events"))
    end_time <- "all"
  }
  surv <- survfit(Surv(OS, EVENT) ~ Ecotype, data = clinical_data)
  
  #file_surv_obj <- paste0(dir_out, "survfit.logrank.rds")
  #saveRDS(surv, file = file_surv_obj)
  theme <- theme_bw() +
    theme(panel.grid = element_blank())
  
  surv_plot <- ggsurvplot(surv, 
                          data = clinical_data,  
                          title = title, 
                          font.title = c(22, "bold", "darkblue"), 
                          pval = TRUE, 
                          conf.int = conf.int,
                          legend="top", 
                          risk.table = risk.table,
                          pval.size = 9,
                          legend.title = "Ecotype",
                          legend.labs = sort(unique(clinical_data$Ecotype)),
                          font.legend = c(25, "bold"),
                          font.x = c(22, "bold"),
                          font.y = c(22, "bold"),
                          font.tickslab = c(21),
                          size = 1.25,
                          palette = palette,
                          ggtheme = theme,
                          xlab = "Overall Survival (Months)", 
                          xlim = c(0, end_time),
                          break.x.by = 60)
  return(surv_plot)
}
