### Results from Predict IO - webtool

#using the output of PredictIO to prepare the plots for our manuscript

library(ggplot2)
library(dplyr)
library(forestplot)

forest_plot <- read.csv("Downloads/forest-plot.csv")
forest_plot = forest_plot %>% mutate(mean=(X_95ci_low+X_95ci_high)/2) %>% arrange(mean) %>%
  mutate(name=paste(study, ' ', "(", primary_tissue, ')'))
forest_plot[18, ]=c("Custom", 'OS', 'COX', "Summary", rep(NA, 5),-1.13,-0.5,0, -0.8147, NA)
forest_plot$mean=as.numeric(forest_plot$mean)
forest_plot$X_95ci_low =as.numeric(forest_plot$X_95ci_low)
forest_plot$X_95ci_high =as.numeric(forest_plot$X_95ci_high)

volcano_plot <- read.csv("Downloads/volcano-plot.csv")
volcano_plot=volcano_plot %>% mutate(OR=(X_95ci_low+X_95ci_high)/2) %>%
  mutate(color=ifelse(pval>0.05 & OR>0, 'grey', 
                      ifelse(pval<= 0.05 & OR>0, 'red', 
                      ifelse(pval<=0.05 & OR<0, 'blue', 
                      ifelse(pval>0.05 & OR<0, 'grey', 'NS')))))
volcano_plot = volcano_plot %>% mutate(color=ifelse(signature=='Custom', 'purple', color))

ggplot(volcano_plot, aes(OR, -log10(pval)))+geom_jitter(aes(color=color), size=3)+
  geom_vline(xintercept = 0, linetype='dashed')+
  geom_hline(yintercept=-log10(0.05), linetype='dashed')+
  theme_minimal()+scale_color_manual(values =c('#2471A3', 'grey', '#673AB7','#B03A2E') )

forestplot(forest_plot, labeltext=name, 
           lower='X_95ci_low', upper='X_95ci_high', 
           boxsize=0.5,
           is.summary=c(rep(FALSE, 17), TRUE), 
           col= fpColors(summary = "darkorange"))


