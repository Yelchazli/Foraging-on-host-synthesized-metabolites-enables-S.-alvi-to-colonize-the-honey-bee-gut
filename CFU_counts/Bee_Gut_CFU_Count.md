# Generate plots for folony forming units (CFU’s) plated from bee guts. Figure 1A, Figure

## environment setup

### Install and load libraries (run coade if missing libraries)

#### These are the required package to reproduce this analysis.

``` r
#required_packages = c("readxl",  "tidyverse",    "RColorBrewer" , "ggbeeswarm","ggpubr", "reshape2" , "gridExtra", "stringr","EnvStats", "multcompView", "stats", "ggh4x")
#need_install = required_packages[!(required_packages) %in% installed.packages()]

#if (length(need_install) > 0){
#  install.packages((need_install))
#}

#lapply(required_packages, require, character.only=T)

#dir_script = dirname(rstudioapi::getSourceEditorContext()$path)

#setwd(dir_script)
```

#Plot CFU counts across diet and colonizations groups

``` r
library(readxl)  
library(tidyverse)
library(RColorBrewer) 
library(ggbeeswarm)
library(ggpubr) 
library(reshape2)  
library(gridExtra) 
library(stringr)
library(EnvStats)
library(multcompView) 
library(stats)
#library(ggh4x)


CFUs <- data.frame(read.csv("All_CFU_Count.csv")) 
CFUs$Colonization <- factor(CFUs$Colonization, levels=c("MD", "Sn", "Gi", "Sn_Gi"))
CFUs$Diet <- factor(CFUs$Diet, levels=c("SW", "P"))
CFUs$Condition <- factor(CFUs$Condition, levels=c("MD_SW", "MD_P", "Sn_SW", "Sn_P","Gi_SW", "Gi_P","Sn_Gi_SW", "Sn_Gi_P"))
CFUs$Experiment <- factor(CFUs$Experiment)

cfu_anova <- aov(log10(CFU)~ Condition, data=CFUs, na.action = na.omit)
tHSD_cfu <- TukeyHSD(cfu_anova, conf.level = 0.95)
#plot(tHSD_cfu , las=1 , col="brown")


generate_label_df <- function(tHSD_cfu, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tHSD_cfu[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Cage=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Cage) , ]
  return(Tukey.labels)
}
LABELS=generate_label_df(tHSD_cfu , "Condition")

tHSD_label <- LABELS %>%
  mutate(yhsd = 1e+08)

compare_means(CFU ~ Diet, group.by = "Colonization", data = CFUs)
```

    ## # A tibble: 3 × 9
    ##   Colonization .y.   group1 group2        p  p.adj p.format p.signif method  
    ##   <fct>        <chr> <chr>  <chr>     <dbl>  <dbl> <chr>    <chr>    <chr>   
    ## 1 Sn           CFU   SW     P      0.000388 0.0012 0.00039  ***      Wilcoxon
    ## 2 Gi           CFU   SW     P      0.0470   0.047  0.04696  *        Wilcoxon
    ## 3 Sn_Gi        CFU   SW     P      0.00782  0.016  0.00782  **       Wilcoxon

``` r
diet_comp <- list( c("Sn_SW", "Sn_P"),c("Gi_SW", "Gi_P"),c("Sn_Gi_SW", "Sn_Gi_P"))

ggplot(data = CFUs, mapping = aes(x = Condition, y = CFU)) +
  geom_quasirandom(fill="grey80",color = "black", height = 0, width = 0.33, shape=21, size = 1) +
       stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon", size=0.25) +
  theme(axis.title.x = element_text(vjust = -0.6, size=8),
    panel.background = element_rect(fill = "white", colour = "grey50") ,   legend.position = "none",
        axis.text.x = element_text(color = "black", size=8), 
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.y = element_text(color = "black", size=8),
        axis.title.y =element_text(color = "black", size=8),
        legend.text = element_text(color = "black", size=6), 
        legend.title = element_text(color = "black", size=6),
        plot.title = element_text(hjust = 0.5, size=8))+
  scale_y_log10( name = "CFU's / Gut", breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)),
                 limits = c(1e+02, 5e+9)) +
  scale_x_discrete(labels=c('MD','MD (P)', 'Sn', 'Sn (P)', 'Gi', 'Gi (P)','Sn+Gi', 'Sn+Gi (P)')) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept=140, col="grey50", linetype = "dashed") +
  #geom_text(data = tHSD_label, aes(x = Cage,  y = yhsd, label = Letters), size=5) +
  stat_compare_means(comparisons = diet_comp, label =  "p.signif", label.y=c(8,8.5,9)) +
  stat_n_text(na.rm = TRUE, size=2, y.pos = 2.5)
```

    ## [1] FALSE

![](Bee_Gut_CFU_Count_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
  ggsave("Figure_S1.pdf" ,device=cairo_pdf, width=4, height= 3, units = "in", dpi = 300)

  
 (ggplot(subset(CFUs, Colonization %in% "Sn"), mapping = aes(x = Condition, y = CFU)) +
  geom_quasirandom(aes(fill=Condition),color = "grey30", height = 0, width = 0.25, shape=21, size = 2) +
     stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon") +
     stat_summary(fun.y = median, geom = "line", color="grey20",width = 1, aes(group = 1))+
  scale_fill_manual(values=c("skyblue","royalblue")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = "none", 
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8),
        axis.title.x = element_text(color = "black",size = 12, face = "bold"), 
        axis.title.y = element_text(color = "black",size = 12, face = "bold")) +
  scale_y_log10( name = "CFU's / Gut", breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)),
                 limits = c(1e+02, 1e+08)) +
  scale_x_discrete(labels=c('Sn', 'Sn (P)')) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept=140, col="grey50", linetype = "dashed") +
  stat_n_text(na.rm = TRUE, size=2, y.pos = 2)+
   facet_wrap(~Experiment, nrow=1))
```

![](Bee_Gut_CFU_Count_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
  CFU_summary <-CFUs %>%
  group_by(Condition) %>%
  summarise(avg = mean(CFU),
            med = median(CFU),
            stdev = sd(CFU))
```

``` r
CFUs_strain <- data.frame(read.csv("CFU_Count_AM.csv"))

CFUs_strain$Colony_count <- as.numeric(CFUs_strain$Colony_count)
CFUs_strain$Dilution <- as.numeric(CFUs_strain$Dilution)


CFUs_strain <- CFUs_strain%>%
  mutate(cfu_count = Colony_count * 700 / 10 * 10 ^ Dilution) %>%
  drop_na(cfu_count) 

Success <- CFUs_strain %>%
  group_by(Species, Colonization) %>%
  summarise(n_bees = n(),
            n_col = sum(cfu_count > 140),
            colonized = n_col / n_bees) %>%
  mutate(colonized = scales::percent(colonized, accuracy = 1),
         ypos = 20)

Success_native <- CFUs_strain %>%
  filter(Native != "None") %>%
  group_by(Species, Native) %>%
  summarise(n_bees = n(),
            n_col = sum(cfu_count > 140),
            colonized = n_col / n_bees) %>%
  mutate(colonized = scales::percent(colonized, accuracy = 1),
         n_md = n_bees - n_col)  

Success_apis <- Success_native %>%
  filter(Species == "Apis") %>%
  select(Native, n_col,n_md) %>%
  column_to_rownames(var = "Native") %>%
  select(n_col,n_md)

fisher.test(Success_apis)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  Success_apis
    ## p-value = 1.352e-06
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  6.164691      Inf
    ## sample estimates:
    ## odds ratio 
    ##        Inf

``` r
cfu_anova <- aov(log10(cfu_count)~ Colonization, data=CFUs_strain, na.action = na.omit)
tHSD_cfu <- TukeyHSD(cfu_anova, conf.level = 0.95)

generate_label_df <- function(tHSD_cfu, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- tHSD_cfu[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Colonization=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Colonization) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS=generate_label_df(tHSD_cfu, "Colonization")


CFUs_strain$Colonization <- as.character(CFUs_strain$Colonization)
CFUs_strain$Colonization <- factor(CFUs_strain$Colonization, levels=c('MD', 'wkb2', 'ESL0304',  'ESL0323',  'ESL0251','ESL0897'))


tHSD_label <- LABELS %>%
  mutate(yhsd = 1e+08)

ggplot(data = CFUs_strain, mapping = aes(x = Colonization, y = cfu_count)) +
  geom_jitter(color="black", height = 0, width = 0.25, size = 2, stroke = 0) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon") +
  #scale_fill_manual(values=c("tan", "skyblue","steelblue","royalblue","orchid","salmon")) +
  #scale_color_manual(values=c("tan", "skyblue","steelblue","royalblue","orchid","salmon")) +
 theme(text=element_text(size=8,  family="Arial"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(color = "black", size=8, angle = 90), axis.text.y = element_text(color = "black", size=8),
        axis.title.x =element_text(color = "black", size=8), axis.title.y =element_text(color = "black", size=8))+
  
  scale_y_log10( name = "CFU's / Gut", breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),limits = c(1e+01, 1e+08)) +
  annotation_logticks(sides = "l") +
  scale_x_discrete(breaks=c('MD',   'wkb2', 'ESL0304',  'ESL0323',  'ESL0251',  'ESL0893',  'ESL0897'), 
                   labels=c('MD',   'wkb2', 'ESL0304',  'ESL0323',  'ESL0251',  'ESL0893',  'ESL0897')) +
  geom_vline(xintercept=4.5, col="grey50", linetype = "dashed") +
  geom_hline(yintercept=140, col="grey50", linetype = "dashed") +
  geom_text(data = tHSD_label, aes(x = Colonization,  y = yhsd, label = Letters), size=3) +
  stat_n_text(na.rm = TRUE, size=3, y.pos = 1) +
  geom_text(data = Success, aes(x = Colonization,  y = ypos, label = colonized), size=3) +
  geom_text(aes(x = 2.5, y = 45), color="black", label = 'Colonization Success', size=3) #+
```

![](Bee_Gut_CFU_Count_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
 # force_panelsizes(rows = unit(3, "in"), cols = unit(3, "in"))

  #ggsave("AM_cfu_plot_2.pdf", device=cairo_pdf, width=3.5, height= 4, units = "in", dpi = 300)


  
CFU_summary_strains <-CFUs_strain %>%
  group_by(Species, Colonization) %>%
  summarise(avg = mean(cfu_count),
            med = median(cfu_count),
            stdev = sd(cfu_count))  
```