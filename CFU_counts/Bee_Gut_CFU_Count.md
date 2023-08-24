Bee Gut CFU Counts
================
Andrew Quinn, University of Lausanne
1/11/2022

# Generate plots for folony forming units (CFU’s) plated from bee guts, as well as weight calculations.

## Figure S1A, Figure S4, Figure 5,

## Environment setup

### Install and load libraries (run code if missing libraries)

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

# Plot CFU counts across diet and colonizations groups

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
library(cowplot)
library(ggfortify)
library(lme4)
library(nlme)

CFUs <- data.frame(read.csv("All_CFU_Count.csv")) 
CFUs$Colonization <- factor(CFUs$Colonization, levels=c("MF", "Sn", "Gi", "Sn_Gi"))
CFUs$Diet <- factor(CFUs$Diet, levels=c("SW", "P"))
CFUs$Condition <- factor(CFUs$Condition, levels=c("MF_SW", "MF_P", "Sn_SW", "Sn_P","Gi_SW", "Gi_P","Sn_Gi_SW", "Sn_Gi_P"))
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

    ## # A tibble: 3 x 9
    ##   Colonization .y.   group1 group2        p  p.adj p.format p.signif method  
    ##   <fct>        <chr> <chr>  <chr>     <dbl>  <dbl> <chr>    <chr>    <chr>   
    ## 1 Sn_Gi        CFU   SW     P      0.0310   0.062  0.03098  *        Wilcoxon
    ## 2 Sn           CFU   SW     P      0.000388 0.0012 0.00039  ***      Wilcoxon
    ## 3 Gi           CFU   SW     P      0.235    0.23   0.23478  ns       Wilcoxon

``` r
diet_comp <- list( c("Sn_SW", "Sn_P"),c("Gi_SW", "Gi_P"),c("Sn_Gi_SW", "Sn_Gi_P"))

ggplot(data = CFUs, mapping = aes(x = Condition, y = CFU)) +
  geom_quasirandom(fill="grey70",color = "black", height = 0, width = 0.25, shape=21, size = 1.5, stroke=0.25) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon") +
  theme(axis.title.x = element_text(vjust = -0.6, size=7),
    panel.background = element_rect(fill = "white", colour = "grey50") ,   legend.position = "none",
        axis.text.x = element_text(color = "black", size=7), 
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.y = element_text(color = "black", size=7),
        axis.title.y =element_text(color = "black", size=7),
        legend.text = element_text(color = "black", size=7), 
        legend.title = element_text(color = "black", size=7),
        plot.title = element_text(hjust = 0.5, size=7),
         aspect.ratio=1)+
  scale_y_log10( name = "CFU's / Gut", breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)),
                 limits = c(1e+01, 5e+9)) +
  scale_x_discrete(labels=c('MF','MF (P)', 'Sn', 'Sn (P)', 'Gi', 'Gi (P)','Sn+Gi', 'Sn+Gi (P)')) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept=70, col="grey50", linetype = "dashed") +
  stat_compare_means(comparisons = diet_comp, label =  "p.signif", label.y=c(8,8.5,9)) +
  stat_n_text(na.rm = TRUE, size=2, y.pos = 1)+
  geom_text(aes(x = 3, y = 120), color="black", label = 'LOD', size=2)
```

![](Bee_Gut_CFU_Count_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
 # ggsave("Figure_S1A.pdf" ,device=cairo_pdf, width=90, height= 120, units = "mm", dpi = 600)

  
 (ggplot(subset(CFUs, Colonization %in% "Sn"), mapping = aes(x = Condition, y = CFU)) +
  geom_quasirandom(fill="grey70",color = "black", height = 0, width = 0.25, shape=21, size = 1.5, stroke=0.25) +
     stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon") +
     stat_summary(fun.y = median, geom = "line", color="grey20",width = 1, aes(group = 1))+
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
                 limits = c(1e+01, 5e+9)) +
  scale_x_discrete(labels=c('Sn', 'Sn (P)')) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept=70, col="grey50", linetype = "dashed") +
  stat_n_text(na.rm = TRUE, size=2, y.pos = 1)+
  geom_text(aes(x = 1, y = 120), color="black", label = 'LOD', size=2)+
  facet_wrap(~Experiment, nrow=1))
```

![](Bee_Gut_CFU_Count_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
  CFU_summary <-CFUs %>%
  group_by(Condition) %>%
  summarise(avg = mean(CFU),
            med = median(CFU),
            stdev = sd(CFU))
```

## Bee weights

``` r
Weights <- CFUs %>%
  select(-CFU) %>%
  mutate(Bee_weight = 1000*Bee_weight,  # convert g to mg
         Gut_weight = 1000*Gut_weight,
         Body_weight=Bee_weight-Gut_weight)

Weights$Condition <- as.character(Weights$Condition)
Weights$Condition <- factor(Weights$Condition, levels=c("MF_SW", "MF_P", "Sn_SW", "Sn_P","Gi_SW", "Gi_P","Sn_Gi_SW", "Sn_Gi_P"))


p1<- ggplot(data = Weights, mapping = aes(x = Condition, y = Bee_weight)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(fill="grey70",color = "black", height = 0, width = 0.15, shape=21, size = 1, stroke=0) +
  theme(text=element_text(size=7,  family="Arial"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(color = "black", size=5), 
        axis.text.y = element_text(color = "black", size=6),
        axis.title.x =element_text(color = "black", size=7), 
       axis.title.y =element_text(color = "black", size=7),
         aspect.ratio=1)+
  scale_y_continuous( name = "Bee weight (mg)",limits=c(0,200),breaks=seq(0,200,25))+
  stat_n_text(na.rm = TRUE, size=1.5, y.pos = 1) 

p2<- ggplot(data = Weights, mapping = aes(x = Condition, y = Gut_weight)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(fill="grey70",color = "black", height = 0, width = 0.15, shape=21, size = 1, stroke=0) +
  theme(text=element_text(size=7,  family="Arial"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(color = "black", size=5), 
        axis.text.y = element_text(color = "black", size=6),
        axis.title.x =element_text(color = "black", size=7), 
       axis.title.y =element_text(color = "black", size=7),
         aspect.ratio=1)+
  scale_y_continuous( name = "Gut weight (mg)",limits=c(0,125),breaks=seq(0,125,25))+
  stat_n_text(na.rm = TRUE, size=1.5, y.pos = 1) 

p3<- ggplot(data = Weights, mapping = aes(x = Condition, y = Body_weight)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(fill="grey70",color = "black", height = 0, width = 0.15, shape=21, size = 1, stroke=0) +
  theme(text=element_text(size=7,  family="Arial"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(color = "black", size=5), 
        axis.text.y = element_text(color = "black", size=6),
        axis.title.x =element_text(color = "black", size=7), 
       axis.title.y =element_text(color = "black", size=7),
         aspect.ratio=1)+
  scale_y_continuous( name = "Body weight (mg)",limits=c(0,150),breaks=seq(0,150,25))+
  stat_n_text(na.rm = TRUE, size=1.5, y.pos = 1) 

plot_grid(p1, p2, p3, nrow=1)
```

![](Bee_Gut_CFU_Count_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# ggsave("FigSx_Weights_v2.pdf" ,device=cairo_pdf, width=210, height= 90, units = "mm", dpi = 600)

Weight_summary_strains <-Weights %>%
  group_by(Condition) %>%
  summarise(across(c(Bee_weight,Gut_weight,Body_weight), 
                   .f = list(mean = mean, sd = sd), na.rm = TRUE))

Diet_weight_dif <- Weights %>%
  group_by(Diet)%>%
  summarise(across(c(Bee_weight,Gut_weight,Body_weight), 
                   .f = list(mean = mean, sd = sd), na.rm = TRUE))
```

## Test for diet and colonization effects on weights

``` r
weight_lmm <- Weights %>%
  select(-Condition) %>%
  mutate(Sn=case_when(str_detect(Colonization,"Sn")~1,
                      TRUE ~ 0),
         Gi=case_when(str_detect(Colonization,"Gi")~1,
                      TRUE ~ 0))%>%
  ungroup() %>%
  select(-Colonization)%>%
  mutate(Diet =as.factor(Diet),
         Sn=as.factor(Sn),
         Gi=as.factor(Gi),
         Experiment=as.factor(Experiment))%>%
  drop_na(Body_weight)

bee_lme_v0 = lme(Bee_weight ~ Diet + Sn*Gi, random = ~1|Experiment, data = weight_lmm, method = "ML")
bee_lme_v1 <- update(bee_lme_v0, weights = varIdent(form = ~ 1 | Experiment))
anova(bee_lme_v0,bee_lme_v1)
```

    ##            Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    ## bee_lme_v0     1  7 2301.279 2326.494 -1143.639                        
    ## bee_lme_v1     2 11 2300.113 2339.736 -1139.056 1 vs 2 9.165795  0.0571

``` r
summary(bee_lme_v1)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: weight_lmm 
    ##        AIC      BIC    logLik
    ##   2300.113 2339.736 -1139.056
    ## 
    ## Random effects:
    ##  Formula: ~1 | Experiment
    ##         (Intercept) Residual
    ## StdDev: 0.001265557 17.91364
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Experiment 
    ##  Parameter estimates:
    ##         1         2         3         4         6 
    ## 1.0000000 0.7279693 0.8417823 1.0707895 0.9826679 
    ## Fixed effects: Bee_weight ~ Diet + Sn * Gi 
    ##                 Value Std.Error  DF  t-value p-value
    ## (Intercept) 105.48077  2.040206 262 51.70105  0.0000
    ## DietP        25.88653  1.947340 262 13.29328  0.0000
    ## Sn1          -2.31804  2.585523 262 -0.89655  0.3708
    ## Gi1          -0.01232  3.079195 262 -0.00400  0.9968
    ## Sn1:Gi1       3.48425  4.033075 262  0.86392  0.3884
    ##  Correlation: 
    ##         (Intr) DietP  Sn1    Gi1   
    ## DietP   -0.458                     
    ## Sn1     -0.613 -0.023              
    ## Gi1     -0.517 -0.014  0.414       
    ## Sn1:Gi1  0.391  0.020 -0.641 -0.764
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.6527930 -0.6494609 -0.1017718  0.5569148  3.6489141 
    ## 
    ## Number of Observations: 271
    ## Number of Groups: 5

``` r
anova(bee_lme_v1)
```

    ##             numDF denDF   F-value p-value
    ## (Intercept)     1   262 14686.768  <.0001
    ## Diet            1   262   176.178  <.0001
    ## Sn              1   262     0.086  0.7690
    ## Gi              1   262     1.031  0.3108
    ## Sn:Gi           1   262     0.746  0.3884

``` r
gut_lme_v0 = lme(Gut_weight ~ Diet + Sn*Gi, random = ~1|Experiment, data = weight_lmm, method = "ML")
gut_lme_v1 <- update(gut_lme_v0, weights = varIdent(form = ~ 1 | Experiment))
anova(gut_lme_v0,gut_lme_v1)
```

    ##            Model df      AIC      BIC   logLik   Test  L.Ratio p-value
    ## gut_lme_v0     1  7 2181.259 2206.474 -1083.63                        
    ## gut_lme_v1     2 11 2171.899 2211.523 -1074.95 1 vs 2 17.35974  0.0016

``` r
summary(gut_lme_v1)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: weight_lmm 
    ##      AIC      BIC   logLik
    ##   2171.9 2211.523 -1074.95
    ## 
    ## Random effects:
    ##  Formula: ~1 | Experiment
    ##         (Intercept) Residual
    ## StdDev:    3.919051 16.16226
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Experiment 
    ##  Parameter estimates:
    ##         1         2         3         4         6 
    ## 1.0000000 0.5587936 0.8240468 0.8581789 0.7498660 
    ## Fixed effects: Gut_weight ~ Diet + Sn * Gi 
    ##                 Value Std.Error  DF   t-value p-value
    ## (Intercept) 26.874654  2.427023 262 11.073095  0.0000
    ## DietP       21.881590  1.478776 262 14.797097  0.0000
    ## Sn1          2.639956  2.040344 262  1.293878  0.1968
    ## Gi1          3.648958  2.775369 262  1.314765  0.1897
    ## Sn1:Gi1     -5.179663  3.422827 262 -1.513271  0.1314
    ##  Correlation: 
    ##         (Intr) DietP  Sn1    Gi1   
    ## DietP   -0.292                     
    ## Sn1     -0.404 -0.015              
    ## Gi1     -0.403 -0.016  0.382       
    ## Sn1:Gi1  0.323  0.021 -0.597 -0.810
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.4030656 -0.6176212 -0.1188048  0.5282569  4.0124448 
    ## 
    ## Number of Observations: 271
    ## Number of Groups: 5

``` r
anova(gut_lme_v1)
```

    ##             numDF denDF  F-value p-value
    ## (Intercept)     1   262 408.5034  <.0001
    ## Diet            1   262 220.0454  <.0001
    ## Sn              1   262   0.2827  0.5954
    ## Gi              1   262   0.0230  0.8796
    ## Sn:Gi           1   262   2.2900  0.1314

``` r
body_lme_v0 = lme(Body_weight ~ Diet + Sn*Gi, random = ~1|Experiment, data = weight_lmm, method = "ML")
body_lme_v1 <- update(body_lme_v0, weights = varIdent(form = ~ 1 | Experiment))
anova(body_lme_v0,body_lme_v1)
```

    ##             Model df      AIC      BIC    logLik   Test L.Ratio p-value
    ## body_lme_v0     1  7 2268.492 2293.707 -1127.246                       
    ## body_lme_v1     2 11 2259.274 2298.898 -1118.637 1 vs 2 17.2177  0.0018

``` r
summary(body_lme_v1)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: weight_lmm 
    ##        AIC      BIC    logLik
    ##   2259.274 2298.898 -1118.637
    ## 
    ## Random effects:
    ##  Formula: ~1 | Experiment
    ##         (Intercept) Residual
    ## StdDev:    2.972602 15.88971
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Experiment 
    ##  Parameter estimates:
    ##         1         2         3         4         6 
    ## 1.0000000 0.6634722 1.1137099 0.8785420 1.0739697 
    ## Fixed effects: Body_weight ~ Diet + Sn * Gi 
    ##                Value Std.Error  DF  t-value p-value
    ## (Intercept) 78.73713  2.300607 262 34.22451  0.0000
    ## DietP        4.09950  1.748108 262  2.34511  0.0198
    ## Sn1         -5.21036  2.365217 262 -2.20291  0.0285
    ## Gi1         -3.34178  3.105791 262 -1.07598  0.2829
    ## Sn1:Gi1      7.65662  3.920184 262  1.95313  0.0519
    ##  Correlation: 
    ##         (Intr) DietP  Sn1    Gi1   
    ## DietP   -0.363                     
    ## Sn1     -0.462 -0.020              
    ## Gi1     -0.442 -0.016  0.326       
    ## Sn1:Gi1  0.345  0.024 -0.567 -0.792
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.1456164 -0.5932805 -0.2044215  0.4279269  3.0485726 
    ## 
    ## Number of Observations: 271
    ## Number of Groups: 5

``` r
anova(body_lme_v1)
```

    ##             numDF denDF   F-value p-value
    ## (Intercept)     1   262 2317.0686  <.0001
    ## Diet            1   262    5.2373  0.0229
    ## Sn              1   262    1.3857  0.2402
    ## Gi              1   262    0.5939  0.4416
    ## Sn:Gi           1   262    3.8147  0.0519

``` r
intervals(body_lme_v1, level = 0.95,which = "fixed")
```

    ## Approximate 95% confidence intervals
    ## 
    ##  Fixed effects:
    ##                    lower      est.      upper
    ## (Intercept) 74.249082161 78.737130 83.2251771
    ## DietP        0.689275213  4.099503  7.5097303
    ## Sn1         -9.824451623 -5.210362 -0.5962722
    ## Gi1         -9.400589458 -3.341780  2.7170295
    ## Sn1:Gi1      0.009087724  7.656624 15.3041604
    ## attr(,"label")
    ## [1] "Fixed effects:"

## CFU counts for Figure 5

``` r
CFUs_strain <- data.frame(read.csv("CFU_Count_AM.csv"))

CFUs_strain$Colony_count <- as.numeric(CFUs_strain$Colony_count)
CFUs_strain$Dilution <- as.numeric(CFUs_strain$Dilution)

CFUs_strain <- CFUs_strain%>%
  mutate(cfu_count = Colony_count * 700 / 10 * 10 ^ Dilution,
         cfu_count = case_when(cfu_count == 70 ~50,
                               TRUE ~ cfu_count)) %>%
  drop_na(cfu_count) 

Success <- CFUs_strain %>%
  group_by(Species, Colonization) %>%
  summarise(n_bees = n(),
            n_col = sum(cfu_count > 70),
            colonized = n_col / n_bees) %>%
  mutate(colonized = scales::percent(colonized, accuracy = 1),
         ypos = 20)

Success_native <- CFUs_strain %>%
  filter(Native != "None") %>%
  group_by(Species, Native) %>%
  summarise(n_bees = n(),
            n_col = sum(cfu_count > 70),
            colonized = n_col / n_bees) %>%
  mutate(colonized = scales::percent(colonized, accuracy = 1),
         n_mf = n_bees - n_col)  

Success_apis <- Success_native %>%
  filter(Species == "Apis") %>%
  select(Native, n_col,n_mf) %>%
  column_to_rownames(var = "Native") %>%
  select(n_col,n_mf)

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
CFUs_strain$Colonization <- factor(CFUs_strain$Colonization, levels=c('MF', 'wkb2', 'ESL0304',  'ESL0323',  'ESL0251','ESL0897'))

tHSD_label <- LABELS %>%
  mutate(yhsd = 1e+08)

ggplot(data = CFUs_strain, mapping = aes(x = Colonization, y = cfu_count)) +
  geom_hline(yintercept=70, col="grey50", linetype = "dashed", size = 0.75) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                  geom = "crossbar",width = 0.5, color="salmon") +
  geom_quasirandom(fill="grey70",color = "black", height = 0, width = 0.15, shape=21, size = 2, stroke=0.25) +
 theme(text=element_text(size=7,  family="Arial"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        axis.text.x = element_text(color = "black", size=7), 
        axis.text.y = element_text(color = "black", size=7),
        axis.title.x =element_text(color = "black", size=7), 
       axis.title.y =element_text(color = "black", size=7),
         aspect.ratio=1)+
  scale_y_log10( name = "CFU's / Gut", breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),limits = c(1e+01, 1e+08)) +
  annotation_logticks(sides = "l") +
  scale_x_discrete(breaks=c('MF',   'wkb2', 'ESL0304',  'ESL0323',  'ESL0251',  'ESL0893',  'ESL0897'), 
                   labels=c('MF',   'wkB2', 'ESL0304',  'ESL0323',  'ESL0251',  'ESL0893',  'ESL0897')) +
  geom_vline(xintercept=4.5, col="grey50", linetype = "dotted") +
  geom_text(data = tHSD_label, aes(x = Colonization,  y = yhsd, label = Letters), size=5) +
  stat_n_text(na.rm = TRUE, size=2, y.pos = 1) +
  geom_text(data = Success, aes(x = Colonization,  y = ypos, label = colonized), size=2) +
  geom_text(aes(x = 3, y = 45), color="black", label = 'Colonization Success', size=2) +
  geom_text(aes(x = 3, y = 120), color="black", label = 'LOD', size=2)
```

![](Bee_Gut_CFU_Count_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
  #ggsave("Figure_5a.pdf", device=cairo_pdf, width=90, height= 120, units = "mm", dpi = 600)

CFU_summary_strains <-CFUs_strain %>%
  group_by(Species, Colonization) %>%
  summarise(avg = mean(cfu_count),
            med = median(cfu_count),
            stdev = sd(cfu_count))  
```
