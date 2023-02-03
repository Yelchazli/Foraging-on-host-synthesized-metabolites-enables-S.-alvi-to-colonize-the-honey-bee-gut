#### These are the required package to reproduce this analysis.

``` r
### env setup

required_packages = c("dplyr", "data.table", "tidyverse" )
need_install = required_packages[!(required_packages) %in% installed.packages()]

if (length(need_install) > 0){
  install.packages((need_install))
}

# Load packages
lapply(required_packages, require, character.only=T)

dir_script = dirname(rstudioapi::getSourceEditorContext()$path)

setwd(dir_script)
#font_import()
#loadfonts(device="win")
```

``` r
rawdata = read.csv("curated_values.csv", sep=",", header = T)


### Split Df into a list with the 4 different targets
rawdatasplit = split(rawdata, rawdata$Target) 


### For each of the dataframes in the list, group by sample_ID and take the one with lowest SD
### This gets rid of rerun samples
rawdatasplit_filtered = 
  lapply(rawdatasplit, function(x){
   x %>%  group_by(Sample_ID ) %>% 
    
    filter(SD %in% c(min(SD)))
})


### Rebind the elements of the list into one dataframe
rawdatafiltered = do.call(rbind.data.frame, rawdatasplit_filtered)

#### Split the df into a list based on the calibration and target
rawdatasplit_cal = split(rawdatafiltered, list(rawdatafiltered$Calibration, rawdatafiltered$Target)) 

### Remove unecessary tables
rm(rawdatasplit,rawdata)

#### Calculating copy number, 200 is the volume of DNA elution, 1.627 is the dilution factor and 1.2 

# efficiency and intercept values from standard curves of 2019 calibration
Eff_cal_2019 <- list(Actin_eff= 1.956, Actin_int=37.1295885699136, UV_eff=2.009,
                     UV_int = 38.32411214,Ga_eff= 1.750, Ga_int=44.97112553,
                     Sa_eff=1.935,Sa_int = 38.7608650411878 )

# efficiency and intercept values from standard curves of 2017 calibration
Eff_cal_2017 <- list(Actin_eff= 2.022, Actin_int=35.09421199, UV_eff=1.945,
                     UV_int = 36.11530377,Ga_eff= 1.882, Ga_int=38.54214641,
                     Sa_eff=1.913,Sa_int = 37.4529874 )


### Enter the efficiency and intercept values for Copynum calculation

Target = rep(c("Actin", "16S", "Ga", "Sa"), each=2)
Calibration =c("C2017", "C2019")
Eff = c(2.022,1.956, 1.945, 2.009, 1.882, 1.750, 1.913, 1.935)
Int = c(35.09421199, 37.1295885699136, 36.11530377, 38.32411214, 38.54214641, 44.97112553, 37.4529874, 38.7608650411878)

## Create Df with vectors
std_curves= data.frame(Target, Calibration, Eff, Int)



## Calculate Copy number
### for each calibration and target copynum is the formula below
setDT(std_curves)
normlaized_split= 
  lapply(rawdatasplit_cal, function(x){
   setDT(x) 
   x[std_curves, on = c("Calibration", "Target"),  Copynum := (i.Eff^(i.Int-Ct))*200*1.627906977][]
})

### Bind back the list into a single df
normalized.df = do.call(rbind.data.frame, normlaized_split)


###########
##############
############### Stopped here 26.10.21

normalized.df$Copynum[is.na(normalized.df$Copynum)] <-  1
normalized.df$Copynum[normalized.df$Copynum <=36000] <-  36000

#### Group by Sample, and calculate normalized copy number by dividing by the actin value and multiplying by the 
# median of actin of the whole data set
normalized.df[,by=.(Sample_ID) ,norm_copynum := (Copynum/Copynum[Target=="Actin"])][,
  norm_copynum := norm_copynum*median(Copynum[Target=="Actin"])][ 
    ,Norm_cell := norm_copynum/4] %>% print()

actin_limit <- (Eff_cal_2017[[1]]^(Eff_cal_2017[[2]] - 28.72))*200*1.627906977

UV_limit <- (Eff_cal_2017[[3]]^(Eff_cal_2017[[4]] - 30.09))*200*1.627906977

Ga_limit <- (Eff_cal_2017[[5]]^(Eff_cal_2017[[6]] - 31.26))*200*1.627906977

Sa_limit <- (Eff_cal_2017[[7]]^(Eff_cal_2017[[8]] - 30.41))*200*1.627906977


# take all values below Calculated limit of detection

normalized.df$norm_copynum[is.na(normalized.df$norm_copynum)] <- 32000
normalized.df$norm_copynum[normalized.df$norm_copynum < 32000] <- 32000

normalized.df$Norm_cell[is.na(normalized.df$Norm_cell)] <- 8000
normalized.df$Norm_cell[normalized.df$Norm_cell < 8000] <- 8000


write.csv(normalized.df, file.path(dir_script, "normalized_data.csv"), row.names = FALSE) # export data in .csv
```
