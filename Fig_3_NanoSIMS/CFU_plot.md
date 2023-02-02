``` r
library(knitr)
hook_output = knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
    # this hook is used only when the linewidth option is not NULL
    if (!is.null(n <- options$linewidth)) {
        x = knitr:::split_lines(x)
        # any lines wider than n should be wrapped
        if (any(nchar(x) > n))
            x = strwrap(x, width = n)
        x = paste(x, collapse = "\n")
    }
    hook_output(x, options)
})
```

#### These are the required package to reproduce this analysis.

``` r
### env setup

required_packages = c("dplyr", "data.table", "ggplot2", "readxl", "ggbeeswarm", "ggpubr",
    "EnvStats", "tidyverse", "patchwork", "ggforce", "wesanderson", "RcolorBrewer",
    "extrafont", "ggthemes", "egg", "scales", "ggh4x", "grid")
need_install = required_packages[!(required_packages) %in% installed.packages()]

if (length(need_install) > 0) {
    install.packages((need_install))
}
```

    ## Error in contrib.url(repos, "source"): essai d'utilisation de CRAN sans fixer un miroir

``` r
# Load packages
lapply(required_packages, require, character.only = T)

dir_script = dirname(rstudioapi::getSourceEditorContext()$path)

setwd(dir_script)
# font_import() loadfonts(device='win')
```

#### Load the data

``` r
cfu_count <- fread("CFU_count.csv") %>%
    drop_na() %>%
    group_by(Time) %>%
    mutate(Time_d = Time/24, cfumax = max(CFU), cfumin = min(CFU), Colonized = case_when(CFU >
        50 ~ 1, TRUE ~ 0))
cfu_summary <- cfu_count %>%
    group_by(Time) %>%
    filter(Gut == "Ileum") %>%
    summarise(avg = mean(CFU), stdev = sd(CFU), Colonization = sum(Colonized)/n() *
        100)



cfu_summary <- cfu_count %>%
    group_by(Time) %>%
    filter(Gut == "Ileum") %>%
    summarise(avg = mean(CFU), stdev = sd(CFU), Colonization = sum(Colonized)/n() *
        100)

inoculum = cfu_count[cfu_count$Time != 0, ]

setDT(cfu_count)
```

##### Figure 3B

``` r
ggplot(cfu_count, aes(x = Time, y = CFU)) + stat_summary(data = cfu_count[Condition !=
    "Inoculum"], fun = mean, geom = "line", color = "gray", alpha = 0.5, stroke = 0) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange",
        color = "gray", alpha = 0.5, stroke = 0) + geom_beeswarm(data = cfu_count[Time ==
    0], size = 1.2, colour = "blue", alpha = 0.5, width = 1, stroke = 0) + geom_quasirandom(data = cfu_count,
    fill = "goldenrod2", height = 0, width = 0.5, size = 1.2, colour = "black", alpha = 0.5,
    stroke = 0) + geom_hline(yintercept = 50, col = "black", linetype = "dashed") +
    theme(text = element_text(size = 8, family = "Arial"), panel.background = element_rect(fill = "white",
        colour = "grey50"), panel.border = element_rect(colour = "grey", fill = NA,
        size = 1), axis.text.x = element_text(color = "black", size = 8), axis.text.y = element_text(color = "black",
        size = 10), axis.title.x = element_text(color = "black", size = 8), axis.title.y = element_text(color = "black",
        size = 10)) + scale_y_log10(name = "CFUs / Ileum", breaks = trans_breaks("log10",
    function(x) 10^x, n = 10), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 8)) + xlab("Time (h)") +
    annotate("text", x = 15, y = 14, size = 2, label = "Colonization Success (%)") +
    geom_text(data = cfu_summary, aes(label = paste0(round(Colonization, 0), "%"),
        y = 6), size = 2)  #+
```

![](CFU_plot_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# force_panelsizes(rows = unit(1.6, 'in'), cols = unit(2.4, 'in'))

# ggsave('Fig.3B.pdf' ,device=cairo_pdf, width=3, height= 2, units = 'in', dpi
# = 300)
```
