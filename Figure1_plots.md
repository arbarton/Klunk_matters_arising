# Matters arising: Generate Fig 1 plots

## Get data from previous runs


```R
# Read in p-values from the permutation
London_p_maf10_e05_permute <- read.csv('~/London_p_maf10_e05_permute.txt',header=TRUE)
colnames(London_p_maf10_e05_permute) <- "London_p_maf10_e05"
London_p_maf10_e01_permute <- read.csv('~/London_p_maf10_e01_permute.txt',header=TRUE)
colnames(London_p_maf10_e01_permute) <- "London_p_maf10_e01"
London_p_maf30_e01_permute <- read.csv('~/London_p_maf30_e01_permute.txt',header=TRUE)
colnames(London_p_maf30_e01_permute) <- "London_p_maf30_e01"

permute_pvalues <- data.frame(cbind(London_p_maf10_e05_permute,London_p_maf10_e01_permute,London_p_maf30_e01_permute))
```


```R
# Read in p-values from the coverage analysis
London_p_maf10_e05_coverage <- read.csv('~/London_p_maf10_e05_coverage.txt',header=TRUE)
colnames(London_p_maf10_e05_coverage) <- "London_p_maf10_e05"
London_p_maf10_e01_coverage <- read.csv('~/London_p_maf10_e01_coverage.txt',header=TRUE)
colnames(London_p_maf10_e01_coverage) <- "London_p_maf10_e01"
London_p_maf30_e01_coverage <- read.csv('~/London_p_maf30_e01_coverage.txt',header=TRUE)
colnames(London_p_maf30_e01_coverage) <- "London_p_maf30_e01"

coverage_pvalues <- data.frame(cbind(London_p_maf10_e05_coverage,London_p_maf10_e01_coverage,London_p_maf30_e01_coverage))
```

## Make Density plots 


```R
library('ggplot2')
```


```R
cdensity_permute <- ggplot(permute_pvalues, aes(London_p_maf10_e01)) +
stat_ecdf(geom = "step") + scale_x_continuous(trans='log10',limits=c(2.745034e-22,0.05),breaks = c(1e-17, 1e-11, 1e-5,0.05),labels=c(1e-17, 1e-11, 1e-5,0.05),name="P-value") + theme_bw(base_size = 24) + theme(axis.ticks = element_line(size = 1),
 panel.grid = element_blank()) + ylab("Cumulative density")  + geom_hline(yintercept=0.05, linetype="dashed", color="darkblue",size=1) +  geom_vline(xintercept=0.05, linetype="dashed", color="darkgreen",size=1) + geom_vline(xintercept=7.89e-12,
color="red",size=1)
cdensity_permute
```

    Warning message:
    “Removed 44 rows containing non-finite values (stat_ecdf).”



    
![png](output_6_1.png)
    



```R
cdensity_coverage <- ggplot(coverage_pvalues, aes(London_p_maf10_e01)) +
stat_ecdf(geom = "step") + scale_x_continuous(trans='log10',limits=c(2.745034e-22,0.05),breaks = c(1e-17, 1e-11, 1e-5,0.05),labels=c(1e-17, 1e-11, 1e-5,0.05),name="P-value") + theme_bw(base_size = 24) + theme(axis.ticks = element_line(size = 1),
 panel.grid = element_blank()) + ylab("Cumulative density")  + geom_hline(yintercept=0.05, linetype="dashed", color="darkblue",size=1) +  geom_vline(xintercept=0.05, linetype="dashed", color="darkgreen",size=1) + geom_vline(xintercept=7.89e-12,
color="red",size=1)
cdensity_coverage
```

    Warning message:
    “Removed 23 rows containing non-finite values (stat_ecdf).”



    
![png](output_7_1.png)
    


## Make q-q plots


```R
# Code for generating q-q plot adopted from Kamil Slowikowski: https://slowkow.com/notes/ggplot2-qqplot/

gg_qqplot <- function(ps, ci = 0.95) {
    n  <- length(ps)
    df <- data.frame(
        observed = -log10(sort(ps)),
        expected = -log10(ppoints(n)),
        clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
        geom_ribbon(
        mapping = aes(x = expected, ymin = clower, ymax = cupper),
        alpha = 0.1
        ) +
        geom_point(aes(expected, observed), shape = 1, size = 3) +
        geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
        xlab(log10Pe) +
        ylab(log10Po)
}

inflation <- function(ps) {
chisq <- qchisq(1 - ps, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda
}
```


```R
qqplot_permute <- gg_qqplot(permute_pvalues$London_p_maf10_e01) +
    theme_bw(base_size = 24) +
    annotate(
        geom = "text",
        x = -Inf,
        y = Inf,
        hjust = -0.15,
        vjust = 1 + 0.15 * 3,
        label = sprintf("λ = %.2f", inflation(permute_pvalues$London_p_maf10_e01)),
        size = 8
        ) +
    theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank()
)
qqplot_permute
```


    
![png](output_10_0.png)
    



```R
qqplot_coverage <- gg_qqplot(coverage_pvalues$London_p_maf10_e01) +
    theme_bw(base_size = 24) +
    annotate(
        geom = "text",
        x = -Inf,
        y = Inf,
        hjust = -0.15,
        vjust = 1 + 0.15 * 3,
        label = sprintf("λ = %.2f", inflation(coverage_pvalues$London_p_maf10_e01)),
        size = 8
        ) +
    theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank()
)
qqplot_coverage
```


    
![png](output_11_0.png)
    


## Make final figure


```R
library(patchwork)
```

    Warning message:
    “package ‘patchwork’ was built under R version 4.2.1”



```R
Figure_1 <- (cdensity_permute + qqplot_permute) / (cdensity_coverage + qqplot_coverage) + plot_annotation(tag_levels = 'A')
Figure_1
```

    Warning message:
    “Removed 44 rows containing non-finite values (stat_ecdf).”
    Warning message:
    “Removed 23 rows containing non-finite values (stat_ecdf).”



    
![png](output_14_1.png)
    



```R

```
