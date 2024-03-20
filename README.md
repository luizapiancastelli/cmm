Clustered Mallows Model
================
Luiza Piancastelli

Configuration requirement: please ensure that you have a C++ compiler
properly set up.

``` r
# Package names
packages <- c("dplyr", "reshape2", "Rcpp", "ggplot2", "gridExtra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


sourceCpp("cmm_cpp.cpp")
source("cmm_functions.R")
```

## Simulating CMM data

Our first code chunk simulates data from the CMM with the Hamming
$d_{oc}$, allocations vector $\boldsymbol = (1,1,2,2,3,3,4,4)$ and
$\theta = 0.8$. This means that the CT is (2,2,2,2), which is a
partition of the total $n=8$ items into $L=4$ ordered groups, each with
two items. The function `partition_from_z` can be used to produce a list
of size $L$ giving which items are in each partition. In this example
our item set is arbitrarily named $\{1, \ldots, n\}$ but a string with
item names can be provided as well.

The function `rcmm` produces $q$ random realizations of a CMM with given
$\boldsymbol{z}$ and $\theta$. Here we ask for $q=100$ simulated ranks.

``` r
z = sort(rep(1:4, 2))
theta = 0.8

partition_from_z(z)
```

    ## [[1]]
    ## [1] 1 2
    ## 
    ## [[2]]
    ## [1] 3 4
    ## 
    ## [[3]]
    ## [1] 5 6
    ## 
    ## [[4]]
    ## [1] 7 8

``` r
Pi = rcmm(100, theta, z, 'h')

head(Pi)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## [1,]    2    5    4    7    3    1    6    8
    ## [2,]    4    7    1    2    5    8    6    3
    ## [3,]    7    1    4    3    5    6    2    8
    ## [4,]    2    8    4    7    5    6    1    3
    ## [5,]    2    4    3    1    6    5    7    8
    ## [6,]    5    6    4    3    2    8    7    1

``` r
dim(Pi)
```

    ## [1] 100   8

``` r
apply(Pi, 2, function(x){prop.table(table(x))})
```

    ##   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## 1 0.19 0.33 0.05 0.06 0.09 0.06 0.11 0.11
    ## 2 0.27 0.15 0.10 0.12 0.11 0.08 0.08 0.09
    ## 3 0.06 0.05 0.22 0.28 0.08 0.14 0.09 0.08
    ## 4 0.06 0.13 0.28 0.15 0.07 0.07 0.16 0.08
    ## 5 0.11 0.06 0.07 0.10 0.31 0.23 0.04 0.08
    ## 6 0.10 0.09 0.12 0.03 0.20 0.24 0.07 0.15
    ## 7 0.12 0.08 0.05 0.12 0.08 0.10 0.25 0.20
    ## 8 0.09 0.11 0.11 0.14 0.06 0.08 0.20 0.21

The final piece of code shows the proportion of times that each position
(row) is occupied by each item (column). As expected since items
$\{1,2\}$ are the first partition, there is a higher prevalence of them
in the first and second positions of the ranks. A similar comment
applies to $\{3, 4\}$ in 3$^{rd}$ and 4$^{th}$ places and so on. Another
important aspect of the CMM is that same group items are exchangeable,
so the probabilities for items with the same allocation are
approximately equal.

## Maximum likelihood estimation

We now apply the maximum likelihood estimation method in sections 6.1
and 6.2 of the paper to the simulated data. The routine `MLE_fit_cmm`
calls the simulated annealing algorithm with some supplied schedule, and
the $\theta$ iterative moment-matching routine that takes tolerance
`eps`.

``` r
anneal = seq(0.01, 2, 0.01) #Annealing sequence for z optimization
eps = 10^(-2)               #Tolerance of algorithm

CT = rep(2, 4) #Clustering Table (CT)
CT
```

    ## [1] 2 2 2 2

``` r
mle = MLE_fit_cmm(CT, Pi, anneal, eps, dist = 'h')

mle$z_mle
```

    ## [1] 1 1 2 2 3 3 4 4

``` r
mle$theta_mle
```

    ## [1] 0.8233959

A graphical summary of the fitted partition can be obtained with
`plot_rank_probs` as illustrated below. In this plot, each bar
corresponds to a ranking position, and colors indicate the prevalence of
clusters therein.

``` r
plot_rank_probs(Pi, mle$z_mle)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Let us also fit the simulated data with the mispecified distance
(Kendall) and investigate model choice. The criterion that we use is
outlined in the manuscript section 4. It is given by the optimized
log-likelihood function value, penalised for complexity of the
partition.

``` r
mle2 = MLE_fit_cmm(CT, Pi, anneal, eps, dist = 'k')

I_criterion(mle, 10000)
```

    ## [1] -998.7394

``` r
I_criterion(mle2, 10000)
```

    ## [1] -1042.304

The maximum value of `I_criterion` indicates the best model, and this
can also be applied to the $CT$ choice.

## Bayesian inference with AEA

Finally, the Bayesian inferential approach is implemented in `cmm_MCMC`.
This function relies on the Approximate Exchange Algorithm routine,
described in Algorithm 2 of the manuscript. It takes a pre-specified
$CT$, the observed data $\underline{\pi}$, and choices of prior
parameters $p(\theta) ~ \mbox{Gamma(shape, rate)}$, number of iterations
and distance.

``` r
mcmc = cmm_MCMC(CT, Pi, 1000, "h", prior_theta= list(shape = 2, rate = 2))
```

    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500
    ## [1] 600
    ## [1] 700
    ## [1] 800
    ## [1] 900
    ## [1] 1000

``` r
mcmc$accepted_theta
```

    ## [1] 0.389

``` r
theta_chain = data.frame("iter" = 1:1000, "value" = mcmc$theta)
density = ggplot(theta_chain, aes(x = value))+geom_density(adjust =5, fill = 'cyan3', alpha = 0.5)+xlim(0,1)+theme_bw()
trace  = ggplot(theta_chain, aes(y = value, x = iter))+geom_line()+theme_bw()

grid.arrange(density, trace, ncol =2)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
