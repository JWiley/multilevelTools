# Function to calculate the mean and deviations from mean

Tiny helper function to calculate the mean and deviations from the mean,
both returned as a list. Works nicely with data.table to calculate a
between and within variable.

## Usage

``` r
meanDeviations(x, na.rm = TRUE)
```

## Arguments

- x:

  A vector, appropriate for the `mean` function.

- na.rm:

  A logical, whether to remove missing or not. Defaults to `TRUE`.

## Value

A list of the mean (first element) and deviations from the mean (second
element).

## Examples

``` r
## simple example showing what it does
meanDeviations(1:10)
#> [[1]]
#> [1] 5.5
#> 
#> [[2]]
#>  [1] -4.5 -3.5 -2.5 -1.5 -0.5  0.5  1.5  2.5  3.5  4.5
#> 

## example use case, applied to a data.table
library(data.table)
d <- as.data.table(iris)
d[, c("BSepal.Length", "WSepal.Length") := meanDeviations(Sepal.Length),
  by = Species]
#>      Sepal.Length Sepal.Width Petal.Length Petal.Width   Species BSepal.Length
#>             <num>       <num>        <num>       <num>    <fctr>         <num>
#>   1:          5.1         3.5          1.4         0.2    setosa         5.006
#>   2:          4.9         3.0          1.4         0.2    setosa         5.006
#>   3:          4.7         3.2          1.3         0.2    setosa         5.006
#>   4:          4.6         3.1          1.5         0.2    setosa         5.006
#>   5:          5.0         3.6          1.4         0.2    setosa         5.006
#>  ---                                                                          
#> 146:          6.7         3.0          5.2         2.3 virginica         6.588
#> 147:          6.3         2.5          5.0         1.9 virginica         6.588
#> 148:          6.5         3.0          5.2         2.0 virginica         6.588
#> 149:          6.2         3.4          5.4         2.3 virginica         6.588
#> 150:          5.9         3.0          5.1         1.8 virginica         6.588
#>      WSepal.Length
#>              <num>
#>   1:         0.094
#>   2:        -0.106
#>   3:        -0.306
#>   4:        -0.406
#>   5:        -0.006
#>  ---              
#> 146:         0.112
#> 147:        -0.288
#> 148:        -0.088
#> 149:        -0.388
#> 150:        -0.688
str(d)
#> Classes ‘data.table’ and 'data.frame':   150 obs. of  7 variables:
#>  $ Sepal.Length : num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
#>  $ Sepal.Width  : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
#>  $ Petal.Length : num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
#>  $ Petal.Width  : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
#>  $ Species      : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ BSepal.Length: num  5.01 5.01 5.01 5.01 5.01 ...
#>  $ WSepal.Length: num  0.094 -0.106 -0.306 -0.406 -0.006 ...
#>  - attr(*, ".internal.selfref")=<externalptr> 

rm(d)
```
