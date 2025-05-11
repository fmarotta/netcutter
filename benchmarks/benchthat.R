library(microbenchmark)
library(ggplot2)
library(R.utils)
library(data.table)

# Generate matrices of various size and density and benchmark our
# implementations.

Density <- seq(0.2, 0.95, 0.15)
MatrixSize <- 36^seq(1, 2.5, 0.5)
SampleSize <- 10^(1:3)
Algorithm <- c("nc_randomize_R", "nc_randomize_fast")

Density <- seq(0.2, 0.95, 0.30)
MatrixSize <- 36^(1:2)
SampleSize <- 10^(1:2)
Algorithm <- c("nc_randomize_R", "nc_randomize_fast")

# we need matrices where each row and column has at least one T... difficult to do with low densities.

set.seed(1)
d <- rbindlist(lapply(Density, function(d) {
  rbindlist(lapply(MatrixSize, function(m) {
    values <- sample(c(T, F), m, prob = c(d, 1-d), replace = T)
    h <- matrix(values, nrow = m / 9, ncol = m / 4)
    v <- matrix(values, nrow = m / 4, ncol = m / 9)
    q <- matrix(values, nrow = m / 6, ncol = m / 6)
    rbindlist(lapply(SampleSize, function(s) {
      rbindlist(lapply(Algorithm, function (a) {
        cat(d, m, s, a, "\n", sep = " ")
        mb <- tryCatch(
          withTimeout(
            setDT(microbenchmark(
              h = match.fun(a)(h, s),
              v = match.fun(a)(v, s),
              q = match.fun(a)(q, s),
              times = 10
            ))[, as.list(summary(time)), by = expr],
            timeout = 10),
          error = function(e) {
            message(e, "\n")
            return(NULL)
          }
        )
      }), idcol = "algorithm", fill = T)
    }), idcol = "sample_size", fill = T)
  }), idcol = "matrix_size", fill = T)
}), idcol = "density", fill = T)
# d$matrix_size <- factor(d$matrix_size)

ggplot(d) +
  aes(
    x = MatrixSize[matrix_size],
    y = Median,
    color = SampleSize[sample_size],
    shape = Algorithm[algorithm],
    linetype = Algorithm[algorithm],
  ) +
  geom_point() +
  geom_line(aes(group = paste(sample_size, algorithm)), orientation = "x") +
  facet_grid(rows = vars(expr), cols = vars(Density[density])) +
  scale_color_viridis_c() +
  NULL

