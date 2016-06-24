set.seed(100)
n.col <- 10; n.row = 4
P.mat <- matrix(runif(n.col*n.row/2), nrow = n.row / 2)
P.mat <- rbind(P.mat[1, ], P.mat[1, ], P.mat[2, ], P.mat[2, ])
Q.mat <- matrix(runif(n.col*n.row), ncol = n.col)

full.correction <- fisherInd(Q.mat, P.mat, 1)


vill.1 <- fisherInd(Q.mat[1:2, 1:2], P.mat[1:2, 1:2], 1)
vill.2 <- fisherInd(Q.mat[3:4, 3:4], P.mat[3:4, 3:4], 1)

geom.mean <- function(x) {(prod(x))^(1/length(x))}

geom.mean(full.correction[1:2])
geom.mean(full.correction[2:3])

geom.mean(vill.1)
geom.mean(vill.2)

geom.mean(full.correction[1:2])
geom.mean(full.correction[2:3])

geom.mean(vill.1) * full.correction
geom.mean(vill.2) * full.correction

geom.mean(full.correction[1:2]) / geom.mean(full.correction[2:3])
geom.mean(vill.1) / geom.mean(vill.2)




set.seed(100)
n.col <- 10; n.row = 6
P.mat <- matrix(runif(n.col*n.row / 2), nrow = n.row / 2)
P.mat <- rbind(P.mat[1, ], P.mat[1, ], P.mat[2, ], P.mat[2, ], P.mat[3, ], P.mat[3, ])
Q.mat <- matrix(runif(n.col*n.row), ncol = n.col)

full.correction <- fisherInd(Q.mat, P.mat, 1)

g.mean.1 <- geom.mean(full.correction[1:2])
g.mean.2 <- geom.mean(full.correction[3:4])
g.mean.3 <- geom.mean(full.correction[5:6])


full.correction <- fisherIndfullmat(Q.mat, P.mat)

geom.mean(as.vector(full.correction[1:2, 1:2] ))
geom.mean(as.vector(full.correction[1:2, 3:4] ))
geom.mean(as.vector(full.correction[1:2, 5:6] ))




vill.1 <- fisherInd(Q.mat[1:2, 1:2], P.mat[1:2, 1:2], 1)
vill.2 <- fisherInd(Q.mat[3:4, 3:4], P.mat[3:4, 3:4], 1)
vill.3 <- fisherInd(Q.mat[3:4, 3:4], P.mat[3:4, 3:4], 1)




matrix(1:3, ncol = 3) %*% matrix(1:6, ncol = 2)


# Must come in as a dataframe

library(data.table)

consol.matrix <- function(x) {
  if (!is.data.table(x)) x <- as.data.table(x)
  x.ret <- x[, .(.N), by = names(x)]
  N.ret <- matrix(x.ret$N, ncol = 1)
  x.ret[, N := NULL]
  list(mat = as.matrix(x.ret), freq = N.ret)
}



set.seed(100)
# n.col <- 100; n.row = 40000
n.row.fact <- 1000
real.rows.factor = 4
n.col <- 100; n.row = real.rows.factor; n.row = n.row * n.row.fact
n.real.rows = n.row / real.rows.factor
P.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )
P.mat <- rbind(P.mat[-1, ], P.mat[1, ])
#P.mat <- rbind(P.mat[1, ], P.mat[1, ], P.mat[2, ], P.mat[2, ])
#P.mat <- matrix(runif(n.col*n.row), nrow = n.row )
# Q.mat <- matrix(runif(n.col*n.row), ncol = n.col)
Q.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )


Q.mat.consol <- consol.matrix(Q.mat)
P.mat.consol <- consol.matrix(P.mat)

#library(Cpp)
#Rcpp::sourceCpp("~/git/PQindex/PQindex/src/fisherInd.cpp")



system.time(
fisherIndfaster.ret <- fisherIndfasterold(Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                #Q_freq = t(Q.mat.consol$freq),
                Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
)

#   user  system elapsed
#  0.383   0.016   0.432

#  arma::mat Q_consol,  arma::mat P_consol,
#                           arma::mat Q_freq,    arma::mat P_freq,
#                           arma::uvec Q_ind,     arma::uvec P_ind)


system.time( slow.fisher <- fisherInd(Q.mat, P.mat, 1) )

summary( slow.fisher - fisherIndfaster.ret )

system.time(
fast.fisher <-
  fisherIndfast(Q = Q.mat, P = P.mat,
                Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                Q_freq = t(Q.mat.consol$freq),
                P_freq = t(P.mat.consol$freq ))
)

summary( slow.fisher - fast.fisher )


system.time(
fast.fisher.no.consol <-
  fisherIndfast(Q = Q.mat, P = P.mat,
                Q_consol = Q.mat,
                P_consol = P.mat,
                Q_freq = matrix(1, ncol = nrow(Q.mat)),
                P_freq = matrix(1, ncol = nrow(P.mat)))
)



slow.fisher ; fast.fisher ; fast.fisher.no.consol
slow.fisher - fast.fisher

          [,1]
[1,] 0.9751347
[2,] 0.9643772
[3,] 1.0079894
[4,] 1.1667398



set.seed(100)
n.col <- 3; n.row = 4
Q.mat <- matrix(runif(n.col*4), ncol = n.col, nrow = n.row, byrow = TRUE )
library(PQindex)
testfn(Q_consol = Q.mat, Q_ind = rep(0:3, 3))
testfn(Q_consol = Q.mat, Q_ind = 0:1)
testfn(Q_consol = Q.mat, Q_ind = 1)




# x must come in as a data.table


consol.matrix <- function(x) {
  x <- copy(x)
  x[, y := do.call(paste, .SD)]
  x[, y := as.numeric(factor(y, unique(y), ordered = T)) ]
  setorder(x, y)
  y <- rle(x[, y])$lengths
  x[, y := NULL]
  # Thanks to http://stackoverflow.com/questions/26160079/fast-concise-way-to-generate-ordered-frequency-count-of-unique-matrix-rows
  list(mat = as.matrix(unique(x, by = NULL)), freq = y)
}


library(data.table)

x <- read.csv(textConnection(
  '0,1,1,0
0,1,1,0
1,0,1,0
0,1,0,1
1,0,0,1'),
  header = FALSE)

set.seed(1)
x20 <- x[sample(nrow(x), 20, TRUE), ]
x1M <- x[sample(nrow(x), 1e6, TRUE), ]

#Editing this a bit:
#  https://gist.github.com/mrdwab/b83017e62915586c8f99

funDt <- function(indt, addRN = TRUE, rnName = "rn") {
  if (!is.data.table(indt)) indt <- as.data.table(indt)
  if (isTRUE(addRN)) indt[, (rnName) := sequence(nrow(indt))]
  byCols <- setdiff(names(indt), rnName)
  indt[, list(
    .N, rowid = get(rnName)[1],
    rows = list(get(rnName))),
    by = byCols]
}

system.time( funDt(x) )
system.time( funDt(x20, addRN = FALSE) )
system.time( funDt(x1M) )
system.time( consol.matrix(as.data.table(x)))
system.time( consol.matrix(as.data.table(x20)))
system.time( consol.matrix(as.data.table(x1M)))







