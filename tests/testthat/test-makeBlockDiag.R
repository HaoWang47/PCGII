test_that("makeBlockDiag works", {
  expect_equal({
    set.seed(1111); sum(makeBlockDiag(blocksize = 4, p = 20, min.beta = .3, max.beta = .9))
  }, {
    set.seed(1111);  reps=20/4
    S=list()
    for (i in 1:reps) {
      bd=matrix(runif(1, .3, .9), 4, 4)
      diag(bd)=runif(1,1,1.25)
      while(!corpcor::is.positive.definite(bd)){diag(bd)=diag(bd)+0.01}
      S[[i]]=bd
    }
    sum(as.matrix(Matrix::bdiag(S)))
  })
})
