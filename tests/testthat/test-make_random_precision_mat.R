test_that("make_random_precision_mat works", {
  expect_equal({
    set.seed(123); sum(make_random_precision_mat(eta=.03, p=20, lower = .2, upper = .5, diag = .2))
  }, {
    set.seed(123)
    {g <- igraph::sample_gnp(n=20, p=.03, directed = FALSE)
    omega=as.matrix(igraph::as_adjacency_matrix(g))
    for(h1 in 1:(20-1)){
      for(h2 in (h1+1):20){
        if(omega[h1,h2]!=0){
          temp=runif(1, .2, .5)*sample(c(-1,1),size=1)
          omega[h1,h2]=temp
          omega[h2,h1]=temp
        }
      }
    }
    diag(omega)=rowSums(abs(omega)) + .2}
    sum(omega)
  })
})
