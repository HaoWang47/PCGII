test_that("make_sf_precision_mat works", {
  expect_equal({
    set.seed(111); sum(make_sf_precision_mat(e=2, power=1, p=100, lower = .2, upper =.5, diag=.2))
  }, {set.seed(111); g <- igraph::sample_pa(n=100, power=1, m=2, directed = FALSE)
  omega=as.matrix(igraph::as_adjacency_matrix(g))
  for(h1 in 1:(100-1)){
    for(h2 in (h1+1):100){
      if(omega[h1,h2]!=0){
        temp=runif(1, .2, .5)*sample(c(-1,1),size=1)
        omega[h1,h2]=temp
        omega[h2,h1]=temp
      }
    }
  }
  diag(omega)=rowSums(abs(omega)) + .2
  sum(omega)})
})
