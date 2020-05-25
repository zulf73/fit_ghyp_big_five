# Analysis of Actual Big Five from
#
# https://openpsychometrics.org/_rawdata/
#
# fit a ghyp on five factors then grab the sigma matrix and use svd
# to change to the eigenbasis. Let's look at the distributions
# in this de-correlated basis.
library(hash)

discrete_pmf_1_5<-function(z)
{
  g<-hash()
  count<-0
  for ( b in z){
    if ( is.na(b) ){
      next
    }
    x<-as.character(b)
    b<-as.numeric(b)
    if (b>=1 & b<=5){
      count<-count+1
      if (! has.key(x,g)){
        g[[x]]<-1
      } else {
        g[[x]]<-g[[x]]+1
      }
    }
  }
  for ( j in 1:5){
    g[[as.character(j)]]<-g[[as.character(j)]]/count 
  }
  c(g[["1"]],g[["2"]],g[["3"]],g[["4"]],g[["5"]])
}


kullback_leibler_divergence_from_unifom_1_5<-function( pmf ){
  p<-0.2
  q<-as.vector(pmf)
  abs(sum(p*log((q+0.01)/p)))
}

to_numeric_matrix<-function( str_mtx ){
  n=dim(str_mtx)[1]
  m=dim(str_mtx)[2]
  G<-base::matrix(0,nrow=n,ncol=m)
  for (r in 1:n){
    for (s in 1:m){
      G[r,s]<-as.numeric(str_mtx[r,s])
    }
  }
  G
}

five_by_five_square_freq_table<-function( x, y)
{
  P<-base::matrix(0, nrow=5, ncol=5)
  n<-dim(K)[1]
  for (j in 1:n){
    P[ x[j], y[j] ]<-P[ x[j], y[j]] + 1
  }
  P
}

discrete_pmf_5_5<-function(x, y){
  P<-five_by_five_square_freq_table(x,y)
  P/sum(P)
}

kullback_leibler_divergence_from_unifom_5_5<-function( pmf ){
  p<-1/25
  q<-as.vector(pmf)
  abs(sum(p*log((q+0.01)/p)))
}

dimension_score<-function( x ){
  y<-discrete_pmf_1_5(x)
  kullback_leibler_divergence_from_unifom_1_5(y)
}

dimension_score_2d<-function( x, y ){
  w<-discrete_pmf_5_5(x, y)
  print( w )
  ans<-kullback_leibler_divergence_from_unifom_5_5(w)
  print(paste('Divergence from uniform: ', ans))
  ans
}

create_dimension_score_2d<-function( mtx )
{
  n=dim(mtx)[2]
  ans<-base::matrix(0,nrow=n,ncol=n)
  for (j in 1:n){
    for (k in 1:(j-1)){
      print("----------------------")
      print(paste("Questions:", j, k))
      entry<-dimension_score_2d(mtx[,j],mtx[,k])
      ans[j,k]<-entry
      ans[k,j]<-entry
    }
  }
  ans
}

feature_selection_from_matrix<-function( mtx, threshold )
{
  ans<-c()
  nr<-dim(mtx)[1]
  nc<-dim(mtx)[2]
  for (j in 1:nr){
    for (k in 1:nc){
      if ( mtx[j,k]>threshold){
        ans<-cbind( ans, c(j,k))
      }
    }
  }
  ans
}

dissimilarity_pair_based<-function( x, y, pair_features){
  # assume values are 1-5
  np<-dim(feature_selected_pairs)[2]
  ans<-0
  for (k in 1:np){
    cp<-feature_selected_pairs[,k]
    #print(paste("pair:",cp[1],cp[2]))
    x1<-x[cp[1]]
    x2<-x[cp[2]]
    y1<-y[cp[1]]
    y2<-y[cp[2]]
    if (is.na(x1) || is.na(x2) ||
        is.na(y1) || is.na(y2)){
      next 
    }
    if ( ( x[cp[1]] == y[cp[1]] ) &&
         ( x[cp[2]] == y[cp[2]] ) ){}
    else {
      ans<-ans+1
    }
  }
  ans
}

library(Rcpp)
library(RcppArmadillo)
library(RcppXPtrUtils)

custom_dist<-function( mtx, fpairs ){
  nr<-dim(mtx)[1]
  do<-rep(0,nr*(nr-1)/2)
  for (i in 1:nr){
    for (j in 1:(i-1)){
      val<-dissimilarity_pair_based( mtx[i,],mtx[j,], fpairs)
      do[nr*(i-1) - i*(i-1)/2 + j-i]<-val
    }
  }
  x<-do
  attrs<-list(Size = nr, Labels = dimnames(mtx)[[1L]], Diag = diag, 
       Upper = do, method = "custom", call = match.call(), 
       class = "dist") 
  attr(x,"Size")<-nr
  attr(x,"class")<-"dist"
  attr(x,"method")<-"custom"
  x
}


construct_sparse_feature_selection_mtx<-function( fpairs ){
  q<-dim(fpairs)[2]
  i<-fpairs[1,]
  j<-fpairs[2,]
  x<-ones(2*q,1)
  base::A<-sparseMatrix( i=i, j=j, x=x)
  base::A 
}

sourceCpp('customDist.cpp')
'
class mm_bin {
  static arma::mat A;
  void construct_sparse_feature_selection_mtx( const arma::mat& fpairs ){
    q=dim(fpairs)[2];
    i=fpairs[1,];
    j=fpairs[2,];
    x=ones(2*q,1);
    mm_bin::=sparseMatrix( i=i, j=j, x=x);
  }

  public:
    mm_bin(){}
    mm_bin(arma::mat& fpairs ){
      A = construct_feature_selection_mtx( fpairs );
    }
    double operator()( const arma::mat& x, const arma::mat& y){
      const arma::mat z = x-y;
      const arma::mat zp=arma::abs(arma::sign(A*z))
      return arma::sum( zp );   
    }
    
}

arma::mat mm_bin=arma::mat();
'
binaryFuncPtr <- cppXPtr('mm_bin::operator()',
  depends = c("RcppArmadillo"))


#dd4<-parDist(K[1:100000,], method="custom", func = euclideanFuncPtr)