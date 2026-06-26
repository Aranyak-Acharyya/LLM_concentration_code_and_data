#CODE FOR OBTAINING TABLE IN REAL DATA ANALYSIS


library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()


library(Matrix)
library(MASS)
library(irlba)
library(igraph)
library(smacof)
library(simcausal)
library(pracma)
library(wordspace)
library(tmvtnorm)

e <- new.env()
e$libs <- c("irlba","Matrix","igraph",
            "simcausal","pracma","wordspace",
            "smacof","MASS","tmvtnorm",
            .libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))


RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

#loading dataset
load("xijk_AA_n10-m2-k5623.RData")


df<-df.aa

print(df)



#declaring parameters
p<-768    #embedding dimension of every vectorized response
N<-10     #total number of models
m<-2      #fixed number of queries
R<-5623   #maximum number of replicates
K<-4     #number of distinct population mean response matrices



#generic function to store an input as a matrix
store_as_matrix<-function(X)
{
  Y<-as.matrix(X)
  
  rownames(Y)<-NULL
  colnames(Y)<-NULL
  
  return(Y)
}


#dissimilarity function between the response matrices of two LLMs
f<-function(x,y)
{
  X<-store_as_matrix(x)
  Y<-store_as_matrix(y)
  
  res<-(m^(-0.5))*(norm(X-Y,type="F"))
  return(res)
}


#function for double centering a dissimilarity matrix
double_centering_func<-function(D)
{
  n<-nrow(D)
  
  H<-diag(n)-(1/n)*matrix(1,n,n)
  
  D2<-D*D
  
  DC<-(-0.5)*H%*%D2%*%H
  
  return(DC)
}



#list of distinct matrices of population mean responses
mu_basic_list<-lapply(1:K,
                     function(i)
                       store_as_matrix(do.call(rbind,lapply(1:m,
                                                            function(j)
                                                              colMeans(do.call(rbind,lapply(1:R, function(k) df$embedding[[(i-1)*m*R+(j-1)*R+k]])))
                       )))
                     )



#an instance depending on n
n<-35

#selecting indices to assign unique population mean response matrices to each LLM
ind_select<-sapply(1:n,
                   function(i)
                     i%%K+1
)



#list of population mean response matrices
mu_true_list<-lapply(1:n, function(i) mu_basic_list[[ind_select[i]]] )
print(mu_true_list)


#population dissimilarity matrix
Delta<-proxy::dist(mu_true_list,method = f,diag=TRUE,upper=TRUE)
Delta<-store_as_matrix(Delta)
print(Delta)

#double-centered population dissimilarity matrix
B<-double_centering_func(Delta)





#finding embedding dimension from scree plot
ee<-eigen(B)$values

x<-1:length(ee)
y<-ee

plot(x,y,type="l")
axis(side = 1, at = x,labels = T)


#setting embedding dimension 
d<-2

#finding population perspectives
B_irlba<-irlba(B,d)
psi<-B_irlba$u%*%(diag(B_irlba$d)^0.5)

print(psi)


#deciding the number of iid replicates
r<-floor(n^(2.75))


sample(1:R,size=r,replace=TRUE)

#list of matrices of sample mean responses
#X_bar_nosample_list<-lapply(1:n,
#                     function(i)
#                       store_as_matrix(do.call(rbind,lapply(1:m,
#                                                            function(j)
#                                                              colMeans(do.call(rbind,lapply(1:r, function(k) df$embedding[[(i%%K)*m*R+(j-1)*R+k]])))
#                       )))
#)


#list of matrices of sample mean responses
X_bar_list<-lapply(1:n,
                   function(i)
                     store_as_matrix(do.call(rbind,lapply(1:m,
                                                          function(j)
                                                            colMeans(do.call(rbind,lapply(sample(1:R,size=r,replace=TRUE), function(k) df$embedding[[(i%%K)*m*R+(j-1)*R+k]])))
                     )))
)





#sample dissimilarity matrix
D<-proxy::dist(X_bar_list,method = f,diag=TRUE,upper=TRUE)
D<-store_as_matrix(D)
print(D)

B_hat<-double_centering_func(D)



#finding sample perspectives by eigendecomposition
B_hat_irlba<-irlba(B_hat,d)
psi_hat<-B_hat_irlba$u%*%(diag(B_hat_irlba$d)^0.5)

err<-procrustes(psi,psi_hat)$d

print(err)





####computing upper bound

#defining function to compute trace of empirical covariance matrix
tr_cov_func<-function(M)
{
 imp_vec<-sapply(1:nrow(M),
        function(i)
          (Norm(M[i,]-colMeans(M),p=2))^2
        )
 
 res<-mean(imp_vec)
  
 return(res)
}




#estimating variability
gamma_basic_mat<-matrix(nrow=K,ncol=m)

for(i in 1:K)
{
  for(j in 1:m)
  {
    
      MM<-store_as_matrix(do.call(rbind,lapply(1:R, function(k) df$embedding[[(i%%K)*m*R+(j-1)*R+k]])))
      
      print(MM)
      
      M<-tr_cov_func(MM)
      
      gamma_basic_mat[i,j]<-M
  }
}

print(gamma_basic_mat)




#defining upper bound as a function 
upper_bd_func<-function(x)
{
  
  lambda_d<-x[1]
  lambda_1<-x[2]
  d<-x[3]
  
  L<-2*sqrt(p)
  
  omega<-1
  
  C<-2*10^(-5)
  
  r<-floor(n^2.75)
  
  Gamma<-max(gamma_basic_mat)
  
  
  
  kappa<-abs(lambda_1)/abs(lambda_d)
  
  print(kappa)
  
  coeff<-((2+sqrt(2))+5*sqrt(2)*d+4*d*sqrt(kappa))/(abs(lambda_d)^0.5)
  
  fnr1<-((2*L*C*omega*sqrt(Gamma))*n)/(sqrt(m*r))
  
  fnr2<-Gamma/r
  
  fnr<-fnr1+fnr2
  
  print(coeff)
  
  res<-coeff*fnr
  
  print(res)
  
  if(fnr<0.5*lambda_d)
  {
    print("Go ahead")
    return(res)  
  }
  if(fnr>=0.5*lambda_d)
  {
    print("do not proceed, condition not met")
    return(res)
  }
  
  print(c(coeff,res,lambda_d))
  
}


lambda_d<-B_irlba$d[d]
lambda_1<-max(B_irlba$d)

print(lambda_d)



upbd<-upper_bd_func(c(lambda_d,lambda_1,d))

print(upbd)


#concatenating important quantities
dec<-c(n,m,d,err,upbd,ifelse(err<upbd,1,0))

print(dec)


#1st iteration
fin_mat<-matrix(,ncol=6)

fin_mat<-rbind(fin_mat,dec)

dff<-data.frame(fin_mat)

print(dff)

save(dff,file="total_real.RData")


#after 1st iteration
load("total_real.RData")

print(dff)

dff<-rbind(store_as_matrix(dff),dec)

save(dff,file="total_real.RData")








  
  
  
  
  
  
  
  

  


