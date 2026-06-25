###CODE FOR OBTAINING THE TABLE IN SIMULATION SECTION

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


#declaring hyperparameters
p<-75
n<-10
m<-2
R<-floor(n^(3.25))
r<-floor(n^(2.25))



#generic function to store an input as a matrix
store_as_matrix<-function(X)
{
  Y<-as.matrix(X)
  
  rownames(Y)<-NULL
  colnames(Y)<-NULL
  
  return(Y)
}


#declaring function for generating random responses
gen_resp<-function(mu,r)
{
  
  
  L<-lapply(1:nrow(mu),
            function(i)
              colMeans(
                rtmvnorm(r,mean = mu[i,], sigma=Sigma_resp,
                         lower=rep(-2,p),upper=rep(2,p),
                         algorithm="rejection")
              ) )
  
  
  
  res<-store_as_matrix(do.call(rbind,L))
  
  
  return(res)
}




#declaring hyperparameters
L_chol<-matrix(0,p,p)
L_chol[lower.tri(L_chol,diag=TRUE)]<-runif(p*(p+1)/2,min=0.15,max=0.20)
Sigma_query<-L_chol%*%t(L_chol)





#list of LLM-wise population mean of responses
K<-3
mu_basic_list<-lapply(1:K, 
                      function(i) 
                        rtmvnorm(m,mean = rep(0,p), sigma=Sigma_query,
                                 lower=rep(-2,p),upper=rep(2,p),
                                 algorithm="rejection")
                        )









#constructing dispersion matrix of responses
L_chol_resp<-matrix(0,p,p)
L_chol_resp[lower.tri(L_chol_resp,diag=TRUE)]<-runif(p*(p+1)/2,min=0.25,max=0.30)
Sigma_resp<-L_chol_resp%*%t(L_chol_resp)
Sigma_resp<-0.005*Sigma_resp


#list of distinct true population mean matrices of responses, obtained empirically
mu_mid_list<-lapply(1:K,
                    function(i)
                      gen_resp(mu_basic_list[[i]],R)
)



#selecting indices for allocating distinct mean responses to individual LLMs
#ind_select<-sample(1:K,size=n,replace=TRUE)
ind_select<-sapply(1:n,
                   function(i)
                     i%%K+1
                   )




#list of true population mean responses for every LLM
mu_true_list<-lapply(1:n, function(i) mu_mid_list[[ind_select[i]]] )
print(mu_true_list)






#dissimilarity function between the response matrices of two LLMs
f<-function(x,y)
{
  X<-store_as_matrix(x)
  Y<-store_as_matrix(y)
  
  res<-(m^(-0.5))*(norm(X-Y,type="F"))
  return(res)
}





#population dissimilarity matrix
Delta<-proxy::dist(mu_true_list,method = f,diag=TRUE,upper=TRUE)
Delta<-store_as_matrix(Delta)
print(Delta)


#function for double centering a dissimilarity matrix
double_centering_func<-function(D)
{
  n<-nrow(D)
  
  H<-diag(n)-(1/n)*matrix(1,n,n)
  
  D2<-D*D
  
  DC<-(-0.5)*H%*%D2%*%H
  
  return(DC)
}


B<-double_centering_func(Delta)

ee<-eigen(B)$values

x<-1:length(ee)
y<-ee

plot(x,y,type="l")
axis(side = 1, at = x,labels = T)



#choosing embedding dimension for perspectives (from scree-plot elbow)
d<-2

B_irlba<-irlba(B,d)

print(B_irlba)

B_irlba$d


lambda_d<-B_irlba$d[d]
lambda_1<-max(B_irlba$d)

print(c(lambda_d,lambda_1,d))


#defining upper bound as a function 
upper_bd_func<-function(x)
{
  
  lambda_d<-x[1]
  lambda_1<-x[2]
  d<-x[3]
  
  L<-2*p
  
  omega<-1
  
  C<-5*10^(-4)
  
  Gamma<-sum(diag(Sigma_resp))
  
  kappa<-abs(lambda_1)/abs(lambda_d)
  
  coeff<-((2+sqrt(2))+5*sqrt(2)*d+4*d*sqrt(kappa))/(abs(lambda_d)^0.5)
  
  fnr<-((2*L*C*omega*sqrt(Gamma))*n)/(sqrt(m*r))+Gamma/r
  
  res<-coeff*fnr
  
  if(res<0.5*lambda_d)
  {
    print("Go ahead")
    return(res)  
  }
  if(res>=0.5*lambda_d)
  {
    print("do not proceed, condition not met")
  }

  

}




upbd<-upper_bd_func(c(lambda_d,lambda_1,d))

print(upbd)


#finding true perspectives
B<-double_centering_func(Delta)
B_irlba<-irlba(B,d)
psi<-B_irlba$u%*%(diag(B_irlba$d)^0.5)









print(psi)

clusterExport(clust,list("upbd","d","psi","f","mu_true_list"))
clusterExport(clust,list("double_centering_func","gen_resp"))



#finding estimation error on MC samples
L<-foreach(trial=1:50,.combine='c') %dopar%
{
  
  #matrix of sample responses
  X_bar_list<-lapply(1:n,
                     function(i)
                       gen_resp(mu_true_list[[i]],r)
                     )
    
    
  
    
  #sample dissimilarity matrix
  D<-proxy::dist(X_bar_list,method = f,diag=TRUE,upper=TRUE)
  
  D<-store_as_matrix(D)
  
  
  
  
  #doubly-centered sample dissimilarity matrix
  B_hat<-double_centering_func(D)
  
  #finding sample perspectives by eigendecomposition
  B_hat_irlba<-irlba(B_hat,d)
  
  psi_hat<-B_hat_irlba$u%*%(diag(B_hat_irlba$d)^0.5)
  
  print(psi_hat)
  
  
  val<-procrustes(psi,psi_hat)$d
  
  val
}


print(L)


print(mean(L))
print(upbd)



prp<-mean(ifelse(L<=upbd,1,0))




#finding minimum membership
freq_table<-table(ind_select)

print(freq_table)

min_freq<-min(freq_table)

min_rel_freq<-min_freq/n




dec<-c(n,d,lambda_d,mean(L),upbd,prp,min_rel_freq)

print(dec)


stopCluster(clust)


#1st iteration
fin_mat<-matrix(,ncol=7)

fin_mat<-rbind(fin_mat,dec)

df<-data.frame(fin_mat)

print(df)

save(df,file="total_improved_K3.RData")




#after 1st iteration


load("total_improved.RData")

print(df)

df<-rbind(store_as_matrix(df),dec)

save(df,file="total_improved.RData")



#plot 
library(ggplot2)
library(reshape2)
library(latex2exp)

ddf<-df[-1,]

print(ddf)



dff<-data.frame(ddf[,1],ddf[,4],ddf[,5])
colnames(dff)<-c("n","err","upbd")



dffm<-reshape2::melt(dff, id.vars = 'n')
print(dffm)
g<-ggplot(dffm, aes(x=n, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("average error and upper bound")) +
  xlab(TeX("number of models (n)")) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = c("red","orange"),
                      labels=unname(TeX(c(
                        "sample $\\min_{W}||\\hat{\\psi} W -\\psi ||$",
                        "upper bound"))))

ggsave(g,file="LLM_conc_plot_simK4.pdf",
       width = 8, height = 4,
       units = "in", dpi = 600)










