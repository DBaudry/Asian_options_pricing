library(randtoolbox)
library(mlmc)

##################Parameters####################

r=0.05 #interest rate
mu=0.05 #drift
sigma=0.25 #vol
T=1 #maturity
k=20 #number of subidivisions
S0=1 #underying price at time 0
K=0.5 #strike price
n=10^4 #number of simulations
alpha=0.05 #for confidence interval
M=2
epsilon = 0.001

u=rnorm(n*k) #we generate n*k random normal distribution


#### Price of the underlying in Black-Scholes model ####

Underlying_BS=function(mu,sigma,t,x,S0){
  u_bs=S0*exp((mu-1/2*sigma^2)*t+sigma*sqrt(t)*x);
  return (u_bs)
}

Underlying_BS(mu,sigma,T,x[1,2],S0)

#### Standard Monte Carlo ####

price_mc=function(mu,r,T,K,k,sigma,S0,n,u){
  t0=Sys.time()
  x=matrix(u,ncol=n,nrow=k)
  z=matrix(S0,ncol=n,nrow=k+1)
  for(i in 2:k+1){
    z[i,] = z[i-1,]*Underlying_BS(mu,sigma,T/k,x[i-1,],1)
  }
  S_asiat=apply(z,2,mean)
  expected_MC=mean((S_asiat-K)*(S_asiat>=K)*exp(-r*T))
  Variance_MC=mean(((S_asiat-K)*(S_asiat>=K)*exp(-r*T)-expected_MC)^2)
  t1=Sys.time()
  return(list(PriceMonte=expected_MC, VarianceMC=Variance_MC,Simulationvector=(S_asiat-K)*(S_asiat>=K)*exp(-r*T),time=t1-t0))
}

price_mc_test<-price_mc(mu,r,T,K,k,sigma,S0,n,u)
price_mc_test$PriceMonte
price_mc_test$VarianceMC
price_mc_test$time
IC_mc=quantile(price_mc_test$Simulationvector,c(alpha/2,1-alpha/2)) #confidence interval

#The variance in the pricing is too high to use the simple Monte Carlo estimator in financial models

###### Variance reduction methods ######

#### Antithetic method ####  

# using the fact that if W is a brownian motion, -W too

Anti_price=function(mu,r,T,K,k,sigma,S0,n,u){
  t0=Sys.time()
  x=matrix(u,ncol=n,nrow=k)
  z=matrix(S0,ncol=n,nrow=k+1)
  y=matrix(S0,ncol=n,nrow=k+1)
  for(i in 2:k+1){
    z[i,]=z[i-1,]*Underlying_BS(mu,sigma,T/k,x[i-1,],1)
    y[i,]=y[i-1,]*Underlying_BS(mu,sigma,T/k,-x[i-1,],1)
  }
  
  S_asiat=apply(z,2,mean)
  S_asiat_anti=apply(y,2,mean)
  q=(S_asiat-K)*(S_asiat>=K)
  s=(S_asiat_anti-K)*(S_asiat_anti>=K)
  j=1/2*(q+s)
  anti_price=mean(j)*exp(-r*T)
  Variance_Anti=mean((j*exp(-r*T)-anti_price)^2)
  t1=Sys.time()
  return(list(AntiPrice=anti_price, Variance_Anti=Variance_Anti,Vect=j*exp(-r*T),time=t1-t0))
}
Anti_price_test=Anti_price(mu,r,T,K,k,sigma,S0,n,u)
Anti_price_test$Variance_Anti
Anti_price_test$AntiPrice
Anti_price_test$time
IC_anti=quantile(Anti_price_test$Vect,c(alpha/2,1-alpha/2)) 


#### Control variable method ####

# Using a control variable Y correlated with the original
# Monte Carlo estimators of X + c*(Y-E(Y))
# We choose the geometric mean as control variate (suggestion from)

price_control=function(r,mu,k,T,K,S0,sigma,n,u){
  t0=Sys.time()
  x=matrix(u,ncol=n,nrow=k)
  z=matrix(S0,ncol=n,nrow=k+1)
  for(i in 2:k+1){
    z[i,]=z[i-1,]*Underlying_BS(mu,sigma,T/k,x[i-1,],1)
  }
  
  var_control=apply(z,2,prod)
  S_asiat=exp(-r*T)*(apply(z,2,mean)-K)*(apply(z,2,mean)>=K)
  var_control=((var_control)^(1/k)-K)*((var_control)^(1/k)>=K)*exp(-r*T)
  model=lm(S_asiat~var_control)
  Vec=S_asiat-model$coefficients[2]*(var_control-mean(var_control))
  Price=mean(Vec)
  variance_control=mean((Vec-Price)^2)
  t1=Sys.time()
  return(list(Price=Price,Vec=Vec,R_squared=summary(model)[8],Variance_control=variance_control,time=t1-t0))
}
price_control_test=price_control(r,mu,k,T,K,S0,sigma,n,u)
price_control_test$Price
price_control_test$R_squared
price_control_test$Variance_control
price_control_test$time
IC_control=quantile(price_control_test$Vec,c(alpha/2,1-alpha/2)) #confidence interval


####Quasi monte carlo #####

QMC=function(mu,sigma,S0,k,K,n,T,r,index){
  t0=Sys.time()
  x=matrix(S0,ncol=n,nrow=k+1)
  if(index==1){
    M=sobol(n,k)
  }
  if(index==2){
    M=halton(n,k)
  }
  
  if(index==3){
    M=torus(n,k)
  }
  
  for(j in 2:k+1){
    x[j,]=x[j-1,]*Underlying_BS(mu,sigma,T/k,qnorm(M[,j-1]),1)
  }
  payoff=exp(-r*T)*(apply(x,2,mean)-K)*(apply(x,2,mean)>=K)
  
  mean_payoff=mean(payoff)
  var_price=mean((payoff - mean_payoff)^2)
  t1=Sys.time()
  return(list(Price_QMC=mean_payoff,Vec=payoff,Var_QMC=var_price,time=t1-t0))
}
Price_Sobol=QMC(mu,sigma,S0,k,K,n,T,r,1)
Price_Sobol$Price_QMC
Price_Sobol$Var_QMC
Price_Sobol$time
IC_Sobol=quantile(Price_Sobol$Vec,c(alpha/2,1-alpha/2))
Price_Halton=QMC(mu,sigma,S0,k,K,n,T,r,2)
Price_Halton$Price_QMC
Price_Halton$Var_QMC
Price_Halton$time
IC_Halton=quantile(Price_Halton$Vec,c(alpha/2,1-alpha/2))
Price_Torus=QMC(mu,sigma,S0,k,K,n,T,r,3)
Price_Torus$Price_QMC
Price_Torus$Var_QMC
Price_Torus$time
IC_Torus=quantile(Price_Torus$Vec,c(alpha/2,1-alpha/2))

plot(cumsum(price_mc_test$Simulationvector)/(1:n),type='l',col='red',xlab='Nombre de Simulations',ylab='Prix',main='Quasi Monte Carlo speed')
lines(cumsum(Price_Sobol$Vec)/(1:n),type='l',col='blue')
lines(cumsum(Price_Halton$Vec)/(1:n),type='l',col='green')
lines(cumsum(Price_Torus$Vec)/(1:n),type='l',col='black')

#### Speed of convergence ####

#plot speed convergence of the three methods (antithetic, control, qmc)
plot(cumsum(price_mc_test$Simulationvector)/(1:n),type='l',col='red',xlab='Nombre de Simulations',ylab='Prix',main=' Monte Carlo, Antithetic, Control')
lines(cumsum(Anti_price_test$Vect)/(1:n),type='l',col='blue')
lines(cumsum(price_control_test$Vec)/(1:n),type='l',col='green')
#legend("bottomright",legend=c(" Prix Monte","Prix Antithétique ","Prix controle"),col=c("red","blue","green"),pch=15)


##### Multi-level Monte-Carlo ####

mlmc_level=function(l,N,mu,r,T,K,sigma,S0,M){
  if(l==0){
    u0=rnorm(N)
    return(as.numeric(price_mc(mu,r,T,K,2,sigma,S0,N,u0)$PriceMonte))}
  u_l = rnorm(N*M^l)
  u_l1 = rnorm(N*M^(l-1))
  P_l = price_mc(mu,r,T,K,M^l,sigma,S0,N,u_l)$PriceMonte
  P_l1 = price_mc(mu,r,T,K,M^(l-1)+1,sigma,S0,N,u_l1)$PriceMonte
  return(as.numeric(P_l-P_l1))
}

mlmc(2, 10, 1000, 0.01, mlmc_level, gamma = 1,mu=mu,r=r,T=T,K=K,sigma=sigma,S0=S0,M=4)

###### Try to implement Multi level monte Carlo method from Giles paper (Oxford, 2006)
## Have to be improved

#### Manual Multi_level ####

Pl=function(mu,r,T,K,M,l,sigma,S0,N){
  u=rnorm(N*M^l)
  x=matrix(u,ncol=N,nrow=M^l)
  z=matrix(S0,ncol=N,nrow=M^l+1)
  if(l==0){
    z[2,]=Underlying_BS(mu,sigma,T,x[1,],S0)
  }
  else{
    for(i in 2:M^l+1){
      z[i,] = z[i-1,]*Underlying_BS(mu,sigma,T/M^l,x[i-1,],1)
    }}
  S_asiat=apply(z,2,mean)
  expected_MC=mean((S_asiat-K)*(S_asiat>=K)*exp(-r*T))
  Variance_MC=mean(((S_asiat-K)*(S_asiat>=K)*exp(-r*T)-expected_MC)^2)
  return(list(PriceMonte=expected_MC,MC_Var=Variance_MC,Simul_Vec=(S_asiat-K)*(S_asiat>=K)*exp(-r*T)))
}

Yl=function(mu,r,T,K,M,l,sigma,S0,N){
  if (l==0){
    return(Pl(mu,r,T,K,M,l,sigma,S0,N))}
  else{
    PriceMonte = Pl(mu,r,T,K,M,l,sigma,S0,N)$PriceMonte-Pl(mu,r,T,K,M,l-1,sigma,S0,N)$PriceMonte
    Simul_Vec = Pl(mu,r,T,K,M,l,sigma,S0,N)$Simul_Vec - Pl(mu,r,T,K,M,l,sigma,S0,N)$Simul_Vec
    MC_Var=mean((Simul_Vec-PriceMonte)^2)
    return(list(PriceMonte=PriceMonte,MC_Var=MC_Var))
  }
}

## Equation (12) in Giles

Nl_calculus=function(M,l,epsilon,V){
  #V Vector containing V0, V1,..., VL
  Nl=trunc(2/epsilon^2*sqrt(V[l+1]/M^l)*sum(sqrt(V*M^l)))
  return(as.numeric(Nl))
}

## Tests : Equation (10) and (11) in Giles

Test = function(Y,M,L,epsilon){
  Test=FALSE
  if(L<2){return(FALSE)}
  if(L>5){return(TRUE)}
  E1=(max(abs(Y[L])/M,abs(Y[L+1]))<1/sqrt(2)*(M-1)*epsilon)
  E2=(abs(Y[L+1]-Y[L]/M)< 1/sqrt(2)*(M^2-1)*epsilon)
  if(E1==TRUE & E2==TRUE){
    Test=TRUE
  }
  return(Test)
}

## MLMC : Check out to create a vector for Yl and their variances

MLMC = function(mu,r,T,K,M,sigma,S0,epsilon){
  L=0
  N_init=10^4
  P0=Yl(mu,r,T,K,M,0,sigma,S0,N_init)
  N0=Nl_calculus(M,0,epsilon,c(P0$MC_Var))
  if(N0>N_init){
    P0=Yl(mu,r,T,K,M,0,sigma,S0,N0)}
  Price_levels=c(P0$PriceMonte)
  Var_levels=c(P0$MC_Var)
  Var_Tests=c(P0$MC_Var)
  N=c(N0)
  while(L<2 | Test(Price_levels,M,L,epsilon)==FALSE){
    L<-L+1
    PL=Yl(mu,r,T,K,M,L,sigma,S0,N_init)
    Var_Tests=c(Var_Tests,PL$MC_Var)
    NL=Nl_calculus(M,L,epsilon,Var_Tests)
    if(NL>N_init){
      PL=Yl(mu,r,T,K,M,L,sigma,S0,NL)}
    Price_levels=c(Price_levels,PL$PriceMonte)
    Var_levels=c(Var_levels,PL$MC_Var)
    N=c(N,NL)
  }
  return(list(levels=Price_levels,Price=sum(Price_levels),Var_list=Var_levels,N=N,L=L))}

m<-MLMC(mu,r,T,K,M,sigma,S0,epsilon)
m$Price
m$L
m$N

price_control_test=price_control(r,mu,k,T,K,S0,sigma,n,u)
price_control_test$Price
price_control_test$R_squared
price_control_test$Variance_control

#### Giles Package ####

Yl1=function(l,N,mu,r,T,K,M,sigma,S0){
  if (l==0){
    return(Pl(mu,r,T,K,M,l,sigma,S0,N))}
  else{
    PriceMonte = Pl(mu,r,T,K,M,l,sigma,S0,N)$PriceMonte-Pl(mu,r,T,K,M,l-1,sigma,S0,N)$PriceMonte
    Simul_Vec = Pl(mu,r,T,K,M,l,sigma,S0,N)$Simul_Vec - Pl(mu,r,T,K,M,l,sigma,S0,N)$Simul_Vec
    return(PriceMonte)
  }
}

mlmc(2, 5, 1000, 0.035, Yl1, gamma = 1,mu=mu,r=r,T=T,K=K,M=M,sigma=sigma,S0=S0)
