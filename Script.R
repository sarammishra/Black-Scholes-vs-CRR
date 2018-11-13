#Project: Black Scholes model/Binomial Tree for European Call Options
#Author: Saran Mishra
#Date: 03/20/2018

# This is my R Script for the Black-Scholes and Binomial Tree Models


# Goal(s): 
# Use a normal distribution table to calculate the price of the call option
# by the Black-Scholes formula.

# Calculate the price of the same call option by implementing binomial
# tree models with the aid of a software package. Experiment with the
# increasing numbers of steps until your numerical result stabilizes to within
# one cent of the result given by the Black-Scholes formula.

#####################################################################################
# Part 1 -- Black - Scholes Model
#####################################################################################

# Parameters

#Strike Price
Sk = 120 

#Expiration Time
Time = 1

#Annual interest rate 

R = 0.07

#Stock volatility

vol = 0.3

#Initial stock price:

## This is the only USER INPUT ##
### Please enter a number that is 7 integers long
#### Example number: 2944231
ID =                      
ID = ID %% 100

So = 100+(ID/10)

#Black Scholes Formula

#Model: 

BlackScholesModel = function(PriceofCallOption){

              dOne = (((log(So/Sk) + ((annualIr+(vol^2)/2)*Time))/(vol* sqrt(Time))))
              
              dTwo = dOne - (vol * sqrt(Time))

              PriceofCallOption = So*pnorm(dOne) - pnorm(dTwo)*(Sk*exp(-annualIr*(Time)))

              return (PriceofCallOption)
}
BlackScholesModel()


#Extra models for final (ignore for this project)

# 
# BSModelforEuroPut  <- function(S, X, t, r, v)
# {
#   d1 <- (log(S/X) + (r + 0.5 * v^2) * t)/(v * sqrt(t))
#   
#   d2 <- d1 - v * sqrt(t)
#   
#   X * exp(-r * t) * pnorm(-d2)- S * pnorm(-d1)
# }
# 
# BSModelforEuroPut(20,23,9,0.1,0.3)

#Where Type is either put or call 

# DeltaforEuroPut <- function(type, S, X, t, r, v)
# {
#   d1 <- (log(S/X) + (r + 0.5 * v^2) * t)/(v * sqrt(t));
#   if(type > 0)	# put
#     pnorm(d1) -1
#   else			# call
#     pnorm(d1);
#   
# }
# 
# DeltaforEuroPut(1,20,23,9,0.1,0.3)
# 
# 
# GammaforEuroPut <- function(type, S, X, t, r, v)
# {
#   d <- (log(S/X) + (r + 0.5 * v^2) * t)/(v * sqrt(t));
#   top <- exp(-d^2/2)
#   bottom <- S * v * sqrt(2*pi*t)
#   top/bottom
# }
# 
# GammaforEuroPut(1,20,23,9,0.1,0.3)

#####################################################################################
#Part two -- utilize binomial tree model to get same call price under risk neutral p
#####################################################################################

#Using the same parameters listed in the Black Scholes model, we will create a binomial tree model


# So we know that using the C-R-R (Cox, Ross, Rubenstien) model that if a price is at T = 0, then
# S can go up (U) or down (D) at S*(U) where U = e^sigma*sqrt(DeltaTime) and S*(D) where D = e^-sigma*sqrt(DeltaTime)
# Note: sigma = vol

#New Variable has to be Time/Number of Periods or Delta(T)

#Calculations 

#dummyvar which will be later replaced by N (a manual input)
NumberofPeriods = 10

DeltaT = Time/NumberofPeriods

#Build Stock Tree


buildStockTree = function(S, sigma, delta_t, N) {
  tree = matrix(0, nrow=N+1, ncol=N+1)
  u = exp(sigma*sqrt(delta_t))
  d = exp(-sigma*sqrt(delta_t))
  for (i in 1:(N+1)) {
    for (j in 1:i) {
      tree[i,j] = S * u^(j-1) * d^((i-1)-(j-1))
    }
  }
  return(tree)
}

buildStockTree(So, vol, DeltaT, NumberofPeriods)

# Find the Q probability or Q = (R-D)/(U-D)

Qprob = function(r, delta_t, sigma) {
  u = exp(sigma*sqrt(delta_t))
  d = exp(-sigma*sqrt(delta_t))
  return((exp(r*delta_t) - d)/(u-d))
}

# Finally, we build the binomial option tree

binomial_option = function(tree, sigma, delta_t, r, X, type) {
  q = Qprob(r, delta_t, sigma)
  option_tree = matrix(0, nrow=nrow(tree), ncol=ncol(tree))
  if(type == 'put') {
    option_tree[nrow(option_tree),] = pmax(X - tree[nrow(tree),], 0)
  } else {
    option_tree[nrow(option_tree),] = pmax(tree[nrow(tree),] - X, 0)
  }
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      option_tree[i, j] = ((1-q)*option_tree[i+1,j] + q*option_tree[i+1,j+1])/exp(r*delta_t)
    }
  }
  return(option_tree)
}


tree = buildStockTree(S=So, sigma=vol, delta_t=1/2, N=10)

binomial_option(tree, sigma=vol, delta_t=1/10, r=0.07, X=Sk, type='call')

# Put everything all together for a binomial tree call option

binomial_optionANDStock = function(type, sigma, T, r, X, S, N) {
  
                                q = Qprob(r=r, delta_t=T/N, sigma=sigma)
                                
                                tree = buildStockTree(S=S, sigma=sigma, delta_t=T/N, N=N)
                                
                                option = binomial_option(tree, sigma=sigma, delta_t=T/N, r=r, X=X, type=type)
                                
                                return(list(q=q, stock=tree, option=option, BinomialTreeCallprice=option[1,1]))
}

BinomialTreemodel = binomial_optionANDStock(type='call', sigma=vol, T=1, r=0.07, X=Sk, S=So, N=72)

BinomialTreemodel
BlackScholesModel()

#####################################################################################
# Fin.
#####################################################################################
# Final Comments: as we notice, the B-S model gave us a call price of 8.89 and the binomial tree model using C-R-R 
# gives us a Call price of 8.88 when n=72
