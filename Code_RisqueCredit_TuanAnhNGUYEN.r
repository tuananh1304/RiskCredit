library(plotly)

rm(list = ls())


#I) Define some useful functions for valuation

#I.1) Function for calculating cumulative payment from cumulative loss
cum_U <- function(L, K_lower = 0, K_upper = 1, notional = 100){
  return( max(0, L - K_lower*notional) - max(0, L - K_upper*notional))
}

#I.2) Function for calculating Delta(U) to estimate the integral with respect to U
delta_U <- function(cum_U){
  return(diff(cum_U))
} 




#I)Hawkes process simulaiton

#initial condition
lambda0 = 1
a = 0.9
kappa = 1
delta = 1

T <- numeric()
T[1] = 0
L = numeric()
L[1] = 0
lambda = numeric()
lambda[1] = lambda0

Maturity <- 100


#Simulate the i arrival time
while(tail(T,1) < Maturity){
  D_i =  1 + kappa * log(runif(1)) / (tail(lambda,1) - a)
  S_i2 <- -1/a * log(runif(1))
  
  if ( D_i < 0){
    S_i <- S_i2
  } else {
    S_i1 <- - 1/kappa * log(D_i)
    S_i <- min(S_i1, S_i2)
  }
  
  T_i <- tail(T,1) + S_i
  T <- c(T, T_i)
  
  #generate loss
  
  L_i <- runif(1, min = 0.25,max = 0.75)
  #luu vao vector
  
  L <- c(L,L_i)
  
  #Save lambda
  
  lambda_i_minus <- (tail(lambda,1) - a) * exp( - kappa * (tail(T,1) - tail(T,2)[1])) + a
  
  lambda_i_plus <- lambda_i_minus +  L_i * delta
  
  lambda <- c(lambda, lambda_i_plus)
}

T <- T[-length(T)]
L <- L[-length(L)]


lambda <- lambda[-length(lambda)]
lambda
length(T)
length(L)
plot.new()

plot_ly() %>%
  add_trace(x = T, y = lambda, type = 'scatter', mode = 'lines', name = 'intensity')
hist(T, breaks = c(0:(Maturity)))
axis(side = 2, labels = TRUE, tick = TRUE )
par(new = TRUE)
plot(T, cumsum(L), xaxt = "n", yaxt = "n", type = "l", col = "red", ylab = "")
axis(side = 4, tick = TRUE, labels = FALSE, las = 1)
mtext("Cumulative Loss", side = 4, col = "red")





#II) Model Valuation
#Explicit formula for expected value of lambda t and N t

#expected value for lambdaT
# lambda0 = 1
# a = 0.9
# kappa = 0.5
# delta = 1
lambda_mc <- function(lambda0 = 1, a = 0.9, kappa = 0.5, delta = 1, t = 10, N_sim = 1000){
  sum_lam <-  0
  sum_N <- 0
  
  #Monte-Carlo Simulation
  for (i in (c(1:N_sim))){
    T <- numeric()
    T[1] = 0
    L = numeric()
    L[1] = 0
    lambda = numeric()
    lambda[1] = lambda0
    #Simulate the i arrival time
    
    while( tail(T,1) < t){
      D_i =  1 + kappa * log(runif(1)) / (tail(lambda,1) - a)
      S_i2 <- -1/a * log(runif(1))
      
      if ( D_i < 0){
        S_i <- S_i2
      } else {
        S_i1 <- - 1/kappa * log(D_i)
        S_i <- min(S_i1, S_i2)
      }
      
      T_i <- tail(T,1) + S_i
      T <- c(T, T_i)
      
      #generate loss
      
      L_i <- runif(1, min = 0.35,max = 0.85)
      #luu vao vector
      
      L <- c(L,L_i)
      
      #luu lambda
      
      lambda_i_minus <- (tail(lambda,1) - a) * exp( - kappa * (tail(T,1) - tail(T,2)[1])) + a
      
      lambda_i_plus <- lambda_i_minus +  L_i * delta
      
      lambda <- c(lambda, lambda_i_plus)
      
    }
    
    T <- T[-length(T)]
    L <- L[-length(L)]
    lambda <- lambda[-length(lambda)]
    lambda_T <- (tail(lambda,1) - a) * exp( - kappa * (t - tail(T,1))) + a
    N_T <- length(T) - 1 
    sum_lam <- sum_lam + lambda_T
    sum_N <- sum_N + N_T
  }
  
  return(c(sum_lam/N_sim, sum_N/N_sim))
  
}

# result after 10000 simulation

mc_result <- sapply(c(1:10), function(x) lambda_mc( t = x, N_sim = 10000))
mc_result
mc_lambda <- mc_result[1,]
elambda <- function(t = 10, kappa = 0.5, c = 0.9, lambda0 = 1, delta = 1 ){
  
  l <- 0.6 # mean of loss
  mu <- delta * l - kappa
  
  mean <- exp( mu * t) * ( kappa * c / mu + lambda0) - kappa * c / mu
  return(mean)
  
}

ex_lambda <- sapply(c(1:10), function(x) elambda(t = x))
ex_lambda

error_lambda <- (mc_lambda - ex_lambda/ex_lambda)


test_lambda <- cbind(mc_lambda, ex_lambda, error_lambda)




#Calculation moment of N_t

eN_t <- function(t = 10, kappa = 0.5, c = 0.9, lambda0 = 1, delta = 1 ){
  
  l <- 0.6 # mean of loss
  mu <- delta * l - kappa
  c_1 <- (kappa * c + mu * lambda0) / (mu ^ 2)
  c_2 <- - kappa * c / mu
  
  mean_N <- c_1 * ( exp(mu * t) - 1) + c_2 * t
  return(mean_N)
  
}



ex_N_t <- sapply(c(1:10), function(x) eN_t(t = x))
ex_N_t

mc_N_t <- mc_result[2,]

error_N <- (mc_N_t - ex_N_t)/ex_N_t
test_N_t <- cbind(mc_N_t, ex_N_t, error_N)

test_N_t


test_table <- cbind(test_lambda, test_N_t)

write.csv(test_table, file = 'valid.csv')

#III) CDX Valuation

cdx_valuation <- function(lambda0 = 0.6, a = 0.5, kappa = 2, delta = 1, r_d = 0.05, n = 100, Maturity = 10, N_sim = 1000){
  T <- numeric()
  T[1] = 0
  L = numeric()
  L[1] = 0
  lambda = numeric()
  lambda[1] = lambda0
  
  day <- seq(0, Maturity, 1/365.25)
  day_end <- day[-1]
  day_start <- day[-length(day)]
  
  
  prem_period <- seq(0, 10, 1/4)
  prem_period <- prem_period[-length(prem_period)]
  
  U <- rep(0, length(day) - 1)
  I <- rep(0, length(prem_period))
  
  U_T <- 0
   
  t <- Maturity
  #Monte-Carlo Simulation
  for (i in (c(1:N_sim))){
    T <- numeric()
    T[1] = 0
    L = numeric()
    L[1] = 0
    lambda = numeric()
    lambda[1] = lambda0
    #Simulate the i arrival time
    
    while( tail(T,1) < t){
      D_i =  1 + kappa * log(runif(1)) / (tail(lambda,1) - a)
      S_i2 <- -1/a * log(runif(1))
      
      if ( D_i < 0){
        S_i <- S_i2
      } else {
        S_i1 <- - 1/kappa * log(D_i)
        S_i <- min(S_i1, S_i2)
      }
      
      T_i <- tail(T,1) + S_i
      T <- c(T, T_i)
      
      #generate loss
      
      L_i <- runif(1, min = 0.35,max = 0.85)
      #luu vao vector
      
      L <- c(L,L_i)
      
      #luu lambda
      
      lambda_i_minus <- (tail(lambda,1) - a) * exp( - kappa * (tail(T,1) - tail(T,2)[1])) + a
      
      lambda_i_plus <- lambda_i_minus +  L_i * delta
      
      lambda <- c(lambda, lambda_i_plus)
      
    }
    
    T <- T[-length(T)]
    L <- L[-length(L)]
    
    #III.1 Calculating the present value of all payment
    
    #Generate the cumulative payment U and delta U
    cum_L <- cumsum(L)
    cummulative_U <- sapply(cum_L, function(x) cum_U(x, notional = n))
    delta_U1 <- delta_U(cummulative_U)
    
    
    vec_integral <- rep(0,length(day)- 1)
    for ( i in c(2:length(T))){
      vec_integral <- vec_integral + (day_end > T[i] & day_start < T[i]) * delta_U1[i-1]  
    }
    
    
    U <- U + vec_integral
    
    U_T <- U_T + tail(cummulative_U,1)
    
    
    #III.2 Calculating the present value of Premium
    #In case of pricing of CDS index, the upfront payment is not required, so I calculate only the spread.
    
    #Calculating the coefficient c so that: S * c = D_t (with S = spread and D_t = Present value of payment)
    
    #For simplicity, the quarterly premium is taken into account.
    
    #Number of default
    
    jump <- rep(1, length(T) -1)
    N_t <- c(0,cumsum(jump))
    I_t <- rep(n, length(N_t)) - N_t
    
    vec_prem <- rep(0, length(prem_period))
    
    for ( i in c(1: length(T))){
      if (i < length(T)){
        vec_prem <- vec_prem + (prem_period >= T[i] & prem_period < T[i+1]) * I_t[i]
      }
      else if ( i == length(T)){
        vec_prem <- vec_prem + (prem_period > T[i]) * I_t[i]
      }
      
    }
    
    I <- I + vec_prem
  }
  
  
  mean_U_T <- U_T / N_sim
  mean_U <- U / N_sim
  mean_I <- I/ N_sim
  
  
  # Function for calculate the present value
  act_factor <- function(t, r = r_d){
    return( exp( -  t * r ))
  }
  
  # Expected Payment Calculation
  
  U_actual <- sapply(day_end, act_factor )
  
  Payment <- act_factor(t = Maturity, r = r_d) * mean_U_T + r_d * sum(U_actual * mean_U)
  
  #Expected coefficient of spread
  
  I_actual <- sapply(prem_period, act_factor)
  
  c <- sum(I_actual * 1/4 * mean_I)
  
  
  #Calculate the spread 
  
  spread <- Payment / c
  
  return(spread)
}



cdx_valuation()

#IV) Sensitivity Testing

# Define the range of parameters to test
m_test <- seq(1, 10, by = 1)
d_test <- seq(0, 2, by = 0.5)

# Initialize a matrix to store the results
results <- matrix(0, nrow = length(m_test), ncol = length(d_test))

# Perform grid testing
for (i in 1:length(m_test)) {
  for (j in 1:length(d_test)) {
    # Set the parameter values
    param1 <- m_test[i]
    param2 <- d_test[j]
    
    # Run your simulation or calculation using the parameter values
    
    # Store the result in the results matrix
    results[i, j] <- cdx_valuation( delta = param2, Maturity = param1, N_sim = 5000) # Store your result here
  }
}

# Print the results
print(results)


write.csv(results, file = 'sensi2.csv')

plot_ly() %>%
  add_trace(x = m_test, y = results[,1], type = 'scatter', mode = 'lines', name = 'delta = 0') %>%
  add_trace(x = m_test, y = results[,2], type = 'scatter', mode = 'lines', name = 'delta = 0.5') %>%
  add_trace(x = m_test, y = results[,3], type = 'scatter', mode = 'lines', name = 'delta = 1') %>%
  add_trace(x = m_test, y = results[,4], type = 'scatter', mode = 'lines', name = 'delta = 1.5') %>%
  add_trace(x = m_test, y = results[,5], type = 'scatter', mode = 'lines', name = 'delta = 2') %>%
  layout(title = 'Annual Spread', xaxis = list(title = 'Delta'), yaxis = list(title = 'Spread'))




#Test sensitivity 2

delta_test <- seq(0,1, by = 0.1)
k_test <- seq(0.5, 2, by = 0.5)

# Initialize a matrix to store the results
results1 <- matrix(0, nrow = length(delta_test), ncol = length(k_test))

# Perform grid testing
for (i in 1:length(delta_test)) {
  for (j in 1:length(k_test)) {
    # Set the parameter values
    param1 <- delta_test[i]
    param2 <- k_test[j]
    
    # Run your simulation or calculation using the parameter values
    
    # Store the result in the results matrix
    results1[i, j] <- cdx_valuation(delta = param1, kappa = param2, N_sim = 5000) # Store your result here
  }
}

# Print the results
print(results1)
write.csv(results1, file ='result2.csv')

plot_ly() %>%
  add_trace(x = delta_test, y = results1[,1], type = 'scatter', mode = 'lines', name = 'kappa = 0.5') %>%
  add_trace(x = delta_test, y = results1[,2], type = 'scatter', mode = 'lines', name = 'kappa = 1') %>%
  add_trace(x = delta_test, y = results1[,3], type = 'scatter', mode = 'lines', name = 'kappa = 1.5') %>%
  add_trace(x = delta_test, y = results1[,4], type = 'scatter', mode = 'lines', name = 'kappa = 2') %>%
  layout(title = 'Annual Spread', xaxis = list(title = 'Delta'), yaxis = list(title = 'Spread'))

cdx_valuation(delta = 0, kappa = 1)
