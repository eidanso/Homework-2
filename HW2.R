rm(list=ls())

library("Matrix")

## Converting decimal to binary
n2binary <- function(n){
  n0 <- n
  i <- 1
  b <- c(NA)
  while (n0 > 0){
    rem <- floor(n0/2)
    b[i] <- n0 - rem*2
    n0 <- rem
    i <- i+1
  }
  
  ### Extending the code to truncate the precision of the frractional part of the number
  
  b <- rev(b)
  n0 <- n - n
  k <- 1
  m <- c(NA)
  while (TRUE){
    rem <- n0 - 2**(-k)
    if (rem >= 0) {
      b[i+k] = 1
      n0 <- rem ;
    }
    else {
      b[i+k] = 0
    }
    k <- k+1;
    if (k > 5) {
      break
    }
  }
  return(b)
}
decimal2binary <- n2binary(300.95)


#### Question 2

#Gaussian
A <- matrix(c(4,2,3,  5,-8, 1,  4, 7, -9), 
            byrow =T,
            ncol = 3)
b <- matrix(c(1, 4, 0))
mod_start <- Sys.time()
for (i in range(2000)){
  x <- solve(A, b)
}
mod_end <- Sys.time()
run_time = mod_start - mod_end
print(run_time)

#### Using LU method to confirm if the run_time is  different than the Gaussian
LU <- expand(lu(A))
L <- as.matrix(LU$L)
U <- as.matrix(LU$U)
P <- as.matrix(LU$P)
LHS <- A
RHS <- P%*%L%*%U
LU.mod_start <- Sys.time()
for (i in range(2000)){
  xLU <- solve(U,solve(L,solve(P,b)))
}
LU.mod_end <- Sys.time()
lu.run_time = LU.mod_start - LU.mod_end
#lu.run_time

sprintf("LU Method is betterr than Gausian Method by %5.10f seconds",run_time - lu.run_time)


#### Q3

#### Jacobi Method
A1 <- matrix(c(3, 1, 1, 0, 1, 5, -1, 2, 1, 0, 3, 1, 0, 1, 1, 4), nrow = 4, ncol = 4, byrow = T)
A2 <- matrix(c(2.5, 1, 1, 0, 1, 4.1, -1, 2, 1, 0, 2.1, 1, 0, 1, 1, 2.1), nrow = 4, ncol = 4, byrow = T)
A3 <- matrix(c(2, 1, 1, 0, 1, 3.5, -1, 2, 1, 0, 2.1, 1, 0, 1, 1, 2.1), nrow = 4, ncol = 4, byrow = T)
b1 <- b2 <- b3 <- matrix(c(1, 4, -2, 1))

jac <- function(A, b, x0, tol = 1e-5, maxN=20){ #### tol = tolerance, maxN = Maximum run

  ## Decomposition of A into D (diagonal) and C (Non-Diagonal)
  D <- diag(diag(A))
  C <- A - D
  n <- length(b)
  x0 <- matrix(x0, ncol = 1)
  DinvC <- -solve(D, C) ## Solve the individual formulas before the for loop
  Dinvb <- solve(D, b)
  
  output <- as.data.frame(matrix(NA, nrow = maxN+1,
                                   ncol = 2 + n))
  output[1,] <- c(0, x0, NA)
  for (i in 1:maxN){
    x1 = DinvC%*%x0 + Dinvb
    output[i+1,] <- c(i, x1, norm(x1-x0) / norm(x0))
    
    ## this function checks if the Convergence criteria is met
    if (norm(x1-x0) <= norm(x0) * tol){
      out <- list('x' = x1, 'info' = output)
      return(out)
    }
    x0 <- x1
  }
  warning("The System did not converge")
  return(out)
}

x0 <- matrix(rep(0, length=nrow(A1)), ncol=1)
maxIterations <- 10000

a1.soln <- jac(A1, b1, x0, 1e-5)
a2.soln <- jac(A2, b1, x0, 1e-5)
a3.soln <- jac(A3, b1, x0, 1e-5)
#print(a1.soln)
#print(a2.soln)



### Q4

gSeidel <- function(A, b, x0 , tol, maxIterations = 1000){
  # Decomposition of A into D (diagonal) and C (Non-Diagonal)
  D <- lower.tri(A, diag = TRUE) * A
  C <- A - D
  n <- length(x0)
  
  D_invB <- solve(D, b)
  D_invC <- -solve(D, C)
  
  # Check for convergence:
  B <- -solve(D) %*% C
  rho <- max(abs(eigen(B)$values)) 
  if (rho >= 1){
    warning("This is an Error Message!!") 
  }
  
  for(i in 1:maxIterations){
    x1 <- D_invC %*% x0 + D_invB 
    if (norm(x1 - x0) < tol * norm(x0)){
      results <- list("x" = x1,"iterations" = i) 
      return(results)
      break
    }
    x0 <- x1
  }
  results <- list("x" = x, "iterations")
  return(results)
}

a1.soln2 <- gSeidel(A1, b1, x0, 1e-5)
a2.soln2 <- gSeidel(A2, b1, x0, 1e-5)
a3.soln3 <- gSeidel(A3, b1, x0, 1e-5)
print(a1.soln2)
print(a2.soln2)
