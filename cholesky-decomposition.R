# Matthew Sudmann-Day
# Barcelona GSE Data Science
# Cholesky Decomposition

my.chol <- function(A) {
  # Initialize the return matrix.
  L <- matrix(0.0, nrow(A), nrow(A))
  
  # Handle the top-left corner as a special case.
  L[1,1] <- sqrt(A[1,1])
  
  # If A is bigger than 1x1, populate a lower triangular matrix for return.
  if (nrow(A) > 1) {
    for (i in 2:nrow(A)) {
      for (j in 1:(i-1))
        L[i,j] <- (A[i,j] - sum(L[i,1:j] * L[j,1:j])) / L[j,j]
      L[i,i] <- sqrt(A[i,i] - sum(L[i,1:(i-1)] ^ 2))
    }
  }
  return(L)
}

my.forward.solve <- function(L, b)
{
  # Initialize the return vector.
  x <- rep(0.0, length(b))
  
  # Populate the first item.
  x[1] <- b[1] / L[1,1]
  
  # If we are dealing with more than a 1x1 situation, perform
  # forward substitution to populate the result vector.
  if (length(b) > 1) {
    for (i in 2:length(b)) {
      x[i] <- (b[i] - sum(L[i,1:(i-1)] * x[1:(i-1)])) / L[i,i]
    }
  }
  return(x)
}

my.backward.solve <- function(U, b)
{
  d <- length(b)

  # Initialize the return vector.
  x <- rep(0.0, d)
  
  # Populate the last item.
  x[d] <- b[d] / U[d,d]

  # If we are dealing with more than a 1x1 situation, perform
  # backward substitution to populate the rest of the result vector.
  if (d > 1) {
    for (i in seq(d-1, 1)) {
      x[i] <- (b[i] - sum(U[i,(i+1):d] * x[(i+1):d])) / U[i,i]
    }
  }
  
  return(x)
}

my.solve <- function(A, b) {
  # Solve for x in Ax=b
  # First, decompose A: A=LL'
  L <- my.chol(A)
  
  # LL'x = b
  # Let y=L'x
  # So solve for y in Ly=b
  y <- my.forward.solve(L, b)
  
  # Now solve for x in L'x=y
  x <- my.backward.solve(t(L), y)

  return(x)
}

test <- function(n)
{
  ev = runif(n, 0, 10)
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  A <- Z
  print ("A")
  print (A)
  
  b <- rnorm(n, 1, 10)

  print("b")
  print(b)
  
  my.x <- my.solve(A, b)
  r.x <- solve(A, b)

  print("x")
  print(cbind(my.x,r.x))
}
test(3)



Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

my.solve(Posdef(3), c(4,2,1))
solve(Posdef(3), c(4,2,1))


Posdef(5)

run()

A <- matrix( c(3,2,1,2,4,1,1,1,5) , 3 , 3 )
b <- c(4,2,1)

# solve using cholesky
L <- t(chol(A))

# solve for y
y <- rep(0,3)
y[1] <- b[1]/L[1,1]
y[2] <- ( b[2] -L[2,1]*y[1] )/L[2,2]
y[3] <- ( b[3] -L[3,1]*y[1] -L[3,2]*y[2] )/L[3,3]

U <- t(L)

# solve for x
x <- rep(0,3)
x[3] <- y[3]/U[3,3]
x[2] <- ( y[2] -U[2,3]*x[3] )/U[2,2]
x[1] <- ( y[1] -U[1,2]*x[2] -U[1,3]*x[3] )/L[1,1]

# compare with exact solution
round( cbind( x , solve(A)%*% b, solve(A,b) ) , 3 )





A <- matrix( c(3,2,1,2,4,1,1,1,5) , 3 , 3 )
b <- c(4,2,1)

# solve using cholesky
L <- t(chol(A))

# solve for y
y <- rep(0,3)
y[1] <- b[1]/L[1,1]
y[2] <- ( b[2] -L[2,1]*y[1] )/L[2,2]
y[3] <- ( b[3] -L[3,1]*y[1] -L[3,2]*y[2] )/L[3,3]

L <- t(L)

# solve for x
x <- rep(0,3)
x[3] <- y[3]/L[3,3]
x[2] <- ( y[2] -L[2,3]*x[3] )/L[2,2]
x[1] <- ( y[1] -L[1,2]*x[2] -L[1,3]*x[3] )/L[1,1]

# compare with exact solution
round( cbind( x , solve(A)%*% b ) , 3 )