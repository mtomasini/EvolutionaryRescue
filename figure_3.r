par(mar=c(5, 4.2, 4, 2))
migr <- seq(1e-4, 0.5, by = 1e-4)

# function to transform negative values of a vector in zeros
positive.value <- function(vector){
  for (i in 1:length(vector)){
    if(vector[i] < 0 ){
      vector[i] <- 0
    }
    else {
      next
    }
  }
  return(vector)
}

# formula for rescue from de novo in phase 1
p.dnm.1 <- function(m, z, s, r, theta, zeta){
  kappa <- 10000
  u <- 1/20000
  p.dnm <- positive.value(1-exp((u*(-((r + m*zeta)*(s^2 + s*z + m*s*(1 - zeta) + 2*m*z*(1 - zeta) - m*s*zeta - s*sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2))*
                                        theta) + ((m*(2*s + z)*zeta - z*(s + z + m*(1 - zeta) + sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2)))*
                                                    (r - m*(1 - zeta) + m*zeta - (r - m*(1 - zeta) + m*zeta)/exp((r + m*zeta)*theta) + m*(1 - zeta)*(r + m*zeta)*theta))/(r + m*zeta))*kappa)/
                                  ((r + m*zeta)*sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2))))
  
  return(p.dnm)
}

# formula for rescue from de novo in phase 2
p.dnm.2 <- function(m, z, s, r, theta, zeta){
  kappa <- 10000
  u <- 1/20000
  term.1 <- (-2*u*z*(r - m*(1 - zeta) + m*zeta + exp((r + m*zeta)*theta)*(r + m*(1 - zeta) + m*zeta))*kappa)
  term.2 <- exp((r + m*zeta)*theta)/(r*(r + m*zeta))
  p.res.inf <- positive.value(1-exp((-2*u*z*(r - m*(1 - zeta) + m*zeta + exp((r + m*zeta)*theta)*(r + m*(1 - zeta) + m*zeta))*kappa) / exp((r + m*zeta)*theta)/(r*(r + m*zeta))))

  return(p.res.inf)
}

# formula to calculate contributions of sgv
p.sgv.1phase <- function(m, z, s, r, theta, zeta){
  kappa <- 10000
  u <- 1/20000
  p.res.sgv <- positive.value(1-(((s*sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2) - u*(s*(s + z) + m*(s + 2*z)*(1 - zeta) - m*s*zeta)*kappa +
                                     s*u*sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2)*kappa)*
                                    (-(u*z*(z + m*(1 - zeta) - m*zeta + sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2))*kappa) +
                                       s*(sqrt((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2) - u*z*kappa + 2*m*u*zeta*kappa)))/
                                   (s^2*((s + z + m*(1 - zeta))^2 + 2*m*(-s - z + m*(1 - zeta))*zeta + m^2*zeta^2))))
  
  return(p.res.sgv)
}


theoretical.formula <- function(m, z, s, r, theta, zeta, sgv){
  
  p.no.res.denovo <- (1-p.dnm.1(m, z, s, r, theta, zeta))*(1-p.dnm.2(m, z, s, r, theta, zeta))
  
  if(sgv == TRUE){
    ptot <- 1 - p.no.res.denovo*(1-p.sgv.1phase(m, z, s, r, theta, zeta))
  } 
  else{
    ptot <- 1 - p.no.res.denovo
  }
  
  return(ptot)
}


# formula for asymmetric migration
theory.denovo <- function(m, s1, s2, r, theta, zeta){
  kappa = 10000
  u = 1/20000
  s2 = -s2
  
  p.res.theta <- positive.value(1-exp(-(u*(2*(2*r + m*(-1 + 2*zeta))*(m*(s1 - 2*s1*zeta + 2*s2*zeta) + s1*(s1 - s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))) -
                                             (2*(2*r + m*(-1 + 2*zeta))*(m*(s1 - 2*s1*zeta + 2*s2*zeta) + s1*(s1 - s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))/exp(((2*r + m*zeta)*theta)/2.) +
                                             (2*r + m*zeta)*(-(m*s1*(m + s1)*(-1 + zeta)) + m*s2^2*zeta + m*s2*(-s1 + m*zeta) + m*(s1 - s1*zeta + s2*zeta)*sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta) +
                                                               2*r*(-(m*(s2 + 2*s1*(-1 + zeta) - 2*s2*zeta)) + s2*(-s1 + s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))*theta)*kappa)/
                                        ((2*r + m*zeta)^2*sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))
  
  p.res.inf <- positive.value(1-exp((2*s1*u*(m - 2*r - exp(r*theta + (m*zeta*theta)/2.)*(m + 2*r) - 2*m*zeta)*kappa)/exp(((2*r + m*zeta)*theta)/2.)/(r*(2*r + m*zeta))))
  
  p.res.sgv <- positive.value(1-(((s1*u*(m + s1 - 2*m*zeta + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))*kappa - s2*(sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta) + u*(s1 - 2*m*zeta)*kappa))*
                                    (2*m*s1*u*kappa + s2^2*u*kappa - 2*m*s1*u*zeta*kappa + s2*(-(s1*u*kappa) + m*u*(-1 + 2*zeta)*kappa + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta)*(-1 + u*kappa))))/
                                   (s2^2*((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))
  
  p.nores.denovo <- (1-p.res.theta)*(1-p.res.inf)
  
  ptot <- 1 - p.nores.denovo
  return(ptot)
}


## OUTPUT FIGURES

# change z
plot(migr, theoretical.formula(migr, 0.01, 1.0, 0.25, 200, 0.5, FALSE), 
     log='x', type="l", lty = 4,
     ylim = c(0, 1.0), lwd = 2,
     xlab = "Migration rate", ylab = "Probability of rescue", cex.lab = 1.5, cex.axis = 1.3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 200, 0.5, FALSE), 
      lwd = 2, lty = 3)
lines(migr, theoretical.formula(migr, 0.05, 1.0, 0.25, 200, 0.5, FALSE),
      lwd = 2, lty = 2)
lines(migr, theoretical.formula(migr, 0.1, 1.0, 0.25, 200, 0.5, FALSE),
      lwd = 2, lty = 1)
grid()
legend("topleft", c("z = 0.1","z = 0.05","z = 0.02"," z = 0.01"), lty = c(1, 2, 3, 4), lwd=c(2,2,2,2), bty ="n", cex = 1.3)

# change s
plot(migr, theoretical.formula(migr, 0.02, 0.1, 0.25, 200, 0.5, TRUE), 
     log='x', type="l",
     ylim = c(0.00, 0.45), lwd = 2,
     xlab = "Migration rate", ylab = "Probability of rescue", cex.lab = 1.5, cex.axis = 1.3)
lines(migr, theoretical.formula(migr, 0.02, 0.25, 0.25, 200, 0.5, TRUE), 
      lwd = 2, lty = 2)
lines(migr, theoretical.formula(migr, 0.02, 0.5, 0.25, 200, 0.5, TRUE),
      lwd = 2, lty = 3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 200, 0.5, TRUE),
      lwd = 2, lty = 4)
grid()
legend("topright", c("s = 0.1","s = 0.25","s = 0.5","s = 1.0"), lty = c(1, 2, 3, 4), bty ="n", cex = 1.3, lwd = c(2,2,2,2))

# change r
plot(migr, theoretical.formula(migr, 0.02, 1.0, 0.1, 200, 0.5, FALSE), 
     log='x', type="l",
     ylim = c(0, 0.45), lwd = 2,
     xlab = "Migration rate", ylab = "Probability of rescue", cex.lab = 1.5, cex.axis = 1.3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 200, 0.5, FALSE), 
      lwd = 2, lty = 2)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.5, 200, 0.5, FALSE),
      lwd = 2, lty = 3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.9, 200, 0.5, FALSE),
      lwd = 2, lty = 4)
grid()
legend("topright", c("r = 0.1","r = 0.25","r = 0.5","r = 0.9"), lty = c(1, 2, 3, 4), bty ="n", cex = 1.3, lwd = c(2,2,2,2))

# change theta
plot(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 50, 0.5, FALSE), 
     log='x', type="l", lty = 4,
     ylim = c(0, 0.45), lwd = 2,
     xlab = "Migration rate", ylab = "Probability of rescue", cex.lab = 1.5, cex.axis = 1.3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 200, 0.5, FALSE), 
      lwd = 2, lty = 3)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 500, 0.5, FALSE),
      lwd = 2, lty = 2)
lines(migr, theoretical.formula(migr, 0.02, 1.0, 0.25, 900, 0.5, FALSE),
      lwd = 2, lty = 1)
grid()
legend("topleft", c(expression(theta == 1000), expression(theta == 500), expression(theta == 100), expression(theta == 50)), 
       lty = c(1, 2, 3, 4), bty ="n", cex = 1.3, lwd = c(2,2,2,2))

