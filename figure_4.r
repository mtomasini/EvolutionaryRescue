options(digits=10)
# COMPARISONS MODEL WITH ASYMMETRIC QUANTITIES
kappa = 10000
u = 1/(2*kappa)

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

# theoretical function for symmetric model
theory.sym <- function(m,s1,s2,r,theta){  
  pr1 <- 1- exp(-(s1*u*(m^2*theta + 2*m*r*theta + ((m + 2*r)*s2*(-s1 + sqrt(m^2 + (s1 - s2)^2) + s2)*theta)/s1 +
                        (s1 + sqrt(m^2 + (s1 - s2)^2) - s2 + (m*s2)/s1)*((4*(1 - exp(-((m + 2*r)*theta)/2.))*r)/(m + 2*r) + m*theta))*kappa)/((m + 2*r)*sqrt(m^2 + (s1 - s2)^2)))
  pr2 <- 1 - exp((-4*(m + r + r/exp(((m + 2*r)*theta)/2.))*s1*u*kappa)/(r*(m + 2*r)))
  pr1 <- positive.value(pr1)
  pn2 <- positive.value(pr2)
  pres <- 1 - (1-pr1)*(1-pr2)
  
  return(pres)
}

# theoretical function for asymmetric carrying capacities
theory.kappa <- function(m,s1,s2,r,theta,beta){    
  pr1 <- 1-exp(-(u*(-2*(m + 2*r)^2*(m*s1 + s2*(-s1 + sqrt(m^2 + (s1 - s2)^2) + s2))*(-1 + beta)*theta*kappa +
                    2*(s1*(s1 + sqrt(m^2 + (s1 - s2)^2) - s2) + m*s2)*(4*r*beta + 2*m*(-1 + 2*beta) + (2*m - 4*(m + r)*beta)/exp(((m + 2*r)*theta)/2.) - m*(m + 2*r)*(-1 + beta)*theta)*kappa))/((m + 2*r)^2*sqrt(m^2 + (s1 - s2)^2)))

  pr2 <- 1-exp((-2*s1*u*(-2*m*(1 - beta)*kappa + 4*exp((m*theta)/2. + r*theta)*(m + r)*(1 - beta)*kappa +
                           2*m*beta*kappa + 4*r*beta*kappa))/exp(((m + 2*r)*theta)/2.)/(r*(m + 2*r)))

  pr1 <- positive.value(pr1) 
  pr2 <- positive.value(pr2)
  pres <- 1 - (1-pr1) * (1-pr2)
  
  return(pres)
}


# theoretical function for asymmetric migration
theory.migr <- function(m, s1, s2, r, theta, zeta){
  pr1 <- 1-exp(-(u*(2*(2*r + m*(-1 + 2*zeta))*(m*(s1 - 2*s1*zeta + 2*s2*zeta) + s1*(s1 - s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))) -
                     (2*(2*r + m*(-1 + 2*zeta))*(m*(s1 - 2*s1*zeta + 2*s2*zeta) + s1*(s1 - s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))/exp(((2*r + m*zeta)*theta)/2.) +
                     (2*r + m*zeta)*(-(m*s1*(m + s1)*(-1 + zeta)) + m*s2^2*zeta + m*s2*(-s1 + m*zeta) + m*(s1 - s1*zeta + s2*zeta)*sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta) +
                                                    2*r*(-(m*(s2 + 2*s1*(-1 + zeta) - 2*s2*zeta)) + s2*(-s1 + s2 + sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta))))*theta)*kappa)/
                ((2*r + m*zeta)^2*sqrt((m + s1 - s2)^2 + 4*m*(-s1 + s2)*zeta)))
  pr2 <- 1-exp((2*s1*u*(m - 2*r - exp(r*theta + (m*zeta*theta)/2.)*(m + 2*r) - 2*m*zeta)*kappa)/exp(((2*r + m*zeta)*theta)/2.)/(r*(2*r + m*zeta)))
  pr1 <- positive.value(pr1) 
  pr2 <- positive.value(pr2)
  pres <- 1 - (1-pr1) * (1-pr2)
  return(pres)
}

test <- function(m, s1, s2, r, theta, zeta){
  exp((2*s1*u*(m - 2*r - exp(r*theta + (m*zeta*theta)/2.)*(m + 2*r) - 2*m*zeta)*kappa)/exp(((2*r + m*zeta)*theta)/2.)/(r*(2*r + m*zeta)))
}

# define parameters
adv <- 0.02
disadv <- -0.5
stress <- 0.5
epoch <- 100
ratio.kappa <- 0.5
ratio.migr <- 0.5
plot.limit <- 0.3

migr <- seq(0.0001,1,by=0.0001)
y1 <- theory.sym(migr, adv, disadv, stress, epoch)
y2 <- theory.kappa(migr, adv, disadv, stress, epoch, ratio.kappa)
y3 <- theory.migr(migr, adv, disadv, stress, epoch, ratio.migr)


## Differences with carrying capacities

png("Differences_Kappa.png", width = 800, height = 600)
par(mar=c(4.5,5,1,1))
plot(migr, y1, type="l",
     ylim=c(0,0.15),log='x', cex.axis = 2.,
     xlab="migration",ylab="Probability of rescue", cex.lab = 2.5,
     lty=1, lwd=2, col="gray50")
lines(migr, theory.kappa(migr, adv, disadv, stress, epoch, 0.1), lty=2, lwd=2)
lines(migr, theory.kappa(migr, adv, disadv, stress, epoch, 0.9), lty=3, lwd=2)
grid()
legend("topleft", c("symmetric", "asym, k2 larger", "asym, k2 smaller"), lty=c(1, 2, 3), lwd=2, col=c("gray50", "black", "black"), bty="n", cex=2.)
dev.off()

## Differences with asymmetric migration
png("Differences_migr.png", width = 800, height = 600)
par(mar=c(4.5,5,1,1))
plot(migr, y1, type="l",
     ylim=c(0,0.35),log='x', cex.axis = 2.,
     xlab="migration",ylab="Probability of rescue", cex.lab=2.5,
     lty=1, lwd=2, col="gray50")
lines(migr, theory.migr(migr, adv, disadv, stress, epoch, 0.1), lty=2, lwd=2)
lines(migr, theory.migr(migr, adv, disadv, stress, epoch, 0.9), lty=3, lwd=2)
grid()
legend("topleft", c("symmetric", "asym, m21 larger", "asym, m21 smaller"), lty=c(1, 2, 3), lwd=2, col=c("gray50", "black", "black"), bty="n",cex=2.5)
dev.off()
