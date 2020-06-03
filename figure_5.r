# FIGURE 5

# NB: x-axis values were modified later using inkscape in order to uniform them
# to scientific notation. Colors and font sizes could be different than the version
# present in the paper.

data.rho <- read.table("190115.txt", header=T)
data.inst <- read.table("180903.txt", header=T)

pars.rho <- (data.rho$Theta == 100 & data.rho$s == 0.1 & data.rho$r == 0.9)
filtered.rho <- data.rho[pars.rho,]
pars.inst <- (data.inst$Theta == 100 & data.inst$s == 0.1 & data.inst$r == 0.9 
              & data.inst$RatioMigration==0.5 & data.inst$RatioCapacities==0.5)
filtered.inst <- data.inst[pars.inst,]

par101 <- (filtered.rho$growth==1.01)
rho101 <- filtered.rho[par101,]
par105 <- (filtered.rho$growth==1.05)
rho105 <- filtered.rho[par105,]
par125 <- (filtered.rho$growth==1.25)
rho125 <- filtered.rho[par125,]
par150 <- (filtered.rho$growth==1.5)
rho150 <- filtered.rho[par150,]

plot(filtered.inst$migration, filtered.inst$Rescue, 
     log="x", ylab="Probability of rescue", xlab="Migration rate",
     xlim=c(0.001, 1), ylim=c(0., 1.0), col="purple", pch=19, 
     cex.axis=1.3, cex.lab=1.3)
points(filtered.inst$migration, filtered.inst$Rescue+filtered.inst$Error, 
       pch="-", col="purple")
points(filtered.inst$migration, filtered.inst$Rescue-filtered.inst$Error, 
       pch="-", col="purple")

points(rho101$migration, rho101$Rescue, pch=19, col="royalblue")
points(rho101$migration, rho101$Rescue+rho101$Error, pch="-", col="royalblue")
points(rho101$migration, rho101$Rescue-rho101$Error, pch="-", col="royalblue")

points(rho105$migration, rho105$Rescue, pch=19, col="tomato3")
points(rho105$migration, rho105$Rescue+rho101$Error, pch="-", col="tomato3")
points(rho105$migration, rho105$Rescue-rho101$Error, pch="-", col="tomato3")

points(rho125$migration, rho125$Rescue, pch=19, col="forestgreen")
points(rho125$migration, rho125$Rescue+rho125$Error, pch="-", col="forestgreen")
points(rho125$migration, rho125$Rescue-rho125$Error, pch="-", col="forestgreen")

points(rho150$migration, rho150$Rescue, pch=19, col="orange")
points(rho150$migration, rho150$Rescue+rho101$Error, pch="-", col="orange")
points(rho150$migration, rho150$Rescue-rho101$Error, pch="-", col="orange")

abline(v=2*(1-1/1.01), col="royalblue", lty=2, lwd=2)
abline(v=2*(1-1/1.05), col="tomato3", lty=2, lwd=2)
abline(v=2*(1-1/1.25), col="forestgreen", lty=2, lwd=2)
abline(v=2*(1-1/1.5), col="orange", lty=2, lwd=2)

grid()

legend("topleft", c("instantaneous growth", expression(rho == 1.5), 
                    expression(rho == 1.25), expression(rho == 1.05),
                    expression(rho == 1.01)),
       col=c("purple", "orange", "forestgreen", "tomato3", "royalblue"),
       pch=19, cex=1.3, bty="n")
 
