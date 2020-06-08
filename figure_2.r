# FIGURE 2

ER <- read.table("~/Documenti/UniBe/EvolutionaryRescue/sim_200204.txt", header = T)

list_theta <- unique(ER$Theta)
list_z <- unique(ER$z)
list_s <- unique(ER$s)
list_r <- unique(ER$r)
list_zeta <- unique(ER$RatioMigration)
list_beta <- unique(ER$RatioCapacities)

ER$selection_rate <- ER$s / ER$z
ER$environment_shape <- ER$r * ER$Theta * (1 - ER$RatioMigration) * (1 - ER$RatioCapacities) / (ER$RatioMigration * ER$RatioCapacities)

par(mfrow=c(4,3))
par(mar=c(5,5,2,1), oma=c(2,2,0,0))
for (z in list_z){
  for (theta in list_theta){
    params <- (ER$Theta == theta & ER$s == 1.0 & ER$r == 0.1 & ER$z == z)
    ER_plot <- ER[params,]
    plot(ER_plot$migration, ER_plot$Rescue,
         log='x', pch = 19,
         xlab = "", ylab = "",
         ylim = c(0, 1.0), cex.axis = 2.0)
    grid()
    m.theory <- seq(0.0001, 0.11, by=0.0001)
    lines(m.theory, theoretical.formula(m.theory, z, 1.0, 0.1, theta, 0.5, T))
    abline(v = 1.0*z/(1.0-z))
    mtext(bquote(paste(s/z == .(ER_plot$selection_rate), ", ", r*theta == .(ER_plot$environment_shape))), cex=1.7)
  }

}
mtext("migration rate",side=1,line=0,outer=TRUE,cex=1.9)
mtext("Probability of rescue",side=2,line=0,outer=TRUE,cex=1.9,las=0)

