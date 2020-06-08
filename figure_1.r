## FIGURE 1

### simulation of population dynamics in deme 1
png("~/Documenti/UniBe/EvolutionaryRescue/figure_1.png", 
    units = "in", width=7, height=7, res=1000)
par(mfrow=c(2:1))
par(mar = c(2,4.1,3,4))

k = 1
theta = 100
time = 1: 3*theta

wt1 = rep(k*0.98,3*theta)
mut1 = rep(0.02,3*theta)
r = 0.98
z = 0.025

for(i in 1:(2*theta))
{
  wt1[theta + i+1] = wt1[theta + i]*r 
}


for(i in 1:(150))
{
  mut1[150 + i+1] = mut1[150 + i]*(1+z) 
}

wt1 = wt1 + rnorm(300,0,0.01)
mut1 = mut1 + rnorm(300,0,0.01)
mut1[mut1<0.03]=0

plot(wt1,ylim=c(0,1),xlim=c(0,300),type="l",lwd=2,col="dodgerblue4",xlab="Time",axes=F,ylab  ="Population size", cex.lab=1.3)
polygon(c(0,theta,theta,0),c(0,0,1,1),col=rgb(0.1,0.8,0.1,0.3))
polygon(c(theta,2*theta,2*theta,theta),c(0,0,1,1),col=rgb(0.8,0.1,0.1,0.3))
polygon(c(2*theta,3*theta,3*theta,2*theta),c(0,0,1,1),col=rgb(0.8,0.1,0.1,0.3))
lines(mut1,col="darkred",lwd=2)
lines(wt1+mut1,col="black",lwd=2)
axis(1,at=c(0,100,200,300),labels=c("",expression(0),expression(theta),expression(2*theta)), cex.axis = 1.3)
axis(2,at=c(0,1),labels=c("0","K"), cex.axis=1.3)
mtext("Deme 1", side = 4, cex=1.3)
legend(10,0.8,col=c("black","darkred","dodgerblue4"),c("total","mutant","wildtype"),lwd=2,horiz=FALSE)


### simulation of population dynamics in deme 2
par(mar = c(4,4.1,1,4))
wt2 = rep(k*0.98,3*theta)
mut2 = rep(0.02,3*theta)

wt2[100:202]=wt2[100:202]*0.9

for(i in 1:(theta))
{
  wt2[2*theta + i+1] = wt2[2*theta + i]*r
}

m = 0.01

for(i in 1:(50))
{
  mut2[150 + i+1] = mut2[150 + i]*r+ m*mut1[150+i]
}

for(i in 1:100)
{
  mut2[200 + i+1] = mut2[200 + i]*(1+z)
}

wt2 = wt2 + rnorm(300,0,0.01)
mut2 = mut2 + rnorm(300,0,0.01)
mut2[mut2<0.03]=0


## plotting
plot(wt2,ylim=c(0,1),xlim=c(0,300),type="l",lwd=2,col="dodgerblue4",xlab="Time",axes=F,ylab  ="Population size", cex.lab=1.3)
polygon(c(0,2*theta,2*theta,0),c(0,0,1,1),col=rgb(0.1,0.8,0.1,0.3))
polygon(c(2*theta,3*theta,3*theta,2*theta),c(0,0,1,1),col=rgb(0.8,0.1,0.1,0.3))
lines(mut2,col="darkred",lwd=2)
lines(wt2+mut2,col="black",lwd=2)
axis(1,at=c(0,100,200,300),labels=c("",expression(0),expression(theta),expression(2*theta)), cex.axis=1.3)
axis(2,at=c(0,1),labels=c("0","K"), cex.axis=1.3)

mtext("Deme 2", side=4, cex=1.3)
legend(10,0.8,fill=c(rgb(0.1,0.8,0.1,0.3),rgb(0.9,0.1,0.1,0.3)),c("original","deteriorated"))
dev.off()