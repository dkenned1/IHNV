

deltat=0.01
TimeSteps=1100

CohortDuration=1

S=numeric(TimeSteps)
I1=numeric(TimeSteps)
I2=numeric(TimeSteps)
SumI1=numeric(TimeSteps)
SumI2=numeric(TimeSteps)


M=100

beta1=0.1
beta2=0.2

mu12=0.01
mu21=0.01

nu1=.1
nu2=.1

alpha1=3
alpha2=6

#S[1]=alpha1/beta1
#I1[1]=M/alpha1
I1[1]=8
I2[1]=0
S[1]=M-I1[1]-I2[1]

SumI1[1]=0
SumI2[1]=0

for (i in 1:(TimeSteps-1))
{
	if ((i %% round(CohortDuration/deltat))!=0)
	{
		S[i+1]= S[i]  + (-beta1*S[i]*I1[i] - beta2*S[i]*I2[i]) *deltat
		I1[i+1]=I1[i] + (beta1*S[i]*I1[i] - mu12*I1[i] + mu21*I2[i] - alpha1*I1[i]) *deltat
		SumI1[i+1]=SumI1[i] + (beta1*S[i]*I1[i]) *deltat
		I2[i+1]=I2[i] + (beta2*S[i]*I2[i] + mu12*I1[i] - mu21*I2[i] - alpha2*I2[i]) *deltat
		SumI2[i+1]=SumI2[i] + (beta1*S[i]*I2[i]) *deltat
	}
	else 
	{
		S[i+1]=M - SumI1[i]*nu1 - SumI2[i]-nu2
		I1[i+1]= SumI1[i]*nu1
		I2[i+1]=SumI2[i]*nu2
		SumI1[i+1]=0
		SumI2[i+1]=0
	}
}


plot(deltat*c(1:TimeSteps), S, type="p", ylim=c(0, max(S, I1, I2)), ylab="Population size", xlab="Cohort")

points(deltat*c(1:TimeSteps), I1, col="blue")
points(deltat*c(1:TimeSteps), I2, col="red")


