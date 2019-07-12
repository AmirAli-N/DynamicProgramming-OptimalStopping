library(MASS)
library(mvtnorm)
############################################################################
y=c() #a list (for each dose) of vectors of observed responses
n_simulation=30#number of simulations
dose=lapply(1:n_simulation, function(x) c()) #vector of optimal doses
J=11 #number of doses
patient=1000 #number of patients
initial_patient=280
true_sigma=10 #observation variance
mu_0=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) #initial hyper-parameter for KG multivariate normal
sigma_0=diag(100, nrow=J, ncol=J) #initial hyperparameter
mu_n=lapply(1:patient, function(x) c())
sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
theta_estimate=lapply(1:n_simulation, function(x) matrix(NA, nrow=patient, ncol=J)) #a matrix of optimal theta estimates
n_p=30 #number of patients in confirmatory phase-III
c1=1000 #cost of sampling one more patient
c2=1000000 #reward of one unit improvement over placebo
start_time=0 #start time of the algorithm
end_time=0 #end time of the algorithm
M=1000 #number of KG iterations
T=100
placebo=0 #placebo mean response for all dose-responses
decision=c() #vector of optimal stopping for patient k
objvalue=c()
alpha=1
#sigma_string="std_dev=sqrt(1)"
#sigma_string="std_dev=sqrt(10)"
sigma_string="std_dev=10"
#############################################################################
true_theta=c(0.0, 0.07, 0.18, 0.47, 1.19, 2.69, 5, 7.31, 8.81, 9.53, 9.82)
curve_st="sigmoid-significant"
target_dose=10
#############################################################################
##read optimal dose allocation and estimation of optimal thetas from a file##
for (reps in 1:n_simulation){
	dose[[reps]]=scan(file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-doses-",toString(reps),".txt", sep=""), sep="\n")
	theta_estimate[[reps]]=read.table(file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-thetas-",toString(reps),".txt", sep=""), sep="\t", nrow=patient)
}
evolution_eq <-function(mu, sigma, new_y, z_j){paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-thetas-",toString(reps),".txt", sep="")
	e_z = rep(0, J)
	e_z[z_j] <- 1
	sigma_tilde <- (sigma %*% e_z)/sqrt(true_sigma^2+sigma[z_j, z_j])
	new_sigma <- sigma - sigma_tilde %*% t(sigma_tilde)
	Var_yF <- true_sigma^2 + sigma[z_j, z_j]
	new_X <- (new_y - mu[z_j])/sqrt(Var_yF)
	new_mu= mu + sigma_tilde * new_X
	return (list("mu"=new_mu, "sigma"=new_sigma))
}
optimal_stopping<-function(k, y, mu, sigma, z_j){
	U0=0 #utility of abandoning
	U1=0 #utility of continuation
	U2=0 #utility of termination
	##########calculate U2#################################
	temp_theta_sample=rmvnorm(M, mu, sigma)
	ED95=c()
	ED95=apply(temp_theta_sample, 1, function(z){
		if (all(z<0)){
			return(NA)
		}else{
			return (min(which(z>=0.95*max(z))))
		}
	})
	theta_ED=unlist(sapply(1:length(ED95), function(nt)	{
		if (is.na(ED95[nt])==FALSE){
			return(temp_theta_sample[nt, ED95[nt]])
		}
	}))
	y_bar=sapply(1:length(theta_ED), function(nt) rnorm(n_p, theta_ED[nt], true_sigma))
	y1_bar=sapply(1:length(theta_ED), function(nt) rnorm(n_p, placebo, true_sigma))
	B=sapply(1:length(theta_ED), function(nt) pnorm(qnorm(0.999999)-sqrt(n_p)*(mean(y_bar[,nt])-mean(y1_bar[,nt]))/sqrt(2*true_sigma^2), lower.tail=FALSE))
	U2=-c1*n_p+c2*mean(B*sapply(1:length(theta_ED), function(nt) theta_ED[nt]-placebo))
	U2=U2-100000
	#########calculate U1#################################
	y_temp=sapply(1:M, function(nt) rnorm(1, temp_theta_sample[nt, z_j], true_sigma))
	U1_nt=sapply(1:M, function(nt) 	{
		res=evolution_eq(mu, sigma, y_temp[nt], z_j)
		mu_nt=res$mu
		sigma_nt=res$sigma
		temp2_theta_sample=rmvnorm(T, mu_nt, sigma_nt)
		ED_95=c()
		ED95=apply(temp2_theta_sample, 1, function(z) 	{
			if (all(z<0)){
				return(NA)
			}else{
				return (min(which(z>=0.95*max(z))))
			}
		})
		theta_ED=unlist(sapply(1:length(ED95), function(nt2) 	{
			if(is.na(ED95[nt2])==FALSE){
				return(temp2_theta_sample[nt2, ED95[nt2]])
			}
		}))
		y_bar=sapply(1:length(theta_ED), function(nt2) rnorm(n_p, theta_ED[nt2], true_sigma))
		y1_bar=sapply(1:length(theta_ED), function(nt2) rnorm(n_p, placebo, true_sigma))
		B=sapply(1:length(theta_ED), function(nt2) pnorm(qnorm(0.999999)-sqrt(n_p)*(mean(y_bar[,nt2])-mean(y1_bar[,nt2]))/sqrt(2*true_sigma^2), lower.tail=FALSE))
		return (-c1*n_p+c2*mean(B*sapply(1:length(theta_ED), function(nt2) theta_ED[nt2]-placebo)))
	})
	U1=-c1+mean(U1_nt)
	U1=U1-100000
	#decision=min(which(c(U0, U1, U2)==max(U0, U1, U2)))-1
	if (U1<=0 & U2<=0){
		gained_value=0
	} else{
		gained_value=U1-U2
	}
	return (list("KG"=gained_value, "termination"=U2, "continuation"=U1))
}
start_time=Sys.time()
save_decision=lapply(1:n_simulation, function(x) c())
save_objvalue=lapply(1:n_simulation, function(x) c())
for (reps in 1:n_simulation){
	set.seed(reps)
	for (k in 1:patient){
		if (k==1){#at first decision epoch
			#sigma_n[[k]]=sigma_0
			for (i in 1:J){
				for(j in 1:J){
					sigma_n[[k]][i,j]<-100*exp(-alpha*(i-j)^2)
				}
			}
			mu_n[[k]]=mu_0
			z_j=dose[[reps]][k]
			y=c(y, rnorm(1, true_theta[z_j], true_sigma))
		} else{#assignment dose is allocated from a file containing optimal allocated doses
			z_j=dose[[reps]][k]
			if(k>=initial_patient){
				res=optimal_stopping(k, y, mu_n[[k-1]], sigma_n[[k-1]], z_j)
				gained_value=res$KG
				termination=res$termination
				continuation=res$continuation
				if	(gained_value<=0){
					if(termination<=0){
						decision[k-initial_patient+1]=0
						objvalue[k-initial_patient+1]=max(0, continuation, termination)
						break
					}else{
						decision[k-initial_patient+1]=2
						objvalue[k-initial_patient+1]=max(0, continuation, termination)
						break
					}
				}else{
					decision[k-initial_patient+1]=1
					objvalue[k-initial_patient+1]=max(0, continuation, termination)
				}
			}
			##restore data after dose allocation
			##observe the true response of the optimal dose and add it to original data
			y_temp=rnorm(1, true_theta[z_j], true_sigma)
			y=c(y, y_temp) #observing and add the true response of the optimal dose
			res=evolution_eq(mu_n[[k-1]], sigma_n[[k-1]], y[length(y)], z_j)
			mu_n[[k]]<-res$mu
			sigma_n[[k]]<-res$sigma
		}
	}
	#############save results###############################
	save_decision[[reps]]=decision
	save_objvalue[[reps]]=objvalue
	#############reset vectors for iterations###############
	y=c() #a list (for each dose) of vectors of observed responses
	decision=c() #vector of optimal stopping for patient k
	mu_n=lapply(1:patient, function(x) c())
	sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
}
end_time=Sys.time()
start_time-end_time

write.table(as.data.frame(t(matrix(unlist(save_objvalue), nrow=28, ncol=200, byrow=TRUE))), file=paste("C:/Results/correct-selection/",toString(patient),"one_step-objvalue-",curve_st,".csv", sep=""), append=FALSE, row.names=FALSE, col.names=FALSE, sep=",", eol="\n")
write.table(as.data.frame(t(matrix(unlist(save_decision), nrow=28, ncol=200, byrow=TRUE))), file=paste("C:/Results/correct-selection/",toString(patient),"one_step-decision-",curve_st,".csv", sep=""), append=FALSE, row.names=FALSE, col.names=FALSE, sep=",", eol="\n")

