library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)
library(plyr)
################################################################################
y=c() #a list (for each dose) of vectors of observed responses
n_simulation=10 #number of simulations
dose=lapply(1:n_simulation, function(x) c()) #vector of optimal doses
J=11 #number of doses
patient=1000 #number of patients
initial_patient=20
true_sigma=1 #observation variance
mu_0=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) #initial hyper-parameter
sigma_0=diag(100, nrow=J, ncol=J) #initial hyperparameter
mu_n=lapply(1:patient, function(x) c())
sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
theta_estimate=lapply(1:n_simulation, function(x) matrix(NA, nrow=patient, ncol=J)) #a matrix of optimal theta estimates
m_0=0
m_n=rep(m_0, patient)
s_0=10
s_n=rep(s_0, patient)
df_star=c()
n_p=30 #number of patients in confirmatory phase-III
c1=1000 #cost of sampling one more patient
c2=1000000 #reward of one unit improvement over placebo
M=1000 #number of KG iterations
n_exp=1000 #number of generated sample paths in KG
start_time=0 #sart time of the algorithm
end_time=0 #end time of the algorithm
decision=c() #vector of optimal stopping for patient k
#################################################################################
true_theta=c(0.0, 0.07, 0.18, 0.47, 1.19, 2.69, 5, 7.31, 8.81, 9.53, 9.82)
curve_st="sigmoid-significant"
target_dose=10
sigma_string="std_dev=10"
##################################################################################
n_grid=20 #number of grid points on each axis
max_m=20
min_m=0
max_s=10
min_s=0
#calculate grid length for each patient across all experiments
length_m=max_m-min_m
length_s=max_s-min_s
#calculate cell length for patient's grid across all experiments
grid_length_m=length_m/n_grid
grid_length_s=length_s/n_grid
seq_m=seq(min_m, max_m, grid_length_m)
seq_s=seq(min_s, max_s, grid_length_s)
###############Vector spaces#############################
#read optimal dose allocation and estimation of optimal thetas from a file
for (reps in 1:n_simulation){
	dose[[reps]]=scan(file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-doses-",toString(reps),".txt", sep=""), sep="\n")
	theta_estimate[[reps]]=read.table(file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-thetas-",toString(reps),".txt", sep=""), sep="\t")
}
evolution_eq <-function(mu, sigma, new_y, z_j){
	e_z = rep(0, J)
	e_z[z_j] <- 1
	sigma_tilde <- (sigma %*% e_z)*(1/sqrt(true_sigma+sigma[z_j, z_j]))
	new_sigma <- sigma - sigma_tilde %*% t(sigma_tilde)
	Var_yF <- true_sigma + sigma[z_j, z_j]
	new_X <- (new_y - mu[z_j])/sqrt(Var_yF)
	new_mu= mu + sigma_tilde * new_X
	return (list("mu"=new_mu, "sigma"=new_sigma))
}
optimal_stopping<-function(k, y, m_n, s_n, mu, sigma, z_star){
#define a list for every experiment containing vectors for all patients
	m_n_i=list()
	s_n_i=list()
	z_n_i=list()
	#create a posterior based on observation up to patient k
	myCluster<-makeCluster(7)#assign six core to the cluster
	registerDoParallel(myCluster)
	#registerDoSNOW(myCluster)
	res_for<-foreach(i=1:n_exp, .packages="mvtnorm", .export=c("reps", "k", "m_n", "s_n", "mu", "sigma", "dose", "patient", "evolution_eq", "J", "M", "n_p", "true_sigma", "z_star"), .verbose=FALSE)%dopar%{
		z_ni=rep(z_star, patient)
		m_ni=rep(m_n, patient)
		s_ni=rep(s_n, patient)
		mu_ni=lapply(1:patient, function(x) mu)
		sigma_ni=lapply(1:patient, function(x) matrix(sigma, nrow=J, ncol=J))
		for(n in k:(patient-1)){
			z_j=dose[[reps]][n-1]
			temp_theta_sample<-rmvnorm(M, mu_ni[[n-1]], sigma_ni[[n-1]])
			temp_theta_estimate<-apply(temp_theta_sample, 2, mean)
			theta_ED=0
			if (all(temp_theta_estimate<=0)){
				theta_ED=temp_theta_estimate[1]
			}else{
				theta_ED=temp_theta_estimate[min(which(temp_theta_estimate>=0.95*max(temp_theta_estimate)))]
			}
			df_star=rnorm(n_p, theta_ED, sqrt(2*true_sigma^2))
			s_ni[n+1]=1/((1/s_ni[n])+(n_p/(2*true_sigma^2)))
			m_ni[n+1]=s_ni[n]*((m_ni[n]/s_ni[n])+(sum(df_star)/(2*true_sigma^2)))
			y_temp=rnorm(1, temp_theta_sample[z_j], true_sigma)
			res=evolution_eq(mu_ni[[n-1]], sigma_ni[[n-1]], y_temp, z_j)
			mu_ni[[n]]<-res$mu
			sigma_ni[[n]]<-res$sigma
			temp_theta_sample<-rmvnorm(M, mu_ni[[n]], sigma_ni[[n]])
			temp_theta_estimate<-apply(temp_theta_sample, 2, mean)
			z_ED=0
			if (all(temp_theta_estimate<=0)){
				z_ED=1
			}else{
				z_ED=min(which(temp_theta_estimate>=0.95*max(temp_theta_estimate)))
			}
			z_ni[n+1]=z_ED
		}
		list(m_ni, s_ni, z_ni)
	}
	m_n_i=append(m_n_i, lapply(res_for, '[[', 1))
	s_n_i=append(s_n_i, lapply(res_for, '[[', 2))
	z_n_i=append(s_n_i, lapply(res_for, '[[', 3))
	
	A=lapply(1:(length(seq_s)-1), function(x) c())
	A_t=lapply(1:(length(seq_m)-1), function(x) A)
	A_n_t=lapply(1:patient, function(x) A_t)
	
	for (i in 1:n_exp){
		if(m_n_i[[i]][k]>seq_m[length(seq_m)]){
			ind_m=length(seq_m)-1
		} else if(m_n_i[[i]][k]<seq_m[1]){
			ind_m=1
		} else{
			ind_m=findInterval(m_n_i[[i]][k], seq_m, rightmost.closed=TRUE)
		}
		ind_s=findInterval(s_n_i[[i]][k], seq_s, rightmost.closed=TRUE)
		A_n_t[[k]][[ind_m]][[ind_s]]=c(A_n_t[[k]][[ind_m]][[ind_s]], i)
	}
		
	res_for<-foreach(i=1:(length(seq_m)-1), .packages="mvtnorm", .export=c("reps", "k", "mu", "sigma", "dose", "patient", "evolution_eq", "J", "M", "n_p", "true_sigma", "seq_m", "seq_s", "A_n_t", "z_star"), .verbose=FALSE)%:%foreach(j=1:(length(seq_s)-1), .packages="mvtnorm", .export=c("reps", "k", "mu", "sigma", "dose", "patient", "evolution_eq", "J", "M", "n_p", "true_sigma", "seq_m", "seq_s", "A_n_t", "z_star"), .verbose=FALSE)%dopar%{
		#create a temporary optimal dose vector
		M_n_t=length(A_n_t[[k]][[i]][[j]])
		if (M_n_t==0){
			z_ni=rep(z_star, patient)
			m_ni=rep((seq_m[i]+seq_m[i+1])/2, patient)
			s_ni=rep((seq_s[j]+seq_s[j+1])/2, patient)
			mu_ni=lapply(1:patient, function(x) mu)
			sigma_ni=lapply(1:patient, function(x) matrix(sigma, nrow=J, ncol=J))
			#for the remaining time(patients) in forward simulation
			for(n in k:(patient-1)){
				z_j=dose[[reps]][n-1]
				temp_theta_sample<-rmvnorm(M, mu_ni[[n-1]], sigma_ni[[n-1]])
				temp_theta_estimate<-apply(temp_theta_sample, 2, mean)
				theta_ED=0
				if (all(temp_theta_estimate<=0)){
					theta_ED=temp_theta_estimate[1]
				}else{
					theta_ED=temp_theta_estimate[min(which(temp_theta_estimate>=0.95*max(temp_theta_estimate)))]
				}
				df_star=rnorm(n_p, theta_ED, sqrt(2*true_sigma^2))
				s_ni[n+1]=1/((1/s_ni[n])+(n_p/(2*true_sigma^2)))
				m_ni[n+1]=s_ni[n]*((m_ni[n]/s_ni[n])+(sum(df_star)/(2*true_sigma^2)))
				y_temp=rnorm(1, temp_theta_sample[z_j], true_sigma)
				res=evolution_eq(mu_ni[[n-1]], sigma_ni[[n-1]], y_temp, z_j)
				mu_ni[[n]]<-res$mu
				sigma_ni[[n]]<-res$sigma
				temp_theta_sample<-rmvnorm(M, mu_ni[[n]], sigma_ni[[n]])
				temp_theta_estimate<-apply(temp_theta_sample, 2, mean)
				z_ED=0
				if (all(temp_theta_estimate<=0)){
					z_ED=1
				}else{
					z_ED=min(which(temp_theta_estimate>=0.95*max(temp_theta_estimate)))
				}
				z_ni[n+1]=z_ED
			}
			#save grid moments, thetas, and optimal dose for all patients in each experiment
			list(m_ni, s_ni, z_ni)
		}
	}
	stopCluster(myCluster)
	m_n_i=compact(append(m_n_i, unlist(lapply(res_for, function(x) lapply(x, '[[', 1)), recursive=FALSE)))
	s_n_i=compact(append(s_n_i, unlist(lapply(res_for, function(x) lapply(x, '[[', 2)), recursive=FALSE)))
	z_n_i=compact(append(z_n_i, unlist(lapply(res_for, function(x) lapply(x, '[[', 3)), recursive=FALSE)))
	
	A=lapply(1:(length(seq_s)-1), function(x) c())
	A_t=lapply(1:(length(seq_m)-1), function(x) A)
	A_n_t=lapply(1:patient, function(x) A_t)
	
	for (i in 1:(length(m_n_i))){
		for(n in k:patient){
			if(m_n_i[[i]][n]>seq_m[length(seq_m)]){
				ind_m=length(seq_m)-1
			} else if(m_n_i[[i]][n]<seq_m[1]){
				ind_m=1
			} else{
				ind_m=findInterval(m_n_i[[i]][n], seq_m, rightmost.closed=TRUE)
			}
			ind_s=findInterval(s_n_i[[i]][n], seq_s, rightmost.closed=TRUE)
			A_n_t[[n]][[ind_m]][[ind_s]]=c(A_n_t[[n]][[ind_m]][[ind_s]], i)
		}
	}

	#Initial values for decisions costs 
	#optimal decision and utility list
	d_n_t=lapply(1:patient, function(x) matrix(NA, nrow=n_grid, ncol=n_grid))
	u_n_t=lapply(1:patient, function(x) matrix(NA, nrow=n_grid, ncol=n_grid))
	ind_n=lapply(1:patient, function(x) c())
	ind_t=lapply(1:((length(seq_m)-1)*100+(length(seq_s)-1)), function(x) ind_n)
	#sweep the grid from the last patient to the current patient and evaluate each cells utilities
	for (n in patient:k){
		for (i in 1:(length(seq_m)-1)){
			for (j in 1:(length(seq_s)-1)){
				U2=0
				U1=0
				M_n_t=length(A_n_t[[n]][[i]][[j]])
				if(M_n_t>0){			
					for (t in A_n_t[[n]][[i]][[j]]){
						m=m_n_i[[t]][n]
						s=s_n_i[[t]][n]
						U2=U2+c2*m*(1-pnorm(qnorm(0.90)-(m*sqrt(n_p))/(sqrt(2*true_sigma^2+s^2))))
						ind_t[[t]][[n]]=c(i, j)
					}
					U2=-c1*n_p+U2/M_n_t
					if (n==patient){
						u_n_t[[n]][i,j]=max(0, U2)
						if (U2 >0){
							d_n_t[[n]][i,j]=2 #terminate
						}else{
							d_n_t[[n]][i,j]=0 #abandon
						}
					} else{
						for (t in A_n_t[[n]][[i]][[j]]){
							i1=ind_t[[t]][[n+1]][1]
							j1=ind_t[[t]][[n+1]][2]
							U1=U1+u_n_t[[n+1]][i1, j1]
						}
						U1=-c1+U1/M_n_t
						u_n_t[[n]][i, j]=max(0, U1, U2)
						d_n_t[[n]][i, j]=which(c(0, U1, U2)==max(0, U1, U2))-1
					}	
				}
			}
		}
	}
	d_n_t[[k]][is.na(d_n_t[[k]])]=-1
	return(list("dmat"=d_n_t[[k]], "umat"=u_n_t[[k]]))
}

start_time=Sys.time()
mat_d=lapply(1:patient, function(x) matrix(0, nrow=n_grid, ncol=n_grid))
alloc_m=lapply(1:patient, function(x) c())
alloc_s=lapply(1:patient, function(x) c())
save_decision=lapply(1:n_simulation, function(x) c())
for(reps in 1:n_simulation){
	set.seed(reps)
	m_n=rep(m_0, patient)
	s_n=rep(s_0, patient)
	for (k in 1:patient){
		if (k==1){#at first decision epoch
			sigma_n[[k]]=sigma_0
			mu_n[[k]]=mu_0
			theta_sample=rmvnorm(M, mu_n[[k]], sigma_n[[k]])
			ED95=c()
			ED95=apply(theta_sample, 1, function(z){
				if (all(z<0)){
					return(NA)
				}else{
					return (min(which(z>=0.95*max(z))))
				}
			})
			theta_ED=unlist(sapply(1:length(ED95), function(nt)	{
				if (is.na(ED95[nt])==FALSE){
					return(theta_sample[nt, ED95[nt]])
				}
			}))
			s_n[k]=var(theta_ED)
			m_n[k]=mean(theta_ED)
			z_j=dose[[reps]][k]
			y=c(y, rnorm(1, true_theta[z_j], true_sigma))		
		} else{#assignment dose is allocated from a file containing optimal allocated doses
			theta_sample=rmvnorm(M, mu_n[[k-1]], sigma_n[[k-1]])
			ED95=c()
			ED95=apply(theta_sample, 1, function(z){
				if (all(z<0)){
					return(NA)
				}else{
					return (min(which(z>=0.95*max(z))))
				}
			})
			theta_ED=unlist(sapply(1:length(ED95), function(nt)	{
				if (is.na(ED95[nt])==FALSE){
					return(theta_sample[nt, ED95[nt]])
				}
			}))
			s_n[k]=var(theta_ED)
			m_n[k]=mean(theta_ED)
			#if(k%%5==0){
			if(k>=initial_patient){
			#if(k>=50 & k<=200){#in optimal decision epochs
				##calculating m,s to see whether stop, continue, or abandon 
				alloc_m[[k]]=c(alloc_m[[k]], m_n[k])
				alloc_s[[k]]=c(alloc_s[[k]], s_n[k])
				##call function optimal_stopping for patient k with already observed responses
				res_opt=optimal_stopping(k, y, m_n[k], s_n[k], mu_n[[k-1]], sigma_n[[k-1]], z_ED)
				##decision grid
				mat_d[[k]]=res_opt$dmat
				##optimal objective grid
				##sequence of grid points along m and s axis
				##find what is the decision based on m,s and the grid
				if (m_n[k]>seq_m[length(seq_m)]){
					m_ind=length(seq_m)-1
				}else if (m_n[k]>seq_m[1]){
					m_ind=1
				} else{
					m_ind=findInterval(m_n[k], seq_m, rightmost.closed=TRUE)
				}
				s_ind=findInterval(s_n[k], seq_s, rightmost.closed=TRUE)
				decision=c(decision, res_opt$dmat[m_ind, s_ind])
				if (res_opt$dmat[m_ind, s_ind]==0 || res_opt$dmat[m_ind, s_ind]==2){
					break
				}
			}
			##observe the true response of the optimal dose and add it to original data
			z_j=dose[[reps]][k] #assignment dose
			y=c(y, rnorm(1, true_theta[z_j], true_sigma)) #observing and add the true response of the optimal dose
			res=evolution_eq(mu_n[[k-1]], sigma_n[[k-1]], y[length(y)], z_j)
			mu_n[[k]]<-res$mu
			sigma_n[[k]]<-res$sigma
		}
	}
	#############save results###############################
	save_decision[[reps]]=decision
	#############reset vectors for iterations###############
	y=c() #a list (for each dose) of vectors of observed responses
	decision=c() #vector of optimal stopping for patient k
	mu_n=lapply(1:patient, function(x) c())
	sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
}
write.table(save_decision, file=paste("C:/Results/correct-selection/gridding-",curve_st,".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
end_time=Sys.time()
start_time-end_time

end_time=Sys.time()
start_time-end_time
