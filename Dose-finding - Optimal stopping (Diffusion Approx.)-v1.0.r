library(MASS)
###############Parameter initialization###################
patient=400#number of patients
initial_patient=280
true_sigma=10#observation variance
n_p=30 #number of patients in confirmatory phase-III
c1=1000 #cost of sampling one more patient
c2=1000000 #reward of one unit improvement over placebo
start_time=0 #start time of the algorithm
end_time=0 #end time of the algorithm
###############Diffision Approximation. Trinomial Tree#############################
m_0=0
nu_0=100
#t_0=(2*true_sigma^2)/(nu_0)
t_0=(true_sigma^2)/(nu_0)
max_m=20
min_m=0
delta_t=0.1
p_max=0.475 #(Arlotto)
#delta_m=(sqrt(2)*true_sigma*sqrt(delta_t))/(sqrt(2*(t_0+initial_patient)*(t_0+initial_patient+1)*p_max)) #(Arlotto)
delta_m=(true_sigma*sqrt(delta_t))/(sqrt(2*(t_0+initial_patient)*(t_0+initial_patient+1)*p_max)) #(Arlotto)
m_interv=seq(min_m, max_m, delta_m)
M=length(m_interv)
t_interv=seq(initial_patient, patient, delta_t)
T=length(t_interv)
B=matrix(NA, nrow=T, ncol=M)
#Final time
start_time=Sys.time()
##Solution of trinomial tree approximation to PDE in the last time interval
for (j in 1:M){
	#first maximand in G(y_T, T)
	B1=0
	#second maximand in G(y_T, T)
	#B2=-c1*n_p+c2*(1-pnorm(qnorm(0.99)-(m_interv[j])/sqrt((2*true_sigma^2+2*true_sigma^2/t_interv[T])/n_p),0, 1, lower.tail = TRUE))*(m_interv[j])
	B2=-c1*n_p+c2*(1-pnorm(qnorm(0.99)-(m_interv[j])/sqrt((2*true_sigma^2+true_sigma^2/t_interv[T])/n_p),0, 1, lower.tail = TRUE))*(m_interv[j])
	#enumerating the B matrix
	B[T, j]=max(B1, B2)
}
#Backward solution of trinomial tree approximation to PDE
for(i in (T-1):1){
	for(j in 1:M){
		#P_d=0.5*(2*true_sigma^2*delta_t/t_interv[i]/delta_m^2)
		#P_d=(2*true_sigma^2*delta_t)/(2*(as.integer(t_interv[i])+t_0)*(t_0+as.integer(t_interv[i])+1)*delta_m^2) #(Arlotto)
		P_d=(true_sigma^2*delta_t)/(2*(as.integer(t_interv[i])+t_0)*(t_0+as.integer(t_interv[i])+1)*delta_m^2) #(Arlotto)
		P_u=P_d
		P_m=1-P_d-P_u
		if(j==1){
			B[i, j]=P_u*B[i+1, j+1]+P_m*B[i+1, j]+P_d*(2*B[i+1, j]-B[i+1, j+1])
		}else if(j==M){
			B[i, j]=P_u*(2*B[i+1, j]-B[i+1, j-1])+P_m*B[i+1, j]+P_d*B[i+1, j-1]
		}else{
			B[i, j]=P_u*B[i+1, j+1]+P_m*B[i+1, j]+P_d*B[i+1, j-1]
		}
	}
}
for(i in 1:T){
	B[i,]=B[i,]-c1*delta_t
}

G=matrix(NA, nrow=T, ncol=M)
difference=lapply(1:T, function(x) c())
for(i in 1:T){
	for(j in 1:M){
		#G[i, j]=max(0, -c1*n_p+c2*(1-pnorm(qnorm(0.99)-(m_interv[j])/sqrt((2*true_sigma^2+2*true_sigma^2/t_interv[i])/n_p),0, 1, lower.tail = TRUE))*(m_interv[j]))
		G[i, j]=max(0, -c1*n_p+c2*(1-pnorm(qnorm(0.99)-(m_interv[j])/sqrt((2*true_sigma^2+true_sigma^2/t_interv[i])/n_p),0, 1, lower.tail = TRUE))*(m_interv[j]))
		if (B[i,j]-G[i,j]>=0){
			difference[[i]]=c(difference[[i]], j)
		}
		# if (abs(B[i,j]-G[i,j])>=0.01){
			# difference[[i]]=c(difference[[i]], j)
		# }
	}
}




bound_u=c(unlist(lapply(difference, function(x) x[length(x)])))
bound_d=c(unlist(lapply(difference, function(x) x[1])))
smooth_b_u=loess.smooth(t_interv[1:length(bound_u)], bound_u*delta_m, span=1, degree=1, family="gaussian", evaluation=T)
smooth_b_d=loess.smooth(t_interv[1:length(bound_u)], bound_d*delta_m, span=1, degree=1, family="gaussian", evaluation=T)

write(smooth_b_u$y, file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/Upper_bound_smooth.txt", sep=""), append=FALSE, sep="\n")
write(smooth_b_d$y, file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/Lower_bound_smooth.txt", sep=""), append=FALSE, sep="\n")
write(smooth_b_d$x, file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/x_axis.txt", sep=""), append=FALSE, sep="\n")

write(bound_u, file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/Upper_bound.txt", sep=""), append=FALSE, sep="\n")
write(bound_d, file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/Lower_bound.txt", sep=""), append=FALSE, sep="\n")

end_time=Sys.time()
start_time-end_time
##############################################################################
#############################diffusion path###################################
##############################################################################
curve_st="sigmoid-significant"
#sigma_string="std_dev=sqrt(1)"
#sigma_string="std_dev=sqrt(10)"
sigma_string="std_dev=10"
n_simulation=30
diffusion_path=lapply(1:n_simulation, function(x) c())
objective_path=lapply(1:n_simulation, function(x) c())
J=11 #number of doses
theta_estimate=lapply(1:n_simulation, function(x) matrix(NA, nrow=patient, ncol=J)) #a matrix of optimal theta estimates
for (reps in 3:n_simulation){
	theta_estimate[[reps]]=read.table(file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/dose-responses/",sigma_string,"/",curve_st,"/",toString(patient),"patients-KGC-thetas-",toString(reps),".txt", sep=""), sep="\t", nrow=patient)
	for (k in 1:patient){
			if (all(theta_estimate[[reps]][k,]<0)){
				diffusion_path[[reps]]=c(diffusion_path[[reps]], theta_estimate[[reps]][k, 1])
			} else{
				diffusion_path[[reps]]=c(diffusion_path[[reps]], theta_estimate[[reps]][k, min(which(theta_estimate[[reps]][k,]>=0.95*max(theta_estimate[[reps]][k,])))])
			}
			if (k>=initial_patient){
				m_path=findInterval(diffusion_path[[reps]][k], m_interv, rightmost.closed=TRUE, all.inside=TRUE)
				t_path=findInterval(k, t_interv, rightmost.closed=TRUE, all.inside=TRUE)
				objective_path[[reps]]=c(objective_path[[reps]], B[t_path, m_path])			
			}
	}
}
write.table(as.data.frame(diffusion_path), file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/correct-selection/",toString(patient),"patients-diff-path-",curve_st,".csv", sep=""), append=FALSE, row.names=FALSE, col.names=FALSE, sep=",", eol="\n")
write.table(as.data.frame(objective_path), file=paste("C:/Users/snasrol/Google Drive/Research-Optimal stopping in dose-finding clinical trials/Codes/Results/09-24-2018/correct-selection/",toString(patient),"patients-obj-path-",curve_st,".csv", sep=""), append=FALSE, row.names=FALSE, col.names=FALSE, sep=",", eol="\n")
