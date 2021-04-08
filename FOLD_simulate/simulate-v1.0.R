
#####################################################################
# simulate-v1.0.R
# This code is for meta-analysis association studies with overlapping samples
# This code simulates splitting, Lin and Sullivan's(LS) approaches
# We only consider overlapping control samples
# Eunji Kim (2016,June)
#####################################################################

simulate <- function(outputfile="Output",nrep=20000,nstudy=5,RRs=rep(1.15,5),ncase=rep(5000,5),ncont=rep(0,5),nshared=5000,t=5E-8) {

	# outputfile: user defined output file names, save power of Lin and Sullivan and splitting
	# nrep: number of simulation repeats
	# nstudy:number of studies involved in meta-analysis
	# RRs,ncase,ncont and nshared are vectors whose sizes are nstudy.
	# RRs: relative risk of disease 
	# ncase: number of case samples per study, i.e) c(1000,2000) represents 1000 for study1 and 2000 for study2
	# ncont: number of study-specific controls per study
	# nshared: number of shared controls per study
	# t: pvalue threshold
	# preval: disease prevalence rate 
	# maf: minor allele frequency(MAF), default value is 0.3

	maf=0.3
	preval=0.001

	# This code uses equal splitting strategy, which equally distributing shared controls 
	nsplit=nshared/nstudy

	# Case/control frequencies
	casemaf=RRs*maf / ((RRs-1)*maf + 1)
	contmaf=(maf-preval*casemaf)/(1-preval)

	# Build Lin and Sullivan's correlation matrix
	R=diag(nstudy)
	for (i in 1:nstudy) {
		for (j in 1:nstudy) {
			if (i == j) {
				next
			} else {
				total=sqrt((ncase[i]+ncont[i]+nshared)*(ncase[j]+ncont[j]+nshared))
				controlpart=nshared*sqrt(ncase[i]*ncase[j]/(ncont[i]+nshared)/(ncont[j]+nshared))
				rij=controlpart/total
				R[i,j]=R[j,i]=rij
			}
		}
	}
	invR=solve(R)

	# This records pvalue of LS and splitting
	splitpvalue=rep(0, nrep)
	LSpvalue=rep(0, nrep)

        # Generate samples 
	for (eachrep in 1:nrep) {
		betas.nonde=rep(0,nstudy)
		stder.nonde=rep(0,nstudy) 
		betas.split=rep(0,nstudy)
		stder.split=rep(0,nstudy) 
		casemac=rep(0,nstudy)
		contmac=rep(0,nstudy)
		splitmac=rep(0,nstudy)

		for (i in 1:nstudy) {
			casemac[i]=rbinom(1, ncase[i]*2, casemaf[i])
			contmac[i]=rbinom(1, ncont[i]*2, contmaf[i])
			splitmac[i]=rbinom(1, 2*nsplit, contmaf[i])
		}

		sharedmac=rep(sum(splitmac),nstudy)

		for (i in 1:nstudy) {
			## non-splitting
			n11=casemac[i]
			n10=2*ncase[i]-casemac[i]
			n01=contmac[i]+sharedmac[i]
			n00=2*ncont[i]+2*nshared-contmac[i]-sharedmac[i]
			betas.nonde[i]=log(n11*n00/n01/n10)
			stder.nonde[i]=sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)

			## splitting
			n11=casemac[i]
			n10=2*ncase[i]-casemac[i]
			n01=contmac[i]+splitmac[i]
			n00=2*ncont[i]+2*nsplit-contmac[i]-splitmac[i]
			betas.split[i]=log(n11*n00/n01/n10)
			stder.split[i]=sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)
		}

		## Calculate pvalue using Lin and Sullivan's(LS) approach
		unit=t(rep(1,nstudy))
		LSinvsigma=solve(diag(stder.nonde)) %*% invR %*% solve(diag(stder.nonde))
		LSinvvar=(unit %*% LSinvsigma %*% t(unit))
		Linbeta=(unit %*% LSinvsigma %*% t(t(betas.nonde)))/LSinvvar
		LSpvalue[eachrep]=2*pnorm(abs(Linbeta/sqrt(1/LSinvvar)),lower.tail=FALSE)

		## Calculate pvalue using splitting approach
		split_invsigma=solve(diag(stder.split)) %*% diag(rep(1,nstudy)) %*% solve(diag(stder.split))
		split_invvar=(unit %*% split_invsigma %*% t(unit))
		splitbeta=(unit %*% split_invsigma %*% t(t(betas.split)))/split_invvar
		splitpvalue[eachrep]=2*pnorm(abs(splitbeta/sqrt(1/split_invvar)),lower.tail=FALSE)
	}
	## Write power to outputfile
	cat("Power of Lin and Sullivan",sum(LSpvalue<t)/nrep,"Power of Splitting",sum(splitpvalue<t)/nrep,"\n",sep=",")
	output.power=paste("Study-specific controls",ncont[1],"shared controls",nshared,"Power of Lin and Sullivan",sum(LSpvalue<t)/nrep,"Power of Splitting",sum(splitpvalue<t)/nrep,sep=",")
	write.table(output.power,file=paste(outputfile,".power",sep=""), quote=F, col.names=F, row.names=F)
}
