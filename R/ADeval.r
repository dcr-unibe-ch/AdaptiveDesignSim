#' ADeval
#'
#' Function to evaluate operating characteristics of trials simulated via ADsim
#'
#' @param simdat matrix or data.frame of multiple runs of ADsim
#' @param zcutoff a matrix (or data frame) with the cut-offs on the z-statstics at interim.
#'  One row per interim analysis and a column for stop for futilty and efficacy (with NA if not planned)
#' @param finalm Measure used for the final analysis, either z (z-statistics) or ci (confidence interval).
#' @param alpha Alpha level used for the final analysis if finalm=="z", otherwise the ci given in simdat is used.
#'
#' @returns A list with the evaluation for each repetition (**rawdata**) and the summarized operating characteristics (**opchar**).
#' **opchar** includes the average estimates in each group, effect (point estimate, **pe**) and (1-alpha) confidence limits (**lci**, **uci**), and
#' - **psig**: the probability for a significant trial (the power or type I error, depending on the assumptions).
#' - **pstop**: the probability for any stop and a stop for futility or efficacy (**pstop_fut**, **pstop_eff**).
#' - the average, minimum and maximum sample size (**avn**, **minn**, **maxn**).
#' - **bias**: the average of estimate minus true value.
#' - **sd**: the standard deviation of the point estimate.
#' - **mse**: the mean squared error, the average of the squared difference between estimate and true value.
#' - **coverage**: the proportion of confidence intervals including the true value.
#'
#' @export
#'
#' @examples
#' #Simulate 1000 trials with one interim analysis
#' simdat <- lapply(1:1000, function(i)
#' 		ADsim(n01 = c(100,100), par01 = c(0.2, 0.4), nia = 1))
#' simdat <- do.call(rbind, simdat)
#'
#' # Define stop for futility stoppin (z<0) and efficacy stopping (z>2)
#' zcutoff<-matrix(c(0,2),1,2)
#' sr<-ADeval(simdat = simdat, zcutoff = zcutoff)
#' head(sr$rawdata)
#' sr$opchar
#'
#'
ADeval<-function(simdat, zcutoff, finalm = c("z","ci"), alpha = 0.025) {

	finalm<-match.arg(finalm)

	#derive n and true values from simdat
	nia<-sum(grepl("ia[0-9]_pe",colnames(simdat)))
	n<-mean(simdat[,"fa_n0"]) + mean(simdat[,"fa_n1"])
	#tp0<-mean(simdat[,"true_p0"])
	#tp1<-mean(simdat[,"true_p1"])
	rdt<-mean(simdat[,"true_pe"])

	#derive direction from pe and z
	spe<-sign(simdat[,"fa_pe"])
	sz<-sign(simdat[,"fa_z"])
	if (all(spe[spe!=0]==sz[spe!=0])) {
		direct<-"higher"
	} else {
		stopifnot((all(spe[spe!=0]!=sz[spe!=0])))
		direct<-"lower"
	}

	resi<-data.frame(simdat)

	fans<-colnames(simdat)[grepl("^fa_",colnames(simdat))]
	fan<-fans[1]
	for (fan in fans) {
		nname<-sub("fa_","obs_",fan)
		resi<-cbind(resi,nname=resi[,fan])
		colnames(resi)[colnames(resi)=="nname"]<-nname
	}


	if (nia>0) {
		sel<-data.frame(matrix(NA,nrow(resi),nia))
		colnames(sel)<-1:nia

		for (i in (1:nia)) {

			zs<-zcutoff[i,]
			seli<-rep(0,nrow(resi))
			if (!is.na(zcutoff[i,1])) {
				seli[resi[,paste0("ia",i,"_z")]<zcutoff[i,1]]<-1
			}
			if (!is.na(zcutoff[i,2])) {
				seli[resi[,paste0("ia",i,"_z")]>zcutoff[i,2]]<-2
			}
			sel[,i]<-seli
		}

		#use rd at interim at the first stopped rd
		resi$tstop<-apply(sel,1,function(x) ifelse(any(x>0),min(which(x>0)),NA))
		#resi$neff<-n

		i<-1
		for (i in (1:nia)) {
			si<-!is.na(resi$tstop) & resi$tstop==i
			
		  	for (fan in fans) {
				nname<-sub("fa_","obs_",fan)
				iname<-sub("fa_",paste0("ia",i,"_"),fan)
				resi[si,nname]<-resi[si,iname]
			}		
		}

		resi$tstop_fut<-apply(sel,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))
		resi$tstop_eff<-apply(sel,1,function(x) ifelse(any(x==2),min(which(x==2)),NA))
		resi$anystop_fut<-!is.na(resi$tstop_fut) & (is.na(resi$tstop_eff) | (resi$tstop_fut<=resi$tstop_eff))
		resi$anystop_eff<-!is.na(resi$tstop_eff) & (is.na(resi$tstop_fut) | (resi$tstop_eff<=resi$tstop_fut))

		resi$anystop<-"Not stopped"
		resi$anystop[resi$anystop_fut]<-"Stopped for futility"
		resi$anystop[resi$anystop_eff]<-"Stopped for efficacy"

	} else {
		resi$anystop<-"Not stopped"

	}


	resi$anystop<-factor(resi$anystop,
		levels=c("Not stopped","Stopped for futility","Stopped for efficacy"))
	
	
	#average summarize: 
	avs<-numeric(0)
	for (fan in fans) {
		oname<-sub("fa_","obs_",fan)
		avs<-append(avs,mean(resi[,oname]))	
	}
	names(avs)<-gsub("fa_","av_",fans)
	
	#power 
	if (finalm == "ci") {
		if (direct=="lower") {
			pwr<-mean(resi$obs_uci<0)
		} else {
			pwr<-mean(resi$obs_lci>0)
		}
	} else {
		pwr<-mean(resi$obs_z > qnorm(1-alpha))
	}

	#stop at any ia
	pstop<-mean(resi$anystop %in% c("Stopped for futility","Stopped for efficacy"))
	pstopfut<-mean(resi$anystop %in% c("Stopped for futility"))
	pstopeff<-mean(resi$anystop %in% c("Stopped for efficacy"))

	#stop at each ia:
	pstopi<-pstopfuti<-pstopeffi<-numeric(0)
	for (i in (1:nia)) {
		pstopfuti<-c(pstopfuti,mean(!is.na(resi$tstop_fut) & resi$tstop_fut==i & resi$anystop_fut))
		pstopeffi<-c(pstopeffi,mean(!is.na(resi$tstop_eff) & resi$tstop_eff==i & resi$anystop_eff))
		pstopi<-c(pstopi,mean(!is.na(resi$tstop) & resi$tstop==i &
			resi$anystop %in% c("Stopped for futility","Stopped for efficacy")))
	}
	nastopi<-c(paste0("pstop_ia",1:nia),paste0("pstop_fut_ia",1:nia),paste0("pstop_eff_ia",1:nia))
	
	#expected total ss 
	av_n<-avs["av_n0"] + avs["av_n1"] 
	min_n<-min(resi$obs_n0 + resi$obs_n1)
	max_n<-max(resi$obs_n0 + resi$obs_n1)
	bias<-mean(resi$obs_pe - rdt)
	sd<-sd(resi$obs_pe)
	mse<-mean((resi$obs_pe - rdt)^2)
	coverage<-mean(resi$obs_lci<rdt & resi$obs_uci>rdt)

	op<-c(rdt,avs,
		pwr,pstop,pstopfut,pstopeff,
			pstopi,pstopfuti,pstopeffi,
			av_n,min_n,max_n,
			bias,sd,mse,coverage)

	names(op)<-c("true_pe",names(avs),
		"psig","pstop","pstop_fut","pstop_eff",nastopi,
		"av_n","min_n","max_n",
		"bias","sd","mse","coverage")

	res<-list(resi,op)
	names(res)<-c("rawdata","opchar")
	return(res)
}

