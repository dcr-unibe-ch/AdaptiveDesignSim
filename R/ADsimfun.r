

#' ADsimfun
#'
#' Simulates trials with a binary endpoint and one or several interim analysis,
#'
#' @param n01 Total sample size in each arm, numeric vector of length 2
#' @param p01 Proportion in each arm, numeric vector of length 2
#' @param nia Number of interim analysis, not including the final analysis
#' @param tia Optional vector with the time point(s) of the interim analyses, as fraction of the total (between 0 and 1)
#' @param alpha Alpha level for final analysis
#' @param direct Direction of effect, "lower" (lower is better) or higher ("higher is better")
#' @param effm Effect measure, "rd" (risk differnce) or "rr" (risk ratio)
#' @param cimethod Confidence interval used for analysis, "Wald" or "score".
#'
#' @returns Numeric vector with true proportions in both arms and point estimate (true_p0, true_p1, true_pe),
#'  and for each interim stage (prefix iax_) and the final stage (prefix fa_),
#'  the simulated number of successes, sample size and observed probabilities for both arms (suffix x0, x1, n0, n1, p0, p1),
#'  point estimate (suffic pe), confidence interval (suffic lci and uci) and z-statistic (z).
#'
#' @export
#'
#' @importFrom stats rbinom qnorm
#'
#' @examples
#' set.seed(1)
#' ADsimfun(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, alpha = 0.025)
#' set.seed(1)
#' ADsimfun(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, tia = 0.4, alpha = 0.025)
#' set.seed(1)
#' ADsimfun(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, alpha=0.025, effm="rr")
#' set.seed(1)
#' ADsimfun(n01 = c(100,100), p01 = c(0.4,0.2), nia = 2, alpha=0.025,
#'   direct="higher", effm="rd")


ADsimfun<-function(n01, p01, nia, tia = NULL, alpha = 0.025, direct = c("lower","higher"), effm = c("rd","rr"),
	cimethod = c("Wald","score")) {


	if (!is.null(tia)) {
		if (length(tia)!=nia) {
			warning("Length of tia does not match nia. Equal spacing is assumed.")
			tia<-NULL
		} 
		else {
			if (min(tia)<=0 | max(tia)>=1) {
				warning("tia should be between 0 and 1. Equal spacing is assumed.")
				tia<-NULL
			}
		}	
	}

	if (is.null(tia)) {
	
		#with ratio 
		sp<-1/(nia+1)
		sp<-seq(sp,1,by=sp)
		sp<-sp[1:(length(sp)-1)]
		pia<-sapply(n01,function(x) round(x*sp))
		
		#with patients
		#sp<-n01/(nia+1)
		#pia<-sapply(1:2,function(x) round(seq(sp[x],l=nia,by=sp[x])))
		#if (nia==1) {
		#	pia<-t(matrix(pia))
		#}
	} else {
		pia<-sapply(n01,function(x) round(x*tia))	
	}
	
	if (nia==1) {
		pia<-t(matrix(pia))
	}
		
	if (length(pia)==0) {
		pia<-NULL
		stopifnot(nia==0)
	}
	
	cimethod<-match.arg(cimethod)
	effm<-match.arg(effm)
	direct<-match.arg(direct)

	clevel<-(1-alpha*2)

	outs0<-rbinom(n=n01[1], size=1, prob=p01[1])
	outs1<-rbinom(n=n01[2], size=1, prob=p01[2])

	nam<-c("x0","x1","n0","n1","p0","p1","pe","lci","uci","z")

	#final
	xs<-c(sum(outs0),sum(outs1))
	ns<-c(length(outs0),length(outs1))
	ps<-xs/ns

	if (cimethod=="Wald") {
		est<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
	}
	if (cimethod=="score") {
		if (effm=="rd") {
			pd<-pairwiseCI::Prop.diff(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
				CImethod="NHS",conf.level = clevel)
			west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
			est<-c(pd$estimate,pd$conf.int,west["z"])
		} else {
			pd<-pairwiseCI::Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
				CImethod="Score",conf.level = clevel)
			west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
			est<-c(log(pd$estimate),log(pd$conf.int),west["z"])
		}

	}

	if (effm=="rd") {
		true_pe<-p01[2]-p01[1]
	} else {
		true_pe<-log(p01[2]/p01[1])
	}

	rif<-c(p01[1],p01[2],true_pe,xs,ns,ps,est)
	names(rif)<-c("true_p0","true_p1","true_pe",paste0("fa_",nam))

	#interim
	if (nia>0) {
		for (ia in 1:nia) {

			ni0<-pia[ia,1]
			ni1<-pia[ia,2]
			
			outi0<-outs0[1:ni0]
			outi1<-outs1[1:ni1]

			xs<-c(sum(outi0),sum(outi1))
			ns<-c(length(outi0),length(outi1))
			ps<-xs/ns

			if (cimethod=="Wald") {
				est<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
			}
			if (cimethod=="score") {
				if (effm=="rd") {
					pd<-pairwiseCI::Prop.diff(x = c(xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
						CImethod="NHS",conf.level = clevel)
					west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
					est<-c(pd$estimate,pd$conf.int,west["z"])
				} else {
					pd<-pairwiseCI::Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
						CImethod="Score",conf.level = clevel)
					west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
					est<-c(log(pd$estimate),log(pd$conf.int),west["z"])
				}
			}

			rii<-c(xs,ns,ps,est)

			names(rii)<-paste0("ia",ia,"_",nam)

			rif<-c(rif,rii)
		}
	}
	return(rif)
}

#Calculates Wald CI

waldest<-function(xs, ns, alpha,  effm = c("rd","rr"), direct = c("lower","higher")) {

	effm<-match.arg(effm)
	direct<-match.arg(direct)

	phat<-xs/ns

	if (effm=="rd") {

		rd<-phat[2]- phat[1]

		#wald ci, individual props
		se <- sqrt((phat[1] * (1 - phat[1])) / ns[1] +
				   (phat[2] * (1 - phat[2])) / ns[2])

		pp<-sum(xs)/sum(ns)
		sepool<-sqrt(pp*(1-pp)*(1/ns[1] + 1/ns[2]))

		lci<-rd - qnorm(1-alpha)*se
		uci<-rd + qnorm(1-alpha)*se

		if (direct=="lower") {
			#z<-(-rd)/sepool
			z<-(-rd)/se
		} else {
			z<-(rd)/se
		}

		est<-c(rd,lci,uci,z)

	} else {

		lrr <- log(phat[2]/phat[1])

		se <- sqrt(1/xs[1] - 1/ns[1] + 1/xs[2] - 1/ns[2])

		lci <- lrr - qnorm(1-alpha)*se
		uci <- lrr + qnorm(1-alpha)*se

		if (direct=="lower") {
			z<-(-lrr)/se
		} else {
			z<-(lrr)/se
		}

		est<-c(lrr,lci,uci,z)
	}


	names(est)<-c("pe","lci","uci","z")
	return(est)

}




