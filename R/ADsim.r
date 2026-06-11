#' ADsim
#'
#' Simulates trials with a binary endpoint and one or several interim analysis,
#'
#' @param n01 Total sample size in control and experimental arm, numeric vector of length 2.
#' @param nia Number of interim analysis, not including the final analysis.
#' @param tia Optional vector with the time point(s) of the interim analyses, as fraction of the total (between 0 and 1).
#' @param par01 Parameters in control and experimental arm (which is used as argument for simfun).
#' @param simfun Simulation function or one of the default options ("binom", "norm").
#' @param estfun Function to calculate point estimate, CIs and z- statistic or one of the
#'	five defaults functions ("rd_wald", "lrr_wald", "rd_score", "lrr_score", "md_t").
#' @param cilevel Level of the two-sided confience interval for the default estfun.
#' @param truefun Function to calculate the true point estimate from par01 (ignored if the default estfun are used).
#' @param direct Direction of the effect, "higher" (default) or "lower" values are better.
#'
#' @details Parameters for simulation are defined with *par01*, the simulation function with *simfun*.
#'	For the default simfun "binom" *par01* must be a numeric vector of length 2 with the probabilities in control 
#'	and experimental group. For the default simfun "norm", *par01* must be a 2x2 data frame 
#'	with rows *mean* and *sd* and first column for control and second for experimental.
#' User-specified *simfun*,  must take *par01* as argument.
#'	
#' For estimation, five defaults can be specified via *estfun* - standard Wald-type approximation 
#'	for the risk difference (*wald_rd*) and risk ratio (*wald_rr*),
#'	the score-based method for the risk difference suggested by Newcombe (*score_rd*), 
#'	the score-based method for the risk ratio according to Koopman (*score_rr*), 
#'	and Student's t-test for the mean difference (*md_t*).
#' 	The effect is allways given as experimental (second element of par01) vs control (first element of par01). 
#'
#' User-specified *estfun* must take the output from *simfun* and *cilevel* as argument,
#'	and produce a named numeric vector with point estimate (*pe*), confidence limits (*lci*, *uci*) and z-statistic (*z*).
#'	If a user-specifid *estfun* is used, a function to calculate the true effect must be given via *truefun*.
#'
#' The direction if the effect in *direct* only affects the sign of the z-statistic.
#'	So that positive values of the z-statistic are always assocaiated with a benefit of the experimental arm.
#'
#' @returns Numeric vector with true point estimate (*true_pe*),
#'  and for each interim stage (prefix *iax_*) and the final stage (prefix *fa_*),
#'	the sample size (suffix *n0* and *n1*), point estimate (suffix *pe*), 
#'	confidence interval (suffix *lci* and *uci*) and z-statistic (*z*)
#'  And depending on the setting further group-specific paremeters
#' (e.g. the simulated number of successes and observed probabilities for both arms).
#'
#' @export
#'
#' @importFrom stats rbinom qnorm uniroot rnorm sd t.test
#'
#' @examples
#' set.seed(1)
#' ADsim(n01 = c(100,100), par01 = c(0.2,0.4), nia = 1)
#'
#' #Earlier interim analysis
#' set.seed(1)
#' ADsim(n01 = c(100,100), par01 = c(0.2,0.4), nia = 1, tia = 0.4)
#'
#' #Difference effect estimation
#' set.seed(1)
#' ADsim(n01 = c(100,100), par01 = c(0.2,0.4), nia = 1, estfun = "lrr_wald")
#'
#' #Assume higher would be better
#' set.seed(1)
#' ADsim(n01 = c(100,100), par01 = c(0.2,0.4), nia = 2, direct="lower")
#'
#' #User defined function that reproduces the default setting
#' efun<-function(out,cilevel) {
#' 	xs<-unlist(lapply(out,sum))
#' 	ns<-unlist(lapply(out,length))
#' 	pe<-diff(xs/ns)
#' 	se<-sqrt(sum(xs/ns*(1-xs/ns)/ns))
#' 	lci<-pe - qnorm(1-(1-cilevel)/2)*se
#' 	uci<-pe + qnorm(1-(1-cilevel)/2)*se
#' 	z<-(-pe)/se
#' 	res<-c(pe,lci,uci,z)
#' 	names(res)<-c("pe","lci","uci","z")
#' 	return(res)
#' }
#' 
#' tfun <- function(x) diff(x)
#' 
#' set.seed(1)
#' ADsim(n01 = c(100,100), par01 = c(0.4,0.2), nia = 1, estfun = efun, truefun = tfun)
#'
#' #Continuous outcome 
#' par01<-data.frame(matrix(c(5,1,5.5,1),2,2))
#' rownames(par01)<-c("mean","sd")
#' ADsim(n01 = c(100,100), par01 = par01, nia = 1, estfun = "md_t", simfun="norm")
ADsim<-function(n01,
	nia,
	tia = NULL,
	par01,
	simfun = "binom",
	estfun = "rd_wald",
	cilevel = 0.95,
	truefun = NULL,
	direct = c("higher","lower")) {

	direct<-match.arg(direct)

	#time point of IA
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

	#number of patients at IA
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

	#function for true point estimate and estimation
	if (is.character(estfun)) {
		if (!(estfun %in% c("rd_wald","lrr_wald","rd_score","lrr_score","md_t"))) {
			stop("estfun has to be a function or 'rd_wald', 'lrr_wald', 'rd_score', 'lrr_score', or 'md_t'")
		}
		if (estfun %in% c("rd_wald","lrr_wald","rd_score","lrr_score")) {
			if (!is.vector(par01) | length(par01)!=2) {
				stop(paste0("par01 should be a vector of length 2 for estfun ",estfun))
			}
			if (is.character(simfun) && simfun=="norm") {
				warning(paste0("simfun ",simfun," does not work with estfun ",estfun," - binom is used"))
				simfun<-"binom"
			}
		}
		if (estfun %in% c("md_t")) {
			if (!is.data.frame(par01) | any(dim(par01)!=c(2,2)) | any(rownames(par01)!=c("mean","sd"))) {
				stop(paste0("par01 should be a 2x2 data frame with row.names 'mean' and 'sd' for estfun ",estfun))
			}
			if (is.character(simfun) && simfun=="binom") {
				warning(paste0("simfun ",simfun," does not work with estfun ",estfun, " - norm is used"))
				simfun<-"norm"
			}
		}
		
		if (is.character(estfun) && estfun=="rd_wald") {	
			estfun<-function(out, cilevel) waldest(out, cilevel, effm="rd")
			truefun<-function(x) return(diff(par01))
		}
		if (is.character(estfun) && estfun=="lrr_wald") {
			estfun<-function(out, cilevel) waldest(out, cilevel, effm="rr")
			truefun<-function(x) return(diff(log(par01)))
		}
		if (is.character(estfun) && estfun=="rd_score") {
			estfun<-function(out, cilevel) scoreest(out, cilevel,effm="rd")
			truefun<-function(x) return(diff(par01))
		}
		if (is.character(estfun) && estfun=="lrr_score") {
			estfun<-function(out, cilevel) scoreest(out, cilevel,effm="rr")
			truefun<-function(x) return(diff(log(par01)))
		}
		if (is.character(estfun) && estfun=="md_t") {
			estfun<-function(out, cilevel) test(out, cilevel, effm="md")
			truefun<-function(x) return(diff(unlist(par01["mean",])))
		}
	}
	
	#function to simulate data:
	if (is.character(simfun)) {
		if (!(simfun %in% c("binom","norm"))) {
			stop("simfun has to be a function or binom, norm")
		}
		if (is.character(simfun) && simfun=="binom") {
			simfun<-function(n01,par01) {
				out0<-rbinom(n=n01[1], size=1, prob=par01[1])
				out1<-rbinom(n=n01[2], size=1, prob=par01[2])
				return(list(out0,out1))
			}
		}
		if (is.character(simfun) && simfun=="norm") {
			simfun<-function(n01,par01) {
				out0<-rnorm(n=n01[1], mean=par01["mean",1], sd=par01["sd",1])
				out1<-rnorm(n=n01[2], mean=par01["mean",2], sd=par01["sd",2])
				return(list(out0,out1))
			}
		}
	}
	


	if (!is.function(estfun)) {
		stop("estfun is not a function")
	}

	if (!is.function(truefun)) {
		stop("truefun is not a function")
	}
	
	if (!is.function(simfun)) {
		stop("simfun is not a function")
	}
	
	
	#target
	true_pe<-truefun(par01)
		
	#simulate
	out<-simfun(n01, par01)
	
	#estimate
	est<-estfun(out=out, cilevel=cilevel)

	if (direct=="lower") {
		est["z"]<-(-1)*est["z"]
	}
	
	rif<-c(true_pe,n01,est)
	names(rif)<-c("true_pe","fa_n0","fa_n1",paste0("fa_",c(names(est))))

	#interim
	if (nia>0) {
		for (ia in 1:nia) {

			ni0<-pia[ia,1]
			ni1<-pia[ia,2]

			outi0<-out[[1]][1:ni0]
			outi1<-out[[2]][1:ni1]
			outi<-list(outi0,outi1)

			est<-estfun(out=outi, cilevel=cilevel)

			if (direct=="lower") {
				est["z"]<-(-1)*est["z"]
			}

			rii<-c(ni0,ni1,est)

			names(rii)<-paste0("ia",ia,"_",c("n0","n1",names(est)))

			rif<-c(rif,rii)
		}
	}
	return(rif)
}


#Calculates Wald CI and Z
#---------------------

waldest<-function(out, cilevel = 0.95,  effm = c("rd","rr")) {
	
	xs<-unlist(lapply(out,sum))
	ns<-unlist(lapply(out,length))
	phat<-xs/ns
	
	effm<-match.arg(effm)

	#one-sided ci level:
	clev<-1-(1-cilevel)/2	

	if (effm=="rd") {

		rd<-phat[2]- phat[1]

		#wald ci, individual props
		se <- sqrt((phat[1] * (1 - phat[1])) / ns[1] +
				   (phat[2] * (1 - phat[2])) / ns[2])

		pp<-sum(xs)/sum(ns)
		sepool<-sqrt(pp*(1-pp)*(1/ns[1] + 1/ns[2]))

		lci<-rd - qnorm(clev)*se
		uci<-rd + qnorm(clev)*se

		z<-rd/se
		est<-c(xs,phat,rd,lci,uci,z)

	} else {

		lrr <- log(phat[2]/phat[1])

		se <- sqrt(1/xs[1] - 1/ns[1] + 1/xs[2] - 1/ns[2])

		lci <- lrr - qnorm(clev)*se
		uci <- lrr + qnorm(clev)*se

		z<-lrr/se

		est<-c(xs,phat,lrr,lci,uci,z)
	}

	names(est)<-c("x0","x1","p0","p1","pe","lci","uci","z")
	return(est)
}


#Calculates score-based CI and Z
#---------------------

scoreest<-function(out, cilevel,  effm = c("rd","rr")) {
	
	xs<-unlist(lapply(out,sum))
	ns<-unlist(lapply(out,length))
	phat<-xs/ns
	
	effm<-match.arg(effm)

	if (effm=="rd") {
		z<-calc_newcombe_z(x1 = xs[2], n1 = ns[2], x2 = xs[1], n2 = ns[1], delta0 = 0)

		pd<-pairwiseCI::Prop.diff(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
			CImethod="NHS",conf.level = cilevel)

		est<-c(xs,phat,pd$estimate,pd$conf.int,z)
	} else {

		z<-calc_koopman_score_z(x1 = xs[2], n1 = ns[2], x2 = xs[1], n2 = ns[1], phi0 = 1)

		pd<-pairwiseCI::Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
			CImethod="Score",conf.level = cilevel)
		est<-c(xs,phat,log(pd$estimate),log(pd$conf.int),z)
	}
	
	names(est)<-c("x0","x1","p0","p1","pe","lci","uci","z")
	return(est)
}

calc_newcombe_z <- function(x1, n1, x2, n2, delta0 = 0) {
  p1_hat <- x1 / n1
  p2_hat <- x2 / n2

  # Function to find the root where Newcombe's system equals zero for a given z
  newcombe_root <- function(z) {
    z <- abs(z)
    # Wilson score intervals for p1 and p2 at a given z level
    # Group 1
    l1 <- (2*x1 + z^2 - z*sqrt(z^2 + 4*x1*(1 - p1_hat))) / (2*(n1 + z^2))
    u1 <- (2*x1 + z^2 + z*sqrt(z^2 + 4*x1*(1 - p1_hat))) / (2*(n1 + z^2))
    # Group 2
    l2 <- (2*x2 + z^2 - z*sqrt(z^2 + 4*x2*(1 - p2_hat))) / (2*(n2 + z^2))
    u2 <- (2*x2 + z^2 + z*sqrt(z^2 + 4*x2*(1 - p2_hat))) / (2*(n2 + z^2))

    # Evaluate Newcombe boundary equation vs null value (delta0)
	# Newcombe (1998): https://doi.org/10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I
	# CI limits for RD: L = theta_hat  - d1, U = theta_hat + d2
	# with d1 = sqrt((p1_hat - l1)^2 + (u2 - p2_hat)^2) and d2 = sqrt((u1 - p1_hat)^2 + (p2_hat - l2)^2)
	# H0: delta0, set equal to L or U, reject if <> CI limit.
    if ((p1_hat - p2_hat) >= delta0) {
      val <- (p1_hat - p2_hat) - sqrt((p1_hat - l1)^2 + (u2 - p2_hat)^2) - delta0
    } else {
      val <- (p1_hat - p2_hat) + sqrt((u1 - p1_hat)^2 + (p2_hat - l2)^2) - delta0
    }
    return(val)
  }

  # Search for the exact Z value that zeros out the Newcombe equation
  res <- uniroot(newcombe_root, interval = c(0, 20))
  sign_multiplier <- sign(p1_hat - p2_hat - delta0)
  return(res$root * sign_multiplier)
}

calc_koopman_score_z <- function(x1, n1, x2, n2, phi0 = 1) {
  p1_hat <- x1 / n1
  p2_hat <- x2 / n2

  # 1. Coefficients for the Koopman quadratic equation
  a <- phi0 * (n1 + n2)
  b <- -(phi0 * (n2 + x1) + n1 + x2)
  c_coeff <- x1 + x2

  # 2. Solve for the restricted MLE p2_tilde using the quadratic root
  p2_tilde <- (-b - sqrt(b^2 - 4 * a * c_coeff)) / (2 * a)
  p1_tilde <- phi0 * p2_tilde

  # 3. Calculate the exact Koopman Score Z-statistic
  numerator <- p1_hat - (phi0 * p2_hat)
  denominator <- sqrt((p1_tilde * (1 - p1_tilde) / n1) + (phi0^2 * (p2_tilde * (1 - p2_tilde) / n2)))

  return(numerator / denominator)
}


#Calculates t CI and z
#---------------------

test<-function(out, cilevel = 0.95, effm = "md") {
	
	effm<-match.arg(effm)
	
	ms<-unlist(lapply(out,mean))
	sds<-unlist(lapply(out,sd))
	ns<-unlist(lapply(out,length))

	if (effm=="md") {
		
		tt<-t.test(x=out[[2]],y=out[[1]], conf.level = cilevel)
		
		md<-ms[2]- ms[1]

		lci<-tt$conf.int[1]
		uci<-tt$conf.int[2]

		#z<-md/se
		z<-tt$statistic
		
		est<-c(ms,sds,md,lci,uci,z)

	} else {

	}

	names(est)<-c("mean0","mean1","sd0","sd1","pe","lci","uci","z")
	return(est)
}


#set.seed(1)
#par01<-data.frame(matrix(c(5,1,4.89,1),2,2))
#rownames(par01)<-c("mean","sd")
#out<-simfun(n01, par01)
#lapply(out,mean)
#lapply(out,sd)
#test(out)
