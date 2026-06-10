#' ADsim
#'
#' Simulates trials with a binary endpoint and one or several interim analysis,
#'
#' @param n01 Total sample size in each arm, numeric vector of length 2
#' @param p01 Proportion in each arm, numeric vector of length 2
#' @param nia Number of interim analysis, not including the final analysis
#' @param tia Optional vector with the time point(s) of the interim analyses, as fraction of the total (between 0 and 1)
#' @param estfun Function to estimate to determine point estimate, CIs and test.
#'	Should take arguments xs (number of success in both groups), ns (sample size in both groups), and cilevel,
#' 	and return a named vector with pe, lci, uci and z.
#'	Four defaults functions can be entered as string: wald_rd, wald_rr, score_rd, score_rr.
#' @param cilevel Level of the two-sided confience interval for the default functions.
#' @param truefun Function to calculate the true point estimate from p01 (ignored if the default estfun are used)
#' @param direct Direction of the effect, "lower" or "higher" values are better.
#'
#' @returns Numeric vector with true proportions in both arms and point estimate (true_p0, true_p1, true_pe),
#'  and for each interim stage (prefix iax_) and the final stage (prefix fa_),
#'  the simulated number of successes, sample size and observed probabilities for both arms (suffix x0, x1, n0, n1, p0, p1),
#'  point estimate (suffic pe), confidence interval (suffic lci and uci) and z-statistic (z).
#'
#' @export
#'
#' @importFrom stats rbinom qnorm uniroot
#'
#' @examples
#' set.seed(1)
#' ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1)
#'
#' #Earlier interim analysis
#' set.seed(1)
#' ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, tia = 0.4)
#'
#' #Difference effect estimation
#' set.seed(1)
#' ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, estfun = "wald_rr")
#'
#' #Assume higher would be better
#' set.seed(1)
#' ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 2, direct="higher")
#
# User defined function that reproduces the default setting
#'efun<-function(xs,ns,cilevel) {
#'	pe<-diff(xs/ns)
#'	se<-sqrt(sum(xs/ns*(1-xs/ns)/ns))
#'	lci<-pe - qnorm(1-(1-cilevel)/2)*se
#'	uci<-pe + qnorm(1-(1-cilevel)/2)*se
#'	z<-(-pe)/se
#'	res<-c(pe,lci,uci,z)
#'	names(res)<-c("pe","lci","uci","z")
#'	return(res)
#'}
#'set.seed(1)
#'ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, estfun = efun, truefun = function(x) diff(x))
#'
ADsim<-function(n01, p01, nia, tia = NULL,
	estfun = "wald_rd", cilevel = 0.95, truefun = NULL,
	direct = c("lower","higher")) {

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
		if (!(estfun %in% c("wald_rd","wald_rr","score_rd","score_rr"))) {
			stop("estfun has to be a function or wald_rd, wald_rr, score_rd, score_rr")
		}
		if (is.character(estfun) && estfun=="wald_rd") {
			estfun<-function(xs, ns, cilevel) waldest(xs, ns, cilevel, effm="rd")
			truefun<-function(x) return(diff(p01))
		}
		if (is.character(estfun) && estfun=="wald_rr") {
			estfun<-function(xs, ns, cilevel) waldest(xs, ns, cilevel, effm="rr")
			truefun<-function(x) return(diff(log(p01)))
		}
		if (is.character(estfun) && estfun=="score_rd") {
			estfun<-function(xs, ns, cilevel) scoreest(xs, ns, cilevel,effm="rd")
			truefun<-function(x) return(diff(p01))
		}
		if (is.character(estfun) && estfun=="score_rr") {
			estfun<-function(xs, ns, cilevel) scoreest(xs, ns, cilevel,effm="rr")
			truefun<-function(x) return(diff(log(p01)))
		}
	}

	if (!is.function(estfun)) {
		stop("estfun is not a function")
	}

	if (!is.function(truefun)) {
		stop("truefun is not a function")
	}

	outs0<-rbinom(n=n01[1], size=1, prob=p01[1])
	outs1<-rbinom(n=n01[2], size=1, prob=p01[2])

	nam<-c("x0","x1","n0","n1","p0","p1")

	#final
	xs<-c(sum(outs0),sum(outs1))
	ns<-c(length(outs0),length(outs1))
	ps<-xs/ns

	true_pe<-truefun(p01)

	est<-estfun(xs=xs, ns=ns, cilevel=cilevel)

	if (direct=="higher") {
		est["z"]<-(-1)*est["z"]
	}

	rif<-c(p01[1],p01[2],true_pe,xs,ns,ps,est)
	names(rif)<-c("true_p0","true_p1","true_pe",paste0("fa_",c(nam,names(est))))

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

			est<-estfun(xs=xs, ns=ns, cilevel=cilevel)

			if (direct=="higher") {
				est["z"]<-(-1)*est["z"]
			}

			rii<-c(xs,ns,ps,est)

			names(rii)<-paste0("ia",ia,"_",c(nam,names(est)))

			rif<-c(rif,rii)
		}
	}
	return(rif)
}


#Calculates Wald CI and Z
#---------------------

waldest<-function(xs, ns, cilevel,  effm = c("rd","rr"), direct = c("lower","higher")) {

	effm<-match.arg(effm)

	#one-sided ci level:
	clev<-1-(1-cilevel)/2

	phat<-xs/ns

	if (effm=="rd") {

		rd<-phat[2]- phat[1]

		#wald ci, individual props
		se <- sqrt((phat[1] * (1 - phat[1])) / ns[1] +
				   (phat[2] * (1 - phat[2])) / ns[2])

		pp<-sum(xs)/sum(ns)
		sepool<-sqrt(pp*(1-pp)*(1/ns[1] + 1/ns[2]))

		lci<-rd - qnorm(clev)*se
		uci<-rd + qnorm(clev)*se

		z<-(-rd)/se
		est<-c(rd,lci,uci,z)

	} else {

		lrr <- log(phat[2]/phat[1])

		se <- sqrt(1/xs[1] - 1/ns[1] + 1/xs[2] - 1/ns[2])

		lci <- lrr - qnorm(clev)*se
		uci <- lrr + qnorm(clev)*se

		z<-(-lrr)/se

		est<-c(lrr,lci,uci,z)
	}

	names(est)<-c("pe","lci","uci","z")
	return(est)
}


#Calculates score-based CI and Z
#---------------------

scoreest<-function(xs, ns, cilevel,  effm = c("rd","rr")) {

	effm<-match.arg(effm)

	if (effm=="rd") {
		z<-calc_newcombe_z(x1 = xs[1], n1 = ns[1], x2 = xs[2], n2 = ns[2], delta0 = 0)

		pd<-pairwiseCI::Prop.diff(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
			CImethod="NHS",conf.level = cilevel)

		est<-c(pd$estimate,pd$conf.int,z)
	} else {

		z<-calc_koopman_score_z(x1 = xs[1], n1 = ns[1], x2 = xs[2], n2 = ns[2], phi0 = 1)

		pd<-pairwiseCI::Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]),
			CImethod="Score",conf.level = cilevel)
		est<-c(log(pd$estimate),log(pd$conf.int),z)
	}

	names(est)<-c("pe","lci","uci","z")
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


#sig<-scoreest(xs=c(20,10), ns=c(946,1000), cilevel=0.95,  effm = "rd")
#nsig<-scoreest(xs=c(20,10), ns=c(947,1000), cilevel=0.95,  effm = "rd")
#sig["uci"]<0 & (1-pnorm(sig["z"]))<0.025
#nsig["uci"]>0 & (1-pnorm(nsig["z"]))>0.025

#sig<-scoreest(xs=c(20,10), ns=c(957,1000), cilevel=0.95,  effm = "rr")
#nsig<-scoreest(xs=c(20,10), ns=c(958,1000), cilevel=0.95,  effm = "rr")
#sig["uci"]<0 & (1-pnorm(sig["z"]))<0.025
#nsig["uci"]>0 & (1-pnorm(nsig["z"]))>0.025

# scoreest(xs=c(13,13), ns=c(50,50), cilevel=0.95)
# calc_newcombe_z(13, 50, 13, 50)
# x1<-13
# n1<-50
# x2<-13
# n2<-50
# delta0=0


