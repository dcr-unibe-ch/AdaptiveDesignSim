#--------------------------
#Shiny app for adaptve designs
#Lukas Bütikofer
#2025-02-15
#--------------------------


#------------
#packages
#------------

library(shiny)
library(here)
library(pairwiseCI)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(DT)
library(shinyWidgets)

#------------
#functions
#------------

ff_ci<-function(est,dig=NULL,fs="pe (lci to uci)") {
  
  if (!is.na(sum(est))) {
    
    if (is.null(dig)) {
      dig1<-ifelse(est[1]<10,2,ifelse(est[1]<100,1,0))
      dig2<-ifelse(est[2]<10,2,ifelse(est[2]<100,1,0))
      dig3<-ifelse(est[3]<10,2,ifelse(est[3]<100,1,0))
    } else {
      dig1<-dig2<-dig3<-dig
    }
    
    fs<-sub("pe",ff(est[1],dig=dig1),fs)
    fs<-sub("lci",ff(est[2],dig=dig2),fs)
    fs<-sub("uci",ff(est[3],dig=dig3),fs)
    fs
  } else {
    fs<-NA
    fs
  }
}

ff<-function(x,dig=2) {
  formatC(x,format="f",digits=dig)
}

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
	
	
	names(est)<-c("rd","lci","uci","z")
	return(est)
	
}	

#xs<-c(200,230)
#ns<-c(500,500)
#
##risk difference
#waldest(xs, ns, alpha=0.025,  effm = "rd")
#prop.test(x=rev(xs), n=rev(ns), correct = FALSE)
##ci,identical, prop.test uses pooled prop for z.
#waldest(xs, ns, alpha=0.025,  effm = "rd", direct="higher")
#
#
##risk ratio
#rw<-waldest(xs, ns, alpha=0.025,  effm = "rr")
#exp(rw)
#tab <- matrix(c(xs[2],(ns-xs)[2],xs[1],(ns-xs)[1]), nrow = 2, byrow = TRUE)
#epitools::riskratio(tab, method="wald",rev="b")
##the same ci 
#rw
#2*(1-pnorm(rw["z"]))
##z not the same, eptab uses chi-squared


simfun<-function(n01, p0, p1, pia, alpha, effm = c("rd","rr"),direct = c("lower","higher"), 
	cimethod = c("Wald","score")) {

	cimethod<-match.arg(cimethod)
	effm<-match.arg(effm)
	direct<-match.arg(direct)
	
	clevel<-(1-alpha*2)
	
	outs0<-rbinom(n=n01[1], size=1, p=p0)
	outs1<-rbinom(n=n01[2], size=1, p=p1)
	
	nam<-c("x0","x1","n0","n1","rd","rd_lci","rd_uci","z")

	#final
	xs<-c(sum(outs0),sum(outs1))
	ns<-c(length(outs0),length(outs1))
	
	if (cimethod=="Wald") {
		est<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
	}		
	if (cimethod=="score") {	
		if (effm=="rd") {
			pd<-Prop.diff(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]), 
				CImethod="NHS",conf.level = clevel)	
			west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)			
			est<-c(pd$estimate,pd$conf.int,west["z"])
		} else {
			pd<-Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]), 
				CImethod="Score",conf.level = clevel)	
			west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)			
			est<-c(log(pd$estimate),log(pd$conf.int),west["z"])
		}
 
	}
		
	rif<-c(xs,ns,est)
	names(rif)<-paste0("fa_",nam)

	#interim
	if (all(!is.na(pia))) {
		for (ia in 1:length(pia)) {

			ni<-pia[ia]
			ni0<-ni1<-ni %/% 2
			left <- ni - (ni0+ni1)
			left <- sample(c(0,left),replace=FALSE)
			ni0<-ni0 + left[1]
			ni1<-ni1 + left[2]
			stopifnot(ni0+ni1==ni)
			
			outi0<-outs0[1:ni0]
			outi1<-outs1[1:ni1]

			xs<-c(sum(outi0),sum(outi1))
			ns<-c(length(outi0),length(outi1))  
			
			if (cimethod=="Wald") {
				est<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)
			} 
			if (cimethod=="score") {
				if (effm=="rd") {
					pd<-Prop.diff(x = c(xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]), 
						CImethod="NHS",conf.level = clevel)
					west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)		
					est<-c(pd$estimate,pd$conf.int,west["z"])
				} else {
					pd<-Prop.ratio(c(x=xs[2],ns[2]-xs[2]), y=c(xs[1],ns[1]-xs[1]), 
						CImethod="Score",conf.level = clevel)	
					west<-waldest(xs=xs, ns=ns, alpha=alpha, effm=effm, direct=direct)			
					est<-c(log(pd$estimate),log(pd$conf.int),west["z"])				
				}
			}
		
			rii<-c(xs,ns,est)
			
			names(rii)<-paste0("ia",ia,"_",nam)
			
			rif<-c(rif,rii)
		}
	}
	return(rif)
}


#set.seed(1)
#simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025)
#set.seed(1)
#simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025, cimethod="score")
#
#set.seed(1)
#simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025, effm="rr")
#set.seed(1)
#simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025, effm="rr", cimethod="score")



#set.seed(1)
#system.time({
#simdat <- lapply(1:10000, function(i) 
#	simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025))
#simdat <- do.call(rbind, simdat)
#})
#user  system elapsed 
#0.63    0.00    0.63  
			
evalfun<-function(simdat,nia, zsc, zsign, n, direct = c("lower","higher")) {
	
	direct<-match.arg(direct)
	
	resi<-simdat
	
	rdt<-unique(resi$rdt)
	
	resi$rdf<-resi$fa_rd
	resi$rdf_lci<-resi$fa_rd_lci
	resi$rdf_uci<-resi$fa_rd_uci
	
	if (nia>0) {
		sel<-data.frame(matrix(NA,nrow(resi),nia))
		colnames(sel)<-1:nia
		
		for (i in (1:nia)) {
			zs<-zsc[[i]]
			seli<-rep(0,nrow(resi))
			if (zsign=="futility") {
				seli[resi[,paste0("ia",i,"_z")]<zs]<-1
			}
			if (zsign=="efficacy") {
				seli[resi[,paste0("ia",i,"_z")]>zs]<-2
			}
			if (zsign=="both") {
				seli[resi[,paste0("ia",i,"_z")]<zs[1]]<-1
				seli[resi[,paste0("ia",i,"_z")]>zs[2]]<-2
			}
			sel[,i]<-seli
		}
	
		#use rd at interim at the first stopped rd
		resi$fstop<-apply(sel,1,function(x) ifelse(any(x>0),min(which(x>0)),NA))
		resi$neff<-n
		
		i<-1
		for (i in (1:nia)) {
		  si<-!is.na(resi$fstop) & resi$fstop==i
		  resi$rdf[si]<-resi[si,paste0("ia",i,"_rd")]
		  resi$rdf_lci[si]<-resi[si,paste0("ia",i,"_rd_lci")]
		  resi$rdf_uci[si]<-resi[si,paste0("ia",i,"_rd_uci")]
		  
		  niobs<-mean(resi[,paste0("ia",i,"_n0")] + resi[,paste0("ia",i,"_n1")])
		  resi$neff[si]<-niobs
		} 
		
		resi$tstop1<-apply(sel,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))
		resi$tstop2<-apply(sel,1,function(x) ifelse(any(x==2),min(which(x==2)),NA))
		resi$anystop1<-!is.na(resi$tstop1) & (is.na(resi$tstop2) | (resi$tstop1<resi$tstop2))
		resi$anystop2<-!is.na(resi$tstop2) & (is.na(resi$tstop1) | (resi$tstop2<resi$tstop1))
		
		resi$anystop<-"Not stopped"
		resi$anystop[resi$anystop1]<-"Stopped for futility"
		resi$anystop[resi$anystop2]<-"Stopped for efficacy"

	} else {
		resi$anystop<-"Not stopped"
		niobs<-mean(resi[,"fa_n0"] + resi[,"fa_n1"])
		resi$neff<-niobs
	}
	
	if (direct=="lower") {
		resi$pwr<-mean(resi$rdf_uci<0)*100
	} else {
		resi$pwr<-mean(resi$rdf_lci>0)*100
	}
	resi$anystop<-factor(resi$anystop,levels=c("Not stopped","Stopped for futility","Stopped for efficacy"))
	resi$pstop<-mean(resi$anystop %in% c("Stopped for futility","Stopped for efficacy"))*100
	resi$pstopfut<-mean(resi$anystop %in% c("Stopped for futility"))*100
	resi$pstopeff<-mean(resi$anystop %in% c("Stopped for efficacy"))*100
	resi$avn<-mean(resi$neff)
	resi$minn<-min(resi$neff)
	resi$maxn<-max(resi$neff)
	resi$bias<-mean(resi$rdf - rdt)
	resi$sd<-sd(resi$rdf)
	resi$mse<-mean((resi$rdf - rdt)^2)
	resi$coverage<-mean(resi$rdf_lci<rdt & resi$rdf_uci>rdt)*100
	
	return(resi)
}

#simdat <- lapply(1:100, function(i) 
#	simfun(n01=c(250,250), p0=0.5, p1=0.5, pia=c(250), alpha=0.025))
#simdat <- do.call(rbind, simdat)
#simdat<-data.frame(p0t=0.5,p1t=0.5,rdt=0,simdat)
#	
#resf<-evalfun(simdat = simdat, nia = 1, zsc = list(0), zsign = "futility", n=500)



#------------
#input first tab
#------------

ui1 <-  function(id) {
	ns <- NS(id)
	fluidPage(
		sidebarLayout(
			sidebarPanel(
				h3("Assumptions"),
				sliderInput(ns("p0"), "Control proportion (%)",
					min = 0, max = 100, value = 50, step = 1
				),
				radioButtons(
					inputId = ns("direction"),        
					label = "Direction", 
					  choices = c("Lower is better","Higher is better"),
					  selected = "Lower is better",
					  inline = TRUE
				),
				radioButtons(
					inputId = ns("effectm"),        
					label = "Effect measure", 
					  choices = c("Risk difference","Risk ratio"),
					  selected = "Risk difference",
					  inline = TRUE
				),
				#numericInput(ns("rd"), "True risk difference (%)",
				#	min = -100, max = 100, value = 0, step = 1
				#),
				#h5("negative risk difference indicates a benefit of the experimental group"),
				uiOutput(ns("dynamic_slider_css")),
				uiOutput(ns("dynamic_effect")),
				uiOutput(ns("dynamic_effect_text")),
				numericInput(ns("n"), "Sample size",
					min = 0, max = Inf, value = 1000, step = 1
				),
				h3("Interim analysis"),
				numericInput(ns("nia"), "Number of interim analyses",
					min = 0, max = Inf, value = 0, step = 1
				),
				uiOutput(ns("dynamic_tp")),
				#radioButtons(
				#	inputId = ns("zsign"),        
				#	label = "Stop for", 
				#	  choices = c("futility","efficacy","both"),
				#	  selected = "futility" ,
				#	  inline = TRUE
				#),
				uiOutput(ns("dynamic_zsign")),
				#h5("a positive Z indicates a benefit of the experimental group"),
				#h5("< is used for futility stopping, > for efficacy stopping"),
				#h5(strong("Cut-off on the Z statistic")),
				uiOutput(ns("dynamic_text_coff")),
				uiOutput(ns("dynamic_inputs")),
				h3("Final analysis"),
				numericInput(ns("finalz"), "Final analyis at Z of",
				             min = -4, max = 4, value = round(qnorm(0.975),3), step = 1
				),
				textOutput(ns("falpha")),
				h3("Simulation"),
				numericInput(ns("reps"), "Number of simulations",
					min = 0, max = Inf, value = 100, step = 1
				),
				numericInput(ns("seed"), "Seed for random number generation (optional)",
				             min = 0, max = Inf, value = NA, step = 1
				),
				br(),
				actionButton(ns("calc_btn"), "Calculate")
			),
			mainPanel(
				h3("Results from all simulated trials:"),
				h5("Ordered by estimated effect size"),
				plotOutput(ns("fplotH1"),height = "900px"),
				br(),
				h3("Operating characteristics:"),
				DTOutput(ns("tab")),
				br()
			)
		)
	)
}


#----------------------------------------------------------------------------------------------#
#Setup first tab
#----------------------------------------------------------------------------------------------#


server1 <- function(id, out_shared) {
	moduleServer(id, function(input, output, session) {
		
		
	#dynamic inputs 
	#--------
		
		#effects 
		#--------
		output$dynamic_slider_css <- renderUI({	
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				minr<- -input$p0
				maxr<-100 - input$p0
				mi<-0
			} else {
				minr <- 0.01
				maxr <-  floor(100/input$p0*100)/100
				mi<-1
			}
			
			 zero_pos <- (mi - minr) / (maxr - minr) * 100
			 
			if (is.null(input$direction) || input$direction=="Lower is better") {
			 colors <- list(neg = "#ADD8E6", pos = "#FFA500")
			} else {
			  colors <- list(neg = "#FFA500", pos = "#ADD8E6")
			}
   
			tagList(
				tags$style(HTML(sprintf("

				#sliderhalf .irs-bar, .irs-bar-edge {
				  background: transparent !important;
				  background-image: none !important;
				  border: none !important;
				}

				#sliderhalf .irs-line {
				  background: linear-gradient(to right, %s 0%%, %s %s%%, %s %s%%, %s 100%%) !important;
				  background-image: linear-gradient(to right, %s 0%%, %s %s%%, %s %s%%, %s 100%%) !important;
				  border: 1px solid #ccc !important;
				  height: 10px !important;
				  top: 25px !important;
				}

				#sliderhalf .irs-handle i:first-child { background: #000 !important; }", 
					colors$neg, colors$neg, zero_pos, colors$pos, zero_pos, colors$pos,
					colors$neg, colors$neg, zero_pos, colors$pos, zero_pos, colors$pos
				)))
			)
		  
		})	
  
		output$dynamic_effect <- renderUI({
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				#numericInput(
				#	inputId = session$ns("rd"),
				#	label = "True risk difference (%)",
				#	min = -input$p0, max = 100 - input$p0, value = 0, step = 1
				#)	
				div(id = "sliderhalf",
				sliderInput(
					inputId = session$ns("rd"), 
					label = "True risk difference (%)",
					min = -input$p0,
					max =  100 - input$p0, 
					value = 0, step = 1
				))		
			} else {			
				#numericInput(
				#	inputId = session$ns("rd"),
				#	label = "True risk ratio",
				#	min = 0, max = Inf, value = 1, step = 0.1
				#)
				div(id = "sliderhalf",
				sliderInput(
					inputId = session$ns("rd"), 
					label = "True risk ratio",
					min = 0.01,
					max =  floor(100/input$p0*100)/100, 
					value = 1, step = 0.01
				))
			}	
		})
		
		output$dynamic_effect_text <- renderUI({
			if (is.null(input$direction) || input$direction=="Lower is better") {
				tdir<-"<"
			} else {
				tdir<-">"
			}
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				tags$p(h5(paste0("A risk difference ",tdir,
					" 0 indicates a benefit of the experimental group")))
			} else {
				tags$p(h5(paste0("A risk ratio ",tdir," 1 indicates a benefit of the experimental group")))
			}
		})
		
		
		#time points
		#--------
		
		output$dynamic_tp <- renderUI({
			
			req(input$n)
			req(input$nia)			
			
			if (input$nia==0) {
				NULL
			} else {
				
				n <- input$n
				nia <- input$nia

				initial_values <- (seq_len(nia) / (nia + 1))*n
				min_val <- 0
				max_val <- n
		
				noUiSliderInput(
					inputId = session$ns("timepoint"),
					label = "Time point of interim analysis (# patients)",
					min = 0,
					max = n,
					value = initial_values,
					step = 1,
					margin = 0,    
					connect = FALSE,
					format = list(decimals = 0),
					tooltip=TRUE
				)
			}	
		})
		
		#futility vs efficacy
		#---------
		output$dynamic_zsign <- renderUI({
		
			req(input$nia)	
			
			if (input$nia==0) {
				NULL
			} else {
				radioButtons(
					inputId = session$ns("zsign"),        
					label = "Stop for", 
					  choices = c("futility","efficacy","both"),
					  selected = "futility" ,
					  inline = TRUE)
			}
		})
		
		
		#text 
		#---------
		output$dynamic_text_coff <- renderUI({
			
			req(input$nia)
			
			if (input$nia==0) {
				NULL
			} else {
				tags$p(h5(strong("Cut-off on the Z statistic")),
				h5("a positive Z indicates a benefit of the experimental group"))
			}
		})
  
		
		#input for each IA:	
		#--------
		
		output$dynamic_inputs <- renderUI({
			
			req(input$nia)
			req(input$zsign)
			
			if (input$nia==0) {
				NULL
			} else {
				
				if (input$zsign=="futility") {
					connect_vec <-  c(TRUE, FALSE)
					#connect_cols<-rgb(1, 0, 0, 0.4)
					initial_values<-0
					slider_css <- "
					#slider-wrapper .noUi-connect:nth-of-type(1) { background: rgba(255, 0, 0, 0.4) !important; }
					"
				} 
				if (input$zsign=="efficacy") {
					connect_vec <-  c(FALSE, TRUE)
					#connect_cols<-rgb(0, 0, 1, 0.4)
					initial_values<-qnorm(0.975)
					slider_css <- "
					#slider-wrapper .noUi-connect:nth-of-type(1) { background: rgba(0, 153, 0, 0.4) !important; }
					"
				}
				if (input$zsign=="both") {
					initial_values<-c(0,1.96)
					connect_vec <- c(TRUE, FALSE, TRUE)
					#connect_cols<-c(rgb(1, 0, 0, 0.4),rgb(0, 0, 1, 0.4))
					slider_css <- "
					#slider-wrapper .noUi-connect:nth-of-type(1) { background: rgba(255, 0, 0, 0.4) !important; }
					#slider-wrapper .noUi-connect:nth-of-type(2) { background: rgba(0, 153, 0, 0.4) !important; }
					"
				}
				
				lapply(seq_len(input$nia), function(i) {
					tagList(
						#sliderInput(
						#	inputId = session$ns(paste0("zval_", i)),
						#	label   = paste0("at interim analysis ", i),
						#	min = -4, max = 4, value = 0, step = 0.05
						#),
						#sliderTextInput(
						#  inputId = session$ns(paste0("zval_", i)),
						#  label = paste0("Z at interim analysis ", i),
						#  choices = round(sort(c(seq(-4,4,by=0.1),
						#	qnorm(c(0.025,0.05,0.1,0.9,0.95,0.975)))),2),
						#  selected = 0
						#),
						#numericInput(
						#  inputId = session$ns(paste0("zval_", i)),
						#  label = paste0("Z at interim analysis ", i),
						#   min = -4, max = 4, 
						#   value = 0, step = 1
						#),
						tags$style(HTML(slider_css)),
						div(id = "slider-wrapper",
						noUiSliderInput(
							inputId = session$ns(paste0("zval_", i)),
							label = HTML(paste0("&nbsp at interim analysis ", i)),
							min = -4,
							max = 4,
							value = initial_values,
							step = 0.01,
							connect = connect_vec, 
							tooltips = TRUE,
							width = "100%",
							#color = connect_cols 
						)),
						plotOutput(
							outputId = session$ns(paste0("zdist_", i)),
							height = 100
						),
						tags$hr()
					)	
				})
			}	
		})
		
		#for a separate plot item
		#output$dynamic_plots <- renderUI({
		#  req(input$nia)
		#
		#  lapply(seq_len(input$nia), function(i) {
		#	plotOutput(
		#	  outputId = session$ns(paste0("zdist_", i)),
		#	  height = 100
		#	)
		#  })
		#})
	
		#standard normal plots
		#-----------------
		
		observe({		
			
			req(input$zsign)
			req(input$nia)
			
			lapply(seq_len(input$nia), function(i) {
				local({
					ii <- i
					output[[paste0("zdist_", ii)]] <- renderPlot({

						id <- paste0("zval_", ii)
						zs<-input[[id]]	
						req(zs)
						
						#zsn<-as.numeric(zs)
						
						xlim<-c(-4,4)
						adx<-0.5
						x <- seq(xlim[1], xlim[2], l = 500)
						y <- dnorm(x)
						ylim<-c(0,max(y))
						colt<-c(rgb(1, 0, 0, 1),rgb(0, 0.6, 0, 1))
						cols<-c(rgb(1, 0, 0, 0.4),rgb(0, 0.6, 0, 0.4))
						
						par(mar=c(2,0,0,0)) 
						plot(0,type="n",axes=FALSE,
							ylim=ylim,xlim=xlim,ylab="",xlab="")
						axis(side=1,pos=0,mgp=c(0,0.6,0),cex.axis=0.8)
						
						lines(x, y)
						
						if (input$zsign %in% c("futility","both")) {
							if (input$zsign=="both") {
								zsi<-zs[1]
							} else {
								zsi<-zs
							}						
							xs <- x[x <= zsi]
							ys <- dnorm(xs)
							polygon(c(xs, rev(xs)),c(rep(0, length(xs)), rev(ys)),
								col = cols[1],border = NA)
							pzs<-paste0(round(pnorm(zsi)*100),"%")	
							text(x = xlim[1]+adx,y=mean(ylim),labels=pzs,
								adj=c(0,0.5),cex=0.8, col=colt[1])
						}
						if (input$zsign %in% c("efficacy","both")) {
							if (input$zsign=="both") {
								zsi<-zs[2]
							} else {
								zsi<-zs
							}	
							xs <- x[x >= zsi]
							ys <- dnorm(xs)
							polygon(c(xs, rev(xs)),
								c(rep(0, length(xs)), rev(ys)),
								col = cols[2],border = NA)
							pzs<-paste0(round((1-pnorm(zsi))*100),"%")	
							text(x = xlim[2]-adx,y=mean(ylim),labels=pzs,
								adj=c(1,0.5),cex=0.8, col=colt[2])	
						}
					  
						for (zi in 1:length(zs)) {
							#abline(v = zs, col = "black", lwd = 1, lty = 2)
							lines(x = c(zs[zi],zs[zi]), y=c(0,dnorm(zs[zi])), lty = 2)
						}
					})
				})
			})
		})
		
		
		#n per group
		#-----------
		
		n01<-reactive({
			n0<-n1<-input$n %/% 2
			left <- input$n - (n0+n1)
			left <- sample(c(0,left),replace=FALSE)
			n0<-n0 + left[1]
			n1<-n1 + left[2]
			stopifnot(n0+n1==input$n)
			return(c(n0,n1))
		})	
		
		rdused<-reactive({
		    if (is.na(input$rd)) {
		        rdused<-0
		    } else {
		        rdused<-input$rd
		    }
		})
		
		#alpha output
		#-----------		
		
		output$falpha <- renderText({
			paste0("Corresponds to an one-sided alpha of ",signif(1-pnorm(input$finalz),2))
		})
		
		
		
		#check: observe 
		#-----------
		
		#observe({
		#req(input$timepoint)
		#print(input$timepoint)
		#req(input$nia)
		#print(input$nia)
		#req(input$zsign)
		#print(input$nia)
		#})
		

		
	#simulation
	#----------
		
		#res<-reactive({
		res<-eventReactive(input$calc_btn, {

		  #browser()
			req(input$nia)
			
			if (!is.na(input$seed)) {
				set.seed(input$seed)
			}
			
			if (is.null(input$direction) || input$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				effm<-"rd"
				rdl<-rdused()/100
			} else {
				effm<-"rr"
				rdl<-log(rdused())
			}

			if (input$finalz==1.96) {
				usefinalz<-qnorm(0.975)
			} else {
				usefinalz<-input$finalz
			}
			
			alpha<-1-pnorm(usefinalz)
			
			res<-vector(length=length(rdl),mode="list")
			names(res)<-rdl
		
			if (input$nia==0) {
				pia<-NA
			} else {
			  req(input$timepoint)
			  pia<-input$timepoint
				stopifnot((length(pia)==input$nia))
			}
			
			for (r in 1:length(rdl)) {

				rdt<-rdl[r]
				
				if (is.null(input$effectm) || input$effectm=="Risk difference") {
					p1<-input$p0/100 + rdt
				} else {
					p1<-input$p0/100 * exp(rdt)
				}

				ri <- lapply(1:input$reps, function(x) 
					simfun(n01 = n01(), p0 = input$p0/100, p1 = p1, pia = pia, 
						alpha = alpha, effm = effm, direct = direct))
				ri <- do.call(rbind, ri)
				
				res[[r]]<-data.frame(p0t=input$p0/100,p1t=p1,rdt=rdt,ri)
			}
			
			return(res)
		})
		
		
	#final rd and CI at a specific threshold
	#----------------------------------
		
		#resf<-reactive({
		resf<-eventReactive(input$calc_btn, {
					
			req(input$nia)
			
			if (is.null(input$direction) || input$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			#z values 
			zsc<-vector(length=input$nia,mode="list")
			if (input$nia>0) {
			  for (i in (1:input$nia)) {	
					id <- paste0("zval_", i)
					zsc[[i]]<-input[[id]]
				}
			}
			
			resf<-res()
			
			for (ri in 1:length(res())) {
				
				resi<-res()[[ri]]
				
				resf[[ri]]<-evalfun(simdat = resi, nia = input$nia, 
					zsc = zsc, zsign = input$zsign, n = input$n,
					direct = direct)

			}
			return(resf)

		})	
		

	#plot 
	#-------------

		output$fplotH1<-renderPlot({
			
			#browser()
			if (is.null(input$direction) || input$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				tfac<-100
				rdl<-rdused()
				xnam<-"Risk difference (%)"
				ref<-0
			} else {
				tfac<-1
				rdl<-log(rdused())
				xnam<-"Risk ratio"
				ref<-1
			}
			
			erdf<-mean(resf()[[length(resf())]]$rdf)*tfac
			ns<-nrow(resf()[[length(resf())]])
			if (direct=="lower") {
				nsig<-sum(resf()[[length(resf())]]$rdf_uci<0)		
			} else {
				nsig<-sum(resf()[[length(resf())]]$rdf_lci>0)		
			}
			
			ids<-(1:ns) / (ns + 1)*100
			psig<-mean(c(ids[nsig],ids[nsig+1]))
			if (nsig==ns) {
				psig<-100  
			}
			
			labnp<-c("Negative trials*","Positive trials*")
			labfav<-c("← Favors Experimental","Favors Control →")
			labci<-"the upper limit of the 95% CI is < "
			if (direct!="lower") {
				psig<-100-psig
				labnp<-rev(labnp)
				labfav<-c("← Favors Control","Favors Experimental →")
				labci<-"the lower limit of the 95% CI is > "
			}
			
			#colors for error bars
			nsim<-nrow(resf()[[length(resf())]])
			map_log <- function(x, xmin=100, xmax=10000, ymin=0.8, ymax=0.05) {
				ymin + (ymax - ymin) *
				(log(x) - log(xmin)) /
				(log(xmax) - log(xmin))
			}		
			cola<-map_log(nsim)
			cols<-c(rgb(0.1,0.1,0.1,cola),rgb(1,0,0,cola),rgb(0,0.6,0,cola))
			
			#color for background 
			colbg<-c(rgb(0.90, 0.49, 0.13, alpha = 1),rgb(0.20, 0.60, 0.86, alpha = 1))
			if (direct!="lower") {
				colbg<-rev(colbg)
			}
			
			alphabg<-0.1
			thjust<-1.1
			txlabpos<-Inf
			
			tsize<-16
			bgpos<-psig
			tylabpos<-psig
			
			if (psig<2.2) {
			  tylabpos<-2.2
			} 
			if (psig>97.8) {
			  tylabpos<-97.8
			} 
			if ((direct=="lower" & nsig/ns>0.9) | (direct=="higher" & nsig/ns<0.1)) {
			  thjust<--0.1
			  txlabpos<-(-Inf)
			}
			if ((direct=="lower" & nsig==0) | (direct=="higher" & nsig==ns)) {
			  tylabpos<-0
				bgpos<-0
				colbg[2]<-NA
			}
			if ((direct=="lower" & nsig==ns) | (direct=="higher" & nsig==0)) {
				tylabpos<-100
				bgpos<-100 
				colbg[1]<-NA	
			}
			
			# Example vertical lines
			vlines <- data.frame(
			  x = c(rdl, erdf),
			  label = c("Truth", "Estimate"),
			  linetype = c("solid", "dashed"),
			  color = c("magenta", "magenta")
			)
			
			
			fpl<-resf()[[length(resf())]] %>%
				arrange(rdf) %>%
				mutate(id = ids) %>%
				mutate(across(starts_with("rdf"), ~ .x*tfac)) %>%
				mutate(rdfs = mean(rdf)) %>%
				ggplot(aes(x = rdf, y = id, group = anystop,color = anystop)) +
					annotate("rect",
						xmin = -Inf, xmax = Inf,
						ymin = bgpos, ymax = Inf,
						fill =   colbg[1],alpha = alphabg) +
					annotate("rect",
						xmin = -Inf, xmax = Inf,
						ymin = -Inf, ymax = bgpos,
						fill = colbg[2],alpha = alphabg) +
					annotate("text", x = txlabpos, y = tylabpos, #y = (100-psig)/2+psig,  
						label = labnp[1],
						hjust = thjust, vjust = -0.4, size = tsize/3, 
						color = colbg[1]) +
					annotate("text", x = txlabpos, y = tylabpos , # y = psig/2,
						label = labnp[2],
						hjust = thjust, vjust = 1.2, size = tsize/3, 
						color =  colbg[2]) +		
				geom_point(size = 1, alpha=cola, show.legend = c(color = TRUE, linetype = FALSE)) +
				geom_errorbarh(aes(xmin = rdf_lci, xmax = rdf_uci), height = 0.2, show.legend = c(color = TRUE, linetype = FALSE)) +
				xlab(xnam) +
				ylab("Proportion of trials (%)") +
			    scale_y_continuous(limits=c(0,100), breaks=seq(0,100,by=20), expand = expansion(mult = c(0, 0))) + 
				#geom_vline(xintercept = rdl, linetype = "solid", colour="magenta") +
			    #geom_vline(xintercept = erdf, linetype = "dashed", colour="magenta") +
				geom_vline(aes(xintercept = erdf, linetype = "Mean estimate"), color = "magenta") +
				geom_vline(aes(xintercept = rdl, linetype = "Truth"), color = "magenta") +
				scale_linetype_manual(name = "Vertical Lines", 
					values = c("Mean estimate" = "dashed", "Truth" = "solid"), drop = FALSE) +
				#scale_color_manual(name = "Legend Title", values = c(anystop = cols)) +
				scale_color_manual(values = cols, drop = FALSE) +
				theme_bw(base_size = tsize) +
				theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10),
					axis.title.x = element_text(margin = margin(t = 20)), 
					legend.title = element_blank(),
					legend.position = "top",
					legend.justification = "right",
					legend.box = "vertical",      
					legend.box.just = "right",
					legend.margin = margin(t = 0, b = 0),
					legend.box.margin = margin(t = -10, b = -10)) + 
				guides(color = guide_legend(nrow=1, order = 2,
						override.aes = list(shape = 16, linetype = "solid",  alpha = 1, color = cols)), 
					linetype = guide_legend(nrow=1, order = 1, override.aes = list(shape = NA ))) +
				annotate("text", x = 0, y = -Inf, label = labfav[1], 
					vjust = 3.5, hjust = 1.1, size = 4, fontface = "italic") +
				annotate("text", x = 0, y = -Inf, label = labfav[2], 
					vjust = 3.5, hjust = -0.1, size = 4, fontface = "italic") +
				coord_cartesian(clip = "off") +
				labs(caption = paste0("*based on whether ",labci , ref))
			
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				fpl
			} else {
				xlims <- c(min(resf()[[length(resf())]]$rdf_lci),max(resf()[[length(resf())]]$rdf_uci))
				original_breaks <- pretty(exp(xlims), n = 6)       
				breaks_log <- log(original_breaks)                 
				fpl +  scale_x_continuous(breaks = log(original_breaks), labels = original_breaks)
			}
			
			#browser()
		}) %>% bindEvent(input$calc_btn)
		
		
		
	#table
	#----------

		output$tab <- renderDT({
			
			#browser()
			#req(input$zsign)
			req(input$nia)
			
			resi<-resf()[[length(resf())]]
			
			
			if (is.null(input$effectm) || input$effectm=="Risk difference") {
				
				tfac<-100
				dig<-1
				
				addp<-"%"
				labeff<-"True risk difference (RD)"
				labest<-"Estimated RD (95% CI)"
				labsd<-"Standard deviation of estimated RD"
				labmse<-"Mean squared error of RD (in %)"
				if (rdused()==0) {
					nap<-"Type I error (one-sided)"
				} else {
					nap<-"Power"
				}
				resi$rdtf<-paste0(ff(tfac*resi$rdt, dig=dig),addp)
				resi$est<-ff_ci(c(apply(resi[,c("rdf","rdf_lci","rdf_uci")]*tfac,2,mean)),
					dig=dig,fs="pe% (lci to uci%)")
				
			} else {
				tfac<-1
				dig<-2
				addp<-""
				labeff<-"True risk ratio (RR)"
				labest<-"Estimated RR (95% CI)"
				labsd<-"Standard deviation of estimated RR"
				labmse<-"Mean squared error of RR"
				
				if (rdused()==1) {
					nap<-"Type I error (one-sided)"
				} else {
					nap<-"Power"
				}
				
				resi$rdtf<-paste0(ff(exp(resi$rdt), dig=dig),addp)
				resi$est<-ff_ci(exp(c(apply(resi[,c("rdf","rdf_lci","rdf_uci")]*tfac,2,mean))),
					dig=dig,fs="pe (lci to uci)")
				
			}
		  
			nas<-c("Stopped","Stopped for futility","Stopped for efficacy")
				
			zf<-""
			if (input$nia>0) {
				
			  req(input$zsign)
			  
				if (input$zsign=="futility") {
					zsign<-"<"
					#nas<-paste0("Stopped for ",input$zsign)
				} 
				if (input$zsign=="efficacy") {
					zsign<-">"
					#nas<-paste0("Stopped for ",input$zsign)
				}
				if (input$zsign=="both") {
					zsign<-c("<",">")
					#nas<-c("Stopped","Stopped for futility","Stopped for efficacy")
				}
							
				#stop at each IA:
				zf<-paste0(zsign,ff(input$zval_1, dig=2),collapse=" & ")
				s1<-!is.na(resi$tstop1) & resi$tstop1==1 & resi$anystop1			
				s2<-!is.na(resi$tstop2) & resi$tstop2==1 & resi$anystop2
				pstopi1<-paste0(ff(mean(s1)*100, dig=1),"%")		
				pstopi2<-paste0(ff(mean(s2)*100, dig=1),"%")
				pstopi<-paste0(ff(mean(s1 | s2)*100, dig=1),"%")
				
				if (input$nia>1) {
					for (i in (2:input$nia)) {
						id <- paste0("zval_", i)
						zi <- input[[id]]
						zf<-paste0(zf,", ",paste0(zsign,ff(input$zval_1, dig=2),collapse=" & "))
						
						s1<-!is.na(resi$tstop1) & resi$tstop1==i & resi$anystop1
						s2<-!is.na(resi$tstop2) & resi$tstop2==i & resi$anystop2
						
						pstopi1<-paste0(pstopi1,", ",ff(mean(s1)*100, dig=1),"%")
						pstopi2<-paste0(pstopi2,", ",ff(mean(s2)*100, dig=1),"%")
						pstopi<-paste0(pstopi,", ",ff(mean(s1 | s2)*100, dig=1),"%")
					}
				}	
			}	
			
			resi$pwrf<-paste0(ff(mean(resi$pwr), dig=1),"%")
			resi$pstopf <- paste0(ff(mean(resi$pstop), dig=1),"%")
			resi$pstopfutf <- paste0(ff(mean(resi$pstopfut), dig=1),"%")
			resi$pstopefff <- paste0(ff(mean(resi$pstopeff), dig=1),"%")
				
			if (input$nia>1) {
				resi$pstopf<-paste0(resi$pstopf," (",pstopi,")")
				resi$pstopfutf<-paste0(resi$pstopfutf," (",pstopi1,")")
				resi$pstopefff<-paste0(resi$pstopefff," (",pstopi2,")")
				
				nas<-paste0(nas," (at each interim analysis)")
			}
					
			d1<-resi %>%
				#mutate(zsfor = zf) %>%
				mutate(pwrfor = pwrf) %>%
				mutate(pstopfor = pstopf) %>%
				mutate(pstopfutfor = pstopfutf) %>%
				mutate(pstopfeffor = pstopefff) %>%
				mutate(avnfor = paste0(round(avn)," (",round(minn),", ",round(maxn),")")) %>%
				mutate(rdfor=rdtf) %>%
				mutate(estfor=est) %>%
				mutate(biasfor = paste0(ff(tfac*bias, dig=dig),addp)) %>%
				mutate(sdfor = paste0(ff(tfac*sd, dig=dig),addp)) %>%
				mutate(msefor = ff(tfac^2*mse, dig=dig^2)) %>%
				mutate(coveragefor = paste0(ff(coverage, dig=1),"%")) %>%
				select(ends_with("for")) %>%
				unique()
				
			#if (input$zsign!="both") {
			#	d1<-d1 %>% select(-pstopfutfor,-pstopfeffor)
			#}
			
			tb<-cbind(c(
				#"Stop if Z statistic",
				nap,
				nas,
				"Sample size (expected, min, max)",
				labeff,
				labest,
				"Bias",
				labsd,
				labmse,
				"Coverage of 95% CI"),
				t(d1))
					
			colnames(tb)<-c(" "," ")
					
			tb %>%
				datatable(
				options = list(searching = FALSE,ordering = FALSE,
					paging = FALSE, info = FALSE),
					rownames = FALSE)
					
			#browser()		
		}) %>% bindEvent(input$calc_btn)
		
		
#share inputs to 2nd tab
#-----------------

    return(reactive({
		
		req(input$nia)
		
		zsc<-vector(length=input$nia,mode="list")
		if (input$nia) {
			for (i in (1:input$nia)) {	
				id <- paste0("zval_", i)
				zsc[[i]]<-input[[id]]
			}
		}
		
		list(
			p0 = input$p0,
			rd = input$rd,
			n = input$n,
			n01 = n01(),
			nia = input$nia,
			timepoint = input$timepoint,
			zsign = input$zsign,
			finalz = input$finalz,
			reps = input$reps,
			zsc = zsc,
			effectm = input$effectm,
			direction = input$direction
		)
    }))
		
		
		
	})
}

#----------------------------------------------------------------------------------------------#
#Second tab, over alternatives
#----------------------------------------------------------------------------------------------#


ui2 <- function(id) {
	ns <- NS(id)
	fluidPage(
		sidebarLayout(
			sidebarPanel(
				h3("Assumptions"),
				#textOutput(ns("p0fw")),
				uiOutput(ns("ass1")),
				br(),
				uiOutput(ns("slider_ui")),
				#h5("negative risk difference indicates a benefit of the experimental group"),
				uiOutput(ns("dynamic_effect_text")),
				numericInput(
					ns("step"),
					"Number of steps:",
					value = 10,
					min   = 1),	
				uiOutput(ns("dynamic_step_size")),
				h3("Interim analysis"),	
				uiOutput(ns("ia1")),
				h3("Final analysis"),	
				uiOutput(ns("an1")),				
				h3("Simulation"),
				numericInput(ns("reps"), "Number of simulations",
					min = 0, max = Inf, value = 100, step = 1
				),
				numericInput(ns("seed"), "Seed for random number generation (optional)",
				             min = 0, max = Inf, value = NA, step = 1
				),
				br(),
				actionButton(ns("calc_btn"), "Calculate")
			),
			mainPanel(
				fluidRow(
				column(width=6,
					h4("Power"),
					plotOutput(ns("ppwr"))),
				column(width=6,
					h4("Proportion stopped"),
					plotOutput(ns("pstop"))),
				column(width=6,
					h4("Bias"),
					plotOutput(ns("pbias"))),
				column(width=6,
					h4("Sample size"),
					plotOutput(ns("pss")))
				)	
			)
		)
	)
}


server2 <- function(id, shared_inputs) {
	moduleServer(id, function(input, output, session) {
		
			
		output$slider_ui <- renderUI({
			
			input1<-shared_inputs()
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				values<- c(-15, 5)
				valuesrr<- c(0.7, 1.1)
			} else {
				values<- c(-5, 15)
				valuesrr<- c(0.9, 1.3)
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				sliderInput(
					inputId=session$ns("rdrange"),
					label = "Range of alternatives:",
					min = - input1$p0,
					max = 100 - input1$p0,
					step = 1,
					value = values)
			} else {
				sliderInput(
					inputId=session$ns("rdrange"),
					label = "Range of alternatives:",
					min = 0.01,
					max = floor(100/input1$p0*100)/100,
					step = 0.01,
					value = valuesrr)
			}	
		})
		
		output$dynamic_effect_text <- renderUI({
		
			input1<-shared_inputs()
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				tdir<-"<"
			} else {
				tdir<-">"
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				#tags$p(h5("A risk difference < 0 indicates a benefit of the experimental group"))
				tags$p(h5(paste0("A risk difference ",tdir,
					" 0 indicates a benefit of the experimental group")))
			} else {
				#tags$p(h5("A risk ratio < 1 indicates a benefit of the experimental group"))
				tags$p(h5(paste0("A risk ratio ",tdir," 1 indicates a benefit of the experimental group")))
			}
		})
		
		#output$dynamic_step_size <- renderUI({
		#
		#input1<-shared_inputs()
		#	
		#	if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
		#		numericInput(
		#			inputId=session$ns("step"),
		#			label = "Step size:",
		#			value = 1,
		#			min   = 0.1)
		#			
		#	} else {
		#		numericInput(
		#			inputId=session$ns("step"),
		#			label = "Step size:",
		#			value = 0.1,
		#			min   = 0.01)			
		#	}
		#})

		
		output$ass1 <- renderUI({
		  input1<-shared_inputs()
		  tagList(
				#hr(),
				p(strong("Control proportion: "), input1$p0, "%"),
				p(strong("Direction:"), input1$direction),
				p(strong("Effect measure: "), input1$effectm),
				p(strong("Sample size: "), input1$n)
		  )
		})
		
		output$ia1 <- renderUI({
			
			input1<-shared_inputs()
			
			req(input1$zsign)
			
			if (input1$zsign=="both") {
				lc<-lapply(input1$zsc,function(x) paste(x,collapse="/"))
				zsc1<-paste(unlist(lc),collapse=", ")
				txt1<-"at Z </>"
			} else {
				zsc1<-paste(unlist(input1$zsc),collapse=", ")
				if (input1$zsign=="futility") {
					txt1<-paste0("at Z < ")
				} else {
					txt1<-paste0("at Z > ")
				}
			}			
			
			tagList(
				p(strong("Number: "), input1$nia),
				p(strong("Time point: "), 
					paste(input1$timepoint,collapse=", ")),
				p(strong("Stopping for: "),input1$zsign), 
				p(strong(txt1),zsc1)
			)
		})
		
		output$an1 <- renderUI({
		  input1<-shared_inputs()
		  tagList(
				paste0("At a Z of ", input1$finalz,
					", which corresponds to a one-sided alpha of ",signif(1-pnorm(input1$finalz),2))
		  )
		})
		
		
		observe({
			input1<-shared_inputs()
		
		#	print(input1$p0)
		#	print(input1$rd)
		#	print(input1$n)
		#	print(input1$nia)
		#	print(input1$n01)
		#	
		#	req(input$rdrange)
		#	print(input$rdrange)
		#	req(input1$zsc)
		#	print(input1$zsc)
		#	print(input1$effectm)
		})

		
		#simulate over range of alternatives
		#-------
		res<-eventReactive(input$calc_btn, {
			
			input1<-shared_inputs()
			
			req(input1$nia)
			req(input$rdrange)
			
			if (!is.na(input$seed)) {
				set.seed(input$seed)
			}
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				effm<-"rd"
				rdl<-seq(min(input$rdrange),max(input$rdrange),l=input$step)/100
			} else {
				effm<-"rr"
				rdl<-seq(log(min(input$rdrange)),log(max(input$rdrange)),l=input$step)
			}
		
			if (input1$finalz==1.96) {
				usefinalz<-qnorm(0.975)
			} else {
				usefinalz<-input1$finalz
			}
			
			alpha<-1-pnorm(usefinalz)
			
			res<-vector(length=length(rdl),mode="list")
			names(res)<-rdl
		
			if (input1$nia==0) {
				pia<-NA
			} else {
			  req(input1$timepoint)
			  pia<-input1$timepoint
				stopifnot((length(pia)==input1$nia))
			}
			
			for (r in 1:length(rdl)) {
		
				rdt<-rdl[r]

				if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
					p1<-input1$p0/100 + rdt
				} else {
					p1<-input1$p0/100 * exp(rdt)
				}
				
				ri <- lapply(1:input$reps, function(x) 
					simfun(n01 = input1$n01, p0 = input1$p0/100, p1 = p1, pia = pia, 
						alpha = alpha, effm = effm, , direct = direct))
				ri <- do.call(rbind, ri)
				
				res[[r]]<-data.frame(p0t=input1$p0/100,p1t=p1,rdt=rdt,ri)
			}
			
			return(res)

		})
				
	#evaluate at threshholds 
	#-------------------
		
		resf<-eventReactive(input$calc_btn, {
				
			input1<-shared_inputs()		
			
			#req(input1$zsign)
			req(input1$nia)
			req(input1$n)
			#req(input1$zsc)
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			outi<-res()
			
			for (ri in 1:length(res())) {
				
				resi<-res()[[ri]]
				
				resii<-evalfun(simdat = resi, nia = input1$nia, 
					zsc = input1$zsc, zsign = input1$zsign, n = input1$n, 
					direct = direct)
					
				ress<-resii %>%
					select(pwr, pstop, pstopfut, pstopeff,
						avn, minn, maxn, bias, sd)
				
				outi[[ri]]<-ress[1,]

			}
			out<-do.call(rbind,outi)
			out<-cbind(rdt=as.numeric(names(res())),out)
			
			return(out)

		})	
		
		#observe({
		#	print(names(res()))	
		#	print(nrow(resf()))
		#	print(head(resf()))
		#}) 


		#plot
		#-------
		output$ppwr<-renderPlot({
			
			input1<-shared_inputs()		
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			tsize<-16
			colbg<-c(rgb(0.90, 0.49, 0.13, alpha = 1),rgb(0.20, 0.60, 0.86, alpha = 1))
			dlab<-c("Control better","Experimental better")
			alphabg<-0.1
			
			if (direct!="lower") {
				colbg<-rev(colbg)
				dlab<-rev(dlab)
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				tfac<-100
				xnam<-"True risk difference (%)"
			} else {
				tfac<-1
				xnam<-"True risk ratio"
			}
			
			ppwr<-resf()  %>%	
				mutate(rdt = rdt*tfac) %>%
				ggplot(aes(x=rdt,y=pwr, group = 1)) +
				annotate("rect",
						xmin = -Inf, xmax = 0,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[2],alpha = alphabg) +
				annotate("rect",
						xmin = 0, xmax = Inf,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[1],alpha = alphabg) +		
				geom_line() +
				geom_hline(yintercept = 80, linetype = "dashed") +
				annotate("text",x = -Inf, y = 80,label = "80%", hjust = -0.2, vjust = -0.5) +
				geom_hline(yintercept = 2.5, linetype = "dashed") +
				annotate("text",x = -Inf, y = 2.5,label = "2.5%", hjust = -0.2, vjust = -0.5) +
				xlab(xnam) +
				ylab("Power (%)") + 
				scale_y_continuous(breaks = seq(0,100,by=20), limits=c(0,100)) +
				coord_cartesian(clip = "off") +
				theme_bw(base_size = tsize)
			
			xsc<-c(min(resf()$rdt),max(resf()$rdt))*tfac
			plow<--(min(xsc)/(xsc[2]-xsc[1]))
			pupp<-max(xsc)/(xsc[2]-xsc[1])
	
			if (plow>0.1) {
				ppwr<-ppwr + 
					annotate("text", x = 0, y = 100,
					label = dlab[2],
					hjust = 1.1, vjust = -0.4, 
					color = colbg[2])
			}
			if (pupp>0.1) {
				ppwr<-ppwr + 
					annotate("text", x = 0, y = 100,
					label = dlab[1],
					hjust = -0.1, vjust = -0.4, 
					color = colbg[1])
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				ppwr		
			} else {
				original_breaks <- pretty(exp(xsc), n = 6)     
				breaks_log <- log(original_breaks)                 
				ppwr +  scale_x_continuous(breaks = log(original_breaks), 
					labels = original_breaks,
					limits = xsc)
			}
			
		})	%>% bindEvent(input$calc_btn)

		
		output$pstop<-renderPlot({
			
			input1<-shared_inputs()	
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			df<-data.frame(
				rdt=rep(resf()$rdt,3),
				pstop=c(resf()$pstop,resf()$pstopfut,resf()$pstopeff),
				type=rep(c("Overall","Futility","Efficacy"),each=nrow(resf())))
				
			df$type<-factor(df$type,levels=c("Overall","Futility","Efficacy"))
			
			cola<-1
			cols<-c(rgb(0.1,0.1,0.1,cola),rgb(1,0,0,cola),rgb(0,0.6,0,cola))
			tsize<-16	
			colbg<-c(rgb(0.90, 0.49, 0.13, alpha = 1),rgb(0.20, 0.60, 0.86, alpha = 1))
			dlab<-c("Control better","Experimental better")
			alphabg<-0.1
			
			if (direct!="lower") {
				colbg<-rev(colbg)
				dlab<-rev(dlab)
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				tfac<-100
				xnam<-"True risk difference (%)"
			} else {
				tfac<-1
				xnam<-"True risk ratio"
			}
			
			pstop<-df  %>%
				mutate(rdt = rdt*tfac) %>%
				ggplot(aes(x=rdt,y=pstop, colour=type, group=type)) +
				annotate("rect",
						xmin = -Inf, xmax = 0,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[2],alpha = alphabg) +
				annotate("rect",
						xmin = 0, xmax = Inf,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[1],alpha = alphabg) +
				geom_line(show.legend=TRUE) +
				xlab(xnam) +
				ylab("Stopped (%)") + 
				scale_color_manual(values = cols, drop = FALSE) +
				scale_y_continuous(breaks = seq(0,100,by=20), limits=c(0,100)) +
				theme_bw(base_size = tsize) +
				theme(legend.title = element_blank(),
					legend.position = "top",
					legend.justification = "right",
					legend.margin = margin(t = 0, b = 0),
					legend.box.margin = margin(t = -10, b = -10))
			
			xsc<-c(min(resf()$rdt),max(resf()$rdt))*tfac
			plow<--(min(xsc)/(xsc[2]-xsc[1]))
			pupp<-max(xsc)/(xsc[2]-xsc[1])
			
			if (plow>0.1) {
				pstop<-pstop + 
					annotate("text", x = 0, y = 100,
					label = dlab[2],
					hjust = 1.1, vjust = -0.4, 
					color = colbg[2])
			}
			if (pupp>0.1) {
				pstop<-pstop + 
					annotate("text", x = 0, y = 100,
					label = dlab[1],
					hjust = -0.1, vjust = -0.4, 
					color = colbg[1])
			}
						
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				pstop		
			} else {
				original_breaks <- pretty(exp(xsc), n = 6)     
				breaks_log <- log(original_breaks)                 
				pstop +  scale_x_continuous(breaks = log(original_breaks), 
					labels = original_breaks,
					limits = xsc)
			}
			
			
		})	%>% bindEvent(input$calc_btn)
		
		
		output$pbias<-renderPlot({
			
			input1<-shared_inputs()		
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			tsize<-16
			colbg<-c(rgb(0.90, 0.49, 0.13, alpha = 1),rgb(0.20, 0.60, 0.86, alpha = 1))
			dlab<-c("Control better","Experimental better")
			alphabg<-0.1
			
			if (direct!="lower") {
				colbg<-rev(colbg)
				dlab<-rev(dlab)
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				tfac<-100
				xnam<-"True risk difference (%)"
			} else {
				tfac<-1
				xnam<-"True risk ratio"
			}
			
			pbias<-resf()  %>%	
				mutate(rdt = rdt*tfac) %>%
				mutate(bias = bias*tfac) %>%
				ggplot(aes(x=rdt,y=bias, group = 1)) +
				annotate("rect",
						xmin = -Inf, xmax = 0,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[2],alpha = alphabg) +
				annotate("rect",
						xmin = 0, xmax = Inf,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[1],alpha = alphabg) +
				geom_line() +
				geom_hline(yintercept = 0, linetype = "dashed") +
				xlab(xnam) +
				ylab("Bias of risk difference (%)") + 
				theme_bw(base_size = tsize) 
			
			xsc<-c(min(resf()$rdt),max(resf()$rdt))*tfac
			plow<--(min(xsc)/(xsc[2]-xsc[1]))
			pupp<-max(xsc)/(xsc[2]-xsc[1])
	
			if (plow>0.1) {
				pbias<-pbias + 
					annotate("text", x = 0, y = max(resf()$bias),
					label = dlab[2],
					hjust = 1.1, vjust = -0.4, 
					color = colbg[2])
			}
			if (pupp>0.1) {
				pbias<-pbias + 
					annotate("text", x = 0, y = max(resf()$bias),
					label = dlab[1],
					hjust = -0.1, vjust = -0.4, 
					color = colbg[1])
			}
			
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				pbias		
			} else {
				original_breaks <- pretty(exp(xsc), n = 6)     
				breaks_log <- log(original_breaks)                 
				pbias +  scale_x_continuous(breaks = log(original_breaks), 
					labels = original_breaks,
					limits = xsc)
			}
			
			#browser()
			
		})	%>% bindEvent(input$calc_btn)
		
		
		output$pss<-renderPlot({
			
			input1<-shared_inputs()	
			
			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}
			
			tsize<-16
			colbg<-c(rgb(0.90, 0.49, 0.13, alpha = 1),rgb(0.20, 0.60, 0.86, alpha = 1))
			dlab<-c("Control better","Experimental better")
			alphabg<-0.1
			
			if (direct!="lower") {
				colbg<-rev(colbg)
				dlab<-rev(dlab)
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				tfac<-100
				xnam<-"True risk difference (%)"
				ebw<-0.2
			} else {
				tfac<-1
				xnam<-"True risk ratio"
				ebw<-0.2/100
			}
			
			pss<-resf()  %>%	
				mutate(rdt = rdt*tfac) %>%
				ggplot(aes(x=rdt,y=avn, group = 1)) +
				annotate("rect",
						xmin = -Inf, xmax = 0,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[2],alpha = alphabg) +
				annotate("rect",
						xmin = 0, xmax = Inf,
						ymin = -Inf, ymax = Inf,
						fill =   colbg[1],alpha = alphabg) +
				geom_line() +
				geom_errorbar(aes(ymin = minn, ymax = maxn), width = ebw, alpha=0.2) +
				scale_y_continuous(limits=c(0,max(resf()$maxn))) +
				xlab(xnam) +
				ylab("Sample size") + 
				theme_bw(base_size = tsize) 
				

			xsc<-c(min(resf()$rdt),max(resf()$rdt))*tfac
			plow<--(min(xsc)/(xsc[2]-xsc[1]))
			pupp<-max(xsc)/(xsc[2]-xsc[1])
	
			if (plow>0.1) {
				pss<-pss + 
					annotate("text", x = 0, y = 0,
					label = dlab[2],
					hjust = 1.1, vjust = 0.4, 
					color = colbg[2])
			}
			if (pupp>0.1) {
				pss<-pss + 
					annotate("text", x = 0, y = 0,
					label = dlab[1],
					hjust = -0.1, vjust = 0.4, 
					color = colbg[1])
			}
			
			if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
				pss		
			} else {
				original_breaks <- pretty(exp(xsc), n = 6)     
				breaks_log <- log(original_breaks)                 
				pss +  scale_x_continuous(breaks = log(original_breaks), 
					labels = original_breaks,
					limits = xsc)
			}
			
			#browser()
			
		})	%>% bindEvent(input$calc_btn)
		
		
	})	
}


#----------------------------------------------------------------------------------------------#
#Combined
#----------------------------------------------------------------------------------------------#

ui <- navbarPage(
  "Adaptive design",
  header = div(
    style = "padding:10px 20px;",
    p("This application allows you to explore operating characteristics 
      for an adaptive design with a binary endpoint and two groups.")
  ),
  tabPanel("Single scenario", ui1("app1")),
  tabPanel("Over alternatives", ui2("app2"))
)

server <- function(input, output, session) {
	
	
	settings_data <- server1("app1")
  
	#server1("app1")
	server2("app2", settings_data)
}


shinyApp(ui, server)
#shinyApp(ui1, server1)

