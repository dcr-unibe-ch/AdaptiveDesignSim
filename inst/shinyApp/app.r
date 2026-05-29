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
library(AdaptiveDesignSim)

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
				#numericInput(ns("n"), "Sample size",
				#	min = 0, max = Inf, value = 1000, step = 1
				#),
				h5(strong("Sample Size")),
				fluidRow(
					column(width = 6,
						numericInput(inputId = ns("n0"),
							label = span("Control", style = "font-weight: normal;"),
							min = 0, max = Inf, value = 500, step = 1 )
					),
					column(width = 6,
						numericInput(inputId = ns("n1"),
							label = span("Intervention", style = "font-weight: normal;"),
							min = 0, max = Inf, value = 500, step = 1)
					)
				),
				h3("Interim analysis"),
				numericInput(ns("nia"), "Number of interim analyses",
					min = 0, max = Inf, value = 0, step = 1
				),
				uiOutput(ns("dynamic_tp_text")),
				uiOutput(ns("dynamic_tp")),
				uiOutput(ns("dynamic_zsign")),
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
		output$dynamic_tp_text <- renderUI({

			req(input$nia)

			if (input$nia==0) {
				NULL
			} else {
				tags$p(h5("as fraction of the total sample size"))
			}
		})

		output$dynamic_tp <- renderUI({

			req(input$nia)

			if (input$nia==0) {
				NULL
			} else {

				nia <- input$nia
				sp<-1/(nia+1)
				sp<-seq(sp,1,by=sp)
				sp<-sp[1:(length(sp)-1)]

				noUiSliderInput(
					inputId = session$ns("timepoint"),
					label = "Time point of interim analysis",
					min = 0,
					max = 1,
					value = sp,
					step = 0.01,
					margin = 0,
					connect = FALSE,
					format = list(decimals = 2),
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
							width = "100%"
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
							text(x = xlim[1]+adx,y=mean(ylim),labels=pzs[1],
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
							text(x = xlim[2]-adx,y=mean(ylim),labels=pzs[1],
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

		#n01<-reactive({
		#	n0<-n1<-input$n %/% 2
		#	left <- input$n - (n0+n1)
		#	left <- sample(c(0,left),replace=FALSE)
		#	n0<-n0 + left[1]
		#	n1<-n1 + left[2]
		#	stopifnot(n0+n1==input$n)
		#	return(c(n0,n1))
		#})

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

		observe({
		#req(input$timepoint)
		#print(input$timepoint)
		#req(input$nia)
		#print(input$nia)
		#req(input$zsign)
		#print(input$nia)
		#req(input$n0)
		#print(input$n0)
		#print(input$n1)
		#print(input$zval_1)
		})



	#simulation
	#----------

		#res<-reactive({
		res<-eventReactive(input$calc_btn, {

		  #browser()
			req(input$nia)
			req(input$n0)
			req(input$n1)

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
				tia<-NULL
			} else {
			  req(input$timepoint)
			  tia<-input$timepoint
			}

			for (r in 1:length(rdl)) {

				rdt<-rdl[r]

				if (is.null(input$effectm) || input$effectm=="Risk difference") {
					p1<-input$p0/100 + rdt
				} else {
					p1<-input$p0/100 * exp(rdt)
				}

				ri <- lapply(1:input$reps, function(x)
					ADsimfun(n01 = c(input$n0,input$n1), p01 = c(input$p0/100, p1),
						nia = input$nia, tia = tia,
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

			#z matrix
			zcut<-matrix(NA, nrow = input$nia, ncol = 2)
			if (input$nia>0) {
			  for (i in (1:input$nia)) {
					id <- paste0("zval_", i)
					zcut[i,]<-input[[id]]
			  }
			  if (input$zsign=="futility") {
			    zcut[,2]<-NA
			  }
			  if (input$zsign=="efficacy") {
			    zcut[,1]<-NA
			  }
			}

			resf<-res()

			for (ri in 1:length(res())) {

				resi<-res()[[ri]]

				#resf[[ri]]<-evalfun(simdat = resi, nia = input$nia,
				#	zsc = zsc, zsign = input$zsign, n = input$n,
				#	direct = direct)

				resf[[ri]]<-ADevalfun(simdat = resi, zcutoff = zcut)

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

			#browser()

			#erdf<-mean(resf()[[length(resf())]]$obs_pe)*tfac
			#ns<-nrow(resf()[[length(resf())]])
			#if (direct=="lower") {
			#	nsig<-sum(resf()[[length(resf())]]$obs_uci<0)
			#} else {
			#	nsig<-sum(resf()[[length(resf())]]$obs_lci>0)
			#}

			nsim<-nrow(resf()[[length(resf())]][["rawdata"]])
			erdf<-resf()[[length(resf())]][["opchar"]]["av_pe"]*tfac
			nsig<-resf()[[length(resf())]][["opchar"]]["psig"]*nsim

			ids<-(1:nsim) / (nsim + 1)*100
			psig<-mean(c(ids[nsig],ids[nsig+1]))
			if (nsig==nsim) {
				psig<-100
			}

			#browser()

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
			if ((direct=="lower" & nsig/nsim>0.9) | (direct=="higher" & nsig/nsim<0.1)) {
			  thjust<--0.1
			  txlabpos<-(-Inf)
			}
			if ((direct=="lower" & nsig==0) | (direct=="higher" & nsig==nsim)) {
			  tylabpos<-0
				bgpos<-0
				colbg[2]<-NA
			}
			if ((direct=="lower" & nsig==nsim) | (direct=="higher" & nsig==0)) {
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


			fpl<-resf()[[length(resf())]][["rawdata"]] %>%
				arrange(obs_pe) %>%
				mutate(id = ids) %>%
				mutate(across(starts_with("obs_"), ~ .x*tfac)) %>%
				mutate(rdfs = mean(obs_pe)) %>%
				ggplot(aes(x = obs_pe, y = id, group = anystop,color = anystop)) +
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
				geom_errorbar(aes(xmin = obs_lci, xmax = obs_uci), orientation =  "y", width = 0.2, show.legend = c(color = TRUE, linetype = FALSE)) +
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
				xlims <- c(min(resf()[[length(resf())]][["rawdata"]][,"obs_lci"]),
					max(resf()[[length(resf())]][["rawdata"]][,"obs_uci"]))
				original_breaks <- pretty(exp(xlims), n = 6)
				breaks_log <- log(original_breaks)
				fpl +  scale_x_continuous(breaks = log(original_breaks), labels = original_breaks)
			}

			#browser()
		}) %>% bindEvent(input$calc_btn)



	#table
	#----------

		output$tab <- renderDT({

			#req(input$zsign)
			req(input$nia)

			resi<-resf()[[length(resf())]][["opchar"]]


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
				rdtf<-paste0(ff(tfac*resi["true_pe"], dig=dig),addp)
				est<-ff_ci(resi[c("av_pe","av_lci","av_uci")]*tfac,
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

				rdtf<-paste0(ff(exp(resi["true_pe"]), dig=dig),addp)
				est<-ff_ci(exp(resi[c("av_pe","av_lci","av_uci")]),
				           dig=dig,fs="pe (lci to uci)")
			}

			nas<-c("Stopped","Stopped for futility","Stopped for efficacy")

			pstopf<-"NA"
			pstopfutf<-"NA"
			pstopefff<-"NA"

			if (input$nia>0) {

			  req(input$zsign)

			  pstopf <- paste0(ff(resi["pstop"]*100, dig=1),"%")
			  pstopfutf <- paste0(ff(resi["pstop_fut"]*100, dig=1),"%")
			  pstopefff <- paste0(ff(resi["pstop_eff"]*100, dig=1),"%")

				if (input$nia>1) {
					for (i in (1:input$nia)) {
					  sep<-ifelse(i==1," (",", ")
						pstopf<-paste0(pstopf,sep,ff(resi[paste0("pstop_ia",i)]*100, dig=1),"%")
						pstopfutf<-paste0(pstopfutf,sep,ff(resi[paste0("pstop_fut_ia",i)]*100, dig=1),"%")
						pstopefff<-paste0(pstopefff,sep,ff(resi[paste0("pstop_eff_ia",i)]*100, dig=1),"%")
					}
				  pstopf<-paste0(pstopf,")")
				  pstopfutf<-paste0(pstopfutf,")")
				  pstopefff<-paste0(pstopefff,")")
				  nas<-paste0(nas," (at each interim analysis)")
				}

			  if (input$zsign=="futility") {
			    pstopefff<-"NA"
			  }
			  if (input$zsign=="efficacy") {
			    pstopfutf<-"NA"
			  }

			}

			avnfor <- paste0(round(resi["avn"])," (",round(resi["minn"]),", ",round(resi["maxn"]),")")

      d1<-c(paste0(ff(resi["psig"]*100, dig=1),"%"),
            pstopf, pstopfutf, pstopefff,
            avnfor,
            rdtf, est,
            paste0(ff(tfac*resi["bias"], dig=dig),addp),
            paste0(ff(tfac*resi["sd"], dig=dig),addp),
            ff(tfac^2*resi["mse"], dig=dig^2),
            paste0(ff(resi["coverage"], dig=1),"%")
          )

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
				d1)

			colnames(tb)<-c(" "," ")

			tb %>%
				datatable(
				options = list(searching = FALSE,ordering = FALSE,
					paging = FALSE, info = FALSE),
					rownames = FALSE)

		}) %>% bindEvent(input$calc_btn)


#share inputs to 2nd tab
#-----------------

    return(reactive({

		req(input$nia)

		zcut<-matrix(NA, nrow = input$nia, ncol = 2)
		if (input$nia>0) {
		  for (i in (1:input$nia)) {
		    id <- paste0("zval_", i)
		    req(input[[id]])
		    zcut[i,]<-input[[id]]
		  }
		  if (input$zsign=="futility") {
		    zcut[,2]<-NA
		  }
		  if (input$zsign=="efficacy") {
		    zcut[,1]<-NA
		  }
		}

		#browser()

		list(
			p0 = input$p0,
			rd = input$rd,
			n01 = c(input$n0,input$n1),
			nia = input$nia,
			timepoint =input$timepoint,
			zsign = input$zsign,
			finalz = input$finalz,
			reps = input$reps,
			zcut = zcut,
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

		output$ass1 <- renderUI({
		  input1<-shared_inputs()
		  tagList(
				#hr(),
				p(strong("Control proportion: "), input1$p0, "%"),
				p(strong("Direction:"), input1$direction),
				p(strong("Effect measure: "), input1$effectm),
				p(strong("Sample size: "), paste0(input1$n01,collapse=", "))
		  )
		})

		output$ia1 <- renderUI({

			input1<-shared_inputs()

			req(input1$zsign)
			
			zcol<-apply(input1$zcut,2,function(x) paste(x,collapse="/"))
			
			if (input1$zsign=="both") {
				txt1<-paste0("at Z < ",zcol[1]," and Z >",zcol[2])
			} else {
				if (input1$zsign=="futility") {
					txt1<-paste0("at Z < ",zcol[1])
				} else {
					txt1<-paste0("at Z > ",zcol[2])
				}
			}

			tagList(
				p(strong("Number: "), input1$nia),
				p(strong("Time point: "),
					paste(input1$timepoint,collapse=", ")),
				p(strong("Stopping for: "),input1$zsign),
				p(txt1)
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

			#print(input1$p0)
			#print(input1$rd)
			#print(input1$nia)
			#print(input1$n01)
			#print(input1$timepoint)
			#print(input1$zcut)
			#print(input1$effectm)
			#print(input1$direction)
			#print(input$rdrange)
			#print(input$calc_btn)
		})


		#simulate over range of alternatives
		#-------
		#res<-reactive({
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
				tia<-NULL
			} else {
				req(input1$timepoint)
				tia<-input1$timepoint
				stopifnot((length(tia)==input1$nia))
			}

			#browser()

			for (r in 1:length(rdl)) {

				rdt<-rdl[r]

				if (is.null(input1$effectm) || input1$effectm=="Risk difference") {
					p1<-input1$p0/100 + rdt
				} else {
					p1<-input1$p0/100 * exp(rdt)
				}

				ri <- lapply(1:input$reps, function(x)
					ADsimfun(n01 = input1$n01, p01 = c(input1$p0/100, p1),
						nia = input1$nia, tia = tia,
						alpha = alpha, effm = effm, direct = direct))

				ri <- do.call(rbind, ri)

				res[[r]]<-data.frame(p0t=input1$p0/100,p1t=p1,rdt=rdt,ri)
			}

			return(res)

		})


		#observe({
		#	req(input$rdrange)
		#	print(head(res()[[1]]))
		#})

	#evaluate at threshholds
	#-------------------

		resf<-eventReactive(input$calc_btn, {

			#browser()

			input1<-shared_inputs()

			#req(input1$zsign)
			req(input1$nia)
			req(input1$n01)
			#req(input1$zsc)

			if (is.null(input1$direction) || input1$direction=="Lower is better") {
				direct<-"lower"
			} else {
				direct<-"higher"
			}

			outi<-res()

			for (ri in 1:length(res())) {

				resi<-res()[[ri]]

				resii<-ADevalfun(simdat = resi, zcutoff = input1$zcut)

				outi[[ri]]<-resii[["opchar"]]

			}
			out<-data.frame(do.call(rbind,outi))

			return(out)

		})

		#observe({
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

			#browser()

			ppwr<-resf()  %>%
				ggplot(aes(x=av_pe,y=psig*100, group = 1)) +
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

			xsc<-c(min(resf()[,"true_pe"]),max(resf()[,"true_pe"]))*tfac
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
				rdt=rep(resf()[,"true_pe"],3),
				pstop=c(resf()[,"pstop"],resf()[,"pstop_fut"],resf()[,"pstop_eff"])*100,
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

			xsc<-c(min(resf()[,"true_pe"]),max(resf()[,"true_pe"]))*tfac
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

		#browser()

			pbias<-resf()  %>%
				mutate(rdt = true_pe*tfac) %>%
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

			xsc<-c(min(resf()[,"true_pe"]),max(resf()[,"true_pe"]))*tfac
			plow<--(min(xsc)/(xsc[2]-xsc[1]))
			pupp<-max(xsc)/(xsc[2]-xsc[1])

			if (plow>0.1) {
				pbias<-pbias +
					annotate("text", x = 0, y = max(resf()[,"bias"]),
					label = dlab[2],
					hjust = 1.1, vjust = -0.4,
					color = colbg[2])
			}
			if (pupp>0.1) {
				pbias<-pbias +
					annotate("text", x = 0, y = max(resf()[,"bias"]),
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
				mutate(rdt = true_pe*tfac) %>%
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


			xsc<-c(min(resf()[,"true_pe"]),max(resf()[,"true_pe"]))*tfac
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

