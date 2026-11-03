* Encoding: windows-1255.
* Kaplan-Meier by clinic. Requires only survival and standard R packages.
* Run the entire file against the active SPSS dataset.
* Unweighted, pooled analysis. Existing FILTER is honored.
* Does not change SPSS FILTER, WEIGHT or SPLIT settings.
* SPSS WEIGHT / SPLIT are not applied to this R analysis.

BEGIN PROGRAM R.
local({
# ================= USER SETTINGS =================
TIME_VAR <- "survtime"
STATUS_VAR <- "status"
GROUP_VAR <- "clinic"
ID_VAR <- "id"                  # NULL disables duplicate-ID check
EVENT_VALUE <- 1
CENSOR_VALUE <- 0
GROUP_ORDER <- NULL             # e.g. c("1", "2"); first is reference
GROUP_LABELS <- NULL            # e.g. c("Clinic 1", "Clinic 2")
TITLE <- "Kaplan-Meier Survival by Clinic"
X_LABEL <- "Time (Days)"
Y_LABEL <- "Survival Probability"
TIME_MAX <- NULL                # NULL = maximum observed follow-up
TIME_BREAKS <- NULL             # e.g. seq(0, 1000, 200)
SHOW_POLISHED_CAPTION <- TRUE   # combine enabled results into one caption
VERBOSE_OUTPUT <- FALSE         # TRUE restores detailed R tables
SHOW_CI <- TRUE
SHOW_SIGNIFICANCE_TEXT <- TRUE  # conditional interpretation of log-rank test
SIGNIFICANCE_ALPHA <- 0.05       # significant when p < this threshold
SHOW_LOGRANK <- TRUE
SHOW_COX_HR <- TRUE            # two groups only; unadjusted Cox model
CONF_LEVEL <- 0.95
OUTPUT_DIR <- ""               # empty = Viewer only; e.g. "C:/Results"
FILE_PREFIX <- "KM_clinic"
WIDTH <- 3                    # exported figure width, inches
HEIGHT <- 3                    # exported figure height, inches
PNG_DPI <- 300                  # 1200 x 1200 pixels at 4 x 4 inches
FONT_PT <- 6                   # base font size in points
CENSOR_SIZE <- 0.8             # censoring cross size
CENSOR_LWD <- 0.8
CURVE_LWD <- 0.8
ALL_LINES_SOLID <- TRUE        # FALSE = different line types by group
# =================================================

if (!requireNamespace("survival", quietly=TRUE)) stop("survival is not available.")
if (!is.numeric(CONF_LEVEL) || length(CONF_LEVEL)!=1 ||
    !is.finite(CONF_LEVEL) || CONF_LEVEL<=0 || CONF_LEVEL>=1)
  stop("CONF_LEVEL must be between 0 and 1.")
quiet_cat <- function(...) { if (VERBOSE_OUTPUT) cat(...) }
quiet_print <- function(...) { if (VERBOSE_OUTPUT) print(...) }
vars <- unique(c(TIME_VAR, STATUS_VAR, GROUP_VAR, ID_VAR))
raw <- spssdata.GetDataFromSPSS(variables=vars, factorMode="none",
                               missingValueToNA=TRUE, keepUserMissing=FALSE)
quiet_cat("\nAnalysis: unweighted; pooled across SPSS SPLIT groups. Active FILTER is honored.\n")
quiet_cat("Rows read:", nrow(raw), "\n")
if (!all(vars %in% names(raw))) stop("Check variable names (including case).")
d <- data.frame(time=raw[[TIME_VAR]], status=raw[[STATUS_VAR]],
                group=raw[[GROUP_VAR]])
keep <- complete.cases(d)
quiet_cat("Rows excluded for missing time/status/clinic:", sum(!keep), "\n")
if (!is.null(ID_VAR)) {
  ids <- raw[[ID_VAR]][keep]
  if (anyNA(ids) || anyDuplicated(ids))
    stop("Missing or duplicate IDs among analysis cases. Expected one row per person.")
}
d <- d[keep,,drop=FALSE]
if (!nrow(d)) stop("No complete analysis cases.")
if (!is.numeric(d$time) || any(!is.finite(d$time)) || any(d$time<0))
  stop("Survival times must be finite, numeric and nonnegative.")
if (EVENT_VALUE==CENSOR_VALUE || !all(d$status %in% c(EVENT_VALUE,CENSOR_VALUE)))
  stop("Unexpected status codes, or identical event/censor codes.")
d$event <- as.integer(d$status==EVENT_VALUE)
observed <- levels(factor(d$group))
if (!is.null(GROUP_ORDER) &&
    (anyDuplicated(GROUP_ORDER) || !setequal(as.character(GROUP_ORDER),observed)))
  stop("GROUP_ORDER must contain each observed clinic code exactly once.")
lev <- if (is.null(GROUP_ORDER)) observed else as.character(GROUP_ORDER)
d$group <- factor(as.character(d$group), levels=lev)
k <- nlevels(d$group)
if (k<2 || k>6) stop("This layout supports 2 to 6 observed clinics.")
labels <- if (is.null(GROUP_LABELS)) paste("Clinic",lev) else GROUP_LABELS
if (length(labels)!=k || anyNA(labels)) stop("Provide one label for each clinic.")
if (!sum(d$event)) stop("No events in analysis cases.")
xmax <- if (is.null(TIME_MAX)) max(d$time) else TIME_MAX
if (length(xmax)!=1 || !is.finite(xmax) || xmax<=0) stop("TIME_MAX must be positive.")
ticks <- if (is.null(TIME_BREAKS)) pretty(c(0,xmax),n=6) else TIME_BREAKS
if (!is.numeric(ticks) || any(!is.finite(ticks))) stop("Invalid TIME_BREAKS.")
ticks <- sort(unique(c(0,ticks[ticks>=0 & ticks<=xmax])))
if (length(ticks)<2) ticks <- c(0,xmax)

fit <- survival::survfit(survival::Surv(time,event)~group, data=d,
                         conf.int=CONF_LEVEL, conf.type="log-log")
quiet_print(fit)
quiet_cat("\nGroup order and event counts:\n")
quiet_print(data.frame(code=lev,label=labels,N=as.vector(table(d$group)),
                 events=as.vector(tapply(d$event,d$group,sum))))
# Number at risk immediately before each time point, including events at t.
# Ordinary right-censoring, no delayed entry; zero beyond observed follow-up.
risk <- t(vapply(lev,function(g)
  vapply(ticks,function(t) sum(d$group==g & d$time>=t),numeric(1)),
  numeric(length(ticks))))
rownames(risk) <- labels
colnames(risk) <- format(ticks,trim=TRUE)
quiet_cat("\nNumber at risk (t-):\n"); quiet_print(risk)

annotations <- character()
p <- NULL
cox_caption <- NULL
if (SHOW_LOGRANK || SHOW_SIGNIFICANCE_TEXT) {
  lr <- survival::survdiff(survival::Surv(time,event)~group,data=d,rho=0)
  p <- pchisq(lr$chisq,df=k-1,lower.tail=FALSE)
  ptext <- if (p<0.001) "p < 0.001" else sprintf("p = %.3f",p)
  if (SHOW_LOGRANK) annotations <- c(annotations,paste("Log-rank:",ptext))
  if (SHOW_SIGNIFICANCE_TEXT) {
    if (length(SIGNIFICANCE_ALPHA)!=1 || !is.numeric(SIGNIFICANCE_ALPHA) ||
        !is.finite(SIGNIFICANCE_ALPHA) || SIGNIFICANCE_ALPHA<=0 || SIGNIFICANCE_ALPHA>=1)
      stop("SIGNIFICANCE_ALPHA must be between 0 and 1.")
    interpretation <- if (!is.finite(p)) {
      "Log-rank significance could not be determined."
    } else if (p < SIGNIFICANCE_ALPHA) {
      "Survival differs significantly between clinics."
    } else {
      "No statistically significant survival difference detected."
    }
    annotations <- c(annotations, interpretation)
    quiet_cat("\n",interpretation," (alpha = ",SIGNIFICANCE_ALPHA,")\n",sep="")
  }
  quiet_cat("\nLog-rank test:\n"); quiet_print(lr)
}
if (SHOW_COX_HR) {
  if (k!=2) stop("SHOW_COX_HR requires exactly two clinics.")
  cw <- character()
  cx <- withCallingHandlers(
    survival::coxph(survival::Surv(time,event)~group,data=d,ties="efron"),
    warning=function(w) { cw <<- c(cw,conditionMessage(w)); invokeRestart("muffleWarning") })
  cs <- summary(cx,conf.int=CONF_LEVEL)
  quiet_cat("\nUnadjusted Cox model (second clinic versus first):\n"); quiet_print(cs)
  ci <- cs$conf.int[1,]
  if (length(cw) || any(!is.finite(ci)) || ci[3]<=0 || ci[4]<=0) {
    cat("Cox HR omitted from figure; model warnings or non-finite estimate.\n")
    if(length(cw)) cat(paste(cw,collapse="\n"),"\n")
  } else {
    cox_caption <- sprintf("Unadjusted Cox proportional hazards model: HR = %.2f, %.0f%% CI [%.2f-%.2f] (%s vs %s)",
      ci[1],100*CONF_LEVEL,ci[3],ci[4],labels[2],labels[1])
    annotations <- c(annotations,sprintf("Unadjusted HR (%s vs %s): %.2f; %.0f%% CI %.2f-%.2f",
      labels[2],labels[1],ci[1],100*CONF_LEVEL,ci[3],ci[4]))
    quiet_cat("\nProportional-hazards diagnostic:\n")
    ph <- tryCatch(survival::cox.zph(cx),error=function(e) conditionMessage(e))
    quiet_print(ph)
    if (inherits(ph,"cox.zph") && is.finite(ph$table["GLOBAL","p"]) &&
        ph$table["GLOBAL","p"] < 0.05)
      cat("Cox diagnostic: evidence against proportional hazards (global P < 0.05). Review the HR interpretation.\n")
    if (is.character(ph)) cat("Cox diagnostic unavailable:",ph,"\n")
  }
}
# Compose from calculated results, never by parsing printed R output.
if (SHOW_POLISHED_CAPTION) {
  group_text <- if(k==2) paste(labels,collapse=" and ") else
    paste0(paste(labels[-k],collapse=", "),", and ",labels[k])
  intro <- paste0("Kaplan-Meier survival curves comparing ",group_text,".")
  details <- character()
  if (SHOW_LOGRANK && !is.null(p) && is.finite(p)) {
    p_caption <- if(p<0.001) "P < 0.001" else
      paste0("P = ",formatC(p,format="f",digits=3,decimal.mark="."))
    details <- c(details,paste("Log-rank test:",p_caption))
  }
  if (!is.null(cox_caption)) details <- c(details,cox_caption)
  finding <- ""
  if (SHOW_SIGNIFICANCE_TEXT && !is.null(p) && is.finite(p)) {
    finding <- if(p<SIGNIFICANCE_ALPHA)
      "The difference in survival was statistically significant" else
      "No statistically significant difference in survival was detected"
    # State the threshold explicitly if the user changes its default.
    if(SIGNIFICANCE_ALPHA!=0.05)
      finding <- paste0(finding," at alpha = ",SIGNIFICANCE_ALPHA)
  }
  ending <- if(length(details)) paste0(" (",paste(details,collapse="; "),")") else ""
  caption <- if(nzchar(finding)) {
    paste0(intro," ",finding,ending,".")
  } else if(length(details)) {
    paste0(intro," ",paste(details,collapse="; "),".")
  } else intro
  # Wrapping uses figure size/font settings; reserve enough panel height.
  annotations <- strwrap(caption,width=max(25,floor((WIDTH-0.4)*72/(FONT_PT*0.55))))
  cat("\n",caption,"\n",sep="")
}
caption_height <- max(0.7,0.22*length(annotations)+0.2)
cols <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9")[seq_len(k)]
ltys <- if (ALL_LINES_SOLID) rep(1L,k) else seq_len(k)
# Convert a right-continuous step function to polygon coordinates.
step_xy <- function(x,y) {
  n <- length(x)
  if (n==1) return(list(x=x,y=y))
  list(x=c(x[1],rep(x[-1],each=2)),
       y=c(rep(y[-n],each=2),y[n]))
}
draw <- function() {
  old <- par(no.readonly=TRUE)
  on.exit({layout(1); par(old)},add=TRUE)
  layout(matrix(1:3,ncol=1),heights=c(4.3,caption_height,0.6+0.30*k))
  par(mar=c(3.3,5.3,3,1),mgp=c(2,0.6,0),tcl=-0.25,
      las=1,family="sans",cex=1,cex.main=1.15)
  plot(NA,xlim=c(0,xmax),ylim=c(0,1),xaxs="i",yaxs="i",
       xlab=X_LABEL,ylab="",axes=FALSE,main=TITLE)
  mtext(Y_LABEL,side=2,line=3.8,las=0)
  abline(h=seq(0,1,.2),v=ticks,col="#E8E8E8",lwd=.7)
  if (SHOW_CI) for (i in seq_len(k)) {
    f <- fit[i]
    # Include time zero and use KM steps, not linear interpolation.
    tt <- c(0,f$time); lo <- c(1,f$lower); hi <- c(1,f$upper)
    a <- step_xy(tt,lo); b <- step_xy(tt,hi)
    xx <- c(a$x,rev(b$x)); yy <- c(a$y,rev(b$y))
    good <- is.finite(xx)&is.finite(yy)
    polygon(xx[good],yy[good],col=adjustcolor(cols[i],alpha.f=.18),border=NA)
  }
  # Explicit base-graphics paths avoid reliance on lines.survfit dispatch
  # in the embedded SPSS session. Draw all curves AFTER all CI ribbons.
  for (i in seq_len(k)) {
    f <- fit[i]
    path <- step_xy(c(0, f$time), c(1, as.numeric(f$surv)))
    graphics::lines.default(path$x, path$y, col=cols[i],
                            lty=ltys[i], lwd=CURVE_LWD)
    censored <- which(f$n.censor > 0)
    if (length(censored))
      graphics::points.default(f$time[censored], f$surv[censored],
                               col=cols[i], pch=3, cex=CENSOR_SIZE, lwd=CENSOR_LWD)
  }
  axis(1,at=ticks); axis(2,at=seq(0,1,.2),labels=paste0(seq(0,100,20),"%"))
  box(col="#777777")
  legend("topright",legend=paste0(labels," (N=",as.vector(table(d$group)),")"),
         col=cols,lty=ltys,lwd=CURVE_LWD,bty="n",cex=.95)
  par(mar=if(SHOW_POLISHED_CAPTION) c(0.2,1,0.2,1) else c(0,5.3,0,1))
  plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
  if(length(annotations)) text(0,seq(.85,.20,length.out=length(annotations)),
                              annotations,adj=c(0,.5),cex=.95)
  par(mar=c(1,5.3,1.6,1),mgp=c(2,.6,0))
  plot(NA,xlim=c(0,xmax),ylim=c(.5,k+.5),xaxs="i",axes=FALSE,xlab="",ylab="")
  title(main="Number at Risk",adj=0,line=.35,cex.main=1)
  abline(v=ticks,col="#E8E8E8",lwd=.7)
  for(i in seq_len(k)) {
    y <- k-i+1
    abline(h=y,col="#F0F0F0")
    mtext(labels[i],side=2,at=y,las=1,line=1.5,col=cols[i],cex=.95)
    text(ticks,rep(y,length(ticks)),risk[i,],cex=.95,xpd=NA)
  }
}
# Render explicitly at full size, then submit the finished PNG to SPSS.
# This bypasses the small automatic Viewer graphics device.
if (!exists("spssRGraphics.Submit",mode="function"))
  stop("SPSS graphics submission function is not available.")
if (any(!is.finite(c(WIDTH,HEIGHT,PNG_DPI,FONT_PT))) ||
    any(c(WIDTH,HEIGHT,PNG_DPI,FONT_PT)<=0)) stop("Invalid figure dimensions/font settings.")
save_plot <- function(path,type) {
  if(type=="pdf") grDevices::pdf(path,width=WIDTH,height=HEIGHT,pointsize=FONT_PT)
  else grDevices::png(path,width=WIDTH,height=HEIGHT,units="in",res=PNG_DPI,
                      pointsize=FONT_PT,bg="white")
  on.exit(grDevices::dev.off(),add=TRUE)
  draw()
}
viewer_file <- tempfile(pattern="SPSS_KM_HiRes_",fileext=".png")
save_plot(viewer_file,"png")
spssRGraphics.Submit(viewer_file)
quiet_cat(sprintf("\nViewer image: %d x %d pixels; %.1f x %.1f inches; %d DPI; font %d pt.\n",
            round(WIDTH*PNG_DPI),round(HEIGHT*PNG_DPI),WIDTH,HEIGHT,PNG_DPI,FONT_PT))
# Leave the temporary image available for the duration of this R session.
if (nzchar(OUTPUT_DIR)) {
  if (!dir.exists(OUTPUT_DIR)) stop("OUTPUT_DIR does not exist; create it first.")
  if (basename(FILE_PREFIX)!=FILE_PREFIX || !nzchar(FILE_PREFIX)) stop("Invalid FILE_PREFIX.")
  pdf_file <- file.path(OUTPUT_DIR,paste0(FILE_PREFIX,".pdf"))
  png_file <- file.path(OUTPUT_DIR,paste0(FILE_PREFIX,".png"))
  if(any(file.exists(c(pdf_file,png_file)))) stop("Output file already exists. Change FILE_PREFIX.")
  save_plot(pdf_file,"pdf")
  if (!file.copy(viewer_file,png_file,overwrite=FALSE)) stop("Could not save PNG.")
  quiet_cat("\nSaved:\n",pdf_file,"\n",png_file,"\n")
}
quiet_cat("\nAnalysis complete. No additional packages installed.\n")
})
END PROGRAM.







