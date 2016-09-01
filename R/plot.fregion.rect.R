#' Visualize rectangular confidence regions
#'
#' @param rect A 'fregion.rect' object, the output from 'fregion.rect()' funciton.
#' @param npc Number of PCs to plot
#' @param grid.size Since this function makes plots, fd objects are evaluated and gird.size determines how find the grid should be.
#' @examples
#' plot(rect)
#' @export

plot.fregion.rect <- function(rect,npc=NULL,returntype=NULL,main=NULL,xlab=NULL,ylab=NULL,grid.size=200){

  if (is.null(rect$harmonics)) {datatype="vector"} else {datatype="fd"}
  if (datatype=="fd") {evalgrid <- make.grid(grid.size, rangevals=rect$harmonics$basis$rangeval)}

  if (is.null(npc)) npc <- dim(rect$table)[1]
  range <- c(1:npc)

  d <- rect$table[1:npc,]
  x <- as.integer(d$PC)
  y <- d$abs.score
  ub <- y + d$score.bound
  lb <- y - d$score.bound
  eps <- 0.05                # for plot
  y2 <- d$cum.var.explained
  y3 <- d$var.explained

  y0 <- abs(d$coefficient)
  y0.u <- y0 + d$MoE
  y0.l <- y0 - d$MoE

  if (substr(rect$type, nchar(rect$type), nchar(rect$type)) == "s") {prefix <- "t-"} else {prefix <- "z-"}

  plot.intervals <- function(){
    par(mar=c(5.1,4.1,4.1,4.1))

    #plot(x, y2,type="n"xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1))
    #polygon(c(min(x),x,max(x)),c(0,y2,0),col='#e0ecf4',lty=0)

    #plot(x, y2,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
    #polygon(c(min(x),x,max(x)),c(min(y),y2,min(y)),col='#e0ecf4',lty=0)
    #axis(4,labels=TRUE)
    #mtext("Variance Explained (Shaded)",side=4,line=3)

    #plot(x, y0,type="b",xaxt="n",yaxt="n",xlab="",ylab="")
    #polygon(c(x,rev(x)),c(y0.u,rev(y0.l)),col='#e0ecf4',lty=0)
    #axis(4,labels=TRUE)
    #mtext("Marginal Margine of Error (Shaded)",side=4,line=3)
    #par(new=TRUE)

    ### Plot As z-Score ###
    #plot(x,y,pch="*",cex=2,ylim=c(min(lb),max(ub)),
    #     main=paste0("Marginal Intervals (as ", prefix,"score) for Each PC"),
    #     xlab="PC",ylab="Absolute Score with Marginal Bounds", xaxt='n')
    #axis(1,at=x)
    #abline(h=0,col="grey")
    #segments(x,lb,x,ub,col="red",lwd=2)
    #segments(x-eps,ub,x+eps,ub,col="red",lwd=2)
    #segments(x-eps,lb,x+eps,lb,col="red",lwd=2)
    #points(x,y,pch="*",cex=2)

    ## Plot as Coefficient ###
    plotmargin <- (max(y0.u) - min(y0.l))*0.1
    if (is.null(main)) main.intervals <- "Marginal Intervals for Each PC" else main.intervals <- main

    plot(x,y0,pch="*",cex=2,ylim=c(min(y0.l)-plotmargin,max(y0.u)+plotmargin),
         main=main.intervals,
         xlab="PC",ylab="Absolute Coefficients with Marginal Intervals", xaxt='n')
    axis(1,at=x)
    abline(h=0,col="grey")
    segments(x,y0.l,x,y0.u,col="red",lwd=2)
    segments(x-eps,y0.u,x+eps,y0.u,col="red",lwd=2)
    segments(x-eps,y0.l,x+eps,y0.l,col="red",lwd=2)
    points(x,y0,pch="*",cex=1.5)
    text(x=x,y=y0.u,labels=round(y2*100,2),col="blue",cex=0.9,pos=3)
    text(x=x,y=y0.l,labels=round(y3*100,2),col="darkgreen",cex=0.9,pos=1)
    #mtext(round(y2,3),side=1,line=2,at=x,col="blue",cex=0.5)
    #mtext(round(y3,3),side=1,line=3,at=x,col="green",cex=0.5)
    #legend(x="topright",legend=c("Cumulative Variance Explained", "Variance Explained"), text.col=c("blue","green"),cex=c(0.8,0.8))
    legend(x="topright",legend=c("Cumulative Variance Explained (%)","Variance Explained (%)"),text.col=c("blue","darkgreen"),cex=0.9)
  }

#   plot1 <- ggplot(rect$table[1:npc,],aes(y=abs.score,x=PC)) +
#       geom_errorbar(aes(ymin=abs.score - score.bound,ymax=abs.score + score.bound),width=0.5,colour='red') +
#       theme_bw() + geom_hline(yintercept=0,col="blue",width=5) +  geom_point(size=10,shape="*") +
#       ggtitle("Marginal Intervals for Each PC")

#   plot2 <- ggplot(rect$table[1:npc,],aes(y=var.explained,x=PC)) +
#     #geom_ribbon(aes(x=PC, ymin=0, ymax=var.explained, colour="red"),alpha=0.3) + theme_bw()
#     geom_area(aes(x=PC, y=var.explained, colour="red"),alpha=0.3) + theme_bw()


  plot.band <- function(pc){
    if (is.null(main)) main.band <- paste0("Confidence Band along PC ", pc) else main.band <- main
    if (is.null(xlab)) xlab.band <- "" else xlab.band=xlab
    if (is.null(ylab)) ylab.band <- "" else ylab.band=ylab

    if (datatype=="fd") {eigenfunction <- as.vector(eval.fd(evalgrid,rect$harmonics[pc]))} else
                        {eigenfunction <- rect$eigenfunctions[,pc]}
    projection <- rect$table$coefficient[pc] * eigenfunction
    marginalprojection.u <- projection + rect$table$MoE[pc] * eigenfunction
    marginalprojection.l <- projection - rect$table$MoE[pc] * eigenfunction
    ylim <- c(min(marginalprojection.u,marginalprojection.l,0),max(marginalprojection.u,marginalprojection.l,0))
    Mat <- cbind(projection,marginalprojection.u,marginalprojection.l)

    if (datatype=="fd") xvalues <- evalgrid else xvalues <- rownames(rect$eigenfunctions)
    matplot(xvalues, Mat, ylim=ylim,
            type="l",xlab=xlab.band,ylab=ylab.band,main=main.band,
            lwd=c(2,2,2), col=c(1,2,2), lty=c(1,3,3))
    abline(h=0,col="grey",lwd=0.5)
    legend(x="topleft",lty=c(1,3),lwd=c(2,2),col=c(1,2),
           legend=c(paste0("Actual Projection on PC ", pc), paste0("Confidence Band along PC ",pc)))
  }

  #grid.newpage()
  #grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size <- "last"))

  if (is.null(returntype)) {
    plot.intervals()
    repeat{
      pc <- readline("Type in the PC to show as band (0 for intervals, ENTER for exit) : ")
      if (pc == "") break else {
        pc <- as.integer(pc)
        if (pc==0) plot.intervals() else plot.band(as.integer(pc))
      }
    }
  } else {
    if (returntype==0) return(plot.intervals()) else
      return(plot.band(as.integer(returntype)))
  }

}


