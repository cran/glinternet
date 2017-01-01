plot.glinternet.cv = function(x, ...){
    #plot cv curve
    bestIndex = which.min(x$cvErr)
        barwidth = 0.25
        xi = 1:length(x$lambda)
        y = x$cvErr
        delta = x$cvErrStd
        plot(xi, y, type="n", xlab="Lambda index", ylab="CV error", xaxt="n", ylim=c(min(y-delta), max(y+delta)))
        segments(xi-barwidth, y+delta, xi+barwidth, y+delta,col="grey")
        segments(xi-barwidth, y-delta, xi+barwidth, y-delta,col="grey")
        segments(xi, y+delta,xi,y-delta,col="grey")
        lines(xi,y,lwd=2)

        abline(v=bestIndex, lty=3)
        axis(1, at=seq(2, length(xi), 2), las=1)
}
