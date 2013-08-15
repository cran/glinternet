plot.glinternet.cv = function(x, ...){

  #plot cv curve
  bestIndex = which.min(x$cvErr)
  if ("ggplot2" %in% rownames(installed.packages())){
    require(ggplot2)
    cvErr=x$cvErr
    cvErrStd=x$cvErrStd
    lambdaIdx=1:length(x$lambda)
    plt = ggplot(data.frame(cvErr, cvErrStd, lambdaIdx), aes(x=lambdaIdx, y=cvErr)) + geom_line(size=1.5) + geom_vline(xintercept=bestIndex, linetype="dashed") + scale_x_discrete("Lambda index", breaks=seq(2,length(x$lambda),2)) + scale_y_continuous("CV error") + theme_bw() + geom_ribbon(aes(ymin=cvErr-cvErrStd, ymax=cvErr+cvErrStd), fill="black", alpha=0.5)
    print(plt)
  }
  else {
    barwidth = 0.25
    x = 1:length(x$lambda)
    y = x$cvErr
    delta = x$cvErrStd
    plot(x, x$cvErr, type="l", lwd=1.5, xlab="Lambda index", ylab="CV error", xaxt="n", ylim=c(min(y-delta), max(y+delta)))
    segments(x-barwidth, y+delta, x+barwidth, y+delta)
    segments(x-barwidth, y-delta, x+barwidth, y-delta)
    abline(v=bestIndex, lty=3)
    axis(1, at=seq(2, length(x), 2), las=1)
  }
}
