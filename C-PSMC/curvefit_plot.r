### rscript <filename>.r
args = commandArgs(trailingOnly=TRUE)



pairs <- args[2:length(args)]


colors <- rep('black',length(pairs))
line_type <- rep('solid',length(pairs))






ymin <- 0
ymax <- 2
xmin <- 0
xmax <- as.numeric(args[1])
method <- 'Average slope difference' ## 'Average slope difference' or 'Pearson correlation'


for (pair in pairs){
	#gs <- paste(pair,'.txt',sep='')
	gs <- pair
	gs <- read.csv(gs)
	names(gs) <- c('g','fit','atomic_time_units','timepoints')
	gs <- gs[,c('g','fit')] ## select two columns (of the four) to plot
	#gscalar <- target_g_factors[which(pairs %in% pair)]
	gscalar <- 1
	gs$g <- gs$g * gscalar
	filename <- paste(pair,'.png',sep='')
	png(filename,height=6,width=12,units='in',res=600)
	#plot(gs,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab='Relative generation time',ylab=method,type='n', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
	plot(gs,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab='Generation time',ylab=method,type='n',xaxt='n')	## type='n' makes the individual points invisible, so they can be traced over with a line
	axis(1, at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5), labels=c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0')) ## side: (1=bottom, 2=left, 3=top, 4=right)
	#axis(1, at=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), labels=c('0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')) ## side: (1=bottom, 2=left, 3=top, 4=right)
	lines(gs, col=colors[which(pairs %in% pair)], type='l', lwd=3, lty=line_type[which(pairs %in% pair)]) ## lwd == line weight
	print(pair)
}






