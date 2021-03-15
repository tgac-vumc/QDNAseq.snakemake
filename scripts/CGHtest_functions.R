Perform_CGHtest <- function(data1,data2,name1,name2, output, teststat = "Chi-square",iterations = 1000, af = 0.1){
        combined <- combine(data1,data2)
        datainfo <- data.frame(chromosomes(combined), bpstart(combined), bpend(combined), nclone(combined), avedist(combined))
        data_combined <- regions(combined)
        for (i in 1:ncol(combined)) {
            regions(combined)[which(regions(combined)[,i] == -2),i] <- -1
            regions(combined)[which(regions(combined)[,i] == 2),i] <- 1
        }
        samples1 <- ncol(data1)
        samples2 <- ncol(data2)
                                        # p values and fdr's for copy number losses
        gfr_loss <- groupfreq(data_combined,sepfile="no",group=c(samples1,samples2),groupnames=c(name1, name2),af=af, lgonly=-1)
        colnames(gfr_loss) <- c(paste(name1,"loss"), paste(name1,"no loss"), paste(name2,"loss"), paste(name1,"no loss"))
        pvs_loss <- pvalstest(data_combined,datainfo,teststat = teststat ,group=c(samples1,samples2 ),groupnames=c(name1, name2),lgonly=-1, af=af, niter=iterations) # ncpus=ncores
        fdrs_loss <- fdrperm(pvs_loss,mtdirection="stepup")

                                        # p values and fdr's for copy number gains
        gfr_gain <- groupfreq(data_combined,sepfile="no",group=c(samples1,samples2),groupnames=c(name1, name2),af=af, lgonly=1)
        colnames(gfr_gain) <- c(paste(name1,"gain"), paste(name1,"no gain"), paste(name2,"gain"), paste(name1,"no gain"))
        pvs_gain <- pvalstest(data_combined,datainfo,teststat = teststat ,group=c(samples1,samples2 ),groupnames=c(name1, name2),lgonly=1, af=af, niter=iterations) # ncpus=ncores
        fdrs_gain <- fdrperm(pvs_gain,mtdirection="stepup")
        save(fdrs_loss, fdrs_gain,datainfo,pvs_loss,pvs_gain,gfr_gain,gfr_loss ,file=output)
 }

freqPlotCompare <- function(x, group1, group2, cohort1_ID, cohort2_ID){
    ###############
    # input parameters
    #y.limit <- c(-0.75, 0.75)
    y.limit <- c(-0.3, 0.55)
    misscol=NA
    build='GRCh37'
    
    ###############
    # Set colors
    
    group1.gaincol='darkgreen' # col2rgb('darkgreen') 0 100 0
    group1.gainFILLcol=rgb(0/255, 100/255, 0/255, alpha=0.4)
    
    #EF.losscol='red4'  # col2rgb('red4') 139 0 0
    group1.losscol='red1'  # col2rgb('red1') 255 0 0
    group1.lossFILLcol=rgb(139/255, 0/255, 0/255, alpha=0.4)
    
    #col2rgb('palegreen1') 154, 255, 154
    group2.gaincol=rgb(154/255, 255/255, 154/255, alpha=0.5)
    #col2rgb('firebrick1') 255, 48, 48
    group2.losscol=rgb(255/255, 48/255, 48/255, alpha=0.5)
    
    ###############
    # coordinate settings
    
    chrom <- chromosomes(x)
    pos <- bpstart(x)
    pos2 <- bpend(x)
    uni.chrom <- unique(chrom)
    chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
    chrom.ends <- integer()
    cumul <- 0
    for (j in uni.chrom) {
        pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]           # each bpstart that is not on chr1, gets the length of the bpstart + cum sum of all previous chromos, but that is not what is says..
        pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
        cumul <- cumul + chrom.lengths[as.character(j)]
        chrom.ends <- c(chrom.ends, cumul)
    }
    names(chrom.ends) <- names(chrom.lengths)
    
    ###############
    # group frequencies
    calls1 <- regions(group1)
    calls2 <- regions(group2)
    group1.loss.freq <- rowMeans(calls1 < 0)
    group1.gain.freq <- rowMeans(calls1 > 0)
    group2.loss.freq <- rowMeans(calls2 < 0)
    group2.gain.freq <- rowMeans(calls2 > 0)

    ###############
	# new vector for vertical line segments
	# not all regions are adjacent; make a dataframe with adjacent regions
	adjacent <- c()
	for (i in 1:length(pos)){
        z <- pos[i+1]-1==pos2[i]
        adjacent <- c(adjacent, z)
	}
	not.adjacent <- which(adjacent==FALSE)
	row.entries <- 1+ not.adjacent
	chrom.entries <- chrom[not.adjacent]
	start.entries <- 1+ pos2[not.adjacent]
	end.entries <- pos[row.entries] -1
	input.rows <- cbind(chrom.entries, start.entries, end.entries, 0, 0, 0, 0)
	allpos <- cbind(chrom,pos, pos2, group1.gain.freq, group1.loss.freq, group2.gain.freq, group2.loss.freq)
	dd <- as.data.frame(rbind(allpos, input.rows))
	dd <- dd[ order(dd[,1], dd[,2]), ]
	new.pos <- as.data.frame(dd)
	rownames(new.pos) <- c(1:nrow(new.pos))
	# New vector for x and y positions to draw vertical lines
	pos3 <- c(new.pos$pos, new.pos$pos2[length(new.pos$pos2)])	
	group1.gain2 <- c(0, new.pos$group1.gain.freq)
	group1.gain3 <- c(new.pos$group1.gain.freq, 0)
	group1.loss2 <- c(0, new.pos$group1.loss.freq)
	group1.loss3 <- c(new.pos$group1.loss.freq, 0)
    ###############
    # the plot

    title <- paste("Distribution of CNAs: ", cohort1_ID, " vs ", cohort2_ID, sep="")
    lab_y <- 'percentage gains and losses'
    plot(NA, xlim=c(0, max(pos2)), ylim=y.limit, type='n', xlab='chromosomes', ylab='', main=title ,xaxs='i', xaxt='n', yaxs='i', yaxt='n') #main=main ??
    if (!is.na(misscol)) {
        rect(0, -1, max(pos2), 1, col=misscol, border=NA)
        rect(pos, -1, pos2, 1, col='white', border=NA)
    }

    segments(pos, group1.gain.freq, pos2, group1.gain.freq, lty=1, lwd=2.5, col=group1.gaincol)         # gains: horizontal segments
    segments(pos3, group1.gain2, pos3, group1.gain3, lty=1, lwd=2.5, col=group1.gaincol)                        # gains: vertical (joining) segments

    segments(pos, -group1.loss.freq, pos2, -group1.loss.freq, lty=1, lwd=2.5, col=group1.losscol)       # losses: horizontal segments
    segments(pos3, -group1.loss2, pos3, -group1.loss3, lty=1, lwd=2.5, col=group1.losscol)                      # losses: vertical (joining) segments

    rect(pos, 0, pos2, group2.gain.freq, col=group2.gaincol, border=NA)
    rect(pos, 0, pos2, -group2.loss.freq, col=group2.losscol, border=NA)

    box()
    abline(h=0)

    if (length(chrom.ends) > 1){
        for (j in names(chrom.ends)[-length(chrom.ends)]){
            abline(v=chrom.ends[j], lty='dashed')
            ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
            axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
            #  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
            #  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
            axis(side=2, at=c(-0.25, 0, 0.25, 0.5), labels=c(' 25 %', '0 %', '25 %', '50 %'), las=1)
            # mtext('gains', side=2, line=3, at=0.4)
            # mtext('losses', side=2, line=3, at=-0.4)
        }
    }
    ### number of data points
    mtext(paste(nrow(x), 'regions'), side=3, line=0, adj=0)
    mtext(paste(ncol(x), 'samples'), side=3, line=0, adj=1)
    mtext(lab_y, side=2, line=3, at=0)
    # LEGEND
    legend('topleft',
           legend = c(paste0(cohort1_ID,"-gains"), paste0(cohort1_ID,'-losses'), paste0(cohort2_ID,"-gains"), paste0(cohort2_ID,'-losses')),
           col = c(group1.gaincol,group1.losscol,group2.gaincol,group2.losscol),
           lwd=2.5, lty=c(1,1,NA,NA),
           pch = c(NA, NA, 15, 15), pt.cex=2.5,
           horiz = FALSE, ncol=2)
            # ncol = 2
            # bty = 'n', xpd = NA, pt.cex = 2, fill='white',
    ###############
 }

plotPFDR <- function(path.to.PFDRdata, cohort1_ID, cohort2_ID){
    load(path.to.PFDRdata)
    ###############
    # log-transform data for plotting
	# p en fdrs for losses
	losses_pvals 	<- fdrs_loss$pvalue
	losses_fdrs 	<- fdrs_loss$fdr	
	# replace nas with 1 (I don't know if this step is necessary)
	losses_match 	<- match(rownames(datainfo),rownames(pvs_loss$info))	
	losses_match[which(!is.na(losses_match))] <- as.vector(losses_pvals[losses_match[which(!is.na(losses_match))]])
	losses_match[which(is.na(losses_match))] <- 1
	losses 			<- losses_match
	# p en fdrs for losses
	gains_pvals 	<- fdrs_gain$pvalue
	gains_fdrs 		<- fdrs_gain$fdr
	# replace nas with 1 (I don't know if this step is necessary)
	gains_match 	<- match(rownames(datainfo),rownames(pvs_gain$info))
	gains_match[which(!is.na(gains_match))] <- as.vector(gains_pvals[gains_match[which(!is.na(gains_match))]])
	gains_match[which(is.na(gains_match))] <- 1
	gains 			<- gains_match
	# convert P values to log scale:
	p_gains_log	<- -log10(gains)
	p_losses_log<-  log10(losses)
	gains[is.infinite(gains)] 	<- 5
	losses[is.infinite(losses)] <- -5
	# convert FDR values to log scale:
	fdr_gains_log	<- -log10(gains_fdrs)
	fdr_losses_log	<-  log10(losses_fdrs)
	fdr_gains_log[is.infinite(fdr_gains_log)] 	<- 5
	fdr_losses_log[is.infinite(fdr_losses_log)] <- -5

    ###############
    # parameter settings : TODO: optional: add as input
	
	misscol=NA	
	build='GRCh37'
	#y.limit (see below)		
	
	# filename & plot_title	
	#filename <- 'profiles/P_FDR_EFvsLR.png'
	#png(file= filename, res=300, width=14, height=7, unit='in')

	# colors
	p.col 		<- "blue" 
	# col2rgb("dodgerblue1") 30, 144, 255
	# col2rgb("dodgerblue1") 0, 0, 255
	p.col.trans <- rgb(0/255, 0/255, 255/255, alpha=0.5)
	fdr.col <- "gold"
	trans.white <- rgb(255/255, 255/255, 255/255, alpha=0.5)
	
    ###############
    # coordinate settings

	# region coordinates
	chrom 	<- fdrs_gain[,1]
	pos 	<- fdrs_gain[,2]
  	pos2 	<- fdrs_gain[,3]
	# chromo info
	uni.chrom <- unique(chrom)
	chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
 	chrom.ends <- integer()
  	cumul <- 0

	# adjust positions, based on chromosome lengths
  	for (j in uni.chrom) {
        pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]		# each bpstart that is not on chr1, gets the length of the bpstart + cum sum of all previous chromos, but that is not what is says..
        pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
        cumul <- cumul + chrom.lengths[as.character(j)]
        chrom.ends <- c(chrom.ends, cumul)
  	}
  	names(chrom.ends) <- names(chrom.lengths)

	# new vector for vertical line segments
	# not all regions are adjacent; make a dataframe with adjacent regions
	adjacent <- c()
	for (i in 1:length(pos)){
	    z <- pos[i+1]-1==pos2[i]
	    adjacent <- c(adjacent, z)
	}
	not.adjacent <- which(adjacent==FALSE)
	row.entries <- 1+ not.adjacent
	chrom.entries <- chrom[not.adjacent]
	start.entries <- 1+ pos2[not.adjacent]
	end.entries <- pos[row.entries] -1
	input.rows <- cbind(chrom.entries, start.entries, end.entries, 0, 0)
	allpos <- cbind(chrom,pos, pos2, fdr_gains_log, fdr_losses_log)
	dd <- as.data.frame(rbind(allpos, input.rows))
	dd <- dd[ order(dd[,1], dd[,2]), ]
	new.pos <- as.data.frame(dd)
	rownames(new.pos) <- c(1:nrow(new.pos))
	# New vector for x and y positions to draw vertical lines
	pos3 <- c(new.pos$pos, new.pos$pos2[length(new.pos$pos2)])	
	FDR.gain2 <- c(0, new.pos$fdr_gains_log)
	FDR.gain3 <- c(new.pos$fdr_gains_log, 0)
	FDR.loss2 <- c(0, new.pos$fdr_gains_log)
    FDR.loss3 <- c(new.pos$fdr_gains_log, 0)

    ###############
    # the actual plot

	y.limit <- c(-1.5,3)
	
	title 	<- paste("P-values and FDRs of CNAs: ", cohort1_ID, " vs ", cohort2_ID, sep="")
	lab_x	<- "chromosomes"
	lab_y	<- "p-value and FDR"
	
 	plot(NA, xlim=c(0, max(pos2)), ylim=y.limit, type='n', xlab='chromosomes', ylab="", main=title ,xaxs='i', xaxt='n', yaxs='i', yaxt='n') #main=main ??
	
  	if (!is.na(misscol)) {
        rect(0, -1, max(pos2), 1, col=misscol, border=NA)
        rect(pos, -1, pos2, 1, col='white', border=NA)
  	}
	# plot p values	
	rect(pos, 0, pos2, p_gains_log, col=p.col.trans, border=NA, lwd=2)
	rect(pos, 0, pos2, p_losses_log, col=p.col.trans, border=NA, lwd=2)
	# plot FDRs
	rect(pos, 0, pos2, fdr_gains_log, col=fdr.col, border=NA, lwd=1.5, lty=2, density=8) 			# grid=T, grid.lty=3, grid.col=fdr.col
	rect(pos, 0, pos2, fdr_losses_log, col=trans.white, border=NA, lwd=1.5, lty=2, density=8) 		# grid=T, grid.lty=3, grid.col=fdr.col)
	segments(pos, fdr_gains_log, pos2, fdr_gains_log, lty=1, lwd=2, col=fdr.col) 					# gains: horizontal segments 
	segments(pos3, FDR.gain2, pos3, FDR.gain3, lty=1, lwd=2, col=fdr.col) 							# gains: vertical (joining) segments 
	segments(pos, fdr_losses_log, pos2, fdr_losses_log, lty=1, lwd=2, col=fdr.col)					# losses: horizontal segments
	segments(pos3, FDR.loss2, pos3, FDR.loss3, lty=1, lwd=2, col=fdr.col) 							# losses: vertical (joining) segments
  	
	box()
	abline(h=0)
	
	if (length(chrom.ends) > 1){
        for (j in names(chrom.ends)){
            print(j)
   	        abline(v=chrom.ends[j], lty='dashed')
	        ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
	        axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
	        axis(2,at=-log10(c(0.1,0.05,0.01,0.001)),labels=c(0.1,0.05,0.01,0.001), las=1)      		#main="test",cex.axis=0.7,cex.lab=0.7)
	        axis(2,at=log10(c(1, 0.1,0.05)),labels=c(1, 0.1,0.05), las=1)     							#main="test",cex.axis=0.7,cex.lab=0.7)
        }
    }
	mtext(lab_y, side=2, line=3, at=0)
	# horizontal lines: significance thresholds P (0.05) and FDR (0.2)
	abline(h=0)
	abline(h=c(-log10(0.05),log10(0.05)),lty=3, lwd=2, col=p.col)
	abline(h=c(-log10(0.1),log10(0.1)),lty=3, lwd=2, col=fdr.col)
	# LEGEND
	legend('topleft', legend = c('p value', 'FDR'), fill=c(p.col.trans,fdr.col), border=c(p.col.trans,fdr.col),	angle=c(NA,45), density=c(NA, 20), horiz = FALSE)
	#lwd=1.5, lty=2, pt.cex=2.5, #ncol=2
}


compareCNAs2cohorts <- function(x, data1, data2, cohort1_ID, cohort2_ID, teststat, output_freqPlotCompare,output_CGHtest,output_plotPFDR,output_combined){
    png(file= output_freqPlotCompare, res=300, width=14, height=7, unit='in')
    freqPlotCompare(x, data1, data2, cohort1_ID, cohort2_ID)
    dev.off()

    Perform_CGHtest(data1,data2,cohort1_ID,cohort2_ID,teststat = teststat,output = output_CGHtest)

    png(file= output_plotPFDR, res=300, width=14, height=7, unit='in')
    plotPFDR(output_CGHtest, cohort1_ID,cohort2_ID)
    dev.off()

    png(file= output_combined, res=300, width=14, height=10, unit='in')
    par(mfrow=c(2,1), mar=c(2, 6.1, 3.1, 2.1))
    freqPlotCompare(x, data1, data2, cohort1_ID, cohort2_ID)
    par(mar=c(4.1, 6.1, 2, 2.1))
    plotPFDR(output_CGHtest, cohort1_ID,cohort2_ID)
    dev.off()
}
    


    
