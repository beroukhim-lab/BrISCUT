library(ismev)
library(extRemes)
library(fitdistrplus)
library(truncdist)
library(segmented)
library(parallel)
library(dplyr)
library(reticulate)

source("plot_BrISCUT_results.R")
source('plot_fig2.R')
source_python('combine_BrISCUT_results.py')
source_python('BrISCUT_preprocessing.py')
source_python('calculate_effect_sizes.py')

clus <- makeCluster(8)

clusterEvalQ(clus, library(ismev))
clusterEvalQ(clus, library(extRemes))
clusterEvalQ(clus, library(fitdistrplus))
clusterEvalQ(clus, library(truncdist))
clusterEvalQ(clus, library(segmented))
clusterEvalQ(clus, library(parallel))
clusterEvalQ(clus, library(dplyr))

genelocs <- read.csv('docs/hg19_genes_refseq_telcentfiltered_020518.txt',sep='\t',header=T)


#'docs/SNP6_hg19_chromosome_locs_200605.txt'
abslocs <- read.table('docs/SNP6_hg19_chromosome_locs_200605.txt',sep='\t',header=T)

split3 <- function(x) {return(strsplit(x,'/')[[1]][3])}
all_tts <- unlist(lapply(list.dirs('telomere/breakpoint_files_200619',recursive=F),split3))
n <- 1000
bgdist <- 'beta'

#Other functions
filter_big_small <- function(df) {
  df <- df[df$percent>=1e-3,]
  df <- df[df$percent<=0.999,]
  df <- df[order(df$percent),]
  return(df)
}

fractionize <- function(x, chromosome, pq, telcent) {
  #abslocs <- read.table('docs/SNP6_hg19_chromosome_locs_200605.txt',sep='\t',header=T)
  #chromosome <- ifelse(is.na(as.integer(arm)),as.integer(substr(arm,1,nchar(arm)-1)),as.integer(arm))
  #pq <- ifelse(is.na(as.integer(arm)),substr(arm,nchar(arm),nchar(arm)+1), 'q')
  p_arm_length <- (abslocs[chromosome,'p_end']-abslocs[chromosome,'p_start']+1)
  q_arm_length <- (abslocs[chromosome,'q_end']-abslocs[chromosome,'q_start']+1)
  if(pq == 'p' && telcent=='tel') {
    return ((x-abslocs[chromosome,'p_start']+1)/p_arm_length) 
  } else if (pq=='q' && telcent=='tel') {
    return (1-((x-abslocs[chromosome,'q_start']+1)/q_arm_length))
  } else if (pq=='q' &&telcent=='cent') {
    return ((x-abslocs[chromosome,'q_start']+1)/q_arm_length) 
  } else if (pq=='p' && telcent=='cent') {
    return (1-((x-abslocs[chromosome,'p_start']+1)/p_arm_length) )
  }
}

unfractionize <- function(y, chromosome, pq, telcent) {
  #abslocs <- read.table('docs/SNP6_hg19_chromosome_locs_200605.txt',sep='\t',header=T)
  #chromosome <- ifelse(is.na(as.integer(arm)),as.integer(substr(arm,1,nchar(arm)-1)),as.integer(arm))
  #pq <- ifelse(is.na(as.integer(arm)),substr(arm,nchar(arm),nchar(arm)+1), 'q')
  p_arm_length <- (abslocs[chromosome,'p_end']-abslocs[chromosome,'p_start']+1)
  q_arm_length <- (abslocs[chromosome,'q_end']-abslocs[chromosome,'q_start']+1)
  if(pq == 'p' && telcent=='tel') {
    return (round(y*p_arm_length-1+abslocs[chromosome,'p_start']))
  } else if (pq=='q' && telcent=='tel') {
    return(round(q_arm_length-1+abslocs[chromosome,'q_start']-(y*q_arm_length)))
    #return (1-((x-abslocs[chromosome,'q_start']+1)/q_arm_length))
  } else if (pq=='q' &&telcent=='cent') {
    return (round(y*q_arm_length-1+abslocs[chromosome,'q_start']))
  } else if (pq=='p' && telcent=='cent') {
    return (round(p_arm_length-1+abslocs[chromosome,'p_start']-(y*p_arm_length)))
  }
}

tablify <- function(df, chromosome, pq, telcent, remove_ones = F) {
  
  if((pq=='p'&&telcent=='tel') || (pq=='q' && telcent=='cent')) {
    tabs <- as.data.frame(table(unfractionize(df$percent,chromosome,pq,telcent)))
    rownames(tabs) <- tabs$Var1
    if(remove_ones) {tabs <- tabs[tabs$Freq>1,]}
    #tabs <-droplevels(tabs)
    
  } else {
    tabs <- as.data.frame(table(unfractionize(df$percent,chromosome,pq,telcent)))
    rownames(tabs) <- tabs$Var1
    if(remove_ones) {tabs <- tabs[tabs$Freq>1,]}
    #tabs <- droplevels(tabs)
    tabs <- tabs[order(-as.integer(tabs$Var1)),]
    
  }
  return(tabs)
}

foop <- function(x) {
  if(!(x %in% c(13,14,15,21,22))) {return(paste(toString(x),'p',sep=''))}
}

fooq <- function(x) {
  if(!(x %in% c(13,14,15,21,22))) {return(paste(toString(x),'q',sep=''))}
}

errorfunc <- function(e) {
  if(!is.null(dev.list())) {dev.off()}
  return(NA)
}
#Other functions

tel <- read.csv('docs/PANCAN_tels_210521.txt',sep='\t')
tel <- filter_big_small(tel)
telemp <- tel$percent

cent <- read.csv('docs/PANCAN_cents_210521.txt',sep='\t')
cent <- filter_big_small(cent)
centemp <- cent$percent

amptel <- tel[tel$amp_del=='amp',]
deltel <- tel[tel$amp_del=='del',]
amptelemp <- amptel$percent
deltelemp <- deltel$percent

amptela <- summary(fitdist(amptelemp,'beta'))$estimate[1]
amptelb <- summary(fitdist(amptelemp,'beta'))$estimate[2]
deltela <- summary(fitdist(deltelemp,'beta'))$estimate[1]
deltelb <- summary(fitdist(deltelemp,'beta'))$estimate[2]

ampcent <- cent[cent$amp_del=='amp',]
delcent <- cent[cent$amp_del=='del',]
ampcentemp <- ampcent$percent
delcentemp <- delcent$percent

ampcenta <- summary(fitdist(ampcentemp,'beta'))$estimate[1]
ampcentb <- summary(fitdist(ampcentemp,'beta'))$estimate[2]
delcenta <- summary(fitdist(delcentemp,'beta'))$estimate[1]
delcentb <- summary(fitdist(delcentemp,'beta'))$estimate[2]

#Main BrISCUT function
run_BrISCUT <- function(arm,tt,direc,telcent,mode, thedate, ci, prefix=NULL,leftright=NULL,df1=NULL,lims=c(0,1), x=NULL, iteration=1, prior_peaks=list()) {

  chromosome <- ifelse(is.na(as.integer(arm)),as.integer(substr(arm,1,nchar(arm)-1)),as.integer(arm))
  pq <- ifelse(is.na(as.integer(arm)),substr(arm,nchar(arm),nchar(arm)+1), 'q')
  p_arm_length <- (abslocs[chromosome,'p_end']-abslocs[chromosome,'p_start']+1)
  q_arm_length <- (abslocs[chromosome,'q_end']-abslocs[chromosome,'q_start']+1)
  
  probefilename = paste('breakpoint_files_',thedate,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'.txt',sep='')

  probes <- read.table(probefilename,sep='\t',header=T,stringsAsFactors = F,fill=T)
  probes <- filter_big_small(probes)
  tabprobes <- tablify(probes,chromosome,pq,telcent,F)
  probelist <- as.numeric(as.character(tabprobes$Var1))
  #print(probelist)
  
  
  if(direc=='amp') {
    if(telcent=='tel') {
      betafit <- fitdist(amptel$percent,'beta')
      alpha = summary(betafit)$estimate[1]
      beta = summary(betafit)$estimate[2]
      emp <- amptelemp
    } 
    else if(telcent=='cent') {
      betafit <- fitdist(ampcent$percent,'beta')
      alpha = summary(betafit)$estimate[1]
      beta = summary(betafit)$estimate[2]
      emp <- ampcentemp
    } 
    else { #telcent is combined
      betafit <- fitdist(ampcomb$percent, 'beta')
      alpha = summary(betafit)$estimate[1]
      beta <- summary(betafit)$estimate[1]
      emp <- ampcombemp
    }
  } else {
    if(telcent=='tel') {
      betafit <- fitdist(deltel$percent,'beta')
      alpha = summary(betafit)$estimate[1]
      beta = summary(betafit)$estimate[2]
      emp <- deltelemp
    } else if (telcent=='cent'){
      betafit <- fitdist(delcent$percent,'beta')
      alpha = summary(betafit)$estimate[1]
      beta = summary(betafit)$estimate[2]
      emp <- delcentemp
    } else {
      betafit <- fitdist(delcomb$percent, 'beta')
      alpha = summary(betafit)$estimate[1]
      beta <- summary(betafit)$estimate[1]
      emp <- delcombemp
    }
  }

  env <- new.env()
  print(paste(arm,tt,direc,telcent))
  if(is.null(df1)) {
 
    filename = paste('breakpoint_files_',thedate,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'.txt',sep='')
    df <- read.table(filename,sep='\t',header=T,stringsAsFactors = F,fill=T)
    df$percent = as.numeric(df$percent)
    dir.create(file.path(getwd(),resultsfolder,tt),showWarnings=F)
    df1 <-df
    df1 <- filter_big_small(df1)

    ## if there's a coverage desert, spread the brekapoints until the next PROBE.
    tabs <- tablify(df1,chromosome,pq,telcent,T)
    
    newc <- df1$percent
    intnewc <- unfractionize(newc,chromosome,pq,telcent)

    for (k in rownames(tabs)) {
      z = as.integer(match(k,intnewc)) #index of first instance
      num_affected <- tabs[k,'Freq']
      if(num_affected+z > length(intnewc)) {
        newseq = seq(from=fractionize(intnewc[z],chromosome,pq,telcent), to=0.999, length.out=num_affected+1)
      }
      else {
        if((pq=='p'&&telcent=='tel') || (pq=='q' && telcent=='cent')) {
          nextprobe <- probelist[match(intnewc[z],probelist)+1]
        } else{
          nextprobe <- probelist[match(intnewc[z],probelist)+1]
        }
        intnewseq = seq(from=intnewc[z], to=nextprobe,length.out=num_affected+1)
        newseq = fractionize(intnewseq, chromosome,pq,telcent)
      }
      newc[z:(z+num_affected-1)] = newseq[1:length(newseq-1)]
    }
    df1$percent <- newc
  }
 
  ## set initial params
  ## MAKE THIS ITERATIVE! USE THE LEFT_BOUNDARY FROM THE PREVIOUS RUN AND UPDATE DF1
  
  decide_selection <- function(dis) {
    if(abs(max(dis)) > abs(min(dis))) return('positive')
    else return('negative')
  }
  
  if(length(rownames(df1))<4) next
  df1$rownamez <- seq(1,length(rownames(df1)))

  
  
  if(is.null(x)) {
    x = qbeta(seq(0,1,length=length(rownames(df1))),alpha,beta)
  }
  
  lowlim = lims[1]
  highlim = lims[2]

## THIS IS NEW 050318, makes the calculation vertical## 
  abs_dis = ptrunc(df1$percent,spec='beta',alpha,beta, a=lowlim,b=highlim)*length(df1$percent) - df1$rownamez
  dis = abs_dis/length(df1$percent)
  df1$dis =dis
  
  
  selection <- decide_selection(dis)
  

  if(selection=='positive') {
    peak_index <- which(dis==max(dis))
  } else { #negative selection
    peak_index <- which(dis==min(dis))
  }

  
  if(iteration==1) {
    prefix = substr(selection,1,1)
  } else {
    prefix <- paste(prefix,'-',substr(leftright,1,1),substr(selection,1,1),sep='')
  }
  
  
  #one-sample
  ## THIS IS NEW 091219, makes the KS to the empirical## 
  curemp <- emp[(emp<=highlim)&(emp>=lowlim)]
  #theks <- ks.test(df1$end,function (x) ptrunc(x,spec='beta',alpha,beta, a=lowlim,b=highlim))
  theks <- ks.test(df1$percent, curemp)
  p_ks <- theks$p
  stat_ks <- theks$statistic[['D']]
  
  print(paste('one-sample ks p-value is ',p_ks,', ks stat is ',stat_ks,sep=''))
  
  comboname <- paste(tt,arm,direc,telcent,prefix,toString(iteration),sep='_')
  saveRDS(list(combo=comboname,one_samp_p=p_ks, ks_stat=stat_ks), file=paste(resultsfolder,'/stats/',comboname,'.rds',sep=''))
  
  
  title <- paste(arm,tt,direc,telcent,selection,prefix,toString(ci),'iteration',toString(iteration),sep=' ')
  
  
  pdf(file=paste(resultsfolder,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'_',prefix,'_',toString(ci),'_iter',toString(iteration),'.pdf',sep=''),width=5,height=7,useDingbats=FALSE)
  par(mfrow=c(2,1),las=1)
  
  plot(df1$percent,df1$rownamez,cex=0.5,xlab='pSCNA Length',ylab='Ranked Tumors',main=paste(title,'\nks p-value = ',toString(signif(p_ks, 3),sep='')),pch=16,xlim=c(lowlim,highlim))
  lines(x,df1$rownamez,cex=0.3,pch=3)
  
  abline(v=df1$percent[peak_index]) #peak
  plot(df1$percent,dis,cex=0.5,pch=16,xlab='pSCNA Length',ylab='Distance from background',xlim=c(lowlim,highlim))
  
  #find number of peak
  peak<- df1$percent[peak_index]
  
  #TOGGLE FOR LINE
  abline(v=peak)
  vector_of_interest <- df1$percent
  
  #stopifnot(p_ks>0.05)
 if(p_ks > 0.05) {
   stop('not significant')
 }

  make_right_dist <- function(rightdf, lims) {
    rightnumrows = length(rightdf$percent)

    rlims = c(min(rightdf$percent),lims[2])
    rightx = qbeta(seq(pbeta(rlims[1],alpha,beta),pbeta(rlims[2],alpha,beta),length=rightnumrows),alpha,beta)
    

    return(list(rlims,rightx))
  }
  
  make_left_dist <- function(leftdf, lims) {
    leftnumrows = length(leftdf$percent)
    llims = c(lims[1],max(leftdf$percent))
    leftx = qbeta(seq(pbeta(llims[1],alpha,beta),pbeta(llims[2],alpha,beta),length=leftnumrows),alpha,beta)

    return(list(llims,leftx))
  }
  
  calc_bounds <-function(p) {
    
      #number of events (telomere-bounded deletions)
    right_numevents <- length(df1[df1$percent>p,]$percent)
    right_many = rep(right_numevents,n) #c(40,40,..) 1000 times. what we wanna do is take the 40 and do 
    right_gevdist <- unlist(lapply(right_many,function (x) min(qbeta(runif(min=pbeta(p,alpha,beta),max=pbeta(highlim,alpha,beta),n=x),alpha,beta))))
    right_bound <- quantile(right_gevdist,c((1-ci)/2, 1-((1-ci)/2)))
    
    left_numevents <- length(df1[df1$percent<p,]$percent)
    left_many = rep(left_numevents,n) #c(40,40,..) 1000 times. what we wanna do is take the 40 and do 
    left_gevdist <- unlist(lapply(left_many,function (x) max(qbeta(runif(max=pbeta(p,alpha,beta),min=pbeta(lowlim,alpha,beta),n=x),alpha,beta))))
    left_bound <- quantile(left_gevdist,c((1-ci)/2, 1-((1-ci)/2)))

    stuff <- list(right_bound,left_bound)
    return(stuff)

  }
  
  calc_left_boundary <- function(p,interval=1e-4) {
    left_boundary <- p
    right_high <- calc_bounds(left_boundary)[[1]][[2]]
    #if left_boundary is less than or eq to right_high, then drift left until it's not
    if (p<=right_high) {
      while (T) {
        #print(left_boundary)
        if (calc_bounds(left_boundary)[[1]][[2]] < p | left_boundary<=0) { #if goes out of bounds
          break
        }
        left_boundary <- left_boundary - interval
        
      }
      #print(left_boundary)
    }
    
    return(left_boundary)
  }
  
  calc_right_boundary <- function(p,interval=1e-4) {
    right_boundary <- p
    left_high <- calc_bounds(right_boundary)[[2]][[1]]
    #if right_boundary is greater than or eq to left_high, then drift right until it's not
    if (p>=left_high) {
      while (T) {
        if (calc_bounds(right_boundary)[[2]][[1]] > p | right_boundary>=1) { #if goes out of bounds
          break
        }
        right_boundary <- right_boundary + interval
        
      }
      #print(right_boundary)
    }
    
    return(right_boundary)
  }
  
  
  calc_both_boundaries <- function(peak) {
    if( (selection=='positive' && telcent == 'tel') || (selection=='negative' && telcent=='cent')) {
      print(paste('original peak', peak,sep=' '))
      right_peak <- peak
      left_peak <- peak
      
      # calculate right-most peak
      left_bound_a <- calc_bounds(right_peak)[[2]]
      right_bound_a <- calc_bounds(right_peak)[[1]]
      #while not at ends, and not in between good ok bounds
      counter<- match(right_peak,vector_of_interest)
      
      while (counter!=1 & 
             counter!=length(vector_of_interest) & 
             vector_of_interest[counter-1]<left_bound_a[[2]]) {     #next one to the left is far enough, so move right
        temp_right_peak <- vector_of_interest[counter+1]
        left_bound_a <- calc_bounds(temp_right_peak)[[2]]
        if (right_peak >=left_bound_a[[2]]) {break}
        right_peak<- temp_right_peak
        counter <- counter+1
      }
      
      # calculate left-most peak
      left_bound_b <- calc_bounds(left_peak)[[2]]
      right_bound_b <- calc_bounds(left_peak)[[1]]
      counter <- match(left_peak,vector_of_interest)
      while (counter!=1 & 
             counter!=length(vector_of_interest) & 
             vector_of_interest[counter+1]<=right_bound_b[[2]]) {# next one to the right is close enough, so move on left
        temp_left_peak <- vector_of_interest[counter-1]
        right_bound_b<- calc_bounds(temp_left_peak)[[1]]
        if (left_peak>right_bound_b[[2]]) {break}
        left_peak <- temp_left_peak
        counter <- counter-1
      }
      
      left_boundary <- calc_left_boundary(left_peak)
      right_boundary <- right_peak
      
      
    } else { 
  
      print(paste('original peak', peak,sep=' '))
      right_peak <- peak
      left_peak <-peak
      
      
      # calculate left-most peak
      left_bound_a <- calc_bounds(left_peak)[[2]]
      right_bound_a <- calc_bounds(left_peak)[[1]]
      counter<- match(left_peak,vector_of_interest)
      
      # NOT to be checking next one to the left is far enough, so move right; now it's
      # checking that the next one to the right is far enough, so move left.
      while (counter!=1 & 
             counter!=length(vector_of_interest) & 
             vector_of_interest[counter+1]>right_bound_a[[1]]) { #next one to the right is far enough, so "true" peak should be left.
        
        temp_left_peak <- vector_of_interest[counter-1]
        right_bound_a <- calc_bounds(temp_left_peak)[[1]]
        if (left_peak <=right_bound_a[[1]]) {break} #oriingal left peak is smaller than new lower right_bound of temp left peak
        left_peak<- temp_left_peak
        counter<-counter-1
      }
      
      
      ##### TODO
      # calculate right-most peak
      left_bound_b <- calc_bounds(right_peak)[[2]]
      right_bound_b <- calc_bounds(right_peak)[[1]]
      counter <- match(right_peak,vector_of_interest)
      while (counter!=1 & 
             counter!=length(vector_of_interest) & 
             vector_of_interest[counter-1]>=left_bound_b[[1]]) {# next one to the left is close enough, so "true" peak should be right.
        temp_right_peak <- vector_of_interest[counter+1]
        left_bound_b <- calc_bounds(temp_right_peak)[[2]]
        if (right_peak<left_bound_b[[1]]) {break} #original right peak is smaller than the new lower left bound of temp right peak
        right_peak <- temp_right_peak
        counter <- counter+1
      }
      
      right_boundary <- calc_right_boundary(right_peak)
      left_boundary <- left_peak
  
    }
    return(c(left_boundary,right_boundary))
  }
  
  boundaries <- calc_both_boundaries(peak)
  
 # putting restrictions here on left_boundary and right_boundary
  if(pq=='p') {
      if (telcent=='tel') { #p arm, simple, telomere, so 0.1-0.4 is exactly as is.
        left_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
      } else { #p arm, simple, centromere, so 0.1-0.4 is actually 0.9-0.6
        left_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
      }

  } else {
      if (telcent=='tel') { #q arm, simple, telomere, 0.1-0.4 is actually 0.9-0.6
        left_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
      } else { #q arm, simple, centromere, so 0.1-0.4 is exactly as is.
        left_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
      }
  }

  if(right_boundary_coord<left_boundary_coord) {
    low <- right_boundary_coord
    high <- left_boundary_coord
    right_boundary_coord <- high
    left_boundary_coord <- low
  }
  
  #print(probelist)
  sorted_probelist<-probelist[order(probelist)]
  left_closest_index <- which(abs(sorted_probelist-left_boundary_coord)==min(abs(sorted_probelist-left_boundary_coord)))
  print(left_closest_index)
  if(left_boundary_coord!= sorted_probelist[left_closest_index]) {
    left_boundary_coord <- sorted_probelist[left_closest_index-1][1]
    } else {
      left_boundary_coord <- sorted_probelist[left_closest_index][1]
  }

  right_closest_index <- which(abs(sorted_probelist-right_boundary_coord)==min(abs(sorted_probelist-right_boundary_coord)))
  print(right_closest_index)
  if(right_boundary_coord!= sorted_probelist[right_closest_index]) {
    right_boundary_coord <- sorted_probelist[right_closest_index+1][1]
  } else {
    right_boundary_coord <- sorted_probelist[right_closest_index][1]
  }
  
  
  if((pq=='p'&&telcent=='tel')||(pq=='q'&&telcent=='cent')) {
    left_boundary <- fractionize(left_boundary_coord,chromosome,pq,telcent)
    right_boundary <- fractionize(right_boundary_coord,chromosome,pq,telcent)
    
  } else {
    left_boundary <- fractionize(right_boundary_coord,chromosome,pq,telcent)
    right_boundary <- fractionize(left_boundary_coord,chromosome,pq,telcent)
  }
  print(prefix)
  print('left boundary')
  print(left_boundary)
  print(left_boundary_coord)
  print('right boundary')
  print(right_boundary)
  print(right_boundary_coord)
  
  abline(v=left_boundary,lty=2)
  abline(v=right_boundary,lty=2)
  
  chrgenelocs <- genelocs[genelocs$Chr == chromosome, ]
  
  filteredgenelocs <- chrgenelocs[chrgenelocs$End >= left_boundary_coord,]
  filteredgenelocs <- filteredgenelocs[filteredgenelocs$Start <= right_boundary_coord,]
  if(length(rownames(filteredgenelocs))==0) {
    rownames(chrgenelocs) <- NULL
    leftsearch <- chrgenelocs$End
    leftdistance <- min(abs(leftsearch-left_boundary_coord))
    leftindex <- which(abs(leftsearch-left_boundary_coord)==leftdistance)
    rightsearch <- chrgenelocs$Start
    rightdistance <- min(abs(rightsearch-right_boundary_coord))
    rightindex <- which(abs(rightsearch-right_boundary_coord)==rightdistance)
    if(leftdistance<rightdistance) {
      filteredgenelocs <- chrgenelocs[leftindex,]
    } else if (rightdistance < leftdistance) {
      filteredgenelocs <- chrgenelocs[rightindex,]
    } else {filteredgenelocs <- chrgenelocs[c(rightindex,leftindex),]
    }
    filteredgenelocs$Gene <- paste('[',filteredgenelocs$Gene,']',sep='')
  }
  filteredgenelocs$Start.1 <- fractionize(filteredgenelocs$Start,chromosome, pq,telcent)
  filteredgenelocs$End.1 <- fractionize(filteredgenelocs$End,chromosome, pq,telcent)
  filteredgenelocs$ks_stat <- stat_ks
  filteredgenelocs$ks_p <- p_ks
  filteredgenelocs$log10_ks_p <- -1*log10(p_ks)
  filteredgenelocs$n_events <- length(rownames(df1))
  filteredgenelocs$Peak.Start <- floor(left_boundary_coord)
  filteredgenelocs$Peak.End <- ceiling(right_boundary_coord)
  filteredgenelocs$Peak.Start.1 <- left_boundary
  filteredgenelocs$Peak.End.1 <- right_boundary
  filteredgenelocs$type <- tt
  filteredgenelocs$arm <- arm
  filteredgenelocs$direction <- direc
  filteredgenelocs$telcent <- telcent
  filteredgenelocs$negpos <- substr(selection,1,1)
  filteredgenelocs$iter <- iteration
  filteredgenelocs$conf <- ci
  
  filteredgenelocs <- filteredgenelocs[order(filteredgenelocs$Start),] 
  
  #NEW 102119, REMOVE BIGGER THAN 0.5 PEAKS
  if (right_boundary-left_boundary > 0.5) {
   # write.table(filteredgenelocs,file=paste(resultsfolder,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'_',prefix,'_',toString(ci),'_iter',toString(iteration),'.txt',sep=''),sep='\t',row.names=FALSE)
    stop('too big!')
    #filteredgenelocs <- filteredgenelocs[0,]
  }
  
  #NEW 210311, stop if boundaries are out of range
  if ((right_boundary>= 1) || (left_boundary<=0)) {
   # write.table(filteredgenelocs,file=paste(resultsfolder,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'_',prefix,'_',toString(ci),'_iter',toString(iteration),'.txt',sep=''),sep='\t',row.names=FALSE)
    stop('out of bounds!')
    #filteredgenelocs <- filteredgenelocs[0,]
  }
  
  
  if(!is.null(dev.list())) {dev.off()}
  

  check_overlaps <- function(prior_peaks, new_peak) {
    statusquo = FALSE
    for(i in prior_peaks) {
      statusquo <- max(new_peak[1],i[1]) <= min(new_peak[2],i[2])
      #print(i)
      #print(statusquo)
      if(statusquo){break}
    }
    return(statusquo)
  }
  
  #080918; got toggled, left and right.
  if(mode=='overlap') {
    rightdf <- df1[vector_of_interest>left_boundary,]
    leftdf <- df1[vector_of_interest<right_boundary,]
    print(prior_peaks)
    print(left_boundary, right_boundary)
    if(check_overlaps(prior_peaks,c(left_boundary,right_boundary))) {
      stop('repeat found!')
    }
  } else {
    rightdf <- df1[vector_of_interest>=right_boundary,]
    leftdf <- df1[vector_of_interest<=left_boundary,]
  }

  
  rightx <- make_right_dist(rightdf,lims)[[2]]
  rlims <- make_right_dist(rightdf,lims)[[1]]
  
  leftx <- make_left_dist(leftdf,lims)[[2]]
  llims <- make_left_dist(leftdf,lims)[[1]]
  
  filteredgenelocs$n_right <- length(rightx)
  filteredgenelocs$n_left <- length(leftx)
  write.table(filteredgenelocs,file=paste(resultsfolder,'/',tt,'/',tt,'_',arm,'_',direc,'_',telcent,'_',prefix,'_',toString(ci),'_iter',toString(iteration),'.txt',sep=''),sep='\t',row.names=FALSE)
  
  
  iteration = iteration +1
  prior_peaks[[prefix]] <- c(left_boundary, right_boundary)
  #print(prior_peaks)
  tryCatch(run_BrISCUT(arm,tt,direc,telcent,mode, thedate, ci,prefix,'right',rightdf,rlims,rightx,iteration,prior_peaks),error=errorfunc)
  tryCatch(run_BrISCUT(arm,tt,direc,telcent,mode, thedate, ci,prefix,'left',leftdf, llims,leftx,iteration,prior_peaks, thedate, ci),error=errorfunc)

}

BrISCUT_pipeline <- function(segfile,ttlist,thedate,ci,infoloc) {
  #RESULTS FOLDER
  resultsfolder = paste('results_',thedate,'_',toString(ci),sep='')
  dir.create(file.path(getwd(),resultsfolder),showWarnings=F)
  dir.create(file.path(getwd(),resultsfolder,'stats'),showWarnings=F)
  
  abslocs <- read.table(infoloc,sep='\t',header=T)
  take_care_arms(segfile,tt,thedate,infoloc)
  all_arms <- 1:22
  all_arms <- c(unlist(lapply(all_arms,foop)),c('13','14','15','21','22'),unlist(lapply(all_arms,fooq)))
  parallelize_tt <- function(t) {
    for (a in all_arms) {
      for (d in c('amp','del')) {
        for (tc in c('cent','tel')) {
          bar <- tryCatch(run_BrISCUT(a,t,d,tc,mode='overlap', thedate, ci),error=errorfunc)
        }
      }
    }
  }
  clusterExport(clus,as.list(unique(c(ls(.GlobalEnv),ls(environment())))),envir=environment())
  parSapply(clus, ttlist, parallelize_tt)
  stopCluster(clus)
  
  files = list.files(path = paste(resultsfolder,'/stats',sep=''))
  newthing = bind_rows(lapply(files, function (x) as.data.frame(readRDS(paste(resultsfolder,'/stats/',x,sep='')))))
  newthing$by = p.adjust(newthing$one_samp_p,method='BY')
  newthing$bonferroni = p.adjust(newthing$one_samp_p,method='bonferroni')
  write.table(newthing,file=paste(resultsfolder,paste("KS_pvalues_",thedate,'_',toString(ci),".txt",sep=''),sep='/'),sep='\t',row.names=FALSE)
  
  #call the function to combine BrISCUT results and process some files for plotting
  selectedgenes = c('CDKN2A','TERT','MYC','BAP1','TERC','TP53','ARID1A','EGFR','PPM1D')
  all_processing(thedate, ci, selectedgenes)
  
  #call the plotting
  plot_BrISCUT_results(thedate, ci, infoloc)
  
  #plotting comparison plots between different tumor types for the same arms or genes
  plot_comparisons(thedate,ci, infoloc)
  
  #calls effect sizes
  calculate_effect_sizes(thedate,ci, infoloc, amptela, amptelb, deltela, deltelb, ampcenta, ampcentb, delcenta, delcentb)
  
}


