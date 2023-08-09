fup.completeness = function(date.inclusion,
                            end.date,
                            last.fup.date,
                            status,
                            cencode =0,
                            death.date = NULL,
                            death = NULL,
                            deathcode = 1,
                            method = "clarkc"){

  dat = cbind.data.frame("date.inclusion" = date.inclusion,
                         "end.date" = end.date,
                         "last.fup.date" = last.fup.date,
                         "status" = status)

  if(method == "clarkc"){

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$theoretical.fup = as.numeric(dat$end.date - dat$date.inclusion)
    dat$theoretical.fup[dat$status !=cencode] = as.numeric(dat$last.fup.date[dat$status !=cencode] - dat$date.inclusion[dat$status !=0])

    clarkc = round(sum(dat$obs.fup)/sum(dat$theoretical.fup)*100,2)

    print(paste("The completness of follow-up is", clarkc, "% according to the method of Clark C") )
    output = list("fupc"= clarkc, "dat"= dat)
    return(output)

  }

  if(method == "percent"){

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$maxfup = as.numeric(dat$end.date - dat$date.inclusion)
    dat$completed.fup = ifelse(dat$status ==cencode & dat$obs.fup < dat$maxfup-1, 0,1)

    complete = round(sum(dat$completed.fup )/nrow(dat)*100,2)

    print(paste("The completness of follow-up is", complete, "% according to the percentage method") )

    return(list("fupc"= complete, "dat"= dat))

  }

  if(method == "mclarkc"){
    if(is.null(death) | is.null(death.date)){
      warning("Vectors death or date.date are not specified; assuming non-censoring in status is death")
      dat$death = status
      dat$death.date = last.fup.date

    } else {dat$death = death
    dat$death.date = death.date}

    if(deathcode == 1){warning("Presuming death code is 1")}

    dat$death.fup = as.numeric(dat$death.date - dat$date.inclusion)

    x = table(dat$death)
    r = sum(table(dat$death)[!rownames(x) == deathcode]) / sum(dat$death.fup)

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$missed.fup =  as.numeric(dat$end.date - dat$last.fup.date)
    dat$missed.fup = replace(dat$missed.fup, dat$missed.fup<0,0)
    dat$missed.fup[dat$death == deathcode] = 0
    dat$theoretical.fup = dat$obs.fup + dat$missed.fup
    dat$corrected.miss = dat$missed.fup - (r/2)*(dat$missed.fup)^2
    dat$corrected.miss = (1-exp(-r*dat$missed.fup))/r
    dat$corrected.theory.fup = dat$obs.fup + dat$corrected.miss


    mclarkc = round(sum(dat$obs.fup)/sum(dat$corrected.theory.fup)*100,2)
    print(paste("The completness of follow-up is", mclarkc, "% according to the modified Clark C method") )

    return(list("fupc" = mclarkc, "dat"= dat ))



  }

  if(method == "FUI"){

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$theoretical.fup = as.numeric(dat$end.date - dat$date.inclusion)
    dat$FUI = dat$obs.fup/dat$theoretical.fup
    dat$FUI[dat$status != cencode] = 1
    dat$FUI = replace(dat$FUI, dat$FUI > 1, 1)

    fui = round(mean(dat$FUI)*100,2)

    print(paste("The mean of follow-up index is", fui, "% according to the FUI method") )
    return(list("fupc" =fui, "dat" = dat))

  }

  if(method == "FPT"){

    nsubject = nrow(dat)
    int.left = rep(0, nrow(dat))
    int.right = rep(0, nrow(dat))
    int.left[dat$status != cencode] = round(as.numeric(dat$last.fup.date[dat$status != cencode] - dat$date.inclusion[dat$status != cencode]),0) - 1
    int.right[dat$status != cencode] = round(as.numeric(dat$last.fup.date[dat$status != cencode] - dat$date.inclusion[dat$status != cencode]),0)
    int.left[dat$status == cencode] = round(as.numeric(dat$last.fup.date[dat$status== cencode] - dat$date.inclusion[dat$status == cencode]),0)
    int.right[dat$status == cencode] = 365*100

    tfldata = cbind.data.frame(first.col = dat$status, int.left, int.right, dat)

    endfl<- round(max(as.numeric(dat$end.date - dat$date.inclusion))/365.25,0)


    # fit interval censored data to obtain survival estimates
    icout<-icfit(Surv(tfldata[,2]/365.25,tfldata[,3]/365.25,type="interval2")~1)
    esurv<-getsurv(c(1:endfl),icout)
    edeath<-rep(0,endfl)
    nnow<-nrow(tfldata)
    ss<-1.0
    for (j in 1:endfl) {
      edeath[j]<-nnow*(ss-esurv[[1]]$S[j])
      nnow<-nnow-edeath[j]
      ss<-esurv[[1]]$S[j]
    }


    newcensor<-ifelse(tfldata$first.col<1,tfldata[,2],tfldata[,3])/365.25
    newsurv<-(tfldata[,2]+tfldata[,3])/2.0/365.25
    # define variables
    nloss<-rep(0,endfl)
    ndeath<-rep(0,endfl)
    ncurrent<-nsubject
    ncurrentl<-nsubject
    ncurnls<-nsubject
    ncurls<-nsubject
    totalobs<-0
    totalmax<-0

    for ( i in 1:endfl) {
      # number of loss to follow-up each interval
      nloss[i]<-sum(newsurv>newcensor & (i-1)<newcensor & newcensor<i)
      # number of death each year
      ndeath[i]<-sum(newsurv<newcensor & (i-1)<newsurv & newsurv<i)
      # adding total observed person years
      totalobs<-totalobs+ncurrent-ndeath[i]-nloss[i]+nloss[i]/2+ndeath[i]/2
      # adding the maximum observed person years if no loss follow-up
      totalmax<-totalmax+ncurnls-edeath[i]/2

      # updating for each interval
      # number of participants remained in each interval
      ncurrent<-ncurrent-ndeath[i]-nloss[i]
      # number of participants remained in each interval if all losses were followed and they are no events;
      ncurrentl<-ncurrentl-ndeath[i]
      # number of participants if events are treated as remaining in the study until the end (simple alternative approach)
      ncurls<-ncurls-nloss[i]
      # number of parcicipants if all losses are followed and they have the same rate of events as the non-losses
      ncurnls<-ncurnls-edeath[i]   }

    py<-round(totalobs/totalmax, 2)*100

    print(paste("The person-time follow-up is", py, "% according to the FPT method (Xue et al, BMC medical research methodology, 2017)") )

    return(list("dat" = tfldata,"fupc" = py))



  }
  if(method == "SPT"){

    nsubject = nrow(dat)
    int.left = rep(0, nrow(dat))
    int.right = rep(0, nrow(dat))
    int.left[dat$status != cencode] = round(as.numeric(dat$last.fup.date[dat$status != cencode] - dat$date.inclusion[dat$status != cencode]),0) - 1
    int.right[dat$status != cencode] = round(as.numeric(dat$last.fup.date[dat$status != cencode] - dat$date.inclusion[dat$status != cencode]),0)
    int.left[dat$status == cencode] = round(as.numeric(dat$last.fup.date[dat$status== cencode] - dat$date.inclusion[dat$status == cencode]),0)
    int.right[dat$status == cencode] = 365*100

    tfldata = cbind.data.frame(first.col = dat$status, int.left, int.right, dat)

    endfl<- round(max(as.numeric(dat$end.date - dat$date.inclusion))/365.25,0)
    # fit interval censored data to obtain survival estimates
    icout<-icfit(Surv(tfldata[,2]/365.25,tfldata[,3]/365.25,type="interval2")~1)
    esurv<-getsurv(c(1:endfl),icout)
    edeath<-rep(0,endfl)
    nnow<-nrow(tfldata)
    ss<-1.0
    for (j in 1:endfl) {
      edeath[j]<-nnow*(ss-esurv[[1]]$S[j])
      nnow<-nnow-edeath[j]
      ss<-esurv[[1]]$S[j]
    }


    newcensor<-ifelse(tfldata$first.col<1,tfldata[,2],tfldata[,3])/365.25
    newsurv<-(tfldata[,2]+tfldata[,3])/2.0/365.25
    # define variables
    nloss<-rep(0,endfl)
    ndeath<-rep(0,endfl)
    ncurrent<-nsubject
    ncurrentl<-nsubject
    ncurnls<-nsubject
    ncurls<-nsubject
    totalobs<-0
    totalmax<-0
    totallb<-0
    satotal<-0

    for ( i in 1:endfl) {
      # number of loss to follow-up each interval
      nloss[i]<-sum(newsurv>newcensor & (i-1)<newcensor & newcensor<i)
      # number of death each year
      ndeath[i]<-sum(newsurv<newcensor & (i-1)<newsurv & newsurv<i)
      # adding total observed person years
      totalobs<-totalobs+ncurrent-ndeath[i]-nloss[i]+nloss[i]/2+ndeath[i]/2
      # adding the maximum observed person years if no loss follow-up
      totalmax<-totalmax+ncurnls-edeath[i]/2
      # simplified method gives full follow-up person years to events
      satotal<-satotal+ncurls-nloss[i]/2
      # updating for each interval
      # number of participants remained in each interval
      ncurrent<-ncurrent-ndeath[i]-nloss[i]
      # number of participants remained in each interval if all losses were followed and they are no events;
      ncurrentl<-ncurrentl-ndeath[i]
      # number of participants if events are treated as remaining in the study until the end (simple alternative approach)
      ncurls<-ncurls-nloss[i]
      # number of parcicipants if all losses are followed and they have the same rate of events as the non-losses
      ncurnls<-ncurnls-edeath[i]   }

    #simplified person time method

  }

  pys<-round(satotal/(nsubject*endfl),2)*100
  print(paste("The person-time follow-up is", pys, "% according to the SPT method (Xue et al, BMC medical research methodology, 2017)") )


  return(list("dat"= tfldata, "fupc"= pys))

}
