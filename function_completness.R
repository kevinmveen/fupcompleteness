
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

    return(clarkc, dat)

  }

  if(method == "percent"){

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$maxfup = as.numeric(dat$end.date - dat$date.inclusion)
    dat$completed.fup = ifelse(dat$status ==cencode & dat$obs.fup < dat$maxfup, 0,1)

    complete = round(sum(dat$completed.fup )/nrow(dat)*100,2)

    print(paste("The completness of follow-up is", complete, "% according to the percentage method") )

    return(complete, dat)

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
    dat$corrected.theory.fup = dat$obs.fup + dat$corrected.miss


    mclarkc = round(sum(dat$obs.fup)/sum(dat$corrected.theory.fup)*100,2)

    return(mclarkc, dat)

    print(paste("The completness of follow-up is", mclarkc, "% according to the modified Clark C method") )

  }

  if(method == "FUI"){

    dat$obs.fup = as.numeric(dat$last.fup.date - dat$date.inclusion)
    dat$theoretical.fup = as.numeric(dat$end.date - dat$date.inclusion)
    dat$FUI = dat$obs.fup/dat$theoretical.fup
    dat$FUI[dat$status != cencode] = 1
    dat$FUI = replace(dat$FUI, dat$FUI > 1, 1)

    ggdensity(dat$FUI)

    fui = round(mean(dat$FUI)*100,2)

    print(paste("The mean of follow-up index is", fui, "% according to the FUI method") )


  }

}






