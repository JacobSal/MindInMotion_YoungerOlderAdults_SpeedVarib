#%% ======================================================================== %%#
#%% CHK OS
print('Adding Functions');
ispc <- function() {
  sys_name <- Sys.info()["sysname"]
  if (sys_name == "Windows") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#%% ======================================================================== %%#
#%% ROUND EVEN FOR CV PROC.
round_even <- function(x,rnd_num) {
  rounded <- round(x)
  if (rounded %% rnd_num != 0) {
    if (x > rounded) {
      rounded <- rounded + 1
    } else {
      rounded <- rounded - 1
    }
    rounded <- round_even(rounded,rnd_num)
  } else{
    return(rounded)
  }
}

#%% ======================================================================== %%#
#%% GET DAT .MAT
get_mat_dat <- function(mat_dat,ext_inds=c(1:11)){
  #%% GET FORMATTING
  tmp = mat_dat[[1]];
  times = tmp[[1]][[1]][[ext_inds[11]]];
  freqs = tmp[[1]][[1]][[ext_inds[10]]];
  ntimes = length(times);
  nfreqs = length(freqs);
  fns = attr(tmp[[1]][[1]],"dimnames");
  print(str(fns[[1]]))
  
  #%% MAKE DATA LISTS
  #-- reduced data list
  tmp_dat_store <- list()
  index_l <- data.frame(
    subj_n=double(),
    cond_n=double(),
    group_n=double(),
    cluster_n=double(),
    index_n=double())
  #%% EXTRACT DATA
  for(i in 1:length(tmp)){
    tt = tmp[[i]][[1]];
    #-- data list
    tmp_dat_store <- append(tmp_dat_store,
                            list(itc_dat=matrix(tt[[ext_inds[9]]],nrow=nfreqs,ncol=ntimes)));
    #-- indexing list
    index_l <- rbind(index_l,data.frame(
      subj_n=as.numeric(tt[[ext_inds[2]]]),
      group_n=as.numeric(tt[[ext_inds[4]]]),
      cond_n=as.numeric(tt[[ext_inds[6]]]),
      cluster_n=as.numeric(tt[[ext_inds[8]]]),
      index_n=i))
  }
  return(list(itc_datl=tmp_dat_store,indexl=index_l,freqs=freqs,times=times))
}


#%% ======================================================================== %%#
#%% GET DAT FOR MCLAPPLY
get_mcl_dat_mat <- function(indexl,itc_datl,nfreqs,ntimes,
                            do_cond_baseline=FALSE){
  subjs=unique(indexl$subj_n);
  conds=unique(indexl$cond_n);
  groups=unique(indexl$group_n);
  clusters=unique(indexl$cluster_n);
  #--
  loop_vals = list();
  cnt = 1;
  for (i in 1:length(clusters)) {
    cli = clusters[i];
    tt <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% cli));
    subjs=unique(tt$subj_n)
    conds=unique(tt$cond_n)
    for(l in 1:length(conds)) {
      ci = conds[l];
      ttc <- filter_at(tt,vars('cond_n'), any_vars(. %in% ci));
      #-- baseline, if one.
      if(do_cond_baseline){
        mask_matrix <- matrix(0,nrow=nfreqs*ntimes,ncol=length(subjs));
        for(k in 1:length(subjs)){
          si = subjs[k];
          tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
          ii = tts$index;
          mask_matrix[,k] = matrix(itc_datl[[ii]],nrow=ntimes*nfreqs,ncol=1);
        }
        cmu = rowMeans(mask_matrix);
      }else{
        cmu = matrix(0,nrow=ntimes*nfreqs,ncol=1);
      }
      #-- loop subjects
      for(k in 1:length(subjs)){
        si = subjs[k];
        tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
        #-- subtract baseline, if one.
        ii = tts$index;
        itcd = matrix(itc_datl[[ii]],nrow=ntimes*nfreqs,ncol=1);
        ttd = itcd-cmu;
        #-- assign to list
        loop_vals <- cbind(loop_vals,
                           list(list(
                             tf_dat=ttd,
                             cli=cli,
                             ci=ci,
                             si=si,
                             freqN=nfreqs,
                             timeN=ntimes)))
        cnt = cnt + 1;
      }
    }
  }
  return(loop_vals)
}

#%% ======================================================================== %%#
#%% GET DAT CSV
# get_mcl_dat <- function(dtbl,clusters){
#   freqs = unique(dtbl$itc_freq);
#   times = unique(dtbl$itc_time);
#   subjs=unique(dtbl$subj_n)
#   conds=unique(dtbl$cond_n)
#   groups=unique(dtbl$group_n)
#   #--
#   loop_vals = list();
#   cnt = 1;
#   for (i in 1:length(clusters)) {
#     cli = clusters[i];
#     tt <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% cli));
#     subjs=unique(tt$subj_n)
#     conds=unique(tt$cond_n)
#     for(l in 1:length(conds)) {
#       ci = conds[l];
#       ttc <- filter_at(tt,vars('cond_n'), any_vars(. %in% ci));
#       # mask_matrix <- matrix(0L,nrow=length(freqs)*length(times),ncol=length(subjs));
#       # for(k in 1:length(subjs)){
#       #   si = subjs[k];
#       #   tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
#       #   mask_matrix[,k] = tts$itc_dat;
#       # }
#       # cmu = rowMeans(mask_matrix);
#       for(k in 1:length(subjs)){
#         si = subjs[k];
#         tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
#         # ttd = tts$itc_dat-cmu;
#         ttd = tts$itc_dat;
#         newl <- list(tf_dat=ttd,cli=cli,ci=ci,si=si,freqN=length(freqs),timeN=length(times))
#         loop_vals <- cbind(loop_vals,list(newl))
#         cnt = cnt + 1;
#       }
#     }
#   }
#   return(loop_vals)
# }

#%% ======================================================================== %%#
# #%% TF GET DATA (csv implement)
# get_tf_dat <- function(cli,ci,si,dtbl){
#   #--
#   dtbl <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% cli));
#   dtbl <- filter_at(dtbl,vars('cond_n'), any_vars(. %in% ci));
#   tt <- filter_at(dtbl,vars('subj_n'), any_vars(. %in% si));
#   #--
#   freqs=unique(dtbl$itc_freq)
#   times=unique(dtbl$itc_time)
#   #-- get data and run
#   y = tt$itc_dat;
#   #--
#   return(list(tf_dat=y,freqs=freqs,times=times))
# }

#%% ======================================================================== %%#
get_tf_dat_mat <- function(cli,ci,si,indexl,itc_dat){
  #--
  indexl <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% cli));
  indexl <- filter_at(indexl,vars('cond_n'), any_vars(. %in% ci));
  tt <- filter_at(indexl,vars('subj_n'), any_vars(. %in% si));
  #-- get data and run
  y = itc_dat[[tt$index]];
  #--
  return(list(tf_dat=y))
}

#%% ======================================================================== %%#
#%% TF PLOT FUNCTION
tf_plot <- function(pwr_dat,times,freqs,title_char,zlim_in,
                    do_cbar=TRUE,do_contour=FALSE,x_lab="Time (ms)",y_lab="Frequency (Hz)"){
  #-- set params
  COLOR_F = 120;
  C_TICK_N = 20;
  #-- set plot props
  zlim_in <- seq(from=zlim_in[1],to=zlim_in[2],length.out=C_TICK_N)
  colMap <- colorRampPalette(c("blue","yellow","red" ))(COLOR_F)
  #--
  if(length(unique(pwr_dat))==1){
    return(grob());
  }
  #--
  pdin = t(matrix(pwr_dat,nrow=length(freqs)));
  h <- levelplot(pdin,
                 row.values=times,
                 column.values=freqs,
                 aspect="fill",
                 col.regions=colMap,
                 colorkey=do_cbar,
                 contour=do_contour,
                 at=zlim_in,
                 xlab=x_lab,
                 ylab=y_lab,
                 main=list(title_char))
  
  return(h)
}

#%% ======================================================================== %%#
#%% FUSEDLASSO CV 
mfusedl2d_cv <- function(lambs,tf_dat,freqN,timeN,nlambs_test=30,block_size=3,K=5){
  #-- lambda vec
  lambda_values <- seq(min(lambs), max(lambs), length.out=nlambs_test)
  #-- get input image
  img = matrix(tf_dat,nrow=freqN);
  
  #%% INITIATE K-FOLD
  m = round_even(timeN,block_size)
  n = round_even(freqN,block_size)
  cm = timeN-m;
  cn = freqN-n;
  errors <- matrix(0, nrow=length(lambda_values), ncol=K)
  #-- Generate K folds as spatial patches
  folds <- sample(rep(1:K, length.out=(n/block_size)*(m/block_size)))
  folds <- matrix(folds, nrow=n/block_size, ncol=m/block_size)
  
  #%% RUN K-FOLD
  for (k in 1:K) {
    #-- Separate training and validation patches
    train_mask <- kronecker((folds != k), matrix(1, block_size, block_size))
    test_mask <- kronecker((folds == k), matrix(1, block_size, block_size))
    
    y_train <- img[1:n,1:m] * train_mask
    y_test <- img[1:n,1:m] * test_mask
    
    fl <- fusedlasso2d(y_train)
    nls <- fl$lambda;
    
    for (i in seq_along(lambda_values)) {
      if(lambda_values[i] < min(nls) || lambda_values[i] > max(nls)){
        errors[i,k] = NA;
      } else {
        #-- Fused lasso estimation
        beta_hat <- coef(fl,lambda=lambda_values[i])$beta
        img_est <- matrix(beta_hat,n,m)
        
        #-- Mean Squared Error on the validation set
        errors[i, k] <- mean((img_est * test_mask - y_test)^2, na.rm=TRUE)
      }
    }
  }
  cve = rowMeans(errors,na.rm=TRUE);
  return(list(cv_errors=cve,lamb_vals=lambda_values))
}

#%% ======================================================================== %%#
#%% FUSED LASSO 2D WRAPPER
mfusedl2d <- function(item,save_dir) {
  #%% PARAMETERS
  block_size = 3;
  nlambs_test = 30;
  K = 10;
  #(05/11/2025) JS, bumping to 10
  tf_dat = item$tf_dat;
  cli = item$cli;
  ci = item$ci;
  si = item$si;
  freqN = item$freqN;
  timeN = item$timeN;
  #--
  message(sprintf("cl%i-c%i-s%i) Starting mfusedl2d.",cli,ci,si));
  
  if(is.null(tf_dat) || is_empty(tf_dat)){
    return()
  }
  
  #%% CALCULATE 1ST FUSEDLASSO
  message(sprintf("cl%i-c%i-s%i) Running fused-lasso alg. ...",cli,ci,si));
  fl <- fusedlasso2d(tf_dat,
                     dim1=freqN,
                     dim2=timeN)
  lambs <- fl$lambda;
  
  #%% CROSS-VALIDATION
  message(sprintf("cl%i-c%i-s%i) Running fused-lasso CV K-fold ...",cli,ci,si));
  out <- mfusedl2d_cv(lambs,tf_dat,freqN,timeN,
                      nlambs_test=nlambs_test,
                      block_size=block_size,
                      K=K);
  #-- get best lambda and estimate from original model
  cv_errors <- out$cv_errors;
  lambda_values <- out$lamb_vals; 
  best_lambda <- lambda_values[which.min(cv_errors)];
  beta_hat <- coef(fl,lambda=best_lambda);
  
  #%% SAVE DATA
  message(sprintf("cl%i-c%i-s%i) saving data ...",cli,ci,si));
  #-- save flasso model (larger data?)
  out_dat = list(cv_errors=cv_errors,
                 lamb_vals=lambda_values,
                 bbeta=beta_hat,
                 blamb=best_lambda,
                 cli=cli,
                 si=si,
                 ci=ci)
  fname = sprintf("cl%i-s%i-c%i_flmoddat.RData",cli,si,ci);
  saveRDS(out_dat,file=file.path(save_dir,fname))
  
  #%% RETURN VALUES?
  message(sprintf("cl%i-c%i-s%i) done ...",cli,ci,si));
  # return(out_dat)
  return(NULL)
}
#(05/09/2025) JS, removing the fl model because its a really big memory
#(05/09/2025) JS, these 2 save methods both work, not entirely sure the difference but I think saveRDS is more efficient memory-wise?
# save() will overwrite any variables with the saved name when load() while saveRDS() values must be extracted upon readRDS(). saveRDS()
# and readRDS() seem to be preferred by most?

#%% ======================================================================== %%#
#%% STAT TESTS
get_stat_vecs <- function(indexl,clust_i,freqs,times,finds,tinds,
                             save_dir){ 
  #%% LOOP THROUGH CONDS & SUBJS
  tt <- filter_at(indexl,vars('cluster_n'), any_vars(. %in% clust_i));
  subjs = unique(tt$subj_n);
  conds = 1:length(unique(tt$cond_n));
  nsubjs = length(subjs)
  nconds = length(conds);
  cond_vals = c(0.25,0.50,0.75,1.0);
  freqso = freqs[finds];
  timeso = times[tinds];
  
  #--
  speed_fl_mod = array(dim=c(length(freqso),length(timeso),nsubjs,nconds))
  subj_id = array(dim=c(nsubjs*nconds))
  group_id = array(dim=c(nsubjs*nconds))
  speed_val = array(dim=c(nsubjs*nconds))
  cnt = 1;
  #--
  for(l in 1:nconds){
    ci = conds[l];
    #-- load data
    fname = sprintf("cl%i-c%i_flmaskagg.RData",clust_i,ci); 
    flmodo = readRDS(file=file.path(save_dir,fname))
    #--
    for(k in 1:nsubjs){
      si = subjs[k];
      tts <- filter_at(tt,vars('subj_n'), any_vars(. %in% si));
      #--
      ttd = matrix(flmodo$mask_mat[,k],nrow=length(freqs),ncol=length(times));
      ttd = ttd[finds,tinds];
      #--
      speed_fl_mod[,,k,l] = matrix(ttd,nrow=length(freqso),ncol=length(timeso));
      subj_id[cnt] = si;
      group_id[cnt] = tts$group_n[1];
      speed_val[cnt] = cond_vals[ci];
      #--
      cnt = cnt + 1;
    }
  }
  return(list(dat_mat=speed_fl_mod,
              subj_vals=subj_id,
              group_vals=group_id,
              speed_vals=speed_val,
              freqs=freqso,
              times=timeso,
              nfreqs=length(freqso),
              ntimes=length(timeso)))
}

#%% ======================================================================== %%#
#%% LOOP DATA
get_stat_dat_mat <- function(dat_mat,nfreqs,ntimes,
                             speed_vals,group_vals,subj_vals){
  #-- loop vals
  cnt = 1;
  loop_vals = list();
  for(i in 1:nfreqs){
    for(j in 1:ntimes){
      y <- c(dat_mat[i,j,,1], dat_mat[i,j,,2], dat_mat[i,j,,3], dat_mat[i,j,,4])
      loop_vals <- cbind(loop_vals,
                         list(list(
                           point_dat=y,
                           speeds=speed_vals,
                           groups=group_vals,
                           subjs=subj_vals,
                           i=i,
                           j=j,
                           cnt=cnt)))
      cnt = cnt + 1;
    }
  }
  return(loop_vals)
}

#%% ======================================================================== %%#
#%% PERFORM STATS
lmer_fl_intstat_ij <- function(loop_dat,DO_INT_MOD=TRUE){
  ngrps = 3;
  ncoeffs = 6;
  naterms = 3;
  nocomps = 3;
  #--
  y = loop_dat$point_dat;
  tspeed = loop_dat$speeds;
  tgrp = as.factor(loop_dat$groups);
  tsubj = as.factor(loop_dat$subjs);
  i = loop_dat$i;
  j = loop_dat$j;
  #--
  estimate <- array(NA,dim=c(ncoeffs));
  pvalue <- array(NA,dim=c(ncoeffs));
  statvalue <- array(NA,dim=c(ncoeffs));
  c_chars <- array("NA",dim=c(ncoeffs));
  anpvalue <- array(NA,dim=c(naterms));
  anstat <- array(NA,dim=c(naterms));
  anchar <- array("NA",dim=c(naterms));
  #--
  cc <- array('NA',dim=c(nocomps));
  ccp <- array(0,dim=c(nocomps));
  cce <- array(0,dim=c(nocomps));
  cccic <- array('NA',dim=c(nocomps));
  cccie <- array(0,dim=c(nocomps));
  ccciu <- array(0,dim=c(nocomps));
  cccil <- array(0,dim=c(nocomps));
  
  #%% MODELS
  tryCatch(
    {
      #-- group-speed interaction model
      if(DO_INT_MOD){
        coeffis = 1:6;
        coeffis_a = 1:3;
        ngrps = 3;
        #--
        fit <- lmer(y ~ tspeed + tgrp + tspeed:tgrp + (1|tsubj))
        # fit <- lm(y ~ tspeed + tgrp + tspeed:tgrp)
        ann <- anova(fit)
        sfit <- summary(fit)
      }else{
        #-- group only model
        coeffis = 1:4;
        coeffis_a = 1:2;
        ngrps = 3;
        #--
        fit <- lmer(y ~ tspeed + tgrp + (1|tsubj))
        # fit <- lm(y ~ tspeed + tgrp)
        ann <- anova(fit)
        sfit <- summary(fit)
      }
      
      #--
      emm = emmeans::emmeans(fit,spec=c('tspeed','tgrp'),level=0.95);
      cfis <- data.frame(emm)
      pwc = emmeans::contrast(emm,method = "pairwise")
      tbp = data.frame(pwc)
      cc = tbp$contrast
      ccp = as.numeric(tbp$p.value)
      cce = as.numeric(tbp$estimate)
      cccic = cfis$tgrp
      cccie = as.numeric(cfis$emmean)
      ccciu = as.numeric(cfis$lower.CL)
      cccil = as.numeric(cfis$upper.CL)
      #--
      estimate[coeffis] <- as.numeric(sfit$coefficients[coeffis, "Estimate"])
      pvalue[coeffis] <- as.numeric(sfit$coefficients[coeffis, "Pr(>|t|)"])
      statvalue[coeffis] <- as.numeric(sfit$coefficients[coeffis, "t value"])
      c_chars[coeffis] <- rownames(sfit$coefficients[coeffis,])
      anpvalue[coeffis_a] <- as.numeric(ann[coeffis_a,"Pr(>F)"])
      anstat[coeffis_a] <- as.numeric(ann[coeffis_a,"F value"])
      anchar[coeffis_a] <- rownames(ann)
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
    }
    )
  return(list(est=estimate,
              pv=pvalue,
              cch=c_chars,
              stat=statvalue,
              apv=anpvalue,
              astat=anstat,
              acch=anchar,
              pwcc = cc,
              pwccp = ccp,
              pwcce = cce,
              cic = cccic,
              cie = cccie,
              ciu = ccciu,
              cil = cccil,
              i=i,
              j=j))
}

#%% ======================================================================== %%#
flstat_int_agg <- function(saves,cli,freqs,times,pinds=1:6,apinds=1:3,pwcinds=1:3){
  nfreqs = length(freqs);
  ntimes = length(times);
  lnpi = length(pinds);
  lnapi = length(apinds);
  lnpwci = length(pwcinds);
  lno = 6;
  lnoa = 3;
  lnop = 3;
  
  #%% EXTRACT DATA
  tmp <- lapply(saves,function(x) x$i)
  is <- as.numeric(tmp)
  tmp <- lapply(saves,function(x) x$j)
  js <- as.numeric(tmp)
  tmp <- lapply(saves,function(x) x$cch)
  cch <- unique(tmp); 
  cch <- cch[[1]];
  tmp <- lapply(saves,function(x) x$acch)
  acch <- unique(tmp); 
  acch <- acch[[1]];
  #--
  tmp <- lapply(saves,function(x) x$pwcc)
  pwcc <- unique(tmp);
  pwcc <- pwcc[[1]];
  tmp <- lapply(saves,function(x) x$cic)
  cic <- unique(tmp);
  cic <- cic[[1]]
  #--
  fdrp <- array(1,dim=c(nfreqs,ntimes,lnpi))
  afdrp <- array(1,dim=c(nfreqs,ntimes,lnapi))
  pwcfdrp <- array(1,dim=c(nfreqs,ntimes,lnpwci))
  #--
  pvo <- array(1,dim=c(nfreqs,ntimes,lno))
  esto <- array(0,dim=c(nfreqs,ntimes,lno))
  stato <- array(0,dim=c(nfreqs,ntimes,lno))
  apvo <- array(1,dim=c(nfreqs,ntimes,lnoa))
  astato <- array(0,dim=c(nfreqs,ntimes,lnoa))
  #--
  pwccpo <- array(1,dim=c(nfreqs,ntimes,lnop))
  pwceo <- array(0,dim=c(nfreqs,ntimes,lnop))
  cieo <- array(0,dim=c(nfreqs,ntimes,lnop))
  ciuo <- array(0,dim=c(nfreqs,ntimes,lnop))
  cilo <- array(0,dim=c(nfreqs,ntimes,lnop))
  #--
  for(cnt in 1:length(saves)){
    ii = as.numeric(saves[[cnt]]$i);
    jj = as.numeric(saves[[cnt]]$j);
    pv = as.numeric(saves[[cnt]]$pv);
    est = as.numeric(saves[[cnt]]$est);
    stt = as.numeric(saves[[cnt]]$stat);
    apv = as.numeric(saves[[cnt]]$apv);
    astt = as.numeric(saves[[cnt]]$astat);
    #--
    pwccp = as.numeric(saves[[cnt]]$pwccp);
    pwce = as.numeric(saves[[cnt]]$pwcce);
    cie = as.numeric(saves[[cnt]]$cie);
    ciu = as.numeric(saves[[cnt]]$ciu);
    cil = as.numeric(saves[[cnt]]$cil);
    #--
    pvo[ii,jj,] = pv;
    esto[ii,jj,] = est;
    stato[ii,jj,] = stt;
    apvo[ii,jj,] = apv;
    astato[ii,jj,] = astt;
    pwccpo[ii,jj,] = pwccp;
    pwceo[ii,jj,] = pwce;
    cieo[ii,jj,] = cie;
    ciuo[ii,jj,] = ciu;
    cilo[ii,jj,] = cil;
  }
  
  #%% FDR CORRECTION
  #(06/20/2025) JS, as per recommendation from Arkaprava, correcting p-values within each image is preferred
  # There isn't a direct connection between model coefficients that might make you assume depedence between them.
  pvot = pvo[,,pinds];
  apvot = apvo[,,apinds];
  #-- correct LMER values (Satterthwaite)
  for(pi in 1:length(pinds)){
    tmp <- p.adjust(matrix(pvo[,,pi],nrow=nfreqs*ntimes),
                     method="fdr",
                     n=nfreqs*ntimes);
    fdrp[,,pi] <- matrix(tmp,nrow=nfreqs,ncol=ntimes);
  }
  
  #-- correct anova values (Satterthwaite F-values)
  for(pi in 1:length(apinds)){
    tmp <- p.adjust(matrix(apvo[,,pi],nrow=nfreqs*ntimes),
                    method="fdr",
                    n=nfreqs*ntimes);
    afdrp[,,pi] <- matrix(tmp,nrow=nfreqs,ncol=ntimes);
  }
  
  #-- correct pwc values 
  for(pi in 1:length(pwcinds)){
    tmp <- p.adjust(matrix(pwccpo[,,pi],nrow=nfreqs*ntimes),
                    method="fdr",
                    n=nfreqs*ntimes);
    pwcfdrp[,,pi] <- matrix(tmp,nrow=nfreqs,ncol=ntimes);
  }
  
  #%% SAVE
  dato = list(fdrp=matrix(fdrp,nrow=nfreqs*ntimes*lnpi),
              rawp=matrix(pvo,nrow=nfreqs*ntimes*lno),
              estimate=matrix(esto,nrow=nfreqs*ntimes*lno),
              stat=matrix(stato,nrow=nfreqs*ntimes*lno),
              afdrp=matrix(afdrp,nrow=nfreqs*ntimes*lnapi),
              arawp=matrix(apvo,nrow=nfreqs*ntimes*lnoa),
              astat=matrix(astato,nrow=nfreqs*ntimes*lnoa),
              freqs=freqs,
              times=times,
              coeff_c=cch,
              acoeff_c=acch,
              pinds=pinds,
              apinds=apinds,
              pwcinds=pwcinds,
              pwc_c=pwcc,
              pwc_est=matrix(pwceo,nrow=nfreqs*ntimes*lnpwci),
              pwcfdrp=matrix(pwcfdrp,nrow=nfreqs*ntimes*lnpwci),
              ci_c=cic,
              ci_err=matrix(cieo,nrow=nfreqs*ntimes*lnop),
              ci_up=matrix(ciu,nrow=nfreqs*ntimes*lnop),
              ci_lw=matrix(cil,nrow=nfreqs*ntimes*lnop));
  return(dato)
}

#%% ======================================================================== %%#
#%% GET ONE COND STAT DATA
get_statdat_onesubj <- function(dat_mat,nfreqs,ntimes,speed_i,cond_i,
                             speed_vals,group_vals,subj_vals){
  #-- loop vals
  indso <- speed_vals == speed_i
  cnt = 1;
  loop_vals = list();
  for(i in 1:nfreqs){
    for(j in 1:ntimes){
      y <- dat_mat[i,j,,cond_i]
      loop_vals <- cbind(loop_vals,
                         list(list(
                           point_dat=y,
                           speeds=speed_vals[indso],
                           groups=group_vals[indso],
                           subjs=subj_vals[indso],
                           i=i,
                           j=j,
                           cnt=cnt)))
      cnt = cnt + 1;
    }
  }
  return(loop_vals)
}

#%% ======================================================================== %%#
#%% RUN ONE COND STAT DATA
flstat_onecond_ij <- function(loop_dat){
  y = loop_dat$point_dat;
  tspeed = loop_dat$speeds;
  tgrp = as.factor(loop_dat$groups);
  tsubj = as.factor(loop_dat$subjs);
  estimate <- vector(mode="double",3);
  pvalue <- vector(mode="double",3);
  i = loop_dat$i;
  j = loop_dat$j;
  cnt = loop_dat$cnt;
  
  #%% MODELS
  #-- group-speed interaction model
  fit <- lm(y ~ tgrp)
  summary_fit <- summary(fit)
  #--
  estimate[1] <- summary_fit$coefficients["tgrp2", "Estimate"]
  pvalue[1] <- summary_fit$coefficients["tgrp2", "Pr(>|t|)"]
  #--
  estimate[2] <- summary_fit$coefficients["tgrp3", "Estimate"]
  pvalue[2] <- summary_fit$coefficients["tgrp3", "Pr(>|t|)"]
  
  return(list(est=estimate,pv=pvalue,i=i,j=j,cnt=cnt))
}