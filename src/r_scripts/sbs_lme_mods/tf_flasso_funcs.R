#%% LOAD FUNCTIONS
print('Adding Functions');
ispc <- function() {
  sys_name <- Sys.info()["sysname"]
  if (sys_name == "Windows") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#%% PAR & PLOT FUNCTIONS
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

#%% GET LOOP VALS ==============================================================
get_mat_dat <- function(mat_dat){
  #%% GET FORMATTING
  tmp = mat_dat[[1]];
  times = tmp[[1]][[1]][[11]];
  freqs = tmp[[1]][[1]][[10]];
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
                            list(itc_dat=matrix(tt[[9]],nrow=nfreqs,ncol=ntimes)));
    #-- indexing list
    index_l <- rbind(index_l,data.frame(
      subj_n=as.numeric(tt[[2]]),
      cond_n=as.numeric(tt[[6]]),
      group_n=as.numeric(tt[[4]]),
      cluster_n=as.numeric(tt[[8]]),
      index_n=i))
  }
  return(list(itc_datl=tmp_dat_store,indexl=index_l,freqs=freqs,times=times))
}

#%%
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

#%% 
get_mcl_dat <- function(dtbl,clusters){
  freqs = unique(dtbl$itc_freq);
  times = unique(dtbl$itc_time);
  subjs=unique(dtbl$subj_n)
  conds=unique(dtbl$cond_n)
  groups=unique(dtbl$group_n)
  #--
  loop_vals = list();
  cnt = 1;
  for (i in 1:length(clusters)) {
    cli = clusters[i];
    tt <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% cli));
    subjs=unique(tt$subj_n)
    conds=unique(tt$cond_n)
    for(l in 1:length(conds)) {
      ci = conds[l];
      ttc <- filter_at(tt,vars('cond_n'), any_vars(. %in% ci));
      # mask_matrix <- matrix(0L,nrow=length(freqs)*length(times),ncol=length(subjs));
      # for(k in 1:length(subjs)){
      #   si = subjs[k];
      #   tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
      #   mask_matrix[,k] = tts$itc_dat;
      # }
      # cmu = rowMeans(mask_matrix);
      for(k in 1:length(subjs)){
        si = subjs[k];
        tts <- filter_at(ttc,vars('subj_n'), any_vars(. %in% si));
        # ttd = tts$itc_dat-cmu;
        ttd = tts$itc_dat;
        newl <- list(tf_dat=ttd,cli=cli,ci=ci,si=si,freqN=length(freqs),timeN=length(times))
        loop_vals <- cbind(loop_vals,list(newl))
        cnt = cnt + 1;
      }
    }
  }
  return(loop_vals)
}

#%% TF GET DATA ================================================================
get_tf_dat <- function(cli,ci,si,dtbl){
  #--
  dtbl <- filter_at(dtbl,vars('cluster_n'), any_vars(. %in% cli));
  dtbl <- filter_at(dtbl,vars('cond_n'), any_vars(. %in% ci));
  tt <- filter_at(dtbl,vars('subj_n'), any_vars(. %in% si));
  #--
  freqs=unique(dtbl$itc_freq)
  times=unique(dtbl$itc_time)
  #-- get data and run
  y = tt$itc_dat;
  #--
  return(list(tf_dat=y,freqs=freqs,times=times))
}

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

#%% TF PLOT FUNCTION ===========================================================
tf_plot <- function(pwr_dat,times,freqs,title_char,zlim_in){
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
                 colorkey=TRUE,
                 contour=FALSE,
                 at=zlim_in,
                 xlab=c('Time'),
                 ylab=c('Frequency (Hz)'),
                 main=list(title_char))
  
  return(h)
}

#%% FUSEDLASSO CV ============================================================
mfusedl2d_cv <- function(lambs,tf_dat,freqN,timeN,nlambs_test=30,block_size=5,K=5){
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

fusedl_valid_plots <- function(tf_dat,times,freqs){
  #%% PLOT RESULTS
  # ggo <- ggplot(data.frame(lambda_values,cv_errors),
  #               aes(x=lambda_values,y=cv_errors),
  #               xlab="Lambda",
  #               ylab="Cross-Validation Error") +
  #   geom_point(pch=19) +
  #   geom_line(col="blue") +
  #   geom_vline(xintercept=best_lambda,col="red",lty=2)
  # #-- save plot
  # fname = sprintf("cl%i-s%i-c%i_cv.png",cli,si,ci); 
  # ggsave(file.path(save_dir,fname),ggo)
  
  #%% TF PLOT
  # zlim_in = range(tf_dat)
  # op <- tf_plot(tf_dat,times,freqs,"Original",zlim_in)
  # #-- plot fusedlasso mask
  # zlim_in = range(beta_hat$beta)
  # tit_in = bquote(lambda==.(sprintf("%.3f",beta_hat$lambda[1])));
  # np <- tf_plot(beta_hat$beta[,1],times,freqs,tit_in,zlim_in)
  # #-- join them
  # pl <- list(op,np)
  # grid.arrange(grid::rectGrob(),grid::rectGrob())
  # ml <- marrangeGrob(pl, nrow=2, ncol=1);
  # #-- save plot
  # fname = sprintf("cl%i-s%i-c%i_tfplot.png",cli,si,ci); 
  # ggsave(file.path(save_dir,fname), ml)
}

#%% FUSED LASSO 2D WRAPPER =====================================================
mfusedl2d <- function(item,save_dir) {
  #%% PARAMETERS
  block_size = 5;
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
  
  #%% PLOT RESULTS
  # message(sprintf("cl%i-c%i-s%i) Plotting results & saving data ...",cli,ci,si));
  # ggo <- ggplot(data.frame(lambda_values,cv_errors),
  #               aes(x=lambda_values,y=cv_errors),
  #               xlab="Lambda",
  #               ylab="Cross-Validation Error") +
  #   geom_point(pch=19) +
  #   geom_line(col="blue") +
  #   geom_vline(xintercept=best_lambda,col="red",lty=2)
  # #-- save plot
  # fname = sprintf("cl%i-s%i-c%i_cv.png",cli,si,ci);
  # ggsave(file.path(save_dir,fname),ggo)
  
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
  #(05/09/2025) JS, removing the fl model because its a really big memory
  fname = sprintf("cl%i-s%i-c%i_flmoddat.RData",cli,si,ci);
  saveRDS(out_dat,file=file.path(save_dir,fname))
  
  #(05/09/2025) JS, these 2 save methods both work, not entirely sure the difference but I think saveRDS is more efficient memory-wise?
  # save() will overwrite any variables with the saved name when load() while saveRDS() values must be extracted upon readRDS(). saveRDS()
  # and readRDS() seem to be preferred by most?
  
  #%% RETURN VALUES?
  message(sprintf("cl%i-c%i-s%i) done ...",cli,ci,si));
  # return(out_dat)
  return(NULL)
  # return(list(flmod=fl,bbeta=beta_hat,blamb=best_lambda,cv_errors=cv_errors,lamb_vals=lambda_values))
}