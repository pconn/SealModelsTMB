# Abundance estimation for multiple species from count data 
# using spatial regression with prior distributions on detection probability at
# each location sampled.  In this version counts are disassociated from species; a
# prior distribution on confusion matrix parameters provides the link to species-specific 
# counts.
# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# Eta -- spatial random effects
# Beta -- spatial regression parameters
# logtau -- precision of spatial random effects

library( RandomFields )
library( TMB )
library( INLA )
library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
library( mvtnorm)
#library(TMBdebug)


# Compile
setwd("C:/Users/paul.conn/git/SealModelsTMB/")
TmbFile = "C:/Users/paul.conn/git/SealModelsTMB/SealModelsTMB/src/st_model_misID_separable_CPIF_knot_DayAvail"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
#source('./OkhotskSeal/R/sim_funcs.R')
#source('./OkhotskSeal/R/sim_data_generic.R')

# Settings
grid_dim = c("x"=30, "y"=30)
n_knots = 64
n_samp_t = 20  #number of samples at each time step
n_cells = grid_dim[1]*grid_dim[2]
prop_sampled=0.2
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_eta = SD_x = 1
SD_delta = 0.2  #innovation SD
iid.sd=0
IID = TRUE
if(IID)iid.sd = 0.1  #simulate data w/ extra overdispersion on log scale?
phi = 0.98 # AR process on random effects
beta0 = 2
Use_REML = FALSE   # 
Spatial_sim_model = "GP_gaussian"
Spatial_model = "SPDE_GMRF"
Alpha = 1  # Smoothness for GMRF, 1 or 2 (1 is faster)
RandomSeed = 12345
n_sim = 10
n_species = 3
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category
t_steps = 20
n_samp = t_steps*n_samp_t  #total number of samples
Prop_sampled=rep(prop_sampled,n_samp)
penalty=0.1

loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
knot_pos_x = knot_pos_y = rep(0,sqrt(n_knots))
knot_incr_x = (grid_dim[1]+1)/(sqrt(n_knots)-1)
knot_incr_y = (grid_dim[2]+1)/(sqrt(n_knots)-1)
for(ix in 2:sqrt(n_knots))knot_pos_x[ix]=knot_pos_x[ix-1]+knot_incr_x
for(iy in 2:sqrt(n_knots))knot_pos_y[iy]=knot_pos_y[iy-1]+knot_incr_y
loc_k = expand.grid("x"=knot_pos_x,"y"=knot_pos_y)


#misID_fix = c(10,10,10)  #which column of confusion matrix to fix to 1.0 on the multinomial logit scale for each species

Thin = rep(0.6,n_species*n_samp)  #keeping this as a vector for use w MVNORM_t in TMB
thin_logit = qlogis(Thin)
Sigma_logit_thin = vector("list",n_species)
for(isp in 1:n_species){
  Sigma_logit_thin[[isp]]=0.05*diag(n_samp)  #recall this is variance, not SD 
}
#convert into block diagonal 
Sigma_logit_thin = as.matrix(bdiag(Sigma_logit_thin))
Sigma_logit_thin = as(Sigma_logit_thin,"dgTMatrix")

#set up misID matrix using parameters on multinomial-logit scale
misIDcols = c(1:n_obs_types)
MisID_pos_rows = c(rep(1,7),rep(2,7),rep(3,7))
MisID_zero_cols = c(3,6,9)
n_misID_par=length(MisID_pos_rows)
MisID_pos_cols = c(1,2,5,6,8,9,10,2,3,4,5,8,9,10,2,3,5,6,7,8,10)
MisID_pars = c(2.5,1.5,-0.5,0,-0.5,0,1,
              -0.5,0,2.5,1.5,-0.5,0,1,
              -0.5,0,-0.5,0,2.5,1.5,1)
MisID_real = matrix(-20,n_species,n_obs_types)
for(ipar in 1:n_misID_par)MisID_real[MisID_pos_rows[ipar],MisID_pos_cols[ipar]]=MisID_pars[ipar]
for(isp in 1:n_species){
  MisID_real[isp,MisID_zero_cols[isp]]=0
  MisID_real[isp,] = exp(MisID_real[isp,])/sum(exp(MisID_real[isp,]))
}
MisID_Sigma = matrix(0.000001,n_misID_par,n_misID_par)
diag(MisID_Sigma) = 0.004

# Results
Results = array(NA,dim=c(n_sim,n_species,6))
# columns are "sim","N_true","N_est","N_SE","Convergence","LogLik"

# Loop through replicates
counter=1
isim=1
#for(isim in 1:n_sim){
  cat(paste("Simulation ",isim,"\n"))
  set.seed( RandomSeed + isim )
  
  N = round(runif(n_species,4000,20000))
  betax = betax_prob = runif(n_species,-.5,.5)       # Impact of covariate (x) on density
  
  # Spatial model
  # Spatial model
  loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
  model_x <- RMgauss(var=SD_x^2, scale=SpatialScale)
  model_eta <- RMgauss(var=SD_eta^2, scale=SpatialScale)
  model_delta <- RMgauss(var=SD_delta^2, scale=SpatialScale)
  
  # Realization from GRF
  # x_s -- spatial variation in covariate
  Eta_s = Delta_s = array(0,dim=c(n_species,n_cells,t_steps))
  Delta_s = array(0,dim=c(n_species,n_cells,t_steps))
  if( Spatial_sim_model=="GP_gaussian"){
    x_s = RFsimulate(model=model_x, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
    for(isp in 1:n_species){
      Eta_s[isp,,1] = RFsimulate(model=model_eta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
      for(it in 2:t_steps){
        Delta_s[isp,,it]=RFsimulate(model=model_delta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
        Eta_s[isp,,it]=phi*Eta_s[isp,,it-1]+Delta_s[isp,,it]
      }
    }
  }
  
  
  X_s = matrix(1,n_cells*t_steps,1)
  X_s[,1]=x_s
  
  # Total abundance
  Ztrue_st = N_st = array(0,dim=c(n_species,n_cells,t_steps))
  total_abundance = matrix(0,n_species,t_steps)
  for(isp in 1:n_species){
    for(it in 1:t_steps){
      Ztrue_st[isp,,it]=exp( betax[isp]*x_s + Eta_s[isp,,it] + rnorm(n_cells,0,iid.sd))
      Ztrue_st[isp,,it]=Ztrue_st[isp,,it]/sum(Ztrue_st[isp,,it])
      N_st[isp,,it]=rmultinom(1,N[isp],Ztrue_st[isp,,it])
      total_abundance[isp,it]=sum(N_st[isp,,it])
    }
  }
  N_s = matrix(0,n_cells*t_steps,n_species)
  for(isp in 1:n_species)N_s[,isp]=as.numeric(N_st[isp,,])
  
  #plot_N_map_xy(N=N_st[1,,1],XY=loc_s,leg.title="Abundance")
  
  # Samping intensity
  R_s = rep(1/n_cells,n_cells)
  
  # Process for locating samples
  Mapping_i = matrix(0,n_samp,2)  #1st column cell, 2nd column time
  Y_s = matrix(0,n_cells,t_steps)
  for(it in 1:t_steps){
    Sampled = sample(1:prod(grid_dim), size=n_samp_t, replace=FALSE, prob=R_s)
    Mapping_i[(it-1)*n_samp_t+c(1:n_samp_t),]=cbind(Sampled,rep(it,n_samp_t))
    Y_s[,it]=ifelse(1:prod(grid_dim) %in% Sampled, 1, 0)
  }
  S_i = (Mapping_i[,2]-1)*n_cells+Mapping_i[,1]
  
  # Sine effect on logit(availability)
  Cur_t = 0.5*c(0:(t_steps-1))/pi
  logit <- function(x) log(x/(1-x))
  for(isp in 1:n_species){
    Sin_eff = sin(Cur_t + runif(1,0,pi))
    Sin_eff = Sin_eff-mean(Sin_eff)  # make effect 
    Sin_eff_samp = rep(Sin_eff,each=n_samp_t)
    Thin[(isp-1)*n_samp+c(1:n_samp)]=plogis(logit(Thin[(isp-1)*n_samp+c(1:n_samp)])+Sin_eff_samp)
  }
  #plot(Sin_eff)
  
  # Counting process - binomial model for counts, 
  C_i_true = matrix(0,n_samp,n_species)  #
  C_i = matrix(0,n_samp,n_obs_types)  #one entry per sample
  for(isp in 1:n_species){
    for(isamp in 1:n_samp){
      C_i_true[isamp,isp] = rbinom(1,N_s[S_i[isamp],isp],prop_sampled*Thin[(isp-1)*n_samp+isamp])
      if(C_i_true[isamp,isp]>0)C_i[isamp,]=C_i[isamp,]+rmultinom(1,C_i_true[isamp,isp],MisID_real[isp,])
    }
  }
  
    # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
    mesh = inla.mesh.create( loc_k )
      
    n_knots = mesh$n
    Dist_sk = matrix(0,n_cells,n_knots)
    for(icell in 1:n_cells){
      for(iknot in 1:n_knots){
        Dist_sk[icell,iknot]=sqrt((loc_s[icell,1]-mesh$loc[iknot,1])^2+(loc_s[icell,2]-mesh$loc[iknot,2])^2)
      }
    }
    Kmap = 0*Dist_sk
    sd_smooth = 5
    for(i in 1:n_cells){
      #Kmap[i,which(Dist_sk[i,]==min(Dist_sk[i,]))]=1  #matrix mapping cells to nearest knot
      #Kmap[i,]=dnorm(0,Dist_sk[i,],sd_smooth)
      #Kmap[i,]=Kmap[i,]/sum(Kmap[i,])
      #Kmap[i,which(Kmap[i,]<0.05)]=0
      third_largest = sort(Dist_sk[i,])[3]
      Which_largest = which(Dist_sk[i,]<=third_largest)
      Kmap[i,Which_largest]=1/Dist_sk[i,Which_largest]
      #Kmap[i,(nrow(loc_k)+1):n_knots]=0       #take out influence of knots on edges
      Kmap[i,]=Kmap[i,]/sum(Kmap[i,])
    }
    
    # Options
    Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?
    
    # Data
    spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
    X_day = matrix(0,n_samp,2)
    X_day[,1] = rep(c(1:t_steps),each=n_samp_t)
    X_day[,2] = (X_day[,1])^2 
    X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
    X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
    X_day = kronecker(diag(n_species),X_day)
    

    Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells*t_steps),"S_i"=S_i-1,"Y_s"=Y_s,"X_s"=X_s, "spde"=spde,"thin_mu_logit"=thin_logit,"Sigma_logit_thin"=Sigma_logit_thin,"X_day"=X_day,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"MisID_zero_cols"=MisID_zero_cols-1,"Kmap"=Kmap)

    # Parameters / Initial values - set close to truth for faster convergence
    Etainput_st = array(0,dim=c(n_species,mesh$n,t_steps))  #spatio-temporal random effects
    Beta_init = matrix(0,n_species,ncol(Data$X_s))
    Beta_init[,1] = betax 
    Params = list("log_N"=log(N*exp(rnorm(n_species,0,0.1))),"Beta"=Beta_init, "logtau_eta"=rep(0,n_species), "logkappa_eta"=rep(0,n_species),"logit_rho"=rep(0,n_species),"Etainput_st"=Etainput_st,"thin_logit_i"=thin_logit,"thin_beta_day"=rep(0,2*n_species),"MisID_pars"=MisID_pars)
    Params$Beta[,1]=beta0 + runif(n_species,-0.1,0.1) #set mean expected abundance close to truth for faster optimization
    if(ncol(Data$X_s)>1)Params$Beta[,2:ncol(Data$X_s)]=betax + runif(n_species,-0.1,0.1)

    # Random
    Random = c( "Etainput_st","thin_logit_i")
    if(Use_REML==TRUE) Random = c(Random,"Beta")

    #Random = NULL
    #Params$Etainput_st=NULL
    #Data$Etainput_st = Etainput_st
    
    
    # Fix parameters
    Map = list()
    #Map[["thin_sigma_log"]]=factor(rep(NA,n_species))
    #Map[["logkappa_eta"]] = factor( rep(1,length(Params[["logkappa_eta"]])) )
    #Map[["logit_rho"]] = factor( rep(2,length(Params[["logit_rho"]])) )
    
      
    # Make object
    #compile( paste0(Version,".cpp") )
    dyn.load( dynlib(TmbFile) )
    Start_time = Sys.time()
    #setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
    Obj$fn( Obj$par )
    #image(Obj$env$spHess(random=TRUE))  #look at covariance structure

        # Run
    #Lower = -Inf
    #Upper = Inf
    start.time = Sys.time()
    Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
    Upper = 50
    #Opt = Optimize(obj=Obj,newtonsteps=2,lower=-50,upper=50,loopnum=1,bias.correct=TRUE)
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
    #Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    Report = Obj$report()
    cat(Sys.time() - start.time)
    # Potentially fix random fields with zero sample or population variance
    #if( any(Report$MargSD_z<0.001) ){
      #Which = which(Report$MargSD_z<0.001)
     # Map[["logtau_z"]] = factor( ifelse(1:2==Which,NA,1:2) )
    #  if(length(Which)==2){
    #    Map[["logkappa_z"]] = factor( c(NA,NA) )
     # }
     # Map[["beta_b"]] = factor(NA)
    #  Params[["beta_b"]] = 0
    #  Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
    #  Params[["etainput_s"]][] = 0
    #  # Re-run
    #  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
    #  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
    #  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    #}
    
    #plot_N_map_xy(N=Report$Z_s[1,1:n_cells],XY=loc_s,leg.title="Abundance")
    #cur_day=10
    #plot_N_map_xy(N=Report$Z_s[1,((cur_day-1)*n_cells+1):(cur_day*n_cells)],XY=loc_s,leg.title="Abundance")

    plot(Report$Thin_trans[1,],ylim=c(0,1))
    points(Thin[1:400],col='red')    
    
    Converge=Opt$convergence
    # SD
    # don't need to do bias-correct for CPIF model
    if(Converge==0){
      Report = Obj$report()
      SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
      for(isp in 1:3)Results[counter,isp,1:5]=c(isim,total_abundance[isp],Report$total_abundance[isp,1], SD$sd[isp],Converge)
      Results[counter,,6]=-Opt$objective
    }
    if(Converge==1){
      Results[counter,,5]=1
    }   
    

    # output
    save( Results, file="Results.RData")

    counter=counter+1
      # Plot stuff
      #par( mfrow=c(1,3) )
      #plot( x=Ztrue_s[s_i], y=c_i)
      #plot( x=Ztrue_s, y=R_s)
      #plot( x=Report$Z_s, y=Report$R_s)
#}

summary((Results[,1,3]-Results[,1,2])/Results[,1,2])
summary((Results[,2,3]-Results[,2,2])/Results[,2,2])
summary((Results[,3,3]-Results[,3,2])/Results[,3,2])

# 
# ####################
# # Read results
# ####################
# TmbFile = "C:/Users/paul.conn/git/pref_sampling/" #TMB_version/inst/executables"
# setwd(TmbFile)
# load( "Results.RData")
# WhichDone = which(Results[,"Convergence"]==0)
# Results=Results[WhichDone,]
# Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
# True.b = Results[,"b"]
# True.b[which(True.b==1)]=0
# True.b[which(True.b==2)]=1
# True.b[which(True.b==3)]=3
# Results[,"b"]=True.b
# Which.est = which(Results[,"EM"]>1)
# Results[,"Bias_b"]=NA
# Results[Which.est,"Bias_b"]=(Results[Which.est,"b_est"]-Results[Which.est,"b"])
# 
# #
# pdf("B_est.pdf")
#   Which_plot = which(Results[,"EM"]>1)
#   plot(Results[Which_plot,"b_est"],xlab="Simulation",ylab=expression(hat(b)))
# dev.off()
# 
# 
# #For figures, limit to results for which |b estimate| < 10 
# Results[,"OnBoundary"]=rep(0,nrow(Results))
# Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
# 
# Results_plot=Results[Results[,"OnBoundary"]==0,]
# Results_plot[,"Est.model"]=rep("Independent",nrow(Results_plot))
# Results_plot[Results_plot[,"EM"]==2,"Est.model"]="Joint"
# 
# library(doBy)
# sd <- function(x)sqrt(var(x))
# Bias_table = summaryBy(Bias_b+RelBias_N~b+EM,data=Results_plot,FUN=c(median,mean,sd))
# 
# #plot bias_N by b_est for EstMod=2
# Cur_data=Results_plot[Results_plot[,"EM"]==2,]
# pdf("Bias_bN.pdf")
# plot(Cur_data[,"Bias_b"],Cur_data[,"RelBias_N"],xlab="Bias (b)",ylab="Relative bias (N)")
# dev.off()
# 
# 
# 
# #plot bias as function of b, estimation method
# #first, rearrange relative bias in form for ggplot
# Bias.df=Results_plot
# Bias.df[,"B"]=Bias.df[,"b"]
# Bias.df[,"Bias"]=Bias.df[,"RelBias_N"]
# 
# 
# library(ggplot2)
# #plot proportion relative bias
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),Bias))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
# bias.plot
# pdf("bias.pdf")
# bias.plot
# dev.off()
# 
# #bias of b parameter
# BiasB.df = Bias.df[which(!is.na(Bias.df[,"Bias_b"])),]
# bias.plot = ggplot(BiasB.df,aes(factor(b),Bias_b))+geom_boxplot() #+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "True b parameter", y="Absolute bias")
# bias.plot
# pdf("biasB.pdf")
# bias.plot
# dev.off()
# 
# #produce plot of bias of species-habitat relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta11))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# pdf("Bias_beta0.pdf")
# bias.plot
# dev.off()
# 
# #produce plot of bias of species-habitat relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta12))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# pdf("Bias_beta1.pdf")
# bias.plot
# dev.off()
# 
# 
# Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
# Results[which(Results[,"Converge"]==2),"Converge"]=1
# Converge_table = summaryBy(Converge~b+EM,data=Results,FUN=sum)
# 
# 
