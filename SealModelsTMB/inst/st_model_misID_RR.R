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
library( mvtnorm)
library( fields)
#library(TMBdebug)


# Compile
setwd("C:/Users/paul.conn/git/SealModelsTMB/")
TmbFile = "C:/Users/paul.conn/git/SealModelsTMB/SealModelsTMB/src/st_model_misID_RR"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('c:/users/paul.conn/git/OkhotskST/OkhotskSeal/R/util_funcs.R')
#source('./OkhotskSeal/R/sim_funcs.R')
#source('./OkhotskSeal/R/sim_data_generic.R')

# Settings
grid_dim = c("x"=20, "y"=20)
n_samp_t = 40  #number of samples at each time step
n_cells = grid_dim[1]*grid_dim[2]
n_knots = 25  #assume square grid starting at 0,0 and ending at max(grid_dim$x)+1,max(grid_dim$y)+1
prop_sampled=0.25
SpatialScale = sqrt(prod(grid_dim))/5  #for covariate generation
matern_range = max(3,sqrt(n_knots)/2)
matern_smoothness = 2 #recall value of 0.5 is exponential; infty in Gaussian
SD_eta = SD_x = 1
SD_delta = 0.2  #innovation SD
phi = 0.98 # AR process on random effects
beta0 = 2
Use_REML = FALSE   
RandomSeed = 123456
n_sim = 10
n_species = 3
n_obs_types = 3*n_species+1 #3 certainty levels for each species + an 'unknown' category
t_steps = 10
n_samp = t_steps*n_samp_t  #total number of samples
Prop_sampled=rep(prop_sampled,n_samp)

#misID_fix = c(10,10,10)  #which column of confusion matrix to fix to 1.0 on the multinomial logit scale for each species

Thin = rep(0.6,n_species*n_samp)  #keeping this as a vector for use w MVNORM_t in TMB
thin_logit = logit(Thin)
Sigma_logit_thin = vector("list",n_species)
for(isp in 1:n_species){
  Sigma_logit_thin[[isp]]=0.05*diag(n_samp)  #recall this is variance, not SD 
}
#convert into block diagonal 
Sigma_logit_thin = as.matrix(bdiag(Sigma_logit_thin))
Sigma_logit_thin = as(Sigma_logit_thin,"dgTMatrix")

#set up misID matrix using parameters on multinomial-logit scale
misIDcols = c(1:n_obs_types)
MisID_pars = matrix(0,n_species,n_obs_types)
MisID_pars[1,] = c(2.5,1.5,1,-20,-0.5,0,-20,-0.5,0,1)
MisID_pars[2,] = c(-20,-0.5,0,2.5,1.5,1,-20,-0.5,0,1)
MisID_pars[3,] = c(-20,-0.5,0,-20,-0.5,0,2.5,1.5,1,1)
MisID_real = MisID_pars
for(isp in 1:n_species)MisID_real[isp,] = exp(MisID_pars[isp,])/sum(exp(MisID_pars[isp,]))
n_misID_par = n_species*n_obs_types
MisID_Sigma = matrix(0.000001,n_misID_par,n_misID_par)
diag(MisID_Sigma) = 0.004

#do some intial spatial distance type calculations
loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
knot_pos_x = knot_pos_y = rep(0,sqrt(n_knots))
knot_incr_x = (grid_dim[1]+1)/(sqrt(n_knots)-1)
knot_incr_y = (grid_dim[2]+1)/(sqrt(n_knots)-1)
for(ix in 2:sqrt(n_knots))knot_pos_x[ix]=knot_pos_x[ix-1]+knot_incr_x
for(iy in 2:sqrt(n_knots))knot_pos_y[iy]=knot_pos_y[iy-1]+knot_incr_y
loc_k = expand.grid("x"=knot_pos_x,"y"=knot_pos_y)
Dist_sk = matrix(0,n_cells,n_knots)
for(icell in 1:n_cells){
  for(iknot in 1:n_knots){
    Dist_sk[icell,iknot]=sqrt((loc_s[icell,1]-loc_k[iknot,1])^2+(loc_s[icell,2]-loc_k[iknot,2])^2)
  }
}
Dist_kk = matrix(0,n_knots,n_knots)
for(iknot1 in 1:n_knots){
  for(iknot2 in 1:n_knots){
    Dist_kk[iknot1,iknot2]=sqrt((loc_k[iknot1,1]-loc_k[iknot2,1])^2+(loc_k[iknot1,2]-loc_k[iknot2,2])^2)
  }
}
Cor_k = matrix(Matern(Dist_kk,range=matern_range,smoothness=matern_smoothness),n_knots,n_knots)
Cross_cor = matrix(Matern(Dist_sk,range=matern_range,smoothness=matern_smoothness),n_cells,n_knots)

# Results
Results_N = array(NA,dim=c(n_sim,n_species,t_steps,3))
Results_par = array(NA,dim=c(n_sim,n_species,5))  #logvar_eta, logrange, logsmooth, logvar_delta, phi_logit
# columns are "N_true","N_est","N_SE"

# Loop through replicates
counter=1
i=1
for(i in 1:n_sim){
    cat(paste("Simulation ",i,"\n"))
    set.seed( RandomSeed + i )

    betax = betax_prob = runif(n_species,-.5,.5)       # Impact of covariate (x) on density

    # Spatial model
    model_x <- RMgauss(var=SD_x^2, scale=SpatialScale)
    #model_eta <- RMgauss(var=SD_eta^2, scale=SpatialScale)
    #model_delta <- RMgauss(var=SD_delta^2, scale=SpatialScale)
    
      
    # Realization from GRF
    # x_s -- spatial variation in covariate
    A = Cross_cor %*% solve(Cor_k)
    Eta_k = Delta_k = array(0,dim=c(n_species,n_knots,t_steps))
    Eta_s = array(0,dim=c(n_species,n_cells,t_steps))
    x_s = RFsimulate(model=model_x, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
    for(isp in 1:n_species){
      Eta_k[isp,,1] = rmvnorm(1,rep(0,n_knots),Cor_k*SD_eta^2)
      Eta_s[isp,,1] = A %*% Eta_k[isp,,1]
      for(it in 2:t_steps){
        Delta_k[isp,,it]= rmvnorm(1,rep(0,n_knots),Cor_k*SD_delta^2)
        Eta_k[isp,,it]=phi*Eta_k[isp,,it-1]+Delta_k[isp,,it]
        Eta_s[isp,,it]=A %*% Eta_k[isp,,it]
      }
    }


    X_s = matrix(1,n_cells*t_steps,2)
    X_s[,2]=x_s
      
    # Total abundance
    Ztrue_st = array(0,dim=c(n_species,n_cells,t_steps))
    total_abundance = matrix(0,n_species,t_steps)
    for(isp in 1:n_species){
      for(it in 1:t_steps){
        Ztrue_st[isp,,it]=exp( beta0 + betax[isp]*x_s + Eta_s[isp,,it] )
        total_abundance[isp,it]=sum(Ztrue_st[isp,,it])
      }
    }
    Ztrue_s = matrix(0,n_cells*t_steps,n_species)
    for(isp in 1:n_species)Ztrue_s[,isp]=as.numeric(Ztrue_st[isp,,])
 
    #plot_N_map_xy(N=Ztrue_st[1,,1],XY=loc_s,leg.title="Abundance")
    
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

    # Counting process - poisson model for counts, 
    C_i_true = matrix(0,n_samp,n_species)  #
    C_i = matrix(0,n_samp,n_obs_types)  #one entry per sample
    for(isp in 1:n_species){
      Cur_lambda = Ztrue_s[S_i,isp]*prop_sampled*Thin[(isp-1)*n_samp+c(1:n_samp)]
      C_i_true[,isp] = rpois(n=n_samp,lambda=Cur_lambda)
      for(isamp in 1:n_samp){
        if(C_i_true[isamp,isp]>0)C_i[isamp,]=C_i[isamp,]+rmultinom(1,C_i_true[isamp,isp],MisID_real[isp,])
      }
    }
 
    # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
    #mesh = inla.mesh.create( loc_s )
      
       
    # Options
    Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?
    
    # Data
    #spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
    Data = list( "Options_vec"=Options_vec, "C_i"=C_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells*t_steps),"S_i"=S_i-1,"Y_s"=Y_s,"X_s"=X_s,"thin_mu_logit"=thin_logit,"Sigma_logit_thin"=Sigma_logit_thin,"MisID_mu"=MisID_pars,"MisID_Sigma"=MisID_Sigma,"Dist_kk"=Dist_kk,"Dist_sk"=Dist_sk)

    # Parameters / Initial values - set close to truth for faster convergence
    Etainput_st = array(0,dim=c(n_species,n_knots,t_steps))  #spatio-temporal random effects; 1st time step holds eta; 2nd-T hold delta
    Etainput_st[,1:n_knots,1]=Eta_k[,,1]
    Etainput_st[,1:n_knots,2:t_steps]=Delta_k[,,2:t_steps]
    Beta_init = matrix(0,n_species,ncol(Data$X_s))
    Beta_init[,1] = log(colSums(Ztrue_s[1:n_cells,])/n_cells-0.5*sqrt(apply(Ztrue_s[1:n_cells,],2,'var')))
    Beta_init[,2] = betax
    Params = list("Beta"=Beta_init, "logvar_eta"=rep(0,n_species), "logrange"=rep(0,n_species),"logsmooth"=rep(0,n_species),"logvar_delta"=rep(-1,n_species),"Etainput_st"=Etainput_st,"thin_logit_i"=thin_logit,"MisID_pars"=MisID_pars,"phi_logit"=rep(4.0,n_species))
    Params$Beta[,1]=beta0 + runif(n_species,-0.1,0.1) #set mean expected abundance close to truth for faster optimization
    if(ncol(Data$X_s)>1)Params$Beta[,2:ncol(Data$X_s)]=betax + runif(n_species,-0.1,0.1)

    # Random
    Random = c( "Etainput_st" )
    if(Use_REML==TRUE) Random = c(Random,"Beta")

    # Fix parameters
    Map = list()
    #Map[["logkappa_z"]] = factor(NA) #factor( rep(1,length(Params[["logkappa_z"]])) )
    #Map[["logtau_z"]] = factor(NA) #factor( rep(1,length(Params[["logkappa_z"]])) )
    
      
    # Make object
    #compile( paste0(Version,".cpp") )
    start.time = Sys.time()
    dyn.load( dynlib(TmbFile) )
    Start_time = Sys.time()
    #setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
    Obj$fn( Obj$par )
    cat(Sys.time() - start.time)

    # Run
    #Lower = -Inf
    #Upper = Inf
    start.time = Sys.time()
    Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
    Upper = 50
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
    Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    Report = Obj$report()
    cat(Sys.time() - start.time)
#     # Potentially fix random fields with zero sample or population variance
#     #if( any(Report$MargSD_z<0.001) ){
#       #Which = which(Report$MargSD_z<0.001)
#      # Map[["logtau_z"]] = factor( ifelse(1:2==Which,NA,1:2) )
#     #  if(length(Which)==2){
#     #    Map[["logkappa_z"]] = factor( c(NA,NA) )
#      # }
#      # Map[["beta_b"]] = factor(NA)
#     #  Params[["beta_b"]] = 0
#     #  Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
#     #  Params[["etainput_s"]][] = 0
#     #  # Re-run
#     #  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
#     #  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
#     #  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
#     #}
#     
     plot_N_map_xy(N=Ztrue_st[2,,1],XY=loc_s,leg.title="Abundance")
     plot_N_map_xy(N=Report$Z_s[2,1:n_cells],XY=loc_s,leg.title="Abundance")
     Report$total_abundance
     N_true = matrix(0,n_species,t_steps)
     for(isp in 1:n_species){
       for(it in 1:t_steps){
         N_true[isp,it] = sum(Ztrue_st[isp,,it])
       }
     }
     N_true

    Converge=Opt$convergence
    # SD
    if(Converge==0){
      Report = Obj$report()
      if( all(c("Etainput_s")%in%names(Map)) ){
        SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
        SD$unbiased$value = c("total_abundance"=Report$total_abundance)
      }else{
        start.time = Sys.time()
        SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )
        cat(Sys.time() - start.time)
      }
      Opt[["run_time"]] = Sys.time()-Start_time
    }
    
    # columns are "N_true","N_est","N_SE"
    Results_N[i,,,1]=N_true
    Results_N[i,,,2]=SD$unbiased$value
    Results_N[i,,,3]=SD$sd
    
    Results_par[i,,1]=Report$logvar_eta
    Results_par[i,,2]=Report$logrange
    Results_par[i,,3]=Report$logsmooth
    Results_par[i,,4]=Report$logvar_delta
    Results_par[i,,5]=Report$phi

    # output
    save( Results, file="Results.RData")

}
# 
# summary((Results[,1,3]-Results[,1,2])/Results[,1,2])
# summary((Results[,2,3]-Results[,2,2])/Results[,2,2])
# summary((Results[,3,3]-Results[,3,2])/Results[,3,2])
# 
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
