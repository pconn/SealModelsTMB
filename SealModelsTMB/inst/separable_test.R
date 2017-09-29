
#generate data for a simple test of the 'SEPARABLE' function in TMB

n_cell_x = 10
n_cells = n_cell_x^2
n_t = 10
Mu = matrix(n_cells,n_t)

D = matrix(0,n_cells,n_cells)
Sigma_s = matrix(0,n_cells,n_cells)  #spatial correlation
lscale = 2

for(icell1 in 1:n_cells){
  for(icell2 in 1:n_cells){
    cur.y1 = icell1%%n_cell_x
    cur.x1 = icell1%/%n_cell_x
    cur.y2 = icell2%%n_cell_x
    cur.x2 = icell2%/%n_cell_x
    D[icell1,icell2] = sqrt((cur.y1-cur.y2)^2+(cur.x1-cur.x2)^2)
    Sigma_s[icell1,icell2]=exp(-(D[icell1,icell2])^2/(2*lscale^2))
  }
}

rho = 0.9
Sigma_t = matrix(0,n_t,n_t)
for(icell1 in 1:n_t){
  for(icell2 in 1:n_t){
    Sigma_t[icell1,icell2] = rho^(abs(icell1-icell2))
  }
}

Sigma = kronecker(Sigma_t,Sigma_s)

library(mvtnorm)
Mu = rmvnorm(1,rep(0,n_cells*n_t),sigma=Sigma)


library( RandomFields )
library( TMB )
library( INLA )
loc_s = expand.grid( "x"=1:n_cell_x, "y"=1:n_cell_x)
mesh = inla.mesh.create( loc_s )
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
Y_st = matrix(Mu,n_cells,n_t)
Data = list( "Y_st"=Y_st, "spde"=spde, "n_x"=n_cell_x,"n_s"=n_cells,"n_t"=n_t,"n_mu"=mesh$n)
Params = list("logkappa"=0,"logtau"=0,"logit_rho"=0, "Mu"=matrix(0,mesh$n,n_t))
Random = c( "Mu" )
Map = list()

setwd("C:/Users/paul.conn/git/SealModelsTMB/")
TmbFile = "C:/Users/paul.conn/git/SealModelsTMB/SealModelsTMB/src/separable_test"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )

start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=500, iter.max=500))         #
Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report()
cat(Sys.time() - start.time)
