// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


/** Precision matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
namespace R_inla_generalized {
using namespace Eigen;
using namespace tmbutils;
using namespace R_inla;

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  	
  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
}

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_aniso_t<Type> spde, Type kappa, matrix<Type> H, int alpha=2){

  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  matrix<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
	  	  
  //Type H_trace = H(0,0)+H(1,1);
  //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices

  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 

  if( alpha==1 ) return kappa_pow2*G0 + G1_aniso;
  if( alpha==2 ) return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}
} // end namespace R_inla

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template <class Type>
Type f(Type x){
  return Type(2)/(Type(1) + exp(-Type(1) * x)) - Type(1);
}

template<class Type>
Type dbern(Type x, Type prob, int give_log=1){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace R_inla_generalized;
  using namespace Eigen;
  using namespace density;
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: compute SE? 
  
  // Data
  DATA_MATRIX( C_i );       	// Matrix of responses (counts) of each species at each sampled location 
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_IVECTOR( S_i ); // Site-time index for each sample
  DATA_IMATRIX( Y_s); //indicator for which sites sampled/not sampled
  DATA_MATRIX( X_s );  //design matrix for fixed effects

  DATA_STRUCT(spde, spde_t);
  
  DATA_VECTOR(thin_mu_logit);
  DATA_SPARSE_MATRIX(Sigma_logit_thin); 
  DATA_MATRIX(X_day);  // design matrix for extra day and day^2 effects
  DATA_ARRAY(MisID_mu);
  DATA_MATRIX(MisID_Sigma);  //variance covariance matrix for MisID_pars on mlogit scale
  
  
  DATA_MATRIX(Kmap);
  //DATA_ARRAY(Etainput_st);
    
  // Parameters 
  PARAMETER_VECTOR(log_N);       //log total abundance for each species
  PARAMETER_MATRIX(Beta);              // fixed effects on density
  PARAMETER_VECTOR(logtau_eta);      //species specific spatial precision
  PARAMETER_VECTOR(logkappa_eta);
  PARAMETER_VECTOR(logit_rho);      //species specific AR1 correlation (logit scale)
  
  // Random effects
  PARAMETER_ARRAY( Etainput_st );           // Spatial variation in abundance
  PARAMETER_VECTOR( thin_logit_i );         // thinning "parameter" for each surveyed location (assumed MVN on logit scale)
  PARAMETER_VECTOR( thin_beta_day );     //extra day and day^2 effects
  
  PARAMETER_ARRAY(MisID_pars);   //confusion matrix estimates (mlogit scale)
  
  // derived sizes
  int n_i = C_i.col(0).size();
  int n_obs_types = C_i.row(0).size();
  int n_s = Y_s.col(0).size();
  int n_t = Y_s.row(0).size();
  int n_st = X_s.col(0).size();
  int n_b = X_s.row(0).size();
  int n_sp = logkappa_eta.size(); //no. of species
  int n_eta = Etainput_st.dim(1);

  // global stuff
  vector<Type> N(n_sp);
  for(int isp=0;isp<n_sp;isp++)N(isp)=exp(log_N(isp));
  MVNORM_t<Type>neg_log_density_misID(MisID_Sigma);
  MVNORM_t<Type>neg_log_density_thin(Sigma_logit_thin);
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  vector<Type> MargSD_eta(n_sp);
  vector<Type> Range_eta(n_sp);
  vector<Type> Rho(n_sp);
  Type cur_sum;
  
  for( int isp=0; isp<n_sp; isp++){
    MargSD_eta(isp) = 1 / sqrt(4.0*M_PI) / exp(logtau_eta(isp)) / exp(logkappa_eta(isp));
    Range_eta(isp) = sqrt(8.0) / exp( logkappa_eta(isp) );
    Rho(isp) = f(logit_rho(isp));  
  }
  

  
  //set up thinning matrix
  matrix<Type> Thin_i(n_sp,n_i);
  matrix<Type> Thin_trans(n_sp,n_i);
  vector<Type> Day_effect(n_i);
  Day_effect = X_day * thin_beta_day;

  for( int i=0;i<n_i;i++){
    for(int isp=0;isp<n_sp;isp++){
      Thin_trans(isp,i)=1/(1+exp(-thin_logit_i(n_i*isp+i)-Day_effect(n_i*isp+i)));
      Thin_i(isp,i)=P_i(i)*A_s(S_i(i))*Thin_trans(isp,i);
    }
  }

  // Transform confusion matrix
  matrix<Type> Psi(n_sp,n_obs_types);
  Type tmp_sum;
  for(int isp=0;isp<n_sp;isp++){
    tmp_sum = 0.0;
    for(int itype=0;itype<n_obs_types;itype++){
      tmp_sum = tmp_sum+exp(MisID_pars(isp,itype));
    }
    for(int itype=0;itype<n_obs_types;itype++){
      Psi(isp,itype)=exp(MisID_pars(isp,itype)) / tmp_sum;
    }
  }
  //std::cout<<Psi<<'\n';

  
  //random effects priors
  SparseMatrix<Type> Q;
  array<Type> Cur_eta(n_eta,n_t);
  GMRF_t<Type> gmrf_Q;
  for(int isp=0;isp<n_sp;isp++){
     Q = Q_spde_generalized(spde, exp(logkappa_eta(isp)), 2);  //default smoothing is alpha = 2
     gmrf_Q = GMRF(Q);
     for(int ieta=0;ieta<n_eta;ieta++){
       for(int it=0;it<n_t;it++){
         Cur_eta(ieta,it)=Etainput_st(isp,ieta,it);
       }
     }
     jnll_comp(1) += SEPARABLE(AR1(Rho(isp)),gmrf_Q)(Cur_eta);
   }

  // Transform random effects
  array<Type> Eta_st( n_sp,n_eta,n_t);
  for(int isp=0;isp<n_sp;isp++){
    for(int ieta=0;ieta<n_eta;ieta++){
      for(int it=0;it<n_t;it++){
        Eta_st(isp,ieta,it) = Etainput_st(isp,ieta,it) / exp(logtau_eta(isp));
      }
    }
  }

  // Predicted densities
  matrix<Type> Pi_s(n_sp,n_st);
  matrix<Type> Z_s(n_sp,n_st);
  matrix<Type> E_count_sp(n_sp,n_i);
  matrix<Type> E_count_obs(n_i,n_obs_types);
  matrix<Type> RE(n_sp,n_st);
  vector<Type> linpredZ_s(n_st);
  vector<Type> Beta_tmp(n_b);
  for(int isp=0;isp<n_sp;isp++){
    Beta_tmp = Beta.row(isp);
    linpredZ_s = X_s * Beta_tmp;
    for(int it=0;it<n_t;it++){
      for(int is=0;is<n_s;is++){
        RE(isp,it*n_s+is)=0.0;
        for(int ieta=0;ieta<n_eta;ieta++){
          RE(isp,it*n_s+is)+=Kmap(is,ieta)*Eta_st(isp,ieta,it);
        }
      }
    }
    for(int ist=0; ist<n_st; ist++){
      Pi_s(isp,ist) = exp( linpredZ_s(ist) + RE(isp,ist) );
    }
    for(int it=0;it<n_t;it++){
      cur_sum = 0;
      for(int is=0;is<n_s;is++){
        cur_sum = cur_sum + Pi_s(isp,it*n_s+is);
      }
      for(int is=0;is<n_s;is++){
        Pi_s(isp,it*n_s+is)=Pi_s(isp,it*n_s+is)/cur_sum;
        Z_s(isp,it*n_s+is)=N(isp)*Pi_s(isp,it*n_s+is);
      }
    }
  }


  // Probability of counts
  E_count_obs = E_count_obs.setZero();
  for(int i=0; i<n_i; i++){
    for(int isp=0;isp<n_sp;isp++){
      E_count_sp(isp,i)=Z_s(isp,S_i(i))*Thin_i(isp,i);
      for(int itype=0;itype<n_obs_types;itype++){
        E_count_obs(i,itype)+=E_count_sp(isp,i)*Psi(isp,itype);
      }
    }
    for(int itype = 0;itype<n_obs_types;itype++){
      if( !isNA(C_i(i,itype)) ) jnll_comp(0) -= dpois( C_i(i,itype), E_count_obs(i,itype), true );
    }
  }

  // Probability of thinning and misID parameters (MVN prior/penalty)
  jnll_comp(2) = neg_log_density_thin(thin_logit_i-thin_mu_logit-Day_effect);
  for(int ispeff = 0; ispeff<(n_sp*2);ispeff++) jnll_comp(2) += pow(thin_beta_day(ispeff),2.0);  //ridge reg
  std::cout<<"thin dens "<<jnll_comp(2)<<'\n';
  jnll_comp(2) += neg_log_density_misID(MisID_pars-MisID_mu);

  matrix<Type> total_abundance(n_sp,n_t);
  for(int isp=0;isp<n_sp;isp++){
    for(int it=0;it<n_t;it++){
      total_abundance(isp,it)=Z_s.block(isp,it*n_s,1,n_s).sum();
    }
  }

  // Total objective
  Type jnll = jnll_comp.sum();

  //std::cout<<"Range "<<Range_eta<<"\n";
  //Reporting
   REPORT( N );
   REPORT( Z_s );
   REPORT( total_abundance );
   REPORT( Range_eta );
   REPORT( MargSD_eta );
   REPORT( Rho );
   REPORT( Beta );
   REPORT( Eta_st );
   REPORT( Psi);
   REPORT( thin_beta_day);
   REPORT( Thin_trans);
   REPORT( jnll_comp );
   REPORT( jnll );

  // Bias correction output
  ADREPORT( N);
  if(Options_vec(0)==1){
    ADREPORT( Beta );
    ADREPORT(Z_s);
  }
  return jnll;
}
