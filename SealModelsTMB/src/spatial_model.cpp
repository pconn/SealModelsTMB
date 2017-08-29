// Spatial model for seal data w/ perfect species ID
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
  DATA_FACTOR( s_i ); // Site for each sample
  DATA_IVECTOR( y_s); //indicator for which sites sampled/not sampled
  DATA_MATRIX( X_s );  //design matrix for fixed effects

  DATA_STRUCT(spde, spde_t);
  
  DATA_VECTOR(thin_mu_logit);
  DATA_MATRIX(Sigma_logit_thin);  //might want to try changing to DATA_SPARSE_MATRIX


  // Parameters 
  PARAMETER_MATRIX(Beta);              // fixed effects on density
  PARAMETER_VECTOR(logtau_z);      //species specific RE precision
  PARAMETER_VECTOR(logkappa_z);
  
  // Random effects
  PARAMETER_MATRIX( Etainput_s );           // Spatial variation in abundance
  PARAMETER_VECTOR( thin_logit_i );         // thinning "parameter" for each surveyed location (assumed MVN on logit scale)

  // derived sizes
  int n_i = C_i.row(0).size();
  int n_s = X_s.col(0).size();
  int n_b = X_s.row(0).size();
  int n_sp = logkappa_z.size(); //no. of species
  int n_eta = Etainput_s.row(0).size();

  // global stuff
  MVNORM_t<Type>neg_log_density_thin(Sigma_logit_thin);
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  vector<Type> MargSD_z(n_sp);
  vector<Type> Range_z(n_sp);
  for( int z=0; z<n_sp; z++){
    MargSD_z(z) = 1 / sqrt(4*M_PI) / exp(logtau_z(z)) / exp(logkappa_z(z));
    Range_z(z) = sqrt(8) / exp( logkappa_z(z) );
  }
  matrix<Type> Thin_i(n_sp,n_i);
  matrix<Type> Pc_i(n_sp,n_i);
  matrix<Type> ThinU_s(n_sp,n_s);
  matrix<Type> Thin_trans(n_sp,n_i);
  for( int is=0;is<n_s;is++){
    for(int isp=0;isp<n_sp;isp++){
      ThinU_s(isp,is)=A_s(is);
    }
  }
  for( int i=0;i<n_i;i++){
    for(int isp=0;isp<n_sp;isp++){
      Thin_trans(isp,i)=1/(1+exp(-thin_logit_i(n_i*isp+i)));
      Pc_i(isp,i)=1-P_i(i)*Thin_trans(isp,i);
      Thin_i(isp,i)=P_i(i)*A_s(s_i(i))*Thin_trans(isp,i);
      ThinU_s(isp,s_i(i))=ThinU_s(isp,s_i(i))*Pc_i(isp,i);
    }
  }
  vector<Type>C_sum(n_sp);
  for(int isp=0;isp<n_sp;isp++){
    C_sum(isp) = C_i.row(isp).sum();  //for finite pop correction on total abundance - (spat only!)
  }
  // Transform random effects
  matrix<Type> Eta_s( n_sp,n_eta);
  for(int isp=0;isp<n_sp;isp++){
    Eta_s.row(isp) = Etainput_s.row(isp) / exp(logtau_z(isp));
  }
  //random effects priors
  SparseMatrix<Type> Q;
  for(int isp=0;isp<n_sp;isp++){
    Q = Q_spde_generalized(spde, exp(logkappa_z(isp)), 2);  //default smoothing is alpha = 2
    jnll_comp(1) += GMRF(Q)(Etainput_s.row(isp));
  }

  // Predicted densities
  matrix<Type> Z_s(n_sp,n_s);
  matrix<Type> E_count(n_sp,n_i);
  matrix<Type> Unsampled_s(n_sp,n_s);
  vector<Type> linpredZ_s(n_s); 
  vector<Type> Beta_tmp(n_b);
  for(int isp=0;isp<n_sp;isp++){
    Beta_tmp = Beta.row(isp);
    linpredZ_s = X_s * Beta_tmp;  
    for(int is=0; is<n_s; is++){
      Z_s(isp,is) = exp( linpredZ_s(is) + Eta_s(isp,is) );
    }
  }
  

  // Probability of counts
  for(int i=0; i<n_i; i++){
    for(int isp=0;isp<n_sp;isp++){
      E_count(isp,i)=Z_s(isp,s_i(i))*Thin_i(isp,i);
      if( !isNA(C_i(isp,i)) ) jnll_comp(0) -= dpois( C_i(isp,i), E_count(isp,i), true );
    }
  }
  //std::cout<<Thin_i<<"\n\n";

  // Probability of thinning parameters (MVN prior/penalty)
  jnll_comp(2) = neg_log_density_thin(thin_logit_i-thin_mu_logit);


  // Posterior predictions of abundance in unsampled areas
  for(int isp=0;isp<n_sp;isp++){
    for(int is=0; is<n_s; is++){
      Unsampled_s(isp,is)=Z_s(isp,is)*ThinU_s(isp,is);
    }
  }

  //Type total_abundance = Z_s.sum();
  vector<Type> total_abundance(n_sp);
  for(int isp=0;isp<n_sp;isp++){
    total_abundance(isp)=Unsampled_s.row(isp).sum() + C_sum(isp);
  }

  // Total objective
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( Z_s );
  REPORT( total_abundance );
  REPORT( Range_z );
  REPORT( MargSD_z );
  REPORT( Beta );
  REPORT( Eta_s );
  REPORT( Unsampled_s );
  REPORT( jnll_comp );
  REPORT( jnll );

  // Bias correction output
  ADREPORT( total_abundance);
  if(Options_vec(0)==1){
    ADREPORT( Beta );
    ADREPORT(Z_s);
  }
  return jnll;
}
