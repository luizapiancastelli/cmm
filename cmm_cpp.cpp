#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
IntegerVector z_pi(IntegerVector pi, IntegerVector z){
  int n = pi.length();
  IntegerVector output(n);
  for(int i=0; i < n; i++){
    output[i] = z[pi[i]-1];
  }
  return(output);
}

// [[Rcpp::export]]
IntegerVector n_z(IntegerVector vec, int L){
  
  int n = vec.length();
  IntegerVector tab(L);
  
  for(int l=1; l <= L; ++l){
    
    for(int i =0; i< n; i++){
      if(vec[i] == l){
        tab[l-1] +=  1;
      }
    }
  }
  return tab; 
}


// [[Rcpp::export]]
IntegerVector full_table_from0(IntegerVector vec, int upper){
  
  int n = vec.length();
  IntegerVector tab(upper+1);
  
  for(int l=0; l <= upper; ++l){
    
    for(int i =0; i< n; i++){
      if(vec[i] == l){
        tab[l] +=  1;
      }
    }
  }
  return tab; 
}

// [[Rcpp::export]]
IntegerVector H_reference(IntegerVector z){
  IntegerVector ref = clone(z);
  std::sort(ref.begin(), ref.end());
  return ref;
}

// [[Rcpp::export]]
int dgo_h_cpp(IntegerVector pi,IntegerVector z){
  
  IntegerVector ref = H_reference(z);
  IntegerVector z_of_pi = z_pi(pi, z);
  
  int distance = sum(z_of_pi != ref);
  
  return distance;
}



// [[Rcpp::export]]
List dgo_h_vector(IntegerVector pi,IntegerVector z){
  int n = z.length();
  int L = unique(z).length();
  
  IntegerVector d_n(n);
  IntegerVector d_l(L);
  d_l = rep(0, L);
  
  IntegerVector ref = H_reference(z);
  IntegerVector z_of_pi = z_pi(pi, z);
  
  int d_n_i;
  for(int i =0; i< n; i++){
    
    if(z_of_pi[i]==ref[i]){
      d_n_i = 0;
    } else {
      d_n_i = 1;
    }
    
    d_n[i] = d_n_i;
    d_l[z_of_pi[i]-1] = d_l[z_of_pi[i]-1] + d_n_i;
    
  }
  
  return List::create(_["d_n"] = d_n, _["d_l"] = d_l);
}

// [[Rcpp::export]]
double CMM_log_q(IntegerVector pi, IntegerVector z, double theta){  //Hamming

  double log_q= -theta*dgo_h_cpp(pi,z);
  return log_q;
  
}


// [[Rcpp::export]]
NumericVector CMM_log_q_matrix(IntegerMatrix Pi,  IntegerVector z, double theta){
  
  int n = Pi.nrow();
  NumericVector log_q(n);
  for(int i =0; i< n; i ++){
    log_q[i] = CMM_log_q(Pi(i, _), z, theta);
  }
  return log_q;
}


// [[Rcpp::export]]
List switch_cpp(IntegerVector pi){
  
  IntegerVector pi_prime = clone(pi);
  IntegerVector switch_ = sample(pi.length(), 2);
  
  pi_prime[switch_[0]-1] = pi[switch_[1]-1];
  pi_prime[switch_[1]-1] = pi[switch_[0]-1];
  
  return List::create(_["prime"]= pi_prime, _["switched"] = switch_);
  
}



// [[Rcpp::export]]
IntegerVector gibbs_CMM_h_cpp(int total_iter, IntegerVector z, double theta){
  
  int n = z.length();

  IntegerVector pi = sample(n, n); //Initial permutation
  int iter =0;
  IntegerVector pi_prime;
  
  while(iter < total_iter){
    
  pi_prime = switch_cpp(pi)["prime"];
  
  double log_alpha = CMM_log_q(pi_prime, z, theta) - CMM_log_q(pi, z, theta);
  double log_u = log(R::runif(0,1));
  
  if(log_u < log_alpha){
    pi = clone(pi_prime);
  }
  
  iter = iter +1;
  
  }
  
  return pi;
  
}

  
// [[Rcpp::export]]
IntegerMatrix rcmm_h_gibbs_cpp(int q,  double theta, IntegerVector z, int gibbs_iter){
  
  IntegerMatrix Pi(q, z.length());
  
  for(int i =0; i< q; i++){
    Pi(i,_) = gibbs_CMM_h_cpp(gibbs_iter, z,  theta);
  }
  return Pi;
}

// [[Rcpp::export]]
double adjust_proposal(double current_sd, double accept_rate, int nprops, double target_ac){
  double adj = (1.0/(1.0*nprops))*(accept_rate - target_ac);
  double new_sd = current_sd*exp(adj);
  return new_sd;
}

// [[Rcpp::export]]
NumericVector adjust_proposal_vec(NumericVector current_sds, NumericVector accept_rates, int nprops, double target_ac){
 int L = current_sds.length();
 NumericVector new_values(L);
 for(int i =0; i < L; i ++){
   new_values[i] = adjust_proposal(current_sds[i], accept_rates[i], nprops, target_ac);
 }
  return new_values;
}

// [[Rcpp::export]]
IntegerVector dgo_h_vector_from_zpi(IntegerVector z_pi){
  int n = z_pi.length();
  int L = unique(z_pi).length();
  
  IntegerVector d_n(n);
  IntegerVector d_l(L);
  d_l = rep(0, L);
  
  IntegerVector ref = H_reference(z_pi);
  
  int d_n_i;
  for(int i =0; i< n; i++){
    
    if(z_pi[i]==ref[i]){
      d_n_i = 0;
    } else {
      d_n_i = 1;
    }
    
    d_n[i] = d_n_i;
    d_l[z_pi[i]-1] = d_l[z_pi[i]-1] + d_n_i;
    
  }
  
  return d_l;
}



// [[Rcpp::export]]
double dgo_lp_cpp(IntegerVector pi,IntegerVector z, double p){
  
  IntegerVector ref = H_reference(z);
  IntegerVector z_of_pi = z_pi(pi, z);
  
  double distance = sum(pow(abs(z_of_pi - ref),p));
  return distance;
}



// [[Rcpp::export]]
double CMM_log_q_lp(IntegerVector pi, IntegerVector z, double theta, double p){  //l-p
  
  double log_q= -theta*dgo_lp_cpp(pi,z,p);
  return log_q;
  
}



// [[Rcpp::export]]
IntegerVector gibbs_CMM_lp_cpp(int total_iter, IntegerVector z, double theta, double p){
  
  int n = z.length();
  
  IntegerVector pi = sample(n, n); //Initial permutation
  int iter =0;
  IntegerVector pi_prime;
  
  while(iter < total_iter){
    
    pi_prime = switch_cpp(pi)["prime"];
    
    double log_alpha = CMM_log_q_lp(pi_prime, z, theta, p) - CMM_log_q_lp(pi, z, theta,p);
    double log_u = log(R::runif(0,1));
    
    if(log_u < log_alpha){
      pi = clone(pi_prime);
    }
    
    iter = iter +1;
    
  }
  
  return pi;
  
}





// [[Rcpp::export]]
IntegerMatrix rcmm_lp_gibbs_cpp(int q,  double theta, IntegerVector z, int gibbs_iter, double p){
  
  IntegerMatrix Pi(q, z.length());
  
  for(int i =0; i< q; i++){
    Pi(i,_) = gibbs_CMM_lp_cpp(gibbs_iter, z,  theta, p);
  }
  return Pi;
}

// [[Rcpp::export]]
IntegerVector which_cpp(IntegerVector z, int l){
  IntegerVector indexes = seq(1, z.length());
  return indexes[z ==l];
} 

// [[Rcpp::export]]
List group_from_z_cpp(IntegerVector z){
  
  int L = unique(z).length();
  int l =1; List group(L);
  
  while(l <= L){
    group[l-1] = which_cpp(z, l);
    l = l+1;
  }
  return group;
}


// [[Rcpp::export]]
IntegerVector in_range(IntegerVector x, int low, int high) {
  return x[Rcpp::Range(low, high)];
}



// [[Rcpp::export]]
List group_from_pi_cpp(IntegerVector pi, IntegerVector z){
  
  int L = unique(z).length();
  int l = 1;
  int l_inf =0; int l_sup = sum(z == l)-1;
  List group(L);
  group[0] = in_range(pi, l_inf, l_sup);
     
  while(l < L){
    l_inf = l_sup + 1;
    l_sup = sum(z <= (l+1)) -1;
    
    group[l] = in_range(pi, l_inf, l_sup);
    
    l = l + 1;
  }
  
  return group;
    
}

// [[Rcpp::export]]
int n_ij(int i, int j, List l1, List l2){
  NumericVector v1 = l1[i-1];
  NumericVector v2 = l2[j-1];
  return intersect(v1,v2).length();
}

// [[Rcpp::export]]
int n_ij_index(int index, IntegerMatrix pairs_mat, List l1, List l2){
  IntegerVector row = pairs_mat(index-1,_);
  int i = row[0];
  int j =  row[1];
  return n_ij( i,j, l1, l2);
}


// [[Rcpp::export]]
IntegerMatrix nij_combs(List l1, List l2){
  
  int L = l1.length();
  IntegerMatrix N( pow(L,2), 3);
  
  int i = 1; int j = 1; int row = 0;
  while(i <= L){
    j=1;
    while(j <= L){
      
      N(row, 0) = i;
      N(row, 1) = j;
      N(row, 2) = n_ij(i, j, l1, l2);
      
      row = row +1;
      j = j+1;
    }
   i = i+1;
  }
  
  return N;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix fun_subset_mat( Rcpp::IntegerMatrix Input_Matrix, 
                                    Rcpp::LogicalVector Input_Log_Vec) { 
  
  // Get the number of rows and columns of "Input_Matrix"
  int Mat_nrows = Input_Matrix.nrow();
  int Mat_ncols = Input_Matrix.ncol();
  
  // Define matrix "Output_Mat" with the aimed for dimensions,
  // i.e. with the number of rows corresponding to the length of "Input_Log_Vec"
  Rcpp::IntegerMatrix Output_Mat(sum(Input_Log_Vec), Mat_ncols);
  
  // Loop over the entries of "Input_Log_Vec" to see whether the
  // corresponding row of "Input_Matrix" should be included in 
  // the output matrix "Output_Mat"
  for (int i = 0, j = 0; i < Mat_nrows; i++) {
    if (Input_Log_Vec[i]) {
      Output_Mat(j, _) = Input_Matrix(i, _);
      j = j+1;
    }
  }
  
  // Return the output matrix "Output_Mat"
  return(Output_Mat);
}

// [[Rcpp::export]]
LogicalVector ij_term_filter(int i, int j, IntegerMatrix N){ //fixed i,j
 
 LogicalVector filt1 = N(_,0) > i;
 LogicalVector filt2 = N(_,1) <= j;
 
 LogicalVector res = filt1 & filt2;
 return res;

}

// [[Rcpp::export]]
LogicalVector ij_equal_filter(int i, int j, IntegerMatrix N){ //fixed i,j
  
  LogicalVector filt1 = N(_,0) == i;
  LogicalVector filt2 = N(_,1) == j;
  
  LogicalVector res = filt1 & filt2;
  return res;
  
}

// [[Rcpp::export]]
int dk_ij_term(int i, int j, IntegerMatrix N){
  
  LogicalVector filter_out = ij_equal_filter(i, j, N);
  
  LogicalVector filter_in = ij_term_filter(i, j, N);
  IntegerMatrix outside = fun_subset_mat(N, filter_out);
  
  IntegerMatrix inside =  fun_subset_mat(N, filter_in);
  int prod = outside(0,2)*sum( inside(_,2) );
  
  return(prod);
  
 //return List::create(_["outside"] = outside, _["inside"] = inside, _["val"] = prod);
  
}


// [[Rcpp::export]]
int dk_cpp(IntegerVector pi, IntegerVector z){
  
  int L = unique(z).length();
  IntegerMatrix N = nij_combs(group_from_pi_cpp(pi, z),group_from_z_cpp(z));
  int value = 0;
  int i = 1; int j = 1;
  while(i <= (L-1)){
    j=1;
    while(j <= L){
      value = value + dk_ij_term( i,  j,  N);
      j =j+1;
    }
    i = i+1;
  }
  return value;
}


// [[Rcpp::export]]
double CMM_log_q_k(IntegerVector pi, IntegerVector z, double theta){  //Hamming
  
  double log_q= -theta*dk_cpp(pi,z);
  return log_q;
  
}

// [[Rcpp::export]]
IntegerVector gibbs_CMM_k_cpp(int total_iter, IntegerVector z, double theta){
  
  int n = z.length();
  
  IntegerVector pi = sample(n, n); //Initial permutation
  int iter =0;
  IntegerVector pi_prime;
  
  while(iter < total_iter){
    
    pi_prime = switch_cpp(pi)["prime"];
    
    double log_alpha = CMM_log_q_k(pi_prime, z, theta) - CMM_log_q_k(pi, z, theta);
    double log_u = log(R::runif(0,1));
    
    if(log_u < log_alpha){
      pi = clone(pi_prime);
    }
    
    iter = iter +1;
    
  }
  
  return pi;
  
}


// [[Rcpp::export]]
IntegerMatrix rcmm_k_gibbs_cpp(int q,  double theta, IntegerVector z, int gibbs_iter){
  
  IntegerMatrix Pi(q, z.length());
  
  for(int i =0; i< q; i++){
    Pi(i,_) = gibbs_CMM_k_cpp(gibbs_iter, z,  theta);
  }
  return Pi;
}
  
