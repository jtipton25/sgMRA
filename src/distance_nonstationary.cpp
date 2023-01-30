#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;

//' Calculate thresheld pairwise distance
//'
//' @param sq_devs An N x N_grid x 2 array of squared deviations between locations and gridcells
//' @param a A N vector of weights in the x direction
//' @param b A N vector of weights in the x direction
//' @param radius The threshold radius
// //' @param n_neighbors The expected number of neighbors based on the MRA grid
//' @return The ()thresheld) pairwise distance
//'
//' @export
//[[Rcpp::export]]
arma::mat distance_nonstationary(arma::cube& sq_devs, arma::vec& a, arma::vec& b) {
    // Rcpp::List distance_near_nonstationary(arma::cube& sq_devs, arma::vec& a, arma::vec& b, double& radius, const int& n_neighbors=68) {
    int N = sq_devs.n_rows;
    int N_grid = sq_devs.n_cols;
    // Rcpp::Rcout << "N = " << N << "\n";
    // Rcpp::Rcout << "N_Grid = " << N_grid << "\n";

    // figure out a smart way of dealing with this
    // int max_points = 2 * n_neighbors * N;
    arma::mat D(N, N_grid);

    // distance loop
    for (int i=0; i<N; i++) {
        for (int j=0; j<N_grid; j++ ) {
            D(i, j) = sqrt(sq_devs(i, j, 0) * a(i) + sq_devs(i, j, 1) * b(i));
        }
    }

    return(D);
}



