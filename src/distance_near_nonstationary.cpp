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
Rcpp::List distance_near_nonstationary(arma::cube& sq_devs, arma::vec& a, arma::vec& b, double& radius) {
    // Rcpp::List distance_near_nonstationary(arma::cube& sq_devs, arma::vec& a, arma::vec& b, double& radius, const int& n_neighbors=68) {
    int N = sq_devs.n_rows;
    int N_grid = sq_devs.n_cols;
    // Rcpp::Rcout << "N = " << N << "\n";
    // Rcpp::Rcout << "N_Grid = " << N_grid << "\n";

    // figure out a smart way of dealing with this
    // int max_points = 2 * n_neighbors * N;
    int max_points = N * N_grid;


    // initialize storage vector
    arma::vec I(max_points);
    arma::vec J(max_points);
    arma::vec V(max_points);
    arma::vec ddista(max_points);
    arma::vec ddistb(max_points);

    // initialize the index
    int idx = 0;

    // distance loop
    for (int i=0; i<N; i++) {
        for (int j=0; j<N_grid; j++ ) {
            // want to keep the explicit zeros
            double value = sqrt(sq_devs(i, j, 0) * a(i) + sq_devs(i, j, 1) * b(i));
            // if ((value < radius) & (value > 0)) {
            if (value < radius) {
                I(idx) = i + 1;
                J(idx) = j + 1;
                V(idx) = value;
                if (value > 0) {
                    ddista(idx) = sq_devs(i, j, 0) / value;
                    ddistb(idx) = sq_devs(i, j, 1) / value;
                } else {
                    ddista(idx) = 0;
                    ddistb(idx) = 0;

                }
                idx += 1;
            }
        }
    }

    // drop the extra data
    I.resize(idx);
    J.resize(idx);
    V.resize(idx);
    ddista.resize(idx);
    ddistb.resize(idx);

    return(Rcpp::List::create(Rcpp::Named("ind") = arma::join_horiz(I, J),
                      Rcpp::Named("V") = V,
                      Rcpp::Named("ddista") = ddista,
                      Rcpp::Named("ddistb") = ddistb));
}



