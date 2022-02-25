// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List distance_near_with_ddist_cpp(arma::mat& locs, arma::mat& locs_grid, double& radius, const int& n_neighbors=68) {
    int N = locs.n_rows;
    int N_grid = locs_grid.n_rows;

    int max_points = 2 * n_neighbors * N_grid;


    // initialize storage vector
    arma::vec I(max_points);
    arma::vec J(max_points);
    arma::vec V(max_points);
    arma::vec ddistx(max_points);
    arma::vec ddisty(max_points);

    // initialize the index
    int idx = 0;

    // distance loop
    for (int i=0; i<N; i++) {
        for (int j=0; j<N_grid; j++ ) {
            double value = sqrt(pow(locs(i, 0) - locs_grid(j, 0), 2) + pow(locs(i, 1) - locs_grid(j, 1), 2));
            if ((value <= radius) & (value > 0)) {
                I(idx) = i + 1;
                J(idx) = j + 1;
                // // fill in the symmetrical distance
                // I(idx+1) = j;
                // J(idx + 1) = i;
                V(idx) = value;
                ddistx(idx) = (locs(i, 0) - locs_grid(j, 0)) / value;
                ddisty(idx) = (locs(i, 1) - locs_grid(j, 1)) / value;
                // V(idx+1) = value;
                idx += 1;
                // idx += 2;
            }
        }
    }

    // drop the extra data
    I.resize(idx);
    J.resize(idx);
    V.resize(idx);
    ddistx.resize(idx);
    ddisty.resize(idx);

    return(Rcpp::List::create(Rcpp::Named("ind") = arma::join_horiz(I, J),
                      Rcpp::Named("V") = V,
                      Rcpp::Named("ddistx") = ddistx,
                      Rcpp::Named("ddisty") = ddisty));
}

