#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Calculate thresheld pairwise distance for a row
//'
//' @param i The row index
//' @param locs An N x 2 matrix of spatial locations
//' @param locs_grid An N_grid x 2 matrix of spatial grids
//' @param radius The thresholded radius
//' @param byrow Perform calculation row/column-wise
//' @return The thresheld pairwise distance
//'
//' @export
//[[Rcpp::export]]
arma::mat distance_near_row_cpp(const double& i, const arma::rowvec& locs, const arma::mat& locs_grid, const double& radius, const bool& byrow=true) {
    // Rcpp::List distance_near_row_cpp(const double& i, const arma::mat& locs, const arma::mat& locs_grid, const double& radius) {
    int N_grid = locs_grid.n_rows;

    // initialize storage vector

    arma::mat out;
    if (byrow) {
        out.resize(5, N_grid);
    } else {
        out.resize(N_grid, 5);
    }

    // initialize the index
    int idx = 0;
    double devs_x;
    double devs_y;
    double value;

    // distance loop
    for (int j=0; j<N_grid; j++ ) {
        devs_x = locs(0) - locs_grid(j, 0);
        devs_y = locs(1) - locs_grid(j, 1);
        value = sqrt(pow(devs_x, 2) + pow(devs_y, 2));
        if ((value < radius) & (value > 0)) {
            if (byrow) {
                out(0, idx) = i + 1;
                out(1, idx) = j + 1;
                out(2, idx) = value;
                out(3, idx) = (devs_x) / value;
                out(4, idx) = (devs_y) / value;
            } else {
                out(idx, 0) = i + 1;
                out(idx, 1) = j + 1;
                out(idx, 2) = value;
                out(idx, 3) = (devs_x) / value;
                out(idx, 4) = (devs_y) / value;

            }
            idx += 1;
        }
    }

    // drop the extra data
    if (byrow) {
        out.resize(5, idx);
    } else{
        out.resize(idx, 5);
    }

    return(out);
}


//' Calculate thresheld pairwise distance for a row using loops for parallelization
//'
//' @param locs An N x 2 matrix of spatial locations
//' @param locs_grid An N_grid x 2 matrix of spatial grids
//' @param radius The thresholded radius
//' @param n_neighbors The expected number of neighbors based on the MRA grid
//' @return The thresheld pairwise distance
//' @param byrow Perform calculation row/column-wise
//'
//' @export
//[[Rcpp::export]]
arma::mat distance_near_loop_cpp(const arma::mat& locs, const arma::mat& locs_grid, const double& radius, const int& n_neighbors=86, const bool& byrow=true) {
    // Rcpp::List distance_near_loop_cpp(const arma::mat& locs, const arma::mat& locs_grid, const double& radius, const int& n_neighbors=86) {

    // https://stackoverflow.com/questions/26514691/dynamically-readjustable-arrays-and-openmp
    // https://stackoverflow.com/questions/59338272/dynamically-add-rows-in-rcpp
    // https://stackoverflow.com/questions/26514691/dynamically-readjustable-arrays-and-openmp
    // https://stackoverflow.com/questions/70374211/rcppparallel-equivalent-of-combine-rbind
    // https://stackoverflow.com/questions/31913437/r-fast-cbind-matrix-using-rcpp
    // https://stackoverflow.com/questions/62534741/produce-sparse-pairwise-distance-matrix-python-avoiding-memory-error
    // https://stackoverflow.com/questions/26001319/efficient-algorithm-to-calculate-pairwise-distances-for-a-large-set-of-locations
    // https://gallery.rcpp.org/articles/parallel-distance-matrix/
    // https://stackoverflow.com/questions/69648919/constructing-distance-matrix-in-parallel-in-c11-using-openmp
    // https://stackoverflow.com/questions/26001319/efficient-algorithm-to-calculate-pairwise-distances-for-a-large-set-of-locations

    int N = locs.n_rows;

    int max_points = 2 * n_neighbors * N;

    // initialize the output list
    arma::mat out;
    if (byrow) {
        out.resize(5, max_points);
    } else {
        out.resize(max_points, 5);
    }

    // run the loop
    // QUESTION: is it faster to subset locs here or inside the distance_near_row_cpp function?

    int idx = 0;
    for (int i=0; i<N; i++) {
        // Add user interrupt check?
        arma::mat tmp = distance_near_row_cpp(i, locs.row(i), locs_grid, radius, byrow);
        int n_tmp;
        if (byrow) {
            n_tmp=tmp.n_cols;
        } else {
            n_tmp=tmp.n_rows;
        }
        // for (int j=0; j<n_tmp; j++) {
        //     if (byrow) {
        //         out.col(idx) = tmp.col(j);
        //     } else {
        //         out.row(idx) = tmp.row(j);
        //     }
        //     idx++;
        // }
        if (byrow) {
            out.cols(idx, idx+n_tmp-1) = tmp;
        } else {
            out.rows(idx, idx+n_tmp-1) = tmp;
        }
        idx += n_tmp;
    }

    // Resize the result
    // figure this out later or do in R
    if (byrow) {
        out.resize(5, idx);
    } else {
        out.resize(idx, 5);
    }

    // return the result
    return(out);
}



//' Calculate thresheld pairwise distance for a row using loops for parallelization
//'
//' @param locs An N x 2 matrix of spatial locations
//' @param locs_grid An N_grid x 2 matrix of spatial grids
//' @param radius The thresholded radius
//' @param n_neighbors The expected number of neighbors based on the MRA grid
//' @param byrow Perform calculation row/column-wise
//' @param nchunks The number chunks to divide the data into
//' @param ncores The number of openmp threads
//' @return The thresheld pairwise distance
//'
//' @export
//[[Rcpp::export]]
arma::field<arma::mat> distance_near_chunk_cpp(const arma::mat& locs,
                                       const arma::mat& locs_grid,
                                       const double& radius,
                                       const int& n_neighbors=86,
                                       const bool& byrow=true,
                                       const bool& joint_index=true,
                                       Rcpp::Nullable<int> nchunks=R_NilValue,
                                       const int& ncores=1) {
        // Rcpp::List distance_near_chunk_cpp(const arma::mat& locs,
        //                            const arma::mat& locs_grid,
        //                            const double& radius,
        //                            const int& n_neighbors=86,
        //                            const bool& byrow=true,
        //                            const bool& joint_index=true,
        //                            Rcpp::Nullable<int> nchunks=R_NilValue,
        //                            const int& ncores=1) {
    // Rcpp::List distance_near_loop_cpp(const arma::mat& locs, const arma::mat& locs_grid, const double& radius, const int& n_neighbors=86) {

    // https://stackoverflow.com/questions/26514691/dynamically-readjustable-arrays-and-openmp
    // https://stackoverflow.com/questions/59338272/dynamically-add-rows-in-rcpp
    // https://stackoverflow.com/questions/26514691/dynamically-readjustable-arrays-and-openmp
    // https://stackoverflow.com/questions/70374211/rcppparallel-equivalent-of-combine-rbind
    // https://stackoverflow.com/questions/31913437/r-fast-cbind-matrix-using-rcpp
    // https://stackoverflow.com/questions/62534741/produce-sparse-pairwise-distance-matrix-python-avoiding-memory-error
    // https://stackoverflow.com/questions/26001319/efficient-algorithm-to-calculate-pairwise-distances-for-a-large-set-of-locations
    // https://gallery.rcpp.org/articles/parallel-distance-matrix/
    // https://stackoverflow.com/questions/69648919/constructing-distance-matrix-in-parallel-in-c11-using-openmp
    // https://stackoverflow.com/questions/26001319/efficient-algorithm-to-calculate-pairwise-distances-for-a-large-set-of-locations


    // define some code chunks
    int N = locs.n_rows;
    int n_chunks;
    if (nchunks.isNotNull()) {
        n_chunks = as<int>(nchunks);
    } else {
        n_chunks = ncores;
    }

    if (n_chunks < ncores) {
        // make there be at least as many chunks as cores
        n_chunks = ncores;
    }

    arma::vec chunk_size(n_chunks, arma::fill::ones);
    chunk_size *= N / n_chunks;
    // Rcpp::Rcout << "Chunk size = " << chunk_size << "\n";
    // Rcpp::Rcout << "N % n_chunks = " <<  N % n_chunks << "\n";
    if (N % n_chunks != 0) {
        // Rcpp::Rcout << "Chunk subvec" << chunk_size.subvec(0, N % n_chunks - 1) << "\n";
        chunk_size.subvec(0, N % n_chunks - 1) += 1;
    }
    // Rcpp::Rcout << "Chunk size = " << chunk_size << "\n";
    // Rcpp::Rcout << "total Chunk size = " << sum(chunk_size) << "\n";

    int max_points_loop = 2 * n_neighbors * N / n_chunks;

    // Rcpp::List out_list(n_chunks);
    // using arma::field rather than Rcpp::List https://stackoverflow.com/questions/29195019/subset-armadillo-field
    arma::field<arma::mat> out_list(n_chunks, 1);

    for (int j=0; j<n_chunks; j++) {
        if (byrow) {
            out_list(j) = arma::mat(5, max_points_loop);
        } else {
            out_list(j) = arma::mat(max_points_loop, 5);
        }
    }

#if defined(_OPENMP)
#pragma omp parallel for num_threads(ncores) if (ncores > 1)
#endif
    for (int j=0; j<n_chunks; j++) {
        // arma::mat out;
        // if (byrow) {
        //     out.resize(5, max_points_loop);
        // } else {
        //     out.resize(max_points_loop, 5);
        // }


        int start_index;
        if (j == 0) {
            start_index = 0;
        } else {
            start_index = sum(chunk_size.subvec(0, j-1));
        }
        // Rcpp::Rcout << "Sum Chunk size for j = " << j << " is " << start_index << "\n";
        // Rcpp::Rcout << "Chunk size for j = "<< j << " is " << chunk_size(j) << "\n";

        int idx = 0;
        for (int i=start_index; i<start_index + chunk_size(j); i++) {
            // Add user interrupt check?
            arma::mat tmp = distance_near_row_cpp(i, locs.row(i), locs_grid, radius, byrow);
            int n_tmp;
            if (byrow) {
                n_tmp=tmp.n_cols;
            } else {
                n_tmp=tmp.n_rows;
            }
            if (joint_index) {
                if (byrow) {
                    out_list(j).cols(idx, idx+n_tmp-1) = tmp;
                } else {
                    out_list(j).rows(idx, idx+n_tmp-1) = tmp;
                }
                idx += n_tmp;
            } else {
                for (int j=0; j<n_tmp; j++) {
                    if (byrow) {
                        out_list(j).col(idx) = tmp.col(j);
                    } else {
                        out_list(j).row(idx) = tmp.row(j);
                    }
                    idx++;
                }
            }
        }
        // Resize the result
        if (byrow) {
            out_list(j).resize(5, idx);
        } else {
            out_list(j).resize(idx, 5);
        }

        // save to a list -- try and make this more efficient
        // out_list(j) = out;
    }

    // return the result
    return(out_list);
}





