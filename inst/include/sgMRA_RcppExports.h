// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_sgMRA_RCPPEXPORTS_H_GEN_
#define RCPP_sgMRA_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace sgMRA {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("sgMRA", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("sgMRA", "_sgMRA_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in sgMRA");
            }
        }
    }

    inline Rcpp::List distance_near_with_ddist_cpp(arma::mat& locs, arma::mat& locs_grid, double& radius, const int& n_neighbors = 68) {
        typedef SEXP(*Ptr_distance_near_with_ddist_cpp)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_distance_near_with_ddist_cpp p_distance_near_with_ddist_cpp = NULL;
        if (p_distance_near_with_ddist_cpp == NULL) {
            validateSignature("Rcpp::List(*distance_near_with_ddist_cpp)(arma::mat&,arma::mat&,double&,const int&)");
            p_distance_near_with_ddist_cpp = (Ptr_distance_near_with_ddist_cpp)R_GetCCallable("sgMRA", "_sgMRA_distance_near_with_ddist_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_distance_near_with_ddist_cpp(Shield<SEXP>(Rcpp::wrap(locs)), Shield<SEXP>(Rcpp::wrap(locs_grid)), Shield<SEXP>(Rcpp::wrap(radius)), Shield<SEXP>(Rcpp::wrap(n_neighbors)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

}

#endif // RCPP_sgMRA_RCPPEXPORTS_H_GEN_