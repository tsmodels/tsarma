#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.armafilter)]]
NumericVector armafilter(NumericVector y, NumericVector epsilon, NumericVector x,
                         NumericVector initstate, NumericVector mu, NumericVector phi,
                         NumericVector theta, IntegerVector model) {
    int timesteps = epsilon.size();
    NumericVector fitted(timesteps);
    NumericVector constant(timesteps);
    for(int j=0;j<model(0);j++){
        fitted(j) = initstate(j);
        constant(j) = mu(0) + x(j);
    }
    for(int i = model(0);i<timesteps;i++) {
        constant(i) = mu(0) + x(i);
        fitted(i) = constant(i);
        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                fitted(i) += phi(j) * (y(i - j - 1) - constant(i - j - 1));
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                fitted(i) += theta(j) * epsilon(i - j - 1);
            }
        }
        epsilon(i) = y(i) - fitted(i);
    }
    return fitted;
}

// [[Rcpp::export(.armafilter2)]]
NumericVector armafilter2(NumericVector y, NumericVector epsilon, NumericVector x,
                         NumericVector mu, NumericVector phi,
                         NumericVector theta, IntegerVector model) {
    int timesteps = epsilon.size();
    NumericVector fitted(timesteps);
    NumericVector constant(timesteps);
    const int sumorder = model(1) + model(2);
    for(int i = 0;i<timesteps;i++) {
        constant(i) = mu(0) + x(i);
        fitted(i) = constant(i);
        if (sumorder > 0) {
            if (i >= model(1)) {
                if (model(1) > 0) {
                    for(int j = 0;j<model(1);j++) {
                        fitted(i) += phi(j) * (y(i - j - 1) - constant(i - j - 1));
                    }
                }
                if (model(2) > 0) {
                    for(int j = 0;j<model(2);j++) {
                        if ((i - j - 1) >= 0) {
                            fitted(i) += theta(j) * epsilon(i - j - 1);
                        }
                    }
                }
            }
        }
        epsilon(i) = y(i) - fitted(i);
    }
    return fitted;
}
