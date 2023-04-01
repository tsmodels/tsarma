/// @file armafun.hpp
#ifndef armafun_hpp
#define armafun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type armafun(objective_function<Type>* obj) {
    DATA_VECTOR(y);
    PARAMETER(mu);
    PARAMETER_VECTOR(phi);
    PARAMETER_VECTOR(theta);
    PARAMETER_VECTOR(xi);
    PARAMETER(sigma);
    PARAMETER_VECTOR(distribution);
    // parameter scaling vector
    DATA_VECTOR(pscale);
    DATA_MATRIX(x);
    // model flags
    DATA_IVECTOR(cmodel);
    const int timesteps = y.rows();
    vector<Type> epsilon(timesteps);
    vector<Type> regressors(timesteps);
    vector<Type> fitted(timesteps);
    vector<Type> std_residuals(timesteps);
    regressors.setZero();
    fitted.setZero();
    epsilon.setZero();
    std_residuals.setZero();
    int m = x.cols();
    int j = 0;
    // re-scale parameters
    int k = 0;
    mu *= pscale(k);
    k += 1;
    for(j = 0;j<cmodel(1);j++) {
        phi(j) *= pscale(j + k);
    }
    if (cmodel(1) == 0) {
        k += 1;
    } else{
        k += cmodel(1);
    }
    for(j = 0;j<cmodel(2);j++){
        theta(j) *= pscale(j + k);
    }
    if (cmodel(2) == 0) {
        k += 1;
    } else{
        k += cmodel(2);
    }
    for(j = 0;j<m;j++){
        xi(j) *= pscale(j + k);
    }
    k += m;
    sigma *= pscale(k);
    k += 1;
    distribution(0) *= pscale(k);
    distribution(1) *= pscale(k + 1);
    distribution(2) *= pscale(k + 2);
    regressors = x * xi;
    vector<Type> const_intercept(timesteps);
    const_intercept.setZero();
    const_intercept.fill(mu);
    const_intercept.array() = const_intercept.array() + regressors.array();
    // initialize y to the mean
    int i = 0;
    Type sumorder = cmodel(1) + cmodel(2);
    for(i = 0;i<timesteps;i++){
        fitted(i) += const_intercept(i);
        if(sumorder > 0) {
            if (i >= cmodel(1)) {
                if(cmodel(1) > 0) {
                    for(j = 0;j<cmodel(1);j++){
                        fitted(i) += phi(j) * (y(i - j - 1) - const_intercept(i - j - 1));
                    }
                }
                if(cmodel(2) > 0) {
                    for(j = 0;j<cmodel(2);j++){
                        if ((i - j - 1) >= 0) {
                            fitted(i) += theta(j) * (y(i - j - 1) - fitted(i - j - 1));
                        }
                    }
                }
            }
        }
        epsilon(i) = y(i) - fitted(i);
        std_residuals(i) = epsilon(i)/sigma;
    }
    vector<Type> ll_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), cmodel(3))/sigma;
    REPORT(epsilon);
    REPORT(fitted);
    REPORT(const_intercept);
    REPORT(std_residuals);
    REPORT(ll_vector);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
