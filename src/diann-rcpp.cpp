#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_hello_world() {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Identity(3, 3);
    // Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);
    // Do not use Random() here to not promote use of a non-R RNG
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Zero(3, 3);
    for (auto i=0; i<m2.rows(); i++)
        for (auto j=0; j<m2.cols(); j++)
            m2(i,j) = R::rnorm(0, 1);

    return m1 + 3 * (m1 + m2);
}

// [[Rcpp::export]]
std::vector<double> col_max(std::vector<double> &quantities, int m, int n) {
    int i, k, shift;
    std::vector<double> res(n);
    for (i = shift = 0; i < n; i++, shift += m) {
        double best = -1000000.0;
        for (k = 0; k < m; k++) if (quantities[shift + k] > best) best = quantities[shift + k];
        res[i] = best;
    }
    return res;
}

// [[Rcpp::export]]
std::vector<double> maxlfq_solve(std::vector<double> &quantities, int peptides, int samples, double margin = -10.0001) {
    int i, j, k;
    
    std::vector<double> ratios, ref(samples);
    Eigen::MatrixXd A(samples, samples);
    Eigen::VectorXd B(samples);
    
    for (i = 0; i < samples; i++) {
        B(i) = 0.0;
        for (j = 0; j < samples; j++) A(i, j) = 0.0;
    }
    
    const double inf = 1000000.0;
    
    for (i = 0; i < samples; i++) {
        double max = -inf;
        double *pi = &quantities[i * peptides];
        for (k = 0; k < peptides; k++) if (pi[k] > max) max = pi[k];
        ref[i] = (max > margin) ? max : (-inf); 
    }
    
    for (i = 0; i < samples; i++) {
        for (j = i + 1; j < samples; j++) {
            ratios.clear();
            double *pi = &quantities[i * peptides], *pj = &quantities[j * peptides];
            for (k = 0; k < peptides; k++) if (pi[k] > margin && pj[k] > margin) ratios.push_back(pi[k] - pj[k]);
            
            if (ratios.size()) {
                double median;
                if (ratios.size() >= 2) {
                    std::sort(ratios.begin(), ratios.end());
                    if (ratios.size() & 1) median = ratios[ratios.size() / 2];
                    else median = 0.5 * (ratios[(ratios.size() / 2) - 1] + ratios[ratios.size() / 2]);
                } else median = ratios[0];
                
                A(i, i) += 1.0;
                A(j, j) += 1.0;
                A(i, j) = A(j, i) = -1.0;
                B(i) += median;
                B(j) -= median;
            }
        }
    }
    
    for (i = 0; i < samples; i++) {
        double reg = 0.0001 * (A(i, i) >= 1.0 ? A(i, i) : 1.0);
        A(i, i) += reg;
        B(i) += ref[i] * reg;
    }
    
    Eigen::VectorXd X = A.llt().solve(B);
    
    std::vector<double> x(samples);
    for (i = 0; i < samples; i++) x[i] = X(i);
    return x;
}
