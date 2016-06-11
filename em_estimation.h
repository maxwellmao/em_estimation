#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>

class EMEstimation{
public:
    int K;

    std::vector<std::vector<double> > data;

    EMEstimation(const int &k)
    {
        this->K=k;
    }

    std::vector<std::vector<double> > posterior_of_z;

    std::vector<double> prior;

    std::vector<double> parameter;

    virtual void update_posterior() = 0;

    virtual double update_prior_and_parameter() = 0;

    virtual void init_model_parameter() = 0;

    void estimation(const int &MAX_ITER=500, const double &tot=10e-6);

    virtual void output_parameter() = 0;
};

class MixExponentialEMEstimation: public EMEstimation{
public:
    MixExponentialEMEstimation(const int &k): EMEstimation(k){}

    void update_posterior();

    double update_prior_and_parameter();

    void init_model_parameter();

    void output_parameter();
//    void estimation(const int &MAX_ITER=500, const double &tot=10e-6);
};

class MixBinomialEMEstimation: public EMEstimation{
public:
    int N;

    MixBinomialEMEstimation(const int &k, const int &n): EMEstimation(k)
    {
        this->N=n;
    }
    
    void update_posterior();

    double update_prior_and_parameter();

    void init_model_parameter();

    void output_parameter();
//    void estimation(const int &MAX_ITER=500, const double &tot=10e-6);
};

class MixPoissonEMEstimation: public MixBinomialEMEstimation{
public:
    MixPoissonEMEstimation(const int &k): MixBinomialEMEstimation(k, 1){}
    void update_posterior();
    
    void output_parameter();
};

class MixGaussianEMEstimation: public EMEstimation{
public:
    unsigned int dim;

    MixGaussianEMEstimation(const int &k, const unsigned int &d): EMEstimation(k), dim(d){}

    std::vector<std::vector<double> > parameter_mu;// mean values
    std::vector<std::vector<std::vector<double> > > parameter_inv_sigma;// covariance matrix
    std::vector<double> parameter_inv_det_sigma;// the inverse of the determinant of sigma

    void update_posterior();

    double update_prior_and_parameter();

    void init_model_parameter();

    void output_parameter();
};
