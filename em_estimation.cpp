#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <math.h>
#include <iostream>
#include <time.h>
#include <assert.h>
#include "em_estimation.h"
#include "vec_mat_operation.h"
#include "common.h"

using namespace std;

double exponential_pdf(const double &x, const double &lambda)
{
    return lambda*exp(-lambda*x);
}

double log_choose(const int &base, const int &chose)
{
    double log_c=0.0;

    if(chose==0 || base==chose)
        return 0.0;

    assert(base>chose);

    int upper=chose;

    if(chose>base-chose)
        upper=base-chose;

    for(int i=1; i<=upper; i++)
    {
        log_c+=log(base-upper+i)-log(i);
    }

    return log_c;
}

double binomial_log_pdf(const double &x, const double &p, const double &n)
{
    assert(p!=0.0 && p!=1.0);
        
    return log_choose(n, x)+x*log(p)+(n-x)*log(1-p);
}

double binomial_pdf(const double &x, const double &p, const double &n)
{
    if(p==0.0)
    {
        return x==0.0? 1.0: 0.0;
    }
    if(p==1.0)
    {
        return x==n? 1.0: 0.0;
    }
        
    return exp(log_choose(n, x)+x*log(p)+(n-x)*log(1-p));
}

double poisson_pdf(const double &x, const double &lambda)
{
    double p=x*log(lambda)-lambda;

    for(int i=1; i<int(x); i++)
    {
        p-=log(i);
    }
    return exp(p);
}

double multivariate_gaussian_pdf(const vector<double> &x, const vector<double> &mu, const vector<vector<double> > &inv_sigma, const double &inv_det_sigma)
{
    vector<double> tmp(x);
    
//    cout << "x: ";
//    output_vec(tmp);
//    cout << "mu: ";
//    output_vec(mu);

    vec_sub(tmp, mu);

    vector<double> first_mult;

    vec_mult_mat(tmp, inv_sigma, first_mult);

//    cout << "Inv det: " << inv_det_sigma << endl;

    double p=exp((-0.5)*vec_dot_product(first_mult, tmp))*sqrt(inv_det_sigma)/(pow(sqrt(2*M_PI), x.size()));

//    cout << "p=" << p << endl;

    return p;
}

void MixExponentialEMEstimation::update_posterior()
{
    // expectation
    for(size_t d_index=0; d_index<this->data.size(); d_index++)
    {
        double cond=0.0;
        for(int index=0; index<this->K; index++)
        {
            cond+=this->prior[index]*exponential_pdf(this->data[d_index][0], this->parameter[index]);
        }
        assert(cond>0);
        for(int index=0; index<this->K; index++)
        {
            posterior_of_z[d_index][index]=this->prior[index]*exponential_pdf(this->data[d_index][0], this->parameter[index])/cond;
        }
    }
}

double MixExponentialEMEstimation::update_prior_and_parameter()
{
    // maximization
    double error=0.0;

    for(int index=0; index<this->K; index++)
    {
        double posterior_sum=0.0;
        
        double rate=0.0;

        for(size_t d_index=0; d_index<this->data.size(); d_index++)
        {
            posterior_sum+=this->posterior_of_z[d_index][index];
            rate+=this->data[d_index][0]*this->posterior_of_z[d_index][index];
        }

        this->prior[index]=posterior_sum/this->data.size();

        assert(rate>0);

        rate=posterior_sum/rate;

        error+=fabs(rate-this->parameter[index]);

//        cout << rate << "-" << this->parameter[index] << "=" << error << "; ";
        
        this->parameter[index]=rate;
    }
    cout << error << endl;
    return error;
}

void MixExponentialEMEstimation::init_model_parameter()
{
    this->prior.clear();

    this->parameter.clear();

    this->posterior_of_z.clear();

    this->prior.resize(this->K, 1.0/this->K);
    
    this->parameter.resize(this->K, 0.0);

    for(int index=0; index<K; index++)
    {
        this->parameter[index]=1.0/this->data[rand() % this->data.size()][0];
    }

    vector<double> tmp(this->K, 1.0/this->K);

    this->posterior_of_z.resize(this->data.size(), tmp);
}

void MixExponentialEMEstimation::output_parameter()
{
}

void EMEstimation::estimation(const int &MAX_ITER, const double &tot)
{
    int iter=0;
    double error=1.0;

//    cout << "Init Parameter: ";
//    this->output_parameter();
//    output_list(this->parameter);

    while(iter<MAX_ITER && error >= tot)
    {
        this->update_posterior();

        error=this->update_prior_and_parameter();

        iter++;

        cout << "Iter: " << iter << ", error: " << error << endl;

        this->output_parameter();
    }
    cout << "Final Prior: ";
    output_list(this->prior);
    
//    this->output_parameter();
    cout << "Final Parameter: ";
    output_list(this->parameter);
}

void MixBinomialEMEstimation::init_model_parameter()
{
    this->prior.clear();

    this->parameter.clear();

    this->posterior_of_z.clear();

    this->prior.resize(this->K, 1.0/this->K);
    
    this->parameter.resize(this->K, 0.0);

    for(int index=0; index<K; index++)
    {
        // make sure no two initialized parameters are same
        double new_para=this->data[rand() % this->data.size()][0]/this->N;

        while(true)
        {
            int j=index-1;
            for(; j>=0; j--)
            {
                if(new_para==this->parameter[j])
                    break;
            }
            if(j<0)
                break;
            
            new_para=this->data[rand() % this->data.size()][0]/this->N;
        }

        this->parameter[index]=new_para;
    }

    vector<double> tmp(this->K, 1.0/this->K);

    this->posterior_of_z.resize(this->data.size(), tmp);
}

void MixBinomialEMEstimation::update_posterior()
{
    // expectation
    for(size_t d_index=0; d_index<this->data.size(); d_index++)
    {
        double cond=0.0;

        for(int index=0; index<this->K; index++)
        {
            double p=binomial_pdf(this->data[d_index][0], this->parameter[index], this->N);
            cond += this->prior[index]*p;
        }
        if(cond>0)
        {
            for(int index=0; index<this->K; index++)
            {
                posterior_of_z[d_index][index]=this->prior[index]*binomial_pdf(this->data[d_index][0], this->parameter[index], this->N)/cond;
            }
        }
        else
        {
            double log_min=1.0; // log-cond

            vector<double> log_prob(this->K, 0.0);

            for(int index=0; index<this->K; index++)
            {
                log_prob[index]=log(this->prior[index])+binomial_log_pdf(this->data[d_index][0], this->parameter[index], this->N);
                if(log_min>log_prob[index])
                    log_min=log_prob[index];
            }
            cond=0.0;

            for(int index=0; index<this->K; index++)
            {
                cond+=exp(log_prob[index]-log_min);
            }

            assert(cond > 0.0);
            
            for(int index=0; index < this->K; index++)
            {
                posterior_of_z[d_index][index]=exp(log_prob[index]-log_min-log(cond));
            }
        }
    }
}

double MixBinomialEMEstimation::update_prior_and_parameter()
{
    // maximization
    double error=0.0;
    
    for(int index=0; index<this->K; index++)
    {
        double posterior_sum=0.0;
        
        double rate=0.0;

        for(size_t d_index=0; d_index<this->data.size(); d_index++)
        {
//            cout << this->posterior_of_z[d_index][index] << endl;

            posterior_sum+=this->posterior_of_z[d_index][index];
            rate+=this->data[d_index][0]*this->posterior_of_z[d_index][index];
        }

        this->prior[index]=posterior_sum/this->data.size();

        assert(rate>0.0);

        rate=rate/posterior_sum/this->N;

        error+=fabs(rate-this->parameter[index]);
        
        this->parameter[index]=rate;
    }
    
    return error;
}

void MixBinomialEMEstimation::output_parameter()
{
    cout << "Parameters:" << endl;

    output_vec(this->parameter);

    cout << "Priors:" << endl;

    output_vec(this->prior);
}

void MixPoissonEMEstimation::update_posterior()
{
    // expectation
    for(size_t d_index=0; d_index<this->data.size(); d_index++)
    {
        double cond=0.0;
        for(int index=0; index<this->K; index++)
        {
            cond+=this->prior[index]*poisson_pdf(this->data[d_index][0], this->parameter[index]);
        }
        assert(cond>0);
        for(int index=0; index<this->K; index++)
        {
            posterior_of_z[d_index][index]=this->prior[index]*poisson_pdf(this->data[d_index][0], this->parameter[index])/cond;
        }
    }
}

void MixPoissonEMEstimation::output_parameter()
{
}

void MixGaussianEMEstimation::update_posterior()
{
    // expectation
    for(size_t d_index=0; d_index<this->data.size(); d_index++)
    {
        double cond=0.0;
        for(int index=0; index<this->K; index++)
        {
            cond+=this->prior[index]*multivariate_gaussian_pdf(this->data[d_index], this->parameter_mu[index], this->parameter_inv_sigma[index], this->parameter_inv_det_sigma[index]);
        }
        assert(cond>0);
        for(int index=0; index<this->K; index++)
        {
            posterior_of_z[d_index][index]=this->prior[index]*multivariate_gaussian_pdf(this->data[d_index], this->parameter_mu[index], this->parameter_inv_sigma[index], this->parameter_inv_det_sigma[index])/cond;
        }
//        output_vec(posterior_of_z[d_index]);
//        cout << vec_fabs_sum(posterior_of_z[d_index]) << endl;;
    }
}

double MixGaussianEMEstimation::update_prior_and_parameter()
{
    // maximization
    double error=0.0;
    
    for(int index=0; index<this->K; index++)
    {
        double posterior_sum=0.0;

        vector<double> mu(this->dim, 0.0);

        vector<vector<double> > sigma(this->dim, mu);

        for(size_t d_index=0; d_index<this->data.size(); d_index++)
        {
            posterior_sum+=this->posterior_of_z[d_index][index];

            vector<double> d_i(this->data[d_index]);

            vec_dot_mult(d_i, this->posterior_of_z[d_index][index]);

            vec_add(mu, d_i);

        }
        this->prior[index]=posterior_sum/this->data.size();

        vec_dot_div(mu, posterior_sum);

        vec_sub(this->parameter_mu[index], mu);

        error+=vec_fabs_sum(this->parameter_mu[index]);

        copy(mu.begin(), mu.end(), this->parameter_mu[index].begin());

        for(size_t d_index=0; d_index<this->data.size(); d_index++)
        {
            vector<double> diff(this->data[d_index]);

            vec_sub(diff, mu);

            vector<vector<double> > sigma_i;
            
            vec_mult_to_mat(diff, diff, sigma_i);

            mat_dot_mult(sigma_i, this->posterior_of_z[d_index][index]);

//            cout << "Sigma" << endl;
//            output_mat(sigma);
//            cout << "Sigma_i" << endl;
//            output_mat(sigma_i);

            mat_add(sigma, sigma_i);
        }
        mat_dot_div(sigma, posterior_sum);

        vector<vector<double> > inv_sigma;

//        cout << "Sigma" << endl;
//        output_mat(sigma);

        this->parameter_inv_det_sigma[index]=fabs(1.0/mat_inverse(sigma, inv_sigma));

//        cout << "Inverse Sigma" << endl;
//        output_mat(inv_sigma);

        mat_sub(this->parameter_inv_sigma[index], inv_sigma);

        error+=mat_fabs_sum(this->parameter_inv_sigma[index]);

        copy(inv_sigma.begin(), inv_sigma.end(), this->parameter_inv_sigma[index].begin());
    }
    return error;
}

void MixGaussianEMEstimation::init_model_parameter()
{
    this->prior.clear();

    this->parameter_mu.clear();

    this->posterior_of_z.clear();

    this->parameter_inv_sigma.clear();

    this->parameter_inv_det_sigma.clear();

    this->prior.resize(this->K, 1.0/this->K);

    this->parameter_mu.reserve(this->K);

    this->parameter_inv_sigma.reserve(this->K);

    vector<double> data_min(this->data[0]);
    vector<double> data_max(this->data[0]);

    for(size_t d=1; d<this->data.size(); d++)
    {
        for(size_t i=0; i<this->dim; i++)
        {
            if(data_min[i]>this->data[d][i])
                data_min[i]=this->data[d][i];

            if(data_max[i]<this->data[d][i])
                data_max[i]=this->data[d][i];
        }
    }

    for(int index=0; index<K; index++)
    {
//        this->parameter_mu.push_back(this->data[rand() % this->data.size()]);
        vector<double> mu(this->dim, 0.0);

        for(size_t d=0; d<this->dim; d++)
        {
            mu[d]=data_min[d]+(data_max[d]-data_min[d])*(double)index/this->K;
        }
        this->parameter_mu.push_back(mu);

//        vector<double> row(this->dim, 0.0);
//
//        vector<vector<double> > sigma(this->dim, row);
//
//        for(unsigned int d=0; d<this->dim; d++)
//        {
//            sigma[d][d]=(double)rand()/RAND_MAX;
//        }
//
//        this->parameter_inv_sigma.push_back(sigma);
//    
//        this->parameter_inv_det_sigma.push_back(fabs(mat_determinant(this->parameter_inv_sigma[index])));
    }

    vector<double> tmp(this->K, 1.0/this->K);

    this->posterior_of_z.resize(this->data.size(), tmp);

    this->parameter_inv_det_sigma.resize(this->K, 0.0);

    this->parameter_inv_sigma.reserve(this->K);

    for(int index=0; index<this->K; index++)
    {
        vector<double> sigma_row(this->dim, 0.0);

        vector<vector<double> > sigma(this->dim, sigma_row);

        double posterior_sum=0.0;

        for(size_t d_index=0; d_index<this->data.size(); d_index++)
        {
            posterior_sum+=this->posterior_of_z[d_index][index];

            vector<double> diff(this->data[d_index]);
    
            vec_sub(diff, this->parameter_mu[index]);
    
            vector<vector<double> > sigma_i;
            
            vec_mult_to_mat(diff, diff, sigma_i);
    
            mat_dot_mult(sigma_i, this->posterior_of_z[d_index][index]);
    
            mat_add(sigma, sigma_i);
        }
        mat_dot_div(sigma, posterior_sum);

        vector<vector<double> > inv_sigma;
        
        this->parameter_inv_det_sigma[index]=fabs(1.0/mat_inverse(sigma, inv_sigma));

        this->parameter_inv_sigma.push_back(inv_sigma);
        
    }
    cout << "Init Parameter: ";
    this->output_parameter();
}

void MixGaussianEMEstimation::output_parameter()
{
    cout << "Prior: ";
    output_vec(this->prior);
    for(int k=0; k < this->K; k++)
    {
        cout << "Cluster " << k << ":" << endl;
        cout << "mu: ";
        output_vec(this->parameter_mu[k]);
        cout << "inv sigma:" << endl;
        output_mat(this->parameter_inv_sigma[k]);
        cout << "inv sigma determinant: " << this->parameter_inv_det_sigma[k] << endl;
        cout << "===============================" << endl;
    }
}

void em_estimation_MoB(const int &K, const unsigned int &N)
{
    MixBinomialEMEstimation model(K, N);

    string line;

    while(getline(cin, line))
    {
        float d=atof(line.c_str());

        vector<double> d_i(1, d);

        model.data.push_back(d_i);
    }
    model.init_model_parameter();

    cerr << "Total number of data points: " << model.data.size() << endl;

    model.estimation();
}

void em_estimation_MoP(const int &K)
{
    MixPoissonEMEstimation model(K);

    string line;

    while(getline(cin, line))
    {
        float d=atof(line.c_str());

        vector<double> d_i(1, d);

        model.data.push_back(d_i);
    }
    model.init_model_parameter();

    cout << "Total number of data points: " << model.data.size() << endl;

    model.estimation();
}

void em_estimation_MoE(const int &K)
{
    MixExponentialEMEstimation model(K);

    string line;

    while(getline(cin, line))
    {
        float d=atof(line.c_str());

        vector<double> d_i(1, d);

        model.data.push_back(d_i);
    }
    model.init_model_parameter();

    cout << "Total number of data points: " << model.data.size() << endl;

    model.estimation();
}

void em_estimation_MoG(const int &K, const unsigned int &dim, const bool &with_label_at_first, const string &save_path)
{
    MixGaussianEMEstimation model(K, dim);

    string line;

    vector<int> data_label;

    while(getline(cin, line))
    {
        vector<string> content;

        string_split(line, "\t", content);

        if((with_label_at_first && content.size()!=dim+1) || (!with_label_at_first && content.size()!=dim))
        {
            string_split(line, " ", content);
        }
        assert((with_label_at_first && content.size()==dim+1) || (!with_label_at_first && content.size()==dim));

        vector<double> d(dim, 0.0);

        if(with_label_at_first)
        {
            data_label.push_back(atoi(content[0].c_str()));
            
            for(size_t i=1; i<content.size(); i++)
            {
                d[i-1]=atof(content[i].c_str());
            }
        }
        else
        {
            for(size_t i=0; i<content.size(); i++)
            {
                d[i]=atof(content[i].c_str());
            }
        }
        model.data.push_back(d);
    }
    cout << "Total number of data points: " << model.data.size() << ", dimension: " << dim << " " << model.data[0].size() << endl;

    model.init_model_parameter();

//    model.update_prior_and_parameter();
    
    model.estimation();

    for(int i=0; i<K; i++)
    {
        cout << "Cluster " << i+1 << endl;
        cout << "Mu: " << endl;
        output_vec(model.parameter_mu[i]);
        cout << "Sigma: " << endl;
        vector<vector<double> > sigma;
        mat_inverse(model.parameter_inv_sigma[i], sigma);
        output_mat(sigma);
    }

    if(save_path.size()>0)
    {
        ofstream fp(save_path.c_str());
        if(!fp)
        {
            cerr << "Cannot save into " << save_path << endl;
            return;
        }
        cerr << "Saving into " << save_path << endl;
        for(size_t index=0; index<model.data.size(); index++)
        {
            fp << data_label[index] << " ";

            for(int k=0; k<K; k++)
            {
                fp << model.posterior_of_z[index][k] << " ";
            }
            fp << endl;
        }
        fp.close();
    }
}

void test_mat_assign()
{
    vector<double> tmp(3, 1.0);

    vector<double> tmp_vec(tmp);

    output_vec(tmp);

    output_vec(tmp_vec);

    tmp_vec[1]=3;

    output_vec(tmp);

    output_vec(tmp_vec);

    vector<vector<double> > mat(3, tmp);
    
    cout << "Original" << endl;

    output_mat(mat);

    vector<vector<double> > tmp_mat;

    tmp_mat.assign(mat.begin(), mat.end());
    
    cout << "Copy" << endl;
    output_mat(tmp_mat);

    mat[0][0]=2;
    tmp_mat[1][1]=3;
    
    cout << "Original" << endl;
    output_mat(mat);
    
    cout << "Copy" << endl;
    output_mat(tmp_mat);
}

void test_mat_inverse()
{
    int dim=0;
    cin >> dim;

    vector<vector<double> > mat;
    for(int i=0; i<dim; i++)
    {
        vector<double> row(dim, 0.0);

        for(int j=0; j<dim; j++)
        {
            double e;
            cin >> e;
            row[j]=e;
        }
        mat.push_back(row);
    }
    vector<vector<double> > inv_mat;

    cout << "Determinant: " << endl;

    cout << mat_inverse(mat, inv_mat) << endl;

    cout << "Inverse matrix: " << endl;

    output_mat(inv_mat);
}

enum MixModel{MixExp, MixPoiss, MixBin, MixGauss, Test};

unsigned int dim=1;

int K=2;

bool with_label_at_first=false;

MixModel model_type=Test;

string save_path="";

void parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++){
        if(argv[i][0] == '-') {
            switch(argv[i][1]){
                case 'g':
                    model_type=MixGauss;
                    dim=atoi(argv[i+1]);
                    i++;
                    break;

                case 'b':
                    model_type=MixBin;
                    dim=atoi(argv[i+1]);
                    i++;
                    break;

                case 'p':
                    model_type=MixPoiss;
                    break;

                case 'e':
                    model_type=MixExp;
                    break;

                case 't':
                    model_type=Test;
                    break;

                case 'k':
                    K=atoi(argv[i+1]);
                    i++;
                    break;

                case 'l':
                    with_label_at_first=true;
                    break;

                case 's':
                    save_path=argv[i+1];
                    i++;
                    break;

                default:
                    break;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    unsigned int seed=time(NULL);

    cerr << "Random seed: " << seed << endl;

    srand(seed);

    parse_args(argc, argv);

    switch(model_type)
    {
        case Test:
            test_mat_inverse();
            test_mat_assign();
            break;

        case MixGauss:
            em_estimation_MoG(K, dim, with_label_at_first, save_path);
            break;

        case MixExp:
            em_estimation_MoE(K);
            break;

        case MixPoiss:
            em_estimation_MoP(K);
            break;

        case MixBin:
            em_estimation_MoB(K, dim);
            break;
    }

    return 1;
}

