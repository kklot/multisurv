#include <TMB.hpp>

#define AGE_MAX 50
#define N_AGE 51
#define N_D 51 // differences, more than needed
#define N_PAR 9
#define N_Q 7

using Eigen::seqN;

template <class T> 
struct Kube {
    vector<T> 
        masterQ = vector<T>(N_AGE * N_Q * N_Q),
        masterP = vector<T>(N_AGE * N_Q * N_Q),
        masterD = vector<T>(N_AGE * N_Q * N_Q * N_D);
    matrix<T> 
        qM = matrix<T>(N_Q, N_Q),
        pM = matrix<T>(N_Q, N_Q), 
        pD = matrix<T>(N_Q, N_Q), 
        est = matrix<T>(N_AGE, N_PAR);
    int len = N_Q * N_Q;
    Kube () {};
    Kube (matrix<T> modelmatrix, vector<T> betas)
    {
        masterP.setZero();
        masterQ.setZero();
        masterD.setZero();
        for (int i = 0; i < N_PAR; i++)
            est.col(i) = modelmatrix * betas({i, i + N_PAR}).matrix();
        est = exp(est.array());
        for (int i = 0; i < N_AGE; i++) {
            qM.setZero(); 
            qM(0, 1) = est(i, 0); // debut
            qM(0, 2) = est(i, 1); // marriage from virgin
            qM(1, 2) = est(i, 2); // marriage from debut
            qM(2, {3,4,5}) = est(i, {3,4,5}); // marriage dissolution
            qM({3,4,5}, 6) = est(i, {6,7,8}); // disso > remarried
            qM(6, {3,4,5}) = est(i, {3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
            qM.diagonal() = T(-1) * qM.rowwise().sum();
            memcpy(&masterQ(0) + i * len, &qM(0), sizeof(T) * len);
            pM = expm(qM);
            memcpy(&masterP(0) + i * len, &pM(0), sizeof(T) * len);
            for (int j = 0; j < N_D; j++) { // prep delta
                matrix<T> tmp = qM * j;
                pD = expm(tmp);
                memcpy(&masterD(0) + i * len * N_D + j * len, &pD(0), sizeof(T) * len);
            }
        }
    };
    matrix<T> operator()(int a, bool diag = true){
        matrix<T> ans = masterP(Eigen::seqN(a * len, len)).reshaped(N_Q, N_Q);
        if (!diag) for (int i = 0; i < N_Q; ++i) ans(i, i) = T(0);
        return ans;
    };
    matrix<T> operator()(int a, int d, bool diag = true){
        matrix<T> ans = masterD(Eigen::seqN(a * len * N_D + d * len, len)).reshaped(N_Q, N_Q);
        return ans;
    }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(afs);
  DATA_IVECTOR(aam);

  DATA_IVECTOR(VV);
  DATA_IVECTOR(VX);
  DATA_IVECTOR(VM);
  DATA_IVECTOR(XX);
  DATA_IVECTOR(XM);
  DATA_IVECTOR(MM);
  DATA_IVECTOR(MJ);
  DATA_IVECTOR(J);
  DATA_IVECTOR(delta);
  DATA_VECTOR(count);
  
  DATA_MATRIX(modelmatrix);

  // priors
  DATA_VECTOR(prior_base);
  DATA_VECTOR(prior_t);

  // Coefs
  PARAMETER_VECTOR(betas);

// base rate | intercept
  vector<Type> intercepts = betas(seqN(0, N_PAR));
  prior -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();
  
  // Soft-constraints remarried to be the same
  prior -= dnorm(betas(6) - betas(7), Type(0), Type(0.001), true);
  prior -= dnorm(betas(6) - betas(8), Type(0), Type(0.001), true);

  // age's coeff | time in the hazard
  vector<Type> beta_t = betas(seqN(N_PAR, N_PAR));
  prior -= dnorm(beta_t, prior_t(0), prior_t(1), true).sum();

  // Soft-constraints remarried to be the same
  prior -= dnorm(betas(N_PAR + 6) - betas(N_PAR + 7), Type(0), Type(0.001), true);
  prior -= dnorm(betas(N_PAR + 6) - betas(N_PAR + 8), Type(0), Type(0.001), true);

  Kube<Type> KM(modelmatrix, betas);

  matrix<Type> o(afs.size(), 7);
  o.setOnes();

  for (int i = 0; i < afs.size(); i++) {
    if (VV[i] > 0)
        for (int j = 0; j <= VV[i]; j++)
            o(i, 0) *= KM(j)(0,0);
    if (VX[i] > 0)
            o(i, 1) *= KM(VX[i])(0,1);
    if (XX[i] > 0)
        for (int j = afs[i]; j <= afs[i] + XX[i]; j++) 
            o(i, 2) *= KM(j)(1,1);
    if (XM[i] > 0)
            o(i, 3) *= KM(XM[i])(1,2);
    if (MM[i] > 0)
        for (int j = aam[i]; j <= aam[i] + MM[i]; j++) 
            o(i, 4) *= KM(j)(2,2);
    if (MJ[i] > 0)
            o(i, 5) *= KM(MJ[i], delta[i])(2,J[i]);
    if (VM[i] > 0)
            o(i, 6) *= KM(VM[i])(0,2);
  }
  REPORT(o);
  o = log(o.array()); // all the 1s disapear here
  vector<Type> oo = o.rowwise().sum();
  oo *= count; // element-wise
  dll += prior - oo.sum();

  REPORT(betas);
  REPORT(KM.est);
  REPORT(KM.masterP);
  REPORT(KM.masterQ);
  REPORT(KM.masterD);
  return dll;
}
