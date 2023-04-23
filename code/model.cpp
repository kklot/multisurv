#include <TMB.hpp>

#define AGE_MAX 50
#define N_AGE 51
#define N_D 51 // differences, more than needed
#define N_PAR 8
#define N_Q 7

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
            qM(1, 2) = est(i, 1); // marriage
            qM(2, {3,4,5}) = est(i, {2,3,4}); // marriage dissolution
            qM({3,4,5}, 6) = est(i, {5,6,7}); // disso > remarried
            qM(6, {3,4,5}) = est(i, {2,3,4}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
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

template <class T> // 3 or 4 need to implement tophi
struct mAR {
    matrix<T> Sigma;
    vector<T> phi;
    mAR () {};
    mAR (int m) : Sigma(m, m), phi(m) {Sigma.setIdentity();};
    vector<T> to_phi (vector<T> x) {
        if (x.size() != 2) Rf_error("not implemented");
        vector<T> psi(2);
        psi[0] = 2. * exp(x[0]) / (1. + exp(x[0])) - 1.;
        psi[1] = 2. * exp(x[1]) / (1. + exp(x[1])) - 1.;
        phi[1] = psi[1];
        phi[0] = psi[0] * (1.0 - phi[1]);
        if (phi[1] == -1) phi[1] += FLT_EPSILON;
        if (phi[1] == 1 - phi.abs()(0) ) phi[1] -= DBL_EPSILON;
        return phi;
    };
    T operator()(vector<T> pacf, vector<T> x, bool zeroing = true){
        T dll = 0.;
        dll += density::MVNORM(Sigma)(pacf);
        phi = to_phi(pacf);
        dll += density::ARk(phi)(x);
        if (zeroing)
        dll -= dnorm(sum(x), T(0), T(0.001) * x.size(), true);
        return dll;
    };
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(afs);
  DATA_IVECTOR(age);
  DATA_IVECTOR(aam);

  DATA_IVECTOR(VV);
  DATA_IVECTOR(VX);
  DATA_IVECTOR(XX);
  DATA_IVECTOR(XM);
  DATA_IVECTOR(MM);
  DATA_IVECTOR(MJ);
  DATA_IVECTOR(J);
  DATA_IVECTOR(delta);
  
  DATA_MATRIX(modelmatrix);

  // priors
  DATA_VECTOR(sd_b);
  DATA_VECTOR(sd_q);

  // Coefs
  PARAMETER_VECTOR(betas); 
  prior -= dnorm(betas, sd_q(0), sd_q(1), true).sum();
  // Soft-constraints remarried to be the same
  prior -= dnorm(betas(5) - betas(6), Type(0), Type(0.001), true);
  prior -= dnorm(betas(5) - betas(7), Type(0), Type(0.001), true);
  prior -= dnorm(betas(N_PAR + 5) - betas(N_PAR + 6), Type(0), Type(0.001), true);
  prior -= dnorm(betas(N_PAR + 5) - betas(N_PAR + 7), Type(0), Type(0.001), true);
  Kube<Type> KM(modelmatrix, betas);

  matrix<Type> o(afs.size(), 6);
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
  }
  REPORT(o);
  o = log(o.array());
  dll += prior - o.sum();

  REPORT(betas);
  REPORT(KM.est);
  REPORT(KM.masterP);
  REPORT(KM.masterQ);
  REPORT(KM.masterD);
  return dll;
}
