#include <TMB.hpp>

#define AGE_MAX 50
#define N_AGE 51
#define N_PAR 8
#define N_Q 7

template <class T> 
struct Kube {
    vector<T> masterQ;
    vector<T> masterP;
    matrix<T> qM;
    matrix<T> pM;
    matrix<T> est;
    int len;
    Kube () {};
    Kube (vector<T> betav, vector<T> qv) : 
        masterQ(N_AGE * N_Q * N_Q), 
        masterP(N_AGE * N_Q * N_Q), 
        qM(N_Q, N_Q), 
        pM(N_Q, N_Q), 
        est(N_AGE, N_PAR),
        len(N_Q * N_Q)
    {
        masterP.setZero();
        masterQ.setZero();
        betav = exp(betav);
        matrix<T> betam = betav.reshaped(N_AGE, N_PAR);
        betam.array().rowwise() *= qv.transpose(); // *= only array
        est = betam; // save for reporting
        for (int i = 0; i < N_AGE; i++) {
            qM.setZero(); 
            qM(0, 1) = betam(i, 0); // debut
            qM(1, 2) = betam(i, 1); // marriage
            qM(2, {3,4,5}) = betam(i, {2,3,4}); // marriage dissolution
            qM({3,4,5}, 6) = betam(i, {5,6,7}); // disso > remarried
            qM(6, {3,4,5}) = betam(i, {2,3,4}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
            qM.diagonal() = T(-1) * qM.rowwise().sum();
            memcpy(&masterQ(0) + i * len, &qM(0), sizeof(T) * len);
            pM = expm(qM);
            memcpy(&masterP(0) + i * len, &pM(0), sizeof(T) * len);
        }
    };
    matrix<T> operator()(int a, bool diag = true){
        matrix<T> ans = masterP(Eigen::seqN(a * len, len)).reshaped(N_Q, N_Q);
        if (!diag) for (int i = 0; i < N_Q; ++i) ans(i, i) = T(0);
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
  DATA_IVECTOR(age);
  DATA_IVECTOR(aam);

  DATA_IVECTOR(u);
  DATA_IVECTOR(ms);
  DATA_IVECTOR(n_X);
  DATA_IVECTOR(n_M);
  DATA_IVECTOR(n_N);

  // priors
  DATA_VECTOR(sd_b);
  DATA_VECTOR(sd_q);

  // Data model
  PARAMETER_VECTOR(betav);
  prior -= dnorm(betav, sd_b(0), sd_b(1), true).sum();
  // constraints IID within a vector
  vector<Type> csum = betav.reshaped(AGE_MAX + 1, N_PAR).colwise().sum();
  prior -= dnorm(csum, Type(0.0), Type(0.001) * (AGE_MAX + 1), true).sum();

  // Baseline intensity
  PARAMETER_VECTOR(lqv); 
  prior -= dnorm(lqv, sd_q(0), sd_q(1), true).sum() + sum(lqv);
  vector<Type> qv = exp(lqv);

  Kube<Type> KM(betav, qv);

  vector<Type> o(afs.size());
  o.setOnes();

  for (int i = 0; i < afs.size(); i++) {
    if (afs[i] == 0) {
        for (int a = 0; a < age[i]; a++) 
            o[i] *= KM(a)(0, 0);
    } else {
        for (int a = 0; a < afs[i]; a++) 
            o[i] *= KM(a)(0, 0);
        o[i] *= KM(afs[i])(0, 1);
        if (n_X[i] == 0 & aam[i] == 0) 
            o[i] *= KM(afs[i])(1, 1);
        if (n_X[i] > 0) 
            for (int a = afs[i]; a < (afs[i] + n_X[i]); a++)
                o[i] *= KM(a)(1, 1); // XX(afs:(afs+n_X))
        if (aam[i] > 0) {
            o[i] *= KM(aam[i])(1, 2); // XM(aam)
            if (u[i] == 1) {
                if (n_M[i] == 0) 
                    o[i] *= KM(aam[i])(2, 2); // MM(aam)
                if (n_M[i] > 0)
                    for (int a = aam[i]; a < age[i]; a++) 
                        o[i] *= KM(a)(2, 2); // MM(aam:age)
                if (n_M[i] == -2908) {
                    int J = ms[i] - 1;
                    if (n_N[i] == 0) 
                        o[i] *= KM(aam[i])(2, J); // MJ(aam)
                    if (n_N[i] > 0) {
                        matrix<Type> oo(2, 2); oo.setOnes();
                        if (n_N[i] > 1)
                            for (int a = 1; a < n_N[i] - 1; a++)
                                oo *= KM(aam[i] + a)({2, J}, {2, J});
                        o[i] *= KM(aam[i])(2, {2, J}) * oo * 
                                KM(aam[i] + n_N[i])({2, J}, J);
                    }
                }
            }
            if (u[i] > 1) {
                int J = ms[i] - 1;
                if (n_N[i] == 3) { // M...J(aam:age+3)
                    o[i] *= (KM(aam[i] + 0)(2, {3,4,5}) *
                            KM(aam[i] + 1)({3,4,5}, 6) *
                            KM(aam[i] + 2)(6, J)).value(); // 77 = RM = RR =? MM
                }
                if (n_N[i] == 4) { // M...J(aam:age+4)
                    if (ms[i] != 7)
                        o[i] *= (KM(aam[i] + 0)(2, {2,3,4,5}) *
                                KM(aam[i] + 1)({2,3,4,5}, {3,4,5,6}) *
                                KM(aam[i] + 2, false)({3,4,5,6}, {J, 6}) *
                                KM(aam[i] + 3)({J, 6}, J)).value();
                    if (ms[i] == 7)
                        o[i] *= (KM(aam[i] + 0)(2, {2,3,4,5}) *
                                KM(aam[i] + 1)({2,3,4,5}, {3,4,5,6}) *
                                KM(aam[i] + 2)({3,4,5,6}, {3,4,5,6}) *
                                KM(aam[i] + 3)({3,4,5,6}, 6)).value();
                }
                if (n_N[i] > 4) { // M...J(aam:age+5+)"
                    matrix<Type> oo(5, 5); oo.setOnes();
                    matrix<Type> qp(2, 2); oo.setOnes();
                    if (n_N[i] > 5) {
                        for (int a = 2; a < (n_N[i] - 4); a++) {
                            oo *= KM(aam[i] + a)({2,3,4,5,6}, {2,3,4,5,6});
                            if (J != 7)
                                qp *= KM(aam[i] + a)({2,J}, {2,J});
                        }
                    }
                    if (ms[i] != 7)
                        o[i] *=
                        (
                            KM(aam[i] + 0)(2, {2,3,4,5}) *
                            KM(aam[i] + 1)({2,3,4,5}, {2,3,4,5,6}) *
                            oo *
                            KM(aam[i] + n_N[i] - 3)({2,3,4,5,6}, {3,4,5,6}) *
                            KM(aam[i] + n_N[i] - 2)({3,4,5,6}, {J, 6}) *
                            KM(aam[i] + n_N[i] - 1)({J, 6}, J) -
                            // substract one union cases
                            KM(aam[i] + 0)(2, {2, J}) *
                            KM(aam[i] + 1)({2, J}, {2, J}) *
                            qp *
                            KM(aam[i] + n_N[i] - 3)({2, J}, J) *
                            KM(aam[i] + n_N[i] - 2)(J, J) *
                            KM(aam[i] + n_N[i] - 1)(J, J)
                        ).value();
                    if (ms[i] == 7)
                        o[i] *=
                        (
                            KM(aam[i] + 0)(2, {2,3,4,5}) *
                            KM(aam[i] + 1)({2,3,4,5}, {2,3,4,5,6}) *
                            oo *
                            KM(aam[i] + 2)({2,3,4,5,6}, {3,4,5,6}) *
                            KM(aam[i] + 3)({3,4,5,6}, {3,4,5,6}) *
                            KM(aam[i] + 4)({3,4,5,6}, J)
                        ).value();
                }
            }
        }
    }
  }
  dll += prior - sum(log(o));
  REPORT(betav);
  REPORT(qv);
  REPORT(KM.est);
  REPORT(KM.masterP);
  REPORT(KM.masterQ);
  return dll;
}
