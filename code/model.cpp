#include <TMB.hpp>
#include <unsupported/Eigen/MatrixFunctions> // matrix exp

template<class Type>
matrix<Type> pM(int x, vector<Type> betav, vector<Type> qv, bool diag = true) {
    Type lp = betav[0] + betav[1] * (x/100);
    qv = qv * exp(lp);
    matrix<Type> qM(7,7);
    qM.setZero();
    qM(0, 1) = qv[0];
    qM(1, 2) = qv[1];
    qM(2, 3) = qv[2];
    qM(2, 4) = qv[3];
    qM(2, 5) = qv[4];
    qM(3, 6) = qv[5];
    qM(4, 6) = qv[6];
    qM(5, 6) = qv[7];
    qM(6, 3) = qv[8];
    qM(6, 4) = qv[9];
    qM(6, 5) = qv[10]; // this might be shortened such as A(Eigen::placeholders::all,{4,2,5,5,3})?
    qM.diagonal() = Type(-1) * qM.rowwise().sum();
    matrix<Type> mexp = qM.exp(); // https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#matrixbase_exp
    if (!diag) for (int i = 0; i < 7; ++i) mexp(i, i) = 0.0;
    return(mexp);
}

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

  // Intensity: change to loggamma?
  PARAMETER_VECTOR(lqv); 
  prior -= dnorm(lqv, sd_q(0), sd_q(1), true).sum();
  vector<Type> qv = exp(lqv);

  vector<Type> o(afs.size());
  o.setOnes();

  for (int i = 0; i < afs.size(); i++) {
    if (afs[i] == 0) {
        for (int a = 0; a < age[i]; a++) 
            o[i] *= pM(a, betav, qv)(0, 0);
    } else {
        for (int a = 0; a < afs[i]; a++) 
            o[i] *= pM(a, betav, qv)(0, 0);
        o[i] *= pM(afs[i], betav, qv)(0, 1);
        if (n_X[i] == 0 & aam[i] == 0) 
            o[i] *= pM(afs[i], betav, qv)(1, 1);
        if (n_X[i] > 0) 
            for (int a = afs[i]; a < (afs[i] + n_X[i]); a++)
                o[i] *= pM(a, betav, qv)(1, 1); // XX(afs:(afs+n_X))
        if (aam[i] > 0) {
            o[i] *= pM(aam[i], betav, qv)(1, 2); // XM(aam)
            if (u[i] == 1) {
                if (n_M[i] == 0) 
                    o[i] *= pM(aam[i], betav, qv)(2, 2); // MM(aam)
                if (n_M[i] > 0)
                    for (int a = aam[i]; a < age[i]; a++) 
                        o[i] *= pM(a, betav, qv)(2, 2); // MM(aam:age)
                if (n_M[i] == -2908) {
                    int J = ms[i] - 1;
                    if (n_N[i] == 0) 
                        o[i] *= pM(aam[i], betav, qv)(2, J); // MJ(aam)
                    if (n_N[i] > 0) {
                        matrix<Type> oo(2, 2); oo.setOnes();
                        if (n_N[i] > 1)
                            for (int a = 1; a < n_N[i] - 1; a++)
                                oo *= pM(aam[i] + a, betav, qv)({2, J}, {2, J});
                        o[i] *= pM(aam[i], betav, qv)(2, {2, J}) * oo * 
                                pM(aam[i] + n_N[i], betav, qv)({2, J}, J);
                    }
                }
            }
            if (u[i] > 1) {
                int J = ms[i] - 1;
                if (n_N[i] == 3) { // M...J(aam:age+3)
                    o[i] *= (pM(aam[i] + 0, betav, qv)(2, {3,4,5}) *
                            pM(aam[i] + 1, betav, qv)({3,4,5}, 6) *
                            pM(aam[i] + 2, betav, qv)(6, J)).value(); // 77 = RM = RR =? MM
                }
                if (n_N[i] == 4) { // M...J(aam:age+4)
                    if (ms[i] != 7)
                        o[i] *= (pM(aam[i] + 0, betav, qv)(2, {2,3,4,5}) *
                                pM(aam[i] + 1, betav, qv)({2,3,4,5}, {3,4,5,6}) *
                                pM(aam[i] + 2, betav, qv, false)({3,4,5,6}, {J, 6}) *
                                pM(aam[i] + 3, betav, qv)({J, 6}, J)).value();
                    if (ms[i] == 7)
                        o[i] *= (pM(aam[i] + 0, betav, qv)(2, {2,3,4,5}) *
                                pM(aam[i] + 1, betav, qv)({2,3,4,5}, {3,4,5,6}) *
                                pM(aam[i] + 2, betav, qv)({3,4,5,6}, {3,4,5,6}) *
                                pM(aam[i] + 3, betav, qv)({3,4,5,6}, 6)).value();
                }
                if (n_N[i] > 4) { // M...J(aam:age+5+)"
                    matrix<Type> oo(5, 5); oo.setOnes();
                    matrix<Type> qp(2, 2); oo.setOnes();
                    if (n_N[i] > 5) {
                        for (int a = 2; a < (n_N[i] - 4); a++) {
                            oo *= pM(aam[i] + a, betav, qv)({2,3,4,5,6}, {2,3,4,5,6});
                            if (J != 7)
                                qp *= pM(aam[i] + a, betav, qv)({2,J}, {2,J});
                        }
                    }
                    if (ms[i] != 7)
                        o[i] *=
                        (
                            pM(aam[i] + 0, betav, qv)(2, {2,3,4,5}) *
                            pM(aam[i] + 1, betav, qv)({2,3,4,5}, {2,3,4,5,6}) *
                            oo *
                            pM(aam[i] + n_N[i] - 3, betav, qv)({2,3,4,5,6}, {3,4,5,6}) *
                            pM(aam[i] + n_N[i] - 2, betav, qv)({3,4,5,6}, {J, 6}) *
                            pM(aam[i] + n_N[i] - 1, betav, qv)({J, 6}, J) -
                            // substract one union cases
                            pM(aam[i] + 0, betav, qv)(2, {2, J}) *
                            pM(aam[i] + 1, betav, qv)({2, J}, {2, J}) *
                            qp *
                            pM(aam[i] + n_N[i] - 3, betav, qv)({2, J}, J) *
                            pM(aam[i] + n_N[i] - 2, betav, qv)(J, J) *
                            pM(aam[i] + n_N[i] - 1, betav, qv)(J, J)
                        ).value();
                    if (ms[i] == 7)
                        o[i] *=
                        (
                            pM(aam[i] + 0, betav, qv)(2, {2,3,4,5}) *
                            pM(aam[i] + 1, betav, qv)({2,3,4,5}, {2,3,4,5,6}) *
                            oo *
                            pM(aam[i] + 2, betav, qv)({2,3,4,5,6}, {3,4,5,6}) *
                            pM(aam[i] + 3, betav, qv)({3,4,5,6}, {3,4,5,6}) *
                            pM(aam[i] + 4, betav, qv)({3,4,5,6}, J)
                        ).value();
                }
            }
        }
    }
  }
  dll += prior + sum(log(o));
  REPORT(betav);
  REPORT(qv);
  return dll;
}
