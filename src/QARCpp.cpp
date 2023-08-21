// Functions for QAR

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

/////////////////
// KUMARASWAMY //
/////////////////

// [[Rcpp::export]]
double dKumaraswamy(double x, double a, double b) {
  return (a * b * pow(x, a - 1) * pow(1 - pow(x, a), b - 1));
}

arma::mat dKumaraswamyMat(arma::mat x, double a, double b) {
  return (a * b * pow(x, a - 1) % pow(1 - pow(x, a), b - 1));
}

// [[Rcpp::export]]
arma::cube dKumaraswamyCube(arma::cube x, double a, double b) {
  return (a * b * pow(x, a - 1) % pow(1 - pow(x, a), b - 1));
}

// [[Rcpp::export]]
double pKumaraswamy(double q, double a, double b) {
  return (1 - pow(1 - pow(q, a), b));
}

// [[Rcpp::export]]
arma::cube qKumaraswamy(arma::cube p, double a, double b) {
  return pow(1 - pow(1 - p, 1 / b), 1 / a);
}

/////////////////
//    OTHER    //
/////////////////

// dnorm ~ log mu = 0
//
// Proportional in x to normal log-density with mu = 0
//
// @param x quantile
// @param sd standard deviation
// @return log-density
// [[Rcpp::export]]
double logdnorm(double x, double sd) {
  return (- pow(x / sd, 2) / 2);
}

// [[Rcpp::export]]
double logdmvnorm(arma::vec x, double mu, double sigma2, arma::mat Rinv, int n) {
  arma::vec vector = x - mu;
  return - (n * log(sigma2) + (vector.t() * Rinv * vector).eval()(0,0) / sigma2) / 2;
}

// Online Update of Mean and Variance
//
// Update mu and Sigma for tunning the Metropolis step
//
// @param x Vector New value
// @param mu Vector Old mean
// @param Sigma Matrix Old variance
// @param n Number of data with x
// @return List with updated mu and Sigma
// [[Rcpp::export]]
Rcpp::List MuSigmaUpdate(arma::vec x, arma::vec mu, arma::mat Sigma, int n) {
  arma::vec u = (x - mu) / n;
  mu += u;
  Sigma = (n - 2) * Sigma / (n - 1) + u * u.t() * n;
  return Rcpp::List::create(mu, Sigma);
}

// [[Rcpp::export]]
double logitInv(double x, double a, double b) {
  return (a + b * exp(x)) / (1 + exp(x));
  //return 1 / (1 + exp(-x));
}

// [[Rcpp::export]]
arma::vec logitInvVec(arma::vec x, double a, double b) {
  return (a + b * exp(x)) / (1 + exp(x));
  //return 1 / (1 + exp(-x));
}

arma::vec transform(arma::vec x, int K) {
  x = exp(x);
  double sum1 = 1 + arma::accu(x(arma::span(0, K - 2)));
  double sum2 = 1 + arma::accu(x(arma::span(K - 1, 2 * K - 3)));
  x(arma::span(0, K - 2)) = x(arma::span(0, K - 2)) / sum1;
  x(arma::span(K - 1, 2 * K - 3)) = x(arma::span(K - 1, 2 * K - 3)) / sum2;
  return x;
}

//////////////////
//    QAR1K1    //
//////////////////

// [[Rcpp::export]]
double fQAR1K1(double x, double y1, double y2, 
               double a1, double b1, double a2, double b2) {
  return(y2 - (y1 * pKumaraswamy(x, a1, b1) + (1 - y1) * pKumaraswamy(x, a2, b2)));
}

// [[Rcpp::export]]
double brentQAR1K1(double y1, double y2, 
                   double a1, double b1, double a2, double b2, 
                   double tol) {
  int n = 1;
  double a = 0;
  double b = 1;
  double d;
  double aux;
  double fa = y2;
  double fb = y2 - 1;
  if (y2 < 0.5) {
    aux = a;
    a = b;
    b = aux;
    aux = fa;
    fa = fb;
    fb = aux;
  }
  double r = a;
  double s;
  double fr = fa;
  double fs;
  bool mflag = true;
  while (n < 10001) {
    if (fa != fr && fb != fr) {
      s = (a*fb*fr)/((fa-fb)*(fa-fr))+(b*fa*fr)/((fb-fr)*(fb-fa))+(r*fa*fb)/((fr-fa)*(fr-fb));
    } else {
      s = b-(fb*(b-a))/(fb-fa);
    }
    if ((((3 * a + b) / 4) > s || s > b) ||
        (mflag && abs(s - b) >= abs(b - r) / 2) ||
        (!mflag && abs(s - b) >= abs(r - d) / 2)
       ) {
      s = (b + a) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = fQAR1K1(s, y1, y2, a1, b1, a2, b2);
    d = r;
    r = b;
    fr = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fa) < abs(fb)) {
      aux = a;
      a = b;
      b = aux;
      aux = fa;
      fa = fb;
      fb = aux;
    }
    if (abs(b - a) < tol || fs == 0) {
      return s;
    } 
    ++n;
  }
  Rcpp::warning("Warning: The rootfinder did not converge \n");
  return s;
}

//double bisecQAR1K1(double y1, double y2, 
//                   double a1, double b1, double a2, double b2, 
//                   double tol) {
//  int n = 0;
//  double a = 0;
//  double b = 1;
//  double r;
//  while ((b - a > tol) && (n < 1000)) {
//    ++n;
//    r = (a + b) / 2;
//    if (fQAR1K1(a, y1, y2, a1, b1, a2, b2) * fQAR1K1(r, y1, y2, a1, b1, a2, b2) > 0) {
//      a = r;
//    } else {
//      b = r;
//    }
//  }
//  r = (a + b) / 2;
//  return r;
//}

// [[Rcpp::export]]
double logLikelihoodQAR1K1(int T, int L, arma::mat Y, 
                           double a1, double b1, double a2, double b2,
                           double tol) {
  
  double tauY;
  double logLikelihood = 0;
  
  for (int t = 0; t < T; ++t) {
    for (int l = 1; l < L; ++l) {
      tauY = brentQAR1K1(Y(t,l-1), Y(t,l), a1, b1, a2, b2, tol);
      logLikelihood -= log(Y(t,l-1) * dKumaraswamy(tauY, a1, b1) + (1 - Y(t,l-1)) * dKumaraswamy(tauY, a2, b2));
    }
  }

  return logLikelihood;
}

// [[Rcpp::export]]
double logfallQAR1K1(int T, int L, arma::vec x, arma::mat Y, double prior, double tol) {
  return logLikelihoodQAR1K1(T, L, Y, exp(x(0)), exp(x(1)), exp(x(2)), exp(x(3)), tol) +
    logdnorm(x(0), prior) +
    logdnorm(x(1), prior) +
    logdnorm(x(2), prior) +
    logdnorm(x(3), prior);
}

// [[Rcpp::export]]
arma::vec rwBmetropolisQAR1K1(arma::vec x, arma::mat sd, int T, int L, arma::mat Y, double prior, double tol) {
  arma::vec y(5); 
  y(arma::span(0, 3)) = x(arma::span(0, 3)) + arma::chol(sd, "lower") * arma::randn(4);
  y(4) = logfallQAR1K1(T, L, y(arma::span(0, 3)), Y, prior, tol);
  while (y(4) == INFINITY) {
    y(arma::span(0, 3)) = x(arma::span(0, 3)) + arma::chol(sd, "lower") * arma::randn(4);
    y(4) = logfallQAR1K1(T, L, y(arma::span(0, 3)), Y, prior, tol);
  }
  double A = y(4) - x(4);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

/////////////////
//   QAR1K2    //
/////////////////

// [[Rcpp::export]]
double fQAR1K2(double x, double y1, double y2, 
               double a1, double b1, double a2, double b2,
               double a3, double b3, double a4, double b4,
               double pi1, double pi2) {
  double mix1 = pi1 * pKumaraswamy(x, a1, b1) + (1 - pi1) * pKumaraswamy(x, a2, b2);
  double mix2 = pi2 * pKumaraswamy(x, a3, b3) + (1 - pi2) * pKumaraswamy(x, a4, b4);
  return(y2 - (y1 * mix1 + (1 - y1) * mix2));
}

// [[Rcpp::export]]
double brentQAR1K2(double y1, double y2, 
                   double a1, double b1, double a2, double b2,
                   double a3, double b3, double a4, double b4,
                   double pi1, double pi2, double tol) {
  int n = 1;
  double a = 0;
  double b = 1;
  double d;
  double aux;
  double fa = y2;
  double fb = y2 - 1;
  if (y2 < 0.5) {
    aux = a;
    a = b;
    b = aux;
    aux = fa;
    fa = fb;
    fb = aux;
  }
  double r = a;
  double s;
  double fr = fa;
  double fs;
  bool mflag = true;
  while (n < 10001) {
    if (fa != fr && fb != fr) {
      s = (a*fb*fr)/((fa-fb)*(fa-fr))+(b*fa*fr)/((fb-fr)*(fb-fa))+(r*fa*fb)/((fr-fa)*(fr-fb));
    } else {
      s = b-(fb*(b-a))/(fb-fa);
    }
    if ((((3 * a + b) / 4) > s || s > b) ||
        (mflag && abs(s - b) >= abs(b - r) / 2) ||
        (!mflag && abs(s - b) >= abs(r - d) / 2)// ||
        //(mflag && abs(b - r) < delta) ||
        //(!mflag && abs(r - d) < delta)
       ) {
      s = (b + a) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = fQAR1K2(s, y1, y2, a1, b1, a2, b2, a3, b3, a4, b4, pi1, pi2);
    d = r;
    r = b;
    fr = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fa) < abs(fb)) {
      aux = a;
      a = b;
      b = aux;
      aux = fa;
      fa = fb;
      fb = aux;
    }
    if (abs(b - a) < tol || fs == 0) {
      return s;
    } 
    ++n;
  }
  Rcpp::warning("Warning: The rootfinder did not converge \n");
  return s;
}

//double bisecQAR1K2(double y1, double y2, 
//                   double a1, double b1, double a2, double b2,
//                   double a3, double b3, double a4, double b4,
//                   double pi1, double pi2, double tol) {
//  int n = 1;
//  double a = 0;
//  double b = 1;
//  double r;
//  double fr;
//  double fa = y2;
//  while (n < 10001) {
//    r = (a + b) / 2;
//    fr = fQAR1K2(r, y1, y2, a1, b1, a2, b2, a3, b3, a4, b4, pi1, pi2);
//    if (b - a < tol || fr == 0) {
//      return r;
//    } 
//    ++n;
//    if (fa * fr > 0) {
//      a = r;
//      fa = fr;
//    } else {
//      b = r;
//    }
//  }
//  Rcpp::warning("Warning: The rootfinder did not converge \n");
//  return r;
//}

// [[Rcpp::export]]
double logLikelihoodQAR1K2(int T, int L, arma::mat Y, 
                           double a1, double b1, double a2, double b2,
                           double a3, double b3, double a4, double b4,
                           double pi1, double pi2, double tol) {
  
  double tauY;
  double mix1;
  double mix2;
  double logLikelihood = 0;
  
  for (int t = 0; t < T; ++t) {
    for (int l = 1; l < L; ++l) {
      tauY = brentQAR1K2(Y(t,l-1), Y(t,l), a1, b1, a2, b2, a3, b3, a4, b4, pi1, pi2, tol);
      mix1 = pi1 * dKumaraswamy(tauY, a1, b1) + (1 - pi1) * dKumaraswamy(tauY, a2, b2);
      mix2 = pi2 * dKumaraswamy(tauY, a3, b3) + (1 - pi2) * dKumaraswamy(tauY, a4, b4);
      logLikelihood -= log(Y(t,l-1) * mix1 + (1 - Y(t,l-1)) * mix2);
    }
  }
  
  return logLikelihood;
}

// [[Rcpp::export]]
double logfallQAR1K2(int T, int L, arma::vec x, arma::mat Y, double prior, double tol) {
  return logLikelihoodQAR1K2(T, L, Y, exp(x(0)), exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), exp(x(6)), exp(x(7)), logitInv(x(8), 0, 0.5), logitInv(x(9), 0, 0.5), tol) +
    logdnorm(x(0), prior) +
    logdnorm(x(1), prior) +
    logdnorm(x(2), prior) +
    logdnorm(x(3), prior) +
    logdnorm(x(4), prior) +
    logdnorm(x(5), prior) +
    logdnorm(x(6), prior) +
    logdnorm(x(7), prior) +
    R::dlogis(x(8), 0, 1, 1) +
    R::dlogis(x(9), 0, 1, 1);
} 

// [[Rcpp::export]]
arma::vec rwBmetropolisQAR1K2(arma::vec x, arma::mat sd, int T, int L, arma::mat Y, double prior, double tol) {
  arma::vec y(11);
  y(arma::span(0, 9)) = x(arma::span(0, 9)) + arma::chol(sd, "lower") * arma::randn(10);
  y(10) = logfallQAR1K2(T, L, y(arma::span(0, 9)), Y, prior, tol);
  while (y(10) == INFINITY) {
    y(arma::span(0, 9)) = x(arma::span(0, 9)) + arma::chol(sd, "lower") * arma::randn(10);
    y(10) = logfallQAR1K2(T, L, y(arma::span(0, 9)), Y, prior, tol);
  }
  double A = y(10) - x(10);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

/////////////////
//    QAR1K    //
/////////////////

// [[Rcpp::export]]
double fQAR1K(double x, double y1, double y2, 
              int K, arma::vec ab, arma::vec lambda, double lambda1K, double lambda2K) {
  double eta1 = 0;
  double eta2 = 0;
  for (int k = 0; k < K - 1; ++k) {
    eta1 += lambda(k) * pKumaraswamy(x, ab(2 * k), ab(2 * k + 1));
    eta2 += lambda(k + K - 1) * pKumaraswamy(x, ab(2 * k), ab(2 * k + 1));
  }
  eta1 += lambda1K * pKumaraswamy(x, ab(2 * K - 2), ab(2 * K - 1));
  eta2 += lambda2K * pKumaraswamy(x, ab(2 * K - 2), ab(2 * K - 1));
  
  return(y2 - (y1 * eta1 + (1 - y1) * eta2));
}

// [[Rcpp::export]]
double brentQAR1K(double y1, double y2, 
                  int K, arma::vec ab, arma::vec lambda, double lambda1K, double lambda2K,
                  double tol) {
  int n = 1;
  double a = 0;
  double b = 1;
  double d;
  double aux;
  double fa = y2;
  double fb = y2 - 1;
  if (y2 < 0.5) {
    aux = a;
    a = b;
    b = aux;
    aux = fa;
    fa = fb;
    fb = aux;
  }
  double r = a;
  double s;
  double fr = fa;
  double fs;
  bool mflag = true;
  while (n < 10001) {
    if (fa != fr && fb != fr) {
      s = (a*fb*fr)/((fa-fb)*(fa-fr))+(b*fa*fr)/((fb-fr)*(fb-fa))+(r*fa*fb)/((fr-fa)*(fr-fb));
    } else {
      s = b-(fb*(b-a))/(fb-fa);
    }
    if ((((3 * a + b) / 4) > s || s > b) ||
        (mflag && abs(s - b) >= abs(b - r) / 2) ||
        (!mflag && abs(s - b) >= abs(r - d) / 2)
    ) {
      s = (b + a) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = fQAR1K(s, y1, y2, K, ab, lambda, lambda1K, lambda2K);
    d = r;
    r = b;
    fr = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fa) < abs(fb)) {
      aux = a;
      a = b;
      b = aux;
      aux = fa;
      fa = fb;
      fb = aux;
    }
    if (abs(b - a) < tol || fs == 0) {
      return s;
    } 
    ++n;
  }
  Rcpp::warning("Warning: The rootfinder did not converge \n");
  return s;
}

//double bisecQAR1K(double y1, double y2, 
//                  int K, arma::vec ab, arma::vec lambda, double lambda1K, double lambda2K,
//                  double tol) {
//  int n = 0;
//  double a = 0;
//  double b = 1;
//  double r;
//  while ((b - a > tol) && (n < 1000)) {
//    ++n;
//    r = (a + b) / 2;
//    if (fQAR1K(a, y1, y2, K, ab, lambda, lambda1K, lambda2K) * fQAR1K(r, y1, y2, K, ab, lambda, lambda1K, lambda2K) > 0) {
//      a = r;
//    } else {
//      b = r;
//    }
//  }
//  r = (a + b) / 2;
//  return r;
//}

// [[Rcpp::export]]
double logLikelihoodQAR1K(int T, int L, arma::mat Y, 
                          int K, arma::vec ab, arma::vec lambda,
                          double tol) {
  
  double lambda1K = 1 - arma::accu(lambda(arma::span(0, K - 2)));
  double lambda2K = 1 - arma::accu(lambda(arma::span(K - 1, 2 * K - 3)));
  
  double tauY;
  double eta1;
  double eta2;
  double logLikelihood = 0;
  
  for (int t = 0; t < T; ++t) {
    for (int l = 1; l < L; ++l) {
      tauY = brentQAR1K(Y(t,l-1), Y(t,l), K, ab, lambda, lambda1K, lambda2K, tol);
      eta1 = 0;
      eta2 = 0;
      for (int k = 0; k < K - 1; ++k) {
        eta1 += lambda(k) * dKumaraswamy(tauY, ab(2 * k), ab(2 * k + 1));
        eta2 += lambda(k + K - 1) * dKumaraswamy(tauY, ab(2 * k), ab(2 * k + 1));
      }
      eta1 += lambda1K * dKumaraswamy(tauY, ab(2 * K - 2), ab(2 * K - 1));
      eta2 += lambda2K * dKumaraswamy(tauY, ab(2 * K - 2), ab(2 * K - 1));
      logLikelihood -= log(Y(t,l-1) * eta1 + (1 - Y(t,l-1)) * eta2);
    }
  }

  return logLikelihood;
}

// [[Rcpp::export]]
double logfallQAR1K(int T, int L, arma::vec x, arma::mat Y, int K, arma::vec ab, double prior, double tol) {
  double suma = logLikelihoodQAR1K(T, L, Y, K, ab, transform(x, K), tol);
  for (int k = 0; k < K - 1; ++k) {
    suma += logdnorm(x(k), prior) + logdnorm(x(k + K - 1), prior);
  }
  return suma;
} 

// [[Rcpp::export]]
arma::vec rwBmetropolisQAR1K(arma::vec x, arma::mat sd, int T, int L, arma::mat Y, int K, arma::vec ab, double prior, double tol) {
  arma::vec y(2 * K - 1);
  y(arma::span(0, 2 * K - 3)) = x(arma::span(0, 2 * K - 3)) + arma::chol(sd, "lower") * arma::randn(2 * K - 2);
  y(2 * K - 2) = logfallQAR1K(T, L, y(arma::span(0, 2 * K - 3)), Y, K, ab, prior, tol);
  while (y(2 * K - 2) == INFINITY) {
    y(arma::span(0, 2 * K - 3)) = x(arma::span(0, 2 * K - 3)) + arma::chol(sd, "lower") * arma::randn(2 * K - 2);
    y(2 * K - 2) = logfallQAR1K(T, L, y(arma::span(0, 2 * K - 3)), Y, K, ab, prior, tol);
  }
  double A = y(2 * K - 2) - x(2 * K - 2);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

//////////////////
//    QAR2K1    //
//////////////////

// [[Rcpp::export]]
double fQAR2K1(double x, double y1, double y2, double y3,
               double a1, double b1, 
               double a2, double b2,
               double a3, double b3,
               double pi) {
  double eta1 = pKumaraswamy(x, a1, b1);
  double eta2 = pKumaraswamy(x, a2, b2);
  double eta3 = pKumaraswamy(x, a3, b3);
  return(y3 - (pi * y2 * eta1 + (1 - pi) * y1 * eta2 + (1 - pi * y2 - (1 - pi) * y1) * eta3));
}

// [[Rcpp::export]]
double brentQAR2K1(double y1, double y2, double y3,
                   double a1, double b1, 
                   double a2, double b2,
                   double a3, double b3,
                   double pi, double tol) {
  int n = 1;
  double a = 0;
  double b = 1;
  double d;
  double aux;
  double fa = y3;
  double fb = y3 - 1;
  if (y3 < 0.5) {
    aux = a;
    a = b;
    b = aux;
    aux = fa;
    fa = fb;
    fb = aux;
  }
  double r = a;
  double s;
  double fr = fa;
  double fs;
  bool mflag = true;
  while (n < 10001) {
    if (fa != fr && fb != fr) {
      s = (a*fb*fr)/((fa-fb)*(fa-fr))+(b*fa*fr)/((fb-fr)*(fb-fa))+(r*fa*fb)/((fr-fa)*(fr-fb));
    } else {
      s = b-(fb*(b-a))/(fb-fa);
    }
    if ((((3 * a + b) / 4) > s || s > b) ||
        (mflag && abs(s - b) >= abs(b - r) / 2) ||
        (!mflag && abs(s - b) >= abs(r - d) / 2)
    ) {
      s = (b + a) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = fQAR2K1(s, y1, y2, y3, a1, b1, a2, b2, a3, b3, pi);
    d = r;
    r = b;
    fr = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fa) < abs(fb)) {
      aux = a;
      a = b;
      b = aux;
      aux = fa;
      fa = fb;
      fb = aux;
    }
    if (abs(b - a) < tol || fs == 0) {
      return s;
    } 
    ++n;
  }
  Rcpp::warning("Warning: The rootfinder did not converge \n");
  return s;
}

//double bisecQAR2K1(double y1, double y2, double y3,
//                   double a1, double b1, 
//                   double a2, double b2,
//                   double a3, double b3,
//                   double pi, double tol) {
//  int n = 0;
//  double a = 0;
//  double b = 1;
//  double r;
//  while ((b - a > tol) && (n < 1000)) {
//    ++n;
//    r = (a + b) / 2;
//    if (fQAR2K1(a, y1, y2, y3, a1, b1, a2, b2, a3, b3, pi) * fQAR2K1(r, y1, y2, y3, a1, b1, a2, b2, a3, b3, pi) > 0) {
//      a = r;
//    } else {
//      b = r;
//    }
//  }
//  r = (a + b) / 2;
//  return r;
//}

// [[Rcpp::export]]
double logLikelihoodQAR2K1(int T, int L, arma::mat Y, 
                           double a1, double b1, 
                           double a2, double b2,
                           double a3, double b3,
                           double pi, double tol) {
  
  double tauY;
  double eta1;
  double eta2;
  double eta3;
  double logLikelihood = 0;
  
  for (int t = 0; t < T; ++t) {
    for (int l = 2; l < L; ++l) {
      tauY = brentQAR2K1(Y(t,l-2), Y(t,l-1), Y(t,l), a1, b1, a2, b2, a3, b3, pi, tol);
      eta1 = dKumaraswamy(tauY, a1, b1);
      eta2 = dKumaraswamy(tauY, a2, b2);
      eta3 = dKumaraswamy(tauY, a3, b3);
      logLikelihood -= log(pi * Y(t,l-1) * eta1 + (1 - pi) * Y(t,l-2) * eta2 + (1 - pi * Y(t,l-1) - (1 - pi) * Y(t,l-2)) * eta3);
    }
  }
  
  return logLikelihood;
}

// [[Rcpp::export]]
double logfallQAR2K1(int T, int L, arma::vec x, arma::mat Y, double prior, double tol) {
  return logLikelihoodQAR2K1(T, L, Y, exp(x(0)), exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), logitInv(x(6), 0, 1), tol) +
    logdnorm(x(0), prior) +
    logdnorm(x(1), prior) +
    logdnorm(x(2), prior) +
    logdnorm(x(3), prior) +
    logdnorm(x(4), prior) +
    logdnorm(x(5), prior) +
    R::dlogis(x(6), 0, 1, 1);
} 

// [[Rcpp::export]]
arma::vec rwBmetropolisQAR2K1(arma::vec x, arma::mat sd, int T, int L, arma::mat Y, double prior, double tol) {
  arma::vec y(8); 
  y(arma::span(0, 6)) = x(arma::span(0, 6)) + arma::chol(sd, "lower") * arma::randn(7);
  y(7) = logfallQAR2K1(T, L, y(arma::span(0, 6)), Y, prior, tol);
  while (y(7) == INFINITY) {
    y(arma::span(0, 6)) = x(arma::span(0, 6)) + arma::chol(sd, "lower") * arma::randn(7);
    y(7) = logfallQAR2K1(T, L, y(arma::span(0, 6)), Y, prior, tol);
  }
  double A = y(7) - x(7);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

//////////////////
//Multiva QAR1k1//
//////////////////

// [[Rcpp::export]]
double logLikelihoodMQAR1K1(int T, int L, arma::cube Y, 
                            double a1, double b1, double a2, double b2,
                            double a3, double b3, double a4, double b4,
                            double rho, double tol) {
  
  arma::cube tauY(T, L, 2);
  for (int t = 0; t < T; ++t) {
    for (int l = 1; l < L; ++l) {
      tauY(t,l,0) = brentQAR1K1(Y(t,l-1,0), Y(t,l,0), a1, b1, a2, b2, tol);
      tauY(t,l,1) = brentQAR1K1(Y(t,l-1,1), Y(t,l,1), a3, b3, a4, b4, tol);
    }
  }
  
  Y.shed_col(L - 1);
  tauY.shed_col(0);
  
  arma::mat logLikelihood = -
    log(Y.slice(0) % dKumaraswamyMat(tauY.slice(0), a1, b1) + (1 - Y.slice(0)) % dKumaraswamyMat(tauY.slice(0), a2, b2)) -
    log(Y.slice(1) % dKumaraswamyMat(tauY.slice(1), a3, b3) + (1 - Y.slice(1)) % dKumaraswamyMat(tauY.slice(1), a4, b4));
  
  arma::mat R(2, 2, arma::fill::eye);
  R(0,1) = rho;
  R(1,0) = rho;
  arma::mat I(2, 2, arma::fill::eye);
  arma::mat Rstar = inv_sympd(R) - I;
  
  arma::mat product(T, L - 1);
  arma::vec vector(2);
  
  for (int t = 0; t < T; ++t) {
    for (int l = 0; l < L - 1; ++l) {
      vector(0) = R::qnorm(tauY(t,l,0), 0, 1, 1, 0);
      vector(1) = R::qnorm(tauY(t,l,1), 0, 1, 1, 0);
      product(t,l) = (vector.t() * Rstar * vector).eval()(0,0);
    }
  }
  
  arma::mat logCopula = - 0.5 * (log_det_sympd(R) + product);
  
  return arma::accu(logLikelihood) + arma::accu(logCopula);
}

// [[Rcpp::export]]
double logfallMQAR1K1(int T, int L, arma::vec x, arma::cube Y, double prior, double tol) {
  return logLikelihoodMQAR1K1(T, L, Y, exp(x(0)), exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), exp(x(6)), exp(x(7)), logitInv(x(8), -1, 1), tol) +
    logdnorm(x(0), prior) +
    logdnorm(x(1), prior) +
    logdnorm(x(2), prior) +
    logdnorm(x(3), prior) +
    logdnorm(x(4), prior) +
    logdnorm(x(5), prior) +
    logdnorm(x(6), prior) +
    logdnorm(x(7), prior) +
    R::dlogis(x(8), 0, 1, 1);
}

// [[Rcpp::export]]
arma::vec rwBmetropolisMQAR1K1(arma::vec x, arma::mat sd, int T, int L, arma::cube Y, double prior, double tol) {
  arma::vec y(10);
  y(arma::span(0, 8)) = x(arma::span(0, 8)) + arma::chol(sd, "lower") * arma::randn(9);
  y(9) = logfallMQAR1K1(T, L, y(arma::span(0, 8)), Y, prior, tol);
  double A = y(9) - x(9);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

//////////////////
//Spatial QAR1k1//
//////////////////

// [[Rcpp::export]]
double logLikelihoodSQAR1K1(int T, int L, int n, arma::cube Y, 
                            double a1, double b1, double a2, double b2,
                            double gamma, arma::mat RPhi, arma::mat I,
                            double tol) {
  
  arma::cube tauY(T, L, n);
  for (int i = 0; i < n; ++i) {
    for (int t = 0; t < T; ++t) {
      for (int l = 1; l < L; ++l) {
        tauY(t,l,i) = brentQAR1K1(Y(t,l-1,i), Y(t,l,i), a1, b1, a2, b2, tol);
      }
    }
  }
  
  Y.shed_col(L - 1);
  tauY.shed_col(0);
  
  arma::cube logLikelihood = - log(Y % dKumaraswamyCube(tauY, a1, b1) + (1 - Y) % dKumaraswamyCube(tauY, a2, b2));
  
  arma::mat R = gamma * RPhi + (1 - gamma) * I;
  arma::mat Rstar = inv_sympd(R) - I;
  
  arma::mat product(T, L - 1);
  arma::vec vector(n);
  
  for (int t = 0; t < T; ++t) {
    for (int l = 0; l < L - 1; ++l) {
      for (int i = 0; i < n; ++i) {
        vector(i) = R::qnorm(tauY(t,l,i), 0, 1, 1, 0);
      }
      product(t,l) = (vector.t() * Rstar * vector).eval()(0,0);
    }
  }
  
  arma::mat logCopula = - 0.5 * (log_det_sympd(R) + product);
  
  return arma::accu(logLikelihood) + arma::accu(logCopula);
}

// [[Rcpp::export]]
double logfallSQAR1K1(int T, int L, int n, arma::vec x, arma::cube Y, arma::mat RPhi, arma::mat I, double prior, double tol) {
  return logLikelihoodSQAR1K1(T, L, n, Y, exp(x(0)), exp(x(1)), exp(x(2)), exp(x(3)), logitInv(x(4), 0, 1), RPhi, I, tol) +
    logdnorm(x(0), prior) +
    logdnorm(x(1), prior) +
    logdnorm(x(2), prior) +
    logdnorm(x(3), prior) +
    R::dlogis(x(4), 0, 1, 1);
}

// [[Rcpp::export]]
arma::vec rwBmetropolisSQAR1K1(arma::vec x, arma::mat sd, int T, int L, int n, arma::cube Y, arma::mat RPhi, arma::mat I, double prior, double tol) {
  arma::vec y(6);
  y(arma::span(0, 4)) = x(arma::span(0, 4)) + arma::chol(sd, "lower") * arma::randn(5);
  y(5) = logfallSQAR1K1(T, L, n, y(arma::span(0, 4)), Y, RPhi, I, prior, tol);
  double A = y(5) - x(5);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

/////////////////////
//Spatial GP QAR1k1//
/////////////////////

// [[Rcpp::export]]
double logLikelihoodSGPQAR1K1(int T, int L, int n, arma::cube Y, 
                              double gamma, arma::vec GP, 
                              arma::mat RPhi, arma::mat I,
                              double tol) {
  
  arma::cube tauY(T, L, n);
  for (int i = 0; i < n; ++i) {
    for (int t = 0; t < T; ++t) {
      for (int l = 1; l < L; ++l) {
        tauY(t,l,i) = brentQAR1K1(Y(t,l-1,i), Y(t,l,i), GP(i), GP(i + n), GP(i + 2 * n), GP(i + 3 * n), tol);
      }
    }
  }
  
  Y.shed_col(L - 1);
  tauY.shed_col(0);
  
  arma::mat logLikelihood(T, L - 1, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    logLikelihood -= log(Y.slice(i) % dKumaraswamyMat(tauY.slice(i), GP(i), GP(i + n)) + (1 - Y.slice(i)) % dKumaraswamyMat(tauY.slice(i), GP(i + 2 * n), GP(i + 3 * n)));
  }
  
  arma::mat R = gamma * RPhi + (1 - gamma) * I;
  arma::mat Rstar = inv_sympd(R) - I;
  
  arma::mat product(T, L - 1);
  arma::vec vector(n);
  
  for (int t = 0; t < T; ++t) {
    for (int l = 0; l < L - 1; ++l) {
      for (int i = 0; i < n; ++i) {
        vector(i) = R::qnorm(tauY(t,l,i), 0, 1, 1, 0);
      }
      product(t,l) = (vector.t() * Rstar * vector).eval()(0,0);
    }
  }
  
  arma::mat logCopula = - 0.5 * (log_det_sympd(R) + product);
  
  return arma::accu(logLikelihood + logCopula);
}

// [[Rcpp::export]]
double logfallSGPQAR1K1(int T, int L, int n, arma::vec x, arma::cube Y, arma::mat RPhi, arma::mat Rinv, arma::mat I, arma::vec prior, double tol) {
  return logLikelihoodSGPQAR1K1(T, L, n, Y, logitInv(x(0), 0, 1), exp(x(arma::span(9, 8 + 4 * n))), RPhi, I, tol) +
    R::dlogis(x(0), 0, 1, 1) +
    logdnorm(x(1), prior(0)) +
    logdnorm(x(2), prior(0)) +
    logdnorm(x(3), prior(0)) +
    logdnorm(x(4), prior(0)) +
    logdnorm(x(5), prior(1)) +
    logdnorm(x(6), prior(1)) +
    logdnorm(x(7), prior(1)) +
    logdnorm(x(8), prior(1)) +
    logdmvnorm(x(arma::span(9, 8 + n)), x(1), exp(x(5)), Rinv, n) +
    logdmvnorm(x(arma::span(9 + n, 8 + 2 * n)), x(2), exp(x(6)), Rinv, n) +
    logdmvnorm(x(arma::span(9 + 2 * n, 8 + 3 * n)), x(3), exp(x(7)), Rinv, n) +
    logdmvnorm(x(arma::span(9 + 3 * n, 8 + 4 * n)), x(4), exp(x(8)), Rinv, n);
}

// [[Rcpp::export]]
arma::vec rwBmetropolisSGPQAR1K1(arma::vec x, arma::mat sd, int T, int L, int n, arma::cube Y, arma::mat RPhi, arma::mat Rinv, arma::mat I, arma::vec prior, double tol) {
  arma::vec y(10 + 4 * n);
  y(arma::span(0, 8 + 4 * n)) = x(arma::span(0, 8 + 4 * n)) + arma::chol(sd, "lower") * arma::randn(9 + 4 * n);
  y(9 + 4 * n) = logfallSGPQAR1K1(T, L, n, y(arma::span(0, 8 + 4 * n)), Y, RPhi, Rinv, I, prior, tol);
  double A = y(9 + 4 * n) - x(9 + 4 * n);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

//////////////////
//    KX2006    //
//////////////////

// [[Rcpp::export]]
double fKX2006(double x, double y1, double y2, 
               double mu, double sigma, double gamma0, double gamma1) {
  double rho = gamma0 + gamma1 * x;
  return(y2 - (mu + sigma * R::qnorm(x, 0 ,1, 1, 0) + (rho < 1 ? rho : 1) * y1));
}

// [[Rcpp::export]]
double brentKX2006(double y1, double y2, 
                   double mu, double sigma, double gamma0, double gamma1,
                   double tol) {
  int n = 1;
  double a = 0;
  double b = 1;
  double d;
  double aux;
  double fa = INFINITY;
  double fb = -INFINITY;
  double r = a;
  double s;
  double fr = fa;
  double fs;
  bool mflag = true;
  bool isInf;
  while (n < 10001) {
    isInf = (abs(fa) == INFINITY || abs(fb) == INFINITY);
    if (fa != fr && fb != fr && !isInf) {
      s = (a*fb*fr)/((fa-fb)*(fa-fr))+(b*fa*fr)/((fb-fr)*(fb-fa))+(r*fa*fb)/((fr-fa)*(fr-fb));
    } else if (!isInf) {
      s = b-(fb*(b-a))/(fb-fa);
    }
    if (isInf ||
        (((3 * a + b) / 4) > s || s > b) ||
        (mflag && abs(s - b) >= abs(b - r) / 2) ||
        (!mflag && abs(s - b) >= abs(r - d) / 2)
    ) {
      s = (b + a) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = fKX2006(s, y1, y2, mu, sigma, gamma0, gamma1);
    d = r;
    r = b;
    fr = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fa) < abs(fb)) {
      aux = a;
      a = b;
      b = aux;
      aux = fa;
      fa = fb;
      fb = aux;
    }
    if (abs(b - a) < tol || fs == 0) {
      return s;
    } 
    ++n;
  }
  Rcpp::warning("Warning: The rootfinder did not converge \n");
  return s;
}

// [[Rcpp::export]]
double logLikelihoodKX2006(int T, int L, arma::mat Y, 
                           double mu, double sigma, double gamma0, double gamma1,
                           double tol) {
  
  double tauY;
  double rho;
  double logLikelihood = 0;
  
  for (int t = 0; t < T; ++t) {
    for (int l = 1; l < L; ++l) {
      tauY = brentKX2006(Y(t,l-1), Y(t,l), mu, sigma, gamma0, gamma1, tol);
      rho = gamma0 + gamma1 * tauY;
      logLikelihood -= log(sigma / R::dnorm(R::qnorm(tauY, 0 ,1, 1, 0), 0, 1, 0) + (rho < 1 ? gamma1 : 0) * Y(t,l-1));
    }
  }
  
  return logLikelihood;
}

// [[Rcpp::export]]
double logfallKX2006(int T, int L, arma::vec x, arma::mat Y, arma::vec prior, double tol) {
  return logLikelihoodKX2006(T, L, Y, x(0), exp(x(1)), logitInv(x(2), 0, 1), exp(x(3)), tol) +
    logdnorm(x(0), prior(0)) +
    logdnorm(x(1), prior(1)) +
    R::dlogis(x(2), 0, 1, 1) +
    logdnorm(x(3), prior(2));
}

// [[Rcpp::export]]
arma::vec rwBmetropolisKX2006(arma::vec x, arma::mat sd, int T, int L, arma::mat Y, arma::vec prior, double tol) {
  arma::vec y(5); 
  y(arma::span(0, 3)) = x(arma::span(0, 3)) + arma::chol(sd, "lower") * arma::randn(4);
  y(4) = logfallKX2006(T, L, y(arma::span(0, 3)), Y, prior, tol);
  while (y(4) == INFINITY) {
    y(arma::span(0, 3)) = x(arma::span(0, 3)) + arma::chol(sd, "lower") * arma::randn(4);
    y(4) = logfallKX2006(T, L, y(arma::span(0, 3)), Y, prior, tol);
  }
  double A = y(4) - x(4);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}

