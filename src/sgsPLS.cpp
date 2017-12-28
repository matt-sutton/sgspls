#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec softthresh(arma::vec Mx, double lambda){
  int n = Mx.size();
  arma::vec out = arma::max(abs(Mx) -lambda, arma::zeros(n));
  out = out%arma::sign(Mx);
  return out;
}


// [[Rcpp::export]]
double lambdazerosubgroup( double lambda, arma::vec Mx, arma::vec subgroupX, double alpha1, double alpha2){

  arma::vec sgUnique = unique(subgroupX);
  int numSG = sgUnique.size();
  int pk = Mx.size();

  NumericVector sgVal(numSG);

  for (int ka = 0; ka < numSG; ka++){

    arma::umat sgind = arma::find( subgroupX == sgUnique(ka) );
    int pka = sgind.size();

    arma::vec res1 = softthresh(Mx.elem(sgind), (1-alpha1-alpha2)*lambda*0.5);
    //Rcpp::Rcout << softthresh(Mx.elem(sgind), (1-alpha1-alpha2)*lambda*0.5) << std::endl;

    double res2 = norm(res1,2) - alpha2*lambda*sqrt(pka);

    sgVal[ka] = pow(std::max(pow(10,-6), res2),2);
  }
  double out = sum(sgVal) - pow(alpha1*lambda,2)*pk;

  // sets the min to set group to zero exactly
  if(std::abs(out) <= pow(10,-6)) {
    out = std::abs(out) - 0.1;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec sgssoftthresh(arma::vec Mx,arma::vec subgroupX, double lambda, double alpha1, double alpha2, double g){
  int pk = Mx.size();
  arma::vec res(pk);
  arma::vec sgUnique = unique(subgroupX);
  int numSG = sgUnique.size();

  for (int ka = 0; ka < numSG; ka++){
    arma::umat sgind = arma::find( subgroupX == sgUnique(ka));
    int pka = sgind.size();

    arma::vec stsubgroup = softthresh(Mx.elem(sgind), (1-alpha1-alpha2)*lambda*0.5);

    if(alpha2*lambda*sqrt(pka) >= norm(stsubgroup)){

      res.elem(sgind) = arma::zeros(pka);

    } else {
      double h = norm(stsubgroup) - alpha2*lambda*sqrt(pka);
      res.elem(sgind) = stsubgroup*((g*h-alpha1*alpha2*pow(lambda,2))/((alpha1*lambda*sqrt(pk)+g)*(alpha2*lambda*sqrt(pka)+h)));
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::vec updataU( arma::vec Mx, arma::vec groupX, arma::vec subgroupX, double lambda, double alpha1, double alpha2){

  arma::vec gUnique = unique(groupX);
  int numG = gUnique.size();
  int p = Mx.size();
  arma::vec out(p);

  for (int k = 0; k < numG; k++){

    arma::umat gind = arma::find( groupX == gUnique(k) );
    int pk = gind.size();

    arma::vec Mxg = Mx.elem(gind);
    arma::vec subgroupXg = subgroupX.elem(gind);

    if ( lambdazerosubgroup(lambda, Mx.elem(gind), subgroupX.elem(gind), alpha1, alpha2) < pow(10,-2) ) {
      out.elem(gind) = arma::zeros(pk);
    } else {

      double g = norm(softthresh(Mxg, (1-alpha1-alpha2)*lambda/2),2) - alpha1*lambda*sqrt(pk);
      out.elem(gind) = sgssoftthresh(Mxg,subgroupXg, lambda,alpha1, alpha2, g);
    }
  }
  return out;
}

