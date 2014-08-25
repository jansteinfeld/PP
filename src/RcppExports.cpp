// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// P_4pl
NumericVector P_4pl(NumericVector delta, double alpha, double theta, double la, double ua);
RcppExport SEXP PP_P_4pl(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP laSEXP, SEXP uaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< double >::type la(laSEXP );
        Rcpp::traits::input_parameter< double >::type ua(uaSEXP );
        NumericVector __result = P_4pl(delta, alpha, theta, la, ua);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// P_4pl4wle
NumericVector P_4pl4wle(NumericVector delta, double alpha, double theta, double la, double ua);
RcppExport SEXP PP_P_4pl4wle(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP laSEXP, SEXP uaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< double >::type la(laSEXP );
        Rcpp::traits::input_parameter< double >::type ua(uaSEXP );
        NumericVector __result = P_4pl4wle(delta, alpha, theta, la, ua);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// r_huber_4pl
double r_huber_4pl(NumericVector delta, double alpha, double theta, double la, double ua, double H);
RcppExport SEXP PP_r_huber_4pl(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP laSEXP, SEXP uaSEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< double >::type la(laSEXP );
        Rcpp::traits::input_parameter< double >::type ua(uaSEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        double __result = r_huber_4pl(delta, alpha, theta, la, ua, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L4pl
NumericMatrix L4pl(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, bool map, NumericVector mu, NumericVector sigma2);
RcppExport SEXP PP_L4pl(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP mapSEXP, SEXP muSEXP, SEXP sigma2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< bool >::type map(mapSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        NumericMatrix __result = L4pl(awm, DELTA, ALPHA, CS, DS, THETA, map, mu, sigma2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L4pl_wle
NumericMatrix L4pl_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA);
RcppExport SEXP PP_L4pl_wle(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        NumericMatrix __result = L4pl_wle(awm, DELTA, ALPHA, CS, DS, THETA);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L4pl_robust
NumericMatrix L4pl_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, double H);
RcppExport SEXP PP_L4pl_robust(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        NumericMatrix __result = L4pl_robust(awm, DELTA, ALPHA, CS, DS, THETA, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// NR_4PL
List NR_4PL(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, String wm, int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H);
RcppExport SEXP PP_NR_4PL(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP wmSEXP, SEXP maxstepsSEXP, SEXP exacSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< String >::type wm(wmSEXP );
        Rcpp::traits::input_parameter< int >::type maxsteps(maxstepsSEXP );
        Rcpp::traits::input_parameter< double >::type exac(exacSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        List __result = NR_4PL(awm, DELTA, ALPHA, CS, DS, THETA, wm, maxsteps, exac, mu, sigma2, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// P_gpcm
double P_gpcm(NumericVector delta, double alpha, double theta, int resp);
RcppExport SEXP PP_P_gpcm(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP respSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< int >::type resp(respSEXP );
        double __result = P_gpcm(delta, alpha, theta, resp);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// r_huber_gpcm
double r_huber_gpcm(NumericVector delta, double alpha, double theta, double H);
RcppExport SEXP PP_r_huber_gpcm(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        double __result = r_huber_gpcm(delta, alpha, theta, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L12gpcm
NumericMatrix L12gpcm(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector THETA, NumericVector mu, NumericVector sigma2, bool map);
RcppExport SEXP PP_L12gpcm(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP THETASEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP mapSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< bool >::type map(mapSEXP );
        NumericMatrix __result = L12gpcm(awm, DELTA, ALPHA, THETA, mu, sigma2, map);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Pcorr1_gpcm
NumericVector Pcorr1_gpcm(NumericVector delta, double alpha, double theta);
RcppExport SEXP PP_Pcorr1_gpcm(SEXP deltaSEXP, SEXP alphaSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< double >::type theta(thetaSEXP );
        NumericVector __result = Pcorr1_gpcm(delta, alpha, theta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L12gpcm_wle
NumericMatrix L12gpcm_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector THETA);
RcppExport SEXP PP_L12gpcm_wle(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP THETASEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        NumericMatrix __result = L12gpcm_wle(awm, DELTA, ALPHA, THETA);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// L12gpcm_robust
NumericMatrix L12gpcm_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector THETA, double H);
RcppExport SEXP PP_L12gpcm_robust(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP THETASEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        NumericMatrix __result = L12gpcm_robust(awm, DELTA, ALPHA, THETA, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// NR_GPCM
List NR_GPCM(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector THETA, String wm, int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H);
RcppExport SEXP PP_NR_GPCM(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP THETASEXP, SEXP wmSEXP, SEXP maxstepsSEXP, SEXP exacSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< String >::type wm(wmSEXP );
        Rcpp::traits::input_parameter< int >::type maxsteps(maxstepsSEXP );
        Rcpp::traits::input_parameter< double >::type exac(exacSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        List __result = NR_GPCM(awm, DELTA, ALPHA, THETA, wm, maxsteps, exac, mu, sigma2, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Lgpcm4pl_mle
NumericMatrix Lgpcm4pl_mle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, CharacterVector model, NumericVector mu, NumericVector sigma2, bool map);
RcppExport SEXP PP_Lgpcm4pl_mle(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP modelSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP mapSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< bool >::type map(mapSEXP );
        NumericMatrix __result = Lgpcm4pl_mle(awm, DELTA, ALPHA, CS, DS, THETA, model, mu, sigma2, map);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Lgpcm4pl_wle
NumericMatrix Lgpcm4pl_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, CharacterVector model);
RcppExport SEXP PP_Lgpcm4pl_wle(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP modelSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP );
        NumericMatrix __result = Lgpcm4pl_wle(awm, DELTA, ALPHA, CS, DS, THETA, model);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Lgpcm4pl_robust
NumericMatrix Lgpcm4pl_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, CharacterVector model, double H);
RcppExport SEXP PP_Lgpcm4pl_robust(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP modelSEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        NumericMatrix __result = Lgpcm4pl_robust(awm, DELTA, ALPHA, CS, DS, THETA, model, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// NR_mixed
List NR_mixed(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector CS, NumericVector DS, NumericVector THETA, CharacterVector model, String wm, int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H);
RcppExport SEXP PP_NR_mixed(SEXP awmSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP CSSEXP, SEXP DSSEXP, SEXP THETASEXP, SEXP modelSEXP, SEXP wmSEXP, SEXP maxstepsSEXP, SEXP exacSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP HSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type CS(CSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type DS(DSSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type THETA(THETASEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP );
        Rcpp::traits::input_parameter< String >::type wm(wmSEXP );
        Rcpp::traits::input_parameter< int >::type maxsteps(maxstepsSEXP );
        Rcpp::traits::input_parameter< double >::type exac(exacSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< double >::type H(HSEXP );
        List __result = NR_mixed(awm, DELTA, ALPHA, CS, DS, THETA, model, wm, maxsteps, exac, mu, sigma2, H);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Likgpcm
NumericVector Likgpcm(IntegerVector awv, NumericMatrix DELTA, NumericVector ALPHA, NumericVector nodes);
RcppExport SEXP PP_Likgpcm(SEXP awvSEXP, SEXP DELTASEXP, SEXP ALPHASEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type awv(awvSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type DELTA(DELTASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ALPHA(ALPHASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type nodes(nodesSEXP );
        NumericVector __result = Likgpcm(awv, DELTA, ALPHA, nodes);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sim_4pl
IntegerMatrix sim_4pl(NumericVector beta, NumericVector alpha, NumericVector lowerA, NumericVector upperA, NumericVector theta);
RcppExport SEXP PP_sim_4pl(SEXP betaSEXP, SEXP alphaSEXP, SEXP lowerASEXP, SEXP upperASEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type lowerA(lowerASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type upperA(upperASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP );
        IntegerMatrix __result = sim_4pl(beta, alpha, lowerA, upperA, theta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ansol
NumericMatrix ansol(IntegerMatrix awm, IntegerVector maxsc);
RcppExport SEXP PP_ansol(SEXP awmSEXP, SEXP maxscSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type awm(awmSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type maxsc(maxscSEXP );
        NumericMatrix __result = ansol(awm, maxsc);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
