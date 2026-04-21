#' PP: Estimation of Person Parameters and Person Fit
#'
#' The PP package provides functions for estimating person parameters for the
#' 1PL, 2PL, 3PL, and 4PL models, as well as the generalized partial credit
#' model (GPCM). Supported estimation methods include maximum likelihood (ML),
#' weighted likelihood (WL; Warm, 1989), maximum a posteriori (MAP), expected
#' a posteriori (EAP), and robust estimation.
#'
#' In addition, the package includes routines for person-fit analysis, including
#' infit, outfit, lz, and lzstar statistics. The implementation is designed for
#' efficient computation and includes compiled code for fast estimation.
#'
#' For an introduction to the package workflow, see the package vignettes.
#'
#' @author Jan Steinfeld and Manuel Reif
#' @seealso [PPass()], [PP_gpcm()], [PP_4pl()], [PPall()], [Pfit()]
#' @references
#' Barton, M. A., & Lord, F. M. (1981). An upper asymptote for the
#' three-parameter logistic item-response model.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#' examinee's ability. In F. M. Lord & M. R. Novick (Eds.), \emph{Statistical
#' theories of mental test scores}. Reading, MA: Addison-Wesley.
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness
#' measurement with polychotomous item response models and standardized indices.
#' \emph{British Journal of Mathematical and Statistical Psychology},
#' \bold{38}(1), 67--86.
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an EM
#' algorithm. \emph{Applied Psychological Measurement}, \bold{16}, 159--176.
#'
#' Samejima, F. (1993). An approximation of the bias function of the maximum
#' likelihood estimate of a latent variable for the general case where the item
#' responses are discrete. \emph{Psychometrika}, \bold{58}, 119--138.
#'
#' Snijders, T. B. (2001). Asymptotic null distribution of person fit statistics
#' with estimated person parameter. \emph{Psychometrika}, \bold{66}(3), 331--342.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. \emph{Psychometrika}, \bold{54}, 427--450.
#'
#' Wright, B. D., & Masters, G. N. (1990). Computation of OUTFIT and INFIT
#' statistics. \emph{Rasch Measurement Transactions}, \bold{3}(4), 84--85.
#'
#' Yen, Y.-C., Ho, R.-G., Liao, W.-W., Chen, L.-J., & Kuo, C.-C. (2012). An
#' empirical evaluation of the slip correction in the four parameter logistic
#' models with computerized adaptive testing. \emph{Applied Psychological
#' Measurement}, \bold{36}, 75--87.
#'
#' @examples
#' set.seed(1522)
#'
#' diffpar <- seq(-3, 3, length = 12)
#' sl <- round(runif(12, 0.5, 1.5), 2)
#' la <- round(runif(12, 0, 0.25), 2)
#' ua <- round(runif(12, 0.8, 1), 2)
#'
#' awm <- matrix(sample(0:1, 10 * 12, replace = TRUE), ncol = 12)
#'
#' res3plmle <- PP_4pl(
#'   respm = awm,
#'   thres = diffpar,
#'   slopes = sl,
#'   lowerA = la,
#'   type = "mle"
#' )
#'
#' res3plwle <- PP_4pl(
#'   respm = awm,
#'   thres = diffpar,
#'   slopes = sl,
#'   lowerA = la,
#'   type = "wle"
#' )
#'
#' res3plmap <- PP_4pl(
#'   respm = awm,
#'   thres = diffpar,
#'   slopes = sl,
#'   lowerA = la,
#'   type = "map"
#' )
#'
#' res3plmlepfit <- Pfit(
#'   respm = awm,
#'   pp = res3plmle,
#'   fitindices = c("infit", "outfit")
#' )
#'
#' out <- PPass(
#'   respdf = data.frame(awm),
#'   thres = diffpar,
#'   items = "all",
#'   mod = c("1PL"),
#'   fitindices = c("lz", "lzstar", "infit", "outfit")
#' )
#'
#' @keywords internal
"_PACKAGE"