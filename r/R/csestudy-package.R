#' csestudy: Time-Series Inference for Cross-Sectional Event Studies
#'
#' R implementation of the methodology of Cohn, Johnson, Liu and Wardlaw
#' (2026), "Past is Prologue: Inference from the Cross Section of Returns
#' Around an Event," *Journal of Financial Economics* 180, 104278
#' \doi{10.1016/j.jfineco.2026.104278}.
#'
#' The package mirrors the Stata (`csestudy.ado` / `csestudy.mata`) and Python
#' (`CSEventStudy`) implementations shipped in the same repository and reproduces
#' their numerical output on the shared sample data.
#'
#' @seealso [csestudy()]
#' @keywords internal
"_PACKAGE"

#' @importFrom stats var pt complete.cases setNames
NULL
