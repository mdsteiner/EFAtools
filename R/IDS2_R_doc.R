#' IDS2_R
#'
#' A matrix containing the bivariate correlations of the IDS-2 subtests. Details can be found in Grieder & Grob (2019).
#'
#'
#' @format A 14 x 14 matrix of bivariate correlations
#' \describe{
#'   \item{id}{(numeric) - A unique identifier (1 through 150000).}
#'   \item{age}{(numeric) - A peson's age in years.}
#'   \item{income}{(numeric) - Income of the household a person lives in.}
#'   \item{weight}{(numeric) - A person's weight in kg.}
#'   \item{height}{(numeric) - A person's height in cm.}
#'   \item{children}{(numeric) - The number of children a person has.}
#'   \item{happiness}{(numeric) - How happy a person is on a scale from 1 to 10.}
#'   \item{fitness}{(numeric) - A person's fitness level rated from 1 to 10.}
#'   \item{food}{(numeric) - How much (CHF) a person spends on food per month.}
#'   \item{alcohol}{(numeric) - How much (CHF) a person spends on alcohol beverages per month.}
#'   \item{tattoos}{(numeric) - The number of tattoos a person has.}
#'   \item{rhine}{(numeric) - How often a person goes swimming in the rhine per month.}
#'   \item{datause}{(numeric) - How many times a day a person checks an app on the mobile phone.}
#'   \item{consultations}{(numeric) - The number of times a person's sees a doctor per year.}
#'   \item{hiking}{(numeric) - How many hours a person spends hiking per year.}
#'   \item{sex}{(character) - A person's sex. Either "male" or "female".}
#'   \item{education}{(character) - The highest completed degree/ education a person has. Levels are "obligatory_school", "apprenticeship", "SEK_II", "SEK_III".}
#'   \item{confession}{(character) - A person's confession. Levels are "catholic", "confessionless", "evangelical-reformed", "muslim", "other".}
#'   \item{fasnacht}{(character) - Whether a person actively (i.e. in a Gugge, Clique or similar) participates at fasnacht. Levels are "yes" and "no".}
#'   \item{eyecor}{(character) - Whether a person needs eye correction. Levels are "yes" and "no".}
#'  }
#' @source Grieder, S., & Grob, A. (2019). Exploratory factor analysis of the intelligence and development scales–2 (IDS-2): Implications for theory and practice. Submitted Manuscript.
"IDS2_R"
