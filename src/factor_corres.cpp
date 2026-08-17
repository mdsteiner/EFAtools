#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <string>
using namespace Rcpp;

// Join the 1-based factor indices into a single key, separated so that the sets
// do not collide for ten or more factors (e.g. {1, 2} -> "1-2", not "12" == {12}).
static std::string join_corres(const Rcpp::CharacterVector& pos) {
  std::string out;
  for (int k = 0; k < pos.size(); k++) {
    if (k > 0) out += "-";
    out += Rcpp::as<std::string>(pos[k]);
  }
  return out;
}

//' Compute number of non-matching indicator-to-factor correspondences
//'
//' @param x numeric matrix. A matrix of pattern coefficients.
//' @param y numeric matrix. A second matrix of coefficients.
//' @param thresh numeric. The threshold to classify a pattern coefficient as substantial.
//'
//' @keywords internal
// [[Rcpp::export(.factor_corres)]]
Rcpp::List factor_corres(NumericMatrix x,
                        NumericMatrix y,
                        double thresh = 0.3) {

 if (x.nrow() != y.nrow() || x.ncol() != y.ncol()) {
   Rcpp::stop("x and y must have the same dimensions.");
 }

 IntegerVector x_corres;
 IntegerVector y_corres;
 StringVector x_corres_cross;
 StringVector y_corres_cross;
 int n_comparable = 0;
 int diff_corres = 0;
 int diff_corres_cross = 0;

 // loop through the rows / indicators to find the corresponding factor
 for (int i = 0; i < x.nrow(); i++) {

   Rcpp::CharacterVector x_pos;
   Rcpp::CharacterVector y_pos;
   int x_primary = 0;
   int y_primary = 0;
   double x_max = -std::numeric_limits<double>::infinity();
   double y_max = -std::numeric_limits<double>::infinity();
   int n_shared = 0;

   // Correspondences compare the same factor positions in both matrices. A cell
   // missing on either side therefore removes that position from both row-wise
   // classifications when na.rm = TRUE reaches this helper.
   for (int jj = 0; jj < x.ncol(); ++jj) {
     const double xv = x(i, jj);
     const double yv = y(i, jj);
     if (std::isnan(xv) || std::isnan(yv)) continue;

     ++n_shared;
     const double ax = std::abs(xv);
     const double ay = std::abs(yv);
     if (ax > x_max) {
       x_max = ax;
       x_primary = jj + 1;
     }
     if (ay > y_max) {
       y_max = ay;
       y_primary = jj + 1;
     }
     if (ax >= thresh) x_pos.push_back(std::to_string(jj + 1));
     if (ay >= thresh) y_pos.push_back(std::to_string(jj + 1));
   }

   if (n_shared == 0) {
     x_corres.push_back(NA_INTEGER);
     y_corres.push_back(NA_INTEGER);
     x_corres_cross.push_back(NA_STRING);
     y_corres_cross.push_back(NA_STRING);
     continue;
   }

   if (x_pos.size() == 0) {
     x_primary = 0;
     x_pos.push_back("0");
   }
   if (y_pos.size() == 0) {
     y_primary = 0;
     y_pos.push_back("0");
   }

   x_corres.push_back(x_primary);
   y_corres.push_back(y_primary);

   x_corres_cross.push_back(join_corres(x_pos));
   y_corres_cross.push_back(join_corres(y_pos));

   ++n_comparable;
   if (x_primary != y_primary) ++diff_corres;
   if (join_corres(x_pos) != join_corres(y_pos)) ++diff_corres_cross;
 }

 if (n_comparable == 0) {
   diff_corres = NA_INTEGER;
   diff_corres_cross = NA_INTEGER;
 }

 return Rcpp::List::create(
   Rcpp::Named("x_corres") = x_corres,
   Rcpp::Named("y_corres") = y_corres,
   Rcpp::Named("diff_corres") = diff_corres,
   Rcpp::Named("x_corres_cross") = x_corres_cross,
   Rcpp::Named("y_corres_cross") = y_corres_cross,
   Rcpp::Named("diff_corres_cross") = diff_corres_cross
 );
}
