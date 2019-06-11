R -e 'Rcpp::compileAttributes()' && R -e 'devtools::document()' && R -e 'devtools::install(quick=TRUE)'
