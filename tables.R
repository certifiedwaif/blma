data_sets <- c("Kakadu", "UScrime", "comCrime", "eyeData")
priors <- c("BIC", "ZE", "liang_g1", "liang_g2", "liang_g_n_approx", "liang_g_n_quad", "robust_bayarri1", "robust_bayarri2", "zellner_siow_gauss_laguerre")
for (data_set in data_sets) {
	for (prior in prior) {
		cat(data_set, prior, "\n")
	}
}
