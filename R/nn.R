.inn_fit <- 
function (d, h, size, decay, seeds, maxit = 2000L) 
{
    cov <- names(.clean_map)
    x <- as.matrix(d[cov])
    center <- colMeans(x)
    scales <- apply(x, 2, sd)
    stopifnot(all(scales > 0))
    x <- scale(x, center, scales)
    a <- make_Dt_weights(d, h, eps_G = 0.001, w_cap = 20)
    w <- a$weights * nrow(d)/sum(a$weights)
    models <- lapply(seeds, function(seed) {
        set.seed(seed)
        o <- nnet::nnet(x = x, y = a$D_t, weights = w, size = size, entropy = TRUE, linout = FALSE, skip = FALSE, 
            decay = decay, maxit = maxit, trace = FALSE, rang = 0.1, abstol = 1e-08, reltol = 1e-09, 
            MaxNWts = 2000)
        stopifnot(o$convergence == 0, all(is.finite(o$wts)), is.finite(o$value))
        o
    })
    structure(list(models = models, covariates = cov, center = center, scale = scales, horizon = h, size = size, 
        decay = decay, seeds = seeds, source_type = "ipcw_nn"), class = "ipcw_nn_source")
}

.inn_prob <- 
function (o, d) 
{
    stopifnot(requireNamespace("nnet", quietly = TRUE))
    x <- scale(as.matrix(d[o$covariates]), o$center, o$scale)
    p <- rowMeans(vapply(o$models, function(m) as.numeric(predict(m, x, type = "raw")), numeric(nrow(d))))
    stopifnot(all(is.finite(p)), all(p >= 0), all(p <= 1))
    p
}

