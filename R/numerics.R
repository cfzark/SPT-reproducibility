# Selected historical numerical routines. No historical runners are loaded.
make_scaler <- 
function (train_data, cont_covariates = c("X3", "X4", "X5", "X6", "X7")) 
{
    scaler <- lapply(cont_covariates, function(v) {
        s <- stats::sd(train_data[[v]], na.rm = TRUE)
        if (is.na(s) || s < 1e-08) 
            s <- 1
        list(mean = mean(train_data[[v]], na.rm = TRUE), sd = s)
    })
    names(scaler) <- cont_covariates
    scaler
}

apply_scaler <- 
function (data, scaler, suffix = "_sc") 
{
    out <- data
    for (v in names(scaler)) {
        out[[paste0(v, suffix)]] <- (out[[v]] - scaler[[v]]$mean)/scaler[[v]]$sd
    }
    out
}

estimate_weights <- 
function (data, t, Y_col = "Y", delta_col = "delta", eps_G = 0.001, w_cap = Inf) 
{
    Y_vals <- data[[Y_col]]
    delta_vals <- data[[delta_col]]
    km_fit <- survival::survfit(survival::Surv(Y_vals, 1 - delta_vals) ~ 1, data = data)
    time_point <- pmin(Y_vals, t)
    Ghat <- vapply(time_point, function(x) {
        s <- summary(km_fit, times = x)$surv
        if (length(s) == 0 || is.na(s)) 
            eps_G
        else max(s, eps_G)
    }, numeric(1))
    indicator <- ifelse((Y_vals > t) | ((Y_vals <= t) & (delta_vals == 1)), 1, 0)
    weights <- indicator/Ghat
    if (is.finite(w_cap)) 
        weights <- pmin(weights, w_cap)
    weights
}

make_Dt_weights <- 
function (data, t, Y_col = "Y", delta_col = "delta", eps_G = 0.001, w_cap = 20) 
{
    list(D_t = as.integer(data[[Y_col]] <= t & data[[delta_col]] == 1), weights = estimate_weights(data, 
        t, Y_col, delta_col, eps_G, w_cap))
}

.guard_km_fit <- 
function (data) 
{
    survival::survfit(survival::Surv(data$Y, 1 - data$delta) ~ 1, data = data)
}

.guard_km_eval <- 
function (fit, times) 
{
    times <- pmax(as.numeric(times), 0)
    out <- rep(1, length(times))
    positive <- times > 0
    if (any(positive)) {
        vals <- as.numeric(summary(fit, times = times[positive], extend = TRUE)$surv)
        if (length(vals) == sum(positive)) 
            out[positive] <- vals
    }
    out[!is.finite(out)] <- 0
    pmin(pmax(out, 0), 1)
}

.frozen_km_weights <- 
function (time, event, h, multiplicity = rep(1, length(time)), eps = 0.001) 
{
    stopifnot(length(time) == length(event), length(time) == length(multiplicity), all(multiplicity >= 
        0), sum(multiplicity) > 0)
    tt <- sort(unique(time))
    idx <- match(time, tt)
    totals <- as.numeric(rowsum(multiplicity, idx, reorder = TRUE))
    cens <- as.numeric(rowsum(multiplicity * (1 - event), idx, reorder = TRUE))
    risk <- rev(cumsum(rev(totals)))
    haz <- ifelse(risk > 0, cens/pmax(risk, 1), 0)
    g <- cumprod(1 - haz)
    lookup <- c(1, g)[findInterval(pmin(time, h), tt) + 1L]
    observed <- time > h | (time <= h & event == 1)
    as.numeric(observed)/pmax(lookup, eps)
}

.usim_auc <- 
function (score, cases, controls) 
{
    stopifnot(length(score) == length(cases), length(score) == length(controls), all(is.finite(score)))
    o <- order(score)
    s <- score[o]
    ca <- cases[o]
    co <- controls[o]
    group <- cumsum(c(TRUE, diff(s) != 0))
    a <- as.numeric(rowsum(ca, group, reorder = FALSE))
    b <- as.numeric(rowsum(co, group, reorder = FALSE))
    sum(a * (cumsum(b) - b/2))/(sum(a) * sum(b))
}

.usim_entropy <- 
function (p) 
{
    a <- b <- numeric(length(p))
    ix <- p > 0
    a[ix] <- -p[ix] * log(p[ix])
    ix <- p < 1
    b[ix] <- -(1 - p[ix]) * log1p(-p[ix])
    a + b
}

.usim_eps <- 
function (u, distribution) 
{
    if (distribution == "loglogistic") 
        qlogis(u)
    else if (distribution == "weibull") 
        log(-log1p(-u)) + 0.577215664901533
    else stop("Unknown distribution")
}

.usim_p <- 
function (x, h, setting, domain, distribution, sigma, b = 0.8) 
{
    z <- (log(h) - log(36) + .usim_g(x, setting, domain, b))/sigma
    if (distribution == "loglogistic") 
        plogis(z)
    else if (distribution == "weibull") 
        -expm1(-exp(z - 0.577215664901533))
    else stop("Unknown distribution")
}

.usim_x <- 
function (n, domain = "source") 
{
    stopifnot(domain %in% c("source", "target"))
    x <- data.frame(X1 = rbinom(n, 1, if (domain == "source") 
        0.4
    else 0.5), X2 = rbinom(n, 1, 0.5))
    for (j in 3:7) x[[paste0("X", j)]] <- runif(n, -sqrt(3), sqrt(3))
    if (domain == "target") {
        u <- (x$X3/sqrt(3) + 1)/2
        d <- 0.5
        y <- 2 * u/(1 - d + sqrt((1 - d)^2 + 4 * d * u))
        x$X3 <- sqrt(3) * (2 * y - 1)
    }
    x
}

.ext_shared_g <- 
function (x, setting, domain, b = 0.8) 
{
    stopifnot(setting %in% 1:4, domain %in% c("source", "target"))
    a <- 0.6 * (x$X1 - 0.4)/sqrt(0.24) - 0.45 * (x$X2 - 0.5)/0.5 + 0.75 * x$X3 - 0.5 * x$X5 + 0.35 * 
        x$X7 + 0.55 * (x$X4^2 - 1)/sqrt(0.8) + 0.45 * sin(pi * x$X6/sqrt(3))
    u <- (x$X1 - 0.4)/sqrt(0.24) * x$X3
    if (domain == "source" || setting == 1) 
        return(a + b * u)
    if (setting == 2) 
        return(a - 0.55 * x$X3 + 0.85 * x$X5 + b * u)
    if (setting == 3) 
        return(a + 0.6 * b * u)
    a
}

.sg_stream <- 
function (n, domain, seed) 
{
    set.seed(seed)
    x <- .usim_x(n, domain)
    u <- runif(n)
    eu <- rexp(n)
    if (.ext_p > 7L) {
        set.seed(seed + 500000L)
        for (j in 8:.ext_p) x[[paste0("X", j)]] <- runif(n, -sqrt(3), sqrt(3))
    }
    list(x = x, u = u, eu = eu)
}

.sg_draw <- 
function (n, setting, domain, distribution, sigma, lambda, seed, prefix) 
{
    z <- .sg_stream(n, domain, seed)
    tt <- exp(log(36) - .usim_g(z$x, setting, domain) + sigma * .usim_eps(z$u, distribution))
    ct <- z$eu/lambda
    stopifnot(all(is.finite(tt)), all(tt > 0), all(is.finite(ct)), all(ct > 0))
    list(data = data.frame(eid = paste0(prefix, ":", seq_len(n)), row_id = seq_len(n), Y = pmin(tt, ct), 
        delta = as.integer(tt <= ct), z$x), event_time = tt, censor_time = ct, true_risk = NULL)
}

.sg_truth <- 
function (z, hs, setting, domain, distribution, sigma) 
{
    z$true_risk <- sapply(hs, function(h) .usim_p(z$data, h, setting, domain, distribution, sigma))
    g <- .usim_g(z$data, setting, domain)
    direct <- sapply(hs, function(h) {
        eta <- (log(h) - log(36) + g)/sigma
        if (distribution == "loglogistic") 
            plogis(eta)
        else pweibull(h, shape = 1/sigma, scale = exp(log(36) - g + sigma * 0.577215664901533))
    })
    stopifnot(all(is.finite(z$true_risk)), all(z$true_risk >= 0), all(z$true_risk <= 1), max(abs(z$true_risk - 
        direct)) < 1e-12)
    z
}

uro_calibrate_teacher <- 
function (source_lp, y, w, clip = 0.001) 
{
    stopifnot(length(source_lp) == length(y), length(w) == length(y), all(is.finite(source_lp)), all(is.finite(y)), 
        all(is.finite(w)), all(w >= 0), sum(w) > 0, sum(w * y) > 0, sum(w * (1 - y)) > 0)
    score <- function(alpha) sum(w * (plogis(source_lp + alpha) - y))/sum(w)
    lo <- -8
    hi <- 8
    while (score(lo) > 0 && abs(lo) < 1024) lo <- lo * 2
    while (score(hi) < 0 && hi < 1024) hi <- hi * 2
    stopifnot(score(lo) <= 0, score(hi) >= 0)
    alpha <- uniroot(score, c(lo, hi), tol = 1e-12)$root
    raw <- plogis(source_lp + alpha)
    q <- pmin(pmax(raw, clip), 1 - clip)
    residual <- score(alpha)
    stopifnot(abs(residual) < 1e-10)
    list(alpha = alpha, q = q, score = residual, clipped_score = sum(w * (q - y))/sum(w), clipped_fraction = mean(raw < 
        clip | raw > 1 - clip), target_risk = sum(w * y)/sum(w), teacher_mean = mean(q), teacher_weighted_mean = sum(w * 
        q)/sum(w))
}

ucv_validation <- 
function (valid, h, p, training_km) 
{
    p <- pmin(pmax(as.matrix(p), 0.001), 0.999)
    stopifnot(nrow(p) == nrow(valid), all(is.finite(p)))
    y <- as.integer(valid$Y <= h & valid$delta == 1)
    observed <- y == 1 | valid$Y > h
    g <- .guard_km_eval(training_km, pmin(valid$Y, h))
    w <- ifelse(observed, 1/pmax(g, 0.001), 0)
    loss <- w * (-y * log(p) - (1 - y) * log1p(-p))
    stopifnot(all(is.finite(loss)), all(w >= 0))
    list(p = p, y = y, w = w, g = g, loss = loss)
}

se_solver <- 
function (B, c, r, y, w, q, tau, S, theta_start, fixed_rho = NULL, fixed_omega = 0, kkt_tol = 2e-06) 
{
    stopifnot(identical(fixed_omega, 0))
    stopifnot(nrow(B) == length(y), length(c) == ncol(B), length(r) == length(y), all(is.finite(B)), 
        all(is.finite(c)), all(is.finite(r)), all(is.finite(S)), all(w >= 0), sum(w) > 0, all(y >= 0 & 
            y <= 1), all(q > 0 & q < 1), tau >= 0)
    n <- nrow(B)
    k <- ncol(B)
    v <- w + tau
    wy <- w * y + tau * q
    Sc <- drop(S %*% c)
    den <- sum(c * Sc)
    stopifnot(den >= -1e-08)
    zero_rho <- abs(den) <= 1e-14 * max(1, sum(abs(c * Sc)))
    zero_omega <- max(abs(r)) <= 1e-12
    if (zero_omega && is.null(fixed_omega)) 
        fixed_omega <- 0
    R <- chol(crossprod(B, B * (v/4)) + S)
    Q <- backsolve(R, diag(k))
    Z <- B %*% Q
    norm <- sum(w)
    free_omega <- is.null(fixed_omega)
    unpack <- function(z) {
        theta <- drop(Q %*% z[seq_len(k)])
        rho <- if (!is.null(fixed_rho)) 
            fixed_rho
        else if (zero_rho) 
            0
        else min(1, max(0, sum(theta * Sc)/den))
        omega <- if (free_omega) 
            z[k + 1L]
        else fixed_omega
        beta <- theta - rho * c
        eta <- drop(B %*% theta) + omega * r
        list(theta = theta, rho = rho, omega = omega, beta = beta, eta = eta)
    }
    value <- function(z) {
        a <- unpack(z)
        (sum(v * .uro_softplus(a$eta) - wy * a$eta) + sum(a$beta * (S %*% a$beta))/2)/norm
    }
    gradient <- function(z) {
        a <- unpack(z)
        d <- v * plogis(a$eta) - wy
        g <- drop(crossprod(Q, crossprod(B, d) + S %*% a$beta))/norm
        if (free_omega) 
            c(g, sum(r * d)/norm)
        else g
    }
    zbase <- drop(R %*% theta_start)
    starts <- if (free_omega) 
        c(0, 0.5, 1)
    else NA_real_
    attempts <- lapply(starts, function(om) {
        initial <- if (free_omega) 
            c(zbase, om)
        else zbase
        opt <- optim(initial, value, gradient, method = "L-BFGS-B", lower = c(rep(-Inf, k), if (free_omega) 0), 
            upper = c(rep(Inf, k), if (free_omega) 1), control = list(maxit = 3000, factr = 1000, pgtol = 1e-10, 
                lmm = 20))
        a <- unpack(opt$par)
        g <- gradient(opt$par)
        residual <- max(abs(c(g[seq_len(k)], if (free_omega) .uro_boxgrad(g[k + 1], a$omega))))
        c(list(par = opt$par, value = value(opt$par), kkt = residual, code = opt$convergence, message = if (is.null(opt$message)) "" else opt$message), 
            a)
    })
    values <- vapply(attempts, `[[`, numeric(1), "value")
    ok <- vapply(attempts, function(a) a$code == 0L && a$kkt <= kkt_tol, logical(1))
    if (!any(ok)) {
        original_attempts <- attempts
        for (a0 in original_attempts) {
            opt <- nlminb(a0$par, value, gradient, lower = c(rep(-Inf, k), if (free_omega) 0), upper = c(rep(Inf, 
                k), if (free_omega) 1), control = list(iter.max = 10000, eval.max = 20000, rel.tol = 1e-12, 
                x.tol = 1e-12))
            aa <- unpack(opt$par)
            gg <- gradient(opt$par)
            residual <- max(abs(c(gg[seq_len(k)], if (free_omega) .uro_boxgrad(gg[k + 1], aa$omega))))
            attempts[[length(attempts) + 1L]] <- c(list(par = opt$par, value = value(opt$par), kkt = residual, 
                code = opt$convergence, message = paste("nlminb same-objective restart:", opt$message)), 
                aa)
        }
        values <- vapply(attempts, `[[`, numeric(1), "value")
        ok <- vapply(attempts, function(a) a$code == 0L && a$kkt <= kkt_tol, logical(1))
        if (!any(ok)) 
            stop("No converged stationary optimizer start, including nlminb restart")
        warning("Accepted numerical repair only after unchanged objective/code-0/KKT checks")
    }
    accepted <- which(ok)[which.min(values[ok])]
    best <- attempts[[accepted]]
    gap <- best$value - min(values)
    stopifnot(gap <= 1e-10)
    if (any(!ok)) 
        warning("Unsuccessful optimizer start retained; accepted code 0, original KKT, objective gap <=1e-10")
    if (best$code != 0L || best$kkt > kkt_tol) 
        stop(sprintf("SPT optimizer failed: code=%d KKT=%.5g %s", best$code, best$kkt, best$message))
    rho_gradient <- -sum(c * (S %*% best$beta))/norm
    rho_kkt <- if (is.null(fixed_rho)) 
        abs(.uro_boxgrad(rho_gradient, best$rho))
    else 0
    stopifnot(rho_kkt < kkt_tol)
    best$diagnostic <- data.frame(accepted_objective_gap = gap, tau = tau, rho = best$rho, omega = best$omega, 
        objective = best$value, kkt = best$kkt, rho_kkt = rho_kkt, optimizer_code = best$code, rho_penalty_curvature = den/norm, 
        zero_rho = zero_rho, zero_omega = zero_omega, residual_rms = sqrt(mean(r^2)), multistart_objective_spread = diff(range(vapply(attempts, 
            `[[`, numeric(1), "value"))), optimizer_attempts = length(attempts), unsuccessful_attempts = sum(vapply(attempts, 
            function(a) a$code != 0L || a$kkt > kkt_tol, logical(1))))
    best$attempts <- do.call(rbind, lapply(seq_along(attempts), function(i) data.frame(attempt = i, objective = attempts[[i]]$value, 
        kkt = attempts[[i]]$kkt, code = attempts[[i]]$code, message = attempts[[i]]$message)))
    best$gradient <- gradient(best$par)
    best
}

usi_single_fit <- 
function (ctx, cal, tau = 0.5, kkt_tol = 2e-06) 
{
    B <- ctx$B
    c <- ctx$c
    r <- ctx$r
    S <- ctx$S
    w <- ctx$w
    y <- ctx$y
    q <- if (tau > 0) 
        cal$q
    else ctx$q
    v <- w + tau
    wy <- w * y + tau * q
    norm <- sum(w)
    k <- ncol(B)
    R <- chol(crossprod(B, B * (v/4)) + S)
    Q <- backsolve(R, diag(k))
    s <- drop(B %*% c) + r
    unpack <- function(z) {
        beta <- drop(Q %*% z[seq_len(k)])
        gamma <- z[k + 1L]
        theta <- beta + gamma * c
        eta <- drop(B %*% beta) + gamma * s
        list(theta = theta, beta = beta, rho = gamma, omega = gamma, gamma = gamma, eta = eta)
    }
    value <- function(z) {
        a <- unpack(z)
        (sum(v * .uro_softplus(a$eta) - wy * a$eta) + sum(a$beta * (S %*% a$beta))/2)/norm
    }
    gradient <- function(z) {
        a <- unpack(z)
        d <- v * plogis(a$eta) - wy
        Sb <- drop(S %*% a$beta)
        c(drop(crossprod(Q, crossprod(B, d) + Sb)), sum(s * d))/norm
    }
    tries <- lapply(c(0, 0.5, 1), function(gamma) {
        initial <- c(drop(R %*% coef(ctx$anchor$fit)), gamma)
        opt <- optim(initial, value, gradient, method = "L-BFGS-B", lower = c(rep(-Inf, k), 0), upper = c(rep(Inf, 
            k), 1), control = list(maxit = 3000, factr = 1000, pgtol = 1e-10, lmm = 20))
        a <- unpack(opt$par)
        g <- gradient(opt$par)
        residual <- max(abs(c(g[seq_len(k)], .uro_boxgrad(g[k + 1L], a$gamma))))
        c(a, list(value = value(opt$par), kkt = residual, code = opt$convergence, par = opt$par))
    })
    values <- vapply(tries, `[[`, numeric(1), "value")
    ok <- vapply(tries, function(a) a$code == 0L && a$kkt <= kkt_tol, logical(1))
    stopifnot(any(ok))
    accepted <- which(ok)[which.min(values[ok])]
    best <- tries[[accepted]]
    gap <- best$value - min(values)
    stopifnot(gap <= 1e-10)
    if (best$code != 0L || best$kkt > kkt_tol) 
        stop("Single-weight optimizer failed KKT or convergence")
    s <- drop(B %*% c) + r
    stopifnot(max(abs(best$eta - (drop(B %*% best$beta) + best$gamma * s))) < 1e-09)
    ref <- (sum(w * (.uro_softplus(best$eta) - y * best$eta)) + tau * sum(.uro_softplus(best$eta) - q * 
        best$eta) + sum(best$beta * (S %*% best$beta))/2)/norm
    stopifnot(abs(ref - best$value) < 1e-10)
    score <- sum(v * plogis(best$eta) - wy)/sum(v)
    stopifnot(abs(score) < 1e-05)
    best$tau <- tau
    best$diagnostic <- data.frame(accepted_objective_gap = gap, tau = tau, gamma = best$gamma, objective = best$value, 
        kkt = best$kkt, optimizer_code = best$code, candidate_intercept_score = score, unsuccessful_attempts = sum(vapply(tries, 
            function(z) z$code != 0L || z$kkt > kkt_tol, logical(1))))
    if (best$diagnostic$unsuccessful_attempts > 0) 
        warning("Single-weight unsuccessful start retained; accepted solution passed KKT")
    best$attempts <- do.call(rbind, lapply(seq_along(tries), function(i) data.frame(attempt = i, objective = tries[[i]]$value, 
        kkt = tries[[i]]$kkt, code = tries[[i]]$code)))
    best
}

paper_penalty <- 
function (fit) 
{
    stopifnot(length(fit$sp) == length(fit$smooth), is.null(fit$full.sp), is.null(attr(fit, "e5_ridge")), 
        all(is.finite(fit$sp)), all(fit$sp >= 0))
    S <- matrix(0, length(coef(fit)), length(coef(fit)))
    for (sm in fit$smooth) {
        stopifnot(length(sm$S) == 1L, sm$first.sp == sm$last.sp)
        ix <- sm$first.para:sm$last.para
        S[ix, ix] <- S[ix, ix] + fit$sp[sm$first.sp] * sm$S[[1]]
    }
    S
}

.gcv_fit <- 
function (formula, d, w) 
{
    args <- list(formula = formula, data = d, family = quasibinomial(), weights = w, method = "GCV.Cp", 
        gamma = 1.4, select = FALSE, control = mgcv::gam.control(maxit = 500))
    valid <- function(f) !is.null(f) && isTRUE(f$converged) && all(is.finite(coef(f))) && is.finite(f$gcv.ubre) && 
        identical(f$method, "GCV") && (is.null(f$outer.info$conv) || identical(f$outer.info$conv, "full convergence"))
    z <- tryCatch(do.call(mgcv::gam, args), error = function(err) {
        warning(paste("GCV initial fit:", conditionMessage(err)))
        NULL
    })
    if (valid(z)) 
        return(z)
    warning("GCV candidate failed original fit; trying same-GCV starts, no REML or fixed-scale substitution")
    for (outer in c("newton", "bfgs")) for (sp0 in c(1, 10, 0.1, 100)) {
        ns <- length(mgcv::interpret.gam(formula)$smooth.spec)
        zz <- tryCatch(do.call(mgcv::gam, c(args, list(optimizer = c("outer", outer), mustart = rep(0.5, 
            nrow(d)), in.out = list(sp = rep(sp0, ns), scale = 1)))), error = function(err) {
            warning(paste("GCV restart:", conditionMessage(err)))
            NULL
        })
        if (valid(zz)) {
            warning(paste("Accepted same-GCV restart", outer, sp0))
            return(zz)
        }
    }
    stop("GCV candidate failed all same-criterion restarts; do not silently select the surviving structure")
}

.gcv_choose <- 
function (d, w, k, role, h, source = FALSE, group = FALSE) 
{
    fits <- lapply(c(FALSE, TRUE), function(ii) .gcv_fit(.gcv_formula(ii, k, source, group = group), 
        d, w))
    scores <- vapply(fits, function(f) unname(f$gcv.ubre), numeric(1))
    chosen <- if (scores[2] < scores[1] - 1e-12 * max(1, abs(scores))) 
        2L
    else 1L
    tab <- data.frame(role = role, horizon = h, n = nrow(d), training_ids_sha256 = digest::digest(d$eid, 
        algo = "sha256"), candidate_interaction = c(FALSE, TRUE), gcv = scores, selected = 1:2 == chosen, 
        edf = vapply(fits, function(f) sum(f$edf), numeric(1)), weight_sum = sum(w), k = k)
    .gcv_trace[[length(.gcv_trace) + 1L]] <<- tab
    list(fit = fits[[chosen]], interaction = chosen == 2L, comparison = tab)
}

.gcv_formula <- 
function (interaction, k = 5, source = FALSE, group = FALSE) 
{
    rhs <- c("X1", "X2", sprintf("s(X%d_sc,k=%d,bs='ts')", 3:.ext_p, k))
    if (interaction) 
        rhs <- c(rhs, if (source) sprintf("s(X3_sc,by=X1,k=%d,bs='ts')", k) else sprintf("s(X3_sc,by=ordered(X1,levels=c(0,1)),k=%d,bs='ts')", 
            k))
    if (group) 
        rhs <- c(rhs, "target_group")
    as.formula(paste("D_t ~", paste(rhs, collapse = " + ")))
}

.cens_initial_retry <- 
function (fit) 
{
    if (isTRUE(fit$converged)) 
        return(fit)
    warning("Initial GAM did not converge; retry identical REML objective from prespecified smoothing starts")
    d <- fit$model
    if ("offset(source_lp)" %in% names(d) && !"source_lp" %in% names(d)) 
        d$source_lp <- d[["offset(source_lp)"]]
    audit <- list(original_converged = fit$converged, original_outer = fit$outer.info$conv, original_scale = fit$scale, 
        original_sp = fit$sp, attempts = list())
    for (outer in c("newton", "bfgs")) for (initial_sp in c(10, 1, 100, 0.1, 1000, 10000)) {
        candidate <- tryCatch(do.call(mgcv::gam, list(formula = formula(fit), data = d, family = quasibinomial(), 
            weights = fit$prior.weights, method = "REML", gamma = 1.4, select = FALSE, mustart = rep(0.5, 
                nrow(d)), in.out = list(sp = rep(initial_sp, length(fit$sp)), scale = 1), optimizer = c("outer", 
                outer), control = mgcv::gam.control(maxit = 500))), error = function(err) {
            warning(paste("Same-objective GAM restart failed:", conditionMessage(err)))
            NULL
        })
        audit$attempts[[length(audit$attempts) + 1L]] <- list(optimizer = outer, initial_sp = initial_sp, 
            converged = if (is.null(candidate)) FALSE else candidate$converged, outer = if (is.null(candidate)) "error" else candidate$outer.info$conv, 
            reml = if (is.null(candidate)) NA_real_ else unname(candidate$gcv.ubre))
        if (!is.null(candidate) && isTRUE(candidate$converged) && identical(candidate$outer.info$conv, 
            "full convergence")) {
            stopifnot(identical(names(coef(candidate)), names(coef(fit))), all(is.finite(coef(candidate))))
            attr(candidate, "same_objective_restart") <- audit
            warning(paste("Accepted same-objective GAM restart:", outer, "initial sp", initial_sp))
            return(candidate)
        }
    }
    stop("All same-objective GAM restarts failed full convergence")
}

.cur_retry_gam <- 
function (fit) 
{
    if (isTRUE(fit$converged)) 
        return(fit)
    z <- tryCatch(.cens_initial_retry(fit), error = function(err) {
        warning(paste("Original GAM retries exhausted:", conditionMessage(err)))
        NULL
    })
    if (!is.null(z)) 
        return(z)
    warning("EXPLICIT FALLBACK: fixed quasi-binomial working scale=1 after original fit/restarts fail; estimator changed")
    d <- fit$model
    if ("offset(source_lp)" %in% names(d) && !"source_lp" %in% names(d)) 
        d$source_lp <- d[["offset(source_lp)"]]
    args <- list(formula = formula(fit), data = d, family = quasibinomial(), weights = fit$prior.weights, 
        method = "REML", gamma = 1.4, select = FALSE, scale = 1, control = mgcv::gam.control(maxit = 500))
    starts <- list(list(optimizer = c("outer", "newton")))
    for (outer in c("bfgs", "newton")) for (sp in c(10, 1, 100, 0.1, 1000, 10000)) starts[[length(starts) + 
        1L]] <- list(optimizer = c("outer", outer), mustart = rep(0.5, nrow(d)), in.out = list(sp = rep(sp, 
        length(fit$sp)), scale = 1))
    trace <- list()
    for (i in seq_along(starts)) {
        z <- tryCatch(do.call(mgcv::gam, c(args, starts[[i]])), error = function(err) {
            warning(paste("Fixed-scale same-objective restart:", conditionMessage(err)))
            NULL
        })
        ok <- !is.null(z) && isTRUE(z$converged) && identical(z$outer.info$conv, "full convergence") && 
            all(is.finite(coef(z)))
        trace[[i]] <- list(attempt = i, accepted = ok, outer = if (is.null(z)) "error" else z$outer.info$conv)
        if (ok) {
            stopifnot(identical(names(coef(z)), names(coef(fit))), z$scale == 1, !z$scale.estimated)
            attr(z, "explicit_dispersion_fallback") <- list(original_scale = fit$scale, new_scale = 1, 
                statistical_objective_changed = TRUE, attempts = trace)
            .cens_fallbacks[[length(.cens_fallbacks) + 1L]] <<- data.frame(formula = paste(deparse(formula(fit)), 
                collapse = ""), n = nrow(d), original_scale = fit$scale, new_scale = 1, converged = z$converged, 
                outer = z$outer.info$conv, fixed_scale_attempt = i)
            return(z)
        }
    }
    stop("All approved fixed-scale same-objective GAM restarts failed full convergence")
}

.paper_folds <- 
function (data, t, k, seed) 
{
    strata <- .state_36m(data, t)
    if (any(table(strata) < k)) 
        stop("A 36-month observability stratum is smaller than k.")
    set.seed(seed)
    fold <- integer(nrow(data))
    for (s in levels(strata)) {
        ix <- sample(which(strata == s))
        fold[ix] <- rep(seq_len(k), length.out = length(ix))
    }
    if (any(tabulate(fold, k) == 0L)) 
        stop("Empty outer fold.")
    fold
}

.score_source <- 
function (d) 
{
    f <- reformulate(paste0("X", 1:.ext_p), response = "survival::Surv(Y,delta)")
    fit <- survival::coxph(f, data = d, ties = "breslow", x = TRUE, y = TRUE)
    beta <- as.numeric(coef(fit))
    stopifnot(all(is.finite(beta)), length(beta) == .ext_p)
    list(beta = beta, baseline = .score_baseline(d$Y, d$delta, drop(.score_x(d) %*% beta)))
}

.score_predict <- 
function (obj, newx, h = c(30, 36, 48)) 
{
    bh <- c(0, cumsum(obj$baseline[, 2]))[findInterval(h, obj$baseline[, 1]) + 1L]
    lp <- drop(as.matrix(newx) %*% obj$beta)
    ans <- -expm1(-outer(exp(lp), bh))
    stopifnot(all(is.finite(ans)), all(ans >= 0 & ans <= 1))
    ans
}

.score_coxkl <- 
function (train, source, seed, cfg) 
{
    x <- .score_x(train)
    center <- colMeans(x)
    scales <- apply(x, 2, sd)
    z <- scale(x, center = center, scale = scales)
    rs <- drop(x %*% source$beta)
    etas <- survkl::generate_eta(n = cfg$coxkl_eta_n, max_eta = cfg$coxkl_eta_max)
    cv <- survkl::cv.coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = etas, 
        nfolds = cfg$inner_folds, criteria = cfg$coxkl_criterion, seed = seed, tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, 
        backtrack = TRUE)
    trace <- cv$internal_stat
    stopifnot(all(is.finite(trace$VVH_Loss)))
    eta <- trace$eta[which.min(trace$VVH_Loss)]
    fit <- survkl::coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = eta, 
        tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, backtrack = TRUE)
    refine <- survkl::coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = eta, 
        tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, backtrack = TRUE, beta_initial = as.numeric(fit$beta))
    change <- max(abs(fit$beta - refine$beta))
    stopifnot(is.finite(change), change < 1e-05)
    beta <- as.numeric(fit$beta)/scales
    list(beta = beta, baseline = .score_baseline(train$Y, train$delta, drop(x %*% beta)), tuning = eta, 
        trace = data.frame(tuning = trace$eta, loss = trace$VVH_Loss), convergence_change = change, center = center, 
        scale = scales)
}

.usim_coxkl_extended <- 
function (train, source, seed, cfg) 
{
    x <- .score_x(train)
    center <- colMeans(x)
    scales <- apply(x, 2, sd)
    z <- scale(x, center = center, scale = scales)
    rs <- drop(x %*% source$beta)
    etas <- sort(unique(c(survkl::generate_eta(n = cfg$coxkl_eta_n, max_eta = cfg$coxkl_eta_max), 10, 
        20, 50)))
    stopifnot(length(etas) == 13, 0 %in% etas, max(etas) == 50)
    cv <- survkl::cv.coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = etas, 
        nfolds = cfg$inner_folds, criteria = cfg$coxkl_criterion, seed = seed, tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, 
        backtrack = TRUE)
    trace <- cv$internal_stat
    stopifnot(all(is.finite(trace$VVH_Loss)))
    eta <- trace$eta[which.min(trace$VVH_Loss)]
    fit <- survkl::coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = eta, 
        tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, backtrack = TRUE)
    refine <- survkl::coxkl(z, train$delta, train$Y, stratum = rep(1, nrow(train)), RS = rs, etas = eta, 
        tol = cfg$coxkl_tol, Mstop = cfg$coxkl_maxit, backtrack = TRUE, beta_initial = as.numeric(fit$beta))
    change <- max(abs(fit$beta - refine$beta))
    stopifnot(is.finite(change), change < 1e-05)
    beta <- as.numeric(fit$beta)/scales
    list(beta = beta, baseline = .score_baseline(train$Y, train$delta, drop(x %*% beta)), tuning = eta, 
        trace = data.frame(tuning = trace$eta, loss = trace$VVH_Loss), convergence_change = change)
}

.score_potl <- 
function (train, source, seed, cfg, progress = function(...) NULL) 
{
    set.seed(seed)
    folds <- integer(nrow(train))
    for (v in 0:1) {
        ix <- which(train$delta == v)
        folds[ix] <- sample(rep(1:cfg$inner_folds, length.out = length(ix)))
    }
    xis <- c(0, 2^seq(cfg$potl_xi_exponents[1], cfg$potl_xi_exponents[2]))
    trace <- list()
    for (i in seq_along(xis)) {
        for (k in 1:cfg$inner_folds) {
            error <- ""
            loss <- tryCatch({
                fit <- .score_potl_fit(train[folds != k, ], source, xis[i], cfg, train[folds == k, ])
                -fit$logPL
            }, error = function(e) {
                error <<- conditionMessage(e)
                Inf
            })
            trace[[length(trace) + 1L]] <- data.frame(tuning = xis[i], inner_fold = k, loss = loss, error = error)
        }
        progress(i, xis[i])
    }
    trace <- do.call(rbind, trace)
    loss <- vapply(xis, function(a) mean(trace$loss[trace$tuning == a]), numeric(1))
    stopifnot(is.finite(loss[1]), sum(is.finite(loss)) >= 2)
    xi <- xis[which.min(loss)]
    fit <- .score_potl_fit(train, source, xi, cfg)
    fit$tuning <- xi
    fit$trace <- trace
    fit$inner_folds <- folds
    fit$convergence_change <- 0
    fit
}

.score_potl_fit <- 
function (train, source, xi, cfg, validation = NULL) 
{
    train <- train[order(train$Y, train$row_id), ]
    x <- .score_x(train)
    center <- colMeans(x)
    scales <- apply(x, 2, sd)
    x <- scale(x, center, scales)
    profile <- !is.null(validation)
    if (!profile) 
        validation <- train[sort(unique(c(1:3, which(train$delta == 1)[1], which(train$delta == 0)[1]))), 
            ]
    validation <- validation[order(validation$Y, validation$row_id), ]
    xv <- scale(.score_x(validation), center, scales)
    end <- max(c(train$Y, validation$Y, max(cfg$horizons))) + 1
    raw <- function(z) cbind(seq_len(nrow(z)) - 1, 0, end, z)
    bh <- c(0, cumsum(source$baseline[, 2]))[findInterval(train$Y, source$baseline[, 1]) + 1L]
    s <- exp(-bh * exp(drop(.score_x(train) %*% source$beta)))
    stopifnot(all(is.finite(s)), all(s >= 0 & s <= 1))
    messages <- capture.output(fit <- .potl_api$TransFit(train$Y, as.integer(train$delta), raw(x), s, 
        rep(1/nrow(train), nrow(train)), validation$Y, as.integer(validation$delta), raw(xv), cbind(c(0, 
            end), c(1, 1)), r = cfg$potl_r, N = 20L, xi = xi, endt = max(cfg$horizons), tol = cfg$potl_tol, 
        maxit = cfg$potl_maxit, verbose = 0L, pll = as.integer(profile)))
    if (any(grepl("Nonconvergence|failed|singular", messages, ignore.case = TRUE))) 
        stop(paste(messages, collapse = "; "))
    stopifnot(all(is.finite(fit$beta)), all(is.finite(fit$lambda)), all(fit$lambda[, 2] >= 0), is.finite(fit$logL))
    if (profile) 
        stopifnot(is.finite(fit$logPL))
    beta <- as.numeric(fit$beta)/scales
    baseline <- fit$lambda
    baseline[, 2] <- baseline[, 2] * exp(-sum(center * beta))
    list(beta = beta, baseline = baseline, logPL = if (profile) fit$logPL else NA_real_, logL = fit$logL, 
        messages = messages)
}

.r100_extra_fit <- 
function (method, s, t, te, hs, seed, e) 
{
    cov <- paste0("X", seq_len(.ext_p))
    .suite_h <<- hs
    clinical <- isTRUE(getOption("spt.paper.clinical", FALSE))
    tuning_hs <- if (clinical && method == "Random-Survival-Forest") c(24,36,48) else hs
    if (method == "TransCox-BIC") {
        set.seed(seed)
        o <- fit_transcox_model(s, t, source_covariates = cov, target_covariates = cov, select_bic = TRUE, 
            backend = "r")
        stopifnot(all(is.finite(c(o$fit$new_beta, o$fit$new_IntH))))
        e$write(as.data.frame(as.table(o$bic$BICmat)), "tuning")
        e$write(data.frame(lambda1 = o$lambda1, lambda2 = o$lambda2, negative_hazard_increments = sum(o$fit$new_IntH < 
            0)), "diagnostics")
        p <- sapply(hs, function(h) predict_transcox_risk_t(o, te, h)$p_event_t)
    }
    else if (method == "CoxTL") {
        o <- .suite_fit_coxtl(s, t, seed)
        z <- attr(o$fit, "native_tuning")
        stopifnot(all(is.finite(z$scores)))
        stopifnot(z$selected == z$size_ratio * z$grid[which.max(rowMeans(z$scores))])
        e$write(data.frame(tuning = rep(z$grid, 5), fold = rep(1:5, each = length(z$grid)), score = as.vector(z$scores)), 
            "tuning")
        e$write(data.frame(row_id = t$row_id, fold = z$folds), "partitions")
        e$write(data.frame(selected_lambda = z$selected, density_code = o$density$code, density_ess = unname(o$density_diagnostic["ess"])), 
            "diagnostics")
        p <- .suite_predict_coxtl(o, te)
    }
    else if (method %in% c("DiscreteKL", "Trans-AFT")) {
        cfg <- jsonlite::fromJSON("config/noncox.json")
        if (method == "DiscreteKL") {
            .nc_load("discretekl")
            cfg$breaks <- c(0, unlist(lapply(1:3, function(j) seq(c(0, hs)[j], hs[j], length.out = 5)[-1])))
            cfg$discrete_administrative_horizon <- hs[3]
            cfg$discrete_interval_months <- 1
            sc <- .nc_disc_source(s, cfg)
            o <- .nc_discretekl(t, sc, seed, cfg)
            e$write(data.frame(boundary = cfg$breaks), "intervals")
        }
        else {
            .nc_load("transaft")
            o <- .nc_transaft(t, s, seed, cfg)
        }
        e$write(o$trace, "tuning")
        e$write(data.frame(index = seq_along(o$inner_foldid), fold = o$inner_foldid), "partitions")
        e$write(data.frame(selected = o$tuning, convergence_change = o$convergence_change, failed_inner_fits = sum(nzchar(o$trace$error))), 
            "diagnostics")
        p <- .nc_predict(o, te[cov], hs)
    }
    else {
        stopifnot(method %in% c("Random-Survival-Forest", "IPCW-XGBoost"))
        inner <- .e7_helpers$.make_stratified_folds(t, if(clinical)36 else hs[2], 3, seed)
        e$write(data.frame(row_id = t$row_id, fold = inner), "partitions")
        grid <- if (method == "Random-Survival-Forest") 
            expand.grid(mtry = unique(c(floor(sqrt(.ext_p)), ceiling(.ext_p/2), .ext_p)), node = c(10, 
                20, 40))
        else expand.grid(depth = c(1, 2, 3), rounds = c(100, 300))
        fitrf <- function(tr, z, ss) ranger::ranger(as.formula(paste("survival::Surv(Y,delta)~", paste(cov, 
            collapse = "+"))), data = tr, num.trees = 500, mtry = z$mtry, min.node.size = z$node, splitrule = "logrank", 
            num.threads = 1, seed = ss)
        predrf <- function(o, d, at=hs) {
            pr <- predict(o, data = d, num.threads = 1)
            sapply(at, function(h) {
                ix <- findInterval(h, pr$unique.death.times)
                if (ix == 0) 
                  rep(0, nrow(d))
                else 1 - pr$survival[, ix]
            })
        }
        fitx <- function(tr, z, h, ss) fit_ipcw_xgb_model(tr, h, covariates = cov, nrounds = z$rounds, 
            max_depth = z$depth, nthread = 1, seed = ss)
        score <- array(NA_real_, c(nrow(grid), 3, if (method == "Random-Survival-Forest") 1 else 3))
        for (k in 1:3) {
            tr <- t[inner != k, ]
            va <- t[inner == k, ]
            km <- .guard_km_fit(tr)
            for (g in seq_len(nrow(grid))) {
                e$stage <- paste(method, "inner", k, "grid", g)
                if (method == "Random-Survival-Forest") {
                  o <- fitrf(tr, grid[g, ], seed + k)
                  pp <- predrf(o, va, tuning_hs)
                  score[g, k, 1] <- mean(sapply(1:3, function(j) mean(ucv_validation(va, tuning_hs[j], matrix(pp[, 
                    j], ncol = 1), km)$loss)))
                }
                else for (j in 1:3) {
                  o <- fitx(tr, grid[g, ], hs[j], seed + k)
                  pp <- predict_ipcw_xgb_risk_t(o, va)$p_event_t
                  score[g, k, j] <- mean(ucv_validation(va, hs[j], matrix(pp, ncol = 1), km)$loss)
                }
            }
        }
        stopifnot(all(is.finite(score)))
        trace <- expand.grid(candidate = seq_len(nrow(grid)), fold = 1:3, horizon_index = seq_len(dim(score)[3]))
        trace$loss <- as.vector(score)
        e$write(trace, "tuning")
        e$write(cbind(candidate = seq_len(nrow(grid)), grid), "grid")
        selected <- vapply(seq_len(dim(score)[3]), function(j) which.min(rowMeans(score[, , j])), integer(1))
        e$write(data.frame(horizon_index = seq_along(selected), candidate = selected), "diagnostics")
        e$stage <- paste(method, "full")
        if (method == "Random-Survival-Forest") 
            p <- predrf(fitrf(t, grid[selected, ], seed), te)
        else p <- sapply(1:3, function(j) predict_ipcw_xgb_risk_t(fitx(t, grid[selected[j], ], hs[j], 
            seed), te)$p_event_t)
    }
    stopifnot(identical(dim(p), c(nrow(te), 3L)), all(is.finite(p)), all(p >= 0 & p <= 1))
    if (method != "IPCW-XGBoost") 
        stopifnot(all(p[, 2] >= p[, 1] - 1e-12), all(p[, 3] >= p[, 2] - 1e-12))
    pmin(pmax(p, 0.001), 0.999)
}

.suite_fit_coxtl <- 
function (source, target, seed) 
{
    cov <- paste0("X", seq_len(.ext_p))
    center <- sapply(target[cov], mean)
    scale <- sapply(target[cov], stats::sd)
    stopifnot(all(is.finite(scale)), all(scale > 0))
    convert <- function(d, r) {
        x <- as.data.frame(scale(as.matrix(d[cov]), center = center, scale = scale))
        x$time <- d$Y
        x$status <- d$delta
        x$R <- factor(rep(r, nrow(d)), levels = c(0, 1))
        x
    }
    s <- convert(source, 1)
    t <- convert(target, 0)
    xs <- cbind(1, as.matrix(s[cov]))
    xt <- cbind(1, as.matrix(t[cov]))
    opt <- nlm(.suite_coxtl_env$density_opt(xs, xt, 1), rep(0, .ext_p + 1), iterlim = 500)
    stopifnot(opt$code %in% 1:2, all(is.finite(opt$estimate)))
    weights <- as.numeric(exp(xs %*% opt$estimate))
    stopifnot(all(is.finite(weights)), all(weights > 0))
    set.seed(seed)
    fit <- .suite_coxtl_env$run_CoxTL(t, s, weights = weights, di = .ext_p, folds_num = 5, lam_set = c(0.001, 
        0.005, 0.01, 0.05, 0.1, 1:10))
    stopifnot(!isTRUE(fit$fail), all(is.finite(stats::coef(fit))))
    list(fit = fit, center = center, scale = scale, density = opt, density_diagnostic = c(min = min(weights), 
        max = max(weights), mean = mean(weights), ess = sum(weights)^2/sum(weights^2)))
}

.suite_predict_coxtl <- 
function (obj, test) 
{
    cov <- paste0("X", seq_len(.ext_p))
    df <- as.data.frame(scale(as.matrix(test[cov]), center = obj$center, scale = obj$scale))
    df$R <- factor(rep(0, nrow(df)), levels = c(0, 1))
    1 - as.matrix(pec::predictSurvProb(obj$fit, newdata = df, times = .suite_h))
}

.nc_disc_source <- 
function (source, cfg) 
{
    s <- .nc_scale(.nc_x(source))
    d <- .nc_disc_data(source, s, cfg)
    list(fit = .nc_disc_fit(d, NULL, 0, cfg), scaling = s, source_n = nrow(source), source_events = sum(d$event), 
        interval_event_counts = tabulate(d$time[d$event == 1], nbins = length(cfg$breaks)-1L))
}

.nc_discretekl <- 
function (train, source, seed, cfg, progress = function(...) NULL) 
{
    set.seed(seed)
    folds <- integer(nrow(train))
    strata <- ifelse(train$delta == 1 & train$Y <= cfg$discrete_administrative_horizon, 1, 0)
    for (v in 0:1) {
        ix <- which(strata == v)
        folds[ix] <- sample(rep(1:cfg$discrete_inner_folds, length.out = length(ix)))
    }
    d <- .nc_disc_data(train, source$scaling, cfg)
    subset_d <- function(ix) list(time = d$time[ix], event = d$event[ix], x = d$x[ix, , drop = FALSE])
    etas <- seq(cfg$discrete_eta_min, cfg$discrete_eta_max, by = cfg$discrete_eta_step)
    traces <- list()
    for (i in seq_along(etas)) {
        for (k in 1:cfg$discrete_inner_folds) {
            error <- ""
            score <- NA_real_
            loss <- tryCatch({
                fit <- .nc_disc_fit(subset_d(folds != k), source$fit, etas[i], cfg)
                score <- fit$score_max
                val <- subset_d(folds == k)
                loss <- -DiscreteKL::DiscLoglik_logit(val$time, val$x, val$event, fit$beta_t, fit$beta_v)
                stopifnot(is.finite(loss))
                as.numeric(loss)
            }, error = function(e) {
                error <<- conditionMessage(e)
                Inf
            })
            traces[[length(traces) + 1L]] <- data.frame(tuning = etas[i], inner_fold = k, loss = loss, 
                score_max = score, error = error)
        }
        if (i%%20 == 0) 
            progress(i, length(etas))
    }
    trace <- do.call(rbind, traces)
    losses <- vapply(etas, function(a) mean(trace$loss[trace$tuning == a]), numeric(1))
    stopifnot(is.finite(losses[1]), sum(is.finite(losses)) >= 2)
    eta <- etas[which.min(losses)]
    fit <- .nc_disc_fit(d, source$fit, eta, cfg)
    list(kind = "discretekl", fit = fit, scaling = source$scaling, interval = 1, breaks = cfg$breaks, 
        tuning = eta, trace = trace, inner_foldid = folds, convergence_change = fit$score_max, training_events = sum(d$event), 
        zero_interval_rows = sum(d$time == 0))
}

.nc_transaft <- 
function (train, source, seed, cfg) 
{
    stopifnot(all(train$Y > 0), all(source$Y > 0))
    scaling <- .nc_scale(.nc_x(train))
    a <- .nc_km_block(train, scaling)
    b <- .nc_km_block(source, scaling)
    n0 <- nrow(train)
    nk <- nrow(source)
    raw <- list(X = rbind(a$x, b$x), y = c(a$y, b$y), status = c(a$status, b$status))
    centered <- .aft_api$KM.Weight.Center(n0, nk, 1, raw)
    stopifnot(all(is.finite(centered$X)), all(is.finite(centered$y)))
    set.seed(seed)
    fit <- .aft_api$Trans.AFT(centered$X, centered$y, n0, nk, 1)
    stopifnot(length(fit$beta) == .ext_p, all(is.finite(c(fit$beta, fit$theta.S))))
    cv <- .aft_api$last_cv
    stopifnot(length(cv$lambda) == 10, all(is.finite(cv$cvm)))
    constant <- cv$lambda.min/sqrt(2 * log(.ext_p)/(n0 + nk))
    intercept <- a$ybar - sum(a$xbar * fit$beta)
    residual <- log(train$Y) - intercept - drop(.nc_z(.nc_x(train), scaling) %*% fit$beta)
    km <- survival::survfit(survival::Surv(residual, train$delta) ~ 1)
    stopifnot(all(is.finite(km$surv)), all(diff(km$surv) <= 1e-12))
    trace <- data.frame(tuning = cv$lambda/sqrt(2 * log(.ext_p)/(n0 + nk)), inner_fold = 0, loss = cv$cvm, 
        loss_se = cv$cvsd, lambda = cv$lambda, error = "")
    list(kind = "transaft", beta = fit$beta, theta_source = fit$theta.S, intercept = intercept, scaling = scaling, 
        residual_km = list(time = km$time, surv = km$surv), tuning = constant, trace = trace, inner_foldid = .aft_api$last_foldid, 
        center_order_target = a$order, center_order_source = b$order, convergence_change = 0, training_events = sum(train$delta), 
        zero_interval_rows = 0, lambda_pool = cv$lambda.min, lambda_correction = constant * sqrt(2 * 
            log(.ext_p)/n0), km_weight_sum_target = sum(a$weight), km_weight_sum_source = sum(b$weight))
}

.nc_predict <- 
function (object, newx, h = c(30, 36, 48)) 
{
    z <- .nc_z(as.matrix(newx), object$scaling)
    if (object$kind == "discretekl") {
        indices <- match(h, object$breaks) - 1L
        stopifnot(!anyNA(indices), all(indices > 0), max(indices) <= length(object$fit$beta_t))
        logit <- outer(drop(z %*% object$fit$beta_v), drop(object$fit$beta_t), "+")
        logsurv <- -pmax(logit, 0) - log1p(exp(-abs(logit)))
        ans <- vapply(indices, function(t) -expm1(rowSums(logsurv[, seq_len(t), drop = FALSE])), numeric(nrow(z)))
    }
    else {
        lp <- object$intercept + drop(z %*% object$beta)
        km <- object$residual_km
        ans <- vapply(h, function(t) 1 - c(1, km$surv)[findInterval(log(t) - lp, km$time) + 1L], numeric(nrow(z)))
    }
    stopifnot(all(is.finite(ans)), all(ans >= 0 & ans <= 1))
    ans
}

.usim_g <- 
function (x, setting, domain, b = 0.8) 
{
    stopifnot(setting %in% 1:4, domain %in% c("source", "target"))
    a <- 0.6 * (x$X1 - 0.4)/sqrt(0.24) - 0.45 * (x$X2 - 0.5)/0.5 + 0.75 * x$X3 - 0.5 * x$X5 + 0.35 * 
        x$X7 + 0.55 * (x$X4^2 - 1)/sqrt(0.8) + 0.45 * sin(pi * x$X6/sqrt(3))
    u <- (x$X1 - 0.4)/sqrt(0.24) * x$X3
    if (domain == "source" || setting == 1) 
        return(a + b * u)
    if (setting == 2) 
        return(a - 0.55 * x$X3 + 0.85 * x$X5 + b * u)
    if (setting == 3) 
        return(a + 0.6 * b * u)
    a
}

.uro_boxgrad <- 
function (g, x) 
{
    if (x <= 1e-08) 
        min(g, 0)
    else if (x >= 1 - 1e-08) 
        max(g, 0)
    else g
}

.uro_softplus <- 
function (x) 
pmax(x, 0) + log1p(exp(-abs(x)))

.score_baseline <- 
function (time, event, lp) 
{
    times <- sort(unique(time[event == 1]))
    stopifnot(length(times) > 0, all(is.finite(lp)))
    increments <- vapply(times, function(t) sum(event[time == t])/sum(exp(lp[time >= t])), numeric(1))
    stopifnot(all(is.finite(increments)), all(increments >= 0))
    cbind(time = times, increment = increments)
}

.score_x <- 
function (d) 
as.matrix(d[paste0("X", 1:.ext_p)])

fit_ipcw_xgb_model <- 
function (train_data, t, covariates = paste0("X", 1:7), eps_G = 0.001, w_cap = 20, nrounds = 150, max_depth = 2, 
    eta = 0.05, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5, nthread = as.integer(getOption("translogistic.xgb_threads", 
        1)), seed = 202605) 
{
    .require_pkg("xgboost")
    aux <- make_Dt_weights(data = train_data, t = t, eps_G = eps_G, w_cap = w_cap)
    x_train <- as.matrix(train_data[, covariates, drop = FALSE])
    storage.mode(x_train) <- "double"
    dtrain <- xgboost::xgb.DMatrix(data = x_train, label = aux$D_t, weight = aux$weights)
    params <- list(objective = "binary:logistic", eval_metric = "logloss", max_depth = max_depth, eta = eta, 
        subsample = subsample, colsample_bytree = colsample_bytree, min_child_weight = min_child_weight, 
        nthread = max(1L, nthread))
    set.seed(seed)
    fit <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
    list(fit = fit, covariates = covariates, model_type = "ipcw_xgb", t = t, nrounds = nrounds, params = params, 
        n_event = sum(aux$D_t == 1, na.rm = TRUE), n_nonevent = sum(aux$D_t == 0, na.rm = TRUE))
}

fit_transcox_model <- 
function (source_data, target_train, source_covariates = paste0("X", 1:7), target_covariates = paste0("X", 
    1:7), select_bic = TRUE, lambda1 = getOption("translogistic.transcox_lambda1", 0.1), lambda2 = getOption("translogistic.transcox_lambda2", 
    0.1), lambda1_vec = getOption("translogistic.transcox_lambda1_vec", c(0.1, 0.5, 1, 2)), lambda2_vec = getOption("translogistic.transcox_lambda2_vec", 
    c(0.1, 0.5, 1, 2)), learning_rate = getOption("translogistic.transcox_learning_rate", 0.004), nsteps = getOption("translogistic.transcox_nsteps", 
    100), backend = getOption("translogistic.transcox_backend", "r"), standardize = TRUE) 
{
    load_transcox_sources()
    cov_tc <- paste0("X", seq_along(target_covariates))
    source_tc <- make_transcox_data(source_data, source_covariates, cov_tc)
    target_tc <- make_transcox_data(target_train, target_covariates, cov_tc)
    if (isTRUE(standardize)) {
        sc <- standardize_transcox_pair(source_tc, target_tc, cov_tc)
        source_tc <- sc$source
        target_tc <- sc$target
    }
    else {
        sc <- list(center = setNames(rep(0, length(cov_tc)), cov_tc), scale = setNames(rep(1, length(cov_tc)), 
            cov_tc), covariates = cov_tc)
    }
    if (sum(source_tc$status == 2) < 1 || sum(target_tc$status == 2) < 1) {
        stop("TransCox requires at least one event in both source and target training data.")
    }
    Cout <- GetAuxSurv(source_tc, cov = cov_tc)
    Pout <- GetPrimaryParam(target_tc, q = Cout$q, estR = Cout$estR)
    bic <- NULL
    if (isTRUE(select_bic)) {
        bic <- SelParam_By_BIC(primData = target_tc, auxData = source_tc, cov = cov_tc, lambda1_vec = lambda1_vec, 
            lambda2_vec = lambda2_vec, learning_rate = learning_rate, nsteps = nsteps, backend = backend, 
            verbose = FALSE)
        lambda1 <- bic$best_la1
        lambda2 <- bic$best_la2
    }
    fit <- runTransCox_one(Pout = Pout, l1 = lambda1, l2 = lambda2, learning_rate = learning_rate, nsteps = nsteps, 
        backend = backend, cov = cov_tc)
    obj <- list(method = if (isTRUE(select_bic)) "transcox_bic" else "transcox_fixed", fit = fit, Pout = Pout, 
        Cout = Cout, bic = bic, lambda1 = lambda1, lambda2 = lambda2, learning_rate = learning_rate, 
        nsteps = nsteps, backend = backend, scaler = sc, source_covariates = source_covariates, target_covariates = target_covariates, 
        covariates = cov_tc)
    class(obj) <- "transcox_model"
    obj
}

predict_ipcw_xgb_risk_t <- 
function (object, newdata) 
{
    x_new <- as.matrix(newdata[, object$covariates, drop = FALSE])
    storage.mode(x_new) <- "double"
    p_event_t <- as.numeric(stats::predict(object$fit, newdata = x_new))
    p_event_t <- clip_probs(p_event_t, eps = 1e-05)
    list(lp = stats::qlogis(p_event_t), p_event_t = p_event_t)
}

predict_transcox_risk_t <- 
function (object, newdata, t) 
{
    x <- as.data.frame(newdata[, object$target_covariates, drop = FALSE])
    names(x) <- object$covariates
    for (v in object$covariates) {
        x[[v]] <- (x[[v]] - object$scaler$center[[v]])/object$scaler$scale[[v]]
    }
    lp <- as.vector(as.matrix(x[, object$covariates, drop = FALSE]) %*% as.numeric(object$fit$new_beta))
    event_times <- as.numeric(object$fit$time)
    increments <- pmax(as.numeric(object$fit$new_IntH), 0)
    H0_t <- sum(increments[event_times <= t], na.rm = TRUE)
    if (!is.finite(H0_t) || H0_t < 0) 
        H0_t <- 0
    p_event_t <- 1 - exp(-H0_t * exp(lp))
    p_event_t <- clip_probs(p_event_t, eps = 1e-05)
    list(lp = lp, p_event_t = p_event_t)
}

.nc_disc_data <- 
function (d, s, cfg) 
{
    stopifnot(all(d$Y > 0), all(d$delta %in% 0:1))
    end <- cfg$discrete_administrative_horizon
    step <- cfg$discrete_interval_months
    event <- as.integer(d$delta == 1 & d$Y <= end)
    time <- ifelse(event == 1, findInterval(d$Y, cfg$breaks, left.open = TRUE), findInterval(pmin(d$Y, 
        end), cfg$breaks) - 1L)
    list(time = as.integer(time), event = event, x = .nc_z(.nc_x(d), s))
}

.nc_disc_fit <- 
function (d, prior, eta, cfg) 
{
    stopifnot(max(d$time) == length(cfg$breaks)-1L)
    fit <- if (is.null(prior)) 
        DiscreteKL::discSurv_logit(d$time, d$x, d$event, tol = cfg$discrete_tol, max_iter = cfg$discrete_max_iter)
    else DiscreteKL::discSurvKL_logit(d$time, d$x, d$event, tol = cfg$discrete_tol, max_iter = cfg$discrete_max_iter, 
        beta_t_tilde = prior$beta_t, beta_v_tilde = prior$beta_v, eta = eta)
    stopifnot(all(is.finite(c(fit$beta_t, fit$beta_v, fit$loglik))))
    o <- order(d$time, decreasing = TRUE)
    t <- matrix(d$time[o])
    x <- d$x[o, , drop = FALSE]
    y <- matrix(d$event[o])
    u <- if (is.null(prior)) 
        DiscreteKL::Update_logit(t, x, y, fit$beta_t, fit$beta_v, max(t), ncol(x), nrow(x))
    else DiscreteKL::UpdateKL_logit(t, x, y, fit$beta_t, fit$beta_v, max(t), ncol(x), nrow(x), prior$beta_t, 
        prior$beta_v, eta)
    fit$score_max <- max(abs(c(u$score_t, u$score_v)))
    stopifnot(is.finite(fit$score_max), fit$score_max < cfg$discrete_max_score)
    fit
}

.nc_scale <- 
function (x) 
{
    ans <- list(center = colMeans(x), scale = apply(x, 2, sd))
    stopifnot(all(is.finite(ans$scale)), all(ans$scale > 0))
    ans
}

.nc_x <- 
function (d) 
as.matrix(d[paste0("X", seq_len(.ext_p))])

.nc_km_block <- 
function (d, scaling) 
{
    o <- order(d$Y, -d$delta, seq_len(nrow(d)))
    d <- d[o, ]
    x <- .nc_z(.nc_x(d), scaling)
    y <- log(d$Y)
    n <- nrow(d)
    w <- numeric(n)
    s <- 1
    for (i in seq_len(n)) {
        w[i] <- s * d$delta[i]/(n - i + 1)
        s <- s * (1 - d$delta[i]/(n - i + 1))
    }
    stopifnot(sum(w) > 0, all(is.finite(w)))
    list(x = x, y = y, status = d$delta, weight = w, xbar = colSums(w * x)/sum(w), ybar = sum(w * y)/sum(w), 
        order = o)
}

.nc_z <- 
function (x, s) 
scale(x, center = s$center, scale = s$scale)

.require_pkg <- 
function (pkg) 
{
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is required. Please install it first.", pkg))
    }
}

make_transcox_data <- 
function (data, covariates, out_covariates = paste0("X", seq_along(covariates))) 
{
    out <- data.frame(time = as.numeric(data$Y), status = ifelse(as.integer(data$delta) == 1L, 2L, 1L), 
        stringsAsFactors = FALSE)
    x <- as.data.frame(data[, covariates, drop = FALSE])
    names(x) <- out_covariates
    cbind(out, x)
}

standardize_transcox_pair <- 
function (source_tc, target_tc, covariates) 
{
    pooled_x <- rbind(source_tc[, covariates, drop = FALSE], target_tc[, covariates, drop = FALSE])
    center <- vapply(pooled_x, mean, numeric(1), na.rm = TRUE)
    scale <- vapply(pooled_x, stats::sd, numeric(1), na.rm = TRUE)
    scale[!is.finite(scale) | scale <= 0] <- 1
    apply_one <- function(df) {
        for (v in covariates) df[[v]] <- (df[[v]] - center[[v]])/scale[[v]]
        df
    }
    list(source = apply_one(source_tc), target = apply_one(target_tc), center = center, scale = scale, 
        covariates = covariates)
}

clip_probs <- 
function (p, eps = 0.001) 
pmin(pmax(p, eps), 1 - eps)

.state_36m <- 
function (data, t = 36) 
{
    factor(ifelse(data$Y <= t & data$delta == 1L, "event_by_t", ifelse(data$Y > t, "observed_event_free_t", 
        "censored_by_t")), levels = c("event_by_t", "observed_event_free_t", "censored_by_t"))
}

