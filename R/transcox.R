# Local TransCox R backend; upstream attribution in THIRD_PARTY.md.
GetAuxSurv <- 
function (auxData, weights = NULL, cov = c("X1", "X2")) 
{
    auxData <- validate_transcox_data(auxData, cov = cov)
    fml <- stats::reformulate(cov, response = "survival::Surv(time, status == 2)")
    res.cox <- survival::coxph(fml, data = auxData, weights = weights, ties = "breslow")
    bhest <- survival::basehaz(res.cox, centered = FALSE)
    estR <- res.cox$coefficients
    q <- data.frame(cumHazards = bhest$hazard, breakPoints = bhest$time)
    return(list(estR = estR, q = q))
}

GetPrimaryParam <- 
function (primData, q, estR) 
{
    primData <- validate_transcox_data(primData, cov = names(estR))
    primData <- primData[order(primData$time), ]
    dQ <- deltaQ(primData, q)
    primData <- primData[order(primData$time), ]
    fullcum <- rep(0, nrow(primData))
    idx0 <- match(dQ$t, primData$time)
    fullcum[idx0] <- dQ$cumQ_upd
    fullcum <- cummax(fullcum)
    primData$fullCumQ = fullcum
    event_idx <- which(primData$status == 2)
    Xinn <- outer(seq_len(nrow(primData)), event_idx, FUN = ">=") * 1
    return(list(primData = primData, Xinn = Xinn, dQ = dQ, estR = estR))
}

SelParam_By_BIC <- 
function (primData, auxData, cov = c("X1", "X2"), statusvar = "status", lambda1_vec = c(0.1, 0.5, seq(1, 
    10, by = 0.5)), lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)), learning_rate = 0.004, nsteps = 100, 
    verbose = TRUE, backend = getOption("translogistic.transcox_backend", "r")) 
{
    Cout <- GetAuxSurv(auxData, cov = cov)
    Pout <- GetPrimaryParam(primData, q = Cout$q, estR = Cout$estR)
    CovData = Pout$primData[, cov, drop = FALSE]
    status = Pout$primData[, statusvar]
    cumH = Pout$primData$fullCumQ
    hazards = Pout$dQ$dQ
    BICmat <- matrix(NA, length(lambda1_vec), length(lambda2_vec))
    total <- length(lambda1_vec) * length(lambda2_vec)
    pb <- if (isTRUE(verbose)) 
        txtProgressBar(min = 0, max = total, style = 2, initial = 0)
    else NULL
    for (i in seq_along(lambda1_vec)) {
        for (j in seq_along(lambda2_vec)) {
            thisi = (i - 1) * length(lambda2_vec) + j
            if (!is.null(pb)) 
                setTxtProgressBar(pb, thisi)
            lambda1 = lambda1_vec[i]
            lambda2 = lambda2_vec[j]
            BICvalue <- tryCatch({
                test <- run_transcox_backend(CovData = as.matrix(CovData), cumH = cumH, hazards = hazards, 
                  status = status, estR = Pout$estR, Xinn = Pout$Xinn, lambda1 = lambda1, lambda2 = lambda2, 
                  learning_rate = learning_rate, nsteps = nsteps, backend = backend)
                names(test) <- c("eta", "xi")
                newBeta = Pout$estR + test$eta
                newHaz = Pout$dQ$dQ + test$xi
                GetBIC(status = status, CovData = CovData, hazards = hazards, newBeta = newBeta, newHaz = newHaz, 
                  eta = test$eta, xi = test$xi, cutoff = 1e-05)
            }, error = function(e) Inf)
            BICmat[i, j] <- BICvalue
        }
    }
    if (!is.null(pb)) 
        close(pb)
    rownames(BICmat) = lambda1_vec
    colnames(BICmat) <- lambda2_vec
    if (all(!is.finite(BICmat))) 
        stop("All TransCox BIC fits failed.")
    idx0 <- which(BICmat == min(BICmat, na.rm = TRUE), arr.ind = TRUE)
    idx0 <- idx0[1, , drop = FALSE]
    b_lambda1 <- lambda1_vec[idx0[1]]
    b_lambda2 <- lambda2_vec[idx0[2]]
    return(list(best_la1 = b_lambda1, best_la2 = b_lambda2, BICmat = BICmat))
}

runTransCox_one <- 
function (Pout, l1 = 1, l2 = 1, learning_rate = 0.004, nsteps = 200, cov = c("X1", "X2"), backend = getOption("translogistic.transcox_backend", 
    "r")) 
{
    CovData <- Pout$primData[, cov, drop = FALSE]
    status <- Pout$primData[, "status"]
    cumH <- Pout$primData$fullCumQ
    hazards <- Pout$dQ$dQ
    test <- run_transcox_backend(CovData = as.matrix(CovData), cumH = cumH, hazards = hazards, status = status, 
        estR = Pout$estR, Xinn = Pout$Xinn, lambda1 = l1, lambda2 = l2, learning_rate = learning_rate, 
        nsteps = nsteps, backend = backend)
    names(test) <- c("eta", "xi")
    new_haz <- as.numeric(Pout$dQ$dQ + test$xi)
    return(list(eta = test$eta, xi = test$xi, new_beta = Pout$estR + test$eta, new_IntH = new_haz, time = Pout$primData[status == 
        2, "time"]))
}

TransCox_R <- 
function (CovData, cumH, hazards, status, estR, Xinn, lambda1, lambda2, learning_rate = 0.004, nsteps = 100, 
    l1_smooth = 1e-06, min_hazard = 1e-08) 
{
    CovData <- as.matrix(CovData)
    Xinn <- as.matrix(Xinn)
    cumH <- as.numeric(cumH)
    hazards <- as.numeric(hazards)
    status <- as.numeric(status)
    estR <- as.numeric(estR)
    p <- length(estR)
    m <- length(hazards)
    if (ncol(CovData) != p) 
        stop("CovData columns do not match length(estR).")
    if (ncol(Xinn) != m) 
        stop("Xinn columns do not match length(hazards).")
    if (length(cumH) != nrow(CovData) || length(status) != nrow(CovData)) {
        stop("cumH/status length must match rows of CovData.")
    }
    hazards0 <- pmax(hazards, min_hazard)
    start <- c(rep(0, p), log(hazards0))
    objective <- function(par) {
        eta <- par[seq_len(p)]
        log_haz <- par[p + seq_len(m)]
        new_haz <- exp(pmin(pmax(log_haz, log(min_hazard)), 30))
        xi <- new_haz - hazards
        beta <- estR + eta
        xb <- as.vector(CovData %*% beta)
        if (any(!is.finite(xb))) 
            return(Inf)
        cum_adj <- cumH + as.vector(Xinn %*% xi)
        if (any(!is.finite(cum_adj))) 
            return(Inf)
        loglik <- sum(xb[status == 2] + log(new_haz)) - sum(cum_adj * exp(pmin(xb, 40)))
        penalty <- lambda1 * sum(sqrt(eta^2 + l1_smooth)) + lambda2 * sum(sqrt(xi^2 + l1_smooth))
        val <- -loglik + penalty
        if (!is.finite(val)) 
            Inf
        else val
    }
    fit <- stats::optim(par = start, fn = objective, method = "BFGS", control = list(maxit = max(1L, 
        as.integer(nsteps))))
    par <- fit$par
    eta <- par[seq_len(p)]
    log_haz <- par[p + seq_len(m)]
    xi <- exp(pmin(pmax(log_haz, log(min_hazard)), 30)) - hazards
    attr(eta, "optim_convergence") <- fit$convergence
    attr(eta, "optim_value") <- fit$value
    list(eta, xi)
}

validate_transcox_data <- 
function (data, cov = c("X1", "X2"), timevar = "time", statusvar = "status") 
{
    needed <- unique(c(timevar, statusvar, cov))
    missing <- setdiff(needed, names(data))
    if (length(missing) > 0) {
        stop("Missing column(s): ", paste(missing, collapse = ", "))
    }
    keep <- stats::complete.cases(data[, needed, drop = FALSE]) & is.finite(data[[timevar]]) & data[[timevar]] > 
        0 & data[[statusvar]] %in% c(1, 2)
    data[keep, , drop = FALSE]
}

GetBIC <- 
function (status, CovData, hazards, newBeta, newHaz, eta, xi, cutoff = 1e-05) 
{
    Logl <- GetLogLike(status = status, CovData = CovData, hazards = hazards, newBeta = newBeta, newHaz = newHaz)
    K = sum(c(abs(eta), abs(xi)) > cutoff)
    n = length(status)
    BIC <- K * log(n) - 2 * Logl
    return(BIC)
}

GetLogLike <- 
function (status, CovData, hazards, newBeta, newHaz) 
{
    XB = as.matrix(CovData) %*% newBeta
    newCumH <- dQtocumQ(newHaz, status)
    haz_event <- pmax(as.numeric(newHaz), .Machine$double.eps)
    LogL <- sum(XB[status == 2] + log(haz_event)) - sum(newCumH * exp(XB))
    return(LogL)
}

deltaQ <- 
function (primData, q) 
{
    obsData <- primData[primData$status == 2, ]
    obsData <- obsData[order(obsData$time), ]
    if (nrow(obsData) == 0) 
        stop("TransCox requires at least one event in the primary data.")
    newQ <- data.frame(dQ = rep(NA, nrow(obsData)), cumQ = rep(NA, nrow(obsData)), t = rep(NA, nrow(obsData)))
    for (i in 1:nrow(obsData)) {
        if (i == 1) {
            newQ$t[i] <- obsData$time[i]
            tmp1 <- q$breakPoints <= obsData$time[i]
            if (all(!tmp1)) {
                idx0 <- 1
            }
            else {
                idx0 <- max(which(tmp1))
            }
            newQ$dQ[i] = newQ$cumQ[i] <- q$cumHazards[idx0]
        }
        else {
            newQ$t[i] <- obsData$time[i]
            tmp1 <- q$breakPoints <= obsData$time[i]
            if (all(!tmp1)) {
                idx0 <- 1
            }
            else {
                idx0 <- max(which(tmp1))
            }
            newQ$cumQ[i] <- q$cumHazards[idx0]
            newQ$dQ[i] <- q$cumHazards[idx0] - newQ$cumQ[(i - 1)]
        }
    }
    newQ$dQ[newQ$dQ == 0] <- 1e-04
    newQ$cumQ_upd <- newQ$cumQ
    for (i in 1:nrow(newQ)) {
        newQ$cumQ_upd[i] <- sum(newQ$dQ[1:i])
    }
    return(newQ)
}

dQtocumQ <- 
function (dQ, status = NULL) 
{
    if (length(status) > 0) {
        cumQ = rep(NA, length(status))
        ncount = 1
        for (i in 1:length(status)) {
            if (status[i] == 2) {
                cumQ[i] = sum(dQ[1:ncount])
                ncount = ncount + 1
                if (ncount > length(dQ)) {
                  ncount = length(dQ)
                }
            }
            else {
                cumQ[i] = sum(dQ[1:ncount])
            }
        }
    }
    else {
        cumQ = rep(NA, length(dQ))
        for (i in 1:length(dQ)) {
            cumQ[i] = sum(dQ[1:i])
        }
    }
    return(cumQ)
}

run_transcox_backend <- function(CovData,cumH,hazards,status,estR,Xinn,lambda1,lambda2,learning_rate=.004,nsteps=100,backend='r') {
 stopifnot(backend=='r')
 TransCox_R(CovData,cumH,hazards,status,estR,Xinn,lambda1,lambda2,learning_rate,nsteps)
}
