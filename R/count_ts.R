#' Apply the algorithms to make decisions for Thompson sampling Negative Binomial (TS-NB) algorithms
#'
#' @param context context at the current decision time
#' @param beta_normalNB the randomly sampled Bayesian estimate
#' @import MASS parallel fastDummies stats
#'
#' @return Intervention option
#' @examples
#' apply_normalNB(matrix(1:10, nrow = 2),matrix(11:20, nrow = 5))
#'
#' @export
#'
apply_normalNB <- function(context, beta_normalNB){
  rewards <- exp(context%*%beta_normalNB)
  arm <- which.max(rewards)
  return(arm)
}

#' Apply the algorithms to make decisions for Thompson sampling Poisson (TS-Poisson) algorithms
#'
#' @param context context at the current decision time
#' @param beta_laplacePoisson the randomly sampled Bayesian estimate
#' @import MASS parallel fastDummies stats
#'
#' @return Intervention option
#' @examples
#' apply_laplacePoisson(matrix(1:10, nrow = 2),matrix(11:20, nrow = 5))
#' @export
#'
apply_laplacePoisson <- function(context, beta_laplacePoisson){
  rewards <- exp(context%*%beta_laplacePoisson)
  arm <- which.max(rewards)
  return(arm)
}

#' Apply the algorithms to make decisions for Linear Thompson sampling (TS) algorithms
#'
#' @param context context at the current decision time
#' @param beta_linearTS the randomly sampled Bayesian estimate
#' @import MASS parallel fastDummies stats
#'
#' @return Intervention option
#' @examples
#' apply_linearTS(matrix(1:10, nrow = 2),matrix(11:20, nrow = 5))
#' @export
apply_linearTS <- function(context, beta_linearTS){
  rewards <- exp(context%*%beta_linearTS)
  arm <- which.max(rewards)
  return(arm)
}

#' Apply the algorithms to make decisions for Thompson sampling Zero-inflated Poisson (TS-ZIP) algorithm
#'
#' @param context context at the current decision time
#' @param beta_ZIP the randomly sampled Bayesian estimate for the Poisson component
#' @param gamma_ZIP the randomly sampled Bayesian estimate for the zero component
#' @import MASS parallel fastDummies stats
#' @examples
#' apply_ZIP(matrix(1:10, nrow = 2),matrix(11:20, nrow = 5),matrix(21:30, nrow = 5))
#' @return Intervention option
#'
#' @export
apply_ZIP <- function(context, beta_ZIP, gamma_ZIP){
  zero <- exp(context%*%gamma_ZIP)/(1+exp(context%*%gamma_ZIP))
  if (is.na(zero[1])){
    zero = 0
  }
  rewards <- exp(context%*%beta_ZIP)*(1-zero)
  arm <- which.max(rewards)
  return(arm)
}

#' Apply the algorithms to make decisions for Thompson sampling Zero-inflated Negative Binomial (TS-ZINB) algorithm
#'
#' @param context context at the current decision time
#' @param beta_ZINB the randomly sampled Bayesian estimate for the Poisson component
#' @param gamma_ZINB the randomly sampled Bayesian estimate for the zero component
#' @import MASS parallel fastDummies stats
#' @examples
#' apply_ZINB(matrix(1:10, nrow = 2),matrix(11:20, nrow = 5),matrix(21:30, nrow = 5))
#' @return Intervention option
#'
#' @export
apply_ZINB <- function(context, beta_ZINB, gamma_ZINB){
  zero <- exp(context%*%gamma_ZINB)/(1+exp(context%*%gamma_ZINB))
  if (is.na(zero[1])){
    zero = 0
  }
  rewards <- exp(context%*%beta_ZINB)*(1-zero)
  arm <- which.max(rewards)
  return(arm)
}

# function of backtrack line search
backtrack <- function(x, dx, f, Df, al=0.2, be=0.6) {
  t <- 1
  g <- Df(x)
  u <- al*(t(g)%*%dx)
  k <- 1
  repeat {
    if (f(x + t*dx) <= f(x) + t*u) break
    t <- be*t
    # print(paste0("backtrack ",k))
    k <- k + 1
  }
  return(t)
}

update_normalNB <- function(Y_normalNB, X_normalNB, alpha){
  fit <- glm(Y_normalNB~-1+X_normalNB,family = poisson("log"))
  mu_normalNB <- coef(fit)
  mu = exp(X_normalNB%*%mu_normalNB)
  len = length(Y_normalNB)
  # r_normalNB = len/sum((Y_normalNB/mu - 1)^2)

  iter = 0
  made.changes = TRUE
  max.iter = 1000
  eps = 0.0001
  r_normalNB = theta.ml(Y_normalNB, mu, limit=5, eps)

  while(made.changes && (iter < max.iter)) {
    loglik_beta <- function(mu_normalNB){
      mu = exp(X_normalNB%*%mu_normalNB)
      -sum(Y_normalNB * log(mu + (Y_normalNB == 0)) - (r_normalNB + Y_normalNB) * log(r_normalNB + mu))
    }
    loglik_r <- function(th, mu, y){
      -sum(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th)
           - (th + y) * log(th + mu))
    }
    Df <- function(mu_normalNB){
      mu = exp(X_normalNB%*%mu_normalNB)
      W <- diag(as.vector(mu/(r_normalNB+mu)), nrow = len)
      -r_normalNB*(t(X_normalNB)%*%W%*%((Y_normalNB - mu)/mu))
    }

    iter = iter + 1
    made.changes = FALSE
    mu = exp(X_normalNB%*%mu_normalNB)
    W <- diag(as.vector(mu/(r_normalNB+mu)), nrow = len)
    I.fish <- r_normalNB*(t(X_normalNB)%*%W%*%X_normalNB)
    gradient <- Df(mu_normalNB)
    d = -ginv(I.fish)%*%gradient
    a = backtrack(mu_normalNB, d, loglik_beta, Df, al=0.2, be=0.8)
    beta0 = mu_normalNB
    mu_normalNB <- mu_normalNB + a*d
    fit_r = optimize(f = loglik_r, interval = c(1e-03,100),mu=exp(X_normalNB%*%mu_normalNB), y=Y_normalNB)
    r = r_normalNB
    r_normalNB = fit_r$minimum

    relative.change.beta = max(abs(mu_normalNB - beta0))
    relative.change.r = abs(r_normalNB - r)
    made.changes = ((relative.change.beta > eps) | (relative.change.r > eps))
  }

  mu = exp(X_normalNB%*%mu_normalNB)
  W <- diag(as.vector(mu/(r_normalNB+mu)), nrow = len)
  I.fish <- r_normalNB*t(X_normalNB)%*%W%*%X_normalNB
  variance <- ginv(I.fish)
  beta_normalNB <- mvrnorm(1, mu_normalNB, alpha * variance)
  return(list("beta_normalNB"=beta_normalNB,"mu_normalNB"=mu_normalNB, "variance"=variance))
}

update_laplacePoisson <- function(Y_laplacePoisson, X_laplacePoisson,alpha){
  iter = 0
  max.iter = 1000
  eps = 0.0001
  made.changes = TRUE
  len <- length(Y_laplacePoisson)
  d <- ncol(X_laplacePoisson)
  mu_laplacePoisson <- rnorm(d, 0, 1)
  #log likelihood
  f = function(b){-sum(Y_laplacePoisson*(X_laplacePoisson%*%b) - exp(X_laplacePoisson%*%b))} #negative log(posterior) ~ log(normal prior*likelihood)
  Df = function(b){-t(X_laplacePoisson)%*%(Y_laplacePoisson - exp(X_laplacePoisson%*%b))} #first derivative

  while(made.changes && (iter < max.iter)) {
    iter = iter + 1
    made.changes = FALSE
    H <- (t(X_laplacePoisson)%*%diag(c(exp(X_laplacePoisson%*%mu_laplacePoisson)), nrow = len)%*%X_laplacePoisson)#heissan matrix

    d <- -ginv(H)%*%Df(mu_laplacePoisson) #direction
    a = backtrack(mu_laplacePoisson, d, f, Df, al=0.2, be=0.8)
    new.beta_t <- mu_laplacePoisson
    mu_laplacePoisson <- mu_laplacePoisson + a*d
    relative.change = max(abs(new.beta_t - mu_laplacePoisson))
    made.changes = (relative.change > eps)
  }

  if (made.changes) {
    warning("Newton's method terminated before convergence")
  }
  H <- (t(X_laplacePoisson)%*%diag(c(exp(X_laplacePoisson%*%mu_laplacePoisson)), nrow = len)%*%X_laplacePoisson)
  variance <- ginv(H)
  beta_laplacePoisson <- mvrnorm(1, mu_laplacePoisson, alpha * variance)
  return(list("beta_laplacePoisson"=beta_laplacePoisson,"mu_laplacePoisson"=mu_laplacePoisson,"variance"=variance))
}

update_ZIP <- function(Y_ZIP, X_ZIP, alpha){
  mubeta_ZIP <- coef(glm(Y_ZIP ~ X_ZIP-1, family = poisson(link="log")))
  len <- length(Y_ZIP)
  made.changes = TRUE
  max.iter = 1000
  iter = 0
  eps = 0.0001
  mugamma_ZIP <- mubeta_ZIP
  # mugamma_ZIP <- rnorm(1,0,1)

  while (made.changes && (iter < max.iter)) {
    iter = iter + 1
    made.changes = FALSE
    # E step

    z_k <- 1/(1+exp(-X_ZIP%*%mugamma_ZIP - exp(X_ZIP%*%mubeta_ZIP))) #vector len * 1
    # z_k <- 1/(1+exp(-mugamma_ZIP - exp(X_ZIP%*%mubeta_ZIP)))
    z_k[Y_ZIP > 0] <- 0  #if y_t!=0, z_t = 0

    # M step for gam_k
    fit2 <- glm(z_k ~ X_ZIP-1, family = quasibinomial(link="logit"))
    # fit2 <- glm(z_k ~ 1, family = quasibinomial(link="logit"))
    gam_k <- mugamma_ZIP
    mugamma_ZIP <- coef(fit2)
    variance2 <- vcov(fit2)


    f = function(b){-sum((1-z_k)*(Y_ZIP*(X_ZIP%*%b) - exp(X_ZIP%*%b)))} #negative log(posterior) ~ log(normal prior*likelihood)
    Df = function(b){-t(X_ZIP)%*%((1-z_k)*(Y_ZIP - exp(X_ZIP%*%b)))} #first derivative


    H <- (t(X_ZIP)%*%diag(c(exp(X_ZIP%*%mubeta_ZIP)), nrow = len)%*%X_ZIP)#heissan matrix

    d <- -ginv(H)%*%Df(mubeta_ZIP) #direction
    a = backtrack(mubeta_ZIP, d, f, Df, al=0.2, be=0.8)
    new.beta_k <- mubeta_ZIP
    mubeta_ZIP <- mubeta_ZIP + a*d

    relative.change.beta = max(abs(new.beta_k - mubeta_ZIP))
    relative.change.gam = max(abs(gam_k - mugamma_ZIP))
    made.changes = ((relative.change.beta > eps) | (relative.change.gam > eps))

  }

  H <- (t(X_ZIP)%*%diag(c((1-z_k)*exp(X_ZIP%*%mubeta_ZIP)), nrow = len)%*%X_ZIP)
  variance1 <- ginv(H)
  beta_ZIP <- mvrnorm(1, mubeta_ZIP, alpha*variance1)
  gamma_ZIP <- mvrnorm(1, mugamma_ZIP, alpha*variance2)

  return(list("beta_ZIP"=beta_ZIP, "gamma_ZIP"=gamma_ZIP,
              "mubeta_ZIP"=mubeta_ZIP, "mugamma_ZIP"=mugamma_ZIP))
}

update_ZINB <- function(Y_ZINB, X_ZINB, alpha) {
  fit <- glm(Y_ZINB ~ X_ZINB-1, family = poisson(link="log"))
  mubeta_ZINB <- coef(fit)
  mu = exp(X_ZINB%*%mubeta_ZINB)
  n = length(Y_ZINB)
  # r_ZINB= n/sum((Y_ZINB/mu - 1)^2)

  made.changes = TRUE
  max.iter = 1000
  iter = 0
  eps = 0.0001
  r_ZINB= theta.ml(Y_ZINB, mu, limit=5, eps)
  mugamma_ZINB <- mubeta_ZINB
  # mugamma_ZINB <- rnorm(1,0,1)

  while(made.changes && (iter < max.iter)) {
    iter = iter + 1
    made.changes = FALSE

    # E step
    z_k <- 1/(1 + exp(-X_ZINB%*%mugamma_ZINB)*(r_ZINB/(r_ZINB + exp(X_ZINB%*%mubeta_ZINB)))^r_ZINB)
    z_k[Y_ZINB > 0] <- 0

    # M step for gam_k

    fit2 <- glm(z_k ~ X_ZINB-1, family = quasibinomial(link="logit"))
    # fit2 <- glm(z_k ~ 1, family = quasibinomial(link="logit"))
    gam_k <- mugamma_ZINB
    mugamma_ZINB <- coef(fit2)
    variance2 <- vcov(fit2)


    loglik_beta <- function(mubeta_ZINB){
      -sum((1-z_k)*(Y_ZINB * log(exp(X_ZINB%*%mubeta_ZINB) + (Y_ZINB == 0)) - (r_ZINB + Y_ZINB) * log(r_ZINB + exp(X_ZINB%*%mubeta_ZINB))))
    }
    loglik_r <- function(th, mu, y){
      -sum((1-z_k)*(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th)
                    - (th + y) * log(th + mu)))
    }
    Df <- function(b){
      -r_ZINB *t(X_ZINB)%*%W%*%((1-z_k)*(Y_ZINB - exp(X_ZINB%*%b))/exp(X_ZINB%*%b))
    }
    mu = exp(X_ZINB%*%mubeta_ZINB)
    W <- diag(as.vector(mu/(r_ZINB+mu)), nrow = n)
    I.fish <-r_ZINB*t(X_ZINB)%*%W%*%X_ZINB
    gradient <- Df(mubeta_ZINB)
    beta_k = mubeta_ZINB
    d = -ginv(I.fish)%*%gradient
    a = backtrack(mubeta_ZINB, d, loglik_beta, Df, al=0.2, be=0.8)
    mubeta_ZINB <- mubeta_ZINB + a*d
    fit_r = optimize(f = loglik_r, interval = c(1e-03,100),mu=exp(X_ZINB%*%mubeta_ZINB), y=Y_ZINB)
    r = r_ZINB
    r_ZINB = fit_r$minimum

    relative.change.beta = max(abs(mubeta_ZINB - beta_k))
    relative.change.gam = max(abs(mugamma_ZINB - gam_k))
    relative.change.r = abs(r_ZINB - r)
    made.changes = ((relative.change.beta > eps) | (relative.change.gam > eps)
                    | (relative.change.r > eps))
  }

  mu = exp(X_ZINB%*%mubeta_ZINB)
  W <- diag(as.vector((1-z_k)*mu/(r_ZINB+mu)), nrow = length(Y_ZINB))
  I.fish <- r_ZINB*t(X_ZINB)%*%W%*%X_ZINB
  variance1 <- ginv(I.fish)

  beta_ZINB <- mvrnorm(1, mubeta_ZINB, alpha*variance1)
  gamma_ZINB <- mvrnorm(1, mugamma_ZINB, alpha*variance2)
  return(list("beta_ZINB"=beta_ZINB, "gamma_ZINB"=gamma_ZINB,
              "mubeta_ZINB"=mubeta_ZINB, "mugamma_ZINB"=mugamma_ZINB))
}

update_linearTS <- function(Y_linearTS, X_linearTS, b_t, B_t, alpha){
  B_t <- B_t + t(X_linearTS)%*%X_linearTS
  b_t <- b_t + t(X_linearTS)%*%Y_linearTS
  mu_linearTS <- ginv(B_t)%*%b_t
  variance <- ginv(B_t)
  beta_linearTS <- mvrnorm(1, mu_linearTS, alpha*variance)
  return(list("beta_linearTS"=beta_linearTS,"b_t"=b_t, "B_t"=B_t, "mu_linearTS" = mu_linearTS))
}

#' Updating parameters in algorithm
#'
#' @param dist tuning parameter that controls which algorithm should be updated, with the options "Negative Binomial", "Poisson", "Linear TS", "ZIP", "ZINB"
#' @param Y_dist History of the observed stochastic outcome at the current decision time
#' @param X_dist History of the observed context at the current decision time
#' @param bt Sum of contexts weighted by the outcome, only for \code{dist} = "Linear TS", default is NULL.
#' @param Bt Outer product of contexts, only for \code{dist} = "Linear TS", default is NULL
#' @param alpha_dist tuning parameter that controls the exploration-exploitation tradeoff. Default is 1.
#' @import MASS parallel fastDummies stats
#'
#' @examples
#' update_algorithm(dist = "Negative Binomial")
#'
#'
#' @return The updated parameter estimates.
#'
#' @export


update_algorithm <- function(dist = c("Negative Binomial", "Poisson", "Linear TS",
                                      "ZIP", "ZINB"), Y_dist = 2, X_dist = 3, alpha_dist = 4, Bt = NULL, bt = NULL){
  if(dist == "Negative Binomial"){
    out <- update_normalNB(Y_normalNB = Y_dist, X_normalNB = X_dist, alpha = alpha_dist)
  }else if(dist == "Poisson"){
    out <- update_laplacePoisson(Y_laplacePoisson = Y_dist, X_laplacePoisson = X_dist, alpha = alpha_dist)
  }else if(dist == "ZIP"){
    out <- update_ZIP(Y_ZIP = Y_dist, X_ZIP = X_dist, alpha = alpha_dist)
  }else if(dist == "ZINB"){
    out <- update_ZINB(Y_ZINB = Y_dist, X_ZINB = X_dist, alpha = alpha_dist)
  }else{
    out <- update_linearTS(Y_linearTS = Y_dist, X_linearTS = X_dist, b_t = bt, B_t = Bt, alpha = alpha_dist)
  }
  return(out)
}

evaluate_policy <- function(beta, gamma, T0 = 1000, num_cov = 4, K = 20, T.init = 20, alpha = 1, gam = 25, distribution = c("Negative Binomial", "Poisson", "Linear TS",
                                                                                                                                                  "ZIP", "ZINB")){
  # initialize the estimator of linear ts
  b_linearTSlog <- rep(0, num_cov)
  B_linearTSlog <- diag(1, nrow = num_cov)

  tryCatch(
    expr = {

      # initialize the cumulative regret for each algorithm
      cumu_regret_normalNB <- 0
      cumu_regret_laplacePoisson <- 0
      cumu_regret_linearTSlog <- 0
      cumu_regret_ZIP <- 0
      cumu_regret_ZINB <- 0

      # save the history of cumulative regret for each algorithm
      cumu_regret_normalNB_hist <- c()
      cumu_regret_laplacePoisson_hist <- c()
      cumu_regret_linearTSlog_hist <- c()
      cumu_regret_ZIP_hist <- c()
      cumu_regret_ZINB_hist <- c()

      # save the history of regret for each algorithm
      regret_normalNB_hist <- c()
      regret_laplacePoisson_hist <- c()
      regret_linearTSlog_hist <- c()
      regret_ZIP_hist <- c()
      regret_ZINB_hist <- c()

      # save the history of beta and gamma for each algorithm
      beta_normalNB_hist <- NULL
      beta_laplacePoisson_hist <- NULL
      beta_linearTSlog_hist <- NULL
      beta_ZIP_hist <- NULL
      beta_ZINB_hist <- NULL
      gamma_ZIP_hist <- NULL
      gamma_ZINB_hist <- NULL

      # rewards observed by each algorithm
      Y_normalNB <- NULL
      Y_laplacePoisson <- NULL
      Y_linearTSlog <- NULL
      Y_ZIP <- NULL
      Y_ZINB <- NULL

      # feature vector collected by each algorithm
      X_normalNB <- NULL
      X_laplacePoisson <- NULL
      X_linearTSlog <- NULL
      X_ZIP <- NULL
      X_ZINB <- NULL


      for (t in c(1:T0)) {
        # generate context / feature
        context_original <- mvrnorm(K, rep(0,num_cov), diag(1,nrow = num_cov))

        # normalize the context
        context <- t(apply(context_original, 1, function(row) row / sqrt(sum(row^2))))

        # initial exploration stage: random selection
        if (t <= T.init) {
          ran_action = sample(seq(1,K),1,replace = T)
          arm_normalNB <- ran_action
          arm_laplacePoisson <- ran_action
          arm_linearTSlog <- ran_action
          arm_ZIP <- ran_action
          arm_ZINB <- ran_action
        }

        # after exploration stage: apply the policy
        if (t > T.init) {
          arm_normalNB <- apply_normalNB(context, beta_normalNB)
          arm_laplacePoisson <- apply_laplacePoisson(context, beta_laplacePoisson)
          arm_linearTSlog <- apply_linearTS(context, beta_linearTSlog)
          arm_ZIP <- apply_ZIP(context, beta_ZIP, gamma_ZIP)
          arm_ZINB <- apply_ZINB(context, beta_ZINB, gamma_ZINB)
        }


        #########when assuming a ZIP or ZINB environment model#########
        if (distribution == "ZIP"){
          # expected rewards for each action
          zero_part <- exp(context%*%gamma)/(1+exp(context%*%gamma))
          poisson_part <- exp(context%*%beta)
          optimal_reward <- max((1-zero_part)*poisson_part)

          # expected rewards of the selected action using different algorithms (for calculating regrets)
          expected_reward_normalNB <- (1-zero_part[arm_normalNB])*poisson_part[arm_normalNB]
          expected_reward_laplacePoisson <- (1-zero_part[arm_laplacePoisson])*poisson_part[arm_laplacePoisson]
          expected_reward_linearTSlog <- (1-zero_part[arm_linearTSlog])*poisson_part[arm_linearTSlog]
          expected_reward_ZIP <- (1-zero_part[arm_ZIP])*poisson_part[arm_ZIP]
          expected_reward_ZINB <- (1-zero_part[arm_ZINB])*poisson_part[arm_ZINB]


          # generate stchastic rewards from the environment
          stochastic_zero_part <- rbinom(K, 1, 1-zero_part)
          z_normalNB <- stochastic_zero_part[arm_normalNB]
          z_laplacePoisson <- stochastic_zero_part[arm_laplacePoisson]
          z_linearTSlog <- stochastic_zero_part[arm_linearTSlog]
          z_ZIP <- stochastic_zero_part[arm_ZIP]
          z_ZINB <- stochastic_zero_part[arm_ZINB]


          # theta <- rgamma(1, shape = gam, rate = gam) # dispersion when assuming a ZINB model
          theta <- 1 # when assuming a ZIP model
          stochastic_reward_normalNB <- z_normalNB*(rpois(1,poisson_part[arm_normalNB]*theta))
          stochastic_reward_laplacePoisson <- z_laplacePoisson*(rpois(1,poisson_part[arm_laplacePoisson]*theta))
          stochastic_reward_linearTSlog <- log(z_linearTSlog*(rpois(1,poisson_part[arm_linearTSlog]*theta))+0.1)
          stochastic_reward_ZIP <- z_ZIP*(rpois(1,poisson_part[arm_ZIP]*theta))
          stochastic_reward_ZINB <- z_ZINB*(rpois(1,poisson_part[arm_ZINB]*theta))
        }

        else if (distribution == "ZINB"){
          # expected rewards for each action
          zero_part <- exp(context%*%gamma)/(1+exp(context%*%gamma))
          poisson_part <- exp(context%*%beta)
          optimal_reward <- max((1-zero_part)*poisson_part)

          # expected rewards of the selected action using different algorithms (for calculating regrets)
          expected_reward_normalNB <- (1-zero_part[arm_normalNB])*poisson_part[arm_normalNB]
          expected_reward_laplacePoisson <- (1-zero_part[arm_laplacePoisson])*poisson_part[arm_laplacePoisson]
          expected_reward_linearTSlog <- (1-zero_part[arm_linearTSlog])*poisson_part[arm_linearTSlog]
          expected_reward_ZIP <- (1-zero_part[arm_ZIP])*poisson_part[arm_ZIP]
          expected_reward_ZINB <- (1-zero_part[arm_ZINB])*poisson_part[arm_ZINB]


          # generate stchastic rewards from the environment
          stochastic_zero_part <- rbinom(K, 1, 1-zero_part)
          z_normalNB <- stochastic_zero_part[arm_normalNB]
          z_laplacePoisson <- stochastic_zero_part[arm_laplacePoisson]
          z_linearTSlog <- stochastic_zero_part[arm_linearTSlog]
          z_ZIP <- stochastic_zero_part[arm_ZIP]
          z_ZINB <- stochastic_zero_part[arm_ZINB]


          # theta <- rgamma(1, shape = gam, rate = gam) # dispersion when assuming a ZINB model
          theta <- 1 # when assuming a ZIP model
          stochastic_reward_normalNB <- z_normalNB*(rpois(1,poisson_part[arm_normalNB]*theta))
          stochastic_reward_laplacePoisson <- z_laplacePoisson*(rpois(1,poisson_part[arm_laplacePoisson]*theta))
          stochastic_reward_linearTSlog <- log(z_linearTSlog*(rpois(1,poisson_part[arm_linearTSlog]*theta))+0.1)
          stochastic_reward_ZIP <- z_ZIP*(rpois(1,poisson_part[arm_ZIP]*theta))
          stochastic_reward_ZINB <- z_ZINB*(rpois(1,poisson_part[arm_ZINB]*theta))
        }


        #########when assuming a Poisson or NB enbironment model#########
        else if (distribution == "Negative Binomial"){
          # expected rewards for each action
          expected_reward_all <- exp(context%*%beta)
          optimal_reward <- max(expected_reward_all)

          # expected rewards of the selected action using different algorithms (for calculating regrets)
          expected_reward_normalNB <- expected_reward_all[arm_normalNB]
          expected_reward_laplacePoisson <- expected_reward_all[arm_laplacePoisson]
          expected_reward_linearTSlog <- expected_reward_all[arm_linearTSlog]
          expected_reward_ZIP <- expected_reward_all[arm_ZIP]
          expected_reward_ZINB <- expected_reward_all[arm_ZINB]

          # generate stochastic rewards from the environment
          theta <- rgamma(1, shape = gam, rate = gam) # dispersion when assuming a NB model
          # theta <- 1 # when assuming a Poisson model
          stochastic_reward_normalNB <- rpois(1,expected_reward_normalNB*theta)
          stochastic_reward_laplacePoisson <- rpois(1,expected_reward_laplacePoisson*theta)
          stochastic_reward_linearTSlog <- log(rpois(1,expected_reward_linearTSlog*theta)+0.1) #log transformation
          stochastic_reward_ZIP <- rpois(1,expected_reward_ZIP*theta)
          stochastic_reward_ZINB <- rpois(1,expected_reward_ZINB*theta)
        }

        else{
          # expected rewards for each action
          expected_reward_all <- exp(context%*%beta)
          optimal_reward <- max(expected_reward_all)

          # expected rewards of the selected action using different algorithms (for calculating regrets)
          expected_reward_normalNB <- expected_reward_all[arm_normalNB]
          expected_reward_laplacePoisson <- expected_reward_all[arm_laplacePoisson]
          expected_reward_linearTSlog <- expected_reward_all[arm_linearTSlog]
          expected_reward_ZIP <- expected_reward_all[arm_ZIP]
          expected_reward_ZINB <- expected_reward_all[arm_ZINB]

          # generate stchastic rewards from the environment
          # theta <- rgamma(1, shape = gam, rate = gam) # dispersion when assuming a NB model
          theta <- 1 # when assuming a Poisson model
          stochastic_reward_normalNB <- rpois(1,expected_reward_normalNB*theta)
          stochastic_reward_laplacePoisson <- rpois(1,expected_reward_laplacePoisson*theta)
          stochastic_reward_linearTSlog <- log(rpois(1,expected_reward_linearTSlog*theta)+0.1) #log transformation
          stochastic_reward_ZIP <- rpois(1,expected_reward_ZIP*theta)
          stochastic_reward_ZINB <- rpois(1,expected_reward_ZINB*theta)
        }

        # calculating regret
        regret_normalNB <- optimal_reward - expected_reward_normalNB
        regret_laplacePoisson <- optimal_reward - expected_reward_laplacePoisson
        regret_linearTSlog <- optimal_reward - expected_reward_linearTSlog
        regret_ZIP <- optimal_reward - expected_reward_ZIP
        regret_ZINB <- optimal_reward - expected_reward_ZINB

        # calculating cumulative regret
        cumu_regret_normalNB <- cumu_regret_normalNB + regret_normalNB
        cumu_regret_laplacePoisson <- cumu_regret_laplacePoisson + regret_laplacePoisson
        cumu_regret_linearTSlog <- cumu_regret_linearTSlog + regret_linearTSlog
        cumu_regret_ZIP <- cumu_regret_ZIP + regret_ZIP
        cumu_regret_ZINB <- cumu_regret_ZINB + regret_ZINB

        # save rewards for updating the estimates
        Y_normalNB <- c(Y_normalNB, stochastic_reward_normalNB)
        Y_laplacePoisson <- c(Y_laplacePoisson, stochastic_reward_laplacePoisson)
        Y_linearTSlog <- stochastic_reward_linearTSlog
        Y_ZIP <- c(Y_ZIP, stochastic_reward_ZIP)
        Y_ZINB <- c(Y_ZINB, stochastic_reward_ZINB)

        # save contexts for updating the estimates
        X_normalNB <- rbind(X_normalNB, context[arm_normalNB,])
        X_laplacePoisson <- rbind(X_laplacePoisson, context[arm_laplacePoisson,])
        X_linearTSlog <- context[arm_linearTSlog,,drop=F]
        X_ZIP <- rbind(X_ZIP, context[arm_ZIP,])
        X_ZINB <- rbind(X_ZINB, context[arm_ZINB,])


        ####updating parameters####
        # linear TS(log)
        para_linearTSlog <- update_linearTS(Y_linearTSlog, X_linearTSlog, b_linearTSlog, B_linearTSlog, alpha)
        b_linearTSlog <- para_linearTSlog$b_t
        B_linearTSlog <- para_linearTSlog$B_t
        beta_linearTSlog <- para_linearTSlog$beta_linearTS
        mu_linearTSlog <- para_linearTSlog$mu_linearTS

        if (t >= T.init) {
          # normalNB
          para_normalNB <- update_normalNB(Y_normalNB, X_normalNB, alpha)
          mu_normalNB <- para_normalNB$mu_normalNB
          beta_normalNB <- para_normalNB$beta_normalNB

          # laplace Poisson
          para_laplacePoisson <- update_laplacePoisson(Y_laplacePoisson, X_laplacePoisson, alpha)
          beta_laplacePoisson <- para_laplacePoisson$beta_laplacePoisson
          mu_laplacePoisson <- para_laplacePoisson$mu_laplacePoisson


          # ZIP
          para_ZIP <- update_ZIP(Y_ZIP, X_ZIP, alpha)
          mubeta_ZIP <- para_ZIP$mubeta_ZIP
          mugamma_ZIP <- para_ZIP$mugamma_ZIP
          beta_ZIP <- para_ZIP$beta_ZIP
          gamma_ZIP <- para_ZIP$gamma_ZIP


          # ZINB
          para_ZINB <- update_ZINB(Y_ZINB, X_ZINB, alpha)
          mubeta_ZINB <- para_ZINB$mubeta_ZINB
          mugamma_ZINB <- para_ZINB$mugamma_ZINB
          beta_ZINB <- para_ZINB$beta_ZINB
          gamma_ZINB <- para_ZINB$gamma_ZINB


          #save estimates history
          beta_normalNB_hist <- rbind(beta_normalNB_hist, t(mu_normalNB))
          beta_laplacePoisson_hist <- rbind(beta_laplacePoisson_hist, t(mu_laplacePoisson))
          beta_linearTSlog_hist <- rbind(beta_linearTSlog_hist, t(mu_linearTSlog))
          beta_ZIP_hist <- rbind(beta_ZIP_hist, t(mubeta_ZIP))
          beta_ZINB_hist <- rbind(beta_ZINB_hist, t(mubeta_ZINB))
          gamma_ZIP_hist <- rbind(gamma_ZIP_hist, t(mugamma_ZIP))
          gamma_ZINB_hist <- rbind(gamma_ZINB_hist, t(mugamma_ZINB))
        }

        #save cumulative regret history
        cumu_regret_normalNB_hist[t] <- cumu_regret_normalNB
        cumu_regret_laplacePoisson_hist[t] <- cumu_regret_laplacePoisson
        cumu_regret_linearTSlog_hist[t] <- cumu_regret_linearTSlog
        cumu_regret_ZIP_hist[t] <- cumu_regret_ZIP
        cumu_regret_ZINB_hist[t] <- cumu_regret_ZINB

        #save instant regret history
        regret_normalNB_hist[t] <- regret_normalNB
        regret_laplacePoisson_hist[t] <- regret_laplacePoisson
        regret_linearTSlog_hist[t] <- regret_linearTSlog
        regret_ZIP_hist[t] <- regret_ZIP
        regret_ZINB_hist[t] <- regret_ZINB
      }
      gen_normalNB <- list("Cum_Regret" = cumu_regret_normalNB_hist, "Regret" = regret_normalNB_hist,
                           "Beta" = beta_normalNB_hist)
      # if(isTRUE(save_file)){
      #   filename_normalNB <- paste0("normalNB_", s, ".RData")
      #   save(cumu_regret_normalNB_hist, regret_normalNB_hist, beta_normalNB_hist, file=filename_normalNB)
      # }

      gen_laplacePoisson <- list("Cum_Regret" = cumu_regret_laplacePoisson_hist, "Regret" = regret_laplacePoisson_hist,
                                 "Beta" = beta_laplacePoisson_hist)
      # if(isTRUE(save_file)){
      #   filename_laplacePoisson <- paste0("laplacePoisson_", s, ".RData")
      #   save(cumu_regret_laplacePoisson_hist, regret_laplacePoisson_hist, beta_laplacePoisson_hist, file=filename_laplacePoisson)
      # }

      gen_linearTSlog <- list("Cum_Regret" = cumu_regret_linearTSlog_hist, "Regret" = regret_linearTSlog_hist,
                              "Beta" = beta_linearTSlog_hist)
      # if(isTRUE(save_file)){
      #   filename_linearTSlog <- paste0("linearTSlog_", s, ".RData")
      #   save(cumu_regret_linearTSlog_hist, regret_linearTSlog_hist, beta_linearTSlog_hist, file=filename_linearTSlog)
      # }

      gen_ZIP <- list("Cum_Regret" = cumu_regret_ZIP_hist, "Regret" = regret_ZIP_hist,
                      "Beta" = beta_ZIP_hist, "Gamma" = gamma_ZIP_hist)
      # if(isTRUE(save_file)){
      #   filename_ZIP <- paste0("ZIP_", s, ".RData")
      #   save(cumu_regret_ZIP_hist, regret_ZIP_hist, beta_ZIP_hist, gamma_ZIP_hist,file = filename_ZIP)
      # }
      gen_ZINB <- list("Cum_Regret" = cumu_regret_ZINB_hist, "Regret" = regret_ZINB_hist,
                       "Beta" = beta_ZINB_hist, "Gamma" = gamma_ZINB_hist)
      # if(isTRUE(save_file)){
      #   filename_ZINB <- paste0("ZINB_", s, ".RData")
      #   save(cumu_regret_ZINB_hist, regret_ZINB_hist, beta_ZINB_hist, gamma_ZINB_hist, file = filename_ZINB)
      # }
      #gen_data = list(gen_normalNB, gen_laplacePoisson, gen_linearTSlog, gen_ZIP, gen_ZINB)
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  )
  return(list(gen_normalNB, gen_laplacePoisson, gen_linearTSlog, gen_ZIP, gen_ZINB))
}


######################################## Plot Function ###########################################
#' Summarize the simulation results and generate the regret plot
#'
#' @param S number of replicates of the experiment (greater than 1). Default is 30.
#' @param num_cov dimension for beta and gamma; we assume that they have the same dimensions for now. Default is 4.
#' @param T0 number of decision times. Default is 1000.
#' @param K number of actions/intervention options. Default is 20.
#' @param T.init length of the initial exploration stage. Default is 20.
#' @param alpha tuning parameter that controls the exploration-exploitation tradeoff. Default is 1.
#' @param gam over dispersion level of the environment model; this is only useful when the environment model is negative binomial or zero-inflated negative binomial. Default is 25.
#' @param dist_env tuning parameter that controls which environment model to use, with the options "Negative Binomial", "Poisson", "Linear TS", "ZIP", "ZINB"
#' @param show_figure A logical flag specifying that the regret plot of the model should be returned if true (default), otherwise, false.
#' @import MASS parallel fastDummies matrixStats ggplot2 stats
#'
#' @return The summary of the simulation results with cumulative regret, regret, and parameters is generated along with the optional
#' output of the regret plot (\code{show_figure} = TRUE).
#'
#' @export
#'
#' @references \itemize{
#' \item Liu, X., Deliu, N., Chakraborty, T., Bell, L., & Chakraborty, B. (2023).
#'       Thompson sampling for zero-inflated count outcomes with an application to the Drink Less mobile health study.
#'       arXiv preprint arXiv:2311.14359. \url{https://arxiv.org/abs/2311.14359}}
#'
#' @examples
#' output_summary(S = 2, num_cov = 2, T.init = 3, T0 = 5, dist_env = "Negative Binomial")

output_summary <- function(S = 30, num_cov = 4, T.init = 20, T0 = 1000, alpha = 1, gam = 25, K = 20,
                           dist_env = c("Negative Binomial", "Poisson", "Linear TS", "ZIP", "ZINB"), show_figure = TRUE){
  # generate the true model parameters and normalize them
  gamm <- rnorm(num_cov, 0, 1)
  gamma_val <- gamm / sqrt(sum(gamm^2))
  bet <- rnorm(num_cov, 0, 1)
  be <- bet/sqrt(sum(bet^2))

  # generate_data = mclapply(1:S, function(s) evaluate_policy(s, beta = be, gamma = gamma_val, T0 = T0, num_cov = num_cov, K = K,
  #                                                           T.init = T.init, alpha = alpha, gam = gam, save_file = save_file, distribution = dist_env),
  #                          mc.cores = mc.cores) #parallel::detectCores() -1)

  generate_data = lapply(1:S, function(s) evaluate_policy(beta = be, gamma = gamma_val, T0 = T0, num_cov = num_cov, K = K,
                                                          T.init = T.init, alpha = alpha, gam = gam, distribution = dist_env)) #parallel::detectCores() -1)
  # matrix for saving the results
  cumu_regret_normalNB_all <- matrix(nrow = S, ncol = T0)
  cumu_regret_laplacePoisson_all <- matrix(nrow = S, ncol = T0)
  cumu_regret_linearTSlog_all <- matrix(nrow = S, ncol = T0)
  cumu_regret_ZIP_all <- matrix(nrow = S, ncol = T0)
  cumu_regret_ZINB_all <- matrix(nrow = S, ncol = T0)

  regret_normalNB_all <- matrix(nrow = S, ncol = T0)
  regret_laplacePoisson_all <- matrix(nrow = S, ncol = T0)
  regret_linearTSlog_all <- matrix(nrow = S, ncol = T0)
  regret_ZIP_all <- matrix(nrow = S, ncol = T0)
  regret_ZINB_all <- matrix(nrow = S, ncol = T0)

  beta_normalNB_all <- list()
  beta_laplacePoisson_all <- list()
  beta_linearTSlog_all <- list()
  beta_ZIP_all <- list()
  beta_ZINB_all <- list()
  gamma_ZIP_all <- list()
  gamma_ZINB_all <- list()

  ###extract the results from the files
  for (s in 1:S) {
    tryCatch(
      expr = {
        # filename_normalNB <- paste0("normalNB_", s, ".RData")
        # filename_laplacePoisson <- paste0("laplacePoisson_", s, ".RData")
        # filename_linearTSlog <- paste0("linearTSlog_", s, ".RData")
        # filename_ZIP <- paste0("ZIP_", s, ".RData")
        # filename_ZINB <- paste0("ZINB_", s, ".RData")
        # load(filename_normalNB)
        # load(filename_laplacePoisson)
        # load(filename_linearTSlog)
        # load(filename_ZIP)
        # load(filename_ZINB)

        #save cumulative regret history
        cumu_regret_normalNB_all[s,] <- generate_data[[s]][[1]]$Cum_Regret
        cumu_regret_laplacePoisson_all[s,] <- generate_data[[s]][[2]]$Cum_Regret #cumu_regret_laplacePoisson_hist
        cumu_regret_linearTSlog_all[s,] <- generate_data[[s]][[3]]$Cum_Regret
        cumu_regret_ZIP_all[s,] <- generate_data[[s]][[4]]$Cum_Regret # cumu_regret_ZIP_hist
        cumu_regret_ZINB_all[s,] <- generate_data[[s]][[5]]$Cum_Regret #cumu_regret_ZINB_hist

        #save instant regret history
        regret_normalNB_all[s,] <- generate_data[[s]][[1]]$Regret #regret_normalNB_hist
        regret_laplacePoisson_all[s,] <- generate_data[[s]][[2]]$Regret #regret_laplacePoisson_hist
        regret_linearTSlog_all[s,] <- generate_data[[s]][[3]]$Regret #regret_linearTSlog_hist
        regret_ZIP_all[s,] <- generate_data[[s]][[4]]$Regret #regret_ZIP_hist
        regret_ZINB_all[s,] <- generate_data[[s]][[5]]$Regret #regret_ZINB_hist

        #save estimates history
        beta_normalNB_all[[s]] <- generate_data[[s]][[1]]$Beta #beta_normalNB_hist
        beta_laplacePoisson_all[[s]] <- generate_data[[s]][[2]]$Beta #beta_laplacePoisson_hist
        beta_linearTSlog_all[[s]] <- generate_data[[s]][[3]]$Beta #beta_linearTSlog_hist
        beta_ZIP_all[[s]] <- generate_data[[s]][[4]]$Beta #beta_ZIP_hist
        beta_ZINB_all[[s]] <- generate_data[[s]][[5]]$Beta #beta_ZINB_hist
        gamma_ZIP_all[[s]] <- generate_data[[s]][[4]]$Gamma #gamma_ZIP_hist
        gamma_ZINB_all[[s]] <- generate_data[[s]][[5]]$Gamma #gamma_ZINB_hist

      },
      error = function(e){
        message('Caught an error!')
        print(e)
      }
    )
  }

  # delete the cumulative regret for the initial stage
  cumu_regret_laplacePoisson_all <- cumu_regret_laplacePoisson_all[,-(1:T.init)]
  cumu_regret_normalNB_all <- cumu_regret_normalNB_all[,-(1:T.init)]
  cumu_regret_ZIP_all <- cumu_regret_ZIP_all[,-(1:T.init)]
  cumu_regret_ZINB_all <- cumu_regret_ZINB_all[,-(1:T.init)]
  cumu_regret_linearTSlog_all <- cumu_regret_linearTSlog_all[,-(1:T.init)]

  #mean cumulative regret
  avgcumu_regret_laplacePoisson <- colMeans(cumu_regret_laplacePoisson_all, na.rm = T)
  avgcumu_regret_normalNB <- colMeans(cumu_regret_normalNB_all, na.rm = T)
  avgcumu_regret_ZIP <- colMeans(cumu_regret_ZIP_all, na.rm = T)
  avgcumu_regret_ZINB <- colMeans(cumu_regret_ZINB_all, na.rm = T)
  avgcumu_regret_linearTSlog <- colMeans(cumu_regret_linearTSlog_all, na.rm = T)

  T1 = T0 - T.init

  # standard errors
  se_laplacePoisson <- colSds(cumu_regret_laplacePoisson_all, na.rm = T) / sqrt(S)
  se_normalNB <- colSds(cumu_regret_normalNB_all, na.rm = T) / sqrt(S)
  se_linearTSlog <- colSds(cumu_regret_linearTSlog_all, na.rm = T) / sqrt(S)
  se_ZIP <- colSds(cumu_regret_ZIP_all, na.rm = T) / sqrt(S)
  se_ZINB <- colSds(cumu_regret_ZINB_all, na.rm = T) / sqrt(S)
  ci_laplacePoisson <- qnorm(0.975) * se_laplacePoisson
  ci_normalNB <- qnorm(0.975) * se_normalNB
  ci_linearTSlog <- qnorm(0.975) * se_linearTSlog
  ci_ZIP <- qnorm(0.975) * se_ZIP
  ci_ZINB <- qnorm(0.975) * se_ZINB

  Trials <- rep(seq(T1), 5)
  Cumulative.regret <- c(avgcumu_regret_laplacePoisson, avgcumu_regret_normalNB,
                         avgcumu_regret_linearTSlog,
                         avgcumu_regret_ZIP, avgcumu_regret_ZINB)
  Algorithm <- c(rep("TS-Poisson", T1), rep("TS-NB", T1) ,
                 rep("Linear TS (log)",T1),
                 rep("TS-ZIP", T1), rep("TS-ZINB", T1))
  data <- data.frame(Trials, Cumulative.regret, Algorithm)

  # calculative the confidence band for each algorithm
  data$lower_laplacePoisson <- avgcumu_regret_laplacePoisson - ci_laplacePoisson
  data$upper_laplacePoisson <- avgcumu_regret_laplacePoisson + ci_laplacePoisson

  data$lower_normalNB <- avgcumu_regret_normalNB - ci_normalNB
  data$upper_normalNB <- avgcumu_regret_normalNB + ci_normalNB

  data$lower_ZIP <- avgcumu_regret_ZIP - ci_ZIP
  data$upper_ZIP <- avgcumu_regret_ZIP + ci_ZIP

  data$lower_linearTSlog <- avgcumu_regret_linearTSlog - ci_linearTSlog
  data$upper_linearTSlog <- avgcumu_regret_linearTSlog + ci_linearTSlog

  data$lower_ZINB <- avgcumu_regret_ZINB - ci_ZINB
  data$upper_ZINB <- avgcumu_regret_ZINB + ci_ZINB

  if(isTRUE(show_figure)){
    figure <- ggplot(data, aes(x=Trials, y=Cumulative.regret, colour=Algorithm)) +
      geom_line(lwd=0.8) +  xlab("t") + ylab("Regret") +
      geom_ribbon(aes(ymin = lower_laplacePoisson, ymax = upper_laplacePoisson), fill="red", alpha = 0.2, color = NA) +
      geom_ribbon(aes(ymin = lower_normalNB, ymax = upper_normalNB), fill="blue", alpha = 0.2, color = NA) +
      geom_ribbon(aes(ymin = lower_ZIP, ymax = upper_ZIP), fill="green", alpha = 0.2, color = NA) +
      geom_ribbon(aes(ymin = lower_ZINB, ymax = upper_ZINB), fill="purple", alpha = 0.2, color = NA) +
      geom_ribbon(aes(ymin = lower_linearTSlog, ymax = upper_linearTSlog), fill = "orange", alpha = 0.5, color = NA) +
      scale_color_manual(values = c(
        "TS-Poisson" = "red",
        "TS-NB" = "blue",
        "TS-ZIP" = "green",
        "TS-ZINB" = "purple",
        "Linear TS (log)" = "orange"
      )) +
      # need to adjust the scale of y axis according to the results
      scale_y_continuous(breaks = c(50, 100,150, 200,250,300, 350, 400)) +
      coord_cartesian(clip = "on", ylim = c(40, 400)) +  # Adjust clipping here, within the main chain
      theme(
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.background = element_blank(),  # Blank background
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(colour = "black", fill = NA) # Add the frame
      )
    val_regret = list(figure, data)
  }else{
    val_regret = data
  }
  return(val_regret)
}
