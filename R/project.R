##' Project forward
##'
##' .. content for \details{} ..
##' @title
##' @param I0 initial incidence. T X 1 matrix.
##' @param R Vector of length p. R values
##' corresponding to each of the p change points
##' specified in tau
##' @param tau a vector of time points at which R
##' changes
##' @param si serial interval distribution. If
##' the length of the vector is less than
##' the nrow(I0) + n_days, it is padded with 0s.
##' @param n_days number of days for which to
##' project forward.
##' @param n_sim number of trajectories desired
##' @return (t + n_days) X n_sim matrix. Each
##' column is a single simulation.
##' @author Sangeeta Bhatia

project <- function(I0, R, tau, si, n_days, n_sim) {
    if (length(R) != length(tau)) {
        stop(
            "Number of R values should be the same
             as the number of changepoints."
        )
    }

    out <- matrix(0, nrow = n_days, ncol = n_sim)
    ## Repeat I0 n_sim times so that we can rbind.
    I0rep <- I0[, rep(seq_len(ncol(I0)), n_sim)]
    I0rep <- matrix(I0rep, nrow = nrow(I0), ncol = n_sim)

    out   <- rbind(I0rep, out)
    start <- nrow(I0) + 1
    end   <- nrow(out)
    if (length(si) < end)
         ws <- rev(c(si, rep(0, end - length(si))))
    else ws <- rev(si)

    tau <- c(0, tau)
    change_at <- diff(tau)

    rvec <- mapply(rep, x = R, times = change_at)
    rvec <- unlist(rvec)
    for (index in start:end) {
          i_t <- out[1:index, ]
          w_t <- utils::tail(ws, index)

          r_t <- rvec[index - start + 1]
          mu  <- (w_t %*% i_t) * r_t
          out[index, ] <- rpois(n_sim, mu)
      }


  out

}
