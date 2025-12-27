


#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#Feb, 2017


##### Robust Mean and Covariance (Algorithm 2) ##########
###
### Input: (W_m, b_m) for m = 1, ..., M
###    Concatenate: u_m <- [W_m, b_m], for all m
###    Initialize: mu <- mean(u_m), C <- cov(u_m)
###    repeat:
###       d_m <- sqrt((u_m - mu)^T C^-1 (u_m-mu))
###       w_m <- phi(d_m)/d_m
###       Update:
###          mu_new <- sum_m(w_m U_m)/sum_m(w_m)
###          C_new <- (sum_m(w_m^2(u_m - mu_new)(u_m - mu_new)^T))/sum_m(w_m^2-1)
###    until: stopping conditions are satisfied
### Output: mu = [w_0, b_0], C_0 = C(1:d, 1:d)
###
############################################################
### Lee, Gyemin, L Stoolman, and C Scott.
### "Transfer Learning for Auto-Gating of Flow Cytometry Data.”
### JMLR(workshop), 2012, 155–66.
############################################################
###
###   Using the GSE package's HuberPairwise function,
###    we estimate the robust mean and covariance.
###    A tuning constant c0 = 1.345 is used in the huber pairwise estimation function
###
############################################################
###  Alqallaf, F.A., Konis, K. P., R. Martin, D., Zamar, R. H.
###  "Scalable Robust Covariance and Correlation Estimates for Data Mining."
###  In Proceedings of the Seventh ACM SIGKDD International Conference
###     on Knowledge Discovery and Data Mining. Edmonton. 2002.
############################################################


#########
######
###

#' Estimate robust hyperplane moments across source models
#'
#' @description
#' Computes robust mean and covariance estimates for the baseline SVM
#' hyperplanes returned by [alg1_baselineClass()]. The routine falls back to
#' empirical moments when robust estimation fails and augments the result with
#' additional derived quantities.
#'
#' @param alg1_result_baselineSVM Matrix of hyperplane parameters from
#'   [alg1_baselineClass()] where each row represents a source model.
#' @param print2screen Logical; print progress updates.
#'
#' @return A list containing robust/simple means, covariance matrices, and
#'   normalized vectors used by downstream algorithms.
#' @export
alg2_rob_meanNCov <- function(alg1_result_baselineSVM, print2screen = F){
  #alg1_result_baselineSVM <- alg1_res$baselineSVM; print2screen = T


  if(is.null(alg1_result_baselineSVM)) print("alg1_result_baselineSVM is Null") else {

    if(nrow(alg1_result_baselineSVM)>1) {

      if(print2screen) print("starting Alg 2 to obtain robust mean and covariance of SVM hyperplane parameters")


      #baselineSVM is a matrix of concatination of the normal vector and the y-intercept (bias)
      baselineSVM <- alg1_result_baselineSVM



      #INITIALIZE: C is covariance and U is mean

      C <- cov(baselineSVM); C
      try.U <- try(CovEM(baselineSVM), silent = T)
      if(!(class(try.U)=="try-error")) {
        res <- CovEM(baselineSVM) #Gaussian MLE of mean and covariance
        U <- getLocation(res); #equal to getting mean as U <- sapply(baselineSVM, mean); U
        names(U) <- colnames(C); U
        colnames(C) <- names(U)
        rownames(C) <- names(U)
      } else {
        U <- as.vector(mode="numeric", colMeans(baselineSVM))
        names(U) <- colnames(C); U
      }

      #C <- getScatter(res);C # different to cov()


      # respmh <- GSE::partial.mahalanobis(baselineSVM, mu=U, S=C)
      # plot(respmh, which="index")
      # getDist(respmh)

      #huber pairwise estimation; using prev. built for now
      #c0 is the tuning constant for the huber function.
      #c0=0 would yield QC. Default is c0=1.345


      resHub <-try(HuberPairwise(as.matrix(baselineSVM), psi=c("huber"), c0=1.345, computePmd=TRUE), silent = T)

      if(!(class(try.U)=="try-error")) {
        #getDist(resHub)
        #plot(resHub, which="index")

        # resHub@R #correlation matrix
        # resHub@mu #the mean weighted by the Huber loss function, a robust loss fx
        # resHub@pmd #partial mahalanobis
        # resHub@S #the cov weighted by the Huber loss function, a robust loss fx
        # resHub@x


        #OUTPUT:
        #U #the simple mean
        U_robust <- resHub@mu; U_robust
        #C #the simple cov
        C_robust <- resHub@S; colnames(C_robust) <- colnames(C);
        rownames(C_robust) <- colnames(C); C_robust

        remove(res)

      } else {

        #OUTPUT:
        #U #the simple mean
        U_robust <- U
        #C #the simple cov
        C_robust <- C;

      }




      alg2_CalcMore_res <- alg2_CalcMore(list(U_simple=U, U_robust=U_robust, C_simple=C, C_robust=C_robust))

      names(alg2_CalcMore_res$U_robust_norm) <- names(baselineSVM)

      if(print2screen) print("Alg 2 has completed!")
    }

    if(nrow(alg1_result_baselineSVM) == 1) {





    }






    return(alg2_CalcMore_res)

  }

}


#' Compute derived statistics for Algorithm 2 results
#'
#' @description
#' Post-processes the outputs from [alg2_rob_meanNCov()] to derive normalized
#' hyperplane parameters, eigenvalues, and helper vectors that guide later RTL
#' steps.
#'
#' @param alg2_res List output from [alg2_rob_meanNCov()].
#'
#' @return The augmented `alg2_res` list with additional normalization values.
#' @export
alg2_CalcMore <- function(alg2_res) {

  #alg2_res = list(U_simple=U, U_robust=U_robust, C_simple=C, C_robust=C_robust)
  tempC <- alg2_res$C_robust[1:(nrow(alg2_res$C_robust)-1),1:(ncol(alg2_res$C_robust)-1)]

  alg2_res$v_0 <- eigen(tempC)$values
  alg2_res$U_robust # <w_0, b_0>
  alg2_res$w_euc_mag <- norm(alg2_res$U_robust[-length(alg2_res$U_robust)], type="2") #euclidean norm

  #w_t,b_t
  alg2_res$U_robust_norm <- alg2_res$U_robust/alg2_res$w_euc_mag



  #orthonormalize v_0 with respect to w_0
  alg2_res$v_t <- orthnormal(rbind(as.numeric(alg2_res$v_0[]),
                                   as.numeric(alg2_res$U_robust_norm[-length(alg2_res$U_robust_norm)])))

  alg2_res$v_t_norm <- alg2_res$v_t[1]#/alg2_res$w_euc_mag

  return(alg2_res)
}




