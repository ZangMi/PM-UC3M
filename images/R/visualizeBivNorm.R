

#' @title Visualize the joint ensity, marginals and conditionals of a bivariate normal
#'
#' @description Generates a plot illustrating the joint/marginal/conditional densities of a bivariate normal.
#'
#' @param data a \code{data.frame} containing \code{X1} and \code{X2} as variables.
#' @param mu1,mu2,sigma1,sigma2,rho parameters defining a bivariate normal.
#' @param nGrid (approximate) number of cuts used to produce the XYZ grid. Computed to match the effective \code{nticks} in \code{\link[graphics]{persp}}.
#' @inheritParams graphics::persp
#' @param basalScatter plot the \code{data} in the XY plane?
#' @param marginals draw the marginals and their means?
#' @param conditionals if different from \code{0}, draws the conditional density of the \code{conditional}-th (\code{1} or \code{2}) variable given the value \code{xc} for the other.
#' @param xc value to condition at.
#' @return
#' Nothing. The function is called to produce a plot.
#' @examples
#' # Generate data
#' library(mvtnorm)
#' n <- 5e2
#' mu1 <- mu2 <- 0
#' mu <- c(mu1, mu2)
#' sigma1 <- sigma2 <- 1
#' rho <- 0.75
#' sigma <- rbind(c(sigma1^2, sigma1 * sigma2 * rho), 
#'                c(sigma1 * sigma2 * rho, sigma2^2))
#' set.seed(4151460)
#' data <- data.frame(rmvnorm(n = n, mean = mu, sigma = sigma))
#' colnames(data) <- c("X1", "X2")
#' 
#' # With marginals and conditionals
#' par(mar = rep(2, 4))
#' visualizeBivNorm(data = data, xc = -2, conditionals = 2)
#'
#' # Dull plot
#' visualizeBivNorm(data = data, conditionals = FALSE, marginals = FALSE)
#' @author Eduardo García-Portugués (\email{edgarcia@est-econ.uc3m.es}).
#' @export
visualizeBivNorm <- function(data, mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 1, 
                             rho = 0.5, nGrid = 5, theta = -30, phi = 20,
                             basalScatter = TRUE, marginals = TRUE, 
                             conditionals = 2, xc = 1) {
  
  # Create and plot initial XY-grid
  sdX1 <- sd(data$X1)
  sdX2 <- sd(data$X2)
  g1 <- pretty(c(min(data$X1), max(data$X1)), nGrid)
  g2 <- pretty(c(min(data$X2), max(data$X2)), nGrid)
  
  # Evaluate joint density
  mu <- c(mu1, mu2)
  sigma <- rbind(c(sigma1^2, sigma1 * sigma2 * rho), 
                 c(sigma1 * sigma2 * rho, sigma2^2))
  x1 <- pretty(g1, 1e2)
  x2 <- pretty(g2, 1e2)
  x1 <- c(g1[1], x1[x1 > g1[1] & x1 < g1[length(g1)]], g1[length(g1)])
  x2 <- c(g2[1], x2[x2 > g2[1] & x2 < g2[length(g2)]], g2[length(g2)])
  lx1 <- length(x1)
  lx2 <- length(x2)
  x12 <- expand.grid(x1 = x1, x2 = x2)
  require(mvtnorm)
  dens <- dmvnorm(x = x12, mean = mu, sigma = sigma)

  # Evaluate marginal densities
  dens1 <- dnorm(x = x1, mean = mu1, sd = sigma1)
  dens2 <- dnorm(x = x2, mean = mu2, sd = sigma2)
  gDens <- seq(0, max(dens, dens1, dens2), length = nGrid)
  
  # Evaluate conditional densities
  if (conditionals == 1) {
    
    muc <- mu1 + sigma1/sigma2 * rho * (xc - mu2)
    sdc <- sqrt(sigma1^2 * (1 - rho^2))
    densc <- dnorm(x = x1, mean = muc, sd = sdc)
    gDens <- seq(0, max(dens, dens1, dens2, densc), length = nGrid)
    
  } else if (conditionals == 2) {
    
    muc <- mu2 + sigma2/sigma1 * rho * (xc - mu1)
    sdc <- sqrt(sigma2^2 * (1 - rho^2))
    densc <- dnorm(x = x2, mean = muc, sd = sdc)
    gDens <- seq(0, max(dens, dens1, dens2, densc), length = nGrid)
    
  }

  # Plot data and regression
  require(plot3D)
  panelFirst <- function(pmat) {
    
    # Plot points projection in the basal plane
    if (basalScatter) {
      
      XY <- trans3D(data$X1, data$X2, rep(0, n), pmat = pmat)
      scatter2D(XY$x, XY$y, pch = 16, cex = 0.5, add = TRUE,
                colkey = FALSE, col = 1)
      
    }
    
    # Grid lines
    for (x in g1) {
      
      lines(trans3D(rep(x, lx2), x2, 
                    dmvnorm(x = cbind(x, x2), mean = mu, sigma = sigma), pmat),
            col = gray(0.5))
      
    }
    for (x in g2) {
      
      lines(trans3D(x1, rep(x, lx1), 
                    dmvnorm(x = cbind(x1, x), mean = mu, sigma = sigma), pmat),
            col = gray(0.5))
      
    }
    
    # Plot marginals
    if (marginals) {
      
      # Densities
      lines(trans3D(x1, rep(x2[lx2], lx1), dens1, pmat), col = 3, lwd = 2)
      lines(trans3D(rep(x1[lx1], lx2), x2, dens2, pmat), col = 3, lwd = 2)
      
      # # Observations
      # points(trans3D(data$X1, rep(x2[lx2], nrow(data)), 0, pmat), 
      #        pch = 16, cex = 0.5, col = 1)
      # points(trans3D(rep(x1[lx1], nrow(data)), data$X2, 0, pmat), 
      #        pch = 16, cex = 0.5, col = 1)
      
      # Means
      lines(trans3D(rep(mu1, 3), c(x2[lx2], x2[lx2], x2[1]), 
                    c(dnorm(mu1, mean = mu1, sd = sigma1), 0, 0), pmat), 
            col = 2)
      lines(trans3D(c(x1[lx1], x1[lx1], x1[1]), rep(mu2, 3), 
                    c(dnorm(mu2, mean = mu2, sd = sigma2), 0, 0), pmat), 
            col = 2)
      
    }
    
    # Plot conditional
    if (conditionals) {
      
      # Density
      if (conditionals == 1) {
        
        lines(trans3D(x1, rep(xc, lx1), densc, pmat), col = "orange", lwd = 2)
        lines(trans3D(x1, rep(xc, lx1), 
                      dmvnorm(x = cbind(x1, xc), mean = mu, sigma = sigma), 
                      pmat), col = "orange", lwd = 2)
        lines(trans3D(x1, rep(xc, lx1), 0, pmat), col = "orange")
        
      } else if (conditionals == 2) {
          
        lines(trans3D(rep(xc, lx2), x2, densc, pmat), col = "orange", lwd = 2)
        lines(trans3D(rep(xc, lx2), x2, 
                      dmvnorm(x = cbind(xc, x2), mean = mu, sigma = sigma), 
                      pmat), col = "orange", lwd = 2)
        lines(trans3D(rep(xc, lx2), x2, 0, pmat), col = "orange")

      }

      # Mean
      if (conditionals == 1) {
        
        dc <- dnorm(x = muc, mean = muc, sd = sdc)
        lines(trans3D(rep(muc, 2), rep(xc, 2), c(0, dc), pmat), col = "orange")
        points(trans3D(muc, xc, 0, pmat), pch = 16, col = "orange", cex = 1)
        
      } else if (conditionals == 2) {
        
        dc <- dnorm(x = muc, mean = muc, sd = sdc)
        lines(trans3D(rep(xc, 2), rep(muc, 2), c(0, dc), pmat), col = "orange")
        points(trans3D(xc, muc, 0, pmat), pch = 16, col = "orange", cex = 1)
        
      }

    }
    
  }
  gridMat <- scatter3D(mu1, mu2, 0, pch = 19, theta = theta, phi = phi,
                       bty = "g", axes = TRUE, colkey = FALSE, col = 2,
                       xlim = range(g1), ylim = range(g2), zlim = range(gDens),
                       panel.first = panelFirst, nticks = nGrid, cex = 1,
                       ticktype = "detailed", xlab = "x1", ylab = "x2", 
                       zlab = "Density")

  M <- mesh(x1, x2)
  surf3D(x = M$x, y = M$y, z = matrix(dens, nrow = lx1, ncol = lx2), 
         col = "lightblue", alpha = 0.2, border = NA, add = TRUE)
  
}


