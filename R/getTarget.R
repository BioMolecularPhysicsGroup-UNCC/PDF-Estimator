getTarget <- function (Ns, target = 70) {
  # getTargets calculates position-dependent threshold values about the mean,
  # according to a beta distribution with parameters k and (n + 1 - k), where
  # k is the position and n is the total number of positions.  These beta
  # distributions represent probability per position for sort order
  # statistics for a uniform distribution.
  
  if (!is.numeric(Ns)) {
    stop("Number of samples must be numeric")
  }
  
  if (!is.numeric(target) || target < 0 || target >= 100) {
    stop("target confidence percentage must be between 1 and 100")
  }
  
  
  above = matrix(0, 1, Ns);          
  below = matrix(0, 1, Ns);        
  
  if (Ns > 100000000) {
    print('maximum Ns = 100000000 = 100M. Must extend the limit');
  }
  Ns1 = Ns + 1;
  # ---------------------------------------------------- set up du - spacing
  np = ceiling(log10(Ns));
  du1000 = 0.001;
  duMIN = du1000 * (0.2) ^ np;                                        # heristic
  nipW = 1000 * (np + 1);       # # of integration points within WINDOW: heristic
  
  gRate = exp(log(du1000 / duMIN) / nipW);               # gRate = growth rate
  duRelative0 = matrix(0, 1, nipW);
  du = duMIN;
  for (i in 1:nipW) {
    duRelative0[i] = du;
    du = du * gRate;
  }
  nipMAX = 1000 + nipW + nipW + 1000;                     # => four sections
  kStart = 1000 + nipW;
  uList = matrix(0, 1, nipMAX);
  duRelative = matrix(0, 1, nipMAX);
  duRelative[1:1000] = du1000;
  duRelative[1001:kStart] = duRelative0[nipW:1];
  duRelative[(kStart + 1):(nipMAX - 1000)] = duRelative0;
  duRelative[(kStart + 1 + nipW):nipMAX] = du1000;
  bd0 = matrix(0, 1, nipMAX);                         # => beta distribution: CDF
  bd1 = matrix(0, 1, nipMAX);                  # => beta distribution: CDF' = PDF
  bd2 = matrix(0, 1, nipMAX);                # => beta distribution: CDF'' = PDF'
  # ------------------------- select the sort order statistic (SOS) k-values
  nHalf = ceiling(Ns / 2);                            # employ refection symmetry
  nksos = nHalf;              # # of k-values for SOS: True when nHalf <= 50
  if (nHalf <= 50) {     # note: 100*(0.025) = 2.5 => floor(100*0.025) = 2
    ksosArray = 1:nHalf;     # k-th SOS index, SOS => sort order statistics
  } else {                       # this part of the code requires nHalf >= 50 
    nksos = 0;
    k075 = floor(0.075 * Ns);
    k = 0;
    dk = 1;
    while (k < k075) {
      for (i16times in 1:16) {
        k = k + dk;
        if (k > k075) {
          break;
        } else {
          nksos = nksos + 1;
        }
      }
      dk = 2 * dk;
    }
    nksos = nksos + 9;
    ksosArray = matrix(0, 1, nksos);
    # --------------------------------------------------- populate k-values
    jk = 0;
    k = 0;
    dk = 1;
    while (k < k075) {
      for (i16times in 1:16) {
        k = k + dk;
        if (k > k075) {
          break;
        } else {
          jk = jk + 1;
          ksosArray[jk] = k;
        }
      }
      dk = 2*dk;
    }
    jk1 = jk + 1;
    jkLast = nksos - 1;
    ff = 0.10;                                          # initialize at 10#
    for (jk in jk1:jkLast) {
      ksosArray[jk] = floor(ff * Ns); 
      ff = ff + 0.05;                                    # increment by 5#
    }
    ksosArray[nksos] = nHalf;   # finish at 50# since nHalf = ceiling(0.50*Ns)
  }
  tag = -matrix(1, 1, nHalf);
  for (jjj in 1:nksos) {
    ksos = ksosArray[jjj];
    tag[ksos] = 1;
    alpha = ksos;
    beta = Ns - ksos + 1;
    alpha1 = alpha - 1;
    beta1 = beta - 1;
    uMode = (alpha - 1) / (alpha + beta - 2);
    uMode = max(1.0e-10, uMode);
    # -------------------------------------------- generate integration points
    k = kStart;
    uList[k] = uMode;
    while (uList[k] < 1) {
      k = k + 1;
      uList[k] = uList[k-1] + duRelative[k];
    }
    kMAX = k;
    uList[kMAX] = 1;
    # --------------
    k = kStart + 1;
    while (uList[k] > 0) {
      k = k - 1;
      uList[k-1] = uList[k] - duRelative[k];
    }
    kMIN = k;
    uList[kMIN] = 1.0e-25;
    u = uList[kMIN:kMAX];
    sMode = kStart - kMIN + 1;
    nip = kMAX - kMIN;
    # ========================================= determine range of integration
    sL = sMode;
    if (Ns < 100) {
      dsL = 1;
    } else {
      dsL = 20;
    }
    while (sL > 1) {
      sL = sL - dsL;
      if (sL < 1) {
        sL = 1;
        break;
      }
      rA = u[sL] / uMode;
      rB = (1 - u[sL]) / (1 - uMode);
      log10TestRatio = alpha1 * log10(rA) + beta1 * log10(rB);
      if (log10TestRatio < -20) {
        break;
      }
      if (log10TestRatio > -15) {
        dsL = 1;
      }
    }
    sR = sMode;
    if (Ns < 100) {
      dsR = 1;
    } else {
      dsR = 20;
    }
    while (sR < nip) {
      sR = sR + dsR;
      if (sR > nip) {
        sR = nip;
        break;
      }
      rA = u[sR] / uMode;
      rB = ( 1 - u[sR] )/( 1 - uMode );
      log10TestRatio = alpha1 * log10(rA) + beta1 * log10(rB);
      if (log10TestRatio < -20) {
        break;
      }
      if (log10TestRatio > -10) {
        dsR = 1;
      }
    }
    # ------------------------------------------------------ right propagation
    ln_rA = log(u[sL] / uMode);
    ln_rB = log((1 - u[sL]) / (1 - uMode));
    ln_bd1 = alpha1 * ln_rA + beta1 * ln_rB;               # using natural logs
    bd1[sL] = exp(ln_bd1);
    ln1 = ln_bd1 - ln_rA;         # => ln1 = alpha2*log(rA) + beta1*log(rB)
    ln2 = ln_bd1 - ln_rB;         # => ln2 = alpha1*log(rA) + beta2*log(rB)
    bd2[sL] = (alpha1 / uMode) * exp(ln1) - (beta1 / (1 - uMode)) * exp(ln2);
    bd0[sL] = 0;
    
    # -------------------------------------------------------------- integrate
    sLess1 = sL - 1;
    sL1 = sL + 1;
    for (s in sL1:sMode) {
      sLess1 = sLess1 + 1;
      du = u[s] - u[sLess1];
      duHalf = du / 2;
      ln_rA = log(u[s]/uMode);
      ln_rB = log((1 - u[s]) / (1 - uMode));
      ln_bd1 = alpha1 * ln_rA + beta1 * ln_rB;            # using natural logs
      bd1[s] = exp(ln_bd1);
      bd2[s] = -bd2[sLess1]  + (bd1[s] - bd1[sLess1]) / duHalf;
      bd0[s] =  bd0[sLess1]  + (bd1[s] + bd1[sLess1]) * duHalf + (bd2[sLess1] - bd2[s]) * du * du / 8;
    } 
    topValue = bd0[sMode];
    # ------------------------------------------------------ left propagation
    ln_rA = log(u[sR] / uMode);
    ln_rB = log((1 - u[sR]) / (1 - uMode));
    ln_bd1 = alpha1 * ln_rA + beta1 * ln_rB;               # using natural logs
    bd1[sR] = exp(ln_bd1);
    ln1 = ln_bd1 - ln_rA;         # => ln1 = alpha2*log(rA) + beta1*log(rB)
    ln2 = ln_bd1 - ln_rB;         # => ln2 = alpha1*log(rA) + beta2*log(rB)
    bd2[sR] = (alpha1 / uMode) * exp(ln1) - (beta1 / (1 - uMode) ) * exp(ln2);
    bd0[sR] = 0;
    # -------------------------------------------------------------- integrate
    sLess1 = sR;
    sMode1 = sMode + 1;
    for (s in sR:sMode1) {
      sLess1 = sLess1 - 1;
      du = u[s] - u[sLess1];
      duHalf = du / 2;
      ln_rA = log(u[sLess1] / uMode);
      ln_rB = log((1 - u[sLess1]) / (1 - uMode));
      ln_bd1 = alpha1 * ln_rA + beta1 * ln_rB;            # using natural logs
      bd1[sLess1] = exp(ln_bd1);
      bd2[sLess1] = -bd2[s]  + (bd1[s] - bd1[sLess1]) / duHalf;
      bd0[sLess1] =  bd0[s]  - (bd1[s] + bd1[sLess1]) * duHalf + (bd2[s] - bd2[sLess1]) * du * du / 8;
    }
    # ----------------------------------------------------- end of integration
    
    # -------------------------------- match left and right integrals at sMode
    shift_bd0 = topValue - bd0[sMode];
    bd0[sMode:sR] = bd0[sMode:sR] + shift_bd0;
    # ---------------------------------------------------------- normalize CDF
    scaleFactor = 1 / bd0[sR];
    bd0 = bd0 * scaleFactor;   
   
    Fbelow = (100 - target) / 200 
    Fabove = 1 - Fbelow
    max_dFbelow = 100;
    max_dFabove = 100;
    sbelow = sL;
    sabove = sR;
    for (s in sL:sR) {
      F = bd0[s];
      dFbelow = abs(F - Fbelow);
      if (dFbelow < max_dFbelow) {
        max_dFbelow = dFbelow;
        sbelow = s;
      }
      dFabove = abs(F - Fabove);
      if (dFabove < max_dFabove) {
        max_dFabove = dFabove;
        sabove = s;
      }
    }
    # --------------------------- define residuals to be w.r.t. mean quantiles
    mu = (1:Ns) / Ns1;                                # => mean SOS positions
    z = (u - mu[ksos]) * sqrt(Ns+2);                   # => scaled residual
    # ------------------------------------------------ quadratic interpolation
    x = Fbelow
    s = sbelow;
    x0 = bd0[s-1];
    x1 = bd0[s];
    x2 = bd0[s+1];
    y0 = z[s-1];
    y1 = z[s];
    y2 = z[s+1];
    c0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
    c1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
    c2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    below[ksos] = y0 * c0 + y1 * c1 + y2 * c2;
    x = Fabove
    s = sabove;
    x0 = bd0[s-1];
    x1 = bd0[s];
    x2 = bd0[s+1];
    y0 = z[s-1];
    y1 = z[s];
    y2 = z[s+1];
    c0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
    c1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
    c2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    above[ksos] = y0 * c0 + y1 * c1 + y2 * c2;
  } 
  # ------------------------------------------ interpolate lemon drop points
  j2 = 3;
  for (k in 1:nHalf) {
    if (tag[k] < 0) {
      x = mu[k];
      while (ksosArray[j2] < k) { 
        j2 = j2 + 1;
        j1 = j2 - 1;
        j0 = j1 - 1;
        k0 = ksosArray[j0];
        k1 = ksosArray[j1];
        k2 = ksosArray[j2];
        dk21 = k2 - k1;
        while ((k1 - k0) < floor(0.92 * dk21)) {
          j0 = j0 - 1;
          k0 = ksosArray[j0];
        }
        # ------------------------------------------------- generate anchor points
        x0 = mu[k0];
        x1 = mu[k1];
        x2 = mu[k2];
        #-----------------
        y0above = above[k0];
        y0below = below[k0];
        #-----------------
        y1above = above[k1];
        y1below = below[k1];
        #-----------------
        y2above = above[k2];
        y2below = below[k2];
      }
      c0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
      c1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
      c2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
      # ---------------------------------------------
      below[k] = c0 * y0below + c1 * y1below + c2 * y2below;
      above[k] = c0 * y0above + c1 * y1above + c2 * y2above;
    }
  }
  # --------------------------------------------------- copy reflected parts
  nHalf1 = nHalf + 1;
  for (k in nHalf1:Ns) {
    above[k] = -below[Ns1-k];
    below[k] = -above[Ns1-k];
  }
  bounds = matrix(c(above, below), nrow = length(below), ncol = 2)
  return (bounds)
}