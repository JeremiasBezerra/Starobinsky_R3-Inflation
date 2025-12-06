int primordial_inflation_potential(
                                   struct primordial * ppm,
                                   double phi,
                                   double * V,
                                   double * dV,
                                   double * ddV
                                   ) {

  double e,de,dde,mu,dmu,ddmu,l,dl,ddl,p,dp,ddp;
  double phip, kappa, alpha;

  switch (ppm->potential) {

    /* V(phi)=polynomial in phi */
  case polynomial:

    kappa = ppm->kappa;
    alpha = ppm->alpha;
    
    phip = phip_val;
    /* -----------------------------------------------------------------------------
 *  φ_p: inflaton value at horizon crossing
 *
 *  The potential-expansion below requires the value of the inflaton field
 *  φ_p(N), evaluated at the e-fold time corresponding to the mode that exits
 *  the Hubble radius. Physically, φ_p is defined implicitly through
 *
 *        N_k = Integrate H dt
 *              from φ_end up to φ_p .
 *
 *  In the literature this relation is often computed using the slow-roll
 *  approximation, which yields closed analytic expressions for φ_p(N).
 *
 *  In our case, however, we avoid slow-roll and evaluate the integral
 *  numerically in an external Python script that solves the background
 *  equations exactly. The result is stored in a table φ_p(N, α), and here we
 *  simply read the appropriate value. Thus this step amounts to replacing the
 *  usual slow-roll estimate by the numerically computed trajectory φ(N).
 *
 *  In summary:
 *    – the code below depends on φ_p;
 *    – φ_p is obtained by inverting the exact e-fold integral;
 *    – it can come either from a slow-roll expression or from a lookup table.
 * ----------------------------------------------------------------------------- */

    /* Now we plug φ = phip into the analytic Taylor-series coefficients
       of the Starobinsky potential around that field value. This allows CLASS
       to compute V, dV, ddV etc., using the local expansion coefficients V0, V1, … defined below. */
       
    ppm->V0 = -(exp(-8 * phip * sqrt(M_PI/3.0)) * (-1 + sqrt(1 + 4 * alpha * (-1 + exp(4 * phip * sqrt(M_PI/3.0))))) * (-1 - 8 * alpha * (-1 + exp(4 * phip * sqrt(M_PI/3.0))) +sqrt(1 + 4 * alpha * (-1 + exp(4 * phip * sqrt(M_PI/3.0))))) *kappa) / (12 * alpha * alpha);
    
    ppm->V1 = -1 * (1.0/6.0) * kappa / (alpha * alpha) * sqrt(16 * M_PI / 3.0) * exp(-2 * phip * sqrt(16 * M_PI / 3.0)) * (1 - sqrt(1 - 4 * alpha * (1 - exp(phip * sqrt(16 * M_PI / 3.0))))) * (1 - 2 * alpha * (4 - exp(phip * sqrt(16 * M_PI / 3.0))) - sqrt(1 - 4 * alpha * (1 - exp(phip * sqrt(16 * M_PI / 3.0)))));
    
    ppm->V2 = (2 * exp(-8 * phip * sqrt(M_PI/3.0)) * (1 - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * (8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) - (8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0)) / sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * sqrt(M_PI/3.0) * kappa) / (3 * alpha * alpha) - (16 * exp(-4 * phip * sqrt(M_PI/3.0)) * (1 - 2 * alpha * (4 - exp(4 * phip * sqrt(M_PI/3.0))) - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * M_PI * kappa) / (9 * alpha * sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) - (16 * exp(-8 * phip * sqrt(M_PI/3.0)) * (1 - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * (1 - 2 * alpha * (4 - exp(4 * phip * sqrt(M_PI/3.0))) - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * M_PI * kappa) / (9 * alpha * alpha);
    
    ppm->V3 = -1 * kappa * ((128 * exp(-8 * phip * sqrt(M_PI/3.0)) * (1 - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * (1 - 2 * alpha * (4 - exp(4 * phip * sqrt(M_PI/3.0))) - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * pow(M_PI, 3.0/2.0)) / (9 * sqrt(3.0) * alpha * alpha) + (2 * exp(-8 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) / (3 * alpha * alpha)) * (- ( (16 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * (8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) - 8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) / sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * sqrt(M_PI/3.0)) / sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0)))) ) + (1 - 2 * alpha * (4 - exp(4 * phip * sqrt(M_PI/3.0))) - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * ( (64 * alpha * alpha * exp(8 * phip * sqrt(M_PI/3.0)) * M_PI) / (3 * pow(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))), 3.0/2.0)) - (32 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * M_PI) / (3 * sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) ) + (1 - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) *
(32.0/3.0 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * M_PI + (64 * alpha * alpha * exp(8 * phip * sqrt(M_PI/3.0)) * M_PI) / (3 * pow(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))), 3.0/2.0)) - (32 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * M_PI) /  (3 * sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) )) - (32 * exp(-8 * phip * sqrt(M_PI/3.0)) * M_PI / (9 * alpha * alpha)) * ((1 - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) *
(8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) - 8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * sqrt(M_PI/3.0) / sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) - (8 * alpha * exp(4 * phip * sqrt(M_PI/3.0)) * (1 - 2 * alpha * (4 - exp(4 * phip * sqrt(M_PI/3.0))) - sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))) * sqrt(M_PI/3.0)) / sqrt(1 - 4 * alpha * (1 - exp(4 * phip * sqrt(M_PI/3.0))))));
    
    ppm->V4 = -(1024.0 * exp(-8.0 * phip * sqrt(M_PI/3.0)) * (1.0 - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * (1.0 - 2.0*alpha*(4.0 - exp(4.0 * phip * sqrt(M_PI/3.0))) - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * M_PI*M_PI * kappa) / (27.0 * alpha*alpha) - (16.0 * exp(-8.0 * phip * sqrt(M_PI/3.0)) * M_PI / (3.0 * alpha*alpha)) * ((-(16.0 * alpha * exp(4.0 * phip * sqrt(M_PI/3.0)) * (((8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0)))) - ((((8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0)))) / sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))) *sqrt(M_PI/3.0)) / (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) + (1.0 - 2.0*alpha*(4.0 - exp(4.0 * phip * sqrt(M_PI/3.0))) - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * ((64.0 * alpha*alpha * exp(8.0 * phip * sqrt(M_PI/3.0)) * M_PI) / (3.0 * pow(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))),3)) - (32.0 * alpha * exp(4.0 * phip * sqrt(M_PI/3.0)) * M_PI) / (3.0 * (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))) + (1.0 - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * (32.0/3.0 * alpha * exp(4.0 * phip * sqrt(M_PI/3.0)) * M_PI + (64.0 * alpha*alpha * exp(8.0 * phip * sqrt(M_PI/3.0)) * M_PI) / (3.0 * pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),3)) - (32.0 * alpha * exp(4.0 * phip * sqrt(M_PI/3.0)) * M_PI) / (3.0 * (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))))) * kappa  + (2.0 * exp(-8.0 * phip * (sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0)) / (3.0 * alpha*alpha)) * ((3.0 * ( (8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0))) - (((8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0)))) / sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * ((64.0*alpha*alpha*exp(8.0 * phip * sqrt(M_PI/3.0))*M_PI)/(3.0 * pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),3)) - (32.0*alpha*exp(4.0 * phip * sqrt(M_PI/3.0))*M_PI)/(3.0*(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))))+ (-(8.0 * alpha * exp(4.0 * phip * sqrt(M_PI/3.0)) * sqrt(3.0*M_PI) * (32.0/3.0*alpha*exp(4.0 * phip * sqrt(M_PI/3.0))*M_PI + (64.0*alpha*alpha*exp(8.0 * phip * sqrt(M_PI/3.0))*M_PI)/(3.0*pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),3)) - (32.0*alpha*exp(4.0 * phip * sqrt(M_PI/3.0))*M_PI)/(3.0*(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))))/(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) + ((1.0 - 2.0*alpha*(4.0 - exp(4.0 * phip * sqrt(M_PI/3.0))) - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))* (-(512.0*pow(alpha,3) * pow((exp(4.0 * phip * sqrt(M_PI/3.0))),1) * pow(M_PI,1.5)) / (sqrt(3.0) * pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),5)) + (256.0*alpha*alpha*exp(8.0 * phip * sqrt(M_PI/3.0))*pow(M_PI,1.5)) / (sqrt(3.0)*pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),3)) - (128.0*alpha*(exp(4.0 * phip * sqrt(M_PI/3.0)))*pow(M_PI,1.5))/(3.0*sqrt(3.0)*(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))))) + ((1.0 - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * ((128.0*alpha*(exp(4.0 * phip * sqrt(M_PI/3.0)))*pow(M_PI,1.5))/(3.0*sqrt(3.0)) - (512.0*pow(alpha,3)*pow((exp(4.0 * phip * sqrt(M_PI/3.0))),1)*pow(M_PI,1.5))/(sqrt(3.0)*pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),5)) + (256.0*alpha*alpha*exp(8.0 * phip * sqrt(M_PI/3.0))*pow(M_PI,1.5))/(sqrt(3.0)*pow((sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))),3)) - (128.0*alpha*(exp(4.0 * phip * sqrt(M_PI/3.0)))*pow(M_PI,1.5))/(3.0*sqrt(3.0)*(sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))))) * kappa + ((128.0 * exp(-8.0 * phip * (sqrt(M_PI/3.0))) * pow(M_PI,1.5) / (3.0 * sqrt(3.0) * alpha*alpha)) * (((1.0 - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * ( (8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0))) - (((8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (sqrt(M_PI/3.0)))) / sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) - 8.0 * alpha * (exp(4.0 * phip * sqrt(M_PI/3.0))) * (1.0 - 2.0*alpha*(4.0 - exp(4.0 * phip * sqrt(M_PI/3.0))) - (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0))))))* (sqrt(M_PI/3.0)) / (sqrt(1.0 - 4.0*alpha*(1.0 - exp(4.0 * phip * sqrt(M_PI/3.0)))))) * kappa));
    
    /* From here, CLASS will use these coefficients to compute V(phi), dV, ddV, etc.
    /* (Rest of the code: compute V, dV, ddV from the series — not shown here.) */
    ...
