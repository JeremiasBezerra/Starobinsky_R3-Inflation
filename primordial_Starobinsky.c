int primordial_inflation_potential(
                                   struct primordial * ppm,
                                   double phi,
                                   double * V,
                                   double * dV,
                                   double * ddV
                                   ) {

  double e,de,dde,mu,dmu,ddmu,l,dl,ddl,p,dp,ddp;
  double phip, kappa;

  switch (ppm->potential) {

    /* V(phi)=polynomial in phi
       (Note: here “polynomial” is just a label — we actually
        substitute the potential corresponding to the Starobinsky model below.) */

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

    ppm->V0 = ppm->kappa*pow((1-exp(-sqrt(16*_PI_/3)*phip)),2);
    ppm->V1 = -2*ppm->kappa*pow(sqrt(16*_PI_/3),1)*(exp(-sqrt(16*_PI_/3)*phip)-1*exp(-2*sqrt(16*_PI_/3)*phip));
    ppm->V2 = -2*ppm->kappa*pow(sqrt(16*_PI_/3),2)*(exp(-sqrt(16*_PI_/3)*phip)-2*exp(-2*sqrt(16*_PI_/3)*phip));
    ppm->V3 = -2*ppm->kappa*pow(sqrt(16*_PI_/3),3)*(exp(-sqrt(16*_PI_/3)*phip)-4*exp(-2*sqrt(16*_PI_/3)*phip));
    ppm->V4 = -2*ppm->kappa*pow(sqrt(16*_PI_/3),4)*(exp(-sqrt(16*_PI_/3)*phip)-8*exp(-2*sqrt(16*_PI_/3)*phip));

    /* From here, CLASS will use these coefficients to compute V(phi), dV, ddV, etc.
    /* (Rest of the code: compute V, dV, ddV from the series — not shown here.) */
    ...
