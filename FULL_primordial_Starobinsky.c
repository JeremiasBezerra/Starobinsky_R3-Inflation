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
  case polynomial:
  {
    static int already_ran = _FALSE_;
    static double N_table[100000];
    static double phi_table[100000];
    static int phi_table_size;

    /* Read pre-computed background trajectory table only once. 
       This table is generated externally (by a Python script without slow-roll approximation),
       containing N (number of efolds) vs φ (inflaton field value) pairs. */
    if (already_ran == _FALSE_) {
      already_ran = _TRUE_;

      FILE *file = fopen("/home/dell/Documentos/COSMO/classy_starobinsky/source/phitrajectory.dat", "r");
      if (file == NULL) {
        fprintf(stderr, "Error opening file phitrajectory.dat\n");
        exit(1);
      }

      phi_table_size = 0;
      while (fscanf(file, "%lf %lf", &N_table[phi_table_size], &phi_table[phi_table_size]) == 2) {
        phi_table_size++;
        if (phi_table_size > 100000) {
          fclose(file);
          fprintf(stderr, "Phi trajectory table exceeds maximum size (100000)\n");
          exit(1);
        }
      }

      fclose(file);
    }

    /* Determine φ_p (phip): the value of the inflaton φ when the mode with
       wavenumber k crosses the horizon (i.e. at a given number of e-folds N_query). */
    const double N_query = ppm->nfolds;
    double phip_val = 0.0;

    if (N_query < N_table[0]) {
      phip_val = phi_table[0];
    }
    else if (N_query > N_table[phi_table_size - 1]) {
      phip_val = phi_table[phi_table_size - 1];
    }
    else {
      for (int i = 1; i < phi_table_size; i++) {
        if (N_table[i] > N_query) {
          /* use the φ value just before N exceeds N_query */
          phip_val = phi_table[i - 1];
          break;
        }
      }
    }

    phip = phip_val;  /* φ at horizon crossing for the given N_query */
  }
  /* ---- End of phip block ---- */

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
