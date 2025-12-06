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
    
{
  const double alpha_query = ppm->alpha;
  const double N_query     = ppm->nfolds;

  /* -------------------------------------------------------------------------
   * R^3 inflation: lookup of φ_p(N, α)
   *
   * φ_p is the inflaton field value at horizon crossing of the mode associated
   * with the comoving wavenumber k(N). Physically, this corresponds to the
   * field value N e-folds before the end of inflation.
   *
   * Instead of using slow-roll approximations such as:
   *       N ≈ ∫ dφ / sqrt(2 ε_V)
   * we solve the exact background equations in an external Python routine.
   * This produces numerical tables of the form:
   *       (N, α, φ)
   * stored in multiple files, each one covering a different range of α.
   *
   * Here we simply read φ from those tables. This ensures full numerical
   * accuracy for the R^3 background evolution without relying on slow-roll.
   *
   * Because reading the tables is expensive, we implement a small cache:
   * if the same (N, α) pair is queried again, we reuse the stored φ_p.
   * ------------------------------------------------------------------------- */

  // --- cache to avoid recomputing φ_p ---
  static int already_ran = _FALSE_;
  static double phip_cached = 0.0;
  static double alpha_cached = 0.0;
  static double N_cached = 0.0;

  // If both N and α are identical to the last query, reuse φ_p
  if (already_ran && N_query == N_cached && alpha_query == alpha_cached) {

    phip = phip_cached;

  } else {

    /* ---------------------------------------------------------------------
     * File selection:
     * We store the numerical φ(N, α) table in 13 different files, each one
     * covering a specific α-range. We determine here which file to load.
     * --------------------------------------------------------------------- */

    int part_to_load = -1;
    if      (alpha_query >= -3e-4 && alpha_query < -1e-4) part_to_load = 0;
    else if (alpha_query >= -1e-4 && alpha_query < -1e-5) part_to_load = 1;
    else if (alpha_query >= -1e-5 && alpha_query < -1e-6) part_to_load = 2;
    else if (alpha_query >= -1e-6 && alpha_query < -1e-7) part_to_load = 3;
    else if (alpha_query >= -1e-7 && alpha_query < -1e-8) part_to_load = 4;
    else if (alpha_query >= -1e-8 && alpha_query < -1e-9) part_to_load = 5;
    else if (alpha_query >= -1e-9 && alpha_query < -1e-10) part_to_load = 6;
    else if (alpha_query >=  1e-10 && alpha_query <  1e-9) part_to_load = 7;
    else if (alpha_query >=  1e-9  && alpha_query <  1e-8) part_to_load = 8;
    else if (alpha_query >=  1e-8  && alpha_query <  1e-7) part_to_load = 9;
    else if (alpha_query >=  1e-7  && alpha_query <  1e-6) part_to_load = 10;
    else if (alpha_query >=  1e-6  && alpha_query <  1e-5) part_to_load = 11;
    else if (alpha_query >=  1e-5  && alpha_query <  1e-4) part_to_load = 12;
    else {
      fprintf(stderr, "Alpha outside tabulated range: %e\n", alpha_query);
      exit(1);
    }

    /* Build filename corresponding to the α-range */
    char filename[512];
    sprintf(filename, "/home/dell/Documentos/COSMO/code/R3_Inversao_Tabela/phik_table_R3_Part%d.dat",
            part_to_load);

    FILE *file = fopen(filename, "r");
    if (!file) {
      fprintf(stderr, "Error opening file %s\n", filename);
      exit(1);
    }

    /* ---------------------------------------------------------------------
     * Step 1 — Scan file to find the closest available α to α_query.
     * The tables contain multiple α-values; we pick the nearest one to
     * ensure smooth behaviour of φ_p.
     * --------------------------------------------------------------------- */

    double best_diff_alpha = 1e99;
    double best_alpha      = 0.0;

    double N_tmp, alpha_tmp, phi_tmp;
    while (fscanf(file, "%lf %lf %lf", &N_tmp, &alpha_tmp, &phi_tmp) == 3) {
      double diff = fabs(alpha_tmp - alpha_query);
      if (diff < best_diff_alpha) {
        best_diff_alpha = diff;
        best_alpha = alpha_tmp;
      }
    }

    rewind(file); // restart file for second scan

    /* ---------------------------------------------------------------------
     * Step 2 — Load only the rows corresponding to the selected α-slice.
     * We dynamically allocate memory and fill temporary buffers (N_buf,
     * phi_buf). They are freed immediately afterwards.
     * --------------------------------------------------------------------- */

    double *N_buf = NULL;
    double *phi_buf = NULL;
    int buf_size = 0;
    int buf_capacity = 1024;

    N_buf   = (double*) malloc(buf_capacity * sizeof(double));
    phi_buf = (double*) malloc(buf_capacity * sizeof(double));

    while (fscanf(file, "%lf %lf %lf", &N_tmp, &alpha_tmp, &phi_tmp) == 3) {
      if (alpha_tmp == best_alpha) {
        if (buf_size >= buf_capacity) {
          buf_capacity *= 2;
          N_buf   = (double*) realloc(N_buf, buf_capacity * sizeof(double));
          phi_buf = (double*) realloc(phi_buf, buf_capacity * sizeof(double));
        }
        N_buf[buf_size]   = N_tmp;
        phi_buf[buf_size] = phi_tmp;
        buf_size++;
      }
    }
    fclose(file);

    if (buf_size == 0) {
      fprintf(stderr, "No rows found for α-neighbour = %e\n", best_alpha);
      exit(1);
    }

    /* ---------------------------------------------------------------------
     * Step 3 — Find the closest N-value to N_query inside the α-slice.
     * Again, we simply pick the nearest neighbour. Interpolation is not
     * required because the tables are dense enough in N.
     * --------------------------------------------------------------------- */

    double phip_val = phi_buf[0];

    if (N_query <= N_buf[0]) {
      phip_val = phi_buf[0];
    }
    else if (N_query >= N_buf[buf_size-1]) {
      phip_val = phi_buf[buf_size-1];
    }
    else {
      for (int i = 1; i < buf_size; i++) {
        if (N_buf[i] > N_query) {
          phip_val = phi_buf[i-1];
          break;
        }
      }
    }

    phip = phip_val;

    /* Free temporary memory buffer */
    free(N_buf);
    free(phi_buf);

    /* Update cache */
    phip_cached   = phip_val;
    N_cached      = N_query;
    alpha_cached  = alpha_query;
    already_ran   = _TRUE_;
  }
}

    /* At this point phip contains the inflaton value φ_p(N, α)
     * obtained without any slow-roll approximation, fully numerical. */

    /* Coeficientes da Serie de Taylor (CLASS internal expansion) */
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
