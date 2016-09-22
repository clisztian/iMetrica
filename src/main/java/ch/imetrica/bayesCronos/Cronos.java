package ch.imetrica.bayesCronos;

import java.io.*;

public class Cronos
{

   int n_obs;
   int n_fsteps;
   int nsims;

   double[] tseries;
   double[] residuals;
   double[] fcasts; 
   double[] fmse; 
   double[] predictiveAvg;
   
   double AICC;
   int model;        //0 = ARMA, 1 = GARCH, 2 = EGARCH, 3 = SV
   int methodtype;   //1 = MLE, 2 = BAYESIAN

   //--- ARMA stuff ------
   int p,q;
   double d; 
   public double[] predictive;
   public double[] ar_params; 
   public double[] ma_params;
   public double innvar;
   int n_params; 
   int[] pindex;            //---index of sorted values
   double[][] bayes_sims;
   

   //---GARCH stuff ------
   double[] garch_ar;
   double[] garch_ma;
   double[] egarch_gs;
   double mu;
   int n_gparams; 

   //--- Stochastic Volatility stuff------
   double[] alpha;
   double[] sv_params;
   double mean, mux, phi, nu;


   //----- HEAVY Parameters -------------------------

 
   double[] w;
   double[] h_alpha; 
   double[] h_lambda;
   double[] h_beta;

   int heavy_retrack;

   HeavyModel heavy; 
 
   //---- MCMC samples for IS SV/FSV -----------

   int seed;         //random seed for Monte Carlo
   int n_samp;       //number of MCMC samples (default 500)     
   int n_iters;      //number of MCMC iterations (default 5000)
   int blocksize;    //blocksize for sample Alpha AR(1) process

   int n_factors;    //number of dynamic factors (k)
   int n_rep;        //number of series

  
   double sigma;                      //precision on alpha AR(1) process
   double[] phi_sim;    //k*n_samp simulations of phi 
   double[] sig_sim;    //k*n_samp simulations of sigma
   double[] mu_sim;     //k*n_samp simulations of mean mu
   double[] alpha_sim;  //k*n_samp simulations of volatility process alpha (length n_obs)


   double[] alpha_avg;  //avg of the alpha process (length n_obs)
   double phi_avg; 
   double mu_avg;
   double sig_avg;
  

   //--------- FOR MULTIVARIATE SVM -------------------------------------------
   double[] m_tseries;   //origanized in n_obs x n_rep (opposite of mdfa)

   //-------- prior parameters for distributions ------------------------------
   double fphi, fsig, fmu;
 
   //--------- histograms for simulated parameters of factor models, length n_fact*n_samp
   double[] fphi_sim; double[] fsig_sim;

   //-------- load matrices corresponding to factors, n_samp * (N x n_facts)
   double[] b_sim;

   //-------- averages of the histograms --------------------
 
   double[][] b_avg; 
   double[][] falpha_avg;
   double[][] factors_avg;   

   double[]   factors;
   double[]   fphi_avg; 
   double[]   fsig_avg;
   double[]   wphi_avg;
   double[]   wmu_avg;
   double[]   wsig_avg;
 
 
   public Cronos(int n, int _model, int _methodtype)
   {
      n_obs = n;
      tseries = new double[n];
      residuals = new double[n];       
   
      model = _model;
      methodtype = _methodtype;

      sv_params = new double[5];   

      //---- initialize ISV/FSV stuff----
      n_factors = 1;
      n_samp = 1000;
      n_iters = 6000;
      blocksize = 50;
      seed = 433;
      n_rep = 1;
      heavy_retrack = 0;
 
      //---------- Initialize Heavy parameters-------------------
      w = new double[3]; h_alpha = new double[3]; h_beta = new double[3]; h_lambda = new double[2];

      w[0] = .1; w[1] = .1; w[2] = .1;
      h_alpha[0] = .1; h_alpha[1] = .1; h_alpha[2] = .1; 
      h_lambda[0] = .1; h_lambda[1] = .1; 
      h_beta[0] = .8; h_beta[1] = .8; h_beta[2] = .8; 

   }

   public void setData(double[] data)
   { 
      n_obs = data.length;
      tseries = new double[n_obs];
      residuals = new double[n_obs]; 
      System.arraycopy(data, 0, this.tseries, 0, n_obs); 

   }


   public void changeBlockSize(int block)
   {blocksize = block;}


   public void setMultiData(double[] data, int nrep, int nobs)
   {
      int i,j;
      n_rep = nrep; n_obs = nobs;

      m_tseries = new double[n_obs*n_rep]; 
  
      for(i=0;i<n_obs;i++)
      {
        for(j=0;j<n_rep;j++)
        {
           m_tseries[i*n_rep + j] = data[j*n_obs + i];
        }
      }      
   }
   
   public void setSVM()
   {
     n_rep = 1;
     phi_sim = new double[n_rep*n_samp]; 
     sig_sim = new double[n_rep*n_samp];
     mu_sim = new double[n_rep*n_samp];
     alpha_sim = new double[n_rep*n_samp*n_obs];
   }


   //---------- setup multivariate FSVM arrays with given variables
   public void setMSVM()
   {

     phi_sim = new double[n_rep*n_samp]; 
     sig_sim = new double[n_rep*n_samp];
     mu_sim = new double[n_rep*n_samp];

     //--------- histograms for simulated parameters of factor models, length n_fact*n_samp
     fphi_sim = new double[n_factors*n_samp]; 
     fsig_sim = new double[n_factors*n_samp]; 
     

     //-------- histograms for alpha process length n_obs, n_fact*n_samp*n_obs
     alpha_sim = new double[n_factors*n_samp*n_obs]; 

     //-------- load matrices corresponding to factors, n_samp * (N x n_facts)
     b_sim = new double[n_rep*n_factors*n_samp]; 
     factors = new double[n_factors*n_samp*n_obs];
   }
     

   public void setDefaultParameterValues()
   {
      fphi=0.4; fsig=0.1; fmu=0.56;
      setInitialParamsISV(0.5, 0.1, 0.3);
   }


   public void setNForecastSteps(int n_t)
   {
     n_fsteps = n_t;
     fcasts = new double[n_obs+n_t];
     fmse = new double[n_t];
   }
    
   public void setNPredictiveSims(int n_f, int _n)
   {
      n_fsteps = n_f;
      nsims = _n;

      fcasts = new double[n_obs+n_fsteps];
      fmse = new double[n_fsteps];

      predictive = new double[nsims*n_fsteps]; 
      predictiveAvg = new double[n_fsteps];
   }



   public void setARMA_Params(int ap, int aq, double ad) 
   {
      p = ap; q = aq; d = ad;
      ar_params = new double[p];
      ma_params = new double[q];      
      n_params = p+q+3;
   }  

   public void setGARCH_Pararms(int ap, int aq)
   {
      p = ap; q = aq; 
      garch_ar = new double[p+1];
      garch_ma = new double[q]; 
      egarch_gs = new double[p+1];
      n_gparams = p+1+q+1; 
   }




   public void computeARIMAModel(int MLE)
   {
      int i,j; double sum = 0.0;
      residuals = new double[n_obs];
      predictive = new double[n_fsteps*nsims];
      predictiveAvg = new double[n_fsteps];
      setNForecastSteps(n_fsteps);
      
      ar_params = new double[p]; ma_params = new double[q]; 
      double[] params = new double[p+3];

      if(MLE == 1)
      {
        computeARIMA(tseries, residuals, n_obs, p, d, q, n_fsteps, nsims, params, ma_params,
                  fcasts, fmse, predictive);
      }
      else 
      {
    
         n_params = p+q+3;
         double[] bayes = new double[n_samp*n_params]; 
         double[] temp = new double[n_samp];
         pindex = new int[n_samp];

         bayes_sims = new double[n_params-2][n_samp]; 

         computeBayesianARIMA(tseries, residuals, n_obs, p, d, q, n_samp, n_iters, n_fsteps, nsims, bayes, params, ma_params,
                  fcasts, fmse, predictive);

         //----- extract sigma coefficient and sort---
         for(j=0;j<n_samp;j++)
         {temp[j] = bayes[n_params*j + 2];}

         sortsims(n_samp, temp, pindex);

         //---- now place in 
         for(j=0;j<n_samp;j++)
         {bayes_sims[0][j] = temp[j];}        

         
         for(i=1;i<n_params-2;i++)
         {       
            for(j=0;j<n_samp;j++)
            {temp[j] = bayes[n_params*j + 2 + i];}            

            sortsims(n_samp, temp, pindex);

            //---- now place in 
            for(j=0;j<n_samp;j++)
            {bayes_sims[i][j] = temp[j];}           
         }
      }

      for(i=0;i<p;i++) {ar_params[i] = params[i];}
 
      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j=0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }

      AICC = params[p];
      mu = params[p+1];
      innvar = params[p+2];

   }

   public void computeGARCHModel(int MLE)
   {
      int i,j; double sum = 0.0;
      residuals = new double[n_obs];
      predictive = new double[n_fsteps*nsims];
      predictiveAvg = new double[n_fsteps];
      setNForecastSteps(n_fsteps);
      
      //---- SPECIAL: COMPUTE h_t + h_{t+steps} --------------------
      fcasts = new double[n_obs+n_fsteps];

      garch_ar = new double[p+1]; garch_ma = new double[q];
      double[] garch = new double[p+3]; 

      if(MLE == 1)
      {
        computeGARCH(tseries, residuals, n_obs, p, q, n_fsteps, nsims, garch, garch_ma,
                  fcasts, fmse, predictive);

      }
      else 
      {
    
         n_gparams = p+1+q+1;
         double[] bayes = new double[n_samp*n_gparams]; 
         double[] temp = new double[n_samp];
         pindex = new int[n_samp];

         bayes_sims = new double[n_gparams][n_samp]; 

         computeBayesianGARCH(tseries, residuals, n_obs, p, q, n_samp, n_iters, n_fsteps, nsims, bayes, garch, garch_ma,
                  fcasts, fmse, predictive);

         //----- extract sigma coefficient and sort---
         for(j=0;j<n_samp;j++)
         {temp[j] = bayes[n_gparams*j];}

         sortsims(n_samp, temp, pindex);

         //---- now place in 
         for(j=0;j<n_samp;j++)
         {bayes_sims[0][j] = temp[j];}        

         
         for(i=1;i<n_gparams;i++)
         {       
            for(j=0;j<n_samp;j++)
            {temp[j] = bayes[n_gparams*j+i];}            

            sortsims(n_samp, temp, pindex);

            //---- now place in 
            for(j=0;j<n_samp;j++)
            {bayes_sims[i][j] = temp[j];}               
         }
      }

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }
  
      for(i=0;i<=p;i++) {garch_ar[i] = garch[i];}


      mu = garch[p+1]; 
      AICC = garch[p+2];
   }

   public void computeEGARCHModel()
   {
      int i,j; double sum;
      residuals = new double[n_obs];
      predictive = new double[n_fsteps*nsims];
      predictiveAvg = new double[n_fsteps];

      setNForecastSteps(n_fsteps);

      //---- SPECIAL: COMPUTE h_t + h_{t+steps} --------------------
      fcasts = new double[n_obs+n_fsteps];

      
      garch_ar = new double[p+1]; garch_ma = new double[q]; egarch_gs = new double[p+1];
      double[] garch = new double[p+3]; 


      computeEGARCH(tseries, residuals, n_obs, p, q, n_fsteps, nsims, garch, garch_ma, egarch_gs,
                  fcasts, fmse, predictive);

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }


      for(i=0;i<=p;i++) {garch_ar[i] = garch[i];}
      mu = garch[p+1]; AICC = garch[p+2];
   }


   public void computeSVModel()
   {
      int i,j; double sum;
      residuals = new double[n_obs];
      predictive = new double[n_fsteps*nsims]; 
      predictiveAvg = new double[n_fsteps];
      setNForecastSteps(n_fsteps);
      fcasts = new double[n_obs+n_fsteps];

      computeSVM(tseries, residuals, n_obs, n_fsteps, nsims, sv_params, fcasts, fmse, predictive);

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }

      //take exponential and divide by 2
      for(i=0;i<fcasts.length;i++)
      {
        fcasts[i] = Math.exp(fcasts[i]/2.0);
      }


      mean = sv_params[0];
      mux = sv_params[1];
      phi = sv_params[2];
      nu = sv_params[3];   
      AICC = sv_params[4];  
   }



   void setTrackReparameter(int h) {heavy_retrack = h;}


   void setHeavyParameters(double[] _w, double[] alpha, double[] beta, double[] lambda)
   {

      int i;
      w = new double[3]; 
      h_alpha = new double[3]; 
      h_lambda = new double[2];
      h_beta = new double[3]; 

      for(i=0;i<3;i++)
      {
        w[i] = _w[i]; h_alpha[i] = alpha[i]; h_beta[i] = beta[i];
      }
      h_lambda[0] = lambda[0]; h_lambda[1] = lambda[1];
   }

   //----------- COMPUTE HEAVY model -----------------------------------
   void computeHeavyModel()
   {
     int k,i;
     heavy = new HeavyModel(n_obs, seed); 
     heavy.setRetrackParameter(heavy_retrack);
     heavy.setForecastDimensions(nsims, n_fsteps);
     

     //--- set forecast dimensions -------------------       
     if(n_rep == 2)
     {
       //--- set (simulation or initialization) parameter values -------------
       heavy.setParameterValues(w[0], w[1], h_alpha[0], h_alpha[1], 
                                h_lambda[0], h_beta[0], h_beta[1]);
     } 
     else if(n_rep == 3)
     {
       heavy.setParameterValuesSV(w[0], w[1], w[2], h_alpha[0], h_alpha[1], h_alpha[2], h_lambda[0], 
                                  h_beta[0], h_beta[1], h_beta[2], h_lambda[1]);  
     }

     //---- now estimate ----------------------
     heavy.setData(n_obs, n_rep, m_tseries);
     heavy.estimateHeavyModel(); 
     heavy.parameterTransform();

     //---- set forecasts and distributions ---------------

     //System.out.println("nsims fsteps = " + nsims + "  " + n_fsteps);
 
     for(k=0;k<nsims;k++)
     {
      for(i=0;i<n_fsteps;i++)
      {
       predictive[n_fsteps*k+i] = heavy.r_forecasts[n_fsteps*k + i];
      }     
     }


     residuals = new double[n_obs];
     falpha_avg = new double[n_rep][n_obs];    //h,mu
     factors_avg = new double[n_rep][n_obs];   //distributions eta, zeta


     fcasts = new double[n_obs+n_fsteps];
     for(i=0;i<n_obs;i++)
     {
        fcasts[i] = heavy.ht[i]; 
        residuals[i] = heavy.llht[i];
        //System.out.println(residuals[i+1]);
     }
     residuals[0] = 0;
     for(i=0;i<n_fsteps;i++)
     {fcasts[n_obs+i] = heavy.h_forecasts[i];}

     for(k=0;k<n_rep;k++)
     {
       for(i=0;i<n_obs;i++)
       {
          falpha_avg[k][i] = heavy.ht[k*n_obs+i]; 
          factors_avg[k][i] = heavy.F_dist[k*n_obs+i];
       }
     }
        
   }



   //------- Functions for ISV/FSV ---------------------------------------------------
   
   public void setNFactors(int k) 
   {
      if(k <= n_rep)
      {n_factors = k;}
      else
      {System.out.println("Number of factors must be less than number of series");}
   }


   //---------- Forecast ARIMA using current parameter selections and data--------
   public void forecastARIMA()
   {
      int i,j; double sum = 0.0;
      predictiveAvg = new double[n_fsteps];            
      predictive = new double[n_fsteps*nsims];

      computeARIMAForecast(tseries, n_obs, p, d, q, n_fsteps,nsims, ar_params, ma_params, innvar, mu,
			   fcasts, fmse, predictive);

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }

   }

   //---------- Forecast ARIMA using current parameter selections and data--------
   public void forecastGARCH()
   {
      int i,j; double sum = 0.0;
      predictiveAvg = new double[n_fsteps];           
      predictive = new double[n_fsteps*nsims];

      fcasts = new double[n_obs+n_fsteps];

      computeGARCHForecast(tseries, n_obs, p, q, n_fsteps, nsims, garch_ar, garch_ma, mu,
                  fcasts, fmse, predictive);

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }

   }
 
   public void forecastEGARCH()
   {
      int i,j; double sum = 0.0;
      predictiveAvg = new double[n_fsteps];      
      predictive = new double[n_fsteps*nsims];
 
      fcasts = new double[n_obs+n_fsteps];

      computeEGARCHForecast(tseries, n_obs, p, q, n_fsteps, nsims, garch_ar, garch_ma, egarch_gs, mu,
                  fcasts, fmse, predictive);

      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }


   }

   public void setNIterations(int nits, int nsamps)
   {
     n_iters = nits; 
     n_samp = nsamps; n_iters = n_samp + 3000;
     
     //adjustSampleStorage();
   }

   public void adjustSampleStorage()
   {
     n_factors = 1;
     phi_sim = new double[n_factors*n_samp];    //k*n_samp simulations of phi 
     sig_sim = new double[n_factors*n_samp];    //k*n_samp simulations of sigma
     mu_sim = new double[n_factors*n_samp];     //k*n_samp simulations of mean mu
     alpha_sim = new double[n_factors*n_samp*n_obs];  //k*n_samp simulations of volatility process alpha (length n_obs)
     pindex = new int[n_samp];

     for(int i=0;i<n_samp;i++) {pindex[i] = i;}

   }


   public void setInitialParamsISV(double _phi, double _sigma, double _mu)
   {phi = _phi; sigma = _sigma; mu = _mu;}


   //---------------------------------------------------------
   //  Assumes shit has been sorted, namely the phi_sim
   //  Uses phi_index to extract other sims
   //---------------------------------------------------------  

   public void computeBayesianAvgARIMA(int start, int end)
   {
     int i, j;
     ar_params = new double[p];
     ma_params = new double[q]; 
     innvar = 0.0;

     if(end > n_samp) {end = n_samp;}
     if(start < 0) {start = 0;}

     System.out.println(start + "  " + end); 
     
     int total = end - start;
      
     for(i = start; i < end; i++)
     {
        innvar = innvar + bayes_sims[0][i]/total;
  
        for(j=0;j<p;j++)
        {ar_params[j] = ar_params[j] + bayes_sims[1+j][i]/total;}

        for(j=0;j<q;j++)
        {ma_params[j] = ma_params[j] + bayes_sims[1+p+j][i]/total;}  
        //System.out.println("indx = " + indx); 
     }
      
     System.out.println(ar_params[0] + "  " + ma_params[0]); 
   }


   public void computeBayesianAvgGARCH(int start, int end)
   {
     int i, j;
     garch_ar = new double[p+1];
     garch_ma = new double[q]; 
     mu = 0.0;

     if(end > n_samp) {end = n_samp;}
     if(start < 0) {start = 0;}

     int total = end - start;
      
     System.out.println(start + "  " + end); 

     for(i = start; i < end; i++)
     {
        mu = mu + bayes_sims[0][i]/total;
  
        for(j=0;j<=p;j++)
        {garch_ar[j] = garch_ar[j] + bayes_sims[1+j][i]/total;}

        for(j=0;j<q;j++)
        {garch_ma[j] = garch_ma[j] + bayes_sims[1+p+1+j][i]/total;}  
     }

     //System.out.println(garch_ar[0] + "  " + garch_ma[0]); 

   } 



   public void computeBayesianAvgISV(int start, int end)
   {

      int indx,t,i; 
      double sumphi = 0.0; 
      double summu = 0.0; 
      double sumsig = 0.0;

      alpha_avg = new double[n_obs];
     
      if(end > n_samp) {end = n_samp;}
      if(start < 0) {start = 0;}

      int total = end - start;
      
      for(i = start; i < end; i++)
      {
         indx = pindex[i];            //get corresponding index
         sumphi = sumphi + phi_sim[i];
    
         summu = summu + mu_sim[i];
         sumsig = sumsig + sig_sim[i];

         for(t=0; t<n_obs; t++)
         {
            alpha_avg[t] = alpha_avg[t] + Math.exp(alpha_sim[n_obs*indx + t]/2.0)/total;
         }
      }

      phi_avg = sumphi/total;
      mu_avg = summu/total;
      sig_avg = sumsig/total;

      alpha = new double[n_obs]; 
      phi = phi_avg; mean = mu_avg; mux = sig_avg;    
      System.arraycopy(alpha_avg, 0, alpha, 0, n_obs);
          
   } 


   //---------------------------------------------------------
   //  Assumes shit has been sorted, namely the phi_sim
   //  Uses phi_index to extract other sims
   //---------------------------------------------------------  

   public void computeBayesianAvgFISV(int start, int end)
   {

      int indx,t,i,k,r; 
      double sumphi = 0.0; 
      double summu = 0.0; 
      double sumsig = 0.0;
      int dim = n_samp*n_obs;

      falpha_avg = new double[n_factors][n_obs];
      fphi_avg = new double[n_factors];     
      fsig_avg = new double[n_factors]; 

      wphi_avg = new double[n_rep];     
      wsig_avg = new double[n_rep]; 
      wmu_avg = new double[n_rep]; 
      b_avg = new double[n_rep][n_factors];
      factors_avg = new double[n_factors][n_obs]; 
 
      if(end > n_samp) {end = n_samp;}
      if(start < 0) {start = 0;}

      int total = end - start;
  
      //System.out.println("start, end = " + start + "  " + end); 

      //----------- Do the factors first ------------------------------
      for(k = 0; k < n_factors; k++)
      {
        sumphi = 0.0; summu = 0.0; sumsig = 0.0;

        for(i = start; i < end; i++)
        {

         indx = pindex[i];            //get corresponding index
         //System.out.println("indx = " + indx);
         sumphi = sumphi + fphi_sim[n_samp*k+indx];  
         sumsig = sumsig + fsig_sim[n_samp*k+indx];

         for(t=0; t<n_obs; t++)
         {
            falpha_avg[k][t] = falpha_avg[k][t] + alpha_sim[dim*k + n_obs*indx + t]/total;
            factors_avg[k][t] = factors_avg[k][t] + factors[dim*k + n_obs*indx + t]/total;
         }

        }

        /* for(t=0; t<n_obs; t++)
         {
            falpha_avg[k][t] =  alpha_sim[dim*k + n_obs*pindex[start] + t];
            factors_avg[k][t] =   factors[dim*k + n_obs*pindex[start] + t];
         }
        */

        fphi_avg[k] = sumphi/total;
        fsig_avg[k] = sumsig/total;
      } 

  
      //----------- Now do the idiosyncratic residuals first ------------------------------
      for(k = 0; k < n_rep; k++)
      {
        sumphi = 0.0; summu = 0.0; sumsig = 0.0;

        for(i = start; i < end; i++)
        {

         indx = pindex[i];            //get corresponding index

         sumphi = sumphi + phi_sim[n_samp*k+indx];  
         summu = summu + mu_sim[n_samp*k+indx];
         sumsig = sumsig + sig_sim[n_samp*k+indx];

         for(r=0;r<n_factors;r++)
         {
          b_avg[k][r] = b_avg[k][r] + b_sim[n_rep*n_factors*indx + n_rep*r + k]/total;  
         }
        }

        wphi_avg[k] = sumphi/total;
        wmu_avg[k] = summu/total;
        wsig_avg[k] = sumsig/total;
      } 
 
   }





   public void computeBayesianISVModel()
   {

     // first compute MLE to get residuals 
     computeSVModel();
     pindex = new int[n_samp]; 
     computeBayesianISV(tseries, n_iters, n_obs, blocksize, n_samp, 
                        phi, sigma, mu, seed, phi_sim, sig_sim, mu_sim, alpha_sim);

    
     sortsims(n_samp, phi_sim, pindex);
     sortsims(n_samp, sig_sim, pindex);
     sortsims(n_samp, mu_sim, pindex);     

 
     AICC = 0.0;

   }



   public void computeFSVModel()
   {
         int i;
         pindex = new int[n_samp];

         //System.out.println("length = " + alpha_sim.length + ", length_series = " + m_tseries.length + 
         //                  ", n_factors = " + n_factors + ", n_rep = " + n_rep + ", b_sim_length = " + b_sim.length);

         computeBayesianFSV(m_tseries, n_iters, n_obs, blocksize, n_factors, n_rep, n_samp, phi, sigma, mu, 
                              fphi, fsig, seed, phi_sim, sig_sim, mu_sim, 
                              fphi_sim, fsig_sim, alpha_sim, b_sim, factors);

         for(i=0;i<alpha_sim.length;i++)
         {alpha_sim[i] = Math.exp(alpha_sim[i]/2.0);}


         //sort first fphi
         double[] temp = new double[n_samp];
         for(i=0;i<n_samp;i++) {temp[i] = fphi_sim[i];} 

         //--- sort according to factor f_phi 
         sortsims(n_samp, temp, pindex);

         //--- compute the Bayesian means to get mean values----------
         computeBayesianAvgFISV(0, n_samp);

         AICC = 0.0;
   }



   //------ Forecastes for the Bayesian model, must need averages ---- 

   public void computeBayesianForecastISVModel()
   {
      int i,j; double sum = 0.0;

      predictive = new double[n_fsteps*nsims];
      predictiveAvg = new double[n_fsteps];
      computeBayesianForecastISV(tseries, residuals, n_obs, n_fsteps, nsims, phi_avg, sig_avg, mu_avg, predictive);     
      
      for(i = 0; i < n_fsteps; i++)
      { 
       sum = 0.0;
       for(j = 0; j < nsims; j++)
       {
         sum = sum + predictive[j*n_fsteps + i];
       }
       predictiveAvg[i] = sum/nsims;
      }


   }







   public native void computeARIMA(double[] data, double[] r, int nobs, int _p, double d, int q, 
                                  int nsteps, int nsims, double[] ar_params, double[] ma_params,
                                  double[] forecasts, double[] fmse, double[] predict);


   public native void computeBayesianARIMA(double[] data, double[] r, int nobs, int _p, double d, int q, int _nsamps, 
                                  int nsteps, int _niters, int nsims, double[] bayes, double[] ar_params, double[] ma_params,
                                  double[] forecasts, double[] fmse, double[] predict);


   public native void computeARIMAForecast(double[] data, int _nobs, int _p, double _d, int _q, int _nsteps, int _nsims,
					       double[] ar_params, double[] ma_params, double _sigma, double _mean, 
					       double[] forecasts, double[] fmse, double[] predictive);



   public native void computeGARCH(double[] data, double[] r, int nobs, int _p, int q, int nsteps, 
                                   int nsims, double[] ar_params, double[] ma_params,
                                   double[] forecasts, double[] fmse, double[] predict);


   public native void computeBayesianGARCH(double[] data, double[] r, int nobs, int _p, int q, int _nsamps, int _niters, 
                                  int nsteps, int nsims, double[] bayes, double[] ar_params, double[] ma_params,
                                  double[] forecasts, double[] fmse, double[] predict);


   public native void computeEGARCH(double[] data, double[] r, int nobs, int _p, int q, int nsteps, 
                                   int nsims, double[] ar_params, double[] ma_params, double[] gs,
                                   double[] forecasts, double[] fmse, double[] predict);


   public native void computeGARCHForecast(double[] data, int _nobs, int _p, int _q, int _nsteps, int _nsims,
					       double[] ar_params, double[] ma_params, double _mean, 
					       double[] forecasts, double[] fmse, double[] predictive);

   public native void computeEGARCHForecast(double[] data, int _nobs, int _p, int _q, int _nsteps, int _nsims,
					       double[] ar_params, double[] ma_params, double[] gs_params, double _mean, 
					       double[] forecasts, double[] fmse, double[] predictive);					       

   public native void computeSVM(double[] data, double[] r, int nobs, int nsteps, 
                                   int nsims, double[] sv_params, double[] forecasts, double[] fmse, double[] predict);


   public native void computeBayesianISV(double[] data, int _iters, int _n, int _m, int _n_samp, 
                                         double _phi, double _sig, double _mu, int _randnum, 
                                         double[] phi_sim, double[] sig_sim, double[] mu_sim, double[] alpha_sim);

   public native void computeBayesianFSV(double[] data, int _iters, int _n, int _m, int _k, int n_rep, int _n_samp, double _phi, double _sig, 
						double _mu, double _fphi, double _fsig, int _randnum, 
						double[] phi_sim, double[] sig_sim, double[] mu_sim, 
						double[] fphi_sim, double[] fsig_sim, double[] falpha_sim, double[] b_sim, double[] factors);
  




   public native void computeBayesianForecastISV(double[] data, double[] r, int nobs, int nsteps, int nsims, double _phi, double _sig, double _mu, double[] predict);


   public native double[] sarima(int N, int burn, double[] params, int[] dim, int n_params, int S, int seed, int rand);


   public native void plotHistogram(double[] d, int n, int ints);

   public static native void sortsims(int n, double[] d, int[] indx);


   static {System.loadLibrary("sv");}
   static {System.loadLibrary("ms");} 
   static {System.loadLibrary("simSarima");}

    public void readData(File file)
    {
          
       String strline; Double D; 
       String[] tokens; String delims = "[ ]+";
       int n_toks;
       double[] values = new double[3000];      
       double val = 0;
       int i = 0;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           D = new Double(tokens[n_toks-1]); //take only the value, no dates
           val = D.doubleValue();

           if(i < 3000) {values[i] = val; i++;}
           else 
           {System.out.println("Maximum times series length is 3000"); break;} 
           
         } 
         
         n_obs = i;
         this.tseries = new double[n_obs]; residuals = new double[n_obs]; 
         System.arraycopy(values, 0, this.tseries, 0, n_obs);    
         
         din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }   
    

}
  