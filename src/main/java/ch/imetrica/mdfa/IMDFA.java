package ch.imetrica.mdfa;

/*!
Copyright (C) 2016 Christian D. Blakely
This file is part of iMetrica, a free-software/open-source application
for interactive graphical econometric analysis - http://imetricablog.com/

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

import java.io.*;


public class IMDFA
{
  
    //---- Data and filter controls -------------
    public int n_obs;
	public int n_rep;
	int K;
	public int L;
	public int Lag;
	int S;
    public int i1;
	public int i2;
	int dd;
	int DD;
	int K1;
	int output;
	int flength;  
    boolean mdfa, mdfa_iter, dfa_iter;
	public boolean iter;
	boolean imdfaReg; 
    public double[] Gamma; 
    //---- extra Filter controls ----------------
    public double lambda;
	public double expweight;
	public double cutoff;
	public double cutoff0;
	double lambda_3;   
    public double smooth;
	public double decay;
	public double cross;
	public double decay2;
	double diff_onestep;
    public double[] w_const;
    double criteria, degrees;
    double tdelay;
    double shift; 
    double diff_stop, diff_band;
    int update_start;
    boolean save_signal = true;
    double max_amp;
    boolean useSD = false;
    boolean b_univ = false;
    boolean H0set = false;
    boolean useH0 = false;
    //------ New idea for periodogram weighting function------------------------------------
    double[] fullRangeData;
    double[] specDensPeriodogram;
    
    //----- generalized spectral density 
    //if true, computes standard DFT, if false uses spectral density function hybrid_specDens
    public boolean estimate_specDens = true; 
    public double[] hybrid_specDens; 
    //---------------------------------------------------------------------------------------   

    
    private int n_comps = 0;

    //--- Min mdfa -----
    double MDFAmin;   
 
    int samp;  //periodogram sampling    

    //---- time series and filter data, n_rep x n_obs or K
    public double[] tseries; 
    public double[] amp_filter; 
    public double[] time_delay;
    public double[] gamma_hat;
    public double[] period_xf;  //--periodograms of x and xf
    public double[] period_hat; //--periodograms of supporting data
    public double[] b;
    public double[] b_old;
    public double[] b_update;
    public double[] h0b0;
    public double[] comb;
    public double[] b_price_filter;         //price filter
    public boolean price_filter = false;
    public int L_price; //length of price filter (must be less than L)
    
    //---- Filtered output and original series 
    public double[] xf; 
    public double[] x;
    public double[] Q_cdev;  // ---- computed in advance the central deviations
    public double[] xf_orig;

    public double f_start; 
    public int n_div; 
    public int n_imfs;
    public double[][] amMap; public double[][] phaseMap; public double[][] imfs;
    public double[][] fmMap; public double[][] ifreqMap; 
    public double[] trend_cycle;
    public double[] spec_dens;
    public int spec; 
    

    
 
    

    //------------------------- Spectral density stuff ---------------------------
    public int arma = 0;
    public int ar_p, ma_q;
    public double[] ar_params;
    public double[] ma_params;
    public double innvar;
    public double[] arg;
    public double[] mod;

    
    //----------------------- Trading parameters ----------------------------------
     public int ahead = 0;           //how many days to look ahead before buying
     public double cost = 0.0;       //the cost to buy/sell
     public int short_sell = 0;      //do we short sell, 0 = no, 1 = yes
     public double risk_free = 0.0;       //a risk free rate 
     public int[] trades;            //trade statistics (succ trades, total trades
     public int[] transactions;      //zero crossings/transactions
     public double[] account; 
     public int trade_obs;           //number of observations - probably different from n_obs    
     public double mse_update;
    //---------- Plotting options --------------
      
     public int update_i1=0;
     public int update_i2=0;
     public int L_update; //--- the update filter length
     public int N_update; //--- the total length of new data (N = number of out-of_sample points + L)
     public double[] update_signal;
     public int K_update; //--- the update filter K
     public double[] amp_update;
     public int nrep_update;
    //------------------------------------------

    public IMDFA(int nobs, int nrep, int _L, int _lag, double _lambda, double _expweight, double _cut, int _i1, int _i2)
    {
        int i; n_obs = nobs; L = _L; i1 = _i1; i2 = _i2; n_rep = nrep;
        lambda = _lambda; expweight = _expweight; cutoff = _cut; Lag = _lag; cutoff0 = 0.0;
        mdfa = true; dd = 0; S=12; iter = false;  smooth=0.0; decay=0.0; cross=0.0; imdfaReg=false; decay2 = 0.0;
        //---- Use mdfa of dfa
        if(n_rep <= 1) mdfa = false;  spec = 0;

        diff_stop = 0.0; diff_band = 0.0;
        shift = 0.0;
        K = (int)n_obs/2; K1 = K+1; // --- frequency         
        flength = n_obs-L+1; //-----  filtered length
        samp = K;
        MDFAmin = 0;
        w_const = new double[n_rep]; w_const[1] = 1.0;
        Gamma = new double[K1]; 
        n_div = 15; f_start = .40;
        trend_cycle = new double[flength];
        arma = 0; ar_p=0; ma_q=0;

        diff_onestep = 0;
        trades = new int[2];
        ar_params = new double[1];
        ma_params = new double[1];

        spec_dens = new double[2*K1];
        for(i = 0; i < K1; i++)
        {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}

        //if(mdfa) 
        //{output = 2*flength + 2*(n_rep+1)*K1 + (n_rep-1)*K1 + (n_rep)*L + 1;} 
        //else  //--- dfa only
        // {output = 2*flength + K1 + K1 + L + 1;} 
        
        //----- set original comb for filter ------
        comb = new double[n_rep];
        for(i=0;i<n_rep;i++)
        {comb[i] = 1.0;}
        comb[0] = 0.0;
         
        computeQcdev();
        set_output_len();
    }


   //---------- set ZPC spect
   public void setSpecDensity(double[] _spec_dens, int _spec)
   {
      spec = _spec;
      if(_spec == 1) 
      {
        spec_dens = new double[2*K1]; 
        System.arraycopy(_spec_dens, 0, spec_dens, 0, 2*K1);
      } 
      else 
      {
        spec_dens = new double[2*K1];
        for(int i = 0; i < K1; i++)
        {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}
      }
   }

   
   public void setSpecDensityParams(double[] ar, double[] ma, double _innvar)
   {
     
      innvar = _innvar;
      ar_p = ar.length; ma_q = ma.length;
      ar_params = new double[ar.length];
      ma_params = new double[ma.length];

      System.arraycopy(ar, 0, ar_params, 0, ar_p);
      System.arraycopy(ma, 0, ma_params, 0, ma_q);

      arma = 1;
   }


   public void setPriceFilter(double[] pr)
   {
      
      L_price = pr.length;
      b_price_filter = new double[pr.length];
   
      System.arraycopy(pr,0,b_price_filter,0,L_price);   
      price_filter = true;
   }
   
   public void useSpectralDensity(boolean c) {if(c) {arma = 1;} else {arma = 0; ar_p=0; ma_q=0;}}
   public void useZPCWeight(int c) 
   {
     spec = c; int i;   
     if(spec == 0)
     {
      spec_dens = new double[2*K1];
      for(i = 0; i < K1; i++)
      {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}
     }
   }
   
   public void setWConstVals(double[] w)
   {w_const = new double[n_rep]; w_const[0] = 0.0; for(int i=1;i<n_rep;i++) {w_const[i] = Gamma[0]*w[i-1];}}
  
   //-------------- Change the parameters of filter data --------------
   public void set_L(int _L) {this.L = _L; this.resetL();} 
   public void set_nobs(int _n) {this.n_obs = _n; this.K = (int)this.n_obs/2; this.K1 = K+1; this.set_output_len();} // --- frequency 
   public void set_nreps(int _n) {this.n_rep = _n; this.set_output_len(); computeQcdev();}
   public void set_lag(int _lag) {this.Lag = _lag;}   
   public void set_S(int _S) {this.S = _S;}
   //------------------------------------------------------------------- 

   public void setRegularization(double _smooth, double _decay, double _decay2, double _cross)
   {this.smooth = _smooth; this.decay = _decay; this.cross = _cross; this.decay2 = _decay2;}

   public void setReg(boolean t) {imdfaReg = t;}
   public void set_dfaIter(boolean m) {iter = m;}
   public void useH0Set(boolean b) {useH0 = b;}
   
   //-------------- Change the parameters of filter output --------------
   public void set_dd(int _dd) {this.dd = _dd;}
   public void set_DD(int _dd) {this.DD = _dd;}
   public void set_shift(double s) {this.shift = s;}
   public void set_bconstraints(int _i1, int _i2) {this.i1 = _i1; this.i2 = _i2; computeQcdev();}
   public void set_lambda(double l) {this.lambda = l;}
   public void set_exp(double l) {this.expweight = l;}
   public void set_cutoff(double l) {this.cutoff = l;}
   public void set_mdfa(boolean m) {this.mdfa = m;}
   public void set_Samps(int _s) {this.samp = _s;}
   public void set_cutoff0(double l) {this.cutoff0 = l;}
   public void setLambda3(double _l) {this.lambda_3 = _l;}
   public void setOnestep(double _l) {this.lambda_3 = -_l;}
   public void setDiffOnestep(double _l) {this.diff_onestep = _l;}
   
   public void setCombValue(int i, double c)
   {
     if(i < n_rep)
     {comb[i] = c;}
   }

   //------------------------------------------------------------------
   
   //------------------ Set time series ------------------------------
   public void set_tseries(double[] _tseries, int N, int R)
   {
      this.n_obs = N; this.n_rep = R; K = N/2; K1 = K+1; flength = n_obs - L + 1;

      if(this.n_rep <= 1) mdfa = false;  
      else {mdfa = true;}
 
      this.set_output_len();

      int le = _tseries.length;
      this.tseries = new double[le];
      System.arraycopy(_tseries, 0, this.tseries, 0, le);           
      
      if(spec == 1 && fullRangeData != null)
      {computeSpecDensPeriodogram();}
      //else if(fullRangeData == null) {System.out.println("full range data null");}
   }
   
   
   //-------------- Set the tseries data in the form for mdfa Stages-------------------
   //   target = series, explanatory = filtered series 
   //----------------------------------------------------------------------------------
   
   public void set_tseries_stage(double[] signal, int N, double[] series)
   { 
      int i;
      this.n_obs = N; this.n_rep = 2; K = N/2; K1 = K+1; flength = this.n_obs - L + 1;
 
      this.set_output_len();

      tseries = new double[N*2];      
      for(i=0;i<N;i++)
      {
        tseries[N-1-i] = series[series.length-1-i];
        tseries[N+N-1-i] = signal[signal.length - 1 - i];
      }
            
      if(spec == 1 && fullRangeData != null)
      {computeSpecDensPeriodogram();}
      //else if(fullRangeData == null) {System.out.println("full range data null");}
   }   
   
   
   public void set_tseries_stageMult(double[] signal, int N, double[] series, int r)
   { 
      int i,k;
      this.n_obs = N; this.n_rep = r + 1; 
      K = N/2; K1 = K+1; flength = this.n_obs - L + 1;
 
      this.set_output_len();
      tseries = new double[N*this.n_rep];      
      for(i=0;i<N;i++)
      {
        tseries[N-1-i] = series[series.length-1-i];
      }
      
      for(k = 0; k < r; k++) 
      {
       for(i=0;i<N;i++)
       {
        tseries[N*(k+2)-1-i] = signal[N*(k+1) - 1 - i];
       }
      }
            
      if(spec == 1 && fullRangeData != null)
      {computeSpecDensPeriodogram();}
      //else if(fullRangeData == null) {System.out.println("full range data null");}
   }  
   
   
   
   
   //------------------ Set symmetric target filter ------------------------------
   public void set_Gamma(double[] _gamma) //---- target symmetric filter ------
   {     
      K1 = _gamma.length;
      this.Gamma = new double[K1];
      System.arraycopy(_gamma, 0, Gamma, 0,K1);
   }

   public void resetL()
   {
        flength = n_obs - L + 1;
        if(mdfa) 
        {
          output = 2*flength + 2*(n_rep+1)*K1 + (n_rep-1)*K1 + (n_rep)*L;
          xf = new double[flength]; x = new double[flength];
          b = new double[(n_rep)*L];
          computeQcdev();
        } 
        else  //--- dfa only
        {
          output = 2*flength + K1 + K1 + L;     
          xf = new double[flength]; 
          x = new double[flength];
          b = new double[L];      
        }     
   }

   //----------------- Update lengths of arrays ---------------------
   public void reset_lengths() 
   {
        int i;
        if(mdfa) 
        {
          amp_filter = new double[(n_rep+1)*K1]; 
          time_delay = new double[(n_rep+1)*K1];
          gamma_hat = new double[(n_rep+1)*K1];
          b = new double[n_rep*L];
          period_hat = new double[(n_rep-1)*K1];
          for(i=0;i<n_rep;i++)
          {comb[i] = 1.0;}
          comb[0] = 0.0;        
        }
        else  //--- dfa only
        {          
          output = 2*flength + K1 + K1 + L;   //System.out.println("gets here = " + output);
          amp_filter = new double[K1]; time_delay = new double[K1]; gamma_hat = new double[K1]; 
          b = new double[L]; 
          //xf = new double[flength]; x = new double[flength];     
        }     
   }


   //----------------- Update lengths of arrays ---------------------
   public void set_output_len() 
   {
        int i;  double[] temp;
        //System.out.println("set_nrep = " + n_rep);
        flength = n_obs-L+1; //-----  filtered length
        if(n_rep == 1){mdfa = false;}
        else{mdfa = true; }

        if(mdfa) 
        {         
          output = 2*flength + 2*(n_rep+1)*K1 + (n_rep-1)*K1 + (n_rep)*L + 1; 
          tseries = new double[n_rep*n_obs];
          amp_filter = new double[(n_rep+1)*K1]; time_delay = new double[(n_rep+1)*K1];
          gamma_hat = new double[(n_rep+1)*K1];
          xf = new double[flength]; 
          x = new double[flength]; 
          b = new double[n_rep*L];
          period_hat = new double[(n_rep-1)*K1];
 
          temp = new double[(int)Math.max(w_const.length,n_rep)]; 
          for(i=0;i<w_const.length;i++) {temp[i] = Gamma[0]*w_const[i];} 
          w_const = new double[n_rep];
          for(i=0;i<n_rep;i++) {w_const[i] = Gamma[0]*temp[i];}
          
          spec_dens = new double[2*K1];
          for(i = 0; i < K1; i++)
          {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}
                   
        }
        else  //--- dfa only
        {
          w_const = new double[2]; w_const[0] = 0.0; w_const[1] = Gamma[0];
          output = 2*flength + K1 + K1 + L; tseries = new double[n_obs];  //System.out.println("gets here = " + output);
          amp_filter = new double[K1]; time_delay = new double[K1]; gamma_hat = new double[K1];
          xf = new double[flength]; x = new double[flength]; b = new double[L];  spec_dens = new double[2*K1];  
          spec_dens = new double[2*K1];
          for(i = 0; i < K1; i++)
          {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}
          
        }     

   }

   /*------------------------------------------------------------------------------
     Computes the filter input of ZPC data. Assumes data is now ZPC filtered output of tseries
   --------------------------------------------------------------------------------*/
   public void computeFilterFromZPC(double[] zpc, int zpc_obs)
   {      
       int i,k,lag2,lag3,lag4; tdelay = 0.0;
       flength = zpc_obs-L+1;
       
        double[] out = new double[output];        
        //System.out.println("length = " + tseries.length + ", nrep = " + n_rep + ", nobs = "+n_obs);          
        //for(i=0;i<n_obs;i++) {System.out.println(tseries[i]);}

        if(mdfa)
        {
               
         out = IMDFAreg(zpc, zpc_obs, L, n_rep, dd, DD, Gamma, shift, lambda, expweight, lambda_3, cutoff0,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, cross, Q_cdev, iter);
         //n_comps++;
         //System.out.println(n_comps);
         xf = new double[flength]; 
         x = new double[flength]; 

         for(i=0;i<flength;i++)  //---- get filtered and original series
         {
           xf[i] = out[i];
           x[i] = out[i+flength]; 
         }
         lag2 = 2*flength; lag3 = 2*flength+(n_rep+1)*K1; 

         for(i=0; i <= n_rep; i++)  //---- get amplitude/phase of filter
         {   
           for(k=0;k<=K;k++)
           {
             amp_filter[K1*i+k] = out[lag2 + K1*i + k];
             time_delay[K1*i+k] = -out[lag3 + K1*i + k];
           }
         }
         
         lag4 = lag2 + 2*(n_rep+1)*K1;
         for(i=0; i <=n_rep; i++)  //---- get frf
         {   
           gamma_hat[K1*i] = amp_filter[K1*i];
           for(k=1;k<=K;k++)
           {gamma_hat[K1*i+k] = amp_filter[K1*i+k]*Math.cos(time_delay[K*i+k]*(Math.PI*k/K));}

         }       
         //----------- frf ------------------
         for(i=0; i < n_rep; i++)  
         {
           for(k=0;k<L;k++)
           {b[L*i + k] = out[lag4 + L*i + k];}
         }
       
         criteria = out[out.length-3];
         degrees = out[out.length-2];      
         MDFAmin = out[out.length-1];
  


       }
       else
       { 
         if(spec == 0) {spec_dens = new double[K1];}
         //System.out.println(cutoff0 + "  " + cutoff);
         //if(iter) {out = IDFAiter(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2);}
         //else {out = IDFAreg(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2, smooth, decay, 1);}          
         out = IDFAreg(zpc, zpc_obs, L, dd, DD, Gamma, shift, lambda, expweight, lambda_3, 
                        cutoff0, cutoff, Lag, i1, i2, smooth, decay,iter,spec_dens,spec);       
         lag2 = 2*flength;
         lag3 = 2*flength+K1;
         lag4 = lag3 + K1;

         xf = new double[flength]; 
         x = new double[flength]; 

        for(i=0;i<flength;i++)  //---- get filtered and original series
        {
           xf[i] = out[i];
           x[i] = out[i+flength]; 
        }
         
        for(k=0;k<=K;k++)
        {amp_filter[k] = out[lag2 + k]; time_delay[k] = -out[lag3 + k];}
                 
        gamma_hat[0] = amp_filter[0];
        for(k=1;k<=K;k++)
        {gamma_hat[k] = amp_filter[k]*Math.cos(time_delay[k]*Math.PI*k/K); tdelay = tdelay + time_delay[k];}
 
        for(k=0;k<L;k++) 
        {b[k] = out[lag4 + k];}

        criteria = out[out.length-3];
        degrees = out[out.length-2];      
        MDFAmin = out[out.length-1];

       }

       computeSampleIns();  
      
   }



   /*------------------------------------------------------------
     - input - We have an asset with n_obs trade observations -- double asset[n_obs]
     - input - We have a signal giving zero crossings on when to buy -- double signal[n_obs]
     - input - transaction costs 
     - input - risk_free rate

     - output - trades, account
   -------------------------------------------------------------*/
   public void tradingFunction(double[] asset, double[] signal, double cost, double risk_free)
   {
       trade_obs = signal.length;
 
       ahead = 0; 
       account = new double[trade_obs];
       double[] _trades = new double[2];
       trades = new int[2];
 
       if(asset.length == trade_obs)
       {
         tradingDiffSingle(asset, signal, account, trade_obs, ahead, cost, short_sell, risk_free, _trades);
         trades[0] = (int)_trades[0]; trades[1] = (int)_trades[1];
       }
       else
       {System.out.println("Number of observations needs to agree");} 

   }


  public void setSpecDensityEstimate(double[] _mod, double[] _arg, int _K, int _n_rep)
  {

    if(_n_rep == n_rep && _K == K)
    {
 
      //this.K = _K;
      //this.K1 = K+1; 
      //this.set_output_len();

      mod = new double[_mod.length];
      arg = new double[_mod.length];

      System.arraycopy(_mod,0,mod,0,mod.length);
      for(int i = 0; i < mod.length;i++) {mod[i] = Math.sqrt(mod[i]);} 
      System.arraycopy(_arg,0,arg,0,arg.length);

      useSpectralDensityEstimate(true);
    }
  }

  public void useSpectralDensityEstimate(boolean tr)
  {
    useSD = tr; 
    computeFilterGeneral(true,false);
  }

  
  public void setH0B0(double[] b0, int _nr, int _L)
  {
    int i,l;
    if(_nr == n_rep && _L == L)
    {
      h0b0 = new double[(n_rep-1)*L];
      //--- scrape off first column,, only last n_rep-1 columns are nonzero
      for(i=1;i<n_rep;i++)
      {
        for(l=0;l<L;l++)
        {
          h0b0[(i-1)*L + l] = b0[L*i + l];
        }
      }   
      H0set = true;
    }
    else
    {System.out.println("Dimensions of the H0 filter must match the current filter");}
  }
  
  public void setIdentityH0()
  {
    
    h0b0 = new double[(n_rep-1)*L];    
    for(int i=0;i<n_rep-1;i++) {h0b0[L*i] = 1.0/(n_rep-1.0);}  
    H0set = true;
  }
  
  //---- correlations identity
  public void setIdentityH0(double[] cor)
  {
    
    h0b0 = new double[(n_rep-1)*L];    
    for(int i=0;i<n_rep-1;i++) {h0b0[L*i] = cor[i];}  
    H0set = true;
  }  
  
  
  public void setIdentityH0(int c)
  {
    if(c < n_rep)
    {
      h0b0 = new double[(n_rep-1)*L]; 
      h0b0[L*c] = 1.0;
      H0set = true;
    }  
  }
  
  public void setTrendH0() //------- sets the prior to the ideal trend
  {
    
    int i; 
    double sum;
    double cutoff2 = cutoff;
    h0b0 = new double[(n_rep-1)*L]; 
    
    h0b0[0] = cutoff2/Math.PI; sum= h0b0[0];
    for(i=1;i<L;i++)
    {h0b0[i] = (1/Math.PI)*Math.sin(cutoff2*i)/((double)i); sum = sum+h0b0[i];} 
    
    for(i=0;i<L;i++)
    {h0b0[i] = h0b0[i]/(sum+(sum-cutoff2/Math.PI));}  

    
  }  
  
   public void setForecastTrend(int lag)
   {

       int i,k,n; double sum = 0.0; double sum2 = 0.0; double sumi = 0.0;
       h0b0 = new double[(n_rep-1)*L];  
   
       for(k=0;k<L;k++)
       {
         sum=0.0;
         for(n=0;n<=K;n++)
         {
          sum = sum + Gamma[n]*Math.cos(-lag*Math.PI*n/K)*Math.cos(Math.PI*n*k/K);
          sumi = sumi + Gamma[n]*Math.sin(-lag*Math.PI*n/K)*Math.sin(Math.PI*n*k/K);
         }     
         if(k==0) {h0b0[0] = sum*sum;}
         else
         {h0b0[k] = 2.0*sum;} //(sum*sum + sumi*sumi);} 
         sum2 = sum2 + h0b0[k];
       }
       for(i=0;i<L;i++)
       {h0b0[i] = h0b0[i]/(sum2-h0b0[0]/2.0);}    
       
       for(i=0;i<L;i++)
       {
         System.out.println("0.0 " + h0b0[i] + " 0.0 0.0");    
       }
       
       
   }  
  
  
  
   /*----------------------------------------------------------
     Computes the filter and filtered data for the given 
     settings of Gamma, L, n_obs, lambda, expweight, cutoff
   -----------------------------------------------------------*/
   public void computeFilter(boolean com)
   {      
       int i,k,lag2,lag3,lag4,lag5; tdelay = 0.0;
       flength = n_obs-L+1;
       if(com)
       {
        double[] out = new double[output];        
        //System.out.println("length = " + tseries.length + ", nrep = " + n_rep + ", nobs = "+n_obs);          
        //for(i=0;i<n_obs;i++) {System.out.println(tseries[i]);}

        if(mdfa)
        {
               
         out = IMDFAreg(tseries, n_obs, L, n_rep, dd, DD, Gamma, shift, lambda, expweight, lambda_3, cutoff0,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, cross, Q_cdev, iter);
         //n_comps++;
         //System.out.println(n_comps);
         xf = new double[flength]; 
         x = new double[flength]; 
         for(i=0;i<flength;i++)  //---- get filtered and original series
         {
           xf[i] = out[i];
           x[i] = out[i+flength]; 
         }
         lag2 = 2*flength; lag3 = 2*flength+(n_rep+1)*K1; 

         for(i=0; i <= n_rep; i++)  //---- get amplitude/phase of filter
         {   
           for(k=0;k<=K;k++)
           {
             amp_filter[K1*i+k] = out[lag2 + K1*i + k];
             time_delay[K1*i+k] = -out[lag3 + K1*i + k];
           }
         }
         
         lag4 = lag2 + 2*(n_rep+1)*K1;
         for(i=0; i <=n_rep; i++)  //---- get frf
         {   
           gamma_hat[K1*i] = amp_filter[K1*i];
           for(k=1;k<=K;k++)
           {gamma_hat[K1*i+k] = amp_filter[K1*i+k]*Math.cos(time_delay[K*i+k]*(Math.PI*k/K));}

         }       
         //----------- frf ------------------
         for(i=0; i < n_rep; i++)  
         {
           for(k=0;k<L;k++)
           {b[L*i + k] = out[lag4 + L*i + k];}
         }

         //---- periodograms------------------
         lag5 = lag4 + n_rep*L;
         for(i=0; i < n_rep-1; i++) 
         {
           for(k=0;k<=K;k++)
           {period_hat[K1*i + k] = out[lag5 + K1*i + k];}
         }       
         criteria = out[out.length-3];
         degrees = out[out.length-2];      
         MDFAmin = out[out.length-1];
  
         for(k=1;k<=K;k++) {tdelay = tdelay + time_delay[K*n_rep + k];}
         tdelay = tdelay/K1;

       }
       else
       { 
         if(spec == 0) {spec_dens = new double[K1];}
         //System.out.println(cutoff0 + "  " + cutoff);
         //if(iter) {out = IDFAiter(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2);}
         //else {out = IDFAreg(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2, smooth, decay, 1);}          
         out = IDFAreg(tseries, n_obs, L, dd, DD, Gamma, shift, lambda, expweight, lambda_3, 
                        cutoff0, cutoff, Lag, i1, i2, smooth, decay,iter,spec_dens,spec);       
         lag2 = 2*flength;
         lag3 = 2*flength+K1;
         lag4 = lag3 + K1;

         xf = new double[flength]; 
         x = new double[flength]; 
        for(i=0;i<flength;i++)  //---- get filtered and original series
        {
           xf[i] = out[i];
           x[i] = out[i+flength]; 
        }
         
        for(k=0;k<=K;k++)
        {amp_filter[k] = out[lag2 + k]; time_delay[k] = -out[lag3 + k];}
                 
        gamma_hat[0] = amp_filter[0];
        for(k=1;k<=K;k++)
        {gamma_hat[k] = amp_filter[k]*Math.cos(time_delay[k]*Math.PI*k/K); tdelay = tdelay + time_delay[k];}
 
        for(k=0;k<L;k++) 
        {b[k] = out[lag4 + k];}

        criteria = out[out.length-3];
        degrees = out[out.length-2];      
        MDFAmin = out[out.length-1];

       }

       computeSampleIns();  
      }
   }
   
   public void setUpdateConstraints(int _i1, int _i2)
   {update_i1 = _i1; update_i2 = _i2;}
   

   public void updateSignalOut(int n, int l, double _lambda, double _alpha, double s, double d, double d2, double c)
   {
   
      int i,j,k; double sum=0; double om;
      int Kn = n/2;
      N_update = n; 
      L_update = l; 
      K_update = Kn;
      
      nrep_update = n_rep;
      
      int nr;
      int end = xf.length;
      int start = n_obs - n;
      update_start = start; 
      int flength_update = n-l+1;      
      double[] tsdata = new double[n_rep*n];   nr = n_rep;
      if(!mdfa) {tsdata = new double[2*n]; nr = 2;}
      update_signal = new double[flength_update];
      double[] cdev; double[] gamma;
 
      if(save_signal)
      {
        xf_orig = new double[xf.length];
        System.arraycopy(xf,0,xf_orig,0,xf.length);
        save_signal = false;
      }
      saveBfilter();
      b_update = new double[l*n_rep];
      amp_update = new double[(Kn+1)*(n_rep+1)];
      
      //---- compute generic gamma on new grid
      gamma = new double[Kn+1];
      for(k=0;k<=Kn;k++)
      {       
         om = (k*Math.PI/Kn);
         if(om < cutoff0) {gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff)
         {gamma[k] = 1.0;}
         else {gamma[k] = 0.0;}
      }
           
      //1) put the time series data in the first row
      for(i=start;i<n_obs;i++)
      {tsdata[i-start] = tseries[i];}
      
      //2) now put the filtered series in
      if(mdfa)
      {     
       cdev = computeQDeviation(l, n_rep, update_i1, update_i2);      
       //------ extract sequences ------     
       for(i=start;i<n_obs;i++)
       {       
        for(j=1;j<n_rep;j++)
        {
         sum = 0.0;
         for(k=0;k<L;k++)
         {sum = sum + b_old[L*j + k]*tseries[n_obs*j + i-k];}         
         tsdata[n*j + i-start] = sum;
        }
       }              
      }
      else
      {     
       cdev = new double[n_rep*l*n_rep*l];
       for(i=start;i<n_obs;i++)
       {
        sum = 0.0; 
        for(k=0;k<L;k++)
        {sum = sum + b_old[k]*tseries[i-k];}
        tsdata[n + i-start] = sum;
       }         
      }   
      double[] out = updateSignal(tsdata, n, l, nr, dd, DD, gamma, shift, _lambda, _alpha, 
                                  lambda_3, cutoff0, cutoff, Lag, update_i1, update_i2, w_const, s, d, d2, c, cdev);   
           
                             
      for(i=0;i<flength_update;i++) 
      {update_signal[i] = out[i];}// System.out.println(out[i]);}
     
      for(i=0;i<flength_update;i++) 
      {xf[end - 1 - i] = update_signal[flength_update-1-i];}
      
      int lag4 = flength_update*2;
      for(i=0; i < n_rep; i++)  
      {
        for(k=0;k<l;k++)
        {b_update[l*i + k] = out[lag4 + l*i + k];}
      }
      
      int lag5 = lag4+n_rep*l;
      for(i=0;i<n_rep;i++)
      {
       for(k=0;k<=Kn;k++)
       {
        amp_update[i*(Kn+1) + k] = out[lag5 + (Kn+1)*i + k];
       }
      }
      for(k=0;k<=Kn;k++)
      {amp_update[n_rep*(Kn+1) + k] = gamma[k];}
       
      mse_update = out[out.length-1]; 
      //System.out.println("new last point = " + xf[end - 1 - i]);
   } 
   
   
   public void updateSignalOut_Univ(int n, int l, double _lambda, double _alpha, double s, double d, double d2, double c)
   {
   
      int i,j,k; double sum=0; double om;
      int Kn = n/2;
      N_update = n; 
      L_update = l; 
      K_update = Kn;
      
      nrep_update = n_rep;
      
      int nr;
      int end = xf.length;
      int start = n_obs - n;
      update_start = start; 
      int flength_update = n-l+1;      
      double[] tsdata = new double[2*n]; nr = 2;
      update_signal = new double[flength_update];
      double[] cdev; double[] gamma;
 
      if(save_signal)
      {
        xf_orig = new double[xf.length];
        System.arraycopy(xf,0,xf_orig,0,xf.length);
        save_signal = false;
      }
      saveBfilter();
      b_update = new double[l];
      amp_update = new double[(Kn+1)*3];
      
      //---- compute generic gamma on new grid
      gamma = new double[Kn+1];
      for(k=0;k<=Kn;k++)
      {       
         om = (k*Math.PI/Kn);
         if(om < cutoff0) {gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff)
         {gamma[k] = 1.0;}
         else {gamma[k] = 0.0;}
      }
           
      //1) put the time series data in the first row
      for(i=start;i<n_obs;i++)
      {tsdata[i-start] = tseries[i];}
      
      //2) now put the filtered series in
           
      cdev = computeQDeviation(l, 2, update_i1, update_i2);      
       //------ extract sequences ------     
      for(i=start;i<n_obs;i++)
      {       
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {         
         for(k=0;k<L;k++)
         {sum = sum + b_old[L*j + k]*tseries[n_obs*j + i-k];}
        }
        tsdata[n + i-start] = sum;
      }              
      
  
      double[] out = updateSignal(tsdata, n, l, nr, dd, DD, gamma, shift, _lambda, _alpha, 
                                  lambda_3, cutoff0, cutoff, Lag, update_i1, update_i2, w_const, s, d, d2, c, cdev);   
           
                             
      for(i=0;i<flength_update;i++) 
      {update_signal[i] = out[i];}// System.out.println(out[i]);}
     
      for(i=0;i<flength_update;i++) 
      {xf[end - 1 - i] = update_signal[flength_update-1-i];}
      
      int lag4 = flength_update*2;
      for(k=0;k<l;k++)
      {b_update[k] = out[lag4 + l + k];}
      
      
      int lag5 = lag4+2*l;
      for(i=0;i<nr;i++)
      {
       for(k=0;k<=Kn;k++)
       {
        amp_update[i*(Kn+1) + k] = out[lag5 + (Kn+1)*i + k];
       }
      }
      for(k=0;k<=Kn;k++)
      {amp_update[2*(Kn+1) + k] = gamma[k];}
       
      mse_update = out[out.length-1]; 
      //System.out.println("new last point = " + xf[end - 1 - i]);
   } 
   
   
  
   //---------- Applies the updated filter to old-filtered data ----------
  
   public void applyUpdatedFilter(boolean upd)
   {
    
      if(b_update != null && upd)
      {     
        //---- first get old filtered data ------------
        int k,i,l,j; double sum;
        int end = xf.length;
        int _start = n_obs - N_update; 
        int start,obs;
        int add = _start - update_start;

        if(add > 0)
        {start = update_start; obs = N_update + add;}
        else
        {start = _start; obs = N_update;}
       
  
        int flength_update = obs - L_update + 1;
        double[] tsdata = new double[n_rep*obs];
        update_signal = new double[flength_update];
        
       if(!b_univ)
       {
        
        if(mdfa)
        {         
         //------ extract sequences ------     
         for(i=start;i<n_obs;i++)
         {       
          for(j=1;j<n_rep;j++)
          {
           sum = 0.0;
           for(k=0;k<L;k++) {sum = sum + b_old[L*j + k]*tseries[n_obs*j + i-k];}         
           tsdata[obs*j + i-start] = sum;
          }
         }              
        }
        else
        {     
          for(i=start;i<n_obs;i++)
          {
           sum = 0.0; 
           for(k=0;k<L;k++) {sum = sum + b_old[k]*tseries[i-k];}
           tsdata[obs + i-start] = sum;
          }         
        }        
  
  
        for(i=L_update-1;i<obs;i++)
        {
          sum = 0.0; 
          for(j=1;j<n_rep;j++)
          {     
            for(l=0;l<L_update;l++)
            {sum = sum + b_update[L_update*j + l]*tsdata[obs*j + i-l]; }   
          }    
          update_signal[i-(L_update-1)] = sum; 
        }
       
       }
       else
       {
         tsdata = new double[obs];
         for(i=start;i<n_obs;i++)
         { 
          sum = 0.0;
          for(j=1;j<n_rep;j++)
          {
           for(k=0;k<L;k++) {sum = sum + b_old[L*j + k]*tseries[n_obs*j + i-k];}                    
          }
          tsdata[i-start] = sum;  //System.out.println(sum);
         }   
         
           
  
         for(i=L_update-1;i<obs;i++)
         {
  
           sum = 0.0;
           for(l=0;l<L_update;l++)
           {sum = sum + b_update[l]*tsdata[i-l];}
           update_signal[i-(L_update-1)] = sum; 
         }
        //---- with new tsdata of length n, apply update filter of length l
       }
       
       for(i=0;i<flength_update;i++) 
       {xf[end - 1 - i] = update_signal[flength_update-1-i];} 
 
      }
   }
  

   public void swapUpdateAndOld(boolean old)
   {
      int i,end,flen;
      if(!save_signal)
      {
        if(old)
        {System.arraycopy(xf_orig,0,xf,0,xf_orig.length);}
        else
        {
           end = xf.length; flen = update_signal.length;
           for(i=0;i<flen;i++) 
           {xf[end - 1 - i] = update_signal[flen-1-i];}
        }
      }    
   }



   
   

   public void computeFilterGeneral(boolean com, boolean print)
   {      
       int i,k,lag2,lag3,lag4,lag5; tdelay = 0.0;
       max_amp = 0.0;       
       //System.out.println("n_obs = " + n_obs);
       if(com)
       {
        save_signal = true;
        flength = n_obs-L+1;
        double[] out = new double[output];        

        if(print)
        {
         System.out.println("length = " + tseries.length + ", nobs = " + n_obs + ", lambda = " + lambda + ", tseries[n_obs-1] = " + tseries[n_obs-1] + 
           ", tseries[tseries.length-1] = " + tseries[tseries.length-1] + ", n_rep = " + n_rep);          
         
         
         System.out.println("cut0 = " + cutoff0 + ", cut = " + cutoff);
         System.out.println("length Gamma = " + Gamma.length + ", expweight = " + expweight + ", regs = " + smooth + " " + decay + " " + decay2 + " " + cross);
         System.out.println("L = " + L + ", lag = " + Lag + ", i1 = " + i1 + ",i2 = " + i2 + " onestep = " + lambda_3);        
         System.out.println(tseries[n_obs-1] + ", " + tseries[n_obs-2] + ", " + tseries[n_obs-3] + ", " + tseries[n_obs-4]);
         System.out.println(tseries[n_obs + n_obs-1] + ", " + tseries[n_obs + n_obs-2] + ", " + tseries[n_obs + n_obs-3] + ", " + tseries[n_obs + n_obs-4]);
         System.out.println("");
        }
         //          
         //if(n_rep >= 4) System.out.println(tseries[n_obs-1] + ", " + tseries[n_obs + n_obs-1] + ", " + tseries[2*n_obs + n_obs-1] + ", " + tseries[3*n_obs + n_obs-1]);
        
        //if(print) System.out.println("MDFA = " + mdfa);
        if(mdfa)
        {
//          System.out.println("spec = " + spec);
//          for(i = 0; i<K1;i++)
//          {System.out.println(spec_dens[i] + " " + Gamma[i]);}
   
         if(useSD)
         {
           //for(i=0;i<100;i++) {System.out.println(tseries[i]);}
           out = SPEC_IMDFAreg(tseries, n_obs, L, n_rep, mod, arg, K, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma);
                       
           //System.out.println("compute with spec dens estimate, nrep = " + n_rep + ", length mod = " + mod.length);            
         }       
        
         else
         {       
          if(useH0 && H0set)
          {
           
           
           
           out = H0_IMDFAreg(tseries, n_obs, L, n_rep, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma, h0b0);
           
           //for(l=0;l<L*(n_rep-1);l++)
           //{System.out.println(h0b0[l]);}
//            if(decay2 > 0)
//            {out = H0_IMDFAreg(tseries, n_obs, L, n_rep, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
//                        cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
//                        ar_p, ma_q, ar_params, ma_params, innvar, arma, h0b0);}
//            else
//            {
//              out = GEN_IMDFAreg(tseries, n_obs, L, n_rep, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
//                        cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
//                        ar_p, ma_q, ar_params, ma_params, innvar, arma);           
//            }
           
          }
          else if((useH0 && !H0set) && (decay2 > 0))
          {
           //setTrendH0();
           setIdentityH0();
           
           //setForecastTrend(-1);
           System.out.println("Using the H0 filter computation - " + h0b0[0] + " " + h0b0[L-1]);
           
           out = H0_IMDFAreg(tseries, n_obs, L, n_rep, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma, h0b0);           
          }
          else
          {
           //for(i=0;i<100;i++) {System.out.println(tseries[i]);}
           out = GEN_IMDFAreg(tseries, n_obs, L, n_rep, dd, 0, Gamma, shift, lambda, expweight, lambda_3, diff_onestep,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma);
          }                                    
         }
         //n_comps++;
         //System.out.println(n_comps);
         xf = new double[flength]; 
         x = new double[flength]; 
         for(i=0;i<flength;i++)  //---- get filtered and original series
         {
           xf[i] = out[i];  
           x[i] = out[i+flength]; 
           
           //System.out.println(x[i] + "   " + xf[i]);
           
         }
         lag2 = 2*flength; lag3 = 2*flength+(n_rep+1)*K1; 

         for(i=0; i <= n_rep; i++)  //---- get amplitude/phase of filter
         {   
           for(k=0;k<=K;k++)
           {
             amp_filter[K1*i+k] = out[lag2 + K1*i + k];
             time_delay[K1*i+k] = -out[lag3 + K1*i + k];
             if(amp_filter[K1*i+k] > max_amp) {max_amp = amp_filter[K1*i+k];}
           }
         }
         
         lag4 = lag2 + 2*(n_rep+1)*K1;
         for(i=0; i <=n_rep; i++)  //---- get frf
         {   
           gamma_hat[K1*i] = amp_filter[K1*i];
           for(k=1;k<=K;k++)
           {gamma_hat[K1*i+k] = amp_filter[K1*i+k]*Math.cos(time_delay[K*i+k]*(Math.PI*k/K));}

         }       
         //----------- frf ------------------
         for(i=0; i < n_rep; i++)  
         {
           for(k=0;k<L;k++)
           {b[L*i + k] = out[lag4 + L*i + k];}
         }

         //---- periodograms------------------
         lag5 = lag4 + n_rep*L;
         for(i=0; i < n_rep-1; i++) 
         {
           for(k=0;k<=K;k++)
           {period_hat[K1*i + k] = out[lag5 + K1*i + k];} // if(i==0 && decay > 0) {System.out.println(period_hat[K1*i+k]);}}
         }       
         criteria = out[out.length-3];
         degrees = out[out.length-2];      
         MDFAmin = out[out.length-1];
  
         for(k=1;k<=K;k++) {tdelay = tdelay + time_delay[K*n_rep + k];}
         tdelay = tdelay/K1;
         //System.out.println("DONE!");
       }
       else
       { 
       
         //System.out.println("spec = " + spec);
         //for(i = 0; i<K1;i++)
         //{System.out.println(spec_dens[i]);}
       
         //if(spec == 0) {spec_dens = new double[K1];}
         //System.out.println(cutoff0 + "  " + cutoff);
         //if(iter) {out = IDFAiter(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2);}
         //else {out = IDFAreg(tseries, n_obs, L, dd, DD, Gamma, lambda, expweight, cutoff, Lag, i1, i2, smooth, decay, 1);}          

         if(useSD)
         {
            out = SPEC_IDFAreg(tseries, n_obs, L, mod, arg, K, dd, 0, Gamma, shift, lambda, expweight, lambda_3, 
                        cutoff0, cutoff, Lag, i1, i2, smooth, decay, decay2, iter,spec_dens,spec,
                        ar_p, ma_q, ar_params, ma_params, innvar, arma);
   

         }
         else
         {
           out = GEN_IDFAreg(tseries, n_obs, L, dd, 0, Gamma, shift, lambda, expweight, lambda_3, 
                        cutoff0, cutoff, Lag, i1, i2, smooth, decay, decay2, iter,spec_dens,spec,
                        ar_p, ma_q, ar_params, ma_params, innvar, arma);
         }               
         lag2 = 2*flength;
         lag3 = 2*flength+K1;
         lag4 = lag3 + K1;

         xf = new double[flength]; 
         x = new double[flength]; 
        for(i=0;i<flength;i++)  //---- get filtered and original series
        {
           xf[i] = out[i];
           x[i] = out[i+flength]; 
        }
         
        for(k=0;k<=K;k++)
        {
          amp_filter[k] = out[lag2 + k]; time_delay[k] = -out[lag3 + k];  
          if(amp_filter[k] > max_amp) {max_amp = amp_filter[k];}
        }
                 
        gamma_hat[0] = amp_filter[0];
        for(k=1;k<=K;k++)
        {gamma_hat[k] = amp_filter[k]*Math.cos(time_delay[k]*Math.PI*k/K); tdelay = tdelay + time_delay[k];}
 
        for(k=0;k<L;k++) 
        {b[k] = out[lag4 + k];}

        criteria = out[out.length-3];
        degrees = out[out.length-2];      
        MDFAmin = out[out.length-1];

       }

       computeSampleIns();  
      }
      //System.out.println("Moving on");
   }   
   






   public void computeQcdev()
   {   
      if(n_rep > 2) {Q_cdev = computeQDeviation(L, n_rep, i1, i2);}
      else {Q_cdev = new double[n_rep*L*n_rep*L];}
   }

   public void saveBfilter()
   {
     int le = b.length; 
     b_old = new double[le];
     System.arraycopy(this.b, 0, this.b_old, 0, le); 
   }

   public void computeRTSE(double[] _tseries, int N, int R)
   {
     
     int i,j,l;
     double sum=0;
     this.n_obs = N; this.n_rep = R; 
      
     int le = _tseries.length;
     this.tseries = new double[le];
     System.arraycopy(_tseries, 0, this.tseries, 0, le); 

     flength = n_obs - L + 1;
     xf = new double[flength];    

    if(mdfa)
    {     
      for(i=L-1;i<n_obs;i++)
      {
      sum = 0.0; 
      for(j=0;j<n_rep;j++)
      {
       for(l=0;l<L;l++)
       { 
         sum = sum + b_old[L*j + l]*tseries[n_obs*j + i-l];       
       }
      }    
      xf[i-(L-1)] = sum;  //System.out.println(sum);
     }
    }
    else
    {
      
      for(i=L-1;i<n_obs;i++)
      {
       sum = 0.0; 
       for(l=0;l<L;l++)
       {sum = sum + b_old[l]*tseries[i-l];}
       xf[i-(L-1)] = sum;  
      }    
      
    }   
   }

   /*----- computes the signal using given data and only returns final out number of signal values
           used for compute in-sample/out-of-sample error 
   ----------------------------------------------------------------------------------------------*/

   public double[] computeFilterSimple(double[] ts, int n, int _K, double[] tGamma)
   { 
        
        int i;
        int _flength = n-L+1; 
        double[] out; 
        double[] _xf,_x;
        _xf = new double[_flength];


        
        
        
        if(mdfa)
        { 
         double[] _spec_dens = new double[2*_K];

          for(i = 0; i < _K; i++)
          {_spec_dens[2*i] = 1.0; _spec_dens[2*i+1] = 0.0;}         
         
         
         if(useSD)
         {
           out = SPEC_IMDFAreg(ts, n, L, n_rep, mod, arg, _K, dd, 0, tGamma, shift, lambda, expweight, lambda_3, cutoff0,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, _spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma);      
         }       
         else
         {   
         
                    
           out = GEN_IMDFAreg(ts, n, L, n_rep, dd, 0, tGamma, shift, lambda, expweight, lambda_3, cutoff0,
                       cutoff, Lag, i1, i2, w_const, smooth, decay, decay2, cross, Q_cdev, iter, _spec_dens, spec,
                       ar_p, ma_q, ar_params, ma_params, innvar, arma);
         }

         _xf = new double[_flength]; 
         _x = new double[_flength]; 
         for(i=0;i<_flength;i++)  //---- get filtered and original series
         {
           _xf[i] = out[i]; 
           _x[i] = out[i+_flength]; 
         }
         //simPanel.plotData2(_x,_xf,_flength);
       }
       return _xf;
   }

   public double[] static_computeRTSE(double[] _tseries, int N, int R)
   {
     
     int i,j,l;
     double sum=0;
     int _flength = N - L + 1;
     double[] _xf = new double[_flength];    

    if(mdfa)
    {     
      for(i=L-1;i<N;i++)
      {
      sum = 0.0; 
      for(j=0;j<R;j++)
      {
       for(l=0;l<L;l++)
       { 
         sum = sum + b_old[L*j + l]*_tseries[N*j + i-l];       
       }
      }    
      _xf[i-(L-1)] = sum;  //System.out.println(sum);
     }
    }
    else
    {
      
      for(i=L-1;i<N;i++)
      {
       sum = 0.0; 
       for(l=0;l<L;l++)
       {sum = sum + b_old[l]*_tseries[i-l];}
       _xf[i-(L-1)] = sum;  
      }    
      
    }
    return _xf;   
   }




   public void setNDiv(double f, int n) {f_start = f; n_div = n;} 
   
//    public void time_freq_decomp()
//    {
//      int m,i; double[] output = new double[30];
//      flength = n_obs-L+1;
//     
//      if(mdfa) {
//       output = timefreqDecomposition(tseries, n_obs, L, n_rep, dd, DD, lambda, expweight, lambda_3, i2, 
//                                     w_const, smooth, decay, cross, f_start, n_div);
//      }
//      else {
//       output = timefreqDecompositiondfa(tseries, n_obs, L, dd, DD, lambda, expweight, lambda_3, i2, 
//                                     w_const, smooth, decay, cross, f_start, n_div);
//      }
// 
//      n_imfs = (int)output[output.length-1];
//      //System.out.println("number of n_imfs = "+n_imfs);
//      amMap = new double[n_imfs][flength]; fmMap = new double[n_imfs][flength];
//      ifreqMap = new double[n_imfs][flength]; phaseMap = new double[n_imfs][flength];
//      trend_cycle = new double[flength]; imfs = new double[n_imfs][flength];
// 
//      for(m=0;m<n_imfs;m++)
//      {
//       for(i=0;i<flength;i++)
//       {
//         amMap[m][i] = output[5*flength*m+i]; 
//         fmMap[m][i] = output[5*flength*m+flength+i];     
//         phaseMap[m][i] = output[5*flength*m+2*flength+i]; 
//         ifreqMap[m][i] = output[5*flength*m+3*flength+i];
//         imfs[m][i] = output[5*flength*m+4*flength+i];
//       }
//      }
// 
//      for(i=0;i<flength;i++)
//      {trend_cycle[i] = output[5*flength*n_imfs+i];}
//   
//    }


   public double[] computeOptimalFilter(double[] sym_sig, double w1, double w2, double w3)
   {
      return findOptimalFilter(tseries, n_obs, K, L, n_rep, dd, DD, cutoff, cutoff0, Gamma,
	smooth, decay, cross, i1, i2, w_const, sym_sig, w1, w2, w3, lambda, Math.abs(expweight));
   } 

   public double[] filterScore(double[] xf, double[] sym, int thresh)
   {
      return computeScore(xf.length, xf, sym, thresh);
   }


   
   public native double[] findOptimalFilter(double[] data, int N, int _K, int _L, 
	int _n_rep, int dd, int DD, double cutoff, double cutoff0, double[] gamma,
	double smooth, double decay, double cross, int i1, int i2, double[] weight_constraint, 
	double[] sym_sig, double w1, double w2, double w3, double lambda_in, double alpha_in);

   public native double[] computeScore(int _N, double[] _xf, double[] _sym, int _th);

   public native double[] IMDFAold(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _lambda, double _expweight, double _cutoff, int _Lag, int _i1, int _i2);

   public native double[] IMDFAiter(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _lambda, double _expweight, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb);

   public native double[] IDFA(double[] data, int _nObs, int _L, int _dd, int _DD, double[] _gamma, 
                               double _lambda, double _expweight, double _cutoff, int _Lag, int _i1, int _i2);

   public native double[] IDFAiter(double[] data, int _nObs, int _L, int _dd, int _DD, double[] _gamma, 
                               double _lambda, double _expweight, double _cutoff, int _Lag, int _i1, int _i2);

   public native double[] IMDFAreg(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb,
                               double _smooth, double _decay, double _cross, double[] _Qcdev, boolean iter);

   public native double[] IDFAreg(double[] data, int _nObs, int _L, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2,
                               double _smooth, double _decay, boolean iter, double[] spec_dens, int spec);
                              
   public native double[] computeQDeviation(int L, int n_rep, int i1, int i2); 

//    public native double[] timefreqDecomposition(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, 
//                                                double _lambda, double _expweight, double _lambda3, int _i2, double[] _comb, 
//                                                double smooth, double decay, double cross, double f_start, int _n_div); 
// 
//    public native double[] timefreqDecompositiondfa(double[] data, int _nObs, int _L,  int _dd, int _DD, 
//                                                double _lambda, double _expweight, double _lambda3, int _i2, double[] _comb, 
//                                                double smooth, double decay, double cross, double f_start, int _n_div);    
// 

   public native void tradingDiffSingle(double[] data, double[] signal, double[] account, int n_obs, int ahead, double cost, 
                                        int _short, double _risk, double[] _trades);

                                               
   //------------------ More general form of ZPC-IMDFA --------------------------------------------------------------
   // Input: spec_dens is length 2*K1 output of concurrent filter from ZPC, holds real and imag coefficients
   // Input: arma p,q estimate of data
                                               
   public native double[] SPEC_IMDFAreg(double[] data, int _nObs, int _L, int _n_rep, double[] mod, double[] arg, int _K, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb,
                               double _smooth, double _decay, double _decay2, double _cross, double[] _Qcdev, boolean iter, double[] spec_dens, int spec,
                               int p, int q, double[] ar, double[] ma, double innvar, int arma);

   public native double[] SPEC_IDFAreg(double[] data, int _nObs, int _L, double[] mod, double[] arg, int _K, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2,
                               double _smooth, double _decay, double _decay2, boolean iter, double[] spec_dens, int spec, 
                               int p, int q, double[] ar, double[] ma, double innvar, int arma);                                               
                                                       

      
   public native double[] GEN_IMDFAreg(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb,
                               double _smooth, double _decay, double _decay2, double _cross, double[] _Qcdev, boolean iter, double[] spec_dens, int spec,
                               int p, int q, double[] ar, double[] ma, double innvar, int arma);

   public native double[] GEN_IDFAreg(double[] data, int _nObs, int _L, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2,
                               double _smooth, double _decay, double _decay2, boolean iter, double[] spec_dens, int spec, 
                               int p, int q, double[] ar, double[] ma, double innvar, int arma);                                               

   public native double[] H0_IMDFAreg(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb,
                               double _smooth, double _decay, double _decay2, double _cross, double[] _Qcdev, boolean iter, double[] spec_dens, int spec,
                               int p, int q, double[] ar, double[] ma, double innvar, int arma, double[] b0); 
 
                               
   public native double[] gridSearch(double[] data, double[] gamma, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, double[] mod, double[] arg, int spc);                                            
                                               
   public native double[] gridSearchCutoff(double[] data, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, double[] mod, double[] arg, int spc);                                            
                                               

   public native double[] gridSearchBandCutoff(double[] data, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, double[] mod, double[] arg, int spc);                                              
                                               
                                               
   public native double[] gridSearchMultiBandCutoff(double[] data, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, double[] mod, double[] arg, int spc);                     
                                               
   public native double[] gridSearchRegularization(double[] data, double[] gamma, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, int var, double[] mod, double[] arg, int spc);    
                              
   public native double[] optimizeRegularization(double[] data, double[] gamma, int _nObs, int _nreps, int _K, int _L, double _shift, double lambda, double alpha, double _cutoff0, double _cutoff, 
                                               int _Lag, int _i1, int _i2, double _smooth, double _decay, double _decay2, double _cross, double[] price, int _opt, int ss, double tc, int var);    
 
                              
   public native double[] updateSignal(double[] data, int _nObs, int _L, int _n_rep, int _dd, int _DD, double[] _gamma, 
                               double _shift, double _lambda, double _expweight, double _lambda3, double _cut0, double _cutoff, int _Lag, int _i1, int _i2, double[] _comb,
                               double _smooth, double _decay, double _decay2, double _cross, double[] _Qcdev);              

                               
   public native static double[] optimalPortfolio(double[] data, int obs, int nreps);                             

   
   static {System.loadLibrary("I_MDFA");}


   public static void main(String args[])
   {
      int n_obs = 300; int n_rep = 1; int L = 28; int Lag = 0; int i1 = 0; int i2 = 0;
      double lambda = 0; double expweight = 0; double cutoff = Math.PI/6.0;
      int K = (int)n_obs/2; int K1 = K+1; 

      int i,j,k,n_toks;
      String[] tokens;
      String strline;
      Double D;
      double val;
      double omegak;
      String delims = "[ ]+";
      //------------ Build symmetric Gamma---------------
      double[] Gamma = new double[K1];
      double[] bc = new double[L+1];
      double sum;

      bc[0] = cutoff/Math.PI; sum= bc[0];
      for(i=1;i<=L;i++)
      {bc[i] = (1/Math.PI)*Math.sin(cutoff*i)/i; sum= sum+bc[i];} 
      for(i=0;i<=L;i++)
      {bc[i] = bc[i]/(sum+(sum-cutoff/Math.PI));}     
      for(k=0; k<=K;k++)
      {
        omegak=k*Math.PI/K;	
        sum=bc[0];
        for(j=1;j<=L;j++) {
          sum = sum + bc[j]*Math.cos(omegak*j)*2.0;
        }
        if(omegak <= cutoff)
        {Gamma[k] = 1.0;} //System.out.println(Gamma[k]);
        else 
        {Gamma[k] = 0.0;}
      }
 
      int n_r = 5; //number of series for mdfa      
      //------------ Get data ---------------
      double[] x1 = new double[n_r*n_obs];
      double[] x = new double[n_obs];
      double[] y = new double[n_r*900]; //all the data
 
      FileInputStream fin;
      DataInputStream din;
      BufferedReader br;
     
      System.loadLibrary("I_MDFA");

      File file = new File("data2.dat");
      int ncount = 0;
      try{

      fin = new FileInputStream(file);
      din = new DataInputStream(fin);
      br = new BufferedReader(new InputStreamReader(din)); 

       while((strline = br.readLine()) != null)
       {
         tokens = strline.split(delims); 
         n_toks = tokens.length; 
         if(n_toks == 0)
         {System.out.println("End of file"); break;}
  
         //D = new Double(tokens[n_toks-1]);
         //val = D.doubleValue();
         for(i=0;i<n_toks;i++) //n_toks should be 5
         {
           D = new Double(tokens[i]);
           val = D.doubleValue(); 
           y[900*i + ncount] = val;
         }
         ncount++;
       }
       din.close();
       fin.close();
       
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}

      for(i=n_obs;i<n_obs+300;i++)
      { 
         for(j=0;j<n_r;j++)
         { x1[n_obs*j + (i-n_obs)] = y[900*j + i]; }
         x[i-n_obs] = y[i]; //System.out.println(x[i-n_obs]);
      }

      L=48;
      n_rep = n_r; lambda = 400; i1 = 0; i2 = 1;
      IMDFA mdfa = new IMDFA(n_obs, n_rep, L, Lag, lambda, expweight, cutoff, i1, i2); 
      mdfa.set_dfaIter(true);   
      mdfa.set_Gamma(Gamma);
      mdfa.set_tseries(x1,n_obs,n_rep);
      mdfa.set_bconstraints(0,0);
      
    
   }
 
  //--------------------- Compute sample Periodogram for x and xf

   public void computeSampleIns()
   {
     //System.out.println("computes periodogram");
     int i,j;
     int n_obsd; double mean_x,mean_xf; 
     double cdenom;  
     double sum = 0.0; double sum2 = 0.0; 
     double sumr, sumi, sumr2, sumi2;
     double[] d_seriesxf;
     double[] d_seriesx;
     int samp1 = samp+1;

     diff_stop = 0.0; diff_band = 0.0;
     int w0 = (int)(cutoff0*samp1/Math.PI);
     int w1 = (int)(cutoff*samp1/Math.PI);
     

     period_xf = new double[2*samp1];
   
     if((dd == 1) && (DD==1))
     { 
       n_obsd =  flength - S - 1; 
       d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd];  

       for(i=S; i < flength-1; i++)
       {
        d_seriesx[i-S] = (x[i+1] - x[i]) - (x[i+1-S] - x[i-S]);  sum = sum + d_seriesx[i-S];
        d_seriesxf[i-S] = (xf[i+1] - xf[i]) - (xf[i+1-S] - xf[i-S]);  sum2 = sum2 + d_seriesxf[i-S];
       }
       mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else if(dd == 1)
     {
      n_obsd = flength - 1; 
      d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd]; 
      for(i=1; i < flength; i++)
      {
        d_seriesx[i-1] = x[i] - x[i-1]; sum = sum + d_seriesx[i-1];
        d_seriesxf[i-1] = xf[i] - xf[i-1]; sum2 = sum2 + d_seriesxf[i-1];
      }
      mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else if(DD == 1)
     {
        n_obsd = flength - S; 
        d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd];    
        for(i=S; i < flength; i++)
        {
          d_seriesx[i-S] = x[i] - x[i-S]; sum = sum + d_seriesx[i-S];
          d_seriesxf[i-S] = xf[i] - xf[i-S]; sum2 = sum2 + d_seriesxf[i-S];
        }
        mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else
     {
       n_obsd = flength; 
       d_seriesx = new double[n_obsd];  
       d_seriesxf = new double[n_obsd];     

       for(i=0; i < flength; i++)
       {
          d_seriesx[i] = x[i]; sum = sum + d_seriesx[i];
          d_seriesxf[i] = xf[i]; sum2 = sum2 + d_seriesxf[i];
       }
       mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;       
     }  
     
     period_xf[0] = mean_x*mean_x/Math.sqrt(2*Math.PI*n_obsd);
     period_xf[samp1] = mean_xf*mean_xf/Math.sqrt(2*Math.PI*n_obsd);

     if(spec == 1) {period_xf[0] = spec_dens[0];}

     if(w0 > 0) {diff_stop = period_xf[samp1];}
     else {diff_band = Math.abs(period_xf[0] - period_xf[samp1]);}
     
     //----------  Testing different transform -----------
     for(j=1;j<samp1;j++)
     {
  
       sumr = 0.0; sumi=0.0;   sumr2 = 0.0; sumi2 =0.0; 
       for(i=0;i<n_obsd;i++)
       {
        sumr = sumr + d_seriesx[i]*Math.cos(Math.PI*(i+1)*j/samp); 
        sumi = sumi + d_seriesx[i]*Math.sin(Math.PI*(i+1)*j/samp);
        sumr2 = sumr2 + d_seriesxf[i]*Math.cos(Math.PI*(i+1)*j/samp); 
        sumi2 = sumi2 + d_seriesxf[i]*Math.sin(Math.PI*(i+1)*j/samp);
       }    

       period_xf[j] = (sumr*sumr + sumi*sumi)/Math.sqrt(2*Math.PI*n_obsd);
       period_xf[j+samp1] = (sumr2*sumr2 + sumi2*sumi2)/Math.sqrt(2*Math.PI*n_obsd);

       if(spec == 1) {period_xf[j] = spec_dens[j];}

       if(dd>0) 
       {
        cdenom = ((1.0 - Math.cos(Math.PI*j/samp))*(1.0 - Math.cos(Math.PI*j/samp)) + Math.sin(Math.PI*j/K)*Math.sin(Math.PI*j/samp));        
        period_xf[j] = period_xf[j]/cdenom; period_xf[j+samp1] = period_xf[j+samp1]/cdenom;
        if(spec == 1) {period_xf[j] = spec_dens[j];}
       }
       if(DD > 0) 
       {
         cdenom = ((1.0 - Math.cos(12*Math.PI*j/samp))*(1.0 - Math.cos(12*Math.PI*j/samp)) + Math.sin(12*Math.PI*j/samp)*Math.sin(12*Math.PI*j/samp));         
         period_xf[j+samp1] = period_xf[j+samp1]/cdenom; period_xf[j] = period_xf[j]/cdenom;
         if(spec == 1) {period_xf[j] = spec_dens[j];}
       }
 
       if((j < w0) || (j > w1))      
       {diff_stop = diff_stop + period_xf[samp1+j];}
       else
       {diff_band = diff_band + Math.abs(period_xf[j] - period_xf[samp1+j]);}      
      }
    } 
    
    
    /*-------------------------------------------------------
       For high-frequency trading, this function sets a much longer 
       data range of high-frequency log-return data than the in-sample. 
       
       The idea is to build a smoothed peridogram of the larger data sample
    
    ---------------------------------------------------------*/
    public void setFullRangeData(double[] full)
    {
    
       if(full.length > n_obs)
       {
          fullRangeData = new double[full.length];          
          System.arraycopy(full, 0, fullRangeData, 0, full.length);
          System.out.println("Periodogram loaded");
       }
       else
       {
         System.out.println("Periodogram not loaded - " + full.length + " < " + n_obs);
       }
    }
    
    /*---------------------------------------------------------
    Here we compute the smoothed periodogram on the targeted frequency set. 
    We can't include future information though */
        
    public void computeSpecDensPeriodogram()
    {
    
       int i,j;
       int n_obsd; 
       double sumr, sumi;    
       double max = 0.0;
       
       spec_dens = new double[2*K1];
       for(i = 0; i < K1; i++)
       {spec_dens[2*i] = 1.0; spec_dens[2*i+1] = 0.0;}       
      
       if(fullRangeData != null)
       {
         n_obsd = fullRangeData.length;
         for(j=0;j<K1;j++)
         {          

           sumr = 0.0; sumi=0.0;   
           for(i=0;i<n_obsd;i++)
           {
             sumr = sumr + fullRangeData[i]*Math.cos(Math.PI*(i+1)*j/K); 
             sumi = sumi + fullRangeData[i]*Math.sin(Math.PI*(i+1)*j/K);     
           }
           
           spec_dens[2*j] = (sumr*sumr + sumi*sumi)/Math.sqrt(Math.PI*n_obsd); 
           spec_dens[2*j+1] = 0.0;
           if(spec_dens[2*j] > max) {max = spec_dens[2*j];}           
         }       
         spec = 1;

         //-------------normalize
         for(j=0;j<K1;j++) {spec_dens[2*j] = spec_dens[2*j]/max;}// System.out.println(spec_dens[2*j]);}      
       }   
    }
    

}

