package ch.imetrica.bayesCronos;

public class HeavyModel
{

  public int n_obs;  //number of observations
  public int n_rep;  //number of data series (2 or 3) 

  public double[] tseries; //A T x n_rep vector of time series 
  public double[] ht;      //A T x n_rep vector of volatilty means, (h_t, u_t, phi_t)   

  public int seed; 
  public int burnin; 

  public int retrack;
  public int n_params;
  

  public int[] A_dim;  //matrix for the data innovations
  public int[] B_dim;  //matrix for the autoregressive components
  public int rt; 

  //---- parameters for the model ----------------------------
  public double alpha, alpha_R, alpha_SR, lambda_1, lambda_2; 
  public double beta, beta_R, beta_SR;
  public double[] W_mean; 
  public double[] alphas; 
  public double[] betas;   
  public double[] lambdas;
  public int n_alpha;
  public double[] parameters;

  //---- Output of heavy_model --------------------------------
  public double[] F_dist;  //empirical distribution of innovation sequence   
  public double[] llht;    //maximum likelihood per observation
  public double[] r_forecasts;  //n_fore series of return forecasts; 
  public double[] h_forecasts;  //K x n_steps forecasts of each mean volatilty series
  public double[] score;        //score of likelihood function
  public double llh; 
  
  //---- Forecast inputs ---------------------------------------
  public int n_fore; 
  public int n_steps; 
  


  public HeavyModel(int _n, int _seed)
  {
    n_obs = _n; 
    n_rep =  2;
    seed = _seed; 
   
    retrack = 0;
    burnin = 100;

    tseries = new double[n_rep*n_obs];
    ht      = new double[n_rep*n_obs];

    //---- default parameters ---------
    n_fore = 3;
    n_steps = 10;

    A_dim = new int[n_rep*n_rep];
    
    A_dim[0] = 0; A_dim[1] = 1; A_dim[2] = 0; A_dim[3] = 1; 
    rt = 0;
  
        

    //[.15 .05 .2 .4 .7 .55]


  }

  public void changeSeed(int s) {seed = s;}

  public void setDimensions(int _K)
  {

    n_rep = _K;

    tseries = new double[n_rep*n_obs];
    ht      = new double[n_rep*n_obs];

    //---- default parameters ---------
    A_dim = new int[n_rep*n_rep];

    if(n_rep == 2)
    {
       A_dim[0] = rt; A_dim[1] = 1; 
       A_dim[2] = 0;  A_dim[3] = 1;
       n_alpha = 2+rt;
    }         
    else if(n_rep == 3)
    {
       A_dim[0] = rt; A_dim[1] = 1; A_dim[2] = 0; 
       A_dim[3] = 0;  A_dim[4] = 1; A_dim[5] = rt;
       A_dim[6] = 0;  A_dim[7] = 0; A_dim[8] = 1;
       n_alpha = 3+2*rt;
    }
    n_params = 2*n_rep + n_alpha;
    parameters = new double[n_params];
  }


  public void setParameterValues(double w1, double w2, double alpha, double alpha_R, double lambda, double beta, double beta_R)  
  {
    
    n_rep = 2; rt = 0;
    W_mean = new double[2]; W_mean[0] = w1; W_mean[1] = w2;

    if(lambda > 0.0) //include rt information
    {
      alphas = new double[3]; rt = 1; 
      alphas[0] = lambda; alphas[1] = alpha; alphas[2] = alpha_R; 
    }
    else //don't include rt information
    {
      alphas = new double[2]; rt = 0; 
      alphas[0] = alpha; alphas[1] = alpha_R;
    }         

    betas = new double[2]; 
    betas[0] = beta; betas[1] = beta_R;      
    
    setDimensions(n_rep);
     
  }


  public void setParameterValuesSV(double w1, double w2, double w3, double alpha, double alpha_R, double alpha_SR,
                                   double lambda, double beta, double beta_R, double beta_SR, double lambda_2)  
  {
    
    n_rep = 3; rt = 0;
    W_mean = new double[3]; W_mean[0] = w1; W_mean[1] = w2; W_mean[2] = w3;

    if(lambda > 0.0) //include rt information
    {
      alphas = new double[5]; rt = 1; 
      alphas[0] = lambda; alphas[1] = alpha; alphas[2] = alpha_R; 
      alphas[3] = lambda_2; alphas[4] = alpha_SR;
    }
    else //don't include rt information
    {
      alphas = new double[3]; rt = 0; 
      alphas[0] = alpha; alphas[1] = alpha_R; alphas[2] = alpha_SR; 
    }         

    betas = new double[3]; 
    betas[0] = beta; betas[1] = beta_R; betas[2] = beta_SR;     

    setDimensions(n_rep);
     
  }


  public void parameterTransform()
  {
    int i;

    lambdas = new double[(n_rep-1)*rt];
    for(i=0; i<n_rep;i++) {W_mean[i] = parameters[i];}

    if(rt == 1)
    {lambdas[0] = parameters[n_rep];}

    for(i=0;i<2;i++)
    {alphas[i] = parameters[n_rep+rt+i];}

    if(rt == 1 && n_rep == 3)
    {
      lambdas[1] = parameters[n_rep+rt+2];
      alphas[2] = parameters[n_rep+rt+3];
    }
    else if(n_rep == 3)
    {alphas[2] = parameters[n_rep+2];}

    for(i=0; i<n_rep;i++) {betas[i] = parameters[n_rep + n_rep + (n_rep-1)*rt + i];}
   
  }

  public void setForecastDimensions(int fore, int steps)
  {
     n_fore = fore; n_steps = steps; 
     r_forecasts = new double[n_fore*n_steps];
     h_forecasts = new double[n_rep*n_steps];
  }


  public void setData(int n, int K, double[] _tseries)
  {
    int i,j;
    n_obs = n; 
    n_rep = K;

    tseries = new double[n_obs*n_rep];
    ht = new double[n_obs*n_rep];

    for(i = 0; i < K; i++)
    {
     for(j = 0;j < n; j++)
     {tseries[K*j+i] = _tseries[K*j+i];}     
    }

    llht = new double[n_obs];
    F_dist = new double[n_obs*K];

  }


  public void setNobs(int n, int K)
  {
    n_obs = n; 
    n_rep = K;

    tseries = new double[n_obs*n_rep];
    ht = new double[n_obs*n_rep];  
  }

  public void setRetrackParameter(int r) {retrack = r;}

  public void estimateHeavyModel()
  { 

      computeHeavyModel(tseries, n_obs, n_rep, n_alpha, n_params, parameters, llht, llh, 
                         ht, r_forecasts, h_forecasts, n_fore, n_steps, F_dist, seed, W_mean, alphas, betas, A_dim, retrack);

      llh = llht[0];
      //System.out.println("min = " + llh);
  }

  
  public double[] simulateHeavy(int m2)
  {
      double[] sim_data = new double[n_obs*n_rep]; ht = new double[n_obs*n_rep];

      //System.out.println("n = " + n_obs + " n_rep = " + n_rep + " m2 = " + m2 + "W = " + W_mean[0] + 
      //" " + W_mean[1] + " alphas = " + alphas[0] + "  " + alphas[1] + ", A_dim = " + A_dim[0]  + "  " + A_dim[1]  + "  " 
      //+ A_dim[2]  + "  " + A_dim[3]); 

      simulateHeavyModel(sim_data, ht, n_obs, n_rep, m2, W_mean, alphas, betas, A_dim, n_alpha, seed);
      return sim_data; 
  }

 
  public native void computeHeavyModel(double[] _data, int n, int k, int nalphas, int nparams, double[] params, 
                                       double[] _llht, double llh, double[] _ht, double[] r_fore, 
                                       double[] h_fore, int nfore, int nsteps, double[] _F, int _seed,
                                       double[] _w, double[] _alphas, double[] _betas, int[] _A, int ret);


  public native void simulateHeavyModel(double[] _data, double[] _ht, int n_obs, int n_rep, int m2, 
                                        double[]_w, double[] _alpha, double[] _beta, int[] _A, int _nalpha, int _seed);


  static {System.loadLibrary("heavy");}


  public static void main(String args[])
  {
  
     int i,k;
     int n_obs = 600; int seed = 234;
     int n_fore = 3; int n_steps = 20;     
     int n_series = 2;
 
     int m2 = 76;
     System.loadLibrary("heavy");
     
     //{.15, .05, .2, .4, .7, .55};
     double w1 = .15; 
     double w2 = .05;
     double lambda = .1;
     double alpha = .2;
     double alpha_R = .4;
     double beta = .7;
     double beta_R = .55;
     

     //---  start the model engine -------------------
     HeavyModel heavy = new HeavyModel(n_obs, seed);  

     //--- set forecast dimensions -------------------       
     heavy.setForecastDimensions(n_fore, n_steps);
     
     //--- set (simulation or initialization) parameter values -------------
     heavy.setParameterValues(w1, w2, alpha, alpha_R, lambda, beta, beta_R); 

     //--- simulate the heavy model with 76? for chisquared process---------
     double[] series = heavy.simulateHeavy(m2); 
     

     //---- now estimate ----------------------

      
     heavy.setParameterValues(w1, w2, alpha, alpha_R, lambda, beta, beta_R); 
     heavy.setData(n_obs, n_series, series);
     heavy.setRetrackParameter(0);
     heavy.estimateHeavyModel(); 
     


     for(i=0;i<n_obs;i++)
     {
       System.out.println(heavy.F_dist[i] + "  " + heavy.F_dist[n_obs + i]);
     }
     
     for(k=0;k<n_series;k++)
     {
      System.out.println("FORECAST_H");
      for(i=0;i<n_steps;i++)
      {
        System.out.println(heavy.h_forecasts[n_steps*k + i]);
      }     
     }
     

     System.out.println("\n min = " + heavy.llh);
     System.out.println("Number of params = " + heavy.parameters.length);

     for(i=0;i<heavy.parameters.length;i++) 
     {System.out.println(heavy.parameters[i]);}

  }

}
