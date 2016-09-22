package ch.imetrica.mdfa;


  /*------------------------------------------------
      Object baraye negah dari ye etele'at e filter      
  --------------------------------------------------*/


public class MDFAFilter
{
  
    //---- Store Data and filter controls -------------
    int n_obs, n_rep, K, L, Lag, S, K1;
    int i1, i2, dd, DD, flength;  
    boolean mdfa, mdfa_iter; 
    
    //---- extra Filter controls ----------------
    double lambda, expweight, lambda_3;   
    double smooth, decay, cross;
    double[] w_const;
    boolean dfa;

    //---- criteria ----------------------------
    double criteria, degrees;
    double MDFAmin;

    //---- data --------------------------------
    double[] xf; 
    double[] b;

    //---- frequency --------------------------
    double cutoff, cutoff0;
    double[] Gamma; 

    public MDFAFilter(int n, int rep, int _L, int _Lag)
    {
       n_obs = n; n_rep = rep; L = _L; Lag = _Lag; 
       w_const = new double[n_rep-1]; 
       xf = new double[n_obs - (L-1)]; 
       K = n_obs/2; K1 = K+1;
       flength = n_obs - (L-1);
       Gamma = new double[K+1]; 
       b = new double[n_rep*L];
       dfa = false;
    }


}