package ch.imetrica.mdfaTradingStrategies;


public class FilterParameters
{

  String[] name;
  int n_obs,L,lag, n_rep;
  double expweight;
 
  double cross;
  double lambda; 
  double smooth; 
  double decay;
  double decay2;
  double cutoff;
  double[] Gamma;
  
  public FilterParameters(String[] nm, int nobs, int nr, int _L, int _lag, double freq, double lam, double exp, double sm, double dec1, double dec2, double cro)
  {
   name = new String[nm.length];
   
   for(int i=0;i<nm.length;i++) {name[i] = nm[i];}
   
   n_rep = nr; lag = _lag; n_obs = nobs; L = _L; cutoff = freq; lambda = lam; expweight = exp; smooth = sm; decay = dec1; decay2 = dec2; cross = cro;
   
   int k,K;
   int K1 = (int)(n_obs/2) + 1; K=K1-1;
   double om; double cutoff0 = 0;
   Gamma = new double[K1];       
   for(k=0; k<K1;k++)
   {       
     om = (k*Math.PI/K);
     if(om < cutoff0) {Gamma[k] = 0.0;}
     else if(om >= 0 && om <= cutoff)
     {Gamma[k] = 1.0;}
     else
     {Gamma[k] = 0.0;}
     
     //System.out.println(Gamma[k]);
   }
   
   //System.out.println(name[0] + ", length of GAMMA " + Gamma.length + ", length of n_obs = " + n_obs + "length of K = " + K1);
   
   
  }

  String getName() {return name[0];}
  
  void printFilterSpec()
  {
  
    System.out.println("Name: " + name[0]);
    System.out.println("N: " + n_obs);
    System.out.println("n_rep: " + n_rep);
    System.out.println("Gamma length: " + Gamma.length);
    System.out.println("Expweight: " + expweight);
    System.out.println("L: " + L);
  
  }
  
  
}