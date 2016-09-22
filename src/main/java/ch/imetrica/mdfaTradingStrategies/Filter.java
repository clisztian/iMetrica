package ch.imetrica.mdfaTradingStrategies;


public class Filter
{

  int L; int n_rep;
  public double[] b_coeffs;

  public Filter(int _L, int rep)
  {L = _L; n_rep = rep; b_coeffs = new double[L*n_rep];}

  public void setCoeffs(double[] b, int _L, int rep)
  {
   b_coeffs = new double[b.length];
   L = _L; n_rep = rep;
   
   System.arraycopy(b,0,b_coeffs,0,b_coeffs.length);
  } 

}