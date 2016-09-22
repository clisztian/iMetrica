package ch.imetrica.mdfaTradingStrategies;


public class StrategyParameters
{

  String name;
  int startTradingInt;
  String insampStart, startTrading, endTrading;
  int n_obs,L,lag, n_rep;
  double expweight;
  double cross;
  double lambda; 
  double smooth; 
  double decay;
  double decay2;
  double cutoff;
  double cutoff0;
  double hybrid_weight,hybrid_weight_diff;
  double time_shift;
  int b0trend;
  int i1,i2;
  int sig_diff;
  int sig_inv;
  int stop_loss;
  int take_profit = 0;
  double meanRevert = -.1;
  
  public StrategyParameters(String filt_name, String insamp, String start, int start_int, String end, int nobs, int nr, int _L, int _lag, double freq, double lam, 
                            double exp, double sm, double dec1, double dec2, double cro,
                            int _i1, int _i2, double phase, double hybrid, double hybrid2, int b0, int sd, int sigi, int st)
  {
   
   startTradingInt = start_int; //0-20
   name = filt_name;
   insampStart = insamp;
   startTrading = start; 
   endTrading = end;
   hybrid_weight = hybrid;
   hybrid_weight_diff = hybrid2; 
   time_shift = phase; 
   cutoff0 = 0.0;
   n_rep = nr; lag = _lag; n_obs = nobs; L = _L; 
   cutoff = freq; lambda = lam; expweight = exp; 
   smooth = sm; decay = dec1; decay2 = dec2; cross = cro;
   i1 = _i1; i2 = _i2; b0trend = b0; //--- mdfa Stages on
   sig_diff = sd; //number of stages stages
   sig_inv = sigi;
   stop_loss = st;
   
  }

  public StrategyParameters(String filt_name, String insamp, String start, double meanR, String end, int nobs, int nr, int _L, int _lag, double freq, double lam, 
                            double exp, double sm, double dec1, double dec2, double cro,
                            int _i1, int _i2, double phase, double hybrid, double hybrid2, int b0, int sd, int sigi, int st)
  {
   
   meanRevert = meanR;
   startTradingInt = 0; //0-20
   name = filt_name;
   insampStart = insamp;
   startTrading = start; 
   endTrading = end;
   hybrid_weight = hybrid;
   hybrid_weight_diff = hybrid2; 
   time_shift = phase; 
   cutoff0 = 0.0;
   n_rep = nr; lag = _lag; n_obs = nobs; L = _L; 
   cutoff = freq; lambda = lam; expweight = exp; 
   smooth = sm; decay = dec1; decay2 = dec2; cross = cro;
   i1 = _i1; i2 = _i2; b0trend = b0; //--- mdfa Stages on
   sig_diff = sd; //number of stages stages
   sig_inv = sigi;
   stop_loss = st;
   
  }  
  
  public StrategyParameters(String filt_name, String insamp, String start, double meanR, String end, int nobs, int nr, int _L, int _lag, double freq, double lam, 
                            double exp, double sm, double dec1, double dec2, double cro,
                            int _i1, int _i2, double phase, double hybrid, double hybrid2, int b0, int sd, int sigi, int st, int tp)
  {
   
   meanRevert = meanR;
   startTradingInt = 0; //0-20
   name = filt_name;
   insampStart = insamp;
   startTrading = start; 
   endTrading = end;
   hybrid_weight = hybrid;
   hybrid_weight_diff = hybrid2; 
   time_shift = phase; 
   cutoff0 = 0.0;
   n_rep = nr; lag = _lag; n_obs = nobs; L = _L; 
   cutoff = freq; lambda = lam; expweight = exp; 
   smooth = sm; decay = dec1; decay2 = dec2; cross = cro;
   i1 = _i1; i2 = _i2; b0trend = b0; //--- mdfa Stages on
   sig_diff = sd; //number of stages stages
   sig_inv = sigi;
   stop_loss = st;
   take_profit = tp;

  }    
  
  
  String getName() {return name;}
  
  public String toString()
  {  
    String filtOut = new String(insampStart + " (" + startTrading + "," + endTrading + ") ");
    filtOut = filtOut + "N = " + n_obs + ", L = " + L + ", nrep = " + n_rep;
    filtOut = filtOut + ", Freq: " + cutoff;
    filtOut = filtOut + "(" + smooth + " " + decay + " " + decay2 + " " + cross+") ";
    filtOut = filtOut + " (" + lambda + " " + expweight + ") ";
    filtOut = filtOut + " (" + lag + ", " + i1 + ", " + i2 + ") ";
    filtOut = filtOut + " TS = " + time_shift + ", OS = " + hybrid_weight + ", OSD = " + hybrid_weight_diff + ", b0 = " + b0trend + ", sigdf = " + sig_diff + ", inv = " + sig_inv + ", stop_loss = " + stop_loss + "\n";  
    
    return filtOut;
  }
  
  
  public String toStringText()
  {    
  
   String filtOut = new String(insampStart + " " + startTrading + " " + startTradingInt + " " + endTrading + " " + n_obs + " " + L + " " + n_rep + " " + 
                               cutoff + " " + smooth + " " + decay + " " + decay2 + " " + cross + " " + lambda + " " + expweight
                                + " " + lag  + " " + i1 + " " + i2 + " " + time_shift + " " + hybrid_weight + " " + hybrid_weight_diff + " " + b0trend + " " + sig_diff
                              + sig_inv + " " + stop_loss + " " + take_profit);
                                
   return filtOut;                             
  
  }
  
  
  public String toStringFunction()
  {
    
     String out = new String(" = new strategyParameters(\""+name+"\",\""+insampStart+"\",\""+startTrading+"\","+startTradingInt+",\""+endTrading+"\",\""+
                              n_obs+","+n_rep+","+L+","+lag+","+cutoff+","+lambda+","+expweight+","+smooth+","+decay+","+decay2+","+cross+
                              ","+i1+","+i2+","+time_shift+","+hybrid_weight+","+hybrid_weight_diff+","+b0trend+","+sig_diff+","+sig_inv+","+stop_loss + ");");
                              
     return out;                          
  }
  
  
}
