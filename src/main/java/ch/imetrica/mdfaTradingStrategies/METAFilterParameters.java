package ch.imetrica.mdfaTradingStrategies;


public class METAFilterParameters
{


  //String[] asset_names;
  boolean longshortPos; 
  boolean longPos;
  boolean shortPos;

  String fnString;
  String which_sym;
  int port_size;
  int window_size;
  int meta_filter; //Sharpe, Rank, TradeRatio
  int weights; //Uniform, Startistic, MaxSharpe
  boolean sig_reverse;   //reverse sign 
  boolean equal_dist;    //equal distribution
  
  
  public METAFilterParameters(boolean ls, boolean l, boolean sp, String assets, String fn, int m, int window, int meta, int w, boolean sigr, boolean equal)
  {

    longshortPos = ls; longPos = l; shortPos = sp;
  
    which_sym = new String(assets);
       
    fnString = new String(fn);
    port_size = m; window_size = window; meta_filter = meta; 
    weights = w; sig_reverse = sigr; equal_dist = equal; 
  
  }
  
  public String toStringFunction()
  {
  
    String out = new String(" = new metaFilterParameters("+longshortPos+","+longPos+","+shortPos+",\""+which_sym +"\","+fnString+","+port_size+","+window_size
                            +","+meta_filter+","+weights+","+sig_reverse+","+equal_dist+");");
  
    return out;
  
  }

}  

