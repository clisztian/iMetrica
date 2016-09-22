package ch.imetrica.mdfaTradingStrategies;

import java.text.DecimalFormat;

//class that holds particular investment strategy for given day

public class JInvestment
{

  String date;              //the date for the given trading day
  String[] assets;          //the assets traded for given day 
  int[] assetsID;         //corresponding number
  String meta_filter;       //Sharpe, rank, trade_ratio, dfa_filter

  boolean[] sig_reverse;     //does the sign change for the signal
  double[] sign; 
  int[] strategy_ind;
  double[] weights;         //the portfolio weights 
  double[] ind_returns;     //individual returns for the day
  double[] metavals; //meta_filter value (reason why asset was chose for that day
  double totreturn;            //the return for the day
  DecimalFormat df = new DecimalFormat("##.####");
  DecimalFormat df4 = new DecimalFormat("##.#####");

  public JInvestment(String d, String[] ass, int[] ids, int[] strat, String meta, double[] rets, double[] w, double[] met, boolean[] sig, double ret)
  {
    date = new String(d); 
    meta_filter = new String(meta);
    
    assets = new String[ass.length];
    System.arraycopy(ass,0,assets,0,ass.length);
    
    assetsID = new int[ids.length];
    System.arraycopy(ids,0,assetsID,0,ids.length);    
    
    strategy_ind = new int[strat.length];
    System.arraycopy(strat,0,strategy_ind,0,strategy_ind.length);
    
    weights = new double[w.length];
    System.arraycopy(w,0,weights,0,w.length);
    
    ind_returns = new double[rets.length];
    System.arraycopy(rets,0,ind_returns,0,ind_returns.length);
    
    metavals = new double[met.length];
    System.arraycopy(met,0,metavals,0,metavals.length);    
    
    sig_reverse = new boolean[sig.length];
    System.arraycopy(sig,0,sig_reverse,0,sig_reverse.length); 
    
    sign = new double[sig.length];
    for(int i = 0; i < sig.length; i++)
    {
     sign[i] = 1.0; 
     if(sig_reverse[i]) {sign[i] = -1.0;}
    }
    
    totreturn = ret;
 
  } 

  public double getReturn() {return totreturn;}
  
  public String getDate() {return date;}
  
  public String toString()
  {
    int i; 
    String invest = new String(date + ": " + df.format(totreturn) + "\n");
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String(assets[i] + "(" + df.format(sign[i]*weights[i]) + "*" + df.format(ind_returns[i]) + "), ");} 
    invest = invest + new String(assets[assets.length-1] + "(" + df.format(sign[i]*weights[assets.length-1]) + "*" + df.format(ind_returns[assets.length-1]) + ")\n");

    invest = invest + " " + meta_filter + ": ";
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String("(" + strategy_ind[i] + "," + df.format(metavals[i]) + "), ");}    
    invest = invest + new String("(" + strategy_ind[assets.length-1] + "," + df.format(metavals[assets.length-1]) + ")\n");
    
    return invest; 
  }


  public String toStringDateReturn()
  {
    String invest = new String(date + ": " + df.format(totreturn) + "\n");  
    return invest; 
  }
  
  public String toStringAsset()
  {
    int i; 
    String invest = new String("");
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String(assets[i] + "(" + df.format(sign[i]*weights[i]) + "*" + df.format(ind_returns[i]) + "), ");} 
    invest = invest + new String(assets[assets.length-1] + "(" + df.format(sign[i]*weights[assets.length-1]) + "*" + df.format(ind_returns[assets.length-1]) + ")\n");

    invest = invest  + " " +  meta_filter + ": ";
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String("(" + strategy_ind[i] + "," + df.format(metavals[i]) + "), ");}    
    invest = invest + new String("(" + strategy_ind[assets.length-1] + "," + df.format(metavals[assets.length-1]) + ")\n");
    
    return invest; 
  }

  public String toStringAssetNo()
  {
    int i; 
    String invest = new String(" ");
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String(assets[i] + "(" + df.format(sign[i]*weights[i]) + "*" + df.format(ind_returns[i]) + "), ");} 
    invest = invest + new String(assets[assets.length-1] + "(" + df.format(sign[i]*weights[assets.length-1]) + "*" + df.format(ind_returns[assets.length-1]) + ")");

    return invest; 
  }  
  
  
  public String toStringMeta()
  {
    int i; 
    String invest = new String("");

    invest = invest + meta_filter + ": ";
    for(i = 0; i < assets.length-1; i++)
    {invest = invest + new String("(" + strategy_ind[i] + "," + df.format(metavals[i]) + "), ");}    
    invest = invest + new String("(" + strategy_ind[assets.length-1] + "," + df.format(metavals[assets.length-1]) + ")\n");
    
    return invest; 
  }
  
  
  
}