package ch.imetrica.mdfaTradingStrategies;

import java.util.ArrayList;


public class FilteredPortfolio
{

  ArrayList<Double[]> returns;
  int n_series; int n_obs;
  
  String filter_name; 
  String[] dates; 
  String[] symbols;

  public FilteredPortfolio(ArrayList<Double[]> r, String[] da, ArrayList<String> nam, String file)
  {
  
   filter_name = new String(file);
   returns = r;
   n_series = r.size();
   if(n_series > 0) n_obs = r.get(0).length;
   else {n_obs = 0;}
   
   symbols = nam.toArray(new String[0]);   
   if(symbols.length != n_series) {System.out.println("Warning: Number of symbols does not match number in portfolio");}
   
   dates = new String[da.length-1];
   for(int i=1; i < da.length; i++)
   {dates[i-1] = da[i];}   
   
   if(n_obs != dates.length)
   {System.out.println("Warning: Number of dates does not match number of trading days: " + n_obs + " vs " + dates.length);} 
  }

  public String getFilterName()
  {return filter_name;}
  
  public Double[] getSeries(int i)
  {  
    if(i < n_series)
    {return returns.get(i);}
    else return null;
  }
  
  public String[] getDates() {return dates;}
  
  public int getNSeries()
  {return n_series;}
  
  public int getNObs()
  {return n_obs;}
  
  public String getSymbolName(int i) {if(symbols.length > 0) {return symbols[i];} else {return "none";}}
  
}