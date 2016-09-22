package ch.imetrica.mdfaTradingStrategies;

import java.io.*;
import java.text.*;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;

public class StrategyFilter 
{

  String fileName;
  
  int n_filters = 0;
  int max_strategies = 0;
  FilteredPortfolio[] port;
  JInvestment[] strategyRets;
  ArrayList<String> names; 
  
  
  //used when combinding strategies to ensure uniquness 
  boolean black_list_on = false;
  int[][] black_list;  
  boolean black_list_set = false;
  int n_blacked,black_start; 
  
  public StrategyFilter(int n) 
  {
    max_strategies = n;
    port = new FilteredPortfolio[n];
    n_filters = 0;
        
  }
  
  public void toggleBlackList(boolean t) {black_list_on = t;}
  
  public void setBlackList(int n, int length, int[][] list)
  {
  
    n_blacked = n;
    black_list = new int[n_blacked][length];
    String[] assts = names.toArray(new String[0]);
  
    for(int k=0;k<n;k++)
    {System.arraycopy(list[k],0,black_list[k],0,length);}
    if(n_blacked > 0) {black_list_set = true;}
    
    for(int i=0; i<(int)(length/6);i++)
    {
     for(int k=0;k<n_blacked;k++)
     {
       System.out.print(assts[black_list[k][length-1-i]] + " ");   
     }
     System.out.println("");
    } 
  }
  
  
  public String[] getSymbols()
  {
    if(port[0].symbols.length > 0) 
    {
      return port[0].symbols;
    }
    else return null; 
  }
  
  public boolean loadFilteredReturns(String file1)
  {
     
     boolean pass = true;
     String[] tokens; String delims = "[ ]+";
     int n_toks; String strline; int i;
     
     names = new ArrayList<String>();
     ArrayList<Double[]> returns_all = new ArrayList<Double[]>();
     String[] dates = new String[10];     
     int count = 0; 
     try
     {
          
           FileInputStream fin = new FileInputStream(new File(file1));
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           names.clear();
           while((strline = br.readLine()) != null)
           {
             
             if(count == 0) //first row are the dates
             {
               dates = strline.split(delims); 
             }  
             else
             {
             
              tokens = strline.split(delims); 
              n_toks = tokens.length; //length should be number of returns + 1
              //System.out.println("n_toks length = " + n_toks);
             
              Double[] returns = new Double[n_toks-1];
              names.add(tokens[0]); //the first value is the asset name
              
              for(i = 1; i < n_toks; i++)
              {returns[i-1] = new Double(tokens[i]);} 
              returns_all.add(returns);   
             }
             count++;
           }
           
           //asset_names.add(names.toArray(new String[0])); 
           
           din.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  

     if(n_filters < max_strategies)
     {
       if(n_filters == 0)  //first one in is okay
       {port[n_filters] = new FilteredPortfolio(returns_all, dates, names, file1); n_filters++;} 
       else if(port[n_filters-1].symbols.length == names.size()) //now check for consistency 
       {
         port[n_filters] = new FilteredPortfolio(returns_all, dates, names, file1);
         n_filters++;
       }
       else {pass = false; System.out.println("Number of symbols must agree for each strategy");
         System.out.println("Portfolio at " + n_filters + " is " + port[n_filters-1].symbols.length + " " + names.size()); 
       }       
     }
     else
     {System.out.println("Maxed out on filtered series"); pass = false;}
     
     return pass; 
     
  }
  
 
  
  /*-----------------------------------------------------------------------
  
  
  
  ------------------------------------------------------------------------*/
  
  
  public JInvestment[] computeStrategyGeneral(int n_basket, int nlags, int filter_method, int normalize, boolean equal, boolean strat_longonly, int nonneg, boolean print,
                                              boolean blackbool, int nblack, int lenghtblack, int[][] black, String name_filter)  
  {
  
    int k,i,j;
    int mindays = 100000; int start; 
    int ndays,nseries,the_one;  
    int M = n_basket;
    double sumw = 0;
    boolean duplicate;
    String meta_stat; 
    String[] dates;
    
    the_one = -1;
    if(filter_method == 0) {meta_stat = "Sharpe";}
    else if(filter_method == 1) {meta_stat = "Rank";}
    else {meta_stat = "TradeRatio";}
    
    if(name_filter.equals("all"))
    {nseries = port[0].getNSeries();}
    else {nseries = 1;}
    
    //System.out.println("nseries = " + nseries + ", name_filter = " + name_filter);
    if(nseries == 1 && !name_filter.equals("all")) //get the one
    {
      //System.out.println("nseries = " + nseries);
      for(i = 0; i < port[0].getNSeries(); i++)
      {
        if(port[0].symbols[i].equals(name_filter))
        {the_one = i; black_list_on = false;  break;}
      }
    }
    else if(nseries == 1 && name_filter.equals("all"))
    {the_one = 0; black_list_on = false;}
    
    if(the_one == -1) //-- revert back to all
    {nseries = port[0].getNSeries();}
    
    
    int[] popoff; int[] heap_index; double[] heap = new double[n_filters];
    int merge_count,best_filter,best_series;
    int min_filter = 0;
  
    //--- first find min ndays
    for(i = 0; i < n_filters; i++)
    {
      if(port[i].getNObs() < mindays) 
      {mindays = port[i].getNObs(); min_filter = i;}         
 /*     if(nseries != port[i].getNSeries())
      {System.out.println("Check number of series in portfolio, seems to not be consistent");}  */        
    }
    
    if(M > nseries)
    {
      //System.out.println("Changing number of portfolio filtered assets to " + nseries);
      M = nseries; n_basket = nseries;
    }
    
    ndays = mindays; 
    dates = port[min_filter].getDates();
    
    double[][] past_returns;
    double portfolio, weightsum; 
    double[][] stats = new double[n_filters][nseries]; 
    double[][] ustats = new double[n_filters][nseries];
    boolean[][] sig_reverse = new boolean[n_filters][nseries];
    
    strategyRets = new JInvestment[ndays-nlags];
    
    int[][] indices = new int[n_filters][nseries];
    Double[][] day_returns = new Double[n_filters][nseries];
    int[][] merged_index = new int[M][2]; //indicates index and strategy
    String[] thegoods = new String[M];
    int[] thegoodsID = new int[M];
    double[] therets = new double[M];
    boolean[] thesigrev = new boolean[M];
    double[] theweights = new double[M];
    int[] thefilters = new int[M];
    double[] thestats = new double[M];    
    
    //where does black list start (make sure black list is at least as big as any tseries

    toggleBlackList(blackbool);
    setBlackList(nblack, lenghtblack, black);    
    

    if(black_list_on) 
    {black_start = black_list[0].length - ndays;}    
    
       
    
    if(black_list_on) System.out.println("black_start = " + black_start + ", " + black_list[0].length + ", n_blacks = " + n_blacked);
    for(i = nlags; i < ndays; i++) //run in real-time 
    {
      for(j = 0; j < n_filters; j++) //for each filtered portfolio
      { 
        for(k=0;k<nseries;k++) //now for each series in the portfolio (these HAVE TO MATCH)
        {       
          Double[] tseries;
          if(nseries == 1)
          {tseries = port[j].returns.get(the_one);}
          else
          {tseries = port[j].returns.get(k);}
          
          start = tseries.length - ndays;
                    
          if(!strat_longonly)
          {
            for(int l=0;l<tseries.length;l++) {tseries[l] = -tseries[l];}
          }
                  
          //if(strat_longonly) 
          //{
           stats[j][k] = computeStatistic(tseries, start + i - nlags, start + i - 1, filter_method);
           if(i < ndays) day_returns[j][k] = tseries[start + i]; sig_reverse[j][k] = false;
          //}
/*          else 
          { 
           stats[j][k] = computeStatistic(tseries, start + i - nlags, start + i - 1,filter_method);
           if(stats[j][k] < 0) 
           {
            if(i < ndays) {day_returns[j][k] = tseries[start + i];} 
            stats[j][k] = Math.abs(stats[j][k]); sig_reverse[j][k] = false;
           }
           else 
           { 
            if(i < ndays) {day_returns[j][k] = -tseries[start + i];}
            sig_reverse[j][k] = true;
           }
          }  */        
        }
        System.arraycopy(stats[j], 0, ustats[j], 0, nseries);
        sort_sims(nseries, stats[j], indices[j]);             
      }
      
      
      if(!equal || M < n_filters)
      {
       //--- now merge the sorted performances --------------
       popoff = new int[n_filters]; 
       heap_index = new int[n_filters];
       merge_count = 0;
       for(j = 0; j < n_filters; j++) {popoff[j] = 0;}
      
       while(merge_count < M)
       {
      
       //scrape the cream off the top layer of the cake yo
        for(j = 0; j < n_filters; j++)
        {heap[j] = stats[j][nseries - 1 - popoff[j]];}
      
        //sort that shit bitch 
        sort_sims(n_filters, heap, heap_index);
      
       //take the cream of the crop, u know u can't stop 
        best_filter = heap_index[n_filters-1];       //get best performer
        best_series = indices[best_filter][nseries - 1 - popoff[best_filter]]; //get best performer's index
        duplicate = false;
        //first make sure not used in any black list
        //System.out.println("Merge count = " + merge_count + " of " + M);
        //System.out.print(assts[best_series] + " in ");
        if(black_list_on)
        {
         for(j=0;j<n_blacked;j++)
         {
          if(black_list[j][black_start + i - 1] == best_series) 
          {duplicate = true;}
         }
         //System.out.println();
        }
        
      
       //--- now make sure no duplicates 
        
        //System.out.print(i);
        if(!duplicate){
        for(k=0;k<merge_count;k++)
        {
          
          if(merged_index[k][1] == best_series)
          {duplicate = true;}
          
//           if(black_list_on && !duplicate) 
//           {
//             System.out.print(assts[best_series] + " in ");
//             for(j=0;j<n_blacked;j++)
//             {
//              System.out.print(assts[black_list[j][black_start + i - 1]] + ", ");
//              if(black_list[j][black_start + i - 1] == best_series) 
//              {duplicate = true; break;}
//             }
//             System.out.println("");
//           }
        }}
        
        if(!duplicate) //if no duplicate, we have a winner 
        {
         merged_index[merge_count][0] = best_filter; 
         merged_index[merge_count][1] = best_series;
         merge_count++;
        }
      
        popoff[best_filter]++;
       } 
      }
      else if(equal && M >= n_filters) //take equal amount of top performers from each strategy
      {

              
        popoff = new int[n_filters]; 
        merge_count = 0;
        for(j = 0; j < n_filters; j++) {popoff[j] = 0;}
      
        best_filter = 0;
        while(merge_count < M)
        {
      
         //scrape the cream off the top layer of the cake yo
//          for(j = 0; j < n_filters; j++)
//          {heap[j] = stats[j][nseries - 1 - popoff[j]];}
      
         //sort that shit bitch 
         //sort_sims(n_filters, heap, heap_index);

         //System.out.println("Merge countequal = " + merge_count + " of " + M);
         best_series = indices[best_filter][nseries - 1 - popoff[best_filter]]; //get best performer's index
         duplicate = false;
         
         //System.out.print(assts[best_series] + " in ");
         if(black_list_on)
         { 
          for(j=0;j<n_blacked;j++)
          {  
           if(black_list[j][black_start + i - 1] == best_series) 
           {duplicate = true; popoff[best_filter]++; break;}
          }
          //System.out.println();
         }
      
      
      
         //--- now make sure no duplicates 
         
         //System.out.print(i);
         if(!duplicate){
         for(k=0;k<merge_count;k++)
         {
          if(merged_index[k][1] == best_series)
          {duplicate = true; popoff[best_filter]++;}

/*          if(black_list_on && !duplicate) 
          {
            System.out.print(assts[best_series] + " in ");
            for(j=0;j<n_blacked;j++)
            {
             System.out.print(assts[black_list[j][black_start + i - 1]] + ", ");
             if(black_list[j][black_start + i - 1] == best_series) 
             {duplicate = true; popoff[best_filter]++; break;}
            }  
            System.out.println("");
          } */         
         }}
      
         if(!duplicate) //if no duplicate, we have a winner 
         {
          merged_index[merge_count][0] = best_filter; 
          merged_index[merge_count][1] = best_series;
          merge_count++; 
          popoff[best_filter]++;
          
          best_filter++;
          if(best_filter >= n_filters) {best_filter = 0;} //reset
         }        
        }      
      
      }
      
  
      if(print)
      {
       for(j=0;j<M;j++)
       {System.out.println(merged_index[j][0] + " " + merged_index[j][1] + " " + ustats[merged_index[j][0]][merged_index[j][1]]);}
       System.out.println();
      }
  

  
  
      if(normalize == 0) //method 1: equal weighting
      {
          //invest in the top M
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];  
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()/M;
             
             thefilters[j] = best_filter;
             theweights[j] = 1.0/M;
             
             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }
             
             
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];
          }
                 
          strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);
          
          //System.out.print("\n");
      }  
      else if(normalize == 1)// method 2: weight according to strength
      {
          weightsum = 0; 
          for(j=0;j<M;j++) 
          {
            best_filter = merged_index[j][0];
            best_series = merged_index[j][1];
            weightsum = weightsum + ustats[best_filter][best_series]; 
          }
          
          portfolio = 0; 
          for(j = 0; j < M; j++)
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()*ustats[best_filter][best_series]/weightsum;
          
             thefilters[j] = best_filter;
             theweights[j] = ustats[best_filter][best_series]/weightsum;
             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];          
          
          
          }// System.out.print(indices[nseries-1-j] + " ");}          
          strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);        
      }
      else //max-sharpe weighting 
      {
          past_returns = new double[nlags][M];
          for(j=0;j<M;j++) //get past returns for the M best assets
          {
            best_filter = merged_index[j][0];
            best_series = merged_index[j][1];  
          
            Double[] tseries = port[best_filter].returns.get(best_series);
            start = tseries.length - ndays;
            
            for(k=0;k<nlags;k++) {past_returns[k][j] = tseries[start + i-nlags+k]; if(print)System.out.println(tseries[start + i-nlags+k]);}
            if(print)System.out.println("\n");
          }
          double[] weights = maximizeSharpe(past_returns, M, nlags, nonneg);
          sumw = 0;
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];        
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()*weights[j]; 
             if(print)System.out.println(day_returns[best_filter][best_series].doubleValue() + ", w = " + weights[j]);
             sumw = sumw + weights[j];
             
             thefilters[j] = best_filter;
             theweights[j] = weights[j];
             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];                     
          }          
         
          strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);       
         
          if(print)System.out.println("\nreturn for day " + i + " is " + strategyRets[i-nlags].toString());
        }
      }
      
      if(print)
      {  
       System.out.println("Best filtered performances:");
       for(j=0;j<M;j++)
       {System.out.println("Filter " + merged_index[j][0] + " at " + merged_index[j][1]);}      
      }
      
      return strategyRets;    
  }
  
  
  
  
  public JInvestment[] computeStrategyPredictive(int n_basket, int nlags, int filter_method, int normalize, boolean equal, boolean strat_longonly, int nonneg, boolean print,
                                              boolean blackbool, int nblack, int lenghtblack, int[][] black, String name_filter)  
  {
  
    int k,i,j;
    int mindays = 100000; int start; 
    int ndays,nseries;  
    int M = n_basket;
    double sumw = 0;
    boolean duplicate;
    String meta_stat; 
    String[] dates;
   
    names.toArray(new String[0]);
   
    int the_one = -1;
    if(filter_method == 0) {meta_stat = "Sharpe";}
    else if(filter_method == 1) {meta_stat = "Rank";}
    else {meta_stat = "TradeRatio";}
    
    if(name_filter.equals("all"))
    {nseries = port[0].getNSeries();}
    else {nseries = 1;}
    
    if(nseries == 1) //get the one
    {
      for(i = 0; i < port[0].getNSeries(); i++)
      {
        if(port[0].symbols[i].equals(name_filter))
        {the_one = i; black_list_on = false; break;}
      }
    }
    
    if(the_one == -1) //-- revert back to all
    {nseries = port[0].getNSeries();}
    
    int[] popoff; int[] heap_index; double[] heap = new double[n_filters];
    int merge_count,best_filter,best_series;
    int min_filter = 0;
  
    //--- first find min ndays
    for(i = 0; i < n_filters; i++)
    {
      if(port[i].getNObs() < mindays) {mindays = port[i].getNObs(); min_filter = i;}
          
      if(nseries != port[i].getNSeries())
      {System.out.println("Check number of series in portfolio, seems to not be consistent");}          
    }
    
    if(M > nseries)
    {
      System.out.println("Changing number of portfolio filtered assets to " + nseries);
      M = nseries; n_basket = nseries;
    }
    
    ndays = mindays; 
    dates = port[min_filter].getDates();
    
    double[][] past_returns;
    double portfolio, weightsum; 
    double[][] stats = new double[n_filters][nseries]; 
    double[][] ustats = new double[n_filters][nseries];
    boolean[][] sig_reverse = new boolean[n_filters][nseries];
    
    strategyRets = new JInvestment[ndays-nlags+1];
    
    int[][] indices = new int[n_filters][nseries];
    Double[][] day_returns = new Double[n_filters][nseries];
    int[][] merged_index = new int[M][2]; //indicates index and strategy
    String[] thegoods = new String[M];
    int[] thegoodsID = new int[M];
    double[] therets = new double[M];
    double[] theweights = new double[M];
    boolean[] thesigrev = new boolean[M];
    int[] thefilters = new int[M];
    double[] thestats = new double[M];    
    
    //where does black list start (make sure black list is at least as big as any tseries

    toggleBlackList(blackbool);
    setBlackList(nblack, lenghtblack, black);    
    

    if(black_list_on && black_list_set) 
    {black_start = black_list[0].length - ndays;}    
    
       
    
  
    for(i = nlags; i <= ndays; i++) //run in real-time 
    {
      for(j = 0; j < n_filters; j++) //for each filtered portfolio
      { 
        for(k=0;k<nseries;k++) //now for each series in the portfolio (these HAVE TO MATCH)
        {       

          Double[] tseries;
          if(nseries == 1)
          {tseries = port[j].returns.get(the_one);}
          else
          {tseries = port[j].returns.get(k);}        
          
          start = tseries.length - ndays;
                  
          if(!strat_longonly)
          {
            for(int l=0;l<tseries.length;l++) {tseries[l] = -tseries[l];}
          }
                  
          //if(strat_longonly) 
          //{
           stats[j][k] = computeStatistic(tseries, start + i - nlags, start + i - 1, filter_method);
           if(i < ndays) day_returns[j][k] = tseries[start + i]; sig_reverse[j][k] = false;
          //}
/*          else 
          { 
           stats[j][k] = computeStatistic(tseries, start + i - nlags, start + i - 1,filter_method);
           if(stats[j][k] < 0) 
           {
            if(i < ndays) {day_returns[j][k] = tseries[start + i];} 
            stats[j][k] = Math.abs(stats[j][k]); sig_reverse[j][k] = false;
           }
           else 
           { 
            if(i < ndays) {day_returns[j][k] = -tseries[start + i];}
            sig_reverse[j][k] = true;
           }
          }  */    
        }
        System.arraycopy(stats[j], 0, ustats[j], 0, nseries);
        sort_sims(nseries, stats[j], indices[j]);             
      }
      
      if(!equal || M < n_filters)
      {
       //--- now merge the sorted performances --------------
       popoff = new int[n_filters]; 
       heap_index = new int[n_filters];
       merge_count = 0;
       for(j = 0; j < n_filters; j++) {popoff[j] = 0;}
      
       while(merge_count < M)
       {
      
       //scrape the cream off the top layer of the cake yo
        for(j = 0; j < n_filters; j++)
        {heap[j] = stats[j][nseries - 1 - popoff[j]];}
      
        //sort that shit bitch 
        sort_sims(n_filters, heap, heap_index);
      
       //take the cream of the crop, u know u can't stop 
        best_filter = heap_index[n_filters-1];       //get best performer
        best_series = indices[best_filter][nseries - 1 - popoff[best_filter]]; //get best performer's index
        duplicate = false;
        //first make sure not used in any black list
        //System.out.println("Merge count = " + merge_count + " of " + M);
        //System.out.print(assts[best_series] + " in ");
        if(black_list_on)
        {
         for(j=0;j<n_blacked;j++)
         {
          if(black_list[j][black_start + i - 1] == best_series) 
          {duplicate = true; break;}
         }
         //System.out.println();
        }
        
      
       //--- now make sure no duplicates 
        
        //System.out.print(i);
        if(!duplicate){
        for(k=0;k<merge_count;k++)
        {
          
          if(merged_index[k][1] == best_series)
          {duplicate = true;}
          
//           if(black_list_on && !duplicate) 
//           {
//             System.out.print(assts[best_series] + " in ");
//             for(j=0;j<n_blacked;j++)
//             {
//              System.out.print(assts[black_list[j][black_start + i - 1]] + ", ");
//              if(black_list[j][black_start + i - 1] == best_series) 
//              {duplicate = true; break;}
//             }
//             System.out.println("");
//           }
        }}
        
        if(!duplicate) //if no duplicate, we have a winner 
        {
         merged_index[merge_count][0] = best_filter; 
         merged_index[merge_count][1] = best_series;
         merge_count++;
        }
      
        popoff[best_filter]++;
       } 
      }
      else if(equal && M >= n_filters) //take equal amount of top performers from each strategy
      {

              
        popoff = new int[n_filters]; 
        merge_count = 0;
        for(j = 0; j < n_filters; j++) {popoff[j] = 0;}
      
        best_filter = 0;
        while(merge_count < M)
        {
      
         //scrape the cream off the top layer of the cake yo
//          for(j = 0; j < n_filters; j++)
//          {heap[j] = stats[j][nseries - 1 - popoff[j]];}
      
         //sort that shit bitch 
         //sort_sims(n_filters, heap, heap_index);

         //System.out.println("Merge countequal = " + merge_count + " of " + M);
         best_series = indices[best_filter][nseries - 1 - popoff[best_filter]]; //get best performer's index
         duplicate = false;
         
         //System.out.print(assts[best_series] + " in ");
         if(black_list_on)
         { 
          for(j=0;j<n_blacked;j++)
          {  
           if(black_list[j][black_start + i - 1] == best_series) 
           {duplicate = true; popoff[best_filter]++; break;}
          }
          //System.out.println();
         }
      
      
      
         //--- now make sure no duplicates 
         
         //System.out.print(i);
         if(!duplicate){
         for(k=0;k<merge_count;k++)
         {
          if(merged_index[k][1] == best_series)
          {duplicate = true; popoff[best_filter]++;}

/*          if(black_list_on && !duplicate) 
          {
            System.out.print(assts[best_series] + " in ");
            for(j=0;j<n_blacked;j++)
            {
             System.out.print(assts[black_list[j][black_start + i - 1]] + ", ");
             if(black_list[j][black_start + i - 1] == best_series) 
             {duplicate = true; popoff[best_filter]++; break;}
            }  
            System.out.println("");
          } */         
         }}
      
         if(!duplicate) //if no duplicate, we have a winner 
         {
          merged_index[merge_count][0] = best_filter; 
          merged_index[merge_count][1] = best_series;
          merge_count++; 
          popoff[best_filter]++;
          
          best_filter++;
          if(best_filter >= n_filters) {best_filter = 0;} //reset
         }        
        }      
      
      }
      
      if(print)
      {
       for(j=0;j<M;j++)
       {System.out.println(merged_index[j][0] + " " + merged_index[j][1] + " " + ustats[merged_index[j][0]][merged_index[j][1]]);}
       System.out.println();
      }
  

  
  
      if(normalize == 0) //method 1: equal weighting
      {
          //invest in the top M
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()/M;
             
             thefilters[j] = best_filter;
             theweights[j] = 1.0/M;

             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }             
             
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];
          }
                 
          if(i < ndays) strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);
          
          //System.out.print("\n");
      }  
      else if(normalize == 1)// method 2: weight according to strength
      {
          weightsum = 0; 
          for(j=0;j<M;j++) 
          {
            best_filter = merged_index[j][0];
            best_series = merged_index[j][1];
            weightsum = weightsum + ustats[best_filter][best_series]; 
          }
          
          portfolio = 0; 
          for(j = 0; j < M; j++)
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()*ustats[best_filter][best_series]/weightsum;
          
             thefilters[j] = best_filter;
             theweights[j] = ustats[best_filter][best_series]/weightsum;

             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }             
             
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];          
          
          
          }// System.out.print(indices[nseries-1-j] + " ");}          
          if(i < ndays) strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);        
      }
      else //max-sharpe weighting 
      {
          past_returns = new double[nlags][M];
          for(j=0;j<M;j++) //get past returns for the M best assets
          {
            best_filter = merged_index[j][0];
            best_series = merged_index[j][1];  
          
            Double[] tseries = port[best_filter].returns.get(best_series);
            start = tseries.length - ndays;
            
            for(k=0;k<nlags;k++) {past_returns[k][j] = tseries[start + i-nlags+k]; if(print)System.out.println(tseries[start + i-nlags+k]);}
            if(print)System.out.println("\n");
          }
          double[] weights = maximizeSharpe(past_returns, M, nlags, nonneg);
          sumw = 0;
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {
             best_filter = merged_index[j][0];
             best_series = merged_index[j][1];        
             portfolio = portfolio + day_returns[best_filter][best_series].doubleValue()*weights[j]; 
             if(print)System.out.println(day_returns[best_filter][best_series].doubleValue() + ", w = " + weights[j]);
             sumw = sumw + weights[j];
             
             thefilters[j] = best_filter;
             theweights[j] = weights[j];

             if(the_one >= 0)
             {
              thegoods[j] = port[best_filter].getSymbolName(the_one);  thegoodsID[j] = the_one;
             }
             else
             {
              thegoods[j] = port[best_filter].getSymbolName(best_series);  thegoodsID[j] = best_series;
             }              
             
             therets[j] = day_returns[best_filter][best_series].doubleValue();
             thesigrev[j] = sig_reverse[best_filter][best_series];
             thestats[j] = ustats[best_filter][best_series];                     
          }          
         
          if(i < ndays) strategyRets[i-nlags] = new JInvestment(dates[i], thegoods, thegoodsID, thefilters, meta_stat, therets, 
                                                    theweights, thestats, thesigrev, portfolio);       
         
          if(i < ndays) {if(print)System.out.println("\nreturn for day " + i + " is " + strategyRets[i-nlags].toString());}
        }
        
        if(i==ndays) //the predicitive state
        {
          //therets = new double[M]; //we don't know then of course
          strategyRets[i-nlags] = new JInvestment("today", thegoods, thegoodsID, thefilters, meta_stat, therets, theweights, thestats, thesigrev, portfolio); 
        
        }
      }
      
      if(print)
      {  
       System.out.println("Best filtered performances:");
       for(j=0;j<M;j++)
       {System.out.println("Filter " + merged_index[j][0] + " at " + merged_index[j][1]);}      
      }
      
      return strategyRets;    
  }    
  
  
 
  
  
  
  
  
  
  public static double[] computeStrategyUnivariate(ArrayList<Double[]> returns, int n_basket, int nlags, int filter_method, 
                                                   int normalize, boolean strat_longonly, int nonneg, boolean print)
  {
    int k,i,j;
    int ndays;  
    int nseries = returns.size(); int M = n_basket;
    double sumw = 0;
    ndays = returns.get(0).length;    
        
    double[][] past_returns;
    double portfolio, weightsum; 
    double[] stats = new double[nseries]; 
    double[] ustats = new double[nseries];
    double[] strategyRets = new double[ndays-nlags];
    int[] indices = new int[nseries];
    Double[] day_returns = new Double[nseries];
    
    if(nseries > 0)
    {       
      for(i = nlags; i < ndays; i++)
      {
        if(print)System.out.println("Day " + i);
        for(k=0;k<nseries;k++)
        {       
          Double[] tseries = returns.get(k);
          if(strat_longonly) 
          {
           stats[k] = computeStatistic(tseries,i-nlags,i-1,filter_method);
           day_returns[k] = tseries[i];
          }
          else 
          { 
           stats[k] = computeStatistic(tseries,i-nlags,i-1,filter_method);
           if(stats[k] < 0) {day_returns[k] = -tseries[i]; stats[k] = Math.abs(stats[k]);}
           else {day_returns[k] = tseries[i];}
          }      
        }
        //--- sort but keep sorted array as well--------
        System.arraycopy(stats, 0, ustats, 0, nseries);
        sort_sims(nseries, stats, indices);
          
          
          
          
        // now compute the weighting strategy for the given day   
        if(normalize == 0) //method 1: equal weighting
        {
          //invest in the top M
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {portfolio = portfolio + day_returns[indices[nseries-1-j]].doubleValue()/M;}// System.out.print(indices[nseries-1-j] + " ");}          
          strategyRets[i-nlags] = portfolio;
          //System.out.print("\n");
        }  
        else if(normalize == 1)// method 2: weight according to strength
        {
          weightsum = 0; 
          for(j=0;j<M;j++) {weightsum = weightsum + stats[nseries-1-j];}
          
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {portfolio = portfolio + day_returns[indices[nseries-1-j]].doubleValue()*stats[nseries-1-j]/weightsum;}          
          strategyRets[i-nlags] = portfolio;          
        }
        else //max-sharpe weighting 
        {
          past_returns = new double[nlags][M];
          for(j=0;j<M;j++) //get past returns for the M best assets
          {
            Double[] tseries = returns.get(indices[nseries-1-j]);
            for(k=0;k<nlags;k++) {past_returns[k][j] = tseries[i-nlags+k]; if(print)System.out.println(tseries[i-nlags+k]);}
            if(print)System.out.println("\n");
          }
          double[] weights = maximizeSharpe(past_returns, M, nlags, nonneg);
          sumw = 0;
          portfolio = 0; 
          for(j = 0; j < M; j++) 
          {
            portfolio = portfolio + day_returns[indices[nseries-1-j]].doubleValue()*weights[j]; 
            if(print)System.out.println(day_returns[indices[nseries-1-j]].doubleValue() + ", w = " + weights[j]);
            sumw = sumw + weights[j];
          }          
          strategyRets[i-nlags] = portfolio;     
          if(print)System.out.println("\nreturn for day " + i + " is " + strategyRets[i-nlags] + ", sum of weights is " + sumw);
        }
      }
    }  
    
    return strategyRets;
  }
  
  
  
  
  
  
  
  public static double computeStatistic(Double[] rets, int start, int end, int method)
  {
  
    //System.out.println("End = " + end + " start = " + start);
    int i;
    double statistic = 0; 
    double[] segment = new double[end-start+1];
    double[] cum_segment = new double[end-start+1];
    int n_pos  = 0; 
    cum_segment[0] = 0;
    for(i=start;i<=end;i++)
    {
     segment[i-start] = rets[i]; 
     if(i-start > 0) {cum_segment[i-start] = cum_segment[i-start-1] + rets[i];}
     if(rets[i] > 0) {n_pos++;}
    }
    
    if(method == 0) //sharpe ratio 
    {
      double[] mstd = mean_std(segment); 
      statistic = Math.sqrt(250)*mstd[0]/mstd[1];
    }
    else if(method == 1) //rank correlation 
    {
      statistic = rankCoefficientSeg(cum_segment,cum_segment.length);  
    }
    else if(method == 2) //trade-success ratio 
    {
      statistic = (double)n_pos/(end-start+1) - .5;    
    }
    
    return statistic; 
  }
      

  public static double rankCoefficientSeg(double[] _account, int n)
  {

   int i; 
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double[] account = new double[n];
   System.arraycopy(_account,0,account,0,n);
   double sum = 0.0;
   double spear;
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i;} 
      
   sort_sims(n, account, index);
   
   for (i=0 ; i<n ; ++i)
   {d[i] = Math.abs(index[i] - rank[i]); d[i] = d[i]*d[i]; sum = sum + 1.0*d[i];} 
   
   spear = 1.0 - (6.0*sum)/(1.0*n*(n*n-1.0));

   return spear; 
  }  

  public static double rankCoefficient(double[] _account, int n)
  {

   int i; 
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double[] account = new double[n];
   
   double[] cumsum = new double[n];
   cumsum[0] = _account[0];
   for(i=1; i < n;i++)
   {
     cumsum[i] = cumsum[i-1] + _account[i];
   }
   
   System.arraycopy(_account,0,account,0,n);
   double sum = 0.0;
   double spear;
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i;} 
      
   sort_sims(n, cumsum, index);
   
   
   
   for (i=0 ; i<n ; ++i)
   {d[i] = Math.abs(index[i] - rank[i]); d[i] = d[i]*d[i]; sum = sum + 1.0*d[i];} 
   
   spear = 1.0 - (6.0*sum)/(1.0*n*(n*n-1.0));

   return spear; 
 }  

 public static double blakelyRatio(int n_days, double[] rets)
  {
      int i; int n = rets.length; 
      double[] ranks = new double[n - n_days]; 
      double ratio = 0;
      double mean_rank = 0;
      for(i=0;i<n-n_days;i++)
      {
        double[] section = new double[n_days];
        System.arraycopy(rets,i,section,0,n_days);
      
        ranks[i] = rankCoefficientSeg(section,n_days);
        mean_rank = mean_rank + ranks[i];     
      }
      mean_rank = mean_rank/(n - n_days);  
      
      double sum = 0; 
      for (i=0; i<n-n_days; i++) 
      { 
       final double v = ranks[i] - mean_rank;  sum += v * v; 
      } 
      
      ratio = Math.sqrt(250)*mean_rank/Math.sqrt(sum/(n-n_days));      
  
      return ratio;   
  }
 
 
 public static void sort_sims(int n, double[] object, int[] indices)
 {
    int l,j,ir,i;
    double tempobj;
    int tempindex=0;

    if (n<2) return;		// doesn't work if only one observation
    // but we don't need to sort then anyway

    if (indices!=null)
      for (i=0 ; i<n ; i++)
	indices[i] = i;     // start off with identity permutation

    l = (n >> 1)+1;
    ir = n;

    for (;;) {
      if (l>1)		// i.e. we are still in hiring phase
	{ 
	  --l;
	  if (indices!=null)
	    tempindex = indices[l-1];
	  tempobj = object[l-1];
	}
      else {
	tempobj=object[ir-1];
	if (indices!=null)
	  tempindex=indices[ir-1];
	object[ir-1]=object[0];
	if (indices!=null)
	  indices[ir-1]=indices[0];
	if (--ir == 1) {	// we are finished
	  object[0]=tempobj;
	  if (indices!=null)
	    indices[0]=tempindex;
	  return;
	}
      }
      i=l;		
      j=l<<1;
      while (j <= ir)
	{
	  if ((j < ir) && (object[j-1]<object[j])) 
	    ++j;	// j is better underling
	  if (tempobj<object[j-1]) {		// demote tempobj
	    object[i-1]=object[j-1];
	    if (indices!=null)
	      indices[i-1]=indices[j-1];
	    j += (i=j);
	  }
	  else j=ir+1;			// finished sifting
	}
      object[i-1]=tempobj;
      if (indices!=null)
	indices[i-1]=tempindex;
    }
  }  
  
  public static double[] cumsum(double[] data, int n)
  {

      double[] cs = new double[n]; double sum; int k;
      sum=0; 
      
      for(k=0;k<n;k++)
      {
        sum = sum+data[k]; cs[k] = sum; 
      }
      return cs;  
  } 
  
  
  public static double[] mean_std( Double[] data ) 
  { 

       double mean = 0; 
       final int n = data.length; 
       
       for ( int i=0; i<n; i++ )  {  mean += data[i]; } 
       mean /= n; 

       double sum = 0; 
       for ( int i=0; i<n; i++ ) { final double v = data[i] - mean;  sum += v * v; 
       } 
       double[] co = {mean,  Math.sqrt( sum / n )};
       return co; 
  }  
 
  public static double[] mean_std( double[] data ) 
  { 

       double mean = 0; 
       final int n = data.length; 
       
       for ( int i=0; i<n; i++ )  {  mean += data[i]; } 
       mean /= n; 

       double sum = 0; 
       for ( int i=0; i<n; i++ ) { final double v = data[i] - mean;  sum += v * v; 
       } 
       double[] co = {mean,  Math.sqrt( sum / n )};
       return co; 
  } 
 
 
   //---- returns are given in each column double[][] data = new double[n_obs][n_basket];
   public static double[] maximizeSharpe(double[][] data, int n_basket, int nobs, int nonneg)
   {
   
      int i,j; 
      double[] means = new double[n_basket];
      double sum=0;
      RealVector sol; 
      double[] w = new double[n_basket]; 
       
      for(i=0;i<n_basket;i++)
      {
       sum=0;
       for(j=0;j<nobs;j++)
       {
        sum = sum + data[j][i];        
       }
       means[i] = sum/nobs;
       //System.out.println(means[i]);
      } 
       
      RealVector m = new ArrayRealVector(means, false);
      Covariance covComp = new Covariance(data);
       
      //LinearConstraint(double[] coefficients, Relationship relationship, double value)
      //public static final Relationship LEQ 
      //RealMatrix rm = covComp.scalarMultiply(10000);
      RealMatrix rm = covComp.getCovarianceMatrix();
//       rm = rm.scalarMultiply(1000000);
//       for(i=0;i<n_basket;i++)
//       {printRow(rm.getRow(i));}
      
       
      try
      { 
        DecompositionSolver solver = new QRDecomposition(rm).getSolver();
        sol = solver.solve(m);
        w = sol.toArray(); 
      }
      catch(SingularMatrixException sme) 
      {
       //System.out.println("Matrix singular: setting weights to uniform"); 
       w = new double[n_basket]; 
       for(i=0;i<n_basket;i++) {w[i] = 1.0/n_basket;}
      }
      
      double sumw = 0;
      for(i=0;i<w.length;i++) 
      {
        if(nonneg == 1)
        {if(w[i] < 0) {w[i] = 1.0/n_basket;}}
        else if(nonneg == 2)
        {w[i] = Math.abs(w[i]);}
        
        sumw = sumw + w[i]; 
      }
       
      for(i=0;i<w.length;i++) {w[i] = w[i]/sumw;}
    
      return w; 
  } 
 

  public static void printRow(double[] r)
  {
    DecimalFormat df2 = new DecimalFormat("##.#########");
    for(int i = 0; i < r.length; i++)
    {System.out.print(df2.format(r[i]) + " ");}
    System.out.print("\n");
  }
 
 
  public static void main(String[] args) throws InterruptedException, ExecutionException 
  {
  
    StrategyFilter sf = new StrategyFilter(3);    
    //strategyFilter.computePortfolioStrategy("./longfilter1.dat"); 
    sf.loadFilteredReturns("hfilteri2m2.dat");
    sf.loadFilteredReturns("longhfilteri2m2.dat");
    sf.loadFilteredReturns("shorthfilteri2m2.dat");
    
    //sf.runStrategy();
    
  }
  
  
}  