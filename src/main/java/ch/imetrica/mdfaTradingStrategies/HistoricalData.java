package ch.imetrica.mdfaTradingStrategies;

import java.io.*;
import java.net.*;
import java.util.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import org.joda.time.DateTime;
import org.joda.time.Period; 
// import com.ib.client.Contract;
// import com.ib.client.ContractDetails;
// import com.ib.client.EClientSocket;
// import com.ib.client.EWrapper;
// import com.ib.client.EWrapperMsgGenerator;
// import com.ib.client.Execution;
// import com.ib.client.Order;
// import com.ib.client.OrderState;
// import com.ib.client.UnderComp;
// import com.ib.client.Util;
// import com.ib.client.CommissionReport;
// import com.ib.client.MarketDataType;
// import com.ib.client.TickType;
// import com.ib.client.TagValue; 
 
public class HistoricalData
{


  String portfolioParamFile;
  String file = "lastSeriesData.dat";
  String targetSymbol;   //the target asset being traded
  String[] expSymbols;   //the explanatory assets

  boolean ln_trans = true;
  int n_holidays;
  BufferedReader in, sin;
  BufferedWriter sout;      
  int obs_freq;
  //---- IQFeed stuff
  String sSymbol, sInterval, sBeginDateTime, sEndDateTime, sBeginFilterTime, sEndFilterTime;
  String sDataDirection, sRequestID, sDatapointsPerSend;
  String sMaxDatapoints = "";
  String sDays;
  String marketOpenTime;
  
  BufferedWriter[] s_outs;
  BufferedReader[] s_ins;   
  boolean daily = false; 
  String dailyTime;
  
  String currentTimeString;
  int timeID;
  long currentTimeID;
  //private EClientSocket client = null; 
  
  boolean zero_marketopen = false;
  int i1,i2,i3;
  int n_obs,n_rep;
  boolean gaussianize_data;
  boolean diff_vol;
  boolean finished_data;
 
  double[] lastSeries;
  double[] price;
  String[] dates;
 
  int n_assets; 
  int[] ishift;
  
  boolean dailyData;
  
  String[] holidays;
  boolean est = false;
  ArrayList<Double> close_series;
  ArrayList<Double> ask_series;
  ArrayList<Double> bid_series;
  ArrayList<Double> exp_series_1;
  ArrayList<Double> exp_series_2;
  ArrayList<Double> exp_series_3;
  ArrayList<String> dates_series1;
  ArrayList<String> dates_series2;
  ArrayList<String> dates_series3;
  ArrayList<String> dates_series;
 
  ArrayList<String> hilo_series;
  ArrayList<String> hi_series;
  ArrayList<String> lo_series;
  ArrayList<Double> asset_series;
  ArrayList<String> asset_dates;
  ArrayList<Double> asset_series1;
  ArrayList<String> asset_dates1;
  ArrayList<Double> asset_series2;
  ArrayList<String> asset_dates2;  
  HashMap<String,Double> spreads;
 
 
  boolean min_5 = false;
  boolean min_30 = false;
  boolean min_60 = false;
  double[] rt;
  double[] vol;
  double[] vold;
  
  double phi1,phi2,mu;
  double[] sigma; 
  double sigma_t,sigma_t1;
  
  String strline; String[] tokens; String delims = "[,]"; int n_toks;  
  
  static final SimpleDateFormat dateF = new SimpleDateFormat("yyyyMMdd HH:mm:ss");  
  
  
  
  
  
  
  
  
  public HistoricalData()
  {
  
    //file = new String(_file);
    ask_series = new ArrayList<Double>();
    bid_series = new ArrayList<Double>();
    close_series = new ArrayList<Double>();
    exp_series_1 = new ArrayList<Double>();
    exp_series_2 = new ArrayList<Double>();
    exp_series_3 = new ArrayList<Double>();
    dates_series1 = new ArrayList<String>();
    dates_series2 = new ArrayList<String>();
    dates_series3 = new ArrayList<String>();
    
     
    
  }
  

  /*
    void onConnect() 
    {
    
        client = new EClientSocket (this);
        if(!client.isConnected())
        {client.eConnect ("localhost", 7496, 3);}    
        if (client.isConnected())
        {
          System.out.println("not yet connected to IQ");
          System.out.println("Connected to TWS server");          
        }
        //System.out.println("connect to IQ"); //tws_panel.dialogPrintln("connect to IQ");

        
        System.out.println(this.getState());
        
    }

    void onDisconnect() {
      if(client != null && client.isConnected()) 
      {client.eDisconnect();}
      else
      {System.out.println("Not currently connected to TWS");} //tws_panel.dialogPrintln("Not currently connected to TWS");}
    }  
  
  */
  
  
  
  
  
  public void setLogTrans(boolean t) {ln_trans = t;}        
    
    
  public void setupHistoricalDownload()
  {
   
        
         portfolioParamFile = "Portfolio.port";
         int i; int M;
         String strline; 
         String[] data_symbols = null;
         obs_freq = 5;
         gaussianize_data = false;
         Integer I;

         String[] tokens; String delims = "[ ]+";
         int n_toks;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(portfolioParamFile);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the data symbols --------------
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           data_symbols = new String[n_toks]; 
           for(i=0;i<n_toks;i++) {data_symbols[i] = tokens[i];}
           
           //--- get start - end time/date
           strline = br.readLine(); 
           sBeginDateTime = strline; 
           
           strline = br.readLine(); 
           sEndDateTime = strline;            
                         
           //--- get hours           
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length;
           sBeginFilterTime = tokens[0]; 
           sEndFilterTime = tokens[1]; 
                   
           //--- get market open
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length;
           marketOpenTime = tokens[0];

           //---observation frequency (in minutes)          
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length;
           I = new Integer(tokens[0]); obs_freq = I.intValue();           
            
           br.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}
      
        
        M = data_symbols.length;
        expSymbols = new String[M];
        for(i=0;i<M;i++)
        {expSymbols[i] = data_symbols[i];}
          
        //---- set number of days
        sInterval = ""+(obs_freq*60);
        n_rep = M; 
        n_assets = M;

        
          
          
          //--- print it all out --------------------------
          System.out.println("New Portfolio uploaded:");
          System.out.print("Target and explanatory symbols: ");
          for(i=0;i<M;i++) {System.out.print(data_symbols[i] + " ");}        
          System.out.println("Observation frequency and number of historical days: " + sInterval);
          System.out.println("Number of explanatory variables: " + n_rep);
          System.out.println("Beginning time for historical data is " + sBeginFilterTime);
          System.out.println("Ending time for historical data is " + sEndFilterTime);
          System.out.println("Number of observations is " + n_obs);

          //generateTimesForex_60();
          //generateTimes24();
//           if((new Integer(sEndFilterTime).intValue()) > 161000)
//           {

           if(min_5)
           {generateTimes5();}
           else if(!min_30 && !min_60)
           {
            if(!est)
            {generateTimesForex_30();}
            else
            {generateTimesForex_30_EST();}             
           }
           else if(min_30 && min_60)
           {generateTimesForex_60();}
           else if(min_60 && !min_30)
           {
            generateTimesForex_60();
            //generateTimes24();
           }
           else
           {generateTimesForex_60();}
           //generateTimes24();}
//           }
//           else
//           {generateTimes();}   
       
        
       
       
 }
                
 
 
  public void setupHistoricalDownload(String asset)
  {
   
        
         //portfolioParamFile = "dailyPortEquity.port";
         portfolioParamFile = "dailyPort.port";
         int i; int M;
         String strline; 
         String[] data_symbols = null;
         obs_freq = 5;
         gaussianize_data = false;
         Integer I;
 
         min_30 = true; min_60 = false; min_5 = false;
         String[] tokens; String delims = "[ ]+";
         try
         {
          
           FileInputStream fin = new FileInputStream(portfolioParamFile);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the data symbols --------------
           strline = br.readLine(); tokens = strline.split(delims); 

           data_symbols = new String[1];
           data_symbols[0] = asset;
           //--- get start - end time/date
           strline = br.readLine(); 
           sBeginDateTime = strline; 
           
           strline = br.readLine(); 
           sEndDateTime = strline;      
                         
           //--- get hours           
           strline = br.readLine(); tokens = strline.split(delims); 
           sBeginFilterTime = tokens[0]; 
           sEndFilterTime = tokens[1]; 
                   
           //--- get market open
           strline = br.readLine(); tokens = strline.split(delims); 
           marketOpenTime = tokens[0];

           //---observation frequency (in minutes)          
           strline = br.readLine(); tokens = strline.split(delims); 
           I = new Integer(tokens[0]); obs_freq = I.intValue();           
              
           br.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}
      
        
        M = data_symbols.length;
        expSymbols = new String[M];
        for(i=0;i<M;i++)
        {expSymbols[i] = data_symbols[i];}
          
        //---- set number of days
        sInterval = ""+(obs_freq*60);
        n_rep = M; 
        n_assets = M;

        
          
          
          //--- print it all out --------------------------
          System.out.println("New Portfolio uploaded:");
          System.out.print("Target and explanatory symbols: ");
          for(i=0;i<M;i++) {System.out.print(data_symbols[i] + " ");}        
          System.out.println("Observation frequency and number of historical days: " + sInterval);
          System.out.println("Number of explanatory variables: " + n_rep);
          System.out.println("Beginning time for historical data is " + sBeginFilterTime);
          System.out.println("Ending time for historical data is " + sEndFilterTime);
          System.out.println("Number of observations is " + n_obs);

          //generateTimesForex_60();
          //generateTimes24();
//           if((new Integer(sEndFilterTime).intValue()) > 161000)
//           {

           if(min_5)
           {generateTimes5();}
           else if(!min_30 && !min_60)
           {
            if(!est)
            {generateTimesForex_30();}
            else
            {generateTimesForex_30_EST();}             
           }
           else if(min_30 && min_60)
           {generateTimesForex_60();}
           else if(min_60 && !min_30)
           {
            generateTimesForex_60();
            //generateTimes24();
           }
           else
           {generateTimesForex_60();}
           //generateTimes24();}
//           }
//           else
//           {generateTimes();}   
       
        
       
       
 } 
 
 
    

  public void downloadHistoricalIQData(String name)
  {
      
      int i,j,m; 
      double[] price;
      String dateStamp;
      int shift = 1;
      DateTime weekend; 
      String[] tokens; String ddelims = "[-]";   
      String[] tokens1; String delims = "[ ]";
      //n_rep = 1;
      
      String[] names = name.split(delims);
      n_rep = names.length;
      
      System.out.println(n_rep + " " + name);
      
      startIQConnect IQ = new startIQConnect();
      IQ.run();     
   

   
   
      asset_series = new ArrayList<Double>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      
      for(m = 0; m < n_rep; m++)
      {
       try {getExplanatoryData(names,m);} catch (Exception e){};
      }
      try {Thread.sleep (300);} catch (Exception e) {};
      
       
      lastSeries = new double[n_rep*n_obs];       
      price = new double[n_obs];     
   
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      String[] dates = new String[n_obs+1];
      
      System.out.println(asset_dates.size() + " " + asset_dates1.size() + " " + asset_dates2.size());
      int min_size = asset_dates.size() - 6;
      if((asset_dates1.size() - 6) <  min_size && n_rep > 1) {min_size = asset_dates1.size() - 6;}
      if((asset_dates2.size() - 6) <  min_size && n_rep > 1) {min_size = asset_dates2.size() - 6;}
      
      
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      
      while(i<=n_obs && min_size >= i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         
         if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
         //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
         
         else
         {
         dateStamp = dates_series.get(dates_series.size() - 1 - i);
         //System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
         //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
                
         //--- dates_series1 -----
         if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
         {
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
         {
           i1++;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
         {
           i1=i1+2;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 4)))
         {
           i1=i1+3;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 5)))
         {
           i1=i1+4;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }         
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 6)))
         {
           i1=i1+5;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }             
         
         else
         {
          close_sp[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1);
          //System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
         }
        
         if(n_rep > 1)
         {
//           if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - shift - i2)))
//           {
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2);
//            i2++;
//           }
          
         if(dateStamp.equals(asset_dates1.get(asset_dates1.size() -  i2 - 1)))
         {
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);
           i2++;
         }
         else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 2)))
         {
           i2++;
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
           i2++;
         }
         else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 3)))
         {
           i2=i2+2;
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
           i2++;
         }
         else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 4)))
         {
           i2=i2+3;
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
           i2++;
         }
         else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 5)))
         {
           i2=i2+4;
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
           i2++;
         }         
         else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 6)))
         {
           i2=i2+5;
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
           i2++;
         }  
          

          else
          {
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2 - 1);
           //System.out.println("Missing data at " + dateStamp + " for 2nd explanatory series");
          }
         
         
//           if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - shift - i3)))
//           {
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3);
//            i3++;
//           }
//           
          
          
         if(dateStamp.equals(asset_dates2.get(asset_dates2.size() -  i3 - 1)))
         {
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);
           i3++;
         }
         else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 2)))
         {
           i3++;
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
           i3++;
         }
         else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 3)))
         {
           i3=i3+2;
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
           i3++;
         }
         else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 4)))
         {
           i3=i3+3;
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
           i3++;
         }
         else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 5)))
         {
           i3=i3+4;
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
           i3++;
         }         
         else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 6)))
         {
           i3=i3+5;
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
           i3++;
         }            
          
          
          else
          {
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3 - 1); 
           //System.out.println("Missing data at " + dateStamp + " for 3rd explanatory series");    
          }        
         }
        
        
         dates[n_obs-i] = dateStamp; 
         i++;
        } 
      }
          
      //--- print to file
      file =  names[0] +".dat";
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
              
        //--- set market opening times to 0
        for(j=1;j<=n_obs;j++)
        {
         if(dates[j]!=null) 
         {
          lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
          if(n_rep > 1)
          {
           lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
           lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
          }
          price[j-1] = log(close_sp[j]);
                 
          if((dates[j]).indexOf(marketOpenTime) != -1)
          {
            for(m=0;m<n_rep;m++)
            {
             lastSeries[m*n_obs + j-1] = 0.0; 
            }
          }

          if(!isHoliday(dates[j])) 
          {          
           out.print(dates[j] + ", " + log(close_sp[j]) + ", ");
           for(m=0;m<n_rep-1;m++)
           {out.print(lastSeries[m*n_obs + j-1] + ", ");}
           out.println(lastSeries[(n_rep-1)*n_obs + j-1]);
           if(j==n_obs || j >= dates.length-1) {break;}
          }
         } 
        }
        out.close(); 
        System.out.println("Data successfully saved in " + file);
         //new_data_set = false;
        
         
     } catch (IOException e) {e.printStackTrace();} 
             
             
       
             
             
  }      
  
  
  
  
  public void downloadHistoricalIQDataFXSpreads(String name)
  {
      
      int i,j,m; 
      String dateStamp;
      String[] times;
      DateTime weekend; 
      String[] tokens; String ddelims = "[-]";   
      String[] tokens1; String delims = "[ ]";
      n_rep = 1;
      double spread;
      
      
      System.out.println(n_rep + " " + name);
      startIQConnect IQ = new startIQConnect(); IQ.run();     
      loadSpreads(name);
     
      String[] names = new String[1];
      names[0] = name;
      
      hilo_series = new ArrayList<String>();
      asset_series = new ArrayList<Double>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      
      System.out.println("Size of spreads hash = " + spreads.size());
      System.out.println(spreads.get("12:00:00"));
      
      n_rep = 3;
      try {getExplanatoryData(names,0);} catch (Exception e){};
      try {Thread.sleep (300);} catch (Exception e) {};
      
       
      lastSeries = new double[n_rep*n_obs];      
      String[] hilo = new String[n_obs + 1];
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      String[] dates = new String[n_obs+1];
      
      //System.out.println(asset_dates.size());
      int min_size = asset_dates.size() - 6;

      
      
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      
      while(i<=n_obs && min_size >= i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         
         if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
         //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
         
         else
         {
         
         if(i == dates_series.size()) {break;}         
         dateStamp = dates_series.get(dates_series.size() - 1 - i);
         //System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
         //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
                
         //--- dates_series1 -----
         if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
         {
           times = dateStamp.split("[ ]+"); 
           
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
         {
           i1++;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint  
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
         {
           i1=i1+2;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 4)))
         {
           i1=i1+3;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 5)))
         {
           i1=i1+4;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           i1++;
         }         
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 6)))
         {
           i1=i1+5;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
           i1++;
         }             
         
         else
         {
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
          //System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
         }
        
        
         dates[n_obs-i] = dateStamp; 
         i++;
        } 
      }
          
      //--- print to file
      file =  name +".dat";
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
              
        //--- set market opening times to 0
        for(j=1;j<=n_obs;j++)
        {
         if(dates[j]!=null) 
         {
          lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
          lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
          lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
         
         
                 
          if((dates[j]).indexOf(marketOpenTime) != -1 && zero_marketopen)
          {
           for(m=0;m<n_rep;m++)
           {
             lastSeries[m*n_obs + j-1] = 0.0; 
           }
          }

          if(!isHoliday(dates[j])) 
          {          
           out.print(dates[j] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
           for(m=0;m<n_rep;m++)
           {out.print(lastSeries[m*n_obs + j-1] + ", ");}
           out.println(hilo[j]);
           if(j==n_obs || j >= dates.length-1) {break;}
          }
         } 
       }
       out.close(); 
       System.out.println("Data successfully saved in " + file);
         //new_data_set = false;
            
      } catch (IOException e) {e.printStackTrace();} 
             
             
     
             
  }   
  
  
  
  public void downloadHistoricalIQDataFXSpreads30(String name) //download at 30min mark
  {
      
      int i,j,m; 
      String dateStamp;
      String[] times;
      DateTime weekend; 
      String[] tokens; String ddelims = "[-]";   
      String[] tokens1; String delims = "[ ]";
      n_rep = 1;
      double spread;
      
      sBeginFilterTime = "000000";
      sInterval = ""+(30*60);
      System.out.println(n_rep + " " + name);
      startIQConnect IQ = new startIQConnect(); IQ.run();     
      loadSpreads(name);
     
      String[] names = new String[1];
      names[0] = name;
      
      asset_series = new ArrayList<Double>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      
      System.out.println("Size of spreads hash = " + spreads.size());
      System.out.println(spreads.get("12:00:00"));
      
      n_rep = 3;
      try {getExplanatoryData(names,0);} catch (Exception e){};
      try {Thread.sleep (300);} catch (Exception e) {};
      
       
      lastSeries = new double[n_rep*n_obs];          
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      String[] dates = new String[n_obs+1];
      
      //System.out.println(asset_dates.size());
      int min_size = asset_dates.size() - 6;

      
      
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      
      while(i<=n_obs && min_size >= i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         
         if(tokens1[1].indexOf(":00:00") != -1) {i1++;}
         else if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
         //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
         
         else
         {
         
         if(i == dates_series.size()) {break;}         
         dateStamp = dates_series.get(dates_series.size() - 1 - i);
         //System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
         //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
                
         //--- dates_series1 -----
         if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
         {
           times = dateStamp.split("[ ]+"); 
           
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
         {
           i1++;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
         {
           i1=i1+2;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 4)))
         {
           i1=i1+3;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 5)))
         {
           i1=i1+4;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           i1++;
         }         
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 6)))
         {
           i1=i1+5;
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
           i1++;
         }             
         
         else
         {
           times = dateStamp.split("[ ]+");
           //System.out.println(times[1]);
           if(spreads.containsKey(times[1]))
           {spread = (spreads.get(times[1])).doubleValue();}
           else 
           {  spread = .0001;}
           
           close_nq[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1); //bid 
           close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
           close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint 
          //System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
         }
        
        
         dates[n_obs-i] = dateStamp; 
         i++;
        } 
      }
          
      //--- print to file
      file =  name +".dat";
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
              
        //--- set market opening times to 0
        for(j=1;j<=n_obs;j++)
        {
         if(dates[j]!=null) 
         {
          lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
          lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
          lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
         
         
                 
          if((dates[j]).indexOf(marketOpenTime) != -1 && zero_marketopen)
          {
           for(m=0;m<n_rep;m++)
           {
             lastSeries[m*n_obs + j-1] = 0.0; 
           }
          }

          if(!isHoliday(dates[j])) 
          {          
           out.print(dates[j] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
           for(m=0;m<n_rep-1;m++)
           {out.print(lastSeries[m*n_obs + j-1] + ", ");}
           out.println(lastSeries[(n_rep-1)*n_obs + j-1]);
           if(j==n_obs || j >= dates.length-1) {break;}
          }
         } 
       }
       out.close(); 
       System.out.println("Data successfully saved in " + file);
         //new_data_set = false;
            
      } catch (IOException e) {e.printStackTrace();} 
             
             
     
             
  }     
  
  
  
  
  public void downloadHistoricalIQDaily30(String name, int freq) //download at 30min mark
  {
      
      int i,j,m;  //int shift = 2;
      String dateStamp;
      String[] tokens; String ddelims = "[-]";   
      String[] tokens1; String delims = "[ ]";
      n_rep = 1;
      double spread;
      
      
      System.out.println("SIZE = " + n_obs);
      
      sBeginFilterTime = "000000";
      sInterval = ""+(freq*60);
      System.out.println(n_rep + " " + name);
      startIQConnect IQ = new startIQConnect(); IQ.run();     
      loadSpreads(name);
     
      String[] names = new String[1];
      names[0] = name;
      
      hilo_series = new ArrayList<String>();
      asset_series = new ArrayList<Double>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      
      System.out.println("Size of spreads hash = " + spreads.size());
      System.out.println(spreads.get("12:00:00"));
      
      n_rep = 3;
      try {getExplanatoryData(names,0);} catch (Exception e){};
      try {Thread.sleep (300);} catch (Exception e) {};
      
      n_obs = asset_dates.size();
      
      System.out.println("SIZE again = " + n_obs);
      
      String[] hilo = new String[n_obs + 1]; 
      lastSeries = new double[n_rep*n_obs];          
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      String[] dates = new String[n_obs+1];
      
      //System.out.println(asset_dates.size());
      int min_size = asset_dates.size() - 6;

      
      
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      
      while(i<=n_obs && min_size >= i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         

           
         dateStamp = asset_dates.get(asset_dates.size() -  i1 - 1);  
         dateStamp.split("[ ]+");   
          
/*         if(spreads.containsKey(times[1]))
         {spread = (spreads.get(times[1])).doubleValue();}
         else 
         {spread = .0001;}  */        
         
          spread = .00015;
          if(name.indexOf("JPY") != -1)
          {spread = .015;}
          
         close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
         close_tf[n_obs-i] = close_nq[n_obs-i] + spread; //ask 
         close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint   
         hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
         
         dates[n_obs-i] = dateStamp; 
         i++; i1++;
         
      }
          
      //--- print to file
      String[] bdate = sBeginDateTime.split("[ ]+");
      file =  name + "-" + bdate[0] +".dat";
      //file =  name + ".dat";
      
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
              
        //--- set market opening times to 0
        for(j=1;j<=n_obs;j++)
        {
         if(dates[j]!=null) 
         {
          lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
          lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
          lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
         
         
                 
          if((dates[j]).indexOf(marketOpenTime) != -1 && zero_marketopen)
          {
           for(m=0;m<n_rep;m++)
           {
             lastSeries[m*n_obs + j-1] = 0.0; 
           }
          }

          if(!isHoliday(dates[j])) 
          {          
           out.print(dates[j] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
           for(m=0;m<n_rep;m++)
           {out.print(lastSeries[m*n_obs + j-1] + ", ");}
           out.println(hilo[j]);
           if(j==n_obs || j >= dates.length-1) {break;}
          }
         } 
       }
       out.close(); 
       System.out.println("Data successfully saved in " + file);
         //new_data_set = false;
            
      } catch (IOException e) {e.printStackTrace();} 
             
             
     
             
  }       
  
  
  
  
  
  
  
  
  public void downloadHistoricalDailyData(String name)
  {
  
      
      int i,j,m;  //int shift = 2;
      String dateStamp;
      DateTime weekend; 
      String[] tokens; String ddelims = "[-]";   
      String[] tokens1; String delims = "[ ]";
      n_rep = 1;
      dailyData = true;
      sBeginFilterTime = "000000";
      sInterval = ""+(30*60);
      System.out.println(n_rep + " " + name);
      startIQConnect IQ = new startIQConnect(); IQ.run();     
      //loadSpreads(name);
     
      String[] names = new String[1];
      names[0] = name;
   
      asset_series = new ArrayList<Double>();
      hilo_series = new ArrayList<String>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      dates_series = new ArrayList<String>();
  
      
      n_rep = 3;
      try {getExplanatoryData(names,0);} catch (Exception e){};
      try {Thread.sleep (300);} catch (Exception e) {};
      
       
      n_obs = asset_dates.size();
  
  
  
      lastSeries = new double[n_rep*n_obs];      
      String[] hilo = new String[n_obs + 1];
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      String[] dates = new String[n_obs+1];
      
      asset_dates.size();

      
      
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      
      while(i<n_obs && asset_dates.size() > i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         
         if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
         //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
         
         else
         {
         
          if(i == asset_dates.size()) {break;}         
          dateStamp = asset_dates.get(asset_dates.size() - 1 - i);
         
          close_nq[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1); //bid 
          close_tf[n_obs-i] = close_nq[n_obs-i] + .00005; //ask 
          close_sp[n_obs-i] = (close_tf[n_obs-i] + close_nq[n_obs-i])/2.0; //midpoint
          hilo[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
          i1++;

          dates[n_obs-i] = dateStamp; 
          i++;
         } 
      }
          
      //--- print to file
      file =  name +".daily.dat";
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
              
        //--- set market opening times to 0
        for(j=1;j<=n_obs;j++)
        {
         if(dates[j]!=null) 
         {
          lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
          lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
          lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
         
         
                 
//           if((dates[j]).indexOf(marketOpenTime) != -1 && zero_marketopen)
//           {
//            for(m=0;m<n_rep;m++)
//            {
//              lastSeries[m*n_obs + j-1] = 0.0; 
//            }
//           }

          if(!isHoliday(dates[j])) 
          {          
           out.print(dates[j] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
           for(m=0;m<n_rep;m++)
           {out.print(lastSeries[m*n_obs + j-1] + ", ");}
           out.println(hilo[j]);
           if(j==n_obs || j >= dates.length-1) {break;}
          }
         } 
       }
       out.close(); 
       System.out.println("Data successfully saved in " + file);
         //new_data_set = false;
            
      } catch (IOException e) {e.printStackTrace();} 
             
  }
  
  
//   public void downloadHistoricalIQDataIBFX(String name)
//   {
//       
//       
//      String durationStr = "1 W";
//      String barSizeSetting = "30 mins";   
//      
// /*     String durationStr = "1 D";
//      String barSizeSetting = "5 mins";    */   
//      
// /*     String durationStr = "1 W";
//      String barSizeSetting = "30 mins";       */ 
//      
//      
//      if(client == null || !client.isConnected())
//      {onConnect();}
//      
//      setLogTrans(false);
//      
//      client.reqCurrentTime();  try {Thread.sleep (150);} catch (Exception e) {}; 
//       
//      HashMap<String,Double> local_spread = new HashMap<String,Double>();
//      HashMap<String,Integer> local_spread_count = new HashMap<String,Integer>();
//      
//       int countv; double val; double spread;
//       int i,j,m; double hilo; //int shift = 2;
//       double[] price;
//       String dateStamp;
//       int shift = 1;
//       DateTime weekend; 
//       String strline; String[] tokens; String ddelims = "[-]"; int n_toks;  
//       String[] tokens1; String delims = "[ ]";
//       //n_rep = 1;
//       String[] time;
//       String[] names = name.split(delims);
//       n_rep = names.length;
//       
//       System.out.println(n_rep + " " + name);
//       n_rep = 3;
// //       startIQConnect IQ = new startIQConnect();
// //       IQ.run();     
// //     
//       String[] forexPair = toIBForex(names[0]);
// 
//       Contract contract = new Contract();
//       contract.m_symbol = forexPair[0];
//       contract.m_exchange = "IDEALPRO"; //"CSFBALGO";
//       contract.m_secType = "CASH";
//       contract.m_currency = forexPair[1];   
//    
//       asset_series = new ArrayList<Double>();
//       hilo_series = new ArrayList<String>();
//       hi_series = new ArrayList<String>();
//       lo_series = new ArrayList<String>();      
//       
//       asset_series1 = new ArrayList<Double>();
//       asset_series2 = new ArrayList<Double>();      
//       asset_dates = new ArrayList<String>();   
//       asset_dates1 = new ArrayList<String>();
//       asset_dates2 = new ArrayList<String>();
//       dates_series = new ArrayList<String>();
//       
// //       for(m = 0; m < n_rep; m++)
// //       {
// //        try {getExplanatoryData(names,m);} catch (Exception e){};
// //       }
// //       try {Thread.sleep (300);} catch (Exception e) {};
//       
//      client.reqHistoricalData(10, contract, currentTimeString, durationStr, barSizeSetting, "MIDPOINT", 1, 1);       
//      try {Thread.sleep (1000);} catch (Exception e) {};
//      client.reqHistoricalData(11, contract, currentTimeString, durationStr, barSizeSetting, "BID", 1, 1);      
//      try {Thread.sleep (1000);} catch (Exception e) {};
//      client.reqHistoricalData(12, contract, currentTimeString, durationStr, barSizeSetting, "ASK", 1, 1); 
//      try {Thread.sleep (1000);} catch (Exception e) {}; 
//       
//        
//      n_obs = dates_series.size();
//      
//        
//       lastSeries = new double[n_rep*n_obs];       
//       price = new double[n_obs];     
//    
//       double[] close_nq = new double[n_obs+1];
//       double[] close_tf = new double[n_obs+1];
//       double[] close_sp = new double[n_obs+1];
//       String[] hilos = new String[n_obs + 1];
//       String[] dates = new String[n_obs+1];
//       
//       System.out.println(asset_dates.size() + " " + asset_dates1.size() + " " + asset_dates2.size());
//       int min_size = asset_dates.size() - 3;
//       if((asset_dates1.size() - 3) <  min_size && n_rep > 1) {min_size = asset_dates1.size() - 3;}
//       if((asset_dates2.size() - 3) <  min_size && n_rep > 1) {min_size = asset_dates2.size() - 3;}
//       
//       
//       i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
//       
//       System.out.println(n_obs + " " + min_size);
//       
//       while(i<=n_obs && min_size >= i1)
//       {
//          tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
//          tokens = tokens1[0].split(ddelims);
//          weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
//          
//          if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
//          //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
//          
//          else
//          {
//          System.out.println(weekend.dayOfWeek().getAsText());
//          dateStamp = dates_series.get(dates_series.size() - 1 - i1);
//          System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
//          //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
//                 
//          //--- dates_series1 -----
//          if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
//          {
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            //hilos[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
//            hilos[n_obs-i] = lo_series.get(lo_series.size() - i1 - 1) + ", " + hi_series.get(hi_series.size() - i1 - 1);
//            i1++;
//          }
//          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
//          {
//            i1++;
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            //hilos[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
//            hilos[n_obs-i] = lo_series.get(lo_series.size() - i1 - 1) + ", " + hi_series.get(hi_series.size() - i1 - 1);
//            i1++;
//          }
//          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
//          {
//            i1=i1+2;
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            //hilos[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
//            hilos[n_obs-i] = lo_series.get(lo_series.size() - i1 - 1) + ", " + hi_series.get(hi_series.size() - i1 - 1);
//            i1++;
//          }
//          else
//          {
//           close_sp[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1);
//           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - 1 - i1 - 1);
//           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - 1 - i1 - 1);
//           //hilos[n_obs-i] = hilo_series.get(hilo_series.size() - i1 - 1);
//           hilos[n_obs-i] = lo_series.get(lo_series.size() - i1 - 1) + ", " + hi_series.get(hi_series.size() - i1 - 1);
//           //System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
//          }
//  
//         
//         
//          dates[n_obs-i] = dateStamp; 
//          i++;
//         } 
//       }
//           
//       //--- print to file
//       file =  forexPair[0] + "" + forexPair[1] +".IB.dat";
//       try
//       {  
//         PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
//         PrintWriter spreadout = new PrintWriter(new FileWriter(new File("spreads5min_" + names[0] + ".dat")));     
//         //--- set market opening times to 0
//         for(j=1;j<n_obs;j++)
//         {
//          if(dates[j]!=null) 
//          {
//           lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
//           if(n_rep > 1)
//           {
//            lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
//            lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
//           }
//           price[j-1] = log(close_sp[j]);
//                  
//           if((dates[j]).indexOf(marketOpenTime) != -1)
//           {
//             for(m=0;m<n_rep;m++)
//             {
//              lastSeries[m*n_obs + j-1] = 0.0; 
//             }
//           }
// 
//           if(!isHoliday(dates[j])) 
//           {          
//            time = dates[j].split("[ ]+");
//            
//            if(!local_spread.containsKey(time[1]))
//            {
//             local_spread.put(time[1], (new Double(close_tf[j] - close_nq[j])));
//             
//             if(time[1].indexOf("16:00") != -1) {local_spread.put(time[1], .00015);}
//             
//             local_spread_count.put(time[1], 1);
//            }
// //            else
// //            {
// //              val = local_spread.get(time[1]);
// //              countv = local_spread_count.get(time[1]);
// //              
// //              val = val + (new Double(close_tf[j] - close_nq[j]));
// //              countv = countv + 1;
// //              
// //              local_spread.put(time[1], val);
// //              local_spread_count.put(time[1], countv);
// //            }
//            
//            out.print(dates[j+1] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
//            for(m=0;m<n_rep;m++)
//            {out.print(lastSeries[m*n_obs + j-1] + ", ");}
//            out.println(hilos[j]);
//            if(j==n_obs || j >= dates.length-1) {break;}
//           }
//          } 
//         }
//         out.close(); 
//         System.out.println("Data successfully saved in " + file);
//          //new_data_set = false;
//         
//         String[] keys = local_spread.keySet().toArray(new String[0]);
//         
//         for(i=0;i<keys.length;i++)
//         {
//           
//           spread = (local_spread.get(keys[i])).doubleValue()/local_spread_count.get(keys[i]);
//           spreadout.println(keys[i] + ", " + spread);
//         }
//         spreadout.close();
//         
//         
//          
//      } catch (IOException e) {e.printStackTrace();} 
//              
//              
//        
//              
//              
//   }      
//     
//   
//   
//   
//   
//   
//   public void downloadHistoricalIQDataIBFX30(String name)
//   {
//       
//       
//      String durationStr = "1 W";
//      String barSizeSetting = "30 mins";   
//      
//      if(client == null || !client.isConnected())
//      {onConnect();}
//      
//      setLogTrans(false);
//      
//      client.reqCurrentTime();  try {Thread.sleep (150);} catch (Exception e) {}; 
//       
//      HashMap<String,Double> local_spread = new HashMap<String,Double>();
//      
//       int i,j,m; double hilo; //int shift = 2;
//       double[] price;
//       String dateStamp;
//       int shift = 1;
//       DateTime weekend; 
//       String strline; String[] tokens; String ddelims = "[-]"; int n_toks;  
//       String[] tokens1; String delims = "[ ]";
//       //n_rep = 1;
//       String[] time;
//       String[] names = name.split(delims);
//       n_rep = names.length;
//       
//       System.out.println(n_rep + " " + name);
//       n_rep = 3;
// //       startIQConnect IQ = new startIQConnect();
// //       IQ.run();     
// //     
//       String[] forexPair = toIBForex(names[0]);
// 
//       Contract contract = new Contract();
//       contract.m_symbol = forexPair[0];
//       contract.m_exchange = "IDEALPRO"; //"CSFBALGO";
//       contract.m_secType = "CASH";
//       contract.m_currency = forexPair[1];   
//    
//       asset_series = new ArrayList<Double>();
//       asset_series1 = new ArrayList<Double>();
//       asset_series2 = new ArrayList<Double>();      
//       asset_dates = new ArrayList<String>();   
//       asset_dates1 = new ArrayList<String>();
//       asset_dates2 = new ArrayList<String>();
//       dates_series = new ArrayList<String>();
//       
// //       for(m = 0; m < n_rep; m++)
// //       {
// //        try {getExplanatoryData(names,m);} catch (Exception e){};
// //       }
// //       try {Thread.sleep (300);} catch (Exception e) {};
//       
//      client.reqHistoricalData(10, contract, currentTimeString, durationStr, barSizeSetting, "MIDPOINT", 1, 1);       
//      try {Thread.sleep (1000);} catch (Exception e) {};
//      client.reqHistoricalData(11, contract, currentTimeString, durationStr, barSizeSetting, "BID", 1, 1);      
//      try {Thread.sleep (1000);} catch (Exception e) {};
//      client.reqHistoricalData(12, contract, currentTimeString, durationStr, barSizeSetting, "ASK", 1, 1); 
//      try {Thread.sleep (1000);} catch (Exception e) {}; 
//       
//        
//      n_obs = dates_series.size();
//      
//        
//       lastSeries = new double[n_rep*n_obs];       
//       price = new double[n_obs];     
//    
//       double[] close_nq = new double[n_obs+1];
//       double[] close_tf = new double[n_obs+1];
//       double[] close_sp = new double[n_obs+1];
//       String[] dates = new String[n_obs+1];
//       
//       System.out.println(asset_dates.size() + " " + asset_dates1.size() + " " + asset_dates2.size());
//       int min_size = asset_dates.size() - 3;
//       if((asset_dates1.size() - 3) <  min_size && n_rep > 1) {min_size = asset_dates1.size() - 3;}
//       if((asset_dates2.size() - 3) <  min_size && n_rep > 1) {min_size = asset_dates2.size() - 3;}
//       
//       
//       i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
//       
//       System.out.println(n_obs + " " + min_size);
//       
//       while(i<=n_obs && min_size >= i1)
//       {
//          tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
//          tokens = tokens1[0].split(ddelims);
//          weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
//          
//          if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
//          //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
//          
//          else
//          {
//          System.out.println(weekend.dayOfWeek().getAsText());
//          dateStamp = dates_series.get(dates_series.size() - 1 - i1);
//          System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
//          //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
//                 
//          //--- dates_series1 -----
//          if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
//          {
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            i1++;
//          }
//          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
//          {
//            i1++;
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            i1++;
//          }
//          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
//          {
//            i1=i1+2;
//            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
//            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
//            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
//            i1++;
//          }
// //          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 4)))
// //          {
// //            i1=i1+3;
// //            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
// //            i1++;
// //          }
// //          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 5)))
// //          {
// //            i1=i1+4;
// //            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
// //            i1++;
// //          }         
// //          else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 6)))
// //          {
// //            i1=i1+5;
// //            close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i1 - 1);
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i1 - 1);
// //            i1++;
// //          }             
//          
//          else
//          {
//           close_sp[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1);
//           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - 1 - i1 - 1);
//           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - 1 - i1 - 1);
//           //System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
//          }
//         
// //          if(n_rep > 1)
// //          {
// // //           if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - shift - i2)))
// // //           {
// // //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2);
// // //            i2++;
// // //           }
// //           
// //          if(dateStamp.equals(asset_dates1.get(asset_dates1.size() -  i2 - 1)))
// //          {
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);
// //            i2++;
// //          }
// //          else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 2)))
// //          {
// //            i2++;
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
// //            i2++;
// //          }
// //          else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 3)))
// //          {
// //            i2=i2+2;
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
// //            i2++;
// //          }
// //          else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 4)))
// //          {
// //            i2=i2+3;
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
// //            i2++;
// //          }
// //          else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 5)))
// //          {
// //            i2=i2+4;
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
// //            i2++;
// //          }         
// //          else if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - i2 - 6)))
// //          {
// //            i2=i2+5;
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - i2 - 1);  
// //            i2++;
// //          }  
// //           
// // 
// //           else
// //           {
// //            close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2 - 1);
// //            //System.out.println("Missing data at " + dateStamp + " for 2nd explanatory series");
// //           }
// //          
// //          
// // //           if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - shift - i3)))
// // //           {
// // //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3);
// // //            i3++;
// // //           }
// // //           
// //           
// //           
// //          if(dateStamp.equals(asset_dates2.get(asset_dates2.size() -  i3 - 1)))
// //          {
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);
// //            i3++;
// //          }
// //          else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 2)))
// //          {
// //            i3++;
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
// //            i3++;
// //          }
// //          else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 3)))
// //          {
// //            i3=i3+2;
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
// //            i3++;
// //          }
// //          else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 4)))
// //          {
// //            i3=i3+3;
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
// //            i3++;
// //          }
// //          else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 5)))
// //          {
// //            i3=i3+4;
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
// //            i3++;
// //          }         
// //          else if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - i3 - 6)))
// //          {
// //            i3=i3+5;
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - i3 - 1);  
// //            i3++;
// //          }            
// //           
// //   
// //           
// //           else
// //           {
// //            close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3 - 1); 
// //            //System.out.println("Missing data at " + dateStamp + " for 3rd explanatory series");    
// //           }        
// //          }
//         
//         
//          dates[n_obs-i] = dateStamp; 
//          i++;
//         } 
//       }
//           
//       //--- print to file
//       file =  forexPair[0] + "" + forexPair[1] +".IB.dat";
//       try
//       {  
//         PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
//         PrintWriter spreadout = new PrintWriter(new FileWriter(new File("spreads5min_" + names[0] + ".dat")));     
//         //--- set market opening times to 0
//         for(j=1;j<=n_obs;j++)
//         {
//          if(dates[j]!=null) 
//          {
//           lastSeries[j-1] = log(close_sp[j]) - log(close_sp[j-1]);
//           if(n_rep > 1)
//           {
//            lastSeries[n_obs + j-1] = log(close_nq[j]) - log(close_nq[j-1]);
//            lastSeries[2*n_obs + j-1] = log(close_tf[j]) - log(close_tf[j-1]);
//           }
//           price[j-1] = log(close_sp[j]);
//                  
//           if((dates[j]).indexOf(marketOpenTime) != -1)
//           {
//             for(m=0;m<n_rep;m++)
//             {
//              lastSeries[m*n_obs + j-1] = 0.0; 
//             }
//           }
// 
//           if(!isHoliday(dates[j])) 
//           {          
//            time = dates[j].split("[ ]+");
//            
//            if(!local_spread.containsKey(time[1]))
//            {
//             local_spread.put(time[1], (new Double(close_tf[j] - close_nq[j])));
//            }
//            
//            out.print(dates[j] + ", " + log(close_sp[j]) + ", " + log(close_nq[j]) + ", " + log(close_tf[j]) + ", ");
//            for(m=0;m<n_rep-1;m++)
//            {out.print(lastSeries[m*n_obs + j-1] + ", ");}
//            out.println(lastSeries[(n_rep-1)*n_obs + j-1]);
//            if(j==n_obs || j >= dates.length-1) {break;}
//           }
//          } 
//         }
//         out.close(); 
//         System.out.println("Data successfully saved in " + file);
//          //new_data_set = false;
//         
//         String[] keys = local_spread.keySet().toArray(new String[0]);
//         
//         for(i=0;i<keys.length;i++)
//         {
//           spreadout.println(keys[i] + ", " + local_spread.get(keys[i]));
//         }
//         spreadout.close();
//         
//         
//          
//      } catch (IOException e) {e.printStackTrace();} 
//              
//              
//        
//              
//              
//   }        
//   
//   
//   
//   
  
  
  
  
 /* 
  
  public void downloadHistoricalIQData24(String name)
  {
      
      int i,j,m; double hilo; //int shift = 2;
      double[] price;
      String dateStamp;
      int shift = 1;
      DateTime weekend; 
      String strline; String[] tokens; String ddelims = "[-]"; int n_toks;  
      String[] tokens1; String delims = "[ ]";
      n_rep = 1;
      
      String[] names = name.split(delims);
      n_rep = names.length;
      
      
      startIQConnect IQ = new startIQConnect();
      IQ.run();     
   

   
   
      asset_series = new ArrayList<Double>();
      asset_series1 = new ArrayList<Double>();
      asset_series2 = new ArrayList<Double>();      
      asset_dates = new ArrayList<String>();   
      asset_dates1 = new ArrayList<String>();
      asset_dates2 = new ArrayList<String>();
      
      for(m = 0; m < n_rep; m++)
      {
       try {getExplanatoryData(names,m);} catch (Exception e){};
      }
      try {Thread.sleep (300);} catch (Exception e) {};
      
       
      lastSeries = new double[n_rep*n_obs];       
      price = new double[n_obs];     
   
      double[] close_nq = new double[n_obs+1];
      double[] close_tf = new double[n_obs+1];
      double[] close_sp = new double[n_obs+1];
      
      ArrayList<Double> close_24 = new ArrayList<Double>();
      ArrayList<String> dates_24 = new ArrayList<String>();
      
      String[] dates = new String[n_obs+1];
      
      System.out.println(asset_dates.size());
      i=0; i1 = 0; m = 0; i2 = 0; i3 = 0;
      while(i<=n_obs && (asset_dates.size() - 6) >= i1)
      {
         tokens1 = (asset_dates.get(asset_dates.size() - 1 - i1)).split(delims);
         tokens = tokens1[0].split(ddelims);
         weekend = new DateTime((new Integer(tokens[0])).intValue(), (new Integer(tokens[1])).intValue(), (new Integer(tokens[2])).intValue(), 14, 0); 
         
         if(weekend.dayOfWeek().getAsText().equals("Saturday") || weekend.dayOfWeek().getAsText().equals("Sunday")) {i1++;}
         //else if(weekend.dayOfWeek().getAsText().equals("Friday")) {i1++;}
         
         else
         {
         dateStamp = dates_series.get(dates_series.size() - 1 - i);
         //System.out.println(dateStamp + " " + asset_dates.get(asset_dates.size() - 1 - i1) + " " + asset_series.get(asset_series.size() - 1- i1));
         //if(i1>1) System.out.println(asset_dates.get(asset_dates.size() - i1) + " " + asset_dates.get(asset_dates.size() - i1 - 2));
                
         //--- dates_series1 -----
         if(dateStamp.equals(asset_dates.get(asset_dates.size() -  i1 - 1)))
         {
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 2)))
         {
           i1++;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 3)))
         {
           i1=i1+2;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 4)))
         {
           i1=i1+3;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 5)))
         {
           i1=i1+4;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }         
         else if(dateStamp.equals(asset_dates.get(asset_dates.size() - i1 - 6)))
         {
           i1=i1+5;
           close_sp[n_obs-i] = asset_series.get(asset_series.size() - i1 - 1);  
           i1++;
         }             
         
         else
         {
          close_sp[n_obs-i] = asset_series.get(asset_series.size() - 1 - i1 - 1);
          System.out.println("Missing data at " + dateStamp + " for 1st explanatory series");
         }
        
         if(n_rep > 1)
         {
          if(dateStamp.equals(asset_dates1.get(asset_dates1.size() - shift - i2)))
          {
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2);
           i2++;
          }
          else
          {
           close_nq[n_obs-i] = asset_series1.get(asset_series1.size() - shift - i2 - 1);
           System.out.println("Missing data at " + dateStamp + " for 2nd explanatory series");
          }
         
         
          if(dateStamp.equals(asset_dates2.get(asset_dates2.size() - shift - i3)))
          {
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3);
           i3++;
          }
          else
          {
           close_tf[n_obs-i] = asset_series2.get(asset_series2.size() - shift - i3 - 1); 
           System.out.println("Missing data at " + dateStamp + " for 3rd explanatory series");    
          }        
         }
        
        
         dates[n_obs-i] = dateStamp; 
         i++;
        } 
      }
      
      
      file =  names[0] +".dat";
      try
      {  
        PrintWriter out = new PrintWriter(new FileWriter(new File(file)));
      
        for(j=1;j<n_obs;j++)
        {
      
          if(dates[j]!=null)
          {
        
           if(dates[j].indexOf(dailyTime) != -1 && !isHoliday(dates[j]))
           {
           
            close_24.add(close_sp[j]);
            
            if(close_24.size() > 1)
            {out.println(dates[j] + ", " + log(close_24.get(close_24.size()-1).doubleValue()) + ", " + (log(close_24.get(close_24.size()-1).doubleValue()) - log(close_24.get(close_24.size()-2).doubleValue())));}
           }
          }
          if(j==n_obs || j >= dates.length-1) {break;}
         } 
         out.close(); System.out.println("Data successfully saved in " + file);
       }
       catch (IOException e) {e.printStackTrace();} 
      
      
  }     
      */
  
 
    
  public void getExplanatoryData(String[] names, int k) throws Exception
  {
                     
                String[] tokens;
		//sSymbol = expSymbols[c];//"@ESU13";
		//int k = 0;
		 sMaxDatapoints = "";
                 sSymbol = names[k];
                 System.out.println("Downloading data for " + sSymbol);
		 sDataDirection = "'1'";
		 sRequestID = "";
		 sDatapointsPerSend = "";
		 
		 String sBeginDate = "20070101";
		 String sEndDate = "20141013";
		 
		 System.out.println("HIT," + sSymbol + "," + sInterval + "," + sBeginDateTime + "," + sEndDateTime + ",," + sBeginFilterTime + "," + sEndFilterTime + "," + sDataDirection + "," + sRequestID + "," + sDatapointsPerSend);
		 
		 if(dailyData)
		 {s_outs[k].write("HDT," + sSymbol + "," + sBeginDate + "," + sEndDate + "," + sMaxDatapoints + "," + sDataDirection + "," + sRequestID + "," + sDatapointsPerSend + "\r\n");}
		 else
		 {s_outs[k].write("HIT," + sSymbol + "," + sInterval + "," + sBeginDateTime + "," + sEndDateTime + ",," + sBeginFilterTime + "," + sEndFilterTime + "," + sDataDirection + "," + sRequestID + "," + sDatapointsPerSend + "\r\n");}
		 
		 
		 s_outs[k].flush();
		 String sLine = s_ins[k].readLine();
		
		 while (sLine != null && sLine.indexOf("!ENDMSG!", 0) == -1)  
		 {
			//System.out.println(sLine+","+ncount);
			if (sLine.equals("E,!SYNTAX_ERROR!,") || sLine.equals("!ENDMSG!"))
			{
				sLine = null;
			}
			else
			{
				sLine = s_ins[k].readLine();
				
				tokens = sLine.split(delims);
				if(tokens.length > 1)
				{
				  if(k==0)
				  {
				   asset_series.add(new Double(tokens[4]));
				   hilo_series.add(tokens[2] + ", " + tokens[1]);
				   asset_dates.add(tokens[0]);
				  }
				  else if(k==1)
				  {
				   asset_series1.add(new Double(tokens[4]));
				   asset_dates1.add(tokens[0]);				  
				  }
				  else if(k==2)
				  {
				   asset_series2.add(new Double(tokens[4]));
				   asset_dates2.add(tokens[0]);				  
				  }
				}
			}
		}
		finished_data = true;
	     	
    }
    
    
/*
   public void historicalData(int reqId, String date, double open, double high, double low, double close, 
                               int volume, int count, double WAP, boolean hasGaps)
   {   
        
     System.out.println(EWrapperMsgGenerator.historicalData(reqId, date, open, high, low, close, volume, count, WAP, hasGaps));  
   
     if(open < 0)
     {
       System.out.println("Historical data complete");
       finished_data = true; 
       return;
     }
     else //if(reqId == es_requestId)
     {
     
       String[] dates = date.split("[ ]+");
       String year = dates[0].substring(0, 4);
       String month = dates[0].substring(4, 6);
       String day = dates[0].substring(6, 8);
   
       if(dates[1].indexOf("17:15:00") != -1) //change to 17:00:00
       {dates[1] = "17:00:00";} 
   
       date = new String(year + "-" + month + "-" + day + " " + dates[1]);
       //System.out.println(date);
       
       if(reqId == 10) //midprice
       {
       
        asset_series.add(new Double(close));  
        asset_dates.add(date);
        dates_series.add(date);
       //String msg = EWrapperMsgGenerator.historicalData(reqId, date, open, high, low, close, volume, count, WAP, hasGaps);  
          //System.out.println(msg);           
       }
       else if(reqId == 11) 
       {
        asset_series1.add(new Double(close));  
        asset_dates1.add(date);
        hi_series.add(""+high);
       
       }
       else if(reqId == 12)
       { 
        asset_series2.add(new Double(close));  
        asset_dates2.add(date);
        lo_series.add(""+low);
       }
       
     }  
   }            
    
    
    
    */
    
    

//   public void getExplanatoryData(String[] names, int k) throws Exception
//   {
//                      
//                 int ncount = 0; 
//                 String[] tokens;
// 		double himlo;
// 		//sSymbol = expSymbols[c];//"@ESU13";
// 		//int k = 0;
// 		 sMaxDatapoints = "";
//                  sSymbol = names[k];
//                  System.out.println("Downloading data for " + sSymbol);
// 		 sDataDirection = "'1'";
// 		 sRequestID = "";
// 		 sDatapointsPerSend = "";
// 		 s_outs[k].write("HIT," + sSymbol + "," + sInterval + "," + sBeginDateTime + "," + sEndDateTime + ",," + sBeginFilterTime + "," + sEndFilterTime + "," + sDataDirection + "," + sRequestID + "," + sDatapointsPerSend + "\r\n");
// 		 s_outs[k].flush();
// 		 String sLine = s_ins[k].readLine();
// 		
// 		
// 		
// 		
// 		
// 		//client.reqHistoricalData(requestIds[m], contract[m], currentTimeString, durationStr, barSizeSetting, "MIDPOINT", 1, 1);
// 		
// 		
// 		
// 		
// 		
// 		 while (sLine != null && sLine.indexOf("!ENDMSG!", 0) == -1)  
// 		 {
// 			System.out.println(sLine+","+ncount);
// 			if (sLine.equals("E,!SYNTAX_ERROR!,") || sLine.equals("!ENDMSG!"))
// 			{
// 				sLine = null;
// 			}
// 			else
// 			{
// 				sLine = s_ins[k].readLine();
// 				
// 				tokens = sLine.split(delims);
// 				if(tokens.length > 1)
// 				{
// 				  if(k==0)
// 				  {
// 				   asset_series.add(new Double(tokens[4]));
// 				   asset_dates.add(tokens[0]);
// 				  }
// 				  else if(k==1)
// 				  {
// 				   asset_series1.add(new Double(tokens[4]));
// 				   asset_dates1.add(tokens[0]);				  
// 				  }
// 				  else if(k==2)
// 				  {
// 				   asset_series2.add(new Double(tokens[4]));
// 				   asset_dates2.add(tokens[0]);				  
// 				  }
// 				  ncount++;
// 				}
// 			}
// 		}
// 		finished_data = true;
// 	     	
//     } 
//  
 
 
 
 
 
 

    

    
    
  public static void main(String args[])
  {  
     if(args.length > 0)
     {
       new String(args[0]);
     }
 
     HistoricalData data = new HistoricalData();
//      data.min_30 = true;
     data.setupHistoricalDownload();
//      
     //data.onConnect();
     data.setLogTrans(false);
     
     try {Thread.sleep(9000);}
      catch (Exception e)
		     {e.printStackTrace();}   
     
     data.setupHistoricalDownload("SGDJPY.FXCM");
     data.downloadHistoricalIQDaily30("SGDJPY.FXCM",30);         
     
     
     
//      data.setLogTrans(true);
//      data.setupHistoricalDownload("VXX");
//      data.downloadHistoricalIQDaily30("VXX",5);       
// 
//      data.setupHistoricalDownload("EEM");
//      data.downloadHistoricalIQDaily30("EEM",5);       
//      
//      data.setupHistoricalDownload("SPY");
//      data.downloadHistoricalIQDaily30("SPY",5);   
//      
//      data.setupHistoricalDownload("QQQ");
//      data.downloadHistoricalIQDaily30("QQQ",5);   
//      
//      data.setupHistoricalDownload("AAPL");
//      data.downloadHistoricalIQDaily30("AAPL",5);   
//      
//      data.setupHistoricalDownload("XLE");
//      data.downloadHistoricalIQDaily30("XLE",5);    
//      
//      data.setupHistoricalDownload("XLV");
//      data.downloadHistoricalIQDaily30("XLV",5);    
//      
//      data.setupHistoricalDownload("GDX");
//      data.downloadHistoricalIQDaily30("GDX",5);     
     
     
/*     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",30);        
     
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",30);       
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",30);        
     
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",30);       */ 
  
/*     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURUAD.FXCM",30);     
  
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",30); */     
  
     
     //data.downloadHistoricalIQDataIBFX("NOKJPY.FXCM");
     //data.downloadHistoricalIQDataIBFX("GBPNOK.FXCM");
     
     //data.setupHistoricalDownload("AUDUSD.FXCM");
     //data.downloadHistoricalIQDaily30("AUDUSD.FXCM");
/*     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",30);    downloadHistoricalIQDaily30
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",30);        
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",30);    
     
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",30);    
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",30);    
     
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",30);  */       

/*     data.setupHistoricalDownload("EURSGD.ABBA");
     data.downloadHistoricalIQDaily30("EURSGD.ABBA",30);     */   

/*          data.setupHistoricalDownload("EURSGD.BNFF");
     data.downloadHistoricalIQDaily30("EURSGD.BNFF",30); 
     
          data.setupHistoricalDownload("EURSGD.MORN");
     data.downloadHistoricalIQDaily30("EURSGD.MORN",30);  
 
          data.setupHistoricalDownload("EURSGD.MIGF");
     data.downloadHistoricalIQDaily30("EURSGD.MIGF",30);   
 
          data.setupHistoricalDownload("EURSGD.SAXO");
     data.downloadHistoricalIQDaily30("EURSGD.SAXO",30);
     
               data.setupHistoricalDownload("EURSGD.FCFX");
     data.downloadHistoricalIQDaily30("EURSGD.FCFX",30);  */ 
 
//                data.setupHistoricalDownload("EURGBP.FCFX");
//      data.downloadHistoricalIQDaily30("EURGBP.FCFX",30); 
//      
//                     data.setupHistoricalDownload("EURUSD.FCFX");
//      data.downloadHistoricalIQDaily30("EURUSD.FCFX",30); 
 
 /*     




     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",30);    
     */
/*     data.setupHistoricalDownload("CADCHF.FXCM");
     data.downloadHistoricalIQDaily30("CADCHF.FXCM",30);        
  
     data.setupHistoricalDownload("NZDCAD.FXCM");
     data.downloadHistoricalIQDaily30("NZDCAD.FXCM",30);      
    
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",30);        
    
     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",30);        
     
     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",30);     */   
   
/*     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",30);               
   
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",30);     
   
        data.setupHistoricalDownload("EURJPY.FXCM");
        data.downloadHistoricalIQDaily30("EURJPY.FXCM",30);       
        
        data.setupHistoricalDownload("EURNZD.FXCM");
        data.downloadHistoricalIQDaily30("EURNZD.FXCM",30);  */        
   
//      data.setupHistoricalDownload("AUDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("AUDCHF.FXCM",30);               
//    
//         data.setupHistoricalDownload("EURNZD.FXCM");
//         data.downloadHistoricalIQDaily30("EURNZD.FXCM",30);    
//      
//         data.setupHistoricalDownload("GBPJPY.FXCM");
//         data.downloadHistoricalIQDaily30("GBPJPY.FXCM",30);        
//   
//         data.setupHistoricalDownload("CHFJPY.FXCM");
//         data.downloadHistoricalIQDaily30("CHFJPY.FXCM",30);          
// 
// //         data.setupHistoricalDownload("USDJPY.FXCM");
// //         data.downloadHistoricalIQDaily30("USDJPY.FXCM",30); 
//         
        
//         data.setupHistoricalDownload("USDJPY.FXCM");
//         data.downloadHistoricalIQDaily30("USDJPY.FXCM",30);   downloadHistoricalIQDaily30
//    
    

/*     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",1);     
    
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",1); 
    
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",1);       
     
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",1);  
     
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",1);     
     
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",1);       
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",1);       

     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",1);
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",1);       
     
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",1);   */ 
     
     
     
     
/*     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",5);       
      
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",5);  
      
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",5);        
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",5);       
     
     data.setupHistoricalDownload("AUDJPY.FXCM");
     data.downloadHistoricalIQDaily30("AUDJPY.FXCM",5);      
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",5);      
     
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",5);       
     
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",5);       
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",5);        
     
     data.setupHistoricalDownload("AUDUSD.FXCM");
     data.downloadHistoricalIQDaily30("AUDUSD.FXCM",5);          
     
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",5);   
     
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",5);   
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",5);      
     
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",5);          
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",5);             
     
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",5);           
     
     data.setupHistoricalDownload("NZDCAD.FXCM");
     data.downloadHistoricalIQDaily30("NZDCAD.FXCM",5);       
     
     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",5);            
     
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",5);  
     
     data.setupHistoricalDownload("CADCHF.FXCM");
     data.downloadHistoricalIQDaily30("CADCHF.FXCM",5);  
     
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",5);       
     
     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",5);    
     
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",5);           
     
     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",5);  
     
     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",5);       
     
     data.setupHistoricalDownload("AUDNZD.FXCM");
     data.downloadHistoricalIQDaily30("AUDNZD.FXCM",5);   */    
     
/*     data.setLogTrans(true);
     
     data.setupHistoricalDownload("VGK");
     data.downloadHistoricalIQDaily30("VGK",1);
      
     data.setupHistoricalDownload("HEDJ");
     data.downloadHistoricalIQDaily30("HEDJ",1);     
      
     data.setupHistoricalDownload("EZU");
     data.downloadHistoricalIQDaily30("EZU",1);        
      
     data.setupHistoricalDownload("EWG");
     data.downloadHistoricalIQDaily30("EWG",1);      
      
     data.setupHistoricalDownload("FEZ");
     data.downloadHistoricalIQDaily30("FEZ",1);      
     
     data.setupHistoricalDownload("EWU");
     data.downloadHistoricalIQDaily30("EWU",1);      
      
     data.setupHistoricalDownload("FAZ");
     data.downloadHistoricalIQDaily30("FAZ",1);      
          
     data.setupHistoricalDownload("DIG");
     data.downloadHistoricalIQDaily30("DIG",1); */       
     
    
     
     
/*     data.setupHistoricalDownload("IWD");
     data.downloadHistoricalIQDaily30("IWD",1);      
     
     data.setupHistoricalDownload("IWF");
     data.downloadHistoricalIQDaily30("IWF",1);       
     
     data.setupHistoricalDownload("JNK");
     data.downloadHistoricalIQDaily30("JNK",1);    */       
     
//      data.setupHistoricalDownload("EURJPY.FXCM");
//      data.downloadHistoricalIQDaily30("EURJPY.FXCM",1);    

//      data.setupHistoricalDownload("USDNOK.FXCM");
//      data.downloadHistoricalIQDaily30("USDNOK.FXCM",1);         
//      
//      data.setupHistoricalDownload("EURSEK.FXCM");
//      data.downloadHistoricalIQDaily30("EURSEK.FXCM",1);        
// 
//      data.setupHistoricalDownload("USDSEK.FXCM");
//      data.downloadHistoricalIQDaily30("USDSEK.FXCM",1);      
//      
//      data.setupHistoricalDownload("NOKJPY.FXCM");
//      data.downloadHistoricalIQDaily30("NOKJPY.FXCM",1);   
    
   
     
     
//      
//      data.setupHistoricalDownload("USDCAD.FXCM");
//      data.downloadHistoricalIQDaily30("USDCAD.FXCM",1);    
//      
//      data.setupHistoricalDownload("GBPCAD.FXCM");
//      data.downloadHistoricalIQDaily30("GBPCAD.FXCM",1);    
// 
//      data.setupHistoricalDownload("EURCAD.FXCM");
//      data.downloadHistoricalIQDaily30("EURCAD.FXCM",1);    
//      
    
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",1);  
    
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",1);  
    
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",1);    
    
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",1);       
      
     data.setupHistoricalDownload("CADCHF.FXCM");
     data.downloadHistoricalIQDaily30("CADCHF.FXCM",1);  
      
     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",1);  
      
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",1);        
     
     data.setupHistoricalDownload("AUDNZD.FXCM");
     data.downloadHistoricalIQDaily30("AUDNZD.FXCM",1);  
  
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",1); 
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",1);  
    
     data.setupHistoricalDownload("EURSGD.COMP");
     data.downloadHistoricalIQDaily30("EURSGD.COMP",1);    
     
    /* 
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",1);       
      
     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",1);  
      
     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",1);        
     
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",1);  
  
     data.setupHistoricalDownload("EURGBP.FXCM");
     data.downloadHistoricalIQDaily30("EURGBP.FXCM",1); 
     
     data.setupHistoricalDownload("NZDCAD.FXCM");
     data.downloadHistoricalIQDaily30("NZDCAD.FXCM",1);      
      
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",1);           
        
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",1);    

     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",1);      
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",1);    
    
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",1);      
    
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",1);    
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",1);    

     data.setupHistoricalDownload("AUDNZD.FXCM");
     data.downloadHistoricalIQDaily30("AUDNZD.FXCM",1);         
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",1);            
     
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",1);      
     
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",1);         
     
     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",1);      
     
     data.setupHistoricalDownload("AUDUSD.FXCM");
     data.downloadHistoricalIQDaily30("AUDUSD.FXCM",1);        
     
     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",1);     
  
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",1);      
  
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",1);    
          
     data.setupHistoricalDownload("AUDJPY.FXCM");
     data.downloadHistoricalIQDaily30("AUDJPY.FXCM",1);  
  
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",1);  
  
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",1); */
  
//      data.setupHistoricalDownload("NZDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("NZDCHF.FXCM",1);   
//      
//      data.setupHistoricalDownload("NZDCAD.FXCM");
//      data.downloadHistoricalIQDaily30("NZDCAD.FXCM",1); 
     
 
/*     data.setupHistoricalDownload("@ME#");
     data.downloadHistoricalIQDaily30("@ME#",1);    

     data.setupHistoricalDownload("@EM#");
     data.downloadHistoricalIQDaily30("@EM#",1);         
  
     data.setupHistoricalDownload("@ME#");
     data.downloadHistoricalIQDaily30("@ME#",1);    

     data.setupHistoricalDownload("LS#");
     data.downloadHistoricalIQDaily30("LS#",1);    

     data.setupHistoricalDownload("@SF#");
     data.downloadHistoricalIQDaily30("@SF#",1);       
     
     data.setupHistoricalDownload("GBE#");
     data.downloadHistoricalIQDaily30("GBE#",1);  */   
     
/*     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",1);    
  
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",1);     
  
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",1);         
     
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",1);         
     
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",1);      
     
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",1);          
     
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",1);      */    
     
     
  
  
/*     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",5);    
 
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",5);    
 
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",5);    

     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",5); 
     
     data.setupHistoricalDownload("EURGBP.FXCM");
     data.downloadHistoricalIQDaily30("EURGBP.FXCM",5);         
 
     data.setupHistoricalDownload("NZDCAD.FXCM");
     data.downloadHistoricalIQDaily30("NZDCAD.FXCM",5);   
  
     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",5);   
  
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",5);        
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",5);    
     
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",5);        
     
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",5);      

     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",5);       
     
     data.setupHistoricalDownload("GBPNZD.FXCM");
     data.downloadHistoricalIQDaily30("GBPNZD.FXCM",5);       
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",5);        
     
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",5); 
         
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",5);     
                
     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",5);     
     
     data.setupHistoricalDownload("AUDNZD.FXCM");
     data.downloadHistoricalIQDaily30("AUDNZD.FXCM",5);        
     
     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",5);      
     
     data.setupHistoricalDownload("AUDUSD.FXCM");
     data.downloadHistoricalIQDaily30("AUDUSD.FXCM",5);               
     
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",5);      
     
     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",5);   
    
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",5);      
     
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",5);       
    
     data.setupHistoricalDownload("AUDJPY.FXCM");
     data.downloadHistoricalIQDaily30("AUDJPY.FXCM",5);   
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",5);   
     
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",5);     
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",5);   */       
     

/*     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",5);  
     
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",5);       
     
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",5);       
     
     data.setupHistoricalDownload("EURCHF.FXCM");
     data.downloadHistoricalIQDaily30("EURCHF.FXCM",5);       
   
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",5);  
     
     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",5);            
   
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",5);  */      
        
     
/*     data.setupHistoricalDownload("AUDCHF.FXCM");
     data.downloadHistoricalIQDaily30("AUDCHF.FXCM",15);      
  
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",15);      
    
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",15);        
        
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",15);                
          
     data.setupHistoricalDownload("CADCHF.FXCM");
     data.downloadHistoricalIQDaily30("CADCHF.FXCM",15);              
           
     data.setupHistoricalDownload("EURSGD.COMP");
     data.downloadHistoricalIQDaily30("EURSGD.COMP",15); 
     
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",15);     
  
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",15);   
  
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",15);     
    
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",15);     
        
     data.setupHistoricalDownload("GBPNZD.FXCM");
     data.downloadHistoricalIQDaily30("GBPNZD.FXCM",15);       
    
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",15);
     
     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",15);       
    
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",15);         
    
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",15);      
      
     data.setupHistoricalDownload("AUDUSD.FXCM");
     data.downloadHistoricalIQDaily30("AUDUSD.FXCM",15);         
     
     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",15);        
     
     data.setupHistoricalDownload("AUDJPY.FXCM");
     data.downloadHistoricalIQDaily30("AUDJPY.FXCM",15);     
     
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",15);         
     
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",15);           
     
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",15);                
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",15);    */         
     
       
//     
    
/*     data.setupHistoricalDownload("@NE#");
     data.downloadHistoricalIQDaily30("@NE#",30);    
 
     data.setupHistoricalDownload("@SF#");
     data.downloadHistoricalIQDaily30("@SF#",30);  
 
     data.setupHistoricalDownload("@RY#");
     data.downloadHistoricalIQDaily30("@RY#",30);     
 
     data.setupHistoricalDownload("@SJY#");
     data.downloadHistoricalIQDaily30("@SJY#",30);   
 
     data.setupHistoricalDownload("@BR#");
     data.downloadHistoricalIQDaily30("@BR#",30);   
     
     data.setupHistoricalDownload("@PSF#");
     data.downloadHistoricalIQDaily30("@PSF#",30);        

     data.setupHistoricalDownload("@PX#");
     data.downloadHistoricalIQDaily30("@PX#",30);  */      
      
//      data.setupHistoricalDownload("USDSGD.FXCM");
//      data.downloadHistoricalIQDaily30("USDSGD.FXCM",30);  
 
//      data.setupHistoricalDownload("USDMXN.FXCM");
//      data.downloadHistoricalIQDaily30("USDMXN.FXCM",30);  
     
//      data.setupHistoricalDownload("XAUUSD.FXCM");
//      data.downloadHistoricalIQDaily30("XAUUSD.FXCM",30); 
//      
//      data.setupHistoricalDownload("XAUEUR.FXCM");
//      data.downloadHistoricalIQDaily30("XAUEUR.FXCM",30);      
//      
//      data.setupHistoricalDownload("XAUAUD.FXCM");
//      data.downloadHistoricalIQDaily30("XAUAUD.FXCM",30);      
     
     
//      data.setupHistoricalDownload("+GC#");
//      data.downloadHistoricalIQDaily30("+GC#",30); 
     
//      data.setupHistoricalDownload("AUDCAD.COMP");
//      data.downloadHistoricalIQDaily30("AUDCAD.COMP",30);      
//      
//      data.setupHistoricalDownload("NZDJPY.FXCM");
//      data.downloadHistoricalIQDaily30("NZDJPY.FXCM",30);    
//      
//      data.setupHistoricalDownload("AUDJPY.FXCM");
//      data.downloadHistoricalIQDaily30("AUDJPY.FXCM",30);        
     
/*     data.setupHistoricalDownload("GBPHKD.COMP");
     data.downloadHistoricalIQDaily30("GBPHKD.COMP",30);    
     
     data.setupHistoricalDownload("USDHKD.COMP");
     data.downloadHistoricalIQDaily30("USDHKD.COMP",30);  */  
     
/*     data.setupHistoricalDownload("GBPNZD.FXCM");
     data.downloadHistoricalIQDaily30("GBPNZD.FXCM",30);    */  
     
/*     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",30);               
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",30);        
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",30);    
     
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",30);    
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",30);    
     
     data.setupHistoricalDownload("USDCAD.FXCM");
     data.downloadHistoricalIQDaily30("USDCAD.FXCM",30);    */ 
     
     
     
/*     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",30);       
 
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",30);    */   
 
/*     data.setupHistoricalDownload("@VX#");
     data.downloadHistoricalIQDaily30("@VX#",5);  

     data.setupHistoricalDownload("@ES#");
     data.downloadHistoricalIQDaily30("@ES#",5);   
 
     data.setupHistoricalDownload("@NQ#");
     data.downloadHistoricalIQDaily30("@NQ#",5);   
 
     data.setupHistoricalDownload("@D1#");
     data.downloadHistoricalIQDaily30("@D1#",5);   
 
     data.setupHistoricalDownload("+NG#");
     data.downloadHistoricalIQDaily30("+NG#",5);    
     
     data.setupHistoricalDownload("+SI#");
     data.downloadHistoricalIQDaily30("+SI#",5);    
     
     data.setupHistoricalDownload("+HG#");
     data.downloadHistoricalIQDaily30("+HG#",5);   */      
//     @VX#
//     

//      data.setupHistoricalDownload("EURJPY.FXCM");
//      data.downloadHistoricalIQDaily30("EURJPY.FXCM",30);        
// 
//      data.setupHistoricalDownload("EURUSD.FXCM");
//      data.downloadHistoricalIQDaily30("EURUSD.FXCM",30);         
//      
//      data.setupHistoricalDownload("GBPUSD.FXCM");
//      data.downloadHistoricalIQDaily30("GBPUSD.FXCM",30);     
//      
//      data.setupHistoricalDownload("USDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("USDCHF.FXCM",30);          
//  
//      data.setupHistoricalDownload("XAUUSDO.ABBA");
//      data.downloadHistoricalIQDaily30("XAUUSDO.ABBA",30);   
//      
//      data.setupHistoricalDownload("XAUCHFO.COMP");
//      data.downloadHistoricalIQDaily30("XAUCHFO.COMP",30);     
//      
//      data.setupHistoricalDownload("XAUEURO.COMP");
//      data.downloadHistoricalIQDaily30("XAUEURO.COMP",30);      
     
//      data.setupHistoricalDownload("EEM");
//      data.downloadHistoricalIQDaily30("EEM",30);    
//      XAUCHFO.COMP    XAUEURO.COMP
//      data.setupHistoricalDownload("QQQ");
//      data.downloadHistoricalIQDaily30("QQQ",30);    
//      

//      data.setupHistoricalDownload("VXX");
//      data.downloadHistoricalIQDaily30("VXX",30);        
//      
//      data.setupHistoricalDownload("XLE");
//      data.downloadHistoricalIQDaily30("XLE",30);    
//      
//      data.setupHistoricalDownload("XLV");
//      data.downloadHistoricalIQDaily30("XLV",30);    
//      
//      data.setupHistoricalDownload("GDX");
//      data.downloadHistoricalIQDaily30("GDX",30);   
/*     
     data.setupHistoricalDownload("BP");
     data.downloadHistoricalIQDaily30("BP",30);    
     
     data.setupHistoricalDownload("LUV");
     data.downloadHistoricalIQDaily30("LUV",30);    
     
     data.setupHistoricalDownload("MET");
     data.downloadHistoricalIQDaily30("MET",30);         
     
     data.setupHistoricalDownload("XOM");
     data.downloadHistoricalIQDaily30("XOM",30);   */        
     
     
 /*    data.setupHistoricalDownload("TXN");
     data.downloadHistoricalIQDaily30("TXN",30); */       
     
//      data.setupHistoricalDownload("MCD");
//      data.downloadHistoricalIQDaily30("MCD",30);    
//      
//      data.setupHistoricalDownload("JNJ");
//      data.downloadHistoricalIQDaily30("JNJ",30);    
//      
//      data.setupHistoricalDownload("IBM");
//      data.downloadHistoricalIQDaily30("IBM",30);           
// //      
//      
     
//      data.setupHistoricalDownload("USDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("USDCHF.FXCM",15);  
         

//      data.setupHistoricalDownload("GBPAUD.FXCM");
//      data.downloadHistoricalIQDaily30("GBPAUD.FXCM",15); 
//      
//      data.setupHistoricalDownload("EURCAD.FXCM");
//      data.downloadHistoricalIQDaily30("EURCAD.FXCM",15);      
//   
//      data.setupHistoricalDownload("EURAUD.FXCM");
//      data.downloadHistoricalIQDaily30("EURAUD.FXCM",15); 

/*     data.setupHistoricalDownload("GBPUSD.COMP");
     data.downloadHistoricalIQDaily30("GBPUSD.COMP",15);      
     
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",15);  
     
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",15);    */ 
     
/*     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",1);  */ 
    
/*     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",1);        

     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",1);        
         
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",1);        
     
     data.setupHistoricalDownload("EURSGD.COMP");
     data.downloadHistoricalIQDaily30("EURSGD.COMP",1);        
     
     data.setupHistoricalDownload("GBPJPY.FXCM");
     data.downloadHistoricalIQDaily30("GBPJPY.FXCM",1);       
    
     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",1);       
    
     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",1);           
    
     data.setupHistoricalDownload("CHFJPY.FXCM");
     data.downloadHistoricalIQDaily30("CHFJPY.FXCM",1);       
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",1);           
     
     data.setupHistoricalDownload("AUDNZD.FXCM");
     data.downloadHistoricalIQDaily30("AUDNZD.FXCM",1);   
         
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",1); */     
    
          
/*   EURNZD.FXCM", "GBPCHF.FXCM", "EURSGD.COMP", "AUDSGD.COMP", "AUDNZD.FXCM", "GBPJPY.FXCM", "CHFJPY.FXCM", 
		             "USDJPY.FXCM", "GBPUSD.FXCM", "AUDCAD.FXCM", "USDCHF.FXCM */  
     
     
//      EURSGD.FXCM", "AUDSGD.FXCM", "CHFJPY.FXCM", "NZDCHF.FXCM", "GBPAUD.FXCM", "EURCAD.FXCM
//      data.setupHistoricalDownload("EURUSD.FXCM");
//      data.downloadHistoricalIQDaily30("EURUSD.FXCM",15);       
     
     
/*     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",15);  */         
     
/*     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",15);        
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",15);    
     
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",15);    
     
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",15);         
           
     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",15);    
     
     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",15); */       
     
     
     
/*     data.setupHistoricalDownload("EURCAD.FXCM");
     data.downloadHistoricalIQDaily30("EURCAD.FXCM",15);  */      
         
         
/*     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",15);         
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",15);      
     
     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",15);     
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",15);        
     
     data.setupHistoricalDownload("EURSGD.COMP");
     data.downloadHistoricalIQDaily30("EURSGD.COMP",15);        
      
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",15);      
     
     data.setupHistoricalDownload("NZDUSD.FXCM");
     data.downloadHistoricalIQDaily30("NZDUSD.FXCM",15);      
     
     data.setupHistoricalDownload("USDJPY.FXCM");
     data.downloadHistoricalIQDaily30("USDJPY.FXCM",15);      */ 
     
     
/*     data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",15); 
     
     data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",15); 
     
     data.setupHistoricalDownload("EURJPY.FXCM");
     data.downloadHistoricalIQDaily30("EURJPY.FXCM",15); 
     
     data.setupHistoricalDownload("EURUSD.FXCM");
     data.downloadHistoricalIQDaily30("EURUSD.FXCM",15);      
 
     data.setupHistoricalDownload("USDCHF.FXCM");
     data.downloadHistoricalIQDaily30("USDCHF.FXCM",15);  
     
     data.setupHistoricalDownload("EURGBP.FXCM");
     data.downloadHistoricalIQDaily30("EURGBP.FXCM",15);   */  
   
/*      data.setupHistoricalDownload("GBPAUD.FXCM");
     data.downloadHistoricalIQDaily30("GBPAUD.FXCM",15);     
     
      data.setupHistoricalDownload("GBPUSD.FXCM");
     data.downloadHistoricalIQDaily30("GBPUSD.FXCM",15);   */     
     
/*     data.setupHistoricalDownload("EURNOK.FXCM");
     data.downloadHistoricalIQDaily30("EURNOK.FXCM",30);    */     
          
/*     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",30);    
     
     data.setupHistoricalDownload("GBPCAD.FXCM");
     data.downloadHistoricalIQDaily30("GBPCAD.FXCM",15);         
     
     data.setupHistoricalDownload("GBPCHF.FXCM");
     data.downloadHistoricalIQDaily30("GBPCHF.FXCM",15); 

     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",15);           
     
     data.setupHistoricalDownload("CADJPY.FXCM");
     data.downloadHistoricalIQDaily30("CADJPY.FXCM",15);     
     
     data.setupHistoricalDownload("EURNZD.FXCM");
     data.downloadHistoricalIQDaily30("EURNZD.FXCM",15);        
     
     data.setupHistoricalDownload("AUDSGD.COMP");
     data.downloadHistoricalIQDaily30("AUDSGD.COMP",15);    */   
     
/*     data.setupHistoricalDownload("EURSGD.COMP");
     data.downloadHistoricalIQDaily30("EURSGD.COMP",15);   */ 

//      data.setupHistoricalDownload("AUDCAD.FXCM");
//      data.downloadHistoricalIQDaily30("AUDCAD.FXCM",15);        
//  
//      data.setupHistoricalDownload("AUDUSD.FXCM");
//      data.downloadHistoricalIQDaily30("AUDUSD.FXCM",15);   
 
     
//      
//      data.setupHistoricalDownload("USDSGD.FXCM");
//      data.downloadHistoricalIQDaily30("USDSGD.FXCM",30);
//      
//      data.setupHistoricalDownload("AUDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("AUDCHF.FXCM",30);
//      
//      data.setupHistoricalDownload("GBPCHF.FXCM");
//      data.downloadHistoricalIQDaily30("GBPCHF.FXCM",15);     
//      
//      data.setupHistoricalDownload("NZDCHF.FXCM");
//      data.downloadHistoricalIQDaily30("NZDCHF.FXCM",30);        
//      
//      data.setupHistoricalDownload("EURSGD.COMP");
//      data.downloadHistoricalIQDaily30("EURSGD.COMP",30);             
//      
//      data.setupHistoricalDownload("AUDSGD.COMP");
//      data.downloadHistoricalIQDaily30("AUDSGD.COMP",15);              
//      
//      data.setupHistoricalDownload("USDSGD.FXCM");
//      data.downloadHistoricalIQDaily30("USDSGD.FXCM",30);  
// 
 
// 
//      data.setupHistoricalDownload("AUDNZD.FXCM");
//      data.downloadHistoricalIQDaily30("AUDNZD.FXCM",15);    
// 
//      data.setupHistoricalDownload("AUDCAD.FXCM");
//      data.downloadHistoricalIQDaily30("AUDCAD.FXCM",15);    
// 
//      data.setupHistoricalDownload("EURUSD.FXCM");
//      data.downloadHistoricalIQDaily30("EURUSD.FXCM",30);    






     
     
     
     
     
/*     data.setupHistoricalDownload("NZDCHF.FXCM");
     data.downloadHistoricalIQDaily30("NZDCHF.FXCM",30);   
     
     data.setupHistoricalDownload("EURGBP.FXCM");
     data.downloadHistoricalIQDaily30("EURGBP.FXCM",30);   
     
     data.setupHistoricalDownload("CADCHF.FXCM");
     data.downloadHistoricalIQDaily30("CADCHF.FXCM",30);   
     
     
     
     data.setupHistoricalDownload("EURAUD.FXCM");
     data.downloadHistoricalIQDaily30("EURAUD.FXCM",30);        
     
     data.setupHistoricalDownload("GBPNZD.FXCM");
     data.downloadHistoricalIQDaily30("GBPNZD.FXCM",30);    
     
     data.setupHistoricalDownload("NZDJPY.FXCM");
     data.downloadHistoricalIQDaily30("NZDJPY.FXCM",30);    
     
     data.setupHistoricalDownload("AUDJPY.FXCM");
     data.downloadHistoricalIQDaily30("AUDJPY.FXCM",30);      
     
     data.setupHistoricalDownload("AUDCAD.FXCM");
     data.downloadHistoricalIQDaily30("AUDCAD.FXCM",30);          
     
     data.setupHistoricalDownload("AUDUSD.FXCM");
     data.downloadHistoricalIQDaily30("AUDUSD.FXCM",30);       */    
     
        
     
     
     
//      data.downloadHistoricalDailyData("CADJPY.FXCM");
//      data.downloadHistoricalDailyData("AUDJPY.FXCM");
//      data.downloadHistoricalDailyData("EURUSD.FXCM");
//      data.downloadHistoricalDailyData("GBPUSD.FXCM");
//      data.downloadHistoricalDailyData("EURGBP.FXCM");    
     
//      data.downloadHistoricalIQDataIBFX("AUDJPY.FXCM");
//      data.downloadHistoricalIQDataIBFX("GBPCHF.FXCM");
//      data.downloadHistoricalIQDataIBFX("AUDUSD.FXCM");
//      data.downloadHistoricalIQDataIBFX("GBPJPY.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURAUD.FXCM");
//      data.downloadHistoricalIQDataIBFX("GBPUSD.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURGBP.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURCHF.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURJPY.FXCM");
//      data.downloadHistoricalIQDataIBFX("USDJPY.FXCM");
//      data.downloadHistoricalIQDataIBFX("USDCHF.FXCM");
//      data.downloadHistoricalIQDataIBFX("AUDCAD.FXCM");
//      data.downloadHistoricalIQDataIBFX("NZDUSD.FXCM");
//      data.downloadHistoricalIQDataIBFX("AUDCHF.FXCM");     
//      data.downloadHistoricalIQDataIBFX("EURCAD.FXCM"); 
//      data.downloadHistoricalIQDataIBFX("GBPCAD.FXCM");
//      data.downloadHistoricalIQDataIBFX("GBPAUD.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURNZD.FXCM");
//      data.downloadHistoricalIQDataIBFX("EURNOK.FXCM");
//      data.downloadHistoricalIQDataIBFX("NZDCAD.FXCM");
//      data.downloadHistoricalIQDataIBFX("GBPHKD.FXCM");
     //data.downloadHistoricalIQDataIBFX("EURUSD.FXCM");
     //data.downloadHistoricalIQDataIBFX("CHFJPY.FXCM");
     //data.downloadHistoricalIQDataIBFX("GBPNZD.FXCM");
//      data.downloadHistoricalIQDataIBFX("USDCAD.FXCM");
//      data.downloadHistoricalIQDataIBFX("CADCHF.FXCM");
//      data.onDisconnect();
//      for(m = 0; m < 6;m++)
//      {data.downloadHistoricalIQData(m);}
   
   
   
  }
  

  
  
  
  /*
    
    
    public void currentTime(long time) 
    {
      try {
        currentTimeString = (dateF.format(new Date(1000*time))).toString() + " EST";
        //currentTimeString = currentTimeString 
        //System.out.println(currentTimeString);
        //tws_panel.dialogPrintln(currentTimeString);
        //tws_panel.uploadHistoricalText.setText(currentTimeString);
        timeID = (int)time;
        currentTimeID = 1000*time;
      }            
      catch (Exception e) {e.printStackTrace ();}
    }
 
 
    public void tickPrice(int tickerId, int field, double price,int canAutoExecute ) {
//         switch (field){
//             case 1:  //bid
//                 System.out.println("tickerID = " + tickerId); tws_panel.dialogPrintln("tickerID = " + tickerId); 
//                 System.out.println("Bid Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("Bid Price = ",String.valueOf(price));
//                 lastBidTick = price; 
//                 break;
//             case 2:  //ask
//                 System.out.println("Ask Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("Ask Price = ",String.valueOf(price));
//                 lastAskTick = price; 
//                 break;
//             case 4:  //last
//                 System.out.println("Last Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("Last Price = ",String.valueOf(price));
//                 lastLastTick = price; 
//                 break;
//             case 6:  //high
//                 System.out.println("High Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("High Price = ",String.valueOf(price));
//                 lastHighTick = price; 
//                 break;
//             case 7:  //low
//                 lastLowTick = price; 
//                 System.out.println("Low Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("Low Price = ",String.valueOf(price));
//                 break;
//             case 9:  //close
//                 lastCloseTick[current_m] = price;
//                 System.out.println("Close Price = "+String.valueOf(price)); tws_panel.dialogPricePrintln("Close Price = ",String.valueOf(price));
//                 break;
//         }
    }

    public void tickSize(int tickerId, int field, int size) {
        switch (field){
            case 0:   //bid
                //System.out.println("Bid Size = "+String.valueOf(size));
                break;
            case 3:   //ask
                //System.out.println("Ask Size = "+String.valueOf(size));
                break;
            case 5:   //last
                //System.out.println("Last Size = "+String.valueOf(size));
                break;
            case 8:   //volume
                //System.out.println("Volume = "+String.valueOf(size));
                break;
        }
    }
 

    public void realtimeBar (int reqId, long time, double open, double high,
            double low, double close, long volume, double wap, int count)
    {


            
    } 
 
 
 
    
    public void bondContractDetails (int reqId, ContractDetails contractDetails)
    {
    }

    public void contractDetails (int reqId, ContractDetails contractDetails)
    {
    }

    public void contractDetailsEnd (int reqId)
    {
    }
    
    public void commissionReport(CommissionReport commissionReport)
    {
    }

    public void tickSnapshotEnd(int reqId)
    {
    }
    
    public void deltaNeutralValidation(int reqId, UnderComp underComp)
    {
    }
    
    public void execDetailsEnd( int reqId)
    {
    }
    
    public void accountDownloadEnd(String accountName)
    {
    }
    
    public void openOrderEnd()
    {
    }
    
    public void tickOptionComputation( int tickerId, int field, double impliedVol,
    		double delta, double optPrice, double pvDividend,
    		double gamma, double vega, double theta, double undPrice)
    {
    }
    
    public void marketDataType(int reqId, int marketDataType)
    {
    }
    
    public void fundamentalData (int reqId, String data)
    {
    }

    public void bondContractDetails (ContractDetails contractDetails)
    {
    }

    public void contractDetails (ContractDetails contractDetails)
    {
    } 
 
    //------ Unimplemented wrapper functions -------------
 
    public void execDetails (int orderId, Contract contract, Execution execution)
    {
    }

    public void managedAccounts (String accountsList)
    {

      
      //for(int i =0; i<account_numbers.length; i++) {System.out.println(account_numbers[i]);}   
    }

    public void openOrder (int orderId, Contract contract, Order order,
            OrderState orderState)
    {

    
    }

    public void orderStatus (int orderId, String status, int filled,
            int remaining, double avgFillPrice, int permId, int parentId,
            double lastFillPrice, int clientId, String whyHeld)
    {
       

       
    
    }

    public void receiveFA (int faDataType, String xml)
    {
    }

    public void scannerData (int reqId, int rank,
            ContractDetails contractDetails, String distance, String benchmark,
            String projection, String legsStr)
    {
    }

    public void scannerDataEnd (int reqId)
    {
    }

    public void scannerParameters (String xml)
    {
    }

    public void tickEFP (int symbolId, int tickType, double basisPoints,
            String formattedBasisPoints, double impliedFuture, int holdDays,
            String futureExpiry, double dividendImpact, double dividendsToExpiry)
    {
    }

    public void tickGeneric (int symbolId, int tickType, double value)
    {
    }

    public void tickOptionComputation (int symbolId, int field,
            double impliedVol, double delta, double modelPrice,
            double pvDividend)
    {
    }

    public void tickString (int symbolId, int tickType, String value)
    {
    }

    public void updateAccountTime (String timeStamp)
    {
    }

    public void updateAccountValue (String key, String value, String currency,
            String accountName)
    {
    }

    public void updateMktDepth (int symbolId, int position, int operation,
            int side, double price, int size)
    {
    }

    public void updateMktDepthL2 (int symbolId, int position,
            String marketMaker, int operation, int side, double price, int size)
    {
    }

    public void updateNewsBulletin (int msgId, int msgType, String message,
            String origExchange)
    {
    }

    public void updatePortfolio (Contract contract, int position,
            double marketPrice, double marketValue, double averageCost,
            double unrealizedPNL, double realizedPNL, String accountName)
    {
    
    }

    public void connectionClosed ()
    {
    }

    public void error (Exception e)
    {
        e.printStackTrace ();
    }

    public void error (String str)
    {
        System.err.println (str); //tws_panel.dialogPrintln(str);
    }

    public void error (int id, int errorCode, String errorMsg)
    {
        System.err.println ("error (id, errorCode, errorMsg): id=" + id + ".  errorCode=" + errorCode + ".  errorMsg=" + errorMsg);
        System.out.println("TWS Message: " + "error (id, errorCode, errorMsg): id=" + id + ".  errorCode=" + errorCode + ".  errorMsg=" + errorMsg);
    }
    
    public void nextValidId (int orderId)
    {
        
    } 
 
   
  
  
  */
  
  
  
  
  
  
  
  
  
  
  class startIQConnect extends Thread
  {
    
  
    
     public void run()
     {
        System.out.println("start up IQ");
		try {

			// Launch IQFeed and Register the app with IQFeed.
			System.out.println("Launching IQConnect.");
			//Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe ‑product SIMON_OTZIGER_6345 -version 2.9.0.13 ‑login 422113 ‑password 76059316");
			System.out.println("Verifying if IQConnect is connected to the server");
			Thread.sleep(9000);
			// verify everything is ready to send commands.
			boolean bConnected = false;
			// connect to the admin port.
			//Socket sockAdmin = new Socket(InetAddress.getByName("192.168.56.10"), 8300);

			System.out.println("Launching IQConnect.");
                        //Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe -product PRABIN_SETH_11606 -version 2.9.0.13 ‑login 438106 ‑password 64568341");
                        //Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe ‑product VIKRAM_KURIYAN_11481 -version 2.9.0.13 ‑login 447691 ‑password 17876547");
                       
                        Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe ‑product INCUBE_GROUP_11864 -version 2.9.0.13 ‑login 449312 ‑password 92480265");
                       
			//Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe -product NAVID_GHATRI_11666 -version 2.9.0.13 ‑login 438955 ‑password 66923996");
			System.out.println("Verifying if IQConnect is connected to the server");
// 			Thread.sleep(12000);
// 			// verify everything is ready to send commands.
			//boolean bConnected = false;
// 			// connect to the admin port.
			Socket sockAdmin = new Socket(InetAddress.getByName("localhost"), 9300);			
			
			BufferedReader bufreadAdmin = new BufferedReader(new InputStreamReader(sockAdmin.getInputStream()));
			BufferedWriter bufwriteAdmin = new BufferedWriter(new OutputStreamWriter(sockAdmin.getOutputStream()));
			String sAdminLine = "";
			// loop while we are still connected to the admin port or until we are connected
			while (((sAdminLine = bufreadAdmin.readLine()) != null) && !bConnected)
			{
				System.out.println(sAdminLine);
				if (sAdminLine.indexOf(",Connected,") > -1)
				{
					System.out.println("IQConnect is connected to the server.");
					bConnected = true;
				}
				else if (sAdminLine.indexOf(",Not Connected,") > -1)
				{
					System.out.println("IQConnect is Not Connected.\r\nSending connect command.");
					bufwriteAdmin.write("S,CONNECT\r\n");
					bufwriteAdmin.flush();
				}
			}
			// cleanup admin port connection
			sockAdmin.shutdownOutput();
			sockAdmin.shutdownInput();
			sockAdmin.close();
			bufreadAdmin.close();
			bufwriteAdmin.close();

			// at this point, we are connected and the feed is ready.
			
			
				
			// creates a socket connection to localhost (IP address 127.0.0.1) on port 9100.
			// This is that port that IQFeed listens on for lookup requests. 
// 			Socket s = new Socket(InetAddress.getByName("192.168.56.10"), 8100);
// 	                //Socket s = new Socket(InetAddress.getByName("localhost"), 9100);
// 			// buffer to incomming data.
// 			sin = new BufferedReader(new InputStreamReader(s.getInputStream()));
// 			// buffer for out going commands.
// 			sout = new BufferedWriter(new OutputStreamWriter(s.getOutputStream()));
// 			// buffer for incomming menu commands, that the user of the application enters.
// 			in = new BufferedReader(new InputStreamReader(System.in));
// 
//                         // Set the lookup port to protocol 5.0 to allow for millisecond times, 
//                         // market center, trade conditions, etc
//                         sout.write("S,SET PROTOCOL,5.0\r\n");
//                         sout.flush();
                        
                        //getFiveMinuteHistoricalData();
                       //n_rep = 1;
                       s_ins = new BufferedReader[n_rep];
                       s_outs = new BufferedWriter[n_rep];

                       Socket[] socks = new Socket[n_rep];
                       
	               System.out.println("Launching IQConnect.");
	               for(int i=0;i<n_rep;i++)
	               {
	                //socks[i] = new Socket(InetAddress.getByName("192.168.56.10"), 8100);
	                socks[i] = new Socket(InetAddress.getByName("localhost"), 9100);
	                s_ins[i] = new BufferedReader(new InputStreamReader(socks[i].getInputStream()));
	                s_outs[i] = new BufferedWriter(new OutputStreamWriter(socks[i].getOutputStream()));

                        s_outs[i].write("S,SET PROTOCOL,5.0\r\n");
                        s_outs[i].flush();
                       }                         
                        
                        
                        
                        
                        
                     }
                     catch (Exception e)
		     {e.printStackTrace();}    
    
    
      } 
 }   
  
  
  //---- Generate 30min times for Forex

  public void generateTimesForex_30()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      //holidays[1] = "07-03";
      //holidays[2] = "07-04";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      

      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":" + emin + ":00");
    
      int endingInt = (new Integer(endTime).intValue());
    

    
 
 
 
 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      
      System.out.println("Ending time! == " + endingInt);
      
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {
       
       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday")))// && !isHoliday(date)) //skip weekends  
       {
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) <= endingInt && HistoricalData.getTime(timeinc) != 0) //get the time
        { 
         //System.out.println(historicalData.getTime(timeinc) + " " + endingInt);
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
        
         
         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         }         
         timeinc = timeinc.plusMinutes(obs_freq);                
       
        
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }    
  
  
  
  public void generateTimesForex_30_EST()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      //holidays[1] = "07-03";
      //holidays[2] = "07-04";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      
      //Moving holidays
//       holidays[5] = "2013-09-02";
//       holidays[6] = "2012-09-03";
//       //holidays[6] = "2011-09-05";
//       
//       holidays[7] = "2013-05-27";
//       holidays[8] = "2012-05-28";
//       holidays[9] = "2011-05-30";
//     
//       holidays[10] = "2013-11-28";
//       holidays[11] = "2012-11-22";
//       holidays[12] = "2011-11-24";
//     
//       holidays[13] = "2013-01-21";
//       holidays[14] = "2012-01-23";
//       holidays[15] = "2011-01-24";    
//     
//       holidays[16] = "2013-02-18";
//       holidays[17] = "2012-02-20";
//       holidays[18] = "2011-02-21";    
//     
//       holidays[19] = "2012-10-29";
//       holidays[20] = "2012-10-30";
//       holidays[21] = "2013-03-29";
//       holidays[22] = "2013-08-22";     
//       holidays[23] = "2013-11-29";
//     
//       holidays[24] = "2014-01-20";
//       holidays[25] = "2014-02-17";
//       holidays[26] = "2014-04-18";
     // holidays[27] = "2014-05-26";    
    
    
      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue()+30; 
      int ehour = (new Integer(sehour)).intValue()+1; int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":" + emin + ":00");
    
      //int endingInt = (new Integer(endTime).intValue()) + 3000;
      int endingInt = 160000;

    
 
 
 
 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      
      System.out.println("Ending time! == " + endingInt);
      
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {
       
       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday")))// && !isHoliday(date)) //skip weekends  
       {
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) <= endingInt && HistoricalData.getTime(timeinc) != 0) //get the time
        { 
         //System.out.println(historicalData.getTime(timeinc) + " " + endingInt);
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
        
         
         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         }         
         timeinc = timeinc.plusMinutes(obs_freq);                
       
        
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }    
  
  
  
  
  
  
  
  
  
  
  
  
  
  public void generateTimesForex_60()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      

      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
//       System.out.println(bhour + ":" + bmin + ":00");
//       System.out.println(ehour + ":" + emin + ":00");
    
      int endingInt = (new Integer(endTime).intValue());
    

    
 
 
 
 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      
      System.out.println("Ending time! == " + endingInt);
      
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {
//        System.out.println(timeinc.toString() + " " + end.toString());
//        System.out.println(timeinc);
       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday")))// && !isHoliday(date)) //skip weekends  
       {
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) <= endingInt)// && historicalData.getTime(timeinc) != 0) //get the time
        { 
         //System.out.println(historicalData.getTime(timeinc) + " " + endingInt);
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
         date = dt_tokens[0];   

         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {
          dates_series.add(date + " " + time); //System.out.println(date + " " + time); 
         
          if(HistoricalData.getTime(timeinc) == endingInt) 
          {break;}      
         }
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time); /*System.out.println(date + " " + time);*/ }
           else {break;}
         }         
         timeinc = timeinc.plusMinutes(obs_freq);                
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }      
  
  
  
  
  public void generateTimesForex_All()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      

      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":" + emin + ":00");
    
      int endingInt = (new Integer(endTime).intValue());
    

    
 
 
 
 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      
      System.out.println("Ending time! == " + endingInt);
      
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {
//        System.out.println(timeinc.toString() + " " + end.toString());
//        System.out.println(timeinc);
       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday")))// && !isHoliday(date)) //skip weekends  
       {
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) <= endingInt)// && historicalData.getTime(timeinc) != 0) //get the time
        { 
         //System.out.println(historicalData.getTime(timeinc) + " " + endingInt);
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
         date = dt_tokens[0];   

         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {
          dates_series.add(date + " " + time); System.out.println(date + " " + time); 
         
          if(HistoricalData.getTime(timeinc) == endingInt) 
          {break;}      
         }
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
           else {break;}
         }         
         timeinc = timeinc.plusMinutes(obs_freq);                
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }        
  
  
  
  
  
  
  
  //------ generates times for 24-FOREX -------------------
  
  public void generateTimes24()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      System.out.println("24 hour FOREX");
    
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      

    
      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      (new Integer(seyear)).intValue(); (new Integer(semonth)).intValue(); (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":0" + emin + ":00");
    
      //int endingInt = (new Integer(endTime).intValue()+1500);
      int endingInt = 233000;

 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin+15, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin+5, 0, 0);     
      DateTime start1 = start.plusDays(1);
      //DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      DateTime end = new DateTime();
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      System.out.println("Ending time! == " + endingInt);
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {

       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday"))) // && !timeinc.dayOfWeek().getAsText().equals("Friday"))// && !isHoliday(date)) //skip weekends  
       {
       
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) != endingInt) //get the time
        { 
        
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
        
         
         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time); System.out.println(date + " " + time); }
         }
         
         timeinc = timeinc.plusMinutes(obs_freq);                
                
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        
        
        if(timeinc.dayOfWeek().getAsText().equals("Saturday"))
        {day = day.plus(nextDay); day = day.plus(nextDay); timeinc = new DateTime(day);}
        
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];            
        timetokens = dt_tokens[1].split(dotdelim); 
        time = timetokens[0];         
         
        dates_series.add(date + " " + time); System.out.println(date + " " + time); 
        
        
        
        timeinc = timeinc.plusMinutes(obs_freq);   

       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }  
   
   
   
  public void generateTimes5()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      System.out.println("5 min FOREX");
    
      n_holidays = 3;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      holidays[1] = "12-24";
      holidays[2] = "12-25";
      

    
      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":0" + emin + ":00");
    
      int endingInt = (new Integer(endTime).intValue()+500);
    

 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin+5, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin+5, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin+5, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin+5, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      System.out.println("Ending time! == " + endingInt);
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {

       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday"))) // && !timeinc.dayOfWeek().getAsText().equals("Friday"))// && !isHoliday(date)) //skip weekends  
       {
       
        //System.out.println("Before the Change! : " + date);
        //dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
        //System.out.println("It's the Big day! : " + date);
 
 
        while(HistoricalData.getTime(timeinc) <= endingInt) //get the time
        { 
        
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
        
         System.out.println(date + " " + time); 
         
         if(!timeinc.dayOfWeek().getAsText().equals("Friday")) //if not friday, all is good
         {dates_series.add(date + " " + time);}
         else                                                  //if Friday, make sure before 5pm
         {
           if(HistoricalData.getTime(timeinc) <= 17*10000)
           {dates_series.add(date + " " + time);}
         }
         
         timeinc = timeinc.plusMinutes(obs_freq);                
       
        
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size();
     
   }     
   
   
  
  
  //--- generates times for 12 hour trading 
   
  public void generateTimes()
  {
  
      String date, time;
      String[] dt_tokens;
      String[] timetokens;
      String Tdelim = "[T]+";
      String dotdelim = "[.]+";
      String delim = "[ ]+";
      n_holidays = 28;
      holidays = new String[n_holidays];
      
      holidays[0] = "01-01";
      holidays[1] = "07-03";
      holidays[2] = "07-04";
      holidays[3] = "12-24";
      holidays[4] = "12-25";
      
      //Moving holidays
      holidays[5] = "2013-09-02";
      holidays[6] = "2012-09-03";
      holidays[6] = "2011-09-05";
      
      holidays[7] = "2013-05-27";
      holidays[8] = "2012-05-28";
      holidays[9] = "2011-05-30";
    
      holidays[10] = "2013-11-28";
      holidays[11] = "2012-11-22";
      holidays[12] = "2011-11-24";
    
      holidays[13] = "2013-01-21";
      holidays[14] = "2012-01-23";
      holidays[15] = "2011-01-24";    
    
      holidays[16] = "2013-02-18";
      holidays[17] = "2012-02-20";
      holidays[18] = "2011-02-21";    
    
      holidays[19] = "2012-10-29";
      holidays[20] = "2012-10-30";
      holidays[21] = "2013-03-29";
      holidays[22] = "2013-08-22";     
      holidays[23] = "2013-11-29"; 
    
      holidays[24] = "2014-01-20";
      holidays[25] = "2014-02-17";
      holidays[26] = "2014-04-18";
      holidays[27] = "2014-05-26";    
    
    
      dates_series = new ArrayList<String>();
    
      String endingDate;
      
      
      String[] BeginDateTime = sBeginDateTime.split(delim);
      String[] EndDateTime = sEndDateTime.split(delim);
      
      
      String beginDate = BeginDateTime[0]; String beginTime = BeginDateTime[1];
      String endDate = EndDateTime[0]; String endTime = EndDateTime[1];
      //endingDate = endDate; 
      //System.out.println(endingDate);
      
      String sbyear = beginDate.substring(0,4); 
      String sbmonth = beginDate.substring(4,6);
      String sbday = beginDate.substring(6,8);
      
      String seyear = endDate.substring(0,4); 
      String semonth = endDate.substring(4,6);
      String seday = endDate.substring(6,8);      
      
      endingDate = seyear+"-"+semonth+"-"+seday;
      
      System.out.println((new Integer(seyear)).intValue() + " " + (new Integer(semonth)).intValue() + " " + (new Integer(seday)).intValue());
      System.out.println((new Integer(sbyear)).intValue() + " " + (new Integer(sbmonth)).intValue() + " " + (new Integer(sbday)).intValue());
      
      int byear,bmonth,bday;
      int eyear,emonth,eday;
      
      byear = (new Integer(sbyear)).intValue(); bmonth = (new Integer(sbmonth)).intValue(); bday = (new Integer(sbday)).intValue();
      eyear = (new Integer(seyear)).intValue(); emonth = (new Integer(semonth)).intValue(); eday = (new Integer(seday)).intValue();

      
      String sbhour = beginTime.substring(0,2); 
      String sbmin = beginTime.substring(2,4);
      
      
      String sehour = endTime.substring(0,2); 
      String semin = endTime.substring(2,4);
       
      int bhour = (new Integer(sbhour)).intValue(); int bmin = (new Integer(sbmin)).intValue(); 
      int ehour = (new Integer(sehour)).intValue(); int emin = (new Integer(semin)).intValue(); 
    
      System.out.println(bhour + ":" + bmin + ":00");
      System.out.println(ehour + ":0" + emin + ":00");
    
      int endingInt = (new Integer(endTime).intValue()+500);
    

 
      DateTime start = new DateTime(byear, bmonth, bday, bhour, bmin+15, 0, 0);
      //DateTime start1 = new DateTime(byear, bmonth, bday+1, bhour, bmin+5, 0, 0);     
      DateTime start1 = start.plusDays(1);
      DateTime end = new DateTime(eyear, emonth, eday, ehour, emin, 0, 0);   
      Period nextDay = new Period(start,start1);
       
 
    
      DateTime timeinc = new DateTime(byear, bmonth, bday, bhour, bmin+15, 0, 0);
      DateTime day = new DateTime(byear, bmonth, bday, bhour, bmin+15, 0, 0);
      
      dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0]; 
      timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
      
      System.out.println("Starting time! == " + time);
      System.out.println(endingDate);
      System.out.println("Ending time! == " + endingInt);
      //while(date.indexOf(endingDate) == -1) //get the day
      while(timeinc.isBefore(end))
      {

       if((!timeinc.dayOfWeek().getAsText().equals("Saturday") && !timeinc.dayOfWeek().getAsText().equals("Sunday")))// && !isHoliday(date)) //skip weekends  
       {
 
        while(HistoricalData.getTime(timeinc) <= endingInt) //get the time
        { 
        
         dt_tokens = (timeinc.toString()).split(Tdelim);   
         timetokens = dt_tokens[1].split(dotdelim); 
         time = timetokens[0]; 
        
         System.out.println(date + " " + time); 
         dates_series.add(date + " " + time);
         
         timeinc = timeinc.plusMinutes(obs_freq);                
       
        
        }
        
        day = day.plus(nextDay);
        timeinc = new DateTime(day);
        dt_tokens = (timeinc.toString()).split(Tdelim); 
        date = dt_tokens[0];         
        //timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0]; 
       }
       else
       {
         day = day.plus(nextDay);
         timeinc = new DateTime(day);
       
         dt_tokens = (timeinc.toString()).split(Tdelim); date = dt_tokens[0];  
         timetokens = dt_tokens[1].split(dotdelim); time = timetokens[0];        
       }
       
       
       //day = day.plus(nextDay);
       //timeinc = new DateTime(day);       
       
     }
     
     n_obs = dates_series.size()-20;
     
   }  
      
   
   
  public String[] toIBForex(String par)
  {
     String[] sym_cur = new String[2];
     String delims = "[.]+";
     String[] tokens = par.split(delims);
     sym_cur[0] = new String(tokens[0].substring(0, 3));
     sym_cur[1] = new String(tokens[0].substring(3, 6));   
     
     return sym_cur;
  }   
   
   
   
   
   
  //---------- Load spreads ------------------------------------------------ 
  
  public void loadSpreads(String fx)
  {
  
    spreads = new HashMap<String, Double>();  

    //read in spread file for given currency pair 
    
    try{  
       
     FileInputStream fin = new FileInputStream(new File("spreads5min_"+fx+".dat"));
     DataInputStream din = new DataInputStream(fin);
     BufferedReader br = new BufferedReader(new InputStreamReader(din));      

     while((strline = br.readLine()) != null)
     {
     
       String[] sp = strline.split("[,]+");
       spreads.put(sp[0], (new Double(sp[1])));
       System.out.println(sp[0] + " " + sp[1]);
       if(sp.length == 0)
       {System.out.println("End of spread file"); break;}
  
     }
     br.close();
    }
    catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
  
    System.out.println("Size of spread hash:" + spreads.size());
 
  }
   

   
   public double log(double p)
   {
     if(ln_trans) { return Math.log(p);}
     else {return p;}
   }
 
  
   
   
   static int getTime(DateTime t)
   {
     int tot;
     int hour = t.getHourOfDay();
     int min = t.getMinuteOfHour();
     
     tot = hour*10000 + (min)*100;
     
     return tot; 
   }  
  
  
   boolean isHoliday(String date)
   {
     boolean holiday = false;
   
     for(int i=0;i<n_holidays;i++)
     {
       if(date.indexOf(holidays[i]) != -1)
       {holiday=true; break;} //System.out.println("It's a holiday: " + holidays[i]); break;}         
     }         
     return holiday;       
   }
  
  
}    

