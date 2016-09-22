package ch.imetrica.dataControl;

import javax.swing.*;

import ch.imetrica.bayesCronos.HeavyModel;
import rcaller.RCaller;
/*----------------------------------------------------------------
 
   Java Interface for HighFrequency and RGL

  
----------------------------------------------------------------*/


public class JHighFreq 
{


  String[] asset_names;
  
  public int n_assets;
  int n_obs, n_days;
  String[] frequency;   //seconds, minutes, hours, days
  int freq_option; 
  
  String[] period; 
  String[] align_secs;
  public String[] kernels;
  String[] make_returns;
  String[] correlate;  
  int n_steps; int n_fore;  

  
  
  JPanel realizedVolPanel;
  public JAsset[] assets;
  RCaller caller;

  
  
  public JHighFreq()
  {
  
     make_returns = new String[2]; 
     make_returns[0] = "makeReturns = FALSE"; make_returns[1] = "makeReturns = TRUE";
  
     frequency = new String[6];
     frequency[0] = "'seconds'"; frequency[1] = "'minutes'"; frequency[2] = "'hours'"; frequency[3] = "'days'";
     frequency[4] = "'weekly'"; frequency[5] = "'monthly'";  
  
     kernels = new String[12];     
     kernels[0] = "'Rectangular'"; kernels[1] = "'Bartlett'"; kernels[2] = "'Second'";              
     kernels[3] = "'Epanechnikov'"; kernels[4] = "'Cubic'"; kernels[5] = "'Fifth'";               
     kernels[6] = "'Sixth'"; kernels[7] = "'Seventh'"; kernels[8] = "'Eighth'";              
     kernels[9] = "'Parzen'";  kernels[10] = "'TukeyHanning'"; kernels[11] =  "'ModifiedTukeyHanning'";

     align_secs = new String[6];
     align_secs[0] = "'seconds'"; align_secs[1] = "'minutes'"; align_secs[2] = "'hours'"; align_secs[3] = "'days'";    
     
     correlate = new String[2]; 
     correlate[0] = "cor = FALSE"; correlate[1] = "cor = TRUE";
     caller = new RCaller();     
     caller.setRscriptExecutable("/usr/bin/Rscript");
     
  }
  
  public int getNDays() {return n_days;}
  
  /*-------------------------------------------------------------------------------
  
    Input: 
    
    instrum - An array of symbol names, (i.e. GOOG.O AAPL.O BAC..) in RDat form
    from,to    - starting date in the form CCYY-MM-DD 
    int freq   - (0, seconds), (1,minutes), (2,hours), (3,days)
    int period - period of freq (i.e. 1 min, 3 min, 1 hour, etc.)
    boolean    - get logreturns of data
    boolean    - get volume of data
  
  ---------------------------------------------------------------------------------*/
  
  public void getMarketDataMult(String[] instrum, String from, String to, int freq, int period, boolean vol)
  {
  
     int i;
     int n_assets = instrum.length;
     String[] aligned_asset = new String[n_assets];
     double[] real_series;
     assets = new JAsset[n_assets];
   
     String closes;
     String cassets = "c("; String aassets = "c(";
     for(i=0;i<n_assets-1;i++) {cassets = cassets + "'" + instrum[i] + "',"; aassets = aassets + instrum[i] + ",";} //wrap name with quotes
     cassets = cassets +  "'" + instrum[n_assets-1] + "')"; 
     aassets = aassets + instrum[n_assets-1] + ")";

     String d_from = "from = '"+from+"'"; String d_to = "to = '"+to+"'"; 

     // 1) upload the symbols
     String getSymbols = "getSymbols.FI("+cassets+", " + d_from +", " + d_to +")";
  
     // 3) set desired frequency
     if(freq == 0)
     {
       for(i=0;i<n_assets;i++)
       {
         aligned_asset[i] = 
           instrum[i] + "<-align.time(to.period("+instrum[i]+"[,5:6],period="+
             frequency[freq]+",k="+period+",name=NULL),"+period+")['T09:30/T16:00',]";     
       }
     }
     else if(freq == 1) //case for seconds, minutes, hours
     {     
       for(i=0;i<n_assets;i++)
       {
        aligned_asset[i] = 
          instrum[i] + "<-align.time(to.period("+instrum[i]+"[,5:6],period="+
           frequency[freq]+",k="+period+",name=NULL),"+period +"*60)['T09:30/T16:00',]";
       }
     }
     else if(freq == 2)
     {
       for(i=0;i<n_assets;i++)
       {
        aligned_asset[i] = 
          instrum[i] + "<-align.time(to.hourly(" + instrum[i] + "[,5:6],name=NULL),3600)['T09:30/T16:00',]";
       }     
     }
     else if(freq == 3)
     {
       for(i=0;i<n_assets;i++)
       {
        aligned_asset[i] = 
          instrum[i] + "<-to.daily("+instrum[i]+"['T09:30/T16:00',5:6],drop.time=TRUE,name=NULL)";
       }
     }
     
  
     String na = "(cbind("; String name_vector = "c("; 

     for(i=0; i < n_assets-1; i++)
     {     
        //---gets closing data ------
        na = na + "Cl("+instrum[i] + "),";
        name_vector = name_vector + "'" + instrum[i] + "',";
     }
     na = na + "Cl("+instrum[n_assets-1]+")))";
     name_vector = name_vector + "'" + instrum[n_assets-1] + "')"; 
     closes = "closes <- na.locf" + na;
     name_vector = "names(closes) <- " + name_vector;  
  

    
     String volume = "volumes<-na.locf(cbind(";
     String name_vector2 = "c("; String instString2 = "";
     String name_vector3 = "c(";
     String name_vector0 = "c(";
     for(i=0;i<n_assets-1;i++)
     {
       instString2 = instString2 + "Vo(" + instrum[i] + "),";
       name_vector2 = name_vector2 + "'assetv"+i+"',";
       name_vector3 = name_vector3 + "'assetOrig"+i+"',";
       name_vector0 = name_vector0 + "'assetLog"+i+"',";
     } 
     instString2 = instString2 + "Vo(" + instrum[n_assets-1] + ")))";
     name_vector2 = name_vector2 + "'assetv"+(n_assets-1)+"')";
     name_vector3 = name_vector3 + "'assetOrig"+(n_assets-1)+"')";
     name_vector0 = name_vector0 + "'assetLog"+(n_assets-1)+"')";
     
     volume = volume+instString2;
     name_vector2 = "names(volumes) <- " + name_vector2;  
     name_vector3 = "names(orig) <- " + name_vector3;
     name_vector0 = "names(closes) <- " + name_vector0;
  
    try
    {

     //RCaller caller = new RCaller();
     //caller.setRscriptExecutable("/usr/bin/Rscript");
     caller.cleanRCode();
     caller.getRCode().addRCode("require (Runiversal)");
     caller.getRCode().addRCode("require (FinancialInstrument)");   
     caller.getRCode().addRCode("require (highfrequency)"); 
     caller.getRCode().addRCode("loadInstruments('/home/lisztian/Data/instruments.rda')");  
     caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/Data/sec',use_identifier='X.RIC',extension='RData')");
     //caller.getRCode().addRCode("loadInstruments('/media/My Passport/Market/instruments.rda')"); 
     //caller.getRCode().addRCode("setSymbolLookup.FI('/media/My Passport/Market/sec',use_identifier='X.RIC',extension='RData')");
     
     //indexTZ(GOOG.O) <- "America/New_York"
     //GOOG.O<-GOOG.O[!is.na(GOOG.O$Trade.Price)]
     
     //----- Call the symbols -------------------------
     //addInst = "getSymbols.FI(" + instString + ",from='" + date1 + "', to='2012-06-19')"; 
     caller.getRCode().addRCode(getSymbols);
     
     for(i=0;i<n_assets;i++)
     {caller.getRCode().addRCode("indexTZ("+instrum[i]+")<-\"America/New_York\"");}
     

     //----- Call the align --------------------------
     for(i=0;i<n_assets;i++) {caller.getRCode().addRCode(aligned_asset[i]);}

      //----- define the matrix ------------------------
      caller.getRCode().addRCode(closes);
      caller.getRCode().addRCode(volume);

      //----- define the name_vector--------------------

      caller.getRCode().addRCode(name_vector);
      caller.getRCode().addRCode(name_vector2);
      

      //---- logreturns ----------------------------------
      caller.getRCode().addRCode("closes <- log(closes)");
      caller.getRCode().addRCode("orig <- closes");
      caller.getRCode().addRCode(name_vector3);     
      caller.getRCode().addRCode("closes <- diff(closes)");
      caller.getRCode().addRCode(name_vector0);
      //---- the (log) closing price-----------------
      
      //caller.getRCode().addRCode("names(orig) <- 'orig'"); 
      
      caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)");        
      caller.getRCode().addRCode("closes[1,] <- 0");
      caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
      caller.getRCode().addRCode("volume_list<-lapply(as.list(volumes), coredata)"); 
      caller.getRCode().addRCode("all_list<-c(orig_list,series_list,volume_list)");
      caller.runAndReturnResult("all_list");
  
  
      for(i=0;i<n_assets;i++)
      {        
        real_series = caller.getParser().getAsDoubleArray("assetOrig"+i);
        assets[i] = new JAsset(real_series.length);  
        assets[i].setName(instrum[i]); 
        
        assets[i].setPrice(real_series);        
        
        real_series = caller.getParser().getAsDoubleArray("assetLog"+i);
        assets[i].setLogReturn(real_series);
        
        if(vol)
        {
         real_series = caller.getParser().getAsDoubleArray("assetv"+i);
         assets[i].setVolume(real_series);              
        }    
      }          
 
    }
    catch(Exception e){System.out.println(e.toString());} 
  
  }
  


  /*-----------------------------------------------------------
   
   Input: instrum - String of instrument names in rda form (i.e. instrum = {GOOG.O, AAPL.O, BAC}
          from, to - String for starting date to ending data in CCYY-MM-DD form
          freq - Compute the realized measures on 0) seconds 1) minutes 2) hours
          kern - choice of kernel used 
          period - the subgrid period (i.e. how many seconds, minutes, hours to align to. 
          lags - number of lags used in the realized measures  

   
  ------------------------------------------------------------*/

  public void getRealizedVolatility(String[] instrum, String from, String to, int freq, int kern, int period, int lags, boolean pr)
  {
  
     int i,j;
     n_assets = instrum.length;
     String[] aligned_asset = new String[n_assets];
     double[] real_series;
     assets = new JAsset[n_assets];
   
     String[] rvs = new String[n_assets];
   
     String cassets = "c("; String aassets = "c(";
     for(i=0;i<n_assets-1;i++) {cassets = cassets + "'" + instrum[i] + "',"; aassets = aassets + instrum[i] + ",";} //wrap name with quotes
     cassets = cassets +  "'" + instrum[n_assets-1] + "')"; 
     aassets = aassets + instrum[n_assets-1] + ")";

     String d_from = "from = '"+from+"'"; String d_to = "to = '"+to+"'"; 

     // 1) upload the symbols
     String getSymbols = "getSymbols.FI("+cassets+", " + d_from +", " + d_to +")";
    
    try
    {

     
     //caller.setRscriptExecutable("/usr/bin/Rscript");
     caller.cleanRCode();
     caller.getRCode().addRCode("require (Runiversal)");
     caller.getRCode().addRCode("require (FinancialInstrument)");   
     caller.getRCode().addRCode("require (highfrequency)"); 
     //caller.getRCode().addRCode("loadInstruments('/home/lisztian/Data/instruments.rda')");  
     //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/Data/sec',use_identifier='X.RIC',extension='RData')");
     caller.getRCode().addRCode("loadInstruments('/media/My Passport/Market/instruments.rda')"); 
     caller.getRCode().addRCode("setSymbolLookup.FI('/media/My Passport/Market/sec',use_identifier='X.RIC',extension='RData')");
     
     //indexTZ(GOOG.O) <- "America/New_York"
     //GOOG.O<-GOOG.O[!is.na(GOOG.O$Trade.Price)]
     
     //appleKern3<-rKernelCov(AAPL.O$Trade.Price,kernel.type = kernels[11], kernel.dofadj = TRUE, period = 5, align.by ="minutes", align.period=5, makeReturns=TRUE);
     
     
     //----- Call the symbols -------------------------
     //addInst = "getSymbols.FI(" + instString + ",from='" + date1 + "', to='2012-06-19')"; 
     caller.getRCode().addRCode(getSymbols);

     for(i=0;i<n_assets;i++)
     {
       String nas = instrum[i] + "<-" + instrum[i] + "[!is.na(" + instrum[i] + "$Trade.Price)]"; //AAPL.O<-AAPL.O[!is.na(AAPL.O$Trade.Price)]
       caller.getRCode().addRCode(nas);
     }
     
     for(i=0;i<n_assets;i++)
     {caller.getRCode().addRCode("indexTZ("+instrum[i]+")<-\"America/New_York\"");}

     String all_list = ""; 
     if(pr){all_list = "all_list<-c(orig_list,";}
     else {all_list = "all_list<-c(";}
     
     //--- get daily data -----------------------
     for(i=0;i<n_assets;i++)
     {
        aligned_asset[i] = 
          "asset"+i+"<-to.daily("+instrum[i]+"['T09:30/T16:00',5],drop.time=TRUE,name=NULL)";       
        
        caller.getRCode().addRCode(aligned_asset[i]);
        caller.getRCode().addRCode("closes"+i+" <- log(Cl(asset"+i+"))");
        
        if(pr && i==0) //get price data too for first asset
        {
           caller.getRCode().addRCode("orig <- closes0[,1]");
           caller.getRCode().addRCode("names(orig) <- 'orig'"); 
           caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)");
        }
        
        caller.getRCode().addRCode("closes"+i+" <- diff(closes"+i+")");
        caller.getRCode().addRCode("closes"+i+"[1,]<-0");
        caller.getRCode().addRCode("names(closes"+i+")<-'assetLog"+i+"'");
        caller.getRCode().addRCode("closesList"+i+"<-lapply(as.list(closes"+i+"), coredata)");
        all_list = all_list+"closesList"+i+",";
     }
       

/*
kernel.type = "rectangular",   # Kernel name (or number)
kernel.param = 1,              # Kernel parameter (usually lags)
kernel.dofadj = TRUE,          # Kernel Degree of freedom adjustment
align.by="seconds",            # Align the tick data to [seconds|minutes|hours]
align.period = 1,              # Align the tick data to this many [seconds|minutes|hours]
cts = TRUE,                    # Calendar Time Sampling is used
makeReturns = FALSE,           # Convert to Returns  
*/
    
     //---- delete nas-------------------------
     for(i=0;i<n_assets;i++)
     {
       //String nas = instrum[i] + "<-" + instrum[i] + "[!is.na(" + instrum[i] + "$Trade.Price)]"; //AAPL.O<-AAPL.O[!is.na(AAPL.O$Trade.Price)]
       String mark = instrum[i] + "<-" + instrum[i] + "['T09:30/T16:00',]";
       //caller.getRCode().addRCode(nas);
       caller.getRCode().addRCode(mark);
       
       String rv = "rv"+i+"<-rKernelCov("+instrum[i]+"$Trade.Price,kernel.type ="+
             kernels[kern]+", kernel.param="+lags+",kernel.dofadj = FALSE, align.by ="+frequency[freq]+
             ", align.period="+period+", cts=TRUE, makeReturns=TRUE)";


       caller.getRCode().addRCode(rv);
       caller.getRCode().addRCode("names(rv"+i+")<-'rv"+i+"'");
       
       rvs[i] = "rv_list"+i;       
       caller.getRCode().addRCode("rv_list"+i+"<-lapply(as.list(rv"+i+"), coredata)");
       
     }
     
     
     for(i=0;i<n_assets-1;i++)
     {
      all_list = all_list + "rv_list"+i+",";
     }
     all_list = all_list + "rv_list"+(n_assets-1)+")";
     caller.getRCode().addRCode(all_list);
  
     caller.runAndReturnResult("all_list");
     //caller.runAndReturnResult("rv_list0");
         

         
     for(i=0;i<n_assets;i++)
     {   
        real_series = caller.getParser().getAsDoubleArray("assetLog"+i);  
        assets[i] = new JAsset(real_series.length);  
        
        assets[i].setLogReturn(real_series);    

        real_series = caller.getParser().getAsDoubleArray("rv"+i);   
        n_days = assets[0].n_days;
        //annualize the RVol
        for(j=0;j<real_series.length;j++) {real_series[j] = Math.sqrt(256.0)*real_series[j];}
        
        assets[i].setRV(real_series);     
     }
     
     if(pr)
     {
       real_series = caller.getParser().getAsDoubleArray("orig");
       assets[0].setPrice(real_series); 
     }     
     
     
     n_days = assets[0].n_days;
 
    }
    catch(Exception e){System.out.println(e.toString());} 
  
  }  
  
  
  public void computeHeavyModel(int nd, int n_f, int n_s)
  {

     
     int seed = 10;
     int n_days = nd;
     HeavyModel heavy = new HeavyModel(n_days, seed);  
     n_fore = n_f; n_steps = n_s;

     int n_rep = 2; 
     //--- initial values------------
     double w1 = .15; 
     double w2 = .05;
     double lambda = .1;
     double alpha = .2;
     double alpha_R = .4;
     double beta = .7;
     double beta_R = .55;


     //--- set forecast dimensions -------------------       
     heavy.setForecastDimensions(n_fore, n_steps);
     
     //--- set (simulation or initialization) parameter values -------------
  
     double[] series = new double[n_days*n_rep];
  
     heavy.setParameterValues(w1, w2, alpha, alpha_R, lambda, beta, beta_R); 
     heavy.setData(n_days, n_rep, series);
     heavy.setRetrackParameter(0);
     heavy.estimateHeavyModel(); 




  }
  
  
  
   public static void main(String args[])
   {
      int i;
      JHighFreq hf = new JHighFreq();
      
      int freq = 1; int period = 5; int lags = 3; int kern = 2;

      String[] instrums = new String[2];
      instrums[0] = "GOOG.O"; instrums[1] = "AAPL.O"; 
      String from = "2012-03-01";
      String to = "2012-06-19";
      
      //hf.getMarketDataMult(instrums, from, to, 2, 0, true);
      hf.getRealizedVolatility(instrums, from, to, freq, kern, period, lags, false);
   
      //simPanel.plotData(hf.assets[0].realized_vol,hf.assets[0].n_days); 

      for(i=0;i<hf.assets[0].n_days;i++)
      {
        System.out.println(hf.assets[0].realized_vol[i] + " " + hf.assets[0].log_return[i]*hf.assets[0].log_return[i]);
      } 

   }
 
 
}
