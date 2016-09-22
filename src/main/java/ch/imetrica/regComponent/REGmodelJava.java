package ch.imetrica.regComponent;

/*!
Copyright (C) 2016 Christian D. Blakely
This file is part of iMetrica, a free-software/open-source application
for interactive graphical econometric analysis - http://imetricablog.com/

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

import java.io.*;
import java.util.ArrayList;

public class REGmodelJava
{


   public int N;
   public int burnin;
   public int seed;
   public int S;
   public int lag;
   public String date_string;               // starting date of data 

   public double[] series;                  // raw data 

   public int n_arimas;                     // number of arima models
   public ArrayList<ARIMAModel> models;     // array of arima model structures
   public ArrayList<ARIMAModel> tv;  

   public boolean check, acf, pacf, hist;   // booleans of the check variable
   public int maxlag;                       // max lag of these 
   
   public boolean transform;                // turn on transform
   public double t_power;                   // transform power

   public boolean regression;               // turn on regression
   
   public String adjust;                    // adjust keyword
   public boolean smooth_est;               // smoothing estimations
 
   public boolean forecast;                 // forecast sets  
   public int maxforecast;                  // maxforecast = 24
   public double fore_prob;                 // probability range for forecast
   
   public boolean easter;                   // easter day regression effects
   public int east_days;                    // number of easter days
   public String east_model_n;                 // to which model does easter go 0 < n < n_models

   public boolean trading_day;              // turn on trading day
   public boolean trade_string;             // turn on 'f' or 't'
   public String tmodel_n;                   // which models are activate with training day

   public int tv_start;
   public boolean tvreg;                    //turn on time varying
   public boolean outlier;
   public String estim; 

   public boolean cmpntreg;                  // component reg values
   public int[] cmpntreg_val;                // values for the regcomponet in string
   public String cmpntreg_string;            // string of values 

   public boolean constant;                  // set constant regression effect true
   public String const_string;               // value of const f or t

   public boolean cmpntlom;                  // use compntlon separately 
   public String lomModel;                   // which model is it attributed

   public String ht_file;                    // file name for h_t series

   double[] regOut;                          // Output by regcmpnt, a ncmpnt x (N + n_fore) array 
   int Ncmpnt;                               // number of total components
   int n_fore;                               // number of forecast observations = 24
   int NFore;                                // total number of observations N + 24; 


   double[] usimCmpnts;                    // all the sig ext components + forecasts
   double[] usimH;                         // all the scaling factors + forecasts
   double[] usimXbeta;                     // all the regression coponents + forecasts
   double[] usimXbetaH;                    // all the regression coponents*H + forecasts
   double[] usimQ;                         // mean squared errors  + forecasts
   double[] usimQbar;                      // mean squared error variances  + forecasts 

   double[] y_frcstm;                      // aggregate data forecasts - mid
   double[] y_frcstl;                      // aggregate data forecasts - low
   double[] y_frcsth;                      // aggregate data forecasts - high
   double[] y_xB_frcst;                    // aggregate data forecasts - regression


   // ---------------- Constructor ----------------------------
   
   //----------------------------------------------------------
   
   public REGmodelJava(int _N, int _burn)
   {
     N = _N; burnin = _burn; seed = 0; lag = 0;
     series = new double[N]; initializeUseless(N);
     maxforecast = 24; n_fore = 24; S = 12; fore_prob = .50;

     check = false; acf = false; pacf = false; hist = false; maxlag = 48;
     n_arimas = 0; estim = "t"; forecast = true; 

     models = new ArrayList<ARIMAModel>();  Ncmpnt = 0;
     regression = false;
     setStartDate(1980,1);  
   } 

   // ---------------- Set the data series properties ----------------------------
   
   public void defaultModel(double[] ser)
   {
      setSeries(ser, ser.length);
      ARIMAModel m = new ARIMAModel(6, 1.0);          // --- airline model with variance 1.0
      models.add(m);
      regression = false;                             // --- no regression 
      n_arimas = models.size();
   }

   public void cleanSlate()
   {
       models.clear(); 
       regression = false;
       n_arimas = 0;   
   }

   public void deleteModel(int i)
   {
     models.remove(i);
     n_arimas = models.size();
   }

   public boolean addTV(ARIMAModel model)
   {
      int i; tv_start = models.size(); boolean res = true;
      if(tv_start < 5)
      {
        for(i=0; i<6; i++) {models.add(model);} tvreg = true;  
      }
      else  
      {System.out.println("Way too many models. Cut down a bit"); res = false;}    
      n_arimas = models.size();  
      return res;
   }
 
   public void deleteTV()
   { 
      int i;
      for(i=tv_start; i < tv_start+6; i++) {models.remove(tv_start);}
      tvreg = false;
      n_arimas = models.size();
   }  

      
   //----------------------------------------------------------    
   public void setSeries(double[] ser, int _N) {setNObs(_N); series = new double[_N]; System.arraycopy(ser, 0, series, 0, N);}
   public void setNObs(int _n) {N = _n; initializeUseless(N);}
   public void setSeasonal(int _n) {S = _n;}
   public void setStartDate(int year, int month)
   {date_string = new String(Integer.toString(year)+","+Integer.toString(month));}


   // ---------------- Add Models to model queue ----------------------------
   
   public void addModel(ARIMAModel s)
   {
      if(n_arimas < 10)
      {models.add(s); n_arimas++;}
      else
      {System.out.println("Currently only 10 SARIMA models total are allowed for modeling your data.");} 
   }

   public void addModel(ARIMAModel s, int i)
   {
      if(i < 10) {models.add(i,s); n_arimas++;}
      else {System.out.println("Currently only 8 SARIMA models total are allowed for modeling your data, quit being so greedy");} 
   }   

   public ARIMAModel getModel(int i)
   {
     return models.get(i);        
   }

   public void replaceModel(int i, ARIMAModel m)
   {
      models.remove(i);
      models.add(i,m); n_arimas = models.size();
   }

   public void deleteAllModels()
   {models.clear(); n_arimas = 0;}


   //---------------- Regression Effect settings------------------------------------------
 
   public void setConstant(boolean t, String ft) {constant = t; const_string = new String(ft); regression = true;} 
  
   public void setCmpntLom(boolean tr, String i) {cmpntlom = tr; lomModel = new String(i); regression = true;}
 
   public void setRegTransforms(int[] _cmptreg)
   {cmpntreg_val = new int[n_arimas]; System.arraycopy(_cmptreg, 0, cmpntreg_val, 0, n_arimas);}

   public void setPowerTransform(boolean tr, double _t) {transform = tr; t_power = _t;}

   public void setTimeVarying(boolean tv) {tvreg = tv; regression = true;}

   public void setHtFile(String file) {ht_file = new String(file);} 
   
   public void setEaster(boolean east, int _d, String modelx) {easter = east; east_days = _d; east_model_n = modelx; regression = true; cmpntreg = east;}

   public void setTradingDay(boolean _t, String d1, String d2, String d3, String d4, String d5, String d6)
   {trading_day = _t; tmodel_n = d1 + " " + d2 + " " + d3+ " " + d4+ " " + d5+ " " + d6 + " "; cmpntreg = _t; regression = true;}

   public void setRegression(boolean _t) {regression = _t;}

   public void setCmpntRegString()
   {
     cmpntreg = true; 
     cmpntreg_string = new String("cmpntreg = "); 

     if(trading_day)
     {cmpntreg_string = cmpntreg_string + tmodel_n;} 
 
     if(easter)
     {cmpntreg_string = cmpntreg_string + east_model_n + " ";}
   }

   //---------------- Set the forecasting and smoothing variables ------------------------
   public void setForecast(boolean _f, int _d, double _p)
   {forecast = _f; maxforecast = _d; fore_prob = _p;}

   public void setEstim(boolean _t) {if(_t) {estim = "t";} else {estim = "f";}}   
   
   public void setCheck(boolean _t) {check = _t;}

   //-----------------Set up the NML file to change model-----------------
   public void setNMLFile()
   {

    String arima_string; int i;
    String reg = new String("&regression ");

    //------ turn on automatic smoothing,estimation,forecasting
    forecast = true; smooth_est = true;

    //------ compute regression string --------------
    if(regression)
    {     
     if(constant) {reg = reg + " const = " + const_string + " ";} else {reg = reg + " const = f ";}
     if(trading_day) {reg = reg + " td = t ";}  else {reg = reg + " td = f ";}
     if(easter) {reg = reg + "easter = " + Integer.toString(east_days) + " ";}
     if(cmpntlom) {reg = reg + lomModel;}
     if(cmpntreg) {setCmpntRegString(); reg = reg + cmpntreg_string;}
     reg = reg + " &end\n";
    }

    try{  
           PrintWriter out = new PrintWriter(new FileWriter("reguSimCmpnt.nml"));

           //---- Create the header for series read-in, useless.dat -----------
           out.println("&series start = " + date_string + " period = " + S); 
           out.println("  title = 'Series'");
           out.println("  file = 'useless.dat'");
           out.println("&end\n");

           if(transform) {out.println("&transform power = " + t_power + " &end\n");}
           if(regression) {out.println(reg);}
 
           //----- Go through all the arima models -----------
           for(i=0; i < n_arimas; i++)
           {
              arima_string = (models.get(i)).arimaToString();
              out.println(arima_string);
           }
           
           out.println("&estimate estim = " + estim + " prtiter = f maxiter = 100 armacorr = f &end\n");

           if(check) {out.println("&check  acf = t  pacf = t  hist = t  maxlag = " + Integer.toString(maxlag) + "&end\n");}
           if(smooth_est) {out.println("&smooth estimate = t &end\n");}
           if(forecast) {out.println("&forecast maxlead = 24 probability = "+fore_prob+" &end\n");}
          
 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
   }


   //-------------------------- Compute the Regcomponent Model with given .nml file --------

   public void computeRegComponent()
   {

     if(n_arimas > 0 || regression)
     {

      int i,t,rows,le; 
      NFore = n_fore+N;
      int f_start; 
      //------ set series ----------
      //setSeries(data,data.length); 
      //------ call regcomponent --- 
      regOut = getRegCmpnt(series,N); 
      le = regOut.length;

      Ncmpnt = (int) regOut[le-1]; 
      rows = Ncmpnt*NFore;      
      f_start = 6*rows;

      usimCmpnts = new double[Ncmpnt*NFore];
      usimH = new double[Ncmpnt*NFore];
      usimXbeta = new double[Ncmpnt*NFore]; 
      usimXbetaH = new double[Ncmpnt*NFore]; 
      usimQ = new double[Ncmpnt*NFore];
      usimQbar = new double[Ncmpnt*NFore];

      y_frcstm = new double[24];
      y_frcstl = new double[24];
      y_frcsth = new double[24];
      y_xB_frcst = new double[24];
      
      for(i = 0; i < Ncmpnt; i++)
      {
        // ------ get data first -----
        for(t=0; t < NFore; t++)
        {
           usimCmpnts[NFore*i + t] = regOut[0*rows + (NFore*i + t)];
           usimH[NFore*i + t]      = regOut[1*rows + (NFore*i + t)];
           usimXbeta[NFore*i + t]  = regOut[2*rows + (NFore*i + t)];
           usimXbetaH[NFore*i + t] = regOut[3*rows + (NFore*i + t)];
           usimQ[NFore*i + t]      = regOut[4*rows + (NFore*i + t)];
           usimQbar[NFore*i + t]   = regOut[5*rows + (NFore*i + t)];
        }
      }

      for(i=0;i<24;i++)
      {
        y_frcstl[i]  = regOut[f_start + 24 + i]; 
        y_frcstm[i]     = regOut[f_start + 48 + i]; 
        y_frcsth[i] = regOut[f_start + 72 + i]; 
        y_xB_frcst[i]  = regOut[f_start + 96 + i];     
      }

     }
     else
     {System.out.println("Cowardly trying to model data with nothing. \nYou need at least one ARIMA or regression component.");} 

     
   }

   // ----------statics and calls to the native functions ------------------------------------
   public native double[] getRegCmpnt(double[] data, int _N);
   public static native void plotData2(double[] data, double[] data2, int _N);
   public static native void plotData(double[] data, int _N);
   public static void plotData2(double[] data, double[] data2)
   {plotData2(data,data2,data.length);}
   public static void plotData(double[] data)
   {plotData(data,data.length);}
   static {System.loadLibrary("REGmodelJava");}
   public static void initializeUseless(int n)
   {
       int i;
       try{  
           PrintWriter out = new PrintWriter(new FileWriter("useless.dat"));
           for(i=0; i < n; i++) {out.println(0.0);}
 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
   }




   public static void main(String args[])
   {

       int n_obs,i;
       double[] series; 
       double[] h;
       double[] macoefs = {-0.737378, -0.627978, -0.430368, -0.360770, -0.219736, -0.180929, -0.088488, -0.071423, -0.020306, -0.016083};
       
       System.loadLibrary("REGmodelJava");

       //-------- Get series -------------------------- //for(i=0;i<n_obs;i++) System.out.println(series[i]); 
       series = readDataFile12("tseries1.dat");
       n_obs = series.length;      

 
       // 1) Setup first model ----------------       
       ARIMAModel seas = new ARIMAModel(1,60.0);      //---- seasonal MA model with var = 60
       seas.setMA(macoefs);                           //---- set MA coefficients
       seas.fixMA(true);                              //---- Make them fixed
 
       // 2) Setup second model ----------------
       ARIMAModel model1 = new ARIMAModel();          //---- an ARIMA (0,2,1) model
       model1.setDimensions(0,2,1,0,0,0); 
       model1.setVar(200.0); double ma[] = {.8};      //set variance and initial ma coeff
       model1.setMA(ma); 

       // 2) Setup third model ----------------
       ARIMAModel model2 = new ARIMAModel(5,200.0);   //---- white noise model = 5
   
       //Start new regComp engine with n_obs 
       REGmodelJava regcmp = new REGmodelJava(n_obs,100); 
       
       regcmp.setPowerTransform(true,0.0);
       //regcmp.setEaster(true, 10, 1);
       //regcmp.setTradingDay(true, 1, 1, 1, 1, 1, 1);
       regcmp.setForecast(true, 24, .50);
       regcmp.setCmpntRegString();
        
       //----- Now add models --------------------------
       regcmp.addModel(seas); 
       regcmp.addModel(model1);
       regcmp.addModel(model2);      
       regcmp.setNMLFile();
       //regcmp.computeRegComponent(series);      //---- add series ----


      
       int NFore = regcmp.NFore;
       double[] data = new double[NFore];
       double[] data2 = new double[NFore];
       
       //------------------------------ -----------------------------------------------
       System.out.println("-----Start new reg-----");
       series = readDataFile12("nr055.dat");
       h = readDataFile12("h_coeff.dat");
       plotData(h);

       n_obs = series.length;       
       System.out.println("Nobs = " + n_obs);
       REGmodelJava reg = new REGmodelJava(n_obs,100);
       
       //---- set regression and other effects -------------------
       reg.setSeries(series, n_obs);
       reg.setRegression(true);
       reg.setForecast(true, 24, .50);

       ARIMAModel m1 = new ARIMAModel();          // --- airline model ------
       m1.setDimensions(0, 1, 1, 0, 1, 1); 
       m1.setVar(.016565);
  
       ARIMAModel m2 = new ARIMAModel();          // --- set AR(2) model with fixed coeffs
       m2.setDimensions(2,0,0,0,0,0);
       m2.setVar(.34488); double[] ar = {.600,.246}; 
       m2.setAR(ar);  
       m2.setAllFixed(true);
       m2.setHtFile(true, "h_coeff.dat");

       //--- add models ----------------------
       reg.addModel(m1); 
       reg.addModel(m2);

       reg.setNMLFile();
       System.out.println("----- Enter -----");
       reg.computeRegComponent();
       System.out.println("----- Done -----");
       NFore = reg.NFore;
       data = new double[NFore];
       data2 = new double[NFore];
 
            
       //for(i=0; i < NFore; i++) {data[i] = regcmp.usimCmpnts[NFore*0 + i];}
       //for(i=0; i < NFore; i++) {data2[i] = regcmp.usimCmpnts[NFore*1 + i];}

       for(i=0; i < NFore; i++) {data[i] = reg.usimCmpnts[NFore*0 + i]; }
       for(i=0; i < NFore; i++) {data2[i] = reg.usimCmpnts[NFore*1 + i]; }

       plotData(series);
       plotData(data); 
       plotData(data2);


    }


   public static double[] readDataFile12(String fileName) 
   {
     String line = ""; int i,n_obs;
     ArrayList<Double> data = new ArrayList<Double>();//consider using ArrayList<int> 
     double[] series; 
       

     try 
     {
        FileReader fr = new FileReader(fileName);
   	    BufferedReader br = new BufferedReader(fr);//Can also use a Scanner to read the file
   	    while((line = br.readLine()) != null) 
        {         
         data.add(Double.valueOf(line).doubleValue());
        }
   	    br.close();
     }
     catch(FileNotFoundException fN) { fN.printStackTrace(); }
     catch(IOException e) { System.out.println(e); }

     n_obs = data.size(); 
     series = new double[n_obs];

     for(i=0;i<n_obs;i++) 
     {series[i] = (double) data.get(i);}

     return series;
   }

}


