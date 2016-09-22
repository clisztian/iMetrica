package ch.imetrica.regComponent;

public class ARIMAModel {

    public int[] dims;           //int array of dimensions 
    public int S; 
    public boolean n_seas,seas;  //does it have nonseasonal and seasonal components
    public double var;           //the set variance

    public boolean fix;          //any of the params fixed? 
    public boolean arfixed;
    public boolean mafixed;
    public boolean sarfixed;
    public boolean smafixed;  

    public double[] arcoefs;     //array of fixed ar params for model
    public double[] macoefs;     //array of fixed ma params for model
    public double[] sarcoefs;    //array of fixed ar params for model
    public double[] smacoefs;    //array of fixed ma params for model 

    public int n_ar, n_ma, n_sar, n_sma; 

    String model_string;        // string version 
    public int model_type;      //0 = trend, 1 = seasonal, 2 = irregular, 3 = rw
    public boolean tvreg;       // not used unless model_type = 4

    public String ht_file;       //string for file name of h_t values 
    boolean hfile;

    public ARIMAModel()
    {
       n_seas = false; seas = false; fix = false; 
       dims = new int[6];  var = 1.0;   model_type = 16;    
      
       arcoefs=null; macoefs=null; sarcoefs=null; smacoefs=null;
       arfixed=false; mafixed=false; sarfixed=false; smafixed=false;
       n_ar=0; n_ma=0; n_sar=0; n_sma=0; hfile = false; tvreg = false;
       model_string = new String("Airline");  model_type = 6; S=12;

    }
 

    public ARIMAModel(int[] dim, double _var)
    {
       var = _var; model_type = 16; dims = new int[6]; S=12;  
       System.arraycopy(dim, 0, dims, 0, 6);   
       if((dim[0] > 0 || dim[2] > 0) || dim[1] > 0)
       {n_seas = true; model_string = new String("ARIMA");}  
       if((dim[3] > 0 || dim[5] > 0) || dim[4] > 0)
       {seas = true; model_string = new String("SARIMA");}  
       
       arcoefs=null; macoefs=null; sarcoefs=null; smacoefs=null;
       arfixed=false; mafixed=false; sarfixed=false; smafixed=false;
       n_ar=0; n_ma=0; n_sar=0; n_sma=0; fix = false;  hfile = false; tvreg = false;
    }
    

    //----- present the type of model -----
    public ARIMAModel(int model, double _var)
    {
       model_type = model;  dims = new int[6]; tvreg = false;
       if(model == 0)  //--- trend ----
       {   
          n_seas = true; seas = false;
          dims[0] = 0; dims[1] = 2; dims[2] = 3;  
          dims[3] = 0; dims[4] = 0; dims[5] = 0;
          model_string = new String("Trend MA"); 
       }
       else if(model == 1)  //--- seasonal ----
       {
          n_seas = true; seas = false;
          dims[0] = 0; dims[1] = 11; dims[2] = 10;  
          dims[3] = 0; dims[4] = 0; dims[5] = 0;  
          model_string = new String("Seasonal MA");               
       } 
       else if(model == 2)  //--- irregular ----
       {
          n_seas = false; seas = false;
          dims[0] = 0; dims[1] = 0; dims[2] = 0;  
          dims[3] = 0; dims[4] = 0; dims[5] = 0; 
          model_string = new String("Irregular"); 
       }
       else if(model == 3)  //--- irregular ----
       {
          n_seas = true; seas = false;
          dims[0] = 0; dims[1] = 1; dims[2] = 0;  
          dims[3] = 0; dims[4] = 0; dims[5] = 0; 
          model_string = new String("Random Walk"); 
       }
       else if(model == 4)  //--- irregular with varying time ----
       {
          n_seas = true; seas = false;
          dims[0] = 0; dims[1] = 1; dims[2] = 0;  
          dims[3] = 0; dims[4] = 0; dims[5] = 0; 
          tvreg = true;
          model_string = new String("Time Varying TD");
       }
       else if(model == 6)  //--- airline ----
       {
          n_seas = true; seas = true;
          dims[0] = 0; dims[1] = 1; dims[2] = 1;  
          dims[3] = 0; dims[4] = 1; dims[5] = 1; 
          model_string = new String("Airline");
       }

       var = _var;  S=12;
       arcoefs=null; macoefs=null; sarcoefs=null; smacoefs=null;
       arfixed=false; mafixed=false; sarfixed=false; smafixed=false;
       n_ar=0; n_ma=0; n_sar=0; n_sma=0; fix = false;  hfile = false;   

    }

    public void setTimeVarying(boolean tv) {tvreg = tv;}

    public void setDimensions(int p, int d, int q, int P, int D, int Q) 
    {
       if((p > 0 || q > 0) || d > 0)
       {n_seas = true; model_string = new String("ARIMA");}  
       if((P > 0 || Q > 0) || D > 0)
       {seas = true;  model_string = new String("SARIMA");} 

       dims[0] = p; dims[1] = d; dims[2] = q;  
       dims[3] = P; dims[4] = D; dims[5] = Q; 
    }


    public void setVar(double _var) {var = _var;}
    public void setAllFixed(boolean _f) {fix = _f;}
    public void fixAR(boolean t) {arfixed = t;}
    public void fixMA(boolean t) {mafixed = t;}
    public void fixSAR(boolean t) {sarfixed = t;}
    public void fixSMA(boolean t) {smafixed = t;}
    //-------------------------------------------------------------
    //   Set coefficient parameters
    //------------------------------------------------------------
    public void setAR(double[] _ar) 
    { 
       n_ar = _ar.length;
       arcoefs = new double[n_ar]; 
       System.arraycopy(_ar, 0, arcoefs, 0, n_ar);  
    }
    public void setMA(double[] _ma) 
    { 
       n_ma = _ma.length;
       macoefs = new double[n_ma]; 
       System.arraycopy(_ma, 0, macoefs, 0, n_ma);  
    }
    //-------------------------------------------------------------
    //   Set coefficient parameters
    //------------------------------------------------------------
    public void setSAR(double[] _ar) 
    { 
       n_sar = _ar.length;
       arcoefs = new double[n_sar]; 
       System.arraycopy(_ar, 0, sarcoefs, 0, n_sar);  
    }
    public void setSMA(double[] _ma) 
    { 
       n_sma = _ma.length;
       macoefs = new double[n_sma]; 
       System.arraycopy(_ma, 0, smacoefs, 0, n_sma);  
    }
 
    public void setHtFile(boolean tr, String file)
    {ht_file = new String(file); hfile = tr;} 


    public String arimaToString()
    {
       int i; 
       String arimaS = new String("&arima   ");

       //---------------------------------------- seasonal --------------------------------- 
       if(model_type == 1) //---------------------------------------------------------------
       {
          arimaS = arimaS + "order = 0 11 10   diffcoefs = 11*-1. "; arimaS = arimaS + " var = " + Double.toString(var) + "\n";  
          if(n_ma > 0) 
          {
            arimaS = arimaS + "  macoefs = "; 
            for(i=0; i < n_ma-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(macoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(macoefs[n_ma-1])+"\n";            
          }
          if(mafixed){arimaS = arimaS + "  mafixed = t\n";} 

          arimaS = arimaS + "&end\n"; 
       }
       //---------------------------------------- trend ---------------------------------
       else if(model_type == 0) //-------------------------------------------------------
       {
          arimaS = arimaS + "order = 0 2 2 "; arimaS = arimaS + " var = " + Double.toString(var) + "\n";  
          if(n_ma > 0) 
          {
            arimaS = arimaS + "  macoefs = "; 
            for(i=0; i < n_ma-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(macoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(macoefs[n_ma-1])+"\n";       
          }
          if(mafixed){arimaS = arimaS + "  mafixed = t\n";} 
          arimaS = arimaS + "&end\n"; 
       }
       //---------------------------------------- irregular ---------------------------------
       else if(model_type == 2) //-----------------------------------------------------------
       {
         arimaS = arimaS + "order = 0 0 0  var = " + Double.toString(var) + " &end\n";
       }
       //---------------------------------------- random walk for time varying regressions --
       else if(model_type == 3) //random-walk 
       {
         arimaS = arimaS + "order = 0 1 0  var = " + Double.toString(var) + " &end\n";
       }
       //---------------------------------------- random walk for time varying regressions --
       else if(model_type == 4) //random-walk time varying
       {
         arimaS = arimaS + "order = 0 1 0  var = " + Double.toString(var) + " tvreg = t &end\n";
       }
       else if(model_type == 5) //white noise 
       {
         arimaS = arimaS + " var = " + Double.toString(var) + " &end\n";
       }
       //---------------------------------------- general (s)arima---------------------------------
       else //-------------------------------------------------------------------------------------
       {

          if(n_seas) {arimaS = arimaS + "order = " + Integer.toString(dims[0]) + " " + Integer.toString(dims[1]) + " " + Integer.toString(dims[2]) + " ";} 
          if(seas) {arimaS = arimaS + " sorder = " + Integer.toString(dims[3]) + " " + Integer.toString(dims[4]) + " " + Integer.toString(dims[5]) + " ";}
          arimaS = arimaS + " var = " + Double.toString(var) + "\n";

          if(n_ma > 0) 
          {
            arimaS = arimaS + "  macoefs = "; 
            for(i=0; i < n_ma-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(macoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(macoefs[n_ma-1])+"\n";    
            if(mafixed){arimaS = arimaS + "  mafixed = t\n";}    
          }
          if(n_ar > 0) 
          {
            arimaS = arimaS + "  arcoefs = "; 
            for(i=0; i < n_ar-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(arcoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(arcoefs[n_ar-1])+"\n";    
            if(arfixed){arimaS = arimaS + "  arfixed = t\n";}    
          } 
          if(n_sma > 0) 
          {
            arimaS = arimaS + "  smacoefs = "; 
            for(i=0; i < n_sma-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(smacoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(smacoefs[n_sma-1])+"\n";    
            if(smafixed){arimaS = arimaS + "  smafixed = t\n";}    
          }
          if(n_sar > 0) 
          {
            arimaS = arimaS + "  sarcoefs = "; 
            for(i=0; i < n_sar-1; i++) 
            {
              if(i == 5) {arimaS = arimaS + "\n            ";} arimaS = arimaS + Double.toString(sarcoefs[i]) + ", "; 
            }
            arimaS = arimaS + Double.toString(sarcoefs[n_sar-1])+"\n";    
            if(sarfixed) {arimaS = arimaS + "  sarfixed = t";} 
          } 
          if(fix) {arimaS = arimaS + "  fix = t ";}
          if(hfile) {arimaS = arimaS + " file = " + "'" + ht_file+ "'" + " \n";}
          arimaS = arimaS + "&end\n";
 

       }
 
       return arimaS;
       
    }


    public static void main(String args[])
    {

       double[] macoefs = {-0.737378, -0.627978, -0.430368, -0.360770, -0.219736, -0.180929, -0.088488, -0.071423, -0.020306, -0.016083};
       ARIMAModel m1 = new ARIMAModel(); 
       

       m1.setDimensions(0, 1, 1, 0, 1, 1); 
       m1.setVar(2.0);
  
       System.out.println(m1.arimaToString());

       ARIMAModel m2 = new ARIMAModel(0, 1.22);
       ARIMAModel m3 = new ARIMAModel(1, 1.22);
       m3.setMA(macoefs);
       m3.fixMA(true);
       m3.setHtFile(true,"uSim_coeffs.dat");
       ARIMAModel m4 = new ARIMAModel(2, 1.22);

  
       System.out.println(m2.arimaToString());
       System.out.println(m3.arimaToString());
       System.out.println(m4.arimaToString());

    }

}
