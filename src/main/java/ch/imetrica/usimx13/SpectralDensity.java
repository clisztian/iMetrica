package ch.imetrica.usimx13;

import java.io.*;

import ch.imetrica.bayesCronos.Cronos;
import rcaller.RCaller;
/*------------------------------------------------------

  Java Class for creating various spectral density estimations

  Input is time series data, n_rep time series that are stationary or require differencing
    
--------------------------------------------------------*/


public class SpectralDensity
{


   int n_rep, N, K, K1; //number of series, length of series, spectral coeffs
   int dd, DD; //differencing operators
   int M; //number of multitaper levels
   double h; //scaling factor 
   int smooth_func;

   double[] tseries;
   public double[] mod_spec;
   public double[] arg_spec; 
   double[] real_spec;
   double[] imag_spec; 
   public double[] period_x;
   double[] orig_spec;
   double[] arma_arg_spec;
   // ---- different spectral density options ------
   boolean arma; 
   boolean multitaper;
   boolean quadratic;
   boolean smooth;
   boolean emd; 
   boolean gaussianize; 
   boolean unit_phase; 
   boolean periodogram;
   boolean arma_phase; 
   boolean gaussianized;
   
   //------ ARMA stuff 
   int ar_order, ma_order; 
   double[] ar_params;
   double[] ma_params;
   double innvar; 


   //---- multitaper stuff
   double W;

   
   

   public SpectralDensity(int n, int nrep)
   {
	 
	 System.loadLibrary("SpectralDens");  
     N = n; n_rep = nrep; K = N/2; K1 = K+1;

     smooth_func = 0; arma = false; multitaper = false; smooth = false; 
     emd = false; gaussianize = false; unit_phase = false; periodogram = true;
     arma_phase = false; gaussianized = false;
     
     M = 5; W = .01;
     mod_spec = new double[K1*n_rep]; 
     arg_spec = new double[K1*n_rep];

     period_x = new double[n_rep*K1];

   }
      
   //---- setting some parameters ----------------------
   public void setNrep(int nrep) {n_rep = nrep;}
   public void setNobs(int n) {N = n; K=N/2; K1 = K+1;}
   public void set_dd(int _dd) {dd = _dd;}
   public void set_DD(int _dd) {DD = _dd;}   
   public void setScale(double _h) {h = _h;}
   public void setMultitaperLevels(int m) {M = m;}
   public void setSpectralConcentration(double w) {W = w;}   

   public void engageARMA(boolean t) {arma = t;}             // create arma model of spectral density
   public void setMultitaper(boolean t) {multitaper = t;}    // turn on/off multitapering
   public void setSmoothing(boolean t) {smooth = t;}         // turn on/off smoothing
   public void setEMD(boolean t) {emd = t;}                  // setEMD filtering of data
   public void setGaussianize(boolean t) {gaussianize = t;}  // Gaussianize the data (using Lambert W)
   public void setUnitPhase(boolean t) {unit_phase = t;}     // use unit phase component (unbiased)
   public void setQuadratic(boolean t) {quadratic = t;}      // use quadratic multitaper method
   public void setSmoothFunction(int s) {smooth_func = s;}
   public void setPeriodogram(boolean t) {periodogram = t;}
   
   @SuppressWarnings("deprecation")
public double[] Gaussianize(double[] series) // Gaussianize nrep series of length n. Series is n_rep x n matrix    
   {
     double[] results = null;
     if(series.length == N)
     {     
        
      try
      {

       RCaller caller = new RCaller();
       caller.setRscriptExecutable("/usr/bin/Rscript");
       caller.cleanRCode();
       
       //---- call libraries --------------
       caller.getRCode().addRCode("require (Runiversal)");
       caller.getRCode().addRCode("require (LambertW)");   
   
       caller.addDoubleArray("x", series); //set series
       caller.getRCode().addRCode("gaussx<-Gaussianize(x,type='hh',method='IGMM')"); //Gaussianize
       caller.getRCode().addRCode("my.all<-list(gaussianized=gaussx)");
       
       caller.runAndReturnResult("my.all");
       results = caller.getParser().getAsDoubleArray("gaussianized");
  
      }
      catch(Exception e){System.out.println(e.toString());}
     }
     
     return results; 
   }


   public void computeARMASpectralDensity(int pc, int qc)
   {
       int n_obs = N;
       int i,j,k; double real, imag;
       double rar_sum = 0.0; double iar_sum = 0.0;
       double rma_sum = 0.0; double ima_sum = 0.0;

       int model = 1; int method = 1; int d = 0; int MLE = 1;
       int n_fsteps = 1; int nsims = 1; double innvar;

       double[] tser = new double[N];
       System.arraycopy(tseries,0,tser,0,n_obs);
       
       arma_arg_spec = new double[K1*n_rep];
       Cronos armam = new Cronos(n_obs, model, method);
       armam.setData(tser);
       armam.setNForecastSteps(n_fsteps);
       armam.setNPredictiveSims(n_fsteps, nsims); 
       armam.setARMA_Params(pc,qc,d); 
       armam.computeARIMAModel(MLE);
 
       innvar = armam.innvar;
       int ar_p = armam.ar_params.length; int ma_q = armam.ma_params.length;
 
 
       for(k=0;k<K1;k++)
       {
        rar_sum = 1.0; iar_sum = 1.0;
        for(i=0;i<ar_p;i++)
        {
         rar_sum = rar_sum - armam.ar_params[i]*Math.cos((i+1.0)*Math.PI*k/K);
         iar_sum = iar_sum - armam.ar_params[i]*Math.sin(-(i+1.0)*Math.PI*k/K);
        }
    
        rma_sum = 1.0; ima_sum = 1.0;
        for(i=0;i<ma_q;i++)
        {    
          rma_sum = rma_sum + armam.ma_params[i]*Math.cos((i+1.0)*Math.PI*k/K);
          ima_sum = ima_sum + armam.ma_params[i]*Math.sin(-(i+1.0)*Math.PI*k/K);
        }
        
        rma_sum = rma_sum*innvar; ima_sum = ima_sum*innvar; 
        rar_sum = 2*rar_sum; iar_sum = 2*rar_sum; 
           
        real = (rma_sum*rar_sum + ima_sum*iar_sum)/(rar_sum*rar_sum + iar_sum*iar_sum);
        imag = (ima_sum*rar_sum - rma_sum*iar_sum)/(rar_sum*rar_sum + iar_sum*iar_sum);
    
        mod_spec[k] = Math.sqrt(real*real + imag*imag);
        arma_arg_spec[k] = Math.atan(imag/real);
           
        mod_spec[K1+k] = mod_spec[k];
        arma_arg_spec[K1+k] = arma_arg_spec[k];           
       }  
 
 
 
       for(j = 2; j < n_rep; j++)
       {
       
          System.arraycopy(tseries,n_obs*j,tser,0,n_obs);
          armam.setData(tser);
          armam.computeARIMAModel(MLE);
         
          //arma.ar_params, arma.ma_params, arma.innvar
         
          for(k=0;k<K1;k++)
          {
           rar_sum = 1.0; iar_sum = 1.0;
           for(i=0;i<ar_p;i++)
           {
            rar_sum = rar_sum - armam.ar_params[i]*Math.cos((i+1.0)*Math.PI*k/K);
            iar_sum = iar_sum - armam.ar_params[i]*Math.sin(-(i+1.0)*Math.PI*k/K);
           }
    
           rma_sum = 1.0; ima_sum = 1.0;
           for(i=0;i<ma_q;i++)
           {    
            rma_sum = rma_sum + armam.ma_params[i]*Math.cos((i+1.0)*Math.PI*k/K);
            ima_sum = ima_sum + armam.ma_params[i]*Math.sin(-(i+1.0)*Math.PI*k/K);
           }
        
           rma_sum = rma_sum*innvar; ima_sum = ima_sum*innvar; 
           rar_sum = 2*rar_sum; iar_sum = 2*rar_sum; 
           
           real = (rma_sum*rar_sum + ima_sum*iar_sum)/(rar_sum*rar_sum + iar_sum*iar_sum);
           imag = (ima_sum*rar_sum - rma_sum*iar_sum)/(rar_sum*rar_sum + iar_sum*iar_sum);
    
           mod_spec[K1*j+k] = Math.sqrt(real*real + imag*imag);
           arma_arg_spec[K1*j+k] = Math.atan(imag/real);
           
          }      
        }
         
        for(j=0;j<n_rep;j++)
        {         
          for(k=0;k<K1;k++)
          {
            if(unit_phase)
            {arg_spec[K1*j+k] = 0.0;}
            else if(arma)
            {arg_spec[K1*j+k] = arma_arg_spec[K1*j+k];}           
          }
        }
      
   }
   
   public void changePhaseInformation(int c)
   {
      int j,k;
   
      if(c == 0) {arma_phase = true; unit_phase = false;}
      else if(c == 1) {arma_phase = false; unit_phase = false;}
      else {unit_phase = true;}
   
        for(j=0;j<n_rep;j++)
        {         
          for(k=0;k<K1;k++)
          {
            if(unit_phase)
            {arg_spec[K1*j+k] = 0.0;}
            else if(arma_phase)
            {arg_spec[K1*j+k] = arma_arg_spec[K1*j+k];}
            else
            {arg_spec[K1*j+k] = orig_spec[K1*j+k];}
          }
        }   
   }
   
   
   public void setARMAModel(boolean t)
   {arma = t;}

   public void computeSpectralDensity(double[] _tseries, int n, int nrep)
   {
      //--- set lenghts
      
      tseries = new double[_tseries.length];
      System.arraycopy(_tseries,0,tseries,0,_tseries.length);
      gaussianized = false;
      int k,i;
      N = n; n_rep = nrep; K = N/2; K1 = K + 1; 
      int K2 = 2*K1;

      double[] out = new double[n_rep*(2*K1)];


      //----- if we gaussianize ------------------
      if(gaussianize && !gaussianized)
      {
        System.out.println("Gaussianize data");
        double[] gseries = new double[n*nrep];
        double[] temp = new double[n];

        for(k=0;k<n_rep;k++)
        {

          System.arraycopy(tseries,k*N,temp,0,N);          
          temp = Gaussianize(temp);
          System.arraycopy(temp,0,gseries,k*N,N); 
        }
        
        System.arraycopy(gseries,0,tseries,0,gseries.length); 
        gaussianized = true;
        
        out = spectral_density(gseries, N, n_rep, smooth, h, smooth_func, multitaper, M, quadratic, W);
      }
      else
      {out = spectral_density(tseries, N, n_rep, smooth, h, smooth_func, multitaper, M, quadratic, W);}
    

      mod_spec = new double[n_rep*K1];
      arg_spec = new double[n_rep*K1];
      orig_spec = new double[n_rep*K1];

      for(i=0;i<n_rep;i++)
      {
        for(k=0;k<K1;k++)
        {
          mod_spec[K1*i + k] = out[K2*i + k];
          arg_spec[K1*i + k] = out[K2*i + K1 + k];
        }
      }

      System.arraycopy(arg_spec,0,orig_spec,0,n_rep*K1);
      
      if(periodogram)
      {computePeriodogram();}

   }

   
   public void computeSpectralDensity(int n, int nrep)
   {
      //--- set lenghts      
      int k,i;
      N = n; n_rep = nrep; K = N/2; K1 = K + 1; 
      int K2 = 2*K1;

      double[] out = new double[n_rep*(2*K1)];

      if(gaussianize && !gaussianized)
      {
        double[] gseries = new double[n*nrep];
        double[] temp = new double[n];

        for(k=0;k<n_rep;k++)
        {

          System.arraycopy(tseries,k*N,temp,0,N);          
          temp = Gaussianize(temp);
          System.arraycopy(temp,0,gseries,k*N,N); 
        }
        
        System.arraycopy(gseries,0,tseries,0,gseries.length); 
        gaussianized = true;
      }      

      //----- if we gaussianize ------------------
      if(!arma)
      {
       out = spectral_density(tseries, N, n_rep, smooth, h, smooth_func, multitaper, M, quadratic, W);
           
       mod_spec = new double[n_rep*K1];
       arg_spec = new double[n_rep*K1];


       for(i=0;i<n_rep;i++)
       {
         for(k=0;k<K1;k++)
         {
          mod_spec[K1*i + k] = out[K2*i + k];
          if(!unit_phase)
          {arg_spec[K1*i + k] = out[K2*i + K1 + k];}
        
          if(arma_phase && arma)
          {arg_spec[K1*i + k] = arma_arg_spec[K1*i + k];}
         }
       }

      }
      if(periodogram)
      {computePeriodogram();}

   }   
   

   public void computePeriodogram()
   {
     int k,i,j;
     double sumr,sumi;
     double[] temp = new double[N]; 
     period_x = new double[K1*n_rep];

     for(k=0;k<n_rep;k++)
     {

       System.arraycopy(tseries,k*N,temp,0,N);

       for(j=0;j<K1;j++)
       {  
         sumr = 0.0; sumi=0.0; 
         for(i=0;i<N;i++)
         {
           sumr = sumr + temp[i]*Math.cos(Math.PI*(i+1)*j/K); 
           sumi = sumi + temp[i]*Math.sin(Math.PI*(i+1)*j/K);
         }    

         period_x[k*K1 + j] = (sumr*sumr + sumi*sumi)/(Math.PI*N);
       }
     }
   }

   //------ Native functions -------------------------------------------
   
   /*----- main spectral density function ------

      Input
      data: a n_rep x n matrix of the time series data
         n: length of time series
      nrep: the total number of time series
  l_smooth: apply smoothing to final spec dens estimate
         h: the localizer for the smoothing
    s_func: type of smoothing function 0-3
      l_mt: apply multitapering to periodogram (using sine basis)
         M: number of levels (or if quadratic tapering, number of sepian functions)
    l_quad: apply multitapering with quadrature (using optimized basis)
         V: the sepian functions M x N matrix
         v: the eigenvalues M vector
         W: spectral concentration 
      

      Output
        output of length nrep x (K1 + K1) (the mod (K1) and arg (K1))  

   --------------------------------------------*/
   public native double[] spectral_density(double[] data, int n, int nrep, boolean l_smooth, double h, int s_func, boolean l_mt, int M,
                                           boolean l_quad, double W); 



   static {System.loadLibrary("SpectralDens");}




   public static void main(String args[])
   {
 
      double val; 

      int n_obs = 385; int n_rep = 3; int K = n_obs/2; int K1 = K+1; 

      double[] tseries = new double[n_obs*n_rep];

      int n_toks;
      String[] tokens;
      String strline;
      Double D;
      String delims = "[ ]+";

      FileInputStream fin;
      DataInputStream din;
      BufferedReader br;
     
      System.loadLibrary("SpectralDens");


      File file = new File("hf_goog.dat");
      int ncount = 0;
      try{

      fin = new FileInputStream(file);
      din = new DataInputStream(fin);
      br = new BufferedReader(new InputStreamReader(din)); 

       while((strline = br.readLine()) != null)
       {
         tokens = strline.split(delims); 
         n_toks = tokens.length; 
         if(n_toks == 0)
         {System.out.println("End of file"); break;}
  
         //D = new Double(tokens[n_toks-1]);
         //val = D.doubleValue();

          D = new Double(tokens[1]);
          val = D.doubleValue(); 
           
          tseries[n_obs*0 + ncount] = val;
          tseries[n_obs*1 + ncount] = val;
   
          D = new Double(tokens[2]);
          val = D.doubleValue(); 
        
          tseries[n_obs*2 + ncount] = val;

         
          ncount++;
        }
       din.close();
       fin.close();       
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}


      System.out.println(ncount + " = " + n_obs);
 
      SpectralDensity sd = new SpectralDensity(n_obs, n_rep);


      //---- start with smoothing ------------
      double h = .1;
      sd.setScale(h);
      sd.setSmoothing(true);
      sd.setSmoothFunction(2);
      sd.setMultitaper(true);

      //sd.setGaussianize(true);

      int J = 15;
      double W = J*1.0/n_obs;
      System.out.println("W = " + W);
      sd.setQuadratic(true);
      sd.setSpectralConcentration(W);
      //---- compute periodogram on -------------
      sd.setPeriodogram(true);


      double[] series = new double[n_obs];
      System.arraycopy(tseries,0,series,0,n_obs);
      //plotData(series,n_obs);
      System.arraycopy(tseries,n_obs,series,0,n_obs);
      //plotData(series,n_obs);
      System.arraycopy(tseries,2*n_obs,series,0,n_obs);      
      //plotData(series,n_obs);

      //----- compute spectral density and period_x
      sd.computeSpectralDensity(tseries, n_obs, n_rep);


      double[] mod = new double[K1];
      double[] spec = new double[K1];
      double[] peri = new double[K1]; 
      System.arraycopy(sd.arg_spec,2*K1,spec,0,K1);
      System.arraycopy(sd.period_x,2*K1,peri,0,K1);
      System.arraycopy(sd.mod_spec,2*K1,mod,0,K1);

      //for(i=0;i<K1;i++) {System.out.println(spec[i]);}



   }
}