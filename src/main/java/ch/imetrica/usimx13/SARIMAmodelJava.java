package ch.imetrica.usimx13;

import java.io.*;
import java.util.Random;

public class SARIMAmodelJava
{

   int N;
   int burnin;
   int seed;
   int S;
   int lag;
   int m_p, m_q, m_P, m_Q, m_d, m_D;
   double[] series;
   double[] d_series;
   double[] sampleIn; 
   double[] sampleFw;
   double[] sampleDn;
   double[] sigDiags;
   double[] arima_params;
   int[] m_dims;
   boolean quarter;
   Polynomial mphi, mPhi, mtheta, mTheta, arpoly, mapoly;

   int n_params;
   double m_innvar; 

   double[] armaps;
   double[] LBV;
   double[] LK;

   double[] trendpoly;
   double[] seaspoly;
   double[] tipoly;
   double[] trendvar;
   double[] seasvar;
   double[] irrvar;
   double[] tivar;   
   int ntrend, nseas, ntip;
  
   int Nfrcst, Nx13out;  
   double[] forecasts;
   double[] x13output;
   double[] resEff;    

   //Optimization stuff
   int nOptParams; // equals complete number of params plus 1
   double[] optParams;

   //--------samples 
   double[] sampleQIf;
   boolean seats; 

   boolean td, outlier, easter, trans;
   int easterDay;

   boolean automdl;    //sets automdl
   boolean fix_model;  //sets fixed parameters
   

   
   public SARIMAmodelJava(int _N, int _burn)
   {
     N = _N; burnin = _burn; seed = 0; lag = 0; nOptParams = 0; 
     series = new double[N];
     m_dims = new int[6];
     sigDiags = new double[4];      
     LBV = new double[4];
     LK = new double[4];
     forecasts = new double[72]; //72 comes from 3x24 24 max forecast
     sampleQIf = new double[600];
     seats = true;
     quarter = false;

     trans = false; 
     td = false; 
     outlier = false;
     easter = false;
     easterDay = 1;

   }

   public double[] sampleSARIMAmodel(int seed)
   {
 

        int i,j,start,total,p,q;
   
        double sqrtinnv, innvar;
        double sum,mean;
   
        double[] samp = new double[N];
        p = arpoly.deg;
        q = mapoly.deg;

        //get the final element, the innovation variance sigma^2
        innvar = arima_params[m_p + m_P + m_q + m_Q];   
        mean = 0.0;   

        start = p+q;
        total = burnin + N + start+1;

        sqrtinnv = Math.sqrt(innvar);  
        //simulate innovation process
        double[] nt = new double[burnin + N + start+1]; 
        double[] xt = new double[burnin + N + start+1];     

        double[] eps = new double[burnin + N + start+1];       
        Random gen = new Random(seed);


        //keep the moving average of innovation

         for(j=0; j < total; j++)
         {eps[j] = mean + sqrtinnv*gen.nextGaussian();}

  
//         //initial values of x
         for(j=0; j <= start; j++)
         {xt[j] = sqrtinnv*gen.nextGaussian();}// System.out.println(xt[j]); }

        //do the moving average part
        sum = 0.0;
        for(i=start; i < total; i++)
        {
           //first sum up the past innovation process
           sum=eps[i]; 
           for(j=1; j <= q; j++)
           {sum = sum + mapoly.coef[j]*eps[i-j]; }  
           nt[i] = sum;
    
           //now sum up the past observations
           sum=0.0;
           for(j=1; j <= p; j++)
           {sum = sum - arpoly.coef[j]*xt[i-j];}
           xt[i] = sum + nt[i]; //cout << i <<"  " << xt[i] << endl; 
        }
   
        //Now give the final N observations 
        for(i=0; i < N; i++)
        {samp[i] = xt[burnin+start+i];} 
        
        return samp;
   }
 

   public void setAutomdl(boolean auto)
   {automdl = auto;}

   public void setFixParameters(boolean fix, int[] dim, double[] params)
   {
     fix_model = fix;
   }


   public void SetObservations(int _N)
   {N = _N; series = new double[N];}
   
   public void SetBurnin(int _b)
   {burnin = _b;}

   public void SetSeriesLength(int _N)
   {N = _N;}

   public void SetSeasonal(int _S)
   {S = _S;}

   public void SetSeed(int s)
   {seed = s;}
   
   public void SetInnVar(double innvar)
   {m_innvar = innvar;}

   public void SetLag(int _lag)
   {lag = _lag;}

   public void setSeats(boolean s)
   {seats = s;}
  
   public void computeSampleIFw(int n, double[] data, int model)
   {
   
     double sum, sum1, sum2;
     double isum1, isum2; 
     double mean_y; 
     int n_obsd, n_obs,length, ardeg, madeg;     
     int i,j,h,lh;
     double[] gamma; double[] temp;
     n_obs = n; int nsamp = 300;
     sampleFw = new double[300];
     sampleDn = new double[300];
     sampleIn = new double[300];
     Polynomial tarpoly;
     Polynomial tmapoly;


     tarpoly = mphi.times(mPhi);
     tmapoly = mtheta.times(mTheta);    
     ardeg = tarpoly.deg; madeg = tmapoly.deg;

     //set differenced data
     if((m_d == 1) && (m_D==1))
     { 
       n_obsd = n_obs - S - 1; 
       d_series = new double[n_obsd];    
       temp = new double[n_obs-1];
       for(i=1; i < n_obs; i++)
       {temp[i-1] = data[i] - data[i-1];}
       for(i=S; i < n_obs-1; i++)
       {d_series[i-S] = temp[i] - temp[i-S];}  
     }
     else if(m_d == 1)
     {
      n_obsd = n_obs - 1; 
      d_series = new double[n_obsd];
      for(i=1; i < n_obs; i++)
      {d_series[i-1] = data[i] - data[i-1];}
     }
     else if(m_D==1)
     {
        n_obsd = n_obs - S; 
        d_series = new double[n_obsd];    
        for(i=S; i < n_obs; i++)
        {d_series[i-S] = data[i] - data[i-S];}
     }
     else
     {
       n_obsd = n_obs; 
       d_series = new double[n_obsd];    
       for(i=0; i < n_obs; i++)
       {d_series[i] = data[i];}
     }  


     length = n_obsd; sum = 0;
     gamma = new double[length];
     for(i=0; i < length; i++)
     {sum = sum + d_series[i];}
     mean_y = sum/length;

     //compute Gamma
     //----compute the sample ACF====
     for(h = 0; h < length; h++)
     {
      sum = 0.0;
      lh = length;
      for(i = h; i < lh; i++)
      {
        sum = sum + (d_series[i] - mean_y)*(d_series[i-h] - mean_y); 
      }
      gamma[h] = sum/length;    
     }       

     //for(i=0;i<=madeg;i++)
     //{System.out.println(tmapoly.coef[i]);} 

     //========================= compute the sample spectrum of Q(I) and Q(f) ==============
     for(i=0; i < nsamp; i++)
     {    
      sum = 0; sum1 = 0; sum2 = 0; 
      isum1 = 0; isum2 = 0; 
                
      // -----------seasonal differencing-------------
      if(model == 2 || model == 6)
      {
	for(j=0;j < 12 ;j++) 
	{sum1 = sum1 + cos(j*(i*Math.PI)/nsamp); isum1 = isum1 + sin(j*(i*Math.PI)/nsamp);}  	
      }
      else if(model == 3) // ------- irregular differencing -------------
      {sum1 = 0; isum1 = 0;} 
      else // ------- trend differencing -------------
      {
	sum1 = 1 - 2*cos((i*Math.PI)/nsamp) + cos(2*(i*Math.PI)/nsamp); 
	isum1 = -2*sin((i*Math.PI)/nsamp) + sin(2*(i*Math.PI)/nsamp);
      }
      sampleDn[i] = sum1*sum1 + isum1*isum1;
      //f2wn = sum3*sum3 + isum3*isum3;
      sum1 = 0; isum1 = 0;
      //---------- f^2_W ------------

      //System.out.println(madeg);
      //System.out.println(ardeg);
      for(j=0; j <= madeg; j++)
      {
       sum2 = sum2 + tmapoly.coef[j]*cos(j*(i*Math.PI)/nsamp);
       isum2 = isum2 + tmapoly.coef[j]*sin(j*(i*Math.PI)/nsamp);
      }
      //System.out.println("i = " + i);
      for(j=0; j <= ardeg; j++)
      {
       //System.out.println(sum1);
       sum1 = sum1 + tarpoly.coef[j]*cos(j*(i*Math.PI)/nsamp);
       isum1 = isum1 + tarpoly.coef[j]*sin(j*(i*Math.PI)/nsamp);
      }

      //System.out.println(sum1 + " " + isum1);
      sampleFw[i] = (sum2*sum2 + isum2*isum2)*m_innvar/(sum1*sum1 + isum1*isum1);
         
      // -------- compute the Q(I) value ---------
      sum = 0;
      for(j=1;j < length;j++)
      {sum = sum + gamma[j]*cos(j*(i*Math.PI)/nsamp);}           
      sampleIn[i] = (2*sum + gamma[0]);  
     
  } 
  //System.out.println("New samples computed");       
}  
   


   public void SetSARIMAparams(double[] params, int[] _dim, int _n_params, int _S)
   {
     int i,total;  
     S = _S;  
     n_params = _n_params;
     arima_params = new double[n_params];   
     
     armaps = new double[n_params];
     Polynomial difd;     
   
     for(i=0; i<6; i++)
     {m_dims[i] = _dim[i]; }
   
     m_p = m_dims[0]; m_d = m_dims[1]; m_q = m_dims[2]; 
     m_P = m_dims[3]; m_D = m_dims[4]; m_Q = m_dims[5];
  
     //Test the consistency
     total = m_p + m_q + m_P + m_Q + 1;
     if(total != n_params)
     {System.out.println("Somethings wrong");}
   
     //System.out.println("Dimensions");
     //System.out.println(m_p + " " + m_q + " " + m_P + " " + m_Q);
     //Copy the parameters
     for(i = 0; i < n_params; i++)
     {arima_params[i] = params[i];}

    
     //get the final element, the innovation variance sigma^2
     m_innvar = arima_params[m_p + m_P + m_q + m_Q];   

     mphi = new Polynomial(m_p);
     mphi.coef[0] = 1.0;
     for(i=1;i<=m_p;i++)
     {mphi.coef[i] = arima_params[i-1];}

     mtheta = new Polynomial(m_q);
     mtheta.coef[0] = 1.0;
     for(i=1;i<=m_q;i++)
     {mtheta.coef[i] = arima_params[m_p + i-1];}
      
     if(m_P != 0)  
     {
      mPhi = new Polynomial(S);
      mPhi.coef[0] = 1.0;      
      mPhi.coef[S] = arima_params[m_p + m_q];
     }
     else
     {mPhi = new Polynomial(0); mPhi.coef[0] = 1.0;}  
     
     if(m_Q != 0)  
     {
      mTheta = new Polynomial(S);
      mTheta.coef[0] = 1.0;      
      mTheta.coef[S] = arima_params[m_p + m_P + m_q];
     }
     else
     {mTheta = new Polynomial(0); mTheta.coef[0] = 1.0;}  
     
     arpoly = mphi.times(mPhi);
     mapoly = mtheta.times(mTheta);     

     if(m_d == 1)
     {
       difd = new Polynomial(1);       
       difd.coef[0] = 1.0;
       difd.coef[1] = -1.0; 
       arpoly = difd.times(arpoly);
     }
        
     if(m_D == 1)
     {
       difd = new Polynomial(12);
       difd.coef[0] = 1.0;
       difd.coef[S] = -1.0;    
       arpoly = difd.times(arpoly);     
     }
   }
  
   public double[] GetSeries()
   {          
      this.series = getSARIMAmodel(N, burnin, arima_params, m_dims, n_params, S, seed);     
      return this.series;
   }

   /*public static void changeModelDimension(int[] dim)
   {
     
     try{  
           PrintWriter out = new PrintWriter(new FileWriter("x13parameters.spc"));

           out.println("series { title = \"Series\" start = 1976.jan file = \"useless.dat\"} ");
           out.println("arima { model = ("+dim[0]+ " 1 " + dim[2]+")("+dim[3]+ " 1 " + dim[5]+")12 }");
           out.println("estimate {  }"); 
           out.println("seats { noadmiss = yes }"); 
           out.close();
        } catch (IOException e) {e.printStackTrace();}

   }*/

   public void changeModelDimension(int[] dim, boolean _trans, boolean _outlier, boolean _easter, boolean _td, int ea)
   {
      trans = _trans; 
      td = _td; 
      outlier = _outlier;
      easter = _easter;
      easterDay = ea;

     try{  
           PrintWriter out = new PrintWriter(new FileWriter("x13parameters.spc"));

           out.println("series { title = \"Series\" start = 1976.jan file = \"useless.dat\"} ");
           out.println("arima { model = ("+dim[0]+ " 1 " + dim[2]+")("+dim[3]+ " 1 " + dim[5]+")12 }");
           if(trans) out.println("transform{ function = auto }");
           if(outlier) out.println("outlier{  }");
           if(td && easter) out.println("regression{ variables = (td    easter["+ea+"]) }");
           else if(td)
           {out.println("regression{ variables = td }");} 
           else if(easter)
           {out.println("regression{ variables = easter["+ea+"] }");}
           out.println("estimate {  }"); 
           out.println("seats { noadmiss = yes }"); 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
   }

   public void changeModelDimension(int[] dim)
   {


     try{  
           PrintWriter out = new PrintWriter(new FileWriter("x13parameters.spc"));

           if(!quarter) {out.println("series { title = \"Series\" start = 1976.jan file = \"useless.dat\"} ");}
           if(quarter) {out.println("series { title = \"Series\" period = 4 start = 1976.1 file = \"useless.dat\"} ");}

           out.println("arima { model = ("+dim[0]+ " 1 " + dim[2]+")("+dim[3]+ " 1 " + dim[5]+")12 }");
           if(trans) out.println("transform{ function = auto }");
           if(outlier) out.println("outlier{ }");
           if(td && easter) out.println("regression{ variables = (td    easter["+easterDay+"]) }");
           else if(td)
           {out.println("regression{ variables = td }");} 
           else if(easter)
           {out.println("regression{ variables = easter["+easterDay+"] }");}
           out.println("estimate {  }"); 
           out.println("seats { noadmiss = yes }"); 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
   }


   public double[] GetDiags(double[] data, int[] dims, double[] params, int n_params, int lag,int obs, int model)
   {

     Nx13out = 4+n_params+4+4+72+900 + 3*obs;
     x13output = new double[Nx13out];     
     
     x13output = getSigExDiagnostics(data,dims,params,n_params, lag, obs, model, seats);
     return x13output;
     //System.out.println("exiting sigExDiag...\n ");
   }

   // Meth = method for optimization; simplex = 0, BFGS = 1
   public double[] computeMinimumModel(int[] dgp, double[] dgpParams, int[] fnd)
   {
      nOptParams = fnd[0] + fnd[2] + fnd[3] + fnd[5] + 2;
      optParams = new double[nOptParams];

      optParams = computeMinimumSarima(dgp, dgpParams, fnd); 
      return optParams;
   } 

   public double[] computeEfficacies(int N, double[] data, int[] tdims, double[] tparams, 
                                     int n_tparams, int[] t_dim, int n_ptparams, int model)
   {
      resEff = new double[n_ptparams+8+1+1200+6];
      resEff = getEfficacies(N, data, tdims, tparams, n_tparams,t_dim, n_ptparams, model, lag);  
      return resEff;
   }

    public void initializeUseless(int n)
    {
       int i;
       try{  
           PrintWriter out = new PrintWriter(new FileWriter("useless.dat"));
           for(i=0; i < n; i++) {out.println(0.0);}
 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
    }



   /*-------------------------------------------------------------
      Native Methods
   --------------------------------------------------------------*/

   
   public native double[] getSigExDiagnostics(double[] data, int[] dims, 
                                              double[] params, int n_params, int lag, 
                                              int obs, int model, boolean seats);

   public native double[] computeMinimumSarima(int[] dgp, double[] dgpParams, int[] fnd);

   public native double[] getEfficacies(int N, double[] data, int[] tdims, double[] tparams, int n_tparams, 
                                                              int[] t_dim, int n_ptparams, int model, int lag);

   public native double[] getSARIMAmodel(int N, int burnin, double[] arima_params, int[] m_dim, 
                                         int n_params, int S, int seed);

   static
   {
     System.loadLibrary("SARIMAmodelJava");
   }


   
   

   static public double cos(double r) {return Math.cos(r);}
   static public double sin(double r) {return Math.sin(r);}
}      
   
   
   



