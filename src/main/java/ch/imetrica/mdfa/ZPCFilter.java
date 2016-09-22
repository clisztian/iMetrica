package ch.imetrica.mdfa;
import java.io.*;

public class ZPCFilter
{
  
    //---- Data and filter controls -------------

    //simPanel sim; 
    int n_obs, K, L, Lag, S, n_rep,dd;
    int K1, output, flength;  
    int order;
    int qlength;
    //---- extra Filter controls ----------------
    double lambda, expweight, cutoff, cutoff0, lambda_3;   
   
    //--- Min mdfa -----
    double zpc_min;   
    double diff_band, diff_stop;
    int normalize = 0;
 
    int samp;  //periodogram sampling    
    public double[] In;
    public double[] period_xf;
    
    //---- time series and filter data, n_rep x n_obs or K
    public double[] Gamma; 
    public double[] tseries;
    public double[] b;
    public double[] b_old;
    public double[] zpc_ar; 
    public double[] zpc_ma;
    public double[] amp; 
    public double[] fRf;
    public double[] time_delay;
    public double[] phase;

    public double[] m_zpc;          //--------multivariate zpc filtered data ---
    public double[] xf_zpc;         //--------filtered with only zpc 
    public double[] xf_hybrid;       //--------IMDFA with zpc-injected gene
  
    public double init_mod, init_mod2, init_step;    // Initial start position and initial stepsize for nelder-mead
    public double[] mods;  
    public double arg; 
    public double b0;
    public double P;
    public int p,q;
    public int i1,i2;

    //---- Filtered output and original series 
    public double[] xf; 
    public double[] x;

    public double modz1;
    public double modp1;
    public double modz2;
    public double modp2;

    ///-------------- grid stuff---------------------
    public double zpc_grid_modz;
    public double zpc_grid_modp;
    public double zpc_grid_min;
    public double[][] gridZPC;


   public ZPCFilter(int _n_obs)
   {
      n_obs = _n_obs; order = 3; 
      K = n_obs/2; K1 = K+1;
      P = .2; dd=0; p=2; q=2; i1 = 0; i2 = 0;

      lambda = 0.0; expweight = 0.0; lambda_3 = 0.0;  
      cutoff = Math.PI/5.0; cutoff0 = Math.PI/12.0; 
      
      //sim = new simPanel(300, 10, 10);

      mods = new double[4]; arg = Math.PI/5.0; b0 = 1.0;
      mods[0] = 1.5; mods[1] = 2.0; mods[2] = 3.0; mods[3] = 4.0;
      samp=150; init_mod = .5; init_mod2 = .5; init_step = .4;
   }
    
   public void set_i1(int i) { if(cutoff0 == 0.0) {i1 = i;}}
   public void set_i2(int i) {i2 = i;}
   public void setNormalize(int n) {normalize = n;}
   public void set_lambda(double l) {this.lambda = l;}
   public void set_exp(double l) {this.expweight = l;}
   public void set_cutoff(double l) {this.cutoff = l;}
   public void set_cutoff0(double l) {this.cutoff0 = l; if(cutoff0 > 0.0) {i1 = 0;}}
   public void set_lambda3(double l) {this.lambda_3 = l;}
   public void set_P(double l) {this.P = l;}
   public void setOrders(int _p, int _q) {p = _p; q = _q;}
   public void set_b0(double l) {this.b0 = l;}
   public void set_Samps(int _s) {this.samp = _s;}
   public void setInitMod(double mod, double mod2) {init_mod = mod; init_mod2 = mod2;}
   public void setInitStep(double step) {init_step = step;}

   public void setIMDFAFilter(double[] _b, int nrep, int _L)
   {
      n_rep = nrep; L = _L;
      b = new double[_b.length];
      System.arraycopy(_b, 0, b, 0, _b.length);
   }

   public void setData(double[] ts, int nrep, int nobs)
   {
      n_rep = nrep; n_obs = nobs; K= nobs/2; K1=K+1;
      tseries = new double[ts.length];
      System.arraycopy(ts, 0, tseries, 0, ts.length);  
   }


   public void computeRTSE(double[] _tseries, int N, int R)
   {
     
     int i,j,l;
     double sum=0; 
     double sum2=0;
     this.n_obs = N; this.n_rep = R; 
      
     int le = _tseries.length;
     this.tseries = new double[le];
     System.arraycopy(_tseries, 0, this.tseries, 0, le); 

     flength = n_obs - (L-1);
 
     xf_zpc = new double[flength];  
     xf_hybrid = new double[flength];  
     double[] xff = new double[flength];

     for(i=L-1;i<n_obs;i++)
     {
       sum = 0.0; 
       for(j=0;j<n_rep;j++)
       {
        for(l=0;l<L;l++)
        { 
         sum = sum + b[L*j + l]*tseries[n_obs*j + i-l];       
        }
       }    
       xff[i-(L-1)] = sum;  //System.out.println(sum);
     }

     //first two initializations
     xf_hybrid[0] = xff[0]; 
     xf_hybrid[1] = xff[1];
     xf_hybrid[2] = xff[2]; 
     xf_hybrid[3] = xff[3];
     xf_hybrid[4] = xff[4];

     xf_zpc[0] = xff[0]; 
     xf_zpc[1] = xff[1];
     xf_zpc[2] = xff[2]; 
     xf_zpc[3] = xff[3];
     xf_zpc[4] = xff[4];

     for(i=p;i<flength;i++)
     {
         sum = 0.0; sum2 = 0.0;         
         for(l=0;l<=q;l++)
         { 
           sum = sum + zpc_ma[l]*xff[i-l];                 
           sum2 = sum2 + zpc_ma[l]*tseries[(L-1) + i - l]; 
         }
         for(l=1;l<=p;l++)
         {
           sum =  -zpc_ar[l]*xf_hybrid[i-l] + sum;
           sum2 = -zpc_ar[l]*xf_zpc[i-l] + sum2;
         }
         xf_hybrid[i] = sum;   
         xf_zpc[i] = sum2;     

         //xf_hybrid[i] = -zpc_ar[1]*xf_hybrid[i-1] - zpc_ar[2]*xf_hybrid[i-2] + sum;   
         //xf_zpc[i] = -zpc_ar[1]*xf_zpc[i-1] - zpc_ar[2]*xf_zpc[i-2] + sum2;     
     }  
  }


  public void injectZPCGene2()
  {

    int i,l,j;
    double sum,sum2; double arsum,arsum2;
    flength = n_obs - (L-1);
 
    xf_zpc = new double[flength];  
    xf_hybrid = new double[flength];  
    double[] xff = new double[flength];


    for(i=L-1;i<n_obs;i++)
    {
      sum = 0.0; 
      for(j=0;j<n_rep;j++)
      {
        for(l=0;l<L;l++)
        { 
         sum = sum + b[L*j + l]*tseries[n_obs*j + i-l];       
        }
      }    
      xff[i-(L-1)] = sum;  //System.out.println(sum);
    }

    //first two initializations
    xf_hybrid[0] = xff[0]; 
    xf_hybrid[1] = xff[1];
    xf_hybrid[2] = xff[2]; 
    xf_hybrid[3] = xff[3];
    xf_hybrid[4] = xff[4];

    xf_zpc[0] = xff[0]; 
    xf_zpc[1] = xff[1];
    xf_zpc[2] = xff[2]; 
    xf_zpc[3] = xff[3];
    xf_zpc[4] = xff[4];

    for(i=p;i<flength;i++)
    {
         sum = 0.0; sum2 = 0.0; 
         arsum = 0.0; arsum2 = 0.0;        
         for(l=0;l<=q;l++)
         { 
           sum = sum + zpc_ma[l]*xff[i-l];                 
           sum2 = sum2 + zpc_ma[l]*tseries[(L-1) + i - l];           
         }

         for(l=1;l<=p;l++)
         {
           arsum =  arsum  - zpc_ar[l]*xf_hybrid[i-l];
           arsum2 = arsum2 - zpc_ar[l]*xf_zpc[i-l];          
         }
         xf_hybrid[i] = arsum+sum;   
         xf_zpc[i] = arsum2+sum2;     
    }  

    this.xf = new double[flength]; this.x = new double[flength]; 

    System.arraycopy(xf_zpc, 0, this.x, 0, flength); System.arraycopy(xf_hybrid, 0, this.xf, 0, flength);
    computeSampleIns();
  
     
  }

   public void injectZPCGene()
   {
      int i,l,j;
      double sum,sum2,arsum,arsum2;
      flength = n_obs - (L-1);
      double[] temp = new double[(n_rep-1)*flength];

      xf_zpc = new double[flength];         //--------filtered with only zpc 
      xf_hybrid = new double[flength];  
      double[] xf_zpc2 = new double[(n_rep)*flength];
      double[] xf_zpc1 = new double[(n_rep)*flength];       

      
      m_zpc = new double[n_rep*n_obs]; 
      
 
     //---- first apply I_MDFA to the raw data -------
     sum = 0.0; sum2 = 0.0; arsum = 0.0; arsum2 = 0.0;  
     for(j=1;j<n_rep;j++)
     {
       for(i=L-1;i<n_obs;i++)
       {
         sum = 0.0; 
         for(l=0;l<L;l++)
         { 
           sum = sum + b[L*j + l]*tseries[n_obs*j + i-l];       
         }
         temp[flength*(j-1) + i-(L-1)] = sum;  //System.out.println(sum);
       }          
     }
    
      //--- now apply ZPC filter to the output of mdfa, get first two
     for(j=1;j<n_rep; j++)
     { 
        xf_zpc2[flength*(j-1)] = temp[flength*(j-1)];
        xf_zpc2[flength*(j-1)+1] = temp[flength*(j-1)+1];
        xf_zpc2[flength*(j-1)+2] = temp[flength*(j-1)+2];
        xf_zpc2[flength*(j-1)+3] = temp[flength*(j-1)+3];
        xf_zpc2[flength*(j-1)+4] = temp[flength*(j-1)+4];
  
        xf_zpc1[flength*(j-1)] = temp[flength*(j-1)];
        xf_zpc1[flength*(j-1)+1] = temp[flength*(j-1)+1];
        xf_zpc1[flength*(j-1)+2] = temp[flength*(j-1)+2];
        xf_zpc1[flength*(j-1)+3] = temp[flength*(j-1)+3];
        xf_zpc1[flength*(j-1)+4] = temp[flength*(j-1)+4];

     } 

     for(j=0;j<n_rep;j++)
     {
        m_zpc[j*n_obs] = tseries[j*n_obs];
        m_zpc[j*n_obs+1] = tseries[j*n_obs+1];
        m_zpc[j*n_obs+2] = tseries[j*n_obs+2];
     }


     //----------- pre filter the data ------------------
     sum = 0.0; sum2 = 0.0; arsum = 0.0; arsum2 = 0.0;  
     for(j=0;j<n_rep;j++)
     { 
       //----- compute filtered zpc ----------------
       for(i=q;i<n_obs;i++)
       {
          sum = 0.0; sum2 = 0.0; arsum = 0.0; arsum2 = 0.0;          
          for(l=0;l<=q;l++)
          {sum = sum + zpc_ma[l]*tseries[j*n_obs + i - l];}  

          for(l=1;l<=p;l++)
          {arsum = -zpc_ar[l]*m_zpc[j*n_obs + i - l] + arsum;}   
          
          m_zpc[j*n_obs + i] = arsum + sum;
       }              
     }



     for(j=1; j < n_rep; j++)
     {       
       for(i=p;i<flength;i++)
       {
         sum = 0.0; sum2 = 0.0; arsum = 0.0; arsum2 = 0.0;          
         for(l=0;l<=q;l++)
         {
           sum2 = sum2 + zpc_ma[l]*tseries[n_obs*(j-1) + (L-1) + i - l];  
           sum = sum + zpc_ma[l]*temp[flength*(j-1) + i - l];            
         }
         
         for(l=1;l<=p;l++)
         {
           arsum =  -zpc_ar[l]*xf_zpc1[flength*(j-1)+i-l] + arsum;
           arsum2 = -zpc_ar[l]*xf_zpc2[flength*(j-1)+i-l] + arsum2;
         }
         xf_zpc1[flength*(j-1)+i] = arsum  + sum;
         xf_zpc2[flength*(j-1)+i] = arsum2 + sum2; 
       }              
     }     

     //----now sum up at each observation--------
     for(i=0;i<flength;i++)
     {
      sum=0.0; sum2 = 0.0;
      for(j=1;j<n_rep;j++)
      { 
        sum2 = sum2 + xf_zpc2[flength*(j-1)+i];
        sum = sum + xf_zpc1[flength*(j-1)+i];
      }
      xf_hybrid[i] = sum;
      xf_zpc[i] = sum2;
     }

    this.xf = new double[flength]; this.x = new double[flength]; 

    System.arraycopy(xf_zpc, 0, this.x, 0, flength); System.arraycopy(xf_hybrid, 0, this.xf, 0, flength);
    computeSampleIns();
      
   }


   //------------------ Set symmetric target filter ------------------------------
   public void set_Gamma(double[] _gamma) //---- target symmetric filter ------
   {     
      K1 = _gamma.length;
      this.Gamma = new double[K1];
      System.arraycopy(_gamma, 0, this.Gamma, 0,K1);
   }

   public void setArg(double l) {arg = l;}

   public void setMods(double mz1, double mz2, double mp1, double mp2)
   {
     mods[0] = mz1; mods[1] = mz2; mods[2] = mp1; mods[3] = mp2;
   }

   public void setPeriodogram(double[] period)
   {
      K1 = period.length;
      this.In = new double[K1];
      System.arraycopy(period, 0, this.In, 0, K1);
   }
 

   public void getModOptimizedZPC()
   {
      int k,i;
      double[] output = optimizeZPC(In, n_obs, p, q, i1, i2, Gamma, lambda, expweight, lambda_3, 
                                    cutoff0, cutoff, P, normalize, init_mod, init_mod2, init_step);

     int n_temp = output.length - 2*K1;
     
      time_delay = new double[K1]; 
      amp = new double[K1];
      phase = new double[K1]; 
      fRf = new double[2*K1];
 
      zpc_ar = new double[p+1];
      zpc_ma = new double[q+1];

      for(k=0;k<K1;k++)
      {
        amp[k] = output[k];
        phase[k] = output[K1+k];
        time_delay[k] = output[2*K1+k];
        fRf[k] = output[n_temp + k];
        fRf[K1+k] = output[n_temp + K1+k];
      }
      for(i=0;i<=p;i++)
      {zpc_ar[i] = output[3*K1+i];}
      for(i=0;i<=q;i++)
      {zpc_ma[i] = output[3*K1+p+1+i];}


      modz1 = output[n_temp - 5];
      modp1 = output[n_temp - 4];
      modz2 = output[n_temp - 3];
      modp2 = output[n_temp - 2];
      zpc_min = output[n_temp-1];
       
   }


   public void computeGridZPC(int combo, int param_set, double fixed_arg, int resolution)
   {
      int k,i; double mod1,mod2;
          
      if(combo == 0) {mod1 = 0; mod2 = 0;}
      else
      {
        if(param_set == 1) 
        {mod1 = modz1; mod2 = modp1;}
        else {mod1 = modz2; mod2 = modp2;}
      }

      double[] output = _computeGridZPC(In, n_obs, p, q, i1, i2, Gamma, lambda, expweight, lambda_3, 
                                    cutoff0, cutoff, P, normalize, init_mod, init_step,
                                    param_set, fixed_arg, mod1, mod2, resolution);

      gridZPC = new double[resolution][resolution]; 

      for(i=0;i<resolution;i++)
      {
        for(k=0;k<resolution;k++)
        {
          gridZPC[i][k] = output[resolution*i+k];
        }
      }
      
      zpc_grid_modz = output[output.length-3];
      zpc_grid_modp = output[output.length-2];
      zpc_grid_min = output[output.length-1];
       
   }




   public native double[] optimizeZPC(double[] sampleIn, int _nObs, int p, int q, int i1, int i2, double[] gamma, double lambda, double expweight,
                                        double lambda3, double _cutoff0, double _cutoff, double _P, int n, double m, double m2, double s);



   public native double[] _computeGridZPC(double[] sampleIn, int _nObs, int p, int q, int i1, int i2, double[] gamma, double lambda, double expweight,
                                        double lambda3, double _cutoff0, double _cutoff, double _P, int n, double m, double s,
                                        int c, double arg, double _modz, double _modp, int N1);

   
   static {System.loadLibrary("zpcFilter");}
   

   public static void main(String args[])
   {
      int n_obs = 300; int n_rep = 1; int L = 28; int Lag = 0; int i1 = 0; int i2 = 0;
      double lambda = 0; double expweight = 0; double cutoff0 = Math.PI/12.0; double cutoff = Math.PI/5.0;
      int K = (int)n_obs/2; int K1 = K+1; 
 
      int i,j,k,n_toks,len;
      String[] tokens;
      String strline;
      Double D;
      double val;
      double omegak;
      String delims = "[ ]+";
      //------------ Build symmetric Gamma---------------
      double[] Gamma = new double[K1];
      double[] bc = new double[L+1];
      double sum;

      bc[0] = cutoff/Math.PI; sum= bc[0];
      for(i=1;i<=L;i++)
      {bc[i] = (1/Math.PI)*Math.sin(cutoff*i)/i; sum= sum+bc[i];} 
      for(i=0;i<=L;i++)
      {bc[i] = bc[i]/(sum+(sum-cutoff/Math.PI));}     
      for(k=0; k<=K;k++)
      {
        omegak=k*Math.PI/K;	
        sum=bc[0];
        for(j=1;j<=L;j++) {
          sum = sum + bc[j]*Math.cos(omegak*j)*2.0;
        }
        if((omegak >= cutoff0) && (omegak <= cutoff))
        {Gamma[k] = 1.0;} //System.out.println(Gamma[k]);
        else 
        {Gamma[k] = 0.0;}
        //System.out.println(Gamma[k]);
      }

      double[] In = new double[K1];  
      int n_r = 5; //number of series for mdfa      
      //------------ Get data ---------------
      double[] x1 = new double[n_r*n_obs];
      double[] x = new double[n_obs];
      double[] y = new double[n_r*900]; //all the data
 
      FileInputStream fin;
      DataInputStream din;
      BufferedReader br;
     
      System.loadLibrary("zpcFilter");

      File file = new File("data2.dat");
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
         for(i=0;i<n_toks;i++) //n_toks should be 5
         {
           D = new Double(tokens[i]);
           val = D.doubleValue(); 
           y[900*i + ncount] = val;
         }
         ncount++;
       }
       din.close();
       fin.close();
       
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}

      for(i=n_obs;i<n_obs+300;i++)
      { 
         for(j=0;j<n_r;j++)
         { x1[n_obs*j + (i-n_obs)] = y[900*j + i]; }
         x[i-n_obs] = y[i]; //System.out.println(x[i-n_obs]);
      }

      L=48;
      len = n_obs - L + 1;
      n_rep = n_r; lambda = 0; i1 = 0; i2 = 0; expweight = 0.0; 
      IMDFA mdfa = new IMDFA(n_obs, n_rep, L, Lag, lambda, expweight, cutoff, i1, i2); 

      mdfa.set_cutoff(cutoff); mdfa.set_cutoff0(cutoff0);
      mdfa.set_dfaIter(false);   
      mdfa.set_Gamma(Gamma);
      mdfa.set_tseries(x1,n_obs,n_rep);
      mdfa.set_bconstraints(0,0);
      
     
      //for(i=0;i<n_obs;i++) {System.out.println(mdfa.tseries[i]);}


  
     mdfa.computeFilter(true);

     //for(k=0;k<len;k++)
     //{System.out.println(mdfa.xf[k]);}

     mdfa.computeSampleIns();   
     for(k=0;k<K1;k++) {In[k] = mdfa.period_xf[k];}
 

     //------- Set zpc filter ---------------  
     ZPCFilter zpc = new ZPCFilter(n_obs);

     zpc.setOrders(2, 2);
     zpc.set_i1(0);
     zpc.set_i2(0);
     zpc.setNormalize(0);
     zpc.set_Gamma(Gamma); 
     zpc.setArg(cutoff);
     zpc.setMods(1.5, 11.0, 2.5, 5.0);
     zpc.setPeriodogram(In);
     zpc.set_lambda(0); 
     zpc.set_exp(0); 
     zpc.set_cutoff(cutoff); 
     zpc.set_cutoff0(cutoff0); 
     zpc.set_lambda3(0.0);
     zpc.set_P(0.2); 
     zpc.set_b0(1.0); 
     zpc.setIMDFAFilter(mdfa.b, n_rep, L);
     zpc.setData(x1, n_rep, n_obs);
     zpc.setInitMod(.5,.5);
     zpc.setInitStep(.4);

     //--- Now optimize --------
     zpc.getModOptimizedZPC();
     zpc.injectZPCGene2();     

     for(k=0;k<len;k++)
     {System.out.println(mdfa.xf[k] + " " + zpc.xf_hybrid[k]);}
 
   }
 
  //--------------------- Compute sample Periodogram for x and xf

   public void computeSampleIns()
   {
     //System.out.println("computes periodogram");
     int i,j;
     int n_obsd; double mean_x,mean_xf; 
     double cdenom;  
     double sum = 0.0; double sum2 = 0.0; 
     double sumr, sumi, sumr2, sumi2;
     double[] d_seriesxf;
     double[] d_seriesx;
     int samp1 = samp+1;
 
     int dd = 0; int DD = 0;

     diff_stop = 0.0; diff_band = 0.0;
     int w0 = (int)(cutoff0*samp1/Math.PI);
     int w1 = (int)(cutoff*samp1/Math.PI);
     

     period_xf = new double[2*samp1];
   
     if((dd == 1) && (DD==1))
     { 
       n_obsd =  flength - S - 1; 
       d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd];  

       for(i=S; i < flength-1; i++)
       {
        d_seriesx[i-S] = (x[i+1] - x[i]) - (x[i+1-S] - x[i-S]);  sum = sum + d_seriesx[i-S];
        d_seriesxf[i-S] = (xf[i+1] - xf[i]) - (xf[i+1-S] - xf[i-S]);  sum2 = sum2 + d_seriesxf[i-S];
       }
       mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else if(dd == 1)
     {
      n_obsd = flength - 1; 
      d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd]; 
      for(i=1; i < flength; i++)
      {
        d_seriesx[i-1] = x[i] - x[i-1]; sum = sum + d_seriesx[i-1];
        d_seriesxf[i-1] = xf[i] - xf[i-1]; sum2 = sum2 + d_seriesxf[i-1];
      }
      mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else if(DD == 1)
     {
        n_obsd = flength - S; 
        d_seriesx = new double[n_obsd];  d_seriesxf = new double[n_obsd];    
        for(i=S; i < flength; i++)
        {
          d_seriesx[i-S] = x[i] - x[i-S]; sum = sum + d_seriesx[i-S];
          d_seriesxf[i-S] = xf[i] - xf[i-S]; sum2 = sum2 + d_seriesxf[i-S];
        }
        mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;
     }
     else
     {
       n_obsd = flength; 
       d_seriesx = new double[n_obsd];  
       d_seriesxf = new double[n_obsd];     

       for(i=0; i < flength; i++)
       {
          d_seriesx[i] = x[i]; sum = sum + d_seriesx[i];
          d_seriesxf[i] = xf[i]; sum2 = sum2 + d_seriesxf[i];
       }
       mean_x = sum/n_obsd; mean_xf = sum2/n_obsd;       
     }  
     
     period_xf[0] = mean_x*mean_x/Math.sqrt(2*Math.PI*n_obsd);
     period_xf[samp1] = mean_xf*mean_xf/Math.sqrt(2*Math.PI*n_obsd);

     if(w0 > 0) {diff_stop = period_xf[samp1];}
     else {diff_band = Math.abs(period_xf[0] - period_xf[samp1]);}
     
     //----------  Testing different transform -----------
     for(j=1;j<samp1;j++)
     {
  
       sumr = 0.0; sumi=0.0;   sumr2 = 0.0; sumi2 =0.0; 
       for(i=0;i<n_obsd;i++)
       {
        sumr = sumr + d_seriesx[i]*Math.cos(Math.PI*(i+1)*j/samp); 
        sumi = sumi + d_seriesx[i]*Math.sin(Math.PI*(i+1)*j/samp);
        sumr2 = sumr2 + d_seriesxf[i]*Math.cos(Math.PI*(i+1)*j/samp); 
        sumi2 = sumi2 + d_seriesxf[i]*Math.sin(Math.PI*(i+1)*j/samp);
       }    

       period_xf[j] = (sumr*sumr + sumi*sumi)/Math.sqrt(2*Math.PI*n_obsd);
       period_xf[j+samp1] = (sumr2*sumr2 + sumi2*sumi2)/Math.sqrt(2*Math.PI*n_obsd);

       if(dd>0) 
       {
        cdenom = ((1.0 - Math.cos(Math.PI*j/samp))*(1.0 - Math.cos(Math.PI*j/samp)) + Math.sin(Math.PI*j/K)*Math.sin(Math.PI*j/samp));        
        period_xf[j] = period_xf[j]/cdenom; period_xf[j+samp1] = period_xf[j+samp1]/cdenom;
       }
       if(DD > 0) 
       {
         cdenom = ((1.0 - Math.cos(12*Math.PI*j/samp))*(1.0 - Math.cos(12*Math.PI*j/samp)) + Math.sin(12*Math.PI*j/samp)*Math.sin(12*Math.PI*j/samp));         
         period_xf[j+samp1] = period_xf[j+samp1]/cdenom; period_xf[j] = period_xf[j]/cdenom;
       
       }
 
       if((j < w0) || (j > w1))      
       {diff_stop = diff_stop + period_xf[samp1+j];}
       else
       {diff_band = diff_band + Math.abs(period_xf[j] - period_xf[samp1+j]);}      

             
      }

    } 
    
  

}
