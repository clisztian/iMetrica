package ch.imetrica.usimx13;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import java.text.*;

public class SARIMAmodelCanvas extends JPanel //implements Runnable
{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//-------SARMA Model stuff------
    SARIMAmodelJava sarima;
    int nObs;
    
    //Thread runner = null; 
    int burnin;
    int S, lag, nx13out;

    int minObs;

    //---- sweep stuff----------------
    double[] means,sdevs;
    double error;

    String[] dates; 
    //----- Sliding spans stuff ---------------------
    int nbackObs = 0;
    int nforeObs = 0;
    int nObsSpan; //= nObs - nbackObs - nforeObs; 
    boolean slidingSpan = false;
    double[] span_data;
    int slideStart = 0; 
    int slideEnd,nfore;


    boolean compute = false;
    public Font mono;
    double[] t_series;  
    double[] real_series;
    double[] sim_params;
    double[] mle_params;
    
    double[] lbv;
    double[] lkhs;
    double[] forecasts;
    double[] sampleQIf;
    double[] seassig;
    double[] trendsig;
    double[] cyclesig;
    double[] sasig;
 
    //------polynomials for building filters -----

    double[] trendpoly;
    double[] seaspoly;
    double[] tipoly;
    double[] trnsdenom;
    double[] trnspoly;

    int n_poly;  //combination of filter coefficients output
    int n_params,seed,model;
    int[] m_dim;    
    int m_p, m_d, m_q, m_P, m_D, m_Q; 
    double m_innvar;

    Polynomial mphi, mPhi, mapoly, arpoly, mtheta, mTheta;

    double[] diagnostics;
    double[] x13output;
    int track_pos_t;
    int track_pos_x;
    boolean plot_tracker=false;
    int tracker=0;
    int nObsex;
    String value;
    double[] mleparams;    
    double[][] records;
    //-------Graphics stuff
    Ellipse2D ellipse; 
    Rectangle2D rectangle;
    Graphics2D g2d;
    Cursor curCursor;
    int height, width;
    double aParam, T0, T1, X0, X1, Y0, Y1;
    double dataMax, dataMin, dataNorm;
    boolean dataSetup;
    boolean plotForecasts; 
    boolean plotSeries;
    boolean plotSeas;
    boolean plotTrend;
    boolean plotSA;
    boolean plotCycle;
    DecimalFormat df,df2;

    int year,month;   

    int pressed_t;
    boolean specPlots; 
    boolean seats;
    Color seasonalColor;
    Color trendColor; 
    Color saColor; 
    Color cycleColor;  
    Color seriesColor;
    Color seriesColor2;
    Color myGray;
    Color forecol; 

    Color highlight; 
    int hlight_indicator; 

    //-----added 5-25-11----------
    boolean sim_series;
    //-----added 6-14-11----------
    boolean trans, td, easter, outlier;
    int eaDays;
    
    JTextField diag1,diag2,diag3,diag4;
    JTextField lb1,lb2,lb3,lb4;
    JTextField aic1,aic2,aic3,aic4;
    JTextField mle1,mle2,mle3,mle4,mle5,mle6,mle7;
  
  
    //---- sweep panel ----------------
    private JLabel aiccLabel;
    private JTextField aiccText;
    private JLabel bicLabel;
    private JTextField bicText;
    private JButton compSweepButton;
    private JScrollBar foreBar;
    private JScrollBar minObsBar;
    private JLabel foreErrorLabel;
    private JTextField foreErrorText;
    private JTextField foreHorizText;
    private JLabel foreLabel;
    private JLabel meanAicLabel;
    private JLabel meanLabel;
    private JLabel meanLabel1;
    private JTextField meanPhiText;
    private JTextField meanPhiText1;
    private JTextField meanphi1Text;
    private JTextField meanphi1Text1;
    private JLabel meanphi2Label;
    private JLabel meanphi2Label1;
    private JTextField meanphi2Text, minObsText;
    private JTextField meanphi2Text1;
    private JLabel sPhiLabel, minObsLabel;
    private JLabel sPhiLabel1;
    private JLabel sphi1Label;
    private JLabel sphi1Label1;
    public JPanel sweepPanel;    
    
  
  
  

    public SARIMAmodelCanvas(int w, int h, int _nObs, int _burnin, int _seed)
    {
       this.height = h; this.width = w; 
       this.nObs = _nObs; this.burnin = _burnin;
       this.seed = _seed;
       sarima = new SARIMAmodelJava(_nObs, _burnin);
   
       
       
       S = 12; minObs = 60;  
       n_poly = 8+27+32+5+36;  //combination of filter coefficients output
       sim_params = new double[3];
       mle_params = new double[3];
       lbv = new double[4];
       lkhs = new double[4];
       forecasts = new double[72];
       sampleQIf = new double[900];
       seassig = new double[nObs];
       trendsig = new double[nObs];
       cyclesig = new double[nObs];
       sasig = new double[nObs];
       t_series = new double[nObs];
       real_series = new double[nObs];
 
       trendpoly = new double[8];
       seaspoly = new double[27];
       tipoly = new double[32];
       trnsdenom = new double[5];
       trnspoly = new double[36];

       value = new String("");
       m_dim = new int[6];
       m_dim[0] = 0; m_dim[1] = 1; m_dim[2] = 1; 
       m_dim[3] = 0; m_dim[4] = 1; m_dim[5] = 1; 
       n_params = 3;

       sim_params[0] = -0.6; sim_params[1] = -0.6; sim_params[2] = 1.0;
       diagnostics = new double[4];
       diagnostics[0] = 0.0; diagnostics[1] = 0.0;
       diagnostics[2] = 0.0; diagnostics[3] = 0.0;
       lag = 0; model = 6;
       nx13out = 4+3+2+4+72+900+3*nObs + n_poly;
       x13output = new double[nx13out];
       mono  = new Font("Monospaced", Font.PLAIN, 12);
       
       setBackground(Color.BLACK);
       //setBackground(Color.WHITE);
       setPreferredSize(new Dimension(w, h));
       dataSetup = false; plotForecasts = false; plotSeries = true;
       plotSeas = false; plotTrend = false; plotSA = false; plotCycle = false;

       seasonalColor = new Color(255,47,47);
       trendColor = new Color(255,139,249); 
       saColor = new Color(241,231,17); 
       cycleColor = new Color(100,229,109);  
       seriesColor = new Color(66,191,213);
       seriesColor2 = new Color(66,218,222);
       forecol = new Color(0,80,82);
       myGray = new Color(67,61,63);
       df = new DecimalFormat("##.##"); 
       df2 = new DecimalFormat("#.###");
       specPlots = false; seats = true;
       highlight = new Color(203,201,188); 
       hlight_indicator = -1;

       trans = false; outlier = false; td = false; eaDays = 1; easter = false; 
       //=====  are we simulating series or uploading =======       
       sim_series = true;
       addMouseMotionListener(new MyArtMouseMotionListener());  
       addMouseListener(new MyArtMouseListener());
       nforeObs = 0; nbackObs = 0; nObsSpan = 0; slideEnd = 0; slideStart = 0;
       curCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);

       mleparams = new double[7]; year = 1980; month = 1;
       
    }

    public void setDiagnosticsFields(JTextField d1, JTextField d2, JTextField d3, JTextField d4)
    {
       diag1 = d1; diag2 = d2; diag3 = d3; diag4 = d4;       
    }
    
    public void setLBFields(JTextField d1, JTextField d2, JTextField d3, JTextField d4)
    {
      lb1=d1; lb2 = d2; lb3 = d3; lb4 = d4;
    }
    
    public void setAICFields(JTextField d1, JTextField d2, JTextField d3, JTextField d4)
    {
      aic1=d1; aic2 = d2; aic3 = d3; aic4 = d4;
    }
    
    public void setMLEFields(JTextField d1, JTextField d2, JTextField d3, JTextField d4, JTextField d5, JTextField d6, JTextField d7)
    {
      mle1 = d1; mle2 = d2; mle3 = d3; mle4 = d4; mle5 = d5; mle6 = d6; mle7 = d7;
    }    
    
    public void setTexts()
    {
      
      int j = 0;
      mleparams = new double[7];
    
      if(m_p >= 1) {mleparams[0] = mle_params[j]; j++;}
      else{mleparams[0] = 0.0;}
      if(m_p == 2) {mleparams[1] = mle_params[j]; j++;}
      else{mleparams[1] = 0.0;}
      if(m_q >= 1) {mleparams[2] = mle_params[j]; j++;}
      else{mleparams[2] = 0.0;}
      if(m_q == 2) {mleparams[3] = mle_params[j]; j++;}
      else{mleparams[3] = 0.0;}
      if(m_P == 1) {mleparams[4] = mle_params[j]; j++;}      
      else{mleparams[4] = 0.0;}
      if(m_Q == 1) {mleparams[5] = mle_params[j]; j++;}       
      else{mleparams[5] = 0.0;}
      mleparams[6] = mle_params[j];
    
      
      mle1.setText("" + df.format(mleparams[0]));
      mle2.setText("" + df.format(mleparams[1]));
      mle3.setText("" + df.format(mleparams[2]));
      mle4.setText("" + df.format(mleparams[3]));
      mle5.setText("" + df.format(mleparams[5]));
      mle6.setText("" + df.format(mleparams[4]));
      mle7.setText("" + df.format(mleparams[6]));
      
      
      diag1.setText("" + df.format(diagnostics[0]));
      diag2.setText("" + df.format(diagnostics[1]));
      diag3.setText("" + df.format(diagnostics[2]));
      diag4.setText("" + df.format(diagnostics[3]));     
      
      lb1.setText("" + df.format(lbv[0]));
      lb2.setText("" + df.format(lbv[1]));
      lb3.setText("" + df.format(lbv[2]));
      lb4.setText("" + df.format(lbv[3]));
 
      aic1.setText("" + df.format(lkhs[0]));
      aic2.setText("" + df.format(lkhs[1]));    
      aic3.setText("" + df.format(lkhs[2]));
      aic4.setText("" + df.format(lkhs[3]));      
       
    }
    
    
    public void changeHighlight(int c)
    {hlight_indicator = c;}

    public void changeModel(int _model){model = _model;}

    public void setSimulationParameters(int[] dim, double[] params)
    { 
        this.m_dim = dim; sarima.changeModelDimension(dim, false, outlier, easter, td, eaDays);
        this.sim_params = params;
    }

    public void setCompute(boolean c) {compute = c;}
    
    //=====  are we simulating series or uploading =======  
    public void setSimSeries(boolean sel) {sim_series = sel;}

    public void setSpecPlots(boolean sel) {specPlots = sel; }

    public void setRealData(double[] _real_series, int _nObs)
    {real_series = _real_series; nObs = _nObs; setMonthlyDates(year,month); minObsBar.setMaximum(nObs-12);}

 
    public void setMonthlyDates(int year, int month)
    {
       
       int i;
       dates = new String[nObs];

       int m,y;

       y = year; m = month;

       for(i=0;i<nObs;i++)
       { 
         if(m > 12) {m=1; y++;} //resent month and increase year
         dates[i] = "" + y + "-" + m;    
         m++;    
       }

    }



    public void changeSlidingSpan(int back, int fore)
    {
       
       if((nObs - (back + fore)) > minObs)  
       {
       
         nbackObs = back;
         nforeObs = fore;
         nObsSpan = nObs - nbackObs - nforeObs;
       
         initializeUseless(nObsSpan);
       
         computeSlidingSpan();
         setTexts();
       }
    }

  
    public void computeSlidingSpan()
    {

      if(slidingSpan)
      {
         
         nObsSpan = nObs - nbackObs - nforeObs;

         //System.out.println("nObsSpan = " + nObsSpan + ", nObs = " + nObs + ", nbackObs = " + nbackObs + ", nforeObs = " + nforeObs);  


          //---- Define shit--------
         int i; int total;  
         nx13out = 4+n_params+4+4+72+900+3*nObsSpan+n_poly;
         int NN = 4+n_params+4+4+72+900+3*nObsSpan;
        
         //----- get tseries span 
         span_data = new double[nObsSpan];

         for(i=0;i<nObsSpan;i++)
         {
           span_data[i] = t_series[nbackObs + i];        
         }

         mle_params = new double[n_params]; x13output = new double[nx13out+n_poly];
         seassig = new double[nObsSpan]; trendsig = new double[nObsSpan];
         cyclesig = new double[nObsSpan]; sasig = new double[nObsSpan];

         sarima.setSeats(seats);
       
         if(compute)
         {x13output = sarima.GetDiags(span_data, this.m_dim, this.sim_params, n_params, this.lag, nObsSpan, model);}
 
         for(i=0;i<4;i++)  
         {diagnostics[i] = x13output[i];}
         for(i=0;i<n_params;i++)
         {mle_params[i] = x13output[4+i];}
         for(i=0;i<4;i++)
         {lbv[i] = x13output[4+n_params+i];}
         for(i=0;i<4;i++)
         {lkhs[i] = x13output[4+n_params+4+i];}
         if(seats){  //don't change forecasts if no seats!
         for(i=0;i<72;i++)
         {forecasts[i] = x13output[4+n_params+4+4+i];}
         }
         for(i=0;i<900;i++)
         {sampleQIf[i] = x13output[4+n_params+4+4+72+i];} 
         if(seats){ //don't change signal extractions
          for(i=0;i<nObsSpan;i++)
          {
           seassig[i] = x13output[i+n_params+4+4+72+4+900]; 
           trendsig[i] = x13output[i+n_params+4+4+72+4+900 + nObsSpan];
           cyclesig[i] = x13output[i+n_params+4+4+72+4+900 + 2*nObsSpan];
          }
         }
         if(seats){ 
          for(i=0;i<8;i++) {trendpoly[i] = x13output[NN+i];}
          for(i=0;i<27;i++) {seaspoly[i] = x13output[NN+8+i];}
          for(i=0;i<32;i++) {tipoly[i] = x13output[NN+8+27+i];}
          for(i=0;i<5;i++) {trnsdenom[i] = x13output[NN+8+27+32+i];}
          for(i=0;i<36;i++) {trnspoly[i] = x13output[NN+8+27+32+5+i];}
         }
    
         m_p = m_dim[0];  m_q = m_dim[2]; m_P = m_dim[3];  m_Q = m_dim[5];
         total = m_p + m_q + m_P + m_Q + 1;
         if(total != n_params) {System.out.println("Somethings wrong");}

         m_innvar = mle_params[m_p + m_P + m_q + m_Q];   

         mphi = new Polynomial(m_p);
         mphi.coef[0] = 1.0;
         for(i=1;i<=m_p;i++)
         {mphi.coef[i] = mle_params[i-1];}

         mtheta = new Polynomial(m_q);
         mtheta.coef[0] = 1.0;
         for(i=1;i<=m_q;i++)
         {mtheta.coef[i] = -1.0*mle_params[m_p + i-1];} //Box notation
      
         if(m_P != 0)  
         {
           mPhi = new Polynomial(S);
           mPhi.coef[0] = 1.0;      
           mPhi.coef[S] = mle_params[m_p + m_q];
         }
         else
         {mPhi = new Polynomial(0); mPhi.coef[0] = 1.0;}  
     
         if(m_Q != 0)  
         {
          mTheta = new Polynomial(S);
          mTheta.coef[0] = 1.0;      
          mTheta.coef[S] = -1.0*mle_params[m_p + m_P + m_q];
         }
         else
         {mTheta = new Polynomial(0); mTheta.coef[0] = 1.0;}       
         arpoly = mphi.times(mPhi);
         mapoly = mtheta.times(mTheta);     
     
         if(specPlots && compute)
         {sarima.computeSampleIFw(nObsSpan, span_data, model);}

      }

    }



    //------ Computes a sample SARIMA model with the given parameters
    public void computeSARIMAmodel(boolean com)
    {
        
        //System.out.println("computes new model = " + com);
      	n_params = this.sim_params.length;
	sarima.SetSARIMAparams(this.sim_params, this.m_dim, n_params, S);  		              

	compute = com;
	
        if(sim_series)
        {t_series = sarima.sampleSARIMAmodel(this.seed);} 

        else{t_series = real_series; nObs = real_series.length;}

        if(com) 
        {seats = true;  dataSetup = true;}
 
        if(!slidingSpan)
        {
         this.computeDiagnostics(model);
         if(specPlots && com)
         {sarima.computeSampleIFw(nObs, t_series, model);}
        }
        else
        {
          computeSlidingSpan();
        }
    }

    public void computeSampleSpec()
    {sarima.computeSampleIFw(nObs, t_series, model);}
    

//     public void computeMinimumKL(int nobs, double[] data, int[] dim, double[] params, int[] est_Dim)
//     {
// 
//      int i,j;
//      int n_params = dim[0] + dim[2] + dim[3] + dim[5] + 1;  
//   
//     
//      sarima.SetSARIMAparams(params, dim, n_params, S);    
//  
//      int estNparams = est_dim[0] + est_dim[2] + est_dim[3] + est_dim[5] + 2;
//      double[] optParams = new double[estNparams];
// 
//      optParams = sarima.computeMinimumModel(dim, params, est_dim);  
//     
//      avocatos =4, chips=3, sausage=4, vodka sauce=3, wine=4, salmon patties=5, burrito=3, 
//      chicken wrap=4, mango=3, berries=3, strawberries=2, oatmeal=5
// 

//     }

    public void setPlotTracker(int t) 
    {
       if(t >= 0)
       {
         setMonthlyDates(year,month);
         plot_tracker = true;
         tracker = t;
         go();
       }
       else 
       {plot_tracker = false;}
    
    }


    public void setDateStart(int y, int m) 
    {year = y; month = m;}

    public void computeDiagnostics(int _model)
    {
       int i; model = _model; int total;  
       nx13out = 4+n_params+4+4+72+900+3*nObs+n_poly;
       int NN = 4+n_params+4+4+72+900+3*nObs;
       
       
       
       mle_params = new double[n_params];
       x13output = new double[nx13out+n_poly];
       seassig = new double[nObs];
       trendsig = new double[nObs];
       cyclesig = new double[nObs];
       sasig = new double[nObs];
       lkhs = new double[4];

       sarima.setSeats(seats);
       
       if(compute)
       {
        x13output = sarima.GetDiags(this.t_series, this.m_dim, this.sim_params, n_params, this.lag, nObs, model);
       }
       for(i=0;i<4;i++)  
       {diagnostics[i] = x13output[i];}
       for(i=0;i<n_params;i++)
       {mle_params[i] = x13output[4+i];}
       for(i=0;i<4;i++)
       {lbv[i] = x13output[4+n_params+i];}
       for(i=0;i<4;i++)
       {lkhs[i] = x13output[4+n_params+4+i];}
      if(seats){  //don't change forecasts if no seats!
       for(i=0;i<72;i++)
       {forecasts[i] = x13output[4+n_params+4+4+i];}
       }
       for(i=0;i<900;i++)
       {sampleQIf[i] = x13output[4+n_params+4+4+72+i];} 
      if(seats){ //don't change signal extractions
       for(i=0;i<nObs;i++)
       {
        seassig[i] = x13output[i+n_params+4+4+72+4+900]; 
        trendsig[i] = x13output[i+n_params+4+4+72+4+900 + nObs];
        cyclesig[i] = x13output[i+n_params+4+4+72+4+900 + 2*nObs];
       }
      }
     if(seats){ 
      for(i=0;i<8;i++) {trendpoly[i] = x13output[NN+i];}
      for(i=0;i<27;i++) {seaspoly[i] = x13output[NN+8+i];}
      for(i=0;i<32;i++) {tipoly[i] = x13output[NN+8+27+i];}
      for(i=0;i<5;i++) {trnsdenom[i] = x13output[NN+8+27+32+i];}
      for(i=0;i<36;i++) {trnspoly[i] = x13output[NN+8+27+32+5+i];}
     }
    
     m_p = m_dim[0];  m_q = m_dim[2]; m_P = m_dim[3];  m_Q = m_dim[5];
     total = m_p + m_q + m_P + m_Q + 1;
     if(total != n_params)
     {System.out.println("Somethings wrong");}
  
     m_innvar = mle_params[m_p + m_P + m_q + m_Q];   

     mphi = new Polynomial(m_p);
     mphi.coef[0] = 1.0;
     for(i=1;i<=m_p;i++)
     {mphi.coef[i] = mle_params[i-1];}

     mtheta = new Polynomial(m_q);
     mtheta.coef[0] = 1.0;
     for(i=1;i<=m_q;i++)
     {mtheta.coef[i] = -1.0*mle_params[m_p + i-1];} //Box notation
      
     if(m_P != 0)  
     {
      mPhi = new Polynomial(S);
      mPhi.coef[0] = 1.0;      
      mPhi.coef[S] = mle_params[m_p + m_q];
     }
     else
     {mPhi = new Polynomial(0); mPhi.coef[0] = 1.0;}  
     
     if(m_Q != 0)  
     {
      mTheta = new Polynomial(S);
      mTheta.coef[0] = 1.0;      
      mTheta.coef[S] = -1.0*mle_params[m_p + m_P + m_q];
     }
     else
     {mTheta = new Polynomial(0); mTheta.coef[0] = 1.0;}  
     
     arpoly = mphi.times(mPhi);
     mapoly = mtheta.times(mTheta);     
     
    }

    public void setSeries(boolean se) {plotSeries = se;}
    public void setForecasts(boolean fr) {plotForecasts = fr;}
    public void setSeas(boolean s) {plotSeas = s;}
    public void setTrend(boolean t) {plotTrend = t;}
    public void setCycle(boolean c) {plotCycle = c;}
    public void setSA(boolean c) {plotSA = c;}    

    public void modelDimensionChange(int[] dim, boolean _trans, boolean _outlier, boolean _easter, boolean _td, int _ea)
    {
      trans = _trans; outlier = _outlier; td = _td; eaDays = _ea; easter = _easter; 
      sarima.changeModelDimension(dim, trans, outlier, easter, td, eaDays);
    }

    public void setLag(int _lag)
    {this.lag = _lag;} 

    
    public void computeDataMax()
    {
         int i;
         dataMax = -1000000;
         dataMin = 10000000;
  
         for(i=0; i < nObs; i++)
         {
            if(t_series[i] > dataMax) dataMax = t_series[i];
            else if(t_series[i] < dataMin) dataMin = t_series[i];   
         }
         dataNorm = Math.abs(dataMax - dataMin);

         setMonthlyDates(year,month);
    }
 
    public void setInnovation(double _innvar)
    {
      int len = this.sim_params.length; 
      this.m_innvar = _innvar;  
      this.sim_params[len-1] = _innvar;
    }  

    public void setNobs(int _nObs)
    {
      this.nObs = _nObs;
      sarima.SetObservations(_nObs);
      t_series = new double[_nObs];
      minObsBar.setMaximum(nObs-12);
    }

    public void setParams(int dims[], double params[], int nparams, double _innvar)
    {
       int j=0;
       this.sim_params = new double[nparams];          //new vector of doubles
        //set innovation var  
       m_p = dims[0]; m_d = dims[1]; m_q = dims[2]; 
       m_P = dims[3]; m_D = dims[4]; m_Q = dims[5]; 
      
      //for(i=0; i < 7; i++) System.out.print(params[i] + " ");
      //System.out.println(" ");
  
      if(m_p >= 1) {sim_params[j] = params[0]; j++;}      
      if(m_p == 2) {sim_params[j] = params[1]; j++;}      
      if(m_q >= 1) {sim_params[j] = params[2]; j++;}      
      if(m_q == 2) {sim_params[j] = params[3]; j++;}      
      if(m_P == 1) {sim_params[j] = params[4]; j++;}            
      if(m_Q == 1) {sim_params[j] = params[5]; j++;}       
       
      this.sim_params[nparams-1] = _innvar;
      this.m_dim = dims;
      this.n_params = nparams;
       //System.out.println(n_params);
      this.m_innvar = _innvar;
      sarima.SetSARIMAparams(this.sim_params, this.m_dim, n_params, S);  
    }
    
    public void setBurnin(int _burnin)
    {      
      this.burnin = _burnin;
      sarima.SetBurnin(burnin);
    }

    public void setSeed(int _seed)
    {this.seed = _seed;}

    
    public void initializeUseless(int n)
    {
       int i;
       try{  
           PrintWriter out = new PrintWriter(new FileWriter("useless.dat"));
           for(i=0; i < n; i++) {out.println(0.0);}
 
           out.close();
        } catch (IOException e) {e.printStackTrace();}
    }    
    
    
    public void go()
    {repaint();}

    public void paintComponent(Graphics g)
    {
        int i,w,h;
        int t0,t1,x0,x1;
		
        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;

        BasicStroke bold = new BasicStroke((float)2.0);

	    super.paintComponent(g);
	    g2d = (Graphics2D)g;

        BasicStroke def = new BasicStroke((float)1.9);
        //Stroke def = g2d.getStroke();




       //--------------------- Draw grid with or without forecast points ------------
        float[] dash1 = {10.0f};
        BasicStroke dashed = new BasicStroke(1.0f, 
                                          BasicStroke.CAP_BUTT, 
                                          BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);



       if(!plotForecasts || slidingSpan) 
       {nObsex = nObs;}
       else 
       {nObsex = nObs+24;}

         
        if(slidingSpan) //---- draw dark blue background
        {
          g2d.setStroke(def);
          t0 = (int)(((double)(nbackObs)/(double)nObsex)*(double)width);
          t1 = (int)(((double)(nbackObs+nObsSpan-1)/(double)nObsex)*(double)width);
          g2d.setPaint(new Color(7,9,34));

          rectangle = new Rectangle(t0, 0, t1-t0, height);
          g2d.draw(rectangle);  
          g2d.fill(rectangle);
           
          g2d.setPaint(new Color(14,38,128));
          g2d.setStroke(def);
          g2d.drawLine(t0, 0, t0, height);
          g2d.drawLine(t1, 0, t1, height);

          slideStart = t0;
          slideEnd = t1;
          
          //System.out.println("line t0 = " + t0 + ", line t1 = " + t1);
        }


        g2d.setStroke(dashed);
        g2d.setPaint(myGray);

       if(!plotForecasts) 
       {   

        int nobsP = (int)Math.floor((double)nObs/12); int p = 12; 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*12)/nObs)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-10);
          g.drawString((String)"" + p, t0, height - 5);
          p = p + 12;
        } 
       }
       else
       {   
        int nobsP = (int)Math.floor((double)nObsex/12); int p = 12; 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*12)/nObsex)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-10);
          g.drawString((String)"" + p, t0, height - 5);
          p = p + 12;
        } 

       }        
        
        for(i=0; i < 9; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }

        g.drawString((String)df.format(dataMax), 5, 15);
        g.drawString((String)df.format(dataMin), 5, height - 5);


        
        

       //----------------------------------------------------------------------------       
       g2d.setStroke(def);
       g2d.setPaint(new Color(110,113,130));
       for(i = 0; i < nObs-1; i++)
       {
	    t0 = (int)(((double)i/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObsex)*(double)width);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((t_series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
       }
	
       if(compute && !slidingSpan)
       {
       

	  //Draw Original Data
        if(plotSeries)
        {
	  g2d.setPaint(seriesColor); g2d.setStroke(def);
	  if(!plotForecasts) 
          { 
           if(hlight_indicator == 0){g2d.setPaint(highlight); g2d.setStroke(bold);}
           for(i = 0; i < nObs-1; i++)
	   {
	    t0 = (int)(((double)i/(double)this.nObs)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)this.nObs)*(double)width);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((t_series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	   }
          }
          else if(plotForecasts)
          {
           g2d.setStroke(def);
           if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}
           for(i = 0; i < nObs-1; i++)
	   {            	
	    t0 = (int)(((double)i/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObsex)*(double)width);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((t_series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	   }	

            i=nObs-1; g2d.setPaint(seriesColor2); g2d.setStroke(def); if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    t0 = (int)(((double)i/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObsex)*(double)width);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[24] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            g2d.setPaint(forecol);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[0] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            g2d.setPaint(forecol);
	    x0 = (int)(((t_series[i] - dataMin)/(dataNorm*1.20))*(double)height);
 	    x1 = (int)(((forecasts[48] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

	   
           //g2d.setPaint(Color.RED); //upper envelope
           for(i = 0; i < 23; i++)
	   {            	
            g2d.setPaint(forecol); g2d.setStroke(def);
	    t0 = (int)(((double)(i+nObs)/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+nObs+1)/(double)nObsex)*(double)width);
	    x0 = (int)(((forecasts[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	    g2d.setPaint(seriesColor2); g2d.setStroke(def); if(hlight_indicator == 1){g2d.setPaint(highlight);}
            x0 = (int)(((forecasts[i+24] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[i+25] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1); 
            g2d.setPaint(forecol); g2d.setStroke(def);
	    x0 = (int)(((forecasts[i+48] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[i+49] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
           }
          }	
        }

        for(i = 0; i < nObs-1; i++)
	{
	    t0 = (int)(((double)i/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObsex)*(double)width);

            if(plotCycle){g2d.setPaint(cycleColor); g2d.setStroke(def);  if(hlight_indicator == 5){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    x0 = (int)(((cyclesig[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((cyclesig[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}   
	    
            if(plotSeas){g2d.setPaint(seasonalColor); g2d.setStroke(def); if(hlight_indicator == 2){g2d.setPaint(highlight);  g2d.setStroke(bold);}
            x0 = (int)((((seassig[i] + t_series[0] - dataMin))/(dataNorm*1.20))*(double)height);
	    x1 = (int)((((seassig[i+1] + t_series[0] - dataMin))/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}
   
            if(plotTrend){g2d.setPaint(trendColor); g2d.setStroke(def); if(hlight_indicator == 3){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    x0 = (int)(((trendsig[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((trendsig[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}

            if(plotSA){g2d.setPaint(saColor); g2d.setStroke(def); if(hlight_indicator == 4){g2d.setPaint(highlight); g2d.setStroke(bold);}
	    x0 = (int)((((t_series[i] - seassig[i]) - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)((((t_series[i+1] - seassig[i+1]) - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}     
    
        }



        
       if(plot_tracker)
       {
        w = 6; h = 6; 
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t-3, (height-23) - track_pos_x, w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);
        
        g.setFont(mono);
        g2d.setPaint(Color.GREEN);
        g.drawString(value, 5, 15);        
        //setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));    
       }        
      }
      else if(compute && slidingSpan) //---plot inside the span only
      {


        if(plotSeries)
        {
	  g2d.setPaint(seriesColor); g2d.setStroke(def);
	  if(!plotForecasts) 
          { 
           if(hlight_indicator == 0){g2d.setPaint(highlight); g2d.setStroke(bold);}
           for(i = 0; i < nObsSpan-1; i++)
	   {
	    t0 = (int)(((double)(nbackObs+i)/(double)this.nObs)*(double)width);
	    t1 = (int)(((double)(nbackObs+i+1)/(double)this.nObs)*(double)width);
	    x0 = (int)(((span_data[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((span_data[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	   }
          }
          else if(plotForecasts)
          {
            g2d.setStroke(def);
            if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}
            for(i = 0; i < nObsSpan-1; i++)
	    {            	
	     t0 = (int)(((double)(nbackObs+i)/(double)nObsex)*(double)width);
	     t1 = (int)(((double)(nbackObs+i+1)/(double)nObsex)*(double)width);
	     x0 = (int)(((span_data[i] - dataMin)/(dataNorm*1.20))*(double)height);
	     x1 = (int)(((span_data[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	     g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	    }	

            i=nObsSpan-1; g2d.setPaint(seriesColor2); g2d.setStroke(def); if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    t0 = (int)(((double)(nbackObs+i)/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(nbackObs+i+1)/(double)nObsex)*(double)width);
	    
	    x0 = (int)(((span_data[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[24] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            /*g2d.setPaint(forecol);
	    x0 = (int)(((span_data[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[0] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            g2d.setPaint(forecol);
	    x0 = (int)(((span_data[i] - dataMin)/(dataNorm*1.20))*(double)height);
 	    x1 = (int)(((forecasts[48] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            */
	   
           //g2d.setPaint(Color.RED); //upper envelope
           for(i = 0; i < 23; i++)
	   {            	
            g2d.setPaint(forecol); g2d.setStroke(def);
	    t0 = (int)(((double)(nbackObs+i+nObsSpan)/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(nbackObs+i+nObsSpan+1)/(double)nObsex)*(double)width);
	    
	    //x0 = (int)(((forecasts[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    //x1 = (int)(((forecasts[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    //g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	    g2d.setPaint(seriesColor2); g2d.setStroke(def); if(hlight_indicator == 1){g2d.setPaint(highlight);}
            x0 = (int)(((forecasts[i+24] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((forecasts[i+25] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1); 
            //g2d.setPaint(forecol); g2d.setStroke(def);
	    //x0 = (int)(((forecasts[i+48] - dataMin)/(dataNorm*1.20))*(double)height);
	    //x1 = (int)(((forecasts[i+49] - dataMin)/(dataNorm*1.20))*(double)height);
	    //g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
           }
          }	
        }

        for(i = 0; i < nObsSpan-1; i++)
	{
	    t0 = (int)(((double)(nbackObs+i)/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(nbackObs+i+1)/(double)nObsex)*(double)width);

            if(plotCycle){g2d.setPaint(cycleColor); g2d.setStroke(def);  if(hlight_indicator == 5){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    x0 = (int)(((cyclesig[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((cyclesig[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}   
	    
            if(plotSeas){g2d.setPaint(seasonalColor); g2d.setStroke(def); if(hlight_indicator == 2){g2d.setPaint(highlight);  g2d.setStroke(bold);}
            x0 = (int)((((seassig[i] + span_data[0] - dataMin))/(dataNorm*1.20))*(double)height);
	    x1 = (int)((((seassig[i+1] + span_data[0] - dataMin))/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}
   
            if(plotTrend){g2d.setPaint(trendColor); g2d.setStroke(def); if(hlight_indicator == 3){g2d.setPaint(highlight);  g2d.setStroke(bold);}
	    x0 = (int)(((trendsig[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((trendsig[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}

            if(plotSA){g2d.setPaint(saColor); g2d.setStroke(def); if(hlight_indicator == 4){g2d.setPaint(highlight); g2d.setStroke(bold);}
	    x0 = (int)((((span_data[i] - seassig[i]) - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)((((span_data[i+1] - seassig[i+1]) - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);}     
    
        }


        g2d.setStroke(dashed);
        g2d.setPaint(myGray);
              
        for(i=0; i < 9; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }

        g.drawString((String)df.format(dataMax), 5, 15);
        g.drawString((String)df.format(dataMin), 5, height - 5);
        
       if(plot_tracker)
       {
        w = 6; h = 6; 
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t-3, (height-23) - track_pos_x, w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);
        
        g.setFont(mono);
        g2d.setPaint(Color.GREEN);
        g.drawString(value, 5, 15);        
        //setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));    
       }        




      }   
        
        
    }


    /*--------------------------------------------------------

       Idea here is to 'sweep' data with given model and compute
       the mean, std, and forecast error up to K points ahead
       of the MLE values and the diagnostics

       //34743
       
    ---------------------------------------------------------*/
    public void sweepWithModel(int _nfore)
    {
      minObs = minObsBar.getValue();
      if(nObs > minObs)
      {

       error = 0.0;
       records = new double[7][nObs+1-minObs];
       nfore = _nfore;
       means = new double[12];
       sdevs = new double[12];

       PrimeThread thread = new PrimeThread();
       
       //(new Thread(this)).start();
       thread.start();
       
    }
   }


    
   class PrimeThread extends Thread 
   {
   public void run()
   {
       int i,j;
       
       error = 0;
         
         for(i=minObs; i <= nObs; i++)
         {
 
          changeSlidingSpan(0, nObs-i);
        
          for(j=0;j<7;j++)
          {
           means[j] = means[j] + mleparams[j];
           records[j][i-minObs] = mleparams[j];
           
          }

          means[7] = means[7] + lbv[0];
          for(j=0;j<3;j++) 
          {means[8+j] = means[8+j] + diagnostics[j];}
                   
          means[11] = means[11] + lbv[2];
          //System.out.println(means[0] + " " + means[1] + " " + means[2] + " " + means[3] + " " + means[4] + " " + means[5]);
          
         
          if(i + nfore < nObs)
          {
           for(j=0;j<nfore;j++)
           {
            error = error + (t_series[i+j] - forecasts[j+24])*(t_series[i+j] - forecasts[j+24]); 
           }
           error = error + Math.sqrt(error)/nfore;
          }
          
          
          repaint();

          try{Thread.sleep(100);}
          catch(Exception e){} 
   
   
       }
       
       //-------- Now compute statistics ---------------------------------------------
       
       for(j=0;j<12;j++)
       {means[j] = means[j]/(nObs-59.0);}  

       //System.out.println(means[7] + " " + means[8]);
       
       
       for(j=0;j<7;j++)
       {           
        for(i=minObs; i <= nObs; i++)
        {
          sdevs[j] = sdevs[j] + (records[j][i-minObs] - means[j])*(records[j][i-minObs] - means[j]);              
        }          
       } 

       for(j=0;j<7;j++)
       {sdevs[j] = Math.sqrt(sdevs[j])/(nObs-minObs-1.0);}// System.out.println("std = " + sdevs[j]);}       
       
       //-------- Done with thread---------------------------------------
       updateSweepStats();
       return;
   }
   
  } 
    
    
     class MyArtMouseMotionListener extends MouseMotionAdapter 
     {
      
      public void mouseDragged(MouseEvent e) 
      {
       
        if(slidingSpan)
        {
         boolean passx = false;
         
         int t1 = (int)(((double)nObsex)*e.getX()/(double)width); 
         int t0 = (int)(((double)nObsex)*pressed_t/(double)width);  
         
         
         //if((e.getX() > slideStart+2) && (e.getX() < slideEnd-2))
         if(curCursor == Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR))
         {                   
           if(t1 - t0 > 0) //going to the right
           {
            if(nforeObs-1 >= 0)
            {
              nforeObs = nforeObs - 1;
              nbackObs = nbackObs + 1;        
              pressed_t = e.getX(); 
              passx = true;
            }
           }
           else if(t1 - t0 < 0)
           {
             if(nbackObs-1 >=0)
             {
               nforeObs = nforeObs + 1;
               nbackObs = nbackObs - 1;
               pressed_t = e.getX();
               passx = true;
             }
           }
           changeSlidingSpan(nbackObs, nforeObs);
           go();
           
         }
         else if(curCursor == Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR))
         {         
           if(t1 - t0 > 0) //going to the right
           {
            if(nObsSpan >= minObs)
            {nbackObs = nbackObs + 1; pressed_t = e.getX(); passx = true;}
           }
           else if(t1 - t0 < 0)
           {
             if(nbackObs-1 >=0)
             {nbackObs = nbackObs - 1; pressed_t = e.getX(); passx = true;}
           }
           changeSlidingSpan(nbackObs, nforeObs);
           go();
           
         }
         else if(curCursor == Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR))
         {
           if(t1 - t0 > 0) //going to the right
           {
            if(nforeObs-1 >= 0)
            {nforeObs = nforeObs - 1; pressed_t = e.getX(); passx = true;}
           }
           else if(t1 - t0 < 0)
           {
             if(nObsSpan >= minObs)
             {nforeObs = nforeObs + 1; pressed_t = e.getX(); passx = true;}
           }
           
           if(passx)
           {changeSlidingSpan(nbackObs, nforeObs);}
           
           go();
           
         }         
        }
       }
      

      public void mouseMoved(MouseEvent e) 
      {

        int j,t0,t1;
        if(plot_tracker)
        {
          if(!slidingSpan)
          {
           curCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
           setCursor(curCursor);  
             
           j = (int)(((double)nObsex)*e.getX()/(double)width);  
           if(j >= nObs) {j=nObs-1;}
           
           track_pos_t = (int)(((double)j/(double)nObsex)*(double)width);       
           //track_pos_t = (int)(((double)j/(double)trade_obs)*(double)width);
           //t1 = (int)(((double)(j+1)/(double)trade_obs)*(double)width);

           if(tracker == 0 && plotSeries) //account
           { 
             track_pos_x = (int)(((t_series[j] - dataMin)/(dataNorm*1.20))*(double)height);                  value = dates[j] + ", " + df.format(t_series[j]);      
           }
           else if(tracker == 1 && plotCycle) //price
           {
             track_pos_x =  (int)(((cyclesig[j] - dataMin)/(dataNorm*1.20))*(double)height);                 value = dates[j] + ", " + df.format(cyclesig[j]);          
           }       
           else if(tracker == 2 && plotTrend) //signal
           {
            track_pos_x = (int)(((trendsig[j] - dataMin)/(dataNorm*1.20))*(double)height);                   value = dates[j] + ", " + df.format(trendsig[j]);
           }
           else if(tracker == 3 && plotSeas) //diff
           {
             track_pos_x =  (int)((((seassig[j] + t_series[0] - dataMin))/(dataNorm*1.20))*(double)height);  value = dates[j] + ", " + df.format(seassig[j]);          
           }
           else if(tracker == 4 && plotSA) //diff
           {
             track_pos_x = (int)((((t_series[j] - seassig[j]) - dataMin)/(dataNorm*1.20))*(double)height);   value = dates[j] + ", " + df.format(t_series[j] - seassig[j]);          
           }           
           //System.out.println(j + " " + track_pos_t + "  " + track_pos_x);
           go();
          }
          else if(slidingSpan)
          {

           j = (int)(((double)nObsex)*e.getX()/(double)width);  
           if(j >= nbackObs+nObsSpan) {j=nObsSpan-1;}
           else if(j < nbackObs) {j = 0;}

           if(j > nObsSpan-1) {j=nObsSpan-1;}

           track_pos_t = (int)(((double)(nbackObs+j)/(double)nObsex)*(double)width);       

           if(tracker == 0 && plotSeries) //account
           { 
             track_pos_x = (int)(((span_data[j] - dataMin)/(dataNorm*1.20))*(double)height); value = dates[j] + ", " + df.format(span_data[j]);      
           }
           else if(tracker == 1 && plotCycle) //price
           {
             track_pos_x =  (int)(((cyclesig[j] - dataMin)/(dataNorm*1.20))*(double)height); value = dates[j] + ", " + df.format(cyclesig[j]);          
           }       
           else if(tracker == 2 && plotTrend) //signal
           {
            track_pos_x = (int)(((trendsig[j] - dataMin)/(dataNorm*1.20))*(double)height); value = dates[j] + ", " + df.format(trendsig[j]);
           }
           else if(tracker == 3 && plotSeas) //diff
           {
             track_pos_x =  (int)((((seassig[j] + span_data[0] - dataMin))/(dataNorm*1.20))*(double)height);  value = dates[j] + ", " + df.format(seassig[j]);          
           }
           else if(tracker == 4 && plotSA) //diff
           {
             track_pos_x = (int)((((span_data[j] - seassig[j]) - dataMin)/(dataNorm*1.20))*(double)height);  value = dates[j] + ", " + df.format(span_data[j] - seassig[j]);          
           }           
           //System.out.println(j + " " + track_pos_t + "  " + track_pos_x);          
          }
 
        }

        t0 = (int)(((double)(nbackObs)/(double)nObsex)*(double)width);
        t1 = (int)(((double)(nbackObs+nObsSpan-1)/(double)nObsex)*(double)width);
        if(e.getX() > t0+2 && e.getX() < t1-2) 
        { 
           curCursor = Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR);
           setCursor(curCursor);                        
        }
        else if(e.getX() >= t0-2 && e.getX() <= t0+2)
        {
           curCursor = Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR);
           setCursor(curCursor);
        }
        else if(e.getX() >= t1-2 && e.getX() <= t1+2)
        {
           curCursor = Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR);
           setCursor(curCursor);
        }
        else
        {
           curCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
           setCursor(curCursor);
        } 
        go();
           
      }
      
     }       

     class MyArtMouseListener extends MouseAdapter 
     {
       
       public void mousePressed(MouseEvent e) 
       { pressed_t = e.getX(); }// System.out.println("pressed = " + pressed_t);}

       public void mouseReleased(MouseEvent e) 
       { }

       public void mouseClicked(MouseEvent e) 
       { }
        
     }
     
     
     
     
     
     
   public void activateSweepInterface()
   {

        sweepPanel = new JPanel();
        sphi1Label = new JLabel();
        meanLabel = new JLabel();
        meanphi1Text = new JTextField();
        meanphi2Label = new JLabel();
        meanphi2Text = new JTextField();
        meanPhiText = new JTextField();
        sPhiLabel = new JLabel();
        meanLabel1 = new JLabel();
        sphi1Label1 = new JLabel();
        meanphi2Label1 = new JLabel();
        meanphi1Text1 = new JTextField();
        meanphi2Text1 = new JTextField();
        meanPhiText1 = new JTextField();
        sPhiLabel1 = new JLabel();
        aiccLabel = new JLabel();
        aiccText = new JTextField();
        bicLabel = new JLabel();
        bicText = new JTextField();
        foreErrorLabel = new JLabel();
        foreErrorText = new JTextField();
        foreBar = new JScrollBar(JScrollBar.HORIZONTAL,1,1,1,24);
        foreLabel = new JLabel();
        foreHorizText = new JTextField();
        compSweepButton = new JButton();
        meanAicLabel = new JLabel();
        minObsLabel = new JLabel();

        minObsBar = new JScrollBar(JScrollBar.HORIZONTAL,60,1,40,nObs-12);
        minObsText = new JTextField();
        minObsText.setText("60"); 
        minObsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        minObsLabel.setText("Starting Observation");

        sphi1Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sphi1Label.setText("phi_1");

        meanLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        meanLabel.setText("Mean, StD");

        meanphi1Text.setText("0");

        meanphi2Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        meanphi2Label.setText("phi_2");

        meanphi2Text.setText("0");

        meanPhiText.setText("0");

        sPhiLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sPhiLabel.setText("Phi");

        meanLabel1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        meanLabel1.setText("Mean, StD");

        sphi1Label1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sphi1Label1.setText("theta_1");

        meanphi2Label1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        meanphi2Label1.setText("theta_2");

        meanphi1Text1.setText("0");

        meanphi2Text1.setText("0");

        meanPhiText1.setText("0");

        sPhiLabel1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sPhiLabel1.setText("Theta");

        aiccLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        aiccLabel.setText("LB-0");

        aiccText.setText("0");

        bicLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        bicLabel.setText("LB-12");

        bicText.setText("0");

        foreErrorLabel.setText("Forecast Error");

        foreErrorText.setText("0");

        foreBar.setOrientation(JScrollBar.HORIZONTAL);

        foreLabel.setText("Forecast Horizon");

        foreHorizText.setText("1");

        compSweepButton.setText("Compute Time Series Sweep");

        meanAicLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        meanAicLabel.setText("Mean");



        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(sweepPanel);
        sweepPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                    .addComponent(meanphi2Label)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(meanphi2Text, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(sphi1Label)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(meanLabel)
                                        .addComponent(meanphi1Text, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE))))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(sPhiLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(meanPhiText, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGap(33, 33, 33)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(meanphi2Label1)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(meanphi2Text1, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(sPhiLabel1)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(meanPhiText1, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(sphi1Label1)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(meanLabel1)
                                    .addComponent(meanphi1Text1, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE))))
                        .addGap(33, 33, 33)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(bicLabel)
                            .addComponent(foreErrorLabel)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(aiccLabel)
                                .addGap(5, 5, 5)))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(aiccText)
                            .addComponent(bicText)
                            .addComponent(foreErrorText, javax.swing.GroupLayout.PREFERRED_SIZE, 70, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(meanAicLabel))
                        .addContainerGap())
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(foreLabel)
                                        .addGap(18, 18, 18)
                                        .addComponent(foreBar, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(foreHorizText, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(minObsLabel)
                                        .addGap(18, 18, 18)
                                        .addComponent(minObsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(minObsText, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                .addGap(129, 129, 129))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(compSweepButton)
                                .addGap(165, 165, 165))))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(meanLabel1)
                            .addComponent(meanAicLabel))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(sphi1Label1)
                            .addComponent(meanphi1Text1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(aiccLabel)
                            .addComponent(aiccText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(meanphi2Label1)
                            .addComponent(meanphi2Text1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(bicLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(bicText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(sPhiLabel1)
                            .addComponent(meanPhiText1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(foreErrorLabel)
                            .addComponent(foreErrorText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(meanLabel)
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(sphi1Label)
                            .addComponent(meanphi1Text, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(meanphi2Label)
                            .addComponent(meanphi2Text, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(sPhiLabel)
                            .addComponent(meanPhiText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(foreBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(foreLabel, javax.swing.GroupLayout.Alignment.TRAILING))
                    .addComponent(foreHorizText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(minObsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(minObsBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(minObsLabel, javax.swing.GroupLayout.Alignment.TRAILING)))
                .addGap(18, 18, 18)
                .addComponent(compSweepButton)
                .addContainerGap())
        );










/*
        GroupLayout layout = new GroupLayout(sweepPanel);
        sweepPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                .addGroup(GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                    .addComponent(meanphi2Label)
                                    .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(meanphi2Text, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE))
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(sphi1Label)
                                    .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                    .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                        .addComponent(meanLabel)
                                        .addComponent(meanphi1Text, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE))))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(sPhiLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(meanPhiText, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE)))
                        .addGap(33, 33, 33)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(GroupLayout.Alignment.TRAILING, layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(meanphi2Label1)
                                    .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(meanphi2Text1, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE))
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(sPhiLabel1)
                                    .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(meanPhiText1, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE)))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(sphi1Label1)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                    .addComponent(meanLabel1)
                                    .addComponent(meanphi1Text1, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE))))
                        .addGap(33, 33, 33)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                            .addComponent(bicLabel)
                            .addComponent(foreErrorLabel)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(aiccLabel)
                                .addGap(5, 5, 5)))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                            .addComponent(aiccText)
                            .addComponent(bicText)
                            .addComponent(foreErrorText, GroupLayout.PREFERRED_SIZE, 70, GroupLayout.PREFERRED_SIZE)
                            .addComponent(meanAicLabel)))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(150, 150, 150)
                        .addComponent(foreLabel)
                        .addGap(18, 18, 18)
                        .addComponent(foreBar, GroupLayout.PREFERRED_SIZE, 88, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(foreHorizText, GroupLayout.PREFERRED_SIZE, 40, GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(177, 177, 177)
                        .addComponent(compSweepButton)))
                .addContainerGap(28, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(55, 55, 55)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(meanLabel1)
                            .addComponent(meanAicLabel))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(sphi1Label1)
                            .addComponent(meanphi1Text1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(aiccLabel)
                            .addComponent(aiccText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(meanphi2Label1)
                            .addComponent(meanphi2Text1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(bicLabel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(bicText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(sPhiLabel1)
                            .addComponent(meanPhiText1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(foreErrorLabel)
                            .addComponent(foreErrorText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(meanLabel)
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(sphi1Label)
                            .addComponent(meanphi1Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(meanphi2Label)
                            .addComponent(meanphi2Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(sPhiLabel)
                            .addComponent(meanPhiText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))))
                .addGap(37, 37, 37)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(foreHorizText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addGroup(GroupLayout.Alignment.TRAILING, layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(foreBar, GroupLayout.Alignment.TRAILING, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(foreLabel, GroupLayout.Alignment.TRAILING)))
                .addGap(18, 18, 18)
                .addComponent(compSweepButton)
                .addGap(29, 29, 29))
        );
     */
          
         AdjustmentListener ad = new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {
            foreHorizText.setText(""+foreBar.getValue());
          }};

         AdjustmentListener ad2 = new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {
            minObsText.setText(""+minObsBar.getValue());
          }};

 

         ActionListener alc = new ActionListener() {
         public void actionPerformed(ActionEvent event) 
         {
         sweepWithModel(foreBar.getValue()); 
         }};

         foreBar.addAdjustmentListener(ad);
         minObsBar.addAdjustmentListener(ad2);
         compSweepButton.addActionListener(alc);

    } 


  public void updateSweepStats()
  {
     
     meanphi1Text.setText(df.format(means[0]) + ", " + df2.format(sdevs[0])); 
     meanphi2Text.setText(df.format(means[1]) + ", " + df2.format(sdevs[1])); 
     meanphi1Text1.setText(df.format(means[2]) + ", " + df2.format(sdevs[2])); 
     meanphi2Text1.setText(df.format(means[3]) + ", " + df2.format(sdevs[3])); 
     meanPhiText.setText(df.format(means[4]) + ", " + df2.format(sdevs[4])); 
     meanPhiText1.setText(df.format(means[5]) + ", " + df2.format(sdevs[5])); 

     aiccText.setText(df.format(means[7]));
     bicText.setText(df.format(means[11]));
     foreErrorText.setText(df.format(error));
    
  }
     
     
     
     
     
     

}

