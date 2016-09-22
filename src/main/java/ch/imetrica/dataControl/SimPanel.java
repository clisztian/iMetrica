package ch.imetrica.dataControl;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.JCheckBox;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.text.*;
import rcaller.RCaller;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.Covariance;

import ch.imetrica.bayesCronos.HeavyModel;


/*----------------------------------------------------------------
 
   Panel for Simulation Controls - 
     Simulation:  
        Two different types of simulations 

   

----------------------------------------------------------------*/


public class SimPanel extends JPanel
{


   /**
	 * 
	 */
   private static final long serialVersionUID = 1L;

   private JFrame frame;

   private SimCanvas theatre; 
   int n_rep;        //---- total number of variables
   private int n_obs;        //---- total number of observations
   public int n_sym;
   int rand;         //---- random data type gaussian, student_t, etc.
   int burnin;       //---- burnin number
   int S;            //---- seasonal pattern 
   
   int MAX_DATA = 550000;
   int Lag;   
   double affect;
   double f1,f2;
   
   


   int max_series = 12;   //---- max_series

   ArrayList<double[]> sim_orig;  //---- original data, original scale---
   private ArrayList<double[]> sim_data;  //---- complete simulated data --------
   ArrayList<double[]> sim_means;   //store the means and standard dev

   int[] seed = {1,1,1,1,1,1};

   JCheckBox[] timePlot;
   JButton[] addMe;
   JButton integrate,integratesv,integrateHeavy;
   JButton addLagTrend;

   double[] tseries;              //----- time series sarima
   double[] cycle_1;              //----- time series cycle 1
   ArrayList<double[]> cycle_2;   //----- time series cycle 2
   double[] trend;                //----- time series trend 
   double[] garch;                //----- time series garch
   double[] real_series;          //----- time series real-uploaded
   double[] volume_series;        //----- time series for volume of asset
   public double[] hf_data; 
   
   //--------- Heavy models------------------------
   HeavyModel heavy;
   ArrayList<double[]> heavy_data;
   int n_heavy = 2;   
   double[] heavy_ts;


   double[] isv;                  //----- time series isv
   ArrayList<double[]> fsv;       //----- time series fsv;
   int sv_rep;
   double phisv, musv, alphasv;
   int acf_plot_no1 = 0; 
   int acf_plot_no2 = 0;  
   int acfs = 0; 

   boolean[] stationary;          //----- indicates whether stationary or not
   int n_cyc;
   boolean[] sim_check;
   JCheckBox[] sim_box;
   DecimalFormat df,df2,df3;
   //------- Arima stuff -------------------------------
   double[] sar_params; double[] params; 
   int[] dim; int n_params; double var; 

   ButtonGroup dist; 
   //----- Random Distribution stuff ------------------
   JRadioButton gaussDist, studDist, levyDist;
   
   //---- Levy Distribution stuff ---------------------
   JScrollBar scaleBar,alphaBar,skewBar; 
   JTextField scaleText,alphaText,skewText;
   JSlider[] seedsofChange; JTextField[] seedText;

   JSlider affectRatio; 
   JSlider LagSlide; 

   //---- GARCH stuff --------------------------------
   JScrollBar b0Bar, b1Bar, a1Bar, thetaBar; 
   JTextField b0Text, b1Text, a1Text, thetaText; 
 
   //---- Cycle stuff --------------------------------
   JScrollBar lagBar, lambdaBar, rhoBar; 
   JTextField lagText, lambdaText, rhoText;

   //---- Trend stuff --------------------------------
   JScrollBar trend1Bar, trend2Bar; 
   JTextField trend1Text, trend2Text;

   //---- SV stuff ----------------------------------
   JScrollBar phisvBar, alphasvBar, musvBar;
   JTextField phisvText, alphasvText, musvText, nrepText;
   JComboBox<String>  nrepBar;
 

   //---- Heavy stuff ---------------------------------
   int m2;
   double w1, w2, h_alpha, h_alpha_R; 
   double h_lambda, h_beta, h_beta_R;

   JScrollBar w1Bar, w2Bar, w3Bar;
   JScrollBar h_alphaBar, h_alphaRBar;
   JScrollBar h_betaBar, h_betaRBar;
   JScrollBar h_lambdaBar; 
   JScrollBar m2Bar; 

   JTextField m2Text;
   JTextField w1Text, w2Text;
   JTextField h_alphaText, h_alphaRText;
   JTextField h_betaText, h_betaRText;
   JTextField h_lambdaText; 

   //sharpe_ratio, port_mean, max_drawdown, total_return, total_rank_port, avg_rank, min_rank   
   JLabel sharpeLabel, meanLabel, maxDrawLabel, totalRetLabel, totalRankLabel, avgRankLabel, minRankLabel;
   JTextField sharpeText, meanText, maxDrawText, totalRetText, totalRankText, avgRankText, minRankText;
   
   
   JLabel m2Label;
   JLabel w1Label, w2Label;
   JLabel h_alphaLabel, h_alphaRLabel;
   JLabel h_betaLabel, h_betaRLabel;
   JLabel h_lambdaLabel; 

   //----- Arima stuff -------------------------------
   JComboBox<String> p,q,P,Q,d,D;
   JLabel pLabel, qLabel, dLabel, DLabel, QLabel, PLabel;
   
   boolean first_in = true;
   double scale,alpha,skew,mu;    
   double b0,b1,a1,theta;  
   double lag,lambda,rho; 
   double trend1,trend2;
   JScrollBar muBar; 
   JTextField muText;
   int n_basket = 0;
   double max_drawdown = 0;
   double[] cmaxx;
   double total_rank_port;
   //---- Projectors ----------
   SmallCanvas[] projectors;
   int projW, projH;    
   JScrollPane timeScrollPane;

   //---- Theaters -----------
   JPanel sarimaPanel, garchPanel, cyclePanel, trendPanel;
   JTextField affectText;
   
   boolean cumsum_on = false;
   boolean sigex_returns_on = false; 
   JCheckBox acfplotBox; 
   boolean acfplot;
   JCheckBox resizeBox;    //---- used for scaling/formating data
   boolean resize;         //---- if off, uses raw data

   JButton deleteAll;
   public boolean[] mix;          //-- boolean array indicating which series to compute target
   private boolean[] represent;    //-- boolean array indicating which series to represent in target 
   double[] weight_mix;    //-- array of weights for computing target
   private double[] target_series;
   private int n_sim_series;
   JCheckBox plotTarget;
   boolean gaussian_data = false;
   int which_gaussian = 0;
   int which_vol = 0;
   boolean diff_vol;

   //----- portfolio sliders and weights
   JSlider[] portsliders;
   JCheckBox[] basketCheck;
   JTextField[] portSlidersText;
   double[] port_weights;
   
   double sharpe_ratio;
   double port_mean; 
   double avg_rank;
   double min_rank;
   int length_rank = 10;
   double total_return;
   
   JSlider[] weightSliders;
   JTextField[] weightSlidersText;
   JCheckBox includeTarget; 
   JCheckBox[] repCheck;
   JCheckBox[] simCheck;
   JCheckBox[] priceCheck;
   JLabel[] weightLabel;

   JLabel phisvLabel, alphasvLabel, musvLabel;
   JCheckBox logBox;    //  ---- take log transforms
   boolean[] l_log;     //  ---- keep track of which series have logs
   boolean simData = true; 
   JScrollBar timeScrollBar;


   //----- Market Data stuff  -------  
   boolean vol_set = false;
   ArrayList<double[]> market; 
   ArrayList<double[]> returnsAll; 
   double[] port_returns;
   double[] cum_port_returns;
   
   int market_obs; 
   
   public double[] logprice_data;
   public int n_price = 1;
   boolean priceDataSet = false;
   
   JButton arimaDialogButton, garchDialogButton; 
   JButton cycleDialogButton, trendDialogButton;
   JButton svDialogButton, heavyDialogButton; 

   JDialog arimaDialog, garchDialog; 
   JDialog cycleDialog, trendDialog;
   JDialog svDialog, heavyDialog; 

   //------ Price spread/cointegration stuff
   double spread_mean = 0.0; 
   double spread_std = 0.0; // the current spread_mean and standard dev
   double coint_coeff = 0.0;

   JButton createSpread;
   JButton deleteSpread;
   JButton gaussianizeTarget;
   JButton removeOutliersSeries;
   JButton normalizeWeightsSeries;
   JButton maxSharpeButton;
   boolean spread_exists = false;
   
   int outlier_rank;
   double outlier_scale;
   
   
  //simulate.target_series, simulate.sim_data, simulate.represent, simulate.n_sym
  public SimPanel(int n, int n_r, int burn, JFrame fra)
  {

      frame = fra; 
      setN_obs(n); n_rep = n_r; burnin = burn;       

      rand = 0;  S = 12;  n_sym = 0;

      setSim_data(new ArrayList<double[]>()); sim_orig = new ArrayList<double[]>();
      sim_means  = new ArrayList<double[]>();
      
      returnsAll = new ArrayList<double[]>();
      mix = new boolean[max_series]; setRepresent(new boolean[max_series]); weight_mix = new double[max_series];

      weightSliders = new JSlider[max_series];
      weightSlidersText = new JTextField[max_series];
      repCheck = new JCheckBox[max_series]; simCheck = new JCheckBox[max_series]; weightLabel = new JLabel[max_series]; priceCheck = new JCheckBox[max_series];
      l_log = new boolean[max_series]; 

      for(int i=0;i<max_series;i++) {mix[i] = false; getRepresent()[i] = false; weight_mix[i] = 1.0; l_log[i] = false;}

      Lag = 0; setN_sim_series(0);
      dim = new int[6]; sar_params = new double[6]; var = 1.0;
      sar_params[2] = -.6; sar_params[5] = -.6;
      n_cyc = 1; sv_rep = 1;
     
      sim_check = new boolean[6];
      tseries = new double[getN_obs()];             
      cycle_1 = new double[getN_obs()];               
      cycle_2 = new ArrayList<double[]>();              
      trend = new double[getN_obs()];                  
      garch = new double[getN_obs()];   
      isv = new double[getN_obs()];
      //fsv = new ArrayList<double[]>();
      f1 = 0.02; f2 = 0.02;

      projectors = new SmallCanvas[6];               
      df = new DecimalFormat("##.##");
      df2 = new DecimalFormat("##.###");
      df3 = new DecimalFormat("##.####");
      projW = 400; projH = 140;
      setTheatre(new SimCanvas(400,260,getN_obs(),max_series));

      timeScrollPane = new JScrollPane();
      timeScrollPane.setPreferredSize(new Dimension(400, 260)); 
      timeScrollPane.setViewportView(getTheatre());
      //timeScrollBar = timeScrollPane.getHorizontalScrollBar();


      stationary = new boolean[6]; buildTargetPanel();
      setupControls();
     
      //start heavy model
      heavy = new HeavyModel(getN_obs(), seed[5]);      

  }
 
   /*----------------------------------------------------------------------
     To change certain options 
   -----------------------------------------------------------------------*/
  public void changeBurnin(int b) {burnin = b; resampleSeries();}
  public void changeSeed(int s) {seed[0] = s; seed[1] = s; seed[2] = s; seed[3] = s; seed[4] = s; seed[5] = s; resampleSeries();}      
  public void setNRep(int nr) { if(nr < max_series) n_rep = nr;}
  public void changeRandomDist(int r) {rand = r; changeARIMA(); resampleSeries();}
  public void changeSeasonal(int _s) {S = _s;}
  public void setNobs(int n) 
  {
    //System.out.println(simData + "  "); 
    if(simData) 
    {
      setN_obs(n); getTheatre().plotACFCross(false); 
      getTheatre().setNobs(getN_obs()); clearSeries(); 
      resampleSeries(); acfplotBox.setSelected(false); 
      acfplot = false;
    }
    simData = true;  
    getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300)); 
    timeScrollPane.setViewportView(getTheatre());    
  } 

  public int getNobs() {return getN_obs();}

  public void setGlobalLag(int _l) {Lag = _l;}

  public void addSeries(int i)
  {
     if(getSim_data().size() < max_series)
     {
      
      if(i==0)    //  ARIMA data 
      {
        if(stationary[i])  //---- if stationary, subtract mean and scale by std
        {standardizeAddSeries(tseries);}
        else
        {getSim_data().add(tseries); sim_orig.add(tseries); getTheatre().addSeries(tseries);double[] ms = {0.0,1.0}; sim_means.add(ms); }
        timePlot[getSim_data().size()-1].setEnabled(true); addToMix(getSim_data().size()-1, true);
      }
      else if(i==1)  // Garch data
      {
        if(stationary[i])  //---- if stationary, subtract mean and scale by std
        {standardizeAddSeries(garch);}
        else
        {getSim_data().add(garch); sim_orig.add(garch); getTheatre().addSeries(garch); double[] ms = {0.0,1.0}; sim_means.add(ms); }
        timePlot[getSim_data().size()-1].setEnabled(true); addToMix(getSim_data().size()-1, true);        
      }
      else if(i==2) //  Cycle data 
      {
        if(getSim_data().size() < (max_series-n_cyc))
        {
         
         standardizeAddSeries(cycle_1); 
         timePlot[getSim_data().size()-1].setEnabled(true);
         addToMix(getSim_data().size()-1, true);
         for(int k=0;k<n_cyc;k++)
         {
          standardizeAddSeries(cycle_2.get(k)); 
          timePlot[getSim_data().size()-1].setEnabled(true);
          addToMix(getSim_data().size()-1, true);
         }
        }
      }
      else if(i==3)  //  Trend data
      { 
         getSim_data().add(trend); sim_orig.add(trend);
         double[] ms = {0.0,1.0}; 
         sim_means.add(ms); getTheatre().addSeries(trend); addToMix(getSim_data().size()-1, true);
         timePlot[getSim_data().size()-1].setEnabled(true);
      }
      else if(i==6)  // real data
      {
         standardizeAddSeries(real_series); 
         timePlot[getSim_data().size()-1].setEnabled(true); 
         if(gaussian_data)
         {timePlot[getSim_data().size()-1].setSelected(true);}
         addToMix(getSim_data().size()-1, true);
         
         if(gaussian_data)
         { 
           computeTarget();
           getTheatre().go();
         }
      }
      else if(i==7)
      {
         standardizeAddSeries(volume_series);
         timePlot[getSim_data().size()-1].setEnabled(true); 
         addToMix(getSim_data().size()-1, true);
      }
      else if(i==4) //fsv data
      { 
        double[] temp = new double[getN_obs()]; 
        if(getSim_data().size() < (max_series-sv_rep))
        {
         for(int k=0;k<sv_rep;k++)
         {            
          for(i=0;i<getN_obs();i++) {temp[i] = isv[getN_obs()*k + i];}
 
          standardizeAddSeries(temp); 
          timePlot[getSim_data().size()-1].setEnabled(true);
          addToMix(getSim_data().size()-1, true);
         }
        } 
      }
      else if(i==5) //heavy data
      {

        double[] temp = new double[getN_obs()]; 
        if(getSim_data().size() < (max_series-3))
        {
          for(int k=0;k<n_heavy;k++)
          {            
           for(i=0;i<getN_obs();i++) {temp[i] = heavy_ts[n_heavy*i + k];}
 
           standardizeAddSeries(temp); 
           timePlot[getSim_data().size()-1].setEnabled(true);
           addToMix(getSim_data().size()-1, true);
          }
        }
      }
    }      
    setN_sim_series(getSim_data().size()); deleteAll.setEnabled(true); acfplotBox.setEnabled(true);
    
  }

  public void addToMix(int i, boolean sel)
  {
     //System.out.println("added " + i + ", for " + sel);
     mix[i] = sel; getRepresent()[i] = sel;
     weight_mix[i] = 1.0;  weightSliders[i].setEnabled(true);
     weightSliders[i].setValue(100); 
     weightSlidersText[i].setText(""+1.0); 
     simCheck[i].setSelected(true); simCheck[i].setEnabled(true);
     repCheck[i].setSelected(true); repCheck[i].setEnabled(true); 
     priceCheck[i].setEnabled(true); priceCheck[i].setSelected(false);
     computeTarget();
  }
   
  public void computeTarget()
  {
     setTarget_series(new double[getN_obs()]);     
     setN_sim_series(getSim_data().size());

     for(int i=0;i<getN_sim_series();i++)
     {
       double[] temp;
       if(!l_log[i])
       {temp = scale(sim_orig.get(i),weight_mix[i]);}
       else
       {temp = scale( logTransform(sim_orig.get(i),i), weight_mix[i]);}
       //sim_data.remove(i); sim_data.add(i,temp);
       //theatre.replaceSeries(i,temp);
       getSim_data().set(i,temp); getTheatre().setSeries(i,temp);

       if(mix[i])
       {           
          setTarget_series(mixSeries(getTarget_series(), getSim_data().get(i), weight_mix[i]));        
       }
       
       if(sigex_returns_on && cumsum_on)
       {setTarget_series(cum_port_returns);} // for(i = 0; i<target_series.length;i++) {System.out.println(target_series[i]);}}
     }
     getTheatre().setTarget(getTarget_series());   
     plotTarget.setEnabled(true);  
  }
  

  public void computePriceSeries()
  {
  
    int i,j; int count = 0; 
    priceDataSet = false;
    
    for(i=0;i<getN_sim_series();i++)
    {
     if(priceCheck[i].isSelected())
     {count++;}    
    }
    
    n_price = count;
    logprice_data = new double[n_price*getN_obs()];
     
    int num = 0;
    for(i=0;i<getN_sim_series();i++)
    {
      if(priceCheck[i].isSelected())
      {
       double[] temp = getSim_data().get(i);

       for(j=0;j<getN_obs();j++)
       {logprice_data[getN_obs()*num + j] = temp[j];}
       num++;
      }
    }
 
    if(n_price > 0) {priceDataSet = true;}
  }
  
  public void clearSeries(int ndays)
  {
  
       setN_obs(ndays); 
       clearSeries(); setNobs(ndays); deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       
       getTheatre().setNobs(ndays);
       getTheatre().setPreferredSize(new Dimension(2*ndays, 300));
       timeScrollPane.setViewportView(getTheatre());  
  }
  
  public void setRealSeries(double[] real)
  {
     if(getN_obs() == real.length)
     {
      real_series = new double[real.length];
      System.arraycopy(real,0,real_series,0,getN_obs());
      addSeries(6);
     }
     else
     {System.out.println("The lengths of the input series and n_obs are not the same. Please check");}
  }
  
  public void setLogPriceData(double[] price)
  {
     if(getN_obs() == price.length)
     {
       logprice_data = new double[getN_obs()];
       System.arraycopy(price,0,logprice_data,0,getN_obs());
     }  
     else
     {System.out.println("The lengths of the input series and n_obs are not the same. Please check");}     
  }
  
  
  public static double[] scale(double[] v, double scale)
  {
    double[] temp = new double[v.length];
    for(int i=0;i<v.length;i++) {temp[i] = v[i]*scale;}
    return temp;
  } 

  public void setTarget(double[] target)
  {
      
     if(getN_obs() == target.length) 
     {
         System.arraycopy(target, 0, getTarget_series(), 0, getN_obs()); 
     }
     else
     {System.out.println("Target data must be same length as explaing data");}

  }
  
 
  public void standardizeSeries() //subtract mean and divide by stan_var
  {
    int i=0; int j=0;
    if(resizeBox.isSelected())
    {
      sim_means.clear(); 
      for(j=0;j<sim_orig.size();j++)
      {     
       double[] temp = sim_orig.get(j);
       double[] mstd = mean_std(temp);
       sim_means.add(mstd);

       for(i=0;i<getN_obs();i++) 
       {temp[i] = (temp[i] - mstd[0]);}// /(mstd[1]*mstd[1]);} 
  
       getSim_data().set(j,temp); getTheatre().setSeries(j,temp);  

       //sim_data.add(temp); 
       //theatre.addSeries(temp);

       weight_mix[j] = 1.0;  
       weightSliders[j].setValue(100); 
       weightSlidersText[j].setText(""+1.0); 
      }    
    }
    else
    {    
     sim_means.clear(); 
     for(j=0;j<sim_orig.size();j++)
     {     
       double[] mstd = {0,1}; sim_means.add(mstd);
       double[] temp = sim_orig.get(j);       
       getSim_data().set(j,temp); getTheatre().setSeries(j,temp); 
     }
    }
    computeTarget();
    getTheatre().go();
  }

  public void standardizeAddSeries(double[] v)
  {
     int diff;
     if(v.length >= getN_obs())
     {
      diff = v.length - getN_obs(); 
      /*if(resize)
      {
      
       double std; double[] ms; double[] temp = new double[n_obs];
       ms = mean_std(v); sim_means.add(ms); 
       std = ms[1]*ms[1];

       for(int i=0;i<n_obs;i++) 
       {temp[i] = (v[diff+i] - ms[0])/std;}// System.out.println(temp[i]);}// System.out.println(temp[i]);}
       //System.out.println("Plotted\n");
       sim_data.add(temp); sim_orig.add(temp); theatre.addSeries(temp);  
      }
      else
      {*/
        double[] ms = {0.0,1.0}; sim_means.add(ms);
        double[] temp = new double[getN_obs()];
        double std = ms[1]*ms[1];
        for(int i=0;i<getN_obs();i++) 
        {temp[i] = (v[diff+i] - ms[0])/std;}// System.out.println(temp[i]);}// System.out.println(temp[i]);}
        //System.out.println("Plotted\n");
        getSim_data().add(temp); sim_orig.add(temp); getTheatre().addSeries(temp);  
      //}
     }
     else
     {System.out.println("New data must be at least length " + getN_obs());}   
  }


  /*----------------------------------------------
    Apply log transforms to all series that are plotted   
  -----------------------------------------------*/
  public void applyLogTransform()
  {
    int i;   

    for(i = 0; i < getSim_data().size(); i++)
    {
        if(timePlot[i].isSelected())
        {
           System.out.println("log transform");
           double[] temp = logTransform(sim_orig.get(i),i);
                 
           //for(j=0;j<n_obs;j++) {System.out.println(temp[j]);}

           getSim_data().set(i,temp); getTheatre().setSeries(i,temp);   
    
           l_log[i] = true;
        }
    }
    computeTarget();
    getTheatre().go();
  }

  public void applyExpTransform()
  {
    int i;  
    for(i = 0; i < getSim_data().size(); i++)
    {
        if(timePlot[i].isSelected() && l_log[i])
        {
           double[] temp = sim_orig.get(i);
           getSim_data().set(i,temp); getTheatre().setSeries(i,temp);      
           l_log[i] = false;
        }
    }
    getTheatre().go();
  }


  public void deleteSeries(int i)
  {     
        int k;
        timePlot[getSim_data().size()-1].setEnabled(false); timePlot[getSim_data().size()-1].setSelected(false);
        getSim_data().remove(i); sim_means.remove(i); sim_orig.remove(i);
        getSim_data().trimToSize(); sim_means.trimToSize(); sim_orig.trimToSize();
        getTheatre().deleteSeries(i);   
        for(k=i;k<getN_sim_series()-1;k++)
        {mix[k] = mix[k+1]; getRepresent()[k] = getRepresent()[k+1];} 
        mix[k] = false; getRepresent()[k] = false; 
        setN_sim_series(getSim_data().size()); 
                   
        computeTarget();
  }

  public void clearSeries()
  {
    getSim_data().clear(); sim_means.clear(); sim_orig.clear();
    setN_sim_series(getSim_data().size());
    setTarget_series(new double[getN_obs()]);
 
    for(int i=0;i<max_series;i++) 
    {
      timePlot[i].setSelected(false); 
      timePlot[i].setEnabled(false); 
      mix[i] = false; getRepresent()[i] = false;
      simCheck[i].setSelected(false); 
      repCheck[i].setSelected(false);
      priceCheck[i].setSelected(false);
      weightSliders[i].setEnabled(false); simCheck[i].setEnabled(false); repCheck[i].setEnabled(false); priceCheck[i].setEnabled(false);
      l_log[i] = false;
    }
 
    plotTarget.setEnabled(false); 
    getTheatre().clearSeries(); 
     
     plotTarget.setEnabled(false);
  }


  public void integrateGarch()
  {
    if(sim_check[1] == true)
    {
       double[] temp = cumsum(garch, garch.length);            
       System.arraycopy(temp, 0, garch, 0, garch.length);
       projectors[1].setTseries(garch, getN_obs(), 1); projectors[1].go();
       stationary[1] = false;
    }  
  }

  public void integrateSV()
  {
    int i;
    if(sim_check[4] == true)
    {
      double[] temp = new double[getN_obs()];
      double[] temp2 = new double[getN_obs()];
      for(int k=0;k<sv_rep;k++)
      { 
       for(i=0;i<getN_obs();i++) {temp[i] = isv[getN_obs()*k + i];}
       temp2 = cumsum(temp, temp.length);      
       System.arraycopy(temp2, 0, isv, getN_obs()*k, getN_obs());
      }
      projectors[4].setTseries(isv, getN_obs(), sv_rep); projectors[4].go();
      stationary[4] = false;
    }  
  }

  public void integrateHEAVY()
  {
    int i,k;
    if(sim_check[5] == true && n_heavy == 2)
    {
      n_heavy = 3;
      double[] temp = new double[getN_obs()];
      double[] temp2 = new double[getN_obs()];
      double[] tempisv = new double[getN_obs()*n_heavy];
 
      for(k=0;k<2;k++)
      { 
        for(i=0;i<getN_obs();i++)
        {
          temp[i] = heavy_ts[2*i];
          tempisv[n_heavy*i+k] = heavy_ts[2*i + k];
        }
      }   
      temp2 = cumsum(temp, temp.length);

      heavy_ts = new double[getN_obs()*n_heavy];
      for(k=0;k<2;k++)
      { 
        for(i=0;i<getN_obs();i++)
        {
          heavy_ts[n_heavy*i + k] = tempisv[n_heavy*i + k];
        }
      }

      for(i=0;i<getN_obs();i++)
      {
          heavy_ts[n_heavy*i + 2] = temp2[i];
      }   

      for(k=0;k<n_heavy;k++) 
      {
          for(i=0;i<getN_obs();i++)
          {tempisv[k*getN_obs() + i] = heavy_ts[n_heavy*i + k];}
      }

      projectors[5].setTseries(tempisv, getN_obs(), n_heavy); projectors[5].go();
      stationary[5] = false;
    }
  }

  public void addTrend()
  {
     int i=0; stationary[0] = false;
     for(i=0;i<Lag;i++) {tseries[i] = tseries[i] + affect*trend[0];}            
     for(i=0;i<getN_obs()-Lag;i++)  {tseries[i+Lag] = tseries[i+Lag] + affect*trend[i];}
     projectors[0].setTseries(tseries, getN_obs(), 1); projectors[0].go();
  }


   /*----------------------------------------------------------------------
     To change parameters n_obs
   -----------------------------------------------------------------------*/
  
  public void setGarchParams(double a, double b, double c, double d)
  {
    b0 = a; b1 = b; a1 = c; theta = d; sampleSeries(1);
  }
  public void setLevy(double a, double b, double c)
  {
    scale=a;alpha=b;skew=c; if(rand == 2) {changeARIMA(); sampleSeries(0);}
  }
  public void setCycleParams(double a, double b, double c)
  {
    lag=a; lambda=b; rho=c;  sampleSeries(2);
  }
  public void setFSVParams(double phi, double alpha, double mu, int nrep)
  {
    phisv = phi; alphasv = alpha; musv = mu; sv_rep = nrep; sampleSeries(4);
  }
  public void setTrendParams(double a, double b)
  {
    trend1 = a; trend2 = b; sampleSeries(3);
  }
  public void setMu(double a)
  {
    mu = a; if(rand == 1) {changeARIMA(); sampleSeries(0);}
  }
  public void setVar(double a)
  {
    var = a; if(rand == 0) {changeARIMA(); sampleSeries(0);}
  }
  public void setHeavyParameters(int _m2, double _w1, double _w2, double _alpha, double _alpha_R, double _beta, double _beta_R, double _lambda)
  {
     m2 = _m2; w1 = _w1; w2 = _w2; h_alpha = _alpha; h_alpha_R = _alpha_R; 
     h_beta = _beta; h_beta_R = _beta_R; h_lambda = _lambda; 

     sampleSeries(5); 
  }


  /*----------------------------------------------------------------------
      Transforms 
  -----------------------------------------------------------------------*/
  
  public double[] logTransform(double[] ts, int x)
  {
    //if(!l_log[x])
    //{
     int i; double min = 1000000;
     double[] logData = new double[ts.length];


     for(i=0; i < ts.length; i++) {if(ts[i] < min) {min = ts[i];}}

     //---- if less than 0----------------------------
     if(min <= 0)
     {for(i=0; i < ts.length; i++) {logData[i] = Math.log(ts[i]-min+10.0) + min-10.0;}}         
     else
     {for(i=0; i < ts.length; i++) {logData[i] = Math.log(ts[i]);}}

     return logData; 
    //}
    //else 
    //{return ts;}  
  }

  public double[] expTransform(double[] ts)
  {
     int i;
     double[] expData = new double[ts.length];
     for(i=0; i < ts.length; i++) {expData[i] = Math.exp(ts[i]);}    
     return expData;
  }

  /*--------------------------------------------------------
    Sets the ACFCross series for two given series
    Default is 
      
  */
 
  public void setACFCross(int s1, int s2)
  {
     if(s1 == 0 && s2 != 0) 
     {
       plotACFCross(getTarget_series(), getSim_data().get(s2-1));
     }
     else if(s1 == 0 && s2 == 0) 
     {
       plotACFCross(getTarget_series(), getTarget_series());
     }  
     else if(s1 != 0 && s2 == 0) 
     {
       plotACFCross(getSim_data().get(s1-1), getTarget_series());
     }    
     else 
     {
       plotACFCross(getSim_data().get(s1-1), getSim_data().get(s2-1));
     }
  }


  public void plotACFCross(double[] v1, double[] v2)
  {
      int i;
      double[] output = crossCorrelation(getN_obs(), v1, v2);
      double[] acf_v1 = new double[getN_obs()];
      double[] acf_v2 = new double[getN_obs()];
      double[] cc12 = new double[getN_obs()];
      double[] cc21 = new double[getN_obs()]; 

      for(i=0; i<getN_obs(); i++)
      { 
        acf_v1[i] = output[i];
        cc12[i] = output[getN_obs() + i];
        cc21[i] = output[2*getN_obs() + i];       
        acf_v2[i] = output[3*getN_obs() + i];        
      }
      getTheatre().setACFCross(acf_v1, acf_v2, cc12, cc21);
   }


  /*-------------------------------------------------------------------------
     Setup controls       
  --------------------------------------------------------------------------*/
  public void setupControls()
  {
   int i;
   sim_box = new JCheckBox[6]; sim_check = new boolean[6];
   sim_box[0] = new JCheckBox("SARIMA",false); sim_box[1] = new JCheckBox("GARCH",false);
   sim_box[2] = new JCheckBox("Cycle",false); sim_box[3] = new JCheckBox("Trend",false);
   sim_box[4] = new JCheckBox("Factor SV",false); sim_box[5] = new JCheckBox("HEAVY",false);

   for(i=0;i<6;i++)
   {
     sim_box[i].setHorizontalTextPosition(JMenuItem.LEFT); 
     sim_box[i].addItemListener(new MyItemListener());
     sim_check[i] = false; 
   }  

   addMe = new JButton[6];
   for(i=0;i<6;i++) {addMe[i] = new JButton("Add"); addMe[i].addActionListener(new MyActionListener());}

   integrate = new JButton("Integrate"); integrate.addActionListener(new MyActionListener());
   addLagTrend = new JButton("Trend"); addLagTrend.addActionListener(new MyActionListener()); addLagTrend.setEnabled(false);
   affectRatio = new JSlider(JSlider.HORIZONTAL, 1, 100, 100); affectRatio.setExtent(1);
   affectRatio.addChangeListener(new MyChangeListener()); affect = affectRatio.getValue()*.01;
  
   integratesv = new JButton("Integrate"); integratesv.addActionListener(new MyActionListener());   
   integrateHeavy = new JButton("Integrate"); integrateHeavy.addActionListener(new MyActionListener()); 

   seedsofChange = new JSlider[6]; seedText = new JTextField[6];
   for(i=0;i<6;i++) 
   {
       seedsofChange[i] = new JSlider(JSlider.HORIZONTAL, 0, 500, 3); seedsofChange[i].setExtent(1);
       seedsofChange[i].addChangeListener(new MyChangeListener()); seedsofChange[i].setPreferredSize(new Dimension(60, 15));  

       seedText[i] = new JTextField(3); seedText[i].setText("" + df.format(seedsofChange[i].getValue()));
   }
   affectText = new JTextField(3); 
   affectText.setText(""+df.format(affect));


   //----- Random Distribution stuff ------------------


   //true_phi=0.8; true_sigma_eta=0.1; true_mu=1.0;      
   phisvBar = new JScrollBar(JScrollBar.HORIZONTAL, 10, 8, 0, 200);
   alphasvBar = new JScrollBar(JScrollBar.HORIZONTAL, 10, 0, 0, 200);
   musvBar = new JScrollBar(JScrollBar.HORIZONTAL, 10, 10, 0, 200);

   String[] items = new String[]{"1", "2", "3", "4"};
   nrepBar = new JComboBox<String>(items);
   nrepBar.addActionListener(new MyActionListener());
   
   phisv = 0.8; alphasv = 0.1; musv = 1.0; 
   phisvText = new JTextField(4); phisvText.setText(""+df.format(phisv));
   alphasvText = new JTextField(4); alphasvText.setText(""+df.format(alphasv));
   musvText = new JTextField(4); musvText.setText(""+df.format(musv));

   phisvLabel = new JLabel("\u03D5"); 
   alphasvLabel = new JLabel("\u03C3");
   musvLabel = new JLabel("\u03BC");

   dist = new ButtonGroup();
   gaussDist = new JRadioButton("Normal",true);      gaussDist.setHorizontalTextPosition(JMenuItem.LEFT);
   studDist = new JRadioButton("Student-t");    studDist.setHorizontalTextPosition(JMenuItem.LEFT);
   levyDist = new JRadioButton("Levy-Stable");  levyDist.setHorizontalTextPosition(JMenuItem.LEFT);
   dist.add(gaussDist); dist.add(studDist);  dist.add(levyDist);
     
   gaussDist.addItemListener(new MyItemListener());
   studDist.addItemListener(new MyItemListener());
   levyDist.addItemListener(new MyItemListener());

   scaleBar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,0,20); 
   alphaBar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,0,20);    
   skewBar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,0,20);    

   b0Bar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,1,200);
   b1Bar = new JScrollBar(JScrollBar.HORIZONTAL,30,1,1,100);  
   a1Bar = new JScrollBar(JScrollBar.HORIZONTAL,60,1,1,100);   
   thetaBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,100);  
 
   lagBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,100);  
   lambdaBar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,0,100);     
   rhoBar = new JScrollBar(JScrollBar.HORIZONTAL,80,1,1,101); 

   trend1Bar = new JScrollBar(JScrollBar.HORIZONTAL,7,1,1,100);  
   trend2Bar = new JScrollBar(JScrollBar.HORIZONTAL,95,1,1,300); 

   muBar = new JScrollBar(JScrollBar.HORIZONTAL,20,1,10,40);

   w1Bar = new JScrollBar(JScrollBar.HORIZONTAL,5,1,0,100);
   w2Bar = new JScrollBar(JScrollBar.HORIZONTAL,15,1,0,100);
   h_alphaBar = new JScrollBar(JScrollBar.HORIZONTAL,20,1,0,100);
   h_alphaRBar = new JScrollBar(JScrollBar.HORIZONTAL,40,1,10,100);
   h_betaBar = new JScrollBar(JScrollBar.HORIZONTAL,60,1,10,100);
   h_betaRBar = new JScrollBar(JScrollBar.HORIZONTAL,50,1,10,100);
   h_lambdaBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,100);
   m2Bar = new JScrollBar(JScrollBar.HORIZONTAL,72,1,1,100);

   w1Label = new JLabel("\u03BC");
   w2Label = new JLabel("\u03BC");
   h_alphaLabel = new JLabel("\u03B1");
   h_alphaRLabel = new JLabel("\u03B1_HF");
   h_betaLabel = new JLabel("\u03D0");
   h_betaRLabel = new JLabel("\u03D0_HF");
   h_lambdaLabel = new JLabel("\u03BB");
   m2Label = new JLabel("HF");

   w1 = .05; w2 = .15; h_alpha = .20; h_alpha_R = .40; 
   h_beta = .60; h_beta_R = .50; h_lambda = 0; m2 = 72; 

   w1Text = new JTextField(4); w1Text.setText(""+df.format(w1));
   w2Text = new JTextField(4); w2Text.setText(""+df.format(w2));
   h_alphaText = new JTextField(4); h_alphaText.setText(""+df.format(h_alpha));
   h_alphaRText = new JTextField(4); h_alphaRText.setText(""+df.format(h_alpha_R)); 
   h_betaText = new JTextField(4); h_betaText.setText(""+df.format(h_beta));
   h_betaRText = new JTextField(4); h_betaRText.setText(""+df.format(h_beta_R));
   h_lambdaText = new JTextField(4); h_lambdaText.setText(""+df.format(h_lambda));
   m2Text = new JTextField(4); m2Text.setText(""+df.format(m2));

   



   
   scaleText = new JTextField(3);  alphaText = new JTextField(3); skewText = new JTextField(3);
   b0Text = new JTextField(3); b1Text = new JTextField(3); a1Text = new JTextField(3); thetaText = new JTextField(3); 
   lagText = new JTextField(3); lambdaText = new JTextField(3); rhoText = new JTextField(3);
   trend1Text = new JTextField(3); trend2Text = new JTextField(3); muText = new JTextField(3);







   timePlot = new JCheckBox[max_series];  
   plotTarget = new JCheckBox("Target:"); plotTarget.setSelected(false); plotTarget.setEnabled(false);
   plotTarget.setHorizontalTextPosition(JMenuItem.LEFT);
   plotTarget.addItemListener(new MyItemListener());

   for(i=0;i<max_series;i++)
   {
    //timePlot[i] = new JCheckBox(" X_"+Integer.toString(i+1)+"(t):");
    timePlot[i] = new JCheckBox(" "+Integer.toString(i+1) + ":");
    timePlot[i].setSelected(false); timePlot[i].setEnabled(false);
    timePlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
    timePlot[i].addItemListener(new MyItemListener());
   }


   logBox = new JCheckBox("Log"); logBox.setHorizontalTextPosition(JMenuItem.LEFT); 
   logBox.setSelected(false); logBox.addItemListener(new MyItemListener());       


   acfplotBox = new JCheckBox("ACF/CC"); acfplotBox.setHorizontalTextPosition(JMenuItem.LEFT); 
   acfplotBox.setSelected(false); acfplotBox.addItemListener(new MyItemListener()); acfplot = false;
   acfplotBox.setEnabled(false);  


   resizeBox = new JCheckBox("Rescale"); resizeBox.setHorizontalTextPosition(JMenuItem.LEFT); 
   resizeBox.setSelected(false); resizeBox.addItemListener(new MyItemListener());
   resize = false;       

   deleteAll = new JButton("Delete"); deleteAll.setEnabled(false); 
   deleteAll.addActionListener(new MyActionListener());

   AdjustmentListener al = new AdjustmentListener()  
   {
        public void adjustmentValueChanged(AdjustmentEvent e) 
        {
         if(e.getAdjustable() == muBar)
         {
           mu = muBar.getValue()*.1; muText.setText(""+df.format(mu)); setMu(mu); 
         } 
         else if(e.getAdjustable() == scaleBar)
         {
           scale = scaleBar.getValue()*.1; scaleText.setText(""+df.format(scale)); setLevy(scale,alpha,skew); 
         }
         else if(e.getAdjustable() == alphaBar)
         { 
           alpha = alphaBar.getValue()*.1; alphaText.setText(""+df.format(alpha)); setLevy(scale,alpha,skew);  
         }
         else if(e.getAdjustable() == skewBar)
         {             
           skew = (skewBar.getValue()-10)*.1; skewText.setText(""+df.format(skew)); setLevy(scale,alpha,skew);           
         }
         else if(e.getAdjustable() == b0Bar)
         {             
           b0 = b0Bar.getValue()*.01; b0Text.setText(""+df.format(b0)); setGarchParams(b0,b1,a1,theta);           
         }
         else if(e.getAdjustable() == b1Bar)
         {             
           b1 = b1Bar.getValue()*.01; b1Text.setText(""+df.format(b1)); setGarchParams(b0,b1,a1,theta); if(gaussian_data) {gaussianizeData(b1,a1);}
         }
         else if(e.getAdjustable() == a1Bar)
         {             
           a1 = a1Bar.getValue()*.01; a1Text.setText(""+df.format(a1)); setGarchParams(b0,b1,a1,theta); if(gaussian_data) {gaussianizeData(b1,a1);}  
         }
         else if(e.getAdjustable() == thetaBar)
         {             
           theta = thetaBar.getValue()*.01; thetaText.setText(""+df.format(theta)); setGarchParams(b0,b1,a1,theta);           
         }
         else if(e.getAdjustable() == lagBar)
         {             
           lag = lagBar.getValue()*.1; lagText.setText(""+df.format(lag)); setCycleParams(lag,lambda,rho);           
         }
         else if(e.getAdjustable() == lambdaBar)
         {             
           lambda = lambdaBar.getValue()*.01; lambdaText.setText(""+df.format(lambda)); setCycleParams(lag,lambda,rho);          
         }
         else if(e.getAdjustable() == rhoBar)
         {             
           rho = rhoBar.getValue()*.01; rhoText.setText(""+df.format(rho)); setCycleParams(lag,lambda,rho);           
         }
         else if(e.getAdjustable() == trend1Bar)
         {             
           trend1 = trend1Bar.getValue()*.01; trend1Text.setText(""+df.format(trend1)); setTrendParams(trend1,trend2);           
         }
         else if(e.getAdjustable() == trend2Bar)
         {             
           trend2 = -1.0*trend2Bar.getValue()*.01; trend2Text.setText(""+df.format(trend2)); setTrendParams(trend1,trend2);             
         }
         else if(e.getAdjustable() == phisvBar)
         {
           phisv = phisvBar.getValue()*.01; phisvText.setText(""+df.format(phisv)); setFSVParams(phisv, alphasv, musv, sv_rep);
         }
         else if(e.getAdjustable() == musvBar)
         {
           musv = musvBar.getValue()*.1; musvText.setText(""+df.format(musv)); setFSVParams(phisv, alphasv, musv, sv_rep);
         }
         else if(e.getAdjustable() == alphasvBar)
         {
           alphasv = alphasvBar.getValue()*.01; alphasvText.setText(""+df.format(alphasv)); setFSVParams(phisv, alphasv, musv, sv_rep);
         }
         else if(e.getAdjustable() == h_alphaBar)
         {
           h_alpha = h_alphaBar.getValue()*.01; h_alphaText.setText(""+df.format(h_alpha));  
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == h_alphaRBar)
         {
           h_alpha_R = h_alphaRBar.getValue()*.01; h_alphaRText.setText(""+df.format(h_alpha_R)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == h_betaBar)
         {
           h_beta = h_betaBar.getValue()*.01; h_betaText.setText(""+df.format(h_beta)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == h_betaRBar)
         {
           h_beta_R = h_betaRBar.getValue()*.01; h_betaRText.setText(""+df.format(h_beta_R)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda); 
         }
         else if(e.getAdjustable() == h_lambdaBar)
         {
           h_lambda = h_lambdaBar.getValue()*.01; h_lambdaText.setText(""+df.format(h_lambda)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == m2Bar)
         {
           m2 = m2Bar.getValue(); m2Text.setText(""+m2); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == w1Bar)
         {
           w1 = w1Bar.getValue()*.01; w1Text.setText(""+df.format(w1)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
         else if(e.getAdjustable() == w2Bar)
         {
           w2 = w2Bar.getValue()*.01; w2Text.setText(""+df.format(w2)); 
           setHeavyParameters(m2, w1, w2, h_alpha, h_alpha_R, h_beta, h_beta_R, h_lambda);
         }
     }
    };
   
   muBar.addAdjustmentListener(al);
   scaleBar.addAdjustmentListener(al);  
   alphaBar.addAdjustmentListener(al); 
   skewBar.addAdjustmentListener(al);    

   b0Bar.addAdjustmentListener(al); 
   b1Bar.addAdjustmentListener(al); 
   a1Bar.addAdjustmentListener(al);   
   thetaBar.addAdjustmentListener(al); 
 
   lagBar.addAdjustmentListener(al); 
   lambdaBar.addAdjustmentListener(al);      
   rhoBar.addAdjustmentListener(al); 

   trend1Bar.addAdjustmentListener(al);   
   trend2Bar.addAdjustmentListener(al); 
   
   phisvBar.addAdjustmentListener(al);
   musvBar.addAdjustmentListener(al);
   alphasvBar.addAdjustmentListener(al);   


   w1Bar.addAdjustmentListener(al); 
   w2Bar.addAdjustmentListener(al); 
   h_alphaBar.addAdjustmentListener(al);
   h_alphaRBar.addAdjustmentListener(al);
   h_betaBar.addAdjustmentListener(al);
   h_betaRBar.addAdjustmentListener(al);
   h_lambdaBar.addAdjustmentListener(al); 
   m2Bar.addAdjustmentListener(al); 



   
   




   initiateParameters();
   setupProjectors(projW, projH);  
   setComboBoxes();
   setupPanels();

   }
   //-------------- Initiate parameters -------------------------
   private void initiateParameters()
   {
           mu = muBar.getValue()*.1; muText.setText(""+df.format(mu)); 
           scale = scaleBar.getValue()*.1; scaleText.setText(""+df.format(scale)); 
           alpha = alphaBar.getValue()*.1; alphaText.setText(""+df.format(alpha));   
           skew = (skewBar.getValue()-10)*.1; skewText.setText(""+df.format(skew));            
           b0 = b0Bar.getValue()*.01; b0Text.setText(""+df.format(b0));         
           b1 = b1Bar.getValue()*.01; b1Text.setText(""+df.format(b1));          
           a1 = b0Bar.getValue()*.01; a1Text.setText(""+df.format(a1));            
           theta = thetaBar.getValue()*.01; thetaText.setText(""+df.format(theta));         
           lag = lagBar.getValue()*.1; lagText.setText(""+df.format(lag));         
           lambda = lambdaBar.getValue()*.01; lambdaText.setText(""+df.format(lambda));           
           rho = rhoBar.getValue()*.01; rhoText.setText(""+df.format(rho));         
           trend1 = trend1Bar.getValue()*.01; trend1Text.setText(""+df.format(trend1));         
           trend2 = -1.0*trend2Bar.getValue()*.01; trend2Text.setText(""+df.format(trend2));                   
   }

   private void buildTargetPanel()
   {
     int i;
     for(i=0; i<max_series; i++)
     {
       weightSliders[i] = new JSlider(JSlider.VERTICAL,0,100,0);      
       weightSliders[i].addChangeListener(new MySliderListener());  
       weightSliders[i].setToolTipText("The weight attached to simulated series "+ Integer.toString(i+1)+" explaing the target series");  
       weightSliders[i].setEnabled(false); weightSliders[i].setPreferredSize(new Dimension(20,100));

       weightSlidersText[i] = new JTextField(3); 
       weightSlidersText[i].setText(""+1.0);

       repCheck[i] = new JCheckBox("",false);
       simCheck[i] = new JCheckBox("",false);
       priceCheck[i] = new JCheckBox("",false);
       
       priceCheck[i].addItemListener(new MyItemListener());
       repCheck[i].addItemListener(new MyItemListener());
       simCheck[i].addItemListener(new MyItemListener());

       weightLabel[i] = new JLabel(" X_"+Integer.toString(i+1)+"(t):");

     }
     
     createSpread = new JButton("Create Spread");
     deleteSpread = new JButton("Delete Spread");
     gaussianizeTarget = new JButton("Gaussianize Target");
     removeOutliersSeries = new JButton("Remove Outliers");
     normalizeWeightsSeries = new JButton("Normalize Weights");
     maxSharpeButton = new JButton("MaxSharpe");
     
     createSpread.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
               computePriceSeriesSpread();      
            }
     }); 
     deleteSpread.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
               if(spread_exists) {deleteSeries(getSim_data().size()-1);}
               spread_exists = false;
            }
     });  
     
     gaussianizeTarget.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
               gaussianizeData(.6, .05);
            }
     });
     
     removeOutliersSeries.addActionListener(new ActionListener() { 
        public void actionPerformed(ActionEvent evt) { 
          //removeOutliers(2,.2);
       }
     }); 
     
     normalizeWeightsSeries.addActionListener(new ActionListener() { 
        public void actionPerformed(ActionEvent evt) { 
           normalizeWeights();
           computePortfolio();
           computeTarget();    
       }
     }); 
     
     maxSharpeButton.addActionListener(new ActionListener() { 
        public void actionPerformed(ActionEvent evt) { 
           maximizeSharpe();  
       }
     });  
     
    sharpeLabel = new JLabel("Sharpe:");  
    meanLabel = new JLabel("Mean:");  totalRankLabel = new JLabel("Rank:");
    maxDrawLabel = new JLabel("Drawdown:"); totalRetLabel = new JLabel("TrRatio:");  
    avgRankLabel = new JLabel("AnyaRatio:"); minRankLabel = new JLabel("minRank:"); 
    
    sharpeText = new JTextField(6);  meanText = new JTextField(6); 
    maxDrawText = new JTextField(6); totalRetText = new JTextField(6); 
    totalRankText = new JTextField(6); avgRankText = new JTextField(6); 
    minRankText = new JTextField(6);     
     
     

     
   }

   private void setupProjectors(int w, int h)
   {     
     for(int i = 0; i < 6; i++)
     {
        projectors[i] = new SmallCanvas(w,h,getN_obs(),i);
     }
   }

   private void setComboBoxes()
   {
         int i;
         p = new JComboBox<String>();
         d = new JComboBox<String>();
         q = new JComboBox<String>();
         P = new JComboBox<String>();
         D = new JComboBox<String>();
         Q = new JComboBox<String>();

         for(i=0;i<3;i++) {p.addItem(Integer.toString(i));} p.setSelectedIndex(0);
         for(i=0;i<3;i++) {q.addItem(Integer.toString(i));} q.setSelectedIndex(1);
         for(i=0;i<2;i++) {P.addItem(Integer.toString(i));} P.setSelectedIndex(0);
         for(i=0;i<2;i++) {Q.addItem(Integer.toString(i));} Q.setSelectedIndex(1);
       
          D.addItem("0"); D.addItem("1"); D.setSelectedIndex(0);     
          d.addItem("0"); d.addItem("1"); d.setSelectedIndex(0);
       

          ActionListener comboListener = new ActionListener()  
          {
           public void actionPerformed(ActionEvent e)
           {
            if((d.getSelectedIndex() > 0) || (D.getSelectedIndex() > 0)) 
            {addLagTrend.setEnabled(false);}
            else
            {addLagTrend.setEnabled(true);}

            setDimensions(p.getSelectedIndex(),d.getSelectedIndex(),q.getSelectedIndex(),
                         P.getSelectedIndex(),D.getSelectedIndex(),Q.getSelectedIndex());
           } 
          };

     p.addActionListener(comboListener); q.addActionListener(comboListener);
     d.addActionListener(comboListener); D.addActionListener(comboListener);
     Q.addActionListener(comboListener); P.addActionListener(comboListener);

     setDimensions(p.getSelectedIndex(),d.getSelectedIndex(),q.getSelectedIndex(),
                   P.getSelectedIndex(),D.getSelectedIndex(),Q.getSelectedIndex());

   }

 
   //----------------- Change Sarima Stuff -------------------------------------
   public void setDimensions(int p, int d, int q, int P, int D, int Q)
   { 
      dim[0] = p; dim[1] = d; dim[2] = q; dim[3] = P; dim[4] = D; dim[5] = Q;
      changeARIMA();
   }

   public void setSarimaParams(double a, double b, double c, double d, double e, double f) 
   {
      
      sar_params[0] = a; sar_params[1] = b; //ar
      sar_params[2] = c; sar_params[3] = d; //ma
      sar_params[4] = e; sar_params[5] = f; //sar,sma
      
      changeARIMA();
      
   }

   //--------- If anything changes besides n_obs, burnin, seed, apply this ----------------------
   public void changeARIMA()
   {
         int i;
         n_params = dim[0] + dim[2] + dim[3] + dim[5] + 1; 
         if(rand == 2) {params = new double[n_params+2];} 
         else {params = new double[n_params];} 

         for(i=0; i < dim[0]; i++) {params[i] = sar_params[i];}
         for(i=0; i < dim[3]; i++) {params[i+dim[0]] = sar_params[i+2+2];}
         for(i=0; i < dim[2]; i++) {params[i+dim[0]+dim[3]] = sar_params[i+2];}
         for(i=0; i < dim[5]; i++) {params[i+dim[0]+dim[3]+dim[2]] = sar_params[i+2+2+1];}
         if(rand==0) {params[n_params-1] = var;}
         else if(rand==1) {params[n_params-1] = mu;} 
         else if(rand==2) {params[n_params-1] = scale; params[n_params] = alpha; params[n_params] = skew;} 

         if(sim_box[0].isSelected()) {sampleSeries(0);}
   }

   /*---------------------------------------------------------------------
      Recomputes all series with new parameters
   ----------------------------------------------------------------------*/
   public void resampleSeries()
   {
      int i,k;
      if(sim_check[0]) //arima
      {        
         
         stationary[0] = true;
         tseries = sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
         projectors[0].setTseries(tseries, getN_obs(), 1); projectors[0].go();
         
         if((dim[1] > 0) || (dim[4] > 0)) {stationary[0] = false;} 
      }
      if(sim_check[1]) //garch
      {
        stationary[1] = true;
        garch = simulateNGARCH(b1, a1, theta, b0);
        projectors[1].setTseries(garch, getN_obs(), 1); projectors[1].go();
      }
      if(sim_check[2]) //cycles
      {
         double[] ts = new double[getN_obs()];
         stationary[2] = true; cycle_2.clear(); 
         double[] temp = simulateCycles(n_cyc, lag, lambda, rho);
         cycle_1 = new double[getN_obs()];  
         for(i=0;i<getN_obs();i++) {cycle_1[i] = temp[i];}
         for(k=0;k<n_cyc;k++) 
         {
          for(i=0;i<getN_obs();i++) {ts[i] = temp[i+(k+1)*getN_obs()];}
          cycle_2.add(ts);
         }
         projectors[2].setTseries(temp, getN_obs(), 2); projectors[2].go();
      }
      if(sim_check[3]) //trend
      {
         stationary[3] = false;
         trend = simulateTrend(trend1, trend2, var);
         projectors[3].setTseries(trend, getN_obs(), 1); projectors[3].go();
      }
      if(sim_check[4]) //factor sv model
      {
         stationary[4] = true; 
         isv = simulateFSV(sv_rep, phisv, alphasv, musv);

         projectors[4].setTseries(isv, getN_obs(), sv_rep); projectors[4].go();       
      }
   }

   /*---------------------------------------------------------------------
      Recomputes specific series with new parameters
   ----------------------------------------------------------------------*/
   public void sampleSeries(int x)
   {
      int i,k;
      if(x==0) //arima
      {        
         stationary[0] = true; if((dim[1] > 0) || (dim[4] > 0)) {stationary[0] = false;}
         tseries = sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
         projectors[0].setTseries(tseries, getN_obs(), 1); projectors[0].go();
          
      }
      else if(x==1) //garch
      {
        stationary[1] = true;
        garch = simulateNGARCH(b1, a1, theta, b0);  
        projectors[1].setTseries(garch, getN_obs(), 1); projectors[1].go();
      }
      else if(x==2) //cycles
      {
         cycle_2.clear(); 
         stationary[2] = true; double[] ts = new double[getN_obs()];
         double[] temp = simulateCycles(n_cyc, lag, lambda, rho); 
         cycle_1 = new double[getN_obs()]; 
         for(i=0;i<getN_obs();i++) {cycle_1[i] = temp[i];} 

         for(k=0;k<n_cyc;k++) 
         {
          for(i=0;i<getN_obs();i++) {ts[i] = temp[i+(k+1)*getN_obs()];}
          cycle_2.add(ts);
         }
         projectors[2].setTseries(temp, getN_obs(), 2); projectors[2].go();
      }
      else if(x==3) //trend
      {
         stationary[3] = false;
         trend = simulateTrend(trend1, trend2, var);
         projectors[3].setTseries(trend, getN_obs(), 1); projectors[3].go();
      }
      else if(x==4) //fsv/isv
      {
         stationary[4] = true; 
         isv = simulateFSV(sv_rep, phisv, alphasv, musv);
         projectors[4].setTseries(isv, getN_obs(), sv_rep); projectors[4].go();       
      }
      else if(x==5) //HEAVY
      {
        stationary[5] = true; n_heavy = 2;
        double[] temp = new double[getN_obs()*2]; 
        heavy_ts = new double[n_heavy*getN_obs()];
        heavy_ts  = simulateHeavy(2, m2, w1, w2, h_alpha, h_alpha_R, h_lambda, h_beta, h_beta_R);

        for(k=0;k<2;k++) 
        {
          for(i=0;i<getN_obs();i++)
          {temp[k*getN_obs() + i] = heavy_ts[2*i + k];}
        }
        projectors[5].setTseries(temp, getN_obs(), 2); projectors[5].go();   
      }
   }



   public void activatePanel(int i, boolean sel)
   {
        if(i==0)
        {
           muBar.setEnabled(sel); alphaBar.setEnabled(sel); 
           skewBar.setEnabled(sel); scaleBar.setEnabled(sel);
            gaussDist.setEnabled(sel);
           levyDist.setEnabled(sel); studDist.setEnabled(sel);
           p.setEnabled(sel); d.setEnabled(sel); q.setEnabled(sel); 
           P.setEnabled(sel); D.setEnabled(sel); Q.setEnabled(sel);
                     
        }
        else if(i==1)
        {
          integrate.setEnabled(sel); 
          b1Bar.setEnabled(sel); b0Bar.setEnabled(sel); a1Bar.setEnabled(sel); thetaBar.setEnabled(sel);
        }
        else if(i==2)
        {
           lagBar.setEnabled(sel); lambdaBar.setEnabled(sel); rhoBar.setEnabled(sel);
        }
        else if(i==3)
        {
           trend1Bar.setEnabled(sel); trend2Bar.setEnabled(sel); 
           if(sim_check[0]) {affectRatio.setEnabled(sel); addLagTrend.setEnabled(sel);}
           
        }
        else if(i==4)
        {
           integratesv.setEnabled(sel);
           alphasvBar.setEnabled(sel); phisvBar.setEnabled(sel); musvBar.setEnabled(sel); nrepBar.setEnabled(sel);
  
        }
        else if(i==5)
        {
          integrateHeavy.setEnabled(sel);
          h_alphaBar.setEnabled(sel); h_alphaRBar.setEnabled(sel); 
          h_betaBar.setEnabled(sel); h_betaRBar.setEnabled(sel);
          h_lambdaBar.setEnabled(sel); m2Bar.setEnabled(sel);
          w1Bar.setEnabled(sel); w2Bar.setEnabled(sel);
        }
        if(!sel) {projectors[i].clear();}
        seedsofChange[i].setEnabled(sel); addMe[i].setEnabled(sel);
        if(sel) {sampleSeries(i);}
   }


   public void updateTime(int i, boolean s)
   {
      getTheatre().setplot(i,s); 
   }



   public void setupPanels()
   {

       GroupLayout paramLayout;
       new JLabel("    "); new JLabel("    "); new JLabel("    ");
       new JLabel("    ");

       JLabel[] seedTitles = new JLabel[6];  JPanel[] seedCon = new JPanel[6];
       for(int i=0;i<6;i++)
       {
         seedTitles[i] = new JLabel(" seed"); seedCon[i] = new JPanel(); //seedCon[i].setPreferredSize(new Dimension(70,20));
         //seedText[i].setPreferredSize(new Dimension(10,10));       
         
         paramLayout = new GroupLayout(seedCon[i]);                
         paramLayout.setAutoCreateGaps(false); paramLayout.setAutoCreateContainerGaps(false);        
         paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()  
           .addComponent(seedTitles[i]).addComponent(seedsofChange[i]));

         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.CENTER)
               .addComponent(seedTitles[i]).addComponent(seedsofChange[i])));
         seedCon[i].setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));
       }
       





     //------------- SET SIM PANELS FIRST ------------------------------------

        /*----------- SET DIALOG BOXES ------------------------*/

        arimaDialog = new JDialog(frame,true);
        garchDialog = new JDialog(frame,true);
        cycleDialog = new JDialog(frame,true);
        trendDialog = new JDialog(frame,true);
        svDialog    = new JDialog(frame,true);
        heavyDialog = new JDialog(frame,true);




       arimaDialogButton = new JButton("Parameters"); arimaDialogButton.setFont(new Font("Ubuntu", 0, 12)); 
       garchDialogButton = new JButton("Parameters"); garchDialogButton.setFont(new Font("Ubuntu", 0, 12)); 
       cycleDialogButton = new JButton("Parameters"); cycleDialogButton.setFont(new Font("Ubuntu", 0, 12)); 
       trendDialogButton = new JButton("Parameters"); trendDialogButton.setFont(new Font("Ubuntu", 0, 12)); 
       svDialogButton = new JButton("Parameters");    svDialogButton.setFont(new Font("Ubuntu", 0, 12)); 
       heavyDialogButton = new JButton("Parameters"); heavyDialogButton.setFont(new Font("Ubuntu", 0, 12)); 


       








       //------------- SARIMA PANEL------------------------------------------
 
 

       sarimaPanel = new JPanel(); 
       sarimaPanel.setBorder(new BevelBorder(BevelBorder.RAISED));
       JPanel simControls = new JPanel();
       JTabbedPane tabSarima = new JTabbedPane((JTabbedPane.TOP));
       //tabSarima.setBorder(new BevelBorder(BevelBorder.RAISED));

        
       //----------- SARIMA SCROLLBARS ------------------------------------
       JPanel modelXPanel = new JPanel();
       JPanel muCon = new JPanel(); 
       JPanel alphaCon = new JPanel(); 
       JPanel skewCon = new JPanel();  
       JPanel scaleCon = new JPanel(); 
       JPanel distCon = new JPanel();   

       muBar.setPreferredSize(new Dimension(150,15));  skewBar.setPreferredSize(new Dimension(150,15));
       alphaBar.setPreferredSize(new Dimension(150,15));  scaleBar.setPreferredSize(new Dimension(150,15));   


       muCon.add(new JLabel("\u03BC")); muCon.add(muBar); muCon.add(muText); muCon.setLayout(new FlowLayout()); 
       alphaCon.add(new JLabel("\u03B1")); alphaCon.add(alphaBar); alphaCon.add(alphaText); alphaCon.setLayout(new FlowLayout());
       skewCon.add(new JLabel("\u03D0")); skewCon.add(skewBar); skewCon.add(skewText);  skewCon.setLayout(new FlowLayout());
       scaleCon.add(new JLabel("c")); scaleCon.add(scaleBar); scaleCon.add(scaleText); scaleCon.setLayout(new FlowLayout());
       distCon.add(gaussDist); distCon.add(studDist); distCon.add(levyDist); distCon.setLayout(new FlowLayout());


        JPanel paramPanel = new JPanel();
        paramLayout = new GroupLayout(paramPanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(muCon).addComponent(alphaCon).addComponent(skewCon));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(muCon).addComponent(alphaCon).addComponent(skewCon)));
        paramPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));

       pLabel = new JLabel("p"); dLabel = new JLabel("d"); qLabel = new JLabel("q");
       PLabel = new JLabel("P"); DLabel = new JLabel("D"); QLabel = new JLabel("Q"); 

       JPanel pX = new JPanel(); pX.setLayout(new FlowLayout()); pX.add(pLabel); pX.add(p); 
       JPanel dX = new JPanel(); dX.setLayout(new FlowLayout()); dX.add(dLabel); dX.add(d);    
       JPanel qX = new JPanel(); qX.setLayout(new FlowLayout()); qX.add(qLabel); qX.add(q); 
       JPanel PX = new JPanel(); PX.setLayout(new FlowLayout()); PX.add(PLabel); PX.add(P);    
       JPanel DX = new JPanel(); DX.setLayout(new FlowLayout()); DX.add(DLabel); DX.add(D); 
       JPanel QX = new JPanel(); QX.setLayout(new FlowLayout()); QX.add(QLabel); QX.add(Q);
   
        paramLayout = new GroupLayout(modelXPanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(addLagTrend).addComponent(pX).addComponent(dX).addComponent(qX).addComponent(PX).addComponent(DX).addComponent(QX));

         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(addLagTrend).addComponent(pX).addComponent(dX).addComponent(qX).addComponent(PX).addComponent(DX).addComponent(QX)));
        modelXPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));
             
        tabSarima.addTab("Model",modelXPanel); 
        tabSarima.addTab("\u03B5_t Dist.",distCon);
        tabSarima.addTab("\u03B5_t Params",paramPanel);


        JPanel topPanel = new JPanel();
        paramLayout = new GroupLayout(topPanel);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[0]).addComponent(addMe[0]).addComponent(arimaDialogButton).addComponent(seedCon[0]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(sim_box[0]).addComponent(addMe[0]).addComponent(arimaDialogButton).addComponent(seedCon[0])));
        topPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));        

        
        paramLayout = new GroupLayout(sarimaPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(topPanel).addComponent(projectors[0])));//.addComponent(tabSarima)));
         
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()          
          .addComponent(topPanel)
          .addComponent(projectors[0]));
          //.addComponent(tabSarima));
         sarimaPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));      



   
       //----------------- GARCH PANEL ------------------------------------

       //---- GARCH --------------------
       garchPanel = new JPanel(); 
       garchPanel.setBorder(new BevelBorder(BevelBorder.RAISED));
  
       JLabel b0Label = new JLabel("b0"); JLabel b1Label = new JLabel("b1"); 
       JLabel a1Label = new JLabel("\u03B1"); JLabel thetaLabel = new JLabel("\u0398");

       JPanel b0Con = new JPanel(); 
       JPanel b1Con = new JPanel();  
       JPanel a1Con = new JPanel(); 
       JPanel thetaCon = new JPanel();   

       b0Bar.setPreferredSize(new Dimension(150,15));  b1Bar.setPreferredSize(new Dimension(150,15));
       a1Bar.setPreferredSize(new Dimension(150,15));  thetaBar.setPreferredSize(new Dimension(150,15));   

       trend1Bar.setPreferredSize(new Dimension(150,15)); 
       trend2Bar.setPreferredSize(new Dimension(150,15));
       lagBar.setPreferredSize(new Dimension(150,15));
       lambdaBar.setPreferredSize(new Dimension(150,15));
       rhoBar.setPreferredSize(new Dimension(150,15));


       b0Con.add(b0Label); b0Con.add(b0Bar); b0Con.add(b0Text); b0Con.setLayout(new FlowLayout()); 
       b1Con.add(b1Label); b1Con.add(b1Bar); b1Con.add(b1Text); b1Con.setLayout(new FlowLayout());
       a1Con.add(a1Label); a1Con.add(a1Bar); a1Con.add(a1Text);  a1Con.setLayout(new FlowLayout());
       thetaCon.add(thetaLabel); thetaCon.add(thetaBar); thetaCon.add(thetaText); thetaCon.setLayout(new FlowLayout());
       nrepText = new JTextField(1);

        /*paramLayout = new GroupLayout(paramPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(muCon).addComponent(alphaCon).addComponent(skewCon).addComponent(scaleCon));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(muCon).addComponent(alphaCon).addComponent(skewCon).addComponent(scaleCon)));
        paramPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));*/

 

       JPanel InnPanel = new JPanel(); 
       paramLayout = new GroupLayout(InnPanel); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(integrate).addComponent(b0Con).addComponent(b1Con));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(integrate).addComponent(b0Con).addComponent(b1Con)));
       InnPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));      

       JPanel InnPanel2 = new JPanel(); 
       paramLayout = new GroupLayout(InnPanel2); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(a1Con).addComponent(thetaCon));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(a1Con).addComponent(thetaCon)));  
       InnPanel2.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));     

       JTabbedPane tabGarch = new JTabbedPane(JTabbedPane.TOP);
       tabGarch.addTab("AR Params",InnPanel);
       tabGarch.addTab("\u03B1^2 Params",InnPanel2);
       //tabGarch.setBorder(new BevelBorder(BevelBorder.RAISED));


        /*JPanel garchParamPanel = new JPanel();
        paramLayout = new GroupLayout(garchParamPanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()          
           .addComponent(b0Con).addComponent(a1Con).addComponent(b1Con).addComponent(thetaCon));    
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(b0Con) .addComponent(b1Con).addComponent(a1Con).addComponent(thetaCon)));
         garchParamPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));      */

        JPanel topPanel2 = new JPanel();
        paramLayout = new GroupLayout(topPanel2);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[1]).addComponent(addMe[1]).addComponent(garchDialogButton).addComponent(seedCon[1]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(sim_box[1]).addComponent(addMe[1]).addComponent(garchDialogButton).addComponent(seedCon[1])));
        topPanel2.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));        

        paramLayout = new GroupLayout(garchPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(topPanel2).addComponent(projectors[1]))); //.addComponent(tabGarch)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
            .addComponent(topPanel2)
            .addComponent(projectors[1]));
            //.addComponent(tabGarch));
        garchPanel.setLayout(paramLayout);











 
       //---------- HEAVY PANEL----------------------------------------------------

       JPanel w1Con = new JPanel(); JPanel w2Con = new JPanel();
       JPanel h_alphaCon = new JPanel(); JPanel h_alphaRCon = new JPanel(); 
       JPanel h_betaCon = new JPanel(); JPanel h_betaRCon = new JPanel(); 
       JPanel h_lambdaCon = new JPanel(); JPanel m2Con = new JPanel(); 

       //------------------------------------- MAKE SLIDERS----------------------------------------


       m2Bar.setPreferredSize(new Dimension(150,15)); w1Bar.setPreferredSize(new Dimension(150,15));
       h_betaBar.setPreferredSize(new Dimension(150,15)); w2Bar.setPreferredSize(new Dimension(150,15));
       h_alphaRBar.setPreferredSize(new Dimension(150,15)); h_alphaBar.setPreferredSize(new Dimension(150,15));
       h_betaRBar.setPreferredSize(new Dimension(150,15)); h_lambdaBar.setPreferredSize(new Dimension(150,15)); 


        paramLayout = new GroupLayout(m2Con); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(m2Label).addComponent(m2Bar).addComponent(m2Text));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(m2Label).addComponent(m2Bar).addComponent(m2Text))); m2Con.setLayout(paramLayout);

        paramLayout = new GroupLayout(w1Con); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(w1Label).addComponent(w1Bar).addComponent(w1Text));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(w1Label).addComponent(w1Bar).addComponent(w1Text))); w1Con.setLayout(paramLayout);

        paramLayout = new GroupLayout(w2Con); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(w2Label).addComponent(w2Bar).addComponent(w2Text));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(w2Label).addComponent(w2Bar).addComponent(w2Text))); w2Con.setLayout(paramLayout);

        paramLayout = new GroupLayout(h_alphaCon); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(h_alphaLabel).addComponent(h_alphaBar).addComponent(h_alphaText));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(h_alphaLabel).addComponent(h_alphaBar).addComponent(h_alphaText))); h_alphaCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(h_alphaRCon); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(h_alphaRLabel).addComponent(h_alphaRBar).addComponent(h_alphaRText));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(h_alphaRLabel).addComponent(h_alphaRBar).addComponent(h_alphaRText))); h_alphaRCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(h_betaCon); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(h_betaLabel).addComponent(h_betaBar).addComponent(h_betaText));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(h_betaLabel).addComponent(h_betaBar).addComponent(h_betaText))); h_betaCon.setLayout(paramLayout);


        paramLayout = new GroupLayout(h_betaRCon); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(h_betaRLabel).addComponent(h_betaRBar).addComponent(h_betaRText));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(h_betaRLabel).addComponent(h_betaRBar).addComponent(h_betaRText))); h_betaRCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(h_lambdaCon); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup().addComponent(h_lambdaLabel).addComponent(h_lambdaBar).addComponent(h_lambdaText));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup().addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(h_lambdaLabel).addComponent(h_lambdaBar).addComponent(h_lambdaText))); h_lambdaCon.setLayout(paramLayout);


       /*
       m2Con.add(m2Label); m2Con.add(m2Bar); m2Con.add(m2Text); m2Con.setLayout(new FlowLayout()); 
       w1Con.add(w1Label); w1Con.add(w1Bar); w1Con.add(w1Text); w1Con.setLayout(new FlowLayout()); 
       w2Con.add(w2Label); w2Con.add(w2Bar); w2Con.add(w2Text); w2Con.setLayout(new FlowLayout()); 
       h_alphaCon.add(h_alphaLabel); h_alphaCon.add(h_alphaBar); h_alphaCon.add(h_alphaText); h_alphaCon.setLayout(new FlowLayout()); 
       h_alphaRCon.add(h_alphaRLabel); h_alphaRCon.add(h_alphaRBar); h_alphaRCon.add(h_alphaRText); h_alphaRCon.setLayout(new FlowLayout()); 
       h_betaCon.add(h_betaLabel); h_betaCon.add(h_betaBar); h_betaCon.add(h_betaText); h_betaCon.setLayout(new FlowLayout()); 
       h_lambdaCon.add(h_betaRLabel); h_betaRCon.add(h_betaRBar); h_betaRCon.add(h_betaRText);  h_betaRCon.setLayout(new FlowLayout()); 
       h_lambdaCon.add(h_lambdaLabel); h_lambdaCon.add(h_lambdaBar); h_lambdaCon.add(h_lambdaText); h_lambdaCon.setLayout(new FlowLayout()); 
       */


       JPanel wPanel = new JPanel(); 
       JPanel RMPanel = new JPanel();
       JPanel betaPanel = new JPanel(); 
       
       
       paramLayout = new GroupLayout(wPanel);                
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(integrateHeavy).addComponent(w1Con).addComponent(w2Con).addComponent(m2Con));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(integrateHeavy).addComponent(w1Con).addComponent(w2Con).addComponent(m2Con)));
        wPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));

   
       paramLayout = new GroupLayout(RMPanel);                
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(h_alphaCon).addComponent(h_alphaRCon).addComponent(h_lambdaCon));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(h_alphaCon).addComponent(h_alphaRCon).addComponent(h_lambdaCon)));
        RMPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));

       
       paramLayout = new GroupLayout(betaPanel);                
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(h_betaCon).addComponent(h_betaRCon));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(h_betaCon).addComponent(h_betaRCon)));
        betaPanel.setLayout(paramLayout); 
       

       /*
       wPanel.add(w1Con); wPanel.add(w2Con); wPanel.add(m2Con); wPanel.setLayout(new FlowLayout());
       RMPanel.add(h_alphaCon); RMPanel.add(h_alphaRCon); 
       RMPanel.add(h_lambdaCon); RMPanel.setLayout(new FlowLayout());
       
       betaPanel.add(h_betaCon); betaPanel.add(h_betaRCon); betaPanel.setLayout(new FlowLayout()); 
       */ 

       JTabbedPane tabHeavy = new JTabbedPane(JTabbedPane.TOP); 
       tabHeavy.addTab("W Params",wPanel);
       tabHeavy.addTab("Alpha Params",RMPanel);
       tabHeavy.addTab("Beta Params", betaPanel);
       //tabHeavy.setBorder(new BevelBorder(BevelBorder.RAISED));
       //tabHeavy.setSize(new Dimension(270,20));

        JPanel h_topPanel = new JPanel();
        paramLayout = new GroupLayout(h_topPanel);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[5]).addComponent(addMe[5]).addComponent(heavyDialogButton).addComponent(seedCon[5]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(sim_box[5]).addComponent(addMe[5]).addComponent(heavyDialogButton).addComponent(seedCon[5])));
        h_topPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));        

        JPanel heavyPanel = new JPanel(); heavyPanel.setBorder(new BevelBorder(BevelBorder.RAISED));
        paramLayout = new GroupLayout(heavyPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(h_topPanel).addComponent(projectors[5])));//.addComponent(tabHeavy)));
         
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()          
          .addComponent(h_topPanel)
          .addComponent(projectors[5]));// 
          //.addComponent(tabHeavy));
         heavyPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));      












       //---- Cycles --------------------
       cyclePanel = new JPanel(); cyclePanel.setBorder(new BevelBorder(BevelBorder.RAISED));
       JLabel lagLabel = new JLabel("lag"); 
       JLabel lambdaLabel = new JLabel("\u03A8"); 
       JLabel rhoLabel = new JLabel("\u03C1");   

       JPanel cycle1Con = new JPanel(); cycle1Con.add(lagLabel); cycle1Con.add(lagBar); cycle1Con.add(lagText); cycle1Con.setLayout(new FlowLayout());
       JPanel cycle2Con = new JPanel(); cycle2Con.add(lambdaLabel); cycle2Con.add(lambdaBar); cycle2Con.add(lambdaText); cycle2Con.setLayout(new FlowLayout()); 
       JPanel cycle3Con = new JPanel(); cycle3Con.add(rhoLabel); cycle3Con.add(rhoBar); cycle3Con.add(rhoText); cycle3Con.setLayout(new FlowLayout());

       JPanel InnPanel3 = new JPanel(); 
       paramLayout = new GroupLayout(InnPanel3); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(cycle1Con).addComponent(cycle2Con).addComponent(cycle3Con));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(cycle1Con).addComponent(cycle2Con).addComponent(cycle3Con)));     
       InnPanel3.setLayout(paramLayout);

 
       /*paramLayout = new GroupLayout(InnPanel3); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(lagLabel).addComponent(lagBar).addComponent(lagText).addComponent(space3)
           .addComponent(lambdaLabel).addComponent(lambdaBar).addComponent(lambdaText).addComponent(space3)
           .addComponent(rhoLabel).addComponent(rhoBar).addComponent(rhoText).addComponent(space3));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(lagLabel).addComponent(lagBar).addComponent(lagText).addComponent(space3)
           .addComponent(lambdaLabel).addComponent(lambdaBar).addComponent(lambdaText).addComponent(space3)
           .addComponent(rhoLabel).addComponent(rhoBar).addComponent(rhoText).addComponent(space3)));
       InnPanel3.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));*/   

        JPanel topPanel3 = new JPanel();
        paramLayout = new GroupLayout(topPanel3);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[2]).addComponent(addMe[2]).addComponent(cycleDialogButton).addComponent(seedCon[2]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(sim_box[2]).addComponent(addMe[2]).addComponent(cycleDialogButton).addComponent(seedCon[2])));
        topPanel3.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));     


        paramLayout = new GroupLayout(cyclePanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()                  
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(topPanel3).addComponent(projectors[2]))); //.addComponent(InnPanel3)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
            .addComponent(topPanel3).addComponent(projectors[2])); // .addComponent(InnPanel3));
        cyclePanel.setLayout(paramLayout);


       //---- Trend --------------------
       trendPanel = new JPanel(); trendPanel.setBorder(new BevelBorder(BevelBorder.RAISED));

       JLabel trend1Label = new JLabel("\u0398_1");  JLabel trend2Label = new JLabel("\u0398_2");
       JPanel trend1Con = new JPanel(); 
       trend1Con.add(trend1Label); trend1Con.add(trend1Bar); trend1Con.add(trend1Text); trend1Con.setLayout(new FlowLayout());
       JPanel trend2Con = new JPanel(); 
       trend2Con.add(trend2Label); trend2Con.add(trend2Bar); trend2Con.add(trend2Text); trend2Con.setLayout(new FlowLayout()); 


       JPanel trendPanel2 = new JPanel(); 
       paramLayout = new GroupLayout(trendPanel2); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(trend1Con).addComponent(trend2Con));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(trend1Con).addComponent(trend2Con)));     
       trendPanel2.setLayout(paramLayout);

 
        JPanel topPanel4 = new JPanel();
        paramLayout = new GroupLayout(topPanel4);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[3]).addComponent(addMe[3]).addComponent(trendDialogButton).addComponent(seedCon[3]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(sim_box[3]).addComponent(addMe[3]).addComponent(trendDialogButton).addComponent(seedCon[3])));
        topPanel4.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));        

      

        paramLayout = new GroupLayout(trendPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()                
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(topPanel4).addComponent(projectors[3]))); //.addComponent(trendPanel2)));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addComponent(topPanel4).addComponent(projectors[3])); //.addComponent(trendPanel2));
        trendPanel.setLayout(paramLayout);




       //----- Add SV Box and Blank box ----------------------

       //---- SV --------------------
       JPanel svPanel = new JPanel(); svPanel.setBorder(new BevelBorder(BevelBorder.RAISED));
       JLabel nrepLabel = new JLabel("NSeries"); 
   
       phisvBar.setPreferredSize(new Dimension(150,15)); 
       musvBar.setPreferredSize(new Dimension(150,15));
       alphasvBar.setPreferredSize(new Dimension(150,15)); 


       JPanel sv0Con = new JPanel(); sv0Con.add(nrepLabel); sv0Con.add(nrepBar); sv0Con.setLayout(new FlowLayout());
       JPanel sv1Con = new JPanel(); sv1Con.add(phisvLabel); sv1Con.add(phisvBar); sv1Con.add(phisvText); sv1Con.setLayout(new FlowLayout());
       JPanel sv2Con = new JPanel(); sv2Con.add(alphasvLabel); sv2Con.add(alphasvBar); sv2Con.add(alphasvText); sv2Con.setLayout(new FlowLayout()); 
       JPanel sv3Con = new JPanel(); sv3Con.add(musvLabel); sv3Con.add(musvBar); sv3Con.add(musvText); sv3Con.setLayout(new FlowLayout());

       JPanel InnPanel4 = new JPanel(); 
       paramLayout = new GroupLayout(InnPanel4); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(integratesv).addComponent(sv0Con).addComponent(sv1Con).addComponent(sv2Con).addComponent(sv3Con));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(integratesv).addComponent(sv0Con).addComponent(sv1Con).addComponent(sv2Con).addComponent(sv3Con)));     
       InnPanel4.setLayout(paramLayout);

 
       /*paramLayout = new GroupLayout(InnPanel3); paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(lagLabel).addComponent(lagBar).addComponent(lagText).addComponent(space3)
           .addComponent(lambdaLabel).addComponent(lambdaBar).addComponent(lambdaText).addComponent(space3)
           .addComponent(rhoLabel).addComponent(rhoBar).addComponent(rhoText).addComponent(space3));        
       paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(lagLabel).addComponent(lagBar).addComponent(lagText).addComponent(space3)
           .addComponent(lambdaLabel).addComponent(lambdaBar).addComponent(lambdaText).addComponent(space3)
           .addComponent(rhoLabel).addComponent(rhoBar).addComponent(rhoText).addComponent(space3)));
       InnPanel3.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));*/   

        JPanel topPanel6 = new JPanel();
        paramLayout = new GroupLayout(topPanel6);                
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(sim_box[4]).addComponent(addMe[4]).addComponent(svDialogButton).addComponent(seedCon[4]));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
               .addComponent(sim_box[4]).addComponent(addMe[4]).addComponent(svDialogButton).addComponent(seedCon[4])));
        topPanel6.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));     


        paramLayout = new GroupLayout(svPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()                  
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(topPanel6).addComponent(projectors[4]))); //.addComponent(InnPanel4)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
            .addComponent(topPanel6).addComponent(projectors[4])); //.addComponent(InnPanel4));
        svPanel.setLayout(paramLayout);







        ///----------------- SET ACTIONS OF DIALOGS-----------------------------------------

        arimaDialog.getContentPane().add(tabSarima); arimaDialog.pack(); arimaDialog.setLocationRelativeTo(frame); 
        arimaDialog.setModal(false); arimaDialog.setVisible(false);

        garchDialog.getContentPane().add(tabGarch); garchDialog.pack(); garchDialog.setLocationRelativeTo(frame); 
        garchDialog.setModal(false); garchDialog.setVisible(false);

        cycleDialog.getContentPane().add(InnPanel3); cycleDialog.pack(); cycleDialog.setLocationRelativeTo(frame); 
        cycleDialog.setModal(false); cycleDialog.setVisible(false);

        trendDialog.getContentPane().add(trendPanel2); trendDialog.pack(); trendDialog.setLocationRelativeTo(frame); 
        trendDialog.setModal(false); trendDialog.setVisible(false);

        svDialog.getContentPane().add(InnPanel4); svDialog.pack(); svDialog.setLocationRelativeTo(frame); 
        svDialog.setModal(false); svDialog.setVisible(false);

        heavyDialog.getContentPane().add(tabHeavy);  heavyDialog.pack(); heavyDialog.setLocationRelativeTo(frame); 
        heavyDialog.setModal(false); heavyDialog.setVisible(false);


        arimaDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            arimaDialog.setModal(false); arimaDialog.setVisible(true);;
            }
        });

       garchDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            garchDialog.setModal(false); garchDialog.setVisible(true);;
            }
        });

       cycleDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            cycleDialog.setModal(false); cycleDialog.setVisible(true);;
            }
        });

       trendDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            trendDialog.setModal(false); trendDialog.setVisible(true);;
            }
        });

       svDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            svDialog.setModal(false); svDialog.setVisible(true);;
            }
        });

       heavyDialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            heavyDialog.setModal(false); heavyDialog.setVisible(true);;
            }
        });







 
       JLabel space5 = new JLabel("  ");
       BevelBorder timeSelBorder = new BevelBorder(BevelBorder.RAISED);
       JPanel timeGrid = new JPanel(); timeGrid.setSize(new Dimension(750, 25));      
       paramLayout = new GroupLayout(timeGrid);
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(acfplotBox).addComponent(resizeBox).addComponent(logBox).addComponent(deleteAll).addComponent(space5).addComponent(plotTarget).addComponent(timePlot[0])
          .addComponent(timePlot[1]).addComponent(timePlot[2])
          .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6])
          .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[10]).addComponent(timePlot[11]));
       paramLayout.setVerticalGroup(
       paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(acfplotBox).addComponent(resizeBox).addComponent(logBox).addComponent(deleteAll).addComponent(space5).addComponent(plotTarget)
           .addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
           .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6]) 
           .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[10]).addComponent(timePlot[11])));    
       timeGrid.setLayout(paramLayout); 
       timeGrid.setBorder(timeSelBorder);  

      Box timePane = Box.createVerticalBox();
      timePane.add(timeScrollPane,BorderLayout.NORTH);
      timePane.add(timeGrid,BorderLayout.SOUTH);

      JPanel sharpeCon = new JPanel(); JPanel meanCon = new JPanel(); JPanel maxDrawCon = new JPanel();
      JPanel totalRetCon = new JPanel(); JPanel totalRankCon = new JPanel(); JPanel minRankCon = new JPanel(); 
      JPanel avgRankCon = new JPanel(); JPanel funcCon = new JPanel();
      
    sharpeCon.add(sharpeLabel); sharpeCon.add(sharpeText); sharpeCon.setLayout(new FlowLayout()); 
    meanCon.add(meanLabel); meanCon.add(meanText); meanCon.setLayout(new FlowLayout());
    maxDrawCon.add(maxDrawLabel); maxDrawCon.add(maxDrawText); maxDrawCon.setLayout(new FlowLayout());
    totalRankCon.add(totalRankLabel); totalRankCon.add(totalRankText); totalRankCon.setLayout(new FlowLayout());
    minRankCon.add(minRankLabel); minRankCon.add(minRankText); minRankCon.setLayout(new FlowLayout());
    avgRankCon.add(avgRankLabel); avgRankCon.add(avgRankText); avgRankCon.setLayout(new FlowLayout());
    totalRetCon.add(totalRetLabel); totalRetCon.add(totalRetText); totalRetCon.setLayout(new FlowLayout());       
    funcCon.add(normalizeWeightsSeries); funcCon.add(maxSharpeButton); funcCon.setLayout(new FlowLayout());  
      
      

      JLabel serLabel =  new JLabel("   Returns");
      JLabel serLabel1 = new JLabel("   Weight in Target");
      JLabel serLabel2 = new JLabel(" Computed in Target");
      JLabel serLabel3 = new JLabel("  Explaining Target");
      JLabel serLabel4 = new JLabel("       Price Series");
 
      JPanel thirdPanel = new JPanel();
      JPanel targetPanel = new JPanel(); 
      BevelBorder targetBorder = new BevelBorder(BevelBorder.RAISED);
      targetPanel.setBorder(targetBorder); targetPanel.setSize(new Dimension(600,300));

       paramLayout = new GroupLayout(targetPanel);
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(serLabel).addComponent(serLabel1).addComponent(serLabel2).addComponent(serLabel3).addComponent(serLabel4))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[0]).addComponent(weightSlidersText[0]).addComponent(weightSliders[0]).addComponent(simCheck[0]).addComponent(repCheck[0]).addComponent(priceCheck[0]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[1]).addComponent(weightSlidersText[1]).addComponent(weightSliders[1]).addComponent(simCheck[1]).addComponent(repCheck[1]).addComponent(priceCheck[1]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[2]).addComponent(weightSlidersText[2]).addComponent(weightSliders[2]).addComponent(simCheck[2]).addComponent(repCheck[2]).addComponent(priceCheck[2]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[3]).addComponent(weightSlidersText[3]).addComponent(weightSliders[3]).addComponent(simCheck[3]).addComponent(repCheck[3]).addComponent(priceCheck[3]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
           .addComponent(weightLabel[4]).addComponent(weightSlidersText[4]).addComponent(weightSliders[4]).addComponent(simCheck[4]).addComponent(repCheck[4]).addComponent(priceCheck[4]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[5]).addComponent(weightSlidersText[5]).addComponent(weightSliders[5]).addComponent(simCheck[5]).addComponent(repCheck[5]).addComponent(priceCheck[5]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[6]).addComponent(weightSlidersText[6]).addComponent(weightSliders[6]).addComponent(simCheck[6]).addComponent(repCheck[6]).addComponent(priceCheck[6]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[7]).addComponent(weightSlidersText[7]).addComponent(weightSliders[7]).addComponent(simCheck[7]).addComponent(repCheck[7]).addComponent(priceCheck[7]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[8]).addComponent(weightSlidersText[8]).addComponent(weightSliders[8]).addComponent(simCheck[8]).addComponent(repCheck[8]).addComponent(priceCheck[8]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[9]).addComponent(weightSlidersText[9]).addComponent(weightSliders[9]).addComponent(simCheck[9]).addComponent(repCheck[9]).addComponent(priceCheck[9]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[10]).addComponent(weightSlidersText[10]).addComponent(weightSliders[10]).addComponent(simCheck[10]).addComponent(repCheck[10]).addComponent(priceCheck[10]))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
            .addComponent(weightLabel[11]).addComponent(weightSlidersText[11]).addComponent(weightSliders[11]).addComponent(simCheck[11]).addComponent(repCheck[11]).addComponent(priceCheck[11]))            
      );     


        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(serLabel).addComponent(weightLabel[0]).addComponent(weightLabel[1]).addComponent(weightLabel[2])
            .addComponent(weightLabel[3]).addComponent(weightLabel[4]).addComponent(weightLabel[5])
            .addComponent(weightLabel[6]).addComponent(weightLabel[7]).addComponent(weightLabel[8])
            .addComponent(weightLabel[9]).addComponent(weightLabel[10]).addComponent(weightLabel[11]))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(serLabel1).addComponent(weightSlidersText[0]).addComponent(weightSlidersText[1]).addComponent(weightSlidersText[2])
            .addComponent(weightSlidersText[3]).addComponent(weightSlidersText[4]).addComponent(weightSlidersText[5])
            .addComponent(weightSlidersText[6]).addComponent(weightSlidersText[7]).addComponent(weightSlidersText[8])
            .addComponent(weightSlidersText[9]).addComponent(weightSlidersText[10]).addComponent(weightSlidersText[11]))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(weightSliders[0]).addComponent(weightSliders[1]).addComponent(weightSliders[2])
            .addComponent(weightSliders[3]).addComponent(weightSliders[4]).addComponent(weightSliders[5])
            .addComponent(weightSliders[6]).addComponent(weightSliders[7]).addComponent(weightSliders[8])
            .addComponent(weightSliders[9]).addComponent(weightSliders[10]).addComponent(weightSliders[11]))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(serLabel2).addComponent(simCheck[0]).addComponent(simCheck[1]).addComponent(simCheck[2])
            .addComponent(simCheck[3]).addComponent(simCheck[4]).addComponent(simCheck[5])
            .addComponent(simCheck[6]).addComponent(simCheck[7]).addComponent(simCheck[8])
            .addComponent(simCheck[9]).addComponent(simCheck[10]).addComponent(simCheck[11]))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(serLabel3).addComponent(repCheck[0]).addComponent(repCheck[1]).addComponent(repCheck[2])
            .addComponent(repCheck[3]).addComponent(repCheck[4]).addComponent(repCheck[5])
            .addComponent(repCheck[6]).addComponent(repCheck[7]).addComponent(repCheck[8])
            .addComponent(repCheck[9]).addComponent(repCheck[10]).addComponent(repCheck[11]))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(serLabel4).addComponent(priceCheck[0]).addComponent(priceCheck[1]).addComponent(priceCheck[2])
            .addComponent(priceCheck[3]).addComponent(priceCheck[4]).addComponent(priceCheck[5])
            .addComponent(priceCheck[6]).addComponent(priceCheck[7]).addComponent(priceCheck[8])
            .addComponent(priceCheck[9]).addComponent(priceCheck[10]).addComponent(priceCheck[11])));
        targetPanel.setLayout(paramLayout);   
     
        thirdPanel.add(targetPanel,BorderLayout.CENTER);

        JPanel spreadPanel = new JPanel();
        paramLayout = new GroupLayout(spreadPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(funcCon).addComponent(sharpeCon).addComponent(totalRetCon).addComponent(maxDrawCon).addComponent(totalRankCon).addComponent(avgRankCon)));      
        paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
           .addComponent(funcCon).addComponent(sharpeCon).addComponent(totalRetCon).addComponent(maxDrawCon).addComponent(totalRankCon).addComponent(avgRankCon));
        spreadPanel.setLayout(paramLayout);
        
        thirdPanel.add(spreadPanel,BorderLayout.WEST);


        JTabbedPane models = new JTabbedPane(JTabbedPane.TOP);      
        paramLayout = new GroupLayout(simControls);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(sarimaPanel).addComponent(garchPanel));   
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
              .addComponent(sarimaPanel).addComponent(garchPanel)));
        simControls.setLayout(paramLayout);

        JPanel simControls2 = new JPanel();
        paramLayout = new GroupLayout(simControls2);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(cyclePanel).addComponent(trendPanel));   
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
              .addComponent(cyclePanel).addComponent(trendPanel)));
        simControls2.setLayout(paramLayout);

        JPanel simControls3 = new JPanel();
        paramLayout = new GroupLayout(simControls3);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(svPanel).addComponent(heavyPanel));   
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
              .addComponent(svPanel).addComponent(heavyPanel)));
        simControls3.setLayout(paramLayout);


        models.addTab("ARMA/GARCH", simControls);
        models.addTab("Cycle/Trend",simControls2);
        models.addTab("SV/Heavy", simControls3);
        models.addTab("Target Series",thirdPanel);
   
        paramLayout = new GroupLayout(this);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()           
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(timePane)
             .addComponent(models)));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
           .addComponent(timePane)
           .addComponent(models));
        this.setLayout(paramLayout);

        for(int i=0;i<6;i++) activatePanel(i, false);

       

   }


   class MyActionListener implements ActionListener {
     public void actionPerformed(ActionEvent e)
     {
        for(int i=0;i<6;i++)
        {
         if(e.getSource() == addMe[i])
         {addSeries(i);}
        } 
        if(e.getSource() == integrate)
        {integrateGarch();}
        if(e.getSource() == addLagTrend)
        {addTrend();}
        if(e.getSource() == deleteAll)
        {
          clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); 
          if(returnsAll.size() > 0) {returnsAll.clear();}
          first_in = true;
        }
        if(e.getSource() == integratesv)
        {integrateSV();}

        if(e.getSource() == integrateHeavy)
        {integrateHEAVY();}
        
        if(e.getSource() == nrepBar)
        {
          sv_rep = nrepBar.getSelectedIndex()+1; 
          setFSVParams(phisv, alphasv, musv, sv_rep);
        }
     }
   }
 
   class MySliderListener implements ChangeListener {
      public void stateChanged(ChangeEvent e) 
      {         
        int v; JSlider s = (JSlider)e.getSource();
        for(int i=0;i<n_rep;i++)
        {
          if(s == weightSliders[i])
          {
            v = s.getValue(); 
            weight_mix[i] = v*.01; 
            weightSlidersText[i].setText(""+df.format(weight_mix[i]));
            if(sigex_returns_on)
            {computePortfolio();}
            computeTarget();
          }
        }
      }
   }

   class MyItemListener implements ItemListener {
      public void itemStateChanged(ItemEvent e)
      {         
         int i; boolean sel; Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
 
         for(i=0;i<6;i++)
         {
           if(source == sim_box[i]) 
           {
             sim_check[i] = sel; activatePanel(i,sel); 
           }
         }

         for(i=0;i<getN_sim_series();i++) 
         {
           if(source == timePlot[i])
           {
             if(acfplot && sel)
             {
               if(acfs == 0)
               {setACFCross(i+1,i+1); acf_plot_no1 = i+1; acf_plot_no2 = i+1; acfs = 1; getTheatre().plot_acf_series(true,false,false,false);} 
               else if(acfs == 1)
               {setACFCross(acf_plot_no2,i+1); acf_plot_no1 = acf_plot_no2; acf_plot_no2 = i+1; acfs = 2; getTheatre().plot_acf_series(false,false,true,true);}                  
               else {timePlot[i].setSelected(false);}
             }
             else if(acfplot && !sel)
             {
              if(acfs == 1)
              {acf_plot_no1 = 0; acf_plot_no2 = 0; acfs = 0; getTheatre().plot_acf_series(false,false,false,false);} 
              else if(acfs == 2)
              {
                if(acf_plot_no1 == i+1) {acf_plot_no1 = acf_plot_no2;}
                setACFCross(acf_plot_no1,acf_plot_no1); acf_plot_no2 = acf_plot_no1; acfs = 1; 
                getTheatre().plot_acf_series(true,false,false,false);
              }                  
             }
             else
             {
               updateTime(i,sel);;
             }             
           }
           else if(source == repCheck[i])
           {getRepresent()[i] = sel;}
           else if(source == simCheck[i])
           {mix[i] = sel; computeTarget();}
           else if(source == priceCheck[i])
           {computePriceSeries();}
           
         }
         if(source == plotTarget)
         {
           if(acfplot && sel)
           {
              if(acfs == 0)
              {setACFCross(0,0); acf_plot_no1 = 0; acf_plot_no2 = 0; acfs = 1; getTheatre().plot_acf_series(true,false,false,false);} 
              else if(acfs == 1)
              {setACFCross(acf_plot_no2,0); acf_plot_no1 = acf_plot_no2; acf_plot_no2 = 0; acfs = 2; getTheatre().plot_acf_series(false,false,true,true);}                  
           }
           else if(acfplot && !sel)
           {
              if(acfs == 1)
              {acf_plot_no1 = 0; acf_plot_no2 = 0; acfs = 0; getTheatre().plot_acf_series(false,false,false,false);} 
              else if(acfs == 2)
              {
                if(acf_plot_no1 == 0) {acf_plot_no1 = acf_plot_no2;}
                setACFCross(acf_plot_no1,acf_plot_no1); acf_plot_no2 = acf_plot_no1; acfs = 1; 
                getTheatre().plot_acf_series(true,false,false,false);
              }                  
           }
           else
           {
            getTheatre().plotTarget(sel);
           }

         }
         else if(source == gaussDist) {rand = 0; sampleSeries(0);}
         else if(source == studDist) {rand = 1; sampleSeries(0);}
         else if(source == levyDist) {rand = 2; sampleSeries(0);}
         else if(source == resizeBox) {resize =  sel; standardizeSeries();}
         else if(source == logBox)
         {
           for(int j=0;j<getSim_data().size();j++)
           {if(timePlot[j].isSelected()) {l_log[j] = sel;}}
           computeTarget();
           //if(sel){applyLogTransform();} 
           //else{applyExpTransform();}
         }
         else if(source == acfplotBox)
         {
              acfplot = sel;  
              if(plotTarget.isSelected())
              {plotTarget.setSelected(false);}
              for(i=0;i<getN_sim_series();i++)
              {
                if(timePlot[i].isSelected())
                {timePlot[i].setSelected(false);}
              }
              getTheatre().plotACFCross(sel);
              if(acfplot) {acfs = 0; acf_plot_no1 = 0; acf_plot_no2 = 0;}                           
         }
        

      }
   }  
         
   class MyChangeListener implements ChangeListener {
     public void stateChanged(ChangeEvent e)
     {
       for(int i=0;i<6;i++) 
       {
         if(seedsofChange[i] == e.getSource())
         {
           seed[i] = seedsofChange[i].getValue(); sampleSeries(i); seedText[i].setText("" + df.format(seed[i]));
         }
       }
       if(affectRatio == e.getSource())
       {
         affect = affectRatio.getValue()*.01; affectText.setText(""+df.format(affect));
         sampleSeries(0); addTrend(); 
       }
 
     }
   }
   /*----------------------------------------------------------------------
     Simulates sarima models 
   -----------------------------------------------------------------------*/

   public double[] simulateAR(double ar, double var)
   {
       int[] dim = {1,0,0,0,0,0};
       int n_params = 2; double[] params = new double[n_params];
       params[0] = ar; params[1] = var;

       return sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
     
   }

   public double[] simulateMA(double ma, double var)
   {
       int[] dim = {0,0,1,0,0,0};
       int n_params = 2; double[] params = new double[n_params];
       params[0] = ma; params[1] = var;

       return sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
     
   }

   public double[] simulateAirline(double ma, double MA, double var)
   {
       int[] dim = {0,1,1,0,1,1};
       int n_params = 3; double[] params = new double[n_params];
       params[0] = ma; params[1] = MA; params[2] = var;

       return sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
     
   }

   public double[] simulateOutliers(double ar, double var)
   {
       int[] dim = {1,0,0,0,0,0}; int _rand = 1; 
       int n_params = 2; double[] params = new double[n_params];
       params[0] = ar; params[1] = var; 

       return sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], _rand);
     
   }


   public double[] simulateTrend(double tr1, double tr2, double var)
   {
       int[] _dim = {2,0,2,0,0,0};
       int _n_params = 5; double[] _params = new double[_n_params];
       _params[0] = -2.0; _params[1] = 1.0; _params[2] = tr1; _params[3] = tr2; _params[4] = var;
 
       return sarima(getN_obs(), burnin, _params, _dim, _n_params, S, seed[3], 0);
   }

   /*----------------------------------------------------------------------
     Simulates a trend + SMA (ma,ar) from a given trend that is lagged

     double[] trend - the basic leading trend/cycle indicator
     double ma - the coefficient for the model 
     double var - the variance 
     int lag - how much the new trend is lagged
     int affect - 0 < affect < 1, the affect of new trend

   -----------------------------------------------------------------------*/

   public double[] simulateLaggedTrendAR(double[] trend, double ar, double var, int lag, double affect)
   {
       
       int i; int[] dim = {1,0,0,0,0,0}; int n_params = 2; double[] params = new double[n_params];
       params[0] = ar; params[1] = var;

       double[] y = sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);
      
       for(i=0;i<lag;i++) {y[i] = y[i] + affect*trend[0];}     
       for(i=0;i<getN_obs()-lag;i++)  {y[i+lag] = y[i+lag] + affect*trend[i];}
       return y;
   }

   public double[] simulateLaggedTrendMA(double[] trend, double ma, double var, int lag, double affect)
   {

       int i; int[] dim = {0,0,1,0,0,0}; int n_params = 2; double[] params = new double[n_params];
       params[0] = ma; params[1] = var;

       double[] y = sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand);   

       for(i=0;i<lag;i++) {y[i] = y[i] + affect*trend[0];}     
       for(i=0;i<getN_obs()-lag;i++)  {y[i+lag] = y[i+lag] + affect*trend[i];} 
       return y;
   }

   public double[] simulateLaggedTrendSMA(double[] trend, double ma, double var, int lag, double affect)
   {
              
       int i; int[] dim = {0,0,0,0,0,1}; int n_params = 2; double[] params = new double[n_params];
       params[0] = ma; params[1] = var;
       double[] y = sarima(getN_obs(), burnin, params, dim, n_params, S, seed[0], rand); 

       for(i=0;i<lag;i++) {y[i] = y[i] + affect*trend[0];}     
       for(i=0;i<getN_obs()-lag;i++)  {y[i+lag] = y[i+lag] + affect*trend[i];} 
       return y;
   }


   /*-------------------------------------------------------------------------
     Returns a vector of length 2*n_nobs, two lagged series of cycles
   --------------------------------------------------------------------------*/

   public double[] simulateCycles(int n_c, double lag, double lambda, double rho)
   {
     return cycles(n_c, getN_obs(), burnin, lag, lambda, rho, seed[2]);
   }

   /*-------------------------------------------------------------------------------------
   *  Simulates a NGarch(1,1) model (p-autoregressive, q-moving average)
   *
   *   e_t = s_t * z_t     z_t ~ N(0,1),   s^2_t = b0 + a_1(e_t-1 - theta s_t-1)^2 + b_1 s^2_t-1
   * 
   *   - (b_1 + theta^2) coefficient on ar volatility 
   *   - a_1 coeff on error process 
   *   - 2 theta a_1  coefficient on correlation s_t-1 e_t-1 
   *   - b0 > 0 mean reverting value 
   *
   *   params = b_1 a_1 theta b0     
   * --------------------------------------------------------------------------------------*/
   public double[] simulateNGARCH(double b1, double a1, double theta, double b0)
   { 
          double[] params = {b1, a1, theta, b0};        
          return nGARCH(getN_obs(), burnin, params, seed[1]);
   }

   /*-----------------------------------------------------------------------------------------
     GARCH model, params = [b1..bp, a0, a1..aq]  , size p+1+q     

   ------------------------------------------------------------------------------------------*/
   public double[] simulateGARCH(int p, int q, double[] params)
   { 
          return GARCH(getN_obs(), burnin, p, q, params, seed[1]);
   }

   public double[] simulateFSV(int n_rep, double phi, double sigma, double mu)
   {       
          return fsv(getN_obs(), n_rep, phi, sigma, mu, f1, f2, seed[4]);
   }

   public double[] simulateHeavy(int n_rep, int m2, double w1, double w2, double alpha, double alpha_R, 
                                 double lambda, double beta, double beta_R)
   {

     heavy.setNobs(getN_obs(), n_rep);
     heavy.changeSeed(seed[5]);
     heavy.setParameterValues(w1, w2, alpha, alpha_R, lambda, beta, beta_R); 

     //--- simulate the heavy model with 76? for chisquared process---------
     return heavy.simulateHeavy(m2); 
   }


    public native double[] crossCorrelation(int N, double[] v1, double[] v2);

    public native double[] fsv(int N, int nrep, double phi, double sigma, double mu, double f1, double f2, int seed);

    public native double[] isv(int N, double phi, double sigma, double mu, int seed);

    public native double[] sarima(int N, int burn, double[] params, int[] dim, int n_params, int S, int seed, int rand);

    public native double[] GARCH(int N, int burn, int p, int q, double[] params, int seed); 

    public native double[] nGARCH(int N, int burn, double[] params, int seed);

    public native double[] cycles(int n_cyc, int N, int burn, double lag, double lambda, double rho, int seed);

    public static native void plotData(double[] stat, int dataSize);
  
    //public native void sortsims(int n, double[] s, int[] in);
  
    public static native void plotData2(double[] stat1, double[] stat2, int dataSize);

    static {System.loadLibrary("simSarima");}
    //static {System.loadLibrary("heavy");}
    static {System.loadLibrary("sv");}
    static {System.loadLibrary("ms");}     

    double[] cumsum(double[] data, int n)
    {

      double[] cs = new double[n]; double sum; int k;
      //sum=Math.abs(data[0]); double min = 1000000;
      sum=0; double min = 1000000;

      for(k=0;k<n;k++)
      {
        sum = sum+data[k]; cs[k] = sum; 
        if(cs[k] < min) {min = cs[k];}
      }

      //System.out.println("mincumsum = " + min);
//       if(min < 0.0)
//       {
//          for(k=0;k<n;k++) {cs[k] = cs[k] - (min-5.0);}
//       }
      return cs;  
    } 

 public void sort_sims(int n, double[] object, int[] indices)
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

    
  public void normalizeWeights()
  {
    int n = returnsAll.size();
    double sum = 0; 
    
    for(int i = 0; i < n; i++)
    {sum = sum + weight_mix[i];}
    
    //now set new weight and new slider value 
    for(int i = 0; i < n; i++)
    {weightSliders[i].setValue((int)(100*(weight_mix[i]/sum)));}       
 
  }
  
    

    public void computePortfolio()
    {
    
      int i,j;
      int n_sections; 
      double rank_coeff;
      avg_rank = 0;
      double[] rets = new double[getN_obs()];
     
      int pdays = 0;
     
      cum_port_returns = new double[getN_obs()];
      n_basket = returnsAll.size();      
      port_returns = new double[getN_obs()];
      
      for(i=0;i<n_basket;i++)
      {
        double[] ret = returnsAll.get(i);
        
        for(j=0;j<getN_obs();j++)
        {
         port_returns[j] = port_returns[j] + weight_mix[i]*ret[j];
        }
      }  
      
      for(i=0;i<getN_obs();i++) 
      {
       if(port_returns[i] > 0) {pdays++;}
      }
      total_return = 1.0*pdays/getN_obs();
      
      
//       sharpe_ratio, port_mean, max_drawdown, total_return, total_rank_port, avg_rank, min_rank
       
      //compute statistics of the portfolio  
      double[] mstd = mean_std(port_returns); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];
      
      port_mean = mstd[0]; 

      cum_port_returns = cumsum(port_returns,getN_obs()); 
      
      max_drawdown = computeDrawdown(cum_port_returns);
      //total_return = cum_port_returns[n_obs-1];
      
      System.arraycopy(cum_port_returns,0,rets,0,cum_port_returns.length);
      total_rank_port = rankCoefficient(rets,cum_port_returns.length);
      
      n_sections = (int)(getN_obs()/length_rank);
      
      min_rank = 1.0;
      for(i=0;i<n_sections;i++)
      {
        double[] section = new double[length_rank];
        System.arraycopy(cum_port_returns,i*length_rank,section,0,length_rank);
      
        rank_coeff = rankCoefficient(section,length_rank);
        avg_rank = avg_rank + rank_coeff;
        
        if(rank_coeff < min_rank)
        {min_rank = rank_coeff;}       
      }
      avg_rank = avg_rank/n_sections;
      
      double bRatio = blakelyRatio(15,cum_port_returns);
      
      //set new values 
      sharpeText.setText(""+df.format(sharpe_ratio));
      meanText.setText(""+df3.format(port_mean));
      maxDrawText.setText(""+df3.format(max_drawdown));
      totalRetText.setText(""+df.format(total_return));
      totalRankText.setText(""+df3.format(total_rank_port));
      minRankText.setText(""+df3.format(min_rank));
      avgRankText.setText(""+df3.format(bRatio));
    
  }
    
  public double blakelyRatio(int n_days, double[] rets)
  {
      int i; int n = rets.length; 
      double[] ranks = new double[n - n_days]; 
      double ratio = 0;
      double mean_rank = 0;
      for(i=0;i<n-n_days;i++)
      {
        double[] section = new double[n_days];
        System.arraycopy(rets,i,section,0,n_days);
      
        ranks[i] = rankCoefficient(section,n_days);
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
   
   
    
    
  public double rankCoefficient(double[] account, int n)
  {

   int i; 
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double sum = 0.0;
   double spear;
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i;} 
      
   sort_sims(n, account, index);
   
   for (i=0 ; i<n ; ++i)
   {d[i] = Math.abs(index[i] - rank[i]); d[i] = d[i]*d[i]; sum = sum + 1.0*d[i];} 
   
   spear = 1.0 - (6.0*sum)/(1.0*n*(n*n-1.0));
      
   return spear;   
  }  
    
  public double computeDrawdown(double[] ret)
  {
     int i;
     double max = -100; 
     double[] cmax = cummax(ret);
     cmaxx = new double[ret.length];
     
     for(i=0;i<ret.length;i++)
     {
      cmaxx[i] = cmax[i] - ret[i]; 
      if(cmaxx[i] > max) {max = cmaxx[i];}
     }
     
     return max; 
  }
  
  public static double[] cummax(double[] ret)
  {
    double[] cummax = new double[ret.length];
    
    
    int n = ret.length;
    double max = ret[0];
    
    for(int i = 1; i < n; i++)
    {
      if(ret[i] > max) {max = ret[i];}
      cummax[i] = max;
    }
    return cummax; 
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

    public static double[] mixSeries(final double[] ts, double[] v, double weight)
    {
       if(ts.length != v.length) 
       {System.out.println("Something wrong with lengths\n"); return null;}

       double[] out = new double[v.length];  
       for(int i=0;i<v.length;i++) {out[i] = ts[i] + v[i];}
       return out;
    }



    /*---------------------------------------------------------------------------------     
       Functions to handle uploading data into simulator   
    -----------------------------------------------------------------------------------*/


   public void readCSVData(File file)
   {
       int j = 0;  double val = 0; int i = 0;
       String strline; 
       Double D; double[][] tempdata; 
       int n_toks; int n_series; int count; 
       String delims = "[ ]+";
       String[] tokens; FileInputStream fin; DataInputStream din; BufferedReader br; String names; 
    
       
       try{ //---------------------------------------------------------------------------
           
            fin = new FileInputStream(file);
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            names = br.readLine(); tokens = names.split(delims); n_series = tokens.length;
            tempdata = new double[n_series][1200];
            System.out.println("tokens = " + n_series);

            count = 0; 
            while((strline = br.readLine()) != null)
            {
              tokens = strline.split(delims); 
              n_toks = tokens.length; 
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              System.out.print(count + " ");
              for(i=1;i<=n_series;i++) 
              {
                D = new Double(tokens[i]);
                val = D.doubleValue();
                System.out.print(val + "   ");

                if(count < 1200)
                {
                  tempdata[i-1][count] = val;
                  
                }
                else {System.out.println("Maximum times series length is 500"); break;} 
              }
              System.out.print("\n"); 
              count++; 
            } //---------------------------------------------------------------------------
            
            //---- Clear series queue and set n_obs ------
            if(count != getN_obs())
            {clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;}           
            setN_obs(count);
             
            //System.out.println("N-series = " + n_series + " " + " n_obs = " + count);

            for(i=0;i<n_series;i++) //For each series add in series queue
            {
               real_series = new double[getN_obs()];   
               for(j=0;j<getN_obs();j++)
               {
                 real_series[j] = tempdata[i][j];
               }            
               addSeries(6);    
            }

            din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }



   public void readCSVDataMarketOneAsset(File file)
   {
       int j = 0;  int i = 0;
       String strline; 
       Double D; 
       int n_toks; int count; 
       String delims = "[,]+";
       String[] tokens; FileInputStream fin; DataInputStream din; BufferedReader br; String names; 

       double time;    
       double bid,ask;
       double bid_size, ask_size;
       double volume; 
     
       double[] tempdata = new double[3];
       

       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
           
            fin = new FileInputStream(file);
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            names = br.readLine(); 
            tokens = names.split(delims); 
            market = new ArrayList<double[]>();       
                   
            count = 0; i=0;
            while(i < MAX_DATA && ((strline = br.readLine()) != null))
            {

              tokens = strline.split(delims); 
              n_toks = tokens.length; 
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[0]);
              time = D.doubleValue();

              //---- get bid -----------------
              D = new Double(tokens[1]);
              bid = D.doubleValue();

              D = new Double(tokens[2]);
              bid_size = D.doubleValue();

              D = new Double(tokens[3]);
              ask = D.doubleValue();

              D = new Double(tokens[4]);
              ask_size = D.doubleValue();

              if(tokens[6].equals("NA"))
              {volume = 0;}
              else 
              {
                D = new Double(tokens[6]);
                volume = D.doubleValue(); 
              } 

              //---- price 
              tempdata[0] = time;  
              tempdata[1] = (bid*ask_size + ask*bid_size)/(ask_size + bid_size);
              tempdata[2] = volume;

              market.add(tempdata);
              i++;

            } //---------------------------------------------------------------------------
            count = i;
            //---- Clear series queue and set n_obs ------
            
            //---- Clear series queue and set n_obs ------
            if(count != getN_obs())
            {clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;}           
            setN_obs(count);
             
            volume_series = new double[getN_obs()];
            real_series = new double[getN_obs()];
   
            for(j=0;j<getN_obs();j++)
            {
              tempdata = market.get(j);
              real_series[j] = tempdata[1];
              volume_series[j] = tempdata[2];
            }            
            addSeries(6);    
            addSeries(7);

            din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }




   public void readCSVDataMarketASSET(File file)
   {
       int j = 0;  int i = 0;
       String strline; 
       Double D; 
       int n_toks; int count; 
       String delims = "[ ]+";
       String[] tokens; FileInputStream fin; DataInputStream din; BufferedReader br; String names; 

       double time;    
       double open;
       double close;
       double high;
       double low;
       double volume; 
     
       double[][] tempdata = new double[6][5000];
       

       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
           
            fin = new FileInputStream(file);
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            names = br.readLine(); 
            tokens = names.split(delims); 
            market = new ArrayList<double[]>();       
                   
            count = 0; i=0;
            strline = br.readLine(); //this is junk

            while(i < MAX_DATA && ((strline = br.readLine()) != null))
            {

              tokens = strline.split(delims); 
              n_toks = tokens.length; System.out.println("Number of toks = "+n_toks);
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[0]);
              time = D.doubleValue();

              //---- get bid -----------------
              D = new Double(tokens[1]);
              open = D.doubleValue();

              D = new Double(tokens[2]);
              high = D.doubleValue();

              D = new Double(tokens[3]);
              low = D.doubleValue();

              D = new Double(tokens[4]);
              close = D.doubleValue();

              if(tokens[5].equals("NA"))
              {volume = 0;}
              else 
              {
                D = new Double(tokens[5]);
                volume = D.doubleValue(); 
              } 

              //---- price 
              tempdata[0][i] = time;  
              tempdata[1][i] = open;
              tempdata[2][i] = high;
              tempdata[3][i] = low;
              tempdata[4][i] = close;
              tempdata[5][i] = volume;

              //System.out.println("close = " + close + ", volume = " + volume);
              
              i++;

            } //---------------------------------------------------------------------------
            count = i;
            //---- Clear series queue and set n_obs ------
            
            if(count != getN_obs())
            {clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;}           
            setN_obs(count);
             
            volume_series = new double[getN_obs()];
            real_series = new double[getN_obs()];
   
            for(j=0;j<getN_obs();j++)
            {
              //tempdata = market.get(j);
              real_series[j] = tempdata[4][j];
              volume_series[j] = tempdata[5][j];
              System.out.println("close = " + real_series[j] + ", volume = " + volume_series[j]);
            }            
            addSeries(6);    
            addSeries(7);

            din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }


   /*------------------------------------------------------------

    0 to.minutes
    1 to.minutes3
    2 to.minutes5
    3 to.minutes10
    4 to.minutes15
    5 to.minutes30
    6 to.hourly
    7 to.daily

    ------------------------------------------------------------*/

   public void getExportedData(ArrayList<double[]> list)
   {
     int i;
     int n_series = list.size();

     if(n_series > 0)
     {

        real_series = list.get(0);
        setN_obs(real_series.length); 


       //---------------- clear series and change n_obs --------------
       clearSeries(); 
       setNobs(getN_obs());
 
       deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       getTheatre().setNobs(getN_obs()); 
 
       getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
       timeScrollPane.setViewportView(getTheatre());        
       simData = false;           
       addSeries(6); 
     
       //----- now get the rest of series -------------------------- 
       for(i=1; i < n_series; i++)
       {
        real_series = list.get(i);     
        addSeries(6); 
       }     
    }
   }

   
   public void getQuandlData(String[] symbs, String from, String to, int freq, boolean fina)
   {
     new String("Quandl.auth(\"1HzkGkvirzCUQLkviAYD\")");
     String collapse = null;
     String subset = null;
     
     System.out.println("Freq = " + freq + ",fina = " + fina);
     
     int i,k,count;

     if(freq == 1 || freq == 3) {collapse = ", collapse = 'weekly')";}
     else if(freq == 2 || freq == 4) {collapse = ", collapse = 'monthly')";}
     else {collapse = ")";}
     
     String addInst,instrums,names;
 
     if(fina)
     {
       subset = "c(1:" + (symbs.length-1) + ",ncol(mydata))";
     }
 
 
 
 
     try
     {

      RCaller caller = new RCaller();
      caller.setRscriptExecutable("/usr/bin/Rscript");
      caller.cleanRCode();
      caller.getRCode().addRCode("require (Runiversal)");     
      caller.getRCode().addRCode("require (quantmod)"); 
      caller.getRCode().addRCode("require (Quandl)"); 
      caller.getRCode().addRCode("Quandl.auth(\"1HzkGkvirzCUQLkviAYD\")");
      
      addInst = "Quandl(";
      
      instrums = "c('";
      for(i = 0; i < symbs.length-1; i++)
      { 
         instrums = instrums + symbs[i] + "','";
      }
      instrums = instrums + symbs[symbs.length-1] + "')";
      addInst = addInst + instrums + ",start_date='" + from + "',end_date='" + to + "',type='xts'" + collapse;
      
         
      
      names = "c('Dates','";
      for(i = 0; i < symbs.length-1; i++)
      { 
         names = names + symbs[i] + "','";
      }
      names = names + symbs[symbs.length-1] + "')";
      
      
      
      
      
      caller.getRCode().addRCode("mydata<-" + addInst);      //get the data from quandl (financial data assumed last)
      caller.getRCode().addRCode("mydata<-na.omit(mydata)"); //omit any na's
      if(fina) {caller.getRCode().addRCode("mydata<-mydata[," + subset + "]");} //if using yahoo finanical data, use only final column
     
      caller.getRCode().addRCode("colnames(mydata)<-" + instrums);
      caller.getRCode().addRCode("series_list<-lapply(as.list(mydata), coredata)");
      caller.getRCode().addRCode("write.table(index(mydata),file='/tmp/quandl_dates.txt')");
//       caller.getRCode().addRCode("mydata.df<-data.frame(Date<-index(mydata),coredata(mydata))");
//       caller.getRCode().addRCode("colnames(mydata.df)<-" + names);
//       
// //       caller.getRCode().addRCode("dates<-index(mydata)");
// //       caller.getRCode().addRCode("dates<-data.frame(dates)");
// //       caller.getRCode().addRCode("colnames(mydata.df)<-c('Dates')");
// //       caller.getRCode().addRCode("mydates<-list(dates=dates)");
//       
//       caller.getRCode().addRCode("all<-as.list(as.data.frame(mydata.df))");
//       caller.getRCode().addRCode("all<-c(mydates,series_list)");
      
      caller.runAndReturnResult("series_list");
      
      
      
      double[] target = caller.getParser().getAsDoubleArray(symbs[0]);
      
      setN_obs(target.length);
      clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(getN_obs()); simData = false;
    
      real_series = new double[target.length];
      System.arraycopy(target,0,real_series,0,target.length);
      addSeries(6);
    
      //--- now set differenced series -----------------------
      real_series = new double[target.length];
      real_series[0] = 0;
      
      if(freq == 3 || freq == 4){
      for(i = 1; i < getN_obs(); i++)
      {real_series[i] = (target[i] - target[i-1]);}
      }
      else{
      for(i = 1; i < getN_obs(); i++)
      {real_series[i] = (target[i] - target[i-1])/target[i-1];}
      }
      addSeries(6);
      
      for(k = 1; k < symbs.length; k++)
      {
        target = caller.getParser().getAsDoubleArray(symbs[k]);
        real_series = new double[target.length];
        real_series[0] = 0;
      
        if(freq == 3 || freq == 4){
        for(i = 1; i < getN_obs(); i++)
        {real_series[i] = (target[i] - target[i-1]);}
        }
        else
        {
         for(i = 1; i < getN_obs(); i++)
         {real_series[i] = (target[i] - target[i-1])/target[i-1];}
        }
        addSeries(6);
      }
      
      
      //--- now get dates ------------------------ 
      String[] dates = new String[getN_obs()];
     
      try{
          
         FileInputStream fin = new FileInputStream(new File("/tmp/quandl_dates.txt"));
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         String strline = br.readLine(); count = 0;
         
         while((strline = br.readLine()) != null)
         {
           String[] tokens = strline.split("[ ]+");  
           dates[count] = tokens[1];
           count++;
         }  
     
         din.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
       
       String[] namesall = new String[symbs.length+1];
       namesall[0] = symbs[0];
       
       for(k = 0; k < symbs.length; k++)
       {namesall[k+1] = "d"+symbs[k];}
       
       getTheatre().setNames(namesall);
       getTheatre().setDates(dates);   

    }
    catch(Exception e){System.out.println(e.toString());}
   }
    

   


   public void getMarketData(String instrum, String date1, boolean time_frame, 
                             String time1, String time2, int freq)
   {

     String addInst;  
     String time_seg;
     String agg_time;
     try
    {

     RCaller caller = new RCaller();
     caller.setRscriptExecutable("/usr/bin/Rscript");
     caller.cleanRCode();
     caller.getRCode().addRCode("require (Runiversal)");
     caller.getRCode().addRCode("require (FinancialInstrument)");                    
     //caller.getRCode().addRCode("loadInstruments('/home/lisztian/Data/instruments.rda')");  
     //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/Data/sec',use_identifier='X.RIC',extension='RData')");
     caller.getRCode().addRCode("loadInstruments('/media/My Passport/Market/instruments.rda')"); 
     caller.getRCode().addRCode("setSymbolLookup.FI('/media/My Passport/Market/sec',use_identifier='X.RIC',extension='RData')");
     
     
     //load data in R
     addInst = "getSymbols('" + instrum + "',from='" + date1 + "',indexTZ='America/New York')"; 
     //addInst = "getSymbols('" + instrum + "',from='" + date1 + "')"; 
     caller.getRCode().addRCode(addInst);

     //filter time
     //System.out.println("time frame = " + time_frame);
     if(time_frame)
     {
       time_seg = "['T13:30/T20:00']";
       caller.getRCode().addRCode("mat <- "+instrum + time_seg);
     }   


     //aggregate time 
     if(freq == 0)
     {agg_time = "to.minutes(mat[,5:6],name=NULL)";}
     else if(freq == 1)
     {agg_time = "to.minutes3(mat[,5:6],name=NULL)";}
     else if(freq == 2)
     {agg_time = "to.minutes5(mat[,5:6],name=NULL)";}
     else if(freq == 3)
     {agg_time = "to.minutes10(mat[,5:6],name=NULL)";}
     else if(freq == 4)
     {agg_time = "to.minutes15(mat[,5:6],name=NULL)";}
     else if(freq == 5)
     {agg_time = "to.minutes30(mat[,5:6],name=NULL)";}
     else if(freq == 6)
     {agg_time = "to.hourly(mat[,5:6],name=NULL)";}
     else if(freq == 7)
     {agg_time = "to.daily(mat[,5:6],name=NULL)";}
     else 
     {agg_time = "to.minutes(mat[,5:6],name=NULL)";}

     caller.getRCode().addRCode("series <- " + agg_time);     
     caller.getRCode().addRCode("series_list<-lapply(as.list(series), coredata)");

     //------ experimental ----------------------------------
     //caller.getRCode().addRCode("series2<-coredata(series)");
     //caller.getRCode().addRCode("series_dat<-list(Close=as.vector(series2[,4]),Volume=as.vector(series2[,5]))");

      //---finally, write out to file marketData.csv-------
     //caller.runAndReturnResult("series_dat");
     //caller.getRCode().addRCode("names(series_dat)");
     caller.runAndReturnResult("series_list");

     real_series = caller.getParser().getAsDoubleArray("Close");
     volume_series = caller.getParser().getAsDoubleArray("Volume");

     setN_obs(real_series.length);
     System.out.println("nobs = " + getN_obs());
     clearSeries(); 
     setNobs(getN_obs());
 
     deleteAll.setEnabled(false); 
     acfplotBox.setEnabled(false); 
     getTheatre().setNobs(getN_obs()); 
 
     getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
     timeScrollPane.setViewportView(getTheatre());        

     simData = false;           
    
     addSeries(6);    
     //addSeries(7);

    }
    catch(Exception e){System.out.println(e.toString());}

   }
  
 



   public void getMarketDataMult(String[] instrum, String date1, String date2, boolean time_frame, 
                             String time1, String time2, int freq, 
                             boolean logreturns, boolean vol, boolean clear)
   {

     int i,j;
     String instString,instString2;
     String addInst;  
     double max;     

     double sum,mean;
     double[] temp;
     int diff;
     int n_inst = instrum.length;
     String[] align = new String[n_inst];
     String[] names = new String[2*n_inst+1];
     String na;
     String closes,volume;
     String name_vector,name_vector2;
     

     instString = "c(";
     for(i=0;i<n_inst-1;i++)
     {
       instString = instString + "'" + instrum[i] + "',";
     }
     instString = instString + "'" + instrum[n_inst-1] + "')";


    //names[0] = "asset0";
    for(i=0; i < n_inst; i++)
    {
     names[i] = "asset" + i;
     
     if(freq == 0)
     {align[i] = names[i] + "<-align.time(to.minutes(" + instrum[i] + "[,5:6],1,name=NULL),60)['T13:30/T20:00',]";}

     if(freq == 1)
     {align[i] = names[i] + "<-align.time(to.minutes3(" + instrum[i] + "[,5:6],name=NULL),3*60)['T13:30/T20:00',]";}

     if(freq == 2)
     {align[i] = names[i] + "<-align.time(to.minutes5(" + instrum[i] + "[,5:6],name=NULL),5*60)['T13:30/T20:00',]";}

     if(freq == 3)
     {align[i] = names[i] + "<-align.time(to.minutes10(" + instrum[i] + "[,5:6],name=NULL),10*60)['T13:30/T20:00',]";}

     if(freq == 4)
     {align[i] = names[i] + "<-align.time(to.minutes15(" + instrum[i] + "[,5:6],name=NULL),15*60)['T13:30/T20:00',]";}

     if(freq == 5)
     {align[i] = names[i] + "<-align.time(to.minutes30(" + instrum[i] + "[,5:6],name=NULL),30*60)['T13:30/T20:00',]";}

     if(freq == 6)
     {align[i] = names[i] + "<-align.time(to.hourly(" + instrum[i] + "[,5:6],name=NULL),60*60)['T13:30/T20:00',]";}

     if(freq == 7)
     {align[i] = names[i] + "<-to.daily(" + instrum[i] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 8)
     {align[i] = names[i] + "<-to.weekly(" + instrum[i] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 9)
     {align[i] = names[i] + "<-to.monthly(" + instrum[i] + "[,5:6],indexAt=’yearmon’,drop.time=FALSE,name=NULL)";}

     names[i+n_inst] = "assetv"+i;
   }
    //apply.daily(x, FUN, ...)
   //getSymbols.FI("GOOG.O",from='2011-10-01', to='2012-06-19')

   //GOOG<-align.time(to.period(GOOG.O[,5:6],period='minutes',k=1,name=NULL),60)
   
     na = "(cbind("; name_vector = "c("; 
     //if(logreturns) {name_vector = name_vector + "'asset0',";}

     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Cl("+names[i] + "),";
        name_vector = name_vector + "'" + names[i] + "',";
     }
     na = na + "Cl("+names[n_inst-1]+")))";
     name_vector = name_vector + "'" + names[n_inst-1] + "')"; 
     closes = "closes <- na.locf" + na;
     name_vector = "names(closes) <- " + name_vector;

     
     
     

     volume = "volumes<-na.locf(cbind(";
     name_vector2 = "c("; instString2 = "";
     for(i=0;i<n_inst-1;i++)
     {
       instString2 = instString2 + "Vo(" + names[i] + "),";
       name_vector2 = name_vector2 + "'assetv"+i+"',";
     } 
     instString2 = instString2 + "Vo(" + names[n_inst-1] + ")))";
     name_vector2 = name_vector2 + "'assetv"+(n_inst-1)+"')";

     volume = volume+instString2;
     name_vector2 = "names(volumes) <- " + name_vector2;


    try
    {

     RCaller caller = new RCaller();
     caller.setRscriptExecutable("/usr/bin/Rscript");
     caller.cleanRCode();
     caller.getRCode().addRCode("require (Runiversal)");
     caller.getRCode().addRCode("require (FinancialInstrument)");   

     //caller.getRCode().addRCode("loadInstruments('/home/lisztian/MarketData/instruments.rda')");  
     //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/MarketData/sec',use_identifier='X.RIC',extension='RData')");
     caller.getRCode().addRCode("loadInstruments('/media/My Passport/Market/instruments.rda')"); 
     caller.getRCode().addRCode("setSymbolLookup.FI('/media/My Passport/Market/sec',use_identifier='X.RIC',extension='RData')");
     
     //indexTZ(GOOG.O) <- "America/New_York"
     //GOOG.O<-GOOG.O[!is.na(GOOG.O$Trade.Price)]
     
     //----- Call the symbols ------------------------- getSymbols.FI(" + instString + ",from='" + date1 + "', to='" + date2 + "')
     addInst = "getSymbols.FI(" + instString + ",from='" + date1 + "', to='" + date2 + "')"; 
     caller.getRCode().addRCode(addInst);

     //----- Call the align --------------------------
     for(i=0;i<n_inst;i++)
     { 
        caller.getRCode().addRCode(align[i]);
     }

      //----- define the matrix ------------------------
      caller.getRCode().addRCode(closes);
      caller.getRCode().addRCode(volume);

      //----- define the name_vector--------------------
      caller.getRCode().addRCode(name_vector);
      caller.getRCode().addRCode(name_vector2);

      if(logreturns)
      {

        //--- somehow, get original series, but only first one, we target first in list
        //caller.getRCode().addRCode("times <- index(closes)");
        caller.getRCode().addRCode("closes <- log(closes)");
        caller.getRCode().addRCode("orig <- closes[,1]");
        caller.getRCode().addRCode("names(orig) <- 'orig'"); 
        caller.getRCode().addRCode("orig[is.na(orig)]<-0");
        caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)");                
        caller.getRCode().addRCode("closes <- diff(closes)");
        caller.getRCode().addRCode("volumes <- diff(log(volumes))");
        caller.getRCode().addRCode("closes[1,] <- 0");
        caller.getRCode().addRCode("volumes[1,] <- 0");
        caller.getRCode().addRCode("closes[is.na(closes)]<-0");  //--- set any final NAs to 0
        caller.getRCode().addRCode("volumes[is.na(volumes)]<-0");  //--- set any final NAs to 0
       
      }
 
      caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
      caller.getRCode().addRCode("volume_list<-lapply(as.list(volumes), coredata)"); 
         
      if(logreturns)
      {caller.getRCode().addRCode("all_list<-c(orig_list,series_list,volume_list)");}
      else
      {caller.getRCode().addRCode("all_list<-c(series_list,volume_list)");}

      caller.runAndReturnResult("all_list");


      //------ Get first series (if logreturns, then first is original series)
      
      if(logreturns) 
      {
        real_series = caller.getParser().getAsDoubleArray("orig");
        logprice_data = new double[real_series.length];
 
        System.arraycopy(real_series,0,logprice_data,0,real_series.length);

      }      
      else{real_series = caller.getParser().getAsDoubleArray(names[0]);}
           
      //System.out.println("Gets here0");
      //---------------- clear series and change n_obs --------------



      if(clear)
      {
       setN_obs(real_series.length); 
       clearSeries(); 
       setNobs(getN_obs());
       deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       getTheatre().setNobs(getN_obs());
       getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
       timeScrollPane.setViewportView(getTheatre());
      }
      // else, keep n_obs fixed, assumes n_obs of new data equal or greater

      simData = false;           
      //addSeries(6); 
     
      int start; 
      if(logreturns) {start = 0;} 
      else {start = 1;}

      //----- now get the rest of series -------------------------- 

     if(logreturns)
     {
      for(i=start; i < n_inst; i++)
      {
       sum = 0.0;
       real_series = caller.getParser().getAsDoubleArray(names[i]);

       for(j=0;j<getN_obs();j++) 
       {sum = sum + real_series[j];}

       mean = sum/getN_obs(); sum = 0.0;

       for(j=0;j<getN_obs();j++) 
       {sum = sum + (real_series[j]-mean)*(real_series[j]-mean);}
       var = Math.sqrt(sum);

       //for(j=0;j<n_obs;j++) 
       //{real_series[j] = real_series[j]/var; arr[j][i] = real_series[j];} 
       addSeries(6); 
      }
     }
     else
     {
      for(i=start; i < n_inst; i++)
      {
       
       real_series = caller.getParser().getAsDoubleArray(names[i]);
       //for(j=0;j<n_obs;j++) {arr[j][i+1] = real_series[j];} 
       addSeries(6); 
      }
     }
       
      if(vol)
      {
        for(i=0; i < n_inst; i++)
        {
         //System.out.println("Gets here " + n_inst);
         temp = caller.getParser().getAsDoubleArray(names[i+n_inst]);
 
         volume_series = new double[getN_obs()];
         diff = temp.length - getN_obs(); max = -1;
 
         for(j=0;j<getN_obs();j++) 
         {
           volume_series[j] = temp[diff+j];           
           if(volume_series[j] > max) {max = volume_series[j];}           
         }
         for(j=0;j<getN_obs();j++) 
         {
           volume_series[j] = volume_series[j]/max; 
           //arr[j][i+n_inst] = volume_series[j]; 
         }

         addSeries(7);
        } 
      }
 

     //if(logreturns)
    // {
       /*String file = "logprice_"+instrum[0]+".dat";       
       PrintWriter out = new PrintWriter(new FileWriter(file));
   
     
     for(i=0;i<n_obs;i++)
     {
        for(j=0;j<sim_orig.size()-1;j++)
        {
          out.print(arr[i][j] + "  "); 
        }
        out.println(arr[i][sim_orig.size()-1]); 
     }
     
     
       //for(i=0;i<n_obs;i++) {out.println(logprice_data[i]);}   
       out.close();*/
      //}
     
   }
     
  
  
   
/*
   for(i=0;i<n_obs;i++)
   {

     for(j=0;j<sim_orig.size();j++)
     {
         sim_orig.get(j);

     }

   }
*/


   /*scaled_design_matrix <- function() {
    getSymbols(c('DXM2','ESM2','STXEM2'),
        from='2012-03-19',
        to='2012-06-15',
        indexTZ='America/New York')

    dxy <- align.time(to.hourly(  DXM2[,5:6], name=NULL), 3600)['T09:00/T17:00',]
    es  <- align.time(to.hourly(  ESM2[,5:6], name=NULL), 3600)['T09:00/T17:00',]
    stxe<- align.time(to.hourly(STXEM2[,5:6], name=NULL), 3600)['T09:00/T17:00',]

    closes <- na.locf(cbind(Cl(es),Cl(es),Cl(dxy),Cl(stxe)))
    names(closes) <- c('y', 'es', 'dxy', 'stxe')

    log_lvls <- log(closes)
    log_rtns <- diff(log_lvls)
    log_rtns[1,] <- 0

    means <- list(dxy=mean(log_rtns$dxy), es=mean(log_rtns$es), stxe=mean(log_rtns$stxe))
    scaled_lret <- scale(log_rtns, scale=TRUE, center=TRUE)

    list(x=scaled_lret, means=means) */

    
    catch(Exception e){System.out.println(e.toString());}

   }


   public void getMarketDataGeneral(String[] instrum, String date1, String date2, boolean time_frame, 
                             String time1, String time2, int freq, 
                             boolean logreturns, boolean vol, boolean high_low, boolean high_low_diff, boolean clear)
   {

     int i;
     String instString,instString2;
     String addInst;  
     int n_inst = instrum.length;
     String[] align = new String[n_inst];
     String[] names = new String[2*n_inst+1];
     String na;
     String closes,volume;
     String name_vector,name_vector2,name_vector3,name_vector4,name_vector5,highs,lows;
     

     instString = "c(";
     for(i=0;i<n_inst-1;i++)
     {
       instString = instString + "'" + instrum[i] + "',";
     }
     instString = instString + "'" + instrum[n_inst-1] + "')";


    //names[0] = "asset0";
    for(i=0; i < n_inst; i++)
    {
     names[i] = "asset" + i;
     
     if(freq == 0)
     {align[i] = names[i] + "<-align.time(to.minutes(" + instrum[i] + "[,5:6],1,name=NULL),60)['T13:30/T20:00',]";}

     if(freq == 1)
     {align[i] = names[i] + "<-align.time(to.minutes3(" + instrum[i] + "[,5:6],name=NULL),3*60)['T13:30/T20:00',]";}

     if(freq == 2)
     {align[i] = names[i] + "<-align.time(to.minutes5(" + instrum[i] + "[,5:6],name=NULL),5*60)['T13:30/T20:00',]";}

     if(freq == 3)
     {align[i] = names[i] + "<-align.time(to.minutes10(" + instrum[i] + "[,5:6],name=NULL),10*60)['T13:30/T20:00',]";}

     if(freq == 4)
     {align[i] = names[i] + "<-align.time(to.minutes15(" + instrum[i] + "[,5:6],name=NULL),15*60)['T13:30/T20:00',]";}

     if(freq == 5)
     {align[i] = names[i] + "<-align.time(to.minutes30(" + instrum[i] + "[,5:6],name=NULL),30*60)['T13:30/T20:00',]";}

     if(freq == 6)
     {align[i] = names[i] + "<-align.time(to.hourly(" + instrum[i] + "[,5:6],name=NULL),60*60)['T13:30/T20:00',]";}

     if(freq == 7)
     {align[i] = names[i] + "<-to.daily(" + instrum[i] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 8)
     {align[i] = names[i] + "<-to.weekly(" + instrum[i] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 9)
     {align[i] = names[i] + "<-to.monthly(" + instrum[i] + "[,5:6],indexAt=’yearmon’,drop.time=FALSE,name=NULL)";}

     names[i+n_inst] = "assetv"+i;
   }

   
   
   //Get Close data -------------
   
     na = "(cbind("; name_vector = "c("; 
     //if(logreturns) {name_vector = name_vector + "'asset0',";}

     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Cl("+names[i] + "),";
        name_vector = name_vector + "'" + names[i] + "',";
     }
     na = na + "Cl("+names[n_inst-1]+")))";
     name_vector = name_vector + "'" + names[n_inst-1] + "')"; 
     closes = "closes <- na.locf" + na;
     name_vector = "names(closes) <- " + name_vector;

   //Get High and Low data-------
     na = "(cbind("; name_vector3 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Hi("+names[i] + "),";
        name_vector3 = name_vector3 + "'assetHi" +i+"',";
     }
     na = na + "Hi("+names[n_inst-1]+")))";
     name_vector3 = name_vector3 + "'assetHi"+(n_inst-1)+"')"; 
     highs = "highs <- na.locf" + na;
     name_vector3 = "names(highs) <- " + name_vector3;  
    
    
     na = "(cbind("; name_vector4 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Lo("+names[i] + "),";
        name_vector4 = name_vector4 + "'assetLo" +i+"',";
     }
     na = na + "Lo("+names[n_inst-1]+")))";
     name_vector4 = name_vector4 + "'assetLo"+(n_inst-1)+"')"; 
     lows = "lows <- na.locf" + na;
     name_vector4 = "names(lows) <- " + name_vector4; 
     
     
 
     //Get volume
     volume = "volumes<-na.locf(cbind(";
     name_vector2 = "c("; instString2 = "";
     name_vector5 = "c("; 
     for(i=0;i<n_inst-1;i++)
     {
       instString2 = instString2 + "Vo(" + names[i] + "),";
       name_vector2 = name_vector2 + "'assetv"+i+"',";
       name_vector5 = name_vector5 + "'assetPrice"+i+"',";
     } 
     instString2 = instString2 + "Vo(" + names[n_inst-1] + ")))";
     name_vector2 = name_vector2 + "'assetv"+(n_inst-1)+"')";
     name_vector5 = name_vector5 + "'assetPrice"+(n_inst-1)+"')";
     volume = volume+instString2;
     name_vector2 = "names(volumes) <- " + name_vector2;
     name_vector5 = "names(orig)<- " + name_vector5;

    try
    {

     RCaller caller = new RCaller();
     caller.setRscriptExecutable("/usr/bin/Rscript");
     caller.cleanRCode();
     caller.getRCode().addRCode("require (Runiversal)");
     caller.getRCode().addRCode("require (FinancialInstrument)");   

     //caller.getRCode().addRCode("loadInstruments('/home/lisztian/MarketData/instruments.rda')");  
     //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/MarketData/sec',use_identifier='X.RIC',extension='RData')");
     caller.getRCode().addRCode("loadInstruments('/media/My Passport/Market/instruments.rda')"); 
     caller.getRCode().addRCode("setSymbolLookup.FI('/media/My Passport/Market/sec',use_identifier='X.RIC',extension='RData')");
     
     //indexTZ(GOOG.O) <- "America/New_York"
     //GOOG.O<-GOOG.O[!is.na(GOOG.O$Trade.Price)]
     
     //----- Call the symbols -------------------------
     addInst = "getSymbols.FI(" + instString + ",from='" + date1 + "', to='" + date2 + "')"; 
     caller.getRCode().addRCode(addInst);

     //----- Call the align --------------------------
     for(i=0;i<n_inst;i++)
     { 
        caller.getRCode().addRCode(align[i]);
     }

      //----- define the matrix ------------------------
      caller.getRCode().addRCode(closes);
      caller.getRCode().addRCode(highs);
      caller.getRCode().addRCode(lows);
      caller.getRCode().addRCode(volume);

      //----- define the name_vector--------------------
      caller.getRCode().addRCode(name_vector);
      caller.getRCode().addRCode(name_vector3);
      caller.getRCode().addRCode(name_vector4);      
      caller.getRCode().addRCode(name_vector2);


      //--- somehow, get original series, but only first one, we target first in list
      //caller.getRCode().addRCode("times <- index(closes)");
      caller.getRCode().addRCode("closes <- log(closes)");
      caller.getRCode().addRCode("highs<-log(highs)");
      caller.getRCode().addRCode("lows<-log(lows)");
      caller.getRCode().addRCode("highs<-highs-lows");
      
      caller.getRCode().addRCode("orig <- closes");
      caller.getRCode().addRCode(name_vector5); 
      caller.getRCode().addRCode("orig[is.na(orig)]<-0");
      caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)"); 
      
      caller.getRCode().addRCode("closes <- diff(closes)");
      caller.getRCode().addRCode("volumes <- diff(volumes)");

      if(high_low_diff)
      {caller.getRCode().addRCode("highs<-diff(highs)"); caller.getRCode().addRCode("highs[1,] <- 0");}
      
      caller.getRCode().addRCode("closes[1,] <- 0");
      caller.getRCode().addRCode("volumes[1,] <- 0");
      caller.getRCode().addRCode("closes[is.na(closes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("volumes[is.na(volumes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("highs[is.na(highs)]<-0");  //--- set any final NAs to 0 
      
      
      caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
      caller.getRCode().addRCode("volume_list<-lapply(as.list(volumes), coredata)"); 
      caller.getRCode().addRCode("highs_list<-lapply(as.list(highs), coredata)");
      
      caller.getRCode().addRCode("all_list<-c(orig_list,series_list,highs_list,volume_list)");


      caller.runAndReturnResult("all_list");


    
      
      real_series = caller.getParser().getAsDoubleArray("assetPrice0");
      if(clear)
      { 
       setN_obs(real_series.length);
       clearSeries(); 
       setNobs(getN_obs());
       deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       getTheatre().setNobs(getN_obs());
       getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
       timeScrollPane.setViewportView(getTheatre());
      }       
      
      //Get Prices first 
      for(i=0;i<n_inst;i++)
      {
         real_series = caller.getParser().getAsDoubleArray("assetPrice"+i);
         setN_obs(real_series.length); 
         addSeries(6);
      }

      //Get log-retuns
      if(logreturns)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("asset"+i);
         addSeries(6);
       }
      }
      
      if(high_low || high_low_diff)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetHi"+i);
         addSeries(6);
       } 
      }
      
      if(vol)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetv"+i);
         addSeries(6);
       }
      }
    }
    catch(Exception e){System.out.println(e.toString());}

   }   
   
   
   public void getHF_CSV_Data(int ndays, int last, int n, boolean volu, boolean high_low, boolean high_low_diff)
   {

       int i = 0; int k = 0;
       String strline; 
       Double D; 
     
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       double volume, end_ask, end_bid, high_ask, high_bid, low_ask, low_bid;
       double log_price1, log_price2;
       int nobs = ndays*n;
 
       int start = last-ndays;      
       double[] log_price = new double[nobs];
       double[] vol = new double[nobs];
       double[] high_low_data = new double[nobs];
       double[] diffhighlow = new double[nobs];

       if(last - ndays > 0) {
       //"volume","bid","bid.size","ask","ask.size","bid2","bid2.size","ask2","ask2.size","bid3","bid3.size","ask3","ask3.size","bid4","bid4.size","ask4","ask4.size","bid5","bid5.size","ask5","ask5.size"
       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
         
         for(k = 0; k < ndays; k++)
         {
              
            fin = new FileInputStream(new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/IF5mtwo/"+(start+k+1)+".front.csv"));
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            br.readLine(); 
            i=0;

            while(i < n && ((strline = br.readLine()) != null))
            {

              tokens = strline.split(delims); 
              n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[1]);
              volume = D.doubleValue();

              //---- get bid -----------------
              D = new Double(tokens[2]);
              high_bid = D.doubleValue(); 

              D = new Double(tokens[3]);
              low_bid = D.doubleValue();

              D = new Double(tokens[4]);
              high_ask = D.doubleValue();

              D = new Double(tokens[5]);
              low_ask = D.doubleValue();

              D = new Double(tokens[6]);
              D.doubleValue();

              D = new Double(tokens[7]);
              end_bid = D.doubleValue();

              D = new Double(tokens[8]);
              D.doubleValue();

              D = new Double(tokens[9]);
              end_ask = D.doubleValue();
 
         
              log_price1 = (end_bid + end_ask)/2.0;
              log_price2 = (high_bid + high_ask)/2.0 - (low_bid + low_ask)/2.0;

              //System.out.println(log_price1 + " " + log_price2);
              
              //---- price 

              high_low_data[n*k+i] = .003*Math.log(log_price2);
              log_price[n*k+i] = Math.log(log_price1);   
              vol[n*k+i] = volume;
           
              i++;
           }
           din.close();
          }
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
         double[] logdiff = new double[ndays*n];
         logdiff[0] = 0.0; int total = ndays*n-1;
         for(i=0;i<total;i++)
         {
           logdiff[i+1] = log_price[i+1] - log_price[i];
           diffhighlow[i+1] = high_low_data[i+1] - high_low_data[i];
           if(Double.isNaN(diffhighlow[i+1]) || Double.isInfinite(diffhighlow[i+1])) 
           {diffhighlow[i+1] = 0.0;}
           //System.out.println(diffhighlow[i+1]);
         }

         setN_obs(nobs);
         clearSeries(); 
         setNobs(getN_obs());
         deleteAll.setEnabled(false); 
         acfplotBox.setEnabled(false); 
         getTheatre().setNobs(getN_obs());
         getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
         timeScrollPane.setViewportView(getTheatre());         
         
         
         real_series = log_price;
         addSeries(6);

         real_series = logdiff;
         addSeries(6); 
         
         real_series = diffhighlow;
         addSeries(6);

     
     
       log_price = new double[nobs];
       vol = new double[nobs];
       high_low_data = new double[nobs];
       diffhighlow = new double[nobs];
     
       //Now do back   
        
       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
         
         for(k = 0; k < ndays; k++)
         {
              
            fin = new FileInputStream(new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/SP_Data/"+(start+k+1)+".csv"));
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            br.readLine(); 
            i=0;

            while(i < n && ((strline = br.readLine()) != null))
            {

              tokens = strline.split(delims); 
              n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[1]);
              volume = D.doubleValue();

              //---- get bid -----------------
              D = new Double(tokens[2]);
              high_bid = D.doubleValue(); 

              D = new Double(tokens[3]);
              low_bid = D.doubleValue();

              D = new Double(tokens[4]);
              high_ask = D.doubleValue();

              D = new Double(tokens[5]);
              low_ask = D.doubleValue();

              D = new Double(tokens[6]);
              D.doubleValue();

              D = new Double(tokens[7]);
              end_bid = D.doubleValue();

              D = new Double(tokens[8]);
              D.doubleValue();

              D = new Double(tokens[9]);
              end_ask = D.doubleValue();
 
         
              log_price1 = (end_bid + end_ask)/2.0;
              log_price2 = (high_bid + high_ask)/2.0 - (low_bid + low_ask)/2.0;

              //---- price 

              high_low_data[n*k+i] = .003*Math.log(log_price2);
              log_price[n*k+i] = Math.log(log_price1);   
              vol[n*k+i] = volume;
           
              i++;
           }
           din.close();
          }
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
         logdiff = new double[ndays*n];
         logdiff[0] = 0.0; 
         for(i=0;i<total;i++)
         {
           logdiff[i+1] = log_price[i+1] - log_price[i];
           diffhighlow[i+1] = high_low_data[i+1] - high_low_data[i];
           if(Double.isNaN(diffhighlow[i+1]) || Double.isInfinite(diffhighlow[i+1])) 
           {diffhighlow[i+1] = 0.0;}
           //System.out.println(diffhighlow[i+1]);
         }

         
         real_series = log_price;
         addSeries(6);

         real_series = logdiff;
         addSeries(6); 
         
         real_series = diffhighlow;
         addSeries(6);         
         
         
         //plotData2(log_ret_data[0],log_ret_data[1], n);
         //plotData(logdiff,nobs);
         //plotData(log_price, nobs);
         //plotData(diffhighlow,nobs);
         //plotData(vol,nobs);
      }
   } 


   public void getSP_CSV_Data(int ndays, int last, int n, boolean volu, boolean high_low, boolean high_low_diff)
   {

       int i = 0; int k = 0;
       String strline; 
       Double D; 
     
       int n_toks; int count; 
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       double high_bid, low_bid, close;
       int nobs = ndays*n;
 
       int start = last-ndays;      
       double[] log_ret_data = new double[nobs];
       double[] log_price = new double[nobs];
       double[] log_ret_russ = new double[nobs];
       double[] log_ret_nasd = new double[nobs];       
       double[] high_low_data = new double[nobs];
       double[] diffhighlow = new double[nobs];

       double[] log_priceNas = new double[nobs];
       double[] log_priceRus = new double[nobs];

       if(last - ndays > 0) {
       //"volume","bid","bid.size","ask","ask.size","bid2","bid2.size","ask2","ask2.size","bid3","bid3.size","ask3","ask3.size","bid4","bid4.size","ask4","ask4.size","bid5","bid5.size","ask5","ask5.size"
       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
         count = 0;
         for(k = 0; k < ndays; k++)
         {
              
            fin = new FileInputStream(new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/SP_Data/day"+(start+k+1)+".csv"));
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            //names = br.readLine(); 
            i=0;

            while(i < n && ((strline = br.readLine()) != null))
            {

              tokens = strline.split(delims); 
              n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
  
              D = new Double(tokens[3]);
              close = D.doubleValue();
              log_price[count] = Math.log(close); 
              
              //---- get highlow -----------------
              D = new Double(tokens[4]);
              high_bid = D.doubleValue(); 
              
              D = new Double(tokens[5]);
              low_bid = D.doubleValue();
              high_low_data[count] = Math.log(high_bid - low_bid);
 
              D = new Double(tokens[8]);
              close = D.doubleValue();
              log_priceNas[count] = Math.log(close);  
           
              D = new Double(tokens[12]);
              close = D.doubleValue();
              log_priceRus[count] = Math.log(close);               
            
              i++; count++;
           }
           din.close();
          }
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
        
         for(i=0;i<nobs-1;i++)
         {
           log_ret_data[i+1] = log_price[i+1] - log_price[i];
           diffhighlow[i+1] = .002*(high_low_data[i+1] - high_low_data[i]);
           if(Double.isNaN(diffhighlow[i+1]) || Double.isInfinite(diffhighlow[i+1])) 
           {diffhighlow[i+1] = 0.0;}
           
           log_ret_russ[i+1] = log_priceRus[i+1] - log_priceRus[i];
           log_ret_nasd[i+1] = log_priceNas[i+1] - log_priceNas[i];           
         }

         setN_obs(nobs);
         clearSeries(); 
         setNobs(getN_obs());
         deleteAll.setEnabled(false); 
         acfplotBox.setEnabled(false); 
         getTheatre().setNobs(getN_obs());
         getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
         timeScrollPane.setViewportView(getTheatre());         
         
         
         real_series = log_price;
         addSeries(6);

         real_series = log_ret_data;
         addSeries(6); 
         
         real_series = diffhighlow;
         addSeries(6);
         
         real_series = log_ret_nasd;
         addSeries(6);
      
         real_series = log_ret_russ;
         addSeries(6);   

      }
      
   }
   

   public void getTWS_CSV_Data(File file)
   {

       int i = 0; 
       String strline; 
       Double D; 
     
       int n_toks=4; int count; 
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       ArrayList<Double> close_series = new ArrayList<Double>();
       ArrayList<Double> highlow_series = new ArrayList<Double>();    
       ArrayList<Double> exp_series_1 = new ArrayList<Double>();  
       new ArrayList<Double>();
       ArrayList<Double> price = new ArrayList<Double>();    
       
       count=0;
       try{
       fin = new FileInputStream(file);
       din = new DataInputStream(fin);
       br = new BufferedReader(new InputStreamReader(din));  

       
       while((strline = br.readLine()) != null)
       {

          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          D = new Double(tokens[1]);
          price.add(D);
          
          D = new Double(tokens[2]);
          close_series.add(D);
 
          if(n_toks > 3)
          {D = new Double(tokens[3]);
          highlow_series.add(D);} 
 
          if(n_toks > 4)
          {D = new Double(tokens[4]);
          exp_series_1.add(D);} 

//           D = new Double(tokens[5]);
//           exp_series_2.add(D);          
           
          count++;
       }
       din.close();
       }catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
        
         setN_obs(count); 
         clearSeries(); 
         setNobs(getN_obs());
         deleteAll.setEnabled(false); 
         acfplotBox.setEnabled(false); 
         getTheatre().setNobs(getN_obs());
         getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
         timeScrollPane.setViewportView(getTheatre());         
         
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = price.get(i).doubleValue();}
         addSeries(6);
         
         for(i=0;i<getN_obs();i++)
         {real_series[i] = close_series.get(i).doubleValue();}
         addSeries(6);
         
         if(n_toks > 3)
         {for(i=0;i<getN_obs();i++)
         {real_series[i] = highlow_series.get(i).doubleValue();}
         addSeries(6);}
        
         if(n_toks > 4)
         {for(i=0;i<getN_obs();i++)
         {real_series[i] = exp_series_1.get(i).doubleValue();}
         addSeries(6);}        

/*         for(i=0;i<n_obs;i++)
         {real_series[i] = exp_series_2.get(i).doubleValue();}
         addSeries(6);     */    
         
      
      
   }   
   
   
   /*---Assumes file format in the following form----
    
    Builds time series based on time x given in argument.
    The series will be built as follows
    
    Price at time x 
    (log) difference at time x
    Observed High Price between x_{t-1} and x
    Observed Low Price between x_{t-1} and x
    (log) difference of high
    (log) difference of low
    
    
    Only the first 4 columns are considered
   
    Date Time, Mid, Bid, Ask, ...... 
   
   
   --------------------------------------------------*/
   
   public void readIntradayBarData(File file, String startDate, String endDate, String timeStamp)
   {
   
       int i = 0; 
       String strline; 
       int indx;
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       double high, low, current_price;
       
       ArrayList<Double> price = new ArrayList<Double>();
       ArrayList<Double> low_price = new ArrayList<Double>();    
       ArrayList<Double> high_price = new ArrayList<Double>();  
       ArrayList<Double> all_prices = new ArrayList<Double>();
       ArrayList<String> all_times = new ArrayList<String>();
       
       boolean ready = false;
       boolean continue_on = false;
       
     
       try
       {
        fin = new FileInputStream(file);
        din = new DataInputStream(fin);
        br = new BufferedReader(new InputStreamReader(din));  
   
        while((strline = br.readLine()) != null)
        {
        
         if(strline.indexOf(startDate) != -1 || continue_on)
         {
          continue_on = true;
          System.out.println(strline);
          
          tokens = strline.split(delims); 
          //collect all the prices  
          all_prices.add(new Double(tokens[1]));
          all_times.add(tokens[0]);
   
          if(tokens[0].indexOf(timeStamp) != -1)
          {
          
            if(ready) //collect the price and look back until previous time x to get high,low
            {
            
             current_price = (new Double(tokens[1])).doubleValue();
             price.add(current_price);

             
             high = current_price;
             low = current_price;             
             
             indx = 1;
             while(all_times.get(all_times.size()-1-indx).indexOf(timeStamp) == -1)
             {
 
              current_price = all_prices.get(all_prices.size()-1-indx);
              
              if(current_price < low) {low = current_price;}
              else if(current_price > high) {high = current_price;}
 
              indx++; 
             }

             System.out.println(low + " " + high);
             all_prices.clear(); 
             all_times.clear();
            
             all_prices.add(new Double(tokens[1]));
             all_times.add(tokens[0]);
                       
             low_price.add(low);
             high_price.add(high);          
                       
            }
          
            
            if(price.size() == 0)
            {ready = true; System.out.println("Ready!");}
            
          }
          
         }
   
        }
        din.close();
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
      
      setN_obs(price.size()); 
      clearSeries(); 
      setNobs(getN_obs());
      deleteAll.setEnabled(false); 
      acfplotBox.setEnabled(false); 
      getTheatre().setNobs(getN_obs());
      getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
      timeScrollPane.setViewportView(getTheatre());         
      
      
      //now load the data into the simPanel
      real_series = new double[getN_obs()];
      
      for(i = 0; i < getN_obs(); i++)
      {real_series[i] = price.get(i);}      
      addSeries(6);
      
      for(i = 0; i < getN_obs(); i++)
      {real_series[i] = low_price.get(i);}
      addSeries(6);
      
      for(i = 0; i < getN_obs(); i++)
      {real_series[i] = high_price.get(i);}      
      addSeries(6);
      
      //now do the differences
      real_series[0] = 0;
      for(i = 1; i < getN_obs(); i++)
      {real_series[i] = price.get(i) - price.get(i-1);}      
      addSeries(6);
      
      real_series[0] = 0;
      for(i = 1; i < getN_obs(); i++)
      {real_series[i] = low_price.get(i) - low_price.get(i-1);}
      addSeries(6);
      
      real_series[0] = 0;
      for(i = 1; i < getN_obs(); i++)
      {real_series[i] = high_price.get(i) - high_price.get(i-1);}      
      addSeries(6); 
 
   
   }
   
   
   
   public void getHK_CSV_Data()
   {

       int i = 0; 
       String strline; 
       Double D; 
     
       int n_toks; int count; 
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       ArrayList<Double> close_series = new ArrayList<Double>();
       ArrayList<Double> highlow_series = new ArrayList<Double>();    
       ArrayList<Double> exp_series_1 = new ArrayList<Double>();  
       ArrayList<Double> exp_series_2 = new ArrayList<Double>();
       ArrayList<Double> price = new ArrayList<Double>();    
       
       count=0;
       try{
       fin = new FileInputStream(new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/IF1306.csv"));
       din = new DataInputStream(fin);
       br = new BufferedReader(new InputStreamReader(din));  

       
       while((strline = br.readLine()) != null)
       {

          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          D = new Double(tokens[2]);
          price.add(D);
          
          D = new Double(tokens[3]);
          close_series.add(D);
 
          D = new Double(tokens[4]);
          highlow_series.add(D); 
 
          D = new Double(tokens[5]);
          exp_series_1.add(D); 

          D = new Double(tokens[6]);
          exp_series_2.add(D);          
           
          count++;
       }
       din.close();
       }catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
        
         setN_obs(count); 
         clearSeries(); 
         setNobs(getN_obs());
         deleteAll.setEnabled(false); 
         acfplotBox.setEnabled(false); 
         getTheatre().setNobs(getN_obs());
         getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
         timeScrollPane.setViewportView(getTheatre());         
         
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = price.get(i).doubleValue();}
         addSeries(6);
         
         for(i=0;i<getN_obs();i++)
         {real_series[i] = close_series.get(i).doubleValue();}
         addSeries(6);
         
         for(i=0;i<getN_obs();i++)
         {real_series[i] = highlow_series.get(i).doubleValue();}
         addSeries(6);
        
         for(i=0;i<getN_obs();i++)
         {real_series[i] = exp_series_1.get(i).doubleValue();}
         addSeries(6);        

         for(i=0;i<getN_obs();i++)
         {real_series[i] = exp_series_2.get(i).doubleValue();}
         addSeries(6);         
         
      
      
   }      
   
   
   
   public void getFutures_CSV_Data(int start, int end, int n, boolean high_low, boolean high_low_diff)
   {

       int i = 0; 
       String strline; 
       Double D; 
     
       int n_toks; int count; 
       String delims = "[,]+";
       String[] tokens; 
        
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       double high_bid, low_bid;
       double close;
       int nobs = n;
 
 
       double[] log_ret_data = new double[nobs];
       double[] log_price = new double[nobs];
       double[] log_ret_russ = new double[nobs];
       double[] log_ret_nasd = new double[nobs];       
       double[] high_low_data = new double[nobs];
       double[] diffhighlow = new double[nobs];

       double[] log_priceNas = new double[nobs];
       double[] log_priceRus = new double[nobs];
   
       //"volume","bid","bid.size","ask","ask.size","bid2","bid2.size","ask2","ask2.size","bid3","bid3.size","ask3","ask3.size","bid4","bid4.size","ask4","ask4.size","bid5","bid5.size","ask5","ask5.size"
       try{ //"","Bid.Price","Bid.Size","Ask.Price","Ask.Size","Trade.Price","Volume"----------------
         
  
            fin = new FileInputStream(new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/Futures5-6-13.csv"));
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            br.readLine(); 
            count = 0; i=0;

            while(i < 1560 && ((strline = br.readLine()) != null))
            {
              
             //5/2/2013 11:15,5/2/2013 11:20,1592.25,1593.00,1593.25,1592.00,"7,671",,2906.50,2908.25,2908.50,2905.50,"1,332",,936.80,937.40,937.50,936.40,588
             if(i>=start && i <= end)
             {
              tokens = strline.split(delims); 
              n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[3]);
              close = D.doubleValue();
              log_price[count] = Math.log(close); 
              
              //---- get highlow -----------------
              D = new Double(tokens[4]);
              high_bid = D.doubleValue(); 
              D = new Double(tokens[5]);
              low_bid = D.doubleValue();
              high_low_data[count] = Math.log(high_bid - low_bid);
 
              D = new Double(tokens[8]);
              close = D.doubleValue();
              log_priceNas[count] = Math.log(close);  
           
              D = new Double(tokens[12]);
              close = D.doubleValue();
              log_priceRus[count] = Math.log(close);               
           
              count++;
             }
             i++;
           }
           din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
        
         for(i=0;i<n-1;i++)
         {
           log_ret_data[i+1] = log_price[i+1] - log_price[i];
           diffhighlow[i+1] = .002*(high_low_data[i+1] - high_low_data[i]);
           if(Double.isNaN(diffhighlow[i+1]) || Double.isInfinite(diffhighlow[i+1])) 
           {diffhighlow[i+1] = 0.0;}
           
           log_ret_russ[i+1] = log_priceRus[i+1] - log_priceRus[i];
           log_ret_nasd[i+1] = log_priceNas[i+1] - log_priceNas[i];           
         }

         setN_obs(nobs);
         clearSeries(); 
         setNobs(getN_obs());
         deleteAll.setEnabled(false); 
         acfplotBox.setEnabled(false); 
         getTheatre().setNobs(getN_obs());
         getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
         timeScrollPane.setViewportView(getTheatre());         
         
         
         real_series = log_price;
         addSeries(6);

         real_series = log_ret_data;
         addSeries(6); 
         
         real_series = diffhighlow;
         addSeries(6);
         
         real_series = log_ret_nasd;
         addSeries(6);

         real_series = log_ret_russ;
         addSeries(6);
         
   }   
   
   
   
  

   /* -------------------------------------------

     The idea here is to take all explanatory components 
     and apply forecastable component analysis to them 

   ---------------------------------------------*/
/*    public void foreCA(int ncomps)
   {
     int n_exp = 0; int i;

     //---- get number of expVariables
     for(i=0;i<max_series;i++)
     {
       if(repCheck[i].isSelected() && repCheck[i].isEnabled()) {n_exp++;}
     }
     double[] expSeries = new double[n_exp*n_obs];

     //----- now get the data into one array
     n_exp = 0;
     for(i=0;i<max_series;i++)
     {
       if(repCheck[i].isSelected() && repCheck[i].isEnabled())
       {
         temp = simData.get(i);
         System.arraycopy(temp,0,expSeries,n_exp*n_obs,n_obs);
         n_exp++;
       }
     }
 

      

     try
     {

      RCaller caller = new RCaller();
      caller.setRscriptExecutable("/usr/bin/Rscript");
      caller.cleanRCode();
      caller.getRCode().addRCode("require (Runiversal)");
      caller.getRCode().addRCode("require (foreCA)");

      caller.addDoubleArray("x", expSeries); //set series
      caller.getRCode().addRCode("gaussx<-Gaussianize(x,type='hh',method='IGMM')"); //Gaussianize
      caller.getRCode().addRCode("my.all<-list(gaussianized=gaussx)");




     }
     catch

   }
   
   
*/   


   public void getMarketDataYahoo(String[] instrum, String date1, boolean time_frame, 
                             String time1, String time2, int freq, 
                             boolean logreturns, boolean vol, boolean clear, boolean high_low, boolean high_low_diff)
   {

     int i;
     String addInst;  
     int n_inst = instrum.length;
     String[] names = new String[2*n_inst];
     String na;
     String closes,volume;
     String name_vector,name_vector2,name_vector3,name_vector4,name_vector5,highs,lows;
     String assetString, instString2;
     String[] transfer = new String[n_inst];
     //---- create matrix of closes------------------------
  
     for(i=0;i<n_inst;i++)
     {names[i] = "asset"+i; transfer[i] = names[i] + "<-" + instrum[i];}



     assetString = "c(";
     closes = "closes<-na.locf(cbind(";
     name_vector = "c("; 
   
   //Get Close data -------------
   
     na = "(cbind("; name_vector = "c("; 
     //if(logreturns) {name_vector = name_vector + "'asset0',";}

     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Cl("+names[i] + "),";
        name_vector = name_vector + "'" + names[i] + "',";
        assetString = assetString + "'"+instrum[i] + "',"; 
     }
     na = na + "Cl("+names[n_inst-1]+")))";
     assetString = assetString + "'" + instrum[n_inst-1] + "')";    
     name_vector = name_vector + "'" + names[n_inst-1] + "')"; 
     closes = "closes <- na.locf" + na;
     name_vector = "names(closes) <- " + name_vector;

   //Get High and Low data-------
     na = "(cbind("; name_vector3 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Hi("+names[i] + "),";
        name_vector3 = name_vector3 + "'assetHi" +i+"',";
     }
     na = na + "Hi("+names[n_inst-1]+")))";
     name_vector3 = name_vector3 + "'assetHi"+(n_inst-1)+"')"; 
     highs = "highs <- na.locf" + na;
     name_vector3 = "names(highs) <- " + name_vector3;  
    
    
     na = "(cbind("; name_vector4 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Lo("+names[i] + "),";
        name_vector4 = name_vector4 + "'assetLo" +i+"',";
     }
     na = na + "Lo("+names[n_inst-1]+")))";
     name_vector4 = name_vector4 + "'assetLo"+(n_inst-1)+"')"; 
     lows = "lows <- na.locf" + na;
     name_vector4 = "names(lows) <- " + name_vector4; 
     
     
 
     //Get volume
     volume = "volumes<-na.locf(cbind(";
     name_vector2 = "c("; instString2 = "";
     name_vector5 = "c("; 
     for(i=0;i<n_inst-1;i++)
     {
       instString2 = instString2 + "Vo(" + names[i] + "),";
       name_vector2 = name_vector2 + "'assetv"+i+"',";
       name_vector5 = name_vector5 + "'assetPrice"+i+"',";
     } 
     instString2 = instString2 + "Vo(" + names[n_inst-1] + ")))";
     name_vector2 = name_vector2 + "'assetv"+(n_inst-1)+"')";
     name_vector5 = name_vector5 + "'assetPrice"+(n_inst-1)+"')";
     volume = volume+instString2;
     name_vector2 = "names(volumes) <- " + name_vector2;
     name_vector5 = "names(orig)<- " + name_vector5;

     


     try
     {

      RCaller caller = new RCaller();
      caller.setRscriptExecutable("/usr/bin/Rscript");
      caller.cleanRCode();
      caller.getRCode().addRCode("require (Runiversal)");
      caller.getRCode().addRCode("require (FinancialInstrument)");                    
      caller.getRCode().addRCode("setDefaults(getSymbols,verbose=TRUE,src='yahoo')");
      //caller.getRCode().addRCode("loadInstruments('/home/lisztian/Data/instruments.rda')");  
      //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/Data/sec',use_identifier='X.RIC',extension='RData')");
 
      //----- Call the symbols -------------------------
      addInst = "getSymbols(" + assetString + ",from='" + date1 + "')"; 
      caller.getRCode().addRCode(addInst);


      for(i=0;i<n_inst;i++)
      {caller.getRCode().addRCode(transfer[i]);}
      

       //----- define the matrix ------------------------
      caller.getRCode().addRCode(closes);
      caller.getRCode().addRCode(highs);
      caller.getRCode().addRCode(lows);
      caller.getRCode().addRCode(volume);

      //----- define the name_vector--------------------
      caller.getRCode().addRCode(name_vector);
      caller.getRCode().addRCode(name_vector3);
      caller.getRCode().addRCode(name_vector4);      
      caller.getRCode().addRCode(name_vector2);


      //--- somehow, get original series, but only first one, we target first in list
      //caller.getRCode().addRCode("times <- index(closes)");
      caller.getRCode().addRCode("closes <- log(closes)");
      caller.getRCode().addRCode("highs<-log(highs)");
      caller.getRCode().addRCode("lows<-log(lows)");
      caller.getRCode().addRCode("highs<-highs-lows");
      
      caller.getRCode().addRCode("orig <- closes");
      caller.getRCode().addRCode(name_vector5); 
      caller.getRCode().addRCode("orig[is.na(orig)]<-0");
      caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)"); 
      
      caller.getRCode().addRCode("closes <- diff(closes)");
      caller.getRCode().addRCode("volumes <- diff(volumes)");

      if(high_low_diff)
      {caller.getRCode().addRCode("highs<-diff(highs)"); caller.getRCode().addRCode("highs[1,] <- 0");}
      
      caller.getRCode().addRCode("closes[1,] <- 0");
      caller.getRCode().addRCode("volumes[1,] <- 0");
      caller.getRCode().addRCode("closes[is.na(closes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("volumes[is.na(volumes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("highs[is.na(highs)]<-0");  //--- set any final NAs to 0 
      
      
      caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
      caller.getRCode().addRCode("volume_list<-lapply(as.list(volumes), coredata)"); 
      caller.getRCode().addRCode("highs_list<-lapply(as.list(highs), coredata)");
      
      caller.getRCode().addRCode("all_list<-c(orig_list,series_list,highs_list,volume_list)");


      caller.runAndReturnResult("all_list");


    
      
      real_series = caller.getParser().getAsDoubleArray("assetPrice0");
      if(clear)
      { 
       setN_obs(real_series.length);
       clearSeries(); 
       setNobs(getN_obs());
       deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       getTheatre().setNobs(getN_obs());
       getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
       timeScrollPane.setViewportView(getTheatre());
      }       
      
      //Get Prices first 
      for(i=0;i<n_inst;i++)
      {
         real_series = caller.getParser().getAsDoubleArray("assetPrice"+i);
         setN_obs(real_series.length); 
         addSeries(6);
      }

      //Get log-retuns
      if(logreturns)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("asset"+i);
         addSeries(6);
       }
      }
      
      if(high_low || high_low_diff)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetHi"+i);
         addSeries(6);
       } 
      }
      
      if(vol)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetv"+i);
         addSeries(6);
       }
      }
    }
    catch(Exception e){System.out.println(e.toString());}

   }   



//-------------- Computes the spread between two log(price) series ------
// log(p) - b log(q) = log(p/c q) 
//-----------------------------------------------------------------------

public void computeSpread(double[] x, double[] y)
{
     
  if(x.length == getN_obs())
  {     
        
      try
      {

       RCaller caller = new RCaller();
       caller.setRscriptExecutable("/usr/bin/Rscript");
       caller.cleanRCode();
       
       //---- call libraries --------------
       caller.getRCode().addRCode("require (Runiversal)");   
//        caller.addDoubleArray("x", x); //set series
//        caller.addDoubleArray("y", y);
       caller.getRCode().addRCode("r<-princomp(~x+y)"); 
       caller.getRCode().addRCode("beta<-r$loadings[1,1]/r$loadings[2,1]");
       caller.getRCode().addRCode("spread<-x - beta*y");
       caller.getRCode().addRCode("my.all<-list(xy_spread=spread,beta=beta)");
       
       caller.runAndReturnResult("my.all");
       real_series = caller.getParser().getAsDoubleArray("xy_spread");
       double[] _coint_coeff = caller.getParser().getAsDoubleArray("beta");  
     
       coint_coeff = _coint_coeff[0];
       addSeries(6);
      }
      catch(Exception e){System.out.println(e.toString());}
  }
}


public void computeMeanStdSpread(double[] x, double[] y, double beta)
{
  int i,n;
  double sum = 0.0;

  n = x.length; 
  double[] spread = new double[n];

  for(i=0;i<n;i++)
  {sum = sum + (x[i] - beta*y[i])*(x[i] - beta*y[i]); spread[i] = x[i] - beta*y[i];}

  sum = Math.sqrt(sum)/n;
  spread_mean = sum;

  sum = 0.0;   
  for(i=0;i<n;i++)
  {sum = sum + (spread[i] - spread_mean)*(spread[i] - spread_mean);}
  spread_std = Math.sqrt(sum)/n;
}


/*-------------------------------------------------------

   The idea here is to consider the first two price series 
   in our data set. If there is less than two, does nothing.
   computes the cointegration constant, the resulting spread mean and standard dev
---------------------------------------------------------*/

public void computePriceSeriesSpread()
{
  int i;
  int nprice = 0; 
  int[] price_loc = new int[2]; 
    
  i = 0;
  while(i < getN_sim_series() && nprice < 2)
  {
    if(priceCheck[i].isSelected()) {price_loc[nprice] = i; nprice++;}    
    i++;
  }  

  if(nprice >= 2)
  {
    double[] tempx = getSim_data().get(price_loc[0]);
    double[] tempy = getSim_data().get(price_loc[1]);
    computeSpread(tempx, tempy);
    computeMeanStdSpread(tempx, tempy, coint_coeff); 
    spread_exists = true;
  }
  
}


  //----- Gaussianize the target data, also creates volatility estimate 
  
//   public void removeOutliers(int rank, double scale)
//   {
//      int i;
//      outlier_scale = scale;
//      outlier_rank = rank;
//   
//      
//      for(i=0;i<n_sim_series;i++)
//      {
//        double[] temp;
//        if(repCheck[i].isSelected())
//        {
//          temp = outlierDetection(sim_orig.get(i)); 
//          sim_data.set(i,temp); theatre.setSeries(i,temp);
// 
//          if(mix[i])
//          {           
//           target_series = mixSeries(target_series, sim_data.get(i), weight_mix[i]);        
//          }
//        }
//      }  
//      theatre.setTarget(target_series);   
//      plotTarget.setEnabled(true);              
//   }
  

//   public double[] outlierDetection(double[] g)
//   {
//        int i;
//        double[] newgdp = new double[g.length];
//        System.arraycopy(g, 0, newgdp, 0, g.length);
//        int[] pindex = new int[g.length];
//        
//        Cronos.sortsims(g.length, newgdp, pindex);
//         
//        double threshold = Math.abs(newgdp[outlier_rank]);
//        
//        for(i=0;i<g.length;i++)
//        {
//          if(Math.abs(g[i]) >= threshold)
//          {
//            g[i] = g[i]*outlier_scale;
//          }
//        }
//        return g;
//   }  
  
  
  
  public void gaussianizeData(double phi1, double phi2)
  {
     int i; double[] rt = new double[getN_obs()];;    
     boolean found = false;
     diff_vol = true;
     double vol_mean = 0;
     for(i = 0; i < getN_sim_series(); i++)
     {
        if(repCheck[i].isSelected() && simCheck[i].isSelected())
        {
          rt = getSim_data().get(i);
          found = true;
          break;        
        }
     }  
          
     if(found)
     {
      int le = getN_obs(); 
      double[] extend_burn = new double[2*le-1];
      double[] sigma = new double[2*le];
      double[] vol = new double[le];
      double mean = 0;
     
      double mu = 2.0;
     
      vol[0] = 0;
      //-- mirror the data
      for(i=0;i<le-1;i++)
      {
       extend_burn[le - 2 - i] = rt[i+1];
       extend_burn[le-1+i] = rt[i];
       
       mean = mean + rt[i];
      }
      mean = mean/(double)le;
     
      extend_burn[2*le-2] = rt[le-1];
     
      System.out.println("mu = " + mu + ", phi1 = " + phi1 + " " + phi2);
     
      sigma[0] = mean;
      for(i=0;i<2*le-1;i++)
      {
       sigma[i+1] = mu*(1.0 - phi1 - phi2) + phi1*sigma[i] + phi2*extend_burn[i]*extend_burn[i];
      } 
   
      for(i=0;i<le;i++)
      {rt[i] = rt[i]/Math.sqrt(sigma[le-1+i]); vol[i] = Math.log(sigma[le-1+i]); vol_mean = vol_mean + vol[i];}
  
      vol_mean = vol_mean/le;
      for(i=0;i<le;i++)
      {vol[i] = vol[i] - vol_mean;}
  
//       if(diff_vol)
//       { 
//         for(i=1;i<le;i++)
//         {vold[i] = 10000*(vol[i] - vol[i-1]);}
//         vold[0] = vold[1];
//       }
  
      
//       if(gaussian_data)
//       {
//         //sim_data.remove(sim_data.size()-1); sim_means.remove(sim_data.size()-1); sim_orig.remove(sim_data.size()-1);
//         //sim_data.remove(sim_data.size()-1); sim_means.remove(sim_data.size()-1); sim_orig.remove(sim_data.size()-1);
//       
//         deleteSeries(sim_data.size()-1);
//         deleteSeries(sim_data.size()-1);            
//       }
  
      real_series = rt;      
      addSeries(6);
      
//       if(!diff_vol)
     
      
//       else
//       {real_series = vold;}
      if(!vol_set)       
      {
       real_series = vol;
       addSeries(6);
       vol_set = true;
      }
      
      gaussian_data = true;
      which_gaussian = getSim_data().size()-2;
      which_vol = getSim_data().size()-1;
         
     }
  
  }



    public double getMin(double[] temp)
    {
       int i;
       double min = 0;

       for(i = 0;i<temp.length;i++)
       {
         if(temp[i] < min) {min = temp[i];}
       }
       return min;
    }

    public void readData(File file)
    {
          
       String strline; Double D; int j,k;
       String[] tokens; String delims = "[ ]+";
       int n_toks;
       double[] values = new double[3000];      
       double val = 0;
       int i = 0;
       int count=0;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
        
           //System.out.println(tokens[0] + " " + tokens[1]);
         
           n_toks = tokens.length; n_rep = n_toks;
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           for(k=0;k<n_rep;k++)
           {
  
            D = new Double(tokens[k]); //take only the value, no dates
            val = D.doubleValue();

            values[i] = val;
            i++;
           }
           
           count++;
         } 
         
            if(count != getN_obs())
            {clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;}
           
            setN_obs(count);  
            
              
            
            for(k=0;k<n_rep; k++)
            {
              real_series = new double[getN_obs()];
            
              for(j=0;j<getN_obs();j++)
              {
                 real_series[j] = values[n_rep*j+k];
              }               
              addSeries(6); 
/*              double[] tempval = cumsum(real_series,real_series.length);
              real_series = tempval;
              addSeries(6);    */                      
            }
            
            
               
            
           din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }

   
   public void readTBillData(File file)
   {
   
     String strline;
     String[] tokens,dates;
     int i = 0;
     int count=0;
     ArrayList<String> _dates = new ArrayList<String>();
     ArrayList<Double> tbill3,tbill6,loan_demand,gdp,M1,inflation,dollar,fed_funds,indust,net_saving,vixx;
     

     tbill3 = new ArrayList<Double>(); tbill6 = new ArrayList<Double>(); loan_demand = new ArrayList<Double>();
     inflation = new ArrayList<Double>(); dollar = new ArrayList<Double>(); fed_funds = new ArrayList<Double>();
     indust = new ArrayList<Double>(); net_saving = new ArrayList<Double>(); vixx = new ArrayList<Double>();
     M1 = new ArrayList<Double>(); gdp = new ArrayList<Double>();
   
     try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {   
     
            tokens = strline.split("[ ]+");
            _dates.add(tokens[0]);

            tbill3.add(new Double(tokens[1])); 
            tbill6.add(new Double(tokens[2])); 
            loan_demand.add(new Double(tokens[3])); 
            gdp.add(new Double(tokens[4]));
            M1.add(new Double(tokens[5])); 
            inflation.add(new Double(tokens[6])); 
            dollar.add(new Double(tokens[7])); 
            fed_funds.add(new Double(tokens[8]));
            indust.add(new Double(tokens[10])); 
            net_saving.add(new Double(tokens[11])); 
            vixx.add(new Double(tokens[12]));
            count++;
         }

         
         
         
         setN_obs(count);
         clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;
                  
         String[] names = {"t3bill","dt3bill","loanD","GDP","M1","Infl","Dollar","fedFund","Indust","Saving","VXX"};
         
         
         //---------- Add target series ------------------------------------
         dates = new String[getN_obs()];
         double[] dtbill = new double[getN_obs()]; dtbill[0] = 0;
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {
            dates[i] = _dates.get(getN_obs() - 1 - i);
            real_series[i] = tbill3.get(getN_obs() - 1 - i);
            
            if(i > 0) {dtbill[i] = (real_series[i] - real_series[i-1])/real_series[i-1];}
         }
         addSeries(6);         
         
         //---------- Add differenced target series ------------------------------------
         real_series = new double[getN_obs()];
         System.arraycopy(dtbill,0,real_series,0,getN_obs());
         addSeries(6);

         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = loan_demand.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = gdp.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = M1.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = inflation.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = dollar.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = fed_funds.get(getN_obs() - 1 - i);}
         addSeries(6);
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = indust.get(getN_obs() - 1 - i);}
         addSeries(6);         
        
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = net_saving.get(getN_obs() - 1 - i);}
         addSeries(6);    
         
         //--------- Add exp series ----------------------------------------------------
         real_series = new double[getN_obs()];
         for(i=0;i<getN_obs();i++)
         {real_series[i] = vixx.get(getN_obs() - 1 - i);}
         addSeries(6);             
         
         getTheatre().setNames(names);
         getTheatre().setDates(dates);     
    
    
    
         din.close();
    }
    catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
  
   }
   
   
   
   
   
   public void readPortfolioData(File file)
   {
   
     String strline;
     String[] names;
     int i = 0;
     int count=0;
     ArrayList<String> _dates = new ArrayList<String>();
     boolean dateF = false;
     names = new String[12];   
     String[] tokens; 
     System.out.println(file.getParent());  
     String path = file.getParent();
       
       
     try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {
           names[count] = strline;              
           count++;
         }  
     
       din.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
       //now open up all files 
       
       for(i = 0; i < count; i++)
       {        
         File nfile = new File("" + path + "/" + names[i]);
         this.readReturnData(nfile);
       }
       
       getTheatre().setNames(names);
       
       File datefile = new File("" + path + "/dates.dat");
       if(datefile.exists())
       {dateF = true;}
       else
       {dateF = false; datefile = new File("" + path + "/" + names[0]);}
       
       try{
          
         FileInputStream fin = new FileInputStream(datefile);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {
           
           if(dateF){ _dates.add(strline);}
           else
           {
             tokens = strline.split("[,]+");
             _dates.add(tokens[0]);
           }
         }  
       
         getTheatre().setDates(_dates.toArray(new String[0]));
       din.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}              
       
       
       
    }
       
   
//    public void printPortfolio()
//    {
//    
//       int i;
//       out.print("Date    ");
//       out.print("Cumulative Returns ");
//       for(i = 0; i < theatre.asset_names.length-1; i++)
//       {
//         out.print(
//    
//    
//    
//    
//    
//    }
   
   
    public void readReturnData(File file)
    {
          
       String strline; Double D; int j;
       String[] tokens; String delims = "[ ]+"; 
       String[] portn;
       String delims2 = "[.]+";
       int n_toks;
       double[] values = new double[600];      
       double val = 0;
       int i = 0;
       int count=0;
       double[] temp = new double[1];
       sigex_returns_on = true;
       cumsum_on = true;
       double mult = 1;
       String fname = file.getName();
       System.out.println("File name is " + fname);
       System.out.println("File path is " + file.getParent());
       portn = fname.split(delims2);
       
       if(portn[0].indexOf("jpy") != -1) {mult = .01;}
       else {mult = 1;}
       
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
        
           //System.out.println(tokens[0] + " " + tokens[1]);
         
           n_toks = tokens.length; n_rep = n_toks;
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           if(n_toks == 2)
           {
            D = new Double(tokens[1]); //take the second value
            val = D.doubleValue();

            values[i] = mult*val;
            i++;
           }
           else
           {
           
            D = new Double(tokens[0]); //take only the value, no dates
            val = D.doubleValue();

            values[i] = val;
            i++;
           }
           
           count++;
         }
         

         
         if(first_in)  //first one in, we'll keep this as the basis for any others that come in
         {
           setN_obs(count);
           clearSeries(); deleteAll.setEnabled(false); acfplotBox.setEnabled(false); getTheatre().setNobs(count); simData = false;
           first_in = false;
         }

         if(count == getN_obs())
         {
          temp = new double[count];
          System.arraycopy(values,0,temp,0,count);
          returnsAll.add(temp);         
          //double[] tempval = cumsum(temp,count);         
         }
         else if(count < getN_obs())
         {
           int descrep = getN_obs() - count;
           //just add zeros
           temp = new double[getN_obs()];
           for(i = 0; i < descrep; i++) {temp[i] = 0;}
           for(i=descrep;i<getN_obs();i++) {temp[i] = values[i-descrep];}
           returnsAll.add(temp);         
           //double[] tempval = cumsum(temp,count);
         }
         else if(count > getN_obs())
         {
           int descrep = count - getN_obs();
           temp = new double[getN_obs()];
           for(i=0;i<getN_obs();i++) {temp[i] = values[i+descrep];}
           returnsAll.add(temp);  
         }
         computePortfolio();   
         
         real_series = new double[getN_obs()];
         for(j=0;j<getN_obs();j++)
         {real_series[j] = temp[j];} 
         addSeries(6); 
         n_rep = returnsAll.size();
         
                  
         weightLabel[returnsAll.size()-1].setText(portn[0]);
 
         din.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

}

public void getIntradayGoogleData(String[] instrum, String[] exchange, String interval, String noDays, boolean clear, boolean logreturns, boolean vol, 
                                  boolean high_low, boolean high_low_diff, boolean high_freq, int period)
{
    //%Y-%m-%d %H:%M:%S

     int i;
     int n_inst = instrum.length;
     String[] names = new String[n_inst];
     String na;
     String closes,volume;
     String name_vector,name_vector2,name_vector3,name_vector4,name_vector5,highs,lows;
     String assetString, instString2;
     //---- create matrix of closes------------------------
  

     assetString = "c(";
     closes = "closes<-na.locf(cbind(";
     name_vector = "c("; 
   
   //Get Close data -------------
   
     na = "(cbind("; name_vector = "c("; 
     //if(logreturns) {name_vector = name_vector + "'asset0',";}

     for(i=0; i < n_inst-1; i++)
     {     
        names[i] = new String("asset"+i);
        //---gets closing data ------
        na = na + "Cl("+names[i] + "),";
        name_vector = name_vector + "'" + names[i] + "',";
        assetString = assetString + "'"+instrum[i] + "',"; 
     }
     names[n_inst-1] = new String("asset"+(n_inst-1));
     na = na + "Cl("+names[n_inst-1]+")))";
     assetString = assetString + "'" + instrum[n_inst-1] + "')";    
     name_vector = name_vector + "'" + names[n_inst-1] + "')"; 
     closes = "closes <- na.locf" + na;
     name_vector = "names(closes) <- " + name_vector;

   //Get High and Low data-------
     na = "(cbind("; name_vector3 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Hi("+names[i] + "),";
        name_vector3 = name_vector3 + "'assetHi" +i+"',";
     }
     na = na + "Hi("+names[n_inst-1]+")))";
     name_vector3 = name_vector3 + "'assetHi"+(n_inst-1)+"')"; 
     highs = "highs <- na.locf" + na;
     name_vector3 = "names(highs) <- " + name_vector3;  
    
    
     na = "(cbind("; name_vector4 = "c("; 
     for(i=0; i < n_inst-1; i++)
     {     
        //---gets closing data ------
        na = na + "Lo("+names[i] + "),";
        name_vector4 = name_vector4 + "'assetLo" +i+"',";
     }
     na = na + "Lo("+names[n_inst-1]+")))";
     name_vector4 = name_vector4 + "'assetLo"+(n_inst-1)+"')"; 
     lows = "lows <- na.locf" + na;
     name_vector4 = "names(lows) <- " + name_vector4; 
     

     Integer hfint = new Integer(interval);
     int hf_interval = hfint.intValue()/period; 
 
     //Get volume
     volume = "volumes<-na.locf(cbind(";
     name_vector2 = "c("; instString2 = "";
     name_vector5 = "c("; 
     for(i=0;i<n_inst-1;i++)
     {
       instString2 = instString2 + "Vo(" + names[i] + "),";
       name_vector2 = name_vector2 + "'assetv"+i+"',";
       name_vector5 = name_vector5 + "'assetPrice"+i+"',";
     } 
     instString2 = instString2 + "Vo(" + names[n_inst-1] + ")))";
     name_vector2 = name_vector2 + "'assetv"+(n_inst-1)+"')";
     name_vector5 = name_vector5 + "'assetPrice"+(n_inst-1)+"')";
     volume = volume+instString2;
     name_vector2 = "names(volumes) <- " + name_vector2;
     name_vector5 = "names(orig)<- " + name_vector5;

     


     try
     {

      RCaller caller = new RCaller();
      caller.setRscriptExecutable("/usr/bin/Rscript");
      caller.cleanRCode();
      caller.getRCode().addRCode("require (Runiversal)");
      caller.getRCode().addRCode("require (FinancialInstrument)");                    
      //caller.getRCode().addRCode("setDefaults(getSymbols,verbose=TRUE,src='yahoo')");
      //caller.getRCode().addRCode("loadInstruments('/home/lisztian/Data/instruments.rda')");  
      //caller.getRCode().addRCode("setSymbolLookup.FI('/home/lisztian/Data/sec',use_identifier='X.RIC',extension='RData')");
 
      caller.getRCode().addRCode("source('/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/intradayGoogle.r')");
 
      //----- Call the symbols -------------------------

      for(i=0;i<n_inst;i++)
      {
        caller.getRCode().addRCode("asset"+i+"<-intradayGoogle('"+instrum[i]+"','"+exchange[i]+"',"+interval+","+noDays+")");
        caller.getRCode().addRCode("asset"+i+"<-asset"+i+"['T09:30/T15:59',]");  
      }
      
      if(high_freq)
      {
        caller.getRCode().addRCode("hfasset<-intradayGoogle('"+instrum[0]+"','"+exchange[0]+"',"+hf_interval+","+noDays+")");
        caller.getRCode().addRCode("hfasset<-log(na.locf(Cl(hfasset)))['T09:30/T15:59',]");
        caller.getRCode().addRCode("names(hfasset)<-'hfasset'");
      }
      
      
       //----- define the matrix ------------------------
      caller.getRCode().addRCode(closes);
      caller.getRCode().addRCode(highs);
      caller.getRCode().addRCode(lows);
      caller.getRCode().addRCode(volume);

      //----- define the name_vector--------------------
      caller.getRCode().addRCode(name_vector);
      caller.getRCode().addRCode(name_vector3);
      caller.getRCode().addRCode(name_vector4);      
      caller.getRCode().addRCode(name_vector2);


      //--- somehow, get original series, but only first one, we target first in list
      //caller.getRCode().addRCode("times <- index(closes)");
      caller.getRCode().addRCode("closes <- log(closes)");
      caller.getRCode().addRCode("highs<-log(highs)");
      caller.getRCode().addRCode("lows<-log(lows)");
      caller.getRCode().addRCode("highs<-highs-lows");
      

      
      caller.getRCode().addRCode("orig <- closes");
      caller.getRCode().addRCode(name_vector5); 
      caller.getRCode().addRCode("orig[is.na(orig)]<-0");
      caller.getRCode().addRCode("orig_list<-lapply(as.list(orig), coredata)"); 
      
      caller.getRCode().addRCode("closes <- diff(closes)");
      caller.getRCode().addRCode("volumes <- diff(volumes)");

      if(high_low_diff)
      {caller.getRCode().addRCode("highs<-diff(highs)"); caller.getRCode().addRCode("highs[1,] <- 0");}
      
      caller.getRCode().addRCode("closes[1,] <- 0");
      caller.getRCode().addRCode("volumes[1,] <- 0");
      caller.getRCode().addRCode("closes[is.na(closes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("volumes[is.na(volumes)]<-0");  //--- set any final NAs to 0
      caller.getRCode().addRCode("highs[is.na(highs)]<-0");  //--- set any final NAs to 0 
      
      
      caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
      caller.getRCode().addRCode("volume_list<-lapply(as.list(volumes), coredata)"); 
      caller.getRCode().addRCode("highs_list<-lapply(as.list(highs), coredata)");

      if(high_freq)
      {
        caller.getRCode().addRCode("hf_list<-lapply(as.list(hfasset), coredata)");
        caller.getRCode().addRCode("all_list<-c(orig_list,series_list,highs_list,volume_list,hf_list)");
      }
      else
      {caller.getRCode().addRCode("all_list<-c(orig_list,series_list,highs_list,volume_list)");}
      
      caller.runAndReturnResult("all_list");


      real_series = caller.getParser().getAsDoubleArray("assetPrice0");
      if(clear)
      { 
       setN_obs(real_series.length);
       clearSeries(); 
       setNobs(getN_obs());
       deleteAll.setEnabled(false); 
       acfplotBox.setEnabled(false); 
       getTheatre().setNobs(getN_obs());
       getTheatre().setPreferredSize(new Dimension(2*getN_obs(), 300));
       timeScrollPane.setViewportView(getTheatre());
      }       
      
      //Get Prices first 
      for(i=0;i<n_inst;i++)
      {
         real_series = caller.getParser().getAsDoubleArray("assetPrice"+i);
         setN_obs(real_series.length); 
         addSeries(6);
      }

      //Get log-retuns
      if(logreturns)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("asset"+i);
         addSeries(6);
       }
      }
      
      if(high_low || high_low_diff)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetHi"+i);
         addSeries(6);
       } 
      }
      
      if(vol)
      {
       for(i=0;i<n_inst;i++)
       {
         real_series = caller.getParser().getAsDoubleArray("assetv"+i);
         addSeries(6);
       }
      }
      
      if(high_freq)
      {
        hf_data = caller.getParser().getAsDoubleArray("hfasset");
        System.out.println("Size of original = " + real_series.length + ", size of hf = " + hf_data.length);
      }
    }
    catch(Exception e){System.out.println(e.toString());}

}      
   
  

   
   
   
   /*
   public void writeX13ParameterFile(int[] sdim, boolean ldefinemodel, boolean ltrans, boolean loutlier, 
                                                 boolean lregression, boolean lsigExdiagnostics, boolean seasonal 
                                                 boolean td, boolean easter, int ea ) 
   {

fileName<-'/tmp/asset.csv'   
download.file(paste("http://www.google.com/finance/getprices?q=",symbol,"&x=",exchange,"&i=",interval,"&p=",noDays,"d&f=d,o,h,l,c,v,t",sep=""), fileName)  
   
unix2POSIXct <- function (time)  structure(time, class = c("POSIXt", "POSIXct"))  
data = read.table(fileName,sep=",",col.names=c("DATE1","CLOSE","HIGH","LOW","OPEN","VOLUME"),fill=TRUE)  
data$DATE = 0   
 
 
for (i in 8:nrow(data))  
 {  
  if(i==8 || substr(as.vector((data$DATE1[i])),1,1) == "a")  
  {  
   tempDate = unix2POSIXct(as.numeric(substr(as.vector((data$DATE1[i])),2,nchar(as.vector((data$DATE1[i]))))))   
   data$DATE[i] = tempDate   
  } else {  
   tempDate1 = tempDate + as.numeric(as.vector(data$DATE1[i]))*interval   
   data$DATE[i] = unix2POSIXct(tempDate1) 
  }   
 }  
data1=as.data.frame(data)  
   
finalData = data.frame(DATE=unix2POSIXct(as.character(data1$DATE[8:nrow(data)])),CLOSE=data1$CLOSE[8:nrow(data)],HIGH=data1$HIGH[8:nrow(data)],LOW=data1$LOW[8:nrow(data)],OPEN=data1$OPEN[8:nrow(data)],VOLUME=data1$VOLUME[8:nrow(data)])  
write.csv(finalData,file=fileName,row.names=FALSE)     
   
asset1<-as.xts(read.zoo(fileName, sep = "," ,header = TRUE, tz = ""))  
   
   
      String check = "check{savelog = (sft)}\n";
      String x11 = "x11{savelog = all}\n";
      String regression = "";
      String outlier = "outlier{types= (AO LS) lsrun = 0}

      String variables = "";
 
      if(seasonal) {variables = " seasonal ";}

      if(td && easter) {
         regression = regression + "regression{ variables = ( " + variables + " td    easter["+ea+"]) \n";
         regression = regression + "            aictest = (td easter) savelog = aictest     \n";        
      }
      else if(td) {
         regression = regression + "regression{ variables = td \n";
         regression = regression + "            aictest = td savelog = aictest     \n";        
      }      
      else if(easter) {
         regression = regression + "regression{ variables = easter["+ea+"] \n";
         regression = regression + "            aictest = easter savelog = aictest     \n";        
      }
      
     try{  
           PrintWriter out = new PrintWriter(new FileWriter("x13parameters.spc"));

           out.println("series { title = \"Series\" start = 1976.jan file = \"useless.dat\"} ");
           if(ldefinemodel) out.println("arima { model = ("+sdim[0]+ " " +sdim[1] + " " + sdim[2]+")("+sdim[3]+ " " +sdim[4] + " " + sdim[5]+")12 }");
           if(ltrans) out.println("transform{ function = auto }");
           if(loutlier) out.println("outlier{  }");
           if(lregression) {out.println(regression);}
           if(lsigExdiagnostics) {out.println("seats { noadmiss = yes }");}
           
           out.close();
        } catch (IOException e) {e.printStackTrace();}   
   
   
   
   }


  */

    
    public static void main(String args[])
    {

      System.loadLibrary("simSarima");
      
      int n_obs; int burnin; int n_rep; 
     
      n_obs = 100; burnin = 100; n_rep = 5; n_obs = n_obs+200;


      SimPanel sim = new SimPanel(n_obs, n_rep, burnin, null);


      //File file = new File("/home/lisztian/Dropbox/iMetrica2013/uSim2012/x13/HFData/402.csv");
      sim.getHF_CSV_Data(15, 20, 50, false, false, false);
   

 


   }
   
   
   public void maximizeSharpe()
   {
   
      int i,j; 
      int n_basket = returnsAll.size();      
      
      if(n_basket > 0)
      {
       double[][] data = new double[getN_obs()][n_basket];
       double[] means = new double[n_basket];
       
       for(i=0;i<n_basket;i++)
       {
        double[] ret = returnsAll.get(i);
        double[] mstd = mean_std(ret); 
        means[i] = mstd[0];
        
        for(j=0;j<getN_obs();j++)
        {data[j][i] = ret[j];}        
       }
       
       RealVector m = new ArrayRealVector(means, false);
       Covariance covComp = new Covariance(data);
       
       DecompositionSolver solver = new LUDecomposition(covComp.getCovarianceMatrix()).getSolver();
       RealVector sol = solver.solve(m);
       
       double[] w = sol.toArray(); 
       double sumw = 0;
       for(i=0;i<w.length;i++) 
       {
        if(w[i] < 0) 
        {w[i] = 0;} 
        sumw = sumw + w[i]; 
       }
       
       for(i=0;i<w.length;i++) 
       {w[i] = w[i]/sumw;}       
       
       for(i = 0; i < w.length; i++)
       {weightSliders[i].setValue((int)(100*(w[i])));}             
       
      } 
   }

public int getN_sim_series() {
	return n_sim_series;
}

public void setN_sim_series(int n_sim_series) {
	this.n_sim_series = n_sim_series;
}

public SimCanvas getTheatre() {
	return theatre;
}

public void setTheatre(SimCanvas theatre) {
	this.theatre = theatre;
}

public ArrayList<double[]> getSim_data() {
	return sim_data;
}

public void setSim_data(ArrayList<double[]> sim_data) {
	this.sim_data = sim_data;
}

public boolean[] getRepresent() {
	return represent;
}

public void setRepresent(boolean[] represent) {
	this.represent = represent;
}

public double[] getTarget_series() {
	return target_series;
}

public void setTarget_series(double[] target_series) {
	this.target_series = target_series;
}

public int getN_obs() {
	return n_obs;
}

public void setN_obs(int n_obs) {
	this.n_obs = n_obs;
}
   

} 
