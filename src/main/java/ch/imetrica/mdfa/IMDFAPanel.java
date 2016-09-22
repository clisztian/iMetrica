package ch.imetrica.mdfa;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.JCheckBox;
import javax.swing.border.*;
import javax.swing.border.LineBorder;
import java.awt.event.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.BorderFactory; 
import javax.swing.border.TitledBorder;
import java.text.*;

import ch.imetrica.bayesCronos.Cronos;
import ch.imetrica.usimx13.Polynomial;
import ch.imetrica.usimx13.SARIMAmodelJava;
import ch.imetrica.usimx13.SpectralDensity;

import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import rcaller.RCaller;
/*----------------------------------------------------------------
 
   Panel for IMDFA Controls - 
     Simulation:  
        Two different types of simulations 
          - DFA uMath.sing uSimX13 data 
          - MDFA uMath.sing simulated GDP series of Wildi 
   

----------------------------------------------------------------*/


public class IMDFAPanel extends JPanel
{


      /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--------- main controls ------------------------------------
	
	  public IMDFAcanvas mdfa_canvas;      // ---------- the mdfa time domain canvas
      public IMDFA mdfa;                  // ---------- the main mdfa engine
  
      FilterCanvas filter_canvas;   // ---------- filter domain canvas
      PeriodCanvas period_canvas;   // ---------- periodogram canvas
      FilterCoefCanvas coef_canvas; // ---------- filter coefficient canvas
      ZPCFilter zpc; 
      ZPCCanvas zpcFreqCanvas;
      private accountCanvas account_canvas;
      SpectralDensityCanvas specDensCanvas;
      SpectralDensity specDens;
      
      //---- Data and filter controls -------------
      public int n_obs;
      public int n_rep;
      public int K;
      public int L;
      public int Lag;
      public int l1;
      public int nsamp;
      public int n_obs_fixed;
      public int n_rep_fixed;
      public int i1, i2, dd, DD, K1, output, flength,S;  
      public int n_tot;

      public int series_max=10; 
      int[] diffInd;                // ----------- indicator for differencing 0,1,2
      double[] expFunc;
      //---- extra Filter controls ----------------
      public double lambda, expweight, cutoff;
      public double smooth, decay, cross, decay2, onestep, onestep_diff;       
      public double cutoff0;
      public boolean reg;
      public boolean useSARIMA;  //--- Use the simulated data from uSimX13
      boolean simulate; //---- simulate IMDFA 
      boolean iter;    //----- iteration style direct filtering
      public int cutN,cutD;   //--- separate denominator and numerator
      public int shiftint;
      public double w2,w3;   //second set of 
 
      boolean return_strategy = false;
      Color[] filter_health;
      double[] sig_update;

      public CrystalBallcanvas crystal_ball;
      public boolean envision = false;
      boolean cPeriod;
      boolean zpc_gene = false;
      boolean diff_sig_trading = false;
 
      int seed;
      JTabbedPane mdfaPlotPane; //Divides pane in time/frequency
      JButton compute;
      JPanel timeGrid;

      public JPanel tradingoptimPanel;
      public JPanel tradingparameterPanel;
      public JPanel tradingstatPanel;
 
      //---------- Scrollbars -------------------------------------
      JScrollBar nObsBar, lambdaBar, lambda3Bar,smoothBar, decayBar, decay2Bar, crossBar,expBar, cutBarD, cutBarN, lagBar, LBar;

	public JScrollBar l1Bar;

	JScrollBar l2Bar;

	JScrollBar nrepBar;

	JScrollBar sampBar;

	JScrollBar deltaBar;

	JScrollBar onestepBar;
      JTextField nText, lambdaText, lambda3Text, smoothText, decayText, decay2Text, crossText, expText, cutTextD, cutTextN, cutText, textLag, LText,l1Text, l2Text, nrepText, sampText, onestepText;
      JTextField mdfaText, deltaText;
      Color myBlue;
      JScrollBar shiftBar; JTextField shiftText; JLabel shiftLabel, deltaLabel;
      JCheckBox updatei1Box,updatei2Box; 
       
      //---------- Radio buttons and CheckBoxes -----------------------
      JRadioButton gammaBut, amplitudeBut, timeDelayBut, periodBut;
      JCheckBox[] freqPlot;
      JCheckBox[] timePlot; 
      JCheckBox[] expVariable;
      JCheckBox plotXf, plotGamma;
      JCheckBox[] periodPlot; // x, xf, x_i
      JCheckBox[] coeffPlot;  //b_1...b_n      
      double[] gamma_zpc;
      double[] gamma_hybrid;
      
      int plotTF; 
      boolean data_uploaded;
      boolean setup;
      
      ButtonGroup freqGroup;
      ButtonGroup simGroup;
      JCheckBox i1Check, i2Check;
      JCheckBox dCheck, DCheck; 
      JCheckBox reComp; 
      JCheckBox iterCheck;
      JRadioButton[] filterPlot;
      ButtonGroup filterPlotGroup;

      JRadioButton simCheck, sarimaCheck, noneCheck;

      boolean computeFilter, reCompFilter;

      //--------- Some labels ----------------------------------
      JLabel nobsLabel,LLabel,lamLabel,expLabel,nrepLabel,diffLabel,cutLabel,sampLabel;
      JLabel lagLabel,cutNLabel,cutDLabel,simulateData, constLabel,l1Label,bandLabel;
      JLabel turnPoint, levelChange;
      JLabel differ, filterLength, mdfaLabel;
     
      //-------- some smaller panels --------------------------
      JPanel simulGrid,freqGrid,typefreqGrid, periodGrid; 
      DecimalFormat df,df2,df4;
      JScrollPane timeScrollPane;
      JScrollBar timeScrollBar;
      
      //------- Supplemental data from X13 -------
      double[] x13data;
      int x13nobs;
      int n_series;
 
      //------- extended series ------------------
      double[] extend;
      double[] symfilt;
      double[] bc;
      double[] fore_extend;
      int L1,L2;
      JCheckBox extFore;
      boolean ext_forecast;
      int n_sym;
      int n_out_samp;
      int prev_out_samp;
      int L_fixed;
  
      
      //------ Filter design controls ------------
      FilterDesign tfilter;
      JCheckBox automatic;
      JRadioButton lowBut, rampBut, bandBut;

	public JRadioButton highBut;
      JRadioButton periodsBut, omegaBut, fracBut;
      JRadioButton x13trend, x13ti, x13seas;
      JCheckBox x13filter;
      JScrollBar num1Bar, num2Bar;  
      JScrollBar den1Bar, den2Bar;
      JScrollBar omega0Bar, omega1Bar, omega2Bar, omega3Bar;
      JScrollBar period1Bar, period2Bar;
      ButtonGroup filterGroup; 
      ButtonGroup omegaGroup;
      boolean low,high,ramp,band;
      boolean autoComp;   

      boolean fixBand;          //---- to fix the bandwidth
      JCheckBox fixBandCheck;   //---- to fix the bandwidth

      double out_samp_error = 0.0;
      int p1,p2;      //------ period parameters  2pi/p1, 2pi/p2
      public double w0;   //------ cutoff omega1, omega2
	  public double w1;
      int num1, num2;  //----- numerators for cutoffs
      int den1, den2;  //----- denominators for cutoffs
      String eof;
 
      JTextField num1Text, num2Text;  
      JTextField den1Text, den2Text;
      JTextField omega0Text, omega1Text;
      JTextField omega2Text, omega3Text;
      JTextField period1Text, period2Text;

      //------------------- X13 filter stuff ----------------------
      double[] trendpoly;
      double[] seaspoly;
      double[] tipoly;
      double[] trnsdenom;
      double[] trnspoly;
      double[] Gamma_x13;
      double[] f_specDens;

      int ntrend, nseas, ntip;
      double trendvar, seasvar, tivar;
      Polynomial mPhi, mphi, mapoly;
      double m_innvar; 
      int x13model;
      boolean useX13filters; 


      //------------ cointegration stuff----------------

      JSlider[] cointSliders;
      //cointRangeModel brm = new cointRangeModel();
      JTextField cointText; 
      double[] w_const;
      double w_sum; 
      public double shift; 

      //---------- out of sample -------------------------
      JScrollBar outof_sample; 
      JTextField outof_sampleText;
      boolean b_outof_sample; int out; 
      ArrayList<double[]> complete_data;
      ArrayList<double[]> subset_data;

      //----Saving stuff-------------------------------------------  
      JCheckBox spec_densBox; 
      JCheckBox plotHist;
      String paramString;      //String of all the parameters in column form 
 
      ArrayList<String> histParams;  //---- should coincide with historical 
      ArrayList<MDFAFilter> savedFilters; //--- queue of older filters
      int n_hist;     

      String curDir;
      JFileChooser fc;
      double[] outSamp_means;   
      

     //-------------------- Stuff for ART display and TimeFrequency Map ------------
     public TimeFreqPlotPanel tfplotPanel;
     public IMFPlotPanel tcimfPanel;
     public artCanvas ARTCanvas;
     public JPanel sigex; 

     public double[] ARTfilter;
     public double[] freq_ints; 
     public double freq_start;
     public int n_div; 
     public int n_imfs; 
     
     private JPanel artPanel;
     private JPanel ARTpanel;

     private JLabel accLabel,freqIntsLabel,maxFreqLabel,mseLabel,reliableLabel,timeLabel;
     private JTextField accText,timeText,reliableText,mseText;       
     private JButton compARTFilter,computeTimeFreq;  
     private JCheckBox contARTCheck;
     private JPanel freqChangePanel;
     private JSlider freqIntsSlider,maxFreqSlider;
     private JPanel reliablePanel;
 
     private JLabel aerrorLabel,turningLabel;
     private JTextField aerrorText,turningText,phasedelayText,rerrorText;
     private JSlider turningSlider;
     private JLabel performLabel,phasedelayLabel,rerrorLabel;
     private JButton restartTimeButton,restartTimeButton1,restartTimeButton2;

     public double lambda_3;   //exponential value for controlling differencing 
     public double modz = .5; 
     public double modp = .5;
 
     private JRadioButton plotAMradioButton, plotIMFradioButton, plotnIMFradioButton;
     private JCheckBox plotDataCheckBox;
     private JPanel plotIMFPanel;

     //-------- zpc interface stuff
     JPanel zpcPanel; 
     JScrollBar alphaBar, hfBar, i1Bar, arg1Bar, modpBar, b0Bar, modzBar, lamzpcBar;
     JLabel amzpcLabel, alphaLabel, customLabel, argz1Label, i1label, zpcLabel, modzLabel;
     JLabel msezpcLabel, b0Label,  modpLabel, hfLabel,lamzpcLabel;
     JTextField alphaText, hfText, lamzpcText, i1Text, arg2Text1, arg1Text1, b0Text;
     JButton optimizeButton, injectButton, prefilterButton; 
     JCheckBox optimizeCheck;
     JTextField arg1Text; 
     JTextField msezpcText, modpText, modzText;
     boolean zpc_mod_optimize;
     JCheckBox[] plot_zpc;
     private int plot_number;
     JCheckBox i1const, i2const, normalizeConst;
     int p_arma,q_arma;
     JComboBox<String> zpcCombo; 


     //------ Time domain diagnostics ------------------------------
     double opt_lambda; double opt_alpha; 
     double w_1,w_2,w_3;
     double[] score;
     boolean timeDiagnostics;

    private JSlider accSlider,relSlider,timeSlider;
    private JLabel d_accLabel,d_dofLabel,d_icLabel,d_relLabel,d_timeLabel,optimLabel,d_mseLabel;
    private JTextField d_accText,d_dofText,d_icText,d_mseText,t_mseText,t_relText,t_excessText;
    private JTextField t_accText,t_corrText,d_relText,d_timeText,t_timeText;
    private JLabel freqdiagLabel,t_ARTLabel,t_corrLabel,t_accLabel,t_excessLabel;
    private JLabel t_mseLabel, t_relLabel, t_timeLabel, timediagLabel;
    private JButton optimizeARTFilter;
    private JPanel optimizeARTPanel,timediagPanel,freqdiagPanel;
    public JPanel diagnosticPanel;
    public JButton zpc_getData;
    public boolean autoCompZPC;
    public JCheckBox continuousUpdate;

    public JComboBox<String> rkhsCombo;
    public JComboBox<String> pBox;
    public JComboBox<String> qBox;
    public JCheckBox armaSD;
    public JCheckBox rkhsBox;
    


    //------------- Trading stuff-------------------

    public double[][] t_price;
    public double[] t_signal;   //original signal
    public double[] d_signal;   //derivative of signal
    public double[] account;
    public double[] target;
    
    public ArrayList<Double> adapt_signal;
    public double[] fixed_signal;   //fixed signal for recomputation
    public double[] fixed_account;  //fixed account for recomputation
    
    public int trade_obs;
    public int days_ahead;
    public boolean short_sell;
    public int succ_trades, total_trades,n_drops;
    public double trading_cost, risk_free,sharpeRatio,avg,risk_reward;
    public double max_drop;
    public boolean priceDataSet = false;
    public boolean trading;  //in trading mode, must have priceDataSet = true
    public int trading_func = 0;
  
    public double tradeDelta = 0.0; //threshold for trades
 

    //------ out_of_sample_tradingSweep ----------------------
    public double out_max_drop = 0.0; 
    public double out_ratio = 0.0; 
    public double out_ROI = 0.0;  
    public double out_n_drops = 0;
    public double sign_correct = 0;

    public boolean sweepMode = false;
    public double opt_value = 0.0; 
    public JTextField opt_valText;
    public int true_out = 0;
    public boolean keep_n_fixed = false;  //when data is added, keeps the number of obs fixed    
    
    
    //----- trading parameters panel------------------------------------------
    private JSlider riskfreeSlider;
    private JTextField riskfreeText,tradingcostText;
    private JCheckBox shortCheck;
    private JSlider tradingcostSlider,tradingfreqSlider;
 
    public JCheckBox diffSigCheck;
    //----- trading optimization panel-----------------------------------------
    private JTextField fin_alphaText;
    private JTextField fin_lambdaText,starobsText,outsampText,otradeRatioText,onumDropsText;
    private TradingOptGrid gridContourPanel;
    private JButton gridSearchButton;
    private JProgressBar gridprogressBar;
    
    @SuppressWarnings("rawtypes")
	private JComboBox optcritCombo;
    public JPanel optimtradePanel,outsampPanel;
    private JButton simannButton,compSweepButton;    
    public JPanel dataMixPanel;
    public ScalePanel optScale;    
    private JLabel otradeRatioLabel,outsampLabel,startObsLabel,omaxDropLabel,ondropsLabel,oROILabel;
    private JTextField oROIText,omaxDropsText;
    public JScrollBar outsampBar,startobsBar;
    public JCheckBox cutOptimizeBox;

    @SuppressWarnings("rawtypes")
	public JComboBox optVariablesCombo;

    //------ adaptive update filtering--------------------------------------------
    private JScrollBar c_updateBar,al_updateBar,d2_updateBar;
    private JTextField al_updateText,lam_updateText,c_updateText,d2_updateText,d_updateText,l_updateText,n_updateText,s_updateText;
    private JCheckBox autoUpdate,plotUpdateBox,shadeBox, uni_updateCheck;   
    private JButton computeUpdate;
    private JLabel c_updateLabel,l_updateLabel,s_updateLabel2,n_updateLabel;
    private JLabel d2_updateLabel,d2_updateLabel1,d_updateLabel,d_updateLabel1;
    public JPanel regPanel2,adaptiveUpdatePanel;
    private JScrollBar l_updateBar,lam_updateBar,n_updateBar, s_updateBar,d_updateBar;
    
    //------- Spectral Density estimate panel ------------------------
    JPanel spectralDensPanel;
    JCheckBox activateEstimationCheck,argBox,autoUpdatesCheck,modBox,smoothCheck,periodogramBox,taperBox;
    JButton estimateARMAButton,setEstimateButton;
    JCheckBox estimateARMACheck,gaussCheck;
    
    @SuppressWarnings("rawtypes")
	JComboBox pComboBox,qComboBox,orthBasisBox,phaseTypeBox,smoothFuncBox;
    JLabel smoothFuncLabel,pLabel,phaseTypeLabel,qLabel,smoothScaleLabel,taperBasisLabel,taperDegreeLabel;
    JPanel taperPanel,gaussPanel,plotPanel,smoothPanel;
   
    JCheckBox[] sd_checkBox;
    JScrollBar smoothScaleBar,taperDegreeBar;
    JTextField smoothScaleText,taperDegreeText; 

    //----- trading statistics panel-------------------------------------------
    private JTextField avgText;
    private JTextField maxdropText;
    private JTextField ndropsText;
    private JTextField ratioText;
    private JTextField roiText;
    private JTextField sharpeText;
    private JPanel statPanel;
    private JTextField succText;
    private JTextField totalText;
    double[] xf_turn_val;      

    public boolean compute_error = false;
    double rank_coeff;
    int tradeFrequency;
    double[] grid_trade;
    int gridRes;
    JCheckBox signalPlot, accountPlot, logreturnPlot, logpricePlot, linesPlot, buysellPlot, filteredPrice;
    JCheckBox outSampOptBox;
    FilterCoefCanvas updateCoeffPanel;
    AmpUpdateCanvas updateAmpPanel;

    private JTextField outavgText;
    private JTextField outmaxdropText;
    private JTextField outndropsText;
    private JTextField outratioText;
    private JTextField outroiText;
    private JTextField outsharpeText;
    private JTextField outsuccText;
    private JTextField outtotalText;
    private JTextField outForeText;
    private JTextField outRRText;
    
    JLabel rrLabel,foreLabel;    
    JPanel outstatPanel;
    public JPanel outsampStatPanel;
    public JPanel trueOutSamplePanel;
    
    JScrollBar true_outBar;
    JTextField true_outText;
    JLabel true_outLabel;

    JCheckBox price_filter;
    double[] price_indicator;

    //---- hf price -------
    double[] hf_price;
    int hf_price_length;
    boolean hf_price_set;
    int hf_period;
    boolean removeOutliers = false;


    static {System.loadLibrary("simSarima");}
    static {System.loadLibrary("sv");}
    static {System.loadLibrary("ms");} 
    
    public IMDFAPanel()
    {
          
          int i;         
        //--------- set initial values ----------------------------
          n_obs = 300; n_rep = 5; L = 28; Lag = 0; i1 = 0; i2 = 0; n_rep_fixed = n_rep;
          num1 = 0; num2 = 1; den1 = 20; den2 = 6; 
          lambda = 0; expweight = 0; w1 = ((1.0*num2)/(1.0*den2))*Math.PI; w0 = 0.0;
          smooth = 0; decay = 0; cross = 0; out = 0; decay2 = 0.0;
          expFunc = new double[K+1]; lambda_3 = 0.0; zpc_mod_optimize = false;
          onestep_diff = 0;
          //------ initialize w_const -------------------
          w_const = new double[9]; w_const[0] = 1.0;
          histParams = new ArrayList<String>(); n_hist = 0;
          complete_data = new ArrayList<double[]>(); 
          subset_data = new ArrayList<double[]>();
          

          filter_health = new Color[125];
          for(i=0;i<125;i++) 
          {filter_health[i] = new Color(Color.HSBtoRGB((float)i/255, (float)100/255, (float)62/255));}

          dd = 0; K = (int)n_obs/2; K1 = K+1; L1 = 200; S=12; nsamp = K; L2 = 0;
          p1 = 300; p2 = 12;
          
          x13model = 0; useX13filters = false;
          
          low=true;high=false;ramp=false;band=false; autoComp = false;
        //--------------------------------------------------------- 
         eof = System.getProperty("line.separator");
         extend = new double[2*L1+n_obs]; fore_extend = new double[36];
         symfilt = new double[n_obs - (L-1)];
         bc = new double[2*L1+1]; Gamma_x13 = new double[K1];
         iter = false; L_fixed = L; 
 
         // ----------------- Initialize ART and time-frequency ----------
         ARTfilter = new double[3]; ARTfilter[0] = 1.0; ARTfilter[1] = 0.0; ARTfilter[2] = 0.0;
         plotTF = 0;  plot_number = 0;
         savedFilters = new ArrayList<MDFAFilter>();


         
         days_ahead=0;
         short_sell=false;
         succ_trades=0; total_trades=0; n_drops=0;
         trading_cost=0; risk_free=0;
         max_drop=0.0;
         priceDataSet = false;

         sharpeRatio=0.0; avg = 0.0; risk_reward = 0;

         int width = 600; int height = 300; 
         seed = 301; cutoff = Math.PI/6;
         //initiate engine
         mdfa = new IMDFA(n_obs, n_rep, L, Lag, lambda, expweight, cutoff, i1, i2);
         mdfa_canvas = new IMDFAcanvas(width+600, height, n_obs, n_rep, L);
         filter_canvas = new FilterCanvas(width, height, K, n_rep);
         period_canvas = new PeriodCanvas(width, height, K, n_rep);
         coef_canvas = new FilterCoefCanvas(width, height, L, n_rep);
         setAccount_canvas(new accountCanvas(width, height, L, n_rep));        
         

         crystal_ball = new CrystalBallcanvas(500, 280, 10, 1, 10); 
         envision = false; 
         
         
         zpc = new ZPCFilter(n_obs); p_arma=2; q_arma=2;
         zpcFreqCanvas = new ZPCCanvas();

         specDens = new SpectralDensity(n_obs, n_rep); 
         specDensCanvas = new SpectralDensityCanvas(width, height, K, n_rep);         
         initSpectralDensityPanel();
         
         df = new DecimalFormat("##.##"); myBlue = new Color(146,196,210);
         df2 = new DecimalFormat("##.###");
         df4 = new DecimalFormat("##.#####");

         timeScrollPane = new JScrollPane();
         timeScrollPane.setPreferredSize(new Dimension(width, height)); 
         timeScrollPane.setViewportView(mdfa_canvas);
         timeScrollBar = timeScrollPane.getHorizontalScrollBar();
         

         timeScrollBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               mdfa_canvas.setCanvasShift(((JScrollBar)e.getSource()).getValue());
            }
         });

         score = new double[7];
         nobsLabel = new JLabel("Obs:"); nobsLabel.setHorizontalTextPosition(JLabel.RIGHT); 
         nobsLabel.setToolTipText("Total number of observations in simulated series.");
         sampLabel = new JLabel("Samples:"); sampLabel.setHorizontalTextPosition(JLabel.RIGHT);  
         LLabel = new JLabel("L:");   LLabel.setHorizontalTextPosition(JLabel.RIGHT);
         LLabel.setToolTipText("Increase the length L of the filter.");
         lamLabel = new JLabel("      \u03BB:"); lamLabel.setHorizontalTextPosition(JLabel.RIGHT);
         lamLabel.setToolTipText("Adjust \u03BB to increase influence on minimizing phase error of filter.");
         expLabel = new JLabel("W(\u03C9):");  expLabel.setHorizontalTextPosition(JLabel.RIGHT);
         expLabel.setToolTipText("Adjust exponent of weight function W(\u03C9) to increase influence on minimizing amplitude error in stop-band.");
         nrepLabel = new JLabel("Series:");  nrepLabel.setHorizontalTextPosition(JLabel.RIGHT);
         nrepLabel.setToolTipText("Total number of representing series in the information supporting the target series Y_t.");
         diffLabel = new JLabel("Differencing|");  diffLabel.setHorizontalTextPosition(JLabel.RIGHT);
         cutLabel = new JLabel("Filter cutoff |");  cutLabel.setHorizontalTextPosition(JLabel.RIGHT);
         lagLabel = new JLabel("   Lag:");   lagLabel.setHorizontalTextPosition(JLabel.RIGHT);
         cutNLabel = new JLabel("Numerator:");  cutNLabel.setHorizontalTextPosition(JLabel.RIGHT);
         cutDLabel = new JLabel("Denominater:");  cutDLabel.setHorizontalTextPosition(JLabel.RIGHT);      
         simulateData = new JLabel("Simulate Data: "); simulateData.setHorizontalTextPosition(JLabel.RIGHT);         
         simulateData.setToolTipText("Different options for simulation data.");
         l1Label = new JLabel("M:"); l1Label.setHorizontalTextPosition(JLabel.RIGHT);         
         l1Label.setToolTipText("If using simulated data, length of approximated symmetric filter.");  
         lagLabel.setToolTipText("Lag for revision or forecasting. Positive lag produces revision, negative produces forecast.");  

         timeDiagnostics = false;
         computeFilter = true;
         reCompFilter = true;
         cPeriod = true; 
         data_uploaded = false;
         n_out_samp = 0;       
         prev_out_samp = 0;

         setup=true;
         //---- setup controls -------        
         tfilter = new FilterDesign(K,L1,w0,w1);

         //---setup optimization panel------------
         gridContourPanel = new TradingOptGrid();
         mdfaPlotPane = new JTabbedPane(JTabbedPane.TOP);
         initDiagnosticPanel();     
         initiateAdaptiveUpdatePanel();
         setupCheckButtons();
         setupFilterAdjusters();
         setupFilterDesign();        
         mdfa.setReg(true);                
         simulateIMDFA(false);
         adjustPlotChecks();
         setFilter();        
         initZPCFilter();
         initZPCPanel(); 
         setupDesign();      
         setup=false;           
         System.loadLibrary("simSarima");


         initTradingParameterDialog();
         initOptimizationDialog();
         initTradingStatDialog();
         initOutTradingStatDialog();
         initiateOutSamplePanel();
         initiateTrueOutSamplePanel();
       

    }

     
    /*-------------------------------------------------------------------
      Input target data targ of length n_obs+2*n_sym 
      n_sym = L1, length of one sided filter
      rep = vector of booleans, which series are providing info for target data
      list = array list of dictionary describing target data
    --------------------------------------------------------------------*/

    public void inputSimulatedData(double[] targ, ArrayList<double[]> list, boolean[] rep, int _n_sym)
    {
        
         int i,l,j,N; int _n_rep = 0; double sum;
         for(i=0;i<series_max;i++)  {if(rep[i]) {_n_rep++;}}
         data_uploaded = false;
         if(_n_rep > 1)
         {n_rep = _n_rep;}
         else 
         {n_rep = 1;}

         //System.out.println("Size = " + list.size());
         //System.out.println("rep = " + rep[0] + " " + rep[1] + " " + rep[2]);  
         //--- Turn off compute filter-------------------  
         computeFilter = false; reCompFilter = false; timeDiagnostics = false;
   
         //---- restart zpc panel 
         zpc_gene = false; 
         filterPlot[1].setSelected(false); filterPlot[2].setSelected(false); 
         filterPlot[1].setEnabled(false); filterPlot[2].setEnabled(false); 
         filterPlot[0].setSelected(true);

         complete_data.clear();
         n_sym = _n_sym; 
         simulate = false;         
            
         
         L1 = n_sym;  
         l1Bar.setValue(L1); l1Text.setText(""+L1);        
         mdfa_canvas.sym_plot = false;       
        
         n_obs = targ.length - n_sym;   
         n_tot = targ.length;
         extend = new double[n_tot];
         
         n_rep_fixed = n_rep;
         n_obs_fixed = n_obs; 
         N = n_obs; 
         
         //System.out.println("n_rep_fixed = " + n_rep_fixed);
         
         n_rep = 1;
         nrepBar.setValue(1);
         activateWConst(i1,false);       
         
         //----------- update data ---------------------------------------------------
         double[] gdp = new double[n_rep*N];
         this.target = new double[n_tot];
         System.arraycopy(targ, 0, extend, 0, n_tot);
         System.arraycopy(targ, 0, this.target, 0, n_tot);
         
         //System.out.println("L1 = " + L1 + " N = " + N + " n_tot =  " + n_tot + " n_rep=  " + n_rep);
         
         
         for(i=0;i<N;i++)
         {gdp[i] = targ[i];}

         _n_rep = 1;
         for(j=1;j<9;j++)
         {
           if(rep[j-1])
           {
            //System.out.println("added " + j);
            double[] temp = list.get(j-1);         
            complete_data.add(temp);
            //for(i=0;i<N;i++)
            //{
            //  gdp[N*_n_rep + i] = temp[i+L1];            
            //}
            _n_rep++;
           }
         }
     
         //----- set outof_sample ----------------------------------------------------------
         b_outof_sample = false; outof_sample.setMaximum(L1); outof_sample.setEnabled(true); n_out_samp = 0;
         if(L1 > 0) 
         {

            System.out.println("SHOULDNT GET HERE");
            b_outof_sample = true;  
            tfilter.setL1(L1); tfilter.recomputeGamma(); setFilter();
	    for(i=L1+(L-1);i<N;i++)
  	    {
      	      sum = 0.0;
              for(l=0;l<L1;l++)
              {sum = sum + bc[l]*extend[i+l];} 
              for(l=1;l<L1;l++)
              {sum = sum + bc[l]*extend[i-l];}
              symfilt[i-(L1+(L-1))] = sum;         
            }
         
      
            mdfa_canvas.setSymSignal(symfilt); 
         }
         
         adjustPlotChecks();
         setTimeSeries(gdp, n_obs, n_rep);
         setLambda(lambda); 

         data_uploaded = true; computeFilter = true; reCompFilter = true;
         mdfa.computeFilterGeneral(computeFilter,false);

         updatePlots(true,true);
 
         timePlot[11].setSelected(false); timePlot[11].setEnabled(false);
         updateData();
         
         //-- reinitialize
         l1Bar.setValue(36);
         
         System.out.println("Initialize fixed signal");
         fixed_account = new double[n_obs_fixed];
         fixed_signal = new double[n_obs_fixed];
         
         //System.out.println("size complete_data = " + complete_data.size());
         
    }

    public void setHFPrice(double[] hf, int per)
    {
      hf_price = new double[hf.length];
      hf_price_length = hf.length;
      System.arraycopy(hf,0,hf_price,0,hf_price_length);
      hf_price_set = true;
      hf_period = per;
      System.out.println("hf_price set");
    }

    //---------- update filter using better optimization out-of-sample -------

    public void updateSignal()
    {

      if(b_outof_sample && !reCompFilter)
      {
         int n_update = n_updateBar.getValue();
         int l_update = l_updateBar.getValue();
         double lam_update = .001*lam_updateBar.getValue();
         double al_update =  .1*al_updateBar.getValue();
         double s =  .001*s_updateBar.getValue();
         double d =  .001*d_updateBar.getValue();
         double d2 =  .001*d2_updateBar.getValue();
         double c =   .001*c_updateBar.getValue();
         
         if(uni_updateCheck.isSelected())
         {mdfa.updateSignalOut_Univ(n_update, l_update, lam_update, al_update, s, d, d2, c);}
         else
         {mdfa.updateSignalOut(n_update, l_update, lam_update, al_update, s, d, d2, c);}
          
         updatePlots(false,false);
      }

    }

    
    public void updateData()
    {
      
      if(data_uploaded && !simulate)
      {
      /* The idea here is to have up to K series stored in the module and be able to select any set of series from panel   */
 
         int i,j,N; int _n_rep = 1; 
         for(i=0;i<series_max;i++)  
         {
           if(expVariable[i].isSelected() && expVariable[i].isEnabled()) 
           {_n_rep++;}
         }
         data_uploaded = false;
         //if(_n_rep > 1) {n_rep = _n_rep+1;}
         //else  {n_rep = 1;}
         n_rep = _n_rep;
         //System.out.println("Update Data");
          
         //--- Turn off compute filter-------------------  
         computeFilter = false; reCompFilter = false; timeDiagnostics = false;
   
         //---- restart zpc panel 
         zpc_gene = false; 
         filterPlot[1].setSelected(false); filterPlot[2].setSelected(false); 
         filterPlot[1].setEnabled(false); filterPlot[2].setEnabled(false); 
         filterPlot[0].setSelected(true);

      
         simulate = false;         
         nrepBar.setValue(n_rep);         
       
         //l1Bar.setValue(L1); l1Text.setText(""+L1);        
         mdfa_canvas.sym_plot = false;       

         activateWConst(i1,false);       
      
         N = n_obs; 
         //----------- update data ---------------------------------------------------
         double[] gdp = new double[n_rep*N];
           
         //System.out.println("L1 = " + L1 + " N = " + N + " n_tot =  " + n_tot + " n_rep=  " + n_rep);

         for(i=0;i<N;i++)
         {gdp[i] = this.target[L2+i];}

         _n_rep = 1;
         for(j=0;j<series_max;j++)
         {
           if(expVariable[j].isSelected() && expVariable[j].isEnabled())
           {
            double[] temp = complete_data.get(j);  
            for(i=0;i<N;i++)
            {
              gdp[N*_n_rep + i] = temp[L2+i];            
            }
            _n_rep++;
           }
         }
     
         //----- set outof_sample ----------------------------------------------------------
         /*b_outof_sample = false; outof_sample.setMaximum(L1); outof_sample.setEnabled(true); n_out_samp = 0;
         if(L1 > 0) 
         {

            System.out.println("SHOULDNT GET HERE");
            b_outof_sample = true;  
            tfilter.setL1(L1); tfilter.recomputeGamma(); setFilter();
	    for(i=L1+(L-1);i<L1+N;i++)
  	    {
      	      sum = 0.0;
              for(l=0;l<L1;l++)
              {sum = sum + bc[l]*extend[i+l];} 
              for(l=1;l<L1;l++)
              {sum = sum + bc[l]*extend[i-l];}
              symfilt[i-(L1+(L-1))] = sum;         
            }
         
            for(i=0;i<36;i++)
            {fore_extend[i] = extend[n_obs+L1+i-1];}

            mdfa_canvas.setForeExt(fore_extend, ext_forecast);
            mdfa_canvas.setSymSignal(symfilt); 
         }*/
         
         //if trading, turn off briefly
         boolean temp_trading = trading;
         trading = false;
         
         setTimeSeries(gdp, n_obs, n_rep);
         trading = temp_trading;
         //setLambda(lambda); 
                 
         data_uploaded = true; computeFilter = true; reCompFilter = true;
         mdfa.computeFilterGeneral(computeFilter,false);
         
         updatePlots(true,true);
         adjustPlotChecksRevise();
         //timePlot[11].setSelected(false); timePlot[11].setEnabled(false);    
 
      }
    
    }
    

  public void toggleH0Filtering(boolean g) {mdfa.useH0Set(g);}
  
  
  public void setH0Filter(File file)
  {

         Double D; int i,k;
         String strline; String[] tokens; String delims = "[ ]+";  
         double[] vals;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
          
           strline = br.readLine(); //get L and nreps
           tokens = strline.split(delims);
          
           int tempL = (new Integer(tokens[0])).intValue();
           int temp_nreps = (new Integer(tokens[1])).intValue();
            
           System.out.println(tempL + "  " + temp_nreps); 
            
           if(tempL < n_obs && temp_nreps == n_rep)
           {
             
             LBar.setValue(tempL); L = tempL;
             
             vals = new double[n_rep*L];
             k=0;
             while((strline = br.readLine()) != null)
             {
                tokens = strline.split(delims);
                
                for(i=0;i<tokens.length;i++)
                {
                 D = new Double(tokens[i]);
                 vals[L*i + k] = D.doubleValue();
                }
                k++;
             }
             mdfa.setH0B0(vals, n_rep, L);   
           }
           System.out.println("H0 Filter set");
           br.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
         catch(IOException ioe){System.out.println("IO out error..." + ioe);}
         
    }     
    
    
    public void adjustSymmetricDataPoints()
    {
       int i,j,l,N; 
       int total = n_obs_fixed - L1 - L2;   
       double sum;    
       
       if(!simulate)
       {
        if((L1 >= 36) && (total > L))
        {
   
         L_fixed = L;
         mdfa_canvas.win_shift = 0;

         b_outof_sample = false; outof_sample.setMaximum(L1); 

         N = n_obs_fixed; 
         n_obs = N - L1 - L2;   
         K = n_obs/2; K1 = K+1;     

         // make new data matrix
         double[] gdp = new double[n_rep*n_obs];
 
         fore_extend = new double[36];
         
         if(n_obs + L2 - L1 - (L-1) > 0)
         {symfilt = new double[n_obs + L2 - L1 - (L-1)];}
         else
         {symfilt = new double[2];}
           
         for(i=0;i<n_obs;i++)
         {gdp[i] = extend[L2+i];}

         int count = 1;
         for(j=0;j<series_max;j++)
         {
           if(expVariable[j].isSelected() && expVariable[j].isEnabled())
           {
            double[] temp = complete_data.get(j);  
            for(i=0;i<n_obs;i++)
            {gdp[n_obs*count + i] = temp[L2+i];}
            count++;
           }  
         }
         
         


         /*for(i=0;i<36;i++)
         {fore_extend[i] = extend[n_obs+i-1];}
         mdfa_canvas.setForeExt(fore_extend, ext_forecast);*/


         if(!reCompFilter)
         {mdfa.computeRTSE(gdp, n_obs, n_rep); updatePlots(true, false);} //-- applies the older filter, nothing recomputed
         else
         {setTimeSeries(gdp, n_obs, n_rep); } //-- recomputes everything with new data 

       
         b_outof_sample = true;   
         if(!useX13filters)
         {
          tfilter.setL1(L1); tfilter.recomputeGamma(); setFilter();
         }
         else
         {computeX13TargetFilter(); tfilter.setL1(L1); tfilter.setGeneralSymmetric(Gamma_x13); setFilter();}

         if(n_obs + L2 - L1 - (L-1) > 0)
         {
	  for(i=L1+(L-1);i<n_obs;i++)
  	  {
      	      sum = 0.0;
              for(l=0;l<L1;l++)
              {sum = sum + bc[l]*extend[i+l];} 
              for(l=1;l<L1;l++)
              {sum = sum + bc[l]*extend[i-l];}
              symfilt[i-(L1+(L-1))] = sum;         
         }
         mdfa_canvas.setSymSignal(symfilt); 
         //timePlot[11].setSelected(true); 
         timePlot[11].setEnabled(true);
         }
         else
         {
           mdfa_canvas.setSymSignal(symfilt); 
           timePlot[11].setEnabled(false);
         }




         timeDiagnostics = false; 
         /*if(n_out_samp == 0) 
         {
           timeDiagnostics = true;       
           filterScore();    
           //System.out.println("Time domain diagnostics available. Select Diagnostics panel under Panels in menu"); 
         }
         else 
         {
           timeDiagnostics = false; 
           System.out.println("To compute time domain diagnostics, set number of out-of-sample to 0.");     
         }*/

       }
       else if(total > L)
       {
         //double[] temp;  
         b_outof_sample = false; 
         timeDiagnostics = false; 
         outof_sample.setMaximum(0);
         N = n_obs_fixed;  
         n_obs = N;   
         K = n_obs/2; K1 = K+1;     
         
         // make new data matrix
         double[] gdp = new double[n_rep*n_obs];
           
         for(i=0;i<n_obs;i++)
         {gdp[i] = extend[i];}

         for(j=1;j<n_rep;j++)
         {
            double[] temp = complete_data.get(j-1);
            for(i=0;i<n_obs;i++)
            {
              gdp[n_obs*j + i] = temp[i]; 
            }
         }
         
         //outof_sample.setMaximum(0); mdfa_canvas.sym_plot = false;
         //timePlot[11].setSelected(false); timePlot[11].setEnabled(false);  

         //keep the symmetric filter, just adjust the position
         mdfa_canvas.win_shift = 36 + (L_fixed - L); 
         //System.out.println("win_shift = " + mdfa_canvas.win_shift);

         
         if(!reCompFilter)
         {
           mdfa.computeRTSE(gdp, n_obs, n_rep); 
           if(plot_number == 1 || plot_number == 2)
           {
             zpc.setData(gdp, n_rep, n_obs); 
             zpc.injectZPCGene2();
             updatePlots(true, false);
           }

         } //-- applies the older filter, nothing recomputed
         else
         {
           setTimeSeries(gdp, n_obs, n_rep); 
         } //-- recomputes everything with new data 

       }
      }

   }

   
   public void getHFPeriodogram(String[] instrum, String date1, String date2, boolean time_frame, 
                             String time1, String time2, int freq)
   {
   
     String addInst;  
     String align = "";
     String names;
     String closes;
     names = "asset";
     
     if(freq == 0)
     {align = names + "<-align.time(to.minutes(" + instrum[0] + "[,5:6],1,name=NULL),60)['T13:30/T20:00',]";}

     if(freq == 1)
     {align = names + "<-align.time(to.minutes3(" + instrum[0] + "[,5:6],name=NULL),3*60)['T13:30/T20:00',]";}

     if(freq == 2)
     {align = names + "<-align.time(to.minutes5(" + instrum[0] + "[,5:6],name=NULL),5*60)['T13:30/T20:00',]";}

     if(freq == 3)
     {align = names + "<-align.time(to.minutes10(" + instrum[0] + "[,5:6],name=NULL),10*60)['T13:30/T20:00',]";}

     if(freq == 4)
     {align = names + "<-align.time(to.minutes15(" + instrum[0] + "[,5:6],name=NULL),15*60)['T13:30/T20:00',]";}

     if(freq == 5)
     {align = names + "<-align.time(to.minutes30(" + instrum[0] + "[,5:6],name=NULL),30*60)['T13:30/T20:00',]";}

     if(freq == 6)
     {align = names + "<-align.time(to.hourly(" + instrum[0] + "[,5:6],name=NULL),60*60)['T13:30/T20:00',]";}

     if(freq == 7)
     {align = names + "<-to.daily(" + instrum[0] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 8)
     {align = names + "<-to.weekly(" + instrum[0] + "[,5:6],drop.time=TRUE,name=NULL)";}

     if(freq == 9)
     {align = names + "<-to.monthly(" + instrum[0] + "[,5:6],indexAt=’yearmon’,drop.time=FALSE,name=NULL)";}

  
     closes = "closes <- na.locf(Cl("+names+"))";
     



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
     
       //indexTZ(GOOG.O) <- "America/New_York"
       //GOOG.O<-GOOG.O[!is.na(GOOG.O$Trade.Price)]
     
       //----- Call the symbols -------------------------
       addInst = "getSymbols.FI('" + instrum[0] + "',from='" + date1 + "', to='" + date2 + "')"; 
       caller.getRCode().addRCode(addInst);
       caller.getRCode().addRCode(align);
       caller.getRCode().addRCode(closes);
       caller.getRCode().addRCode("closes <- log(closes)");
       caller.getRCode().addRCode("closes <- diff(closes)");
       caller.getRCode().addRCode("closes[1,] <- 0");
       caller.getRCode().addRCode("closes[is.na(closes)]<-0");  //--- set any final NAs to 0
       caller.getRCode().addRCode("series_list<-lapply(as.list(closes), coredata)");
       caller.runAndReturnResult("series_list");
       mdfa.setFullRangeData(caller.getParser().getAsDoubleArray("Close"));
 
    }     
    catch(Exception e){System.out.println(e.toString());}

    //now compute the periodo-weight----
    mdfa.computeSpecDensPeriodogram();
    
    //now recompute the filter-------
    if(reCompFilter)
    {computeFilterNow();}    
  }

  public void togglePeriodoWeight(boolean o)
  {
    if(o){mdfa.useZPCWeight(1);}
    else {mdfa.useZPCWeight(0);}
    computeFilterNow();
  }
   
   /*-------------------------------------------------------------------
      
      Function that applies new uploaded filter to current data loaded
    
    
    ------------------------------------------------------------------*/
   
   public void applyNewFilter()
   {
       int i,j; 
       if(!simulate)
       {   
         L_fixed = L;
         mdfa_canvas.win_shift = 0;

         //N = n_obs_fixed; 
         //n_obs = N - 2*L1;   
         //K = n_obs/2; K1 = K+1;     

         // make new data matrix
         double[] gdp = new double[n_rep*n_obs];
 
         //fore_extend = new double[36];
         //symfilt = new double[n_obs - (L-1)];
           
         for(i=0;i<n_obs;i++)
         {gdp[i] = extend[L2+i];}

         int count = 1;
         for(j=0;j<series_max;j++)
         {
           if(expVariable[j].isSelected() && expVariable[j].isEnabled())
           {
            double[] temp = complete_data.get(j);  
            for(i=0;i<n_obs;i++)
            {gdp[n_obs*count + i] = temp[L2+i];}
            count++;
           }  
         }
       
         if(!reCompFilter)
         {mdfa.computeRTSE(gdp, n_obs, n_rep); updatePlots(true, false);} //-- applies the older filter, nothing recomputed
         else
         {setTimeSeries(gdp, n_obs, n_rep); } //-- recomputes everything with new data 
      }

   }
   
   

   
   
   
   
   
   
   
   
   
  public void computeOptimalFilter()
  {
     if(reCompFilter)
     {
      double[] opt_params = mdfa.computeOptimalFilter(symfilt, w_1, w_2, w_3);
      opt_lambda = opt_params[0]; opt_alpha = opt_params[1];

      if(computeFilter)
      {computeFilter = false;} 
      lambdaBar.setValue((int)opt_lambda*10); 
      //setLambda(opt_lambda);
      computeFilter = true;
      expBar.setValue((int)opt_alpha*10); 
     }
  } 

  public void filterScore()
  {
 
      if(b_outof_sample)
      {
        if(symfilt.length == mdfa_canvas.xf.length-L1)
        {
          score = new double[7];
          score = mdfa.filterScore(mdfa_canvas.xf, symfilt, 8);
        }
      }
      else
      {
        System.out.println("In-sample partition must first be created using the 'M' scrollbar to build approximation of symmetric filter."); 
      }
      updateDiagnosticPanel();
  }





    public void addOutofSample()
    {
        
      if(!simulate && b_outof_sample) 
      { 
         int i,j,N,N_fixed; 
         N_fixed = n_obs_fixed;
         n_out_samp = outof_sample.getValue(); 
         outof_sampleText.setText(""+n_out_samp);
          
         //n_obs = N_fixed - L1 - L2 + n_out_samp;   
         n_obs = N_fixed - L1 - L2;  
         N = n_obs; K = n_obs/2; K1 = K+1;     
         
         //if(n_out_samp > 0) {timeDiagnostics = false;} 
         //else {timeDiagnostics = true; filterScore();}
         
         double[] gdp = new double[n_rep*N];
         
         //for(i=0;i<N;i++)
         //{gdp[i] = extend[L2+i];}
         for(i=0;i<N;i++)
         {gdp[i] = extend[n_out_samp+i];}
         
         int count = 1; 
         for(j=0;j<n_rep;j++)
         {
           if(expVariable[j].isSelected())
           {
            double[] temp = complete_data.get(j);
            for(i=0;i<N;i++)
            {
              gdp[N*count + i] = temp[n_out_samp+i]; 
            }
            count++;
           }
         }




         if(!reCompFilter)
         {
            
            mdfa.computeRTSE(gdp, n_obs, n_rep); 
            
            if(continuousUpdate.isSelected() && plotUpdateBox.isSelected())
            {updateSignal();}
            else
            {mdfa.applyUpdatedFilter(plotUpdateBox.isSelected()); updatePlots(true, false);}

            if(zpc_gene)
            {zpc.computeRTSE(gdp, n_obs, n_rep); updatePlots(true, false);}

            if(!continuousUpdate.isSelected() && !plotUpdateBox.isSelected())
            {updatePlots(true, false);}
             //-- applies the older filter, nothing recomputed
         }
         else
         {setTimeSeries(gdp, n_obs, n_rep); } //-- recomputes everything with new data 


     }
    }

  

  

    
    public void adjustPlotChecks()
    {
      for(int i=0;i<n_rep;i++)
      {
        timePlot[i].setSelected(false); timePlot[i].setEnabled(true);
        coeffPlot[i].setSelected(false); coeffPlot[i].setEnabled(true);
        periodPlot[i+1].setSelected(false); periodPlot[i+1].setEnabled(true);
        freqPlot[i].setSelected(false); freqPlot[i].setEnabled(true);
        expVariable[i].setSelected(false); expVariable[i].setEnabled(false);
      }
      for(int i=n_rep;i<=10;i++)
      {
        timePlot[i].setSelected(false); timePlot[i].setEnabled(false);
        coeffPlot[i].setSelected(false); coeffPlot[i].setEnabled(false);
        periodPlot[i+1].setSelected(false); periodPlot[i+1].setEnabled(false);
        freqPlot[i].setSelected(false); freqPlot[i].setEnabled(false);
        expVariable[i].setSelected(false); expVariable[i].setEnabled(false);
      }
      coeffPlot[n_rep].setSelected(true); coeffPlot[n_rep].setEnabled(true);
      timePlot[0].setSelected(true); freqPlot[10].setSelected(false); freqPlot[10].setEnabled(true);
      expVariable[n_rep_fixed-1].setSelected(false); expVariable[n_rep_fixed-1].setEnabled(false);
      //expVariable[0].setSelected(true);

      for(int i=0;i<n_rep_fixed;i++)
      {expVariable[i].setSelected(false); expVariable[i].setEnabled(true);}
      for(int i=n_rep_fixed;i<=10;i++)
      {expVariable[i].setSelected(false); expVariable[i].setEnabled(false);}
    }

    public void adjustPlotChecksRevise()
    {
      for(int i=0;i<n_rep;i++)
      {
        timePlot[i].setSelected(false); timePlot[i].setEnabled(true);
        coeffPlot[i].setSelected(false); coeffPlot[i].setEnabled(true);
        periodPlot[i+1].setSelected(false); periodPlot[i+1].setEnabled(true);
        freqPlot[i].setSelected(false); freqPlot[i].setEnabled(true);
        
      }
      for(int i=n_rep;i<=10;i++)
      {
        timePlot[i].setSelected(false); timePlot[i].setEnabled(false);
        coeffPlot[i].setSelected(false); coeffPlot[i].setEnabled(false);
        periodPlot[i+1].setSelected(false); periodPlot[i+1].setEnabled(false);
        freqPlot[i].setSelected(false); freqPlot[i].setEnabled(false);

      }
      coeffPlot[n_rep].setSelected(true); coeffPlot[n_rep].setEnabled(true);
      timePlot[0].setSelected(true); freqPlot[10].setSelected(false); freqPlot[10].setEnabled(true);
    }
    
    
    
    
    
    public void setSimulate(boolean sim)
    {
       simulate = sim; 
       if(sim) 
       {
         mdfa_canvas.win_shift = 0;
         n_rep = 5; mdfa.set_nreps(n_rep); L1 = 200; l1Bar.setValue(L1);  
         b_outof_sample = false;  outof_sample.setEnabled(false); outof_sampleText.setText("0"); 
         mdfa_canvas.setNRep(n_rep); filter_canvas.setNRep(n_rep); period_canvas.setNRep(n_rep);  
         activateWConst(i1,false);
         simulateIMDFA(false); 
         adjustPlotChecks();
       }
    }


    public void computeExpWeightFunc()
    {
       int i;
       int omega_Gamma = (int)(w1*K/Math.PI);  //omega_Gamma<-as.integer(cutoff*K/pi)
       int omega_Gamma0 = (int)(w0*K/Math.PI);

       double scale2 = 1.0*omega_Gamma0;
       double scale0 = 1.0*(K-omega_Gamma);
       double scale = Math.max(scale0,scale2);

       expFunc = new double[K1];
       for(i=0;i<=omega_Gamma0;i++) expFunc[i] = Math.pow((scale*(omega_Gamma0 - i)/scale2)*Math.PI/K + 1, expweight/10);
       for(i=omega_Gamma0+1;i<omega_Gamma;i++) expFunc[i] = 1.0; 
       for(i=omega_Gamma; i<K1; i++) expFunc[i] = Math.pow((scale*(i-omega_Gamma)/scale0)*Math.PI/K + 1, expweight/10);
    }


    public void saveFilteredData(int k)
    {
      String file = "iMetric-filter-" + Integer.toString(L) + 
                     "-"+Double.toString(smooth)+"-"
                        +Double.toString(decay)+"-"+Double.toString(cross)+"-"+k+".dat";
      try{  
            
           PrintWriter out = new PrintWriter(new FileWriter(file));

           out.println(mdfa.xf.length);
           for(int i=0; i < mdfa.xf.length; i++) {out.println(mdfa.xf[i]);}
           out.println(paramString);

           out.close(); System.out.println("Data successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();} 
    }

    public void saveParameters(boolean f)
    {
       paramString = new String("");
       paramString = Integer.toString(L)+"\n"+Integer.toString(Lag)+"\n"+Integer.toString(i1)+"\n"+Integer.toString(i2)+"\n"
                    +df.format(lambda)+"\n"+df.format(expweight)+"\n"
                    +df.format(smooth)+"\n"+df.format(decay)+"\n"+df.format(cross)+"\n"
                    +df.format(mdfa.criteria)+"\n"+df.format(mdfa.degrees);


       if(f)
       {
       String file = "iMetric-params-" + Integer.toString(L) + 
                     "-"+df.format(smooth)+"-"
                        +df.format(decay)+"-"+df.format(cross)+".dat";
        try{  
            
           PrintWriter out = new PrintWriter(new FileWriter(file));
           out.println(paramString);
           out.close(); System.out.println("Data successfully saved in " + file);
          } catch (IOException e) {e.printStackTrace();} 
       }

    }

    public void parameterSnapshot()
    {

       paramString = new String("");
       paramString = Integer.toString(L)+"\n"+Integer.toString(Lag)+"\n"+Integer.toString(i1)+"\n"+Integer.toString(i2)+"\n"
                    +df.format(lambda)+"\n"+df.format(expweight)+"\n"
                    +df.format(smooth)+"\n"+df.format(decay)+"\n"+df.format(cross)+"\n"
                    +df.format(mdfa.criteria)+"\n"+df.format(mdfa.degrees);

       mdfa_canvas.setParamSnapshot(paramString); 
    }


    public void uploadParameters(File file)
    {
         int count = 0;
         String strline; 
         double[] parameters = new double[10];
         Double D; double val;
         //default parameters
         parameters[1] = 0.56; parameters[8] = 28;
         String[] tokens; String delims = "[ ]+";
         int n_toks;
         
         computeFilter = false;
         try
         {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           strline = br.readLine();
           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks <= 2)
           {
             D = new Double(tokens[0]); w0 = D.doubleValue();
             D = new Double(tokens[1]); w1 = D.doubleValue();  
             
             if(w0 > 0.0) //bandpass
             {bandBut.setSelected(true);}
             else{lowBut.setSelected(true);}
      
             omega0Bar.setValue((int)(100*w0)); omega1Bar.setValue((int)(100*w1));
             tfilter.setBand(w0,w1);                
           }
           else if(n_toks == 4)
           {
             D = new Double(tokens[0]); w0 = D.doubleValue();
             D = new Double(tokens[1]); w1 = D.doubleValue();  
             D = new Double(tokens[2]); w2 = D.doubleValue();
             D = new Double(tokens[3]); w3 = D.doubleValue(); 
             
             highBut.setSelected(true);
      
             omega0Bar.setValue((int)(100*w0)); omega1Bar.setValue((int)(100*w1));
             omega2Bar.setValue((int)(100*w2)); omega3Bar.setValue((int)(100*w3));
             tfilter.setBand(w0,w1); tfilter.setBand2(w2,w3);                       
           }
          
           count = 2;
           while((strline = br.readLine()) != null)
           {    
              
              D = new Double(strline); 
              val = D.doubleValue();
    
              parameters[count] = val;
              System.out.println(val);
              count++;
           }
           din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
         catch(IOException ioe){System.out.println("IO out error..." + ioe);}
        
 
         //------ set parameters -------------
         
         if(parameters[2] < 4000.0 && parameters[3] < 100.0)
         {
          lambda = parameters[2]; expweight = parameters[3];
          lambdaBar.setValue((int)(10*lambda)); expBar.setValue((int)(10*expweight));
         }
         
         if(parameters[4] >= 0.0 && parameters[5] >= 0.0)
         {
          smooth = parameters[4]; decay = parameters[5]; 
          smoothBar.setValue((int)(1000*smooth)); decayBar.setValue((int)(1000*decay));
         }

         if(parameters[6] >= 0.0 && parameters[7] >= 0.0)
         {
          decay2 = parameters[6]; cross = parameters[7];
          decay2Bar.setValue((int)(1000*decay2)); crossBar.setValue((int)(1000*cross));
         } 
         
         if((int)parameters[8] < n_obs && (int)parameters[8] >= 5)
         {
          L = (int)parameters[8];
          LBar.setValue(L);
         }
         
              
         Lag = (int)parameters[9];         
         lagBar.setValue(Lag+36);
                  
         
         computeFilter = true;
         computeFilterNow();        
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


 public  void rankCoefficient(double[] account, int n)
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
      
   rank_coeff = spear;   
  }  




  public void setTradingInterface(boolean r) 
  {
    if(priceDataSet) 
    {trading = r;}
    else
    {System.out.println("The Financial Trading platform must first have asset price data loaded");}
  }


  public void updateTradingInterface()
  {
    if(trading)
    {
      updatePlots(false,false);
    }
  }

  public void postTradingStatistics()
  {
    if(trading)
    {

      //succ_trades, total_trades,n_drops;
      //public double trading_cost, risk_free;
      //public double max_drop
      computeSharpeRatio(250);
      
      roiText.setText(""+df.format(account[account.length-1]));
      succText.setText(""+succ_trades);
      totalText.setText(""+total_trades);
      maxdropText.setText(""+df.format(max_drop));  
      ndropsText.setText(""+df.format(sign_correct));
      ratioText.setText(""+df.format((double)succ_trades/(double)total_trades));
      sharpeText.setText(""+df.format(sharpeRatio));
      avgText.setText(""+df.format(rank_coeff));
    }
  }


  public double[] integrate(double[] data, double s)
  {
      int n = data.length;
      double[] cs = new double[n]; double sum; int k;
      sum=Math.abs(s); double min = 1000000;

      for(k=0;k<n;k++)
      {
        sum = sum+data[k]; cs[k] = sum; 
        if(cs[k] < min) {min = cs[k];}
      }

      if(min < 0.0)
      {
        for(k=0;k<n;k++) {data[k] = data[k] - min;}
      }
      return cs;  
  } 

  public void toggleDiffSigTrading(boolean t) 
  {diff_sig_trading = t; updatePlots(false,false);}

  public void setTradingFunction(int c)
  {trading_func = c;}


  public void insampleTradingStrategy(double[] _price, double[] sig, int n)
  {  
    int i;
    xf_turn_val = new double[n];
    succ_trades = 0; total_trades = 0;
    trade_obs = n;
    account = new double[n];
    double[] logret = new double[n];
    double[] d_signal = new double[n];
    
    for(i=1; i<n;i++)
    { 
      logret[i] = _price[i] - _price[i-1];      
      d_signal[i] = sig[i] - sig[i-1];
    }  
    
    
    for(i = 1; i < n-1; i++)
    {
     if(Math.abs(d_signal[i]) > tradeDelta && Math.abs(sig[i]) > tradeDelta)
     {
     
      if(d_signal[i] > 0 && d_signal[i-1] < 0) {xf_turn_val[i] = 1;}
      else if(d_signal[i] < 0 && d_signal[i-1] > 0) {xf_turn_val[i] = -1;} 
      if(d_signal[i] > 0)  // going long on strategy
      {account[i] = account[i-1] + logret[i+1] - trading_cost;}
      else
      {account[i] = account[i-1] - logret[i+1] - trading_cost;}
     }
     else
     {account[i] = account[i-1];}
     
     if((account[i] - account[i-1]) >= 0) {succ_trades=succ_trades+1;}     
    }
    account[n-1] = account[n-2];
    total_trades = n;
 
    dropDown();
  
    double[] temp_acc = new double[trade_obs];
    System.arraycopy(account, 0, temp_acc, 0, trade_obs);
    rankCoefficient(temp_acc, trade_obs);    
  }
  
  public void setReturnStrategyOn()
  {return_strategy = true; System.out.println("Switched to Macro-Strategy mode");}
  
  
  
  public void insampleTradingDiff(double[] _price, double[] sig, int n)
  {
  
   if(return_strategy)
   {insampleTradingStrategy(_price,sig,n);}
   else
   {
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
    int sig_corr = 0;
   
    succ_trades = 0; total_trades = 0;
   
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    
    trade_obs = n;

    for(i=0;i<trade_obs;i++)
    {
      prix[i] = _price[i]; 
      t_signal[i] = sig[i];
      
    
      if(i > 0)
      {
        if(diff_sig_trading)
        {t_signal[i] = sig[i] - sig[i-1];}
        
        //logdf=prix[i]-prix[i-1];
        if(prix[i] < 0 && t_signal[i-1] < 0) {sig_corr++;}
        else if(prix[i] > 0 && t_signal[i-1] > 0) {sig_corr++;}
      }  
    }// System.out.println(prix[i] + "  " + t_signal[i]);}

    sign_correct = (double)sig_corr/(trade_obs-1.0); 


    // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
    int in_transaction = 0;  
    // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
    int out_transaction = 0; 
    account[0] = 0.0;    
    profit = 0; 
    //--- find where to start, first negative value
    i=0;
    price_bought = 0.0;
    price_borrowed = 0.0;

    while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = i;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {


        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
          price_sold = prix[i+1];
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

	xf_turn_val[i+1] = 1;
	price_bought = prix[i+1];
	in_transaction = 1;	

      }
      else if(in_transaction == 1 && t_signal[i+1] < -tradeDelta) //if in transaction and signal goes below, sell
      {

	 price_sold = prix[i+1];	 
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0; 
      }
      if(t_signal[i] > 0.0 && t_signal[i+1] < -tradeDelta) //went from positive to negative slope, short sell
      {
         if(short_sell)
         {   
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;
         }
      }
   }
   //---- close any transaction-------
   if(in_transaction == 1)
   {
      price_sold = prix[trade_obs-1]; 
      profit = price_sold - price_bought;
	 
      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
      amount = account[trade_obs-1]; 
      in_transaction = 0; 
   }
   else if(out_transaction == 1)
   {

      price_sold = prix[trade_obs-1];
      profit = price_borrowed - price_sold;

      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      //printf("total = %u\n", *total_trades);
      account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
      amount = account[trade_obs-1];
      out_transaction = 0;
   }
   else{account[trade_obs-1] = account[trade_obs-2];}  //else  //not in transaction, last worth is previous observation worth

   //----- with new account data, compute dropDown();
   dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
   rankCoefficient(temp_acc, trade_obs);

   //simPanel.plotData(t_signal,n);
   //simPanel.plotData(account, account.length);
   //simPanel.plotData(prix, prix.length);
 }
 
}


/*------------------------------------------------------------------

  This function computes trades according to the log-price of the data



------------------------------------------------------------------*/
public void changeTradingThreshold(double t) {tradeDelta = t; updatePlots(false,false);}

public void insampleTradingLogPrice(double[] _price, double[] sig, int n)
{
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
       
    d_signal = new double[n];

    succ_trades = 0; total_trades = 0;
   
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    double[] diff_sig = new double[n]; 
    diff_sig[0] = 0.0;
    account = new double[n];
    
    trade_obs = n;

    d_signal[0] = 0.0;
    for(i=0;i<trade_obs;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i];}// System.out.println(prix[i] + "  " + t_signal[i]);}

    for(i=0;i<trade_obs-1;i++)
    {diff_sig[i+1] = t_signal[i+1] - t_signal[i]; d_signal[i+1] = diff_sig[i+1];}       


    // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
    int in_transaction = 0;  
    // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
    int out_transaction = 0; 
    account[0] = 0.0;   
    profit = 0; 
    //--- find where to start, first negative value
    i=0;
    price_bought = 0.0;
    price_borrowed = 0.0;

    while(diff_sig[i] >= 0 && i < diff_sig.length-1) {i++;}  
    start = i;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(diff_sig[i+1] > tradeDelta && diff_sig[i] < 0.0) //new point positive, we see momentum, buy if above threshold
      {
       
        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
          price_sold = prix[i+1];
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

	xf_turn_val[i+1] = 1;
	price_bought = prix[i+1];
	in_transaction = 1;	

      }
      else if(in_transaction == 1 && diff_sig[i+1] < -tradeDelta) //if in transaction and signal goes below, sell
      {

	 price_sold = prix[i+1];	 
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0; 
      }
      if(diff_sig[i] > 0.0 && diff_sig[i+1] < -tradeDelta) //went from positive to negative slope, short sell
      {
         if(short_sell)
         {   
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;
         }
      }
   }

   account[trade_obs-1] = account[trade_obs-2];

   //----- with new account data, compute dropDown();
   dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
   rankCoefficient(temp_acc, trade_obs);

   //simPanel.plotData(diff_sig,n);
   //simPanel.plotData(account, account.length);
   //simPanel.plotData(prix, prix.length);
   //simPanel.plotData(d_signal,n);

}



/*-----------------------------------------------------------------------------------------

  First attempt at out-of-sample multiscale filter

  The idea is to have a signal computed for the low frequency. When out of sample data comes, and the filter
  signals a transaction, the trading interface then enters a higher frequency scale, where the higher 
  frequency filter domain is entered and the higher frequency filter attempts to find a better price during 
  the longer frequency duration. If no better price was found, then it simply obtain 

  In this first (simple) version, we don't filter at the higher frequency, but simply just try to find a 
  better price that at the given triggered price

  INPUT: low-freq price  length n
         low-freq signal length n
 
         high-freq price length 

-------------------------------------------------------------------------------------------*/


public void multiscaleFilterSimple(double[] _price, double[] sig, int n, double[] _hfprice, int hf_length, int period)
{

    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
   
    succ_trades = 0; total_trades = 0;
   
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    
    trade_obs = n;
    price_bought = 0.0;
    price_borrowed = 0.0;

    for(i=0;i<trade_obs;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i];}// System.out.println(prix[i] + "  " + t_signal[i]);}

       
    // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
    int in_transaction = 0;  
    // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
    int out_transaction = 0; 
    account[0] = 0.0;   
    profit = 0; 
    //--- find where to start, first negative value
    i=0;
    price_bought = 0.0;
    price_borrowed = 0.0;
    price_sold = 0.0;

    if(t_signal[0] < 0) 
    {
      out_transaction = 1;
      price_borrowed = prix[0];
    }
    else if(t_signal[0] > 0)
    {
      in_transaction = 1;
      price_bought = prix[0];
    }

    start = 1;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {

        //get best price (find price < prix[i+1])
        price_sold = getBetterPrice(prix[i+1],i+1,_hfprice,period,0);

        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 

          //price_sold = prix[i+1];
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

	xf_turn_val[i+1] = 1;
	//price_bought = prix[i+1];
        price_bought = price_sold;
	in_transaction = 1;	

      }
      else if(in_transaction == 1 && t_signal[i+1] < -tradeDelta) //if in transaction and signal goes below, sell
      {

	 //price_sold = prix[i+1];	 
         //find a higher price to sell at
         price_sold = getBetterPrice(prix[i+1],i+1,_hfprice,period,1);
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0; 
      }
      if(t_signal[i] > 0.0 && t_signal[i+1] < -tradeDelta) //went from positive to negative slope, short sell
      {
         if(short_sell)
         {   
	  xf_turn_val[i+1] = -1;
	  //price_borrowed = prix[i+1];
          price_borrowed = price_sold;
          out_transaction = 1;
         }
      }
   }
   //---- close any transaction-------
   if(in_transaction == 1)
   {
      price_sold = prix[trade_obs-1]; 
      profit = price_sold - price_bought;
	 
      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
      amount = account[trade_obs-1]; 
      in_transaction = 0; 
   }
   else if(out_transaction == 1)
   {

      price_sold = prix[trade_obs-1];
      profit = price_borrowed - price_sold;

      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      //printf("total = %u\n", *total_trades);
      account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;
      amount = account[trade_obs-1];
      out_transaction = 0;
   }
   //else  //not in transaction, last worth is previous observation worth
   else {account[trade_obs-1] = account[trade_obs-2];}

   //----- with new account data, compute dropDown();
   dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
   rankCoefficient(temp_acc, trade_obs);
}


//This function assumes the hfprice has period observations between each lower frequency observation

public double getBetterPrice(double p, int start, double[] hfprice, int period, int dir)
{
  int i = 0;
  int hf_start = start*period; // at the higher frequency, the starting point is period*start
  double price = p;
  if(dir == 1) //find a higher price to sell at
  {
    for(i=0;i<period;i++)
    {
      if(hfprice[hf_start + i] > p) {price = hfprice[hf_start + i]; System.out.println("higer price found at " + i + ": " + hfprice[hf_start + i] + " > " + p); break;}
    }
    if(price == p) {price = hfprice[hf_start + period];} //if could not find better price, stuck with next price
  }
  else if(dir == 0) //find a lower price to buy
  {
    for(i=0;i<period;i++)
    {
      if(hfprice[hf_start + i] < p) {price = hfprice[hf_start + i]; System.out.println("lower price found at " + i + ": " + hfprice[hf_start + i] + " < " + p); break;}      
    }
    if(price == p) {price = hfprice[hf_start + period];} //if could not find better price, stuck with next price
  }

  //System.out.println(count + " times a better price was found in the 5 min interval\n\n");
  return price;
}



/*--------------------------------------------------------------------------

   Idea for this trading function is as follows:
   
   1) We have an indicator which is a smooth real-time signal of the price
   2) The signal is the real-time signal of the log-return
   
   When the signal indicates a transaction, we check if it corresponds with the momentum 
   of the smooth indicator. If they agree, then we proceed with the transaction 
 
----------------------------------------------------------------------------*/


public void trading_function_Duplex(double[] _price, double[] ind, double[] sig, int n)
{

    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
   
    succ_trades = 0; total_trades = 0;
   
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    
    trade_obs = n;

    for(i=0;i<trade_obs;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i];}// System.out.println(prix[i] + "  " + t_signal[i]);}

       


    // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
    int in_transaction = 0;  
    // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
    int out_transaction = 0; 
    account[0] = 0.0;   
    profit = 0; 
    //--- find where to start, first negative value
    i=0;
    price_bought = 0.0;
    price_borrowed = 0.0;

    while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = i;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {
       if((ind[i+1] - ind[i]) > tradeDelta)
       {

        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
          price_sold = prix[i+1];
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

	xf_turn_val[i+1] = 1;
	price_bought = prix[i+1];
	in_transaction = 1;	
       }
      }
      else if(in_transaction == 1 && t_signal[i+1] < -tradeDelta) //if in transaction and signal goes below, sell
      {
       if((ind[i+1] - ind[i]) < -tradeDelta)
       {
	 price_sold = prix[i+1];	 
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0;
	}
      }
      if(t_signal[i] > 0.0 && t_signal[i+1] < -tradeDelta) //went from positive to negative slope, short sell
      {
       if((ind[i+1] - ind[i]) < -tradeDelta)
       {
         if(short_sell)
         {   
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;
         }
       }
      }
   }


   
 
   
   account[trade_obs-1] = account[trade_obs-2];

   //----- with new account data, compute dropDown();
   dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
   rankCoefficient(temp_acc, trade_obs);

}


/*------------------------------------------------------------------

   In this duplex trading function, the signal generated from the log return determines 
   the action on the indicator, which is a timely but somewhat noisy price signal 

-------------------------------------------------------------------*/


public void trading_function_Duplex_Price(double[] _price, double[] sig, double[] logsig, int n)
{

    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
       
    d_signal = new double[n];

    succ_trades = 0; total_trades = 0;
   
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    double[] diff_sig = new double[n]; 
    diff_sig[0] = 0.0;
    account = new double[n];
    
    trade_obs = n;

    d_signal[0] = 0.0;
    for(i=0;i<trade_obs;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i];}// System.out.println(prix[i] + "  " + t_signal[i]);}

    for(i=0;i<trade_obs-1;i++)
    {diff_sig[i+1] = t_signal[i+1] - t_signal[i]; d_signal[i+1] = diff_sig[i+1];}       


    // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
    int in_transaction = 0;  
    // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
    int out_transaction = 0; 
    account[0] = 0.0;   
    profit = 0; 
    //--- find where to start, first negative value
    i=0;
    price_bought = 0.0;
    price_borrowed = 0.0;

    while(diff_sig[i] >= 0 && i < diff_sig.length-1) {i++;}  
    start = i;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {
      account[i] = amount;
      if(diff_sig[i+1] > 0.0 && diff_sig[i] < 0.0) //new point positive, we see momentum, buy if above threshold
      {
       if(logsig[i+1] > tradeDelta) //duplex layer, if signal is positive (we trade)
       {
       
        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
          price_sold = prix[i+1];
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

	xf_turn_val[i+1] = 1;
	price_bought = prix[i+1];
	in_transaction = 1;	
       }
      }
      else if(in_transaction == 1 && diff_sig[i+1] < 0.0) //if in transaction and signal goes below, sell
      {
       if(logsig[i+1] < -tradeDelta)
       {
	 price_sold = prix[i+1];	 
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0; 
	}
      }
      if(diff_sig[i] > 0.0 && diff_sig[i+1] < 0.0) //went from positive to negative slope, short sell
      {
       if(logsig[i+1] < -tradeDelta)
       {
         if(short_sell)
         {   
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;
         }
       }
      }
   }

   account[trade_obs-1] = account[trade_obs-2];

   //----- with new account data, compute dropDown();
   dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
   rankCoefficient(temp_acc, trade_obs);

}



//--------- computes largest drop_down given an account, also computes mean and std of returns -----
  public void dropDown()
  {
   int i;
   max_drop = 0.0; 
   double drop;
   double mean, std;   
   int ticks = 0;


   n_drops = 0; mean = 0.0; std = 0.0;
   for(i=0;i<trade_obs-1;i++)
   {
     drop = account[i+1] - account[i]; 
     //printf("%lf\n",drop);
     if(drop < 0.0)
     {
       n_drops=n_drops+1;
       if(drop < max_drop)
       {max_drop = drop;}// printf("%lf\n",max_drop);}
     }

     if(drop != 0.0)
     { 
       ticks++;
       mean = mean + drop;
     }
   }

   mean = mean/ticks; std = 0.0;
   for(i=0;i<trade_obs-1;i++)
   {
     drop = account[i+1] - account[i];
     if(drop != 0)
     {std = std + (drop - mean)*(drop - mean);}
   }
   std = Math.sqrt(std);
 

 
   //return max_drop;
  }


  public void outsampleStats(int n_in)
  {
   int i;

   out_max_drop = 0.0; out_ratio = 0.0; out_ROI = 0.0;  out_n_drops = 0; 
   double drop;
      
   int ticks = 0;
   int start = trade_obs - n_in; // ---- number of out-of-sample

   //System.out.println("trade_obs = " + trade_obs + ", n_out = " + n_in);
   
   //1)  Compute number of drops -------------------------------
   for(i=start;i<trade_obs-1;i++)
   {
     drop = account[i+1] - account[i]; 
     //printf("%lf\n",drop);
     if(drop < 0.0)
     {
       out_n_drops=out_n_drops+1;
       if(drop < out_max_drop)
       {out_max_drop = drop;}// printf("%lf\n",max_drop);}
     }
     if(drop != 0.0)
     {ticks = ticks+1;}

   }

   //2) Total account increase/decrease----
   out_ROI = account[trade_obs-1] - account[start];
   if(ticks == 0) {out_ratio = 1.0;}
   else {out_ratio =  1.0*(ticks - out_n_drops)/(1.0*ticks);} 
    
   //System.out.println("outROI = " + out_ROI);
    
  }






  
 void setTradingFrequency(int c)
 {
   tradeFrequency = c;
 }
   

 public void setPriceData(double[] _p, int npr)
 {
   
    //---- no matter what, price data length must be n_obs-----
    if((int)(_p.length/npr) == n_obs)
    {
      priceDataSet = true;   
      t_price = new double[npr][n_obs];
      System.arraycopy(_p, 0, t_price[0], 0, n_obs);
      System.out.println("Portfolio for MDFA set");
    }
    else
    {System.out.println("Length of price data must be same as all other data in MDFA analysis");}
 }



 public void gridSearchTrading() 
 {


  int i,j,k;   
  double _lambda=0.0; 
  double _alpha=0.0; 

  gridRes = 50;
 
  double return_val = 0.0;
  double max = -10000000;
  double value,penalty=0.0;
  gridprogressBar.setMaximum(gridRes);  

  double t_lambda; double t_alpha; double d_lambda = 1; double d_alpha = .5;
  int choice = optcritCombo.getSelectedIndex();  

 
  int fl = n_obs - (L-1);
  double[] _price = new double[fl];
  double[] _sig = new double[fl];            

  grid_trade = new double[gridRes*gridRes];

  for(i=0;i<fl;i++)
  {_price[i] = t_price[0][L2+(L-1)+i];}
           
            
  t_lambda = 0.0;
  for(i=0;i<gridRes;i++)
  {
    t_alpha = 0.0;
    for(j=0;j<gridRes;j++)
    { 
        //compute new filter with new lambda, alpha
        
        mdfa.set_lambda(t_lambda); mdfa.set_exp(t_alpha); 
        mdfa.computeFilterGeneral(true,false);  
         
        for(k=0;k<fl;k++)  {_sig[k] = mdfa.xf[k];}

        if(trading_func == 0)
        {insampleTradingDiff(_price, _sig, fl);}
        else if(trading_func == 1)
        {insampleTradingLogPrice(_price, _sig, fl);}

        value = (account[account.length-1]);

        if(choice == 0)
        {penalty = value;}                     //1) maximize total worth

        else if(choice == 1)
        {penalty = -max_drop*max_drop;}          //2) minimize square of maximum drop 
    
        else if(choice == 2)
        {penalty = -max_drop*max_drop + value;}  //3) maximize total worth while minimizing max-drop   
        
        else if(choice == 3)
        {penalty = succ_trades/(1.0*total_trades);} //4) maximize ratio of successful trades
        
        else if(choice == 4)
        {penalty = (total_trades - n_drops)/(1.0*total_trades)*value;}

        else if(choice == 5)
        {penalty = rank_coeff;}      

        return_val = penalty;
        grid_trade[i*gridRes+j] = return_val;   

        if(return_val > max)  //keep new values ---------
        {
          max = return_val; _lambda = t_lambda; _alpha = t_alpha; 
          System.out.println(max); 
        }
        t_alpha = t_alpha + d_alpha;

        //gridprogressBar (i*50+j)
        System.out.println(i + "  " + j);
        gridprogressBar.setValue(i*gridRes+j);
      }
      t_lambda = t_lambda + d_lambda;
   }


   System.out.println("minimum found at " + _lambda + "  " + _alpha);
   fin_lambdaText.setText(""+df.format(_lambda)); 
   fin_alphaText.setText(""+df.format(_alpha)); 
   //gridContourPanel.updateFM(grid_trade, 50, 50);
   gridContourPanel.setData(grid_trade, gridRes, gridRes);

   optScale.updateColor(gridContourPanel.colorArray);
   opt_valText.setText(df.format(opt_value));

   lambdaBar.setValue((int)(_lambda*10));
   expBar.setValue((int)(_alpha*10)); 



   //setLambdaExp(_lambda,_alpha); 
 
                 

 }  

 

 public void gridSearchTradingFast()
 {

  int i,ss;   
  double alpha;
  gridRes = 50;
  gridprogressBar.setMaximum(gridRes);  

  int choice = optcritCombo.getSelectedIndex();  

 
  if(shortCheck.isSelected()) {ss = 1;}
  else {ss = 0;}
  

  int fl = n_obs - (L-1);
  double[] _price = new double[fl];
  grid_trade = new double[gridRes*gridRes];

  for(i=0;i<fl;i++)
  {_price[i] = t_price[0][L2+(L-1)+i];}

  int spc = 0; 
  if(autoUpdatesCheck.isSelected()) {spc=1;} 
 
  double[] output = mdfa.gridSearch(mdfa.tseries, tfilter.Gamma, n_obs, n_rep, mdfa.K, L, shift, lambda, expweight, 
                             w0, w1, mdfa.Lag, i1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost,mdfa.mod,mdfa.arg,spc);

  lambda = output[0]; 
  alpha = output[1];
  opt_value = output[2];
 
  //System.out.println("fl = " + fl + ", fl2 = " + output[3]);
  
  t_signal = new double[fl]; 
  grid_trade = new double[output.length-4-fl];
  for(i=0;i<output.length-4-fl;i++)
  {grid_trade[i] = output[4+i];}
  for(i=0;i<fl;i++) {t_signal[i] = output[i+4+grid_trade.length];}  
  
  
  
  //System.out.println("minimum found at " + lambda + "  " + alpha);
  fin_lambdaText.setText(""+df.format(lambda)); 
  fin_alphaText.setText(""+df.format(alpha)); 

  //gridRes = (int)Math.sqrt(output.length-3);
  gridContourPanel.setCutoffMode(false,grid_trade);
  gridContourPanel.setData(grid_trade, gridRes, gridRes);
  optScale.updateColor(gridContourPanel.colorArray);
  opt_valText.setText(df.format(opt_value));

  lambdaBar.setValue((int)(lambda*10));
  expBar.setValue((int)(alpha*10));  

 } 


 public void gridSearchTradingFast_Cutoff()
 {

  int i,ss;   
  int choice = optcritCombo.getSelectedIndex();  

  if(shortCheck.isSelected()) {ss = 1;}
  else {ss = 0;}
  

  int fl = n_obs - (L-1);
  double[] _price = new double[fl];
           

  for(i=0;i<fl;i++)
  {_price[i] = t_price[0][L2+(L-1)+i];}

  if(bandBut.isSelected())
  {
    double[] output = mdfa.gridSearchBandCutoff(mdfa.tseries, n_obs, n_rep, mdfa.K, L, shift, lambda, expweight, 
                             0.0, w1-w0, mdfa.Lag, i1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost, mdfa.mod, mdfa.arg, 0);
  
    
    w0 = output[0]; w1 = output[1]; opt_value = output[2]; 
    if(fl != (int)output[3]) {System.out.println("somehting is wrong with length");} 

    t_signal = new double[fl];   
    System.out.println("w0 = " + w0 + " w1 = " + w1);
    grid_trade = new double[output.length-4-fl];
    for(i=0;i<output.length-4-fl;i++)
    {grid_trade[i] = output[4+i];}

    for(i=0;i<fl;i++) {t_signal[i] = output[i+4+grid_trade.length];}
    
    computeFilter = false;
    
    omega1Text.setText(""+df.format(w1));
    omega1Bar.setValue((int)(w1*100));
    omega0Text.setText(""+df.format(w0));
    omega0Bar.setValue((int)(w0*100));    
    
    computeFilter = true;
    w0 = output[0]; w1 = output[1]; cutoff = w1; cutoff0 = w0;
    tfilter.setBand(w0,w1);
    setFilter(); computeFilterNow();    
    
    gridContourPanel.setCutoffMode(true,grid_trade); 
    
    
    
  }
  else
  {
    double[] output = mdfa.gridSearchCutoff(mdfa.tseries, n_obs, n_rep, mdfa.K, L, shift, lambda, expweight, 
                             w0, w1, mdfa.Lag, i1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost, mdfa.mod, mdfa.arg, 0);

    
    t_signal = new double[fl];  
    opt_value = output[1];
 
    grid_trade = new double[output.length-3-fl];
    for(i=0;i<output.length-3-fl;i++)
    {grid_trade[i] = output[4+i];}

    for(i=0;i<fl;i++) {t_signal[i] = output[i+3+grid_trade.length];}
 
    w1 = output[0]; 
    System.out.println("w0 = " + w0 + " w1 = " + w1);
    computeFilter = false;
    
    omega1Text.setText(""+df.format(w1));
    omega1Bar.setValue((int)(w1*100));
     
    computeFilter = true;
    w1 = output[0];
    tfilter.setBand(w0,w1);
    setFilter(); computeFilterNow();    
    
    gridContourPanel.setCutoffMode(true,grid_trade);  
  }
 
 } 
 
 
 public void out_of_sample_Optimization()
 {
 
    //extract the entire sequence of data
  int N,_n_rep,i,j,ss;
  int choice = optcritCombo.getSelectedIndex(); 
  int variable = optVariablesCombo.getSelectedIndex(); 
  int reg_var;
  gridRes = 50;
  true_out = true_outBar.getValue();    
  
  //System.out.println("true_out = " + true_out);
  
  int price_length = t_price[0].length-(L-1) - true_out;
  
  if(shortCheck.isSelected()) {ss = 1;}
  else {ss = 0;}

  int spc = 0; 
  if(autoUpdatesCheck.isSelected()) {spc=1;}   
  
  int n_out  = l1Bar.getValue() - true_out;
  if(n_out >= 36 && reComp.isSelected())
  {
      
    double[] _price = new double[price_length];
           
    for(i=0;i<price_length;i++)
    {_price[i] = t_price[0][L-1+i];} 
 
    N = n_obs_fixed - true_out;
    double[] gdp = new double[N*n_rep];
    
    System.out.println("n_out = " + price_length + " N = " + N + " n_rep=  " + n_rep);


    //------------- GET DATA -------------------------------------------------------
    for(i=0;i<N;i++)
    {gdp[i] = this.target[i];}

    _n_rep = 1;
    for(j=0;j<series_max;j++)
    {
       if(expVariable[j].isSelected() && expVariable[j].isEnabled())
       {
          double[] temp = complete_data.get(j); 
          for(i=0;i<N;i++)
          {
            gdp[N*_n_rep + i] = temp[i]; //System.out.println(temp[i]);           
          }
          _n_rep++;
       }
    }    
    //--------------------------------------------------------------------------------
    
   if(variable == 0)
   {

      double[] output = mdfa.gridSearch(gdp, tfilter.Gamma, N, n_rep, K, L, shift, lambda, expweight, 
                             w0, w1, n_out, 1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost, mdfa.mod, mdfa.arg, spc);

      lambda = output[0]; expweight = output[1]; opt_value = output[2];
  
      grid_trade = new double[output.length-4-n_out];
      for(i=0;i<output.length-4-n_out;i++) {grid_trade[i] = output[4+i];}
      
      fin_lambdaText.setText(""+df.format(lambda)); fin_alphaText.setText(""+df.format(expweight)); 

      gridContourPanel.setCutoffMode(false,grid_trade); 
      gridContourPanel.setData(grid_trade, gridRes, gridRes);
      optScale.updateColor(gridContourPanel.colorArray);
      opt_valText.setText(df.format(opt_value));

      lambdaBar.setValue((int)(lambda*10));
      expBar.setValue((int)(expweight*10));  

   }
   else if(variable == 1) //cutoff optimization
   {
    
    double[] output = mdfa.gridSearchBandCutoff(gdp, N, n_rep, K, L, shift, lambda, expweight, 
                             w0, w1, n_out, 1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost, mdfa.mod, mdfa.arg, spc);
  

    gridContourPanel.setMinOmega(w1, 2.0*Math.PI/3.0);
    w0 = output[0]; w1 = output[1]; opt_value = output[2]; 
    if(n_out != (int)output[3]) {System.out.println("something is wrong with length");} 

    //t_signal = new double[n_out];   
    System.out.println("w0 = " + w0 + " w1 = " + w1 + ", value = " + opt_value);
    grid_trade = new double[output.length-4-n_out];
    for(i=0;i<output.length-4-n_out;i++) {grid_trade[i] = output[4+i];}
 
    computeFilter = false;
    omega1Text.setText(""+df.format(w1));
    omega1Bar.setValue((int)(w1*100));
    omega0Text.setText(""+df.format(w0));
    omega0Bar.setValue((int)(w0*100));    
    
    computeFilter = true;
    w0 = output[0]; w1 = output[1]; cutoff = w1; cutoff0 = w0;
    tfilter.setBand(w0,w1);
    setFilter(); computeFilterNow();    
    
    gridContourPanel.setCutoffMode(true,grid_trade); 
  }
  else if(variable >= 2 && variable < 6) //regularization
  {
    reg_var = variable - 2; 
    double[] output = mdfa.gridSearchRegularization(gdp, tfilter.Gamma, N, n_rep, K, L, shift, lambda, expweight, 
                             w0, w1, n_out, 1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost,reg_var, mdfa.mod, mdfa.arg, spc);

    gridContourPanel.setMinOmega(0.0, 1.0);

    if(reg_var == 0)
    {smooth = output[0]; smoothBar.setValue((int)(smooth*1000));}
    else if(reg_var == 1)
    {decay = output[1]; decayBar.setValue((int)(decay*1000));}
    else if(reg_var == 2)
    {decay2 = output[2]; decay2Bar.setValue((int)(decay2*1000));}
    else if(reg_var == 3)
    {cross = output[3]; crossBar.setValue((int)(cross*1000));}

    grid_trade = new double[output.length-4-n_out];
    for(i=0;i<output.length-4-n_out;i++) {grid_trade[i] = output[4+i];}

    gridContourPanel.setCutoffMode(true,grid_trade); 
    //computeFilter = true;
    //computeFilterNow(); 

  }
  else if(variable == 6) //multibandpass filter optimization 
  {
  
  
      double[] output = mdfa.gridSearchMultiBandCutoff(gdp, N, n_rep, K, L, shift, lambda, expweight, 
                             w0, w1, n_out, 1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost, mdfa.mod, mdfa.arg, spc);
  
      w2 = output[0]; w3 = output[1]; opt_value = output[2]; 
      if(n_out != (int)output[3]) {System.out.println("something is wrong with length");} 

      gridContourPanel.setMinOmega(0.0, 1.0);
      //t_signal = new double[n_out];   
      System.out.println("w2 = " + w2 + " w3 = " + w3 + ", value = " + opt_value);
      grid_trade = new double[output.length-4-n_out];
      for(i=0;i<output.length-4-n_out;i++) {grid_trade[i] = output[4+i];}
 
      computeFilter = false;
      omega2Text.setText(""+df.format(w2)); omega2Bar.setValue((int)(w2*100));
      omega3Text.setText(""+df.format(w3)); omega3Bar.setValue((int)(w3*100));    
    
      computeFilter = true;
      tfilter.setBand2(w2,w3);
      setFilter(); computeFilterNow();     
      gridContourPanel.setCutoffMode(true,grid_trade);    
  }
  else if(variable == 7)
  {
  
    double[] output = mdfa.optimizeRegularization(gdp, tfilter.Gamma, N, n_rep, K, L, shift, lambda, expweight, 
                             w0, w1, n_out, 1, i2, smooth, decay, decay2, cross, _price, choice, ss, trading_cost,0);
  
    
    computeFilter = false;
    smooth = output[0]; smoothBar.setValue((int)(smooth*1000));
    decay = output[1]; decayBar.setValue((int)(decay*1000));
    decay2 = output[2]; decay2Bar.setValue((int)(decay2*1000));
    cross = output[3]; crossBar.setValue((int)(cross*1000));
    computeFilter = true;
    computeFilterNow();
  
  }

 }
 else
 {System.out.println("Check reComputeFilter and number of out_of_sample points first");}
 
 
 }
 
 
 
 
 



  private void initTradingParameterDialog() 
  {
    tradingparameterPanel = new JPanel();
    // Variables declaration - do not modify//GEN-BEGIN:variables
    JLabel jLabel1,jLabel13,jLabel2,jLabel3,jLabel4;
    JPanel jPanel1 = new JPanel();

    deltaLabel = new JLabel("Trading Threshold"); deltaLabel.setHorizontalTextPosition(JLabel.LEFT);
    shortCheck = new JCheckBox(); riskfreeText = new JTextField(); riskfreeSlider = new JSlider(JSlider.HORIZONTAL,0,100,0);
    jLabel1 = new JLabel(); jLabel13 = new JLabel(); jLabel4 = new JLabel();jLabel3 = new JLabel();
    jLabel2 = new JLabel(); tradingfreqSlider = new JSlider(JSlider.HORIZONTAL,0,100,50);  tradingcostText = new JTextField();
    tradingcostSlider = new JSlider(JSlider.HORIZONTAL,0,10,0); jPanel1.setBorder(BorderFactory.createTitledBorder("Trading Parameters"));
    shortCheck.setFont(new Font("Ubuntu", 0, 12)); shortCheck.setText("Short Sell");
    diffSigCheck = new JCheckBox("Diff Signal"); diffSigCheck.setEnabled(true); diffSigCheck.setSelected(false);
    
    
    shortCheck.addItemListener(new ItemListener() {
       public void itemStateChanged(ItemEvent e)
       {
         boolean sel; 
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         short_sell = sel;
         updatePlots(false,false);
       }
     }
    );
    
    diffSigCheck.addItemListener(new ItemListener() {
       public void itemStateChanged(ItemEvent e)
       {
         boolean sel; 
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         toggleDiffSigTrading(sel);
         updatePlots(false,false);
       }
     }
    );

    riskfreeText.setText("0");
    riskfreeSlider.setFont(new Font("Ubuntu", 0, 3)); // NOI18N
    jLabel1.setFont(new Font("Ubuntu", 0, 12)); // NOI18N     
    jLabel1.setText("Trading Frequency:");
    jLabel13.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
    jLabel13.setText("Risk Free Rate:");

    jLabel4.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
    jLabel4.setText("Trading Costs:");
    jLabel3.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
    jLabel3.setText("More");
    jLabel2.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
    jLabel2.setText("Less");
    tradingfreqSlider.setFont(new Font("Ubuntu", 0, 3)); // NOI18N
    tradingcostText.setText("0");
    tradingcostSlider.setFont(new Font("Ubuntu", 0, 3)); // NOI18N

    tradingfreqSlider.addChangeListener(new MyChangeListener2());
    tradingcostSlider.addChangeListener(new MyChangeListener2());
    riskfreeSlider.addChangeListener(new MyChangeListener2());


        GroupLayout jPanel1Layout = new GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel1)
                    .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                        .addComponent(jLabel13)
                        .addComponent(jLabel4)
                        .addComponent(deltaLabel))
                    .addComponent(shortCheck)
                    .addComponent(diffSigCheck))
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel2)
                        .addGap(6, 6, 6)
                        .addComponent(tradingfreqSlider, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE))
                    .addGroup(GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 43, Short.MAX_VALUE)
                        .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(tradingcostSlider, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE)
                            .addComponent(riskfreeSlider, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE)
                            .addComponent(deltaBar, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE))))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel3)
                    .addComponent(riskfreeText, GroupLayout.DEFAULT_SIZE, 45, Short.MAX_VALUE)
                    .addComponent(tradingcostText)
                    .addComponent(deltaText))
                .addContainerGap(54, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(tradingfreqSlider, GroupLayout.PREFERRED_SIZE, 33, GroupLayout.PREFERRED_SIZE)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                                .addComponent(jLabel2)
                                .addComponent(jLabel1))
                            .addComponent(jLabel3))))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                    .addComponent(tradingcostText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(tradingcostSlider, GroupLayout.PREFERRED_SIZE, 33, GroupLayout.PREFERRED_SIZE)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel4)
                        .addGap(9, 9, 9)))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                        .addComponent(riskfreeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(riskfreeSlider, GroupLayout.PREFERRED_SIZE, 33, GroupLayout.PREFERRED_SIZE))
                    .addComponent(jLabel13))
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                        .addComponent(deltaText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(deltaBar, GroupLayout.PREFERRED_SIZE, 33, GroupLayout.PREFERRED_SIZE))
                    .addComponent(deltaLabel))    
                .addGap(18, 18, 18)
                .addComponent(shortCheck)
                .addComponent(diffSigCheck)
                .addContainerGap(27, Short.MAX_VALUE))
        );

        GroupLayout layout = new GroupLayout(tradingparameterPanel);
        tradingparameterPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel1, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel1, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
  }// </editor-fold>//GEN-END:initComponents
  
 
  @SuppressWarnings({ "unchecked", "rawtypes" })
private void initOptimizationDialog() 
  {
        tradingoptimPanel = new JPanel();
        JLabel jLabel1,jLabel2,jLabel3;


        optVariablesCombo = new JComboBox();
        optVariablesCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        optVariablesCombo.setModel(new DefaultComboBoxModel(new String[] 
        { "Customization", "Bandpass Frequency", "Regularization-Smooth", "Regularization-Decay", "Regularization-Decay2", "Regularization-Cross", "MultiBand Pass", "Regularization" }));


        cutOptimizeBox = new JCheckBox("Optimize Bandwidth Only");
        cutOptimizeBox.setSelected(false);
        
        outSampOptBox = new JCheckBox("Out-of-sample Optimization");
        outSampOptBox.setSelected(false);

        JPanel cutOptimizeBoxSuit = new JPanel();
        cutOptimizeBoxSuit.add(outSampOptBox);
        cutOptimizeBoxSuit.add(optVariablesCombo);
        cutOptimizeBoxSuit.setLayout(new GridLayout(1,2,2,0)); 
        
        
        
        optScale = new ScalePanel();
        optimtradePanel = new JPanel();
        simannButton = new JButton();
        gridSearchButton = new JButton();
        jLabel1 = new JLabel();
        optcritCombo = new JComboBox();
        
        gridprogressBar = new JProgressBar(0,50*50);
        jLabel2 = new JLabel();
        jLabel3 = new JLabel();
        fin_alphaText = new JTextField();
        fin_lambdaText = new JTextField();

        gridprogressBar.setValue(0);
        gridprogressBar.setStringPainted(true);

        optimtradePanel.setBorder(BorderFactory.createTitledBorder("Financial Optimization Panel"));

        simannButton.setText("Simulated Annealing");
        simannButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent evt) {                
              gridSearchTrading();              
            }
        });

        gridSearchButton.setText("Grid Search");

        gridSearchButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent evt) {                
              
              if(outSampOptBox.isSelected())
              {
                out_of_sample_Optimization();
              }
              else
              {
                if(optVariablesCombo.getSelectedIndex() == 0)
                {gridSearchTradingFast();}
                else
                {gridSearchTradingFast_Cutoff();}
              }              
           }
        });

        jLabel1.setText("In-Sample Financial Optimization Criterion");

        optcritCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        optcritCombo.setModel(new DefaultComboBoxModel(new String[] { "Maximize Return", "Minimize Max Loss", "Min-Max Method", "Maximize Trade Success Ratio", "Max Return*Succ Ratio", "Maximize Rank Coefficient" }));

        gridContourPanel.setBackground(new Color(1, 1, 1));
        gridContourPanel.setForeground(new Color(1, 1, 1));

        GroupLayout gridContourPanelLayout = new GroupLayout(gridContourPanel);
        gridContourPanel.setLayout(gridContourPanelLayout);
        gridContourPanelLayout.setHorizontalGroup(
            gridContourPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
        gridContourPanelLayout.setVerticalGroup(
            gridContourPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 210, Short.MAX_VALUE)
        );

        optScale.setBackground(new java.awt.Color(1, 1, 1));

        javax.swing.GroupLayout optScaleLayout = new javax.swing.GroupLayout(optScale);
        optScale.setLayout(optScaleLayout);
        optScaleLayout.setHorizontalGroup(
            optScaleLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 20, Short.MAX_VALUE)
        );
        optScaleLayout.setVerticalGroup(
            optScaleLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        opt_valText = new JTextField();
        opt_valText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        opt_valText.setText("0.0");       


        jLabel2.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel2.setText("\u03B1");

        jLabel3.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel3.setText("\u03BB");

        fin_alphaText.setText("0.0");

        fin_lambdaText.setText("0.0");

        GroupLayout optimtradePanelLayout = new GroupLayout(optimtradePanel);
        optimtradePanel.setLayout(optimtradePanelLayout);
        optimtradePanelLayout.setHorizontalGroup(
            optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optimtradePanelLayout.createSequentialGroup()
                .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(optimtradePanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(cutOptimizeBoxSuit, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(optimtradePanelLayout.createSequentialGroup()
                                .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(optimtradePanelLayout.createSequentialGroup()
                                        .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                            .addComponent(fin_alphaText, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(jLabel2, javax.swing.GroupLayout.Alignment.LEADING))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(gridContourPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(optScale, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(opt_valText, javax.swing.GroupLayout.PREFERRED_SIZE, 38, javax.swing.GroupLayout.PREFERRED_SIZE))
                                    .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addGroup(optimtradePanelLayout.createSequentialGroup()
                                            .addComponent(simannButton, javax.swing.GroupLayout.PREFERRED_SIZE, 190, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addGap(37, 37, 37)
                                            .addComponent(gridSearchButton, javax.swing.GroupLayout.PREFERRED_SIZE, 190, javax.swing.GroupLayout.PREFERRED_SIZE))
                                        .addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 299, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(optcritCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 299, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                .addGap(0, 48, Short.MAX_VALUE))))
                    .addGroup(optimtradePanelLayout.createSequentialGroup()
                        .addGap(168, 168, 168)
                        .addComponent(jLabel3)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fin_lambdaText, javax.swing.GroupLayout.PREFERRED_SIZE, 55, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
       optimtradePanelLayout.setVerticalGroup(
            optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optimtradePanelLayout.createSequentialGroup()
                .addGap(30, 30, 30)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(optcritCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(37, 37, 37)
                .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(simannButton)
                    .addComponent(gridSearchButton))
                .addGap(18, 18, 18)
                .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(gridContourPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(optimtradePanelLayout.createSequentialGroup()
                        .addGap(2, 2, 2)
                        .addComponent(opt_valText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(64, 64, 64)
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fin_alphaText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(optScale, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optimtradePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(fin_lambdaText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 18, Short.MAX_VALUE)
                .addComponent(cutOptimizeBoxSuit, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        GroupLayout layout = new GroupLayout(tradingoptimPanel);
        tradingoptimPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(optimtradePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(optimtradePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents


  public void initOutTradingStatDialog()
  {

        JLabel jLabel10,jLabel11,jLabel12,jLabel5,jLabel6,jLabel7,jLabel8,jLabel9;


        outsampStatPanel = new JPanel();
        outstatPanel = new javax.swing.JPanel();
        jLabel7 = new javax.swing.JLabel();
        jLabel8 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jLabel9 = new javax.swing.JLabel();
        outsharpeText = new javax.swing.JTextField();
        outavgText = new javax.swing.JTextField();
        jLabel5 = new javax.swing.JLabel();
        jLabel10 = new javax.swing.JLabel();
        outroiText = new javax.swing.JTextField();
        outratioText = new javax.swing.JTextField();
        jLabel11 = new javax.swing.JLabel();
        outsuccText = new javax.swing.JTextField();
        outtotalText = new javax.swing.JTextField();
        outForeText = new JTextField("0, 0");
        outndropsText = new javax.swing.JTextField();
        outmaxdropText = new javax.swing.JTextField();
        jLabel12 = new javax.swing.JLabel();

        outstatPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Out-of-Sample Trading Statistics"));
        foreLabel = new JLabel("ForeMSE");
        foreLabel.setFont(new java.awt.Font("Ubuntu", 0, 12));
        rrLabel = new JLabel("Risk/Reward");
        rrLabel.setFont(new java.awt.Font("Ubuntu", 0, 12));
        outRRText = new JTextField("0");
        jLabel7.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel7.setText("Max Drop:");

        jLabel8.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel8.setText("Sign Corr:");

        jLabel6.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel6.setText("Sharpe Ratio:");

        jLabel9.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel9.setText("Rank Coeff:");

        outsharpeText.setText("0.0");

        outavgText.setText("0.0");

        jLabel5.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel5.setText("ROI:");

        jLabel10.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel10.setText("Succ %");

        outroiText.setText("0.0");

        outratioText.setText("0.0");

        jLabel11.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel11.setText("Succ Trades:");

        outsuccText.setText("0.0");

        outtotalText.setText("0.0");

        outndropsText.setText("0.0");

        outmaxdropText.setText("0.0");

        jLabel12.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel12.setText("Num Trades:");

        javax.swing.GroupLayout statPanelLayout = new javax.swing.GroupLayout(outstatPanel);
        outstatPanel.setLayout(statPanelLayout);
        statPanelLayout.setHorizontalGroup(
            statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(statPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(rrLabel)
                    .addComponent(jLabel7)
                    .addComponent(jLabel8)
                    .addComponent(jLabel6)
                    .addComponent(jLabel5))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(outRRText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outmaxdropText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outndropsText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outsharpeText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outroiText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(26, 26, 26)
                .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel9)
                    .addComponent(jLabel10)
                    .addComponent(jLabel11)
                    .addComponent(jLabel12)
                    .addComponent(foreLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(outtotalText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outsuccText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outavgText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outratioText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(outForeText, javax.swing.GroupLayout.PREFERRED_SIZE, 74, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(33, Short.MAX_VALUE))
        );
        statPanelLayout.setVerticalGroup(
            statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(statPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(statPanelLayout.createSequentialGroup()
                        .addComponent(outtotalText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel11)
                            .addComponent(outsuccText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel10)
                            .addComponent(outratioText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(foreLabel)
                            .addComponent(outForeText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel9)
                            .addComponent(outavgText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                   
                   .addGroup(statPanelLayout.createSequentialGroup()
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel5)
                            .addComponent(outroiText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel12))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel6)
                            .addComponent(outsharpeText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel7)
                            .addComponent(outmaxdropText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel8)
                            .addComponent(outndropsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(statPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(rrLabel)
                            .addComponent(outRRText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                             
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(outsampStatPanel);
        outsampStatPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(outstatPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(outstatPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
  }// </editor-fold>
    // Variables declaration - do not modify




  public void initTradingStatDialog() 
  {
        JLabel jLabel10, jLabel11,jLabel12,jLabel5,jLabel6,jLabel7,jLabel8,jLabel9;


        tradingstatPanel = new JPanel();
        statPanel = new JPanel();
        jLabel7 = new JLabel(); jLabel8 = new JLabel();
        jLabel6 = new JLabel();
        jLabel9 = new JLabel();
        sharpeText = new JTextField();
        avgText = new JTextField();
        jLabel5 = new JLabel();
        jLabel10 = new JLabel();
        roiText = new JTextField();
        ratioText = new JTextField();
        jLabel11 = new JLabel();
        succText = new JTextField();
        totalText = new JTextField();
        ndropsText = new JTextField();
        maxdropText = new JTextField();
        jLabel12 = new JLabel();

        statPanel.setBorder(BorderFactory.createTitledBorder("Trading Statistics"));

        jLabel7.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel7.setText("Max Drop:");

        jLabel8.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel8.setText("Sign Corr:");

        jLabel6.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel6.setText("Sharpe Ratio:");

        jLabel9.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel9.setText("Risk Coeff:");

        sharpeText.setText("0.0");

        avgText.setText("0.0");

        jLabel5.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel5.setText("ROI:");

        jLabel10.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel10.setText("Succ %");

        roiText.setText("0.0");

        ratioText.setText("0.0");

        jLabel11.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel11.setText("Succ Trades:");

        succText.setText("0.0");

        totalText.setText("0.0");

        ndropsText.setText("0.0");

        maxdropText.setText("0.0");

        jLabel12.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        jLabel12.setText("Num Trades:");

        GroupLayout statPanelLayout = new GroupLayout(statPanel);
        statPanel.setLayout(statPanelLayout);
        statPanelLayout.setHorizontalGroup(
            statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(statPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel7)
                    .addComponent(jLabel8)
                    .addComponent(jLabel6)
                    .addComponent(jLabel5))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(maxdropText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(ndropsText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(sharpeText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(roiText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE))
                .addGap(26, 26, 26)
                .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel9)
                    .addComponent(jLabel10)
                    .addComponent(jLabel11)
                    .addComponent(jLabel12))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(totalText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(succText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(avgText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE)
                    .addComponent(ratioText, GroupLayout.PREFERRED_SIZE, 74, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(33, Short.MAX_VALUE))
        );
        statPanelLayout.setVerticalGroup(
            statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(statPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(statPanelLayout.createSequentialGroup()
                        .addComponent(totalText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel11)
                            .addComponent(succText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel10)
                            .addComponent(ratioText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel9)
                            .addComponent(avgText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                    .addGroup(statPanelLayout.createSequentialGroup()
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel5)
                            .addComponent(roiText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel12))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel6)
                            .addComponent(sharpeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel7)
                            .addComponent(maxdropText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(statPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel8)
                            .addComponent(ndropsText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        GroupLayout layout = new GroupLayout(tradingstatPanel);
        tradingstatPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(statPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(statPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
    }



    public void loadFilterCoeffs(File file)
    {

         Double D; int i,k;
         String strline; String[] tokens; String delims = "[ ]+";  
         double[] vals;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
          
           strline = br.readLine(); //get L and nreps
           tokens = strline.split(delims);
          
           int tempL = (new Integer(tokens[0])).intValue();
           int temp_nreps = (new Integer(tokens[1])).intValue();
            
           System.out.println(tempL + "  " + temp_nreps); 
            
           if(tempL < n_obs && temp_nreps == n_rep)
           {
             
             LBar.setValue(tempL); L = tempL;
             
             vals = new double[n_rep*L];
             k=0;
             while((strline = br.readLine()) != null)
             {
                tokens = strline.split(delims);
                
                for(i=0;i<tokens.length;i++)
                {
                 D = new Double(tokens[i]);
                 vals[L*i + k] = D.doubleValue();
                }
                k++;
             }
             
             if(k != L) {System.out.println("Problems in filter length, k = "+k);}
             else
             {
              mdfa.b = new double[n_rep*L];
              for(i=0; i < n_rep; i++)  
              {
               for(k=0;k<L;k++)
               {mdfa.b[L*i + k] = vals[L*i + k];}
              }
              
              mdfa.saveBfilter();            //save filter              
              reComp.setSelected(false);     //turn off recompute
              applyNewFilter();
                            
             }
           }
           br.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
         catch(IOException ioe){System.out.println("IO out error..." + ioe);}
   
    }


    public void loadPriceFilterCoeffs(File file)
    {
    
         Double D; int k;
         String strline; String[] tokens; String delims = "[ ]+";  
         double[] vals;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
          
           strline = br.readLine(); //get L and nreps
           tokens = strline.split(delims);
          
          
           int tempL = (new Integer(tokens[0])).intValue();
           int temp_nreps = (new Integer(tokens[1])).intValue();
            
           System.out.println(tempL + "  " + temp_nreps); 
            
           if(tempL < L && temp_nreps == 1)
           {
                          
             vals = new double[tempL];
             k=0;
             while((strline = br.readLine()) != null)
             {
                tokens = strline.split(delims);
                D = new Double(tokens[0]);
                vals[k] = D.doubleValue();
                k++;
             }
             
             mdfa.setPriceFilter(vals);
           }
           br.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
         catch(IOException ioe){System.out.println("IO out error..." + ioe);}
   
    }






    public void saveFilter()
    {
        MDFAFilter filter = new MDFAFilter(n_obs, n_rep, L, Lag);

        filter.w_const = new double[w_const.length]; 
        filter.lambda = lambda;  filter.expweight=expweight; filter.lambda_3=lambda_3;   
        filter.smooth=smooth; filter.decay=decay; filter.cross=cross;
        filter.i1 = i1; filter.i2 = i2; filter.dd = dd; filter.DD = DD;
        filter.criteria=mdfa.criteria; filter.degrees=mdfa.degrees; filter.MDFAmin = mdfa.MDFAmin;
        filter.b = new double[mdfa.b.length]; 
        filter.dfa = zpc_gene;
        System.arraycopy(mdfa.xf, 0, filter.xf, 0, n_obs-(L-1)); 
        System.arraycopy(mdfa.Gamma, 0, filter.Gamma, 0, K1);

        //---- frequency --------------------------
        filter.cutoff=mdfa.cutoff; filter.cutoff0=mdfa.cutoff0;
        System.arraycopy(w_const, 0, filter.w_const, 0, w_const.length);        
        System.arraycopy(mdfa.b, 0, filter.b,0,mdfa.b.length);  

        savedFilters.add(filter);
         
    }

    public void clearFilters()
    {savedFilters.clear();}

    public void set_iter(boolean m)
    {mdfa.set_dfaIter(m); iter = m; mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);}

    /*public void set_reg(boolean m)
    {
      mdfa.setReg(m); mdfa.reset_lengths(); reg = m; iter = !m; 
      mdfa.computeFilterGeneral(computeFilter);  updatePlots(false,true);
      smoothBar.setEnabled(m); decayBar.setEnabled(m); crossBar.setEnabled(m); iterCheck.setSelected(false);
    }*/

    //---------- Various Signal Extraction Polynomials --------------------------

    public void setTrendPoly(double[] trend, double le, double var)
    { 
      ntrend = (int)le; trendpoly = new double[ntrend]; trendvar = var; 
      System.arraycopy(trend, 0, trendpoly, 0, ntrend);
      trendpoly[0] = 1.0;
    }
    public void setSeasPoly(double[] trend, double le, double var)
    { 
      nseas = (int)le; seaspoly = new double[nseas]; seasvar = var; 
      System.arraycopy(trend, 0, seaspoly, 0, nseas);
      seaspoly[0] = 1.0;
    }
    public void setTIPoly(double[] trend, double le, double var)
    { 
      ntip = (int)le; tipoly = new double[ntip]; tivar = var; 
      System.arraycopy(trend, 0, tipoly, 0, ntip);
      tipoly[0] = 1.0;
    } 
    
    public void setModelPoly(Polynomial _mphi, Polynomial _mPhi, Polynomial _mapoly, double _m)
    {mphi = _mphi; mPhi = _mPhi; mapoly = _mapoly; m_innvar = _m;}

    //------- Compute the X13 target filters ---------------------------------
    //    model = 0 (trend), model = 1 (seasonal), model = 2 (trend-irregular)
    //------------------------------------------------------------------------     

    public void computeX13TargetFilter()
    {
        double sum, sum1, sum2, sum3;
        double isum, isum1, isum2, isum3; 
        int i,j,madeg;
        int nsamp = K;
        double sampleFw;  //--model spectral density
        double sampleDn;  //--noise differencing
        double sampleSn;  //--signal polynomial 
        double seasD, nseasD;
        double f2wn;
        double sigvar=0.0;
        double ar; 
        double maxval = -1000;
        madeg = mapoly.deg;

        double[] ar_params = new double[mphi.deg];
        double[] ma_params = new double[madeg];

        if(useSARIMA)
        {
           seasD=0.0; nseasD=0.0;
           Gamma_x13 = new double[K1];
           f_specDens = new double[K1];

           if(x13model == 0) {sigvar = trendvar; }//System.out.println(sigvar + "  " + ntrend);} 
           else if(x13model == 1) {sigvar = seasvar;}// System.out.println(sigvar + "  " + nseas);} 
           else {sigvar = tivar;}// System.out.println(sigvar + "  " + ntip);} 

          for(i=0; i < K1; i++)
          {    
           sum = 0; sum1 = 0; sum2 = 0; sum3 = 0; isum = 0; isum1 = 0; isum2 = 0; isum3 = 0;

           //------------ AR components -----------------------           
           for(j=0;j <= mphi.deg; j++) 
           {
             sum3 = sum3 + mphi.coef[j]*Math.cos(j*(i*Math.PI)/nsamp); 
             isum3 = isum3 + mphi.coef[j]*Math.sin(j*(i*Math.PI)/nsamp);
           }     
           for(j=0;j <= mPhi.deg; j++) 
           {
	     sum1 = sum1 + mPhi.coef[j]*Math.cos(j*(i*Math.PI)/nsamp);
	     isum1 = isum1 + mPhi.coef[j]*Math.sin(j*(i*Math.PI)/nsamp);
           }
                
           ar = (sum3*sum3 + isum3*isum3)*(sum1*sum1 + isum1*isum1);

           //sum3 = 0; isum3 = 0;
           // ------------  seasonal differencing-------------
           if(x13model == 0 || x13model == 2)
           {
	    for(j=0; j < 12; j++) 
	    {sum = sum + Math.cos(j*(i*Math.PI)/nsamp); isum = isum + Math.sin(j*(i*Math.PI)/nsamp);} 

            f2wn = sum3*sum3 + isum3*isum3; 

	    //sum3 = 1 - 2*Math.cos((i*Math.PI)/nsamp) + Math.cos(2*(i*Math.PI)/nsamp); 
	    //isum3 = -2*Math.sin((i*Math.PI)/nsamp) + Math.sin(2*(i*Math.PI)/nsamp);

 	    //ar = (sum*sum + isum*isum)*(sum3*sum3 + isum3*isum3);
            // (1-B)(1-B^12) = 1 - B^12 - B + B^12
           }
           else // ------- trend differencing -------------
           {
	    sum = 1 - 2*Math.cos((i*Math.PI)/nsamp) + Math.cos(2*(i*Math.PI)/nsamp); 
	    isum = -2*Math.sin((i*Math.PI)/nsamp) + Math.sin(2*(i*Math.PI)/nsamp);

            f2wn = sum1*sum1 + isum1*isum1; 

           }
           sampleDn = sum*sum + isum*isum;

           //---------- f^2_W_MA ------------
           sum2 = 0; isum2 = 0;              
           for(j=0; j <= madeg; j++)
           {            
            sum2 = sum2 + mapoly.coef[j]*Math.cos(j*(i*Math.PI)/nsamp);
            isum2 = isum2 + mapoly.coef[j]*Math.sin(j*(i*Math.PI)/nsamp);
           } 
           sampleFw = (sum2*sum2 + isum2*isum2)*m_innvar/f2wn; 


           seasD = 1 - Math.cos((i*Math.PI)/nsamp) - Math.cos((12.0*i*Math.PI)/nsamp) + Math.cos((13.0*i*Math.PI)/nsamp); 
           nseasD = Math.sin((i*Math.PI)/nsamp) - Math.sin((12.0*i*Math.PI)/nsamp) + Math.sin((13.0*i*Math.PI)/nsamp); 

           //----- spectral density-------------------------------
           f_specDens[i] = (sum2*sum2 + isum2*isum2)*m_innvar/(ar*(seasD*seasD + nseasD*nseasD));
           //-----------------------------------------------------
           

           sum2 = 0; isum2 = 0;   
           if(x13model == 0)  //trend
           {         
             for(j=0; j < ntrend; j++)
             {sum2 = sum2 + trendpoly[j]*Math.cos(j*(i*Math.PI)/nsamp); isum2 = isum2 + trendpoly[j]*Math.sin(j*(i*Math.PI)/nsamp);}  
           }
           else if(x13model == 1)  //seasonal
           {         
             for(j=0; j < nseas; j++)
             {sum2 = sum2 + seaspoly[j]*Math.cos(j*(i*Math.PI)/nsamp); isum2 = isum2 + seaspoly[j]*Math.sin(j*(i*Math.PI)/nsamp);}  
           }
           else if(x13model == 2) //trend-irr
           {         
             for(j=0; j < ntip; j++)
             {sum2 = sum2 + tipoly[j]*Math.cos(j*(i*Math.PI)/nsamp); isum2 = isum2 + tipoly[j]*Math.sin(j*(i*Math.PI)/nsamp);}  
           }
           sampleSn = sum2*sum2 + isum2*isum2;
           Gamma_x13[i] = sampleSn*sampleDn*sigvar/sampleFw;

           if(Gamma_x13[i] > maxval) {maxval = Gamma_x13[i];}
          }
          for(i=0;i<K1;i++) {Gamma_x13[i] = Gamma_x13[i]/maxval;}

          f_specDens[0] = 500;  f_specDens[K] = 500;

          if(spec_densBox.isSelected())
          {
            for(i=1;i<=madeg;i++)
            {ma_params[i-1] = mapoly.coef[i];}

            for(i=1;i<=mphi.deg;i++)
            {ar_params[i-1] = mphi.coef[i];}

            mdfa.setSpecDensityParams(ar_params, ma_params, m_innvar);
          }
        }
        else
        {System.out.println("SARIMA data import must be turned on to activate computation of this target filter");}

     
    }


    public void computeARMASpectralDensity(int p, int q)
    {

       f_specDens = new double[K1];
       int model = 1; int method = 1; int d = 0;
       int n_fsteps = 1; int nsims = 1;
       int MLE = 1; 

       double[] tser = new double[n_obs];
       System.arraycopy(mdfa.tseries,0,tser,0,n_obs);
       
       Cronos arma = new Cronos(n_obs, model, method);
       arma.setData(mdfa.tseries);
       arma.setNForecastSteps(n_fsteps);
       arma.setNPredictiveSims(n_fsteps, nsims); 
       arma.setARMA_Params(p,q,d); 
       arma.computeARIMAModel(MLE);

       System.out.println("Recomputing ARMA estimate");
       mdfa.setSpecDensityParams(arma.ar_params, arma.ma_params, arma.innvar);

       /*
       for(i=0;i<K1;i++)
       {

        sum2 = 1.0; isum2 = 0.0;
        for(j=1; j <= q; j++)
        {            
            sum2 = sum2 + arma.ma_params[j-1]*Math.cos(j*(i*Math.PI)/K);
            isum2 = isum2 + arma.ma_params[j-1]*Math.sin(j*(i*Math.PI)/K);
        } 
         
        sum_ma = (sum2*sum2 + isum2*isum2);

        sum2 = 1.0; isum2 = 0.0;
        for(j=1;j <= p; j++) 
        {
	  sum2 = sum2 + arma.ar_params[j-1]*Math.cos(j*(i*Math.PI)/K);
	  isum2 = isum2 + arma.ar_params[j-1]*Math.sin(j*(i*Math.PI)/K);
        }
        sum_ar = (sum2*sum2 + isum2*isum2);         
       
        f_specDens[i] = sum_ma*arma.innvar/sum_ar;
      }

      mdfa.setSpecDensity(f_specDens, 1);
      */
    }




    //----------------------------------------------------------------------------
   

    public void inputX13data(double[] _x)
    {
      x13nobs = _x.length; x13data = new double[x13nobs];
      System.arraycopy(_x, 0, x13data, 0, x13nobs);
    } 
    

    public void setupCheckButtons()
    {
      int i;
      /*----------------------------------------------
        Radio button group for chooMath.sing freq domain plots
      --------------------------------------------------*/ 
      gammaBut = new JRadioButton("\u0393(\u03C9):"); gammaBut.setHorizontalTextPosition(JMenuItem.LEFT);
      gammaBut.setSelected(true); gammaBut.setToolTipText("Plot the different frequency response functions of the filters.");
      amplitudeBut = new JRadioButton("A(\u03C9):"); amplitudeBut.setHorizontalTextPosition(JMenuItem.LEFT);
      amplitudeBut.setToolTipText("Plot the different amplitude functions of the filters.");
      timeDelayBut = new JRadioButton("\u03C6(\u03C9)/\u03C9:"); timeDelayBut.setHorizontalTextPosition(JMenuItem.LEFT);
      periodBut = new JRadioButton("I_N(\u03C9):"); timeDelayBut.setHorizontalTextPosition(JMenuItem.LEFT);  //--periodogram 
      timeDelayBut.setToolTipText("Plot the different phase delay functions of the filters.");
      periodBut.setToolTipText("Plot the different periodograms of the filters.");

      freqGroup = new ButtonGroup();
      freqGroup.add(gammaBut);
      freqGroup.add(amplitudeBut);
      freqGroup.add(timeDelayBut);
      freqGroup.add(periodBut);

      /*---------------------------------------------------
        Simulation type button group 
      ----------------------------------------------------*/
      sarimaCheck = new JRadioButton("SARIMA:"); sarimaCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      sarimaCheck.setToolTipText("Use the time series data imported from the SARIMA modeling and simulation module.");
      simCheck = new JRadioButton("GDP:"); simCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      simCheck.setToolTipText("Use up to 5 multi-time series data simulated from various AR and MA components.");
      noneCheck = new JRadioButton("None:");  noneCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      noneCheck.setToolTipText("Use up to 5 multi-time series data simulated from various AR and MA components.");
      
      simGroup = new ButtonGroup();
      simGroup.add(sarimaCheck);
      simGroup.add(simCheck);
      simGroup.add(noneCheck);
      simCheck.setSelected(true); simulate = true;  useSARIMA = false;




 
      //---------------------------------------------------
      
      /*---------------------------------------------------
      x13filter
      -----------------------------------------------------*/

      ButtonGroup x13Group = new ButtonGroup();
      x13filter = new JCheckBox("X13-Filters"); x13filter.setToolTipText("Build the X-13-ARIMA-SEATS target filters");
      x13trend = new JRadioButton("X13-Trend"); x13trend.setToolTipText("Import the X-13-ARIMA-SEATS trend signal target filter");      
      x13seas = new JRadioButton("X13-Seas");  x13seas.setToolTipText("Import the X-13-ARIMA-SEATS trend signal target filter");
      x13ti = new JRadioButton("X13-TI"); x13ti.setToolTipText("Import the X-13-ARIMA-SEATS trend-irregular signal target filter");
      x13Group.add(x13trend);  x13Group.add(x13seas);  x13Group.add(x13ti);
      x13trend.setSelected(true); x13model = 0;
      x13filter.setHorizontalTextPosition(JMenuItem.LEFT); x13filter.addItemListener(new MyItemListener());
      x13trend.setHorizontalTextPosition(JMenuItem.LEFT);  x13trend.addItemListener(new MyItemListener()); 
      x13seas.setHorizontalTextPosition(JMenuItem.LEFT);   x13seas.addItemListener(new MyItemListener());
      x13ti.setHorizontalTextPosition(JMenuItem.LEFT);     x13ti.addItemListener(new MyItemListener());
      setEnableX13(false); 

      plotHist = new JCheckBox("History: "); plotHist.setHorizontalTextPosition(JMenuItem.LEFT); plotHist.setSelected(false); 
      plotHist.setToolTipText("Plot the saved historical plots if any."); plotHist.addItemListener(new MyItemListener());

      /*------------------------------------------------------------
        Radio button group for chooMath.sing freq domain plots
      -------------------------------------------------------------*/ 
      expVariable = new JCheckBox[11];
      freqPlot = new JCheckBox[11];
      timePlot = new JCheckBox[12];
      periodPlot = new JCheckBox[12];
      coeffPlot = new JCheckBox[11];
      cointSliders = new JSlider[12];  // --- co-integration sliders for i1 at frequency 0

      plotXf = new JCheckBox(" \u0176(t):"); plotXf.setSelected(false); plotXf.setHorizontalTextPosition(JMenuItem.LEFT);
      plotGamma = new JCheckBox("\u0393(\u03C9):"); plotGamma.setSelected(true); updateFreq(0,true);
      plotGamma.setHorizontalTextPosition(JMenuItem.LEFT);
      
      timePlot[0] = new JCheckBox(" X(t):"); timePlot[0].setSelected(true); updateTime(1,true); timePlot[0].setHorizontalTextPosition(JMenuItem.LEFT); 
      timePlot[1] = new JCheckBox(" W_1(t):"); timePlot[1].setSelected(false); timePlot[1].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[2] = new JCheckBox(" W_2(t):"); timePlot[2].setSelected(false); timePlot[2].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[3] = new JCheckBox(" W_3(t):"); timePlot[3].setSelected(false); timePlot[3].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[4] = new JCheckBox(" W_4(t):"); timePlot[4].setSelected(false); timePlot[4].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[5] = new JCheckBox(" W_5(t):"); timePlot[5].setSelected(false); timePlot[5].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[6] = new JCheckBox(" W_6(t):"); timePlot[6].setSelected(false); timePlot[6].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[7] = new JCheckBox(" W_7(t):"); timePlot[7].setSelected(false); timePlot[7].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[8] = new JCheckBox(" W_8(t):"); timePlot[8].setSelected(false); timePlot[8].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[9] = new JCheckBox(" W_9(t):"); timePlot[9].setSelected(false); timePlot[9].setHorizontalTextPosition(JMenuItem.LEFT);     
      timePlot[10] = new JCheckBox(" W_10(t):"); timePlot[10].setSelected(false); timePlot[10].setHorizontalTextPosition(JMenuItem.LEFT);
      timePlot[11] = new JCheckBox(" Y(t):"); timePlot[11].setSelected(false); timePlot[11].setHorizontalTextPosition(JMenuItem.LEFT);

      freqPlot[0] = new JCheckBox(" \u0393_1(\u03C9):"); freqPlot[0].setSelected(false); freqPlot[0].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[1] = new JCheckBox(" \u0393_2(\u03C9):"); freqPlot[1].setSelected(false); freqPlot[1].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[2] = new JCheckBox(" \u0393_3(\u03C9):"); freqPlot[2].setSelected(false); freqPlot[2].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[3] = new JCheckBox(" \u0393_4(\u03C9):"); freqPlot[3].setSelected(false); freqPlot[3].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[4] = new JCheckBox(" \u0393_5(\u03C9):"); freqPlot[4].setSelected(false); freqPlot[4].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[5] = new JCheckBox(" \u0393_6(\u03C9):"); freqPlot[5].setSelected(false); freqPlot[5].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[6] = new JCheckBox(" \u0393_7(\u03C9):"); freqPlot[6].setSelected(false); freqPlot[6].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[7] = new JCheckBox(" \u0393_8(\u03C9):"); freqPlot[7].setSelected(false); freqPlot[7].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[8] = new JCheckBox(" \u0393_9(\u03C9):"); freqPlot[8].setSelected(false); freqPlot[8].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[9] = new JCheckBox(" \u0393_10(\u03C9):"); freqPlot[9].setSelected(false); freqPlot[9].setHorizontalTextPosition(JMenuItem.LEFT);
      freqPlot[10] = new JCheckBox(" \u0393_X(\u03C9):"); freqPlot[10].setSelected(false); freqPlot[10].setHorizontalTextPosition(JMenuItem.LEFT);    

      periodPlot[0] = new JCheckBox(" Y(t):"); periodPlot[0].setSelected(true); periodPlot[0].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[1] = new JCheckBox(" X(t):"); periodPlot[1].setSelected(false); periodPlot[1].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[2] = new JCheckBox(" W_1(t):"); periodPlot[2].setSelected(false); periodPlot[2].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[3] = new JCheckBox(" W_2(t):"); periodPlot[3].setSelected(false); periodPlot[3].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[4] = new JCheckBox(" W_3(t):"); periodPlot[4].setSelected(false); periodPlot[4].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[5] = new JCheckBox(" W_4(t):"); periodPlot[5].setSelected(false); periodPlot[5].setHorizontalTextPosition(JMenuItem.LEFT); 
      periodPlot[6] = new JCheckBox(" W_5(t):"); periodPlot[6].setSelected(false); periodPlot[6].setHorizontalTextPosition(JMenuItem.LEFT);  
      periodPlot[7] = new JCheckBox(" W_6(t):"); periodPlot[7].setSelected(false); periodPlot[7].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[8] = new JCheckBox(" W_7(t):"); periodPlot[8].setSelected(false); periodPlot[8].setHorizontalTextPosition(JMenuItem.LEFT);
      periodPlot[9] = new JCheckBox(" W_8(t):"); periodPlot[9].setSelected(false); periodPlot[9].setHorizontalTextPosition(JMenuItem.LEFT); 
      periodPlot[10] = new JCheckBox(" W_9(t):"); periodPlot[10].setSelected(false); periodPlot[10].setHorizontalTextPosition(JMenuItem.LEFT);  
      periodPlot[11] = new JCheckBox(" W_10(t):"); periodPlot[11].setSelected(false); periodPlot[11].setHorizontalTextPosition(JMenuItem.LEFT);  

      coeffPlot[0] = new JCheckBox(" sym_b:"); coeffPlot[0].setSelected(false); coeffPlot[0].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[1] = new JCheckBox(" b_1:"); coeffPlot[1].setSelected(true); coeffPlot[1].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[2] = new JCheckBox(" b_2:"); coeffPlot[2].setSelected(false); coeffPlot[2].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[3] = new JCheckBox(" b_3:"); coeffPlot[3].setSelected(false); coeffPlot[3].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[4] = new JCheckBox(" b_4:"); coeffPlot[4].setSelected(false); coeffPlot[4].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[5] = new JCheckBox(" b_5:"); coeffPlot[5].setSelected(false); coeffPlot[5].setHorizontalTextPosition(JMenuItem.LEFT); 
      coeffPlot[6] = new JCheckBox(" b_6:"); coeffPlot[6].setSelected(false); coeffPlot[6].setHorizontalTextPosition(JMenuItem.LEFT); 
      coeffPlot[7] = new JCheckBox(" b_7:"); coeffPlot[7].setSelected(false); coeffPlot[7].setHorizontalTextPosition(JMenuItem.LEFT);
      coeffPlot[8] = new JCheckBox(" b_8:"); coeffPlot[8].setSelected(false); coeffPlot[8].setHorizontalTextPosition(JMenuItem.LEFT); 
      coeffPlot[9] = new JCheckBox(" b_9:"); coeffPlot[9].setSelected(false); coeffPlot[9].setHorizontalTextPosition(JMenuItem.LEFT); 
      coeffPlot[10] = new JCheckBox(" b_10:"); coeffPlot[10].setSelected(false); coeffPlot[10].setHorizontalTextPosition(JMenuItem.LEFT); 
     
      for(i=0;i<11;i++) 
      {
        expVariable[i] = new JCheckBox(" W_" + (i+1) + "(t):"); expVariable[i].setSelected(false); expVariable[i].setEnabled(false);
        expVariable[i].setHorizontalTextPosition(JMenuItem.LEFT);
        expVariable[i].addItemListener(new MyItemListener2()); 
      }
     
      //----------------- Make the sliders -------------------------
      cointText = new JTextField(6); cointText.setText("1.00"); 
      cointText.setToolTipText("If i1 set, the value of the sum of the all the filter coefficients. When 1.0, new filter will be computed"); 
      //brm.setMinimum(0); brm.setMaximum(100); brm.setValue(0); brm.setExtent(1);
      
      for(i=0;i<12;i++)
      {
        cointSliders[i] = new JSlider(JSlider.VERTICAL,0,100,0); 
        cointSliders[i].setPreferredSize(new Dimension(20,50));      
        cointSliders[i].addChangeListener(new MyChangeListener());  
        cointSliders[i].setToolTipText("The value of the constraint at frequency 0 for filter " + Integer.toString(i+1));    
      }
      //brm.addChangeListener(new MyChangeListener());
      for(i=0;i<12;i++)
      {cointSliders[i].setEnabled(false);}

      plot_zpc = new JCheckBox[4];
      for(i=0;i<4;i++) 
      {       
        plot_zpc[i] = new JCheckBox();
        plot_zpc[i].setSelected(false); 
        plot_zpc[i].setHorizontalTextPosition(JMenuItem.LEFT);
        plot_zpc[i].addItemListener(new MyPlotListener()); 
      }
      plot_zpc[0].setSelected(true);  plot_zpc[1].setSelected(true);
      plot_zpc[0].setText("\u0393(\u03C9)"); plot_zpc[1].setText("A(\u03C9)"); 
      plot_zpc[2].setText("\u03A6(\u03C9)"); plot_zpc[3].setText("\u03A6(\u03C9)/\u03C9"); 
 




      rkhsCombo = new JComboBox<String>();
      rkhsCombo.addItem("Beta Kernel");
      rkhsCombo.addItem("Henderson Kernel"); 
      pBox = new JComboBox<String>();
      qBox = new JComboBox<String>();
      for(i=0;i<10;i++)
      {pBox.addItem(""+i); qBox.addItem(""+i);}
      
      rkhsCombo.addActionListener(new MyComboActionListener());
      pBox.addActionListener(new MyComboActionListener());
      qBox.addActionListener(new MyComboActionListener());

      rkhsBox = new JCheckBox("RKHS:");
      armaSD = new JCheckBox("ARMA:");
      rkhsBox.setSelected(false); 
      armaSD.setSelected(false);
      rkhsBox.setHorizontalTextPosition(JMenuItem.LEFT);
      armaSD.setHorizontalTextPosition(JMenuItem.LEFT);

      rkhsBox.addItemListener(new MyItemListener()); 
      armaSD.addItemListener(new MyItemListener());

      rkhsBox.setToolTipText("Use symmetric reproducing kernel Hilbert Space target filter.");   
      armaSD.setToolTipText("Use optimized ARMA spectral density for data.");      

      rkhsCombo.setPreferredSize(new Dimension(70, 15));
      pBox.setPreferredSize(new Dimension(40, 15));
      qBox.setPreferredSize(new Dimension(40, 15));


      i1Check = new JCheckBox("i1:"); i1Check.setHorizontalTextPosition(JMenuItem.LEFT); i1Check.setSelected(false);
      i2Check = new JCheckBox("i2:"); i2Check.setHorizontalTextPosition(JMenuItem.LEFT); i2Check.setSelected(false);
      i1Check.setToolTipText("Apply constraints on the filter coefficients that they must sum to \u0393(0). Use when data not zero mean.");   
      i2Check.setToolTipText("Apply constraints on the filter coefficients that first moment of the filter coefficients sum to 0. Use when detecting turning points.");      

      dCheck = new JCheckBox("d:");  dCheck.setSelected(false); dCheck.setHorizontalTextPosition(JMenuItem.LEFT);   
      DCheck = new JCheckBox("D:");  DCheck.setSelected(false); DCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      dCheck.setToolTipText("Accounts for first-order differencing of the data. Apply when data is nonstationary."); 
      DCheck.setToolTipText("Accounts for first-order seasonal differencing of the data. Apply when data is seasonal.");          


      reComp = new JCheckBox("   Recompute Filter:"); reComp.setSelected(true); reComp.setHorizontalTextPosition(JMenuItem.LEFT);  
      iterCheck = new JCheckBox("   Iterative Method:"); iterCheck.setSelected(false); iterCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      reComp.setToolTipText("When off, only additional points can be added to time series to see effects of filter on out-of-sample data."); 
      iterCheck.setToolTipText("Activates an iterative technique to direct filter computation that typically yields better separation of amplitude and phase error effects."); 
    
      amplitudeBut.addItemListener(new MyItemListener());
      gammaBut.addItemListener(new MyItemListener());
      timeDelayBut.addItemListener(new MyItemListener());
   
      dCheck.addItemListener(new MyItemListener()); 
      DCheck.addItemListener(new MyItemListener());

      plotXf.addItemListener(new MyItemListener()); 
      plotGamma.addItemListener(new MyItemListener());

      for(i=0; i <= 10; i++)
      {freqPlot[i].addItemListener(new MyItemListener()); coeffPlot[i].addItemListener(new MyItemListener());}
  
      for(i=0; i < 11; i++)
      {periodPlot[i].addItemListener(new MyItemListener()); timePlot[i].addItemListener(new MyItemListener());}
      timePlot[11].addItemListener(new MyItemListener());

      accountPlot = new JCheckBox("Account Returns:");  accountPlot.setHorizontalTextPosition(JMenuItem.LEFT); accountPlot.setSelected(true);
      logpricePlot = new JCheckBox("Log-Price Series:"); logpricePlot.setHorizontalTextPosition(JMenuItem.LEFT); logpricePlot.setSelected(false);
      logreturnPlot = new JCheckBox("Target Series:"); logreturnPlot.setHorizontalTextPosition(JMenuItem.LEFT); logreturnPlot.setSelected(false);
      signalPlot = new JCheckBox("Trading Signal:"); signalPlot.setHorizontalTextPosition(JMenuItem.LEFT); signalPlot.setSelected(false);
      linesPlot = new JCheckBox("Derivative Sig:"); linesPlot.setHorizontalTextPosition(JMenuItem.LEFT); linesPlot.setSelected(false);
      buysellPlot = new JCheckBox("Buy Indicators:"); buysellPlot.setHorizontalTextPosition(JMenuItem.LEFT); buysellPlot.setSelected(false);
      filteredPrice = new JCheckBox("Price filter:"); filteredPrice.setHorizontalTextPosition(JMenuItem.LEFT); filteredPrice.setSelected(false); 
      filteredPrice.setEnabled(false);
      
      linesPlot.setEnabled(false);
      accountPlot.addItemListener(new MyItemListener());
      logpricePlot.addItemListener(new MyItemListener());
      logreturnPlot.addItemListener(new MyItemListener());
      signalPlot.addItemListener(new MyItemListener());
      linesPlot.addItemListener(new MyItemListener());
      buysellPlot.addItemListener(new MyItemListener());
      filteredPrice.addItemListener(new MyItemListener());
      
      
      
      i1Check.addItemListener(new MyItemListener()); 
      i2Check.addItemListener(new MyItemListener());
      reComp.addItemListener(new MyItemListener());
      iterCheck.addItemListener(new MyItemListener());
      
      sarimaCheck.addItemListener(new MyItemListener());
      simCheck.addItemListener(new MyItemListener());
      noneCheck.addItemListener(new MyItemListener());  


      //----- Add stuff for ARTpanel-------------------------------------
      freqIntsSlider = new JSlider();
      maxFreqSlider = new JSlider();

      freqIntsSlider.setFont(new Font("Arial", 0, 10));  freqIntsSlider.setMajorTickSpacing(1);
      freqIntsSlider.setMaximum(20); freqIntsSlider.setMinimum(10); freqIntsSlider.setOrientation(JSlider.VERTICAL);
      freqIntsSlider.setValue(13); freqIntsSlider.addChangeListener(new MyOtherChangeListener());

      maxFreqSlider.setFont(new Font("Arial", 0, 10)); maxFreqSlider.setMajorTickSpacing(1);
      maxFreqSlider.setMinimum(10); maxFreqSlider.setMaximum(100); maxFreqSlider.setOrientation(JSlider.VERTICAL); maxFreqSlider.setSnapToTicks(true);
      maxFreqSlider.setToolTipText("Change maximum frequency for extracting trend/cycle"); maxFreqSlider.setValue(40);
      maxFreqSlider.addChangeListener(new MyOtherChangeListener());

      turningSlider = new JSlider();
      turningSlider.setFont(new Font("Arial", 0, 10));  turningSlider.setMajorTickSpacing(1);
      turningSlider.setMaximum(101); turningSlider.setMinimum(0); turningSlider.setOrientation(JSlider.HORIZONTAL);
      turningSlider.setValue(0); turningSlider.addChangeListener(new MyOtherChangeListener());       


      n_div = 13; freq_start = .40; mdfa.setNDiv(freq_start, n_div); computeFreqDivision();

      computeTimeFreq = new JButton();
      contARTCheck = new JCheckBox();
      compARTFilter = new JButton();

      computeTimeFreq.setText("Compute"); computeTimeFreq.setFont(new Font("Arial", 0, 12));
      computeTimeFreq.setToolTipText("Compute the time-freq decomposition and plot the maps");
      computeTimeFreq.addActionListener(new MyTimeActionListener());

      contARTCheck.setFont(new Font("Arial", 0, 12)); contARTCheck.setSelected(false);
      contARTCheck.setText("Continuous"); contARTCheck.addItemListener(new MyItemListener());
      contARTCheck.setToolTipText("Continuously compute the filter using trilemma settings upon mouse mouvments in ART panel");
      compARTFilter.setText("Compute Filter");
      compARTFilter.addActionListener(new MyTimeActionListener());

      //----- Initialize the rest of TimeFreqPlotPanel and such----------
      initTMPanelComponents(); initARTPanelComponents();
      tfplotPanel = new TimeFreqPlotPanel();

      //---- Initialize zpc_gene checkboxes 
      filterPlot = new JRadioButton[3]; 
      filterPlotGroup = new ButtonGroup();

      for(i=0;i<3;i++) 
      {
        filterPlot[i] = new JRadioButton();
        filterPlot[i].setFont(new Font("Arial", 0, 12));
        filterPlot[i].setEnabled(false); 
        filterPlotGroup.add(filterPlot[i]);
        filterPlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
      }
      filterPlot[0].setText("I-MDFA"); filterPlot[0].setToolTipText("Filter using only I-MDFA filter (default)");
      filterPlot[1].setText("ZPC");    filterPlot[1].setToolTipText("Filter using only ZPC filter");
      filterPlot[2].setText("Hybrid"); filterPlot[2].setToolTipText("Inject the ZPC Gene into the I-MDFA filter. Only available after optimizing zpc filter.");
      filterPlot[0].setSelected(true);
      filterPlot[0].addItemListener(new PlotItemListener()); 
      filterPlot[1].addItemListener(new PlotItemListener());
      filterPlot[2].addItemListener(new PlotItemListener());
      plot_number = 0; filterPlot[0].setEnabled(true);
   }


      class PlotItemListener implements ItemListener
      {
          public void itemStateChanged(ItemEvent e)
          {
            JRadioButton source = (JRadioButton)e.getSource();
            for(int i=0;i<3;i++)
            { 
              if((source == filterPlot[i]) && filterPlot[i].isSelected()) 
              {
               if((i==1) || (i==2))
               {
                 //computeZPCFilter();
                 injectZPCGeneSimple();
                 plot_number = i;             
                 mdfa_canvas.plotHybrid(i); 
                 filter_canvas.plotHybrid(i); 
                 period_canvas.plotHybrid(i);
                 updateTrading(plot_number);
               }
               else
               {
                 plot_number = 0;      
                 mdfa_canvas.plotHybrid(0); 
                 filter_canvas.plotHybrid(0); 
                 period_canvas.plotHybrid(0);
                 updateTrading(plot_number);
               }
              }
            }
            //if(b_outof_sample) {filterScore();}
            
          }
      }

     class MyTimeActionListener implements ActionListener
     {
      
      public void actionPerformed(ActionEvent e)
      {
         double[] xtemp = new double[mdfa.flength];
         if(computeTimeFreq == e.getSource())
         {
             
             //mdfa.time_freq_decomp(); n_imfs = mdfa.n_imfs;
             for(int i = 0; i<mdfa.flength; i++) xtemp[i] = mdfa.tseries[L-1 + i];

             tfplotPanel.setSlider(0);
             tfplotPanel.setNDiv(freq_ints);

             if(plotTF == 0) {tfplotPanel.setTimeFreqData(mdfa.amMap, mdfa.flength, n_imfs, 0);}
             else if(plotTF == 1) {tfplotPanel.setTimeFreqData(mdfa.fmMap, mdfa.flength, n_imfs, 1);}
             else if(plotTF == 2) {tfplotPanel.setTimeFreqData(mdfa.phaseMap, mdfa.flength, n_imfs, 2);}
             else if(plotTF == 3) {tfplotPanel.setTimeFreqData(mdfa.ifreqMap, mdfa.flength, n_imfs, 3);}        
     
             tcimfPanel.setIMFData(mdfa.amMap, mdfa.fmMap, n_imfs, mdfa.flength);
             tcimfPanel.setTrend(mdfa.trend_cycle);
             tcimfPanel.setIMF(0,K+1,freq_ints);
             tcimfPanel.setData(xtemp);
             //computeTimeFreq.setEnabled(false);
         }
         else if(compARTFilter == e.getSource())
         {
             computeFilter = true;
             mdfa.computeFilterGeneral(computeFilter,false);
             updatePlots(false,true);
             compute.setEnabled(false); 
             computeTimeFreq.setEnabled(true);
         }
       }
    }  





    //------------------------------------------------------------------------
    class MyOtherChangeListener implements ChangeListener 
    {
       public void stateChanged(ChangeEvent e) 
       {   
         JSlider s = (JSlider)e.getSource(); int v;
         if(s == freqIntsSlider) 
         { 
            n_div = s.getValue(); 
            mdfa.setNDiv(freq_start, n_div);
            computeFreqDivision();
            tfplotPanel.setNDiv(freq_ints);
            computeTimeFreq.setEnabled(true);
         }  
         else if(s == maxFreqSlider) 
         {
            freq_start = s.getValue()*.01;
            mdfa.setNDiv(freq_start, n_div);
            computeFreqDivision();
            tfplotPanel.setNDiv(freq_ints);
            computeTimeFreq.setEnabled(true);
         } 
         else if(s == turningSlider)
         {
           if(!s.getValueIsAdjusting()) 
           {
             v = s.getValue();
             lambda_3 = v*.01; 
             turningText.setText(df.format(lambda_3)); 
             setTurningLambda3(); 
           }          
         }       
       }
    } //------------------------------------------------------------------------

    class MyChangeListener implements ChangeListener 
    {
       int n_max; int v;
       public void stateChanged(ChangeEvent e) 
       {   
         JSlider s = (JSlider)e.getSource();
         if(!s.getValueIsAdjusting()) 
         {
            //System.out.println("adjusting");
            n_max = getMaxAvailable();
            v = s.getValue();
            if(n_max < 0) 
            {s.setValue(v + n_max);} 
            else if(n_max==0)
            {setWConstVals(true);}
            
         }         
       }
    }

    
    class MyChangeListener2 implements ChangeListener 
    {
       int n_max; int v; double nval; double omega0val;
       public void stateChanged(ChangeEvent e) 
       {   
         JSlider s = (JSlider)e.getSource();
         if(s == tradingfreqSlider)       
         {
          if(!s.getValueIsAdjusting()) 
          {
            omega0val = Math.max(w0,(314.0/8.0));
            v = s.getValue(); 
            nval = omega0val*(100.0 - (double)v)/100.0 + (314.0/2.0)*((double)v)/100.0; 
            omega1Bar.setValue((int)nval);
            if(!autoComp) {computeFilterNow();}
          }
         }  
         if(s == tradingcostSlider)       
         {
          if(!s.getValueIsAdjusting()) 
          {            
            v = s.getValue();
            trading_cost = .0001*v;    
            updatePlots(false,false);        
            tradingcostText.setText(""+trading_cost); 
          }
         }
         if(s == riskfreeSlider)       
         {
          if(!s.getValueIsAdjusting()) 
          {
            v = s.getValue();
            risk_free = .0001*v;
            riskfreeText.setText(""+df.format(risk_free));  
            updatePlots(false,false);          
          }
         }
       }
    }



    class MyPlotListener implements ItemListener 
    {
        public void itemStateChanged(ItemEvent e)
        {
         int i; boolean sel; 
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}

         for(i=0;i<4;i++)
         {
          if(e.getSource() == plot_zpc[i])
          {zpcFreqCanvas.setPlotDraw(sel,i);}
         }
        }
    }

    class MyComboActionListener implements ActionListener
    {
       int pc,qc;
       @SuppressWarnings("rawtypes")
	   public void actionPerformed(ActionEvent e)
       {

        if(e.getSource() == pBox || e.getSource() == qBox)
        {                         
         pc = ((JComboBox)e.getSource()).getSelectedIndex();
         qc = ((JComboBox)e.getSource()).getSelectedIndex();
         if(armaSD.isSelected())
         {
           computeARMASpectralDensity(pc, qc);
           mdfa.computeFilterGeneral(computeFilter,false);
           updatePlots(false,true);
           compute.setEnabled(false); 
         }
        }
        else if(e.getSource() == rkhsCombo)
        {
         if(rkhsBox.isSelected())
         {
          pc = ((JComboBox)e.getSource()).getSelectedIndex();
         
          if(pc == 0)
          {tfilter.setReproducingKernelFilter(3);}
          else if(pc == 1)
          {tfilter.setReproducingKernelFilter(9);}

          w0 = 0.0; w1 = tfilter.w1;
          mdfa.set_cutoff(tfilter.w1); 
          mdfa.set_cutoff0(0.0);
          set_Gamma(tfilter.Gamma); setBC(tfilter.bsym);      
          period_canvas.setCutoffs(0, tfilter.w1);
         }
        }
      }
   }



    public void sliderChangeMax(JSlider s, int n_max)
    {
      for(int i=0;i<n_rep-1;i++)
      {
         if(s != cointSliders[i]) {cointSliders[i].setMaximum(min(100,cointSliders[i].getValue() + n_max));}
      }
    }

    public int getMaxAvailable()
    {
      int sum = 0; int n; w_sum = 0.0;
      for(int i=0;i<n_rep-1;i++) 
      {
       n = cointSliders[i].getValue();
       sum = sum + n; w_sum = w_sum + .01*n;
      }
      cointText.setText(" "+ df.format(w_sum));
      return 100-sum;
    } 

    public void setWConstVals(boolean comp)
    {      
      int i,n; w_sum = 0.0;
      for(i=0;i<n_rep-1;i++) 
      {
         n = cointSliders[i].getValue(); w_const[i] = .01*n;
      }
      mdfa.setWConstVals(w_const); 
      if(comp && (i1==1))
      {
        mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);
      } 
    } 


    public void activateWConst(int a, boolean comp)
    {    
       int i; boolean a1=false; if(a==1){a1 = true;}
       for(i=0;i<n_rep-1;i++)
       {cointSliders[i].setEnabled(a1); cointSliders[i].setValue(0);}
       for(i=n_rep-1;i<12;i++)
       {cointSliders[i].setEnabled(false);}     
       cointSliders[0].setValue(100);
       setWConstVals(false); 
    }


    //-------------------Testing inputs from keyboarded values----------
    public void test_nobs(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i >=120 && i < 2000) {nObsBar.setValue(i); n_obs = i; nText.setText(""+n_obs);}
                                          else nText.setText(""+n_obs);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); nText.setText(""+n_obs);}
    }

    public void test_lambda(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 400.0) {lambdaBar.setValue((int)(e*10)); lambda = e; 
                                                                    lambdaText.setText(""+lambda);}
                                            else lambdaText.setText(df.format(lambda));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); lambdaText.setText(""+lambda);}
    }

    public void test_delta(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 1.0) {deltaBar.setValue((int)(e*1000)); tradeDelta = e; 
                                                                    deltaText.setText(""+tradeDelta);}
                                            else deltaText.setText(""+tradeDelta);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); deltaText.setText(""+tradeDelta);}
    }    
    
    public void test_smooth(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 1.0) {smoothBar.setValue((int)(e*1000)); smooth = e; 
                                                                    smoothText.setText(""+smooth);}
                                            else smoothText.setText(""+smooth);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); smoothText.setText(""+smooth);}
    }
    public void test_decay(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 1.0) {decayBar.setValue((int)(e*1000)); decay = e; 
                                                                    decayText.setText(""+decay);}
                                            else decayText.setText(""+decay);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); decayText.setText(""+decay);}
    }
    public void test_decay2(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 1.0) {decay2Bar.setValue((int)(e*1000)); decay2 = e; 
                                                                    decay2Text.setText(""+decay2);}
                                            else decay2Text.setText(""+decay2);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); decay2Text.setText(""+decay2);}
    }    
    public void test_cross(String s)
    {   
      double e; 
      try{e = Double.parseDouble(s.trim()); if(e >=0 && e < 1.0) {crossBar.setValue((int)(e*1000)); cross = e; 
                                                                    crossText.setText(""+cross);}
                                            else smoothText.setText(""+cross);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); crossText.setText(""+cross);}
    }



    public void test_L(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i >=10 && i < n_obs) {LBar.setValue(i); L = i; LText.setText(""+L);}
                                          else LText.setText(""+L);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); LText.setText(""+L);}
    }

    public void test_l1(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i >=100 && i <= 500) {l1Bar.setValue(i); L1 = i; l1Text.setText(""+L1);}
                                          else l1Text.setText(""+L1);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); l1Text.setText(""+L1);}
    }

    public void test_l2(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i >=100 && i <= 300) {l2Bar.setValue(i); L2 = i; l2Text.setText(""+L2);}
                                          else l2Text.setText(""+L2);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); l2Text.setText(""+L2);}
    }
    
    
    
    public void test_exp(String s)
    {   
      double e; int expint;
      try{e = Double.parseDouble(s.trim()); if(e >=0.0 && e <= 40.0) {expint=(int)(10*e); expBar.setValue(expint); 
                                                                  expweight = e; expText.setText(""+df.format(expweight));}
                                          else expText.setText(""+df.format(expweight));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); expText.setText(""+df.format(expweight));}
    }


    public void test_omega0(String s)
    {   
      double e; int expint;  omega0Bar.getValue();
      try{e = Double.parseDouble(s.trim()); if(e >=0.0 && e <= w1) {expint=(int)(100*e); omega0Bar.setValue(expint); 
                                                                  w0 = e; omega0Text.setText(""+df.format(w0));}
                                          else omega0Text.setText(""+df.format(w0));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); omega0Text.setText(""+df.format(w0));}
    }
    public void test_omega1(String s)
    {   
      double e; int expint;
      try{e = Double.parseDouble(s.trim()); if(e >=w0 && e <= Math.PI) {expint=(int)(100*e); omega1Bar.setValue(expint); 
                                                                  w1 = e; omega1Text.setText(""+df.format(w1));}
                                          else omega1Text.setText(""+df.format(w1));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); omega1Text.setText(""+df.format(w1));}
    }
    public void test_omega2(String s)
    {   
      double e; int expint;  omega2Bar.getValue();
      try{e = Double.parseDouble(s.trim()); if(e >=w1 && e <= w3) {expint=(int)(100*e); omega2Bar.setValue(expint); 
                                                                  w2 = e; omega2Text.setText(""+df.format(w0));}
                                          else omega2Text.setText(""+df.format(w2));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); omega2Text.setText(""+df.format(w2));}
    }
    public void test_omega3(String s)
    {   
      double e; int expint;
      try{e = Double.parseDouble(s.trim()); if(e >=w2 && e <= Math.PI) {expint=(int)(100*e); omega3Bar.setValue(expint); 
                                                                  w3 = e; omega3Text.setText(""+df.format(w3));}
                                          else omega3Text.setText(""+df.format(w3));}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); omega3Text.setText(""+df.format(w3));}
    }    
    public void test_num1(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if((1.0*i)/(1.0*den1) < (1.0*num2)/(1.0*den2)) {num1Bar.setValue(i); num1 = i; num1Text.setText(""+i);}
                                          else num1Text.setText(""+num1);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); num1Text.setText(""+num1);}
    }
    public void test_num2(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if((1.0*num1)/(1.0*den1) < (1.0*i)/(1.0*den2)) {num2Bar.setValue(i); num2 = i; num2Text.setText(""+i);}
                                          else num2Text.setText(""+num2);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); num2Text.setText(""+num2);}
    }
    public void test_den1(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i>0 && ((1.0*num1)/(1.0*i) < (1.0*num2)/(1.0*den2))) {den1Bar.setValue(i); den1 = i; den1Text.setText(""+i);}
                                          else den1Text.setText(""+den1);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); den1Text.setText(""+den1);}
    }
    public void test_den2(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if((1.0*num1)/(1.0*den1) < (1.0*num2)/(1.0*i)) {den2Bar.setValue(i); den2 = i; den2Text.setText(""+i);}
                                          else den2Text.setText(""+den2);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); den2Text.setText(""+den2);}
    }
    public void test_nrep(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i > 0 && i <= 5) {nrepBar.setValue(i); n_rep = i; nrepText.setText(""+i);}
                                          else nrepText.setText(""+n_rep);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); nrepText.setText(""+n_rep);}
    }

    //-----------------------------------------------------------------

    public void setupFilterAdjusters()
    {

       nObsBar = new JScrollBar(JScrollBar.HORIZONTAL,300,2,120,2000);
       nObsBar.setUnitIncrement(1);      
       nText = new JTextField(3); nText.setColumns(3);
       nText.setText(""+n_obs);     
       nText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_nobs(nText.getText());}} );  
       nObsBar.setEnabled(simulate);

       sampBar = new JScrollBar(JScrollBar.HORIZONTAL,K,2,K,2000);
       sampBar.setUnitIncrement(2);       
       sampBar.setMinimum(K);
       sampText = new JTextField(3); sampText.setColumns(3);
       sampText.setText(""+nsamp);      
       sampBar.setEnabled(simulate);

       nrepBar = new JScrollBar(JScrollBar.HORIZONTAL,5,1,1,6);
       nrepBar.setUnitIncrement(1);      
       nrepText = new JTextField(1); nrepText.setColumns(3);
       nrepText.setText(""+n_rep);
       nrepText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_nrep(nrepText.getText());}} );  
       nrepBar.setEnabled(simulate);

       LBar = new JScrollBar(JScrollBar.HORIZONTAL,L,2,1,100);
       LBar.setUnitIncrement(1);  
       LBar.setMinimum(1);
       LBar.setMaximum(n_obs-1);
       LText = new JTextField(2); LText.setColumns(3);
       LText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_L(LText.getText());}} );  
       LText.setText(""+L);
       
       l1Bar = new JScrollBar(JScrollBar.HORIZONTAL,L1,2,0,500);
       l1Bar.setUnitIncrement(1);      
       l1Text = new JTextField(2); l1Text.setColumns(3);
       l1Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_l1(l1Text.getText());}} );  
       l1Text.setText(""+L1);
       
       l2Bar = new JScrollBar(JScrollBar.HORIZONTAL,L1,2,0,300);
       l2Bar.setUnitIncrement(1);      
       l2Text = new JTextField(2); l2Text.setColumns(3);
       l2Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_l2(l2Text.getText());}} );  
       l2Text.setText(""+L2);       
                       
       lambdaBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,4000);
       lambdaBar.setValue(0);
       lambdaBar.setUnitIncrement(5);
       lambdaText = new JTextField(3); lambdaText.setColumns(3);
       lambdaText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_lambda(lambdaText.getText());}} );  
       lambdaText.setText("" + lambda);

       shiftBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,400);
       shiftBar.setValue(200);
       shiftBar.setUnitIncrement(1);
       shiftText = new JTextField(3); shiftText.setColumns(3);
       shiftText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_lambda(shiftText.getText());}} );  
       shiftText.setText("" + 0.0);

       smoothBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       smoothBar.setValue(0);
       smoothBar.setUnitIncrement(1);
       smoothText = new JTextField(3); smoothText.setColumns(3);
       smoothText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_smooth(smoothText.getText());}} );  
       smoothText.setText("" + smooth);
 
       decayBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       decayBar.setValue(0);
       decayBar.setUnitIncrement(1);
       decayText = new JTextField(3); decayText.setColumns(3);
       decayText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_decay(decayText.getText());}} );  
       decayText.setText("" + decay);

       decay2Bar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       decay2Bar.setValue(0);
       decay2Bar.setUnitIncrement(1);
       decay2Text = new JTextField(3); decay2Text.setColumns(3);
       decay2Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_decay2(decay2Text.getText());}} );  
       decay2Text.setText("" + decay2);


       crossBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       crossBar.setValue(0);
       crossBar.setUnitIncrement(1);
       crossText = new JTextField(3); crossText.setColumns(3);
       crossText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_cross(crossText.getText());}} );  
       crossText.setText("" + cross);

       deltaBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       deltaBar.setValue(0);
       deltaBar.setUnitIncrement(1);
       deltaText = new JTextField(3); crossText.setColumns(3);
       deltaText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_delta(deltaText.getText());}} );  
       deltaText.setText("" + tradeDelta);       
           
       expBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,400);
       expBar.setValue(0);
       expBar.setUnitIncrement(5);
       expText = new JTextField(3); expText.setColumns(3);
       expText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_exp(expText.getText());}} );  
       expText.setText("" + expweight);
                  
       lagBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,72);
       lagBar.setValue(36);
       lagBar.setUnitIncrement(1);
       textLag = new JTextField(2);
       Lag = lagBar.getValue()-36;
       textLag.setText("" + Lag);

       onestepBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,1300);
       onestepBar.setValue(0);
       onestepBar.setUnitIncrement(1);
       onestepText = new JTextField(3);
       onestepText.setText("" + 0);       
       
       lambda3Bar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,700);
       lambda3Bar.setValue(0);
       lambda3Bar.setUnitIncrement(1);
       lambda3Text = new JTextField(3);
       lambda3Text.setText("" + 0);             
       
       
       outof_sample = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,100); 
       outof_sample.setEnabled(true); outof_sample.setValue(0); outof_sample.setUnitIncrement(1);
       outof_sample.setMaximum(0);
       outof_sampleText = new JTextField(3); outof_sampleText.setText("" + 0);
       b_outof_sample = true;     

      
       compute = new JButton("Apply");
       compute.setActionCommand ("compute");
       compute.addActionListener(new MyActionListener());
       compute.setEnabled(false);
       compute.setToolTipText("Recompute the filter with given settings of target filter."+eof+"Only applicable when new target filter has been set.");

       mdfaText = new JTextField(4);
       mdfaText.setText("" + df.format(mdfa.MDFAmin)); 
       mdfaText.setToolTipText("The pseudo-information criterion and effective degrees of freedom for the computed filter"); 
       //----- Make adjustement listener -----------

       AdjustmentListener al = new AdjustmentListener()  {
        public void adjustmentValueChanged(AdjustmentEvent e) {

            int lambdaint; int expint; int smoothint; int decayint, decay2int, crossint, deltaint;
           // computeFilter = true;  //System.out.println(e.getAdjustable());           

             if(e.getAdjustable() == sampBar)
             {
                 //System.out.println("samp change");
                 if(cPeriod){
                 nsamp = sampBar.getValue();
                 sampText.setText(""+nsamp);
                 mdfa.set_Samps(nsamp); zpc.set_Samps(nsamp);
                 mdfa.computeSampleIns();
                 period_canvas.setPeriodogramXf(mdfa.period_xf, nsamp); 
                 period_canvas.go();}
             }
             if(e.getAdjustable() == outof_sample)
             {
                addOutofSample();
             }
             else if(e.getAdjustable() == nrepBar)
             {      
                      
               if(simulate)
               {
                
                n_rep = nrepBar.getValue();
                nrepText.setText(""+ n_rep);
                mdfa.set_nreps(n_rep); 
                mdfa_canvas.setNRep(n_rep);  
                filter_canvas.setNRep(n_rep);  
                activateWConst(i1,false);
                simulateIMDFA(false);                         
               }
               else
               { 
                 //nrepBar.setValue(n_rep);
                 nrepText.setText(""+ n_rep);
                 mdfa.set_nreps(n_rep); 
                 mdfa_canvas.setNRep(n_rep);  
                 filter_canvas.setNRep(n_rep);
                 activateWConst(i1,false);
               }
             }
            else if(e.getAdjustable() == lambdaBar)
            { 
                
                lambdaint = lambdaBar.getValue();
                lambda = lambdaint*.1;
                lambdaText.setText(df.format(lambda)); 
                setLambda(lambda);  
                
             }
            else if(e.getAdjustable() == shiftBar)
            { 
             
                shiftint = shiftBar.getValue();
                shift = (shiftint-200)*.1;
                shiftText.setText(df.format(shift)); 
                setShift(shift);        
            } 
            else if(e.getAdjustable() == smoothBar)
            { 
             
                smoothint = smoothBar.getValue();
                smooth = smoothint*.001;
                smoothText.setText(df2.format(smooth)); 
                setSmooth(smooth);  
                
             }
            else if(e.getAdjustable() == deltaBar)
            {
               deltaint = deltaBar.getValue();
               tradeDelta = deltaint*.00001; 
               deltaText.setText(df4.format(tradeDelta));
               changeTradingThreshold(tradeDelta);
               getAccount_canvas().setTradeThreshold(tradeDelta);
               //System.out.println("tradeDelta = " + tradeDelta);           
            }
            else if(e.getAdjustable() == decayBar)
            { 
             
                decayint = decayBar.getValue();
                decay = decayint*.001;
                decayText.setText(df2.format(decay)); 
                setDecay(decay);  
                
             }
            else if(e.getAdjustable() == decay2Bar)
            { 
             
                decay2int = decay2Bar.getValue();
                decay2 = decay2int*.001;
                decay2Text.setText(df2.format(decay2)); 
                setDecay2(decay2);  
                
             }
            else if(e.getAdjustable() == crossBar)
            { 
             
                crossint = crossBar.getValue();
                cross = crossint*.001;
                crossText.setText(df2.format(cross)); 
                setCross(cross);  
                
             }
             else if(e.getAdjustable() == expBar)
             {
                expint = expBar.getValue();
                expweight = .1*expint;
                expText.setText(""+df2.format(expweight)); 
                setExp(expweight);              
             }
             else if(e.getAdjustable() == LBar)
             {
                //System.out.println(computeFilter);
                L = LBar.getValue();               
                LText.setText(""+L); 
                setL(L);              
             }
             else if(e.getAdjustable() == lagBar)
             {
               Lag = lagBar.getValue()-36;
               textLag.setText("" + Lag);
               setLag(Lag);
             }
             else if(e.getAdjustable() == onestepBar)
             {
               onestep = onestepBar.getValue()*.01;
               onestepText.setText("" + df2.format(onestep));
               setOnestep();
             }
             else if(e.getAdjustable() == lambda3Bar)
             {
               onestep_diff = lambda3Bar.getValue()*.01;
               lambda3Text.setText("" + df2.format(onestep_diff));
               setDiffOnestep();
             }             
             
             else if(e.getAdjustable() == l1Bar)
             { 
               if((L1 == 0) && (l1Bar.getValue() > 0))
               {
                   L1 = 36; l1Bar.setValue(L1);             
               }
               else if((L1 > 35) && (l1Bar.getValue() < 36))
               {
                   L1 = 0; l1Bar.setValue(L1);
               }
               else
               {
                   L1 = l1Bar.getValue();
                   l1Text.setText(""+L1); 

                   if(!simulate && data_uploaded) {adjustSymmetricDataPoints();}
                   else
                   {
                    tfilter.setL1(L1); tfilter.recomputeGamma(); setFilter();               
                    simulateIMDFA(true);
                   }         
               }
                 
             }
             else if(e.getAdjustable() == l2Bar)
             { 
               L2 = l2Bar.getValue();
               l2Text.setText(""+L2); 

               if(!simulate && data_uploaded) {adjustSymmetricDataPoints();}                 
             }             
           }
         };
        
        
         sampBar.addAdjustmentListener(al);
         //nObsBar.addAdjustmentListener(al);
         lagBar.addAdjustmentListener(al);
         outof_sample.addAdjustmentListener(al); 
         expBar.addAdjustmentListener(al);
         onestepBar.addAdjustmentListener(al);
         deltaBar.addAdjustmentListener(al);
         smoothBar.addAdjustmentListener(al);
         decayBar.addAdjustmentListener(al);
         decay2Bar.addAdjustmentListener(al);       
         crossBar.addAdjustmentListener(al);
         lambdaBar.addAdjustmentListener(al);
         lambda3Bar.addAdjustmentListener(al);
         l1Bar.addAdjustmentListener(al);
         LBar.addAdjustmentListener(al);
         nrepBar.addAdjustmentListener(al);
         shiftBar.addAdjustmentListener(al);
    }


   public void setupFilterDesign()
   {
      
      lowBut = new JRadioButton("Low-Pass:"); lowBut.setHorizontalTextPosition(JMenuItem.LEFT); lowBut.setSelected(true);
      highBut = new JRadioButton("Multi-Pass:"); highBut.setHorizontalTextPosition(JMenuItem.LEFT); 
      rampBut = new JRadioButton("Ramp-Pass:"); rampBut.setHorizontalTextPosition(JMenuItem.LEFT); 
      bandBut = new JRadioButton("Band-Pass:"); bandBut.setHorizontalTextPosition(JMenuItem.LEFT); 

      lowBut.setToolTipText("Low-pass filter design.");
      highBut.setToolTipText("Multi-pass filter design.");
      rampBut.setToolTipText("Low-pass filter design with ramping effect to the stop-band.");
      bandBut.setToolTipText("Band-pass filter design.");


       new JLabel(); 
       JLabel omegaLabel = new JLabel();
       JLabel fracLabel = new JLabel();    
 
      omegaLabel.setToolTipText("Set the frequency ranges for \u03C9_1 and \u03C9_2 directly by entering decimal values less than \u03C0 or uMath.sing the scrolls.");
      fracLabel.setToolTipText("Set the frequency ranges for \u03C9_1 and \u03C9_2 by uMath.sing the fractional format n \u03C0/d. The left two boxes determine \u03C9_1 and the left two boxes determine \u03C9_2");

      automatic = new JCheckBox("Auto:"); automatic.setSelected(false); automatic.setHorizontalTextPosition(JMenuItem.LEFT); 
      automatic.setToolTipText("Sets automatic recomputation of filter when any frequency adjustments are made.");
      automatic.addItemListener(new MyItemListener());

      //periodsBut = new JRadioButton("Periods p_1,p_2:"); periodsBut.setHorizontalTextPosition(JMenuItem.LEFT); 
      omegaBut = new JRadioButton("(\u03C9_0, \u03C9_1):"); omegaBut.setHorizontalTextPosition(JMenuItem.LEFT); omegaBut.setSelected(true);
      fracBut = new JRadioButton("(\u03C9_2, \u03C9_3):"); fracBut.setHorizontalTextPosition(JMenuItem.LEFT); 
      omegaBut.setToolTipText("Set the frequency ranges for \u03C9_0 and \u03C9_1 directly by entering decimal values less than \u03C0 or uMath.sing the scrolls.");
      fracBut.setToolTipText("Set the frequency ranges for \u03C9_2 and \u03C9_3 by uMath.sing the fractional format n \u03C0/d. The left two boxes determine \u03C9_1 and the left two boxes determine \u03C9_2"); 

      filterGroup = new ButtonGroup();
      filterGroup.add(lowBut); lowBut.addItemListener(new MyItemListener());
      filterGroup.add(rampBut); rampBut.addItemListener(new MyItemListener());
      filterGroup.add(bandBut); bandBut.addItemListener(new MyItemListener());
      filterGroup.add(highBut); highBut.addItemListener(new MyItemListener());

      omegaGroup = new ButtonGroup();
      //omegaGroup.add(periodsBut);
      omegaGroup.add(omegaBut);
      omegaGroup.add(fracBut);


       num1Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,0,100);
       num1Bar.setValue(num1);
       num1Bar.setUnitIncrement(1);

       num2Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,1,100);
       num2Bar.setValue(num2);
       num2Bar.setUnitIncrement(1);

       den1Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,1,100);
       den1Bar.setValue(den1);
       den1Bar.setUnitIncrement(1);

       den2Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,1,100);
       den2Bar.setValue(den2);
       den2Bar.setUnitIncrement(1);
       
       num1Text = new JTextField(2); num2Text = new JTextField(2);
       den1Text = new JTextField(2); den2Text = new JTextField(2);
       num1Text.setText(""+num1); num2Text.setText(""+num2);
       den1Text.setText(""+den1); den2Text.setText(""+den2); 

       num1Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_num1(num1Text.getText());}} ); 
       num2Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_num2(num2Text.getText());}} ); 
       den1Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_den1(den1Text.getText());}} ); 
       den2Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_den2(den2Text.getText());}} ); 

       omega0Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,0,314);
       omega0Bar.setValue(0); omega0Bar.setUnitIncrement(1);
      
       omega1Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,1,314);
       omega1Bar.setValue(314/6); omega1Bar.setUnitIncrement(1);

       omega2Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,0,314);
       omega2Bar.setValue(0); omega2Bar.setUnitIncrement(1);
      
       omega3Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,2,1,314);
       omega3Bar.setValue(314/6); omega3Bar.setUnitIncrement(1);



       extFore = new JCheckBox("Plot Future"); extFore.setHorizontalTextPosition(JMenuItem.LEFT);
       extFore.setToolTipText("Plot the future simulated data when forecasting computed.");
       extFore.setSelected(false); 
       extFore.addItemListener(new MyItemListener());
       ext_forecast = false;
 
       spec_densBox = new JCheckBox("Spectral Density"); spec_densBox.setHorizontalTextPosition(JMenuItem.LEFT);
       spec_densBox.addItemListener(new MyItemListener()); 
       /*period1Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,10,1,300);
       period1Bar.setValue(p1); period1Bar.setUnitIncrement(1);
      
       period2Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,10,1,300);
       period2Bar.setValue(p2); period2Bar.setUnitIncrement(1);
       */

       omega0Text = new JTextField(3); omega1Text = new JTextField(3); 
       //period1Text = new JTextField(2); period2Text = new JTextField(2); 
       omega0Text.setText(""+df.format(w0)); omega1Text.setText(""+df.format(w1));
       omega0Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_omega0(omega0Text.getText());}} ); 
       omega1Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_omega1(omega1Text.getText());}} ); 

           
       omega2Text = new JTextField(3); omega3Text = new JTextField(3); 
       //period1Text = new JTextField(2); period2Text = new JTextField(2); 
       omega2Text.setText(""+df.format(w2)); omega3Text.setText(""+df.format(w3));
       omega2Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_omega2(omega2Text.getText());}} ); 
       omega3Text.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_omega3(omega3Text.getText());}} );            


      fixBand = false;         
      fixBandCheck = new JCheckBox("Fix Bandpass width");  fixBandCheck.setHorizontalTextPosition(JMenuItem.LEFT);
      fixBandCheck.setSelected(false); fixBandCheck.addItemListener(new MyItemListener());

       AdjustmentListener al = new AdjustmentListener()  {
        public void adjustmentValueChanged(AdjustmentEvent e) {
          
          int curVal; int dVal; double temp;
          if(e.getAdjustable() == omega0Bar)
          {               
               if(.01*omega0Bar.getValue() < w1)
               {
                 curVal = omega0Bar.getValue();
                 dVal = curVal - (int)(100*w0);

                 w0 = .01*curVal;               
                 omega0Text.setText(""+df.format(w0));
 
                 if((band || high) && fixBand) //adjust upper as well
                 {
                    curVal = omega1Bar.getValue(); 
                    temp =  .01*(curVal + dVal);
                    if(temp < 3.14)
                    {w1 = .01*(curVal + dVal);
                    omega1Bar.setValue(curVal + dVal);
                    omega1Text.setText(""+df.format(w1));}
                 }
                 tfilter.setBand(w0,w1);
                 if(autoComp){computeFilterNow();}
                 cutoff0 = w0; cutoff = w1;
               }
          } 
          else if(e.getAdjustable() == omega1Bar)
          {
               if(.01*omega1Bar.getValue() > w0)
               {

                 curVal = omega1Bar.getValue();
                 dVal = curVal - (int)(100*w1);
                 w1 = .01*curVal;                                   
                 omega1Text.setText(""+df.format(w1));

                 if((band || high) && fixBand) //adjust upper as well
                 {
                    curVal = omega0Bar.getValue(); temp = .01*(curVal + dVal);
                    if(temp > 0.0)
                    {w0 = temp;
                    omega0Bar.setValue(curVal + dVal);
                    omega0Text.setText(""+df.format(w0));}
                 }
                 tfilter.setBand(w0,w1); 
                 cutoff0 = w0; cutoff = w1;
                 if(autoComp){computeFilterNow();}
               }  
          } 
          else if(e.getAdjustable() == omega2Bar)
          {               
               if(.01*omega2Bar.getValue() < w3)
               {
                 curVal = omega2Bar.getValue();
                 dVal = curVal - (int)(100*w2);

                 w2 = .01*curVal;               
                 omega2Text.setText(""+df.format(w2));
 
                 if(high && fixBand) //adjust upper as well
                 {
                    curVal = omega3Bar.getValue(); 
                    temp =  .01*(curVal + dVal);
                    if(temp < 3.14)
                    {w3 = .01*(curVal + dVal);
                    omega3Bar.setValue(curVal + dVal);
                    omega3Text.setText(""+df.format(w3));}
                 }
                 tfilter.setBand2(w2,w3);
                 if(autoComp){computeFilterNow();}
                 //cutoff0 = w0; cutoff = w1;
               }
          } 
          else if(e.getAdjustable() == omega3Bar)
          {
               if(.01*omega3Bar.getValue() > w2)
               {

                 curVal = omega3Bar.getValue();
                 dVal = curVal - (int)(100*w3);
                 w3 = .01*curVal;                                   
                 omega3Text.setText(""+df.format(w3));

                 if(high && fixBand) //adjust upper as well
                 {
                    curVal = omega2Bar.getValue(); temp = .01*(curVal + dVal);
                    if(temp > 0.0)
                    {w2 = temp;
                    omega2Bar.setValue(curVal + dVal);
                    omega2Text.setText(""+df.format(w2));}
                 }
                 tfilter.setBand2(w2,w3); 
                 //cutoff0 = w0; cutoff = w1;
                 if(autoComp){computeFilterNow();}
               }  
          }
/*          else if(e.getAdjustable() == num1Bar)
          {
              if((1.0*num1Bar.getValue())/(1.0*den1) < (1.0*num2)/(1.0*den2))
              {
                 num1 = num1Bar.getValue();
                 num1Text.setText(""+num1);
                 w0 = 1.0*num1*Math.PI/(1.0*den1);
                 omega0Text.setText(""+df.format(w0));
                 tfilter.setBand(w0,w1);
                 if(autoComp){computeFilterNow();}
                 cutoff0 = w0; cutoff = w1;
              }
          } 
          else if(e.getAdjustable() == den1Bar)
          {
              if((1.0*num1/(1.0*den1Bar.getValue())) < (1.0*num2)/(1.0*den2))
              {
                 den1 = den1Bar.getValue();
                 den1Text.setText(""+den1);
                 w0 = 1.0*num1*Math.PI/(1.0*den1);
                 omega0Text.setText(""+df.format(w0));
                 tfilter.setBand(w0,w1); 
                 if(autoComp){computeFilterNow();}
                 cutoff0 = w0; cutoff = w1;
              }
          } 
          else if(e.getAdjustable() == num2Bar)
          {
  
              if((1.0*num1)/(1.0*den1) < (1.0*num2Bar.getValue())/(1.0*den2))
              {
                 num2 = num2Bar.getValue();
                 num2Text.setText(""+num2);
                 w1 = 1.0*num2*Math.PI/(1.0*den2);
                 omega1Text.setText(""+df.format(w1));
                 tfilter.setBand(w0,w1);
                 if(autoComp){computeFilterNow();}
                 cutoff0 = w0; cutoff = w1;
              }
          } 
          else if(e.getAdjustable() == den2Bar)
          {
              if((1.0*num1)/(1.0*den1) < (1.0*num2)/(1.0*den2Bar.getValue()))
              {
                 den2 = den2Bar.getValue();
                 den2Text.setText(""+den2);
                 w1 = 1.0*num2*Math.PI/(1.0*den2);
                 omega1Text.setText(""+df.format(w1));
                 tfilter.setBand(w0,w1);
                 if(autoComp){computeFilterNow();}
                 cutoff0 = w0; cutoff = w1;
              }
          }*/
          if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
          setFilter(); if(!autoComp){compute.setEnabled(true);}
       }
    };


         //num1Bar.addAdjustmentListener(al);
         //num2Bar.addAdjustmentListener(al);
         //den1Bar.addAdjustmentListener(al);
         //den2Bar.addAdjustmentListener(al);
         omega0Bar.addAdjustmentListener(al);
         omega1Bar.addAdjustmentListener(al);
         omega2Bar.addAdjustmentListener(al);
         omega3Bar.addAdjustmentListener(al);


   }


   public void computeFilterNow()
   {
     if(!sweepMode)
     {
      computeFilter = true;
      mdfa.computeFilterGeneral(computeFilter,false);
      updatePlots(false,true);
      compute.setEnabled(false);
     }
     else
     {
       out_of_sample_tradingSweep(getAccount_canvas().nbackObs, getAccount_canvas().n_insamp, getAccount_canvas().n_outsamp);      
     }
   }


    public void checkSamp()
    {                 
         //System.out.println(K + "  " + nsamp);
         if(K > nsamp) //periodogram e jadid, mikhahim ke mimune ye periodogram hadd e aghal K baashe
         { 
           //System.out.println(K + "  " + nsamp);
           nsamp = K;
           cPeriod = false; sampBar.setValue(K);  cPeriod = true;
           sampText.setText(""+nsamp);
           mdfa.set_Samps(nsamp); zpc.set_Samps(nsamp);
           //mdfa.computeSampleIns();
           
           //period_canvas.go();
         }
    }

    public void reset()  // ---- resets GDP data to 300 pts with 5 series
    {
      n_obs = 300; n_rep = 5; K = n_obs/2; K1 = K+1;
      nObsBar.setValue(n_obs); 
      nText.setText(""+ n_obs);   
      mdfa.set_nobs(n_obs); mdfa_canvas.setNobs(n_obs); filter_canvas.setK(K);
      tfilter.setK(K); tfilter.recomputeGamma();
      simulateIMDFA(true);    
    }
  


    public void changeSeed(int _s) {seed = _s; simulateIMDFA(true);}

    public void plotAmplitude()
    {
     filter_canvas.setGamma_hat(mdfa.amp_filter,K,n_rep);

     if(filterPlot[1].isSelected() || filterPlot[2].isSelected())
     {
      for(int i=0; i < n_rep; i++)  //---- get frf
      {   
        for(int k=0;k<=K;k++)
        {
           gamma_zpc[K1*i+k] = zpc.amp[k];  
           gamma_hybrid[K1*i+k] = zpc.amp[k]*mdfa.amp_filter[K1*i+k];                  
        }
      }
     filter_canvas.setGamma_zpc(gamma_zpc);
     filter_canvas.setGamma_hybrid(gamma_hybrid);
     }
    }
 

    public void plotTimeDelay()
    {
    
     filter_canvas.setGamma_hat(mdfa.time_delay,K,n_rep); 
     if(filterPlot[1].isSelected() || filterPlot[2].isSelected())
     {
      for(int i=0; i < n_rep; i++)  //---- get frf
      {   
        for(int k=0;k<=K;k++)
        {
           gamma_zpc[K1*i+k] = zpc.time_delay[k];    
           gamma_hybrid[K1*i+k] = zpc.time_delay[k]*mdfa.time_delay[K1*i+k];             
        }
      }
      filter_canvas.setGamma_zpc(gamma_zpc);
      filter_canvas.setGamma_hybrid(gamma_hybrid);
     }
    }   
 
    public void plotGamma()
    {
      filter_canvas.setGamma_hat(mdfa.gamma_hat,K,n_rep);

      if(filterPlot[1].isSelected() || filterPlot[2].isSelected())
      {
       for(int i=0; i < n_rep; i++)  //---- get frf
       {   
        for(int k=0;k<=K;k++)
        {
           gamma_zpc[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k]);    
           gamma_hybrid[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k])*mdfa.gamma_hat[K1*i+k];             
        }
       }
       filter_canvas.setGamma_zpc(gamma_zpc);
       filter_canvas.setGamma_hybrid(gamma_hybrid);
      }
    }   

    public void updateFreq(int i, boolean sel)
    {filter_canvas.setPlots(i,sel); if(i==11){filter_canvas.setPlots(n_rep+1,sel);}}

    public void updatePeriod(int i, boolean sel)
    {period_canvas.setPlots(i,sel);}    

    public void updateCoeffs(int i, boolean sel)
    {coef_canvas.setPlots(i,sel);}

    public void updateTime(int i, boolean sel)
    {mdfa_canvas.setPlots(i,sel);}

    //------------------------------------------------------------
    public void computeTargetGamma(int l1)
    {
       int i,k; bc = new double[L1+1];
       double[] Gamma = new double[K1];
       
       double sum;
       bc[0] = cutoff/Math.PI; sum= bc[0];
       for(i=1;i<=l1;i++)
       {bc[i] = (1/Math.PI)*Math.sin(cutoff*i)/i; sum= sum+bc[i];} 
       for(i=0;i<=l1;i++)
       {bc[i] = bc[i]/(sum+(sum-cutoff/Math.PI));}     
       for(k=0; k<=K;k++)
       {       
        if((k*Math.PI/K) <= cutoff) {Gamma[k] = 1.0;}
        else {Gamma[k] = 0.0;}

       }           
       
       mdfa.set_Gamma(Gamma);  //--- set the new gammas
       //for(i=0;i<K1;i++) {System.out.println(Gamma[i]);}
       //mdfa.computeFilterGeneral(computeFilter);   //--- compute the new filter with target gamma   
       filter_canvas.setGamma(Gamma);
       //updatePlots(false);
    }

    public void reComputeSymetricSignal()
    {
      int i,l; double sum;
      
      if(simulate)
      {
       symfilt = new double[n_obs - (L-1)];
       for(i=L1+(L-1);i<n_obs;i++)
       {
      	sum = 0.0;
        for(l=0;l<L1;l++)
        {sum = sum + bc[l]*extend[i+l];} 
        for(l=1;l<L1;l++)
        {sum = sum + bc[l]*extend[i-l];}
        symfilt[i-(L1+(L-1))] = sum;         
       }
       mdfa_canvas.setSymSignal(symfilt); 
      }      
    }


    public void setL(int _L)
    {
      L = _L; 
      mdfa.set_L(L); 
      mdfa_canvas.setL(L);
  
      if(L1 < 36 || simulate)
      {
        mdfa.computeFilterGeneral(computeFilter,false);  
        mdfa_canvas.setSymShift(L+L1); reComputeSymetricSignal();
        updatePlots(false,true);
        if(!simulate) {mdfa_canvas.win_shift = 36 + (L_fixed - L);}
      } 
      else if(!simulate)
      {
        adjustSymmetricDataPoints();
      }
    }

    public void setl1(int _L)
    {
      L1 = _L; tfilter.setL1(l1); tfilter.recomputeGamma(); 
      if(useX13filters) {tfilter.setL1(L1); tfilter.setGeneralSymmetric(Gamma_x13);}
    }
 
    public void setCutoff(double d0, double _d)
    {
       cutoff = _d; cutoff0 = d0; 
       if(cutoff0 > 0.0) {i1const.setSelected(false); i1const.setEnabled(false);}
       else {i1const.setEnabled(true);}
       mdfa.set_cutoff(_d); mdfa.set_cutoff0(d0); 
       computeTargetGamma(L1); reComputeSymetricSignal(); 
       period_canvas.setCutoffs(d0, cutoff);
    }

    public void setTurningLambda3()
    {
        mdfa.setLambda3(lambda_3);
        mdfa.computeFilterGeneral(computeFilter,false); 
        updatePlots(false,true);
    }

    public void setDiffOnestep()
    {
        mdfa.setDiffOnestep(onestep_diff);
        mdfa.computeFilterGeneral(computeFilter,false); 
        updatePlots(false,true);
    }    
    
    public void setOnestep()
    {
        mdfa.setOnestep(onestep);
        mdfa.computeFilterGeneral(computeFilter,false); 
        updatePlots(false,true);
    }    
    

    public void setFilter()
    {
      if(!useX13filters)
      {
        mdfa.set_cutoff(w1); 
        if(bandBut.isSelected()){mdfa.set_cutoff0(w0);}
        else if(lowBut.isSelected() || rampBut.isSelected()) {w0 = 0.0; mdfa.set_cutoff0(w0);}
        else if(highBut.isSelected()) {mdfa.set_cutoff0(w0); mdfa.set_cutoff(w3);}
        
        set_Gamma(tfilter.Gamma); setBC(tfilter.bsym);      
        if(rampBut.isSelected() || lowBut.isSelected()) {period_canvas.setCutoffs(0, tfilter.w1);}
        else if(highBut.isSelected()) {period_canvas.setCutoffs2(tfilter.w0, tfilter.w1, tfilter.w2, tfilter.w3);}
        else{period_canvas.setCutoffs(tfilter.w0, tfilter.w1);}
        if(b_outof_sample){mdfa_canvas.setSymShift(L+L1);}
      }
      else
      { 
        mdfa.set_cutoff(3.14);  mdfa.set_cutoff0(0.0);     
        set_Gamma(Gamma_x13); setBC(tfilter.bsym);
        period_canvas.setCutoffs(0.0, 3.14); mdfa_canvas.setSymShift(L+L1);
      }     
 
    }

    //----------------------------------------------------------

    //------------------------------------------------------------

    public void setNobs(int n)
    {
             
         if(simulate)
         {                          
           if(reCompFilter)
           {
              n_obs = n; K = n_obs/2; K1 = K+1;
              mdfa.set_nobs(n_obs); mdfa_canvas.setNobs(n_obs); filter_canvas.setK(K);           
              checkSamp(); cPeriod = false; sampBar.setMinimum(K); cPeriod = true;
	      simulateIMDFA(true); tfilter.setK(K); tfilter.recomputeGamma(); setFilter();
                 
           } 
           else
           {
              n_obs = n; mdfa_canvas.setNobs(n_obs);
              checkSamp(); cPeriod = false; sampBar.setMinimum(K); cPeriod = true;
              simulateIMDFA(true);                            
           }
         }
         else
         {
           n_obs =n; K = n_obs/2; K1 = K+1; checkSamp(); 
           cPeriod = false; sampBar.setMinimum(K); cPeriod = true;
           mdfa.set_nobs(n_obs); mdfa_canvas.setNobs(n_obs); filter_canvas.setK(K); 
           tfilter.setK(K); 
           if(useX13filters) {tfilter.setL1(L1); tfilter.setGeneralSymmetric(Gamma_x13);}
           else {tfilter.recomputeGamma();}
           setFilter();
         }
         tfplotPanel.setNobs();  
         
         //System.out.println("NOBS = " + n_obs);
         mdfa_canvas.setPreferredSize(new Dimension(3*n_obs, 400)); 
         timeScrollPane.setViewportView(mdfa_canvas);                 
   } 




    public void setTimeSeries(double[] _ts, int _n, int _nr)
    {  
       
      n_obs = _n; n_rep = _nr; K = n_obs/2; K1 = K+1; updateTexts();
      mdfa.set_tseries(_ts,_n,_nr); mdfa_canvas.setNobs(n_obs); filter_canvas.setK(K);
      

      checkSamp(); sampBar.setMinimum(K);
	
      if(!useSARIMA)
      {tfilter.setK(K); tfilter.recomputeGamma(); setFilter();}
      else if(useSARIMA && useX13filters)
      {computeX13TargetFilter(); tfilter.setGeneralSymmetric(Gamma_x13); setFilter(); }
      else
      {tfilter.setK(K); tfilter.recomputeGamma(); setFilter();}      

      if(activateEstimationCheck.isSelected() && autoUpdatesCheck.isSelected())
      {computeSpectralDensityEstimateSeries();}
      else
      {mdfa.computeFilterGeneral(reCompFilter,false);}      
      
      

      //if(filterPlot[1].isSelected() || filterPlot[2].isSelected()) {computeZPCFilter();}
      
      updatePlots(true,true); 

      //if(!setup) {setZPC();}
    }

    public void setData(double[] _ts, int _n, int _r, int d, int D)
    {
  
           //----- update sizes -------------------
           n_rep = _r; n_obs = _n;  
           dd = d; DD = D;             
           setTimeSeries(_ts, _n, n_rep);
  
           //----- none of the simulated stuff-------
           noneCheck.setSelected(true); //---- disable the simulate IMDFA stuff -------
           simulate = false; useSARIMA = false;
           nrepBar.setEnabled(simulate);                 
    }


    public void setLag(int _lag)
    {
      if(reCompFilter)
      {
       Lag = _lag; 
       mdfa.set_lag(Lag); 
       mdfa.computeFilterGeneral(computeFilter,false);
       mdfa_canvas.setLag(Lag);
       updatePlots(false,true);       
      }
    }

    public void setX13Data(double[] _ts, int _n)
    {
       if(useSARIMA)
       {   
          //----- update sizes -------------------
          n_rep = 1; n_obs = _n; updateTexts();      
          //---- difference data and set into system-----                              
          setTimeSeries(_ts, _n, n_rep);
       } 
    }

    public void setUseSarima(boolean sar)
    {useSARIMA = sar;}

    public void setSmooth(double _l) 
    {
     smooth = _l; mdfa.setRegularization(smooth, decay, decay2, cross);
     mdfa.computeFilterGeneral(computeFilter,false); 
     updatePlots(false,true);
    }
    public void setDecay(double _l) 
    {
     decay = _l; mdfa.setRegularization(smooth, decay, decay2, cross);
     mdfa.computeFilterGeneral(computeFilter,false); 
     updatePlots(false,true);
    }
    public void setDecay2(double _l) 
    {
     decay2 = _l; mdfa.setRegularization(smooth, decay, decay2, cross);
     mdfa.computeFilterGeneral(computeFilter,false); 
     updatePlots(false,true);
    }
    public void setCross(double _l) 
    {
     cross = _l; mdfa.setRegularization(smooth, decay, decay2, cross); 
     mdfa.computeFilterGeneral(computeFilter,false); 
     updatePlots(false,true);
    }

    public void setLambda(double _l) 
    {
     lambda = _l; mdfa.set_lambda(lambda); 
     mdfa.computeFilterGeneral(computeFilter,false);  
     updatePlots(false,true);
    }
  
    public void setLambdaExp(double _l, double _e)
    {
       lambda = _l;    mdfa.set_lambda(lambda);
       expweight = _e; mdfa.set_exp(expweight);        

       expText.setText(""+df.format(expweight)); 
       lambdaText.setText(""+df.format(lambda));

       mdfa.computeFilterGeneral(computeFilter,false);  
       updatePlots(false,true);
    }

    public void setShift(double _s)
    {
      shift = _s; mdfa.set_shift(shift);
      mdfa.computeFilterGeneral(computeFilter,false);  
      updatePlots(false,true);
    }


    public void setExp(double _l) 
    {
     expweight = _l; mdfa.set_exp(expweight); 
     mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);
    }

    public void setDD(int _dd, int _DD) 
    {dd = _dd; DD = _DD; mdfa.set_dd(dd); mdfa.set_DD(DD); mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);}

    public void setDDNoCompute(int _dd, int _DD) 
    {dd = _dd; DD = _DD; mdfa.set_dd(dd); mdfa.set_DD(DD); mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);}

    public void setBConstraints(int _i1, int _i2)
    {i1 = _i1; i2 = _i2; mdfa.set_bconstraints(_i1, _i2); mdfa.computeFilterGeneral(computeFilter,false); updatePlots(false,true);}


    //----------------------update the plots----------
    public void updatePlots(boolean tsp, boolean freqe)
    {     
      //-----if update of time series needed----
      int fl,i,L_price,l; double sum;

     //update plots of coefficients
     if(autoUpdate.isSelected() && b_outof_sample)
     {
       
       for(i=0;i<=mdfa.n_rep;i++) {updateCoeffPanel.b_plots[i] = true;} 
       updateCoeffPanel.b_plots[1] = false;
       if(uni_updateCheck.isSelected())
       {
         updateCoeffPanel.b_plots[1] = true;
         updateCoeffPanel.setBCoeffs(mdfa.b_update, mdfa.L_update, 1);
         updateAmpPanel.setGammas(mdfa.amp_update, mdfa.K_update, 2);
       }
       else
       {
        updateCoeffPanel.setBCoeffs(mdfa.b_update, mdfa.L_update, mdfa.n_rep);
        updateAmpPanel.setGammas(mdfa.amp_update, mdfa.K_update, mdfa.nrep_update);
       }
       updateAmpPanel.setMSEValue(mdfa.mse_update);
     }      
      
     if(mdfaPlotPane.getSelectedIndex() == 5 && sweepMode) // --- if in sweepmode of financial trading
     {
       out_of_sample_tradingSweep(getAccount_canvas().nbackObs, getAccount_canvas().n_insamp, getAccount_canvas().n_outsamp);    
     }
     else if(computeFilter)
     {
      if(tsp)
      { 
        mdfa_canvas.setTseries(mdfa.tseries, n_obs, n_rep);
        period_canvas.setPeriodograms(mdfa.period_xf, mdfa.period_hat, K, n_rep);
      }
      mdfa_canvas.setRTSignal(mdfa.xf, n_obs);      

      if(freqe)
      {
       
       if(gammaBut.isSelected())
       {filter_canvas.setGamma_hat(mdfa.gamma_hat, K, n_rep);}
       if(amplitudeBut.isSelected())
       {filter_canvas.setGamma_hat(mdfa.amp_filter, K, n_rep);}
       if(timeDelayBut.isSelected())
       {filter_canvas.setGamma_hat(mdfa.time_delay, K, n_rep);}
       filter_canvas.go();
 
       if(!tsp)
       {period_canvas.setPeriodogramXf(mdfa.period_xf, nsamp);}
      }
      //System.out.println(mdfa.b.length + "  " + L + "  " + n_rep);
      reComputeSymetricSignal(); coef_canvas.setBCoeffs(mdfa.b,L,n_rep);
      mdfa_canvas.go(); period_canvas.go();
      mdfaText.setText("IC/DF = " + df.format(Math.abs(mdfa.criteria)) + ", " + df.format(Math.abs(mdfa.degrees)));
      changeHealthColor();
       
  
      if(filterPlot[1].isSelected() || filterPlot[2].isSelected())
      {
         if(reCompFilter && autoCompZPC) {updateZPC();}
 
         injectZPCGeneSimple(); 
         if(mdfaPlotPane.getSelectedIndex() == 0) {mdfa_canvas.plotHybrid(plot_number);}
         else if(mdfaPlotPane.getSelectedIndex() == 1) {filter_canvas.plotHybrid(plot_number);}
         else if(mdfaPlotPane.getSelectedIndex() == 2) {period_canvas.plotHybrid(plot_number);}
         else {mdfa_canvas.plotHybrid(plot_number);}
      }
      //----- Update all the diagnostics and score the filter if available-------------
      //if(b_outof_sample && (n_out_samp == 0))
      //{filterScore();}
      //updateDiagnosticPanel();

      //------ set new price data, and new signal -------------
      if(trading)
      {
        int hf_fl = 0;
        fl = n_obs - (L-1);
        double[] _price = new double[fl];
        double[] _sig = new double[fl];            
        double[] logreturns = new double[fl];
        double[] _hf_price = null;

        for(i=0;i<fl;i++)
        {_price[i] = t_price[0][n_out_samp+(L-1)+i];}
           
        if(hf_price_set && trading_func == 4) //adjust the hf_price to match the lf price
        {
          hf_fl = fl*hf_period;
          _hf_price = new double[hf_fl];
         for(i=0;i<hf_fl;i++) 
         {
           _hf_price[i] = hf_price[hf_period*(L-1) + i]; 
           if((hf_period*(L-1) + i + 1) >= hf_price_length) break;
         }
        }
// 

        //System.out.println("lengths = " + mdfa.xf.length + "  " +  zpc.xf_hybrid.length + "  " + zpc.xf_zpc.length);
         
           
        for(i=0;i<fl;i++)
        {
          if(filterPlot[0].isSelected())
          {_sig[i] = mdfa.xf[i];}
          else if(filterPlot[1].isSelected())
          {_sig[i] = zpc.xf_zpc[i];}
          else if(filterPlot[2].isSelected())
          {_sig[i] = zpc.xf_hybrid[i];}
          
          logreturns[i] = mdfa.tseries[(L-1)+i];
        }       
        
        //compute price filter
        if(trading_func == 2 && mdfa.price_filter) 
        {
          price_indicator = new double[fl];
          L_price = mdfa.L_price;
                     
          for(i = 0; i < fl; i++)
          {
           sum=0.0;
           for(l=0; l < L_price; l++)
           {
            sum = sum + mdfa.b_price_filter[l]*t_price[0][L2+(L-1) - l + i];
           }
           price_indicator[i] = sum;
          }
          
        }

        //System.out.println(""+continuousUpdate.isSelected() + " " + plotUpdateBox.isSelected());
        if(reCompFilter)
        {
          if(n_out_samp == 0)
          {
            for(i=0;i<fl;i++) {fixed_signal[fixed_signal.length-i-1] = _sig[_sig.length-i-1];}          
          }
          else if(n_out_samp > 0)
          {
            for(i=0;i<fl-1;i++) {_sig[fl-i-2] = fixed_signal[fixed_signal.length-i-1];} //applies the old signal for      
            fixed_signal[fixed_signal.length-1] = _sig[fl-1];
            
            for(i=0;i<fl;i++) {fixed_signal[fixed_signal.length-i-1] = _sig[fl-i-1];}
            //System.out.println(_sig[fl-1]);
          }        
        }
        else if(continuousUpdate.isSelected() && plotUpdateBox.isSelected())
        {
          if(n_out_samp == 0)
          {
            adapt_signal = new ArrayList<Double>(); 
            for(i=0;i<fl;i++) {fixed_signal[fixed_signal.length-i-1] = _sig[_sig.length-i-1];}
            prev_out_samp = 0;
          }
          else if(n_out_samp > 0)
          {
           if(prev_out_samp == n_out_samp-1)
           {
            adapt_signal.add(_sig[_sig.length-1]);
//             for(i=0;i<fl-1;i++) {_sig[fl-i-2] = fixed_signal[fixed_signal.length-i-1];} //applies the old signal for      
//             fixed_signal[fixed_signal.length-1] = _sig[fl-1];
//             
//             for(i=0;i<fl;i++) {fixed_signal[fixed_signal.length-i-1] = _sig[fl-i-1];}
            //System.out.println(_sig[_sig.length-1]);
            
            for(i=0;i<n_out_samp;i++)
            {_sig[_sig.length-i-1] = adapt_signal.get(adapt_signal.size()-i-1);} // System.out.println(_sig[_sig.length-i-1]);}
            //System.out.println("\n");
            prev_out_samp++;
           }           
          }   
        }
        
        
        
        if(trading_func == 0)
        {insampleTradingDiff(_price, _sig, fl);}
        else if(trading_func == 1)
        {insampleTradingLogPrice(_price, _sig, fl);} 
        else if(trading_func == 2)
        {trading_function_Duplex(_price, price_indicator, _sig, fl);}
        else if(trading_func == 3)
        {trading_function_Duplex_Price(_price, price_indicator, _sig, fl);}
        else if(trading_func == 4)
        {multiscaleFilterSimple(_price, _sig, fl, _hf_price, hf_fl, hf_period); getAccount_canvas().setHFPrice(_hf_price);}
        

        postTradingStatistics();
 

        
        
        
        getAccount_canvas().setAccount(account);
        getAccount_canvas().setSignal(_sig);
        getAccount_canvas().setPrice(_price);
        getAccount_canvas().setLogReturns(logreturns); 
        getAccount_canvas().setBuySellLines(xf_turn_val);
        if(trading_func == 2 || trading_func == 3) {getAccount_canvas().setFilteredPrice(price_indicator); filteredPrice.setEnabled(true);}

        //if(trading_func == 4 && hf_price_set)
        //{account_canvas.setHFPrice(_hf_price);}

        getAccount_canvas().go();
        //t_signal = null;

        computeOutOfSampleStats();

        if(compute_error)
        {computeOutOfSampleError();}
        
      }

     }
     else if(b_outof_sample && !reCompFilter)
     {mdfa_canvas.setRTSignal(mdfa.xf, n_obs);}
     

     
     
    }

    public void changeHealthColor()
    {
       int index = 0;
       double degreesMax = 1.0*n_rep*L;
       double degreesMin = degreesMax/3.0; 
       double amp = mdfa.max_amp;
       double degrees = mdfa.degrees;
       Color col;       

       double health = degrees*(1.0 + (Math.abs(1.0 - 1.0/(amp+.00001))));       

       if(health >= degreesMax) {col = filter_health[0];}
       else if(health < degreesMax/3.0) {col = filter_health[124];}
       else
       {
         index = (int)(124.0*((degreesMax - health)/(degreesMax - degreesMin)));
       }        
       col = filter_health[index];
    
       mdfaText.setBackground(col);

       //Minimum = degreesMax/3; Maximum = degreesMax
    }

    public void updateTrading(int c)
    {
      int fl,i; //System.out.println("c = " + c);
      if(trading)
      {
        
        fl = n_obs - (L-1);
        double[] _price = new double[fl];
        double[] _sig = new double[fl];            
        double[] logreturns = new double[fl];
        

        for(i=0;i<fl;i++)
        {_price[i] = t_price[0][L2+(L-1)+i];}
           
        //System.out.println("lengths = " + mdfa.xf.length + "  " +  zpc.xf_hybrid.length + "  " + zpc.xf_zpc.length);
         
           
        for(i=0;i<fl;i++)
        {
          if(c==0)
          {_sig[i] = mdfa.xf[i];}
          else if(c==1)
          {_sig[i] = zpc.xf_hybrid[i];}
          else if(c==2)
          {_sig[i] = zpc.xf_zpc[i];}
          
          logreturns[i] = mdfa.tseries[(L-1)+i];
        }       
 
       
        //System.out.println("NEW PRICE COMPARISON");
        //for(i=0;i<fl;i++)
        //{System.out.println(_price[i] + "  " + int_price[i]);}
        //System.out.println(); 
  

        if(trading_func == 0)
        {insampleTradingDiff(_price, _sig, fl);}
        else if(trading_func == 1)
        {insampleTradingLogPrice(_price, _sig, fl);} 
        //updateTradingInterface();
        postTradingStatistics();
 
        getAccount_canvas().setAccount(account);
        getAccount_canvas().setSignal(_sig);
        getAccount_canvas().setPrice(_price);
        getAccount_canvas().setLogReturns(logreturns); 
        getAccount_canvas().setBuySellLines(xf_turn_val);
        
        if(trading_func >= 1) {getAccount_canvas().setDerivativeSignal(d_signal);}
        getAccount_canvas().go();

        computeOutOfSampleStats();

      }    
    
    
    }
    
    
    public void set_Gamma(double[] _gamma)
    {
      mdfa.set_Gamma(_gamma); 
      filter_canvas.setGamma(_gamma); 
      zpcFreqCanvas.setGamma(_gamma); 
    } 

    public void setBC(double[] _bc)
    {
      int le = _bc.length;
      bc = new double[le]; 
      System.arraycopy(_bc, 0, bc, 0, le); 
      //coef_canvas.setSymB(bc,le);
    }

    public void setSarimaCheck(boolean x)
    {   
        int i;
        simulate = !x;  useSARIMA = x; timePlot[11].setEnabled(simulate); mdfa_canvas.setPlots(7, simulate);
        setX13Data(x13data, x13nobs); setEnableX13(useSARIMA);  
 
        double mean = 0.0; 
        for(i = 0; i<x13nobs;i++) {mean = mean + x13data[i];} 
        mean = mean/x13nobs;
        for(i = 0; i<x13nobs;i++) {x13data[i] = x13data[i] - mean;} 

        boolean[] rep = new boolean[series_max];
        for(i=0;i<series_max;i++) rep[i] = false;
        rep[0] = true;

        ArrayList<double[]> list = new ArrayList<double[]>();
        list.add(x13data);

        inputSimulatedData(x13data, list, rep, 0);         
        setEnableX13(x);
    }


    public void simulateIMDFA(boolean comp)
    {
      if(simulate) 
      {         

         //---------  set-up simulation parameters ------------------------
         int i,j,l; int N = n_obs; int len = L1; int burnin = 150; int S = 12; 
         int len1 = N+2*len+1; double sum;
         double[] params;
         int[] dim = new int[6];
         dim[0] = 1; dim[1] = 0; dim[2] = 0; dim[3] = 0; dim[4] = 0; dim[5] = 0;
         int n_params = dim[0] + dim[2] + dim[3] + dim[5] + 1;  
         params = new double[n_params];
         params[0] = -.8; params[1] = 1.0;
         //----------------------------------------------------------------

         extend = new double[2*L1+n_obs]; fore_extend = new double[36];
         symfilt = new double[n_obs - (L-1)];
                   
         double[] z1; double[] z = new double[4];
         double[] y; double[] y1; double[] gdp;
         double[] y2; double[] y3;
         SARIMAmodelJava model = new SARIMAmodelJava(len1+1, burnin);     
         model.SetSARIMAparams(params, dim, n_params, S);
         y = model.sampleSARIMAmodel(seed);
         params[0] = -.99;
         model.SetSARIMAparams(params, dim, n_params, S);
         y1 = model.sampleSARIMAmodel(seed+256);
         y2 = model.sampleSARIMAmodel(seed+301);
         y3 = model.sampleSARIMAmodel(seed+7216);

         z1 = cumsum(y,len1); //integrate series y 
         gdp = new double[n_rep*N];
 
         for(i=1; i<len1; i++)
         {   
           z[0] = y[i];
           z[1] = (z1[i] + 10*y1[i]) - (z1[i-1] + 10*y1[i-1]);
           z[2] = (z1[i] + 10*y2[i]) - (z1[i-1] + 10*y2[i-1]); 
           z[3] = (z1[i] + 10*y3[i]) - (z1[i-1] + 10*y3[i-1]);            
           extend[i-1] =  z[0] + z[1] + z[2] + z[3];

           if((i > L1) && (i <= (n_obs + L1)))
           {
            gdp[i-L1-1] =  z[0] + z[1] + z[2] + z[3];
            for(j = 1; j < n_rep; j++)
            {gdp[N*j + (i-L1-1)] = z[j-1];}
           }
         }

	 for(i=L1+(L-1);i< N;i++)
  	 {
      	   sum = 0.0;
           for(l=0;l<L1;l++)
           {sum = sum + bc[l]*extend[i+l];} 
           for(l=1;l<L1;l++)
           {sum = sum + bc[l]*extend[i-l];}
           symfilt[i-(L1+(L-1))] = sum;         
         }

         mdfa_canvas.setSymShift(L+L1);
         for(i=0;i<36;i++)
         {fore_extend[i] = extend[n_obs+L1+i-1];}

          mdfa_canvas.setForeExt(fore_extend, ext_forecast);

         /*for(i=0; i<N; i++)
         {   
           z[0] = y[len+i];
           z[1] = (z1[len+i] + 10*y1[len+i]) - (z1[len+i-1] + 10*y1[len+i-1]);
           z[2] = (z1[len+i] + 10*y2[len+i]) - (z1[len+i-1] + 10*y2[len+i-1]); 
           z[3] = (z1[len+i] + 10*y3[len+i]) - (z1[len+i-1] + 10*y3[len+i-1]);            
           gdp[i] =  z[0] + z[1] + z[2] + z[3];
           //System.out.println(gdp[i]);
           for(j = 1; j < n_rep; j++)
           {gdp[N*j + i] = z[j-1];}
         }*/
         mdfa_canvas.setSymSignal(symfilt); 
         
         if(!reCompFilter && comp)
         {mdfa.computeRTSE(gdp, n_obs, n_rep); updatePlots(true, false);} //-- applies the older filter, nothing recomputed
         else
         {setTimeSeries(gdp, n_obs, n_rep); } //-- recomputes everything with new data       
      }  
    }


   /*---------------------------------------------------------------------------------
     Reads multivariate data of length back_samp + in_samp + for_samp
     n_r total series, first being the target 
   -----------------------------------------------------------------------------------*/

    public void readMultiDataFile(File file, int n_r, int in_samp, int back_samp, int for_samp)
    {

       String strline; Double D; String[] tokens; String delims = "[ ]+"; int n_toks; 
       ArrayList<double[]> values = new ArrayList<double[]>();
       double sum;        
       extend = new double[in_samp+back_samp+for_samp];
       double[] row = new double[n_r];   
       double[] gdp = new double[n_rep*in_samp];
       int i = 0; int j,l; int pos=0;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); n_toks = tokens.length; 
           if(n_toks != n_r)
           {System.out.println("End of file: added values = "+values.size()); break;}

           for(i=0;i<n_toks;i++) 
           {  
             D = new Double(tokens[i]); row[i] = D.doubleValue();           
           }

           extend[pos] = row[0]; 
           if((pos > back_samp) && (pos <= (in_samp + back_samp)))
           {
            for(j = 0; j < n_toks; j++)
            {gdp[in_samp*j + (pos-back_samp-1)] = row[j];}                      
           }
           pos++;
           // for(j=0;j<n_r;j++) {System.out.print(row[j] + " ");}
           // System.out.print("\n");
         } 
         din.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}

           
         n_obs = in_samp; n_rep = n_r; nObsBar.setValue(n_obs); nText.setText(""+ n_obs);   
         nrepBar.setValue(n_rep); nrepText.setText(""+n_rep);

         simulate = false;  useSARIMA = false; timePlot[10].setEnabled(false); mdfa_canvas.setPlots(7, false);
         nObsBar.setEnabled(false); nrepBar.setEnabled(false);  setEnableX13(false);  noneCheck.setSelected(true);

         if((back_samp >= L1) && (for_samp >= L1))
         {
          for(i=L1+(L-1);i<in_samp;i++)
  	  {
      	   sum = 0.0;
           for(l=0;l<L1;l++)
           {sum = sum + bc[l]*extend[i+l];} 
           for(l=1;l<L1;l++)
           {sum = sum + bc[l]*extend[i-l];}
           symfilt[i-(L1+(L-1))] = sum;         
          }
         

         } 

         mdfa.set_nobs(n_obs); mdfa_canvas.setNobs(n_obs); filter_canvas.setK(K); tfilter.setK(K); tfilter.recomputeGamma();
         setTimeSeries(gdp, n_obs, n_rep);
 
    }


    public void setEnableX13(boolean x)
    {
       x13filter.setEnabled(x); x13trend.setEnabled(x); x13seas.setEnabled(x); x13ti.setEnabled(x); 
    }

 
    public void setupDesign()
    {
       int i;
       JPanel paramPanel1, paramPanel2;
       JPanel resPanel;
       JPanel tpPanel, filterPanel;
       
       GroupLayout paramLayout;
       new JLabel("|       |");      
       JLabel bLabel2 = new JLabel("      "); 
       JLabel bLabel3 = new JLabel(" ");
       JLabel pLabel = new JLabel("p:");
       JLabel qLabel = new JLabel("q:");

       /*----------------------------------------------------------------------------
          Regularization stuff 
       -------------------------------------------------------------------------------*/
       JLabel smoothLabel = new JLabel("Smooth:"); 
       JLabel decayLabel = new JLabel("Decay:"); 
       JLabel decay2Label = new JLabel("Decay2:"); 
       JLabel crossLabel = new JLabel("Cross:"); 
       JLabel shiftLabel = new JLabel("i2 Shift:");
       JLabel onestepLabel = new JLabel("Hybrid:");
       JLabel lambda3Label = new JLabel("Turning:");
       
       lambda3Label.setHorizontalTextPosition(JLabel.LEFT); 
       smoothLabel.setHorizontalTextPosition(JLabel.LEFT); 
       decayLabel.setHorizontalTextPosition(JLabel.LEFT); 
       decay2Label.setHorizontalTextPosition(JLabel.LEFT); 
       crossLabel.setHorizontalTextPosition(JLabel.LEFT); 
       shiftLabel.setHorizontalTextPosition(JLabel.LEFT); 
       onestepLabel.setHorizontalTextPosition(JLabel.LEFT); 
       smoothLabel.setToolTipText("Adjust and apply smoothness regularization to coefficients");
       decayLabel.setToolTipText("Adjust and apply decay regularization to coefficients");
       crossLabel.setToolTipText("Adjust and apply cross-component regularization to coefficients");
       shiftLabel.setToolTipText("Adjust time-shift of i2 constraint. Works only when i2 constraint is on");

       JPanel smoothCon = new JPanel(); 
       JPanel decayCon = new JPanel(); 
       JPanel decay2Con = new JPanel(); 
       JPanel crossCon = new JPanel();  
       JPanel lagCon = new JPanel(); 
       JPanel shiftCon = new JPanel();
       JPanel onestepCon = new JPanel();
       JPanel lambda3Con = new JPanel();
       
       lagBar.setPreferredSize(new Dimension(80, 15));
       crossBar.setPreferredSize(new Dimension(80, 15));
       smoothBar.setPreferredSize(new Dimension(80, 15));
       decayBar.setPreferredSize(new Dimension(80, 15));
       decay2Bar.setPreferredSize(new Dimension(80, 15));
       lambdaBar.setPreferredSize(new Dimension(100, 15));
       expBar.setPreferredSize(new Dimension(100, 15));
       shiftBar.setPreferredSize(new Dimension(80,15));
       onestepBar.setPreferredSize(new Dimension(80,15));
       textLag.setPreferredSize(new Dimension(20,15));
       lambda3Bar.setPreferredSize(new Dimension(80,15));
       /*smoothCon.add(smoothLabel); smoothCon.add(smoothBar); smoothCon.add(smoothText); 
       decayCon.add(decayLabel);  decayCon.add(decayBar);  decayCon.add(decayText);
       crossCon.add(crossLabel);  crossCon.add(crossBar);  crossCon.add(crossText);
       lagCon.add(lagLabel); lagCon.add(lagBar); lagCon.add(textLag);
       smoothCon.setLayout(new GridLayout(1,3,3,0)); 
       decayCon.setLayout(new GridLayout(1,3,3,0)); 
       crossCon.setLayout(new GridLayout(1,3,3,0)); 
       lagCon.setLayout(new GridLayout(1,3,3,0));*/

        JLabel outofSampLabel = new JLabel("OutSampN:");
        JPanel outsampPanel = new JPanel();
        paramLayout = new GroupLayout(outsampPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()      
          .addComponent(outofSampLabel) .addComponent(outof_sample).addComponent(outof_sampleText));                                                                         
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup()
           .addComponent(outofSampLabel).addComponent(outof_sample).addComponent(outof_sampleText)));
        outsampPanel.setLayout(paramLayout);       
       
        paramLayout = new GroupLayout(lagCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
        .addComponent(lagLabel).addComponent(lagBar).addComponent(textLag));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(lagLabel).addComponent(lagBar).addComponent(textLag)));
        lagCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(onestepCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
        .addComponent(onestepLabel).addComponent(onestepBar).addComponent(onestepText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(onestepLabel).addComponent(onestepBar).addComponent(onestepText)));
        onestepCon.setLayout(paramLayout);
        
        
        paramLayout = new GroupLayout(smoothCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(smoothLabel).addComponent(smoothBar).addComponent(smoothText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(smoothLabel).addComponent(smoothBar).addComponent(smoothText)));
        smoothCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(decayCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(decayLabel).addComponent(decayBar).addComponent(decayText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(decayLabel).addComponent(decayBar).addComponent(decayText)));
        decayCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(decay2Con); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(decay2Label).addComponent(decay2Bar).addComponent(decay2Text));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(decay2Label).addComponent(decay2Bar).addComponent(decay2Text)));
        decay2Con.setLayout(paramLayout);

        paramLayout = new GroupLayout(crossCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(crossLabel).addComponent(crossBar).addComponent(crossText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(crossLabel).addComponent(crossBar).addComponent(crossText)));
        crossCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(shiftCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(shiftLabel).addComponent(shiftBar).addComponent(shiftText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(shiftLabel).addComponent(shiftBar).addComponent(shiftText)));
        shiftCon.setLayout(paramLayout);

        paramLayout = new GroupLayout(lambda3Con); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(lambda3Label).addComponent(lambda3Bar).addComponent(lambda3Text));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(lambda3Label).addComponent(lambda3Bar).addComponent(lambda3Text)));
        lambda3Con.setLayout(paramLayout);        
        

       JPanel regPanel = new JPanel(); 
       TitledBorder regBorder = new TitledBorder(new LineBorder(myBlue),"Regularization");

       JPanel reg2Panel = new JPanel(); 
       TitledBorder reg2Border = new TitledBorder(new LineBorder(myBlue),"Advanced Filter Ingredients");       
        
       regBorder.setTitleColor(Color.LIGHT_GRAY);
       regPanel.setBorder(regBorder);
       paramLayout = new GroupLayout(regPanel);
       paramLayout.setAutoCreateGaps(false); paramLayout.setAutoCreateContainerGaps(false);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(smoothCon).addComponent(decayCon).addComponent(decay2Con).addComponent(crossCon));
       paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(smoothCon).addComponent(decayCon).addComponent(decay2Con).addComponent(crossCon)));
       regPanel.setLayout(paramLayout);

       reg2Border.setTitleColor(Color.LIGHT_GRAY);
       reg2Panel.setBorder(reg2Border);
       paramLayout = new GroupLayout(reg2Panel);
       paramLayout.setAutoCreateGaps(false); paramLayout.setAutoCreateContainerGaps(false);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(lambda3Con).addComponent(lagCon).addComponent(shiftCon).addComponent(onestepCon));
       paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(lambda3Con).addComponent(lagCon).addComponent(shiftCon).addComponent(onestepCon)));
       reg2Panel.setLayout(paramLayout);       
       
       
      JPanel plotHybridPanel = new JPanel();
      
      paramLayout = new GroupLayout(plotHybridPanel);
      paramLayout.setAutoCreateGaps(false);
      paramLayout.setAutoCreateContainerGaps(false);        
      paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
             .addComponent(filterPlot[0]) 
             .addComponent(filterPlot[1])
             .addComponent(filterPlot[2]));

      paramLayout.setVerticalGroup(
      paramLayout.createSequentialGroup()
        .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(filterPlot[0]) 
             .addComponent(filterPlot[1])
             .addComponent(filterPlot[2])));
      plotHybridPanel.setLayout(paramLayout);

       /*------------------------------------------------------------------------------
          Panel for cointegration factors at freq 0
       -------------------------------------------------------------------------------*/
       JPanel cointPanel = new JPanel();
       BevelBorder cointBorder = new BevelBorder(BevelBorder.RAISED);
       cointPanel.setSize(new Dimension(700, 60));
       cointPanel.setBorder(cointBorder);

       cointPanel.setLayout(new GridLayout(1,10,2,0)); 
       cointPanel.add(cointText);
       for(i=0;i<10;i++) {cointPanel.add(cointSliders[i]);}
       
       /*
        paramLayout = new GroupLayout(cointPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(cointSliders[0]).addComponent(cointSliders[1]).addComponent(cointSliders[2])
          .addComponent(cointSliders[3]).addComponent(cointSliders[4]).addComponent(cointSliders[5]).addComponent(cointSliders[6])
          .addComponent(cointSliders[7]).addComponent(cointSliders[8]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(cointSliders[0]).addComponent(cointSliders[1]).addComponent(cointSliders[2])
          .addComponent(cointSliders[3]).addComponent(cointSliders[4]).addComponent(cointSliders[5]).addComponent(cointSliders[6])
          .addComponent(cointSliders[7]).addComponent(cointSliders[8])));    
        cointPanel.setLayout(paramLayout); 
        */




       new JLabel("(\u03C9_1, \u03C9_2):");
       new JLabel("(N_1 \u03C0/D_1), (N_2 \u03C0/D_2):");    
 
       turnPoint = new JLabel("     Turning-Point"); turnPoint.setHorizontalTextPosition(JLabel.RIGHT); 
       turnPoint.setToolTipText("Adjust filter constraints with i2 and lambda to produce better turning point estimation.");
       levelChange = new JLabel("     Level-Change "); levelChange.setHorizontalTextPosition(JLabel.RIGHT); 
       differ = new JLabel("Differencing"); differ.setHorizontalTextPosition(JLabel.RIGHT);  
       filterLength = new JLabel("Filter Lag/Length |"); filterLength.setHorizontalTextPosition(JLabel.RIGHT); 
       mdfaLabel = new JLabel("Min:");       

       bandLabel = new JLabel("Order "); bandLabel.setToolTipText("Adjust lengths of computed concurrent filter and/or symmetric filter.");
       new JLabel("SymOrder");
        
       JPanel compPan = new JPanel();
       paramPanel2 = new JPanel();

       //---------- Build borders around areas of selections --------------------------------------
       //TitledBorder timeSelBorder = new TitledBorder(new LineBorder(myBlue),"Time Domain Plots");
       //timeSelBorder.setTitleColor(Color.LIGHT_GRAY);

       dataMixPanel = new JPanel();
       BevelBorder dataMixBorder = new BevelBorder(BevelBorder.RAISED);
       JLabel dataMixLabel = new JLabel("Choose Explanatory Series:");
       dataMixPanel.setSize(new Dimension(800, 80));      
        paramLayout = new GroupLayout(dataMixPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(dataMixLabel).addComponent(expVariable[0]).addComponent(expVariable[1]).addComponent(expVariable[2])
          .addComponent(expVariable[3]).addComponent(expVariable[4]).addComponent(expVariable[5]).addComponent(expVariable[6])
          .addComponent(expVariable[7]).addComponent(expVariable[8]).addComponent(expVariable[9]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(dataMixLabel).addComponent(expVariable[0]).addComponent(expVariable[1]).addComponent(expVariable[2])
          .addComponent(expVariable[3]).addComponent(expVariable[4]).addComponent(expVariable[5]).addComponent(expVariable[6])
          .addComponent(expVariable[7]).addComponent(expVariable[8]).addComponent(expVariable[9])));
        dataMixPanel.setLayout(paramLayout); 
        dataMixPanel.setBorder(dataMixBorder);         
       


       BevelBorder perSelBorder = new BevelBorder(BevelBorder.RAISED);
       JLabel pdp = new JLabel("Periodograms I_N(\u03C9):  ");
       periodGrid = new JPanel(); periodGrid.setSize(new Dimension(700, 50));
       JPanel coeffGrid = new JPanel(); coeffGrid.setSize(new Dimension(700, 50));
       JLabel coeffLabel = new JLabel("Filter Coefficients:  ");

       JPanel accountGrid = new JPanel(); accountGrid.setSize(new Dimension(700, 50));
       JLabel accountLabel = new JLabel("Trading Plots: ");


       
        paramLayout = new GroupLayout(periodGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(pdp).addComponent(periodPlot[0]).addComponent(periodPlot[1]).addComponent(periodPlot[2]).addComponent(periodPlot[3])
          .addComponent(periodPlot[4]).addComponent(periodPlot[5]).addComponent(periodPlot[6])
          .addComponent(periodPlot[7]).addComponent(periodPlot[8]).addComponent(periodPlot[9]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(pdp).addComponent(periodPlot[0]).addComponent(periodPlot[1]).addComponent(periodPlot[2]).addComponent(periodPlot[3])
           .addComponent(periodPlot[4]).addComponent(periodPlot[5]).addComponent(periodPlot[6])   
           .addComponent(periodPlot[7]).addComponent(periodPlot[8]).addComponent(periodPlot[9]))); 
        periodGrid.setLayout(paramLayout); 
        periodGrid.setBorder(perSelBorder); 


     
        paramLayout = new GroupLayout(coeffGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(coeffLabel).addComponent(coeffPlot[0]).addComponent(coeffPlot[1]).addComponent(coeffPlot[2]).addComponent(coeffPlot[3])
          .addComponent(coeffPlot[4]).addComponent(coeffPlot[5]).addComponent(coeffPlot[6])
          .addComponent(coeffPlot[7]).addComponent(coeffPlot[8]).addComponent(coeffPlot[9]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(coeffLabel).addComponent(coeffPlot[0]).addComponent(coeffPlot[1]).addComponent(coeffPlot[2]).addComponent(coeffPlot[3])
           .addComponent(coeffPlot[4]).addComponent(coeffPlot[5]).addComponent(coeffPlot[6])
           .addComponent(coeffPlot[7]).addComponent(coeffPlot[8]).addComponent(coeffPlot[9])));    
        coeffGrid.setLayout(paramLayout); 
        coeffGrid.setBorder(perSelBorder); 


        paramLayout = new GroupLayout(accountGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(accountLabel).addComponent(accountPlot).addComponent(signalPlot).addComponent(logpricePlot).addComponent(logreturnPlot).addComponent(linesPlot).addComponent(buysellPlot).addComponent(filteredPrice));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(accountLabel).addComponent(accountPlot).addComponent(signalPlot).addComponent(logpricePlot).addComponent(logreturnPlot).addComponent(linesPlot).addComponent(buysellPlot).addComponent(filteredPrice)));    
        accountGrid.setLayout(paramLayout); 
        accountGrid.setBorder(perSelBorder); 




       BevelBorder sampSelBorder = new BevelBorder(BevelBorder.RAISED);    
       JPanel sampGrid = new JPanel(); sampGrid.setSize(new Dimension(200, 50));
        paramLayout = new GroupLayout(sampGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(sampLabel).addComponent(sampBar).addComponent(sampText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(sampLabel).addComponent(sampBar).addComponent(sampText)));    
        sampGrid.setLayout(paramLayout);   
        sampGrid.setBorder(sampSelBorder);  


        
       BevelBorder timeSelBorder = new BevelBorder(BevelBorder.RAISED);
       new JLabel("Time Domain Plots  ");
       timeGrid = new JPanel(); timeGrid.setSize(new Dimension(750, 50));      
        paramLayout = new GroupLayout(timeGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(plotHist).addComponent(plotXf).addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
          .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6])
          .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[11])
          );
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(plotHist).addComponent(plotXf).addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
           .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6]) 
           .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[11])
          ));    
        timeGrid.setLayout(paramLayout); 
        timeGrid.setBorder(timeSelBorder);  

       BevelBorder diagBorder = new BevelBorder(BevelBorder.RAISED);      
       JPanel diagGrid = new JPanel(); diagGrid.setSize(new Dimension(60, 50));      
        paramLayout = new GroupLayout(diagGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(mdfaText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(mdfaText)));    
        diagGrid.setLayout(paramLayout); 
        diagGrid.setBorder(diagBorder); 





      //---------------------------------------------------------------------------------

       //---------- Add the freq series domain checkboxes ----------
       TitledBorder freqSelBorder = new TitledBorder(new LineBorder(myBlue),"Filter Characteristic");
       freqSelBorder.setTitleColor(Color.LIGHT_GRAY);
       JLabel fdp = new JLabel("Plot");
        typefreqGrid = new JPanel(); 
        typefreqGrid.setPreferredSize(new Dimension(200, 60));
        paramLayout = new GroupLayout(typefreqGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(fdp).addComponent(gammaBut).addComponent(amplitudeBut).addComponent(timeDelayBut));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(fdp).addComponent(gammaBut).addComponent(amplitudeBut).addComponent(timeDelayBut)));    
        typefreqGrid.setLayout(paramLayout);   
        typefreqGrid.setBorder(timeSelBorder);  

       //---------- Add the freq series domain checkboxes ----------
       TitledBorder seriesSelBorder = new TitledBorder(new LineBorder(myBlue),"Filters");
       seriesSelBorder.setTitleColor(Color.LIGHT_GRAY);
     
        freqGrid = new JPanel(); freqGrid.setPreferredSize(new Dimension(500, 60));
        paramLayout = new GroupLayout(freqGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(plotGamma).addComponent(freqPlot[0]).addComponent(freqPlot[1]).addComponent(freqPlot[2])
          .addComponent(freqPlot[3]).addComponent(freqPlot[4]).addComponent(freqPlot[5])
          .addComponent(freqPlot[6]).addComponent(freqPlot[7]).addComponent(freqPlot[8]).addComponent(freqPlot[10])
          );
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(plotGamma).addComponent(freqPlot[0]).addComponent(freqPlot[1]).addComponent(freqPlot[2])
           .addComponent(freqPlot[3]).addComponent(freqPlot[4]).addComponent(freqPlot[5])   
           .addComponent(freqPlot[6]).addComponent(freqPlot[7]).addComponent(freqPlot[8]).addComponent(freqPlot[10])
           ));    
        freqGrid.setLayout(paramLayout);    
        freqGrid.setBorder(timeSelBorder);  

       
       /*---------- Add the simulation type checkboxes ----------
             N:  []   nrep: [] 
             sarimaCheck[]  simCheck[]  noneCheck[]
       ----------------------------------------------------------*/
       //GroupLayout paramLayout;
       TitledBorder simSelBorder = new TitledBorder(new LineBorder(myBlue),"Simulation Data Type");
       simSelBorder.setTitleColor(Color.LIGHT_GRAY);
       simulGrid = new JPanel();
       simulGrid.setBorder(simSelBorder); 

       //.addComponent(nrepLabel)  .addComponent(sarimaCheck) .addComponent(nrepBar) .addComponent(nrepText)

       /* paramLayout = new GroupLayout(simulGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup() 
             .addComponent(nobsLabel).addComponent(simCheck).addComponent(extFore))
          .addGroup(paramLayout.createParallelGroup() 
             .addComponent(nObsBar)))
          .addGroup(paramLayout.createParallelGroup() 
             .addComponent(nText).addComponent(noneCheck)));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(nobsLabel).addComponent(nObsBar).addComponent(nText))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(nrepLabel).addComponent(nrepBar).addComponent(nrepText))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(simCheck).addComponent(sarimaCheck).addComponent(noneCheck)) 
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(extFore)));
        simulGrid.setLayout(paramLayout);*/

        
  
     
        
       /*---------- Add the filter parameters -------------------

         Turning-Point:  lambda[]  i2[]   Filter Lag/Length  Lag[]  L[]
         Level-Change:   exp[]     i1[]   Cuttoff:  cutN[] cutD[] cutoff()
         Differencing:  d[]   D[]         Gamma: 
 
       ----------------------------------------------------------*/
      
       TitledBorder filterParamBorder = new TitledBorder(new LineBorder(myBlue),"iMDFA Filter Parameters");
       filterParamBorder.setTitleColor(Color.LIGHT_GRAY);
       paramPanel1 = new JPanel();
       //paramPanel1.setBorder(filterParamBorder);            

      
       tpPanel = new JPanel();
       paramLayout = new GroupLayout(tpPanel);
       paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(plotHybridPanel))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(iterCheck) .addComponent(reComp))
          .addGroup(paramLayout.createParallelGroup() 
                 .addComponent(turnPoint).addComponent(levelChange))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(i2Check).addComponent(i1Check))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(dCheck).addComponent(DCheck))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(lamLabel).addComponent(expLabel))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(lambdaBar).addComponent(expBar))
          .addGroup(paramLayout.createParallelGroup()
                 .addComponent(lambdaText).addComponent(expText)));
        paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(paramLayout.createSequentialGroup()
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(plotHybridPanel).addComponent(iterCheck) .addComponent(turnPoint) .addComponent(i2Check) .addComponent(dCheck) .addComponent(lamLabel) 
                    .addComponent(lambdaBar) .addComponent(lambdaText))
               .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                  .addComponent(reComp) .addComponent(levelChange) .addComponent(i1Check) .addComponent(DCheck) .addComponent(expLabel) 
                   .addComponent(expBar) .addComponent(expText) ) ) ));
        tpPanel.setLayout(paramLayout);

        resPanel = new JPanel();
        paramLayout = new GroupLayout(resPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(bandLabel)         
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(LLabel)
             .addComponent(l1Label))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(LBar)
             .addComponent(l1Bar))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(LText)
             .addComponent(l1Text)));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(bandLabel)
           .addComponent(LLabel)
           .addComponent(LBar)
           .addComponent(LText))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(l1Label)
           .addComponent(l1Bar)
           .addComponent(l1Text)));
        resPanel.setLayout(paramLayout);
                                 


       //
        //Box iterCheck2 = Box.createVerticalBox();
        //iterCheck2.add(iterCheck); iterCheck2.add(reComp);




        new JPanel();
          
          
          
          
          
          
        paramLayout = new GroupLayout(paramPanel1);
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()      
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
             .addComponent(tpPanel).addComponent(reg2Panel))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING) 
              .addComponent(resPanel).addComponent(outsampPanel)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup()
           .addComponent(tpPanel)
           .addComponent(resPanel))
         .addGroup(paramLayout.createParallelGroup()
          .addComponent(reg2Panel)
          .addComponent(outsampPanel)));
        paramPanel1.setLayout(paramLayout);
       
        paramLayout = new GroupLayout(paramPanel2);
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()           
          .addComponent(paramPanel1));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup()
           .addComponent(paramPanel1)));
       paramPanel2.setLayout(paramLayout);


      //------------------ Make Filter Design Panel ---------------------------------------------
      filterPanel = new JPanel(); new JPanel(); new JPanel();
      TitledBorder filterDesBorder = new TitledBorder(new LineBorder(myBlue),"MDFA Target Filter Selection");
      filterDesBorder.setTitleColor(Color.LIGHT_GRAY);
      filterPanel.setBorder(filterDesBorder);           

      JPanel filtType = new JPanel();  
      JPanel filtXPanel = new JPanel(); 
      JPanel x13filtPanel = new JPanel();
      JPanel rkhsfiltPanel = new JPanel();
       
 




      
      TitledBorder filterDesBorderX13 = new TitledBorder(new LineBorder(myBlue),"X13 Target Filter Selection");
      filterDesBorderX13.setTitleColor(Color.LIGHT_GRAY);
      x13filtPanel.setBorder(filterDesBorderX13);          

      TitledBorder filterDesBorderARMA = new TitledBorder(new LineBorder(myBlue),"RKHS Target Filter Selection");
      filterDesBorderARMA.setTitleColor(Color.LIGHT_GRAY);
      


        paramLayout = new GroupLayout(filtType);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(lowBut) 
             .addComponent(rampBut)
             .addComponent(bandBut)
             .addComponent(highBut)));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
             .addComponent(lowBut) 
             .addComponent(rampBut)
             .addComponent(bandBut)
             .addComponent(highBut));
        filtType.setLayout(paramLayout);

 
        paramLayout = new GroupLayout(filtXPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(omegaBut)
             .addComponent(fracBut)
             .addComponent(fixBandCheck))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)             
             .addComponent(omega0Text) 
             .addComponent(omega2Text)
             .addComponent(den1Text))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)           
             .addComponent(omega1Text) 
             .addComponent(omega3Text)
             .addComponent(den2Text))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)             
             .addComponent(omega0Bar)
             .addComponent(omega2Bar)
             .addComponent(den1Bar))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(omega1Bar)
             .addComponent(omega3Bar)
             .addComponent(den2Bar)));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(omegaBut) 
             .addComponent(omega0Text)
             .addComponent(omega1Text)            
             .addComponent(omega0Bar)
             .addComponent(omega1Bar))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(fracBut) 
             .addComponent(omega2Text)
             .addComponent(omega3Text) 
             .addComponent(omega2Bar)
             .addComponent(omega3Bar))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(fixBandCheck)
             .addComponent(den1Text)
             .addComponent(den2Text) 
             .addComponent(den1Bar)
             .addComponent(den2Bar)));
        filtXPanel.setLayout(paramLayout);


        //------ Set layout for x13filter choice -----
        paramLayout = new GroupLayout(x13filtPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addComponent(x13filter).addComponent(spec_densBox).addComponent(bLabel2).addComponent(x13trend).addComponent(x13seas).addComponent(x13ti));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(x13filter).addComponent(spec_densBox).addComponent(bLabel2).addComponent(x13trend).addComponent(x13seas).addComponent(x13ti)));
         x13filtPanel.setLayout(paramLayout);

         //------ Set layout for arma/rkhs choice -----
        paramLayout = new GroupLayout(rkhsfiltPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);      
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addComponent(rkhsBox).addComponent(rkhsCombo).addComponent(bLabel3).addComponent(armaSD).addComponent(pLabel).addComponent(pBox).addComponent(qLabel).addComponent(qBox));      
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(rkhsBox).addComponent(rkhsCombo).addComponent(bLabel3).addComponent(armaSD).addComponent(pLabel).addComponent(pBox).addComponent(qLabel).addComponent(qBox)));
         rkhsfiltPanel.setLayout(paramLayout);  

        rkhsfiltPanel.setBorder(filterDesBorderARMA); 

        Box filterPanel4 = Box.createHorizontalBox();
        filterPanel4.add(x13filtPanel); filterPanel4.add(rkhsfiltPanel);      


        //------ Set layout for automatic/compute buttons------
        paramLayout = new GroupLayout(compPan);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addComponent(automatic).addComponent(compute));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(automatic).addComponent(compute)));
        compPan.setLayout(paramLayout);


        //------ Set layout for symmetric filter choice-----
        paramLayout = new GroupLayout(filterPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(filtType)
             .addComponent(compPan))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(filtXPanel)
             .addComponent(filterPanel4)));         
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(filtType).addComponent(filtXPanel)) 
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(compPan).addComponent(filterPanel4)));
        filterPanel.setLayout(paramLayout);

       /* //------ Set layout for x13 and compute filter choice-----
        paramLayout = new GroupLayout(filterPanel2);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
            .addComponent(x13filtPanel).addComponent(compPan));
         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)              
             .addComponent(x13filtPanel).addComponent(compPan)));
        filterPanel2.setLayout(paramLayout);

        //------ Set final layout --------------------------------
        paramLayout = new GroupLayout(filterPanel1);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(filterPanel)
             .addComponent(filterPanel2)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()            
             .addComponent(filterPanel).addComponent(filterPanel2));
        filterPanel1.setLayout(paramLayout);*/




      //------------------------------------------------------------------------------------------  
     
      paramPanel1.setPreferredSize(new Dimension(450, 100));
      filterPanel.setPreferredSize(new Dimension(450, 100)); 

      Box timePane = Box.createVerticalBox();
      Box timePane2 = Box.createHorizontalBox();
      Box freqPane = Box.createVerticalBox();   
      Box freqPane2 = Box.createHorizontalBox();
      Box periodPane2 = Box.createHorizontalBox();  
      Box periodPane = Box.createVerticalBox();  
      Box coeffPane = Box.createVerticalBox();
      Box accountPane = Box.createVerticalBox();

      timePane2.add(timeGrid);
      timePane2.add(diagGrid); 

        JPanel allParam = new JPanel();
        paramLayout = new GroupLayout(allParam);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(paramPanel1)
             .addComponent(regPanel)));
        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()            
             .addComponent(paramPanel1).addComponent(regPanel));
        allParam.setLayout(paramLayout);


      //timePane.add(mdfa_canvas,BorderLayout.NORTH);
      timePane.add(timeScrollPane,BorderLayout.NORTH);
      timePane.add(timePane2,BorderLayout.SOUTH);

      freqPane2.add(freqGrid);
      freqPane2.add(typefreqGrid); 

      periodPane2.add(periodGrid);
      periodPane2.add(sampGrid);

      freqPane.add(filter_canvas,BorderLayout.NORTH);
      freqPane.add(freqPane2,BorderLayout.SOUTH);      
 
      periodPane.add(period_canvas,BorderLayout.NORTH);
      periodPane.add(periodPane2,BorderLayout.SOUTH);

      coeffPane.add(coef_canvas,BorderLayout.NORTH);
      coeffPane.add(coeffGrid,BorderLayout.SOUTH);

      accountPane.add(getAccount_canvas(), BorderLayout.NORTH); 
      accountPane.add(accountGrid, BorderLayout.SOUTH);

      //--------------------------------------------------

      //--------- Design for layout ---------------
      
      mdfaPlotPane.addTab("Time Domain", timePane);
      mdfaPlotPane.addTab("Filters", freqPane);
      mdfaPlotPane.addTab("Periodograms", periodPane);
      mdfaPlotPane.addTab("Coefficients", coeffPane); 
      mdfaPlotPane.addTab("Time-Frequency Domain", tfplotPanel);
      mdfaPlotPane.addTab("Spectral Density Canvas", specDensCanvas);
      mdfaPlotPane.addTab("Financial Trading Canvas", accountPane);
   



         //-------------Tabbed Pane for Parametric/NonParametric -----
          JTabbedPane filtTabbedPane = new JTabbedPane(JTabbedPane.TOP);
          /*filtTabbedPane.addTab("Real-Time Filter Design", paramPanel1,"Set parameters and constraints of real-time filter to approximated the target filter for given series.");
          filtTabbedPane.addTab("Target Filter Design", filterPanel,"Set the type of filter and frequency bands for the target filter of interest.");*/

          filtTabbedPane.addTab("Real-Time Filter Design", allParam);
          filtTabbedPane.addTab("Target Filter Design", filterPanel);
          filtTabbedPane.addTab("ZPC Filter Design",zpcPanel);
          filtTabbedPane.addTab("ART Filter Design",ARTpanel);
          filtTabbedPane.addTab("Time-Frequency",sigex);
          filtTabbedPane.addTab("Spectral Density Design", spectralDensPanel);
          
          filtTabbedPane.addChangeListener(new ChangeListener() {
              public void stateChanged(ChangeEvent evt) 
              { 
               
               JTabbedPane pane = (JTabbedPane)evt.getSource();
               int sel = pane.getSelectedIndex(); 
               if(sel == 2) 
               { 
                  computeZPCFilter(); 
               }
              }                    
           });
         filtTabbedPane.setSelectedIndex(0);


      Box controlsPane = Box.createHorizontalBox();
      //controlsPane.add(simulGrid);
      controlsPane.add(filtTabbedPane);
      

      this.add(mdfaPlotPane);
      //this.add(paramPanel1);
      this.add(controlsPane);
      this.setLayout(new GridLayout(2, 1));

        paramLayout = new GroupLayout(this);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()           
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(mdfaPlotPane)
             .addComponent(controlsPane)));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
           .addComponent(mdfaPlotPane)
           .addComponent(controlsPane));
        this.setLayout(paramLayout);

     

     /*
     System.out.println("last Group");
      GroupLayout epLayout = new GroupLayout(this);
       epLayout.setAutoCreateGaps(true);
       epLayout.setAutoCreateContainerGaps(true);
       epLayout.setHorizontalGroup(epLayout.createSequentialGroup()
          .addGroup(epLayout.createParallelGroup() 
           .addComponent(mdfaPlotPane)
              .addGroup(epLayout.createSequentialGroup()
                .addComponent(paramPanel1)
                .addComponent(simulGrid))));                
       epLayout.setVerticalGroup(epLayout.createSequentialGroup()
          .addComponent(mdfaPlotPane)
            .addGroup(epLayout.createParallelGroup()
              .addComponent(paramPanel1)
              .addComponent(simulGrid)));   
        
       
       this.setLayout(epLayout);
      */

      //System.out.println("end");

    }

    double[] differenceData(double[] data, int n, int nr, int d, int D)
    {
         int i; int n_obsd;
         double[] d_series; double[] temp;
         //set differenced data
         if((d == 1) && (D==1))
         {
           n_obsd = n - S - 1;  d_series = new double[n_obsd]; temp = new double[n-1];
           for(i=1; i < n; i++) {temp[i-1] = data[i] - data[i-1];}         
           for(i=S; i < n; i++) {d_series[i-S] = temp[i] - temp[i-S];}  
         }
         else if(d == 1)
         {
           n_obsd = n - 1; d_series = new double[n_obsd];
           for(i=1; i < n; i++) {d_series[i-1] = data[i] - data[i-1];}
         }
         else if(D==1)
         {
           n_obsd = n - S; d_series = new double[n_obsd];
           for(i=1; i < n; i++) {d_series[i-S] = data[i] - data[i-S];}
         }        
         else {return data;}

         return d_series;
    }

    double[] cumsum(double[] data, int n)
    {
       double[] cs = new double[n]; double sum; int k;
       sum = data[0]; cs[0] = data[0];
       for(k=1;k<n;k++)
       {
        sum = sum+data[k]; cs[k] = sum;
       }
       return cs;  
    } 

    public void updateTexts()
    {nText.setText(""+ n_obs); nrepText.setText(""+ n_rep);}


    public void setFilterEnabled(boolean c)
    {
       nrepBar.setEnabled(c); LBar.setEnabled(c);      
       l1Bar.setEnabled(c); lambdaBar.setEnabled(c);
       expBar.setEnabled(c); smoothBar.setEnabled(c);
       i1Check.setEnabled(c);  decayBar.setEnabled(c); decay2Bar.setEnabled(c);
       i2Check.setEnabled(c);  dCheck.setEnabled(c);  
       DCheck.setEnabled(c); crossBar.setEnabled(c);
       lagBar.setEnabled(c); compute.setEnabled(c);
       injectButton.setEnabled(c); ARTCanvas.setEnabled(c); prefilterButton.setEnabled(c);

       omega0Bar.setEnabled(c);  num2Bar.setEnabled(c); 
       den2Bar.setEnabled(c);    automatic.setEnabled(c); 
       omega1Bar.setEnabled(c);  optimizeButton.setEnabled(c); 
       injectButton.setEnabled(c); prefilterButton.setEnabled(c); 
       num1Bar.setEnabled(c);      compARTFilter.setEnabled(c);
       computeTimeFreq.setEnabled(c); optimizeARTFilter.setEnabled(c);
       den1Bar.setEnabled(c);

       for(int i=0;i<series_max;i++) {expVariable[i].setEnabled(c);}
  
    }


      class MyActionListener implements ActionListener {
        public void actionPerformed(ActionEvent e)
        {
          if("compute".equals(e.getActionCommand()))
          {
             computeFilter = true;// System.out.println("Computing filter");
             mdfa.computeFilterGeneral(computeFilter,false); 
             updatePlots(false,true);
             compute.setEnabled(false); 
          }
        }
      }


      class MyItemListener2 implements ItemListener {
        public void itemStateChanged(ItemEvent e)
        {
         updateData();
        }
      }

      class MyItemListener implements ItemListener {
        public void itemStateChanged(ItemEvent e)
        {
         int i; boolean sel; //computeFilter = true;
         Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}

         if(source == periodPlot[0]) {updatePeriod(0,sel);}
         else if(source == periodPlot[1]) {updatePeriod(1,sel);}
         for(i=0;i<n_rep;i++)
         {          
          if(source == timePlot[i]) {updateTime(i+1,sel);}
          else if(source == periodPlot[i+2]) {updatePeriod(i+2,sel);}         
         }
         for(i=0;i<=9;i++)
         {
          if(source == freqPlot[i])  {updateFreq(i+1,sel);} 
          if(source == coeffPlot[i]) {updateCoeffs(i,sel);}
         }
         if(source == freqPlot[10]) {updateFreq(11,sel);}
         /*if(source == extFore) 
         {
           ext_forecast = sel;  mdfa_canvas.setForeExt(fore_extend, ext_forecast);
         } */
         if(source == timePlot[11]) {if(simulate || b_outof_sample){updateTime(12,sel);}}
         if(source == plotXf) {updateTime(0,sel);}
         else if(source == plotGamma) {updateFreq(0,sel);}
         else if(source == amplitudeBut) {plotAmplitude();}
         else if(source == gammaBut) {plotGamma();}
         else if(source == timeDelayBut) {plotTimeDelay();}
         else if(source == i1Check) 
         {
           if(sel){i1 = 1;} else{i1=0;} 
           setBConstraints(i1, i2);
           activateWConst(i1,true);
         }
         else if(source == i2Check) {if(sel){i2 = 1;} else{i2=0;} setBConstraints(i1, i2);}
         else if(source == dCheck) {if(sel){dd = 1;} else{dd = 0;} setDDNoCompute(dd,DD);} 
         else if(source == DCheck) {if(sel){DD = 1;} else{DD = 0;} setDDNoCompute(dd,DD);}
         else if(source == simCheck) 
         {
           simulate = sel;  
           nObsBar.setEnabled(simulate); nrepBar.setEnabled(simulate);
           simulateIMDFA(sel); timePlot[10].setEnabled(simulate); setEnableX13(sel);  
         } 
         else if(source == noneCheck) 
         { 
           simulate = false;  useSARIMA = false; timePlot[10].setEnabled(false); mdfa_canvas.setPlots(7, false);
           nObsBar.setEnabled(false); nrepBar.setEnabled(false);  setEnableX13(false);  
         }
         else if(source == reComp)
         {
           reCompFilter = sel;
           if(reCompFilter)
           {
             setFilterEnabled(true); simulateIMDFA(true); addOutofSample();
           }
           else
           {
             mdfa.saveBfilter(); 
             setFilterEnabled(false);            
           }       
         }
         else if(source == iterCheck)
         {set_iter(sel);} 

         else if(source == lowBut)
         {
           omega0Bar.setValue(0); w0=0.0; omega0Bar.setEnabled(false); 
           omega1Bar.setEnabled(true); 
           num1Bar.setEnabled(false); 
           den1Bar.setEnabled(false);
           
           tfilter.setFilterType(0);  reComputeSymetricSignal(); low=true;high=false;ramp=false;band=false; setFilter(); 
           if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
           if(trading_func == 1)
           {getAccount_canvas().shiftPrice(false);}
         }
         else if(source == rampBut)
         {
           omega0Bar.setEnabled(true); 
           omega1Bar.setEnabled(true); 
           num1Bar.setEnabled(true); 
           den1Bar.setEnabled(true);

          tfilter.setFilterType(2);  reComputeSymetricSignal(); low=false;high=false;ramp=true;band=false; setFilter(); 
          if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
         } 
         else if(source == bandBut)
         {
           omega0Bar.setEnabled(true); 
           omega1Bar.setEnabled(true); 
           num1Bar.setEnabled(true); 
           den1Bar.setEnabled(true);
           tfilter.setFilterType(3);  reComputeSymetricSignal(); low=false;high=false;ramp=false;band=true; setFilter(); 
           if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
           //signal to accountCanvas to 
           if(trading_func == 1)
           {getAccount_canvas().shiftPrice(true);}
         }
         else if(source == highBut)
         {
      
           omega1Bar.setEnabled(true); 
           omega0Bar.setEnabled(true); 
           omega2Bar.setEnabled(true); 
           omega3Bar.setEnabled(true);

           tfilter.setFilterType(1);  reComputeSymetricSignal(); low=false;high=true;ramp=false;band=false; setFilter(); 
          if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
         }    
         else if(source == omegaBut)
         {

                  num1Bar.setEnabled(false); num2Bar.setEnabled(false);
                  den1Bar.setEnabled(false); den2Bar.setEnabled(false);
                  omega0Bar.setEnabled(true); omega1Bar.setEnabled(true);
                              
                  w0 = omega0Bar.getValue();  w1 = omega1Bar.getValue();
                  omega0Text.setText(""+df.format(w0)); omega1Text.setText(""+df.format(w1));
                  tfilter.setBand(w0, w1); setFilter(); 
                  if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
         }
         else if(source == fracBut)
         {

                  num1Bar.setEnabled(true); num2Bar.setEnabled(true);
                  den1Bar.setEnabled(true); den2Bar.setEnabled(true);
                  omega0Bar.setEnabled(false); omega1Bar.setEnabled(false);
                 
                 num1 = num1Bar.getValue(); den1 = den1Bar.getValue();
                 num2 = num2Bar.getValue(); den2 = den2Bar.getValue();
                 w0 = 1.0*num1*Math.PI/(1.0*den1); w1 = 1.0*num2*Math.PI/(1.0*den2);
                 tfilter.setBand(w0, w1); setFilter();
                 num1Text.setText(""+num1); num2Text.setText(""+num2); 
                 den1Text.setText(""+den1); den2Text.setText(""+den2);
                 if(useX13filters) {useX13filters = false; x13filter.setSelected(false);} 
                
          }
          else if(source == fixBandCheck)
          {fixBand = sel;}
          else if(source == automatic)
          {autoComp = sel; contARTCheck.setSelected(sel);}
          else if(source == contARTCheck) {autoComp = sel; automatic.setSelected(sel);}
          else if(source == x13filter)
          { 
             //expBar.setEnabled(sel);
            if(sel)
            {useX13filters = true; computeX13TargetFilter(); setFilter(); if(!autoComp){compute.setEnabled(true);}}
            else 
            {useX13filters = false;}                     
          }
          else if(source == x13trend) {x13model = 0; computeX13TargetFilter(); setFilter(); if(!autoComp){compute.setEnabled(true);}}
          else if(source == x13seas)  {x13model = 1; computeX13TargetFilter(); setFilter(); if(!autoComp){compute.setEnabled(true);}}
          else if(source == x13ti)    {x13model = 2; computeX13TargetFilter(); setFilter(); if(!autoComp){compute.setEnabled(true);}}
          else if(source == plotHist) {mdfa_canvas.plotHistorical(sel);}
          else if(source == spec_densBox) 
          {           
               computeX13TargetFilter();
               mdfa.useSpectralDensity(sel);
               mdfa.computeFilterGeneral(reCompFilter,false);
          }
          else if(source == rkhsBox)
          {
           if(sel)
           {
            if(useX13filters) {useX13filters = false; x13filter.setSelected(false);}
            lowBut.setSelected(false); rampBut.setSelected(false);
            bandBut.setSelected(false); highBut.setSelected(false);
           }    
           else
           {
             lowBut.setSelected(true); 
           }

          }
          else if(source == armaSD)
          {
              
          }
          else if(source == accountPlot)
          {getAccount_canvas().plot_account = sel; getAccount_canvas().go();}          
          else if(source == logpricePlot)
          {getAccount_canvas().plot_logprice = sel; getAccount_canvas().go();}
          else if(source == logreturnPlot)
          {getAccount_canvas().plot_logreturn = sel; getAccount_canvas().go();}      
          else if(source == signalPlot)
          {getAccount_canvas().plot_signal = sel; getAccount_canvas().go();}
          else if(source == linesPlot)
          {getAccount_canvas().plot_lines = sel; getAccount_canvas().go();}
          else if(source == buysellPlot)
          {getAccount_canvas().plot_buy = sel; getAccount_canvas().go();}     
          else if(source == filteredPrice)
          {getAccount_canvas().plot_indicator = sel; getAccount_canvas().go();}
       }
      }  
          

  public int min(int u, int v)
  {if(u < v) return u; else return v;}

 



  /*-----------------------------------------------------------------------------------------
         Stuff concerning the time-frequency and ART canvas

     artCanvas -- canvas to draw mouse controlled panel for accuracy reliability and timeliness
     timefreqpanel -- panel for controling ART and time-frequency map 
     timefreqCanvas -- panel and canvas for time frequency map and slider for selecting frequencies
     IMF plot       -- small panel for plotting IMFs and trend/cycle    

  ------------------------------------------------------------------------------------------*/

  /*---------------------------------------------------

     x,y \in [0,1]x[0,1] computes the corresonding lambda and exp 
     maps to ART   A = 3(1 - beta_r - beta_t), R = beta_r, T = beta_t         

  -----------------------------------------------------*/
  public void setART(double w, double h,boolean comp)
  {
     
      if(w > -100.0) {lambda = w*12.0;} else {w = lambda/12.0;}
      if(h > -100.0) {expweight = h*10.0;}    else {h = expweight/10.0;}

      ARTfilter[0] = ((1.0 - w - h) + 1.0)/2.0;   // --- accuracy
      ARTfilter[1] = h;                           // --- reliability
      ARTfilter[2] = w;                           // --- timeliness

      computeFilter = false;
      setLambda(lambda); //lambdaBar.setValue((int)(10*lambda));
      setExp(-expweight); //expBar.setValue((int)(10*expweight));
      computeFilter = true;
      if(comp){computeFilterNow();}  
      accText.setText(df.format(ARTfilter[0])); reliableText.setText(df.format(ARTfilter[1])); timeText.setText(df.format(ARTfilter[2]));
      //---- TO DO---------
      //if(autoComp){mseText.setText(df.format(mean_sqr_error));  
      aerrorText.setText(df.format(Math.log(mdfa.diff_band))); 
      rerrorText.setText(df.format(Math.log(mdfa.diff_stop)));
      mseText.setText(df.format(mdfa.MDFAmin));
      phasedelayText.setText(df.format(mdfa.tdelay));
      
  }


  
 

  public void computeFreqDivision()
  {
    double delta; int i;   
    freq_ints = new double[n_div+1];

    //------------cut the domain [0,Pi] into divisions 
  delta = (Math.PI-freq_start)/n_div;  
  freq_ints[0] = freq_start;
  freq_ints[1] = 0.0;
   
  for(i=2;i<=n_div;i++)
  {freq_ints[i] = freq_ints[i-1] + delta; delta = delta+.02;}
  
  for(i=1;i<=n_div;i++)
  {freq_ints[i] = freq_ints[i]*(Math.PI-freq_start)/freq_ints[n_div]+freq_start;}
  freq_ints[n_div] = freq_ints[n_div]-.01;
 
  }

 


    





















  

  /*---------------------------------------------------------------------

      Canvases and Panels galore 

    -----------------------------------------------------------------------*/

  public void setPlotTF(int c)
  {plotTF = c;}



  class artCanvas extends JPanel 
  {
         /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//-----------stuff for art canvas ---------------
     Graphics2D g2D;
     public int art_width, art_height;
     public double x = 5, y = 5, w = 10, h = 10;
     public int x1, y1, x2, y2;

     public Ellipse2D ellipse;
     public Ellipse2D selectedShape;
     public Cursor curCursor;
     public Color fillColor;

     public boolean mouseDown = false; 
     public boolean contComp = false;
 
     public artCanvas() 
     {
       
       setBackground(Color.black);
       addMouseListener(new MyArtMouseListener());
       addMouseMotionListener(new MyArtMouseMotionListener());
       
       fillColor = new Color(0, 128, 128);
     }

     public void setContComp(boolean t) {contComp = t;}
     public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  
       g2D = (Graphics2D)g;
       super.paintComponent(g);
       setBackground(Color.black);
       //g2D.setPaint(Color.white);
       //ellipse = new Ellipse2D.Double(x, y, w, h);
       //g2D.draw(ellipse);
       Dimension D =  this.getSize(); 
       this.art_width = D.width; this.art_height = D.height;
       

       //if(mouseDown)
       //{
        //g2D.setPaint(fillColor); // a dull blue-green
        //g2D.fill(ellipse);
       //}
       if (curCursor != null) setCursor(curCursor);
     }

     class MyArtMouseListener extends MouseAdapter 
     {
       
       public void mousePressed(MouseEvent e) 
       {        

           x = e.getX(); y = e.getY();
           setART((double)(x/art_width), (double)(y/art_height),false); 
           go();
       }

       public void mouseReleased(MouseEvent e) 
       {
         //mouseDown = false; 
         x = e.getX(); y = e.getY();
         go();                 
         setART((double)(x/art_width), (double)(y/art_height),true);  
                      
       }

       public void mouseClicked(MouseEvent e) 
       { }
        
     }
      

     class MyArtMouseMotionListener extends MouseMotionAdapter 
     {
      
      public void mouseDragged(MouseEvent e) 
      {
       boolean passx = false; boolean passy = false;
    
        if((e.getX() >= 1) && (e.getX() < art_width - 1)) passx = true;
        if((e.getY() >= 1) && (e.getY() < art_height - 1)) passy = true;
        
         if(passx && passy)
         {
          mouseDown = true; go();    
          x = e.getX(); y = e.getY();
          go();                 
          setART((double)(x/art_width), (double)(y/art_height),false);        
         }
       
        
      }

      public void mouseMoved(MouseEvent e) {

            curCursor = Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
            go();
      }
     }
  }
  

  class TimeFreqPlotPanel extends JPanel 
  {

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	// Variables declaration ---------------------------------------------------
    private JRadioButton amButton,fmButton,ifreqButton,phaseButton,saharButton,padideButton;
    public JSlider freqPickSlider;
    private JPanel plotChoicePanel;
    private JLabel plotLabel,timeLabel;
    private ButtonGroup selectplotGroup;
    private ButtonGroup colorGroup;
    //---- Custom panels -----------------
    private ScalePanel heatScalePanel;
    private FrequencyScalePanel freqScalePanel;
    private TimeScalePanel timeScalePanel;
    private TimeFrequencyMap timefreqMapCanvas;

    public TimeFreqPlotPanel() { initComponents();}

    private void initComponents() 
    {
       
        freqScalePanel = new FrequencyScalePanel();
        timeScalePanel = new TimeScalePanel();
        timefreqMapCanvas = new TimeFrequencyMap();        
        heatScalePanel = new ScalePanel();

        timeLabel = new JLabel(); plotLabel = new JLabel();
        freqPickSlider = new JSlider(); plotChoicePanel = new JPanel();
        ifreqButton = new JRadioButton(); amButton = new JRadioButton(); 
        phaseButton = new JRadioButton(); fmButton = new JRadioButton();
        
        saharButton = new JRadioButton();
        padideButton = new JRadioButton(); padideButton.setSelected(true);
        JLabel colorLabel = new JLabel("Colors"); colorLabel.setFont(new Font("Arial",0,10));


        timeLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        timeLabel.setText("Time");

        //freqScalePanel.setBackground(new Color(20, 20, 20));
        freqScalePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        freqScalePanel.setToolTipText("Frequency bar");

        GroupLayout freqScalePanelLayout = new GroupLayout(freqScalePanel);
        freqScalePanel.setLayout(freqScalePanelLayout);
        freqScalePanelLayout.setHorizontalGroup(
            freqScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 8, Short.MAX_VALUE)
        );
        freqScalePanelLayout.setVerticalGroup(
            freqScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        
        timeScalePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        timeScalePanel.setToolTipText("Time scale");

        GroupLayout timeScalePanelLayout = new GroupLayout(timeScalePanel);
        timeScalePanel.setLayout(timeScalePanelLayout); timeScalePanel.setL(mdfa.L);
        timeScalePanelLayout.setHorizontalGroup(
            timeScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 775, Short.MAX_VALUE)
        );
        timeScalePanelLayout.setVerticalGroup(
            timeScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 11, Short.MAX_VALUE)
        );

        
        timefreqMapCanvas.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        timefreqMapCanvas.setBackground(new Color(0, 0, 0));
        GroupLayout timefreqMapCanvasLayout = new GroupLayout(timefreqMapCanvas);
        timefreqMapCanvas.setLayout(timefreqMapCanvasLayout);
        timefreqMapCanvasLayout.setHorizontalGroup(
            timefreqMapCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 775, Short.MAX_VALUE)
        );
        timefreqMapCanvasLayout.setVerticalGroup(
            timefreqMapCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        freqPickSlider.setFont(new Font("Arial", 0, 3)); // NOI18N
        freqPickSlider.setForeground(new Color(41, 41, 41));
        freqPickSlider.setMaximum(150);
        freqPickSlider.setOrientation(JSlider.VERTICAL);
        freqPickSlider.setToolTipText("Choose frequency to plot the corresponding Instrinsic mode below");
        freqPickSlider.setValue(0);
        freqPickSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                freqPickSliderStateChanged(evt);
            }
        });

        heatScalePanel.setBackground(new Color(0, 0, 0));
        heatScalePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        heatScalePanel.setToolTipText("Heat map scale of time-frequency values");

        GroupLayout heatScalePanelLayout = new GroupLayout(heatScalePanel);
        heatScalePanel.setLayout(heatScalePanelLayout);
        heatScalePanelLayout.setHorizontalGroup(
            heatScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 56, Short.MAX_VALUE)
        );
        heatScalePanelLayout.setVerticalGroup(
            heatScalePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        plotChoicePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        selectplotGroup = new ButtonGroup(); colorGroup = new ButtonGroup();
        ifreqButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        ifreqButton.setText("Phase");
        ifreqButton.setToolTipText("Plot the time-frequency map of the instantaneous frequency");
        ifreqButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                ifreqButtonActionPerformed(evt);
            }
        }); selectplotGroup.add(ifreqButton);

        amButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        amButton.setText("AM"); amButton.setSelected(true);
        amButton.setToolTipText("Plot the Amplitude-Modulated components");
        amButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                amButtonActionPerformed(evt);
            }
        }); selectplotGroup.add(amButton);

        phaseButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        phaseButton.setText("iFreq");
        phaseButton.setToolTipText("Plot the map of each instrinsic mode function instantaneous frequency");
        phaseButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                phaseButtonActionPerformed(evt);
            }
        }); selectplotGroup.add(phaseButton);

        fmButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        fmButton.setText("FM");
        fmButton.setToolTipText("Plot the FM Modulated Components");
        fmButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                fmButtonActionPerformed(evt);
            }
        }); selectplotGroup.add(fmButton);

        saharButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        saharButton.setText("Sahar");
        saharButton.setToolTipText("Use the Sahar color scheme, Sahar meaning sunrise in Persian");
        saharButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                saharButtonActionPerformed(evt);
            }
        }); colorGroup.add(saharButton);

        padideButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        padideButton.setText("Padideh");
        padideButton.setToolTipText("Use the Padideh color scheme, Padideh meaning phenomenon in Persian");
        padideButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                padidehButtonActionPerformed(evt);
            }
        }); colorGroup.add(padideButton);



        plotLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        plotLabel.setText("Plot Map"); 
        plotLabel.setForeground(new Color(67, 144, 239));
        colorLabel.setForeground(new Color(67, 144, 239));
 

        /*GroupLayout plotChoicePanelLayout = new GroupLayout(plotChoicePanel);
        plotChoicePanel.setLayout(plotChoicePanelLayout);
        plotChoicePanelLayout.setHorizontalGroup(
            plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotChoicePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(plotChoicePanelLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addGroup(plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(plotLabel)
                            .addComponent(phaseButton))
                        .addGap(21, 21, 21))
                    .addGroup(plotChoicePanelLayout.createSequentialGroup()
                        .addGroup(plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(fmButton)
                            .addComponent(amButton)
                            .addComponent(ifreqButton)
                            .addComponent(colorLabel)
                            .addComponent(saharButton)
                            .addComponent(padideButton)
                        .addGap(0, 0, Short.MAX_VALUE))))
        ));
        plotChoicePanelLayout.setVerticalGroup(
            plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotChoicePanelLayout.createSequentialGroup()
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(plotLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(amButton)
                .addGap(4, 4, 4)
                .addComponent(fmButton, GroupLayout.PREFERRED_SIZE, 20, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(phaseButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ifreqButton)
                .addGap(91, 91, 91)
                .addComponent(colorLabel)
                .addComponent(saharButton)
                .addComponent(padideButton))
        );*/

        GroupLayout plotChoicePanelLayout = new GroupLayout(plotChoicePanel);
        plotChoicePanel.setLayout(plotChoicePanelLayout);
        plotChoicePanelLayout.setHorizontalGroup(
            plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotChoicePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(plotLabel)
                    .addComponent(amButton)
                    .addComponent(fmButton)
                    .addComponent(phaseButton)
                    .addComponent(ifreqButton)
                    .addComponent(colorLabel)
                    .addComponent(saharButton)
                    .addComponent(padideButton))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        plotChoicePanelLayout.setVerticalGroup(
            plotChoicePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotChoicePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(plotLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(amButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(fmButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(phaseButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ifreqButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(colorLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(saharButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(padideButton)
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        GroupLayout layout = new GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(timeScalePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(20, 20, 20)
                        .addComponent(heatScalePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 51, Short.MAX_VALUE)
                        .addComponent(freqScalePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(timefreqMapCanvas, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(freqPickSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(3, 3, 3)
                .addComponent(plotChoicePanel, GroupLayout.PREFERRED_SIZE, 95, GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
            .addGroup(GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addComponent(timeLabel)
                .addGap(523, 523, 523))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(plotChoicePanel, GroupLayout.PREFERRED_SIZE, 200, GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(freqPickSlider, GroupLayout.DEFAULT_SIZE, 290, Short.MAX_VALUE)
                            .addComponent(freqScalePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(heatScalePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(timefreqMapCanvas, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(timeScalePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(timeLabel)))
                .addContainerGap())
        );
    }

   //---------- The listeners --------------
    private void ifreqButtonActionPerformed(ActionEvent evt)   {setTimeFreqData(mdfa.ifreqMap, mdfa.flength, n_imfs, 3); plotTF = 3;}
    private void amButtonActionPerformed(ActionEvent evt)      {setTimeFreqData(mdfa.amMap, mdfa.flength, n_imfs, 0); plotTF = 0;}
    private void phaseButtonActionPerformed(ActionEvent evt)   {setTimeFreqData(mdfa.phaseMap, mdfa.flength, n_imfs, 2); plotTF = 2;}
    private void fmButtonActionPerformed(ActionEvent evt)      {setTimeFreqData(mdfa.fmMap, mdfa.flength, n_imfs, 1); plotTF = 1;}
    private void saharButtonActionPerformed(ActionEvent evt)   
    {timefreqMapCanvas.changeColor(1);  heatScalePanel.updateColor(timefreqMapCanvas.colorArray);}                             
    private void padidehButtonActionPerformed(ActionEvent evt)  
    {timefreqMapCanvas.changeColor(0); heatScalePanel.updateColor(timefreqMapCanvas.colorArray);} 
    private void freqPickSliderStateChanged(ChangeEvent evt) {
        JSlider s = (JSlider)evt.getSource();
        changeIMFfreq(s.getValue());        
    }

    public void changeIMFfreq(int val)
    {
       timefreqMapCanvas.setLine(val,K+1);
       tcimfPanel.setIMF(val,K+1,freq_ints); 
    }

    public void setTimeFreqData(double[][] data, int N, int n_imfs, int _sel)
    {
       timefreqMapCanvas.updateFM(data, N, n_imfs,  _sel);
       heatScalePanel.updateColor(timefreqMapCanvas.colorArray);
       timefreqMapCanvas.setL(mdfa.L); tcimfPanel.setL(mdfa.L);
       heatScalePanel.setExtremum(timefreqMapCanvas.dataRangeMin,timefreqMapCanvas.dataRangeMax);
    }

    public void setNobs() 
    {
        timeScalePanel.setNObs(mdfa.flength); 
        timeScalePanel.setL(mdfa.L); 
        freqPickSlider.setMaximum(mdfa.K+1);
        timefreqMapCanvas.setL(mdfa.L);
        tcimfPanel.setL(mdfa.L);

    }
    public void setNDiv(double[] f) {freqScalePanel.setNDiv(f); timefreqMapCanvas.setNDiv(f);}
    public void setSlider(int v) {freqPickSlider.setValue(v);}

  }



    private void initTMPanelComponents() 
    {
        sigex = new JPanel();

        ButtonGroup plotGroup = new ButtonGroup();
        freqChangePanel = new JPanel(); freqIntsLabel = new JLabel(); maxFreqLabel = new JLabel(); 
        tcimfPanel = new IMFPlotPanel(); 

        plotIMFPanel = new JPanel();
        plotnIMFradioButton = new JRadioButton();
        plotAMradioButton = new JRadioButton();
        plotDataCheckBox = new JCheckBox();
        plotIMFradioButton = new JRadioButton();
        
        plotnIMFradioButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        plotnIMFradioButton.setText("nIMF"); plotGroup.add(plotnIMFradioButton);
        plotnIMFradioButton.setToolTipText("Plot the detrended IMFs associated with frequnecy given from Time-Frequency Panel");
        plotnIMFradioButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                plotButtonActionPerformed(evt);
            }
        });
        plotnIMFradioButton.setSelected(true);


        plotAMradioButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        plotAMradioButton.setText("AM"); plotGroup.add(plotAMradioButton); 
        plotAMradioButton.setToolTipText("Plot the AM component associated with the detrended IMF at the frequnecy given from the Time-Frequency Panel");
        plotAMradioButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                plotButton1ActionPerformed(evt);
            }
        });

        plotDataCheckBox.setFont(new Font("Arial", 0, 12)); // NOI18N
        plotDataCheckBox.setText("Data"); plotDataCheckBox.setSelected(true);
        plotDataCheckBox.setToolTipText("Plot the original time series data with the IMFs");
        plotDataCheckBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                plotButton2ActionPerformed(evt);
            }
        });
        

        plotIMFradioButton.setFont(new Font("Arial", 0, 12)); // NOI18N
        plotIMFradioButton.setText("IMF"); plotGroup.add(plotIMFradioButton);
        plotIMFradioButton.setToolTipText("Plot the IMF plus residual trend associated with frequnecy given from Time-Frequency Panel");
        plotIMFradioButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                plotButton3ActionPerformed(evt);
            }
        });


        freqChangePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        freqIntsLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        freqIntsLabel.setText("Intervals"); 
        freqIntsLabel.setToolTipText("Number of frequency intervals in the computation for the realTime-Frequency decomposition");

        maxFreqLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        maxFreqLabel.setText("Max \u03C9");
        maxFreqLabel.setToolTipText("The maximum frequency used for extracting the trend/cycles in the time series data");

     
        GroupLayout freqChangePanelLayout = new GroupLayout(freqChangePanel);
        freqChangePanel.setLayout(freqChangePanelLayout);
        freqChangePanelLayout.setHorizontalGroup(
            freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(freqChangePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                    .addComponent(freqIntsLabel)
                    .addComponent(freqIntsSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(maxFreqSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(maxFreqLabel))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        freqChangePanelLayout.setVerticalGroup(
            freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(freqChangePanelLayout.createSequentialGroup()
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(freqIntsLabel)
                    .addComponent(maxFreqLabel))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqChangePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(maxFreqSlider, GroupLayout.PREFERRED_SIZE, 141, GroupLayout.PREFERRED_SIZE)
                    .addComponent(freqIntsSlider, GroupLayout.PREFERRED_SIZE, 141, GroupLayout.PREFERRED_SIZE)))
        );

        plotIMFPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        GroupLayout plotIMFPanelLayout = new GroupLayout(plotIMFPanel);
        plotIMFPanel.setLayout(plotIMFPanelLayout);
        plotIMFPanelLayout.setHorizontalGroup(
            plotIMFPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, plotIMFPanelLayout.createSequentialGroup()
                .addContainerGap(10, Short.MAX_VALUE)
                .addGroup(plotIMFPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(plotDataCheckBox)
                    .addComponent(plotAMradioButton)
                    .addComponent(plotnIMFradioButton, GroupLayout.PREFERRED_SIZE, 98, GroupLayout.PREFERRED_SIZE)
                    .addComponent(computeTimeFreq)
                    .addComponent(plotIMFradioButton)))
        );
        plotIMFPanelLayout.setVerticalGroup(
            plotIMFPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, plotIMFPanelLayout.createSequentialGroup()
                .addGap(20, 20, 20)
                .addComponent(plotDataCheckBox)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(plotAMradioButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotIMFradioButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotnIMFradioButton)
                .addGap(18, 18, 18)
                .addComponent(computeTimeFreq)
                .addContainerGap())
        );

        GroupLayout layout = new GroupLayout(sigex);
        sigex.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(freqChangePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 805, Short.MAX_VALUE)
                .addComponent(plotIMFPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
            .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(156, 156, 156)
                    .addComponent(tcimfPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addContainerGap(0, Short.MAX_VALUE)))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addComponent(freqChangePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(plotIMFPanel, GroupLayout.DEFAULT_SIZE, 0, Short.MAX_VALUE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(tcimfPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addContainerGap()))

        );
    }
    private void plotButtonActionPerformed(ActionEvent evt) 
    {tcimfPanel.plotType(0); tcimfPanel.setIMFData(mdfa.amMap, mdfa.fmMap, mdfa.n_imfs, mdfa.flength);}
    private void plotButton1ActionPerformed(ActionEvent evt) 
    {tcimfPanel.plotType(1); tcimfPanel.setIMFData(mdfa.amMap, mdfa.fmMap, mdfa.n_imfs, mdfa.flength);}
    private void plotButton3ActionPerformed(ActionEvent evt) 
    {tcimfPanel.plotType(2); tcimfPanel.setIMFData(mdfa.imfs, mdfa.fmMap, mdfa.n_imfs, mdfa.flength);}
    private void plotButton2ActionPerformed(ActionEvent evt) 
    {tcimfPanel.plotData();}

    private void initARTPanelComponents() 
    {
        ARTpanel = new JPanel();

        reliablePanel = new JPanel(); reliableText = new JTextField(); reliableLabel = new JLabel();
        timeLabel = new JLabel(); timeText = new JTextField(); accLabel = new JLabel(); accText = new JTextField();
        mseLabel = new JLabel(); mseText = new JTextField(); artPanel = new JPanel();
        ARTCanvas = new artCanvas();

        reliablePanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        reliableText.setFont(new Font("Arial", 0, 12)); // NOI18N
        reliableText.setText("0.0"); reliableText.setColumns(4);

        reliableLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        reliableLabel.setText("Reliability:");

        timeLabel.setFont(new Font("Arial", 0, 12)); timeLabel.setText("Timeliness:");
        timeText.setFont(new Font("Arial", 0, 12)); timeText.setText("0.0"); timeText.setColumns(4);
        accLabel.setFont(new Font("Arial", 0, 12)); accLabel.setText("Accuracy:");
        accText.setFont(new Font("Arial", 0, 12)); accText.setText("1.0");  accText.setColumns(4);
        mseLabel.setFont(new Font("Arial", 0, 12)); mseLabel.setText("MSE:");
        mseText.setFont(new Font("Arial", 0, 12)); mseText.setText("0.0"); mseText.setColumns(6);

        turningLabel = new JLabel(); turningText = new JTextField();

        aerrorLabel = new JLabel(); rerrorLabel = new JLabel(); phasedelayLabel = new JLabel();
        aerrorText = new JTextField(); rerrorText = new JTextField(); phasedelayText = new JTextField();
        performLabel = new JLabel(); restartTimeButton = new JButton(); restartTimeButton1 = new JButton();
        restartTimeButton2 = new JButton();

        aerrorText.setColumns(6); rerrorText.setColumns(6); phasedelayText.setColumns(6); 

        turningLabel.setText("Turning");

        turningText.setText("0.0"); turningText.setColumns(4);

        aerrorLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        aerrorLabel.setText("Accuracy Error ");
        aerrorLabel.setToolTipText("A measurement of the log-error associated in the pass-band of filter fRf");

        rerrorLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        rerrorLabel.setText("Reliabiliy Error");
        rerrorLabel.setToolTipText("A measurment of the log-error in the stop-band of the filter fRf");

        phasedelayLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        phasedelayLabel.setText("Phase Delay Int");
        phasedelayLabel.setToolTipText("A measurement of the phase delay integral in the pass-band");

        mseLabel.setFont(new Font("Arial", 0, 12)); // NOI18N
        mseLabel.setText("MS Error");
        mseLabel.setToolTipText(" A measure of the error contributed by all the comonents on the frequency interval");

        aerrorText.setFont(new Font("Arial", 0, 12)); // NOI18N
        aerrorText.setText("0.0");

        rerrorText.setFont(new Font("Arial", 0, 12)); // NOI18N
        rerrorText.setText("0.0");

        phasedelayText.setFont(new Font("Arial", 0, 12)); // NOI18N
        phasedelayText.setText("0.0");

        performLabel.setText("ART Performance");


        GroupLayout ARTCanvasLayout = new GroupLayout(ARTCanvas);
        ARTCanvas.setLayout(ARTCanvasLayout);
        ARTCanvasLayout.setHorizontalGroup(
            ARTCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 270, Short.MAX_VALUE)
        );
        ARTCanvasLayout.setVerticalGroup(
            ARTCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        GroupLayout artPanelLayout = new GroupLayout(artPanel);
        artPanel.setLayout(artPanelLayout);
        artPanelLayout.setHorizontalGroup(
            artPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(artPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(ARTCanvas, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        artPanelLayout.setVerticalGroup(
            artPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(artPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(ARTCanvas, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );




        GroupLayout reliablePanelLayout = new GroupLayout(reliablePanel);
        reliablePanel.setLayout(reliablePanelLayout);
        reliablePanelLayout.setHorizontalGroup(
            reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(reliablePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(reliablePanelLayout.createSequentialGroup()
                        .addComponent(aerrorLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(aerrorText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                    .addGroup(reliablePanelLayout.createSequentialGroup()
                        .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(reliablePanelLayout.createSequentialGroup()
                                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                    .addComponent(rerrorLabel)
                                    .addComponent(phasedelayLabel)
                                    .addComponent(mseLabel))
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                    .addComponent(phasedelayText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(rerrorText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(mseText, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                            .addComponent(performLabel))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        reliablePanelLayout.setVerticalGroup(
            reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(reliablePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(performLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(aerrorLabel)
                    .addComponent(aerrorText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(rerrorText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(rerrorLabel))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(phasedelayLabel)
                    .addComponent(phasedelayText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGap(12, 12, 12)
                .addGroup(reliablePanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(mseLabel)
                    .addComponent(mseText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGap(0, 17, Short.MAX_VALUE))
        );

        restartTimeButton.setText("Reset");
        restartTimeButton.setToolTipText("Reset all values so Accuracy = 1.0");

        restartTimeButton1.setText("Reset");
        restartTimeButton1.setToolTipText("Reset Reliability = 0.0 ");

        restartTimeButton2.setText("Reset");
        restartTimeButton2.setToolTipText("Reset Timeliness = 0.0");

        restartTimeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                restartActionPerformed(evt);
            }
        });
        restartTimeButton1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                restart1ActionPerformed(evt);
            }
        });
        restartTimeButton2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                restart2ActionPerformed(evt);
            }
        });

        GroupLayout layout = new GroupLayout(ARTpanel);
        ARTpanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(artPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(accLabel)
                            .addComponent(reliableLabel)
                            .addComponent(timeLabel))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(timeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(restartTimeButton2))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(reliableText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(restartTimeButton1))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(accText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(restartTimeButton))))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(turningLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(turningSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(turningText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 189, Short.MAX_VALUE)
                .addComponent(reliablePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                    .addGroup(GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(artPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(reliablePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                        .addGap(19, 19, 19)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                            .addComponent(turningText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                                    .addComponent(accLabel)
                                    .addComponent(accText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(restartTimeButton))
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                                    .addComponent(reliableLabel)
                                    .addComponent(reliableText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(restartTimeButton1))
                                .addGap(6, 6, 6)
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                                    .addComponent(timeLabel)
                                    .addComponent(timeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(restartTimeButton2))
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGap(18, 18, 18)
                                        .addComponent(turningLabel))
                                    .addGroup(layout.createSequentialGroup()
                                        .addGap(1, 1, 1)
                                        .addComponent(turningSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))))))
                .addContainerGap())
        );
    }
    private void restartActionPerformed(ActionEvent evt) {setART(0.0, 0.0,true);}
    private void restart1ActionPerformed(ActionEvent evt) {setART(0.0, -1000.0,true);}
    private void restart2ActionPerformed(ActionEvent evt) {setART(-1000.0, 0.0,true);}









    /*--------------------------------------------------------------------------------------------------
       ZPC Filter Panel 
    ----------------------------------------------------------------------------------------------------*/

    private void initZPCPanel() 
    {
        zpcPanel = new JPanel(); 
  
        alphaBar = new JScrollBar(); amzpcLabel = new JLabel();
        alphaLabel = new JLabel(); alphaText = new JTextField(); lamzpcText = new JTextField();
        customLabel = new JLabel(); lamzpcBar = new JScrollBar(); hfText = new JTextField();
        hfLabel = new JLabel(); hfBar = new JScrollBar(); i1Text = new JTextField();
        i1Bar = new JScrollBar(); i1label = new JLabel(); argz1Label = new JLabel(); arg1Text = new JTextField();
        zpcLabel = new JLabel(); arg1Bar = new JScrollBar(); arg2Text1 = new JTextField();
        modzBar = new JScrollBar(); modpLabel = new JLabel(); arg1Text1 = new JTextField();
        modzLabel = new JLabel(); modpBar = new JScrollBar(); modzText = new JTextField(); modpText = new JTextField();
        optimizeButton = new JButton(); optimizeCheck = new JCheckBox();
        msezpcLabel = new JLabel(); msezpcText = new JTextField();
        b0Bar = new JScrollBar(); b0Label = new JLabel(); b0Text = new JTextField();
        injectButton = new JButton(); lamzpcLabel = new JLabel(); prefilterButton = new JButton();
     


        i1const = new JCheckBox("i1:"); i1const.setHorizontalTextPosition(JMenuItem.LEFT); i1const.setSelected(false);
        i2const = new JCheckBox("i2:"); i2const.setHorizontalTextPosition(JMenuItem.LEFT); i2const.setSelected(false);
        i1const.setToolTipText("Apply normalization on ZPC constant so that \u0393(0).");   
        i2const.setToolTipText("Inject with another ZPC structure to ensure that the time delay vanishes at 0. Control with i2 pole scrollbar above.");      
        normalizeConst = new JCheckBox("Normalize"); normalizeConst.setHorizontalTextPosition(JMenuItem.LEFT); i1const.setSelected(false);
        normalizeConst.setToolTipText("Apply normalization on ZPC constant for non i1 process so that \u0393(\u03C9_0) is close to 1."); 


        JLabel zpcComboLabel = new JLabel("ZPC Order");
        zpcComboLabel.setFont(new Font("Ubuntu", 0, 12));
        zpcCombo = new JComboBox<String>();
        zpcCombo.addItem("ARMA (2,2)");
        zpcCombo.addItem("ARMA (4,4)");
        zpcCombo.addActionListener(new ActionListener() {
           public void actionPerformed(ActionEvent evt) 
           {
              if(zpcCombo.getSelectedIndex() == 0)
              {p_arma = 2; q_arma = 2; zpc.setOrders(p_arma, q_arma);}
              else if(zpcCombo.getSelectedIndex() == 1)
              {p_arma = 4; q_arma = 4; zpc.setOrders(p_arma, q_arma);}
           }
        });

        zpcFreqCanvas.setBackground(new Color(0, 0, 0));
        zpcFreqCanvas.setPreferredSize(new Dimension(330,184));
        zpc_getData = new JButton("Reset Data");
        zpc_getData.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) { getData(); }
        });


        /*GroupLayout zpcFreqCanvasLayout = new GroupLayout(zpcFreqCanvas);
        zpcFreqCanvas.setLayout(zpcFreqCanvasLayout);
        zpcFreqCanvasLayout.setHorizontalGroup(
            zpcFreqCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 326, Short.MAX_VALUE)
        );
        zpcFreqCanvasLayout.setVerticalGroup(
            zpcFreqCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 184, Short.MAX_VALUE)
        );*/

        lamzpcLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        lamzpcLabel.setText("\u03BB");
        lamzpcLabel.setToolTipText("Set the lambda to emphasize better time delay properties of filter.");

        alphaLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        alphaLabel.setText("\u03B1");
        alphaLabel.setToolTipText("Set the smoothness properties by emphasizing more weight in stop-band.");

        alphaText.setText("0");

        lamzpcText.setText("0");

        customLabel.setFont(new Font("Ubuntu", 0, 14)); // NOI18N
        customLabel.setText("ZPC  Customization");

        hfText.setText("0");

        hfLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        hfLabel.setText("H-F \u0394");
        hfLabel.setToolTipText("Emphasizes control over dampening differences in high-frequency noise.");

        i1Text.setText("0");

        i1label.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        i1label.setText("i1 Pole");
        i1label.setToolTipText("Emphasize control over level-change properties of filter.");

        argz1Label.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        argz1Label.setText("argZP");
        argz1Label.setToolTipText("Change the complex argument of the first zero/pole. If nonzero, uses the complex conjugate.");

        

        zpcLabel.setFont(new Font("Ubuntu", 0, 14));  zpcLabel.setText("Zero-Pole Arg/Mod"); arg2Text1.setText("0");

        modpLabel.setFont(new Font("Ubuntu", 0, 12));   modpLabel.setText("modP"); 
        modpLabel.setToolTipText("Change the modulus of the pole. If real roots, changes 1st of 2nd real root.");

  

        modzLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        modzLabel.setText("modZ");
        modzLabel.setToolTipText("Change the modulus of the zero root. If real roots, changes 1st of 2nd real root.");

        optimizeButton.setFont(new Font("Ubuntu", 0, 14)); // NOI18N
        optimizeButton.setText("Optimize Filter");
        optimizeButton.setToolTipText("Optimize over both arguments and moduli of zeros/poles.");
        optimizeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                optimizeActionPerformed(evt);
            }
        });


        optimizeCheck.setFont(new Font("Ubuntu", 0, 14)); // NOI18N
        optimizeCheck.setText("Automatic Optimization");
        optimizeCheck.setToolTipText("Automatically optimize filter for any given changes to the filter parameters. Recommended only for fast machines.");

        msezpcLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        msezpcLabel.setText("MSE");

        msezpcText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        msezpcText.setText("0");

        b0Label.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        b0Label.setText("Step size");

        b0Text.setText("0.1");

        injectButton.setFont(new Font("Ubuntu", 0, 14)); // NOI18N
        injectButton.setText("Inject ZPC-Gene");
        injectButton.setToolTipText("Apply the ZPC-ARMA filter to the I-MDFA filter for a  hybrid filter with given ZPC properties.");
        injectButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                injectActionPerformed(evt);
            }
        });
       
        prefilterButton.setFont(new Font("Ubuntu", 0, 14)); // NOI18N
        prefilterButton.setText("Prefilter for IMDFA");
        prefilterButton.setToolTipText("Applies the IMDFA directly on ZPC-filter output");
        prefilterButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                complimentActionPerformed(evt);
            }
        });   
 
        //-------------------------- set up scrollbars --------------


       alphaBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,300);
       alphaBar.setValue(0); alphaBar.setUnitIncrement(1); 
       alphaText.setColumns(4); alphaText.setText(df.format(0.0));
       alphaBar.setPreferredSize(new Dimension(130, 15));
      
       hfBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,100);
       hfBar.setValue(0); hfBar.setUnitIncrement(10); 
       hfText.setColumns(4); hfText.setText(df.format(0.0));
       hfBar.setPreferredSize(new Dimension(130, 15));

       i1Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,12,1,1000);
       i1Bar.setValue(200); i1Bar.setUnitIncrement(1); 
       i1Text.setColumns(4); i1Text.setText(df.format(0.0));
       i1Bar.setPreferredSize(new Dimension(130, 15));  

       modzBar = new JScrollBar(JScrollBar.HORIZONTAL,500,10,1,999);
       modzBar.setValue(500); modzBar.setUnitIncrement(1); 
       modzText.setColumns(4);  modzText.setText(df.format(0.5));
       modzBar.setPreferredSize(new Dimension(130, 15));

       modpBar = new JScrollBar(JScrollBar.HORIZONTAL,500,10,1,999);
       modpBar.setValue(500); modzBar.setUnitIncrement(1); 
       modpText.setColumns(4); modpText.setText(df.format(0.5));
       modpBar.setPreferredSize(new Dimension(130, 15));
 

       b0Bar = new JScrollBar(JScrollBar.HORIZONTAL,1,10,1,50);
       b0Bar.setValue(10); b0Bar.setUnitIncrement(1); 
       b0Text.setColumns(4);  b0Text.setText(df.format(0.1));
       b0Bar.setPreferredSize(new Dimension(130, 15));

       lamzpcBar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,1000);
       lamzpcBar.setValue(0); lamzpcBar.setUnitIncrement(1); 
       lamzpcText.setColumns(4); lamzpcText.setText(df.format(0.0));
       lamzpcBar.setPreferredSize(new Dimension(130, 15));
  
       AdjustmentListener adj = new AdjustmentListener()  {
        public void adjustmentValueChanged(AdjustmentEvent e) {
             
             if(e.getAdjustable() == alphaBar)
             {
                zpc.set_exp((double)alphaBar.getValue()*.1);
                alphaText.setText(df.format((double)alphaBar.getValue()*.1));
                if(autoCompZPC) {zpc.getModOptimizedZPC(); injectZPCGeneSimple(); setZPC(); updatePlots(false,true);}
             }
             else if(e.getAdjustable() == hfBar)
             {   
                zpc.set_lambda3((double)hfBar.getValue()*.01);
                hfText.setText(df.format((double)hfBar.getValue()*.01));
                if(autoCompZPC) {zpc.getModOptimizedZPC(); injectZPCGeneSimple(); setZPC(); updatePlots(false,true);}
             }
             else if(e.getAdjustable() == i1Bar)
             {
                zpc.set_P((double)i1Bar.getValue()*.001);
                i1Text.setText(df.format((double)i1Bar.getValue()*.001));
                if(autoCompZPC) {zpc.getModOptimizedZPC(); injectZPCGeneSimple(); setZPC(); updatePlots(false,true);}

             }

             else if(e.getAdjustable() == modzBar)
             {

               modz = (double)modzBar.getValue()*.001;
               //modp = (double)modpBar.getValue()*.001;               
               modzText.setText(df2.format((double)modzBar.getValue()*.001));
               zpc.setInitMod(modz,modp);
               if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC(); updatePlots(false,true);}   
             }
             else if(e.getAdjustable() == modpBar)
             {

               //modz = (double)modzBar.getValue()*.001;
               modp = (double)modpBar.getValue()*.001;               
               modpText.setText(df2.format((double)modpBar.getValue()*.001));
               zpc.setInitMod(modz,modp);
               if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC(); updatePlots(false,true);}      

             }
             else if(e.getAdjustable() == b0Bar)
             {
               zpc.setInitStep((double)b0Bar.getValue()*.01);  
               b0Text.setText(df2.format((double)b0Bar.getValue()*.01));              
               if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC(); updatePlots(false,true);}  
             }
             else if(e.getAdjustable() == lamzpcBar)
             {
                zpc.set_lambda((double)lamzpcBar.getValue()*.10);
                lamzpcText.setText(df.format((double)lamzpcBar.getValue()*.10));
                if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC(); updatePlots(false,true);}
             }
         }
        };

        ItemListener zpcItem = new ItemListener() { 
          public void itemStateChanged(ItemEvent e)
          {
           boolean sel; 
           Object source = e.getItemSelectable();
           if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
           else{sel = true;}

           if(source == optimizeCheck) 
           {
             autoCompZPC = sel; if(sel) {zpc.getModOptimizedZPC(); setZPC();}
           }
           else if(source == i1const)
           {
             if(sel) {zpc.set_i1(1); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
             else {zpc.set_i1(0); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
           } 
           else if(source == i2const)
           {
              if(sel) {i1Bar.setEnabled(true);} 
              else {i1Bar.setEnabled(false);}
              if(sel) {zpc.set_i2(1); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
              else {zpc.set_i2(0); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
           }
           else if(source == normalizeConst)
           {
              if(sel) {zpc.setNormalize(1); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
              else {zpc.setNormalize(0); if(autoCompZPC) {zpc.getModOptimizedZPC(); setZPC();}}
           }           
          }
        };


        alphaBar.addAdjustmentListener(adj);
        lamzpcBar.addAdjustmentListener(adj);
        hfBar.addAdjustmentListener(adj);
        i1Bar.addAdjustmentListener(adj);
        
        modzBar.addAdjustmentListener(adj);
        modpBar.addAdjustmentListener(adj);
        b0Bar.addAdjustmentListener(adj);
        optimizeCheck.addItemListener(zpcItem);
        i1const.addItemListener(zpcItem);  
        i2const.addItemListener(zpcItem);
        normalizeConst.addItemListener(zpcItem);


        GroupLayout paramLayout;
        JPanel zpcParamPanel = new JPanel();
        JPanel zpcCustomPanel = new JPanel();
        JPanel zpcplotPanel = new JPanel();
        JPanel optimPanel = new JPanel();
         
        zpcParamPanel = new JPanel();
        //zpcParamPanel.setBorder(new TitledBorder(new LineBorder(myBlue),"ZPC Parameters"));
        paramLayout = new GroupLayout(zpcParamPanel);

        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup()
            .addComponent(zpcComboLabel).addComponent(modzLabel).addComponent(modpLabel).addComponent(b0Label))
          .addGroup(paramLayout.createParallelGroup()
            .addComponent(zpcCombo).addComponent(modzBar).addComponent(modpBar).addComponent(b0Bar))
          .addGroup(paramLayout.createParallelGroup()
            .addComponent(modzText).addComponent(modpText).addComponent(b0Text)));
        paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(paramLayout.createSequentialGroup()
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(zpcComboLabel).addComponent(zpcCombo))
                 .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(modzLabel).addComponent(modzBar).addComponent(modzText)) 
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(modpLabel).addComponent(modpBar).addComponent(modpText))
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(b0Label).addComponent(b0Bar).addComponent(b0Text))) ));
        zpcParamPanel.setLayout(paramLayout);


        zpcCustomPanel = new JPanel(); 
        //zpcCustomPanel.setBorder(new TitledBorder(new LineBorder(myBlue),"Customization"));
        paramLayout = new GroupLayout(zpcCustomPanel);

        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup()
            .addComponent(lamzpcLabel).addComponent(alphaLabel).addComponent(hfLabel).addComponent(i1label))
          .addGroup(paramLayout.createParallelGroup()
            .addComponent(lamzpcBar).addComponent(alphaBar).addComponent(hfBar).addComponent(i1Bar))
          .addGroup(paramLayout.createParallelGroup()
             .addComponent(lamzpcText).addComponent(alphaText).addComponent(hfText).addComponent(i1Text)));
        paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(paramLayout.createSequentialGroup()
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lamzpcLabel).addComponent(lamzpcBar).addComponent(lamzpcText))
                 .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(alphaLabel).addComponent(alphaBar).addComponent(alphaText)) 
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(hfLabel).addComponent(hfBar).addComponent(hfText))
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(i1label).addComponent(i1Bar).addComponent(i1Text)))));
        zpcCustomPanel.setLayout(paramLayout);


        JPanel zpcConstraintPanel = new JPanel(); 
        //zpcplotPanel.setBorder(new TitledBorder(new LineBorder(myBlue),"ZPC Plots"));
        paramLayout = new GroupLayout(zpcConstraintPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
             .addComponent(i1const) 
             .addComponent(i2const)
             .addComponent(normalizeConst));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createSequentialGroup()
               .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(i1const) 
             .addComponent(i2const)
             .addComponent(normalizeConst))));
        zpcConstraintPanel.setLayout(paramLayout);






        zpcplotPanel = new JPanel(); 
        //zpcplotPanel.setBorder(new TitledBorder(new LineBorder(myBlue),"ZPC Plots"));
        paramLayout = new GroupLayout(zpcplotPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
             .addComponent(plot_zpc[0]) 
             .addComponent(plot_zpc[1])
             .addComponent(plot_zpc[2])
             .addComponent(plot_zpc[3]));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createSequentialGroup()
               .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(plot_zpc[0]) 
             .addComponent(plot_zpc[1])
             .addComponent(plot_zpc[2])
             .addComponent(plot_zpc[3]))));
        zpcplotPanel.setLayout(paramLayout);

       
        JPanel mseComp = new JPanel();
        mseComp.add(msezpcLabel); mseComp.add(msezpcText); mseComp.setLayout(new GridLayout(1,2,1,1));

        optimPanel = new JPanel();
        optimPanel.setBorder(new TitledBorder(new LineBorder(myBlue),"Optimization"));
        paramLayout = new GroupLayout(optimPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
          .addComponent(optimizeButton)
          .addComponent(injectButton)
          .addComponent(prefilterButton)
          //.addComponent(zpc_getData)
          .addComponent(optimizeCheck)
          .addComponent(mseComp)));

        paramLayout.setVerticalGroup(paramLayout.createSequentialGroup() 
          .addComponent(optimizeButton)
          .addComponent(injectButton)
          .addComponent(prefilterButton)
          //.addComponent(zpc_getData)
          .addComponent(optimizeCheck)
          .addComponent(mseComp));
        optimPanel.setLayout(paramLayout);
         
        Box zpcControl1 = Box.createHorizontalBox(); 
        zpcControl1.add(zpcParamPanel);
        zpcControl1.add(zpcCustomPanel);

        Box zpcControl2 = Box.createHorizontalBox(); 
        zpcControl2.add(zpcplotPanel);
        zpcControl2.add(zpcConstraintPanel);

        JPanel zpcControl = new JPanel(); 
        paramLayout = new GroupLayout(zpcControl);
        zpcControl.setLayout(paramLayout);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(zpcControl1).addComponent(zpcControl2)));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
           .addComponent(zpcControl1).addComponent(zpcControl2));

        
        
        paramLayout = new GroupLayout(zpcPanel);
        zpcPanel.setLayout(paramLayout);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
             .addComponent(zpcFreqCanvas).addComponent(zpcControl).addComponent(optimPanel));
        paramLayout.setVerticalGroup( paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(zpcFreqCanvas).addComponent(zpcControl).addComponent(optimPanel)));
 


    }

    private void optimizeActionPerformed(ActionEvent evt) {zpc.getModOptimizedZPC(); setZPC();}
    private void injectActionPerformed(ActionEvent evt) {zpc_gene = true; injectZPCGene();}
    private void complimentActionPerformed(ActionEvent evt) {zpc_gene = true; computeZPCtoIMDFA();}
 
    //------------------ Convolve current ZPC filter with the I-MDFA filter -----
 
    public void injectZPCGene()
    {

       filterPlot[1].setEnabled(true); filterPlot[2].setEnabled(true);
       gamma_zpc = new double[(n_rep+1)*K1];
       gamma_hybrid = new double[(n_rep+1)*K1];
       //update b,tseries,and hybrid filters
       zpc.setIMDFAFilter(mdfa.b, mdfa.n_rep, mdfa.L);
       zpc.setData(mdfa.tseries, mdfa.n_rep, mdfa.n_obs);
       zpc.injectZPCGene2();
      
       //---update timedomain canvas
       mdfa_canvas.setHybridZPC(zpc.xf_hybrid);
       mdfa_canvas.setZPC(zpc.xf_zpc);

 
     if(reCompFilter)
     {
        //System.out.println("n_rep = " + n_rep); 
       //---update frequency canvas -------------------------------     
       for(int i=0; i < n_rep; i++)  //---- get frf
       {   
         for(int k=0;k<=K;k++)
         {
           if(amplitudeBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k];  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*mdfa.amp_filter[K1*i+k];                  
           }
           else if(gammaBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k]);  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k])*mdfa.gamma_hat[K1*i+k];             
           }
           else if(timeDelayBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.time_delay[k];  
              gamma_hybrid[K1*i+k] = zpc.time_delay[k]*mdfa.time_delay[K1*i+k];       
           }                  
         }
       }  
       filter_canvas.setGamma_zpc(gamma_zpc);
       filter_canvas.setGamma_hybrid(gamma_hybrid);


       //---update periodogram canvas--------------------------------
       period_canvas.setPeriodogramXf(mdfa.period_xf, nsamp);
       period_canvas.setPeriodogram_zpc(zpc.period_xf);
       period_canvas.setPeriodogram_hybrid(zpc.period_xf);
     }
   }
  
 
   public void computeZPCtoIMDFA()
   {
       
       zpc.getModOptimizedZPC(); 
       setZPC();
       
       mdfa.setSpecDensity(zpc.fRf, 1);
       
       zpc.setIMDFAFilter(mdfa.b, n_rep, L);
       zpc.injectZPCGene(); 

       //mdfa.computeFilterFromZPC(zpc.m_zpc, mdfa.n_obs);
       mdfa.computeFilterGeneral(true,false);
       
       mdfa_canvas.setHybridZPC(mdfa.xf);
       mdfa_canvas.setZPC(zpc.xf_zpc);
       filterPlot[1].setEnabled(true); filterPlot[2].setEnabled(true);

     if(reCompFilter)
     {
      if((plot_number == 1) || (plot_number == 2))
      {

       gamma_zpc = new double[(mdfa.n_rep+1)*mdfa.K1];
       gamma_hybrid = new double[(mdfa.n_rep+1)*mdfa.K1];
       //---update frequency canvas -------------------------------  

       //System.out.println("n_rep = " + n_rep);   
       for(int i=0; i < n_rep; i++)  //---- get frf
       {   
         for(int k=0;k<=mdfa.K;k++)
         {
           if(amplitudeBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k];  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*mdfa.amp_filter[K1*i+k];                  
           }
           else if(gammaBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k]);  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k])*mdfa.gamma_hat[K1*i+k];             
           }
           else if(timeDelayBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.time_delay[k];  
              gamma_hybrid[K1*i+k] = zpc.time_delay[k]*mdfa.time_delay[K1*i+k];       
           }                  
         }
       }  

       filter_canvas.setGamma_zpc(gamma_zpc);
       filter_canvas.setGamma_hybrid(gamma_hybrid);


       //---update periodogram canvas--------------------------------
       period_canvas.setPeriodogramXf(mdfa.period_xf, nsamp);
       period_canvas.setPeriodogram_zpc(zpc.period_xf);
       period_canvas.setPeriodogram_hybrid(zpc.period_xf);
      }
     }






   }
  

   public void injectZPCGeneSimple()
   {
       
       //zpc.setIMDFAFilter(mdfa.b, mdfa.n_rep, mdfa.L);
       //zpc.setData(mdfa.tseries, mdfa.n_rep, mdfa.n_obs);
       zpc.setIMDFAFilter(mdfa.b, n_rep, L);
       zpc.injectZPCGene2(); 
   
       //---update timedomain canvas
       mdfa_canvas.setHybridZPC(zpc.xf_hybrid);
       mdfa_canvas.setZPC(zpc.xf_zpc);

     if(reCompFilter)
     {
      if((plot_number == 1) || (plot_number == 2))
      {

       gamma_zpc = new double[(mdfa.n_rep+1)*mdfa.K1];
       gamma_hybrid = new double[(mdfa.n_rep+1)*mdfa.K1];
       //---update frequency canvas -------------------------------  

       //System.out.println("n_rep = " + n_rep);   
       for(int i=0; i < n_rep; i++)  //---- get frf
       {   
         for(int k=0;k<=mdfa.K;k++)
         {
           if(amplitudeBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k];  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*mdfa.amp_filter[K1*i+k];                  
           }
           else if(gammaBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k]);  
              gamma_hybrid[K1*i+k] = zpc.amp[k]*Math.cos(zpc.phase[k])*mdfa.gamma_hat[K1*i+k];             
           }
           else if(timeDelayBut.isSelected())
           {
              gamma_zpc[K1*i+k] = zpc.time_delay[k];  
              gamma_hybrid[K1*i+k] = zpc.time_delay[k]*mdfa.time_delay[K1*i+k];       
           }                  
         }
       }  

       filter_canvas.setGamma_zpc(gamma_zpc);
       filter_canvas.setGamma_hybrid(gamma_hybrid);


       //---update periodogram canvas--------------------------------
       period_canvas.setPeriodogramXf(mdfa.period_xf, nsamp);
       period_canvas.setPeriodogram_zpc(zpc.period_xf);
       period_canvas.setPeriodogram_hybrid(zpc.period_xf);
      }
     }
   }

   //------ With MDFA set, compute the ZPC filter---------

   public void computeZPCFilter()
   {
     double[] In = new double[K1]; for(int k = 0; k < K1; k++) {In[k] = mdfa.period_xf[k];}
     zpc.set_Gamma(mdfa.Gamma); 
     zpcFreqCanvas.setGamma(mdfa.Gamma);
    
     zpc.setOrders(p_arma, q_arma);
     if(i1const.isSelected()) {zpc.set_i1(1);} else {zpc.set_i1(0);}
     if(i2const.isSelected()) {zpc.set_i2(1);} else {zpc.set_i2(0);}
     if(normalizeConst.isSelected()) {zpc.setNormalize(1);} else {zpc.setNormalize(0);} 
  
     zpc.setPeriodogram(In);
     zpc.set_cutoff(mdfa.cutoff);  
     zpc.set_cutoff0(mdfa.cutoff0); 
     zpc.setIMDFAFilter(mdfa.b, mdfa.n_rep, mdfa.L);
     zpc.setData(mdfa.tseries, mdfa.n_rep, mdfa.n_obs);
     zpc.getModOptimizedZPC();
     setZPC();
   }


   public void initZPCFilter()
   {
     //System.arraycopy(mdfa.period_xf,0,zpc.In,0,K1);  
     //zpcFilter zpc = new zpcFilter(n_obs);

     double[] In = new double[K1]; 
     for(int k = 0; k < K1; k++) {In[k] = mdfa.period_xf[k];}
     zpc.set_Gamma(mdfa.Gamma); 
     zpc.setPeriodogram(In);
     zpc.set_cutoff(mdfa.cutoff); 
     zpc.set_cutoff0(mdfa.cutoff0); 

     zpc.setOrders(p_arma, q_arma);
     zpc.set_i1(0);
     zpc.set_i2(1);
     zpc.setNormalize(0);

     zpc.set_lambda(0);
     zpc.set_lambda3(0);
     zpc.set_exp(0); 
     zpc.setInitMod(.5,.5);
     zpc.setInitStep(.4);
     zpc.set_P(.2);
       
     zpc.setIMDFAFilter(mdfa.b, n_rep, L);
     zpc.setData(mdfa.tseries, n_rep, n_obs);
     zpc.getModOptimizedZPC();
     setZPC();
   }

   public void getData()
   {
      zpc.setIMDFAFilter(mdfa.b, n_rep, L);
      zpc.setData(mdfa.tseries, n_rep, n_obs);
   }


   public void updateZPC()
   {
      double[] In = new double[mdfa.K1]; 
      for(int k = 0; k < K1; k++) {In[k] = mdfa.period_xf[k];}
   
     zpc.set_Gamma(mdfa.Gamma); 
     zpc.setPeriodogram(In);
     zpc.set_cutoff(mdfa.cutoff); 
     zpc.set_cutoff0(mdfa.cutoff0); 

     zpc.setOrders(p_arma, q_arma);
     if(i1const.isSelected()) {zpc.set_i1(1);} else {zpc.set_i1(0);}
     if(i2const.isSelected()) {zpc.set_i2(1);} else {zpc.set_i2(0);}
     if(normalizeConst.isSelected()) {zpc.setNormalize(1);} else {zpc.setNormalize(0);} 

     zpc.set_lambda((double)lamzpcBar.getValue()*.10);
     zpc.set_lambda3((double)hfBar.getValue()*.01);
     zpc.set_exp((double)alphaBar.getValue()*.1); 
     zpc.setInitMod(modz,modp);
     zpc.setInitStep((double)b0Bar.getValue()*.01);
     zpc.set_P((double)i1Bar.getValue()*.001);
     

     zpc.setIMDFAFilter(mdfa.b, n_rep, L);
     zpc.setData(mdfa.tseries, n_rep, n_obs);
     zpc.getModOptimizedZPC();


      

   }


   //------ sets the frequency functions on the zpc canvas------
   public void setZPC()
   {
     zpcFreqCanvas.setZPCFilter(zpc.amp, zpc.phase, zpc.time_delay);
     if(zpc_mod_optimize) msezpcText.setText(df.format(zpc.zpc_min));
   }












   //---------------------- Diagnostics Panel----------------------------------

   private void initDiagnosticPanel() 
   {
        Font f1 = new Font("Ubuntu", 0, 12);
        diagnosticPanel = new JPanel();

        freqdiagPanel = new JPanel(); freqdiagLabel = new JLabel();
        d_mseLabel = new JLabel();  d_mseText = new JTextField(); d_icLabel = new JLabel();
        d_icText = new JTextField(6); d_dofText = new JTextField(6); d_dofLabel = new JLabel();
        d_accText = new JTextField(6); d_accLabel = new JLabel(); d_relLabel = new JLabel();
        d_relText = new JTextField(6); d_timeText = new JTextField(6); d_timeLabel = new JLabel();
        timediagPanel = new JPanel(); timediagLabel = new JLabel();
        t_mseLabel = new JLabel(); t_mseText = new JTextField(6); t_corrLabel = new JLabel();
        t_corrText = new JTextField(6); t_excessText = new JTextField(6); t_excessLabel = new JLabel();
        t_accText = new JTextField(6); t_accLabel = new JLabel(); t_relLabel = new JLabel();
        t_relText = new JTextField(6); t_timeText = new JTextField(6); t_timeLabel = new JLabel();
        optimizeARTPanel = new JPanel(); optimizeARTFilter = new JButton();
        accSlider = new JSlider(); relSlider = new JSlider(); timeSlider = new JSlider();
        optimLabel = new JLabel(); t_ARTLabel = new JLabel();

        freqdiagPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        freqdiagLabel.setFont(f1); // NOI18N
        freqdiagLabel.setText("Frequency Domain Diagnostics");

        d_mseLabel.setFont(f1); d_mseLabel.setText("MS Error");
        d_mseText.setFont(f1); d_mseText.setText("0");
        d_icLabel.setFont(f1); d_icLabel.setText("Info Criterion");
        d_icText.setFont(f1);  d_icText.setText("0");
        d_dofText.setFont(f1);  d_dofText.setText("0");
        d_dofLabel.setFont(f1);  d_dofLabel.setText("Effective DoF");
        d_accText.setFont(f1);  d_accText.setText("0");
        d_accLabel.setFont(f1); d_accLabel.setText("Accuracy");
        d_accLabel.setToolTipText("The frequency domain estimate of the accuracy of the filter. Based on periodogram. Same diagnostic given in ART Panel.");
        d_relLabel.setFont(f1); // NOI18N
        d_relLabel.setText("Reliability"); d_relLabel.setToolTipText("The frequency domain estimate of the reliability of the filter. Based on periodogram. Same diagnostic given in ART Panel.");
        d_relText.setFont(f1);  d_relText.setText("0");
        d_timeText.setFont(f1); d_timeText.setText("0");

        d_timeLabel.setFont(f1); d_timeLabel.setText("Timeliness");
        d_timeLabel.setToolTipText("The frequency domain estimate of the timeliness of the filter. Based on integral of phase delay. Same diagnostic given in ART Panel.");

        GroupLayout freqdiagPanelLayout = new GroupLayout(freqdiagPanel);
        freqdiagPanel.setLayout(freqdiagPanelLayout);
        freqdiagPanelLayout.setHorizontalGroup(
            freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(freqdiagPanelLayout.createSequentialGroup()
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(freqdiagPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_accLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_accText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_relLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_relText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_timeLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_timeText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_dofLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_dofText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_icLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_icText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, freqdiagPanelLayout.createSequentialGroup()
                                .addComponent(d_mseLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(d_mseText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))))
                    .addGroup(freqdiagPanelLayout.createSequentialGroup()
                        .addGap(12, 12, 12)
                        .addComponent(freqdiagLabel)))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        freqdiagPanelLayout.setVerticalGroup(
            freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(freqdiagPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(freqdiagLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_mseLabel)
                    .addComponent(d_mseText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_icLabel)
                    .addComponent(d_icText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_dofLabel)
                    .addComponent(d_dofText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_accLabel)
                    .addComponent(d_accText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_relLabel)
                    .addComponent(d_relText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(freqdiagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(d_timeLabel)
                    .addComponent(d_timeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(57, Short.MAX_VALUE))
        );

        timediagPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        timediagLabel.setFont(f1);  timediagLabel.setText("Time Domain Diagnostics");

        t_mseLabel.setFont(f1); t_mseLabel.setText("|| Y_t - Y_t ||_2");
        t_mseLabel.setToolTipText("The L_2 error of the computed filter and the symmetric filter. ");

        t_mseText.setFont(f1); t_mseText.setText("0");

        t_corrLabel.setFont(f1);  t_corrLabel.setText("Cov(Y_t,Y_t)");
        t_corrLabel.setToolTipText("Measure of the cross correlation between computed filter and symmetric filter.");

        t_corrText.setFont(f1);  t_corrText.setText("0");

        t_excessText.setFont(f1);  t_excessText.setText("0");

        t_excessLabel.setFont(f1);  t_excessLabel.setText("Excess TP"); t_excessLabel.setToolTipText("Excess number of turning points relative to symmetric filter. Meaures reliability of filter.");

        t_accText.setFont(f1); t_accText.setText("0");

        t_accLabel.setFont(f1);  t_accLabel.setText("Accuracy"); t_accLabel.setToolTipText("Estimate of the accuracy of filter relative to symmetric filter estimate. Minimum value is 0.0, maximum value is 1.0. ");

        t_relLabel.setFont(f1); t_relLabel.setText("Reliability"); t_relLabel.setToolTipText("Estimate of smoothness or reliability of filter. Minimum value is -1.0, maximum value is 1.0. ");

        t_relText.setFont(f1);  t_relText.setText("0");

        t_timeText.setFont(f1); t_timeText.setText("0");

        t_timeLabel.setFont(f1);  t_timeLabel.setText("Timeliness");
        t_timeLabel.setToolTipText("Estimate of the timeliness diagnostic on the time domain. Takes into account turning point detection. Minimum value is 0.0, maximum value is 1.0. ");

        GroupLayout timediagPanelLayout = new GroupLayout(timediagPanel);
        timediagPanel.setLayout(timediagPanelLayout);
        timediagPanelLayout.setHorizontalGroup(
            timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(timediagPanelLayout.createSequentialGroup()
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(timediagPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_accLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_accText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_relLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_relText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_timeLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_timeText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_excessLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_excessText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_corrLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_corrText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))
                            .addGroup(GroupLayout.Alignment.TRAILING, timediagPanelLayout.createSequentialGroup()
                                .addComponent(t_mseLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(t_mseText, GroupLayout.PREFERRED_SIZE, 62, GroupLayout.PREFERRED_SIZE))))
                    .addGroup(timediagPanelLayout.createSequentialGroup()
                        .addGap(12, 12, 12)
                        .addComponent(timediagLabel)))
                .addContainerGap(38, Short.MAX_VALUE))
        );
        timediagPanelLayout.setVerticalGroup(
            timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(timediagPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(timediagLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_mseLabel)
                    .addComponent(t_mseText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_corrLabel)
                    .addComponent(t_corrText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_excessLabel)
                    .addComponent(t_excessText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_accLabel)
                    .addComponent(t_accText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_relLabel)
                    .addComponent(t_relText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(timediagPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(t_timeLabel)
                    .addComponent(t_timeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        optimizeARTPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        optimizeARTFilter.setText("Optimize Filter");
        optimizeARTFilter.setToolTipText("Attempt to optimize filter with respect to in-sample symmetric filter using the ART priorities. May take up to one minute depending on processor speed.");
        optimizeARTFilter.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                optimizeARTFilterNow(evt);
            }
        });

        accSlider.setFont(f1); // NOI18N
        accSlider.setOrientation(JSlider.VERTICAL);
        accSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                accSliderStateChanged(evt);
            }
        });

        relSlider.setFont(f1); // NOI18N
        relSlider.setOrientation(JSlider.VERTICAL);
        relSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                relSliderStateChanged(evt);
            }
        });

        timeSlider.setFont(f1); // NOI18N
        timeSlider.setOrientation(JSlider.VERTICAL);
        timeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                timeSliderStateChanged(evt);
            }
        });

        optimLabel.setFont(f1); // NOI18N
        optimLabel.setText("Optimization of Filter");

        t_ARTLabel.setFont(new Font("Ubuntu", 1, 15)); // NOI18N
        t_ARTLabel.setText(" A       R       T");

        GroupLayout optimizeARTPanelLayout = new GroupLayout(optimizeARTPanel);
        optimizeARTPanel.setLayout(optimizeARTPanelLayout);
        optimizeARTPanelLayout.setHorizontalGroup(
            optimizeARTPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(optimizeARTPanelLayout.createSequentialGroup()
                .addGap(18, 18, 18)
                .addGroup(optimizeARTPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(optimLabel, GroupLayout.PREFERRED_SIZE, 126, GroupLayout.PREFERRED_SIZE)
                    .addGroup(optimizeARTPanelLayout.createSequentialGroup()
                        .addGap(6, 6, 6)
                        .addGroup(optimizeARTPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(optimizeARTPanelLayout.createSequentialGroup()
                                .addComponent(accSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(relSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(timeSlider, GroupLayout.PREFERRED_SIZE, 38, GroupLayout.PREFERRED_SIZE))
                            .addComponent(optimizeARTFilter)
                            .addComponent(t_ARTLabel, GroupLayout.PREFERRED_SIZE, 114, GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(35, Short.MAX_VALUE))
        );
        optimizeARTPanelLayout.setVerticalGroup(
            optimizeARTPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, optimizeARTPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(optimLabel, GroupLayout.PREFERRED_SIZE, 15, GroupLayout.PREFERRED_SIZE)
                .addGap(22, 22, 22)
                .addComponent(t_ARTLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 20, Short.MAX_VALUE)
                .addGroup(optimizeARTPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(accSlider, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 92, GroupLayout.PREFERRED_SIZE)
                    .addComponent(relSlider, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 92, GroupLayout.PREFERRED_SIZE)
                    .addComponent(timeSlider, GroupLayout.Alignment.TRAILING, GroupLayout.PREFERRED_SIZE, 92, GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(optimizeARTFilter)
                .addGap(54, 54, 54))
        );

        GroupLayout layout = new GroupLayout(diagnosticPanel);
        diagnosticPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(freqdiagPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(timediagPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(optimizeARTPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                )
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addComponent(optimizeARTPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(freqdiagPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(timediagPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );
    }

    private void optimizeARTFilterNow(ActionEvent evt) {
        computeOptimalFilter();
    }

    private void accSliderStateChanged(ChangeEvent evt) {
        JSlider s = (JSlider)evt.getSource();
        w_1 = (double)s.getValue()*1.0;
    }

    private void relSliderStateChanged(ChangeEvent evt) {
       JSlider s = (JSlider)evt.getSource();
       w_2 = (double)s.getValue()*1.0;
    }

    private void timeSliderStateChanged(ChangeEvent evt) {
        JSlider s = (JSlider)evt.getSource();
        w_3 = (double)s.getValue()*1.0;
    }

    private void updateDiagnosticPanel()
    {
      d_mseText.setText(df.format(mdfa.MDFAmin)); 
      d_icText.setText(df.format(Math.abs(mdfa.criteria))); 
      d_dofText.setText(df.format(Math.abs(mdfa.degrees)));
        
      d_accText.setText(df.format(Math.log(mdfa.diff_band))); 
      d_relText.setText(df.format(Math.log(mdfa.diff_stop))); 
      d_timeText.setText(df.format(mdfa.tdelay));
           
      if(timeDiagnostics)
      {
        t_mseText.setText(df.format(score[6]));  
        t_corrText.setText(df.format(score[4])); 
        t_excessText.setText(df.format(score[5]));
 
        t_accText.setText(df.format(score[0]));
        t_relText.setText(df.format(score[1])); 
        t_timeText.setText(df.format(score[2]));
      }
      else 
      {
        t_mseText.setText("N/A");  
        t_corrText.setText("N/A"); 
        t_excessText.setText("N/A");
 
        t_accText.setText("N/A");
        t_relText.setText("N/A"); 
        t_timeText.setText("N/A");
      }
    }



   /* ------------------------------------------------------------------------------
      
      The idea here is to start with a subset of observations, compute the in-sample filter, 
      then add M out-of-sample values and compute the performance of the filter out-of-sample,
      then average the statistics. 

      1) Computes the filter and constructs signal starting at start for n_insample points
      2) Applies filter on n_insamp + n_outsamp points to get signal out-of_sample
      3) Computes trading statistics based on in-sample + out_sample

      


      
      --------------------------------------------------------------------------------  */
   public void out_of_sample_tradingSweep(int start, int n_in, int n_out)
   {
      int i,j;        
      int n_insamp = n_in+(L-1); 
      int n_outsamp = n_out;
      int n_total = n_insamp+n_outsamp;
       
      if((L-1) + start + n_in + n_out <= n_obs_fixed)
      {

       //---------------------------------------------       
       // make new data matrix with in_sample data only
       double[] gdp = new double[n_rep*n_insamp];
    
       for(i=0;i<n_insamp;i++) {gdp[i] = extend[start+i];}

       int count = 1;
       for(j=0;j<series_max;j++)
       {
          if(expVariable[j].isSelected() && expVariable[j].isEnabled())
          {
            double[] temp = complete_data.get(j);  
            for(i=0;i<n_insamp;i++)
            {gdp[n_insamp*count + i] = temp[start+i];}
            count++;
          }  
       }
  
       //----- compute everything in-sample and save filter --------------
       //setTimeSeries(gdp, n_insamp, n_rep);
       mdfa.set_tseries(gdp,n_insamp,n_rep);
       mdfa.computeFilterGeneral(true,false);
       mdfa.saveBfilter(); 
       
       
       // ------------------------ now get new series with added out of sample
       gdp = new double[n_rep*n_total];
    
       for(i=0;i<n_total;i++) {gdp[i] = extend[start+i];}

       count = 1;
       for(j=0;j<series_max;j++)
       {
          if(expVariable[j].isSelected() && expVariable[j].isEnabled())
          {
            double[] temp = complete_data.get(j);  
            for(i=0;i<n_total;i++)
            {gdp[n_total*count + i] = temp[start+i];}
            count++;
          }  
       }
  
       //--------------now apply the filter out of sample---------- 
       double[] _sig = mdfa.static_computeRTSE(gdp, n_total, n_rep);  //-- applies the older filter, nothing recomputed
   
       //------------- now extract and compute the trading stuff------------
  
       int fl = n_total - (L-1);
       double[] _price = new double[fl];         
       double[] logreturns = new double[fl];
        
                   
       for(i=0;i<fl;i++)
       {
          _price[i] = t_price[0][L2+(L-1)+start+i];          
          logreturns[i] = mdfa.tseries[(L-1)+i];
       }       
 
       //-----out-of-sample stats go here------------------

       if(trading_func == 0)
       {insampleTradingDiff(_price, _sig, fl);}
       else if(trading_func == 1)
       {insampleTradingLogPrice(_price, _sig, fl);}

       postTradingStatistics();
       if(n_outsamp > 0)
       {outsampleStats(n_outsamp); printOutSampStats();}


 
       getAccount_canvas().setAccount(account);
       getAccount_canvas().setSignal(_sig);
       getAccount_canvas().setPrice(_price);
       getAccount_canvas().setLogReturns(logreturns); 
       getAccount_canvas().setBuySellLines(xf_turn_val);
       if(trading_func == 1) {getAccount_canvas().setDerivativeSignal(d_signal);}
       getAccount_canvas().go();
    
    }

   }
   
   
   
   /*-----------------------------------------------------------------------------
     The idea here is to first compute the out-of-sample signal for the given number 
     of L1 points, then compute the signal in-sample with the additional L1 points. 
     We then take the L^2 error between out-of-sample and in-sample and use it as a metric.
     
      
   ------------------------------------------------------------------------------*/
   public void computeOutOfSampleError()
   {
      int count,i,j,np,N_fixed,N;
      double diff;
      double[] gdp; double sum; 

      //System.out.println("Gets in compute error function");
      
      if(L1 > 0 && reComp.isSelected())
      {
        
         //step 1) compute out of sample signal for given parameterization
         mdfa.saveBfilter(); 
         np = L1 + (L-1) - true_out;
         gdp = new double[n_rep*np];
    
         for(i=0;i<np;i++) {gdp[i] = extend[extend.length-true_out-np+i];}

         count = 1;
         for(j=0;j<series_max;j++)
         {
          if(expVariable[j].isSelected() && expVariable[j].isEnabled())
          {
            double[] temp = complete_data.get(j);  
            for(i=0;i<np;i++)
            {gdp[np*count + i] = temp[temp.length-true_out-np+i];}
            count++;
          }  
         }       
         double[] out_sample_sig = mdfa.static_computeRTSE(gdp, np, n_rep);
   
         //step 2) compute in-sample signal for given parameterization
         
         N_fixed = n_obs_fixed; N = N_fixed;
         //System.out.println("N_fixed = " + N_fixed);
         gdp = new double[n_rep*N];
         
         for(i=0;i<N;i++)
         {gdp[i] = extend[i];}
         
         count = 1; 
         for(j=0;j<series_max;j++)
         {
           if(expVariable[j].isSelected() && expVariable[j].isEnabled())
           {
            double[] temp = complete_data.get(j);
        
            for(i=0;i<N;i++)
            {
              gdp[N*count + i] = temp[i]; 
            }
            count++;
           }
         }   
          
         int Ktemp = N_fixed/2+1; 
         tfilter.setK(Ktemp); tfilter.recomputeGamma();
   
         double[] in_sample_sig = mdfa.computeFilterSimple(gdp, N_fixed, Ktemp, tfilter.Gamma);
   
         //step 3) compute error
         double[] in = new double[L1]; 
         double[] out = new double[L1]; 
         double[] tse = new double[L1];
             
         
         sum = 0;
         for(i=0; i < L1; i++)
         {
            diff = (out_sample_sig[out_sample_sig.length-1-i] - in_sample_sig[in_sample_sig.length-1-i]);
            in[L1-1-i] = out_sample_sig[out_sample_sig.length-1-i];
            out[L1-1-i] = in_sample_sig[in_sample_sig.length-1-i];
            tse[L1-1-i] = extend[extend.length-1-i];
            sum = sum + diff*diff;
         }
         out_samp_error = 10000.0*Math.sqrt(sum)/(double)L1; 
         mdfaText.setText("Out-samp error = " + df4.format(out_samp_error));
         
         if(envision)
         {
          crystal_ball.setSeries(tse, in, out);
          crystal_ball.go();
         } 
         //--- put Gamma back to normal, just in case -----
         tfilter.setK(K); 
         tfilter.recomputeGamma();
         
      }  
   }
   

   public void computeOutOfSampleStats()
   {   
      int count,i,j,np,l; int L_price; int hf_fl = 0;
      double[] gdp; double sum; double[] indicator = new double[1];

      double[] _hf_price = null;
      true_out = true_outBar.getValue();  //System.out.println("True_out = " + true_out);  
      if(true_out > L1-10) {true_out = 0;} //System.out.println(true_out);
      
      if(L1 > 0 && reComp.isSelected()) //need out-of-sample data for this
      {

       mdfa.saveBfilter(); 
       np = L1 + (L-1) - true_out;
    
       // ------------------------ now get new series with added out of sample
       gdp = new double[n_rep*np];
    
       for(i=0;i<np;i++) {gdp[i] = extend[extend.length-true_out-np+i];}

       count = 1;
       for(j=0;j<series_max;j++)
       {
          if(expVariable[j].isSelected() && expVariable[j].isEnabled())
          {
            double[] temp = complete_data.get(j);  
            for(i=0;i<np;i++)
            {gdp[np*count + i] = temp[temp.length-true_out-np+i];}
            count++;
          }  
       }
  
       //--------------now apply the filter out of sample---------- 
       double[] _sig = mdfa.static_computeRTSE(gdp, np, n_rep);  //-- applies the older filter, nothing recomputed
   
       //compute 1-step ahead forecast error (assumes first row is target)
       double[] fore_error = forecastError(_sig, gdp, np); 
   
   
       //------------- now extract and compute the trading stuff------------
       int start;
       int fl = np - (L-1);
       double[] _price = new double[fl];
        
                   
       for(i=0;i<fl;i++)
       {_price[i] = t_price[0][t_price[0].length-true_out-fl+i];}       
 
       //-----out-of-sample stats go here------------------

        if(hf_price_set && trading_func == 4) //adjust the hf_price to match the lf price
        {
          start = (t_price[0].length-true_out-fl)*hf_period;
          hf_fl = hf_price_length - start;
          _hf_price = new double[hf_fl];
         for(i=0;i<hf_fl;i++) {_hf_price[i] = hf_price[start + i];}
        } 




       //compute price filter
       if(trading_func == 2 && mdfa.price_filter) 
       {
          indicator = new double[fl];
          L_price = mdfa.L_price;
           
          for(i = 0; i < fl; i++)
          {
           sum=0.0;
           for(l=0; l < L_price; l++)
           {
            sum = sum + mdfa.b_price_filter[l]*t_price[0][t_price[0].length-true_out-fl - l + i];
           }
           indicator[i] = sum;
          }
       }
               
       
       
       
       
       
       if(trading_func == 0)
       {insampleTradingDiff(_price, _sig, fl);}
       else if(trading_func == 1)
       {insampleTradingLogPrice(_price, _sig, fl);}  
       else if(trading_func == 2)
       {trading_function_Duplex(_price, indicator, _sig, fl);}
       else if(trading_func == 3)
       {trading_function_Duplex_Price(_price, indicator, _sig, fl);} 
       else if(trading_func == 4)
       {multiscaleFilterSimple(_price, _sig, fl, _hf_price, hf_fl, hf_period);}       

       computeSharpeRatio(250);

       outavgText.setText(""+df.format(rank_coeff));
       outmaxdropText.setText(""+df.format(max_drop));  
       outndropsText.setText(""+df.format(sign_correct));
       outratioText.setText(""+df.format((double)succ_trades/(double)total_trades));
       outroiText.setText(""+df4.format(account[account.length-1]));
       outsharpeText.setText(""+df.format(sharpeRatio));
       outsuccText.setText(""+succ_trades);
       outtotalText.setText(""+total_trades);
       outForeText.setText(""+df.format(fore_error[0]) + ", " + df.format(fore_error[1]));
       outRRText.setText(""+df.format(risk_reward));
       //System.out.println(fl);
       
      }
   }

   
   
   public double[] forecastError(double[] signal, double[] target, int np)
   {
      int i;
      double mse = 0; double sd = 0; 
      double[] msv = new double[2];
      
      System.out.println("\nPrinting Errors.........");
      for(i = 0; i < signal.length-1; i++)
      {
        mse = mse + Math.abs(target[np - 1 - i] - signal[signal.length - 2 - i]);  
        System.out.println(target[np - 1 - i] + ", " + signal[signal.length - 2 - i] + ", " + Math.abs(target[np - 1 - i] - signal[signal.length - 2 - i])); 
      }      
      mse = mse/(double)(signal.length-1.0);
      
      for(i = 0; i < signal.length-1; i++)
      {
        sd = sd + (Math.abs(target[np - 1 - i] - signal[signal.length - 2 - i]) - mse)*(Math.abs(target[np - 1 - i] - signal[signal.length - 2 - i]) - mse);   
      } 
      sd = Math.sqrt(sd)/(double)(signal.length-1.0);
    
      msv[0] = mse; msv[1] = sd;
      
      System.out.println("MSE and STD = " + msv[0] + " " + msv[1]);
      
      
      return msv;
   }
   
   
   
   
   
   public void activateSweepMode()
   {
      int i;
      
      int fl = n_obs - (L-1);
      double[] _price = new double[fl]; 
      
      getAccount_canvas().sweepMode(true);
      
      for(i=0;i<fl;i++)
      {
        _price[i] = t_price[0][L2+(L-1)+i];                            
      }
      getAccount_canvas().setTotalPrice(_price); 
      sweepMode = true;
      
      //----- Disable some shit---------
      computeFilter = false;
      mdfaPlotPane.setEnabledAt(0,!sweepMode);
      mdfaPlotPane.setEnabledAt(1,!sweepMode);
      mdfaPlotPane.setEnabledAt(2,!sweepMode);
      mdfaPlotPane.setEnabledAt(3,!sweepMode);
      mdfaPlotPane.setEnabledAt(4,!sweepMode);
      l1Bar.setEnabled(!sweepMode);
      LBar.setEnabled(!sweepMode);
      
      
   }
   
   public void deactivateSweepMode()
   {
   
      sweepMode = false;
      computeFilter = true;
      getAccount_canvas().sweepMode(sweepMode); 
      out_of_sample_tradingSweep(0, n_obs, 0);
           
      mdfaPlotPane.setEnabledAt(0,!sweepMode);
      mdfaPlotPane.setEnabledAt(1,!sweepMode);
      mdfaPlotPane.setEnabledAt(2,!sweepMode);
      mdfaPlotPane.setEnabledAt(3,!sweepMode);
      mdfaPlotPane.setEnabledAt(4,!sweepMode);
      
      l1Bar.setEnabled(!sweepMode);
      LBar.setEnabled(!sweepMode);
      
   }
    

   public void printOutSampStats()
   {
     omaxDropsText.setText(df.format(out_max_drop));
     onumDropsText.setText(df.format(out_n_drops));
     oROIText.setText(df4.format(out_ROI));
     otradeRatioText.setText(df.format(out_ratio));
   }

    public void computeSharpeRatio(int ann)
    {
      int i; double sum; double mean; double stand;
      int length = account.length;
      double[] diff_account = new double[length];
      double reward = 0; double risk = 0; int n_reward = 0; int n_risk = 0;
      
      diff_account[0] = 0.0;
      
      for(i=0;i<length-1;i++)
      {
       diff_account[i+1] = account[i+1] - account[i]; 
       
       if(diff_account[i+1] > 0) {reward = reward + diff_account[i+1]; n_reward++; }
       else if(diff_account[i+1] < 0) {risk = risk - diff_account[i+1]; n_risk++; }
      }
      
      reward = reward/(double)n_reward;
      risk = risk/(double)n_risk;
      
      sum = 0;
      for(i=0; i < length; i++)
      {sum = sum + diff_account[i];}
      mean = sum/(1.0*length); 
 
      sum = 0;
      for(i=0; i < length; i++)
      {sum = sum + (diff_account[i] - mean)*(diff_account[i] - mean);}     
      stand = Math.sqrt(sum)/(1.0*length);
      
      //System.out.println("sharpe mean = " + mean + ", std = " + stand);
      
      sharpeRatio = Math.sqrt(ann)*mean/stand; 
      risk_reward = risk/reward;
    }   


    
   public void initiateTrueOutSamplePanel()
   {
   
     trueOutSamplePanel = new JPanel(); 
     
     true_outBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,151);
     true_outText = new JTextField("0");
     true_outLabel = new JLabel("True Out-of-Sample Points");
     
     GroupLayout paramLayout = new GroupLayout(trueOutSamplePanel); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(true_outLabel).addComponent(true_outBar).addComponent(true_outText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(true_outLabel).addComponent(true_outBar).addComponent(true_outText)));
        trueOutSamplePanel.setLayout(paramLayout);    
      
     AdjustmentListener ad = new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {
            true_outBar.getValue();           
            true_outText.setText(""+true_outBar.getValue());
          }};    
     
     true_outBar.addAdjustmentListener(ad);
   }
   
   
   //------ out-of-sample sweep canvas -----------------
   public void initiateOutSamplePanel()
   {


        outsampPanel = new JPanel(); JPanel jPanel1 = new JPanel();
        omaxDropLabel = new javax.swing.JLabel();
        omaxDropsText = new javax.swing.JTextField();
        onumDropsText = new javax.swing.JTextField();
        ondropsLabel = new javax.swing.JLabel();
        oROILabel = new javax.swing.JLabel();
        oROIText = new javax.swing.JTextField();
        otradeRatioLabel = new javax.swing.JLabel();
        otradeRatioText = new javax.swing.JTextField();
        
        outsampLabel = new javax.swing.JLabel();
        outsampText = new javax.swing.JTextField();
        compSweepButton = new javax.swing.JButton();
        
        startObsLabel = new javax.swing.JLabel();
        starobsText = new javax.swing.JTextField();

        omaxDropLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        omaxDropLabel.setText("Max-Drop");

        omaxDropsText.setText("0");

        onumDropsText.setText("0");

        ondropsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        ondropsLabel.setText("Number Drops");

        oROILabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        oROILabel.setText("ROI");

        oROIText.setText("0");

        otradeRatioLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        otradeRatioLabel.setText("Trade Ratio");

        otradeRatioText.setText("0");

        outsampBar = new JScrollBar(JScrollBar.HORIZONTAL,10,1,1,40);

        outsampLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        outsampLabel.setText("Out-of-Sample Observations");
        outsampLabel.setToolTipText("Choose number of out-of-sample points.");

        outsampText.setText("10");

        compSweepButton.setText("Compute Out-of-Sample Sweep");
        compSweepButton.setToolTipText("Computes the out-of-sample trading statistics with given signal. Averages statistics and displays them in panel.");

        startobsBar = new JScrollBar(JScrollBar.HORIZONTAL,60,1,40,n_obs-60);

        startObsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        startObsLabel.setText("Starting Observation");
        startObsLabel.setToolTipText("Choose starting observation for the out-of-sample sweep.");

        starobsText.setText("60");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(31, 31, 31)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(compSweepButton)
                    .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(jPanel1Layout.createSequentialGroup()
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(ondropsLabel)
                                .addComponent(omaxDropLabel))
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(jPanel1Layout.createSequentialGroup()
                                    .addComponent(omaxDropsText, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGap(18, 18, 18)
                                    .addComponent(oROILabel))
                                .addGroup(jPanel1Layout.createSequentialGroup()
                                    .addComponent(onumDropsText, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGap(18, 18, 18)
                                    .addComponent(otradeRatioLabel))))
                        .addGroup(jPanel1Layout.createSequentialGroup()
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(startObsLabel)
                                .addComponent(outsampLabel))
                            .addGap(18, 18, 18)
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(outsampBar, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(startobsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(outsampText, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(otradeRatioText, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(starobsText, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(oROIText, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(76, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(48, 48, 48)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(omaxDropLabel)
                    .addComponent(omaxDropsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(oROILabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(oROIText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(onumDropsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(otradeRatioLabel)
                        .addComponent(otradeRatioText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(ondropsLabel))
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(38, 38, 38)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(outsampBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(outsampLabel))
                        .addGap(18, 18, 18)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(startobsBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(startObsLabel)))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(32, 32, 32)
                        .addComponent(outsampText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(starobsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(21, 21, 21)
                .addComponent(compSweepButton)
                .addGap(31, 31, 31))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(outsampPanel);
        outsampPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 480, Short.MAX_VALUE)
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(0, 21, Short.MAX_VALUE)
                    .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGap(0, 21, Short.MAX_VALUE)))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 283, Short.MAX_VALUE)
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(0, 0, Short.MAX_VALUE)
                    .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGap(0, 0, Short.MAX_VALUE)))
        );


         AdjustmentListener ad = new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {
            int nout = outsampBar.getValue(); 
            getAccount_canvas().setNOutSample(nout);
            out_of_sample_tradingSweep(getAccount_canvas().nbackObs, getAccount_canvas().n_insamp, nout);
            outsampText.setText(""+outsampBar.getValue());
          }};

         AdjustmentListener ad2 = new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {            
            starobsText.setText(""+startobsBar.getValue());
          }};

         ActionListener alc = new ActionListener() {
         public void actionPerformed(ActionEvent event) 
         {
           sweepOutSample();        
         }};

         outsampBar.addAdjustmentListener(ad);
         startobsBar.addAdjustmentListener(ad2);
         compSweepButton.addActionListener(alc);     

    }



    public void sweepOutSample()
    {
      int startObs = startobsBar.getValue();

      if(n_obs > startObs)
      {

       outSamp_means = new double[4];
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
       int startObs = startobsBar.getValue();
       int n_out = outsampBar.getValue();
         
       for(i=startObs; i < n_obs-n_out; i++)
       {
         
          out_of_sample_tradingSweep(0, i, n_out);

          outSamp_means[0] = outSamp_means[0] + out_max_drop;
          outSamp_means[1] = outSamp_means[1] + out_n_drops;
          outSamp_means[2] = outSamp_means[2] + out_ROI;
          outSamp_means[3] = outSamp_means[3] + out_ratio;
          getAccount_canvas().go();

          try{Thread.sleep(50);}
          catch(Exception e){}  
       }
       
       //-------- Now compute statistics ---------------------------------------------
       
       for(j=0;j<4;j++)
       {outSamp_means[j] = outSamp_means[j]/(n_obs-n_out-startObs-1);}  

       omaxDropsText.setText(df.format(outSamp_means[0]));
       onumDropsText.setText(df.format(outSamp_means[1]));
       oROIText.setText(df.format(outSamp_means[2]));
       otradeRatioText.setText(df.format(outSamp_means[3]));

       //-------- Done with thread---------------------------------------      
       return;
    }   
  } 





















   //------ account canvas -----------------------

   public class accountCanvas extends JPanel 
   {  

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int trade_obs;
    int pressed_t;
    int total_obs;

    int nbackObs = 0;
    int nforeObs = 0;
    int nObsSpan;  

    int n_insamp, n_outsamp;

    int slideStart = 0; 
    int slideEnd;    
    Cursor curCursor;
    double tradeDelta = 0.0;
    
    double[] account;
    double[] signal;
    double[] logprice;
    double[] logreturn;
    double[] buysell;   
    double[] totalprice; //sweepMode only
    double[] d_signal; //derivative of signal
    double[] filtered_price;
    double[] hf_price;    

    boolean plot_hfprice = false;
    boolean sweepMode = false;
    boolean plot_account=true;
    boolean plot_signal=false;
    boolean plot_logprice=false;
    boolean plot_logreturn=false;
    boolean accountSet = false;
    boolean logreturnSet = false;
    boolean logpriceSet = false;
    boolean buysellSet = false;
    boolean signalSet = false;
    boolean plot_lines = false;
    boolean plot_buy = false;
    boolean plot_indicator = false;
    boolean filtered_priceSet = false;
    boolean shift_mean_signal = false;
    
    public Font mono;
    public Cursor curCursor2;
    public Color fillColor; 
    public String value;
    public DecimalFormat df;
    
    boolean plot_tracker = false;
    int tracker = 0;  

    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;
    Ellipse2D ellipse;
    Rectangle2D rectangle;
    
    int track_pos_t,track_pos_x;
    int minObs = 60;
    
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2,myGray4; 
    Color plotGray;

    double totalpriceMin,totalpriceMax,totalpriceNorm;
    double priceMin,priceMax,priceNorm;
    double returnMin, returnMax, returnNorm;  
    double dsigMax, dsigMin, dsigNorm;   
    double sigMax,sigMin, sigNorm;
       
    int green = 158; int blue = 224;
    Color dLines;
    Color dLinesLoss;
    Color dLinesBuy;
    Color dLinesSell;
    Color myGray3;


    public accountCanvas(int w, int h, int _L, int _nrep)
    {
      trade_obs = 300-_L;
      
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      setBackground(Color.BLACK);
      setPreferredSize(new Dimension(w, h));    
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
      myGray = new Color(107,101,133);
      myGray2 = new Color(60,61,69);
      myGray3 = new Color(157,151,183);
      myGray4 = new Color(150,101,153);
      plotGray = new Color(92,172,238);
      df = new DecimalFormat("##.###"); 
      dLines = new Color(0,0,140);
      dLinesLoss = new Color(140,0,0);
      dLinesBuy = new Color(0,140,0);
      dLinesSell = new Color(140,0,140);
      value = new String("");
      n_outsamp = 10;

      addMouseMotionListener(new MyArtMouseMotionListener());  
      addMouseListener(new MyArtMouseListener());
      nforeObs = 0; nbackObs = 0; slideEnd = 0; slideStart = 0;
      curCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);  

      mono  = new Font("Monospaced", Font.PLAIN, 12);
      addMouseMotionListener(new MyArtMouseMotionListener());    
    
    } 

    public void shiftPrice(boolean t)
    {shift_mean_signal = t;}
    
    public void setPlotDim()  //only considers dimensions of account 
    {
        int i;
        dataMax = -100000.0; dataMin = 100000.0; 
        priceMax = -100000.0; priceMin = 100000.0; 
        returnMax = -100000.0; returnMin = 100000.0; 
        dsigMax = -10000.0; dsigMin = 100000.0; sigMax = -10000.0; sigMin = 100000.0;
        totalpriceMax = -10000.0; totalpriceMin = 100000.0;
 
        for(i=0;i<trade_obs;i++)
        {
             if(account[i] < dataMin) dataMin = account[i];
             else if(account[i] > dataMax) dataMax = account[i]; 

             if(logprice[i] < priceMin) priceMin = logprice[i];
             else if(logprice[i] > priceMax) priceMax = logprice[i]; 

             if(logreturn[i] < returnMin) returnMin = logreturn[i];
             else if(logreturn[i] > returnMax) returnMax = logreturn[i]; 
             
             if(trading_func == 1)
             { 
              if(d_signal[i] < dsigMin) dsigMin = d_signal[i];
              else if(d_signal[i] > dsigMax) dsigMax = d_signal[i];
             }
             
             if(trading_func == 1 && shift_mean_signal)
             {
                 if(signal[i] < sigMin) sigMin = signal[i];
                 else if(signal[i] > sigMax) sigMax = signal[i];
             }
                          
        }           
        dataNorm = Math.abs(dataMax - dataMin);
        priceNorm = Math.abs(priceMax - priceMin);
        returnNorm = Math.abs(returnMax - returnMin);
        dsigNorm = Math.abs(dsigMax - dsigMin);
        if(shift_mean_signal) {sigNorm = Math.abs(sigMax - sigMin);}
        
        if(!sweepMode) {totalpriceNorm = priceNorm; totalpriceMax = priceMax; totalpriceMin = priceMin;}
        else if(sweepMode) 
        {
           //System.out.println("total_obs = " + total_obs);
           for(i=0;i<totalprice.length;i++)
           {
             if(totalprice[i] < totalpriceMin) totalpriceMin = totalprice[i];
             else if(totalprice[i] > totalpriceMax) totalpriceMax = totalprice[i];
           }
           nObsSpan = trade_obs - n_outsamp; totalpriceNorm = Math.abs(totalpriceMax - totalpriceMin);
           //System.out.println("nObsSpan = " + nObsSpan);
        }
        
        if(!shift_mean_signal)
        {sigMin = returnMin; sigMax = returnMax; sigNorm = returnNorm;}
        
    }
 
   public void setTradeThreshold(double t) {tradeDelta = .001*t; go();}


   public void setTotalPrice(double[] _p)
   {
     if(sweepMode)
     {
       total_obs = _p.length;
       totalprice = new double[total_obs];
       System.arraycopy(_p, 0, totalprice, 0, total_obs);
       go(); 
     }
   
   }
   
   public void sweepMode(boolean f)
   {
     sweepMode = f; if(!sweepMode) {go();}
   }
  
 
    public void setPlotTracker(int t) 
    {
       if(t >= 0)
       {
         plot_tracker = true;
         tracker = t;
         go();
       }
       else 
       {plot_tracker = false; go();}    
    }
 

    public void setAccount(double[] acc)
    {
      trade_obs = acc.length;
      account = new double[trade_obs];
      System.arraycopy(acc, 0, account, 0, trade_obs);
      accountSet = true; 
    }
    


    public void setSignal(double[] acc)
    {
      if(trade_obs != acc.length)
      {System.out.println("Dimensions of account and signal must be the same");}
      else
      {
       signal = new double[trade_obs];
       System.arraycopy(acc, 0, signal, 0, trade_obs);
       signalSet = true; 
      }
    }

    public void setDerivativeSignal(double[] acc)
    {
      if(trade_obs != acc.length)
      {System.out.println("Dimensions of account and signal must be the same, trade_obs = " + trade_obs + " " + acc.length);}
      else
      {
       d_signal = new double[trade_obs];
       System.arraycopy(acc, 0, d_signal, 0, trade_obs);
       signalSet = true; 
      }
    }

    public void setPrice(double[] acc)
    {
      if(trade_obs != acc.length)
      {System.out.println("Dimensions of account and logprice must be the same");}
      else
      {
       logprice = new double[trade_obs];
       System.arraycopy(acc, 0, logprice, 0, trade_obs);
       logpriceSet = true; 
      }
    }

    public void setHFPrice(double[] hfacc)
    {
      //System.out.println("set's the hf price in accountCanvas");
      hf_price = new double[hfacc.length];
      System.arraycopy(hfacc,0,hf_price,0,hf_price.length);
      plot_hfprice = true;
    }
    
    public void setFilteredPrice(double[] pr)
    {
      if(trade_obs != pr.length)
      {System.out.println("Dimensions of account and logprice must be the same");}
      else
      {
       filtered_price = new double[trade_obs];
       System.arraycopy(pr, 0, filtered_price, 0, trade_obs);
       filtered_priceSet = true; 
      }         
    }

    public void setNOutSample(int n) {n_outsamp = n;}
    
    public void setLogReturns(double[] acc)
    {
      if(trade_obs != acc.length)
      {System.out.println("Dimensions of account and logreturns must be the same");}
      else
      {
       logreturn = new double[trade_obs];
       System.arraycopy(acc, 0, logreturn, 0, trade_obs);
       logreturnSet = true; 
      }
    }

    public void setBuySellLines(double[] acc)
    {
      if(trade_obs != acc.length)
      {System.out.println("Dimensions of account and logreturns must be the same");}
      else
      {
       buysell = new double[trade_obs];
       System.arraycopy(acc, 0, buysell, 0, trade_obs);
       buysellSet = true; 
      }
    }

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,w,h,obs;
     int t0, t1, x0, x1;
     
     if(accountSet)
     {
      setPlotDim();
     }
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
    
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
     BasicStroke def = new BasicStroke((float)1.9);    
     Color myPurple = new Color(179, 34, 250);
     //height = (int)height*8/10;

     new BasicStroke((float)1.7); 
     Stroke thin = new BasicStroke((float)1.2);   
   
     g2d.setStroke(thin);
        
     //System.out.println(accountSet + "  " + plot_logprice + "  " + plot_logreturn + "  " + plot_signal);

    
     if(!sweepMode)     
     {obs = trade_obs;}
     else
     {
       obs = total_obs;
       //---- plot total price -----

       g2d.setStroke(def);
       t0 = (int)(((double)(nbackObs)/(double)total_obs)*(double)width);
       t1 = (int)(((double)(nbackObs+nObsSpan-1)/(double)total_obs)*(double)width);
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
     

       g2d.setPaint(myGray2); 
       for(j = 0; j < obs-1; j++)
       {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)obs)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)obs)*(double)width);
	  x0 = (int)(((totalprice[j] - totalpriceMin)/totalpriceNorm)*(double)height);
	  x1 = (int)(((totalprice[j + 1] - totalpriceMin)/totalpriceNorm)*(double)height);
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
       }   
         
     }
     
     g2d.setPaint(myGray3);
     g.drawString((String)df.format(dataMax), 5, 15);
     g.drawString((String)df.format(dataMin), 5, height - 24);    
     
    
      if(plot_logprice)
      {     
        g2d.setStroke(thin);
        g2d.setPaint(myGray); 
        for(j = 0; j < trade_obs-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	  t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	  x0 = (int)(((logprice[j] - totalpriceMin)/totalpriceNorm)*(double)height);
	  x1 = (int)(((logprice[j + 1] - totalpriceMin)/totalpriceNorm)*(double)height);
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
        }      
       
        if(plot_hfprice)
        {
          g2d.setPaint(myGray4);
          for(j = 0; j < hf_price.length-1; j++)
          {
           //System.out.println(gamma_hat[N*k+j]);
	   t0 = (int)(((double)(nbackObs+j)/(double)hf_price.length)*(double)width);
	   t1 = (int)(((double)(nbackObs+j+1)/(double)hf_price.length)*(double)width);
	   x0 = (int)(((hf_price[j] - totalpriceMin)/totalpriceNorm)*(double)height);
	   x1 = (int)(((hf_price[j + 1] - totalpriceMin)/totalpriceNorm)*(double)height);
	   g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
	   //System.out.println(hf_price[j]);
          }
        }

      }
      if(plot_indicator)
      {
      
        g2d.setStroke(thin);
        g2d.setPaint(myPurple); 
        for(j = 0; j < trade_obs-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	  t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	  x0 = (int)(((filtered_price[j] - totalpriceMin)/totalpriceNorm)*(double)height);
	  x1 = (int)(((filtered_price[j + 1] - totalpriceMin)/totalpriceNorm)*(double)height);
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
        }    
     
      }
      if(plot_logreturn)
      {     
        g2d.setPaint(new Color(137,116,164));
        for(j = 0; j < trade_obs-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	  t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	  x0 = (int)(((logreturn[j] - returnMin)/returnNorm)*(double)height);
	  x1 = (int)(((logreturn[j + 1] - returnMin)/returnNorm)*(double)height);
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
        }                   
     }
     if(plot_signal)
     {     

        g2d.setPaint(new Color(16,204,50));
        for(j = 0; j < trade_obs-1; j++)
        {
          //System.out.println(signal[j]);
	  t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	  t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	  x0 = (int)(((signal[j] - sigMin)/sigNorm)*(double)height);
	  x1 = (int)(((signal[j + 1] - sigMin)/sigNorm)*(double)height);
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
        } 

        g2d.setPaint(new Color(90,90,100));
        g2d.setStroke(dashed);
        x0 = (int)(((0.0 - sigMin)/sigNorm)*(double)height);
        g2d.drawLine(0, height - x0, width, height - x0);

        if(tradeDelta > 0)
        {
           g2d.setPaint(new Color(160,90,160));
           g2d.setStroke(dashed);
           x0 = (int)(((tradeDelta - returnMin)/returnNorm)*(double)height);
           g2d.drawLine(0, height - x0, width, height - x0);           

           g2d.setPaint(new Color(160,90,160));
           g2d.setStroke(dashed);
           x0 = (int)(((-tradeDelta - returnMin)/returnNorm)*(double)height);
           g2d.drawLine(0, height - x0, width, height - x0); 
        }
        
//         for(j = 0; j < trade_obs-1; j++)
//         {
//          if(account[j+1] != account[j])
//          {System.out.println(df4.format(signal[j]) + " " + df4.format(account[j]) + ", change signs here");}
//          else
//          {System.out.println(df4.format(signal[j]) + " " + df4.format(account[j]));}
//         } 
     }
     if(plot_lines)
     { 
           g2d.setStroke(thin);
           g2d.setPaint(new Color(10,180,150));
           for(j = 0; j < trade_obs-1; j++)
           {
             //System.out.println(signal[j]);
	     t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	     t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	     x0 = (int)(((d_signal[j] - dsigMin)/dsigNorm)*(double)height);
	     x1 = (int)(((d_signal[j + 1] - dsigMin)/dsigNorm)*(double)height);
	     g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
           } 

           g2d.setPaint(new Color(90,90,100));
           g2d.setStroke(dashed);
           x0 = (int)(((0.0 - dsigMin)/dsigNorm)*(double)height);
           g2d.drawLine(0, height - x0, width, height - x0);
           
           g2d.setPaint(new Color(160,90,160));
           g2d.setStroke(dashed);
           x0 = (int)(((tradeDelta - dsigMin)/dsigNorm)*(double)height);
           g2d.drawLine(0, height - x0, width, height - x0);           

           g2d.setPaint(new Color(160,90,160));
           g2d.setStroke(dashed);
           x0 = (int)(((-tradeDelta - dsigMin)/dsigNorm)*(double)height);
           g2d.drawLine(0, height - x0, width, height - x0);            

     } 

     if(accountSet && plot_account)
     {     
        g2d.setStroke(thin);
        for(j = 0; j < trade_obs-1; j++)
        { 
          g2d.setStroke(thin);
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	  t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
	  x0 = (int)(((account[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((account[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.setPaint(new Color((int)((x1/(double)height)*255), green,blue));
	  g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          //System.out.println((int)(x1/(double)height)*255);

          //---- Draw derivative of signal-----

   
        }                    
     } 

     if(plot_buy)
     {

         g2d.setStroke(thin); g2d.setStroke(dashed);
         for(j = 0; j < trade_obs-1; j++)
         { 
           //System.out.println(gamma_hat[N*k+j]);
	   t0 = (int)(((double)(nbackObs+j)/(double)obs)*(double)width);
	   t1 = (int)(((double)(nbackObs+j+1)/(double)obs)*(double)width);
           if(buysell[j+1] >= 1.0)     
           {g2d.setPaint(dLinesBuy); g2d.drawLine(t1, 0, t1, height);}
           else if(buysell[j+1] <= -1.0)
           {g2d.setPaint(dLinesSell); g2d.drawLine(t1, 0, t1, height);}
         }  
     }

     if(plot_tracker)
     {
        w = 4; h = 4; 
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t-2, (height-2) - track_pos_x, w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);

        
        g.setFont(mono);
        g2d.setPaint(Color.GREEN);
        g.drawString(value, 5, 15);           
        //setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));    
     }
     
    
 
   }

 

  
   
     class MyArtMouseMotionListener extends MouseMotionAdapter 
     {
      
      public void mouseDragged(MouseEvent e) 
      {

        if(sweepMode)
        {

         boolean passx = false;
         
         int t1 = (int)(((double)total_obs)*e.getX()/(double)width); 
         int t0 = (int)(((double)total_obs)*pressed_t/(double)width);  
                  
         //if((e.getX() > slideStart+2) && (e.getX() < slideEnd-2))
         if(curCursor2 == Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR))
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
           if(passx)
           {
            n_insamp = total_obs-nforeObs-nbackObs-n_outsamp;  
            out_of_sample_tradingSweep(nbackObs, n_insamp, n_outsamp);
            go();
           }
         }
         else if(curCursor2 == Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR))
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

           if(passx)
           {             
             n_insamp = total_obs-nforeObs-nbackObs-n_outsamp;                   
             out_of_sample_tradingSweep(nbackObs, n_insamp, n_outsamp);
           }
           go();          
         }
         else if(curCursor2 == Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR))
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
           {           
             n_insamp = total_obs-nforeObs-nbackObs-n_outsamp; 
             out_of_sample_tradingSweep(nbackObs, n_insamp, n_outsamp);
           }           
           go();
         }         
        }      
      }


      public void mouseMoved(MouseEvent e) 
      {

        int j; int t0,t1;
        if(plot_tracker && account != null)
        {
            
           j = (int)(((double)trade_obs)*e.getX()/(double)width);  
        
           track_pos_t = (int)(((double)j/(double)trade_obs)*(double)width);
           //t1 = (int)(((double)(j+1)/(double)trade_obs)*(double)width);
           if(tracker == 0 && plot_account) //account
           { 
             track_pos_x = (int)(((account[j] - dataMin)/dataNorm)*(double)height); value = df.format(account[j]);
             //x1 = (int)(((account[j+1] - dataMin)/dataNorm)*(double)height);
           }
           else if(tracker == 1 && plot_logprice) //price
           {
             track_pos_x = (int)(((logprice[j] - priceMin)/priceNorm)*(double)height); value = df.format(logprice[j]);         
             //x1 = (int)(((logprice[j+1] - priceMin)/priceNorm)*(double)height);   
           }       
           else if(tracker == 2 && plot_signal) //signal
           {
            track_pos_x = (int)(((signal[j] - returnMin)/returnNorm)*(double)height); value = ""+signal[j];
            //x1 = (int)(((signal[j+1] - returnMin)/returnNorm)*(double)height);
           }
           else if(tracker == 3 && plot_logreturn) //diff
           {
             track_pos_x = (int)(((logreturn[j] - returnMin)/returnNorm)*(double)height);  value = df.format(logreturn[j]);     
             //x1 = (int)(((logreturn[j+1] - returnMin)/returnNorm)*(double)height);       
           }
           //System.out.println(j + " " + track_pos_t + "  " + track_pos_x);
           go();           
        }    
        
       if(sweepMode)
       {
        t0 = (int)(((double)(nbackObs)/(double)total_obs)*(double)width);
        t1 = (int)(((double)(nbackObs+nObsSpan-1)/(double)total_obs)*(double)width);
        if(e.getX() > t0+2 && e.getX() < t1-2) 
        { 
           curCursor2 = Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR);
           setCursor(curCursor2);                        
        }
        else if(e.getX() >= t0-2 && e.getX() <= t0+2)
        {
           curCursor2 = Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR);
           setCursor(curCursor2);
        }
        else if(e.getX() >= t1-2 && e.getX() <= t1+2)
        {
           curCursor2 = Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR);
           setCursor(curCursor2);
        }
        else
        {
           curCursor2 = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
           setCursor(curCursor2);
        } 
        go();
       }
        
        
        
        
      }
      
     }   
   

     class MyArtMouseListener extends MouseAdapter 
     {       
       public void mousePressed(MouseEvent e) 
       { pressed_t = e.getX();}

       public void mouseReleased(MouseEvent e) { }
       public void mouseClicked(MouseEvent e) { }
        
     }

   }















 



  public void initiateAdaptiveUpdatePanel()
  {


        adaptiveUpdatePanel = new JPanel();
        continuousUpdate = new JCheckBox("Update:"); 
        regPanel2 = new JPanel();
        d_updateLabel = new JLabel();
        d_updateBar = new JScrollBar();
        d_updateText = new JTextField();
        d2_updateText = new JTextField();
        d2_updateLabel = new JLabel();
        d2_updateBar = new JScrollBar();
        s_updateText = new JTextField();
        s_updateLabel2 = new JLabel();
        s_updateBar = new JScrollBar();
        c_updateBar = new JScrollBar();
        c_updateText = new JTextField();
        c_updateLabel = new JLabel();
        lam_updateBar = new JScrollBar();
        al_updateBar = new JScrollBar();
        al_updateText = new JTextField();
        d2_updateLabel1 = new JLabel();
        d_updateLabel1 = new JLabel();
        lam_updateText = new JTextField();
        n_updateLabel = new JLabel();
        n_updateText = new JTextField();
        n_updateBar = new JScrollBar();
        l_updateBar = new JScrollBar();
        l_updateText = new JTextField();
        l_updateLabel = new JLabel();
        computeUpdate = new JButton();
        autoUpdate = new JCheckBox();
        shadeBox = new JCheckBox();
        uni_updateCheck = new JCheckBox();
        plotUpdateBox = new JCheckBox();
        updatei1Box = new JCheckBox();
        updatei2Box = new JCheckBox();

        regPanel2.setBorder(BorderFactory.createTitledBorder("Regularization"));

        d_updateLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d_updateLabel.setText("decay");

        d_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        d_updateBar.setMinimum(0); d_updateBar.setMaximum(1000); d_updateBar.setValue(0);
        d_updateBar.setUnitIncrement(1);


        d_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d_updateText.setText("0");

        d2_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d2_updateText.setText("0");

        d2_updateLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d2_updateLabel.setText("decay2");

        d2_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        d2_updateBar.setMinimum(0); d2_updateBar.setMaximum(1000); d2_updateBar.setValue(0);


        s_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        s_updateText.setText("0");

        s_updateLabel2.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        s_updateLabel2.setText("smooth");

        s_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        s_updateBar.setMinimum(0); s_updateBar.setMaximum(1000); s_updateBar.setValue(0);

        c_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        c_updateBar.setMinimum(0); c_updateBar.setMaximum(1000); c_updateBar.setValue(0);

        c_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        c_updateText.setText("0");

        c_updateLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        c_updateLabel.setText("cross");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(regPanel2);
        regPanel2.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(c_updateLabel)
                    .addComponent(s_updateLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(c_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(s_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(c_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(s_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(d_updateLabel)
                    .addComponent(d2_updateLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(d2_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(d_updateBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(d_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(d2_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap(14, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(d_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED))
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(d_updateLabel)
                                    .addComponent(d_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(13, 13, 13)))
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(d2_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(d2_updateLabel)
                                .addComponent(d2_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(s_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(s_updateLabel2)
                                .addComponent(s_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGap(16, 16, 16)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addComponent(c_updateLabel)
                                .addGap(3, 3, 3))
                            .addComponent(c_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(c_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
        );

        lam_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        al_updateBar.setOrientation(JScrollBar.HORIZONTAL);
 
        
        lam_updateBar.setMinimum(0); lam_updateBar.setMaximum(4000); lam_updateBar.setUnitIncrement(5); lam_updateBar.setValue(0); 
        al_updateBar.setMinimum(0); al_updateBar.setMaximum(400); al_updateBar.setUnitIncrement(5); al_updateBar.setValue(0); 

        al_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        al_updateText.setText("0");

        d2_updateLabel1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d2_updateLabel1.setText("\u03C9");

        d_updateLabel1.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        d_updateLabel1.setText("\u03BB");

        lam_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lam_updateText.setText("0");

        n_updateLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        n_updateLabel.setText("Obs");

        n_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        n_updateText.setText("20");

        n_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        n_updateBar.setMinimum(10); n_updateBar.setMaximum(n_obs); 
        n_updateBar.setValue(20); n_updateBar.setUnitIncrement(1);
       
        l_updateBar.setOrientation(JScrollBar.HORIZONTAL);
        l_updateBar.setMinimum(5); l_updateBar.setMaximum(n_obs-L); 
        l_updateBar.setValue(10); l_updateBar.setUnitIncrement(1);

        l_updateText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        l_updateText.setText("10");

        l_updateLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        l_updateLabel.setText("L");

        computeUpdate.setText("Adaptive Update");

        autoUpdate.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        autoUpdate.setText("Auto Update");
        autoUpdate.setToolTipText("Turn on/off automatic updating when parameters change.");

        uni_updateCheck.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        uni_updateCheck.setText("Univarate");
        uni_updateCheck.setSelected(false);
      
        
        shadeBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        shadeBox.setText("Shade Region");
        shadeBox.setToolTipText("Shade the region in which adaptive update is computed.");

        plotUpdateBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        plotUpdateBox.setText("Plot Update");
        plotUpdateBox.setToolTipText("Plot the updated signal");
        plotUpdateBox.setSelected(true);

        
        updatei1Box.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        updatei1Box.setText("i1");
        updatei1Box.setToolTipText("Set i1 constraint");
        updatei1Box.setSelected(false);

        updatei2Box.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        updatei2Box.setText("i2");
        updatei2Box.setToolTipText("Set i2 constraint");
        updatei2Box.setSelected(false);
        
        continuousUpdate.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        continuousUpdate.setText("Adaptive Update");
        continuousUpdate.setToolTipText("Continuous updating with new data");
        continuousUpdate.setSelected(false);

 
        updateCoeffPanel = new FilterCoefCanvas(600, 280, 10, 2); 
        updateCoeffPanel.setUpdateCoeff(true);
        
        updateAmpPanel = new AmpUpdateCanvas(600,280,50,n_rep);
        
        //updateCoeffPanel.setBackground(new java.awt.Color(1, 1, 1));
        updateCoeffPanel.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.LOWERED));
        updateAmpPanel.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.LOWERED));
        
        javax.swing.GroupLayout updateCoeffPanelLayout = new javax.swing.GroupLayout(updateCoeffPanel);
        updateCoeffPanel.setLayout(updateCoeffPanelLayout);
        updateCoeffPanelLayout.setHorizontalGroup(
            updateCoeffPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 496, Short.MAX_VALUE)
        );
        updateCoeffPanelLayout.setVerticalGroup(
            updateCoeffPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        
        javax.swing.GroupLayout updateAmpPanelLayout = new javax.swing.GroupLayout(updateAmpPanel);
        updateAmpPanel.setLayout( updateAmpPanelLayout);
         updateAmpPanelLayout.setHorizontalGroup(
             updateAmpPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 496, Short.MAX_VALUE)
        );
         updateAmpPanelLayout.setVerticalGroup(
             updateAmpPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        
        


        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(adaptiveUpdatePanel);
        adaptiveUpdatePanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(34, 34, 34)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(n_updateLabel)
                            .addComponent(l_updateLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(n_updateBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(l_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(n_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(l_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(d_updateLabel1)
                            .addComponent(d2_updateLabel1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(lam_updateBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(al_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(lam_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(al_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(updatei1Box, javax.swing.GroupLayout.PREFERRED_SIZE, 60, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(updatei2Box, javax.swing.GroupLayout.PREFERRED_SIZE, 60, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(continuousUpdate)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(autoUpdate)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(uni_updateCheck)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(plotUpdateBox))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(regPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(18, 18, 18)
                .addComponent(updateCoeffPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addComponent(updateAmpPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(31, 31, 31)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(n_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(47, 47, 47))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(n_updateLabel)
                                    .addComponent(n_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(22, 22, 22)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(l_updateLabel)
                                    .addComponent(l_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(l_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(d_updateLabel1)
                                    .addComponent(lam_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(lam_updateText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(35, 35, 35)
                                .addComponent(al_updateBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(al_updateText, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(d2_updateLabel1, javax.swing.GroupLayout.Alignment.LEADING))
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(updatei1Box, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(updatei2Box, javax.swing.GroupLayout.Alignment.LEADING)) 
                                )
                        .addGap(31, 31, 31)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(continuousUpdate)
                            .addComponent(autoUpdate)
                            .addComponent(uni_updateCheck)
                            .addComponent(plotUpdateBox))
                        .addGap(18, 18, 18)
                        .addComponent(regPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        //.addContainerGap()
                        .addComponent(updateCoeffPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(updateAmpPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap())
        );



        AdjustmentListener al5 = new AdjustmentListener()  {
          public void adjustmentValueChanged(AdjustmentEvent e) 
          {
        
             if(e.getAdjustable() == n_updateBar)
             {n_updateText.setText(""+n_updateBar.getValue());}
             else if(e.getAdjustable() == l_updateBar)
             {l_updateText.setText(""+l_updateBar.getValue());}           
             else if(e.getAdjustable() == lam_updateBar)
             {lam_updateText.setText(df.format(.1*lam_updateBar.getValue()));} 
             else if(e.getAdjustable() == al_updateBar)
             {al_updateText.setText(df.format(.1*al_updateBar.getValue()));} 
             else if(e.getAdjustable() == s_updateBar)
             {s_updateText.setText(df2.format(.001*s_updateBar.getValue()));} 
             else if(e.getAdjustable() == d_updateBar)
             {d_updateText.setText(df2.format(.001*d_updateBar.getValue()));}     
             else if(e.getAdjustable() == d2_updateBar)
             {d2_updateText.setText(df2.format(.001*d2_updateBar.getValue()));} 
             else if(e.getAdjustable() == c_updateBar)
             {c_updateText.setText(df2.format(.001*c_updateBar.getValue()));}    

             if(autoUpdate.isSelected())
             {updateSignal();}
             
             mdfa_canvas.shadeRegion(shadeBox.isSelected(),n_updateBar.getValue() - l_updateBar.getValue()+1);

          }
        };

        n_updateBar.addAdjustmentListener(al5); d_updateBar.addAdjustmentListener(al5);
        l_updateBar.addAdjustmentListener(al5); s_updateBar.addAdjustmentListener(al5);
        lam_updateBar.addAdjustmentListener(al5); al_updateBar.addAdjustmentListener(al5);
        c_updateBar.addAdjustmentListener(al5); d2_updateBar.addAdjustmentListener(al5);

        plotUpdateBox.addItemListener(new ItemListener() {
          public void itemStateChanged(ItemEvent e)
          {
            boolean sel; e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED) {sel = false;}       
            else {sel = true;}
             
            mdfa.swapUpdateAndOld(!sel);
            updatePlots(false,false);
          }
        });
        
        uni_updateCheck.addItemListener(new ItemListener() {
          public void itemStateChanged(ItemEvent e)
          {
            boolean sel; e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED) {sel = false;}       
            else {sel = true;}
             
            mdfa.b_univ = sel;
            updateSignal();
          }
        });        
        
        updatei1Box.addItemListener(new ItemListener() {
          public void itemStateChanged(ItemEvent e)
          {
            boolean sel; e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED) {sel = false;}       
            else {sel = true;}
             
            //--- shades the final n-l+1 points
            if(sel) {mdfa.setUpdateConstraints(1,mdfa.update_i2);}
            else {mdfa.setUpdateConstraints(0,mdfa.update_i2);}
            updateSignal(); 
          }  
        });      

        updatei2Box.addItemListener(new ItemListener() {
          public void itemStateChanged(ItemEvent e)
          {
            boolean sel; e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED) {sel = false;}       
            else {sel = true;}
             
            //--- shades the final n-l+1 points
            if(sel) {mdfa.setUpdateConstraints(mdfa.update_i1,1);}
            else {mdfa.setUpdateConstraints(mdfa.update_i1,0);}
            updateSignal(); 
          }  
        });        
          
        computeUpdate.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent evt) {                
              updateSignal();            
            }
        });

    }
    
    
    
    
  public void computeSpectralDensityEstimate()
  {

   if(activateEstimationCheck.isSelected())
   {
    specDens.computeSpectralDensity(n_obs, n_rep);
    specDensCanvas.setDimensions(K, n_rep);
    specDensCanvas.setPeriodogram(specDens.period_x, n_rep, K);
    specDensCanvas.setMod(specDens.mod_spec, n_rep, K);
    specDensCanvas.setArg(specDens.arg_spec, n_rep, K);
   
    if(autoUpdatesCheck.isSelected())
    {
      mdfa.setSpecDensityEstimate(specDens.mod_spec, specDens.arg_spec, K, n_rep);
      updatePlots(true,true);
    }
   }
  }    
    
  public void computeSpectralDensityEstimateSeries()
  {

   if(activateEstimationCheck.isSelected())
   {
    specDens.computeSpectralDensity(mdfa.tseries, n_obs, n_rep);
    specDensCanvas.setDimensions(K, n_rep);
    specDensCanvas.setPeriodogram(specDens.period_x, n_rep, K);
    specDensCanvas.setMod(specDens.mod_spec, n_rep, K);
    specDensCanvas.setArg(specDens.arg_spec, n_rep, K);
   
    if(autoUpdatesCheck.isSelected())
    {
      mdfa.setSpecDensityEstimate(specDens.mod_spec, specDens.arg_spec, K, n_rep);
      updatePlots(true,true);
    }
   }
  }        
    
  @SuppressWarnings({ "unchecked", "rawtypes" })
private void initSpectralDensityPanel()
  {
        int i;
        specDens.computeSpectralDensity(mdfa.tseries, n_obs, n_rep);
        spectralDensPanel = new JPanel();
        plotPanel = new JPanel();   
        periodogramBox = new JCheckBox();
        modBox = new JCheckBox();
        argBox = new JCheckBox();
        smoothPanel = new JPanel();
        smoothCheck = new JCheckBox();
        smoothFuncBox = new JComboBox();
        smoothFuncLabel = new JLabel();
        smoothScaleLabel = new JLabel();
        smoothScaleBar = new JScrollBar(JScrollBar.HORIZONTAL,1,10,1,200);
        smoothScaleBar.setValue(10);  smoothScaleBar.setUnitIncrement(1);

        smoothScaleText = new JTextField();
        taperPanel = new JPanel(); taperBox = new JCheckBox(); orthBasisBox = new JComboBox();
        taperBasisLabel = new JLabel();
       
        taperDegreeBar = new JScrollBar(JScrollBar.HORIZONTAL,1,1,1,30);
        taperDegreeBar.setValue(10); taperDegreeBar.setUnitIncrement(1);

        taperDegreeLabel = new JLabel();
        taperDegreeText = new JTextField();
        gaussPanel = new JPanel(); gaussCheck = new JCheckBox();
        pComboBox = new JComboBox();
        pLabel = new JLabel(); qLabel = new JLabel();
        qComboBox = new JComboBox();
        estimateARMACheck = new JCheckBox();
        estimateARMAButton = new JButton();
        phaseTypeBox = new JComboBox();
        phaseTypeLabel = new JLabel();
        activateEstimationCheck = new JCheckBox();
        setEstimateButton = new JButton();
        autoUpdatesCheck = new JCheckBox();

        plotPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        sd_checkBox = new JCheckBox[8];

        for(i = 1; i<= 8; i++) 
        {
          sd_checkBox[i-1] = new JCheckBox("Series " + i);
          sd_checkBox[i-1].setHorizontalTextPosition(SwingConstants.LEADING);
        }

        periodogramBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        periodogramBox.setText("Periodogram");
        periodogramBox.setHorizontalTextPosition(SwingConstants.LEADING);

        modBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        modBox.setText("Mod");
        modBox.setHorizontalTextPosition(SwingConstants.LEADING);

        argBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        argBox.setText("Arg");
        argBox.setHorizontalTextPosition(SwingConstants.LEADING);




        ItemListener il = new ItemListener() { 
          public void itemStateChanged(ItemEvent e)
          {
           int i; boolean sel; //computeFilter = true;
           Object source = e.getItemSelectable();
           if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
           else{sel = true;}

           if(source == smoothCheck)
           {specDens.setSmoothing(sel); enableSmooth(sel); if(activateEstimationCheck.isSelected()) {computeSpectralDensityEstimate();}}
           else if(source == taperBox)
           {specDens.setMultitaper(sel); enableTapering(sel); if(activateEstimationCheck.isSelected()) {computeSpectralDensityEstimate();}}            
           else if(source == gaussCheck)
           {specDens.setGaussianize(sel); enableGaussianize(sel); if(sel && activateEstimationCheck.isSelected()) {computeSpectralDensityEstimateSeries();}}
           else if(source == estimateARMACheck)
           {specDens.setARMAModel(sel);}
           else if(source == periodogramBox)
           {specDensCanvas.setPlotPer(sel);}
           else if(source == modBox)
           {specDensCanvas.setPlotMod(sel);}
           else if(source == argBox)
           {specDensCanvas.setPlotArg(sel);}
           else if(source == activateEstimationCheck)
           {
              computeSpectralDensityEstimateSeries();
              if(!sel)
              {mdfa.useSpectralDensityEstimate(false);}
           }       
           for(i=0;i<8;i++) 
           { 
            if(source == sd_checkBox[i]) {specDensCanvas.setPlot(sel, i);}
           }
                       
          }
        };
           
       ItemListener il2 = new ItemListener() { 
          public void itemStateChanged(ItemEvent e)
          {
           int i; boolean sel; //computeFilter = true;
           Object source = e.getItemSelectable();
           if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
           else{sel = true;}

           if(source == periodogramBox)
           {specDensCanvas.setPlotPer(sel);}
           else if(source == modBox)
           {specDensCanvas.setPlotMod(sel);}
           else if(source == argBox)
           {specDensCanvas.setPlotArg(sel);}     
           for(i=0;i<8;i++) 
           {if(source == sd_checkBox[i]) {specDensCanvas.setPlot(sel, i);}}

          }
        };           
           
           
        smoothCheck.addItemListener(il);
        taperBox.addItemListener(il);
        gaussCheck.addItemListener(il);
        estimateARMACheck.addItemListener(il);
        periodogramBox.addItemListener(il2);
        modBox.addItemListener(il2);
        argBox.addItemListener(il2);
 
        for(i = 0; i < 8; i++) 
        {sd_checkBox[i].addItemListener(il2);}

        AdjustmentListener al = new AdjustmentListener()  {
         public void adjustmentValueChanged(AdjustmentEvent e) {        
          
             double smooth; int tap;

             if(e.getAdjustable() == smoothScaleBar)
             {
                smooth = .01*smoothScaleBar.getValue();
                specDens.setScale(smooth);
                smoothScaleText.setText(df.format(smooth));

             }
             else if(e.getAdjustable() == taperDegreeBar)
             {
                tap = taperDegreeBar.getValue();

                specDens.setSpectralConcentration(1.0*tap/(double)n_obs);
                specDens.setMultitaperLevels(tap);
                taperDegreeText.setText(""+tap);
             }
             computeSpectralDensityEstimateSeries();
         }
        };

        smoothScaleBar.addAdjustmentListener(al);
        taperDegreeBar.addAdjustmentListener(al);


      ActionListener mc = new ActionListener() 
      {

       int pc,qc;
       public void actionPerformed(ActionEvent e)
       {

        if(e.getSource() == qComboBox)
        {                         
         pc = ((JComboBox)e.getSource()).getSelectedIndex();
         qc = ((JComboBox)e.getSource()).getSelectedIndex();
         if(estimateARMACheck.isSelected())
         {
           specDens.computeARMASpectralDensity(pc, qc);
           computeSpectralDensityEstimate();
         }
         
        }
        else if(e.getSource() == smoothFuncBox)
        {specDens.setSmoothFunction(smoothFuncBox.getSelectedIndex()); computeSpectralDensityEstimate();}
        else if(e.getSource() == orthBasisBox)
        {
          specDens.setMultitaper(true);
          if(orthBasisBox.getSelectedIndex() == 1)
          {specDens.setQuadratic(true);}
          else{specDens.setQuadratic(false);}
          computeSpectralDensityEstimate();
        }
        else if(e.getSource() == phaseTypeBox)
        {
          specDens.changePhaseInformation(phaseTypeBox.getSelectedIndex());
          computeSpectralDensityEstimate();
        }
       }
      }; 
 
        smoothFuncBox.addActionListener(mc);
        orthBasisBox.addActionListener(mc);
        pComboBox.addActionListener(mc);
        qComboBox.addActionListener(mc);





        GroupLayout plotPanelLayout = new GroupLayout(plotPanel);
        plotPanel.setLayout(plotPanelLayout);
        plotPanelLayout.setHorizontalGroup(
            plotPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotPanelLayout.createSequentialGroup()
                //.addGap(32, 32, 32)
                .addComponent(sd_checkBox[0])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[1])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[2])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[3])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[4])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[5])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[6])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sd_checkBox[7])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                //.addGap(107, 107, 107)
                .addComponent(periodogramBox)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(modBox)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(argBox)
                .addContainerGap(69, Short.MAX_VALUE))
        );
        plotPanelLayout.setVerticalGroup(
            plotPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, plotPanelLayout.createSequentialGroup()
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(plotPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(sd_checkBox[0])
                    .addComponent(sd_checkBox[1])
                    .addComponent(sd_checkBox[2])
                    .addComponent(sd_checkBox[3])
                    .addComponent(sd_checkBox[4])
                    .addComponent(sd_checkBox[5])
                    .addComponent(sd_checkBox[6])
                    .addComponent(sd_checkBox[7])
                    .addComponent(periodogramBox)
                    .addComponent(modBox)
                    .addComponent(argBox)))
                //.addGap(314, 314, 314))
        );

        smoothPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        smoothCheck.setText("Smooth Spectral Density");

        smoothFuncBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        smoothFuncBox.setModel(new DefaultComboBoxModel(new String[] { "Wendland1", "Wendland2", "Phi1", "Phi2" }));

        smoothFuncLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        smoothFuncLabel.setText("Smoothing Function:");

        smoothScaleLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        smoothScaleLabel.setText("Smoothing Scale:");

        smoothScaleBar.setOrientation(JScrollBar.HORIZONTAL);

        smoothScaleText.setText("0.1");

        GroupLayout smoothPanelLayout = new GroupLayout(smoothPanel);
        smoothPanel.setLayout(smoothPanelLayout);
        smoothPanelLayout.setHorizontalGroup(
            smoothPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(smoothPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(smoothPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(smoothPanelLayout.createSequentialGroup()
                        .addComponent(smoothScaleBar, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(smoothScaleText, GroupLayout.PREFERRED_SIZE, 51, GroupLayout.PREFERRED_SIZE))
                    .addComponent(smoothCheck)
                    .addGroup(smoothPanelLayout.createSequentialGroup()
                        .addComponent(smoothFuncLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(smoothFuncBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                    .addComponent(smoothScaleLabel)))
                //.addContainerGap(50, Short.MAX_VALUE))
        );
        smoothPanelLayout.setVerticalGroup(
            smoothPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(smoothPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(smoothCheck)
                .addGap(18, 18, 18)
                .addGroup(smoothPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(smoothFuncLabel)
                    .addComponent(smoothFuncBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(smoothScaleLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(smoothPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(smoothScaleBar, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(smoothScaleText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                //.addContainerGap(41, Short.MAX_VALUE))
        );

        taperPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        taperBox.setText("Quadratic Tapering");

        orthBasisBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        orthBasisBox.setModel(new DefaultComboBoxModel(new String[] { "Sine", "Slepian" }));

        taperBasisLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        taperBasisLabel.setText("Orthonormal Basis:");

        taperDegreeBar.setOrientation(JScrollBar.HORIZONTAL);

        taperDegreeLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        taperDegreeLabel.setText("Degree of Tapering:");

        taperDegreeText.setText("10");

        GroupLayout jPanel1Layout = new GroupLayout(taperPanel);
        taperPanel.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(taperBox)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(taperBasisLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(orthBasisBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(taperDegreeBar, GroupLayout.PREFERRED_SIZE, 163, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(taperDegreeText, GroupLayout.PREFERRED_SIZE, 51, GroupLayout.PREFERRED_SIZE))
                    .addComponent(taperDegreeLabel))
                .addContainerGap(83, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(taperBox)
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(taperBasisLabel)
                    .addComponent(orthBasisBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                
                .addComponent(taperDegreeLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(taperDegreeBar, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(taperDegreeText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        gaussPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        gaussCheck.setText("Gaussianize Data");

        pComboBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        pComboBox.setModel(new DefaultComboBoxModel(new String[] { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" }));

        pLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        pLabel.setText("p:");

        qLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        qLabel.setText("q:");

        qComboBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        qComboBox.setModel(new DefaultComboBoxModel(new String[] { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" }));

        estimateARMACheck.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        estimateARMACheck.setText("Estimate ARMA Spectral Density");

        estimateARMAButton.setText("Estimate");

        phaseTypeBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        phaseTypeBox.setModel(new DefaultComboBoxModel(new String[] { "ARMA Phase", "DFT Phase", "Unit Phase" }));

        phaseTypeLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        phaseTypeLabel.setText("Phase:");

        GroupLayout jPanel3Layout = new GroupLayout(gaussPanel);
        gaussPanel.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel3Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(gaussCheck)
                    .addComponent(estimateARMACheck)
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addGroup(jPanel3Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel3Layout.createSequentialGroup()
                                .addComponent(pLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(pComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(qLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(qComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                            .addGroup(jPanel3Layout.createSequentialGroup()
                                .addComponent(phaseTypeLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(phaseTypeBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                        .addGap(18, 18, 18)
                        .addComponent(estimateARMAButton)))
                .addContainerGap(25, Short.MAX_VALUE))
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(gaussCheck)
                .addGap(22, 22, 22)
                .addComponent(estimateARMACheck)
                .addGap(18, 18, 18)
                .addGroup(jPanel3Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(pComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(pLabel)
                    .addComponent(qLabel)
                    .addComponent(qComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(jPanel3Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(phaseTypeBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(phaseTypeLabel)
                    .addComponent(estimateARMAButton))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        activateEstimationCheck.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        activateEstimationCheck.setText("Activate Estimation");

        setEstimateButton.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        setEstimateButton.setText("Set Estimate");

        autoUpdatesCheck.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        autoUpdatesCheck.setText("Automatic Updating");

        GroupLayout layout = new GroupLayout(spectralDensPanel);
        spectralDensPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(plotPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(smoothPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(taperPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(gaussPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(activateEstimationCheck)
                            .addComponent(autoUpdatesCheck)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(20, 20, 20)
                                .addComponent(setEstimateButton)))))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(plotPanel, GroupLayout.PREFERRED_SIZE, 52, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(smoothPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(activateEstimationCheck)
                        .addGap(18, 18, 18)
                        .addComponent(autoUpdatesCheck)
                        .addGap(18, 18, 18)
                        .addComponent(setEstimateButton)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(gaussPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(taperPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
    

    
      enableSmooth(false); enableTapering(false); enableGaussianize(false);
    
  }    
    public void enableSmooth(boolean t)
    {smoothFuncBox.setEnabled(t); smoothScaleBar.setEnabled(t);}

    public void enableTapering(boolean t)
    {orthBasisBox.setEnabled(t); taperDegreeBar.setEnabled(t);}

    public void enableGaussianize(boolean t)
    {
      estimateARMAButton.setEnabled(t); pComboBox.setEnabled(t); qComboBox.setEnabled(t); 
      //phaseTypeBox.setEnabled(t);    
    }


	public accountCanvas getAccount_canvas() {
		return account_canvas;
	}


	public void setAccount_canvas(accountCanvas account_canvas) {
		this.account_canvas = account_canvas;
	}    
    


}



