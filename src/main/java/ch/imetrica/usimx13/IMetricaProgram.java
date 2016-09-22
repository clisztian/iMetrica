package ch.imetrica.usimx13;

import java.io.*;
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
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import javax.swing.filechooser.FileFilter;
import java.util.ArrayList;
import java.util.Random;

import ch.imetrica.bayesCronos.BayesCronos;
import ch.imetrica.dataControl.JHighFreq;
import ch.imetrica.dataControl.SimPanel;
import ch.imetrica.emd.EMDamCanvas;
import ch.imetrica.emd.EMDamfm;
import ch.imetrica.mdfa.IMDFAPanel;
import ch.imetrica.mdfaTradingStrategies.EvolutionPanel;
import ch.imetrica.mdfaTradingStrategies.MDFAStrategyPanel;
import ch.imetrica.regComponent.REGmodelPanel;




public class IMetricaProgram extends JFrame implements ActionListener,AdjustmentListener,ItemListener
{

        /**
	 * 
	 */
	    private static final long serialVersionUID = 1L;
	    private IMDFAPanel mdfa;
        private SARIMAmodelCanvas smc;
        private SARIMAspectrumCanvas spec;
        private SARIMApseudoCanvas pseudo;
        private REGmodelPanel reg;          
        private BayesCronos BayesCronosPan; 
        private MDFAStrategyPanel mdfaStrat; 
        private EvolutionPanel evolve; 
         
        //---- Canvas for EMD ---------
        EMDamfm emd; 
        EMDamCanvas emdAM;
        EMDfmcanvas emdFM;
        EMDPhasecanvas emdPhase;
        EMDfmcanvas emdIF;
        EMDIMFcanvas emdIMF;

        EMDMap emdm;
        EMDMapKey emdKey;
        EMDMapScale emdScale;
        EMD3dScale hScale;
        EMDFreqScale freqScale;

        JPanel amfmComps, phaseComps;
        JCheckBox[] imfs;
        JCheckBox[] ams;
        JCheckBox[] fms;
        JCheckBox[] ifs;
        JCheckBox[] phase;
        JRadioButton[] mps2d;
        JPanel jplRadio;

        boolean emdPanelx, sarimaPanelx, mdfaPanelx, regPanelx, simPanelx, bayesPanelx, mdfaStratx, evolveX;
        //-----------------------------
        
	    private JScrollBar innvarParameter, nObsBar, seedBar, burnBar;
	    private JTextField text1, text2, text3, text4;
	    private JTextField text5, text6, text7, text8, textLag;
        private JTextField text9, text10, text11, text12, text13, text0;
        private JTabbedPane tabbedPane;
        private JTabbedPane tabbedContPane;
        private JTabbedPane paramTabbedPane;
        private JTabbedPane emdPlotPane;

        public JButton uSimX13compute;

               
        private JComboBox<String> p,q,P,Q;
        public JComboBox<String> pp,pq,pP,pQ;
        public JComboBox<String> tp,tq,tP,tQ;


	    private int nObs, seed, burnin, lag, n_params;
        private double innvar;
        private int m_p, m_q, m_P, m_Q;
        private int p_p, p_d, p_q, p_P, p_D, p_Q;
        private int[] dims;
        private double[] params; 
        private double[] mleparams;
        
        int minute_df = 15;
        private double phi1, phi2, Phi, theta1, theta2, Theta;
        private int iphi1, iphi2, iPhi, itheta1, itheta2, iTheta, i_innvar;
        private JScrollBar phi1bar, phi2bar, Phibar;
        private JScrollBar theta1bar, theta2bar, Thetabar;
        private JScrollBar lagBar;
        private DecimalFormat df;
        //x13 estimates
        private JTextField textphi1, textphi2, textPhi, textInnvar, textSigvar;
        private JTextField textTheta, texttheta1, texttheta2;
        private JTextField textLk1, textLk2, textLk3, textLk4;
        private JTextField textLB12, textLB24, textLB0, textLB1;

        private boolean estimate;
        private JScrollBar slidingSpan;
	private JTextField slidingSpanText;
	//--------- Check boxes for different plots -------
        private JCheckBox seriesBox;
        private JCheckBox foreBox;
        private JCheckBox seasBox;
        private JCheckBox trendBox;
        private JCheckBox saBox;
        private JCheckBox cycleBox;
        private JRadioButton TrendModel; 
        private JRadioButton SeasonalModel; 
        private JRadioButton TrendIrregModel; 
        private JRadioButton IrregModel; 
       
        //----------------- Controls and data for the pseudo panel --------------
        public int model;
        private boolean pseudoCntrl; 
        private boolean specCntrl;
        private double[] real_series;
        private double[] orig_series;
        public JScrollBar series_scroll;
        public JTextField series_name;
        public int nObsFixed;
        
        
        
        //--------------- group series stuff for uploading metafiles------------
        private double[][] group_series;
        private int[] gs_lengths;
        private String[] list_files;
        private int n_series, i_val; 
        
        boolean sim_series;
        boolean iter;
        private JCheckBoxMenuItem sim_check;

        //------ true model stuff
        double[] true_params;
        int[] true_dim;
        int n_true;

        //------ pseudo-true model stuff
        double[] pseudo_params;
        int[] pseudo_dim;
        int n_pseudo;
        
	public JTextField ptext1, ptext2, ptext3, ptext4;
	public JTextField ptext5, ptext6, ptext7, ptext8, ptextLag;
        public JTextField ptext9, ptext10, ptext11, ptext12, ptext13, ptext0;
        public JTextField ptextefft, ptexteffs, ptexteffi, ptexteffti, ptextMinKL;

        public JTextField ptextphi1, ptextphi2, ptextPhi, ptextInnvar;
        public JTextField ptextTheta, ptexttheta1, ptexttheta2;
        public JTextField ptextLk1, ptextLk2, ptextLk3, ptextLk4;
        public JTextField ptextLB12, ptextLB24;

    JTextField endDay,endMonth,endYear;
    //JLabel yearLabel,startLabel,instrumLabel,endLabel;
    JButton enterCompute;
    JTextField instrumText,startDay,periodText;
    @SuppressWarnings("rawtypes")
	JComboBox kernelCombo,lagCombo; 
    JTextField startMonth,startYear;
    @SuppressWarnings("rawtypes")
	JComboBox timeCombo;
    JPanel periodCombo;
    JScrollBar periodBar;
    JCheckBox mdfaBox;
    JCheckBoxMenuItem computeErrorCheck;

        public JLabel estLabel1, estLabel2;
          JLabel modelLabel;
          JLabel pLabel;
          JLabel qLabel;
          JLabel QLabel;
          JLabel PLabel;
          JLabel phi1Label;
          JLabel phi2Label;
          JLabel theta1Label;
          JLabel theta2Label;
          JLabel ThetaLabel;

          JLabel PhiLabel;
          
          JLabel diagLabel;
          JLabel tlabel;
          JLabel slabel;
          JLabel ilabel;
          JLabel tilabel;

        public JMenuBar menuBar;
        public JFileChooser fc;
        public JButton openButton;

        //---- New spectral plotting options added 12-1-10
        public JPanel sdPlots;
        public boolean plotIn, plotFw, plotFu, plotDn, plotGF, plotGI, plotG;        
        public JCheckBox plotInBox, plotFwBox, plotFuBox, 
               plotDnBox, plotIGBox, plotFGBox, plotGBox;

        public Color myBlue;


        //----- EMD Stuff ---------
        public boolean[] plots_agg;
        public boolean[][] plots_emd;
        public int n_imfs;
        public boolean c_emd; 
        public int selMap;
        public JButton reCompEMD;
        

        //-------- Regression Stuff -------------------
        public boolean trans, td, easter, outlier;
        public JCheckBox tdBox, easterBox, outlierBox, transBox;
        public JComboBox<String> easterDay; 
        public int easterDays;

        //-------- Menu stuff ------------------------
        public JMenuItem seriesItem;
        public JMenuItem seriesFItem;
        public JMenuItem seasItem;
        public JMenuItem trendItem;
        public JMenuItem saItem;
        public JMenuItem cycleItem;  

        public JMenuItem trendPoly;
        public JMenuItem seasPoly;
        public JMenuItem tiPoly;
        public JMenuItem trnsPoly;    
        public JMenuItem simulMenu;
        public JMenuItem diagnosticMenu;
        public JMenuItem marketMenu;
        public JMenuItem uploadFilter;
        public JMenu panelMenu;
        public JCheckBoxMenuItem autoComp;
       
        public JCheckBoxMenuItem slideTrade; 
        public JMenuItem uploadPrice;
        public JMenuItem tradeparamDialog; 
        public JMenuItem tradeoptimDialog; 
        public JMenuItem tradestatDialog;
        public JCheckBoxMenuItem tradingMode;
        
        public JMenuItem track_account,track_price,track_logreturns,track_signal,track_none;
        public JMenuItem track_series,track_cycle,track_seas,track_trend,track_SA,track_noneUsim;

        public String curDir; 
        public FileNameExtensionFilter filter;

        //------- Input Data --------------------------
        public double[][] userData;
         
        //------- Spectral density polynomials --------
        public double[] signalPoly; 
        public double[] f_w_poly; 

        //----- uSIMX13-S----------------------
        public JMenu fileMenu;
        public JMenuItem slideSpanMenu;
        //------- REG menu stuff ---------------
        public JMenu regMenu;
        
        public JMenuItem getComps;
        public JMenuItem getData;           // --- upload data
        public JMenuItem[] compItems;       // --- number of compItems
        public JCheckBoxMenuItem toFile;    // --- Check box for printing to File or uSim
        public int n_regCmpnts;             // --- keeps track of number of regComps
        public boolean regFore;             // --- print with forecasts
        public boolean toFilex;             
       

        //----- I-MDFA menu stuff ------------
        public JMenu mdfaMenu;
        public JCheckBoxMenuItem useSarima;
        public JMenuItem saveFilterDes,saveFilterEncryp;
        public JMenuItem savePlotFile; 
        public JMenuItem savePlot;
        public JMenuItem clearPlots;
        
        public JMenuItem[] dfaItems;
        public int n_hist;
        ExtensionFilter dfaFiles;
  
        JHighFreq hfreq;
        JDialog hfreqDialogPanel; 
        public JMenuItem hfreqDialogPanelMenu;
        public JPanel realizedVolPanel;

        //------- Sim Panel Stuff --------------------------------------
        public SimPanel simulate;  
        public JMenu simMenu;
        public JMenuItem saveData;
        public JMenuItem saveTarget;
        public JMenuItem exportMDFA;
        public JMenuItem exportREG;
        public JMenuItem exportSARIMA;
        public JMenuItem exportBayes;
        public JMenuItem sim_mdfa_check;  
        public JMenuItem addFilteredData;
        public JMenuItem dataMixDialog;
        public JMenuItem loadH0Filter;
        public JCheckBoxMenuItem turnOnH0;
        
        
        public JMenuItem csvData;
        public JMenuItem uploadData,uploadReturnData,uploadPortfolioData,uploadTBillData;
        public JMenuItem marketData;

        public JMenu bayesMenu;
         
        public JMenuItem bayes_export_data, bayes_save_single, bayes_update;
        public JCheckBoxMenuItem plot_vol;
        public JRadioButtonMenuItem[] plot_bseries;
        
        
        public JMenuItem sweepSpanMenu;
        public JCheckBoxMenuItem[] seriesItems;

        double sim_phi1, sim_phi2, sim_theta1;
        double sim_theta2, sim_Phi, sim_Theta;
        
        boolean tsim_data;
        double[] sim_data;
        int n_sym;
        Box hBox22;  
        JDialog simulatePanel;
        JDialog dfaDiagnosticPanel;
        JDialog marketDataPanel;

        JDialog tradeoptimDialogPanel;
        JDialog tradeparamDialogPanel;
        JDialog tradestatDialogPanel;
        JDialog tradesoutstatDialogPanel;
        JDialog dataMixDialogPanel;
        JPanel marketChoosePanel;
        JMenu tradingMenu;
        JDialog outSampstatDialogPanel;
        JMenuItem outSampstatDialog;
        JDialog adaptiveUpdateDialogPanel;
        JMenuItem adaptiveUpdateMenu;
        JMenuItem tradeoutstatDialog; 
        JDialog envisionDialogPanel;
        JMenuItem envisionDialog;
        
    private JLabel endLabel;
    private JLabel yearEndLabel,monthEndLabel,dayEndLabel;
    private JLabel yearLabel,startLabel,monthLabel,mainLabel,instrumLabel,hoursLabel,dayLabel,frequencyLabel;
    private JTextField dayText;
    private JComboBox<String> frequencyCombo,hoursCombo;
    private JTextField yearText,monthText,instrumText2,yearEndText,monthEndText,dayEndText;
    private JButton enterMarket; 
    private JCheckBox logReturns, volumeData, fromYahoo, clearData, high_lowData, high_lowDataReturn;
    private JCheckBoxMenuItem tradingInterface;
    private JCheckBoxMenuItem slideSpanCheckMenu;
    // Variables declaration - do not modify
    @SuppressWarnings("rawtypes")
	private JComboBox googDaysBox;
    private JLabel googDaysLabel,googExchLabel;
    private JTextField googExchText;
    @SuppressWarnings("rawtypes")
	private JComboBox googFreqBox;
    private JLabel googFreqLabel;
    private JCheckBox googHiLoDiffBox,googHiLowBox,googLogReturnsBox;
    private JLabel googSymbolsLabel;
    private JTextField googSymbolsText;
    private JCheckBox googVolBox;
    private JPanel googleIntradayPanel;
    private JButton googIntraday;

    // Variables declaration - do not modify                     
    private JPanel quandlPanel;
    private JLabel fromQLable;
    private JTextField fromQText;
    private JButton goQuandlButton;
    private JLabel quandlLable;
    private JLabel searchQLabel;
    private JTextField searchQText;
    private JLabel toQLabel;
    private JTextField toQText;        
    @SuppressWarnings("rawtypes")
	private JComboBox quandlComboBox;
    private JCheckBox quandlFinance;
    private JLabel quandlFreqLabel; 
    
    public JMenuItem quandlPanelMenu;
    public JMenuItem googIntradayMenuItem;
    public JDialog googIntradayDialog,quandlDialog;

    public JScrollBar nbackObsBar,nforeObsBar;
    private JLabel nbackObsLabel,nforeObsLabel;
    private JTextField nbackObsText,nforeObsText;
    private JCheckBox slideSpanCheck;
    private JPanel slideSpanPanel;
    public JDialog slideSpanDialog;
    public JDialog sweepSpanDialog;
    public JDialog trueOutSampleDialogPanel;
    public JMenuItem trueOutSamplePanelMenu;
    
    public JMenuItem loadFilterCoeffsMenu;
    public JMenuItem saveFilterCoeffsMenu;
    public JMenuItem loadPriceFilterCoeffs;
 
    public JCheckBoxMenuItem periodoWeightBox;
    public JMenuItem tradeLogDiff,tradeLogPrice, tradeDuplex, tradeDuplex2, tradeMultiFreq;
    public JMenuItem IF5mItem;
    public JLabel googHigherFreqPeriodLabel;
    public JCheckBox googHigherFreq;
    @SuppressWarnings("rawtypes")
	public JComboBox googHigherFreqPeriod;
    

        //---stuff to fill market dialog box 
        


	public IMetricaProgram(String title)
	{

           /*-----------------------------------------------------------------
             Initialize Global Parameters 
           -------------------------------------------------------------------*/

	       super(title); 
           df = new DecimalFormat("##.##"); myBlue = new Color(146,196,210);
           new DecimalFormat("###");
           //-----initial values-----
	       nObs = 144; seed = 100;
           burnin = 100; lag = 0; 
           estimate = true;
           td=false; easter=false; outlier=false; easterDays = 1;
           initializeUseless(nObs);  regFore = false; toFilex = true;
           n_hist = 0;
           //Initialize parameters
           phi1 = 0.0; phi2 = 0.0; Phi = 0.0;
           theta1 = -0.60; theta2 = 0.0; Theta = -0.60; innvar = 1.00;
           smc = new SARIMAmodelCanvas(700, 300, nObs, burnin, seed);
           spec = new SARIMAspectrumCanvas(700, 300, 900);
           sim_phi1=0.0; sim_phi2=0.0; sim_theta1=-.6; sim_theta2=0.0; sim_Phi=0.0; sim_Theta=-.6;
         
           params = new double[7];
           mleparams = new double[7];

           params[0] = phi1; params[1] = phi2;
           params[2] = theta1; params[3] = theta2;
           params[4] = Phi; params[5] = Theta; params[6] = innvar;

           dims = new int[6];
           dims[0] = 0; dims[1] = 1; dims[2] = 1;
           dims[3] = 0; dims[4] = 1; dims[5] = 1; 
           m_p = 0; m_q = 1; m_P = 0; m_Q = 1;  
           n_params = m_p + m_q + m_P + m_Q + 1; 
           smc.setParams(dims, params, n_params, innvar);
           smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);
           smc.computeSARIMAmodel(estimate); 
           //smc.setMLEparams();
           smc.computeDataMax();
           Font f = new Font("Dialog", Font.PLAIN, 20);


           tsim_data = false;
           sim_data = new double[nObs];
           n_sym = 100;
          

          //  ---- Initialize the dimensions -----------
          pseudo = new SARIMApseudoCanvas(700, 200, nObs, burnin);         
          true_dim = new int[6]; pseudo_dim = new int[6];
          true_dim = dims; 

          //------- Initialize the emd plots----------------
          plots_agg = new boolean[5];
          plots_emd = new boolean[5][4];
          for(int i=0;i<5;i++)
          {
           plots_agg[i] = false; plots_emd[i][0] = false;   
           plots_emd[i][1] = false; plots_emd[i][2] = false; plots_emd[i][3] = false;
          }
          emdPanelx = false; sarimaPanelx = false; c_emd = false; selMap = 0;
          mdfaPanelx = false; regPanelx = false; simPanelx =true;
          //--------------------------------------------------
                   
          p_p = 0; p_q = 1; p_P = 0; p_Q = 1; p_d = 1; p_D = 1;
          pseudo_dim[0] = p_p; pseudo_dim[1] = p_d; pseudo_dim[2] = p_q; 
          pseudo_dim[3] = p_P; pseudo_dim[4] = p_D; pseudo_dim[5] = p_Q;  
        
          n_true = true_dim[0] + true_dim[2] + true_dim[3] + true_dim[5] + 1;
          n_pseudo = pseudo_dim[0] + pseudo_dim[2] + pseudo_dim[3] + pseudo_dim[5] + 1;

          true_params = params;
          pseudo_params = new double[7];

          pseudo.setPseudoParameters(pseudo_dim);
          pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
          pseudo.setModel(6);
          pseudo.setSeed(seed); 
          
          pseudoCntrl = false; specCntrl = false;
          smc.setSpecPlots(specCntrl);

          dfaFiles = new ExtensionFilter("Data files", new String[] {".sig", ".dat", ".coeff"});
          curDir=System.getProperty("user.dir") + File.separator;
          fc = new JFileChooser(curDir);
          fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
          fc.addChoosableFileFilter(dfaFiles);  
          //filter=new FileNameExtensionFilter("uSim files","dta","dat","usim","cfs");
          //fc.addChoosableFileFilter(filter);

          //String[] ex_files = new String[] {"dat","dta","usim","cfs"};
          //fc.addChoosableFileFilter(new ExtFilter(ex_files, "uSim files (*.dta, *.dat, *.usim, *.cfs)"));
          
          //------- Added May 25 2011 -------------------------------------
          sim_series = true; 
          n_series = 0;
          
          //-------Added May 31st 2011--------------------------------------
          //series_label = new JLabel("Series:");

          series_name = new JTextField(10);
          //----------------------------------------------------------------           



          modelLabel = new JLabel("Model  ");
          modelLabel.setToolTipText("Choose SARIMA model dimensions.");
          pLabel = new JLabel("  p:");
          pLabel.setToolTipText("Nonseasonal autoregressive order");
          qLabel = new JLabel("  q:");
          qLabel.setToolTipText("Nonseasonal moving-average order");
          QLabel = new JLabel("  Q:");
          QLabel.setToolTipText("Seasonal moving-average order");
          PLabel = new JLabel("  P:");
          PLabel.setToolTipText("Seasonal autoregressive order");
          phi1Label = new JLabel("\u03D5_1"); phi1Label.setToolTipText("Set value in simulation for \u03D5_1.");
          phi2Label = new JLabel("\u03D5_2"); phi2Label.setToolTipText("Set value in simulation for \u03D5_2."); 
          theta1Label = new JLabel("\u03D1_1"); theta1Label.setToolTipText("Set value in simulation for \u03D1_1."); 
          theta2Label = new JLabel("\u03D1_2"); theta2Label.setToolTipText("Set value in simulation for \u03D1_2.");  
          ThetaLabel = new JLabel("\u0398");    
          ThetaLabel.setFont(f);                ThetaLabel.setToolTipText("Set value in simulation for \u0398."); 
          PhiLabel = new JLabel("\u03A6");      PhiLabel.setToolTipText("Set value in simulation for \u03A6.");  
          PhiLabel.setFont(f);
          diagLabel = new JLabel("GOF Signal Extraction Diagnostics");
          diagLabel.setToolTipText("Diagnostics for assessing model goodness-of-fit based on signal-noise decomposition.");
          tlabel = new JLabel("Trend:");
          tlabel.setToolTipText("GOF diagnostic where signal-noise given by trend and seasonal+irregular.");
          slabel = new JLabel("Seasonal:");
          slabel.setToolTipText("GOF diagnostic where signal to noise given by seasonal to trend+irregular.");
          ilabel = new JLabel("Irregular:");
          ilabel.setToolTipText("GOF diagnostic where signal to noise given by irregular/cycle to seasonal+trend.");
          tilabel = new JLabel("Trend-Irregular:");
          tilabel.setToolTipText("GOF diagnostic where signal to noise given by irregular+trend to seasonal.");

          textphi1 = new JTextField(3); textphi2 = new JTextField(3); 
          textPhi = new JTextField(3); textInnvar = new JTextField(3);
          textTheta = new JTextField(3); texttheta1 = new JTextField(3); 
          texttheta2 = new JTextField(3);textLk1 = new JTextField(3); 
          textLk2 = new JTextField(3); textLk3 = new JTextField(3); 
          textLk4 = new JTextField(3); textLB12 = new JTextField(3); 
          textLB24 = new JTextField(3); textSigvar = new JTextField(3);
          textLB0 = new JTextField(3);  textLB1 = new JTextField(3); 
          //===========================================================
          
          // 3.5 ============= Configure MLE panel ===============
          ptextphi1 = new JTextField(3); ptextphi2 = new JTextField(3); ptextPhi = new JTextField(3); 
          ptextInnvar = new JTextField(3);
          ptextTheta = new JTextField(3); ptexttheta1 = new JTextField(3); ptexttheta2 = new JTextField(3);
          ptextLk1 = new JTextField(3); ptextLk2 = new JTextField(3); ptextLk3 = new JTextField(3); 
          ptextLk4 = new JTextField(3);
          ptextLB12 = new JTextField(3); ptextLB24 = new JTextField(3); ptextLk1 = new JTextField(4); 
          ptextLk2 = new JTextField(4);
          ptextLk3 = new JTextField(4); ptextLk4 = new JTextField(4); ptextMinKL = new JTextField(3);  
          reCompEMD = new JButton("Recompute");
           
	   Container cp = getContentPane();
         
           //Box sarimaPane = Box.createVerticalBox();

           JPanel sarimaPane = new JPanel();
           JPanel emdPane = new JPanel();
           
          /*--------------------------------------------------------------------
            Top Panel - Top panel for Menu, Observations, SARIMA Model Param
          ----------------------------------------------------------------------*/
        
           Box.createVerticalBox();  
                
           //-------Setup Controls and add to topBox-----
           //setUpMenu(topBox);
           setUpMenu(this);
           JPanel controls = new JPanel();  
           TitledBorder contBorder = new TitledBorder(new LineBorder(myBlue), "Observations");
           //contBorder.setTitleColor(Color.LIGHT_GRAY);
           contBorder.setTitleColor(Color.BLACK);
           controls.setBorder(contBorder);
           controls.setLayout(setUpControlLayout(controls));
           initSlidingSpanComponents();
           smc.activateSweepInterface();

           // ------------------ MDFA ---------------------------------
           mdfa = new IMDFAPanel();
           reg = new REGmodelPanel(smc.t_series, this, 700, 300);
           simulate = new SimPanel(300,10,100,this);
 
 
           //------------------ Bayesian ------------------------------
           BayesCronosPan = new BayesCronos(this);
           mdfaStrat = new MDFAStrategyPanel(this,this);
           mdfaStrat.setPreferredSize(new Dimension(700,300));
           mdfaStrat.initStrategyPanel();
           
           evolve = new EvolutionPanel(this,this);
           evolve.setPreferredSize(new Dimension(700,300));
           
           //---- setup sarima model params panel --------
           JPanel paramPanel = new JPanel(); 
           TitledBorder paramBorder = new TitledBorder(new LineBorder(myBlue), "Model Simulation Parameters");
           //paramBorder.setTitleColor(Color.LIGHT_GRAY);
           paramBorder.setTitleColor(Color.BLACK);
           paramPanel.setBorder(paramBorder);
           paramPanel.setLayout(setUpParameterLayout(paramPanel));
           //---- Put in horizontal box
           hBox22 = Box.createHorizontalBox();       
           hBox22.add(controls);
           hBox22.add(paramPanel);           
           //---- add to top box 
           //topBox.add(hBox22);
           //---- add top box to main content pane in the North
           //cp.add(topBox, BorderLayout.NORTH);     
           //---------------------------------------------------------------
            
           simulatePanel = new JDialog(this,true);
           simulatePanel.getContentPane().add(hBox22);
           simulatePanel.pack();
           simulatePanel.setLocationRelativeTo(this);
           simulatePanel.setModal(false);
           simulatePanel.setVisible(false);




           /*----------------------------------------------------------------
             Plotting component for the modeling tab
           -----------------------------------------------------------------*/
           
           //------------- create plot box vertical ----------------
           Box plotsBox = Box.createVerticalBox();  
         
           //-------------Add plotting tabbed pane as center------------       	   
           tabbedPane = new JTabbedPane(JTabbedPane.TOP);
           tabbedPane.addTab("Time Domain", smc);
           tabbedPane.addTab("Spectral Domain", spec);
           tabbedPane.addChangeListener(new ChangeListener() {
              public void stateChanged(ChangeEvent evt) 
              { 
               
               JTabbedPane pane = (JTabbedPane)evt.getSource();
               int sel = pane.getSelectedIndex(); 
               if(sel == 1){specCntrl = true; smc.setSpecPlots(specCntrl); smc.computeSampleSpec();} 
               else{specCntrl = false;} 
               setSpecControls();
               smc.setSpecPlots(specCntrl); 
              }});
           
           //----- Add Tabbed Pane for Plots and Pane for plotting options
           plotsBox.add(tabbedPane);
           plotsBox.add(createPlottingBox()); //plotsBox.add(hBoxTop);
           //sarimaPane.add(plotsBox,BorderLayout.NORTH);
           //-----------------------------------------------------------------
         
           //----- Create Boxes for bottom MLE Panel------
           Box baseBox = Box.createVerticalBox();
           Box hBox1 = Box.createHorizontalBox();
           Box hBox2 = Box.createHorizontalBox();

           //----- create top two horiz. boxes-----
           hBox1.add(createModelPanel()); //modelPanel);
           hBox1.add(createMLEPanel()); //mlePanel);

           //----- add to base box-----
           baseBox.add(hBox1);
             
           //----- create botton two horiz. boxes-----
           hBox2.add(createDiagPanel()); //diagPanel);
           hBox2.add(createX13Panel()); //x13estPanel);
           baseBox.add(hBox2);
 

       //--------------- Pseudo Panel Stuff -------------
        
          //---- Horizontal and vertical boxes -----
          Box phBox1 = Box.createHorizontalBox();
          Box phBox2 = Box.createHorizontalBox();
          Box phBox3 = Box.createHorizontalBox();  
          Box pseudoBox = Box.createVerticalBox();

          //----- Add the two panels to make top half ----- 
          phBox1.add(createPModelPanel()); //pmodelPanel);
          phBox1.add(createPValPanel()); //pvalPanel); 
          //----- Add the two panels to make top half -----
          phBox2.add(createPDiagPanel()); //pdiagPanel);
          phBox2.add(createPeffPanel()); //effPanel);
          //----- Add the bottom stat panels -----
          setPseudoStatistics();
          phBox3.add(createPLBPanel()); //lbPanel);
          phBox3.add(createCritPanel()); //criteria); 
          //---- Add all three horizboxes -------
          pseudoBox.add(phBox1);
          pseudoBox.add(phBox2);
          pseudoBox.add(phBox3); 

          //-----create Tabbed pane----
          tabbedContPane = new JTabbedPane(JTabbedPane.TOP);
          tabbedContPane.addTab("ML Estimation", baseBox);
          tabbedContPane.addTab("Pseudo Estimation", pseudoBox);
          tabbedContPane.setSelectedIndex(0);
          //----- Add listener --------------------
          tabbedContPane.addChangeListener(new ChangeListener() {
              public void stateChanged(ChangeEvent evt) 
              {                
               JTabbedPane pane = (JTabbedPane)evt.getSource();
               int sel = pane.getSelectedIndex(); 
               //if(sel == 1){pseudoCntrl = false;}  ??
               //else{pseudoCntrl = true;}           ??     
               matchSelectedIndices(sel);               
               if(sel == 0)
               {
                 pseudoCntrl = false;
                 smc.computeSARIMAmodel(estimate);
                 smc.computeDataMax();   
                 setDiagnostics();
                 smc.go(); spec.go(); 
               }
               else if(sel == 1)
               {
                  pseudoCntrl = true; //System.out.println("Pseudo comp"); 
                  pseudo.changeModelDimension();
                  //for(int i = 0; i < 7; i++) {System.out.print(true_params[i] + " ");}
                  //System.out.println("");
                  smc.computeSARIMAmodel(false);
                  smc.computeDataMax();
                  pseudo.setData(nObs, smc.t_series);                 
                  pseudo.computeEfficacies();                 
                  setPseudoStatistics();                   
                  smc.go();                 
                  spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); 
                  spec.go();
               }
              }
          });
          //====================================================================================
 
          sarimaPane.add(tabbedContPane, BorderLayout.SOUTH);

        
          GroupLayout sarimaLayout = new GroupLayout(sarimaPane);
          sarimaLayout.setAutoCreateGaps(true); sarimaLayout.setAutoCreateContainerGaps(true);
          sarimaLayout.setHorizontalGroup(sarimaLayout.createSequentialGroup()
          .addGroup(sarimaLayout.createParallelGroup() 
            .addComponent(plotsBox).addComponent(tabbedContPane)));
             
          sarimaLayout.setVerticalGroup(sarimaLayout.createSequentialGroup()
           .addComponent(plotsBox)
           .addComponent(tabbedContPane));  
          sarimaPane.setLayout(sarimaLayout);


          setupEMDPanels(emdPane);
          fileMenu.setEnabled(false); regMenu.setEnabled(false); mdfaMenu.setEnabled(false); simMenu.setEnabled(true);
          //-------------Tabbed Pane for Parametric/NonParametric -----
          paramTabbedPane = new JTabbedPane(JTabbedPane.TOP);
          paramTabbedPane.addTab("Data Control", simulate);
          paramTabbedPane.addTab("Multivariate DFA", mdfa);
          paramTabbedPane.addTab("uSimX13-S Modeling", sarimaPane);
          paramTabbedPane.addTab("EMD", emdPane);         
          paramTabbedPane.addTab("State Space Modeling", reg);
          paramTabbedPane.addTab("BayesCronos", BayesCronosPan);
          paramTabbedPane.addTab("MDFA Strategy", mdfaStrat);
          paramTabbedPane.addTab("Evolution", evolve);
          
          paramTabbedPane.addChangeListener(new ChangeListener() {
              public void stateChanged(ChangeEvent evt) 
              { 
               
               JTabbedPane pane = (JTabbedPane)evt.getSource();
               int sel = pane.getSelectedIndex(); 
               if(sel == 3) 
               {
                 mdfaPanelx = false; bayesPanelx  = false;
                 emdPanelx = true; sarimaPanelx = false; pseudoCntrl = false; regPanelx = false; simPanelx =false;
                 if(!c_emd){computeEMD(smc.t_series,nObs); c_emd = true;}
               }
               else if(sel == 1) 
               {
                 bayesPanelx  = false;
                 mdfaPanelx = true;  emdPanelx = false; sarimaPanelx = false; setX13SeatsFilters(); regPanelx = false; simPanelx =false;
                 nObs = mdfa.n_obs; nObsBar.setValue(nObs); text2.setText(""+ nObs); initializeUseless(nObs);
                 
               }
               else if(sel == 4)
               {
                  bayesPanelx  = false;
                  regPanelx = true; mdfaPanelx = false; emdPanelx = false; sarimaPanelx = false; simPanelx =false;
                  setRegCmpntStuff(); 
               }
               else if(sel == 0)
               {
                 bayesPanelx  = false;
                 mdfaPanelx = false; emdPanelx = false; sarimaPanelx = true; regPanelx = false; simPanelx = true;                  
               }
               else if(sel == 2)
               {
                 bayesPanelx  = false;
                 mdfaPanelx = false; emdPanelx = false; sarimaPanelx = true; regPanelx = false; simPanelx =false;
                 smc.computeSARIMAmodel(estimate); smc.computeDataMax();   setDiagnostics();
                 smc.go(); spec.go();               
               }
               else if(sel == 5)
               {
                 bayesPanelx  = true;
                 mdfaPanelx = false; emdPanelx = false; sarimaPanelx = false; regPanelx = false; simPanelx = false;                  
               
               }  
               else if(sel == 6)
               {
                 bayesPanelx  = false; mdfaStratx = true; evolveX = false;
                 mdfaPanelx = false; emdPanelx = false; sarimaPanelx = false; regPanelx = false; simPanelx = false;   
               
               }
               else if(sel == 7)
               {
                 bayesPanelx  = false; mdfaStratx = false; evolveX = true;
                 mdfaPanelx = false; emdPanelx = false; sarimaPanelx = false; regPanelx = false; simPanelx = false; 
                 
                 //---- Pass the current MDFA settings to strategy tester -----
                 evolve.setMDFAParameters(mdfa, turnOnH0.isSelected());
                 
               }
               regMenu.setEnabled(regPanelx); mdfaMenu.setEnabled(mdfaPanelx); //simMenu.setEnabled(simPanelx);
               fileMenu.setEnabled(sarimaPanelx); diagnosticMenu.setEnabled(mdfaPanelx);
               bayesMenu.setEnabled(bayesPanelx);
              }
                                  
           });
         paramTabbedPane.setSelectedIndex(0);
         cp.add(paramTabbedPane,BorderLayout.SOUTH);
           
         mdfa.inputX13data(smc.t_series);


        dfaDiagnosticPanel = new JDialog(this,true);
        dfaDiagnosticPanel.getContentPane().add(mdfa.diagnosticPanel);
        dfaDiagnosticPanel.pack();
        dfaDiagnosticPanel.setLocationRelativeTo(this);
        dfaDiagnosticPanel.setModal(false);
        dfaDiagnosticPanel.setVisible(false);

        slideSpanDialog = new JDialog(this,"Sliding Window Controls",true);
        slideSpanDialog.getContentPane().add(slideSpanPanel);
        slideSpanDialog.pack();
        slideSpanDialog.setLocationRelativeTo(this);
        slideSpanDialog.setModal(false);
        slideSpanDialog.setVisible(false);

        sweepSpanDialog = new JDialog(this,"Sweep Time Series Controls",true);
        sweepSpanDialog.getContentPane().add(smc.sweepPanel);
        sweepSpanDialog.pack();
        sweepSpanDialog.setLocationRelativeTo(this);
        sweepSpanDialog.setModal(false);
        sweepSpanDialog.setVisible(false);


        createMarketDialog();
        marketDataPanel = new JDialog(this,true);        
        marketDataPanel.getContentPane().add(marketChoosePanel);
        marketDataPanel.pack();
        marketDataPanel.setLocationRelativeTo(this);
        marketDataPanel.setModal(false);
        marketDataPanel.setVisible(false);

        createQuandlDialog();
        quandlDialog = new JDialog(this,true);
        quandlDialog.getContentPane().add(quandlPanel);
        quandlDialog.pack();
        quandlDialog.setLocationRelativeTo(this);
        quandlDialog.setModal(false);
        quandlDialog.setVisible(false);        
        
        
        createGoogleIntradayDialog();
        googIntradayDialog = new JDialog(this,true);
        googIntradayDialog.getContentPane().add(googleIntradayPanel);
        googIntradayDialog.pack();
        googIntradayDialog.setLocationRelativeTo(this);
        googIntradayDialog.setModal(false);
        googIntradayDialog.setVisible(false);
        
        hfreq = new JHighFreq();
        initRealizedVolControl();
        hfreqDialogPanel = new JDialog(this,true);
        hfreqDialogPanel.getContentPane().add(realizedVolPanel);
        hfreqDialogPanel.pack();
        hfreqDialogPanel.setLocationRelativeTo(this);
        hfreqDialogPanel.setModal(false);
        hfreqDialogPanel.setVisible(false);

        tradeoptimDialogPanel = new JDialog(this,true);
        tradeoptimDialogPanel.getContentPane().add(mdfa.tradingoptimPanel);
        tradeoptimDialogPanel.pack();
        tradeoptimDialogPanel.setLocationRelativeTo(this);
        tradeoptimDialogPanel.setModal(false);
        tradeoptimDialogPanel.setVisible(false);

        tradestatDialogPanel = new JDialog(this,true);
        tradestatDialogPanel.getContentPane().add(mdfa.tradingstatPanel);
        tradestatDialogPanel.pack();
        tradestatDialogPanel.setLocationRelativeTo(this);
        tradestatDialogPanel.setModal(false);
        tradestatDialogPanel.setVisible(false);

        tradesoutstatDialogPanel = new JDialog(this,true);
        tradesoutstatDialogPanel.getContentPane().add(mdfa.outsampStatPanel);
        tradesoutstatDialogPanel.pack();
        tradesoutstatDialogPanel.setLocationRelativeTo(this);
        tradesoutstatDialogPanel.setModal(false);
        tradesoutstatDialogPanel.setVisible(false);

        envisionDialogPanel = new JDialog(this,true);
        envisionDialogPanel.getContentPane().add(mdfa.crystal_ball);
        envisionDialogPanel.pack();
        envisionDialogPanel.setLocationRelativeTo(this);
        envisionDialogPanel.setModal(false);
        envisionDialogPanel.setVisible(false);        
        
        trueOutSampleDialogPanel = new JDialog(this,true);
        trueOutSampleDialogPanel.getContentPane().add(mdfa.trueOutSamplePanel);
        trueOutSampleDialogPanel.pack();
        trueOutSampleDialogPanel.setLocationRelativeTo(this);
        trueOutSampleDialogPanel.setModal(false);
        trueOutSampleDialogPanel.setVisible(false);
        
        tradeparamDialogPanel = new JDialog(this,true);
        tradeparamDialogPanel.getContentPane().add(mdfa.tradingparameterPanel);
        tradeparamDialogPanel.pack();
        tradeparamDialogPanel.setLocationRelativeTo(this);
        tradeparamDialogPanel.setModal(false);
        tradeparamDialogPanel.setVisible(false);

        outSampstatDialogPanel = new JDialog(this,true);
        outSampstatDialogPanel.getContentPane().add(mdfa.outsampPanel);
        outSampstatDialogPanel.pack();
        outSampstatDialogPanel.setLocationRelativeTo(this);
        outSampstatDialogPanel.setModal(false);
        outSampstatDialogPanel.setVisible(false);
        
        dataMixDialogPanel = new JDialog(this,true);
        dataMixDialogPanel.getContentPane().add(mdfa.dataMixPanel);
        dataMixDialogPanel.pack();
        dataMixDialogPanel.setLocationRelativeTo(this);
        dataMixDialogPanel.setModal(false);
        dataMixDialogPanel.setVisible(false);
        
        adaptiveUpdateDialogPanel = new JDialog(this,"Adaptive Filtering", true);
        adaptiveUpdateDialogPanel.getContentPane().add(mdfa.adaptiveUpdatePanel);
        adaptiveUpdateDialogPanel.pack();
        adaptiveUpdateDialogPanel.setLocationRelativeTo(this);
        adaptiveUpdateDialogPanel.setModal(false);
        adaptiveUpdateDialogPanel.setVisible(false);
        
        
        estimate = false; specCntrl = false; setSpecControls();
        
      
      //set texts        
      smc.setDiagnosticsFields(text10,text11,text12,text13);
      smc.setLBFields(textLB0, textLB1, textLB12, textLB24);
      smc.setAICFields(textLk1, textLk2, textLk3, textLk4);
      smc.setMLEFields(textphi1, textphi2, texttheta1, texttheta2, textTheta, textPhi, textInnvar);
    }





    public void setSpecControls()
    {
      plotIGBox.setEnabled(specCntrl); plotFwBox.setEnabled(specCntrl); plotGBox.setEnabled(specCntrl);
      plotInBox.setEnabled(specCntrl); plotDnBox.setEnabled(specCntrl); plotFGBox.setEnabled(specCntrl);
    }

    
    public void setRegCmpntStuff()
    {
       int i; n_regCmpnts = reg.getNcmpnts();
   
       for(i=0; i < n_regCmpnts; i++) 
       {
          compItems[i].setEnabled(true);
       }
       for(i = n_regCmpnts; i < 10; i++)
       {
          compItems[i].setEnabled(false); 
       }
 
       regFore = reg.isForecasting(); 
            
    }


    //--------- Use X13seasts filters in IMDFA only if useSARIMA is on -------    
    //------------------------------------------------------------------------
    public void setX13SeatsFilters()
    {
      if(mdfa.useSARIMA)
      {  
       mdfa.setTrendPoly(smc.trendpoly, smc.trendpoly[0],smc.trendpoly[7]);
       mdfa.setSeasPoly(smc.seaspoly, smc.seaspoly[0],smc.seaspoly[26]);
       mdfa.setTIPoly(smc.tipoly, smc.tipoly[0],smc.tipoly[31]);
       mdfa.setModelPoly(smc.mphi, smc.mPhi, smc.mapoly, smc.m_innvar);
      }
    }

   public void enableX13Simulator(boolean sim)
    {
       
       nObsBar.setEnabled(sim); burnBar.setEnabled(sim);
       seedBar.setEnabled(sim); theta1bar.setEnabled(sim);
       theta2bar.setEnabled(sim); Thetabar.setEnabled(sim);
       phi1bar.setEnabled(sim); phi2bar.setEnabled(sim); 
       Phibar.setEnabled(sim);  innvarParameter.setEnabled(sim);
       sim_check.setSelected(sim); 
    }


    public void setMLEs()
    {
      int j = 0;
    
      if(m_p >= 1) {mleparams[0] = smc.mle_params[j]; j++;}
      else{mleparams[0] = 0.0;}
      if(m_p == 2) {mleparams[1] = smc.mle_params[j]; j++;}
      else{mleparams[1] = 0.0;}
      if(m_q >= 1) {mleparams[2] = smc.mle_params[j]; j++;}
      else{mleparams[2] = 0.0;}
      if(m_q == 2) {mleparams[3] = smc.mle_params[j]; j++;}
      else{mleparams[3] = 0.0;}
      if(m_P == 1) {mleparams[4] = smc.mle_params[j]; j++;}      
      else{mleparams[4] = 0.0;}
      if(m_Q == 1) {mleparams[5] = smc.mle_params[j]; j++;}       
      else{mleparams[5] = 0.0;}
      mleparams[6] = smc.mle_params[j];

      textphi1.setText("" + df.format(mleparams[0]));
      textphi2.setText("" + df.format(mleparams[1]));
      texttheta1.setText("" + df.format(mleparams[2]));
      texttheta2.setText("" + df.format(mleparams[3]));
      textTheta.setText("" + df.format(mleparams[5]));
      textPhi.setText("" + df.format(mleparams[4]));
      textInnvar.setText("" + df.format(mleparams[6]));
 
      textLB0.setText("" + df.format(smc.lbv[0]));
      textLB1.setText("" + df.format(smc.lbv[1]));
      textLB12.setText("" + df.format(smc.lbv[2]));
      textLB24.setText("" + df.format(smc.lbv[3]));
 
      textLk1.setText("" + df.format(smc.lkhs[0]));
      textLk2.setText("" + df.format(smc.lkhs[1]));    
      textLk3.setText("" + df.format(smc.lkhs[2]));
      textLk4.setText("" + df.format(smc.lkhs[3]));

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

    public void setDiagnostics()
    {
       //DecimalFormat df = new DecimalFormat("##.###");
       text10.setText("" + df.format(smc.diagnostics[0]));
       text11.setText("" + df.format(smc.diagnostics[1]));
       text12.setText("" + df.format(smc.diagnostics[2]));
       text13.setText("" + df.format(smc.diagnostics[3]));
      
       //System.out.println(smc.sampleQIf[0]);
       spec.setSpectrumData(smc.sampleQIf, 900);
       setMLEs();

       //If x13 on mdfa, send over filters
       if(estimate) setX13SeatsFilters();

    }


    public void setPseudoStatistics()
    {
   
      int j = 0;   
      if(p_p >= 1) {pseudo_params[0] = pseudo.pseudo_params[j]; j++;}
      else{pseudo_params[0] = 0.0;}
      if(p_p == 2) {pseudo_params[1] = pseudo.pseudo_params[j]; j++;}
      else{pseudo_params[1] = 0.0;}
      if(p_q >= 1) {pseudo_params[2] = pseudo.pseudo_params[j]; j++;}
      else{pseudo_params[2] = 0.0;}
      if(p_q == 2) {pseudo_params[3] = pseudo.pseudo_params[j]; j++;}
      else{pseudo_params[3] = 0.0;}
      if(p_P == 1) {pseudo_params[4] = pseudo.pseudo_params[j]; j++;}      
      else{pseudo_params[4] = 0.0;}
      if(p_Q == 1) {pseudo_params[5] = pseudo.pseudo_params[j]; j++;}       
      else{pseudo_params[5] = 0.0;}
      pseudo_params[6] = pseudo.pseudo_params[j];
      
 
      ptextphi1.setText("" + df.format(-pseudo_params[0])); 
      ptextphi2.setText("" + df.format(-pseudo_params[1])); 
      ptexttheta1.setText("" + df.format(-pseudo_params[2])); 
      ptexttheta2.setText("" + df.format(-pseudo_params[3])); 
      ptextTheta.setText("" + df.format(-pseudo_params[5]));
      ptextPhi.setText("" + df.format(-pseudo_params[4]));
      ptextInnvar.setText("" + df.format(pseudo_params[6]));
      ptextMinKL.setText("" + df.format(pseudo.minKL));
      
      ptext10.setText("" + df.format(pseudo.sigex_diags[0]));
      ptext11.setText("" + df.format(pseudo.sigex_diags[1]));
      ptext12.setText("" + df.format(pseudo.sigex_diags[2]));
      ptext13.setText("" + df.format(pseudo.sigex_diags[3]));
      
      //System.out.println(pseudo.efficacies[0] + ", " + pseudo.efficacies[1] + ", " + pseudo.efficacies[2]);
      ptextefft.setText("" + df.format(pseudo.efficacies[0]));
      ptexteffs.setText("" + df.format(pseudo.efficacies[1]));
      ptexteffi.setText("" + df.format(pseudo.efficacies[2]));
      ptexteffti.setText("" + df.format(pseudo.efficacies[3]));
      
      ptextLk1.setText("" + df.format(pseudo.lk[0]));
      ptextLk2.setText("" + df.format(pseudo.lk[1]));
      ptextLk3.setText("" + df.format(pseudo.lk[2]));
      ptextLk4.setText("" + df.format(pseudo.lk[3]));

      ptextLB12.setText("" + df.format(pseudo.lbv[0]));
      ptextLB24.setText("" + df.format(pseudo.lbv[1]));   
      

     }

     public static boolean isSubset(int[] true_d, int[] pseudo_d)
     {
        boolean isSub = true;  
        if((pseudo_d[0] < true_d[0]) || (pseudo_d[2] < true_d[2]) 
        || (pseudo_d[3] < true_d[3]) || (pseudo_d[5] < true_d[5]))
        {isSub = false;}  
        return isSub;
     }


    //-------------------Testing inputs from keyboarded values----------
    public void test_nobs(String s)
    {   
      int i;
      try{i = Integer.parseInt(s.trim()); if(i >=120 && i < 5000) {nObsBar.setValue(i); nObs = i; text2.setText(""+nObs);}
                                          else text2.setText(""+nObs);}
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); text2.setText(""+nObs);}
    }




    public void adjustmentValueChanged(AdjustmentEvent e)
    { 
      int i;
      if(sim_series)
      { 
       if(e.getAdjustable() == nObsBar)
       {
           nObs = nObsBar.getValue();
           text2.setText(""+ nObs);
	   smc.setNobs(nObs);  pseudo.setNobs(nObs);
           initializeUseless(nObs);
  
           if(!pseudoCntrl)
           {
            smc.computeSARIMAmodel(estimate);
            smc.computeDataMax();   
            setDiagnostics();
            smc.go(); spec.go();
           }
           else
           {
             smc.computeSARIMAmodel(false);
             smc.computeDataMax();
             pseudo.setData(nObs, smc.t_series);
             pseudo.computeEfficacies();
             setPseudoStatistics();  
             smc.go(); 
             spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true, pseudo.sampG); spec.go();  
           }
           computeEMD(smc.t_series, nObs); updateScale();
         
           
           if(mdfaPanelx)
           {
              if(mdfa.useSARIMA) {mdfa.setX13Data(smc.t_series, nObs);}
              else {mdfa.setNobs(nObs);}
           }
        
          if(simPanelx){simulate.setNobs(nObs);}
       }
       else if(e.getAdjustable() == seedBar)
       {
           seed = seedBar.getValue();
           text3.setText(""+ seed);
	   smc.setSeed(seed); pseudo.setSeed(seed);

 
           if(!pseudoCntrl)
           {
            smc.computeSARIMAmodel(estimate); 
            smc.computeDataMax(); setDiagnostics(); 
            smc.go(); spec.go();
           }
           else
           {
             smc.computeSARIMAmodel(false);
             smc.computeDataMax();
             pseudo.setData(nObs, smc.t_series);
             pseudo.computeEfficacies();
             setPseudoStatistics();  
             smc.go(); 
             spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
           } 
           
           computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs); mdfa.changeSeed(seed);}
           if(simPanelx){simulate.changeSeed(seed);}
       }
       else if(e.getAdjustable() == innvarParameter)
       {
           i_innvar = innvarParameter.getValue();
           innvar = i_innvar*.01;
           text1.setText(""+ df.format(innvar));
           params[6] = innvar;
           true_params = params; 
	   smc.setInnovation(innvar);
           pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);

           
           if(!pseudoCntrl) 
           {
            smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
            smc.go(); spec.go();
           }
           else
           {   
             smc.computeSARIMAmodel(false);
             smc.computeDataMax();
             pseudo.setData(nObs, smc.t_series);
             pseudo.computeEfficacies();
             setPseudoStatistics();  
             smc.go(); 
             spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
           } 
          
           computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == burnBar)
       {
           burnin = burnBar.getValue();
           text0.setText(""+ burnin);
	   smc.setBurnin(burnin); pseudo.setBurnin(burnin);
        
           if(!pseudoCntrl)
           {
            smc.computeSARIMAmodel(estimate);
            smc.computeDataMax();   
            setDiagnostics();
            smc.go(); spec.go();
           }
           else
           {
             smc.computeSARIMAmodel(false);
             smc.computeDataMax();
             pseudo.setData(nObs, smc.t_series);
             pseudo.computeEfficacies();
             setPseudoStatistics();  
             smc.go(); 
             spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
           }   
          
           computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == phi1bar)
       {
           iphi1 = phi1bar.getValue();
           phi1 = iphi1*.01; sim_phi1 = phi1;
           text4.setText("" + df.format(-phi1));
           m_p = p.getSelectedIndex();
  
          if(!pseudoCntrl) 
          {
           if(m_p > 0)
           { 
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[0] = m_p;
              params[0] = phi1;
              smc.setParams(dims, params, n_params, innvar);
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
              smc.go(); spec.go();
           }
           else //m_p = 0
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[0] = m_p; phi1 = 0.0;
              params[0] = phi1; 
           }
          }
          else
          {
             m_p = tp.getSelectedIndex();
             if(m_p > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[0] = m_p;
               true_params[0] = phi1;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[0] = m_p; phi1 = 0.0;
              true_params[0] = phi1; 
             }
           }
          
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}         
       }
       else if(e.getAdjustable() == phi2bar)
       {
          iphi2 = phi2bar.getValue();
          phi2 = iphi2*.01; sim_phi2 = phi2; 
          text5.setText("" +  df.format(-phi2));
          m_p = p.getSelectedIndex();
        
          if(!pseudoCntrl)
          {
           if(m_p == 2)
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[0] = m_p;
              params[1] = phi2;
              smc.setParams(dims, params, n_params, innvar);               
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
              smc.go(); spec.go();
           }
           else //m_p = 0
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[0] = m_p; //phi2 = 0.0;
              //params[1] = phi2; 
           }
          }
          else
          {
             m_p = tp.getSelectedIndex();
             if(m_p > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[0] = m_p;
               true_params[1] = phi2;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[0] = m_p; phi2 = 0.0;
              true_params[1] = phi2; 
             }
          }
         
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == Phibar)
       {
          iPhi = Phibar.getValue();
          Phi = iPhi*.01; sim_Phi = Phi;
          text6.setText("" +  df.format(-Phi));
          m_P = P.getSelectedIndex();

    
          if(!pseudoCntrl)
          {
           if(m_P == 1)
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[3] = m_P;
              params[4] = Phi;
              smc.setParams(dims, params, n_params, innvar);
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
              smc.go(); spec.go();
           }
           else
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[3] = m_P; //Phi = 0.0;
              //params[2] = Phi; 
           }
          }
          else
          {
             m_P = tP.getSelectedIndex();
             if(m_P > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[3] = m_P;
               true_params[4] = Phi;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[3] = m_P; Phi = 0.0;
              true_params[4] = Phi; 
             }
          } 
         
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == theta1bar)
       {
           itheta1 = theta1bar.getValue();
           theta1 = itheta1*.01; sim_theta1 = theta1;
           text7.setText("" +  df.format(-theta1));
           m_q = q.getSelectedIndex();

        
           if(!pseudoCntrl)
           {
            if(m_q > 0)
            { 
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[2] = m_q;
              params[2] = theta1;
              smc.setParams(dims, params, n_params, innvar);
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
              smc.go(); spec.go();
            }
            else //m_p = 0
            {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[2] = m_q; //theta1 = 0.0;
              //params[3] = theta1; 
            }
          }
          else
          {
             m_q = tq.getSelectedIndex();
             if(m_q > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[2] = m_q;
               true_params[2] = theta1;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[2] = m_q; theta1 = 0.0;
              true_params[2] = theta1; 
             }
          } 
         
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == theta2bar)
       {
           itheta2 = theta2bar.getValue();
           theta2 = itheta2*.01;  sim_theta2 = theta2;
           text8.setText("" +  df.format(-theta2));
           m_q = q.getSelectedIndex();
      
           if(!pseudoCntrl)
           { 
            if(m_q > 1)
            { 
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[2] = m_q;
              params[3] = theta2;
              smc.setParams(dims, params, n_params, innvar);
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics(); 
              smc.go(); spec.go();
            }
            else //m_p = 0
            {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[2] = m_q; //theta2 = 0.0;
              //params[4] = theta2; 
            }
          }
          else
          {
             m_q = tq.getSelectedIndex();
             if(m_q > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[2] = m_q;
               true_params[3] = theta2;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[2] = m_q; theta2 = 0.0;
              true_params[3] = theta2; 
             }
          } 
         
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       else if(e.getAdjustable() == Thetabar)
       {
           iTheta = Thetabar.getValue();
           Theta = iTheta*.01;  sim_Theta = Theta; 
           text9.setText("" +  df.format(-Theta));
           m_Q = Q.getSelectedIndex();
        
           if(!pseudoCntrl)
           { 
           if(m_Q == 1)
           { 
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[5] = m_Q;
              params[5] = Theta;
              //System.out.println(Theta);
              smc.setParams(dims, params, n_params, innvar);
              smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics();  
              smc.go(); spec.go();
           }
           else //m_p = 0
           {
              n_params = m_p + m_P + m_q + m_Q + 1;
              dims[5] = m_Q; //Theta = 0.0;
              //params[5] = Theta; 
           }
          }
          else
          {
             m_Q = tQ.getSelectedIndex();
             if(m_Q > 0)
             {
               n_true = m_p + m_P + m_q + m_Q + 1;
               true_dim[5] = m_Q;
               true_params[5] = Theta;
               pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
               smc.setParams(true_dim, true_params, n_true, innvar);
               smc.computeSARIMAmodel(false);
               smc.computeDataMax();
               pseudo.setData(nObs, smc.t_series);
               pseudo.computeEfficacies();
               setPseudoStatistics();  
               smc.go(); 
               spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go();
             }
             else //m_p = 0
             {
              n_true = m_p + m_P + m_q + m_Q + 1;
              true_dim[5] = m_Q; Theta = 0.0;
              true_params[5] = Theta; 
             }
          } 
         
          computeEMD(smc.t_series, nObs); if(mdfaPanelx){mdfa.setX13Data(smc.t_series, nObs);}
       }
       if(simPanelx){simulate.setSarimaParams(sim_phi1, sim_phi2, sim_theta1, sim_theta2, sim_Phi, sim_Theta);}
      } 
      if(e.getSource() == lagBar)
      {
        //System.out.println(lag);
        lag = lagBar.getValue();
        smc.setLag(lag);         
        textLag.setText("" + lag);
        if(sarimaPanelx)
        {        
         if(!pseudoCntrl)
         {smc.computeDiagnostics(smc.model);
         setDiagnostics();}
         else
         {pseudo.computeEfficacies();
         setPseudoStatistics();}
        }
        if(mdfaPanelx) {mdfa.setLag(lag);}
        if(simPanelx) {simulate.setGlobalLag(lag);}
       }
       if(specCntrl){setSpectralStuff();}
     
       if(e.getSource() == series_scroll)
       {
          
          if(!sim_series)
          {  
             //System.out.print("adjust series, number = ");

             i_val = series_scroll.getValue();
             //System.out.println(i_val);
             if(i_val >= 0 && i_val < n_series)
             { 
                series_name.setText(list_files[i_val]+" "+i_val+"/"+n_series);                
                //change real_series 
                nObs = gs_lengths[i_val];
                real_series = new double[nObs];
                
                for(i=0; i < nObs; i++)
                {real_series[i] = group_series[i_val][i];}     

                //--------baraaye hallaaji karde shodan khoub e------------                
                initializeUseless(nObs);  
                smc.setNobs(nObs);                                
                smc.setRealData(real_series, nObs);
                smc.computeSARIMAmodel(estimate);
                smc.computeDataMax(); setDiagnostics(); 
                if(!specCntrl) smc.go();     
                else{spec.go();}  
                //---------------------------------------------------

                if(c_emd){computeEMD(real_series, nObs);}
                if(mdfaPanelx){mdfa.setX13Data(real_series, nObs);}
                reg.inputX13data(real_series); BayesCronosPan.setTimeSeries(real_series, 1, nObs);      
                
                orig_series = new double[real_series.length]; System.arraycopy(real_series,0,orig_series,0,nObs);
             }
          }
       }
       mdfa.inputX13data(smc.t_series); 

    }
    //================= Action performed sequences ====================
    @SuppressWarnings("rawtypes")
	public void actionPerformed(ActionEvent e)
    {

      if(e.getSource() == p  || e.getSource() == tp)
      {                         
         m_p = ((JComboBox)e.getSource()).getSelectedIndex();
         dims[0] = m_p;
         if(m_p == 0)
         {phi1 = 0.0; phi2 = 0.0; params[0] = phi1; params[1] = phi2;}
         if(m_p == 1)
         {phi1 = .01*phi1bar.getValue(); params[0] = phi1;}
         else if(m_p == 2)
         {phi1 = .01*phi1bar.getValue(); phi2 = .01*phi2bar.getValue(); params[0] = phi1; params[1] = phi2;}         
         n_params = m_p + m_P + m_q + m_Q + 1; 
         smc.setParams(dims, params, n_params, innvar);      

         n_true = n_params; true_dim = dims; 
         true_params = params; 
         pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
          
         if(!pseudoCntrl)
         {smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);}   
      
         if(!sim_series){
         smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics();  
         smc.go(); spec.go();} 

       
      }
      else if(e.getSource() == P  || e.getSource() == tP)
      {                         
         m_P = ((JComboBox)e.getSource()).getSelectedIndex();
         dims[3] = m_P;
         if(m_P == 0)
         {Phi = 0.0; params[4] = Phi;}
         if(m_P == 1)
         {Phi = .01*Phibar.getValue(); params[4] = Phi;}
         n_params = m_p + m_P + m_q + m_Q + 1;
         smc.setParams(dims, params, n_params, innvar);

         n_true = n_params; true_dim = dims; 
         true_params = params; 
         pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
         smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);       
 
         if(!pseudoCntrl)
         {smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);}      

         if(!sim_series){
         smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics();  
         smc.go(); spec.go();}
 
      }
      else if(e.getSource() == q || e.getSource() == tq)
      {                         
         m_q = ((JComboBox)e.getSource()).getSelectedIndex();
         //System.out.println("m_q = " + m_q);
         dims[2] = m_q;
         if(m_q == 0)
         {theta1 = 0.0; theta2 = 0.0; params[2] = theta1; params[3] = theta2;}
         if(m_q == 1)
         {theta1 = .01*theta1bar.getValue(); params[2] = theta1;}
         else if(m_q == 2)
         {theta1 = .01*theta1bar.getValue(); theta2 = .01*theta2bar.getValue(); params[2] = theta1; params[3] = theta2;}         
         n_params = m_p + m_P + m_q + m_Q + 1;

         smc.setParams(dims, params, n_params, innvar);      
         n_true = n_params; true_dim = dims; 
         true_params = params; 
         
         pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
         
         if(!pseudoCntrl)
         {smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);}          

         if(!sim_series){
         smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics();  
         smc.go(); spec.go();}      
 

      }
      else if(e.getSource() == Q || e.getSource() == tQ)
      {                         
         m_Q = ((JComboBox)e.getSource()).getSelectedIndex();
         dims[5] = m_Q;
         if(m_Q == 0)
         {Theta = 0.0; params[5] = Theta;}
         if(m_Q == 1)
         {Theta = .01*Thetabar.getValue(); params[5] = Theta;} 
         n_params = m_p + m_P + m_q + m_Q + 1;
         smc.setParams(dims, params, n_params, innvar);      

         n_true = n_params; true_dim = dims; 
         true_params = params; 
         pseudo.setTrueParameters(n_true, true_dim, true_params, innvar);
                    
         if(!pseudoCntrl)
         {smc.modelDimensionChange(dims, trans,outlier,easter, td, easterDays);}   

         if(!sim_series){
         smc.computeSARIMAmodel(estimate); smc.computeDataMax(); setDiagnostics();  
         smc.go(); spec.go();}
            
      }

      // ----- pseudo controls -----------------------------
      else if(e.getSource() == pp)
      {    
          p_p = pp.getSelectedIndex(); pseudo_dim[0] = p_p; 
          n_pseudo = p_p + p_P + p_q + p_Q + 1;
          smc.modelDimensionChange(pseudo_dim, trans,outlier,easter, td, easterDays); 
          pseudo.setPseudoParameters(pseudo_dim);
          pseudo.computeEfficacies();
          pseudo.setSpectrumDataLimits();
          setPseudoStatistics();
          spec.go();
         
          if(isSubset(true_dim, pseudo_dim))
          {System.out.println("Warning: True-model is submodel of Pseudo-model");}
          if(!sim_series)
          {smc.computeSARIMAmodel(false);}  

      }
          
      else if(e.getSource() == pP)
      {                         
          p_P = pP.getSelectedIndex(); pseudo_dim[3] = p_P; 
          n_pseudo = p_p + p_P + p_q + p_Q + 1;
          smc.modelDimensionChange(pseudo_dim, trans,outlier,easter, td, easterDays);
          pseudo.setPseudoParameters(pseudo_dim);
          pseudo.computeEfficacies();
          pseudo.setSpectrumDataLimits();
          setPseudoStatistics();
          spec.go();
         
          if(isSubset(true_dim, pseudo_dim))
          {System.out.println("Warning: True-model is submodel of Pseudo-model");}

          if(!sim_series)
          {smc.computeSARIMAmodel(false);}  
      }
      else if(e.getSource() == pq)
      {                         
          p_q = pq.getSelectedIndex(); pseudo_dim[2] = p_q; 
          n_pseudo = p_p + p_P + p_q + p_Q + 1;
          smc.modelDimensionChange(pseudo_dim, trans,outlier,easter, td, easterDays);
          pseudo.setPseudoParameters(pseudo_dim);
          pseudo.computeEfficacies();
          pseudo.setSpectrumDataLimits();
          setPseudoStatistics();
          spec.go();
         
          if(isSubset(true_dim, pseudo_dim))
          {System.out.println("Warning: True-model is submodel of Pseudo-model");}

          if(!sim_series)
          {smc.computeSARIMAmodel(false);}  
      }
      else if(e.getSource() == pQ)
      {                         
          p_Q = pQ.getSelectedIndex(); pseudo_dim[0] = p_Q; 
          n_pseudo = p_p + p_P + p_q + p_Q + 1;
          smc.modelDimensionChange(pseudo_dim, trans,outlier,easter, td, easterDays);
          pseudo.setPseudoParameters(pseudo_dim);
          pseudo.computeEfficacies();
          pseudo.setSpectrumDataLimits();
          setPseudoStatistics();
          spec.go();
         
          if(isSubset(true_dim, pseudo_dim))
          {System.out.println("Warning: True-model is submodel of Pseudo-model");}

          if(!sim_series)
          {smc.computeSARIMAmodel(false);}  
      }   
      //-----------------combo box for easter-----------------------------
      else if(e.getSource() == easterDay) 
      { 
         int indx = ((JComboBox)e.getSource()).getSelectedIndex();
         if(indx == 0)
         {easterDays = 1;}
         else if(indx == 1)
         {easterDays = 8;} 
         else if(indx == 2)
         {easterDays = 15;}        
      }

    }
    
    //--------- For the itemstate changed checkboxes ------
    public void itemStateChanged(ItemEvent e) 
    {
      boolean sel = true;
      Object source = e.getItemSelectable();

      if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
      if(source == seriesBox){smc.setSeries(sel);}         
      else if(source == foreBox) {smc.setForecasts(sel);}
      else if(source == seasBox) {smc.setSeas(sel);}
      else if(source == trendBox) {smc.setTrend(sel);}
      else if(source == saBox) {smc.setSA(sel);}
      else if(source == cycleBox) {smc.setCycle(sel);}
      else if(source == TrendModel) 
      {
         if(!pseudoCntrl)
         {smc.changeModel(6); smc.computeDiagnostics(6); setDiagnostics(); spec.go();}
         else
         {
          pseudo.setModel(6); pseudo.computeEfficacies(); setPseudoStatistics();
          spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go(); 
         }
 
      }
      else if(source == SeasonalModel)
      {
         if(!pseudoCntrl)
         {smc.changeModel(1); smc.computeDiagnostics(1); setDiagnostics(); spec.go();}
         else
         {
          pseudo.setModel(1); pseudo.computeEfficacies(); setPseudoStatistics();
          spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go(); 
         }    
      }
      else if(source == IrregModel) 
      {
        if(!pseudoCntrl)
        {smc.changeModel(3); smc.computeDiagnostics(3); setDiagnostics(); spec.go();}
         else
         {
          pseudo.setModel(3); pseudo.computeEfficacies(); setPseudoStatistics();
          spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go(); 
         }
      }
      else if(source == TrendIrregModel) 
      {
        if(!pseudoCntrl)
        {smc.changeModel(2);  smc.computeDiagnostics(2); setDiagnostics(); spec.go();}
         else
         {
          pseudo.setModel(2); pseudo.computeEfficacies(); 
          setPseudoStatistics();
          spec.setSpectralData(pseudo.sampQf_pseudo, pseudo.sampQf_true,pseudo.sampG); spec.go(); 
         }
      }
      else if(source == plotInBox)
      {plotIn = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}                    
      else if(source == plotFwBox)
      {plotFw = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}
      else if(source == plotDnBox)
      {plotDn = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}
      else if(source == plotFGBox)
      {plotGF = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}
      else if(source == plotIGBox)
      {plotGI = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}
      else if(source == plotGBox)
      {plotG = sel; spec.setSpectralPlots(plotIn, plotFw, plotDn, plotGF, plotGI, plotG);  if(sel){setSpectralStuff();}}       
      
      else if(source == sim_check)
      {sim_series = sel; smc.setSimSeries(sel); enableX13Simulator(sel);}

      else if(source == sim_mdfa_check)
      {mdfa.setSimulate(sel);}

      else if(source == transBox)
      {trans = sel; x13kardan();}
      else if(source == tdBox)
      {td = sel;  x13kardan();  }
      else if(source == easterBox)
      {easter = sel; x13kardan();  }
      else if(source == outlierBox)
      {outlier = sel; x13kardan();  }
      
      smc.go(); spec.go();
    }


    public void x13kardan()
    {
       //1) ------------ Change x13cpp file--------------------------
       smc.modelDimensionChange(dims, trans, outlier,easter, td, easterDays);

       //2) ------------ Compute new model---------------------------
       smc.computeSARIMAmodel(estimate); 
       smc.computeDataMax(); 

       setDiagnostics();

    }


    public void setSpectralStuff()
    {spec.setSpectralFuncs(smc.sarima.sampleIn, smc.sarima.sampleFw, smc.sarima.sampleDn);}

    public void readData(File file)
    {
          
       String strline; Double D;
       String[] tokens; String delims = "[ ]+";
       int n_toks;
       double[] values = new double[1500];     
       
       double val = 0;
       int i = 0;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           

           D = new Double(tokens[n_toks-1]); //take only the value, no dates
           val = D.doubleValue();

           if(i < 1500) {values[i] = val; i++;}
           else 
           {System.out.println("Maximum times series length is 1500"); break;} 
           
         } 
         din.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}

        // ------ Update the inputed data--------------------------------------
        nObs = i; 
        real_series = new double[nObs]; 
        for(i=0; i < nObs; i++)
        {real_series[i] = values[i];}//System.out.println(real_series[i]);}
        normalizeRealData();

        sim_series = false; initializeUseless(nObs);  smc.setNobs(nObs); 
        smc.setSimSeries(sim_series); 
        sim_check.setSelected(sim_series);
        smc.setRealData(real_series, nObs);
        smc.computeSARIMAmodel(estimate);
        smc.computeDataMax(); setDiagnostics(); smc.go(); if(specCntrl){spec.go();}        

        n_series = 1;
        group_series = new double[n_series][500];   //--------- data from series -----------
        gs_lengths = new int[n_series];             //--------- lengths of series ----------
        list_files = new String[n_series];          // ----- file names 
       
        gs_lengths[0] = nObs;
        for(i=0; i < gs_lengths[0]; i++)
        {group_series[0][i] = real_series[i];}
        list_files[0] = file.getName();        

        series_scroll.setMaximum(0);
        series_scroll.setValue(0);
        series_name.setText(file.getName()+" "+1+"/"+n_series);  
      
        //--- input into I_MDFA
        mdfa.inputX13data(real_series);
        BayesCronosPan.setTimeSeries(real_series, 1, nObs);    
        reg.inputX13data(real_series);
        orig_series = new double[real_series.length]; System.arraycopy(real_series,0,orig_series,0,nObs);

    }


    public void normalizeRealData()
    {
       int i;
       double p, norma;  
       double min = 10000000000000000000000000.0;
     
       for(i=0;i<nObs;i++)
       {  
          if(real_series[i] < min) {min = real_series[i];}
       }
    
       p = Math.floor(Math.log10(min)); 
       if(p > 3.0)
       {
          norma = Math.pow(10,p-1);           
          for(i=0;i<nObs;i++) {real_series[i] = real_series[i]/norma;}    
       }
    } 

    public void saveFilterCoeffs()
    {
       int i,j;
       FileFilter ft = new FileNameExtensionFilter("Coefficient files", "cfs");
       fc.addChoosableFileFilter(ft);
       
       int returnVal = fc.showSaveDialog(this);
       int L = mdfa.L; int n_rep = mdfa.n_rep;

       if(returnVal == javax.swing.JFileChooser.APPROVE_OPTION && mdfaPanelx)
       {
         File saved_file = fc.getSelectedFile();
         String file_name = saved_file.toString();
         
         try
         {
           PrintWriter out = new PrintWriter(new FileWriter(file_name));  

           //--- print L then n_reps ------
           out.println(mdfa.L + " " + mdfa.n_rep);

          //          for(i=0; i < n_rep; i++)  
        // {
         //  for(k=0;k<L;k++)
        //   {b[L*i + k] = out[lag4 + L*i + k];}
        // }

           for(i=0;i<L;i++)
           {
             for(j=0;j<n_rep-1;j++)
             {
               out.print(mdfa.mdfa.b[L*j+i] + " "); 
             }
             out.println(mdfa.mdfa.b[L*(n_rep-1)+i]);
           }  
           
           out.close(); System.out.println("Filter successfully saved in " + file_name);
         }
         catch (IOException e) {e.printStackTrace();}
       } 
    }


  

    public void saveFilterFile()
    {
       //open filter dialog
       FileFilter ft = new FileNameExtensionFilter("Filter files", "flt");
       fc.addChoosableFileFilter(ft);
       
       int returnVal = fc.showSaveDialog(this);
       DecimalFormat df3 = new DecimalFormat("##.###");
       
       if(returnVal == javax.swing.JFileChooser.APPROVE_OPTION && mdfaPanelx)
       {
         File saved_file = fc.getSelectedFile();
         String file_name = saved_file.toString();
         
         try
         {
           PrintWriter out = new PrintWriter(new FileWriter(file_name));           
           String allText = "";
           
           //---- Filter ---
           if(!mdfa.highBut.isSelected()) 
           {
             allText = allText + df.format(mdfa.w0) + " " + df.format(mdfa.w1) + "\n";
           } 
           else
           {allText = allText + df.format(mdfa.w0) + " " + df.format(mdfa.w1) + " " + df.format(mdfa.w2) + " " + df.format(mdfa.w3) + "\n";}

           allText = allText + df.format(mdfa.lambda) + "\n";
           allText = allText + df.format(mdfa.expweight) + "\n";
           allText = allText + df3.format(mdfa.smooth) + "\n";
           allText = allText + df3.format(mdfa.decay) + "\n";           
           allText = allText + df3.format(mdfa.decay2) + "\n";
           allText = allText + df3.format(mdfa.cross) + "\n";
           allText = allText + df.format(mdfa.L) + "\n";
           allText = allText + Integer.parseInt(""+mdfa.Lag) + "\n";
           
           out.print(allText);           
           out.close(); System.out.println("Filter successfully saved in " + file_name);
           
         }
         catch (IOException e) {e.printStackTrace();}
         
         
       }   
    }
    
    public void saveFilterFileEncrypt()
    {
       //open filter dialog
       int i;
       Random generator = new Random((int)System.currentTimeMillis()/100000);
       
       FileFilter ft = new FileNameExtensionFilter("Filter files", "flt");
       fc.addChoosableFileFilter(ft);
       
       int returnVal = fc.showSaveDialog(this);
       DecimalFormat df3 = new DecimalFormat("##.###");
       
       if(returnVal == javax.swing.JFileChooser.APPROVE_OPTION && mdfaPanelx)
       {
         File saved_file = fc.getSelectedFile();
         String file_name = saved_file.toString();
         
         try
         {
           PrintWriter out = new PrintWriter(new FileWriter(file_name));           
           String allText = "";
           
                   
           //---- Filter ------------------ 3
           allText = allText + df3.format(mdfa.w0) + " ";
           for(i=1;i<=2;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";} 
           allText = allText + df3.format(mdfa.w1) + " ";
           for(i=4;i<=6;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";}
           allText = allText + df3.format(generator.nextDouble()) + "\n";
           
           
           //---------- lambda-----------0
           allText = allText + df3.format(mdfa.lambda) + " ";
           for(i=1;i<=6;i++) {allText = allText + df3.format(generator.nextDouble()*10.0) + " ";}
           allText = allText + df3.format(generator.nextDouble()*10.0) + "\n";
        
           //---------- exweight----------1
           allText = allText + df3.format(generator.nextDouble()*10.0) + " " + df3.format(mdfa.expweight) + " ";
           for(i=2;i<=6;i++) {allText = allText + df3.format(generator.nextDouble()*10.0) + " ";}
           allText = allText + df3.format(generator.nextDouble()*10.0) + "\n";        
        
           //-------- smooth ------------ 2
           allText = allText + df3.format(generator.nextDouble()) + " " + df3.format(generator.nextDouble()) + " " + df3.format(mdfa.smooth) + " ";
           for(i=3;i<=6;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";}
           allText = allText + df3.format(generator.nextDouble()) + "\n";  
           
           //-------- decay ------------ 5
           for(i=0;i<=4;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";}
           allText = allText + df3.format(mdfa.decay) + " ";
           allText = allText + df3.format(generator.nextDouble()) + " " + df3.format(generator.nextDouble()) + "\n";

           //-------- decay ------------ 6
           for(i=0;i<=5;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";}
           allText = allText + df3.format(mdfa.decay2) + " ";
           allText = allText + df3.format(generator.nextDouble()) + "\n";
           
           //-------- decay ------------ 7
           for(i=0;i<=6;i++) {allText = allText + df3.format(generator.nextDouble()) + " ";}
           allText = allText + df3.format(mdfa.cross) + "\n";
        
                    
           //--------------------------- 2  1373349600 1373522400 1373608800 259200 172800
           int shift = generator.nextInt(10);
           allText = allText + shift;
           for(i=0;i<shift;i++) {allText = allText + generator.nextInt(10);}
           
           String str = ""+System.currentTimeMillis();
           //System.out.println(str);
           char[] charArray = str.toCharArray();
           
           for(i=2;i<charArray.length;i++)
           {allText=allText+charArray[i];}
           
           //allText = allText + System.currentTimeMillis());
           for(i=0;i<5;i++) {allText = allText + generator.nextInt(10);}
           allText = allText + "\n";
           
           //---- L---------------------- 1
           allText = allText + generator.nextInt(100) + " " + mdfa.L + " "; 
           for(i=2;i<=6;i++) {allText = allText + generator.nextInt(100) + " ";} 
           allText = allText + generator.nextInt(100) + "\n";

           //---- display Lag ----------- 6          
           for(i=0;i<=5;i++) {allText = allText + generator.nextInt(5) + " ";} 
           allText = allText + mdfa.Lag + " " + generator.nextInt(5) + "\n";           
           
           //byte[] buffer = allText.getBytes();
           
           out.print(allText);           
           //output = new BufferedOutputStream(new FileOutputStream(file_name));
           //output.write(buffer);
           //output.close();
           out.close(); 
           System.out.println("Filter successfully saved in " + file_name);
           
         }
         catch (IOException e) {e.printStackTrace();}
         
         
       }   
    }

    //==========  Read Meta Data file in =====

    public void readMetaData(File file)
    {
  

       String strline; 
       String[] _list_files = new String[50];
       File file2; 
       Double D;
       int j = 0;  double val = 0; int i = 0;

       int n_toks; 
       String delims = "[ ]+";
       String[] tokens;
       FileInputStream fin;
       DataInputStream din;
       BufferedReader br;
     
       //file = new File("/home/lisztian/uSimDev/x13/Foreign/Foreign.dta");
       String pathFile = file.getAbsolutePath();
       String name = file.getName();
       String dir_path = pathFile.replaceAll(name, "");     
       //System.out.println(dir_path);

       //---------- Get all the files first-----------------------
       try{
          
         fin = new FileInputStream(file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {
            if(i < 50) {_list_files[i] = strline; i++;}
            else
            {System.out.println("Maximum number of files is 50. Please reduce amount.");}
            
         }
         din.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found...." + fe);}
       catch(IOException ioe){System.out.println("ooops... File not found" + ioe);}

       n_series = i;
       group_series = new double[n_series][500];   //--------- data from series -----------
       gs_lengths = new int[n_series];             //--------- lengths of series ----------
       list_files = new String[n_series];          // ----- file names --------------------

       for(i=0; i < n_series; i++)
       {  
         //--- open i_th file for reading ----
         list_files[i] = _list_files[i]; j = 0;        
         try{
           
            name = dir_path.concat(list_files[i]);

            file2 = new File(name);
            fin = new FileInputStream(file2);
            din = new DataInputStream(fin);
            br = new BufferedReader(new InputStreamReader(din));
  
            while((strline = br.readLine()) != null)
            {
              tokens = strline.split(delims); 
              n_toks = tokens.length; 
              if(n_toks == 0)
              {System.out.println("End of file"); break;}
  
              D = new Double(tokens[n_toks-1]);
              val = D.doubleValue();
              if(j < 500) {group_series[i][j] = val; j++;}
              else 
              {System.out.println("Maximum times series length is 500"); break;} 
           
            } 
            din.close();
          }
          catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
          catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}
            
          gs_lengths[i] = j;
        }
 
        // --- now setup all the files -----
        updateGroupData();
    }


    public void updateGroupData()
    {
        int i,n;
        double p, norma;  
        double min;
        
        //System.out.println("print no_series");
        //System.out.println(n_series);
        series_scroll.setMinimum(0); 
        series_scroll.setMaximum(n_series-1);
        series_scroll.setValue(0);

        //----- normalize data--------------------
        for(n = 0; n < n_series; n++)
        {
           min = 10000000000000000000000000.0;
           for(i = 0; i < gs_lengths[n]; i++) 
           {if(group_series[n][i] < min) {min = group_series[n][i];} }
       
           p = Math.floor(Math.log10(min)); 
           if(p > 3.0)
           {
              norma = Math.pow(10,p-1);           
              for(i=0;i<gs_lengths[n];i++) {group_series[n][i] = group_series[n][i]/norma;}    
           }
        }
        //-----------------------------------------
        
        
        //------ Begin with first series in sequence---------
        nObs = gs_lengths[0]; 
        real_series = new double[nObs];
        for(i=0; i < nObs; i++)
        {real_series[i] = group_series[0][i];}     

        //--------baraaye hallaaji karde shodan khoub e------------
        sim_series = false; 
        initializeUseless(nObs);  
        smc.setNobs(nObs); 
        smc.setSimSeries(sim_series); 
        sim_check.setSelected(sim_series);
        smc.setRealData(real_series, nObs);
        smc.computeSARIMAmodel(estimate);
        smc.computeDataMax(); setDiagnostics(); smc.go();     
     
        //---------------------------------------------------
        series_name.setText(list_files[0]+" "+1+"/"+n_series);
  
        mdfa.inputX13data(real_series);
        reg.inputX13data(real_series);
        BayesCronosPan.setTimeSeries(real_series, 1, nObs);    
        
        orig_series = new double[real_series.length]; System.arraycopy(real_series,0,orig_series,0,nObs);

    }

    public void inputFromSimulator(double[] v)
    {
       
       tsim_data = true; 
       nObs = v.length; //System.out.println("nObs = " + nObs + " v.length = " + v.length);
       double[] sym_data = new double[v.length];
       double[] sim_data = new double[nObs];
          
       System.arraycopy(v, 0, sym_data, 0, v.length);
       System.arraycopy(v, 0, sim_data, 0, nObs);

        sim_series = false; 
        initializeUseless(nObs);  
        smc.setNobs(nObs); 
        smc.setSimSeries(sim_series); 
        sim_check.setSelected(sim_series);
        smc.setRealData(sim_data, nObs);
        x13kardan();
        smc.computeSARIMAmodel(estimate);
        smc.computeDataMax(); setDiagnostics(); smc.go();         
        
        //---------------------------------------------------
        series_name.setText("Simuated 1");
  
        mdfa.inputX13data(sim_data);
        reg.inputX13data(sim_data);      
        BayesCronosPan.setTimeSeries(sim_data, 1, nObs);  

    }

    public static void main(String args[])
    {


       Toolkit.getDefaultToolkit().setDynamicLayout(true);
       System.setProperty("sun.awt.noerasebackground", "true");
       //JFrame.setDefaultLookAndFeelDecorated(true);
       JDialog.setDefaultLookAndFeelDecorated(true);

       try {
            UIManager.setLookAndFeel("de.muntjak.tinylookandfeel.TinyLookAndFeel");
       } catch(Exception ex) {
                 ex.printStackTrace();
       }

       //try{UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());} 
       //try{UIManager.setLookAndFeel("javax.swing.plaf.nimbus.NimbusLookAndFeel");}
       //catch (UnsupportedLookAndFeelException e) {}
       //catch (ClassNotFoundException e) {}
       //catch (InstantiationException e) {}
       //catch (IllegalAccessException e) {}


       JFrame f = new IMetricaProgram("iMetrica - Copyright 2009-2016 by C. Blakely");
       f.setDefaultCloseOperation(EXIT_ON_CLOSE);
       f.pack();
       f.setSize(new Dimension(1160, 775));
       f.setVisible(true);
    }



   public GroupLayout setUpModelLayout(JPanel pmodelPanel)
   {
     // 2) ================= Parameter control panel====================

       JLabel pLabel = new JLabel(" p:");
       JLabel qLabel = new JLabel(" q:");
       JLabel QLabel = new JLabel(" Q:");
       JLabel PLabel = new JLabel(" P:");
       JLabel pLabel2 = new JLabel(" p:");
       JLabel qLabel2 = new JLabel(" q:");
       JLabel QLabel2 = new JLabel(" Q:");
       JLabel PLabel2 = new JLabel(" P:");


       GroupLayout pLayout = new GroupLayout(pmodelPanel);
       pLayout.setAutoCreateGaps(true);
       pLayout.setAutoCreateContainerGaps(true);
       pLayout.setHorizontalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup() 
                 .addComponent(estLabel1)
                 .addComponent(estLabel2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(pLabel)
                 .addComponent(pLabel2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(tp)
                 .addComponent(pp))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(qLabel)
                 .addComponent(qLabel2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(tq)
                 .addComponent(pq))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(PLabel)
                 .addComponent(PLabel2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(tP)
                 .addComponent(pP))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(QLabel)
                 .addComponent(QLabel2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(tQ)
                 .addComponent(pQ)));
              


       pLayout.setVerticalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(pLayout.createSequentialGroup()
                  .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(estLabel1)
                   .addComponent(pLabel)
                   .addComponent(tp)
                   .addComponent(qLabel)
                   .addComponent(tq)
                   .addComponent(PLabel) 
                   .addComponent(tP)
                   .addComponent(QLabel)
                   .addComponent(tQ))                 
                 .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE) 
                   .addComponent(estLabel2)
                   .addComponent(pLabel2)
                   .addComponent(pp)
                   .addComponent(qLabel2)
                   .addComponent(pq)
                   .addComponent(PLabel2) 
                   .addComponent(pP)
                   .addComponent(QLabel2)
                   .addComponent(pQ)))));         

     return pLayout;
   }



   public GroupLayout setUpMLEValuesLayout(JPanel vPanel)
   {
     // 2) ================= Parameter control panel====================
      Font f = new Font("Dialog", Font.PLAIN, 20);
      JLabel lphi1 = new JLabel("\u03D5_1");
      JLabel lphi2 = new JLabel("\u03D5_2");
      JLabel ltheta1  = new JLabel("\u03D1_1"); 
      JLabel ltheta2  = new JLabel("\u03D1_2");
      JLabel lTheta = new JLabel("\u0398");
      lTheta.setFont(f);
      JLabel lPhi = new JLabel("\u03A6");
      lPhi.setFont(f);
      JLabel linnvar = new JLabel("Inn.Var.");
      JLabel lsigvar = new JLabel("Sig.Var.");
      

       GroupLayout pLayout = new GroupLayout(vPanel);
       pLayout.setAutoCreateGaps(true);
       pLayout.setAutoCreateContainerGaps(true);
       pLayout.setHorizontalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup() 
                 .addComponent(lphi1)
                 .addComponent(lphi2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(textphi1)
                 .addComponent(textphi2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ltheta1)
                 .addComponent(ltheta2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(texttheta1)
                 .addComponent(texttheta2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(lPhi)
                 .addComponent(lTheta))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(textPhi)
                 .addComponent(textTheta))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(linnvar)
                 .addComponent(lsigvar))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(textInnvar)
                 .addComponent(textSigvar)));
              
       pLayout.setVerticalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(pLayout.createSequentialGroup()
                  .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lphi1)
                   .addComponent(textphi1)
                   .addComponent(ltheta1)
                   .addComponent(texttheta1)
                   .addComponent(lPhi)
                   .addComponent(textPhi) 
                   .addComponent(linnvar)
                   .addComponent(textInnvar))                 
                 .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE) 
                   .addComponent(lphi2)
                   .addComponent(textphi2)
                   .addComponent(ltheta2)
                   .addComponent(texttheta2)
                   .addComponent(lTheta)
                   .addComponent(textTheta) 
                   .addComponent(lsigvar)
                   .addComponent(textSigvar)))));         

                   return pLayout;
   }

   public GroupLayout setUpPseudoValuesLayout(JPanel vPanel)
   {
     // 2) ================= Parameter control panel====================
      Font f = new Font("Dialog", Font.PLAIN, 20);
      JLabel lphi1 = new JLabel("\u03D5_1");
      JLabel lphi2 = new JLabel("\u03D5_2");
      JLabel ltheta1  = new JLabel("\u03D1_1"); 
      JLabel ltheta2  = new JLabel("\u03D1_2");
      JLabel lTheta = new JLabel("\u0398");
      lTheta.setFont(f);
      JLabel lPhi = new JLabel("\u03A6");
      lPhi.setFont(f);
      JLabel linnvar = new JLabel("Inn.Var."); linnvar.setToolTipText("ML estimation of innovation variance.");
      JLabel lminKL = new JLabel("MinKL."); lminKL.setToolTipText("Minimum Value of Likelihood.");
      

       GroupLayout pLayout = new GroupLayout(vPanel);
       pLayout.setAutoCreateGaps(true);
       pLayout.setAutoCreateContainerGaps(true);
       pLayout.setHorizontalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup() 
                 .addComponent(lphi1)
                 .addComponent(lphi2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ptextphi1)
                 .addComponent(ptextphi2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ltheta1)
                 .addComponent(ltheta2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ptexttheta1)
                 .addComponent(ptexttheta2))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(lPhi)
                 .addComponent(lTheta))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ptextPhi)
                 .addComponent(ptextTheta))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(linnvar)
                 .addComponent(lminKL))
              .addGroup(pLayout.createParallelGroup()
                 .addComponent(ptextInnvar)
                 .addComponent(ptextMinKL)));
              
       pLayout.setVerticalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(pLayout.createSequentialGroup()
                  .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lphi1)
                   .addComponent(ptextphi1)
                   .addComponent(ltheta1)
                   .addComponent(ptexttheta1)
                   .addComponent(lPhi)
                   .addComponent(ptextPhi) 
                   .addComponent(linnvar)
                   .addComponent(ptextInnvar))                 
                 .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE) 
                   .addComponent(lphi2)
                   .addComponent(ptextphi2)
                   .addComponent(ltheta2)
                   .addComponent(ptexttheta2)
                   .addComponent(lTheta)
                   .addComponent(ptextTheta) 
                   .addComponent(lminKL)
                   .addComponent(ptextMinKL)))));         

                   return pLayout;
   }


   public GroupLayout setUpControlLayout(JPanel vPanel)
   {

       JLabel obsL = new JLabel("Obs:"); obsL.setToolTipText("Number of observation for simulation.");
       JLabel burnL = new JLabel("Burnin:"); burnL.setToolTipText("Number of burnin points for simulation of model realization.");
       JLabel seedL = new JLabel("Seed:"); seedL.setToolTipText("Seed for random number generator of innovation process.");
       JLabel lagL = new JLabel("Lag:"); lagL.setToolTipText("Lag used for evaluating different diagnostics. Lags 0,1,12, most commonly used.");


       nObsBar = new JScrollBar(JScrollBar.HORIZONTAL,144,10,120,15360);
       nObsBar.setUnitIncrement(1);      
       text2 = new JTextField(3); 
       text2.setText(""+nObs);
       nObsBar.addAdjustmentListener(this);
       text2.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_nobs(text2.getText());}} );  
        
       burnBar = new JScrollBar(JScrollBar.HORIZONTAL,100,10,100,600);
       burnBar.setUnitIncrement(10);
       burnBar.addAdjustmentListener(this);
       text0 = new JTextField(3);
      
       text0.setText("" + burnin);       
       seedBar = new JScrollBar(JScrollBar.HORIZONTAL,500,10,0,1000);
       seedBar.setUnitIncrement(1);

       text3 = new JTextField(3);
       text3.setText(""+seed);
       seedBar.addAdjustmentListener(this);

       lagBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,24);
       lagBar.setUnitIncrement(1);
       textLag = new JTextField(2);
       textLag.setText(""+lag);
       lagBar.addAdjustmentListener(this);


     // 2) ================= Parameter control panel====================
       GroupLayout contLayout = new GroupLayout(vPanel);
       contLayout.setAutoCreateGaps(true);
       contLayout.setAutoCreateContainerGaps(true);
       contLayout.setHorizontalGroup(contLayout.createSequentialGroup()
              .addGroup(contLayout.createParallelGroup() 
                 .addComponent(obsL)
                 .addComponent(seedL))
              .addGroup(contLayout.createParallelGroup()
                 .addComponent(nObsBar)
                 .addComponent(seedBar))
              .addGroup(contLayout.createParallelGroup()
                 .addComponent(text2)
                 .addComponent(text3))
              .addGroup(contLayout.createParallelGroup()
                 .addComponent(burnL)
                 .addComponent(lagL))
              .addGroup(contLayout.createParallelGroup()
                 .addComponent(burnBar)
                 .addComponent(lagBar))
              .addGroup(contLayout.createParallelGroup()
                 .addComponent(text0)
                 .addComponent(textLag)));


       contLayout.setVerticalGroup(contLayout.createSequentialGroup()
              .addGroup(contLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(contLayout.createSequentialGroup()
                  .addGroup(contLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(obsL)
                   .addComponent(nObsBar)
                   .addComponent(text2)
                   .addComponent(burnL)
                   .addComponent(burnBar)
                   .addComponent(text0)) 
                 .addGroup(contLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                  .addComponent(seedL)
                  .addComponent(seedBar)
                  .addComponent(text3)
                  .addComponent(lagL)
                  .addComponent(lagBar)
                  .addComponent(textLag)))));  

                  return contLayout;  
   }

   public void matchSelectedIndices(int sel)
   {
     int _p, _q, _Q, _P;
     if(sel == 1)
     {
      _p = p.getSelectedIndex();
      _q = q.getSelectedIndex();
      _P = P.getSelectedIndex();
      _Q = Q.getSelectedIndex();
      
      tp.setSelectedIndex(_p);
      tq.setSelectedIndex(_q);
      tP.setSelectedIndex(_P);
      tQ.setSelectedIndex(_Q); 
     }
     else
     {
      p.setSelectedIndex(tp.getSelectedIndex());
      q.setSelectedIndex(tq.getSelectedIndex());
      P.setSelectedIndex(tP.getSelectedIndex());
      Q.setSelectedIndex(tQ.getSelectedIndex());        
     }
   }

   /*------------------------------------------------------------------------------------
       Save the series to files in a uSimData directory 



   ------------------------------------------------------------------------------------*/
   public void saveSeries(double[] series, String file) 
   {
       int i; int le = series.length;
       if(le != nObs)
       {System.out.println("You might have length problems, just sayin'");}

       try{  
           PrintWriter out = new PrintWriter(new FileWriter(file));
           for(i=0; i < le; i++) {out.println(series[i]);}
 
           out.close(); System.out.println("Series successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();}
   }
   public void saveForecast(double[] series, double[] forecast, String file)
   {
       int i; int le = series.length; int lef = forecast.length/3; 
       int p = 1; // select which forecast upper = 0, mid = 1, low = 2   

       if(le != nObs)
       {System.out.println("You might have length problems, just sayin'");}

       try{  
           PrintWriter out = new PrintWriter(new FileWriter(file));
           for(i=0; i < le; i++) {out.println(series[i]);}
           for(i=0; i < lef; i++) {out.println(forecast[p*24 + i]);}
           out.close(); System.out.println("Series successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();}

   }

   public void saveSeasAdj(double[] series, double[] seas, String file)
   {
       int i; int le = series.length;
       if(le != nObs)
       {System.out.println("You might have length problems, just sayin'");}

       try{  
           PrintWriter out = new PrintWriter(new FileWriter(file));
           for(i=0; i < le; i++) {out.println(series[i] - seas[i]);}
           out.close(); System.out.println("Series successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();}
   }

   public void saveSignalCoeffs(double[] poly, String file)
   {
       int i; int le = poly.length;

       try{  
           PrintWriter out = new PrintWriter(new FileWriter(file));
           for(i=0; i < le; i++) {out.println(poly[i]);}
 
           out.close(); System.out.println("Data successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();} 
   }
   /*------------------------------------------------------------------------------------
   ------------------------------------------------------------------------------------*/





   public void setUpMenu(JFrame frame)
   {

      ActionListener aListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
           FileFilter type1 = new ExtensionFilter("Data files", new String[] {".dat"});
           //filter=new FileNameExtensionFilter("series files","dat","usim");
           fc.addChoosableFileFilter(type1);          
          
           int returnVal = fc.showOpenDialog(IMetricaProgram.this);

           if (returnVal == JFileChooser.APPROVE_OPTION) 
           {
                File file = fc.getSelectedFile();
                System.out.println("Opening: " + file.getName() + "." + "\n");
                readData(file);
           }
           else 
           {System.out.println("Open command cancelled by user." + "\n");}
        
      }
      };

      ActionListener getDataActionListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
          if(event.getSource() == getComps)
          {setRegCmpntStuff();}
          else if(event.getSource() == getData)
          {
           FileFilter type1 = new ExtensionFilter("Data files", new String[] {".dat"});
           //filter=new FileNameExtensionFilter("series files","dat","usim");
           fc.addChoosableFileFilter(type1);          
          
           int returnVal = fc.showOpenDialog(IMetricaProgram.this);

           if (returnVal == JFileChooser.APPROVE_OPTION) 
           {
                File file = fc.getSelectedFile();
                System.out.println("Opening: " + file.getName() + "." + "\n");
                reg.readData(file); 
                reg.reInitializePanel();
           }
           else 
           {System.out.println("Open command cancelled by user." + "\n");}   
          }
          else if(event.getSource() == bayes_save_single)
          {
             BayesCronosPan.print_to_file();
          }

          else if(event.getSource() == bayes_update)
          {updateBayesianPlot();}

          else if(event.getSource() == bayes_export_data) //save plotted data to data control
          {simulate.getExportedData(BayesCronosPan.getPlottedData());}
          
      }
      };

      ActionListener anListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
           FileFilter type1 = new ExtensionFilter("Metafiles",
            new String[] { ".dta"});


           //filter=new FileNameExtensionFilter("metafiles","dta");
           fc.addChoosableFileFilter(type1);
           int returnVal = fc.showOpenDialog(IMetricaProgram.this);

           if (returnVal == JFileChooser.APPROVE_OPTION) 
           {
                File file = fc.getSelectedFile();
                System.out.println("Opening: " + file.getName() + "." + "\n");
                readMetaData(file);
           }
           else 
           {System.out.println("Open command cancelled by user." + "\n");}
        
      }
      };


      ActionListener polyListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
         if(event.getSource() == trendPoly) 
         {saveSignalCoeffs(smc.trendpoly, "iMetric-trend.sig");}
         else if(event.getSource() == tiPoly) 
         {saveSignalCoeffs(smc.tipoly, "iMetric-trendirr.sig");}        
         else if(event.getSource() == seasPoly) 
         {saveSignalCoeffs(smc.seaspoly, "iMetric-seas.sig");}  
         else if(event.getSource() == trnsPoly) 
         {saveSignalCoeffs(smc.trnspoly, "iMetric-trns.sig"); saveSignalCoeffs(smc.trnsdenom, "iMetric-trnsdenom.sig"); }       
      }
      };

      ActionListener compActionListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {

        for(int i=0; i < n_regCmpnts; i++)
        {
          if(event.getSource() == compItems[i])
          {
             if(toFilex == true)
             {
               saveSeries(reg.GetComponent(i,regFore), "iMetric-regSeries" + Integer.toString(i+1) + ".dat");
               reg.changeHighlight(-1);
             }
          }
        }
   
      }
      };


      ActionListener dfaActionListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
        int i; File file;
        for(i=0; i < n_hist; i++)
        {
          if(event.getSource() == dfaItems[i])
          {mdfa.mdfa_canvas.saveFromQueue(i); mdfa.mdfa_canvas.changeHighlight(-1);}
        }
        if(event.getSource() == exportMDFA)
        {
            //System.out.println("target n_obs = " + simulate.target_series.length);
            mdfa.inputSimulatedData(simulate.getTarget_series(), simulate.getSim_data(), simulate.getRepresent(), simulate.n_sym);
            //enableX13Simulator(false); sim_check.setSelected(false); 
            sim_mdfa_check.setSelected(false);
            
            if(tradingInterface.isSelected() && simulate.logprice_data != null)
            {
              mdfa.setPriceData(simulate.logprice_data,simulate.n_price);
              tradingMode.setEnabled(true);
              tradingMenu.setEnabled(true);
              
              if(googHigherFreq.isSelected() && simulate.hf_data != null)
              {mdfa.setHFPrice(simulate.hf_data, 5); tradeMultiFreq.setEnabled(true);} 
            }
          
            mdfa.l1Bar.setValue(36);
        }
        else if(event.getSource() == exportREG)
        {
            inputFromSimulator(simulate.getTarget_series()); enableX13Simulator(false);
        }
        else if(event.getSource() == exportSARIMA)
        {
            inputFromSimulator(simulate.getTarget_series()); enableX13Simulator(false);
        }
        else if(event.getSource() == exportBayes)
        {
          BayesCronosPan.inputMultipleData(simulate.getSim_data(), simulate.getRepresent(), simulate.getNobs());
                
        }
        else if(event.getSource() == saveTarget) 
        {
           file = new File("iMetrica-simtarget.dat");            
           try //before u buy
           {  
             PrintWriter out = new PrintWriter(new FileWriter(file));
             for(i=0; i < simulate.getTarget_series().length; i++) {out.println(simulate.getTarget_series()[i]);}
             out.close(); System.out.println("Series successfully saved in " + file);
           } 
           catch (IOException e) {e.printStackTrace();}
        }
        else if(event.getSource() == saveData)
        {

           file = new File("iMetrica-simdata-.dat");            
           try //before u buy
           {  
             PrintWriter out = new PrintWriter(new FileWriter(file));
             for(i=0; i < simulate.getTarget_series().length; i++) {out.println(simulate.getTarget_series()[i]);}
          
             for(int k=0;k<simulate.getN_sim_series();k++)
             {
               if(simulate.getRepresent()[k])
               { 
                  out.println("series " + k);
                  double[] temp = simulate.getSim_data().get(k);                 
                  for(i=0; i < temp.length; i++) {out.println(temp[i]);}
               }               
             }
             out.close(); System.out.println("Series successfully saved in " + file);
           } 
           catch (IOException e) {e.printStackTrace();}
        }
 
      }
      };

      ItemListener seriesItemListener = new ItemListener() {
       public void itemStateChanged(ItemEvent e)
       {
         int i; boolean sel;  Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;} 
         else {sel = true;}

         for(i=0;i<simulate.getN_sim_series();i++)
         {
            if(source == seriesItems[i])
            {
              simulate.addToMix(i,sel); simulate.getTheatre().changeHighlight(-1);
            }
         }
         if(source == useSarima)
         {mdfa.setSarimaCheck(sel);}  
         if(source == autoComp)
         {
           estimate = sel; 
           if(estimate) {x13kardan();} 
           smc.setCompute(estimate); 
           smc.go();
           activateSlidingSpan(sel);
         }
       }
      };



      ActionListener seriesListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {
         
         if(event.getSource() == seriesItem) 
         {saveSeries(smc.t_series,"iMetric-tseries.dat");}
         else if(event.getSource() == seriesFItem)
         {saveForecast(smc.t_series,smc.forecasts,"iMetric-forecasted.dat");}
         else if(event.getSource() == seasItem)
         {saveSeries(smc.seassig,"iMetric-seasonal.dat");}
         else if(event.getSource() == trendItem)
         {saveSeries(smc.trendsig,"iMetric-trend.dat");}
         else if(event.getSource() == saItem)
         {saveSeasAdj(smc.t_series,smc.seassig,"iMetric-seasAdj.dat");}
         else if(event.getSource() == cycleItem)
         {saveSeries(smc.cyclesig,"iMetric-regAdj.dat");}
         
          else if(event.getSource() == track_series)
          {smc.setPlotTracker(0);}
          else if(event.getSource() == track_cycle)
          {smc.setPlotTracker(1);}
          else if(event.getSource() == track_trend)
          {smc.setPlotTracker(2);}
          else if(event.getSource() == track_seas)
          {smc.setPlotTracker(3);}
          else if(event.getSource() == track_SA)
          {smc.setPlotTracker(4);}     
          else if(event.getSource() == track_noneUsim)
          {smc.setPlotTracker(-1);} 
         
      }
      };

      MouseListener menuSlistener = new MouseListener() {
      public void mouseEntered(MouseEvent event) {
        
         boolean enable = false;
         if(event.getSource() == seriesItem) 
         {smc.changeHighlight(0); smc.go();}
         else if(event.getSource() == seriesFItem)
         {smc.changeHighlight(1); smc.go();}
         else if(event.getSource() == seasItem)
         {smc.changeHighlight(2); smc.go();}
         else if(event.getSource() == trendItem)
         {smc.changeHighlight(3); smc.go();}
         else if(event.getSource() == saItem)
         {smc.changeHighlight(4); smc.go();}
         else if(event.getSource() == cycleItem)
         {smc.changeHighlight(5); smc.go();}   
         else if(event.getSource() == simMenu)
         {
                 
           for(int i=0;i<simulate.getN_sim_series();i++) {seriesItems[i].setEnabled(simulate.mix[i]);}
           for(int i=simulate.getN_sim_series(); i<10;i++) {seriesItems[i].setEnabled(false);}

           if(simulate.getN_sim_series() > 0)  
           {enable = true;}  

           saveData.setEnabled(enable); exportSARIMA.setEnabled(enable);
           exportMDFA.setEnabled(enable); saveTarget.setEnabled(enable);
           exportREG.setEnabled(enable);  exportBayes.setEnabled(enable);
                   
         }
         for(int i=0;i<n_regCmpnts;i++)
         {
           if(event.getSource() == compItems[i]) 
           {reg.changeHighlight(i);}
         }

         for(int i=0;i<n_hist;i++)
         {
           if(event.getSource() == dfaItems[i]) 
           {mdfa.mdfa_canvas.changeHighlight(i);}
         }  

         if(simPanelx)
         {
          for(int i=0;i<simulate.getN_sim_series();i++)
          {
           if(event.getSource() == seriesItems[i])
           {simulate.getTheatre().changeHighlight(i);}
          }
         }
 
      }

      public void mouseExited(MouseEvent event) {
         smc.changeHighlight(-1); smc.go();  
         reg.changeHighlight(-1); //reg.go();
         mdfa.mdfa_canvas.changeHighlight(-1); 
         if(simPanelx) {simulate.getTheatre().changeHighlight(-1);}
         
      }
      public void mouseReleased(MouseEvent e) {
      }
      public void mouseClicked(MouseEvent e) {
      }
      public void mousePressed(MouseEvent e) {
      }
      };

      ItemListener fileListener = new ItemListener() {
        public void itemStateChanged(ItemEvent e)
        {
         boolean sel;
         Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
 
         if(source == toFile) {toFilex = sel;}
       } 
      };

      ActionListener exitListener = new ActionListener() {
      public void actionPerformed(ActionEvent event) 
      {System.exit(0);}};

 
      //------------- Action Listener for the mdfa stuff -----------

      ItemListener bayesItemListener = new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          int i;
          JRadioButtonMenuItem source = (JRadioButtonMenuItem)e.getItemSelectable();
          boolean sel; 
          if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
          else{sel = true;}
 
          for(i=0; i<14; i++)
          {
             if(source == plot_bseries[i] && sel)
             {BayesCronosPan.highlight_plot(i);} 
          }
        }      
      };
      
      
      ActionListener mdfaActionListener = new ActionListener() 
      {
        File file; int val;
        public void actionPerformed(ActionEvent event)
        {
          if(event.getSource() == saveFilterDes)
          {saveFilterFile();}//    mdfa.saveParameters(true);}
          else if(event.getSource() == saveFilterEncryp)
          {saveFilterFileEncrypt();}
          else if(event.getSource() == saveFilterCoeffsMenu)
          {saveFilterCoeffs();}
          else if(event.getSource() == loadPriceFilterCoeffs)
          {
             val = fc.showOpenDialog(IMetricaProgram.this);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               mdfa.loadPriceFilterCoeffs(file);
               tradeDuplex.setEnabled(true); tradeDuplex2.setEnabled(true);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }             
          else if(event.getSource() == loadFilterCoeffsMenu)
          {
             val = fc.showOpenDialog(IMetricaProgram.this);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               mdfa.loadFilterCoeffs(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }          
          else if(event.getSource() == savePlotFile)
          {mdfa.saveParameters(false); mdfa.saveFilteredData(0); }
          else if(event.getSource() == savePlot)
          {
             mdfa.mdfa_canvas.saveToHistorial(); 
             mdfa.parameterSnapshot(); 
             dfaItems[n_hist].setEnabled(true); n_hist++;
             mdfa.saveFilter();
          }
          else if(event.getSource() == clearPlots)
          {
             mdfa.mdfa_canvas.clearHistorical(); 
             for(int i=0;i<n_hist;i++) 
             {dfaItems[i].setEnabled(false);} n_hist=0; 
             mdfa.clearFilters();
          } 
          else if(event.getSource() == uploadData)
          {
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readData(file);  
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally 
              mdfa.setReturnStrategyOn();
            } 
            else {System.out.println("Open command cancelled by user.");}          
          }
          else if(event.getSource() == uploadReturnData)
          {
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readReturnData(file);  
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally          
            } 
            else {System.out.println("Open command cancelled by user.");}          
          } 
          else if(event.getSource() == uploadPortfolioData)
          {
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readPortfolioData(file);  
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally          
            } 
            else {System.out.println("Open command cancelled by user.");}          
          }        
          else if(event.getSource() == uploadTBillData)
          {
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readTBillData(file);  
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally          
            } 
            else {System.out.println("Open command cancelled by user.");}          
          }             
          else if(event.getSource() == uploadFilter)
          {
             val = fc.showOpenDialog(IMetricaProgram.this);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               mdfa.uploadParameters(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }
          else if(event.getSource() == csvData)
          {
            FileFilter type1 = new ExtensionFilter("CSV Data files", new String[] {".csv"});
            fc.addChoosableFileFilter(type1);          
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readCSVData(file);
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally
            } 
            else {System.out.println("Open command cancelled by user.");}          
          }
          else if(event.getSource() == marketData)
          {
            FileFilter type1 = new ExtensionFilter("CSV Market Data files", new String[] {".csv"});
            fc.addChoosableFileFilter(type1);          
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              simulate.readCSVDataMarketASSET(file);
              nObsBar.setValue(simulate.getNobs());  //-- change number of Nobs globally
            } 
            else {System.out.println("Open command cancelled by user.");}          
          }
          else if(event.getSource() == marketMenu)
          {
            marketDataPanel.setModal(false);
            marketDataPanel.setVisible(true);
          }
          else if(event.getSource() == googIntradayMenuItem)
          {
            googIntradayDialog.setModal(false);
            googIntradayDialog.setVisible(true);
          }
          else if(event.getSource() == IF5mItem)
          {
            FileFilter type1 = new ExtensionFilter("IQFeed Data", new String[] {".dat"});
            fc.addChoosableFileFilter(type1);          
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              file = fc.getSelectedFile();
              
          
            //simulate.getHF_CSV_Data(20, 390, 50, false, true, false);
            //simulate.getFutures_CSV_Data(781, 1560, true, true);
             //simulate.getSP_CSV_Data(4, 36, 78, false, true, true);
             //simulate.getTWS_CSV_Data();
              //simulate.getTWS_CSV_Data(file);
              
              simulate.readIntradayBarData(file, "2013-02-01", "2015-05-11", "16:00:00");
              
           }  
             //simulate.getHK_CSV_Data();
          }
          else if(event.getSource() == quandlPanelMenu)
          {
            quandlDialog.setModal(false);
            quandlDialog.setVisible(true);
          }
          else if(event.getSource() == hfreqDialogPanelMenu)
          {
            hfreqDialogPanel.setModal(false); 
            hfreqDialogPanel.setVisible(true);
          }
          else if(event.getSource() == tradeoptimDialog)
          {
            tradeoptimDialogPanel.setModal(false);
            tradeoptimDialogPanel.setVisible(true);
          }
          else if(event.getSource() == tradestatDialog)
          {
            tradestatDialogPanel.setModal(false);
            tradestatDialogPanel.setVisible(true);
          }
          else if(event.getSource() == tradeoutstatDialog)
          {
            tradesoutstatDialogPanel.setModal(false);
            tradesoutstatDialogPanel.setVisible(true);
          }        
          else if(event.getSource() == tradeparamDialog)
          {
            tradeparamDialogPanel.setModal(false);
            tradeparamDialogPanel.setVisible(true);
          }
          else if(event.getSource() == outSampstatDialog)
          {
            outSampstatDialogPanel.setModal(false);
            outSampstatDialogPanel.setVisible(true);
          }          
          else if(event.getSource() == trueOutSamplePanelMenu)
          {
             trueOutSampleDialogPanel.setModal(false);
             trueOutSampleDialogPanel.setVisible(true);        
          }
          else if(event.getSource() == dataMixDialog)
          {
            dataMixDialogPanel.setModal(false);
            dataMixDialogPanel.setVisible(true);
          }
          else if(event.getSource() == simulMenu)
          {
            simulatePanel.setModal(false);
            simulatePanel.setVisible(true);
          }
          else if(event.getSource() == diagnosticMenu)
          {
            dfaDiagnosticPanel.setModal(false);
            dfaDiagnosticPanel.setVisible(true);
          }
          else if(event.getSource() == slideSpanMenu)
          {
            slideSpanDialog.setModal(false);
            slideSpanDialog.setVisible(true);
          }
          else if(event.getSource() == sweepSpanMenu)
          {
            sweepSpanDialog.setModal(false);
            sweepSpanDialog.setVisible(true);
          }
          else if(event.getSource() == adaptiveUpdateMenu)
          {
             adaptiveUpdateDialogPanel.setModal(false);
             adaptiveUpdateDialogPanel.setVisible(true);
          }
          else if(event.getSource() == envisionDialog)
          {
             envisionDialogPanel.setModal(false);
             envisionDialogPanel.setVisible(true);      
             if(computeErrorCheck.isSelected()) mdfa.envision = true;
          }
          else if(event.getSource() == addFilteredData)
          {
             //add filtered data to M-DFA module
             int nobs,nreps,shift,j,i;

             nobs = mdfa.n_obs; nreps = mdfa.n_rep;
             double[] temp; double[] targ;

             //get filtere data------
             int xf_length = mdfa.mdfa.xf.length;
 
             ArrayList<double[]> listdata = new ArrayList<double[]>();

             temp = new double[xf_length];
             for(i=0;i<xf_length;i++)
             {
               temp[i] = mdfa.mdfa.xf[i];
             }
                
             shift = nobs - xf_length;
             targ = new double[xf_length];

            
             for(i=shift;i<nobs;i++)
             {
                 targ[i-shift] = mdfa.mdfa.tseries[i];
             }
   
             //----------- add target to explanatory if only one series --------
             if(nreps == 1)
             {
               double[] series = new double[xf_length];
               for(i=shift;i<nobs;i++)
               {series[i-shift] = mdfa.mdfa.tseries[i];}
               listdata.add(series);
             }

             //---- add the old explanatory series--------------
             for(j=1;j<nreps;j++)
             {
               double[] series = new double[xf_length];
               for(i=shift;i<nobs;i++)
               {
                 series[i-shift] = mdfa.mdfa.tseries[nobs*j+i];
               }
               listdata.add(series);
             }
   
             //---- add the filtered series--------------
             listdata.add(temp);

             //---- now update number of reps -----------            
             nreps = listdata.size();

             //---- now add series ----------------------
             boolean[] reps = new boolean[mdfa.series_max];
             for(j=0;j<mdfa.series_max;j++) {reps[j] = false;} 
             for(j=0;j<nreps;j++) {reps[j] = true;}

             mdfa.inputSimulatedData(targ, listdata, reps, 0);
             sim_mdfa_check.setSelected(false);

             mdfa.l1Bar.setValue(36);
             mdfa.l1Bar.setValue(0);
 
             //---- now if price data, change price data
             if(mdfa.priceDataSet)
             {  
               int nn = targ.length;
               double[] prc = new double[nn*simulate.n_price];
 
               for(i=0;i<simulate.n_price;i++)
               {
                for(j=0;j<nn;j++)
                { 
                  prc[nn*i + j] = simulate.logprice_data[simulate.getN_obs()*i + (mdfa.L-1) + j];
                }
               }
               mdfa.setPriceData(prc,simulate.n_price);
             }            
          }
          else if(event.getSource() == loadH0Filter)
          {
            FileFilter type1 = new ExtensionFilter("IQFeed Data", new String[] {".dat"});
            fc.addChoosableFileFilter(type1);          
            val = fc.showOpenDialog(IMetricaProgram.this);
            if (val == JFileChooser.APPROVE_OPTION) 
            {
              mdfa.setH0Filter(fc.getSelectedFile());
            }            
          }
          else if(event.getSource() == track_account)
          {mdfa.getAccount_canvas().setPlotTracker(0);}
          else if(event.getSource() == track_price)
          {mdfa.getAccount_canvas().setPlotTracker(1);}
          else if(event.getSource() == track_signal)
          {mdfa.getAccount_canvas().setPlotTracker(2);}
          else if(event.getSource() == track_logreturns)
          {mdfa.getAccount_canvas().setPlotTracker(3);}
          else if(event.getSource() == track_none)
          {mdfa.getAccount_canvas().setPlotTracker(-1);}
          else if(event.getSource() == tradeLogDiff)
          {mdfa.setTradingFunction(0);}
          else if(event.getSource() == tradeLogPrice)
          {mdfa.setTradingFunction(1);}
          else if(event.getSource() == tradeDuplex)
          {mdfa.setTradingFunction(2);}
          else if(event.getSource() == tradeDuplex2)
          {mdfa.setTradingFunction(3);}
          else if(event.getSource() == tradeMultiFreq)
          {mdfa.setTradingFunction(4);}          
          
        }
      };


     //------------------  Create menu bar -----------------------

     int i;
     menuBar = new JMenuBar();
     menuBar.setBorder(new BevelBorder(BevelBorder.RAISED));
     menuBar.setBorderPainted(true);
     //container.add(menuBar, BorderLayout.NORTH);


     //-------Data Control and Simulation Menu------------------------------------------------    
     //
     //---------------------------------------------------------------------------------------

     simMenu = new JMenu("Data Input/Export",true); simMenu.addActionListener(dfaActionListener);
     
     tradingInterface = new JCheckBoxMenuItem("Financial Trading On/Off");
     tradingInterface.addItemListener(new ItemListener()  {
        public void itemStateChanged(ItemEvent e)
        {
         boolean sel;
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         logReturns.setSelected(sel);
         
        }
       }
     );
 

     simMenu.add(tradingInterface); 

     panelMenu = new JMenu("Extra Panels",true); 
     simulMenu = new JMenuItem("Observations Panel"); simulMenu.addActionListener(mdfaActionListener);
     //panelMenu.add(simulMenu);

     diagnosticMenu = new JMenuItem("DFA Diagnostics Panel"); diagnosticMenu.addActionListener(mdfaActionListener);
     panelMenu.add(diagnosticMenu);     
     diagnosticMenu.setEnabled(false);

     csvData =  new JMenuItem("Load CSV Data"); csvData.setToolTipText("Upload multivariate data from CSV file");    
     uploadData = new JMenuItem("Load Data"); uploadData.setToolTipText("Upload multivariate data from file");
     uploadReturnData = new JMenuItem("Load Return Data"); uploadReturnData.setToolTipText("Upload univariate return data from file");
     uploadPortfolioData = new JMenuItem("Load Portfolio File"); uploadPortfolioData.setToolTipText("Upload meta files from file");
     uploadTBillData = new JMenuItem("Load TBill File"); uploadTBillData.setToolTipText("Upload TBill file from file");
     marketData = new JMenuItem("Load Market CSV Data"); marketData.setToolTipText("Upload multivariate market data from CSV file");      
     saveData = new JMenuItem("Save Data"); saveData.setToolTipText("Save data to file");
     marketMenu = new JMenuItem("Load Market Data"); marketMenu.setToolTipText("Upload Market Data Menu");         
     hfreqDialogPanelMenu = new JMenuItem("Load Realized Volatilty");
     hfreqDialogPanelMenu.setToolTipText("Upload Log-Returns and Estimated Realized Volatilty");     
     
     googIntradayMenuItem = new JMenuItem("Load Google Intraday Data");
     IF5mItem = new JMenuItem("Load IF5m Data");
      
     quandlPanelMenu = new JMenuItem("Download Data from QuandL");
     
      
     exportMDFA = new JMenuItem("Export Data to IMDFA module");
     exportREG = new JMenuItem("Export Data to State Space module");
     exportSARIMA = new JMenuItem("Export Data to X13 module");
     exportBayes =  new JMenuItem("Export Data to BayesCronos module");
     saveTarget = new JMenuItem("Construct Target"); 

     simMenu.add(csvData); csvData.addActionListener(mdfaActionListener);
     simMenu.add(marketData); marketData.addActionListener(mdfaActionListener);
     simMenu.add(marketMenu); marketMenu.addActionListener(mdfaActionListener);
     simMenu.add(googIntradayMenuItem); googIntradayMenuItem.addActionListener(mdfaActionListener);
     simMenu.add(quandlPanelMenu); quandlPanelMenu.addActionListener(mdfaActionListener);
     simMenu.add(IF5mItem); IF5mItem.addActionListener(mdfaActionListener);
     simMenu.add(hfreqDialogPanelMenu); hfreqDialogPanelMenu.addActionListener(mdfaActionListener);
     simMenu.add(uploadData); uploadData.addActionListener(mdfaActionListener);
     simMenu.add(uploadReturnData); uploadReturnData.addActionListener(mdfaActionListener);
     simMenu.add(uploadPortfolioData); uploadPortfolioData.addActionListener(mdfaActionListener);
     simMenu.add(uploadTBillData); uploadTBillData.addActionListener(mdfaActionListener);
     simMenu.addSeparator();

     simMenu.add(saveData);      saveData.setEnabled(false);     saveData.addActionListener(dfaActionListener);
     simMenu.add(exportMDFA);    exportMDFA.setEnabled(false);   exportMDFA.addActionListener(dfaActionListener);
     simMenu.add(exportREG);     exportREG.setEnabled(false);    exportREG.addActionListener(dfaActionListener);
     simMenu.add(exportSARIMA);  exportSARIMA.setEnabled(false); exportSARIMA.addActionListener(dfaActionListener);
     simMenu.add(exportBayes);   exportBayes.setEnabled(false);  exportBayes.addActionListener(dfaActionListener);
 
     simMenu.addSeparator();

     saveTarget.setEnabled(false); saveTarget.addActionListener(dfaActionListener);
     simMenu.add(saveTarget);

 

     seriesItems = new JCheckBoxMenuItem[12];        
     for(i=0;i<12;i++)
     {
       seriesItems[i] = new JCheckBoxMenuItem("Series " + Integer.toString(i+1),true);     
       simMenu.add(seriesItems[i]); 
       seriesItems[i].addItemListener(seriesItemListener); 
       seriesItems[i].addMouseListener(menuSlistener);  
       seriesItems[i].setEnabled(false); 
     }
     menuBar.add(simMenu);
     simMenu.setEnabled(true);
     simMenu.addMouseListener(menuSlistener);

     JMenuItem exit = new JMenuItem("Exit");
     simMenu.add(exit);
     exit.addActionListener(exitListener);

     
     
     //----------------------------- I_MDFA Menu -----------------------
 
     //------------------------------------------------------------------

     mdfaMenu = new JMenu("MDFA"); 
      

     useSarima = new JCheckBoxMenuItem("Use SARIMA data");
     useSarima.setToolTipText("Use the data from the SARIMA modeling module");
     useSarima.addItemListener(seriesItemListener);
  
     mdfaMenu.add(useSarima);
     mdfaMenu.addSeparator();
     
     saveFilterCoeffsMenu = new JMenuItem("Save Coefficients");      saveFilterCoeffsMenu.addActionListener(mdfaActionListener);
     loadFilterCoeffsMenu = new JMenuItem("Open Filter Coefficient File");  loadFilterCoeffsMenu.addActionListener(mdfaActionListener);
     loadPriceFilterCoeffs = new JMenuItem("Open Price Filter Coefficient File"); loadPriceFilterCoeffs.addActionListener(mdfaActionListener);
     saveFilterDes = new JMenuItem("Save Filter Parameters"); saveFilterDes.addActionListener(mdfaActionListener);
     saveFilterEncryp = new JMenuItem("Save Filter Parameters Encrytped"); saveFilterEncryp.addActionListener(mdfaActionListener);
     uploadFilter = new JMenuItem("Open Filter File");        uploadFilter.addActionListener(mdfaActionListener);
     savePlotFile = new JMenuItem("Save Filterd Data");       savePlotFile.addActionListener(mdfaActionListener);
     savePlot = new JMenuItem("Save Filter Plot");            savePlot.addActionListener(mdfaActionListener);

     adaptiveUpdateMenu = new JMenuItem("Adaptive Updating");
     adaptiveUpdateMenu.setToolTipText("Launch the adaptive filter updating control panel");
     adaptiveUpdateMenu.addActionListener(mdfaActionListener);
     
     addFilteredData = new JMenuItem("Add Filtered Data");  
     addFilteredData.setToolTipText("Add filtered data output X(t) to the explanatory set W_i(t), i=1..K (adaptive filter method)");  
     addFilteredData.addActionListener(mdfaActionListener);
   
     loadH0Filter = new JMenuItem("Load H0 Filter");  
     loadH0Filter.addActionListener(mdfaActionListener);
     
     turnOnH0 = new JCheckBoxMenuItem("Use H0 Filter");
     turnOnH0.addItemListener(new ItemListener() {
       public void itemStateChanged(ItemEvent e)
        {
         boolean sel; //computeFilter = true;
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         
         mdfa.toggleH0Filtering(sel);
         mdfa.mdfa.computeFilterGeneral(true,false);  
         mdfa.updatePlots(true,true);
        }
     });    
        
     dataMixDialog = new JMenuItem("Explanatory Data Mix Menu");
     dataMixDialog.setToolTipText("Open menu to choose explanatory data series");
     dataMixDialog.setEnabled(true);   
     dataMixDialog.addActionListener(mdfaActionListener);    

     clearPlots = new JMenuItem("Clear Saved Filters");       clearPlots.addActionListener(mdfaActionListener);

     mdfaMenu.add(saveFilterCoeffsMenu);
     mdfaMenu.add(loadFilterCoeffsMenu);
     mdfaMenu.add(loadPriceFilterCoeffs);
     mdfaMenu.add(uploadFilter);
     mdfaMenu.add(saveFilterDes);
     mdfaMenu.add(saveFilterEncryp);
     mdfaMenu.add(savePlotFile);
     mdfaMenu.addSeparator();
     mdfaMenu.add(savePlot); 
     mdfaMenu.add(adaptiveUpdateMenu);
     mdfaMenu.add(addFilteredData);
     mdfaMenu.add(loadH0Filter);
     mdfaMenu.add(turnOnH0);
     mdfaMenu.add(dataMixDialog);
     mdfaMenu.add(clearPlots);
     mdfaMenu.addSeparator();

     periodoWeightBox = new JCheckBoxMenuItem("HF Periodogram Weight On/Off");
     periodoWeightBox.setSelected(true);
     periodoWeightBox.setEnabled(false);
     periodoWeightBox.addItemListener(new ItemListener()  {
        public void itemStateChanged(ItemEvent e)
        {
         boolean sel;
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         mdfa.togglePeriodoWeight(sel);         
        }
       }
     );

     mdfaMenu.add(periodoWeightBox); 
 
     

     dfaItems = new JMenuItem[10];
     for(i=0; i < 10; i++)
     {
       dfaItems[i] = new JMenuItem("Filter " + Integer.toString(i+1));     
       mdfaMenu.add(dfaItems[i]); 
       dfaItems[i].addActionListener(dfaActionListener); 
       dfaItems[i].addMouseListener(menuSlistener);
       dfaItems[i].setEnabled(false);
     }
     
     envisionDialog = new JMenuItem("Show Cystal Ball");
     envisionDialog.setEnabled(false);
     envisionDialog.addActionListener(mdfaActionListener);
 
 
     computeErrorCheck = new JCheckBoxMenuItem("Compute Out-of-sample Error");
     computeErrorCheck.setSelected(false);
     computeErrorCheck.addItemListener(new ItemListener()
          {  
            public void itemStateChanged(ItemEvent e)
            {
             e.getItemSelectable();
             boolean sel; 
             if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
             else{sel = true;}
             mdfa.compute_error = sel; envisionDialog.setEnabled(sel);
            }
     });    
     
  
     mdfaMenu.addSeparator();
     sim_mdfa_check = new JCheckBoxMenuItem("IMDFA Simulator:");
     sim_mdfa_check.setSelected(true);
     sim_mdfa_check.addItemListener(this);
     mdfaMenu.add(sim_mdfa_check);    
     mdfaMenu.add(computeErrorCheck);
     mdfaMenu.add(envisionDialog);

     //------ TRADING MENU -------------------------------------------------------
     tradingMenu = new JMenu("Financial Trading",false);

     tradingMode = new JCheckBoxMenuItem("Turn On/Off Trading Platform");
     tradingMode.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e)
        {
         boolean sel;
         e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
         mdfa.setTradingInterface(sel);
         tradestatDialog.setEnabled(sel);
         tradeoutstatDialog.setEnabled(sel);
         tradeoptimDialog.setEnabled(sel);
         tradeparamDialog.setEnabled(sel);
         slideTrade.setEnabled(sel);       
         trueOutSamplePanelMenu.setEnabled(sel);
         
        }
       }
     );
    
     slideTrade = new JCheckBoxMenuItem("Sliding Window Analysis");
     slideTrade.setEnabled(false);
     slideTrade.setToolTipText("Enter sliding window analysis");
     slideTrade.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e)
        {
          if(e.getStateChange() == ItemEvent.DESELECTED){mdfa.deactivateSweepMode();outSampstatDialog.setEnabled(false);}
          else{mdfa.activateSweepMode(); outSampstatDialog.setEnabled(true);}
        }
       }
     );

     tradingMenu.setToolTipText("Turn the trading platform on in the MDFA module. Must have price/portfolio data loaded");
     uploadPrice = new JMenuItem("Upload Portfolio");
     
     tradeparamDialog = new JMenuItem("Trading Parameters");
     tradeparamDialog.setToolTipText("Open Trading Parameter Dialog Panel");
     tradeparamDialog.setEnabled(false);  
     tradeoptimDialog = new JMenuItem("Trading Optimization");
     tradeoptimDialog.setToolTipText("Open Trading Optimization Dialog Panel");
     tradeoptimDialog.setEnabled(false);
     tradestatDialog = new JMenuItem("Trading Statistics");
     tradestatDialog.setToolTipText("Open Trading Statistics Dialog Panel");
     tradestatDialog.setEnabled(false);


     tradeoutstatDialog = new JMenuItem("Out-of-sample Trading Statistics");
     tradeoutstatDialog.setToolTipText("Open Out-of-Sample Trading Statistics Panel");
     tradeoutstatDialog.setEnabled(false);     
     tradeoutstatDialog.addActionListener(mdfaActionListener);   

     trueOutSamplePanelMenu = new JMenuItem("True Out-of-sample points");
     trueOutSamplePanelMenu.setEnabled(false);     
     trueOutSamplePanelMenu.addActionListener(mdfaActionListener); 
               
     tradeparamDialog.addActionListener(mdfaActionListener); 
     tradeoptimDialog.addActionListener(mdfaActionListener); 
     tradestatDialog.addActionListener(mdfaActionListener);       
     tradingMenu.add(tradingMode);  
     tradingMenu.add(uploadPrice);
     tradingMenu.addSeparator();
     tradingMenu.add(tradeparamDialog);
     tradingMenu.add(tradeoptimDialog);
     tradingMenu.add(tradestatDialog);
     tradingMenu.add(tradeoutstatDialog);
     tradingMenu.add(trueOutSamplePanelMenu);
     tradingMode.setEnabled(false);

     tradingMenu.addSeparator();

     JMenu tradeFuncMenu = new JMenu("Financial Trading Function");
     tradeLogDiff = new JMenuItem("Trade Log-Returns");
     tradeLogPrice = new JMenuItem("Trade Log-Price");
     tradeDuplex = new JMenuItem("Trade Duplex"); tradeDuplex.setEnabled(false);
     tradeDuplex2 = new JMenuItem("Trade Duplex-Price"); tradeDuplex2.setEnabled(false);
     tradeMultiFreq = new JMenuItem("Trade MultiFrequency"); tradeMultiFreq.setEnabled(false);     
     tradeLogDiff.addActionListener(mdfaActionListener);
     tradeLogPrice.addActionListener(mdfaActionListener);  
     tradeDuplex.addActionListener(mdfaActionListener);
     tradeDuplex2.addActionListener(mdfaActionListener);
     tradeMultiFreq.addActionListener(mdfaActionListener);     
     tradeFuncMenu.add(tradeLogDiff);
     tradeFuncMenu.add(tradeLogPrice);
     tradeFuncMenu.add(tradeDuplex);
     tradeFuncMenu.add(tradeDuplex2);     
     tradeFuncMenu.add(tradeMultiFreq);
     tradingMenu.add(tradeFuncMenu);
     tradingMenu.addSeparator();

     JMenu data_tracker = new JMenu("Track Data");     
     track_account = new JMenuItem("Track Cursor-Account");
     track_price = new JMenuItem("Track Cursor-logPrice");
     track_logreturns = new JMenuItem("Track Cursor-logReturns");
     track_signal = new JMenuItem("Track Cursor-Signal");
     track_none = new JMenuItem("Track Cursor-None");
     
     track_none.addActionListener(mdfaActionListener);
     track_account.addActionListener(mdfaActionListener); 
     track_price.addActionListener(mdfaActionListener); 
     track_logreturns.addActionListener(mdfaActionListener); 
     track_signal.addActionListener(mdfaActionListener);
     data_tracker.add(track_none);
     data_tracker.add(track_account);
     data_tracker.add(track_price);
     data_tracker.add(track_logreturns);
     data_tracker.add(track_signal);
     tradingMenu.add(data_tracker);
     tradingMenu.setEnabled(false);
     tradingMenu.addSeparator();
     tradingMenu.add(slideTrade);
     
     outSampstatDialog = new JMenuItem("Out-of-Sample trading statistics interface");
     outSampstatDialog.setToolTipText("Open Out-of-Sample Trading Statistics Dialog Panel");
     outSampstatDialog.setEnabled(false);     
     outSampstatDialog.addActionListener(mdfaActionListener);   
     tradingMenu.add(outSampstatDialog);




     //-------uSimX13S Menu------------------------------------------------    
     //
     //---------------------------------------------------------------------------------------


     fileMenu = new JMenu("uSimX13-SEATS",false);
     //menuBar.add(fileMenu);

          

     slideSpanMenu = new JMenuItem("Sliding Span Control Panel");
     slideSpanMenu.setToolTipText("Open the sliding span control panel. uSimX13 computation must be turned on");
     slideSpanMenu.addActionListener(mdfaActionListener);
 
     sweepSpanMenu = new JMenuItem("Sweep Time Series Control Panel");
     sweepSpanMenu.setToolTipText("Open the sweep span control panel. uSimX13 computation must be turned on");
     sweepSpanMenu.addActionListener(mdfaActionListener);     


     autoComp = new JCheckBoxMenuItem("Activate uSimX13 Engine",false);
     autoComp.setToolTipText("Turn on automatic uSimX13 computation"); 
     autoComp.addItemListener(seriesItemListener);
     
     
     JMenuItem open = new JMenuItem("Open Data File"); 
     open.setToolTipText("Open file to import a single seasonal time series. Must be in column format.");

     open.addActionListener(aListener);  

     JMenuItem openMeta = new JMenuItem("Open MetaFile");
     openMeta.setToolTipText("Open .dta file to import multiple seasonal time series. Must be .dta file");
     openMeta.addActionListener(anListener);

  

     JMenu save_series = new JMenu("Save Series",true);
     JMenu save_polys = new JMenu("Save Signal Extraction Polynomials");

 
     seriesItem = new JMenuItem("Series"); seriesItem.addActionListener(seriesListener); 
     seriesItem.addMouseListener(menuSlistener); 
     seriesFItem = new JMenuItem("Series+Forecast"); seriesFItem.addActionListener(seriesListener);
     seriesFItem.addMouseListener(menuSlistener);
     seasItem = new JMenuItem("Seasonal");  seasItem.addActionListener(seriesListener);
     seasItem.addMouseListener(menuSlistener);
     trendItem = new JMenuItem("Trend"); trendItem.addActionListener(seriesListener);
     trendItem.addMouseListener(menuSlistener);
     saItem = new JMenuItem("Seasonally Adjusted Series"); saItem.addActionListener(seriesListener);
     saItem.addMouseListener(menuSlistener);
     cycleItem = new JMenuItem("Regression Adjusted Series"); cycleItem.addActionListener(seriesListener);
     cycleItem.addMouseListener(menuSlistener);

     save_series.add(seriesItem); 
     save_series.add(seriesFItem);
     save_series.add(seasItem);
     save_series.add(trendItem);
     save_series.add(saItem);
     save_series.add(cycleItem);

     trendPoly = new JMenuItem("Trend"); trendPoly.addActionListener(polyListener);
     seasPoly = new JMenuItem("Seasonal"); seasPoly.addActionListener(polyListener);
     tiPoly = new JMenuItem("Trend-Irregular"); tiPoly.addActionListener(polyListener);
     trnsPoly = new JMenuItem("Trans/Cyle"); trnsPoly.addActionListener(polyListener);

     
     JMenu track_series_menu = new JMenu("Track Data",true);
     track_series = new JMenuItem("Series");
     track_cycle = new JMenuItem("RegAdjusted"); 
     track_seas = new JMenuItem("Seasonal");
     track_trend = new JMenuItem("Trend");
     track_SA = new JMenuItem("Seasonally Adjusted");
     track_noneUsim = new JMenuItem("None");
     
     track_series.addActionListener(seriesListener);
     track_cycle.addActionListener(seriesListener);
     track_seas.addActionListener(seriesListener);
     track_trend.addActionListener(seriesListener);
     track_SA.addActionListener(seriesListener);
     track_noneUsim.addActionListener(seriesListener);
     
     track_series_menu.add(track_series);
     track_series_menu.add(track_cycle);
     track_series_menu.add(track_seas);
     track_series_menu.add(track_trend);
     track_series_menu.add(track_SA);
     track_series_menu.add(track_noneUsim);
     
 
     save_polys.add(trendPoly);
     save_polys.add(seasPoly);
     save_polys.add(tiPoly);
     save_polys.add(trnsPoly);
 

     sim_check = new JCheckBoxMenuItem("SARIMA Simulator Activate");
     sim_check.setSelected(true);
     sim_check.addItemListener(this);


     slideSpanCheckMenu = new JCheckBoxMenuItem("Sliding Span Activate");
     slideSpanCheckMenu.setSelected(false);
     slideSpanCheckMenu.addItemListener(new ItemListener() {
         public void itemStateChanged(ItemEvent e)
         {
            boolean sel; 
            e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
            else{sel = true;}

            nforeObsBar.setValue(0); 
            nbackObsBar.setValue(0);
            nforeObsText.setText("0"); 
            nbackObsText.setText("0");
            smc.nbackObs = 0; smc.nObsSpan = 0; smc.nforeObs = 0;

            smc.slidingSpan = sel;
            initializeUseless(nObs-smc.nbackObs-smc.nforeObs);
            smc.computeSARIMAmodel(estimate);            
            
            setDiagnostics(); 
            smc.go(); if(specCntrl){spec.go();}

            TrendModel.setEnabled(!sel); 
            SeasonalModel.setEnabled(!sel); 
            TrendIrregModel.setEnabled(!sel); 
            IrregModel.setEnabled(!sel); 
            sweepSpanMenu.setEnabled(sel);

         }
        }
       ); 

    
   
     simulMenu = new JMenuItem("Simulator Panel"); simulMenu.addActionListener(mdfaActionListener);

     
    //6985 5761
     fileMenu.add(autoComp);
     
     slideSpanCheckMenu.setEnabled(false);
     fileMenu.add(slideSpanCheckMenu);
     sweepSpanMenu.setEnabled(false);
     fileMenu.add(sweepSpanMenu);
     
     
     fileMenu.addSeparator();
     fileMenu.add(open);
     fileMenu.add(openMeta);

     fileMenu.addSeparator(); 

     fileMenu.add(track_series_menu);
     fileMenu.add(save_series);
     fileMenu.add(save_polys);  


     


     fileMenu.addSeparator();

     fileMenu.add(sim_check);  
     fileMenu.add(simulMenu);

   

     menuBar.add(fileMenu);
     


     //-------------------- State Space Menu -------------------------

     regMenu = new JMenu("State Space Menu");
     menuBar.add(regMenu); 

     getData = new JMenuItem("Open Data File");  
     getData.addActionListener(getDataActionListener);
     regMenu.add(getData);      
     regMenu.addSeparator();
     toFile = new JCheckBoxMenuItem("Output to File");
     toFile.setSelected(true); toFile.addItemListener(fileListener);
     regMenu.add(toFile);
     regMenu.addSeparator();

     getComps = new JMenuItem("Update Components");
     getComps.addActionListener(getDataActionListener);
     regMenu.add(getComps);

     compItems = new JMenuItem[10];
     for(i=0; i < 10; i++)
     {
       compItems[i] = new JMenuItem("Component " + Integer.toString(i+1));     
       regMenu.add(compItems[i]); 
       compItems[i].addActionListener(compActionListener); 
       compItems[i].addMouseListener(menuSlistener);
       compItems[i].setEnabled(false);
     }


     bayesMenu = new JMenu("BayesCronos Menu");
     menuBar.add(bayesMenu);
     
     plot_bseries = new JRadioButtonMenuItem[14];
 
     plot_bseries[0] = new JRadioButtonMenuItem("none");
     plot_bseries[0].setEnabled(true); plot_bseries[0].setSelected(true);
  
     for(i=0;i<5;i++) 
     {
       plot_bseries[i+1] = new JRadioButtonMenuItem("y_"+(i+1)+"(t)");
       plot_bseries[i+1].setEnabled(false); plot_bseries[i+1].setSelected(false); 
     }
     for(i=0;i<3;i++) 
     {
       plot_bseries[i+6] = new JRadioButtonMenuItem("h_"+(i+1)+"(t)");
       plot_bseries[i+6].setEnabled(false); plot_bseries[i+6].setSelected(false); 
     }   
     for(i=0;i<3;i++) 
     {
       plot_bseries[i+9] = new JRadioButtonMenuItem("f_"+(i+1)+"(t)");
       plot_bseries[i+9].setEnabled(false); plot_bseries[i+9].setSelected(false); 
     }
     plot_bseries[12] = new JRadioButtonMenuItem("r(t)");
     plot_bseries[12].setEnabled(false); plot_bseries[12].setSelected(false);

     plot_bseries[13] = new JRadioButtonMenuItem("e(t)");
     plot_bseries[13].setEnabled(false); plot_bseries[13].setSelected(false);  

     /*
     plot_series = new JRadioButtonMenuItem("Plot Series/Forecasts"); plot_series.setSelected(true);
     plot_residuals = new JRadioButtonMenuItem("Plot Residuals");
     plot_vol = new JCheckBoxMenuItem("Plot Volatilty"); plot_vol.setSelected(false);
     */

     bayes_save_single = new JMenuItem("Save Single Series");
     bayes_export_data = new JMenuItem("Export Plotted Series");
     bayes_update = new JMenuItem("Update Plots");

     bayes_save_single.setToolTipText("Save a single time series to a file. Select an available series in the list below and then click this menu item");
     bayes_update.setToolTipText("Update the radio buttons for the plots that are plotted on the plot canvas");
     bayes_export_data.setToolTipText("Export all the time series data on the plot canvas to the Data Control Module");

     bayes_save_single.addActionListener(getDataActionListener);
     bayes_export_data.addActionListener(getDataActionListener);
     bayes_update.addActionListener(getDataActionListener);

     bayesMenu.add(bayes_update);
     bayesMenu.add(bayes_save_single);
     bayesMenu.add(bayes_export_data);
     bayesMenu.addSeparator();

     ButtonGroup group = new ButtonGroup();
 
     for(i=0;i<14;i++)
     {
        group.add(plot_bseries[i]); 
        plot_bseries[i].addItemListener(bayesItemListener);
        bayesMenu.add(plot_bseries[i]);
     }
     
     
     //------------------------ Other Menu -------------------------------------
     // 
     //------------------------------------------------------------------------
     menuBar.add(mdfaMenu);
     menuBar.add(tradingMenu);
     menuBar.add(panelMenu);
     JMenu optionsMenu = new JMenu("Other");
     menuBar.add(optionsMenu);   
     optionsMenu.add(new JMenuItem("Go to Guide"));
     optionsMenu.add(new JMenuItem("About"));
     optionsMenu.addSeparator();     
   
     optionsMenu.add(new JMenuItem("Copy"));
     optionsMenu.add(new JMenuItem("Paste"));
     optionsMenu.addSeparator();


 
     frame.setJMenuBar(menuBar);

  }


  //----------------------- Method to set up EMD panels -----------


  //----------------------------------------------------------------

  public void updateBayesianPlot()
  {

   int i,nreps,n_factors;
   nreps = BayesCronosPan.n_rep;
   n_factors = BayesCronosPan.getN_factors();

   for(i=0;i<nreps;i++) 
   {
     if(BayesCronosPan.seriesPlot[i].isSelected())
     {plot_bseries[i+1].setEnabled(true);}
   }
   for(i=0;i<n_factors;i++) 
   {
     if(BayesCronosPan.alphasPlot[i].isSelected())
     {plot_bseries[6+i].setEnabled(true);}
   }  
   for(i=0;i<n_factors;i++) 
   {
     if(BayesCronosPan.factorPlot[i].isSelected())
     {plot_bseries[9+i].setEnabled(true);}
   }

   if(BayesCronosPan.residualPlot.isSelected())
   {plot_bseries[12].setEnabled(true);}

   if(BayesCronosPan.predictivePlot.isSelected())
   {plot_bseries[13].setEnabled(true);}   

  }
  


  public void setupEMDPanels(JPanel emdP)
  {
 
     int i; 
     //---- setup emd panel and border------------------------
     JPanel emdPanel = new JPanel();  
     TitledBorder emdBorder = new TitledBorder(new LineBorder(myBlue),"Select Empirical AM-FM Components");
     emdBorder.setTitleColor(Color.BLACK);
     TitledBorder emdBorder1 = new TitledBorder(new LineBorder(myBlue),"Time Series and EMD Components");
     emdBorder1.setTitleColor(Color.BLACK);
     TitledBorder emdBorder2 = new TitledBorder(new LineBorder(myBlue),"AM Components");
     emdBorder2.setTitleColor(Color.BLACK);
     TitledBorder emdBorder3 = new TitledBorder(new LineBorder(myBlue),"FM and Instantaneous Frequency Components"); 
     emdBorder3.setTitleColor(Color.BLACK);

     TitledBorder emdBorder4 = new TitledBorder(new LineBorder(myBlue),"Time-Frequency Plot"); 
     emdBorder4.setTitleColor(Color.BLACK);     

     //---- setup plot canvas-------------------------------------
     emd = new EMDamfm(800, 120, nObs);  
     emdAM = new EMDamCanvas(800, 100, nObs);
     emdFM = new EMDfmcanvas(800, 100, nObs, true);
     
     emdIMF = new EMDIMFcanvas(800, 120, nObs);
     emdIF = new EMDfmcanvas(800, 100, nObs, false);
     emdPhase = new EMDPhasecanvas(800, 100, nObs);

     emdPanel.setBorder(emdBorder);    
     //emd.setBorder(emdBorder1);
     //emdAM.setBorder(emdBorder2);
     //emdFM.setBorder(emdBorder3);

     //---- setup labels -------------------------------------
     String lbls1[] = new String[5];
     lbls1[0] = "Instrinstic Modes"; lbls1[1] = "AM Components"; 
     lbls1[2] = "FM Components"; lbls1[3] = "iFrequencies";
     lbls1[4] = "Phase Functions";
     JLabel lbls[] = new JLabel[5];
     for(i=0;i<5;i++) lbls[i] = new JLabel(lbls1[i]);
     //-----------------------------------------------------------
 
     mps2d = new JRadioButton[4];
     mps2d[0] = new JRadioButton("AM"); mps2d[0].setSelected(true); mps2d[0].setToolTipText("Amplitude Modulated Components");   
     mps2d[1] = new JRadioButton("FM"); mps2d[1].setToolTipText("Frequency Modulated Components");   
     mps2d[2] = new JRadioButton("IF"); mps2d[2].setToolTipText("Instantaneous Frequency derived from FM components");
     mps2d[3] = new JRadioButton("\u03D5"); mps2d[3].setToolTipText("Phase derived from FM components");
     ButtonGroup group = new ButtonGroup();
     group.add(mps2d[0]);group.add(mps2d[1]);group.add(mps2d[2]);group.add(mps2d[3]);
 
     // Put the radio buttons in a column in a panel
     jplRadio = new JPanel();
     jplRadio.setPreferredSize(new Dimension(50, 120));
     jplRadio.setLayout(new GridLayout(0, 1));
     jplRadio.add(mps2d[0]); jplRadio.add(mps2d[1]); jplRadio.add(mps2d[2]);
     jplRadio.add(mps2d[3]);

     JPanel plot2D = createEMDMap();
     plot2D.setBorder(emdBorder4);

     //---- create check boxes for control panel ----
     imfs = new JCheckBox[5];
     ams = new JCheckBox[5];
     fms = new JCheckBox[5];
     ifs = new JCheckBox[5]; 
     phase = new JCheckBox[5];
  
     for(i=0; i < 5; i++)
     {
       imfs[i] = new JCheckBox(""); imfs[i].setSelected(false);      
       ams[i] = new JCheckBox(""); ams[i].setSelected(false);
       fms[i] = new JCheckBox(""); fms[i].setSelected(false);     
       ifs[i] = new JCheckBox(""); ifs[i].setSelected(false);     
       phase[i] = new JCheckBox(""); phase[i].setSelected(false);   
     }
     //-----------------------------------------------------------

     //----- Setup recompute button -----------------------------
     ActionListener recompListener = new ActionListener() {
     public void actionPerformed(ActionEvent event) 
     {computeEMD(smc.t_series,nObs);}};
     reCompEMD.addActionListener(recompListener);



      //------ Create Item Listener ---------------------------
      ItemListener iListener = new ItemListener() {
        public void itemStateChanged(ItemEvent e)
        {
         int i; boolean sel;
         Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}

         for(i=0;i<5;i++)
         {
          if(source == imfs[i]) {plots_agg[i] = sel; updateAgg(i,sel); updateIMFs(i,sel);}
          else if(source == ams[i]) {plots_emd[i][0] = sel; updateAM(i,sel);}
          else if(source == fms[i]) {plots_emd[i][1] = sel; updateFM(i,sel);} 
          else if(source == ifs[i]) {plots_emd[i][2] = sel; updateIF(i,sel);}
          else if(source == phase[i]) {plots_emd[i][3] = sel; updatePhase(i,sel);}
         }
         if(source == mps2d[0]) {updateMap(0); selMap = 0;}  
         else if(source == mps2d[1]) {updateMap(1); selMap = 1;}
         else if(source == mps2d[2]) {updateMap(2); selMap = 2;}
         else if(source == mps2d[3]) {updateMap(3); selMap = 3;}
        }  
      };    
      //-----------------------------------------------------------
 
      //------------- Add ItemListener ---------------------------
      for(i=0; i < 5; i++)
      {
       imfs[i].addItemListener(iListener);
       ams[i].addItemListener(iListener); 
       fms[i].addItemListener(iListener);
       ifs[i].addItemListener(iListener);
       phase[i].addItemListener(iListener);
      }
       mps2d[0].addItemListener(iListener);
       mps2d[1].addItemListener(iListener);
       mps2d[2].addItemListener(iListener);
       mps2d[3].addItemListener(iListener);


     //------------- Create Layout ------------------------------
     GroupLayout emdLayout = new GroupLayout(emdPanel);
     emdLayout.setAutoCreateGaps(true);
     emdLayout.setAutoCreateContainerGaps(true);
     emdLayout.setHorizontalGroup(emdLayout.createSequentialGroup()
              .addGroup(emdLayout.createParallelGroup() 
                 .addComponent(lbls[0])
                 .addComponent(lbls[1])
                 .addComponent(lbls[2])
                 .addComponent(lbls[3])
                 .addComponent(lbls[4]))
              .addGroup(emdLayout.createParallelGroup()
                 .addComponent(imfs[0])
                 .addComponent(ams[0])
                 .addComponent(fms[0])
                 .addComponent(ifs[0])
                 .addComponent(phase[0]))
              .addGroup(emdLayout.createParallelGroup()
                 .addComponent(imfs[1])
                 .addComponent(ams[1])
                 .addComponent(fms[1])
                 .addComponent(ifs[1])
                 .addComponent(phase[1]))
              .addGroup(emdLayout.createParallelGroup()
                 .addComponent(imfs[2])
                 .addComponent(ams[2])
                 .addComponent(fms[2])
                 .addComponent(ifs[2])
                 .addComponent(phase[2]))
              .addGroup(emdLayout.createParallelGroup()
                 .addComponent(imfs[3])
                 .addComponent(ams[3])
                 .addComponent(fms[3])
                 .addComponent(ifs[3])
                 .addComponent(phase[3]))
              .addGroup(emdLayout.createParallelGroup()
                 .addComponent(imfs[4])
                 .addComponent(ams[4])
                 .addComponent(fms[4])
                 .addComponent(ifs[4])
                 .addComponent(phase[4])));


       emdLayout.setVerticalGroup(emdLayout.createSequentialGroup()
              .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(emdLayout.createSequentialGroup()
                  .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lbls[0])
                   .addComponent(imfs[0])
                   .addComponent(imfs[1])
                   .addComponent(imfs[2])
                   .addComponent(imfs[3])
                   .addComponent(imfs[4])) 
                 .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                   .addComponent(lbls[1])
                   .addComponent(ams[0])
                   .addComponent(ams[1])
                   .addComponent(ams[2])
                   .addComponent(ams[3])
                   .addComponent(ams[4]))
                 .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                   .addComponent(lbls[2])
                   .addComponent(fms[0])
                   .addComponent(fms[1])
                   .addComponent(fms[2])
                   .addComponent(fms[3])
                   .addComponent(fms[4]))
                 .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                   .addComponent(lbls[3])
                   .addComponent(ifs[0])
                   .addComponent(ifs[1])
                   .addComponent(ifs[2])
                   .addComponent(ifs[3])
                   .addComponent(ifs[4]))
                 .addGroup(emdLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                   .addComponent(lbls[4])
                   .addComponent(phase[0])
                   .addComponent(phase[1])
                   .addComponent(phase[2])
                   .addComponent(phase[3])
                   .addComponent(phase[4])))));         
       emdPanel.setLayout(emdLayout);    

       /*------------------------------------------------------------
          AMFMComp Layout
       -------------------------------------------------------------*/
       amfmComps = new JPanel(); phaseComps = new JPanel();
       GroupLayout amfmLayout = new GroupLayout(amfmComps);
       amfmLayout.setAutoCreateGaps(true);
       amfmLayout.setAutoCreateContainerGaps(true);
       amfmLayout.setHorizontalGroup(amfmLayout.createSequentialGroup()
          .addGroup(amfmLayout.createParallelGroup() 
           .addComponent(emd)
           .addComponent(emdAM)
           .addComponent(emdFM)));
       amfmLayout.setVerticalGroup(amfmLayout.createSequentialGroup()
          .addComponent(emd)
          .addComponent(emdAM)
          .addComponent(emdFM));
       amfmComps.setLayout(amfmLayout);

       /*------------------------------------------------------------
          AMFMComp Layout
       -------------------------------------------------------------*/
       GroupLayout phaseLayout = new GroupLayout(phaseComps);
       phaseLayout.setAutoCreateGaps(true);
       phaseLayout.setAutoCreateContainerGaps(true);
       phaseLayout.setHorizontalGroup(phaseLayout.createSequentialGroup()
          .addGroup(phaseLayout.createParallelGroup() 
           .addComponent(emdIMF)
           .addComponent(emdPhase)
           .addComponent(emdIF)));
       phaseLayout.setVerticalGroup(phaseLayout.createSequentialGroup()
          .addComponent(emdIMF)
          .addComponent(emdPhase)
          .addComponent(emdIF));
       phaseComps.setLayout(phaseLayout);
       //--------------------------------------------------------------

       emdPlotPane = new JTabbedPane(JTabbedPane.TOP);
       emdPlotPane.addTab("AM-FM Components", amfmComps);
       emdPlotPane.addTab("Phase-Instantaneous Freq.", phaseComps);
      

       GroupLayout epLayout = new GroupLayout(emdP);
       epLayout.setAutoCreateGaps(true);
       epLayout.setAutoCreateContainerGaps(true);
       epLayout.setHorizontalGroup(epLayout.createSequentialGroup()
          .addGroup(epLayout.createParallelGroup() 
           .addComponent(emdPlotPane)
              .addGroup(epLayout.createSequentialGroup()
                .addComponent(emdPanel)
                .addComponent(plot2D))));                
       epLayout.setVerticalGroup(epLayout.createSequentialGroup()
          .addComponent(emdPlotPane)
            .addGroup(epLayout.createParallelGroup()
              .addComponent(emdPanel)
              .addComponent(plot2D)));   
       emdP.setLayout(epLayout);
       //emdP.add(emdP);
                         
  }



  //----------Methods to update the Aggregate,AM,FM,Inst_f plots -----
  public void updateAgg(int i, boolean sel)
  {
     emd.setAggregateBool(i, sel);

     if(sel)
     {emd.addAggregate(i);}
     else
     {emd.subAggregate(i);}  
  }

  public void updateAM(int i, boolean sel)
  {emdAM.updateAM(i, sel);}

  public void updateFM(int i, boolean sel)
  {emdFM.updateFM(i,sel);}

  public void updateIF(int i, boolean sel)
  {emdIF.updateIF(i,sel);}
 
  public void updatePhase(int i, boolean sel)
  {emdPhase.updatePhase(i,sel);}  

  public void updateIMFs(int i, boolean sel)
  {emdIMF.updateIMFs(i,sel);} 
 
  public void updateMap(int i)
  {
     if(i==0) emdm.updateFM(emd.getAm(), nObs, emd.getN_imfs(), i);
     else if(i==1) emdm.updateFM(emd.getFm(), nObs, emd.getN_imfs(), i);
     else if(i==2) emdm.updateFM(emd.getInst_f(), nObs, emd.getN_imfs(), i);
     else if(i==3) emdm.updateFM(emd.getPhase(), nObs, emd.getN_imfs(), i);
     updateKey();
  }

  public void updateKey()
  {emdKey.updateColor(emdm.colorArray); hScale.updateData(emdm.dataRangeMax, emdm.dataRangeMin);}

  public void updateScale()
  {emdScale.updateN(nObs);}

  public void computeEMD(double[] series, int N)
  {
    if(emdPanelx) 
    {
      emd.AM_FM_decomp(series, N); 
      emdAM.setAM(emd.getAm(), N, emd.getN_imfs()); 
      emdFM.setFM(emd.getFm(), emd.getInst_f(), N, emd.getN_imfs());
      emdIF.setFM(emd.getFm(), emd.getInst_f(), N, emd.getN_imfs());
      emdPhase.setPhase(emd.getPhase(),N,emd.getN_imfs());   
      emdIMF.setIMFs(emd.getFm(), emd.getAm(), N, emd.getN_imfs());
      updateMap(selMap); freqScale.updateN_IMFS(emd.getN_imfs());
    }
  }
 

  /*-----------------------------------------------------------
    Create box for SARIMA model parameters- 
    -----------------------------------------------------------*/

  public void slidingSpanReset()
  {
    slidingSpan.setMaximum(nObsFixed/2);
    slidingSpan.setValue(0);
    slidingSpanText.setText("0");
  
  }
        
    
  public GroupLayout setUpParameterLayout(JPanel paramPanel)
  {
 
      
       text4 = new JTextField(3);
       text4.setText("" + phi1); 
       phi1bar = new JScrollBar(JScrollBar.HORIZONTAL,-10,1,-200,200);
       phi1bar.setValue(-10);      
       phi1bar.setUnitIncrement(1);
       phi1bar.addAdjustmentListener(this);
       
        //================ phi2 bar =======================       
       text5 = new JTextField(3);
       text5.setText("" + phi2);
       phi2bar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,-200,200);
       phi2bar.setValue(0);
       phi2bar.setUnitIncrement(1);
       phi2bar.addAdjustmentListener(this);

       //================ Phi bar =======================
       text6 = new JTextField(3);     
       text6.setText("" + Phi);
       Phibar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,-200,200);
       Phibar.setValue(0);
       Phibar.setUnitIncrement(1);     
       Phibar.addAdjustmentListener(this);

       //================ theta_1 bar =======================      
       text7 = new JTextField(3);     
       text7.setText("" + theta1);
       theta1bar = new JScrollBar(JScrollBar.HORIZONTAL,-60,1,-200,200);
       theta1bar.setValue(-60);
       theta1bar.setUnitIncrement(1);     
       theta1bar.addAdjustmentListener(this);

       //================ theta_2 bar =======================
       text8 = new JTextField(3);
       text8.setText("" + theta2);
       theta2bar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,-200,200);
       theta2bar.setValue(0);
       theta2bar.setUnitIncrement(1);
       theta2bar.addAdjustmentListener(this);

       //================ Theta bar ======================
       text9 = new JTextField(3);
       text9.setText("" + Theta);    
       Thetabar = new JScrollBar(JScrollBar.HORIZONTAL,-60,1,-200,200);
       Thetabar.setValue(-60);
       Thetabar.setUnitIncrement(1);
       Thetabar.addAdjustmentListener(this);

       JLabel innvarLabel = new JLabel("InnVar.: "); 
       innvarParameter = new JScrollBar(JScrollBar.HORIZONTAL,100,1,10,600);
       innvarParameter.setUnitIncrement(1); 
       text1 = new JTextField(3); 
       text1.setText("" + df.format(innvar));
       innvarParameter.addAdjustmentListener(this);

       JLabel tphi1,tphi2,ttheta1,ttheta2,tTheta,tPhi; 
       tphi1= new JLabel("AR_1");  tphi1.setToolTipText("Set value in simulation for nonseasonal \u03D5_1."); 
       tphi2= new JLabel("AR_2");  tphi2.setToolTipText("Set value in simulation for nonseasonal \u03D5_2."); 
       ttheta1= new JLabel("MA_1");   ttheta1.setToolTipText("Set value in simulation for nonseasonal \u03D1_1.");  
       ttheta2= new JLabel("MA_2");   ttheta2.setToolTipText("Set value in simulation for nonseasonal \u03D1_2."); 
       tTheta=new JLabel("S_MA");     tTheta.setToolTipText("Set value in simulation for seasonal \u0398."); 
       tPhi = new JLabel("S_AR");     tPhi.setToolTipText("Set value in simulation for seasonal \u03A6."); 
       //JLabel regLabel = new JLabel("Reg.:");       
       //regBar = new JScrollBar(JScrollBar.HORIZONTAL,500,10,0,1000);
       //regBar.setUnitIncrement(1);
       ////regBar.addAdjustmentListener(this)
  
       slidingSpanText = new JTextField(3);
       
  
       slidingSpan = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,nObs/2);
       slidingSpan.setUnitIncrement(1);
       slidingSpan.addAdjustmentListener(this);
       slidingSpanText.setText("0"); 
    
  
       series_scroll = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,50);
       series_scroll.setUnitIncrement(1);
       series_scroll.addAdjustmentListener(this);
       series_name.setText("0 files");
       JLabel files = new JLabel("Loaded files:");

       GroupLayout paramLayout = new GroupLayout(paramPanel);
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup() 
                .addComponent(tphi1) .addComponent(tphi2)) // .addComponent(phi1Label) .addComponent(phi2Label))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(phi1bar) .addComponent(phi2bar))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text4) .addComponent(text5))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(ttheta1) .addComponent(ttheta2)) //.addComponent(theta1Label) .addComponent(ttheta2Label))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(theta1bar) .addComponent(theta2bar))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text7) .addComponent(text8))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(tPhi) .addComponent(tTheta)) //.addComponent(PhiLabel) .addComponent(ThetaLabel))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(Phibar) .addComponent(Thetabar))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text6) .addComponent(text9))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(innvarLabel) .addComponent(files))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(innvarParameter) .addComponent(series_scroll))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text1) .addComponent(series_name)));


       paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                //.addGroup(paramLayout.createSequentialGroup()
                 // .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(tphi1) //.addComponent(phi1Label)
                   .addComponent(phi1bar)
                   .addComponent(text4)
                   .addComponent(ttheta1) //.addComponent(theta1Label)
                   .addComponent(theta1bar)
                   .addComponent(text7)
                   .addComponent(tPhi) //.addComponent(PhiLabel)
                   .addComponent(Phibar)
                   .addComponent(text6) 
                   .addComponent(innvarLabel)
                   .addComponent(innvarParameter)
                   .addComponent(text1)) 
                 //.addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)   
                .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)   
                  .addComponent(tphi2) //.addComponent(phi2Label)
                  .addComponent(phi2bar)
                  .addComponent(text5)
                  .addComponent(ttheta2) //.addComponent(ttheta2)
                  .addComponent(theta2bar)
                  .addComponent(text8)
                  .addComponent(tTheta) //.addComponent(tTheta)
                  .addComponent(Thetabar)
                  .addComponent(text9)
                  .addComponent(files)
                  .addComponent(series_scroll)
                  .addComponent(series_name)));

            return paramLayout;    
    }

    /*-------------------------------------------------------------
      Create Plottting Box
    ---------------------------------------------------------------*/

    public Box createPlottingBox()
    {
  
      //---------- Build the check boxes and group for plotting time domain
      seriesBox = new JCheckBox("Series"); seriesBox.setToolTipText("Plot the raw time series.");
      seriesBox.setSelected(true); seriesBox.setHorizontalTextPosition(JMenuItem.LEFT);
      foreBox = new JCheckBox("Forecasts"); foreBox.setToolTipText("Compute and plot the 24-step-ahead forecast.");
      foreBox.setSelected(false);  foreBox.setHorizontalTextPosition(JMenuItem.LEFT);
      seasBox = new JCheckBox("Seasonal"); seasBox.setToolTipText("Compute and plot the x13-ARIMA-SEATS seasonal-component.");
      seasBox.setSelected(false);   seasBox.setHorizontalTextPosition(JMenuItem.LEFT);
      trendBox = new JCheckBox("Trend"); trendBox.setToolTipText("Compute and plot the x13-ARIMA-SEATS trend-component.");
      trendBox.setSelected(false);  trendBox.setHorizontalTextPosition(JMenuItem.LEFT);
      saBox = new JCheckBox("Seas-Adjusted"); saBox.setToolTipText("Compute and plot the x13-ARIMA-SEATS seasonally-adjusted component.");
      saBox.setSelected(false);  saBox.setHorizontalTextPosition(JMenuItem.LEFT);
      cycleBox = new JCheckBox("Reg-Adjusted"); cycleBox.setToolTipText("Compute and plot the regression-adjusted component.");
      cycleBox.setSelected(false); cycleBox.setHorizontalTextPosition(JMenuItem.LEFT);

      seriesBox.addItemListener(this);
      foreBox.addItemListener(this);
      seasBox.addItemListener(this);
      saBox.addItemListener(this);
      cycleBox.addItemListener(this);   
      trendBox.addItemListener(this);


      //---------- Build the check boxes and group for plotting spectral
      plotInBox = new JCheckBox("I_n"); plotInBox.setToolTipText("Plot the periodogram of the (differenced) time series data."); 
      plotInBox.setSelected(false); plotInBox.setHorizontalTextPosition(JMenuItem.LEFT);
      plotFwBox = new JCheckBox("F_W,\u03C8"); plotFwBox.setToolTipText("Plot the spectral density of the estimated SARIMA model."); 
      plotFwBox.setSelected(false); plotFwBox.setHorizontalTextPosition(JMenuItem.LEFT);
      plotDnBox = new JCheckBox("D_n"); plotDnBox.setToolTipText("Plot the spectral density of differencing operator."); 
      plotDnBox.setSelected(false); plotDnBox.setHorizontalTextPosition(JMenuItem.LEFT);
      plotIGBox = new JCheckBox("I_n g_\u03C8"); plotIGBox.setToolTipText("Plot periodogram times signal-noise/F_W-squared ratio. Choose signal in panel to right."); plotIGBox.setHorizontalTextPosition(JMenuItem.LEFT);
      plotIGBox.setSelected(true);
      plotFGBox = new JCheckBox("F_W g_\u03C8"); plotFGBox.setToolTipText("Plot spectral density times signal-noise/F_W-squared ratio. Choose signal in panel to right."); plotFGBox.setHorizontalTextPosition(JMenuItem.LEFT);
      plotFGBox.setSelected(true);
      plotGBox = new JCheckBox("g_\u03C8"); plotGBox.setToolTipText("Plot signal-noise/F_W-squared ratio. Choose signal in panel to right.");
      plotGBox.setSelected(false); plotGBox.setHorizontalTextPosition(JMenuItem.LEFT);

      plotInBox.addItemListener(this);
      plotFwBox.addItemListener(this);
      plotDnBox.addItemListener(this);
      plotIGBox.addItemListener(this);
      plotFGBox.addItemListener(this);
      plotGBox.addItemListener(this);

      sdPlots = new JPanel(new GridLayout(2, 3));
      sdPlots.add(plotInBox);
      sdPlots.add(plotFwBox);
      sdPlots.add(plotDnBox);
      sdPlots.add(plotGBox); 
      sdPlots.add(plotIGBox); 
      sdPlots.add(plotFGBox);
   

      TitledBorder jrbBorder = new TitledBorder(new LineBorder(myBlue), "Spectral Plots");
      jrbBorder.setTitleColor(Color.BLACK);
      sdPlots.setBorder(jrbBorder);

      plotIn = false; plotFw = false; plotFu = false;
      plotG = false; plotGI = true; plotGF = true;

      JPanel plots = new JPanel(new GridLayout(2, 3));
      plots.add(seriesBox);
      plots.add(foreBox);
      plots.add(seasBox);
      plots.add(trendBox);
      plots.add(saBox);
      plots.add(cycleBox);

      TrendModel = new JRadioButton("Trend");  TrendModel.setHorizontalTextPosition(JMenuItem.LEFT);
      TrendModel.setToolTipText("Decomposition of SARIMA model into trend signal and seasonal+irregular noise model.");
      TrendModel.setSelected(true);
      SeasonalModel = new JRadioButton("Seasonal"); SeasonalModel.setHorizontalTextPosition(JMenuItem.LEFT);
      SeasonalModel.setToolTipText("Decomposition of SARIMA model into seasonal signal and trend+irregular noise model.");
      TrendIrregModel = new JRadioButton("Trend-Irr");  TrendIrregModel.setHorizontalTextPosition(JMenuItem.LEFT);
      TrendIrregModel.setToolTipText("Decomposition of SARIMA model into trend+irregular signal and seasonal noise model.");
      IrregModel = new JRadioButton("Irregular");  IrregModel.setHorizontalTextPosition(JMenuItem.LEFT);
      IrregModel.setToolTipText("Decomposition of SARIMA model into irregular/cycle signal and trend-seasonal noise model.");

      ButtonGroup group = new ButtonGroup();
      group.add(TrendModel);
      group.add(SeasonalModel);
      group.add(IrregModel);
      group.add(TrendIrregModel);

      TrendModel.addItemListener(this);
      SeasonalModel.addItemListener(this);
      TrendIrregModel.addItemListener(this);
      IrregModel.addItemListener(this);

      JPanel specPlots = new JPanel(new GridLayout(2, 2));
      specPlots.add(TrendModel);
      specPlots.add(SeasonalModel);
      specPlots.add(IrregModel); 
      specPlots.add(TrendIrregModel);


      TitledBorder specBorder = new TitledBorder(new LineBorder(myBlue), "Signal Component");
      specBorder.setTitleColor(Color.BLACK);
      specPlots.setBorder(specBorder);


      TitledBorder plotBorder = new TitledBorder(new LineBorder(myBlue), "Series Components");
      plotBorder.setTitleColor(Color.BLACK);
      plots.setBorder(plotBorder);

   
      Box hBoxTop = Box.createHorizontalBox();
      hBoxTop.add(plots);
      hBoxTop.add(sdPlots);
      hBoxTop.add(specPlots); 

      return hBoxTop;

  }


  /*-------------------------------------------------------------
      Create Plottting Box
  ---------------------------------------------------------------*/


  public JPanel createModelPanel()
  {
       //------Creat Model Selection Panel-------------------------
       JPanel modelPanel = new JPanel();
       TitledBorder modelBorder = new TitledBorder(new LineBorder(myBlue), "Model");
       modelBorder.setTitleColor(Color.BLACK);
       modelPanel.setBorder(modelBorder);

       p = new JComboBox<String>();
       p.addItem("0");
       p.addItem("1");
       p.addItem("2");
       p.addActionListener(this);
       q = new JComboBox<String>();
       q.addItem("0");
       q.addItem("1");
       q.addItem("2");
       q.setSelectedIndex(1);
       q.addActionListener(this);

       P = new JComboBox<String>();
       P.addItem("0");
       P.addItem("1");
       P.addActionListener(this);
       Q = new JComboBox<String>();
       Q.addItem("0");
       Q.addItem("1");
       Q.setSelectedIndex(1);
       Q.addActionListener(this);

       //------------------------------ Set up regression components ---------------
       easterDay = new JComboBox<String>(); easterDay.setToolTipText("Check and adjust for Easter-Day effects in data using x-day Easter regressors."); 
       easterDay.addItem("1"); easterDay.addItem("8"); easterDay.addItem("15"); easterDay.setSelectedIndex(0); 
       easterDay.addActionListener(this);
       tdBox = new JCheckBox("TD:"); tdBox.setHorizontalTextPosition(JLabel.LEFT);
       outlierBox = new JCheckBox("Outliers:"); 
       outlierBox.setHorizontalTextPosition(JLabel.LEFT); 

       easterBox = new JCheckBox("Easter:"); easterBox.setHorizontalTextPosition(JLabel.LEFT); 
       transBox = new JCheckBox("BC-Trans:"); transBox.addItemListener(this); transBox.setHorizontalTextPosition(JLabel.LEFT);
       transBox.setToolTipText("Set Box-Cox Transformation"); 
       easterBox.setToolTipText("Check and adjust for Easter-Day effects in data using x-day Eastor regressors."); 
       tdBox.setToolTipText("Check and adjust for Trading-Day effects in data.");
       outlierBox.setToolTipText("Check and adjust for outliers in data.");

       tdBox.addItemListener(this); outlierBox.addItemListener(this); easterBox.addItemListener(this);
 
       //--------------------------------------------------------------------------------------------------------
   
      
       //modelPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 15, 5));
       //modelPanel.add(pLabel);
       //modelPanel.add(p);
       //modelPanel.add(qLabel);
       //modelPanel.add(q);
       //modelPanel.add(PLabel);
       //modelPanel.add(P);
       //modelPanel.add(QLabel);
       //modelPanel.add(Q);       
       //----------------------------------------------------------

       modelPanel.setLayout(setUpModelLayout2(modelPanel)); 

       return modelPanel;
  } 


  /*----------------------------------------------------------------
    Setup panel for model order and regression  
  -----------------------------------------------------------------*/

  public GroupLayout setUpModelLayout2(JPanel pmodelPanel)
  {

       JLabel armalab = new JLabel("ARMA Orders:");
       JLabel reglab = new JLabel("Regression:"); reglab.setToolTipText("Use regressors to adjust for certain deterministic effects in data.");
       JLabel eastlab = new JLabel("Days:"); 
       //JLabel spacelab = new JLabel("");

       JPanel botP = new JPanel();
       GroupLayout bot = new GroupLayout(botP);
       bot.setAutoCreateGaps(true);
       bot.setAutoCreateContainerGaps(true);  
       bot.setHorizontalGroup(bot.createSequentialGroup()
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(reglab))
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(transBox))
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(outlierBox))
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(tdBox))
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(easterBox))
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(eastlab))  
              .addGroup(bot.createParallelGroup()                  
                 .addComponent(easterDay)));
       bot.setVerticalGroup(bot.createSequentialGroup()
              .addGroup(bot.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(bot.createSequentialGroup()
                  .addGroup(bot.createParallelGroup(GroupLayout.Alignment.BASELINE)     
                   .addComponent(reglab) 
                   .addComponent(transBox)           
                   .addComponent(outlierBox)
                   .addComponent(tdBox)
                   .addComponent(easterBox)
                   .addComponent(eastlab)
                   .addComponent(easterDay)))));


       JPanel topP = new JPanel();
       GroupLayout top = new GroupLayout(topP);
       top.setAutoCreateGaps(true);
       top.setAutoCreateContainerGaps(true);  
       top.setHorizontalGroup(top.createSequentialGroup()
              .addGroup(top.createParallelGroup()                  
                 .addComponent(armalab))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(pLabel))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(p))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(qLabel))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(q))  
              .addGroup(top.createParallelGroup()                  
                 .addComponent(PLabel))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(P))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(QLabel))
              .addGroup(top.createParallelGroup()                  
                 .addComponent(Q)));
       top.setVerticalGroup(top.createSequentialGroup()
              .addGroup(top.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(top.createSequentialGroup()
                  .addGroup(top.createParallelGroup(GroupLayout.Alignment.BASELINE)     
                   .addComponent(armalab)            
                   .addComponent(p)
                   .addComponent(pLabel)
                   .addComponent(q)
                   .addComponent(qLabel)
                   .addComponent(P)
                   .addComponent(PLabel)
                   .addComponent(Q)
                   .addComponent(QLabel)))));

       GroupLayout pLayout = new GroupLayout(pmodelPanel);
       pLayout.setAutoCreateGaps(true);
       pLayout.setAutoCreateContainerGaps(true);
       pLayout.setHorizontalGroup(pLayout.createSequentialGroup()
             .addGroup(pLayout.createParallelGroup() 
                 .addComponent(topP)
                 .addComponent(botP)));
       pLayout.setVerticalGroup(pLayout.createSequentialGroup()
              .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(pLayout.createSequentialGroup()
                  .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)   
                     .addComponent(topP))
                  .addGroup(pLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)  
                     .addComponent(botP))))); 



     return pLayout;

  }

  /*--------------------------------------------------------------------
    MLE Panel
  ----------------------------------------------------------------------*/
  public JPanel createMLEPanel()
  {

      JPanel mlePanel = new JPanel();
      TitledBorder mleBorder = new TitledBorder(new LineBorder(myBlue), "MLE Values");
      mleBorder.setTitleColor(Color.BLACK);
      mlePanel.setBorder(mleBorder);

      mlePanel.setLayout(setUpMLEValuesLayout(mlePanel));
      mlePanel.add(new JLabel("\u03D5_1"));
      mlePanel.add(textphi1);
      mlePanel.add(new JLabel("\u03D1_1"));
      mlePanel.add(texttheta1);
      mlePanel.add(new JLabel("\u0398"));
      mlePanel.add(textTheta);
      mlePanel.add(new JLabel("Inn.Var."));
      mlePanel.add(textInnvar);

      mlePanel.add(new JLabel("\u03D5_2"));
      mlePanel.add(textphi2);      
      mlePanel.add(new JLabel("\u03D1_2"));
      mlePanel.add(texttheta2);
      mlePanel.add(new JLabel("\u03A6"));
      mlePanel.add(textPhi);
      mlePanel.add(new JLabel("Sig.Var."));
      mlePanel.add(textSigvar);
    
      return mlePanel;
  }

  public JPanel createDiagPanel()
  {

       //-------Create Diagnostics Box------------------------
       JPanel diagPanel = new JPanel();
       TitledBorder diagBorder = new TitledBorder(new LineBorder(myBlue), "GoF Signal Extraction and Likelihood Statistics");
       diagBorder.setTitleColor(Color.BLACK);
       diagPanel.setBorder(diagBorder);
       
       JLabel trendL = new JLabel(" Trend:"); JLabel seasL = new JLabel("Seasonal:");
       JLabel irrL = new JLabel("Irregular:"); JLabel trendirrL = new JLabel("Trend-Irr.:");
       JLabel aicL = new JLabel("      AIC:"); JLabel aiccL = new JLabel("      AICC:");
       JLabel hnqL = new JLabel(" HNQuin:"); JLabel bicL = new JLabel("         BIC:");


       text10 = new JTextField(5);
       text10.setText("" + df.format(smc.diagnostics[0]));   
       text11 = new JTextField(5);      
       text11.setText("" + df.format(smc.diagnostics[1])); 
       text12 = new JTextField(5); 
       text12.setText("" + df.format(smc.diagnostics[2]));   
       text13 = new JTextField(5);
       text13.setText("" + df.format(smc.diagnostics[3]));
    
       textLk1.setText("" + df.format(smc.lkhs[0]));
       textLk2.setText("" + df.format(smc.lkhs[1]));
       textLk3.setText("" + df.format(smc.lkhs[2]));
       textLk4.setText("" + df.format(smc.lkhs[3]));
      

       GroupLayout paramLayout = new GroupLayout(diagPanel);
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup() 
                 .addComponent(trendL)
                 .addComponent(aicL))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text10)
                 .addComponent(textLk1))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(seasL)
                 .addComponent(aiccL))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text11)
                 .addComponent(textLk2))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(irrL)
                 .addComponent(hnqL))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text12)
                 .addComponent(textLk3))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(trendirrL)
                 .addComponent(bicL))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(text13)
                 .addComponent(textLk4)));

       paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(paramLayout.createSequentialGroup()
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(trendL)
                   .addComponent(text10)
                   .addComponent(seasL)
                   .addComponent(text11)
                   .addComponent(irrL)
                   .addComponent(text12)
                   .addComponent(trendirrL)               
                   .addComponent(text13)) 
                 .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)      
                   .addComponent(aicL)
                   .addComponent(textLk1)
                   .addComponent(aiccL)
                   .addComponent(textLk2)
                   .addComponent(hnqL)
                   .addComponent(textLk3)
                   .addComponent(bicL)               
                   .addComponent(textLk4)))));


       diagPanel.setLayout(paramLayout);

        
       return diagPanel;
  }

  public JPanel createX13Panel()
  {
       JPanel x13estPanel = new JPanel();
       TitledBorder x13estBorder = new TitledBorder(new LineBorder(myBlue), "Ljung-Box Stats");
       x13estBorder.setTitleColor(Color.BLACK);
       x13estPanel.setBorder(x13estBorder);
     
       JLabel lb0 = new JLabel("LB-Lag 0:");
       JLabel lb1 = new JLabel("LB-Lag 1:");
       JLabel lb12 = new JLabel("LB-Lag 12:");
       JLabel lb24 = new JLabel("LB-Lag 16:");


       GroupLayout paramLayout = new GroupLayout(x13estPanel);
       paramLayout.setAutoCreateGaps(true);
       paramLayout.setAutoCreateContainerGaps(true);
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup() 
                 .addComponent(lb0).addComponent(lb12))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(textLB0).addComponent(textLB12))
              .addGroup(paramLayout.createParallelGroup() 
                 .addComponent(lb1).addComponent(lb24))
              .addGroup(paramLayout.createParallelGroup()
                 .addComponent(textLB1).addComponent(textLB24)));

       paramLayout.setVerticalGroup(paramLayout.createSequentialGroup()
              .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                .addGroup(paramLayout.createSequentialGroup()
                  .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lb0).addComponent(textLB0).addComponent(lb1).addComponent(textLB1))
                 .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                   .addComponent(lb12).addComponent(textLB12).addComponent(lb24).addComponent(textLB24)))));

       x13estPanel.setLayout(paramLayout);
       textLB0.setText("" + df.format(smc.lbv[0]));
       textLB1.setText("" + df.format(smc.lbv[1])); 
       textLB12.setText("" + df.format(smc.lbv[2]));
       textLB24.setText("" + df.format(smc.lbv[3]));
 
       return x13estPanel;
  }

  public JPanel createPModelPanel()
  {

       // ------- Define all the pseudo-panels 
       JPanel pmodelPanel = new JPanel();
       TitledBorder pmodelBorder = new TitledBorder(new LineBorder(myBlue), "Model");
       pmodelBorder.setTitleColor(Color.BLACK);
       pmodelPanel.setBorder(pmodelBorder);

       //============= model selection tools ========================
       tp = new JComboBox<String>();
       tp.addItem("0");
       tp.addItem("1");
       tp.addItem("2");
       tp.addActionListener(this);
       tq = new JComboBox<String>();
       tq.addItem("0");
       tq.addItem("1");
       tq.addItem("2");
       tq.setSelectedIndex(1);
       tq.addActionListener(this);

       tP = new JComboBox<String>();
       tP.addItem("0");
       tP.addItem("1");
       tP.addActionListener(this);
       tQ = new JComboBox<String>();
       tQ.addItem("0");
       tQ.addItem("1");
       tQ.setSelectedIndex(1);
       tQ.addActionListener(this);     

       //============= model selection tools ========================
       pp = new JComboBox<String>();
       pp.addItem("0");
       pp.addItem("1");
       pp.addItem("2");
       pp.addActionListener(this);
       pq = new JComboBox<String>();
       pq.addItem("0");
       pq.addItem("1");
       pq.addItem("2");
       pq.setSelectedIndex(1);
       pq.addActionListener(this);

       pP = new JComboBox<String>();
       pP.addItem("0");
       pP.addItem("1");
       pP.addActionListener(this);
       pQ = new JComboBox<String>();
       pQ.addItem("0");
       pQ.addItem("1");
       pQ.setSelectedIndex(1);
       pQ.addActionListener(this);

       estLabel1 = new JLabel("Simulated Model:  ");
       estLabel2 = new JLabel("Estimated Model:  ");
       pmodelPanel.setLayout(this.setUpModelLayout(pmodelPanel)); 

       return pmodelPanel;
   }

   public JPanel createPValPanel()
   {
       
       JPanel pvalPanel = new JPanel();
       TitledBorder pvalBorder = new TitledBorder(new LineBorder(myBlue), "Pseudo-True Values");
       pvalBorder.setTitleColor(Color.BLACK);
       pvalPanel.setBorder(pvalBorder);
       pvalPanel.setLayout(setUpPseudoValuesLayout(pvalPanel));
        
       pvalPanel.add(new JLabel("\u03D5_1"));
       pvalPanel.add(ptextphi1);
       pvalPanel.add(new JLabel("\u03D1_1"));
       pvalPanel.add(ptexttheta1);
       pvalPanel.add(new JLabel("\u0398"));
       pvalPanel.add(ptextTheta);
       pvalPanel.add(new JLabel("InnVar"));
       pvalPanel.add(ptextInnvar);

       pvalPanel.add(new JLabel("\u03D5_2"));
       pvalPanel.add(ptextphi2);      
       pvalPanel.add(new JLabel("\u03D1_2"));
       pvalPanel.add(ptexttheta2);
       pvalPanel.add(new JLabel("\u03A6"));
       pvalPanel.add(ptextPhi);
       pvalPanel.add(new JLabel("minKL"));
       pvalPanel.add(ptextMinKL);

       return pvalPanel;
   }

   public JPanel createPDiagPanel()
   {

       JPanel pdiagPanel = new JPanel();
       TitledBorder pdiagBorder = new TitledBorder(new LineBorder(myBlue), "Pseudo-True Diagnostics");
       pdiagBorder.setTitleColor(Color.BLACK);
       pdiagPanel.setBorder(pdiagBorder);

       ptext10 = new JTextField(4);  
       ptext11 = new JTextField(4);         
       ptext12 = new JTextField(4);       
       ptext13 = new JTextField(4);
       
       //----------------------------------------
      
       pdiagPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 3));
       pdiagPanel.add(new JLabel(" Trend:"));
       pdiagPanel.add(ptext10);
       pdiagPanel.add(Box.createHorizontalGlue());
       pdiagPanel.add(new JLabel("Seasonal:"));
       pdiagPanel.add(ptext11);
       pdiagPanel.add(new JLabel("Irregular:"));
       pdiagPanel.add(ilabel);
       pdiagPanel.add(ptext12);
       pdiagPanel.add(Box.createHorizontalGlue());
       pdiagPanel.add(new JLabel("Trend-Irr.:"));
       pdiagPanel.add(ptext13);
 
       return pdiagPanel;
 
   }

   public JPanel createPeffPanel()
   {

       JPanel effPanel = new JPanel();
       TitledBorder effBorder = new TitledBorder(new LineBorder(myBlue), "Signal Extraction Efficacies");
       effBorder.setTitleColor(Color.BLACK);
       effPanel.setBorder(effBorder);

       ptextefft = new JTextField(4);         
       ptexteffs = new JTextField(4);           
       ptexteffi = new JTextField(4);          
       ptexteffti = new JTextField(4);

       effPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 3));
       effPanel.add(tlabel);
       effPanel.add(ptextefft);
       effPanel.add(Box.createHorizontalGlue());
       effPanel.add(slabel);
       effPanel.add(ptexteffs);
       effPanel.add(Box.createHorizontalGlue());
       effPanel.add(ilabel);
       effPanel.add(ptexteffi);
       effPanel.add(Box.createHorizontalGlue());
       effPanel.add(tilabel);
       effPanel.add(ptexteffti);
  
       return effPanel;
   }

   public JPanel createPLBPanel()
   {
       JPanel lbPanel = new JPanel();
       TitledBorder lbBorder = new TitledBorder(new LineBorder(myBlue), "Ljung-Box Stats");
       lbBorder.setTitleColor(Color.BLACK);
       lbPanel.setBorder(lbBorder);
       // --------- Add two bottom boxes ------------
       lbPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 5));
       lbPanel.add(new JLabel("LB-Lag 12"));
       lbPanel.add(ptextLB12);
       lbPanel.add(new JLabel("LB-Lag 24"));
       lbPanel.add(ptextLB24);
    
       return lbPanel;

   }

   public JPanel createCritPanel()
   {
       JPanel criteria = new JPanel();
       TitledBorder critBorder = new TitledBorder(new LineBorder(myBlue), "Likelihood Stats");
       critBorder.setTitleColor(Color.BLACK);
       criteria.setBorder(critBorder); 

       criteria.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 5));
       criteria.add(new JLabel("AIC"));
       criteria.add(ptextLk1);
       criteria.add(new JLabel("AICC"));
       criteria.add(ptextLk2);
       criteria.add(new JLabel("HNQuin"));
       criteria.add(ptextLk3);
       criteria.add(new JLabel("BIC"));
       criteria.add(ptextLk4);

       return criteria;

   }

   
   @SuppressWarnings({ "rawtypes", "unchecked" })
private void createQuandlDialog()
   {
   
        quandlPanel = new JPanel();
     

        quandlLable = new JLabel();
        searchQLabel = new JLabel();
        searchQText = new JTextField();
        fromQLable = new JLabel();
        fromQText = new JTextField();
        toQText = new JTextField();
        toQLabel = new JLabel();
        goQuandlButton = new JButton();
        quandlComboBox = new JComboBox();
        quandlFinance = new JCheckBox("YahooFinance");
        quandlFreqLabel = new JLabel("Data Frequency");
        
        quandlLable.setText("Quandl Data Download");

        searchQLabel.setText("Quandl Search");

        searchQText.setText("FRED/ICSA MOODY/WAAAYLD YAHOO/INDEX_GSPC");

        fromQLable.setText("From");

        fromQText.setText("2000-01-01");

        toQText.setText("2014-09-20");

        toQLabel.setText("To");

        goQuandlButton.setText("Go Quandl");

        quandlComboBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        quandlComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Default", "Weekly", "Monthly", "Weekly - no log", "Monthly - no log"}));

        quandlFreqLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        quandlFreqLabel.setText("Data Frequency");

        quandlFinance.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        quandlFinance.setText("YahooFinance");
        quandlFinance.setHorizontalTextPosition(javax.swing.SwingConstants.LEADING);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(quandlPanel);
        quandlPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(quandlLable)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(searchQLabel)
                            .addComponent(fromQLable))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(searchQText)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(goQuandlButton, javax.swing.GroupLayout.PREFERRED_SIZE, 132, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(fromQText, javax.swing.GroupLayout.PREFERRED_SIZE, 109, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(18, 18, 18)
                                        .addComponent(toQLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(toQText, javax.swing.GroupLayout.PREFERRED_SIZE, 109, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(18, 18, 18)
                                        .addComponent(quandlFreqLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(quandlComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(quandlFinance)))
                                .addGap(0, 41, Short.MAX_VALUE)))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(quandlLable)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(searchQLabel)
                    .addComponent(searchQText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fromQLable)
                    .addComponent(fromQText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(toQLabel)
                    .addComponent(toQText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(quandlComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(quandlFreqLabel)
                    .addComponent(quandlFinance))
                .addGap(18, 18, 18)
                .addComponent(goQuandlButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
                  


      goQuandlButton.addActionListener(new ActionListener() 
      {
            public void actionPerformed(ActionEvent evt) 
            {
                
              String[] symbols; 
              String delims = "[ ]+";

              symbols = searchQText.getText().split(delims); 
              simulate.getQuandlData(symbols, fromQText.getText(), toQText.getText(), quandlComboBox.getSelectedIndex(), quandlFinance.isSelected());
   
            }
      });        
   
   }
   

   @SuppressWarnings({ "rawtypes", "unchecked" })
private void createGoogleIntradayDialog()
   {

        googleIntradayPanel = new JPanel();
 
        googSymbolsLabel = new JLabel();
        googSymbolsText = new JTextField();
        googExchLabel = new JLabel();
        googExchText = new JTextField();
        googFreqBox = new JComboBox();
        googFreqLabel = new JLabel();
        googDaysLabel = new JLabel();
        googDaysBox = new JComboBox();
        googLogReturnsBox = new JCheckBox();
        googVolBox = new JCheckBox();
        googHiLowBox = new JCheckBox();
        googHiLoDiffBox = new JCheckBox();
        googHigherFreq = new JCheckBox();
        googHigherFreqPeriod = new JComboBox();
        googHigherFreqPeriodLabel = new JLabel();

        googIntraday = new JButton("Download Market Data");

        googSymbolsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googSymbolsLabel.setText("Symbol(s)");

        googSymbolsText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googSymbolsText.setText("GOOG AAPL");

        googExchLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googExchLabel.setText("Exchange(s)");

        googExchText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googExchText.setText("NASDAQ NASDAQ");

        googFreqBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googFreqBox.setModel(new DefaultComboBoxModel(new String[] { "1", "30", "60", "300", "600", "900", "1200", "1800", "3600" }));
        googFreqBox.setSelectedIndex(2);

        googHigherFreqPeriodLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googHigherFreqPeriodLabel.setText("Period:");
        
        googHigherFreqPeriod.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googHigherFreqPeriod.setModel(new DefaultComboBoxModel(new String[] { "5", "10", "15", "20" }));
        googHigherFreqPeriod.setSelectedIndex(0);        
        
        googFreqLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googFreqLabel.setText("Frequency (in sec)");

        googDaysLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googDaysLabel.setText("Number Days");

        googDaysBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googDaysBox.setModel(new DefaultComboBoxModel(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20" }));
        googDaysBox.setSelectedIndex(1);

        googLogReturnsBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googLogReturnsBox.setText("Get LogReturns");

        googVolBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googVolBox.setText("Get Volume");

        googHiLowBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googHiLowBox.setText("Get High-Low");

        googHiLoDiffBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googHiLoDiffBox.setText("Get High-Low Difference");

        googHigherFreq.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        googHigherFreq.setText("Get Higher-Frequency Price Data");
        
        GroupLayout layout = new GroupLayout(googleIntradayPanel);
        googleIntradayPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(24, 24, 24)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(googLogReturnsBox)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(googVolBox)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(googHiLowBox)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(googHiLoDiffBox))
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(googSymbolsLabel)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(googSymbolsText, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGroup(layout.createSequentialGroup()
                                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(googExchLabel)
                                        .addComponent(googFreqLabel))
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addGroup(layout.createSequentialGroup()
                                            .addComponent(googFreqBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                            .addComponent(googDaysLabel)
                                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                            .addComponent(googDaysBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                        .addComponent(googExchText, javax.swing.GroupLayout.PREFERRED_SIZE, 250, javax.swing.GroupLayout.PREFERRED_SIZE))))))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(133, 133, 133)
                        .addComponent(googIntraday))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(googHigherFreq)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(googHigherFreqPeriodLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(googHigherFreqPeriod)))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(36, 36, 36)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(googSymbolsLabel)
                    .addComponent(googSymbolsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(googExchLabel)
                    .addComponent(googExchText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(googFreqBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(googFreqLabel)
                    .addComponent(googDaysBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(googDaysLabel))
                .addGap(29, 29, 29)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(googLogReturnsBox)
                    .addComponent(googVolBox)
                    .addComponent(googHiLowBox)
                    .addComponent(googHiLoDiffBox))
                .addGap(18, 18, 18)
                .addComponent(googIntraday)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(googHigherFreq)
                    .addComponent(googHigherFreqPeriodLabel)
                    .addComponent(googHigherFreqPeriod))
                .addContainerGap(40, Short.MAX_VALUE))
        );
       
        googIntraday.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                
              String[] symbols; String[] exchange;
              String interval; String noDays; String delims = "[ ]+";
 
              String sym = googSymbolsText.getText();
              String exch = googExchText.getText();

              symbols = sym.split(delims); 
              exchange = exch.split(delims);

              int symbtoks = symbols.length;
              int exchtoks = exchange.length;

              int period = 5; 
              
              //if(googHigherFreq.isSelected()) {period = (int)googHigherFreqPeriod.getSelectedItem();}
              
              interval = (String)googFreqBox.getSelectedItem();             
              noDays = (String)googDaysBox.getSelectedItem();

              if(symbtoks == exchtoks)
              {  
                simulate.getIntradayGoogleData(symbols, exchange, interval, noDays, true, googLogReturnsBox.isSelected(), 
                  googVolBox.isSelected(), googHiLowBox.isSelected(), googHiLoDiffBox.isSelected(), googHigherFreq.isSelected(), period);
              }
 
            }
        });


   }




    private void createMarketDialog() 
    {

        marketChoosePanel = new JPanel();

        high_lowData = new JCheckBox();
        high_lowDataReturn = new JCheckBox();
        volumeData = new JCheckBox();
        enterMarket = new JButton("Download Market Data");
        instrumLabel = new JLabel();
        instrumText = new JTextField();
        hoursLabel = new JLabel();
        startLabel = new JLabel();
        endLabel = new JLabel();
        frequencyLabel = new JLabel();
        hoursCombo = new JComboBox<String>(new String[] { "US Market Hours", "24 Hours" });
        frequencyCombo = new JComboBox<String>(new String[] { "Minute", "3-Minute", "5-Minute", "10-Minute", "15-Minute", "30-Minute", "Hourly", "Daily", "Weekly", "Monthly" });
        frequencyCombo.setSelectedIndex(4); 
        yearText = new JTextField();
        yearLabel = new JLabel();
        monthLabel = new JLabel();
        monthText = new JTextField();
        dayLabel = new JLabel();
        dayText = new JTextField();
        mainLabel = new JLabel();

        yearEndText = new JTextField();
        yearEndLabel = new JLabel();
        monthEndLabel = new JLabel();
        monthEndText = new JTextField();
        dayEndLabel = new JLabel();
        dayEndText = new JTextField();        
        
        fromYahoo = new JCheckBox();

        //setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

        instrumLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        instrumLabel.setText("Symbol(s)");

        instrumText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        instrumText.setText("GOOG AAPL");

        hoursLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        hoursLabel.setText("Hours");

        startLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        startLabel.setText("Start Date");

        endLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        endLabel.setText("End Date");        
        
        frequencyLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        frequencyLabel.setText("Frequency");

        hoursCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        logReturns = new JCheckBox("Log-Returns");        
        clearData = new JCheckBox("Clear Data");

        frequencyCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        yearText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        yearText.setText("2013");

        yearLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        yearLabel.setText("Year");

        monthLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        monthLabel.setText("Month");

        monthText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        monthText.setText("01");

        dayLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        dayLabel.setText("Day");

        dayText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        dayText.setText("04");

        yearEndText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        yearEndText.setText("2013");

        yearEndLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        yearEndLabel.setText("Year");

        monthEndLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        monthEndLabel.setText("Month");

        monthEndText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        monthEndText.setText("02");

        dayEndLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        dayEndLabel.setText("Day");

        dayEndText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        dayEndText.setText("01");        
        
        
        mainLabel.setText("Choose Market Instrument");

     
        ActionListener marketListener = new ActionListener() {
         public void actionPerformed(ActionEvent event) 
         {
           String symbols;
           String year,month,day;
           String delims = "[ ]+";
           String[] tokens;
           String date1,date2;
           String time1,time2;
           int hoursIndex, freqIndex; 

           if(event.getSource() == enterMarket)
           {

             symbols = instrumText.getText();
             year    = yearText.getText();
             month   = monthText.getText(); 
             day     = dayText.getText();
             time1 = "09:30";  
             time2 = "16:00";
             tokens = symbols.split(delims); 
             hoursIndex = hoursCombo.getSelectedIndex();
             freqIndex  = frequencyCombo.getSelectedIndex();
             date1 = year + "-" + month + "-" + day;
             
             year    = yearEndText.getText();
             month   = monthEndText.getText(); 
             day     = dayEndText.getText();             
             
             date2 = year + "-" + month + "-" + day;
             
             //simulate.getMarketData(tokens[0], date1, hoursIndex==0, time1, time2, freqIndex);

             if(!fromYahoo.isSelected())
             {
                          
                simulate.getMarketDataGeneral(tokens, date1, date2, hoursIndex==0, time1, time2, 
                          freqIndex, logReturns.isSelected(), volumeData.isSelected(), high_lowData.isSelected(), 
                           high_lowDataReturn.isSelected(), clearData.isSelected());                
                          
             }
             else if(fromYahoo.isSelected())
             {
                simulate.getMarketDataYahoo(tokens, date1, hoursIndex==0, time1, time2, 
                          freqIndex, logReturns.isSelected(), volumeData.isSelected(), clearData.isSelected(), high_lowData.isSelected(), 
                           high_lowDataReturn.isSelected());
             }
             else if(!clearData.isSelected())
             {
               mdfa.getHFPeriodogram(tokens, date1, date2, hoursIndex==0, time1, time2, freqIndex);
               periodoWeightBox.setEnabled(true);
             }
            }   
         }
       };
       enterMarket.addActionListener(marketListener);

        GroupLayout layout = new GroupLayout(marketChoosePanel);
        marketChoosePanel.setLayout(layout);


        logReturns.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        logReturns.setText("Log Returns");
        logReturns.setToolTipText("Uploads log returns of asset price ");        

        clearData.setFont(new Font("Ubuntu",0,12));
        clearData.setText("New Data Set");
        clearData.setToolTipText("Clears all data in the Data Control Module for new data");

        volumeData.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        volumeData.setText("Get Volume Data");
        volumeData.setText("Volume Data");

        fromYahoo.setFont(new Font("Ubuntu",0,12)); 
        fromYahoo.setText("Yahoo Source"); 
        fromYahoo.setToolTipText("Downloads the instrument data from Yahoo (only in Daily prices)");        

        high_lowData.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        high_lowData.setText("Get High-Low Data");
        high_lowData.setToolTipText("Downloads the high-low of instrument data");  
        
        high_lowDataReturn.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        high_lowDataReturn.setText("Get High-Low Return Data");        
        high_lowDataReturn.setToolTipText("Downloads the high-low return of instrument data");  

        marketChoosePanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(instrumLabel)
                            .addComponent(hoursLabel)
                            .addComponent(frequencyLabel))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(hoursCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(instrumText, GroupLayout.PREFERRED_SIZE, 283, GroupLayout.PREFERRED_SIZE)
                            .addComponent(frequencyCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                            .addComponent(mainLabel)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(clearData)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(logReturns)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(volumeData)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(fromYahoo))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(high_lowData)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(high_lowDataReturn))))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(startLabel)
                                .addGap(18, 18, 18)
                                .addComponent(yearLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(yearText)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(monthLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(monthText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                            .addGroup(layout.createSequentialGroup()
                                .addGap(107, 107, 107)
                                .addComponent(enterMarket))
                           .addGroup(layout.createSequentialGroup()
                                .addComponent(endLabel)
                                .addGap(18, 18, 18)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(yearEndLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(yearEndText)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(monthEndLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(monthEndText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(dayLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(dayText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(dayEndLabel)
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(dayEndText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))))
                .addContainerGap(98, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(13, 13, 13)
                .addComponent(mainLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(instrumLabel)
                    .addComponent(instrumText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(dayText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(dayLabel))
                    .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(startLabel)
                        .addComponent(yearText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(yearLabel)
                        .addComponent(monthText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(monthLabel)))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(dayEndText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(dayEndLabel)
                    .addComponent(monthEndText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(monthEndLabel)
                    .addComponent(yearEndText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(yearEndLabel)
                    .addComponent(endLabel))                        
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(hoursLabel, GroupLayout.PREFERRED_SIZE, 27, GroupLayout.PREFERRED_SIZE)
                    .addComponent(hoursCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(frequencyLabel)
                    .addComponent(frequencyCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(clearData)
                    .addComponent(logReturns)
                    .addComponent(volumeData) 
                    .addComponent(fromYahoo))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(high_lowData)
                    .addComponent(high_lowDataReturn))                    
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(enterMarket)
                .addContainerGap(21, Short.MAX_VALUE))
        );



        pack();
    }



  @SuppressWarnings({ "rawtypes", "unchecked" })
public void initRealizedVolControl() 
  {

    JLabel jLabel10,jLabel11,jLabel13,jLabel5,jLabel7,jLabel9,jLabel4,jLabel6,jLabel12;

        realizedVolPanel = new JPanel();
        JLabel instrumLabel = new JLabel();
        instrumText2 = new JTextField();
        JLabel startLabel = new JLabel();
        startYear = new JTextField();
        JLabel yearLabel = new JLabel();
        jLabel4 = new JLabel();
        startMonth = new JTextField();
        jLabel5 = new JLabel();
        startDay = new JTextField();
        jLabel6 = new JLabel();
        jLabel7 = new JLabel();
        endYear = new JTextField();
        JLabel endLabel = new JLabel();
        endMonth = new JTextField();
        jLabel9 = new JLabel();
        endDay = new JTextField();
        kernelCombo = new JComboBox();
        jLabel10 = new JLabel();
        timeCombo = new JComboBox();
        jLabel11 = new JLabel();
        
        periodBar = new JScrollBar(JScrollBar.HORIZONTAL,1,1,1,30);
        periodText = new JTextField(); periodText.setText("1"); 
        periodBar.setPreferredSize(new Dimension(80, 15));
        periodBar.addAdjustmentListener( new AdjustmentListener()  {
           public void adjustmentValueChanged(AdjustmentEvent e) {
             periodText.setText(""+periodBar.getValue());
        }});
        
        mdfaBox = new JCheckBox("Prepare Data for MDFA Module"); mdfaBox.setSelected(false); mdfaBox.setHorizontalTextPosition(JLabel.LEFT);
        mdfaBox.setToolTipText("Prepares the realized volatility measure for the MDFA module");
        
        periodCombo = new JPanel();
        GroupLayout paramLayout = new GroupLayout(periodCombo); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(periodBar).addComponent(periodText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(periodBar).addComponent(periodText)));
        periodCombo.setLayout(paramLayout);        
        
        
        
        jLabel12 = new JLabel();
        lagCombo = new JComboBox();
        jLabel13 = new JLabel();
        enterCompute = new JButton();

        instrumLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        instrumLabel.setText("Choose Instrument");

        instrumText2.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        instrumText2.setText("GOOG.O");

        startLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        startLabel.setText("Start Date:");

        startYear.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        startYear.setText("2012");

        yearLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        yearLabel.setText("ccyy");

        jLabel4.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel4.setText("mm");

        startMonth.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        startMonth.setText("03");

        jLabel5.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel5.setText("dd");

        startDay.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        startDay.setText("01");

        jLabel6.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel6.setText("dd");

        jLabel7.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel7.setText("ccyy");

        endYear.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        endYear.setText("2012");

        endLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        endLabel.setText("End Date:");

        endMonth.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        endMonth.setText("06");

        jLabel9.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel9.setText("mm");

        endDay.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        endDay.setText("19");

        kernelCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        kernelCombo.setModel(new DefaultComboBoxModel(new String[] { hfreq.kernels[0], hfreq.kernels[1],hfreq.kernels[2], hfreq.kernels[3], hfreq.kernels[4],
                                   hfreq.kernels[5], hfreq.kernels[6], hfreq.kernels[7], hfreq.kernels[8], hfreq.kernels[9], hfreq.kernels[10], hfreq.kernels[11] }));

        jLabel10.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel10.setText("Kernel:");

        timeCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        timeCombo.setModel(new DefaultComboBoxModel(new String[] { "seconds", "minutes", "hours"}));

        jLabel11.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel11.setText("Time Scale:");


        jLabel12.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel12.setText("Period:");

        lagCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lagCombo.setModel(new DefaultComboBoxModel(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"}));

        jLabel13.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        jLabel13.setText("Lags:");

        enterCompute.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        enterCompute.setText("Compute Realized Volatility");

  
        
       GroupLayout layout = new javax.swing.GroupLayout(realizedVolPanel);
        realizedVolPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(startLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(yearLabel)
                                        .addGap(3, 3, 3)
                                        .addComponent(startYear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jLabel4)
                                        .addGap(9, 9, 9)
                                        .addComponent(startMonth, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(6, 6, 6)
                                        .addComponent(jLabel5)
                                        .addGap(8, 8, 8)
                                        .addComponent(startDay, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(endLabel))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(instrumLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(instrumText2, javax.swing.GroupLayout.PREFERRED_SIZE, 201, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(enterCompute)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jLabel7)
                                        .addGap(3, 3, 3)
                                        .addComponent(endYear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jLabel9)
                                        .addGap(9, 9, 9)
                                        .addComponent(endMonth, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(6, 6, 6)
                                        .addComponent(jLabel6)
                                        .addGap(8, 8, 8)
                                        .addComponent(endDay, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE))))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jLabel10)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(kernelCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(jLabel11)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(timeCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jLabel12)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(periodBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(periodText, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jLabel13)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(lagCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGap(38, 38, 38))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(mdfaBox)
                        .addContainerGap())))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(35, 35, 35)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(instrumLabel)
                            .addComponent(instrumText2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(enterCompute)))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(startLabel)
                    .addComponent(startYear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(yearLabel)
                    .addComponent(jLabel4)
                    .addComponent(jLabel5)
                    .addComponent(startMonth, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(startDay, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(endLabel)
                    .addComponent(endYear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel7)
                    .addComponent(jLabel9)
                    .addComponent(jLabel6)
                    .addComponent(endMonth, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(endDay, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(30, 30, 30)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(kernelCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel10)
                    .addComponent(timeCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel11)
                    .addComponent(jLabel12)
                    .addComponent(lagCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel13)
                    .addComponent(periodBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(periodText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(mdfaBox)
                .addGap(17, 17, 17))
        );
        
        
        
        
        



       ActionListener enterComputeListener = new ActionListener() {
         public void actionPerformed(ActionEvent event) 
         {
           String symbols;
           String year,month,day;
           int i,ndays; 
        
           String delims = "[ ]+";
           String[] tokens;
           String startDate,endDate;
       
           if(event.getSource() == enterCompute)
           {

             symbols = instrumText2.getText();
             year    = startYear.getText();
             month   = startMonth.getText(); 
             day     = startDay.getText();
             startDate = year + "-" + month + "-" + day;         

             year    = endYear.getText();
             month   = endMonth.getText(); 
             day     = endDay.getText();
             endDate = year + "-" + month + "-" + day;   


             tokens = symbols.split(delims); 
             
             
             //instrums = new String[n_toks];
             //for(i=0;i<n_toks;i++) {instrums[i] = tokens[i];}

             //----- get the daily log-returns

             
            hfreq.getRealizedVolatility(tokens, startDate, endDate, timeCombo.getSelectedIndex(), 
               kernelCombo.getSelectedIndex(), periodBar.getValue(), lagCombo.getSelectedIndex()+1,mdfaBox.isSelected());
               
                 
            ndays = hfreq.getNDays();
            simulate.clearSeries(ndays);
            
            if(mdfaBox.isSelected()) //set price info first
            {simulate.setLogPriceData(hfreq.assets[0].price);}
            
            for(i=0;i<hfreq.n_assets;i++)
            {
            
              if(mdfaBox.isSelected())  //transform realized measure to    
              {
                 double[] realizedRet = getRealizedReturns(hfreq.assets[i].log_return,hfreq.assets[i].realized_vol);
                 simulate.setRealSeries(hfreq.assets[i].log_return);
                 simulate.setRealSeries(realizedRet);
              }
              else
              {
               simulate.setRealSeries(hfreq.assets[i].log_return);
               simulate.setRealSeries(hfreq.assets[i].realized_vol);            
              }
            }           
           }                           
          }      
       };

       enterCompute.addActionListener(enterComputeListener);
   }
 
   
   public double[] getRealizedReturns(double[] logret, double[] realized)
   {
      int i; int n = logret.length; 
      double[] real_ret = new double[n];
      if(logret.length == realized.length)
      {      
        for(i=0;i<n;i++)
        {
          real_ret[i] = Math.sqrt(realized[i]);
          if(logret[i] < 0) {real_ret[i] = -1.0*real_ret[i];}
        }                
      }
      return real_ret;   
   }












   /*-----------------------------------------------------

   ------------------------------------------------------*/
   public JPanel createEMDMap()
   {
 
     JPanel plots2x = new JPanel();
     emdm = new  EMDMap(500,120,nObs,120);
     emdKey =  new EMDMapKey(30,120,120);
     emdScale = new EMDMapScale(500,15);
     freqScale = new EMDFreqScale(15,120);
     hScale = new EMD3dScale(15,120);
         
     //-------- Set Borders ---------------------------------------
     Border raisedetched = BorderFactory.createRaisedBevelBorder();
     BorderFactory.createLoweredBevelBorder();
     emdm.setBorder(raisedetched);
     emdKey.setBorder(raisedetched);
     emdScale.setBorder(raisedetched);
     jplRadio.setBorder(raisedetched);

     GroupLayout layout = new GroupLayout(plots2x);
     layout.setAutoCreateGaps(true);
     layout.setAutoCreateContainerGaps(true);

     layout.setHorizontalGroup(
        layout.createSequentialGroup()
         .addComponent(jplRadio,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
         .addComponent(freqScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
         .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING,false)
           .addComponent(emdm,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
           .addComponent(emdScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE))
         .addComponent(hScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
         .addComponent(emdKey,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
     );
 
     layout.setVerticalGroup(
        layout.createSequentialGroup()
         .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(jplRadio,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
           .addComponent(freqScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
           .addComponent(emdm,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
           .addComponent(hScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
           .addComponent(emdKey,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE))
        .addComponent(emdScale,GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE,
          GroupLayout.PREFERRED_SIZE)
     );
     plots2x.setLayout(layout);

     return plots2x;
   }


   public boolean modelOkay(int[] dim)
   {
      int t = dim[0] + dim[2] + dim[3] + dim[5]; 
      if(t > 0) 
      {return true;}
      else 
      {System.out.println("Model choice insufficient, choose again"); return false;}
   }


    private void initSlidingSpanComponents() 
    {

        slideSpanCheck = new javax.swing.JCheckBox();
        nbackObsBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,1000);
        nbackObsLabel = new javax.swing.JLabel();
        nbackObsText = new javax.swing.JTextField();
        nforeObsText = new javax.swing.JTextField();
        nforeObsLabel = new javax.swing.JLabel();
        nforeObsBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,1000);

        

        slideSpanCheck.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        slideSpanCheck.setText("Sliding Span Activate");
        slideSpanCheck.setEnabled(false);

        nbackObsBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nbackObsBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        nbackObsBar.setEnabled(false);

        nbackObsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nbackObsLabel.setText("No. Back Points");

        nbackObsText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nbackObsText.setText("0");

        nforeObsText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nforeObsText.setText("0");

        nforeObsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nforeObsLabel.setText("No. Lead Points");

        nforeObsBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nforeObsBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        nforeObsBar.setEnabled(false);
        nforeObsBar.setMinimum(0);
        nforeObsBar.setMaximum(100);
 
        nforeObsBar.setMinimum(0);
        nforeObsBar.setMaximum(100);       

        slideSpanPanel = new JPanel();
        slideSpanPanel.setBorder(BorderFactory.createTitledBorder("Sliding Span Interface"));

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(slideSpanPanel);
        slideSpanPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(slideSpanCheck)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(nbackObsLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nbackObsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 198, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(nbackObsText, javax.swing.GroupLayout.PREFERRED_SIZE, 29, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(nforeObsLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nforeObsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 198, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(nforeObsText, javax.swing.GroupLayout.PREFERRED_SIZE, 29, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(33, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(nbackObsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nforeObsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(slideSpanCheck)
                        .addGap(27, 27, 27)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(nbackObsLabel)
                            .addComponent(nbackObsBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(27, 27, 27)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(nforeObsLabel)
                            .addComponent(nforeObsBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(39, Short.MAX_VALUE))
        );


        slideSpanCheckMenu.addItemListener(new ItemListener() {
         public void itemStateChanged(ItemEvent e)
         {
            boolean sel; 
            e.getItemSelectable();
            if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
            else{sel = true;}

            nforeObsBar.setValue(0); 
            nbackObsBar.setValue(0);
            nforeObsText.setText("0"); 
            nbackObsText.setText("0");
            smc.nbackObs = 0; smc.nObsSpan = 0; smc.nforeObs = 0;

            smc.slidingSpan = sel;
            initializeUseless(nObs-smc.nbackObs-smc.nforeObs);
            smc.computeSARIMAmodel(estimate);            
            
            setDiagnostics(); 
            smc.go(); if(specCntrl){spec.go();}

            TrendModel.setEnabled(!sel); 
            SeasonalModel.setEnabled(!sel); 
            TrendIrregModel.setEnabled(!sel); 
            IrregModel.setEnabled(!sel); 

         }
        }
       ); 
 

 
       AdjustmentListener al = new AdjustmentListener()  {
        public void adjustmentValueChanged(AdjustmentEvent e) {
        
             if(e.getAdjustable() == nbackObsBar)
             {
                if(smc.nObs - nbackObsBar.getValue() - smc.nforeObs > 60)
                {
                   smc.nbackObs = nbackObsBar.getValue();
                   nbackObsText.setText(""+smc.nbackObs);
                   initializeUseless(nObs-smc.nbackObs-smc.nforeObs);
                   smc.changeSlidingSpan(smc.nbackObs, smc.nforeObs);
                   setDiagnostics(); 
                   smc.go(); if(specCntrl){spec.go();}                   
                }
             }
             if(e.getAdjustable() == nforeObsBar)
             {
                if(smc.nObs - smc.nbackObs - nforeObsBar.getValue() > 60)
                {
                   smc.nforeObs = nforeObsBar.getValue();
                   nforeObsText.setText(""+smc.nforeObs);
                   initializeUseless(nObs-smc.nbackObs-smc.nforeObs);
                   smc.changeSlidingSpan(smc.nbackObs, smc.nforeObs);
                   setDiagnostics(); 
                   smc.go(); if(specCntrl){spec.go();}                   
                }
             }
          }
       };

       nbackObsBar.addAdjustmentListener(al);
       nforeObsBar.addAdjustmentListener(al);

    
       
   }


   public void activateSlidingSpan(boolean sel)
   {
     
     nforeObsBar.setValue(0);   nbackObsBar.setValue(0); 
     nforeObsText.setText("0"); nbackObsText.setText("0");
        
     slideSpanCheck.setEnabled(sel); nforeObsBar.setEnabled(sel); nbackObsBar.setEnabled(sel);
     slideSpanCheckMenu.setEnabled(sel);
     
     
   }
















   
   public static void addTab(JTabbedPane tabPane, String text, Component comp){ 
    int tabPlacement = tabPane.getTabPlacement(); 
    switch(tabPlacement){ 
        case JTabbedPane.LEFT: 
        case JTabbedPane.RIGHT: 
            tabPane.addTab(null, new VerticalTextIcon(text, tabPlacement==JTabbedPane.RIGHT), comp); 
            return; 
        default: 
            tabPane.addTab(text, null, comp); 
    } 
   }

  public class ExtensionFilter extends FileFilter {
    private String extensions[];

    private String description;

    public ExtensionFilter(String description, String extension) {
      this(description, new String[] { extension });
    }

    public ExtensionFilter(String description, String extensions[]) {
      this.description = description;
      this.extensions = (String[]) extensions.clone();
    }

    public boolean accept(File file) {
      if (file.isDirectory()) {
        return true;
      }
      int count = extensions.length;
      String path = file.getAbsolutePath();
      for (int i = 0; i < count; i++) {
        String ext = extensions[i];
        if (path.endsWith(ext)
            && (path.charAt(path.length() - ext.length()) == '.')) {
          return true;
        }
      }
      return false;
    }

    public String getDescription() {
      return (description == null ? extensions[0] : description);
    }
  } 

  


 
  
 
  
  
  
  
  
  
  public static double computeDrawdown(double[] ret)
  {
     int i;
     double max = -100; 
     double[] cmax = cummax(ret);
     double[] cmaxx = new double[ret.length];
     
     for(i=0;i<ret.length;i++)
     {
      cmaxx[i] = cmax[i] - ret[i]; 
      if(cmaxx[i] > max) {max = cmaxx[i];}
     }
     
     return max; 
  }
  
  public static double[] cumsum(double[] data, int n)
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
  
  public static double[] cumsum(double[] data)
    {
      int n = data.length;
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

  public static double[] cumsumReal(double[] data)
  {
      int n = data.length;
      double[] cs = new double[n]; double sum; int k;
      double[] account = new double[n];
      account[0] = 100000;
      //sum=Math.abs(data[0]); double min = 1000000;
      sum=data[0]; double min = 1000000;

      
      for(k=1;k<n;k++)
      {
        sum = sum + account[k-1]*data[k]; cs[k] = sum; 
        account[k] = account[k-1] + account[k-1]*data[k];
        if(cs[k] < min) {min = cs[k];}
      }
      return account;  
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
  
  public static double[] mean_std( Double[] data ) 
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
  
    

  class EMDMap extends JComponent
  {


  	private static final long serialVersionUID = 1L;
  	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      
      BasicStroke dashed;
      DecimalFormat df;
      float[] dash1;
      Color myGray;
      Color[] colorArray,colorAM,colorFM,colorIF,colorPh;
      int colorCount;

      int xCanvasPanelCount;  // determines the resolution of the canvas
      int yCanvasPanelCount;

      int xDataCount;
      int yDataCount;
      int dataType;

      double[][] contourData;
      double dataRangeMin;
      double dataRangeMax;


      //---- EMD stuff -------
      int N,n_imfs;
      double[][] fm; //

      public EMDMap(int w, int h, int _N, int co)
      {
       this.height = h; this.width = w; 
       dataRangeMax = -100000000.0; dataRangeMin = 10000000.0;
       setPreferredSize(new Dimension(w, h));

       N=144;
       n_imfs = 4;
       fm = new double[n_imfs][N];
       colorCount = co;
       colorArray = new Color[co];

       colorAM = new Color[co]; 
       colorFM = new Color[co]; 
       colorIF = new Color[co];
       colorPh = new Color[co];
       setupColors(); 
       dataType = 0;
       
      }

      public void setResolution(int PanelCountx, int PanelCounty)
      {xCanvasPanelCount = PanelCountx; yCanvasPanelCount = PanelCounty;}

      public void updateFM(double[][] _fm, int _N, int _n_imfs, int sel)
      {
        int i,m;
        N = _N; n_imfs = _n_imfs;

        if(sel == 0)
        {fm = sqr(_fm,n_imfs,N);}
        else if(sel == 3)
        {fm = icos(_fm,n_imfs,N);} //get phase
        else{fm = _fm;}
    
        dataRangeMax = -100000000.0; dataRangeMin = 10000000.0;
        for(m=0;m<n_imfs;m++)
        {
           for(i=0; i < N; i++)
           {
            if(fm[m][i] < dataRangeMin) dataRangeMin = fm[m][i];
            else if(fm[m][i] > dataRangeMax) dataRangeMax = fm[m][i];
           } 
        }
        computeColorInterp(20);
   
        if(sel == 0) System.arraycopy(colorAM, 0, colorArray, 0, colorAM.length);  
        else if(sel == 1) System.arraycopy(colorFM, 0, colorArray, 0, colorFM.length);
        else if(sel == 2) System.arraycopy(colorIF, 0, colorArray, 0, colorIF.length);
        else System.arraycopy(colorPh, 0, colorArray, 0, colorIF.length);    
        dataType = sel;
        go();
      }
      
      /*--------------------------------------------
        Compute the color map by interpolation between
        each FM mode
      ----------------------------------------------*/
    
      public void computeColorInterp(int hres)
      {
         int m,i,j; 
         int fade = 20;
         double base = 0.0;
         yDataCount = (n_imfs-1)*hres+fade;
         xDataCount = N;
   
         contourData = new double[yDataCount][xDataCount];
         double p1,p2,t; 
        
         for(i=0;i<N;i++)
         {
          
           for(m=0;m<n_imfs-1;m++)
           { 
             p1 = fm[m][i]; p2 = fm[m+1][i]; //interpolate between these two pnts   

             for(j=0;j<hres;j++)
             {
                t = (double)j/(hres-1);
                contourData[m*hres + j][i] =  p1 + t*(p2 - p1);
             }
           }
           //--- do additional last row, taper to zero
           p1 = fm[n_imfs-1][i]; p2 = base; 
           for(j=0;j<fade;j++)
           {
                t = (double)j/(fade-1);
                contourData[(n_imfs-1)*hres + j][i] =  p1 + t*(p2 - p1);
           }
         }
      }
            

      public void setupColors()
      {      
          int i,color;
          float hinc = (float)(0.9/colorCount);
          float h = hinc;
          float saturation = (float)(0.8);
          float intensity  = (float)(0.5);
          int bo = 25;
          float iinc = (float)intensity/bo; float ia = intensity;
          int co4 = colorCount/4;
          float inc4 = 0.5f/co4; 
          float inc3 = 0.75f/co4;

        
            int co2 = colorCount/2; float inc2 = 0.7f/co2; 
            float r = (float)0.0; 
            float b = 0.0f;   
            float g = 0.0f;      


           //-------------Phase Colors---------------
           b=0f; g=0f; r=0f;
           for(i = 0; i < co4; i++)
           {colorPh[i] = new Color(0f, 0f, b+0.20f); b+=inc3;}
           for(i = 0; i < co4; i++)
           {colorPh[co4+i] = new Color(0f, g, b+0.20f ); g+=inc4;}
           for(i = 0; i < co4; i++)
           {colorPh[2*co4+i] = new Color(r, g, b+0.20f ); r+=inc4;}        
           for(i = 0; i < co4; i++)
           {colorPh[3*co4+i] = new Color(r, g, b+0.20f ); r+=inc2;}     
            

          
          //-------------FM Colors------------------------
          for(i = 0; i < colorCount/2-bo; i++)
          {
             color = Color.HSBtoRGB(h, saturation, intensity);
             h +=  hinc;
             colorFM[i] = new Color(color);
          }
          for(i = colorCount/2-bo; i < colorCount/2; i++)
          {
                 ia = ia - iinc; h += hinc;
                 color = Color.HSBtoRGB(h, saturation, ia);
                 colorFM[i] = new Color(color);
          } 
          for(i = colorCount/2; i < colorCount/2+bo; i++)
          {
                 ia = ia + iinc; h += hinc;
                 color = Color.HSBtoRGB(h, saturation, ia);
                 colorFM[i] = new Color(color);
          } 
          for(i = colorCount/2+bo; i < colorCount; i++)
          {
                 color = Color.HSBtoRGB(h, saturation, intensity);
                 h +=  hinc;
                 colorFM[i] = new Color(color);
          }             
           
           //-------------IF Colors---------------
        
           hinc = (float)(1.0/colorCount); h = hinc;
           saturation = (float)(0.7);
           intensity  = (float)(0.6);
           for(i = 0; i < colorCount; i++)
           {
             color = Color.HSBtoRGB(h, saturation, intensity);
             h +=  hinc;
             colorIF[i] = new Color(color);
           }


           //-------------AM Colors---------------
           b=0f; g=0f; r=0f;
           for(i = 0; i < co4; i++)
           {colorAM[i] = new Color(r, 0f, 0f); r+=inc4;}
           for(i = 0; i < co4; i++)
           {colorAM[co4+i] = new Color(r, 0f, 0f ); r+=inc4;}
           for(i = 0; i < co4; i++)
           {colorAM[2*co4+i] = new Color(r-inc4, 0f, b); b+=inc4;}        
           for(i = 0; i < co4; i++)
           {colorAM[3*co4+i] = new Color(r-inc4, 0f, b ); b+=inc4;}        
      }   



     /*------------------------------------------------------------
       Using bilinear interpolation of the input data; Data is
       assumed to be at the nodes of a grid. The rectangle color value
       is taken to be the value at the center of the rectangle of
       a "canvas" grid decomposed into xCanvasPanelCount X yCanvasPanelCount 
       panels 
     ------------------------------------------------------------------------*/

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {

      if(contourData == null) return;
      if(Math.abs(dataRangeMax- dataRangeMin) < 1.0e-10) return;

      Dimension D =  this.getSize(); 
      setResolution(3*D.width/4,1*D.height/2);

      //System.out.println(D.width + " " + D.height);
      int boxWidth  = D.width/xCanvasPanelCount;
      int boxHeight = D.height/yCanvasPanelCount;
      if(boxWidth  <= 1)
      {
       xCanvasPanelCount = D.width/2;
       boxWidth  =  D.width/xCanvasPanelCount;
      }
      if(boxHeight <= 1)
      {
      yCanvasPanelCount = D.height/2;
      boxHeight = D.height/yCanvasPanelCount;
      }
      int xOffset = (D.width  - xCanvasPanelCount*boxWidth)/2;
      int yOffset = (D.height - yCanvasPanelCount*boxHeight)/2;
      int xCanvasLocation;
      int yCanvasLocation;

      double canvasHx; double canvasHy;
      double dataHx;   double dataHy;
      canvasHx = 1.0/((double)(xCanvasPanelCount));
      canvasHy = 1.0/((double)(yCanvasPanelCount));
      dataHx   = 1.0/((double)(xDataCount-1));
      dataHy   = 1.0/((double)(yDataCount-1));

      double xLocation;
      double yLocation;

      int xDataIndex;
      int yDataIndex;

      int i; int j;
      double colorLocation;
      double dataValue;
      int colorIndex;

      double p; double q;
      double v1; double v2; double v3; double v4;
      double phi0; double phi1;

      for(i = 0; i < xCanvasPanelCount; i++)
      {
         xCanvasLocation  = xOffset + i*boxWidth;
         xLocation        = canvasHx/2.0 + ((double)i)*canvasHx;
         xDataIndex = (int)(xLocation/dataHx);
         for(j = 0; j < yCanvasPanelCount; j++)
         {
           yCanvasLocation  = yOffset + j*boxHeight;
           yLocation        = canvasHy/2.0 + ((double)j)*canvasHy;
           yDataIndex = (int)(yLocation/dataHy);
       
           //do bilinear interpolation
    
           p = (xLocation - xDataIndex*dataHx)/dataHx;
           q = (yLocation - yDataIndex*dataHy)/dataHy;
           v1 = contourData[(yDataCount-1) - yDataIndex][xDataIndex]; //*yDataCount];
           v2 = contourData[(yDataCount-1) - yDataIndex][xDataIndex+1];//*yDataCount];
           v3 = contourData[(yDataCount-1)  - (yDataIndex + 1)][(xDataIndex+1)];//*yDataCount];
           v4 = contourData[(yDataCount-1) - (yDataIndex + 1)][xDataIndex];//*yDataCount];
           phi0 = v1 + p*(v2 - v1);
           phi1 = v4 + p*(v3 - v4);
           dataValue = phi0 + q*(phi1 - phi0);

           colorLocation =
             ((colorCount-1.0)*(dataValue - dataRangeMin))/(dataRangeMax - dataRangeMin);
           colorIndex = (int)colorLocation;
           if(colorIndex >= colorCount) colorIndex = colorCount-1;
           if(colorIndex < 0) colorIndex = 0;
           g.setColor(colorArray[colorIndex]);
           g.fillRect(xCanvasLocation,yCanvasLocation,boxWidth,boxHeight);
         }
      }
      g.setColor(Color.black);
    }

    public double[][] sqr(double[][] _fm, int _n_imfs, int _N)
    {
      int m,i;
      double[][] res = new double[_n_imfs][_N];
      for(m=0;m<_n_imfs;m++)
      {
         for(i=0;i<_N;i++) res[m][i] = _fm[m][i]*_fm[m][i];
      }
      return res;  
    }

    public double[][] icos(double[][] _fm, int _n_imfs, int _N)
    {
      int m,i;
      double[][] res = new double[_n_imfs][_N];
      for(m=0;m<_n_imfs;m++)
      {
         for(i=0;i<_N;i++) res[m][i] = Math.acos(_fm[m][i]);
      }
      return res;      
    }  
  }




  class EMDfmcanvas extends JPanel
  {

      /**
  	 * 
  	 */
  	private static final long serialVersionUID = 1L;
  	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      double dataMax, dataMin, dataNorm;
      double dataMaxIF, dataMinIF, dataNormIF;
      BasicStroke dashed;
      DecimalFormat df;
      float[] dash1;
      Color myGray;


      //---- EMD stuff -------
      int N,n_imfs;
      double[][] fm;
      double[][] inst_f;
      boolean[][] plot_fm;
      boolean FM;


      public EMDfmcanvas(int w, int h, int _N, boolean _FM)
      {
       this.height = h; this.width = w; FM = _FM;
       dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
       dash1 = new float[1]; dash1[0] = 10.0f;
       dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                            10.0f, dash1, 0.0f); 
       myGray = new Color(67,71,73);

       N = _N;  n_imfs = 6;
       fm = new double[6][N]; 
       inst_f = new double[6][N]; 
       plot_fm = new boolean[6][4];
       setBackground(Color.BLACK);
       df = new DecimalFormat("##.##"); 
       //setPreferredSize(new Dimension(w, h));
      }

      public void setFM(double[][] _fm, double[][] _instf, int _n, int nmfs)
      {
         int i,m;
         N = _n; n_imfs = nmfs;

         fm = _fm;
         inst_f = _instf;
         dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
         dataMaxIF = -100000000.0; dataMinIF = 10000000.0; dataNormIF = 0.0;
         for(m=0; m < nmfs; m++) 
         {
           for(i=0;i<N;i++)
           {
             if(fm[m][i] < dataMin) dataMin = fm[m][i];
             else if(fm[m][i] > dataMax) dataMax = fm[m][i];
             if(inst_f[m][i] < dataMinIF) dataMinIF = inst_f[m][i];
             else if(inst_f[m][i] > dataMaxIF) dataMaxIF = inst_f[m][i];
           }
         }
         dataNorm = Math.abs(dataMax - dataMin);
         dataNormIF = Math.abs(dataMaxIF - dataMinIF);
         go();
      }

      public void updateFM(int i, boolean w) {plot_fm[i][1] = w; go();}
      public void updateIF(int i, boolean w) {plot_fm[i][2] = w; go();}

      public void setNObs(int n) {N = n;}    
      public void setPlotOptions(boolean[][] _plot_fm) {plot_fm = _plot_fm;}
      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {     
    
        int i,j,colx,blue,green;
        int t0, t1, x0, x1;
        super.paintComponent(g);
        g2d = (Graphics2D)g; 
        colx = (int)(240/n_imfs);      

        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;

        green = 0;
        Color grad = new Color(73,green,253);

        blue = 0;
        Color grad1 = new Color(blue,162,255);


          //Draw dashed lines 
          
          g2d.setStroke(dashed);
          g2d.setPaint(myGray);

          for(i=0; i < 9; i++)
          {
              x0 = (int)(((double)i/(double)8)*(double)height);
              g2d.drawLine(0, x0, width, x0);
          }
          g.drawString((String)df.format(dataMax), 5, 15);
          g.drawString((String)df.format(dataMin), 5, height - 5);
          int nobsP = (int)Math.floor((double)N/12); 
          for(i=1; i <= nobsP; i++) 
          {
            t0 = (int)(((double)(i*12)/N)*(double)width);
  	  g2d.drawLine(t0, 0, t0, height-20);
          }


        g2d.setStroke(new BasicStroke(1.0f));
        

        for(i=0; i < n_imfs; i++)
        {

          if(plot_fm[i][1] && FM)
          {
            g2d.setPaint(grad); 
            for(j = 0; j < N-1; j++)
  	  {
  	    t0 = (int)(((double)j/(double)N)*(double)width);
  	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
  	    x0 = (int)(((fm[i][j] - dataMin)/dataNorm)*(double)height);
  	    x1 = (int)(((fm[i][j+1] - dataMin)/dataNorm)*(double)height);
  	    g2d.drawLine(t0, height - x0, t1, height - x1);
  	  }
            green = green + colx;
            grad = new Color(73,green,253);                      
          }
          else if(plot_fm[i][2] && !FM)
          {
            g2d.setPaint(grad1); 
            for(j = 0; j < N-1; j++)
  	  {
  	    t0 = (int)(((double)j/(double)N)*(double)width);
  	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
  	    x0 = (int)(((inst_f[i][j] - dataMinIF)/(dataNormIF*1.10))*(double)height);
  	    x1 = (int)(((inst_f[i][j+1] - dataMinIF)/(dataNormIF*1.10))*(double)height);
  	    g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
  	  }
            blue = blue + colx;
            grad1 = new Color(blue,162,255);             
          }
        }
      }
  }



 
  class EMDMapKey extends JComponent
  {


	private static final long serialVersionUID = 1L;
	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      int co;
      Color[] colorArray;

      public EMDMapKey(int w, int h, int _co)
      { 
        this.height = h; this.width = w; 
        setPreferredSize(new Dimension(w, h));
        co = _co;
        colorArray = new Color[co];           
      }

      public void updateColor(Color[] c)
      {System.arraycopy(c, 0, colorArray, 0, c.length); go();}

      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {
          int i;
          Dimension D =  this.getSize(); 
          int boxHeight = D.height/co;

          for(i=0; i<co; i++)
          {
            g.setColor(colorArray[i]);
            g.fillRect(0,D.height-i,D.width,boxHeight);  
          }
      }

  }



  class EMDPhasecanvas extends JPanel
  {


  	private static final long serialVersionUID = 1L;
  	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      double dataMax, dataMin, dataNorm;
      BasicStroke dashed;
      DecimalFormat df;
      float[] dash1;
      Color myGray;

      //---- EMD stuff -------
      int N,n_imfs;
      double[][] phase;
      boolean[][] plot_fm;

      public EMDPhasecanvas(int w, int h, int _N)
      {
       this.height = h; this.width = w;
       dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
       dash1 = new float[1]; dash1[0] = 10.0f;
       dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                            10.0f, dash1, 0.0f); 
       myGray = new Color(67,71,73);

       N = _N;  n_imfs = 6;
       phase = new double[6][N]; 
       plot_fm = new boolean[6][4];
       setBackground(Color.BLACK);
       df = new DecimalFormat("##.##"); 
       //setPreferredSize(new Dimension(w, h));
      }

      public void setPhase(double[][] _fm, int _n, int nmfs)
      {
         int i,m;
         N = _n; n_imfs = nmfs;

         phase = _fm;
         
         dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
         for(m=0; m < nmfs; m++) 
         {
           for(i=0;i<N;i++)
           {
             if(phase[m][i] < dataMin) dataMin = phase[m][i];
             else if(phase[m][i] > dataMax) dataMax = phase[m][i];
           }
         }
         dataNorm = Math.abs(dataMax - dataMin);
         go();
      }

      public void updatePhase(int i, boolean w) {plot_fm[i][3] = w; go();}
      

      public void setNObs(int n) {N = n;}    
      public void setPlotOptions(boolean[][] _plot_fm) {plot_fm = _plot_fm;}
      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {     
    
        int i,j,colx,blue,green;
        int t0, t1, x0, x1;
        super.paintComponent(g);
        g2d = (Graphics2D)g; 
        colx = (int)(240/n_imfs);      

        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;

        green = 0;
        Color grad = new Color(73,green,253);

        blue = 0;
        new Color(blue,162,255);


          //Draw dashed lines 
          
          g2d.setStroke(dashed);
          g2d.setPaint(myGray);

          for(i=0; i < 9; i++)
          {
              x0 = (int)(((double)i/(double)8)*(double)height);
              g2d.drawLine(0, x0, width, x0);
          }
          g.drawString((String)df.format(dataMax), 5, 15);
          g.drawString((String)df.format(dataMin), 5, height - 5);
          int nobsP = (int)Math.floor((double)N/12); 
          for(i=1; i <= nobsP; i++) 
          {
            t0 = (int)(((double)(i*12)/N)*(double)width);
  	  g2d.drawLine(t0, 0, t0, height-20);
          }


        g2d.setStroke(new BasicStroke(1.0f));
        
        for(i=0; i < n_imfs; i++)
        {
          if(plot_fm[i][3])
          {
            g2d.setPaint(grad); 
            for(j = 0; j < N-1; j++)
  	  {
              //System.out.println(phase[i][j]);
  	    t0 = (int)(((double)j/(double)N)*(double)width);
  	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
  	    x0 = (int)(((phase[i][j] - dataMin)/dataNorm)*(double)height);
  	    x1 = (int)(((phase[i][j+1] - dataMin)/dataNorm)*(double)height);
  	    g2d.drawLine(t0, height - x0, t1, height - x1);
  	  }
            green = green + colx;
            grad = new Color(73,green,253);                      
          }
        }

      }
  }


  class EMD3dScale extends JPanel
  {

      /**
  	 * 
  	 */
  	private static final long serialVersionUID = 1L;
  	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      double dataMax, dataMin; 
      Color myGray;
      Font font;
      
      public EMD3dScale(int w, int h)
      { 
        this.height = h; this.width = w; 
        setPreferredSize(new Dimension(w, h));
        setBackground(Color.BLACK);        
        myGray = new Color(67,71,73);  
        dataMax = -10000000.0; dataMin = 100000000.0;
        font = new Font("serif", Font.PLAIN, 9);
      }

      public void updateData(double _max, double _min) 
      {dataMax = _max; dataMin = _min; go();}
      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {
          super.paintComponent(g);
          g2d = (Graphics2D)g; 
          Dimension D =  this.getSize(); 
   
          double mid = (dataMax - dataMin)/2.0;
          g2d.setPaint(myGray);
          g2d.setFont(font);
   
          g.drawString(""+dataMin, 0 , D.height-3);
          g.drawString(""+mid,0 , D.height/2);
          g.drawString(""+dataMax, 0, 3);
      }
  }



  class EMDIMFcanvas extends JPanel
  {

   
  	private static final long serialVersionUID = 1L;
  	//------ Other stuff
      Graphics2D g2d;
      int height, width;
      double dataMax, dataMin, dataNorm;
      BasicStroke dashed;
      float[] dash1;
      Color myGray;
      DecimalFormat df;
      //---- EMD stuff -------
      int N,n_imfs;
      double[][] imfs;
      boolean[] plot_imfs;

      public EMDIMFcanvas(int w, int h, int _N)
      {
       int i;
       this.height = h; this.width = w; 
       dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
       dash1 = new float[1]; dash1[0] = 10.0f;
       dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                            10.0f, dash1, 0.0f); 
       myGray = new Color(67,71,73);
      
       N = _N; n_imfs = 6;
       imfs = new double[6][N]; plot_imfs = new boolean[6];
       for(i=0; i < 6; i++) plot_imfs[i] = false;

       setBackground(Color.BLACK); df = new DecimalFormat("##.##"); 
       //setPreferredSize(new Dimension(w, h));
      }

      public void setNObs(int n) {N = n;}

      public void setIMFs(double[][] _am, double[][] _fm, int _n, int nmfs)
      {
         int i,m;
         N = _n; n_imfs = nmfs;
         imfs = new double[n_imfs][N];
    
         for(m=0; m < nmfs; m++) 
         {
           for(i=0;i<N;i++)
           {imfs[m][i] = _am[m][i]*_fm[m][i];}
         }       
        
         dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
         for(m=0; m < nmfs; m++) 
         {
           for(i=0;i<N;i++)
           {
             if(imfs[m][i] < dataMin) dataMin = imfs[m][i];
             else if(imfs[m][i] > dataMax) dataMax = imfs[m][i];
           }
         }
         dataNorm = Math.abs(dataMax - dataMin);
         go();
      }
          
      public void updateIMFs(int i, boolean w) {plot_imfs[i] = w; go();}
    
      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {     
        int i,j,colx,red;
        super.paintComponent(g);
        int t0, t1, x0, x1;
        g2d = (Graphics2D)g;

        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;  

          //Draw dashed lines 
          g2d.setStroke(dashed);
          g2d.setPaint(myGray);

          for(i=0; i < 8; i++)
          {
              x0 = (int)(((double)i/(double)8)*(double)height);
              g2d.drawLine(0, x0, width, x0);
          }
          g.drawString((String)df.format(dataMax), 5, 15);
          g.drawString((String)df.format(dataMin), 5, height - 5);
          int nobsP = (int)Math.floor((double)N/12); 
          for(i=1; i <= nobsP; i++) 
          {
            t0 = (int)(((double)(i*12)/N)*(double)width);
  	  g2d.drawLine(t0, 0, t0, height-20);
          }
   
        g2d.setStroke(new BasicStroke(1.0f));
        colx = (int)(240/n_imfs);      
        red = 253;
        Color grad = new Color(165,red,211);
        g2d.setPaint(grad);

      
        for(i=0; i < n_imfs; i++)
        {
          if(plot_imfs[i])
          {
            for(j = 0; j < N-1; j++)
  	  {
  	    t0 = (int)(((double)j/(double)N)*(double)width);
  	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
  	    x0 = (int)(((imfs[i][j] - dataMin)/(dataNorm*1.20))*(double)height);
  	    x1 = (int)(((imfs[i][j+1] - dataMin)/(dataNorm*1.20))*(double)height);
  	    g2d.drawLine(t0, height - x0, t1, height - x1);
  	  }
            red = red - colx;
            grad = new Color(165,red,211);
            g2d.setPaint(grad);    
          }
        }

      }
  }


  class EMDMapScale extends JPanel
  {



  	private static final long serialVersionUID = 1L;
  	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      int N; 
      Color myGray;
      Font font;
      
      public EMDMapScale(int w, int h)
      { 
        this.height = h; this.width = w; 
        setPreferredSize(new Dimension(w, h));
        setBackground(Color.BLACK);        
        myGray = new Color(67,71,73); N = 144; 
        font = new Font("serif", Font.PLAIN, 10);
      }

      public void updateN(int _N) {N = _N; go();}

      public void go() {repaint();}

      public void paintComponent(Graphics g)
      {
          int i; int t0; int p;
          super.paintComponent(g);
          g2d = (Graphics2D)g; 
          Dimension D =  this.getSize(); 

          //Draw dashed lines 
          
          g2d.setPaint(myGray); p = 12; g2d.setFont(font);
          int nobsP = (int)Math.floor((double)N/12); 
          for(i=1; i <= nobsP; i++) 
          {
            t0 = (int)(((double)(i*12)/N)*(double)D.width);
  	  g2d.drawLine(t0, 0, t0, 5);
            g.drawString((String)"" + p, t0-3, height - 3); p = p + 12;
          }

      }

  }


  class EMDFreqScale extends JPanel
  {

      /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//------ Plot stuff
      Graphics2D g2d;
      int height, width;
      int N,n_imfs; 
      Color myGray;
      Font font;
      
      public EMDFreqScale(int w, int h)
      { 
        this.height = h; this.width = w; 
        setPreferredSize(new Dimension(w, h));
        setBackground(Color.BLACK);        
        myGray = new Color(67,71,73); N = 144; n_imfs = 3;
        font = new Font("serif", Font.PLAIN, 9);
      }

      public void updateN_IMFS(int _N) {n_imfs = _N; go();}
      public void go() {repaint();}
      public void paintComponent(Graphics g)
      {
          int i; int p;
          super.paintComponent(g);
          g2d = (Graphics2D)g; 
          Dimension D =  this.getSize(); 

          g2d.setPaint(myGray);
          g2d.setFont(font);
   
          int mDist = D.height/n_imfs; p = mDist/2;
          for(i=1; i < n_imfs; i++) 
          {g.drawString("2\u03C0 N/" + 2*i, 0, p); p+=mDist;}

          

      }

  }


}
