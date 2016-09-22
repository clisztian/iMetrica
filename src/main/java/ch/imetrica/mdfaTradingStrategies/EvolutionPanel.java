package ch.imetrica.mdfaTradingStrategies;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.net.Socket;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.joda.time.DateTime;
import org.joda.time.Days;

import ch.imetrica.mdfa.IMDFAPanel;
import ch.imetrica.mdfaTradingStrategies.HistoricalData;
import ch.imetrica.mdfaTradingStrategies.MDFAStrategyEvolution;

public class EvolutionPanel extends JPanel
{
   
   
   /**
	 * 
	 */
   private static final long serialVersionUID = 1L;
   int n_obs, n_assets; 
   int i1,i2,M,nlag; 
   int meta_filter = 0;
   int port_normalize = 0;
   double hybrid_weight=0; 
   double hybrid_weight_diff = 0;
   int b0trend = 0; int sig_diff = 0;
   double time_shift=0;
   double smooth, decay, decay2, cross; 
   double lambda, expweight; 
   boolean filter_file_set = false; 
   boolean params_set = false;
   boolean historicalData_set = false;
   String[] portfolio; 
   boolean portfolio_file_set = false;
   boolean strategy_computed = false;
   String[] data_files;
   String[] sp_portfolio; 
   boolean sp_files_set = false;
   boolean longOnly = false; 
   boolean equal_dist = false;
   boolean all_returns_computed = false;
   boolean black_list_computed = false;
   int l_filter=0;
   int k_file = 0; 
   int minute_df = 15;
   //Task task; 
   int n_filters =0; 
   boolean saving_perf = false;
   double avg_rank,rank_coeff;
   boolean realrets; 
   JCheckBox realretsCheck; 
   double tradingCost = 0;
   JButton exportResult; 
   String final_time = "16:00";
   int sig_inv;
   JButton exportStrategy;
   ArrayList<METAFilterParameters> metaToFile; 
   JScrollBar timeScrollBar;
   int L, lag;
   double cutoff0, cutoff; 
   File filter_file; 
   JFrame frame;
   int n_threads = 10;
   int stop_loss_int;
   double stop_loss;
   int n_saved_perf = 0;
   JPanel dataConfigPanel;
   JPanel strategyConfigPanel;
   JPanel metaFilterConfigPanel;
   DecimalFormat df,df2,df3,df4; 
   File filterParamFile; 
   JDialog dataConfigDialog, strategyConfigDialog, metaFilterConfigDialog, stratSelectDialog;
   JButton computeStratButton, dataConfigButton, metaFilterConfigButton, strategyConfigButton, stratSelectButton;
   JCheckBox[] jCheckBox; 
   JScrollBar ncoresBar;
   JLabel ncoresLabel;
   JTextField ncoresText;
   JProgressBar progressBarFiles, progressBar;
   JPanel theScreenPanel;
   ArrayList<String> expvar;
   double sharpe_ratio,max_drawdown,min_rank;   
   int n_sections;
   boolean exp_var = false;
   boolean pairs_trading = false;
   JInvestment[] performances;
   ArrayList<JInvestment[]> portfolio_invest;
   JFileChooser fc; 
   //Data panel 
   JPanel assetFilePanel,downloadUniversePanel;
   JRadioButton assetUnivButton,downloadUnivCheck;
   JButton connectIQButton,downloadIQButton;
   JScrollBar fileNameScroll,fileNameScroll1;
   JTextField fileNameText,fileNameText1;
   JLabel jLabel1,noYearsLabel,startTimeLabel;
   JButton loadMetaButton,loadMetaButton1,loadMetaButton2;
   JTextField metafileAssetText,metafileAssetText1,metafileAssetText2;
   @SuppressWarnings("rawtypes")
   JComboBox noYearsCombo;
   JCheckBox forexData,data24,data30,data60,futuresData,data5,zeroButton;
   JTextField tradeTimeText;
   String tradeTime24;
   JComboBox<String> asset_plots;
   JButton add_plots;
   JButton clear_asset_plots;
   String[] symbols_string;
   JButton save_stratButton;
   JButton loadPerformances;  


   JLabel tradingCostsLabel;
   JScrollBar tradingCostBar;
   JTextField tradingCostText;
   String filter_stock_choice; 
   //Strategy panel 
   JButton addStrategyButt,deleteStrategyButt,saveFiltersButton;
   JLabel hybridLabel;
   JScrollBar hybridScroll;
   JTextField hybridText;
   JCheckBox i1stratCheck,i2stratCheck;
   @SuppressWarnings("rawtypes")
   JComboBox insampStartCombo;
   @SuppressWarnings("rawtypes")
   JComboBox historicalDataStart;
   JCheckBox mdfaPanelParamsBox,multivarBox;
   JPanel mdfaSettingsPanel;
   JScrollBar nobsBar,stopLossBar;
   JLabel nobsLabel;
   JTextField nobsText;
   JButton paramFileButton;
   JLabel paramFileLabel;
   JTextField paramFileText;
   JLabel stopLossLabel;
   JTextField stopLossText;
   JLabel timeShiftLabel;
   JScrollBar timeShiftScroll;
   JTextField timeShiftText;
   @SuppressWarnings("rawtypes")
   JComboBox tradeRuleCombo;
   JLabel tradeRuleLabel;
   @SuppressWarnings("rawtypes")
   JComboBox tradingEndCombo;
   JLabel tradingEndLabel,tradingInSampLabel;
   @SuppressWarnings("rawtypes")
   JComboBox tradingStartCombo;
   JLabel tradingStartLabel;   
   String curDir;   
   JScrollPane imetrica_scrollPane;                   
    JCheckBox equalDistCheck;
    JLabel lagDaysLabel;
    JScrollBar lagDaysScroll;
    JTextField lagDaysText;
    JCheckBox longOnlyBox;
    @SuppressWarnings("rawtypes")
	JComboBox metaFilterCombo;
    JLabel metaFilterLabel;
    JPanel metaFilterPanel;
    @SuppressWarnings("rawtypes")
	JComboBox portWeightsCombo;
    JLabel portWeightsLabel;
    JLabel sizePortfolioLabel;
    JScrollBar sizePortfolioScroll;
    JTextField sizePortfolioText;      
    JTabbedPane rt_plotPanelPane;
    EvolutionCanvas mdfaEvolutionCanvas;
    StrategyCanvas stratCanvas;     
    JCheckBox[] stratCheck;
    JCheckBox stratLongCheck,stratShortCheck,stratBothCheck;
    JButton saveStrategyButton;
    JCheckBox aggregatePerformance;
    JRadioButton uniformWeightsCheck,maxSharpeWeightsCheck;
    JButton uploadFiltersButton;
    JCheckBox blackListCheck;
    int blength = 0;
    int nblack = 0;
    int[][] black_list;    
    JComboBox<String> filter_stock; 
     
     
    IMDFAPanel mdfa;
    boolean turnOnH0 = false;
    JPanel stratSelectPanel;
    ArrayList<String> stock_names;
    ArrayList<double[]> stock_returns_all;    
    BufferedWriter[] s_outs;
    BufferedReader[] s_ins;
      
    //--- the computational core   
   
    private HistoricalData data; 
    MDFAStrategyEvolution evolution;  //contains the list of 
    StrategyFilter strategy_suite; 
    ArrayList<StrategyParameters> strategies; 
    private Component parent;
    
   
   public EvolutionPanel(JFrame _f, Component _p)
   {
   
     //--- set-up auxillary panels-----
     this.frame = _f; 
     parent = _p;
     i1 = 0;
     //sthis.setPreferredSize(new Dimension(1000,700));
     curDir=System.getProperty("user.dir") + File.separator;
     fc = new JFileChooser(curDir);
     fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);     
     df = new DecimalFormat("##.##"); 
     df2 = new DecimalFormat("##.###");
     df3 = new DecimalFormat("##.####");   
     df4 = new DecimalFormat("##.#####");     
     filter_stock_choice = new String("all");
     M = 5; nlag = 10; n_saved_perf = 0;
     setupDataPanel();
     setupStrategyPanel();
     setupMetaFilterPanel();
    
     metaToFile = new ArrayList<METAFilterParameters>();
     black_list = new int[0][0];
     portfolio_invest = new ArrayList<JInvestment[]>();
     strategies = new ArrayList<StrategyParameters>();
     //setup Evolution Components
     initComponents();
   
   }
   

 
   
    public void initComponents() 
    {
    
      int width = 1000; int height = 700; int i; GroupLayout paramLayout;
      stratCanvas = new StrategyCanvas(width, height);
      mdfaEvolutionCanvas = new EvolutionCanvas(width, height, 400, 2);  
        
      imetrica_scrollPane = new JScrollPane();
      imetrica_scrollPane.setPreferredSize(new Dimension(width, height)); 
      imetrica_scrollPane.setVerticalScrollBarPolicy(
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
      imetrica_scrollPane.setViewportView(stratCanvas);
      
         
      
        
        
      rt_plotPanelPane = new JTabbedPane(JTabbedPane.TOP);  
      rt_plotPanelPane.add(mdfaEvolutionCanvas,"Evolution Performance");
      rt_plotPanelPane.add(stratCanvas, "Strategies");
      //rt_plotPanelPane.setPreferredSize(new Dimension(800,500));  
        
      JPanel controlsPane = new JPanel();
      
       paramLayout = new GroupLayout(controlsPane); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
          .addComponent(rt_plotPanelPane));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()  
          .addComponent(rt_plotPanelPane));
       controlsPane.setLayout(paramLayout);  
      //controlsPane.add(rt_plotPanelPane);          
      //controlsPane.setPreferredSize(new Dimension(1000,700));    
        
        dataConfigButton = new JButton();
        strategyConfigButton = new JButton();
        metaFilterConfigButton = new JButton();
        computeStratButton = new JButton();
        stratSelectButton = new JButton();
        ncoresText = new JTextField();
        ncoresBar = new JScrollBar();
        ncoresLabel = new JLabel();
        tradingCostsLabel = new JLabel("Costs");
        tradingCostBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,30);
        tradingCost = 0.0;
        progressBar = new JProgressBar(0,100);
        progressBarFiles = new JProgressBar(0,100);
        progressBar.setStringPainted(true);
        progressBarFiles.setStringPainted(true);
        
        tradingCostText = new JTextField(5);
        
        jCheckBox = new JCheckBox[42];
        for(i = 0; i < 42; i++) {jCheckBox[i] = new JCheckBox(""); }
        
        JLabel asset_plotsLabel = new JLabel("   Plot:");
        asset_plots = new JComboBox<String>();
        add_plots = new JButton("Add Plot");
        clear_asset_plots = new JButton("Clear Plots");
        add_plots.setEnabled(false); clear_asset_plots.setEnabled(false);
        asset_plots.setEnabled(false);
      
        dataConfigDialog = new JDialog(frame,true);
        strategyConfigDialog = new JDialog(frame,true);
        metaFilterConfigDialog = new JDialog(frame,true);
        stratSelectDialog = new JDialog(frame,true);

        ncoresBar.setMaximum(30);
        ncoresBar.setMinimum(1);
        ncoresBar.setOrientation(JScrollBar.HORIZONTAL);
        ncoresBar.setValue(10);


        //---- setup the dialog panels
        dataConfigDialog.getContentPane().add(dataConfigPanel); dataConfigDialog.pack(); 
        dataConfigDialog.setLocationRelativeTo(frame); 
        dataConfigDialog.setModal(false); dataConfigDialog.setVisible(false);        
        
        strategyConfigDialog.getContentPane().add(strategyConfigPanel); strategyConfigDialog.pack(); 
        strategyConfigDialog.setLocationRelativeTo(frame); 
        strategyConfigDialog.setModal(false); 
        strategyConfigDialog.setVisible(false);
        
        metaFilterConfigDialog.getContentPane().add(metaFilterConfigPanel); metaFilterConfigDialog.pack(); 
        metaFilterConfigDialog.setLocationRelativeTo(frame); 
        metaFilterConfigDialog.setModal(false); metaFilterConfigDialog.setVisible(false);        
        
        stratSelectDialog.getContentPane().add(stratSelectPanel); stratSelectDialog.pack(); 
        stratSelectDialog.setLocationRelativeTo(frame); 
        stratSelectDialog.setModal(false); stratSelectDialog.setVisible(false);        
                
               
        uniformWeightsCheck = new JRadioButton("Uniform"); uniformWeightsCheck.setSelected(true); uniformWeightsCheck.setEnabled(false); 
        maxSharpeWeightsCheck = new JRadioButton("MaxSharpe"); maxSharpeWeightsCheck.setSelected(false); maxSharpeWeightsCheck.setEnabled(false);
        aggregatePerformance = new JCheckBox("Aggregate"); aggregatePerformance.setSelected(false); aggregatePerformance.setEnabled(false);
        blackListCheck = new JCheckBox("BlackList"); blackListCheck.setSelected(false); blackListCheck.setEnabled(false);
        realretsCheck = new JCheckBox("real"); realretsCheck.setSelected(false); realretsCheck.setEnabled(false);
        
        ButtonGroup weightGroup = new ButtonGroup();
        weightGroup.add(uniformWeightsCheck);
        weightGroup.add(maxSharpeWeightsCheck);
                

//         GroupLayout theScreenPanelLayout = new GroupLayout(controlsPane);
//         controlsPane.setLayout(theScreenPanelLayout);
//         theScreenPanelLayout.setHorizontalGroup(
//             theScreenPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
//             .addGap(0, 0, Short.MAX_VALUE)
//         );
//         theScreenPanelLayout.setVerticalGroup(
//             theScreenPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
//             .addGap(0, 609, Short.MAX_VALUE)
//         );

        dataConfigButton.setText("Data Configuration");

        strategyConfigButton.setText("Strategy Configuration");

        metaFilterConfigButton.setText("Meta-Filter Configuration");

        computeStratButton.setText("Compute Strategy");

        stratSelectButton.setText("Select Strategies");
        
        ncoresText = new JTextField(3);
        ncoresText.setText("10");

        ncoresBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        ncoresBar.setMaximum(20);
        ncoresBar.setMinimum(1);
        ncoresBar.setOrientation(JScrollBar.HORIZONTAL);
        ncoresBar.setValue(10);

        ncoresLabel.setText("#Threads");

        
        
/*       JPanel coresCon = new JPanel(); 
       ncoresBar.setPreferredSize(new Dimension(80, 15));        
        paramLayout = new GroupLayout(coresCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
        .addComponent(ncoresLabel).addComponent(ncoresBar).addComponent(ncoresText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(ncoresLabel).addComponent(ncoresBar).addComponent(ncoresText)));
        coresCon.setLayout(paramLayout);   */    
        
 
       JPanel costCon = new JPanel(); 
       ncoresBar.setPreferredSize(new Dimension(80, 15));        
        paramLayout = new GroupLayout(costCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
        .addComponent(tradingCostsLabel).addComponent(tradingCostBar).addComponent(tradingCostText).addComponent(ncoresLabel).addComponent(ncoresBar).addComponent(ncoresText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(tradingCostsLabel).addComponent(tradingCostBar).addComponent(tradingCostText).addComponent(ncoresLabel).addComponent(ncoresBar).addComponent(ncoresText)));
        costCon.setLayout(paramLayout);   
        
        
        
        JPanel buttonCon = new JPanel();
        paramLayout = new GroupLayout(buttonCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(dataConfigButton).addComponent(strategyConfigButton).addComponent(metaFilterConfigButton).addComponent(computeStratButton).addComponent(stratSelectButton).addComponent(costCon));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(dataConfigButton).addComponent(strategyConfigButton).addComponent(metaFilterConfigButton).addComponent(computeStratButton).addComponent(stratSelectButton).addComponent(costCon)));
        buttonCon.setLayout(paramLayout);
        
        JPanel plotChecksCon = new JPanel();
        paramLayout = new GroupLayout(plotChecksCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup() 
                    .addComponent(jCheckBox[0]).addComponent(jCheckBox[1]).addComponent(jCheckBox[2]).addComponent(jCheckBox[3])
                    .addComponent(jCheckBox[4]).addComponent(jCheckBox[5]).addComponent(jCheckBox[6])
                    .addComponent(jCheckBox[7]).addComponent(jCheckBox[8]).addComponent(jCheckBox[9])
                    .addComponent(jCheckBox[10]).addComponent(jCheckBox[11])
                    .addComponent(jCheckBox[12]).addComponent(jCheckBox[13]).addComponent(jCheckBox[14])
                    .addComponent(jCheckBox[15]).addComponent(jCheckBox[16]).addComponent(aggregatePerformance)
                    .addComponent(uniformWeightsCheck).addComponent(maxSharpeWeightsCheck).addComponent(blackListCheck).addComponent(realretsCheck)
                    .addComponent(asset_plotsLabel).addComponent(asset_plots).addComponent(add_plots)
                    .addComponent(clear_asset_plots));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(jCheckBox[0]).addComponent(jCheckBox[1]).addComponent(jCheckBox[2]).addComponent(jCheckBox[3])
                    .addComponent(jCheckBox[4]).addComponent(jCheckBox[5]).addComponent(jCheckBox[6])
                    .addComponent(jCheckBox[7]).addComponent(jCheckBox[8]).addComponent(jCheckBox[9])
                    .addComponent(jCheckBox[10]).addComponent(jCheckBox[11])
                    .addComponent(jCheckBox[12]).addComponent(jCheckBox[13]).addComponent(jCheckBox[14])
                    .addComponent(jCheckBox[15]).addComponent(jCheckBox[16]).addComponent(aggregatePerformance)
                    .addComponent(uniformWeightsCheck).addComponent(maxSharpeWeightsCheck).addComponent(blackListCheck).addComponent(realretsCheck)
                    .addComponent(asset_plotsLabel).addComponent(asset_plots).addComponent(add_plots)
                    .addComponent(clear_asset_plots)));
        plotChecksCon.setLayout(paramLayout); 
        
        
        JPanel progressCon = new JPanel();
        paramLayout = new GroupLayout(progressCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup() 
          .addComponent(progressBar).addComponent(progressBarFiles));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(progressBar).addComponent(progressBarFiles)));
        progressCon.setLayout(paramLayout); 
          
          
        paramLayout = new GroupLayout(this);  
        this.setLayout(paramLayout);
         paramLayout.setAutoCreateGaps(false);
         paramLayout.setAutoCreateContainerGaps(false);        
         paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(controlsPane)
              .addComponent(plotChecksCon)
              .addComponent(buttonCon)
              .addComponent(progressCon)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
              .addComponent(controlsPane)
              .addComponent(plotChecksCon)
              .addComponent(buttonCon)
              .addComponent(progressCon));
        

    
      ItemListener MyItemListener = new ItemListener() 
      {
        public void itemStateChanged(ItemEvent e)
        {         
         int i; boolean sel; Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
 
         for(i=0;i<16;i++)
         {     
           if(source == jCheckBox[i]) 
           {mdfaEvolutionCanvas.setplot(i,sel);}
         }
         if(source == aggregatePerformance)
         {if(n_saved_perf > 0) {plotPerformances(); mdfaEvolutionCanvas.plotTarget(sel);}}
         else if(source == uniformWeightsCheck) {if(aggregatePerformance.isSelected() && n_saved_perf > 0) {plotPerformances(); mdfaEvolutionCanvas.plotTarget(sel);}}
         else if(source == maxSharpeWeightsCheck) {if(aggregatePerformance.isSelected() && n_saved_perf > 0) {plotPerformances(); mdfaEvolutionCanvas.plotTarget(sel);}}
         else if(source == blackListCheck) {if(black_list_computed) recomputeStrategy();}
         else if(source == realretsCheck) { realrets = sel; plotPerformances();}
       }  
      };      
    
      for(i=0;i<42;i++) {jCheckBox[i].setSelected(false); jCheckBox[i].setEnabled(false); jCheckBox[i].addItemListener(MyItemListener);}
      uniformWeightsCheck.addItemListener(MyItemListener);
      aggregatePerformance.addItemListener(MyItemListener);
      maxSharpeWeightsCheck.addItemListener(MyItemListener);
      blackListCheck.addItemListener(MyItemListener);
      realretsCheck.addItemListener(MyItemListener);
    
      ActionListener buttonActionListener = new ActionListener() 
      {
        public void actionPerformed(ActionEvent event)
        {
          if(event.getSource() == dataConfigButton)
          {dataConfigDialog.setModal(false); dataConfigDialog.setVisible(true);}
          if(event.getSource() == strategyConfigButton)
          {strategyConfigDialog.setModal(false); strategyConfigDialog.setVisible(true);}    
          if(event.getSource() == metaFilterConfigButton)
          {metaFilterConfigDialog.setModal(false); metaFilterConfigDialog.setVisible(true);} 
          if(event.getSource() == stratSelectButton)
          {stratSelectDialog.setModal(false); stratSelectDialog.setVisible(true);}           
          if(event.getSource() == computeStratButton)
          { 
           try{computeStrategy();} catch(InterruptedException ie){} catch(ExecutionException ee){}
          }
          if(event.getSource() == add_plots)
          {
            if(asset_plots.isEnabled()) addAssetPlot(asset_plots.getSelectedIndex()); 
          }          
          if(event.getSource() == clear_asset_plots)
          {mdfaEvolutionCanvas.clearAssetPlots();}
          if(event.getSource() == exportStrategy)
          {buildStrategyToFile();} 
          if(event.getSource() == exportResult)
          {printAggregateStrategy();}
        }
      };  
      
      exportResult.addActionListener(buttonActionListener);
      exportStrategy.addActionListener(buttonActionListener);
      add_plots.addActionListener(buttonActionListener);
      clear_asset_plots.addActionListener(buttonActionListener);
      dataConfigButton.addActionListener(buttonActionListener);
      strategyConfigButton.addActionListener(buttonActionListener);
      metaFilterConfigButton.addActionListener(buttonActionListener);
      computeStratButton.addActionListener(buttonActionListener);
      stratSelectButton.addActionListener(buttonActionListener); 
       
      ncoresBar.addAdjustmentListener(new AdjustmentListener() {
           public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               numberOfThreads(((JScrollBar)e.getSource()).getValue());
               ncoresText.setText(""+n_threads);
            }
       });    
    
      tradingCostBar.addAdjustmentListener(new AdjustmentListener() {
           public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               tradingCost = ((JScrollBar)e.getSource()).getValue()*.00001;
               tradingCostText.setText(""+df4.format(tradingCost));
               System.out.println("avg trading costs = " + tradingCost); 
            }
       });        
    
    
    }                  
    
  /*--------------------------------------------------------------------
  
    This function is used when the filtering strategies are complete and ready 
    to be used in live trading 
    
    -- Builds three separatate functions, one for the data, one for setting the filtering strategies
       and one for setting the meta-filtering strategies
    
    
   --------------------------------------------------------------------*/ 
  public void buildStrategyToFile()
  {
   if(strategies.size() > 0 && metaToFile.size() > 0)
   {
      
     FileFilter ft = new FileNameExtensionFilter("Strategy files", "strat");
     fc.addChoosableFileFilter(ft);       
     int returnVal = fc.showSaveDialog(frame);
      
     if(returnVal == javax.swing.JFileChooser.APPROVE_OPTION)
     {
       File saved_file = fc.getSelectedFile();
       String file_name = saved_file.toString();

       try
       {
       PrintWriter out = new PrintWriter(new FileWriter(file_name));   
   
   
       //---- set the asset universe first 
       String portFile = new String("public void setAssetUniverse()\n{\n");
       portFile = portFile + "String[] universe = {";
       int count = 0;
       for(int i = 0; i < stock_names.size()-1;i++)
       {
         portFile = portFile + "\"" + stock_names.get(i) + "\",";
         count++;
      
         if(count == 10)
         {portFile = portFile + "\n"; count = 0;}
       }
       portFile = portFile + "\"" + stock_names.get(stock_names.size()-1) + "\"};\n";
    
       portFile = portFile + "portfolio = new String[universe.length];\n";
       portFile = portFile + "System.arraycopy(universe,0,portfolio,0,universe.length);\n\n";
       portFile = portFile + "}\n\n";
    
    
       String filterFile = new String("public void setFilterUniverse()\n{\n");
       filterFile = filterFile + "filterUniverse = new strategyParameters[" + strategies.size() + "];\n\n"; 
    
       for(int i = 0; i < strategies.size(); i++)
       {
        StrategyParameters sp = strategies.get(i);
        filterFile = filterFile + "filterUniverse["+i+"] = new strategyParameters(\""+sp.name+"\",\""+sp.insampStart+"\",\""+sp.startTrading+"\","+
        sp.startTradingInt+",\""+sp.endTrading+"\","+sp.n_obs+","+sp.n_rep+","+sp.L+","+sp.lag+","+sp.cutoff+","+sp.lambda+","+sp.expweight
        +","+sp.smooth+","+sp.decay+","+sp.decay2+","+sp.cross+","+sp.i1+","+sp.i2+","+sp.time_shift+","+sp.hybrid_weight + "," + sp.hybrid_weight_diff + "," + sp.b0trend + 
        ","+sp.sig_diff+","+sp.sig_inv+","+sp.stop_loss+");\n";
       } 
       filterFile = filterFile + "\n}\n\n";
           
       String metaFile = new String("public void setMetaUniverse()\n{\n");
    
       metaFile = metaFile + "metaUniverse = new metaFilterParameters["+metaToFile.size()+"];\n\n"; 
       for(int i = 0; i < metaToFile.size(); i++)
       {   
        METAFilterParameters metaFilter = metaToFile.get(i);
        metaFile = metaFile + "metaUniverse["+i+"]" + metaFilter.toStringFunction() + "\n";
       }
       metaFile = metaFile + "\n}\n\n";
    
   
       //now print them out to File
       out.println(portFile);
       out.println(filterFile);
       out.println(metaFile);
       
       out.close(); System.out.println("Strategy successfully saved in " + file_name);
      }
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
     }  
   }
   else
   {
     System.out.println("Must first set the filtering and meta strategies in the filtering panel");
   }
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
  
  public static double computeDrawdown(double[] ret)
  {
     int i;
     double max = -100000; 
     double[] cmax = cummax(ret);
     double[] cmaxx = new double[ret.length];
     
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
  
  public static double[] cumsum(double[] data, int n)
  {

    double[] cs = new double[n]; double sum; int k;
    
    sum=0; double min = 1000000;

    for(k=0;k<n;k++)
    {
      sum = sum+data[k]; cs[k] = sum; 
      if(cs[k] < min) {min = cs[k];}
    }


    return cs;  
  } 
  
  
  public void printAggregateStrategy()
  {
    if(n_saved_perf > 0)
    {
    
      int i,j,k; 
      int n_basket = n_saved_perf+1;
      int min_obs = mdfaEvolutionCanvas.min_obs;
      double[][] data = new double[min_obs][n_basket];
      double[] w = new double[n_basket];
      double[] target = new double[min_obs];
      String[] invests = new String[min_obs];
      //fill with current strategy first
      for(i=0;i<min_obs;i++)
      { 
        invests[min_obs-1-i] = new String(performances[performances.length - 1 - i].getDate());
        data[min_obs - 1 - i][0] = performances[performances.length - 1 - i].getReturn(); 
        invests[min_obs-1-i] = invests[min_obs-1-i] + performances[performances.length - 1 - i].toStringAssetNo();
      }
      
      for(k=0;k<n_saved_perf;k++)
      { 
       JInvestment[] temp = portfolio_invest.get(k);
       for(i=0;i<min_obs;i++)
       { 
        data[min_obs - 1 - i][k+1] = temp[temp.length - 1 - i].getReturn();
        invests[min_obs-1-i] = invests[min_obs-1-i] + " " + temp[temp.length - 1 - i].toStringAssetNo();
       } 
      }

      double[] means = new double[n_basket];
      double sum=0;
      RealVector sol; 
       
      for(i=0;i<n_basket;i++)
      {
       sum=0;
       for(j=0;j<min_obs;j++)
       {sum = sum + data[j][i];}
       means[i] = sum/min_obs;
      } 
       
      RealVector m = new ArrayRealVector(means, false);
      Covariance covComp = new Covariance(data);
      RealMatrix rm = covComp.getCovarianceMatrix();  
  
      if(uniformWeightsCheck.isSelected())
      {
        for(i=0;i<n_basket;i++) {w[i] = 1.0/n_basket;}
      }
      else if(maxSharpeWeightsCheck.isSelected())
      {
      
       try
       { 
        DecompositionSolver solver = new QRDecomposition(rm).getSolver();
        sol = solver.solve(m);
        w = sol.toArray(); 
       }
       catch(SingularMatrixException sme) 
       {
         System.out.println("Matrix singular: setting weights to uniform"); 
         w = new double[n_basket]; 
         for(i=0;i<n_basket;i++) {w[i] = 1.0/n_basket;}
       }
      
       double sumw = 0;
       for(i=0;i<w.length;i++) 
       {
        if(w[i] < 0) {w[i] = 1.0/n_basket;}       
        sumw = sumw + w[i]; 
       }
       for(i=0;i<w.length;i++) {w[i] = w[i]/sumw;}
      }   
      
      for(i=0;i<min_obs;i++)
      {
        sum = 0;
        for(k=0;k<n_basket;k++)
        {sum = sum + data[i][k]*w[k];}
        target[i] = sum;  
        if(target[i] > 0) {}
      }
      
      double[] mstd = mean_std(target); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];    
      double[] cum_port_returns = cumsum(target,target.length); 
      max_drawdown = computeDrawdown(cum_port_returns);
      try
      {
        PrintWriter out = new PrintWriter(new FileWriter("aggregateReturns.dat"));        
        for(i=0;i<min_obs;i++)
        {
          out.println(invests[i] + " " + target[i] + " " + cum_port_returns[i]);  
        }  
        out.close(); 
      
       System.out.println("Strategy successfully saved in aggregateReturns.dat");
      }
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
      

    }
  }
  
  
  public static class Result 
  {
        private final String perform_stats;
        private ArrayList<String> performance;
        private ArrayList<Double> rets,longrets,shortrets;
        public Result(String stats, ArrayList<String> per, ArrayList<Double> r, ArrayList<Double> r2,ArrayList<Double> r3) 
        {
            this.perform_stats = stats;
            this.performance = per; 
            this.rets = r; 
            this.longrets = r2; 
            this.shortrets = r3; 
        }
  }  
  

  public static Result compute(MDFAStrategyEvolution strategy) throws InterruptedException 
  {
        
     boolean output = strategy.startStrategy();
     if(output)
     {
      double ratio_trades = strategy.trade_succ_ratio/(double)strategy.returns.size();
      double sharpe = Math.sqrt(250)*strategy.mean_perf/strategy.standard_deviation;
      String perf = new String(" n_obs = " + strategy.n_obs + ", avg_rank = " + strategy.rank_coeff + ", maxdraw = " + strategy.maxdraw + ", mean = " + strategy.mean_perf + ", success_ratio = " + ratio_trades + ", sharpe " + sharpe);        
      
      return new Result(perf,strategy.date_returns,strategy.returns,strategy.longreturns,strategy.shortreturns);
     }
     else
     {
       String perf = new String("");
       System.out.println("\n\nNO MAN'S LAND!");
       return new Result(perf,strategy.date_returns,strategy.returns,strategy.longreturns,strategy.shortreturns);
     }      
  }  
  
  
  public static String[] forkDataFilePairs(File file, int m_days, int nobs, int cushion, ArrayList<String> expvar)
  {  
  
      int i,m;
          
     ArrayList<String> data_line = new ArrayList<String>();
     ArrayList<String> data_line_exp1 = new ArrayList<String>();
     new ArrayList<String>();
     
     String lines;
     double prevprice,logprice;
     
       int n_dates; 
       String strline; 
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
       ArrayList<String>forkedFiles = new ArrayList<String>();
        
       prevprice = 0;  
     //---- Get the target data file and data  
     try{  
       
     FileInputStream fin = new FileInputStream(file);
     DataInputStream din = new DataInputStream(fin);
     BufferedReader br = new BufferedReader(new InputStreamReader(din));      

     while((strline = br.readLine()) != null)
     {

       data_line.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
     }
     br.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
     
     //---- Now get each explanatory variable-------------------------------
     if(expvar.size() > 0)
     {
      try{  
       
      FileInputStream fin = new FileInputStream(new File(expvar.get(0)));
      DataInputStream din = new DataInputStream(fin);
      BufferedReader br = new BufferedReader(new InputStreamReader(din));      

       while((strline = br.readLine()) != null)
       {

       data_line_exp1.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
       }
       br.close();
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
      
     } 
     
   
     
     int n_lines = data_line.size();
     int n_linesexp1 = data_line_exp1.size();


     System.out.println("N of lines = " + n_lines + " " + n_linesexp1);
     int local_lines = (int)n_lines/m_days;
  
  
     try{
  
     int start = 0; m = 0; 
     while(start < n_lines)
     { 
       File file1 = new File(""+file.getName()+"_"+m);
       PrintWriter out = new PrintWriter(new FileWriter(file1));    
       
       n_dates = 0;
       for(i = start; i < start + local_lines; i++)
       {
        if(i < n_lines && i >= 0) 
        {
         //out.println(data_line.get(i)); 
         
         //2013-10-07 10:30:00, 3.1041381473977774, 0.0015714449299597533
         lines = new String(data_line.get(i));
         
         String[] asset1 = lines.split("[,]+");
         String[] asset2 = (data_line_exp1.get(i)).split("[,]+");
         logprice = (new Double(asset1[1])).doubleValue() - (new Double(asset2[1])).doubleValue();
         
         lines = asset1[0] + ", " + logprice + ", " + (logprice - prevprice);
         out.println(lines);
         n_dates++;
         prevprice = logprice;
         
        } 
        else {break;}        
       }
       out.close();
       start = start + local_lines - nobs - cushion;
       if(start < 0) {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       else if(n_dates < (nobs + cushion))
       {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       
       forkedFiles.add(""+file.getName()+"_"+m);
       m++;
     } 
    }
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
    
    return forkedFiles.toArray(new String[0]);
  }
    
  public static String[] forkDataFile(File file, int m_days, int nobs, int cushion)
  {
     
     int i,m;
          
     ArrayList<String> data_line = new ArrayList<String>();
     //int cushion = 51;
     int n_dates; 
       String strline; 
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
       ArrayList<String>forkedFiles = new ArrayList<String>();
        
       try{  
       
     FileInputStream fin = new FileInputStream(file);
     DataInputStream din = new DataInputStream(fin);
     BufferedReader br = new BufferedReader(new InputStreamReader(din));      

     while((strline = br.readLine()) != null)
     {

       data_line.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
     }

     int n_lines = data_line.size();
     int local_lines = (int)n_lines/m_days;
     
     
     int start = 0; m = 0; 
     while(start < n_lines)
     { 
       File file1 = new File(""+file.getName()+"_"+m);
       PrintWriter out = new PrintWriter(new FileWriter(file1));    
       
       n_dates = 0;
       for(i = start; i < start + local_lines; i++)
       {
        if(i < n_lines && i >= 0) 
        {out.println(data_line.get(i)); n_dates++;}
        else {break;}        
       }
       out.close();
       start = start + local_lines - nobs - cushion;
       if(start < 0) {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       else if(n_dates < (nobs + cushion))
       {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       
       forkedFiles.add(""+file.getName()+"_"+m);
       m++;
     } 
    }
    catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
    
    return forkedFiles.toArray(new String[0]);
  }
  
  

  
  
  public static String[] forkDataFileExp(File file, int m_days, int nobs, int cushion, ArrayList<String> expvar)
  {  
  
      int i,m;
          
     ArrayList<String> data_line = new ArrayList<String>();
     ArrayList<String> data_line_exp1 = new ArrayList<String>();
     ArrayList<String> data_line_exp2 = new ArrayList<String>();
     
     int nexp = 0;
     String lines;
     
     
       int n_dates; 
       String strline; 
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
       ArrayList<String>forkedFiles = new ArrayList<String>();
        
       //---- Get the target data file and data  
     try{  
       
     FileInputStream fin = new FileInputStream(file);
     DataInputStream din = new DataInputStream(fin);
     BufferedReader br = new BufferedReader(new InputStreamReader(din));      

     while((strline = br.readLine()) != null)
     {

       data_line.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
       
     }
     br.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
     
     //---- Now get each explanatory variable-------------------------------
     if(expvar.size() > 0)
     {
      try{  
       
      FileInputStream fin = new FileInputStream(new File(expvar.get(0)));
      DataInputStream din = new DataInputStream(fin);
      BufferedReader br = new BufferedReader(new InputStreamReader(din));      

       while((strline = br.readLine()) != null)
       {

       data_line_exp1.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
       }
       br.close();
       nexp = 1; 
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
      
     } 
     
      //---- Now get each explanatory variable 2 --------------------------------------------
     if(expvar.size() > 1)
     {
      try{  
       
      FileInputStream fin = new FileInputStream(new File(expvar.get(1)));
      DataInputStream din = new DataInputStream(fin);
      BufferedReader br = new BufferedReader(new InputStreamReader(din));      

       while((strline = br.readLine()) != null)
       {

       data_line_exp2.add(strline);
       
       tokens = strline.split(delims); 
       n_toks = tokens.length; 
        
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
       if(n_toks == 0)
       {System.out.println("End of file"); break;}
  
       }
       br.close();
       nexp = 2;
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
     }      
     
     int n_lines = data_line.size();
     int n_linesexp1 = data_line_exp1.size();
     int n_linesexp2 = data_line_exp2.size();

     System.out.println("N of lines = " + n_lines + " " + n_linesexp1 + " " + n_linesexp2);
     int local_lines = (int)n_lines/m_days;
  
  
     try{
  
     int start = 0; m = 0; 
     while(start < n_lines)
     { 
       File file1 = new File(""+file.getName()+"_"+m);
       PrintWriter out = new PrintWriter(new FileWriter(file1));    
       
       n_dates = 0;
       for(i = start; i < start + local_lines; i++)
       {
        if(i < n_lines && i >= 0) 
        {
         //out.println(data_line.get(i)); 
         
         //2013-10-07 10:30:00, 3.1041381473977774, 0.0015714449299597533
         lines = new String(data_line.get(i));
         if(nexp >=1) //only get the first
         {
            String[] line_exp1 = (data_line_exp1.get(i)).split("[,]+");
            lines = lines + ", " + line_exp1[2];
            
            if(nexp == 2)
            {
               line_exp1 = (data_line_exp2.get(i)).split("[,]+");
               lines = lines + ", " + line_exp1[2];
            }
            
         }
         out.println(lines);
         n_dates++;
         
        } 
        else {break;}        
       }
       out.close();
       start = start + local_lines - nobs - cushion;
       if(start < 0) {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       else if(n_dates < (nobs + cushion))
       {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       
       forkedFiles.add(""+file.getName()+"_"+m);
       m++;
     } 
    }
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
    
    return forkedFiles.toArray(new String[0]);
  }
  
  
  public void computeStrategy() throws InterruptedException, ExecutionException 
  {
  
    if(strategies.size() > 0 && data_files.length > 0) //should be ready to go
    {
     progressBar.setValue(0); progressBarFiles.setValue(0);
     progressBar.setMaximum(strategies.size());
     progressBarFiles.setMaximum(data_files.length);

    
     
     int k,j,i,obs;
     
     String delimsp = "[.]+";
     String delims = "[,]+";
     String[] tokens; 
     String[] date_tokens;     
     String asset; 
     int n_trade_days=0;
     n_assets = data_files.length;
     n_filters = strategies.size();
     
     String[] datafiles; 
     ArrayList<String> total_perf = new ArrayList<String>();
     ArrayList<Double> longreturns_all = new ArrayList<Double>();
     ArrayList<Double> shortreturns_all = new ArrayList<Double>();
     ArrayList<Double> returns_all = new ArrayList<Double>();      
     
     List<MDFAStrategyEvolution> objects = new ArrayList<MDFAStrategyEvolution>();
     
     for(l_filter=0;l_filter<n_filters;l_filter++) //apply the strategy 
     {
                   
      try  //for each strategy, open three different performance files
      {

        
        
        PrintWriter full_out = new PrintWriter(new FileWriter("returns_full_"+l_filter+".dat")); 
        PrintWriter long_out = new PrintWriter(new FileWriter("returns_long_"+l_filter+".dat"));
        PrintWriter short_out = new PrintWriter(new FileWriter("returns_short_"+l_filter+".dat"));
        PrintWriter asset_out = new PrintWriter(new FileWriter("asset_returns.dat"));
        
        int cusion = computeMaxTradeObs(new File(data_files[0]));
        obs = (strategies.get(l_filter)).n_obs;

        for(k_file=0;k_file<n_assets;k_file++)
        {
          if(pairs_trading) {datafiles = forkDataFilePairs(new File(data_files[k_file]), n_threads, obs, cusion, expvar);}
          else if(exp_var) {datafiles = forkDataFileExp(new File(data_files[k_file]), n_threads, obs, cusion, expvar);}
          else {datafiles = forkDataFile(new File(data_files[k_file]), n_threads, obs, cusion);}
            
          tokens = data_files[k_file].split(delimsp); //string the .dat and keep name
          asset = tokens[0];
     
          total_perf.clear();
          longreturns_all.clear();
          shortreturns_all.clear();
          returns_all.clear(); 
          objects.clear();
        

          for(j = 0; j < datafiles.length; j++) //parallelization occurs on spit files 
          { 
            
            MDFAStrategyEvolution strategy = new MDFAStrategyEvolution(10, "filter_trend_es_3.params"); 
            strategy.setLogTrans(forexData.isSelected()); 
            strategy.forex24 = data24.isSelected();
            if(n_threads == 1) {strategy.togglePrint(true);}
            
            strategy.setFinalTradeTime(final_time);
            strategy.setObservationFrequency(minute_df);      
     
            strategy.setStrategyParameters(strategies.get(l_filter), datafiles[j], tradingCost);
            objects.add(strategy);
                                   
            
          } 
        
          List<Callable<Result>> tasks = new ArrayList<Callable<Result>>();
      
          for (final MDFAStrategyEvolution object : objects) 
          {
            Callable<Result> c = new Callable<Result>() {
            @Override
            public Result call() throws Exception 
            {return compute(object);}
            };
            tasks.add(c);
          }

          ExecutorService exec = Executors.newCachedThreadPool();

          try 
          {
           List<Future<Result>> results = exec.invokeAll(tasks);
           for (Future<Result> fr : results) 
           {
             if(fr.get().perform_stats.equals(""))
             {
               System.out.println("Empty file...");
             }
             else
             {
              System.out.println(fr.get().perform_stats);
              total_perf.addAll(fr.get().performance);
              returns_all.addAll(fr.get().rets);
              longreturns_all.addAll(fr.get().longrets);
              shortreturns_all.addAll(fr.get().shortrets);
             }
           }
          } 
          finally 
          {System.out.println("shutting down exec"); }   
     
          System.out.println("Printing headers");
          //-------------------print headers
          if(k_file==0) 
          {
           n_trade_days = returns_all.size(); 
           full_out.print("Date "); long_out.print("Date "); short_out.print("Date "); 
           
           for(i = 0; i < total_perf.size()-1; i++)
           {
              date_tokens = (total_perf.get(i)).split(delims);
              full_out.print(date_tokens[0] + " "); long_out.print(date_tokens[0] + " "); short_out.print(date_tokens[0] + " ");
              asset_out.print(date_tokens[0] + " ");
           }
           date_tokens = (total_perf.get(total_perf.size()-1)).split(delims);
           full_out.println(date_tokens[0]); long_out.println(date_tokens[0]); short_out.println(date_tokens[0]); asset_out.println(date_tokens[0]);
          
           if(n_trade_days == returns_all.size())
           {full_out.print(asset + " "); long_out.print(asset + " "); short_out.print(asset + " "); asset_out.print(asset + " ");}    
           System.out.println(date_tokens[0]);
          }
          else //only keep those with same number of trading days, to be consistent  
          {
           if(n_trade_days == returns_all.size())
           {full_out.print(asset + " "); long_out.print(asset + " "); short_out.print(asset + " "); asset_out.print(asset + " ");}
           System.out.println(asset + " ");
          }
          
          System.out.println("PERFORMANCE LENGTHS:" + n_trade_days + " " + total_perf.size() + " " + returns_all.size() + " " + longreturns_all.size() + " "
            + shortreturns_all.size());
          
          if(n_trade_days == returns_all.size())
          {
           //now place the performances in their respective files 
           for(i = 0; i < total_perf.size()-1; i++)
           {
           
            System.out.println(total_perf.get(i) + " " + returns_all.get(i) + " " + longreturns_all.get(i) + " " + shortreturns_all.get(i));
            
            full_out.print(returns_all.get(i) + " "); 
            long_out.print(longreturns_all.get(i) + " ");
            short_out.print(shortreturns_all.get(i) + " ");
            
            date_tokens = (total_perf.get(i)).split(delims);
            asset_out.print(date_tokens[1] + " ");
           }
           //print the final column with new line
           
           i = total_perf.size()-1;
           full_out.println(returns_all.get(i)); 
           long_out.println(longreturns_all.get(i));
           short_out.println(shortreturns_all.get(i));
           
           date_tokens = (total_perf.get(i)).split(delims);
           asset_out.println(date_tokens[1]);
           
          }
     
          progressBarFiles.setValue(k_file+1);
          Rectangle progressRectF = progressBarFiles.getBounds(); 
          progressRectF.x = 0;
          progressRectF.y = 0;
          progressBarFiles.paintImmediately( progressRectF );        
     
     
        }
        full_out.close();
        long_out.close();
        short_out.close();
        asset_out.close();
        all_returns_computed = true;
      }  
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
      
        progressBar.setValue(l_filter);
        Rectangle progressRect = progressBar.getBounds(); 
        progressRect.x = 0;
        progressRect.y = 0;
        progressBar.paintImmediately( progressRect );       
      
     }     
    

     //---- now that the returns have been computed, load them into the strategy suite
     if(all_returns_computed) 
     {
      strategy_suite = new StrategyFilter(n_filters*3);
      for(l_filter = 0; l_filter < n_filters; l_filter++)
      {
        strategy_suite.loadFilteredReturns("returns_full_"+l_filter+".dat");
        strategy_suite.loadFilteredReturns("returns_long_"+l_filter+".dat");
        strategy_suite.loadFilteredReturns("returns_short_"+l_filter+".dat");
      }
            
      performances = strategy_suite.computeStrategyGeneral(M, nlag, meta_filter, port_normalize, equal_dist, !longOnly, 1, false, blackListCheck.isSelected(), 
      nblack, blength, black_list, filter_stock_choice);
      plotPerformances();      
      
      for(k=0; k < n_filters; k++) {stratCheck[k].setEnabled(true); stratCheck[k].setSelected(true);}
      stratShortCheck.setEnabled(true); stratShortCheck.setSelected(true);
      stratLongCheck.setEnabled(true); stratLongCheck.setSelected(true);
      stratBothCheck.setEnabled(true); stratBothCheck.setSelected(true);
      
      strategy_computed = true;
      
      //System.out.println("Printing out assets..");
      readAssetReturnFile();
      asset_plots.setEnabled(true);
      asset_plots.removeAllItems();
      //asset_plots = new JComboBox<String>();
      for(k=0;k<n_assets;k++)
      {asset_plots.addItem(stock_names.get(k)); filter_stock.addItem(stock_names.get(k));} // System.out.println(stock_names.get(k));}
      add_plots.setEnabled(true); clear_asset_plots.setEnabled(true);

      
     }

    }
  
  }
   
   
  public void loadSavedPerformances()
  {
      int k; boolean pass = true; 
      n_filters = 0;
      strategy_suite = new StrategyFilter(36);
      //for(l_filter = 0; l_filter < n_filters; l_filter++)
      l_filter = 0;
      while(l_filter < 12)
      {
        if((new File("returns_full_"+l_filter+".dat")).exists())
        {
         pass = strategy_suite.loadFilteredReturns("returns_full_"+l_filter+".dat");
         pass = strategy_suite.loadFilteredReturns("returns_long_"+l_filter+".dat");
         pass = strategy_suite.loadFilteredReturns("returns_short_"+l_filter+".dat");
         
         if(!pass) {break;}
         else {l_filter++;}
        }
        else
        {break;}
      }
      n_filters = l_filter;
      
      if(n_filters > 0)
      {
       performances = strategy_suite.computeStrategyGeneral(M, nlag, meta_filter, port_normalize, equal_dist, !longOnly, 1, false, blackListCheck.isSelected(), 
       nblack, blength, black_list, filter_stock_choice);
       plotPerformances();      
      
       for(k=0; k < n_filters; k++) {stratCheck[k].setEnabled(true); stratCheck[k].setSelected(true);}
       stratShortCheck.setEnabled(true); stratShortCheck.setSelected(true);
       stratLongCheck.setEnabled(true); stratLongCheck.setSelected(true);
       stratBothCheck.setEnabled(true); stratBothCheck.setSelected(true);
      
       strategy_computed = true;
      
       //System.out.println("Printing out assets..");
       readAssetReturnFile();
       asset_plots.setEnabled(true);
       asset_plots.removeAllItems();
       //asset_plots = new JComboBox<String>();
       //n_assets = stock_names.size();
       for(k=0;k<stock_names.size();k++)
       {asset_plots.addItem(stock_names.get(k)); filter_stock.addItem(stock_names.get(k));}// System.out.println(stock_names.get(k));}
       add_plots.setEnabled(true); clear_asset_plots.setEnabled(true);
       
       all_returns_computed = true;
       strategy_computed = true;
       
      }
      else
      {
        System.out.println("Couldn't find any performance files");
      }
  }
   
  
  
  
  public void computeBlackList() //compute once every time new strategy saved
  {
    int i,j,count,ndiff;
    if(portfolio_invest.size() == 1) //first one 
    {    
      JInvestment[] jinv = portfolio_invest.get(0);
      blength = strategy_suite.port[0].getNObs();
      nblack = jinv[0].assets.length;
      black_list = new int[nblack][blength];
      
      count = 0;
      for(i = jinv.length; i > 0; i--)
      {
         for(j = 0; j < nblack; j++)
         {black_list[j][blength-1-count] = jinv[i-1].assetsID[j];}
         count++;
      }
      ndiff = blength - count; //this SHOULD be positive
      for(i = 0; i < ndiff; i++)
      {
         for(j = 0; j < nblack; j++)
         {black_list[j][i] = -1;}
      }
      black_list_computed = true;
    }
    else if(portfolio_invest.size() > 1) //just add on to previous
    {
    
      if(black_list_computed)
      {
        JInvestment[] jinv = portfolio_invest.get(portfolio_invest.size()-1);
        
        blength = black_list[0].length;
        int[][] black_temp = new int[nblack + jinv[0].assets.length][blength];
  
        //--copy old black list 
        for(i = 0; i < blength; i++)
        {
         for(j = 0; j < nblack; j++)
         {black_temp[j][i] = black_list[j][i];}
        }
  
  
        count = 0;
        for(i = jinv.length; i > 0; i--)
        {
         for(j = 0; j < jinv[0].assetsID.length; j++)
         {black_temp[nblack+j][blength-1-count] = jinv[i-1].assetsID[j];}
         count++;
        }
        
        ndiff = blength - count; //this SHOULD be positive
        for(i = 0; i < ndiff; i++)
        {
         for(j = 0; j < jinv[0].assetsID.length + nblack; j++)
         {black_temp[j][i] = -1;}
        }
        
        nblack = nblack + jinv[0].assets.length; 
        black_list = black_temp;    
        black_list_computed = true;  
      }
    }       
  } 
   
  
  public void uploadFromFilterFile(File file)
  {
    
	 String[] tokens; String delims = "[ ]+";
     String strline;  
     strategies = new ArrayList<StrategyParameters>();
    
     try
     {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           //names = new ArrayList<String>();
           while((strline = br.readLine()) != null)
           {
             tokens = strline.split(delims);
             System.out.println(strline);
             System.out.println(tokens.length);
                          
             if(tokens.length == 24)
             {
              StrategyParameters sp = new StrategyParameters("Filter_"+strategies.size(), tokens[0], tokens[1], (new Integer(tokens[2])).intValue(),
                          tokens[3], (new Integer(tokens[4])).intValue(), (new Integer(tokens[6])).intValue(),
                          (new Integer(tokens[5])).intValue(), (new Integer(tokens[14])).intValue(),                           
                          (new Double(tokens[7])).doubleValue(), (new Double(tokens[12])).doubleValue(), 
                          (new Double(tokens[13])).doubleValue(),                          
                          (new Double(tokens[8])).doubleValue(),  (new Double(tokens[9])).doubleValue(), 
                          (new Double(tokens[10])).doubleValue(),  (new Double(tokens[11])).doubleValue(), 
                          (new Integer(tokens[15])).intValue(), (new Integer(tokens[16])).intValue(),
                          (new Double(tokens[17])).doubleValue(),  (new Double(tokens[18])).doubleValue(),
                          (new Double(tokens[19])).doubleValue(), (new Integer(tokens[20])).intValue(), (new Integer(tokens[21])).intValue(),
                          (new Integer(tokens[22])).intValue(), (new Integer(tokens[23])).intValue());
                          
  
             strategies.add(sp);    
             stratCanvas.addTextNewLine(strategies.get(strategies.size()-1).toString(), Color.GREEN);
             }           
           }
    
           din.close();
       }               
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
  }
  
  
  public int computeMaxTradeObs(File file)
  {
  
       String strline; 
       String prev_date,prev_time; 
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
       String date_delims = "[ ]+"; 
       int n_intervals = 0; 
       int max_count = -10;
       new ArrayList<String>();
 
     prev_date = new String(" ");  prev_time = new String(" ");
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
       
       
       String[] datetime = tokens[0].split(date_delims);
       if(datetime[0].equals(prev_date))
       {
        n_intervals++;
        minute_df = minute_diff(datetime[1],prev_time);        
       }
       else
       {
          if(n_intervals > max_count) {max_count = n_intervals;}
          n_intervals = 0; prev_date = datetime[0]; 
          final_time = new String(prev_time);
       }
       prev_time = datetime[1];
      } 
      br.close();
       //n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
     
     return max_count; 
  }
    
    
    
  public int minute_diff(String t1, String t2)
  {
     String[] time1 = t1.split("[:]+");
     String[] time2 = t2.split("[:]+");
     
     int min1 = (new Integer(time1[1])).intValue();
     int min2 = (new Integer(time2[1])).intValue();
     
     int min_diff = Math.abs(min1 - min2);
     if(min_diff == 0) min_diff = 60;
     
     return min_diff;
  }  
  
 
  
 
  
  
  
  
  public void saveStrategy()
  {
  
    if(strategy_computed)
    {
     int k; String fnString = "\""; 
     portfolio_invest.add(performances);
     n_saved_perf++; 
    
     ArrayList<Integer> fnInd = new ArrayList<Integer>();
     System.out.println("Saving strategy. We now have " + n_saved_perf + " performance(s) saved");
     
     FileFilter ft = new FileNameExtensionFilter("Strategy files", "strat");
     fc.addChoosableFileFilter(ft);       
     int returnVal = fc.showSaveDialog(frame);
      
     if(returnVal == javax.swing.JFileChooser.APPROVE_OPTION)
     {
       File saved_file = fc.getSelectedFile();
       String file_name = saved_file.toString();
  
      //--- Print filter strategies first 
      //save strategy to file 
      //M, nlag, meta_filter, port_normalize, equal_dist
      try
      {
       PrintWriter out = new PrintWriter(new FileWriter(file_name));

       out.println("Filtering strategies: Long-Short (" + stratBothCheck.isSelected() + "), Long (" + stratLongCheck.isSelected()
                    + "), Short (" + stratShortCheck.isSelected() + ")\n");
                   
       out.print("Filter number: ");            
       for(k=0;k<strategies.size();k++)
       {
        if(stratCheck[k].isSelected())
        {
          out.println(strategies.get(k).toString());
        }
       }
       if(strategies.size() > 0)
       {
         for(k=0;k<stratCheck.length;k++) 
         if(stratCheck[k].isSelected()) 
         {
           out.print(k+" ");
           fnInd.add(new Integer(k));          
         }  
       }
      
       out.println("\nMeta-Filtering: Portfolio size (" + M + "), Lag (" + nlag + "), Filter (" + metaFilterCombo.getSelectedItem() 
                  + "), Weights (" + portWeightsCombo.getSelectedItem() + "), Signal Reversal (" + longOnlyBox.isSelected() 
                  + "), Equal Distribution (" + equalDistCheck.isSelected() + ")\n");
      
       out.println("Assets:");
       for(k=0;k<stock_names.size();k++)
       {out.println(stock_names.get(k));}
      
       out.close(); System.out.println("Strategy successfully saved in " + file_name);
      }
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
       
    }
  
    //--- now save to ArrayList ------ 
    for(int i=0;i<fnInd.size()-1;i++)
    {fnString = fnString + fnInd.get(i)+",";}
    fnString = fnString + fnInd.get(fnInd.size()-1)+"\"";    
    
                                 
    METAFilterParameters meta = new METAFilterParameters(stratBothCheck.isSelected(),stratLongCheck.isSelected(),stratShortCheck.isSelected(), 
                                (String)filter_stock.getSelectedItem(), fnString, M, nlag, metaFilterCombo.getSelectedIndex(), portWeightsCombo.getSelectedIndex(),
                                longOnlyBox.isSelected(), equalDistCheck.isSelected());
  
  
    metaToFile.add(meta);
     
    saving_perf = true;
    plotPerformances();
    
    aggregatePerformance.setEnabled(true);
    uniformWeightsCheck.setEnabled(true);
    maxSharpeWeightsCheck.setEnabled(true);
    
    blackListCheck.setEnabled(true); 
    realretsCheck.setEnabled(true);
    computeBlackList();
    
   }
   else
   {System.out.println("No performances computed or uploaded into system yet");}
  }
   
  public void deletePortfolioStrategy()
  {
    portfolio_invest.clear();
  } 
   
   
  public void readAssetReturnFile()
  {
  
       
       
     String[] tokens; String delims = "[ ]+";
     int n_toks; String strline; int i;
     
     stock_names = new ArrayList<String>();
     stock_returns_all = new ArrayList<double[]>();
     int count = 0; 
     try
     {
          
           FileInputStream fin = new FileInputStream(new File("asset_returns.dat"));
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           //names = new ArrayList<String>();
           while((strline = br.readLine()) != null)
           {
             
             if(count == 0) //first row are the dates
             {
               strline.split(delims); 
             }  
             else
             {
             
              tokens = strline.split(delims); 
              n_toks = tokens.length; //length should be number of returns + 1
              //System.out.println("n_toks length = " + n_toks);
             
              double[] returns = new double[n_toks-1];
              stock_names.add(tokens[0]); //the first value is the asset name
              //System.out.println(tokens[0]);
              for(i = 1; i < n_toks; i++)
              {returns[i-1] = (new Double(tokens[i])).doubleValue();} 
              stock_returns_all.add(returns);   
             }
             count++;
           }
           
           //asset_names.add(names.toArray(new String[0])); 
           
           din.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
       
       
  
  }
   
  
  public void printStrategyParameters()
  {
     if(strategies.size() > 0)
     {
       String file_name = new String("saved_strategies.txt");
       
       try
       {
        PrintWriter out = new PrintWriter(new FileWriter(new File(file_name)));
       
        for(int i=0; i < strategies.size(); i++)
        {out.println(strategies.get(i).toStringText());}
       
        out.close(); System.out.println("Strategies successfully saved in " + file_name);
       }
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
       
     }
  }
  
  
   
   
  public void recomputeStrategyReturns()
  {
     int l;  
     if(all_returns_computed && strategy_computed) 
     {
      strategy_suite = new StrategyFilter(n_filters*3);
      
      for(l = 0; l < n_filters; l++)
      {
       if(stratCheck[l].isSelected())
       {
        //System.out.println("loading strategy " + l);
        if(stratBothCheck.isSelected())
			strategy_suite.loadFilteredReturns("returns_full_"+l+".dat");
        if(stratLongCheck.isSelected())
			strategy_suite.loadFilteredReturns("returns_long_"+l+".dat");
        if(stratShortCheck.isSelected())
			strategy_suite.loadFilteredReturns("returns_short_"+l+".dat");
       } 
      }
            
      if(strategy_suite.n_filters == 0) //load in default 
      {
        System.out.println("At least one strategy must be loaded idiot... Loading default strategy");
        stratBothCheck.setSelected(true);
        stratCheck[0].setSelected(true);
      }      
            
      performances = strategy_suite.computeStrategyGeneral(M, nlag, meta_filter, port_normalize, equal_dist, !longOnly, 1, false, blackListCheck.isSelected(), 
      nblack, blength, black_list,filter_stock_choice);
      plotPerformances();      
      strategy_computed = true;
     }  
  
  }
   
   
   
  public void recomputeStrategy()
  {
    performances = strategy_suite.computeStrategyGeneral(M, nlag, meta_filter, port_normalize, equal_dist, !longOnly, 1, false, blackListCheck.isSelected(), 
    nblack, blength, black_list,filter_stock_choice);
    plotPerformances();
  }
   
   
  //------- Send the data to the EvolutionCanvas ------------------------------- 
  public void plotPerformances()
  {
    int i; int totlength = performances.length;
    int n_trading_days = totlength;
    double[] rets = new double[performances.length];
    String[] info = new String[performances.length];
    String[] asset_info = new String[performances.length];
    String[] meta_info = new String[performances.length];
    int npos = 0;
   
    for(i = 0; i < totlength; i++)
    {
      rets[i] = performances[i].getReturn();
      info[i] = performances[i].toStringDateReturn();
      asset_info[i] = performances[i].toStringAsset();
      meta_info[i] = performances[i].toStringMeta();
      
      //System.out.println(info[i]);
      if(rets[i] > 0) {npos++;}
    }
    
    
    
    double[] mstd = mean_std(rets); 
    sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];
    
    double[] cum_port_returns = cumsum(rets,n_trading_days); 
    //double[] stock_returns = cumsum(stock_rets,n_trading_days);
      
    max_drawdown = computeDrawdown(cum_port_returns);
    
    double bRatio = (double)npos/totlength;

    
      //mdfaAnalysisCanvas.setDates(dates);
      //mdfaAnalysisCanvas.setStockReturns(stock_returns);
      if(mdfaEvolutionCanvas.sim_series.size() == 0 || saving_perf)
      {
       mdfaEvolutionCanvas.addSeries(cum_port_returns, new String(""+df2.format(sharpe_ratio)+", " +df2.format(max_drawdown)+", "+df.format(bRatio))); 
       
       if(saving_perf)
       {mdfaEvolutionCanvas.saveAsset();}
       
       saving_perf = false;
      } 
      else mdfaEvolutionCanvas.setSeries(n_saved_perf,cum_port_returns,new String(""+df2.format(sharpe_ratio)+", " +df2.format(max_drawdown)+", "+df.format(bRatio)));
      mdfaEvolutionCanvas.setDates(info);
      mdfaEvolutionCanvas.setAsset(asset_info);
      mdfaEvolutionCanvas.setMeta(meta_info);
      for(i=0;i<mdfaEvolutionCanvas.sim_series.size();i++) {jCheckBox[i].setEnabled(true);}// jCheckBox[i].setSelected(true);}     
  
      if(aggregatePerformance.isSelected() && n_saved_perf > 0)
      {fuseStrategies();}
  
  }  
    

  public void fuseStrategies()
  {
    if(n_saved_perf > 0)
    {
    
      int i,j,k; int npos = 0;
      int n_basket = n_saved_perf+1;
      int min_obs = mdfaEvolutionCanvas.min_obs;
      double[][] data = new double[min_obs][n_basket];
      double[] w = new double[n_basket];
      double[] target = new double[min_obs];
      //fill with current strategy first
      for(i=0;i<min_obs;i++)
      { 
        data[min_obs - 1 - i][0] = performances[performances.length - 1 - i].getReturn();        
      }
      
      for(k=0;k<n_saved_perf;k++)
      { 
       JInvestment[] temp = portfolio_invest.get(k);
       for(i=0;i<min_obs;i++)
       { 
        data[min_obs - 1 - i][k+1] = temp[temp.length - 1 - i].getReturn();
       } 
      }

      double[] means = new double[n_basket];
      double sum=0;
      RealVector sol; 
       
      for(i=0;i<n_basket;i++)
      {
       sum=0;
       for(j=0;j<min_obs;j++)
       {sum = sum + data[j][i];}
       means[i] = sum/min_obs;
      } 
       
      RealVector m = new ArrayRealVector(means, false);
      Covariance covComp = new Covariance(data);
      RealMatrix rm = covComp.getCovarianceMatrix();  
  
      if(uniformWeightsCheck.isSelected())
      {
        for(i=0;i<n_basket;i++) {w[i] = 1.0/n_basket;}
      }
      else if(maxSharpeWeightsCheck.isSelected())
      {
      
       try
       { 
        DecompositionSolver solver = new QRDecomposition(rm).getSolver();
        sol = solver.solve(m);
        w = sol.toArray(); 
       }
       catch(SingularMatrixException sme) 
       {
         System.out.println("Matrix singular: setting weights to uniform"); 
         w = new double[n_basket]; 
         for(i=0;i<n_basket;i++) {w[i] = 1.0/n_basket;}
       }
      
       double sumw = 0;
       for(i=0;i<w.length;i++) 
       {
        if(w[i] < 0) {w[i] = 1.0/n_basket;}       
        sumw = sumw + w[i]; 
       }
       for(i=0;i<w.length;i++) {w[i] = w[i]/sumw;}
      }   
      
      for(i=0;i<min_obs;i++)
      {
        sum = 0;
        for(k=0;k<n_basket;k++)
        {sum = sum + data[i][k]*w[k];}
        target[i] = sum;  
        if(target[i] > 0) {npos++;}
      }
      
      double[] mstd = mean_std(target); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];    
      double[] cum_port_returns = cumsum(target,min_obs); 
      

      
      max_drawdown = computeDrawdown(cum_port_returns);
     
      if(realrets) 
      {cum_port_returns = cumsum(target,target.length);}     
      
      double bRatio = (double)npos/min_obs;      
      
      mdfaEvolutionCanvas.addAggregate(cum_port_returns, new String(""+df2.format(sharpe_ratio)+", " +df2.format(max_drawdown)+", "+df.format(bRatio)));
      
    }
  }  
    
  public void addAssetPlot(int c)
  {
    if(c < stock_returns_all.size())
    {mdfaEvolutionCanvas.setStockReturns(cumsum(stock_returns_all.get(c),stock_returns_all.get(c).length));}
  }  
    
    
    
   @SuppressWarnings({ "unchecked", "rawtypes" })
public void setupDataPanel()
   {
   
        dataConfigPanel = new JPanel();
        
        assetFilePanel = new javax.swing.JPanel();
        loadMetaButton = new javax.swing.JButton();
        metafileAssetText = new javax.swing.JTextField();
        assetUnivButton = new javax.swing.JRadioButton();
        fileNameScroll = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,10);
        fileNameScroll.setValue(0); fileNameScroll.setUnitIncrement(1);
        fileNameScroll1 = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,10);
        fileNameScroll1.setValue(0); fileNameScroll1.setUnitIncrement(1);        
        
        tradeTimeText = new JTextField("11:00:00");
        tradeTime24 = new String("11:00:00");
        
        fileNameText = new javax.swing.JTextField();
        metafileAssetText1 = new javax.swing.JTextField();
        loadMetaButton1 = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        downloadUniversePanel = new javax.swing.JPanel();
        downloadUnivCheck = new javax.swing.JRadioButton();
        metafileAssetText2 = new javax.swing.JTextField();
        loadMetaButton2 = new javax.swing.JButton();
       
        fileNameText1 = new javax.swing.JTextField();
        connectIQButton = new javax.swing.JButton();
        downloadIQButton = new javax.swing.JButton();
        noYearsCombo = new javax.swing.JComboBox();
        noYearsLabel = new javax.swing.JLabel();
        startTimeLabel = new javax.swing.JLabel();

        assetFilePanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Upload Asset Universe ", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, null, new java.awt.Color(0, 255, 240)));
        assetFilePanel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N

        loadMetaButton.setText("Choose File");

        assetUnivButton.setText("Upload asset universe from meta file");
        assetUnivButton.setSelected(true);
  

        loadMetaButton1.setText("Choose File");

        jLabel1.setText("Additional Explanatory Variable");

        javax.swing.GroupLayout assetFilePanelLayout = new javax.swing.GroupLayout(assetFilePanel);
        assetFilePanel.setLayout(assetFilePanelLayout);
        assetFilePanelLayout.setHorizontalGroup(
            assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(assetFilePanelLayout.createSequentialGroup()
                .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(assetFilePanelLayout.createSequentialGroup()
                        .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(loadMetaButton, javax.swing.GroupLayout.PREFERRED_SIZE, 109, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(fileNameScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 103, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(fileNameText)
                            .addComponent(metafileAssetText)))
                    .addGroup(assetFilePanelLayout.createSequentialGroup()
                        .addComponent(loadMetaButton1, javax.swing.GroupLayout.PREFERRED_SIZE, 112, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(metafileAssetText1))
                    .addComponent(assetUnivButton)
                    .addComponent(jLabel1))
                .addGap(17, 17, 17))
        );
        assetFilePanelLayout.setVerticalGroup(
            assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(assetFilePanelLayout.createSequentialGroup()
                .addContainerGap(14, Short.MAX_VALUE)
                .addComponent(assetUnivButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(loadMetaButton)
                    .addComponent(metafileAssetText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(assetFilePanelLayout.createSequentialGroup()
                        .addGap(8, 8, 8)
                        .addComponent(fileNameText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(assetFilePanelLayout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fileNameScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 23, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jLabel1)
                .addGap(8, 8, 8)
                .addGroup(assetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(loadMetaButton1)
                    .addComponent(metafileAssetText1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(27, 27, 27))
        );

        downloadUniversePanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Download Asset Universe ", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, null, new java.awt.Color(0, 252, 255)));

        downloadUnivCheck.setText("Download asset universe from IQFeed");

        loadMetaButton2.setText("Choose File");

  

        connectIQButton.setText("Connect to IQFeed");

        downloadIQButton.setText("Download");

        noYearsCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9"}));
        
        
        noYearsLabel.setText("Number of Years");
        forexData = new JCheckBox("Forex");
        forexData.setToolTipText("Downloads 24-hour Forex Data 1am to 23pm"); 
        
        futuresData = new JCheckBox("Futures");
        futuresData.setToolTipText("Downloads 16-hour Futures Data from 3am to 16pm");         
        
        data5 = new JCheckBox("5min");
        data30 = new JCheckBox("30min");
        data60 = new JCheckBox("1hr");
        data24 = new JCheckBox("24hr");
        data24.setToolTipText("Gets daily data at 11am");
        zeroButton = new JCheckBox("Zero Open");
        zeroButton.setToolTipText("If market gap at open, zero the log-return value to make continuous");
        
        forexData.setSelected(false); 
        data5.setSelected(false);
        data24.setSelected(false);
        zeroButton.setSelected(false);
        data30.setSelected(false);
        data60.setSelected(false);
        futuresData.setSelected(false);
        
        javax.swing.GroupLayout downloadUniversePanelLayout = new javax.swing.GroupLayout(downloadUniversePanel);
        downloadUniversePanel.setLayout(downloadUniversePanelLayout);
        downloadUniversePanelLayout.setHorizontalGroup(
            downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                .addComponent(downloadUnivCheck)
                .addGap(0, 0, Short.MAX_VALUE))
            .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(fileNameScroll1, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(loadMetaButton2, javax.swing.GroupLayout.Alignment.LEADING))
                    .addComponent(noYearsLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                        .addComponent(noYearsCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(forexData)
                        .addComponent(futuresData)
                        .addComponent(data5)
                        .addComponent(data30)
                        .addComponent(data60)
                        .addComponent(data24)
                        .addComponent(zeroButton)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(fileNameText1)
                    .addComponent(metafileAssetText2)))
            .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                .addComponent(connectIQButton, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(downloadIQButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        downloadUniversePanelLayout.setVerticalGroup(
            downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(downloadUnivCheck)
                .addGap(18, 18, 18)
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(loadMetaButton2)
                    .addComponent(metafileAssetText2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                        .addGap(8, 8, 8)
                        .addComponent(fileNameText1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(downloadUniversePanelLayout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fileNameScroll1, javax.swing.GroupLayout.PREFERRED_SIZE, 23, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(noYearsCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(noYearsLabel)
                    .addComponent(forexData)
                        .addComponent(futuresData)
                        .addComponent(data5)
                        .addComponent(data30)
                        .addComponent(data60)
                        .addComponent(data24)
                        .addComponent(zeroButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(downloadUniversePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(connectIQButton)
                    .addComponent(downloadIQButton))
                .addContainerGap(14, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(dataConfigPanel);
        dataConfigPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(downloadUniversePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(assetFilePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(assetFilePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(downloadUniversePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        
       boolean sel = true; 
       loadMetaButton.setEnabled(sel); fileNameScroll.setEnabled(sel);
       downloadIQButton.setEnabled(!sel); loadMetaButton2.setEnabled(!sel);
       fileNameScroll1.setEnabled(!sel); connectIQButton.setEnabled(!sel); 
       noYearsCombo.setEnabled(!sel);  
        
       class MyItemListener implements ItemListener {
        public void itemStateChanged(ItemEvent e)
        {
         boolean sel; //computeFilter = true;
         Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}        
        
         if(source == assetUnivButton)
         {
           loadMetaButton.setEnabled(sel);
           fileNameScroll.setEnabled(sel);
           downloadIQButton.setEnabled(!sel);  
           loadMetaButton2.setEnabled(!sel);
           fileNameScroll1.setEnabled(!sel);
           connectIQButton.setEnabled(!sel);
           noYearsCombo.setEnabled(!sel);
           forexData.setEnabled(!sel);
           futuresData.setEnabled(!sel);
           data24.setEnabled(!sel);
           data30.setEnabled(!sel);
           zeroButton.setEnabled(!sel);
         }
         else if(source == downloadUnivCheck)
         {  
           loadMetaButton2.setEnabled(!sel);
           fileNameScroll.setEnabled(!sel);
           downloadIQButton.setEnabled(sel);  
           loadMetaButton2.setEnabled(sel);
           fileNameScroll1.setEnabled(sel);
           connectIQButton.setEnabled(sel);
           noYearsCombo.setEnabled(sel);
           forexData.setEnabled(sel);
           futuresData.setEnabled(sel);
           data24.setEnabled(sel);
           data30.setEnabled(sel);
           zeroButton.setEnabled(sel);
         } 
        }
       }; 
         
       ButtonGroup filterGroup = new ButtonGroup();
       filterGroup.add(assetUnivButton); assetUnivButton.addItemListener(new MyItemListener());
       filterGroup.add(downloadUnivCheck); downloadUnivCheck.addItemListener(new MyItemListener());        
        
        //-------------Listeners------------------------------
       ActionListener buttonActionListener2 = new ActionListener() 
       {
        File file; int val;
        public void actionPerformed(ActionEvent event)
        {
          if(event.getSource() == loadMetaButton2)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               setPortfolioFile(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }  
          if(event.getSource() == loadMetaButton)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               setSPPortfolioFile(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          } 
          if(event.getSource() == connectIQButton)
          {
            startIQConnect IQ = new startIQConnect();
            IQ.run();        
          }
          if(event.getSource() == downloadIQButton)
          {
            downloadDataIQFeed();           
          }
          
        }
       };

       connectIQButton.addActionListener(buttonActionListener2);
       downloadIQButton.addActionListener(buttonActionListener2);
       
       loadMetaButton.addActionListener(buttonActionListener2);
       loadMetaButton2.addActionListener(buttonActionListener2);
       
       metafileAssetText2.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_metaparameterFile(metafileAssetText2.getText());}} );
       metafileAssetText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_spparameterFile(metafileAssetText.getText());}} );        
        
       fileNameText1.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_manualDataFile(fileNameText1.getText());}} );      
        
       fileNameText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_manualDataFile(fileNameText.getText());}} );        
        
       tradeTimeText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_time(tradeTimeText.getText());}} ); 
        
       metafileAssetText1.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_manualExpDataFile(metafileAssetText1.getText());}}); 
        
        
       fileNameScroll.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               fileNameText.setText(portfolio[(((JScrollBar)e.getSource()).getValue())]);   
            }
         }); 
        
       fileNameScroll1.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               fileNameText1.setText(portfolio[(((JScrollBar)e.getSource()).getValue())]);   
            }
         });        
        
        
    
  }
   
 
  @SuppressWarnings({ "unchecked", "rawtypes" })
public void setupMetaFilterPanel()
  {

        metaFilterConfigPanel = new JPanel();
        metaFilterPanel = new javax.swing.JPanel();
        sizePortfolioLabel = new javax.swing.JLabel();
        sizePortfolioScroll = new javax.swing.JScrollBar();
        sizePortfolioText = new javax.swing.JTextField();
        metaFilterLabel = new javax.swing.JLabel();
        metaFilterCombo = new javax.swing.JComboBox();
        lagDaysText = new javax.swing.JTextField();
        
        loadPerformances = new JButton("Load Saved Performances");
        saveStrategyButton = new JButton("Save Strategy");
        
        lagDaysLabel = new javax.swing.JLabel();
        portWeightsLabel = new javax.swing.JLabel();
        portWeightsCombo = new javax.swing.JComboBox();
        longOnlyBox = new javax.swing.JCheckBox();
        equalDistCheck = new javax.swing.JCheckBox();

        metaFilterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Meta-Filter Strategy", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, null, new java.awt.Color(0, 252, 255)));

        sizePortfolioLabel.setText("Size of Evolution Portfolio");

        sizePortfolioScroll.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sizePortfolioScroll.setOrientation(javax.swing.JScrollBar.HORIZONTAL);

        sizePortfolioText.setText("1");

        metaFilterLabel.setText("Meta-Filter");

        metaFilterCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        metaFilterCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Sharpe", "Rank", "TradeRatio" }));

        
        stratLongCheck = new JCheckBox("Long Only"); stratLongCheck.setEnabled(false); stratLongCheck.setSelected(false);
        stratShortCheck = new JCheckBox("Short Only"); stratShortCheck.setEnabled(false); stratShortCheck.setSelected(false);
        stratBothCheck = new JCheckBox("Long-Short"); stratBothCheck.setEnabled(false); stratBothCheck.setSelected(false);
        
        stratCheck = new JCheckBox[18];
        for(int k = 0; k < 18; k++) {stratCheck[k] = new JCheckBox(""+k); stratCheck[k].setEnabled(false); stratCheck[k].setSelected(false);}        
        exportStrategy = new JButton("Evolutionize"); 
        exportResult = new JButton("Export Result");
        
        
    
        JPanel stratCon = new JPanel(); GroupLayout paramLayout;
        paramLayout = new GroupLayout(stratCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(stratBothCheck).addComponent(stratLongCheck).addComponent(stratShortCheck));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(stratBothCheck).addComponent(stratLongCheck).addComponent(stratShortCheck)));
        stratCon.setLayout(paramLayout);    
    
        
        JPanel stratCon2 = new JPanel();
        paramLayout = new GroupLayout(stratCon2); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(stratCheck[0]).addComponent(stratCheck[1]).addComponent(stratCheck[2]).addComponent(stratCheck[3]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(stratCheck[0]).addComponent(stratCheck[1]).addComponent(stratCheck[2]).addComponent(stratCheck[3])));
        stratCon2.setLayout(paramLayout);   
        
    
        JPanel stratCon3 = new JPanel();
        paramLayout = new GroupLayout(stratCon3); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(stratCheck[4]).addComponent(stratCheck[5]).addComponent(stratCheck[6]).addComponent(stratCheck[7]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(stratCheck[4]).addComponent(stratCheck[5]).addComponent(stratCheck[6]).addComponent(stratCheck[7])));
        stratCon3.setLayout(paramLayout);       
   
        JPanel stratCon4 = new JPanel();
        paramLayout = new GroupLayout(stratCon4); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(stratCheck[8]).addComponent(stratCheck[9]).addComponent(stratCheck[10]).addComponent(stratCheck[11]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(stratCheck[8]).addComponent(stratCheck[9]).addComponent(stratCheck[10]).addComponent(stratCheck[11])));
        stratCon4.setLayout(paramLayout);       
   
        stratSelectPanel = new JPanel();
        paramLayout = new GroupLayout(stratSelectPanel); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(stratCon)
              .addComponent(stratCon2)
              .addComponent(stratCon3)
              .addComponent(stratCon4)
              .addComponent(loadPerformances)
              .addComponent(saveStrategyButton)
              .addComponent(exportStrategy))
              .addComponent(exportResult));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
              .addComponent(stratCon)
              .addComponent(stratCon2)
              .addComponent(stratCon3)
              .addComponent(stratCon4)
              .addComponent(loadPerformances)
              .addComponent(saveStrategyButton)
              .addComponent(exportStrategy)
              .addComponent(exportResult));
        stratSelectPanel.setLayout(paramLayout);           
        Border raisedetched = BorderFactory.createRaisedBevelBorder();
        stratSelectPanel.setBorder(raisedetched);
           
        
        lagDaysText.setText("10");  
        lagDaysScroll = new JScrollBar(JScrollBar.HORIZONTAL,10,5,10,150);
        lagDaysScroll.setValue(10); nlag = 10;
        lagDaysScroll.setUnitIncrement(1);
        
        sizePortfolioScroll = new JScrollBar(JScrollBar.HORIZONTAL,5,1,1,51);
        sizePortfolioScroll.setValue(1); 
        sizePortfolioScroll.setUnitIncrement(1);
        
        JLabel indiv = new JLabel("Select Symbols:");
        filter_stock = new JComboBox<String>(); 
        filter_stock.addItem("all");
        
        lagDaysLabel.setText("Number of Lag Days");
        portWeightsLabel.setText("Portfolio Weights");

        portWeightsCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        portWeightsCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Uniform", "Statistic", "Max Sharpe" }));

        longOnlyBox.setText("Negative Signal Reversed");

        equalDistCheck.setText("Equal Distribution in Strategies");

        javax.swing.GroupLayout metaFilterPanelLayout = new javax.swing.GroupLayout(metaFilterPanel);
        metaFilterPanel.setLayout(metaFilterPanelLayout);
        metaFilterPanelLayout.setHorizontalGroup(
            metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(metaFilterPanelLayout.createSequentialGroup()
                .addComponent(lagDaysLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 184, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(lagDaysScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 87, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(lagDaysText))
            .addGroup(metaFilterPanelLayout.createSequentialGroup()
                .addComponent(sizePortfolioLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sizePortfolioScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 87, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sizePortfolioText))
            .addGroup(metaFilterPanelLayout.createSequentialGroup()
                .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(metaFilterPanelLayout.createSequentialGroup()
                        .addComponent(metaFilterLabel)
                        .addGap(18, 18, 18)
                        .addComponent(metaFilterCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(metaFilterPanelLayout.createSequentialGroup()
                        .addComponent(portWeightsLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(portWeightsCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(longOnlyBox)
                    .addComponent(equalDistCheck)
                    .addComponent(indiv)
                    .addComponent(filter_stock))
                .addGap(0, 0, Short.MAX_VALUE))
        );
        metaFilterPanelLayout.setVerticalGroup(
            metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(metaFilterPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(metaFilterLabel)
                    .addComponent(metaFilterCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(metaFilterPanelLayout.createSequentialGroup()
                        .addGap(4, 4, 4)
                        .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(sizePortfolioScroll, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(sizePortfolioLabel)))
                    .addComponent(sizePortfolioText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(lagDaysText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(metaFilterPanelLayout.createSequentialGroup()
                        .addGap(4, 4, 4)
                        .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(metaFilterPanelLayout.createSequentialGroup()
                                .addGap(2, 2, 2)
                                .addComponent(lagDaysLabel))
                            .addComponent(lagDaysScroll, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addGap(18, 18, 18)
                .addGroup(metaFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(portWeightsLabel)
                    .addComponent(portWeightsCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(longOnlyBox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(equalDistCheck)
                .addComponent(indiv)
                .addComponent(filter_stock)
                .addContainerGap(149, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(metaFilterConfigPanel);
        metaFilterConfigPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(metaFilterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(metaFilterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        
        lagDaysScroll.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {  
               setLagDays(((JScrollBar)e.getSource()).getValue());
               lagDaysText.setText(""+nlag);
            }
        });
        
        sizePortfolioScroll.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {  
               setSizePortfolio(((JScrollBar)e.getSource()).getValue());
               sizePortfolioText.setText(""+M);
            }
        });
        
       ActionListener comboActionListener = new ActionListener() 
       {
        public void actionPerformed(ActionEvent event)
        {
          if(event.getSource() == metaFilterCombo)
          {meta_filter = metaFilterCombo.getSelectedIndex(); if(strategy_computed) {recomputeStrategy();}}
          else if(event.getSource() == portWeightsCombo)
          {port_normalize = portWeightsCombo.getSelectedIndex(); if(strategy_computed) {recomputeStrategy();}}
          else if(event.getSource() == filter_stock)
          {filter_stock_choice = (String)filter_stock.getSelectedItem(); if(strategy_computed) {recomputeStrategy();}}
          
          else if(event.getSource() == loadPerformances)
          {loadSavedPerformances();}
          else if(event.getSource() == saveStrategyButton)
          {saveStrategy();}
        }  
       };
       
       ItemListener mitemListener = new ItemListener()
       {
         public void itemStateChanged(ItemEvent e)
         {         
          boolean sel; Object source = e.getItemSelectable();
          if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
          else{sel = true;}
        
          if(source == longOnlyBox) {longOnly = sel;}
          else if(source == equalDistCheck) {equal_dist = sel;}
          if(strategy_computed) {recomputeStrategy();}
         }
       };
       
       saveStrategyButton.addActionListener(comboActionListener);
       loadPerformances.addActionListener(comboActionListener);
       longOnlyBox.addItemListener(mitemListener); 
       equalDistCheck.addItemListener(mitemListener);
       metaFilterCombo.addActionListener(comboActionListener);
       portWeightsCombo.addActionListener(comboActionListener);
       filter_stock.addActionListener(comboActionListener);
 
 
      ItemListener checkItemListener = new ItemListener() {
       public void itemStateChanged(ItemEvent e)
       {
         int i; Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){} 
         else {}

         for(i=0;i<stratCheck.length;i++)
         {
          if(source == stratCheck[i]) {recomputeStrategyReturns();}          
         }
         if(source == stratBothCheck)
         {recomputeStrategyReturns();}
         else if(source == stratLongCheck)
         {recomputeStrategyReturns();}
         else if(source == stratShortCheck)
         {recomputeStrategyReturns();}
       }
      };
      
      for(int k = 0; k<18;k++) {stratCheck[k].addItemListener(checkItemListener);}
      stratLongCheck.addItemListener(checkItemListener);
      stratShortCheck.addItemListener(checkItemListener);
      stratBothCheck.addItemListener(checkItemListener); 
 
 
 
 
 
 
   }
    
    //------------- Setup strategyPanel -------------------------------------
    
   public void setLagDays(int i) 
   {
     nlag = i;
     if(strategy_computed)
     {recomputeStrategy();}     
   }  
    
   public void setSizePortfolio(int i)
   {
     M = i;
     if(strategy_computed)
     {recomputeStrategy();} 
   }
    
    
   @SuppressWarnings({ "rawtypes", "unchecked" })
public void setupStrategyPanel()
   {
   
        strategyConfigPanel = new JPanel();    
   
        mdfaSettingsPanel = new javax.swing.JPanel();
        insampStartCombo = new javax.swing.JComboBox();
        tradingInSampLabel = new javax.swing.JLabel();
        tradingStartLabel = new javax.swing.JLabel();
        tradingStartCombo = new javax.swing.JComboBox();
        tradingEndCombo = new javax.swing.JComboBox();
        tradingEndLabel = new javax.swing.JLabel();
        nobsLabel = new javax.swing.JLabel();
        nobsBar = new javax.swing.JScrollBar();
        nobsText = new javax.swing.JTextField();
        tradeRuleLabel = new javax.swing.JLabel();
        tradeRuleCombo = new javax.swing.JComboBox();
        paramFileLabel = new javax.swing.JLabel();
        paramFileButton = new javax.swing.JButton();
        paramFileText = new javax.swing.JTextField();
        multivarBox = new javax.swing.JCheckBox();
        stopLossBar = new javax.swing.JScrollBar();
        stopLossText = new javax.swing.JTextField();
        stopLossLabel = new javax.swing.JLabel();
        i1stratCheck = new javax.swing.JCheckBox();
        i2stratCheck = new javax.swing.JCheckBox();
        timeShiftText = new javax.swing.JTextField();
        timeShiftLabel = new javax.swing.JLabel();
        timeShiftScroll = new javax.swing.JScrollBar();
        hybridText = new javax.swing.JTextField();
        hybridLabel = new javax.swing.JLabel();
        hybridScroll = new javax.swing.JScrollBar();
        addStrategyButt = new javax.swing.JButton();
        mdfaPanelParamsBox = new javax.swing.JCheckBox();
        deleteStrategyButt = new javax.swing.JButton();

        uploadFiltersButton = new JButton("Upload From File");
        saveFiltersButton = new JButton("Save Filters"); 
        
        nobsBar.setMaximum(550);
        nobsBar.setMinimum(150);
        nobsBar.setOrientation(JScrollBar.HORIZONTAL);
        nobsBar.setValue(400);        
        
        
        mdfaSettingsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "MDFA Trading Parameters", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, null, new java.awt.Color(0, 255, 243)));

        insampStartCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        insampStartCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "01:30", "01:45", "02:00", "02:15", "02:30", "02:45", "03:00", "03:15", "03:30", "03:45", "04:00", "04:15", "04:30", "04:45", "05:00", "05:15", "05:30",  "05:45",  "06:00",  "06:15",  "06:30",  "06:45", "07:00", "07:15", "07:30", "07:45", "08:00", "08:15", "08:30", "08:45", "09:00", "09:15", "09:30" }));

        tradingInSampLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradingInSampLabel.setText("In-sample Start Time");

        tradingStartLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradingStartLabel.setText("Trading Start Time");

        tradingStartCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradingStartCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "01:30", "01:45", "02:00", "02:15", "02:30", "02:45", "03:00", "03:15", "03:30", "03:45", "04:00", "04:15", "03:30", "03:45", "04:00", "04:15", "03:30", "03:45", "04:00", "04:15", "04:30", "04:45", "05:00", "05:15", "05:30",  "05:45",  "06:00",  "06:15",  "06:30",  "06:45", "07:00", "07:15", "07:30", "07:45", "08:00", "08:15", "08:30", "08:45", "09:00", "09:15", "09:30", "09:45", "10:00", "10:15", "10:30", "10:45", "11:00", "11:15", "11:30", "11:45", "12:00" }));

        tradingEndCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradingEndCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] {"10:00", "10:15", "10:30", "10:45","11:00", "11:15", "11:30", "11:45","12:00", "12:15", "12:30", "12:45", "13:00", "13:15", "13:30", "13:45", "14:00", "14:15", "14:30", "14:45","15:00", "15:15", "15:30", "15:45", "16:00", "16:15", "16:30", "16:45","17:00", "17:15", "17:30", "17:45", "18:00", "18:15", "18:30", "18:45", "19:00", "19:15", "19:30", "19:45", "20:00", "20:15", "20:30", "20:45", "21:00", "21:15", "21:30", "21:45", "22:00", "22:15","22:30", "22:45", "23:00", "23:15", "23:30" }));

        tradingEndLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradingEndLabel.setText("Trading End Time");

        nobsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nobsLabel.setText("Observations");

        nobsBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nobsBar.setMaximum(500);
        nobsBar.setMinimum(150);
        nobsBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        nobsBar.setValue(400);

        nobsText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        nobsText.setText("400");

        tradeRuleLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradeRuleLabel.setText("Trading Rule");

        tradeRuleCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        tradeRuleCombo.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Binary", "Stop-Loss Binary", "Signal Strength" }));

        paramFileLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        paramFileLabel.setText("Parameter FIle");

        paramFileButton.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        paramFileButton.setText("Choose File");

        multivarBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        multivarBox.setSelected(true);
        multivarBox.setText("Multivariate");

        stopLossBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18
        stopLossBar.setMaximum(211);
        stopLossBar.setMinimum(1);
        stopLossBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        stopLossBar.setValue(10);        
        
        
        
        stopLossText.setColumns(5);
        stopLossText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        stopLossText.setText("0");

        stopLossLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        stopLossLabel.setText("Stop-Loss");

        i1stratCheck.setText("i1");

        i2stratCheck.setText("i2");

        timeShiftText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        timeShiftText.setText("10");

        timeShiftLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        timeShiftLabel.setText("Time Shift");

        timeShiftScroll = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,400);
        timeShiftScroll.setValue(200);
        timeShiftScroll.setUnitIncrement(1);
        timeShiftText = new JTextField(3); timeShiftText.setColumns(3);
        timeShiftText.setText("" + 0.0);        
        
        
        hybridText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        hybridText.setText("10");

        hybridLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        hybridLabel.setText("Hybrid");

        hybridScroll.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        hybridScroll.setMaximum(20);
        hybridScroll.setMinimum(1);
        hybridScroll.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        hybridScroll.setValue(10);

    
       hybridScroll = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,500);
       hybridScroll.setValue(0);
       hybridScroll.setUnitIncrement(1);
       hybridText = new JTextField(3);
       hybridText.setText("" + 0);       
       


        
        addStrategyButt.setText("Add Strategy");

        mdfaPanelParamsBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        mdfaPanelParamsBox.setText("Use MDFA Panel Parameters");
        mdfaPanelParamsBox.setSelected(false);
        
        deleteStrategyButt.setText("Delete Strategy");

        javax.swing.GroupLayout mdfaSettingsPanelLayout = new javax.swing.GroupLayout(mdfaSettingsPanel);
        mdfaSettingsPanel.setLayout(mdfaSettingsPanelLayout);
        mdfaSettingsPanelLayout.setHorizontalGroup(
            mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(nobsLabel)
                                .addGap(6, 6, 6)
                                .addComponent(nobsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(nobsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addComponent(paramFileText, javax.swing.GroupLayout.PREFERRED_SIZE, 142, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGap(52, 52, 52))
                                .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addComponent(paramFileLabel)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(paramFileButton))
                                .addComponent(multivarBox)
                                .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(tradingInSampLabel)
                                        .addComponent(tradingStartLabel)
                                        .addComponent(tradingEndLabel))
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(tradingStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(insampStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(tradingEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                        .addGap(18, 18, 18)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(tradeRuleLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(tradeRuleCombo, javax.swing.GroupLayout.PREFERRED_SIZE, 108, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(stopLossLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(stopLossBar, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(stopLossText, javax.swing.GroupLayout.PREFERRED_SIZE, 24, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                        .addComponent(i1stratCheck)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(i2stratCheck))
                                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(timeShiftLabel)
                                            .addComponent(hybridLabel))
                                        .addGap(6, 6, 6)
                                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(hybridScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(timeShiftScroll, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE))))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(hybridText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(timeShiftText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addComponent(mdfaPanelParamsBox))
                        .addContainerGap(38, Short.MAX_VALUE))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(addStrategyButt, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(deleteStrategyButt)
                        .addGap(19, 19, 19)
                        .addComponent(saveFiltersButton)
                        .addGap(19, 19, 19)
                        .addComponent(uploadFiltersButton)
                        .addGap(19, 19, 19))
                        
           
        )));
        mdfaSettingsPanelLayout.setVerticalGroup(
            mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(mdfaPanelParamsBox)
                    .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(nobsLabel)
                                .addComponent(nobsBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGap(4, 4, 4))
                        .addComponent(nobsText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(i1stratCheck)
                            .addComponent(i2stratCheck))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(timeShiftText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(hybridText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(timeShiftScroll, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(timeShiftLabel))
                                .addGap(10, 10, 10)
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                        .addComponent(hybridLabel)
                                        .addGap(2, 2, 2))
                                    .addComponent(hybridScroll, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(tradeRuleLabel)
                            .addComponent(tradeRuleCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(6, 6, 6)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(stopLossText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(stopLossLabel)
                                .addComponent(stopLossBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(multivarBox)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(insampStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(tradingInSampLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(tradingStartLabel)
                            .addComponent(tradingStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(tradingEndLabel)
                            .addComponent(tradingEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(paramFileLabel)
                            .addComponent(paramFileButton))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramFileText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(19, 19, 19)
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(deleteStrategyButt)
                    .addComponent(addStrategyButt)
                    .addComponent(saveFiltersButton)
                    .addComponent(uploadFiltersButton)))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(strategyConfigPanel);
        strategyConfigPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(mdfaSettingsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(mdfaSettingsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    
    
       nobsBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setNobs(((JScrollBar)e.getSource()).getValue());
               nobsText.setText(""+n_obs);
            }
         });
         
       stopLossBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setStopLoss((((JScrollBar)e.getSource()).getValue())*.0001);
               stopLossText.setText(""+((JScrollBar)e.getSource()).getValue());
            }
         });  
         

       nobsText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_nobs(nobsText.getText());}} );  
       
       
       paramFileText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_parameterFile(paramFileText.getText());}} );
              
    
       hybridScroll.addAdjustmentListener(new AdjustmentListener()  {
          public void adjustmentValueChanged(AdjustmentEvent e) 
          {
           hybrid_weight = hybridScroll.getValue()*.01;
           hybridText.setText("" + hybrid_weight);
       }}); 
       
       timeShiftScroll.addAdjustmentListener(new AdjustmentListener() {
          public void adjustmentValueChanged(AdjustmentEvent e)
          {
             time_shift = (timeShiftScroll.getValue()-200)*.1;
             timeShiftText.setText(df.format(time_shift));
          }   
       });    
       
    
       ActionListener buttonActionListener = new ActionListener() 
       {
        File file; int val;
        public void actionPerformed(ActionEvent event)
        {
          if(event.getSource() == paramFileButton)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               setFilterFile(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }  
          else if(event.getSource() == addStrategyButt)
          {enterNewStrategy();}
          else if(event.getSource() == deleteStrategyButt)
          {deleteLatestStrategy();}
          else if(event.getSource() == uploadFiltersButton)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               uploadFromFilterFile(fc.getSelectedFile());    
             }
             else {System.out.println("Open command cancelled by user.");}                 
          } 
          else if(event.getSource() == saveFiltersButton)
          {printStrategyParameters();}
        }
       };

       uploadFiltersButton.addActionListener(buttonActionListener);
       paramFileButton.addActionListener(buttonActionListener);
       addStrategyButt.addActionListener(buttonActionListener);
       deleteStrategyButt.addActionListener(buttonActionListener);
       saveFiltersButton.addActionListener(buttonActionListener);
  }                    

  public void test_metaparameterFile(String s)
  {
     File file = new File(s);
     setPortfolioFile(file);
  }
  
  public void test_spparameterFile(String s)
  {
     File file = new File(s);
     setSPPortfolioFile(file);
  }  
  
  public void test_manualDataFile(String s)
  {
   
    portfolio = s.split("[ ]+"); 
    n_assets = portfolio.length;
    data_files = new String[n_assets];
    fileNameScroll1.setMaximum(n_assets);
    fileNameScroll1.setValue(0);  
    portfolio_file_set = true;
     
    sizePortfolioScroll.setValue(1); 
    sizePortfolioScroll.setMaximum(n_assets+1);
     
  }
  
  public void test_time(String time)
  {
    tradeTime24 = new String(time);
    System.out.println("Daily trade time set to " + tradeTime24);
  }
    
  public void test_manualExpDataFile(String s)
  {
    pairs_trading = false;
    String[] files = s.split("[ ]+");
    expvar = new ArrayList<String>();
    
    for(int i = 0; i < files.length; i++)
    {expvar.add(files[i]);}
    
    exp_var = true;
    if(files.length > 0)
    System.out.println("Explanatory data files added " + s);
    if(files.length == 1)
    {System.out.println("Additional file " + s + " set for pairs trading"); pairs_trading = true;}
    
  }
  
  public void test_parameterFile(String s)
  {
     File file = new File(s);
     setFilterFile(file);  
  }


   
  public void test_nobs(String s)
  {   
       int i;
       try{i = Integer.parseInt(s.trim()); if(i >=150 && i < 500) {nobsBar.setValue(i); n_obs = i; nobsText.setText(""+n_obs);}
                                          else nobsText.setText(""+n_obs);}
       catch (NumberFormatException nfe)
       {System.out.println("NumberFormatException: " + nfe.getMessage()); nobsText.setText(""+n_obs);}
  }  
   
  public void setNobs(int n) {n_obs = n;} 
  public void setStopLoss(double s) {stop_loss = s;}  
  public void numberOfThreads(int th) {n_threads = th;}
   
   
  /*-----------------------------------------------------------
    
    Takes either the filter settings from the imdfaPanel, 
    or uses a filter parameter file and settings given in panel
    
  ------------------------------------------------------------*/
  public void enterNewStrategy()
  {
    
    if(mdfaPanelParamsBox.isSelected()) //get parameters from imdfaPanel
    {
      L = mdfa.L; 
      cutoff = mdfa.cutoff;
      i1 =  mdfa.i1; i2 =  mdfa.i2;
      smooth = mdfa.smooth; decay = mdfa.decay; decay2 = mdfa.decay2; cross = mdfa.cross;
      lambda = mdfa.lambda; expweight = mdfa.expweight; 
      time_shift = mdfa.shift; hybrid_weight = mdfa.onestep; hybrid_weight_diff = mdfa.onestep_diff;
      lag = mdfa.Lag; b0trend = 0; sig_diff = 0; sig_inv = 0; stop_loss_int = 0;
      
      if(turnOnH0) {b0trend = 1;}   
      params_set = true;
      
      if(mdfa.diffSigCheck.isSelected()) {sig_diff = 1;}
      
    }
    else //get from file and panel 
    {
      if(i1stratCheck.isSelected()) {i1 = 1;} else {i1 = 0;}
      if(i2stratCheck.isSelected()) {i2 = 1;} else {i2 = 0;}
      sig_inv = 0; stop_loss_int = 0;
      if(filter_file_set) {readFilterParamFile(); params_set = true;}
      else {System.out.println("Filter file not yet selected"); params_set = false;}
    }
    
    if(params_set)
    {
     StrategyParameters sp = new StrategyParameters("Filter_"+strategies.size(), (String)insampStartCombo.getSelectedItem(), 
                          (String)tradingStartCombo.getSelectedItem(), tradingStartCombo.getSelectedIndex(), (String)tradingEndCombo.getSelectedItem(), 
                          nobsBar.getValue(), 2, L, lag, cutoff, lambda, expweight, smooth, decay, decay2, cross,
                          i1, i2, time_shift, hybrid_weight,hybrid_weight_diff,b0trend,sig_diff, sig_inv, stop_loss_int);    
  
     strategies.add(sp);
    }
    //Repaint - Post to terminal 
    stratCanvas.addTextNewLine(strategies.get(strategies.size()-1).toString(), Color.GREEN);
    
  }
  
  public void deleteLatestStrategy()
  {
   if(strategies.size() > 0)
   {
    strategies.remove(strategies.size()-1);
    //--- Repaint    
    stratCanvas.clearText();
    
    for(int i = 0; i < strategies.size(); i++)
    {
     stratCanvas.addTextNewLine(strategies.get(i).toString(), Color.GREEN);
    }
   } 
  }
  
  public void displayFiles()
  {
    if(data_files.length > 0)
    {
      for(int i = 0; i<data_files.length; i++)
      {stratCanvas.addSameLineText(data_files[i],Color.BLUE);}
    }  
  }
  
  
  public void setFilterFile(File file)
  {filter_file = file; filter_file_set = true; paramFileText.setText(filter_file.getName());}    
   
   
  public void readFilterParamFile()
  {
         String strline; 
         //file_valid = true;
         double[] parameters = new double[10];
         Double D; Integer I;
         //default parameters
         parameters[1] = 0.52; parameters[8] = 28;
         String[] tokens; String delims = "[ ]+";
         int n_toks;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(filter_file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           strline = br.readLine();
           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks <= 2)
           {
             D = new Double(tokens[0]); cutoff0 = D.doubleValue();
             D = new Double(tokens[1]); cutoff = D.doubleValue();                 
           }
          
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); lambda = D.doubleValue();

           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); expweight = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); smooth = D.doubleValue();           
    
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); decay = D.doubleValue();    
   
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); decay2 = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           D = new Double(tokens[0]); cross = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           I = new Integer(tokens[0]); L = I.intValue();           
           
           strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
           I = new Integer(tokens[0]); lag = I.intValue();

           br.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}
 
  }     
   
   
  public void setPortfolioFile(File f)
  {
     String strline;
     ArrayList<String> names = new ArrayList<String>();
     try
     {
          
       FileInputStream fin = new FileInputStream(f);
       DataInputStream din = new DataInputStream(fin);
       BufferedReader br = new BufferedReader(new InputStreamReader(din));
  
       while((strline = br.readLine()) != null)    
       {names.add(strline);}
       
       br.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
     catch(IOException ioe){System.out.println("IO out error..." + ioe);}     
  
     portfolio = names.toArray(new String[0]);
     n_assets = portfolio.length;
     data_files = new String[n_assets];
     fileNameScroll1.setMaximum(n_assets);
     fileNameScroll1.setValue(0);  
     portfolio_file_set = true;
     
     sizePortfolioScroll.setValue(1); 
     sizePortfolioScroll.setMaximum(n_assets+1);
     
  }
   
   
   
  public void setSPPortfolioFile(File f)
  {
     String strline; 
     ArrayList<String> names = new ArrayList<String>();
     try
     {
          
       FileInputStream fin = new FileInputStream(f);
       DataInputStream din = new DataInputStream(fin);
       BufferedReader br = new BufferedReader(new InputStreamReader(din));
  
       while((strline = br.readLine()) != null)    
       {names.add(strline);}
       
       br.close();
     }
     catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
     catch(IOException ioe){System.out.println("IO out error..." + ioe);}     
  
     portfolio = names.toArray(new String[0]);
     n_assets = portfolio.length;
  
     data_files = new String[n_assets];
     System.arraycopy(portfolio,0,data_files,0,n_assets);
  
     fileNameScroll.setMaximum(n_assets);
     fileNameScroll.setValue(0);
     sp_files_set = true;
        
     sizePortfolioScroll.setValue(1); 
     sizePortfolioScroll.setMaximum(n_assets+1);     
     
  }   
   
   
   
   
   void createPortfolioFile()
   {
   
     int i;
     File file = new File("Portfolio.port");
     String last_trading_day;
     String first_trading_day;
     //--- Get latest date string first------------------------
     DateTime today = new DateTime();

     if(today.getHourOfDay() < 16) 
     {today = today.minus(Days.days(1));}
     
     DateTime four_months = today.minusYears(noYearsCombo.getSelectedIndex()+1);
     
     if(today.dayOfWeek().getAsText().equals("Sunday")) {today = today.minusDays(2);}
     else if(today.dayOfWeek().getAsText().equals("Saturday")) {today = today.minusDays(1);}
     
     int iday = today.getDayOfMonth();
     int iyear = today.getYear();
     int imonth = today.getMonthOfYear();
     
     if(iday < 10 && imonth < 10) 
     {last_trading_day = ""+iyear+"0"+imonth+"0"+iday;}
     else if(iday < 10 && imonth >= 10)
     {last_trading_day = ""+iyear+""+imonth+"0"+iday;}
     else if(iday >= 10 && imonth < 10)
     {last_trading_day = ""+iyear+"0"+imonth+""+iday;}     
     else
     {last_trading_day = ""+iyear+""+imonth+""+iday;}
     
     iday = four_months.getDayOfMonth();
     iyear = four_months.getYear();
     imonth = four_months.getMonthOfYear();
     
    
     
     if(iday < 10 && imonth < 10) 
     {first_trading_day = ""+iyear+"0"+imonth+"0"+iday;}
     else if(iday < 10 && imonth >= 10)
     {first_trading_day = ""+iyear+""+imonth+"0"+iday;}
     else if(iday >= 10 && imonth < 10)
     {first_trading_day = ""+iyear+"0"+imonth+""+iday;}     
     else
     {first_trading_day = ""+iyear+""+imonth+""+iday;}     
     
     
     if(portfolio_file_set){
     //----------------------------------------------------------
     try
     {
      PrintWriter out = new PrintWriter(new FileWriter(file));
     
      for(i=0;i<n_assets-1;i++)
      {out.print(portfolio[i] + " ");}
      out.println(portfolio[n_assets-1]); 
     
      
     
      
      //starting date/time  
      if(forexData.isSelected())
      {
       if(data5.isSelected()) //24-hour 5 minute
       {
//         out.println(first_trading_day + " 012500"); 
//         out.println(last_trading_day + " 231000");
//         out.println("012500 231000"); 
//         out.println("01:30");
//         out.println("5");
//         out.close();
        
        out.println(first_trading_day + " 000000"); 
        out.println(last_trading_day + " 235500");
        out.println("000000 235500"); 
        out.println("00:00:00");
        out.println("5");
        out.close();
        
        System.out.println("Data5 is selected");
       }
       else if(data30.isSelected() && !data60.isSelected()) //24-hour 30 minute
       {
/*        out.println(first_trading_day + " 013000"); 
        out.println(last_trading_day + " 233000");
        out.println("010000 230000"); 
        out.println("01:30");
        out.println("30");
        out.close();  */    
        
        out.println(first_trading_day + " 000000"); 
        out.println(last_trading_day + " 233000");
        out.println("000000 233000"); 
        out.println("00:00:00");
        out.println("30");
        out.close();           
        
        
        
        System.out.println("Data30 is selected");
       }
       else if(data60.isSelected() && !data30.isSelected()) //24-hour 60 minute
       {
        out.println(first_trading_day + " 000000"); 
        out.println(last_trading_day + " 230000");
        out.println("000000 230000"); 
        out.println("00:00:00");
        out.println("60");
        out.close();
        System.out.println("Data60 is selected");
       }
       else if(data60.isSelected() && data30.isSelected())
       {
               
        out.println(first_trading_day + " 003000"); 
        out.println(last_trading_day + " 233000");
        out.println("003000 233000"); 
        out.println("00:30:00");
        out.println("60");
        out.close(); 
       
       }
       else  //default is 24-hour 15 minute
       {
//         out.println(first_trading_day + " 011500");        
//         out.println(last_trading_day + " 231000");
//         out.println("011500 231000"); 
//         out.println("01:30");
//         out.println("15");
//         out.close();       
        
        out.println(first_trading_day + " 000000");        
        out.println(last_trading_day + " 234500");
        out.println("000000 234500"); 
        out.println("00:00:00");
        out.println("15");
        out.close();            
        
        
        
        System.out.println("Data15 is selected");
       }
     }  
     else if(futuresData.isSelected())
     {
       if(data5.isSelected()) //24-hour 5 minute
       {
        out.println(first_trading_day + " 032500"); 
        out.println(last_trading_day + " 161000");
        out.println("032500 161000"); 
        out.println("03:30");
        out.println("5");
        out.close();
       }     
       else if(data30.isSelected()) //24-hour 30 minute
       {
        out.println(first_trading_day + " 033000"); 
        out.println(last_trading_day + " 160000");
        out.println("033000 160000"); 
        out.println("03:30");
        out.println("30");
        out.close();                 
       }
       else if(data60.isSelected()) //24-hour 60 minute
       {
        out.println(first_trading_day + " 040000"); 
        out.println(last_trading_day + " 160000");
        out.println("040000 160000"); 
        out.println("04:00");
        out.println("60");
        out.close();
       }     
       else
       {
         out.println(first_trading_day + " 031500"); 
         out.println(last_trading_day + " 161000");
         out.println("031500 161000"); 
         out.println("03:30");
         out.println("15");
         out.close();
       }  
     }
     else //default equity/pairs data 8h30 to 16h15
     {
     
     
       if(data5.isSelected()) //24-hour 5 minute
       {
        out.println(first_trading_day + " 082500"); 
        out.println(last_trading_day + " 161000");
        out.println("082500 161000"); 
        out.println("08:30");
        out.println("5");
        out.close();
       }     
       else if(data30.isSelected()) 
       {
                
        out.println(first_trading_day + " 080000"); 
        out.println(last_trading_day + " 153000");
        out.println("080000 153000"); 
        out.println("08:30");
        out.println("30");
        out.close();    
        
       }
       else if(data60.isSelected()) //24-hour 60 minute
       {
        out.println(first_trading_day + " 090000"); 
        out.println(last_trading_day + " 160000");
        out.println("090000 160000"); 
        out.println("09:00");
        out.println("60");
        out.close();
       }       
       else
       {
        out.println(first_trading_day + " 081500"); 
        //out.println("20120606 082500"); 
        out.println(last_trading_day + " 161000");
        out.println("081500 161000"); 
        out.println("08:30");
        out.println("15");
        out.close();
       }
     }
      
      
    }
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
    }
   }
   
   
   public void downloadDataIQFeed()
   {
      
     createPortfolioFile();
      
     if(data5.isSelected())
     {data.min_60 = false; data.min_30 = false; data.min_5 = true;} 
     else if(data30.isSelected())
     {data.min_30 = true; data.min_60 = false; data.min_5 = false;}
     else if(data60.isSelected())
     {data.min_60 = true; data.min_30 = false; data.min_5 = false;}
     
     if(forexData.isSelected()) {data.est = false;}
     else {data.est = true;}
     
     if(data24.isSelected()) {data.daily = true; data.dailyTime = new String(tradeTime24);}
     
     data.zero_marketopen = zeroButton.isSelected();
     data.setupHistoricalDownload();
     if(forexData.isSelected())
     {data.setLogTrans(false);}
         
     for(int m = 0; m < n_assets; m++)
     {
      if(data.daily)
      { } //data.downloadHistoricalIQData24(portfolio[m]);}
      else if(forexData.isSelected())
      {
       if(data30.isSelected() && data60.isSelected())
       {
        data.downloadHistoricalIQDataFXSpreads30(portfolio[m]);
        //data.downloadHistoricalIQDataIBFX30(portfolio[m]);
        //mdfaStrat.ib_data = true; 
       }
       else 
       {
        data.downloadHistoricalIQDataFXSpreads(portfolio[m]);
        //data.downloadHistoricalIQDataIBFX(portfolio[m]);
        //mdfaStrat.ib_data = true;
       }         
      }       
      else
      {data.downloadHistoricalIQData(portfolio[m]);}  
      String[] names = portfolio[m].split("[ ]+");
      data_files[m] = new String(names[0] + ".dat");     
     } 
     historicalData_set = true;

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
    
  
  public static void sort_sims(int n, double[] object, int[] indices)
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
       
   
   
   
    class startIQConnect extends Thread
    {
    
     public void run()
     {
        System.out.println("start up IQ");
		try {
                       //
			// Launch IQFeed and Register the app with IQFeed.
// 			System.out.println("Launching IQConnect.");
// 			//Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe ‑product SIMON_OTZIGER_6345 -version 2.9.0.13 ‑login 422113 ‑password 76059316");
// 			System.out.println("Verifying if IQConnect is connected to the server");
// 			Thread.sleep(12000);
// 			// verify everything is ready to send commands.
// 			boolean bConnected = false;
// 			// connect to the admin port.
// 			Socket sockAdmin = new Socket(InetAddress.getByName("192.168.56.10"), 8300);
			
			System.out.println("Launching IQConnect.");
	
                        //Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe -product JEANPIERRE_ROBALO_11529 -version 2.9.0.13 ‑login 431660 ‑password 70433790");
                        Runtime.getRuntime().exec("/usr/bin/wine IQConnect.exe -product PRABIN_SETH_11606 -version 2.9.0.13 ‑login 438106 ‑password 64568341");
			System.out.println("Verifying if IQConnect is connected to the server");
			Thread.sleep(12000);
			//verify everything is ready to send commands.
			boolean bConnected = false;
			//connect to the admin port.
			Socket sockAdmin = new Socket(InetAddress.getByName("localhost"), 9300);			
			
			
			BufferedReader bufreadAdmin = new BufferedReader(new InputStreamReader(sockAdmin.getInputStream()));
			BufferedWriter bufwriteAdmin = new BufferedWriter(new OutputStreamWriter(sockAdmin.getOutputStream()));
			String sAdminLine = "";
			// loop while we are still connected to the admin port or until we are connected
			while (((sAdminLine = bufreadAdmin.readLine()) != null) && !bConnected)
			{
				System.out.println(sAdminLine);
				if (sAdminLine.indexOf(",Connected,") > -1)
				{
					System.out.println("IQConnect is connected to the server.");
					bConnected = true;
				}
				else if (sAdminLine.indexOf(",Not Connected,") > -1)
				{
					System.out.println("IQConnect is Not Connected.\r\nSending connect command.");
					bufwriteAdmin.write("S,CONNECT\r\n");
					bufwriteAdmin.flush();
				}
			}
			// cleanup admin port connection
			sockAdmin.shutdownOutput();
			sockAdmin.shutdownInput();
			sockAdmin.close();
			bufreadAdmin.close();
			bufwriteAdmin.close();

			// at this point, we are connected and the feed is ready.
			int total_num_assets = 8;
			
// 			for(int i = 0; i < n_basket; i++)
// 			{total_num_assets = total_num_assets + ppms[i].name.length;}
// 			
// 			total_num_assets = ppms.length*3; 
			s_outs = new BufferedWriter[total_num_assets];
			s_ins = new BufferedReader[total_num_assets];
			Socket[] socks = new Socket[total_num_assets];
			
			for(int i = 0; i < total_num_assets; i++)
                       {
                         socks[i] = new Socket(InetAddress.getByName("localhost"), 9100);
	                 //socks[i] = new Socket(InetAddress.getByName("192.168.56.10"), 8100);
	                 s_ins[i] = new BufferedReader(new InputStreamReader(socks[i].getInputStream()));
	                 s_outs[i] = new BufferedWriter(new OutputStreamWriter(socks[i].getOutputStream()));
	                 // Set the lookup port to protocol 5.0 to allow for millisecond times, 
                         // market center, trade conditions, etc
                         s_outs[i].write("S,SET PROTOCOL,5.0\r\n");
                         s_outs[i].flush();   
	               }
                        
                     }
                     catch (Exception e)
		     {e.printStackTrace();}    
    
    
      } 
   }


	public void setMDFAParameters(IMDFAPanel mdfa, boolean h0) {
		
		this.mdfa = mdfa;	
		turnOnH0 = h0;
	}          
   
     
}

