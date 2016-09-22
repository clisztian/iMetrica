package ch.imetrica.mdfaTradingStrategies;

import java.io.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.JCheckBox;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.BorderFactory; 
import java.text.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/*----------------------------------------------------------------
 
   Panel for IMDFA Controls - 
     Simulation:  
        Two different types of simulations 
          - DFA uMath.sing uSimX13 data 
          - MDFA uMath.sing simulated GDP series of Wildi 
   

----------------------------------------------------------------*/


public class MDFAStrategyPanel extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	mdfaStrategyCanvas mdfaAnalysisCanvas;      
    private Component parent;
    JButton computeButton,deleteHistButton,deleteVarRatioButton,historicalFileButton;
    JLabel historicalFileLabel;
    JTextField historicalFileText;
    @SuppressWarnings("rawtypes")
	JComboBox insampStartCombo;
    JPanel plotPanel;
    JScrollBar kdiffBar;
    JLabel kdiffLabel;
    JTextField kdiffText;
    JPanel mdfaSettingsPanel;
    JLabel mdfaStratLabel;
    JScrollBar ncoresBar;
    JLabel ncoresLabel;
    JTextField ncoresText;
    JScrollBar nobsBar, samplesBar;
    JLabel nobsLabel,sharpeLabel;
    JTextField nobsText;
    JButton paramFileButton,paramFileButton1;
    JLabel paramFileLabel,paramFileLabel1;
    JTextField paramFileText, paramFileText1, samplesText;
    JCheckBox[] plotHist,plotVar;
    JCheckBox plotStock;
    JProgressBar stratProgressBar;
    @SuppressWarnings("rawtypes")
	JComboBox tradeRuleCombo,varRatioEndCombo,tradingEndCombo,tradingStartCombo;
    JLabel tradeRuleLabel,tradingEndLabel,tradingInSampLabel;
    JLabel tradingStartLabel;
    JButton varRatioButton;
    JLabel varRatioEndLabel;
    JPanel varRatioPanel;
    JLabel varRatioPanelLabel;
    JProgressBar varRatioProgressBar;
    @SuppressWarnings("rawtypes")
	JComboBox varRatioStartCombo;
    JLabel varRatioStartLabel, windowSizeLabel;
    JLabel avgRankLabel,maxDrawLabel,minRankLabel,sharpeRatioText;
    JTextField avgRankText,maxDrawText,minRankText,sharpeText; 
    DecimalFormat df3,df;
    JFileChooser fc; 
    String curDir;
    JCheckBox multivarBox;
    JCheckBox autoVar; 
    boolean auto; 
    int max_series; 
    JTextField stopLossText;
    JScrollBar stopLossBar;
    JLabel stopLossLabel;
    
    JCheckBox[] rollingCheck;
    
    boolean readyTradeAnalysis = false;
    JScrollBar recompDayBar;
    JLabel recompDayLabel;
    JTextField recompDayText;
    String final_time = "16:15";
    double tradingCost = 0; 
    ArrayList<StrategyParameters> strategies; 
    // MDFA strategy selection
    VarRatio vr; int insamp_start_int;
    ArrayList<String> total_perf; 
    ArrayList<Double> rets; 
    String[] dates;
    double[] strat_returns; 
    String endingTime; 
    String[] startingHeures = new String[5];
    String insamp_startingHeure; 
    int start_seq; int nobs; 
    boolean hist_set = false; 
    boolean filter_file_set = false;
    File hist_data; File filter_file; 
    String[] filter_files;
    int n_threads;
    ArrayList<double[]> return_collection; 
    int n_sections,n_trading_days; 
    int length_rank = 15;
    double sharpe_ratio;
    double port_mean;
    double max_drawdown;
    double total_return;
    double min_rank;
    double bRatio;        
    int Kdiff, n_samples;
    double rank_coeff; 
    double avg_rank;
    ArrayList<Double> returns_all;
    ArrayList<Double> stock_perf;
    ArrayList<MDFATrade> trade_perf;
    boolean multivar,meta_mode;
    double stop_loss; 
    int minute_df = 15;
    int recompute_day = 0;
    JCheckBox shortBox; 
    JCheckBox firm_closeBox; 
    boolean short_sell = true;
    boolean ib_data = false;
    String ib_filename;
    ArrayList<String> rolling_ind;
    int rolling_length = 30;
    JDialog tradeAnalysisDialog;
    MDFAStrategyTradeCanvas tradeAnalysis;
    JButton launchTradeAnalysis;
    MDFATrade[] tradesAll_adj;
    double stop_loss_thresh = 1;
    double take_profit = 0;
    JPanel tradeAnalysisPanel;
    static boolean dailyStrat = false;
    private JLabel adjustSLLabel;
    private JTextField adjustSLText;
    private JScrollBar adjustSLBar;
    private JScrollBar takeProfitBar;
    private JLabel takeProfitLabel;
    private JTextField takeProfitText;     
    
    
    
    
  public MDFAStrategyPanel(Component _p, JFrame f)
  {
  
  
     //--- initialize some defaults
     nobs = 400; parent = _p;
     insamp_startingHeure = "08:30"; start_seq = 4; 
     startingHeures[0] = "08:30"; startingHeures[1] = "08:45";
     startingHeures[2] = "09:00"; startingHeures[3] = "09:15";
     startingHeures[4] = "09:30"; 
     n_threads = 10; 
     return_collection = new ArrayList<double[]>(); 
     Kdiff = 5;
     n_samples = 20; insamp_start_int = 0;
     df = new DecimalFormat("##.##"); 
     df3 = new DecimalFormat("##.###");
     multivar = true;
     curDir=System.getProperty("user.dir") + File.separator;
     fc = new JFileChooser(curDir);
     fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
     meta_mode = false;      
     
     tradeAnalysis = new MDFAStrategyTradeCanvas(1200,800);
     initTradeComponents();
     
     tradeAnalysisDialog = new JDialog(f,"Trade Analysis",true);
     tradeAnalysisDialog.getContentPane().add(tradeAnalysisPanel);
     tradeAnalysisDialog.pack();
     tradeAnalysisDialog.setLocationRelativeTo(f);
     tradeAnalysisDialog.setModal(false);
     tradeAnalysisDialog.setVisible(false);     
     
     
     rolling_ind = new ArrayList<String>();
  }   
    
    
  public void setInSampStartingTime(int c)
  {insamp_startingHeure = startingHeures[c]; insamp_start_int = c;}
    
  public void setTradingCost(int basis) {tradingCost = 0.00001*basis;}  
    
  public void setStartingTradeTime(int c)
  {start_seq = c + 4 - insamp_start_int;}
  
  public void setEndingTradeTime(int c)
  {
     if(c == 0) {endingTime = "15:00";}
     else if(c == 1) {endingTime = "15:15";}    
     else if(c == 2) {endingTime = "15:30";}   
     else if(c == 3) {endingTime = "15:45";}        
     else if(c == 4) {endingTime = "16:00";} 
     
  }
  
  public void setNobs(int n) {nobs = n;}
    
  public void setHistoricalDataFile(File hist_data_name)
  {
   hist_data = hist_data_name; hist_set = true; historicalFileText.setText(hist_data.getName());
   
   if(ib_data)
   {
     String[] ibname = hist_data.getName().split("[.]+");
     ib_filename = new String(ibname[0] + ".IB.dat");
   }
  }
    
  public void setFilterFile(File file)
  {
     filter_file = file;  paramFileText.setText(filter_file.getName());
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
                          (new Double(tokens[19])).doubleValue(), (new Integer(tokens[20]).intValue()), (new Integer(tokens[21]).intValue()),
                          (new Integer(tokens[22])).intValue(), (new Integer(tokens[23]).intValue()));
             strategies.add(sp);  
             filter_file_set = true;
             }           
           }
    
           br.close();
   }               
   catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
   catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}      
   
  }    
    
  public void setMetaFilterFile(File file)
  {
  
    
     String strline;
     ArrayList<String> names = new ArrayList<String>();   
       
     System.out.println(file.getParent());  
     file.getParent();
       
       
     try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {
           names.add(strline);              
         }  
     
         din.close(); br.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         
  
    meta_mode = true;
    filter_files = names.toArray(new String[0]);
    filter_file_set = true;
    paramFileText1.setText(file.getName());
  }
    
    
  public void numberOfThreads(int th) {n_threads = th;}  
    
  public void setKDifferences(int k) {Kdiff = k;}
  public void setNSamples(int s) {n_samples = s; rolling_length = s;}
  public void multivarMDFA(boolean f) {multivar = f;} 
  public void setStopLoss(double s) {stop_loss = s; System.out.println(stop_loss);} 
  public void setRecompDay(int r) {recompute_day = r; System.out.println(recompute_day);}
  public void shortSell(boolean f) {short_sell = f; System.out.println("shorting = " + short_sell);}
    
  public void computeVarRatio()
  {
   if(hist_set)
   {
   
     vr = new VarRatio();

     String begin_time = (String)varRatioStartCombo.getSelectedItem();  
     String end_time = (String)varRatioEndCombo.getSelectedItem();
     
     vr.computeVarRatio(hist_data, Kdiff, begin_time, end_time, n_samples);
     
     double[] means = mean_std( vr.varianceRatio.toArray(new Double[0]));
     
     if(auto)
     {
     
      if(mdfaAnalysisCanvas.var_series.size() == 0)
      {
        mdfaAnalysisCanvas.addVarSeries(vr.varianceRatio.toArray(new Double[0]), new String("K=" + Kdiff + ", window=" + n_samples + 
       ", time=("+begin_time+", "+end_time+"), mean="+df3.format(means[0])+",std="+df3.format(means[1])));
      }
      else
      {
      mdfaAnalysisCanvas.setVarSeries(mdfaAnalysisCanvas.var_series.size()-1, vr.varianceRatio.toArray(new Double[0]), new String("K=" + Kdiff + ", window=" + n_samples + 
       ", time=("+begin_time+", "+end_time+"), mean="+df3.format(means[0])+",std="+df3.format(means[1])));
      }
     }  
     else
     {
      mdfaAnalysisCanvas.addVarSeries(vr.varianceRatio.toArray(new Double[0]), new String("K=" + Kdiff + ", window=" + n_samples + 
       ", time=("+begin_time+", "+end_time+"), mean="+df3.format(means[0])+",std="+df3.format(means[1])));
     } 
       

       
       
     for(int i=0;i<mdfaAnalysisCanvas.var_series.size();i++) {plotVar[i].setEnabled(true); plotVar[i].setSelected(true);}       
      
   }
  }
    
    
  public void computeMDFA_EvolutionStrategy() throws InterruptedException, ExecutionException 
  {
  
    if(filter_file_set && hist_set)
    {
     int j;
     String[] datafiles; 
     total_perf = new ArrayList<String>();
     returns_all = new ArrayList<Double>();      
     stock_perf = new ArrayList<Double>();
     trade_perf = new ArrayList<MDFATrade>();
     
     for(int l = 0; l < strategies.size(); l++)
     {
     
      total_perf.clear(); returns_all.clear(); stock_perf.clear(); trade_perf.clear();
      List<MDFAStrategyEvolution> objects = new ArrayList<MDFAStrategyEvolution>(); 
      
      int cusion=20;
      if(!dailyStrat)
      {cusion = computeMaxTradeObs(hist_data); nobs = (strategies.get(l)).n_obs;}
      else
      {      
       final_time = "23:00"; minute_df = 15;
       nobs = (strategies.get(l)).n_obs;
       
      }
      datafiles = MDFAStrategyEvolution.forkDataFile(hist_data, n_threads, nobs, cusion);
  
      for(j = 0; j < datafiles.length-1; j++) //parallelization occurs on spit files 
      { 
            System.out.println(j);
            MDFAStrategyEvolution strategy = new MDFAStrategyEvolution(10, "filter_trend_es_3.params"); 
            if(n_threads == 1) {strategy.togglePrint(true);}

            if(ib_data && ib_filename != null)
            {strategy.ib_data = true; strategy.ib_data_file = new String(ib_filename);}
            
            strategy.setRecomputeDay(recompute_day);
            strategy.setFinalTradeTime("23:00");
            //strategy.min30trade = true;
            if(final_time.indexOf("16:30") != -1)
            {strategy.setFinalTradeTime("23:30");}
            if(final_time.indexOf("17:00") != -1)
            {strategy.setFinalTradeTime("23:00");}            
            else strategy.setFinalTradeTime(final_time);
            strategy.setObservationFrequency(minute_df);
            
            
            
            //System.out.println(strategy.stop_loss_thresh);
            strategy.setStrategyParameters(strategies.get(l), datafiles[j], tradingCost);
            strategy.let_shop = !firm_closeBox.isSelected();
            //System.out.println("stop loss = " + stop_loss);
            if(stop_loss > 0) {strategy.stop_loss = true; strategy.stop_loss_thresh = stop_loss;}
            //else {strategy.stop_loss = false; strategy.stop_loss_thresh = stop_loss;}
            
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
          // some other exectuors you could try to see the different behaviours
          //ExecutorService exec = Executors.newFixedThreadPool(8);
          // ExecutorService exec = Executors.newSingleThreadExecutor();
          try 
          {
           System.currentTimeMillis();
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
              stock_perf.addAll(fr.get().asset);
              trade_perf.addAll(fr.get().mdfaTrades);
             }
           }
          } 
          finally 
          {System.out.println("shutting down exec"); }   
  
          parseAndPlotReturns(); 
          
       }   
    }
    else{System.out.println("Must specify historical data and filter file first");} 
    System.out.println("Done with performance on filter " + filter_file.getName());
  
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
       
     }
     catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
     
     return max_count; 
  }
    
    
    
  public static int minute_diff(String t1, String t2)
  {
     String[] time1 = t1.split("[:]+");
     String[] time2 = t2.split("[:]+");
     
     int min1 = (new Integer(time1[1])).intValue();
     int min2 = (new Integer(time2[1])).intValue();
     
     int min_diff = Math.abs(min1 - min2);
     if(min_diff == 0) min_diff = 60;
     
     return min_diff;
  }  
/*    
  public void computeMDFAStrategy() throws InterruptedException, ExecutionException 
  {  
    if(filter_file_set && hist_set)
    {
     //--- Get set time frames ------------------------------------
     setInSampStartingTime(insampStartCombo.getSelectedIndex());
     setStartingTradeTime(tradingStartCombo.getSelectedIndex());
     setEndingTradeTime(tradingEndCombo.getSelectedIndex());
     
     //create strategy objects

     total_perf = new ArrayList<String>();
     returns_all = new ArrayList<Double>();
     stock_perf = new ArrayList<Double>();
     
     //int obs = 478;
     String[] data_files = mdfaStrategyFork.forkDataFile(hist_data, n_threads, nobs);
       
     List<mdfaStrategyFork> objects = new ArrayList<mdfaStrategyFork>();
     int adaptive_recomp = 81;
     int recomp = 0;
     int close = 0;

  
     for(int j = 0; j < data_files.length-1; j++)
     { 
      
        mdfaStrategyFork strategy = new mdfaStrategyFork(10, "filter_trend_es_3.params"); 
        if(multivar) {strategy.n_rep = 4;} else {strategy.n_rep = 2;}
        strategy.i1 = 0;
        if(tradeRuleCombo.getSelectedIndex() == 1) {strategy.stop_loss = true; strategy.stop_loss_thresh = stop_loss;} 
        else {strategy.stop_loss = false;}
        
        strategy.start_seq = this.start_seq;
        strategy.close_time = close;
        strategy.setTradingParameters(this.short_sell, 0);
        strategy.turnOnVolatilityThresh(false);
        strategy.turnOnGaussianize(false);
        strategy.saveVolatility(true);
        strategy.setAdaptiveFiltering(false);
        strategy.setAdaptiveParameterFile("child_filter.params");
        strategy.setAdaptiveRecompPulse(31);
        strategy.setLagLookback(false,false);
        //strategy.setTradingParameters(true, 0);
        strategy.setStartIndex(0);
        strategy.setRecompPulse(recomp);     
        strategy.setParameterFile(filter_file.getName()); 
        strategy.readFilterParamFile();
        strategy.setMDFAParameters();    
        strategy.dataFiles[0] = data_files[j]; 
        strategy.setNobs(nobs);            
        strategy.setStartingTime(insamp_startingHeure);
        strategy.setEndingTime((String)tradingEndCombo.getSelectedItem()); 
        //System.out.println(strategy.endingTime + " " + strategy.trade_obs);
        objects.add(strategy);  
     }
    
    
     List<Callable<Result>> tasks = new ArrayList<Callable<Result>>();
      
     for (final mdfaStrategyFork object : objects) 
     {
        Callable<Result> c = new Callable<Result>() {
        @Override
          public Result call() throws Exception 
          {return compute(object);}
        };
        tasks.add(c);
     }

     ExecutorService exec = Executors.newCachedThreadPool();
        // some other exectuors you could try to see the different behaviours
        // ExecutorService exec = Executors.newFixedThreadPool(3);
        // ExecutorService exec = Executors.newSingleThreadExecutor();
     try 
     {
       long start = System.currentTimeMillis();
       List<Future<Result>> results = exec.invokeAll(tasks);
       int sum = 0;
       for (Future<Result> fr : results) 
       {
           //sum += fr.get().wait;
//            System.out.println(String.format("Task waited %d ms",
//            fr.get().wait));
          System.out.println(fr.get().perform_stats);
          total_perf.addAll(fr.get().performance);
          returns_all.addAll(fr.get().rets);
          stock_perf.addAll(fr.get().asset);
       }
//             long elapsed = System.currentTimeMillis() - start;
//             System.out.println(String.format("Elapsed time: %d ms", elapsed));
//             System.out.println(String.format("... but compute tasks waited for total of %d ms; speed-up of %.2fx", sum, sum / (elapsed * 1d)));
      } 
      finally 
      {exec.shutdown();}   
      
      
      parseAndPlotReturns(); 
    }
    else{System.out.println("Must specify historical data and filter file first");} 
    System.out.println("Done with performance on filter " + filter_file.getName());
  }    
    */
    
  private static class Result 
  {
        private final String perform_stats;
        private ArrayList<String> performance;
        private ArrayList<Double> rets;
        private ArrayList<Double> asset;
        private ArrayList<MDFATrade> mdfaTrades;
        
        public Result(String stats, ArrayList<String> per, ArrayList<Double> r, ArrayList<Double> a, ArrayList<MDFATrade> trades) 
        {
            this.perform_stats = stats;
            this.performance = per; 
            this.rets = r; 
            this.asset = a; 
            this.mdfaTrades = trades;
        }
  }  
  

  public static Result compute(MDFAStrategyEvolution strategy) throws InterruptedException 
  {
        
     if(!dailyStrat) strategy.startStrategy();
     else strategy.startStrategyDaily();
     
     double ratio_trades = 0; String perf = "";
     if(!dailyStrat) {ratio_trades = strategy.trade_succ_ratio/(double)strategy.returns.size();}

     double sharpe = Math.sqrt(250)*strategy.mean_perf/strategy.standard_deviation;
     if(!dailyStrat) {perf = new String(" n_obs = " + strategy.n_obs + ", avg_rank = " + strategy.rank_coeff + ", maxdraw = " + strategy.maxdraw + ", mean = " + strategy.mean_perf + ", success_ratio = " + ratio_trades + ", sharpe " + sharpe + ", ulcer = " 
     + strategy.ulcer_index);}
     else
     {
       perf = new String(" n_obs = " + strategy.n_obs  + ", mean = " + strategy.mean_perf + ", success_ratio = " + strategy.win_ratio + ", sharpe " + sharpe + ", ulcer = " 
     + strategy.ulcer_index);
     }
     
 
     return new Result(perf,strategy.date_returns,strategy.returns,strategy.dailyoutret, strategy.mdfaTrades);
  }   
    
    
  /*------------------------------------------------------------
    Here we consider all trades will unrealized pnl > stop_loss
    and then set pnl and loss to stop_loss
  -------------------------------------------------------------*/
  
  public void applyStopLoss()
  {
    int i;
    int n_trades = trade_perf.size();
    double pnl = 0;
    ArrayList<String> perf_dates = new ArrayList<String>();
    ArrayList<Double> perf_returns = new ArrayList<Double>();
    
    //stop_loss_tresh given in .00000 format 
    //create a tradesAll array
    if(n_trades > 0 && readyTradeAnalysis)
    {
    
      tradesAll_adj = new MDFATrade[n_trades];    
    
      //adjust all trades larger than stop_loss
      for(i = 0; i < n_trades; i++)
      {
        
        tradesAll_adj[i] = new MDFATrade(trade_perf.get(i).date, trade_perf.get(i).startTransTime, trade_perf.get(i).endTransTime, trade_perf.get(i).drawdown,
                                         trade_perf.get(i).drawup, trade_perf.get(i).pl);
        
                                    
        if(Math.abs(tradesAll_adj[i].drawdown) > stop_loss_thresh)
        {
          tradesAll_adj[i].drawdown = -stop_loss_thresh;
          tradesAll_adj[i].pl = -stop_loss_thresh;
        }
      
        if(!perf_dates.contains(tradesAll_adj[i].date)) //first date entry
        {
         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
         }
         
         pnl = tradesAll_adj[i].pl;
         perf_dates.add(tradesAll_adj[i].date);
         //System.out.println(perf_dates.get(perf_dates.size()-1));
        }
        else //already contains the date, so add on pnl
        {pnl = pnl + tradesAll_adj[i].pl;}
        
      }
      //add the final day 
      perf_returns.add(pnl);
      
      
//       for(i=0;i<perf_dates.size();i++)
//       {     
//         System.out.println(perf_dates.get(perf_dates.size() - 1 - i) + " " + perf_returns.get(perf_returns.size() - 1 - i));
//       }
      
      
      //----- now create a new cumsum plot and trade analysis plot
//       System.out.println(n_trading_days + " " + perf_returns.size() + " " + perf_dates.size() + " " + mdfaAnalysisCanvas.dates.length);
//       if(n_trading_days != perf_returns.size())
//       {System.out.println("There is problem with sizing");}
      n_trading_days = Math.min(mdfaAnalysisCanvas.dates.length, perf_returns.size());
      
      String[] dates = new String[n_trading_days];
      double[] strat_returns = new double[n_trading_days];
      for(i = 0; i < n_trading_days; i++)
      {
        dates[n_trading_days - i - 1]         = perf_dates.get(perf_dates.size() - i - 1);  
        strat_returns[n_trading_days - i - 1] = (double)perf_returns.get(perf_returns.size() - i - 1);     
      }   
   

      double[] mstd = mean_std(strat_returns); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];
      port_mean = mstd[0]; 
      double[] cum_port_returns = cumsum(strat_returns,n_trading_days);    
      max_drawdown = computeDrawdown(cum_port_returns); 
      n_sections = (int)(n_trading_days/length_rank);
      min_rank = 1.0;
      double bRatio = 0;
      
      sharpeText.setText(""+df.format(sharpe_ratio));
      maxDrawText.setText(""+df3.format(max_drawdown));
      minRankText.setText(""+df3.format(min_rank));
      avgRankText.setText(""+df3.format(bRatio));     
     
      if(tradeRuleCombo.getSelectedIndex() == 1) {Integer.toString(stopLossBar.getValue());}
      String time_int = ", "+(String)tradingStartCombo.getSelectedItem() + "-" + (String)tradingEndCombo.getSelectedItem() + ", ";        
      String[] filtab = (filter_file.getName()).split("[.]+");
      
      //mdfaAnalysisCanvas.setDates(dates);
      
      if(mdfaAnalysisCanvas.sim_series.size() == 1)
      {
       mdfaAnalysisCanvas.addSeries(cum_port_returns,new String(""+df3.format(sharpe_ratio)+", " +df3.format(max_drawdown)+", "+df.format(bRatio)+
                              ", " + stop_loss_thresh + ", " + time_int + nobs + ", " + filtab[0]));
      }
      else
      {
        mdfaAnalysisCanvas.setSeries(cum_port_returns,new String(""+df3.format(sharpe_ratio)+", " +df3.format(max_drawdown)+", "+df.format(bRatio)+
                              ", " + stop_loss_thresh + ", " + time_int + nobs + ", " + filtab[0]));
      }
      
      for(i=0;i<mdfaAnalysisCanvas.sim_series.size();i++) {plotHist[i].setEnabled(true); plotHist[i].setSelected(true);}
     
   
      //------ now set new trade analysis plot --
      tradeAnalysis.keepFrameFixed(true);
      tradeAnalysis.setStopLossLine(stop_loss_thresh);
      tradeAnalysis.setTrades(tradesAll_adj);
      tradeAnalysis.go();
      launchTradeAnalysis.setEnabled(true);
      //----------------------------------------------------------------------------

    }
  }
  
    
  public void applyTakeProfit()
  {
    int i;
    int n_trades = trade_perf.size();
    double pnl = 0;
    ArrayList<String> perf_dates = new ArrayList<String>();
    ArrayList<Double> perf_returns = new ArrayList<Double>();
    
    //stop_loss_tresh given in .00000 format 
    //create a tradesAll array
    if(n_trades > 0 && readyTradeAnalysis)
    {
    
      tradesAll_adj = new MDFATrade[n_trades];    
    
      //adjust all trades larger than stop_loss
      for(i = 0; i < n_trades; i++)
      {
        
        tradesAll_adj[i] = new MDFATrade(trade_perf.get(i).date, trade_perf.get(i).startTransTime, trade_perf.get(i).endTransTime, trade_perf.get(i).drawdown,
                                         trade_perf.get(i).drawup, trade_perf.get(i).pl);
        
                                    
        if(tradesAll_adj[i].drawup > take_profit)
        {
          tradesAll_adj[i].drawup = take_profit;
          tradesAll_adj[i].pl = take_profit;
        }
      
        if(!perf_dates.contains(tradesAll_adj[i].date)) //first date entry
        {
         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
         }
         
         pnl = tradesAll_adj[i].pl;
         perf_dates.add(tradesAll_adj[i].date);
         //System.out.println(perf_dates.get(perf_dates.size()-1));
        }
        else //already contains the date, so add on pnl
        {pnl = pnl + tradesAll_adj[i].pl;}
        
      }
      //add the final day 
      perf_returns.add(pnl);
      
      
//       for(i=0;i<perf_dates.size();i++)
//       {     
//         System.out.println(perf_dates.get(perf_dates.size() - 1 - i) + " " + perf_returns.get(perf_returns.size() - 1 - i));
//       }
      
      
      //----- now create a new cumsum plot and trade analysis plot
//       System.out.println(n_trading_days + " " + perf_returns.size() + " " + perf_dates.size() + " " + mdfaAnalysisCanvas.dates.length);
//       if(n_trading_days != perf_returns.size())
//       {System.out.println("There is problem with sizing");}
      n_trading_days = Math.min(mdfaAnalysisCanvas.dates.length, perf_returns.size());
      
      String[] dates = new String[n_trading_days];
      double[] strat_returns = new double[n_trading_days];
      for(i = 0; i < n_trading_days; i++)
      {
        dates[n_trading_days - i - 1]         = perf_dates.get(perf_dates.size() - i - 1);  
        strat_returns[n_trading_days - i - 1] = (double)perf_returns.get(perf_returns.size() - i - 1);     
      }   
   

      double[] mstd = mean_std(strat_returns); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];
      port_mean = mstd[0]; 
      double[] cum_port_returns = cumsum(strat_returns,n_trading_days);    
      max_drawdown = computeDrawdown(cum_port_returns); 
      n_sections = (int)(n_trading_days/length_rank);
      min_rank = 1.0;
      double bRatio = 0;
      
      sharpeText.setText(""+df.format(sharpe_ratio));
      maxDrawText.setText(""+df3.format(max_drawdown));
      minRankText.setText(""+df3.format(min_rank));
      avgRankText.setText(""+df3.format(bRatio));     
     
      if(tradeRuleCombo.getSelectedIndex() == 1) {Integer.toString(stopLossBar.getValue());}
      String time_int = ", "+(String)tradingStartCombo.getSelectedItem() + "-" + (String)tradingEndCombo.getSelectedItem() + ", ";        
      String[] filtab = (filter_file.getName()).split("[.]+");
      
      //mdfaAnalysisCanvas.setDates(dates);
      
      if(mdfaAnalysisCanvas.sim_series.size() == 1)
      {
       mdfaAnalysisCanvas.addSeries(cum_port_returns,new String(""+df3.format(sharpe_ratio)+", " +df3.format(max_drawdown)+", "+df.format(bRatio)+
                              ", " + stop_loss_thresh + ", " + time_int + nobs + ", " + filtab[0]));
      }
      else
      {
        mdfaAnalysisCanvas.setSeries(cum_port_returns,new String(""+df3.format(sharpe_ratio)+", " +df3.format(max_drawdown)+", "+df.format(bRatio)+
                              ", " + stop_loss_thresh + ", " + time_int + nobs + ", " + filtab[0]));
      }
      
      for(i=0;i<mdfaAnalysisCanvas.sim_series.size();i++) {plotHist[i].setEnabled(true); plotHist[i].setSelected(true);}
     
   
      //------ now set new trade analysis plot --
      tradeAnalysis.keepFrameFixed(true);
      tradeAnalysis.setStopLossLine(stop_loss_thresh);
      tradeAnalysis.setTrades(tradesAll_adj);
      tradeAnalysis.go();
      launchTradeAnalysis.setEnabled(true);
      //----------------------------------------------------------------------------

    }
  }
      
    
    
  public void parseAndPlotReturns()
  {
     n_trading_days = total_perf.size();
     String[] dates = new String[n_trading_days];
     double[] strat_returns = new double[n_trading_days];
     double[] stock_rets = new double[n_trading_days];
     
     
     String delims = "[,]+";
     String[] tokens;   
     
     for(int i = 0; i < n_trading_days; i++)
     {
       tokens = (total_perf.get(i)).split(delims); 
       dates[i] = tokens[0];  
       strat_returns[i] = (double)returns_all.get(i);
       stock_rets[i] = (double)stock_perf.get(i);
     }
    
     double maxDD = 0;
     //---- trade analysis plot - first set the trade canvas and then plot---------
     MDFATrade[] tradesAll = new MDFATrade[trade_perf.size()];
     for(int i = 0; i < trade_perf.size(); i++)
     {
      tradesAll[i] = trade_perf.get(i);
      if(tradesAll[i].drawdown < maxDD) {maxDD = tradesAll[i].drawdown;}
     } 
     
     readyTradeAnalysis = false; 
     System.out.println("Setting up tradeAnalysis...");
     adjustSLBar.setMaximum((int)(10000*Math.abs(maxDD)));
     adjustSLBar.setValue((int)(10000*Math.abs(maxDD)) - 2);
      
     
     tradeAnalysis.setTrades(tradesAll);
     tradeAnalysis.go();
     launchTradeAnalysis.setEnabled(true);
     //----------------------------------------------------------------------------
     
     
     
      double[] mstd = mean_std(strat_returns); 
      sharpe_ratio = Math.sqrt(250)*mstd[0]/mstd[1];
      
      port_mean = mstd[0]; 

      double[] cum_port_returns = cumsum(strat_returns,n_trading_days); 
      double[] stock_returns = cumsum(stock_rets,n_trading_days);
      
      max_drawdown = computeDrawdown(cum_port_returns);
      
      n_sections = (int)(n_trading_days/length_rank);
      
      min_rank = 1.0;
      for(int i=0;i<n_sections;i++)
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
      maxDrawText.setText(""+df3.format(max_drawdown));
      minRankText.setText(""+df3.format(min_rank));
      avgRankText.setText(""+df3.format(bRatio));     
     
      String sl = "00";
      if(tradeRuleCombo.getSelectedIndex() == 1) {sl = Integer.toString(stopLossBar.getValue());}
      String time_int = ", "+(String)tradingStartCombo.getSelectedItem() + "-" + (String)tradingEndCombo.getSelectedItem() + ", ";        
      String[] filtab = (filter_file.getName()).split("[.]+");
      
      mdfaAnalysisCanvas.setDates(dates);
      mdfaAnalysisCanvas.setStockReturns(stock_returns);
      mdfaAnalysisCanvas.addSeries(cum_port_returns,new String(""+df3.format(sharpe_ratio)+", " +df3.format(max_drawdown)+", "+df.format(bRatio)+
                              ", " + sl + time_int + nobs + ", " + filtab[0]));
  
      for(int i=0;i<mdfaAnalysisCanvas.sim_series.size();i++) {plotHist[i].setEnabled(true); plotHist[i].setSelected(true);}
      plotStock.setEnabled(true); plotStock.setSelected(false);
  
      readyTradeAnalysis = true;
  }
  
  //compute the rolling indicators (uses the first performance trajectory)
  public void computeRollingIndicators()
  {
   if(mdfaAnalysisCanvas.sim_series.size() > 0)
   {
    int i;
    rolling_ind.clear();
    double[] xt = mdfaAnalysisCanvas.sim_series.get(0);
    double[] rets = new double[xt.length];
    
    //compute returns first
    rets[0] = 0;
    for(i=1;i<xt.length;i++) {rets[i] = xt[i] - xt[i-1];}     
    for(i=0;i<rolling_length;i++) {rolling_ind.add("0 0 0 0");}
  
    for(i = rolling_length; i < xt.length; i++)
    {
       double[] mstd = mean_std(rets,i-rolling_length,i); 
       
       double sh = Math.sqrt(250)*mstd[0]/mstd[1];
       double rc = rankCoefficient(rets,i-rolling_length,i); 
       double ui = ulcerIndex(rets,i-rolling_length,i);     
       double kp = kellyIndex(rets,i-rolling_length,i); 
                    
       rolling_ind.add(sh + " " + rc + " " + ui + " " + kp);     
    }
   
    mdfaAnalysisCanvas.setRollingIndicators(rolling_ind.toArray(new String[0]));
    
    for(i=0;i<4;i++) {rollingCheck[i].setEnabled(true);}
    
   }
  }  
    
  public void clearReturnSeries()
  {
    for(int i = 0; i < plotHist.length; i++) {plotHist[i].setSelected(false); plotHist[i].setEnabled(false);}      
    for(int i = 0; i < 4; i++) {rollingCheck[i].setSelected(false); rollingCheck[i].setEnabled(false); mdfaAnalysisCanvas.plot_rolling[i] = false;}
    mdfaAnalysisCanvas.clearSeries();
  }

  public void clearVarSeries()
  {
    for(int i = 0; i < plotVar.length; i++) {plotVar[i].setSelected(false); plotVar[i].setEnabled(false);}  
    mdfaAnalysisCanvas.clearVarSeries();
  }    
    
    
    
    
  public void initTradeComponents() 
  {

        tradeAnalysisPanel = new JPanel();
        adjustSLBar = new javax.swing.JScrollBar();
        adjustSLLabel = new javax.swing.JLabel();
        adjustSLText = new javax.swing.JTextField();
       
        takeProfitLabel = new javax.swing.JLabel();
        takeProfitBar = new javax.swing.JScrollBar();
        takeProfitText = new javax.swing.JTextField();

        adjustSLBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        adjustSLBar.setMinimum(5); 
        adjustSLBar.setMaximum(300); 
        adjustSLBar.setValue(250); 
        adjustSLBar.setUnitIncrement(1);
         
        
        adjustSLLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        adjustSLLabel.setText("Adjust Stop-Loss");

        adjustSLText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        adjustSLText.setText(".0200");

        tradeAnalysis.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));

        javax.swing.GroupLayout tradeAnalysisLayout = new javax.swing.GroupLayout(tradeAnalysis);
        tradeAnalysis.setLayout(tradeAnalysisLayout);
        tradeAnalysisLayout.setHorizontalGroup(
            tradeAnalysisLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        tradeAnalysisLayout.setVerticalGroup(
            tradeAnalysisLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 574, Short.MAX_VALUE)
        );

        takeProfitLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        takeProfitLabel.setText("Adjust Take Profit");

        takeProfitBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        takeProfitBar.setMinimum(0); 
        takeProfitBar.setMaximum(600); 
        takeProfitBar.setValue(1); 
        takeProfitBar.setUnitIncrement(1);
        
        
        
        takeProfitText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        takeProfitText.setText(".0000");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(tradeAnalysisPanel);
        tradeAnalysisPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(tradeAnalysis, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 418, Short.MAX_VALUE)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(adjustSLLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 106, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(adjustSLBar, javax.swing.GroupLayout.PREFERRED_SIZE, 353, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(adjustSLText, javax.swing.GroupLayout.PREFERRED_SIZE, 61, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addComponent(takeProfitLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 117, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(takeProfitBar, javax.swing.GroupLayout.PREFERRED_SIZE, 353, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(takeProfitText, javax.swing.GroupLayout.PREFERRED_SIZE, 61, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tradeAnalysis, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(adjustSLText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(adjustSLBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(adjustSLLabel, javax.swing.GroupLayout.Alignment.TRAILING)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(takeProfitText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(takeProfitBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(takeProfitLabel, javax.swing.GroupLayout.Alignment.TRAILING)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        
        
        
        
        AdjustmentListener tradeListener = new AdjustmentListener()  {
         public void adjustmentValueChanged(AdjustmentEvent e) {        
         
             //double take_profit = 0;
             if(e.getAdjustable() == adjustSLBar)
             {
               stop_loss_thresh = .0001*adjustSLBar.getValue();
               adjustSLText.setText(""+stop_loss_thresh);
               applyStopLoss();
             }
             else if(e.getAdjustable() == takeProfitBar)
             {
               take_profit = .0001*takeProfitBar.getValue();
               takeProfitText.setText(""+take_profit);
               tradeAnalysis.setTakeProfit(take_profit);    
               applyTakeProfit();
             }
         }
        };
        
        adjustSLBar.addAdjustmentListener(tradeListener);
        takeProfitBar.addAdjustmentListener(tradeListener);
         
    }                       


          
 
    
    
    
    
    
    
    
    
    @SuppressWarnings({ "rawtypes", "unchecked" })
	public void initStrategyPanel()
    {

        mdfaAnalysisCanvas = new mdfaStrategyCanvas(200, 400, 300, 1);
        mdfaSettingsPanel = new JPanel();
        mdfaStratLabel = new JLabel();
        insampStartCombo = new JComboBox();
        tradingInSampLabel = new JLabel();
        tradingStartLabel = new JLabel();
        tradingStartCombo = new JComboBox();
        tradingEndCombo = new JComboBox();
        tradingEndLabel = new JLabel();
        nobsLabel = new JLabel();
        nobsBar = new JScrollBar();
        nobsText = new JTextField();
        tradeRuleLabel = new JLabel();
        tradeRuleCombo = new JComboBox();
        paramFileLabel = new JLabel();
        paramFileButton = new JButton();
        paramFileText = new JTextField();
        paramFileLabel1 = new JLabel();
        paramFileButton1 = new JButton();
        paramFileText1 = new JTextField();        
        historicalFileLabel = new JLabel();
        historicalFileButton = new JButton();
        historicalFileText = new JTextField();
        ncoresText = new JTextField();
        ncoresBar = new JScrollBar();
        ncoresLabel = new JLabel();
        stratProgressBar = new JProgressBar();
        computeButton = new JButton();
        varRatioPanel = new JPanel();
        kdiffLabel = new JLabel();
        kdiffText = new JTextField();
        kdiffBar = new JScrollBar();
        varRatioPanelLabel = new JLabel();
        varRatioEndCombo = new JComboBox();
        varRatioEndLabel = new JLabel();
        varRatioStartCombo = new JComboBox();
        varRatioStartLabel = new JLabel();
        varRatioProgressBar = new JProgressBar();
        varRatioButton = new JButton();
        plotPanel = new JPanel();
        deleteHistButton = new JButton();
        deleteVarRatioButton = new JButton();
        avgRankLabel = new JLabel(); maxDrawLabel = new JLabel();
        avgRankText = new JTextField(); maxDrawText = new JTextField();
        minRankLabel = new JLabel(); minRankText = new JTextField();
        sharpeText = new JTextField(); sharpeLabel = new JLabel();
        windowSizeLabel = new JLabel(); samplesText = new JTextField();
        samplesBar = new JScrollBar(); stopLossBar = new JScrollBar();
        stopLossLabel = new JLabel(); stopLossText = new JTextField();
        shortBox = new JCheckBox("IBData"); shortBox.setSelected(true);
        firm_closeBox = new JCheckBox("Firm Close"); firm_closeBox.setSelected(false);
        paramFileText.setText("./stuff_time.txt");
        historicalFileText.setText("./EURUSD.FXCM.dat");
        
        rollingCheck = new JCheckBox[4];
        rollingCheck[0] = new JCheckBox("Sharpe");
        rollingCheck[1] = new JCheckBox("Rank");
        rollingCheck[2] = new JCheckBox("Ulcer");
        rollingCheck[3] = new JCheckBox("Kelly");
        for(int i = 0; i < 4; i++) {rollingCheck[i].setSelected(false); rollingCheck[i].setEnabled(false);}
        
        max_series = 20; 
        plotStock = new JCheckBox("Underlying");
        plotHist = new JCheckBox[max_series]; plotVar = new JCheckBox[max_series];
        for(int i = 0; i < max_series; i++) 
        {
          plotHist[i] = new JCheckBox(""+i,false); plotHist[i].setSelected(false); plotHist[i].setEnabled(false);
          plotVar[i] = new JCheckBox(""+i,false);  plotVar[i].setSelected(false); plotVar[i].setEnabled(false);
        }
        plotStock.setSelected(false); plotStock.setEnabled(false); 
        mdfaAnalysisCanvas.setBackground(new Color(1, 1, 1));
        mdfaAnalysisCanvas.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        mdfaAnalysisCanvas.setForeground(new Color(1, 1, 1));

        JPanel rollingPanel = new JPanel();
        GroupLayout rollingIndLayout = new GroupLayout(rollingPanel);
        rollingPanel.setLayout(rollingIndLayout);
        rollingIndLayout.setHorizontalGroup(
            rollingIndLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(rollingIndLayout.createSequentialGroup()
                .addContainerGap()        
                 .addComponent(rollingCheck[0]).addComponent(rollingCheck[1]).addComponent(rollingCheck[2]).addComponent(rollingCheck[3])));
        rollingIndLayout.setVerticalGroup(
            rollingIndLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(rollingIndLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                .addComponent(rollingCheck[0]).addComponent(rollingCheck[1]).addComponent(rollingCheck[2]).addComponent(rollingCheck[3])));
        rollingPanel.setLayout(rollingIndLayout);       
        
        
        
        
//         GroupLayout mdfaAnalysisCanvasLayout = new GroupLayout(mdfaAnalysisCanvas);
//         mdfaAnalysisCanvas.setLayout(mdfaAnalysisCanvasLayout);
//         mdfaAnalysisCanvasLayout.setHorizontalGroup(
//             mdfaAnalysisCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
//             .addGap(0, 0, Short.MAX_VALUE)
//         );
//         mdfaAnalysisCanvasLayout.setVerticalGroup(
//             mdfaAnalysisCanvasLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
//             .addGap(0, 435, Short.MAX_VALUE)
//         );

        mdfaSettingsPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        //mdfaStratLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        mdfaStratLabel.setText("MDFA Strategy Settings");

        //insampStartCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        insampStartCombo.setModel(new DefaultComboBoxModel(new String[] { "08:30", "08:45", "09:00", "09:15", "09:30" }));

        //tradingInSampLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradingInSampLabel.setText("In-sample Start Time");

        //tradingStartLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradingStartLabel.setText("Trading Start Time");

        //tradingStartCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradingStartCombo.setModel(new DefaultComboBoxModel(new String[] { "09:30", "09:45", "10:00", "10:15", "10:30", "10:45", "11:00", "11:15", "11:30", "11:45", "12:00" }));

        //tradingEndCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradingEndCombo.setModel(new DefaultComboBoxModel(new String[] {"13:00", "13:15", "13:30", "13:45", "14:00", "14:15", "14:30", "14:45","15:00", "15:15", "15:30", "15:45", "16:00" }));
        tradingEndCombo.setSelectedIndex(11);
        //tradingEndLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradingEndLabel.setText("Trading End Time");

        //nobsLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        nobsLabel.setText("TradingCost");

        //nobsBar.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        nobsBar.setMaximum(100);
        nobsBar.setMinimum(0);
        nobsBar.setOrientation(JScrollBar.HORIZONTAL);
        nobsBar.setValue(0);

        //nobsText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        nobsText.setText("0");

        //tradeRuleLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradeRuleLabel.setText("Trading Rule");

        //tradeRuleCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        tradeRuleCombo.setModel(new DefaultComboBoxModel(new String[] { "Binary", "Stop-Loss Binary", "Signal Strength" }));

        //paramFileLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        paramFileLabel.setText("Parameter FIle");

        //paramFileButton.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        paramFileButton.setText("Choose File");

        //paramFileLabel1.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        paramFileLabel1.setText("MetaParam File");

        //paramFileButton1.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        paramFileButton1.setText("Choose File");        
        
        
        //historicalFileLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        historicalFileLabel.setText("Historical Data");

        //historicalFileButton.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        historicalFileButton.setText("Choose File");

        //ncoresText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        ncoresText.setText("10");

        //ncoresBar.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        ncoresBar.setMaximum(30);
        ncoresBar.setMinimum(1);
        ncoresBar.setOrientation(JScrollBar.HORIZONTAL);
        ncoresBar.setValue(10);

        //ncoresLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        ncoresLabel.setText("Number of Threads");

        computeButton.setText("Compute Strategy");

        //sharpeLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sharpeLabel.setText("Sharpe");

        sharpeText.setColumns(4);
        //sharpeText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sharpeText.setForeground(new java.awt.Color(0, 255, 12));
        sharpeText.setText("0");

        //avgRankLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        avgRankLabel.setText("RankRatio");

        avgRankText.setColumns(4);
        avgRankText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        avgRankText.setForeground(new java.awt.Color(0, 255, 47));
        avgRankText.setText("0");

        maxDrawText.setColumns(4);
        //maxDrawText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        maxDrawText.setForeground(new java.awt.Color(0, 255, 13));
        maxDrawText.setText("0");

        //maxDrawLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        maxDrawLabel.setText("MaxDrawDown");

        //minRankLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        minRankLabel.setText("Min Ratio");

        minRankText.setColumns(4);
        //minRankText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        minRankText.setForeground(new java.awt.Color(15, 255, 0));
        minRankText.setText("0");

        multivarBox = new JCheckBox("Multivariate", true);
        autoVar = new JCheckBox("Auto",false);
        
  
  
  
  
  
  
  
  
  
  
  

        //stopLossBar.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        stopLossBar.setMaximum(211);
        stopLossBar.setMinimum(0);
        stopLossBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        stopLossBar.setValue(0);

        stopLossText.setColumns(7);
        //stopLossText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        stopLossText.setText("10");

        //stopLossLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        stopLossLabel.setText("StpLs");

        recompDayLabel = new JLabel("recomp Days");
        recompDayBar = new JScrollBar(JScrollBar.HORIZONTAL,0,1,0,20);
        recompDayText = new JTextField("0");
        recompDayText.setColumns(7);
        launchTradeAnalysis = new JButton("Trades");
        launchTradeAnalysis.setEnabled(false);
        
        JPanel recompCon;
        
        new JPanel(); recompCon = new JPanel();
        
//         GroupLayout paramLayout = new GroupLayout(stopCon); paramLayout.setAutoCreateGaps(true);
//         paramLayout.setAutoCreateContainerGaps(true);        
//         paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
//          .addComponent(stopLossLabel).addComponent(stopLossBar).addComponent(stopLossText));
//         paramLayout.setVerticalGroup(
//          paramLayout.createSequentialGroup()
//          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
//             .addComponent(stopLossLabel).addComponent(stopLossBar).addComponent(stopLossText)));
//         stopCon.setLayout(paramLayout);
        
        GroupLayout paramLayout = new GroupLayout(recompCon); paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
         .addComponent(recompDayLabel).addComponent(recompDayBar).addComponent(recompDayText));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
            .addComponent(recompDayLabel).addComponent(recompDayBar).addComponent(recompDayText)));
        recompCon.setLayout(paramLayout);
        
        
        
        
        
        JPanel checks = new JPanel();
        checks.add(multivarBox); checks.add(shortBox); checks.add(firm_closeBox); checks.setLayout(new FlowLayout());
        
        javax.swing.GroupLayout mdfaSettingsPanelLayout = new javax.swing.GroupLayout(mdfaSettingsPanel);
        mdfaSettingsPanel.setLayout(mdfaSettingsPanelLayout);
        mdfaSettingsPanelLayout.setHorizontalGroup(
            mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(mdfaStratLabel)
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(nobsLabel)
                                .addGap(6, 6, 6)
                                .addComponent(nobsBar, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(nobsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addComponent(tradeRuleLabel)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(tradeRuleCombo, 0, 1, Short.MAX_VALUE))
                                .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(tradingInSampLabel)
                                        .addComponent(tradingStartLabel)
                                        .addComponent(tradingEndLabel))
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(tradingStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(insampStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(tradingEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(checks)
                                .addGap(20, 20, 20))
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(stopLossLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(stopLossBar, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(stopLossText, javax.swing.GroupLayout.PREFERRED_SIZE, 24, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(paramFileLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramFileButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramFileText, javax.swing.GroupLayout.PREFERRED_SIZE, 142, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(108, 108, 108)
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(historicalFileLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(historicalFileButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(historicalFileText))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(sharpeLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(maxDrawText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(ncoresLabel)
                                .addGap(6, 6, 6)
                                .addComponent(ncoresBar, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(minRankLabel)
                                .addComponent(maxDrawLabel)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addComponent(ncoresText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, mdfaSettingsPanelLayout.createSequentialGroup()
                                .addGap(0, 0, Short.MAX_VALUE)
                                .addComponent(minRankText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(sharpeText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                    .addComponent(avgRankLabel)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(avgRankText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addComponent(recompCon)
                            .addComponent(computeButton))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(paramFileLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramFileButton1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramFileText1, javax.swing.GroupLayout.PREFERRED_SIZE, 142, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        mdfaSettingsPanelLayout.setVerticalGroup(
            mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addComponent(mdfaStratLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(nobsLabel)
                                            .addComponent(nobsBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                        .addGap(10, 10, 10))
                                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, mdfaSettingsPanelLayout.createSequentialGroup()
                                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                            .addComponent(nobsText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(checks))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)))
                                .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(tradeRuleLabel)
                                    .addComponent(tradeRuleCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(stopLossLabel)
                                    .addComponent(stopLossBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
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
                                    .addComponent(tradingEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                                .addGap(24, 24, 24)
                                .addComponent(stopLossText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGap(1, 1, 1)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(paramFileLabel)
                            .addComponent(paramFileButton)
                            .addComponent(paramFileText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(mdfaSettingsPanelLayout.createSequentialGroup()
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(historicalFileLabel)
                            .addComponent(historicalFileButton)
                            .addComponent(historicalFileText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(ncoresText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(ncoresLabel)
                                .addComponent(ncoresBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(computeButton)
                        .addComponent(recompCon)
                        
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(sharpeLabel)
                            .addComponent(sharpeText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(maxDrawLabel)
                            .addComponent(maxDrawText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(avgRankLabel)
                            .addComponent(avgRankText, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(minRankLabel)
                                .addComponent(minRankText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mdfaSettingsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(paramFileLabel1)
                            .addComponent(paramFileButton1)
                            .addComponent(paramFileText1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addContainerGap())
        );        
    
        

        
        varRatioPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        //kdiffLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        kdiffLabel.setText("K");

        //kdiffText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        kdiffText.setText("10");

        //kdiffBar.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        kdiffBar.setMaximum(30);
        kdiffBar.setMinimum(1);
        kdiffBar.setOrientation(JScrollBar.HORIZONTAL);
        kdiffBar.setValue(10);

        //varRatioPanelLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        varRatioPanelLabel.setText("Variance Ratio Settings");

        //varRatioEndCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        varRatioEndCombo.setModel(new DefaultComboBoxModel(new String[] { "10:00", "10:15", "10:30", "10:45",  "11:00", "11:15", "11:30", "11:45", "12:00", "12:15", "12:30", "12:45", "13:00", "13:15", "13:30", "13:45", "14:00", "14:15", "14:30", "14:45", "15:00", "15:15", "15:30", "15:45", "16:00" }));

        //varRatioEndLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        varRatioEndLabel.setText("VarRatio End Time");

        //varRatioStartCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        varRatioStartCombo.setModel(new DefaultComboBoxModel(new String[] { "03:00", "03:15", "03:30", "03:45", "04:00", "04:15", "04:30", "04:45", "05:00", "05:15", "05:30", "05:45", "06:00", "06:15", "06:30", "06:45", "07:00", "07:15", "07:30", "07:45", "08:00", "08:15", "08:30", "08:45", "09:00", "09:15", "09:30", "09:45", "10:00", "10:15", "10:30", "10:45", "11:00", "11:15", "11:30", "11:45", "12:00", "12:15", "12:30" }));

        //varRatioStartLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        varRatioStartLabel.setText("VarRatio Start Time");

        varRatioButton.setText("Compute VarRatio");
        
        JPanel varRatioButtonPanel = new JPanel(); 
        varRatioButtonPanel.add(varRatioButton); varRatioButtonPanel.add(autoVar);
        varRatioButtonPanel.setLayout(new FlowLayout());
        
        //windowSizeLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        windowSizeLabel.setText("M-Samples");
        
        //samplesText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        samplesText.setText("10");

        //samplesBar.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        samplesBar.setMaximum(80);
        samplesBar.setMinimum(10);
        samplesBar.setOrientation(JScrollBar.HORIZONTAL);
        samplesBar.setValue(20);        
        
        
        javax.swing.GroupLayout varRatioPanelLayout = new javax.swing.GroupLayout(varRatioPanel);
        varRatioPanel.setLayout(varRatioPanelLayout);
        varRatioPanelLayout.setHorizontalGroup(
            varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(varRatioPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(varRatioPanelLayout.createSequentialGroup()
                        .addComponent(kdiffLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(kdiffBar, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(kdiffText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(varRatioPanelLabel)
                    .addGroup(varRatioPanelLayout.createSequentialGroup()
                        .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(varRatioStartLabel)
                            .addComponent(varRatioEndLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(varRatioStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(varRatioEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addComponent(varRatioButtonPanel)
                    .addGroup(varRatioPanelLayout.createSequentialGroup()
                      .addComponent(rollingPanel))
                    .addGroup(varRatioPanelLayout.createSequentialGroup()
                        .addComponent(windowSizeLabel)
                        .addGap(24, 24, 24)
                        .addComponent(samplesBar, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(samplesText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(119, Short.MAX_VALUE))
        );
        varRatioPanelLayout.setVerticalGroup(
            varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(varRatioPanelLayout.createSequentialGroup()
                .addGap(8, 8, 8)
                .addComponent(varRatioPanelLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(kdiffText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(kdiffLabel)
                        .addComponent(kdiffBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(samplesText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(varRatioPanelLayout.createSequentialGroup()
                        .addComponent(windowSizeLabel)
                        .addGap(2, 2, 2))
                    .addComponent(samplesBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(varRatioStartLabel)
                    .addComponent(varRatioStartCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(varRatioPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(varRatioEndLabel)
                    .addComponent(varRatioEndCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(varRatioButtonPanel)
                .addComponent(rollingPanel)
                .addGap(95, 95, 95))
        );        
        



        plotPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        //deleteHistButton.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        deleteHistButton.setText("Delete Historicals");

        //deleteVarRatioButton.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        deleteVarRatioButton.setText("Delete VarRatios");


        GroupLayout plotPanelLayout = new GroupLayout(plotPanel);
        plotPanel.setLayout(plotPanelLayout);
        plotPanelLayout.setHorizontalGroup(
            plotPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(deleteHistButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[0])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[1])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[2])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[3])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[4])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[5])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[6])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[7])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[8])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[9])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[10])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[11])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[12])     
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[13])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[14])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotHist[15])
                .addComponent(launchTradeAnalysis)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, 139, Short.MAX_VALUE)
                .addComponent(deleteVarRatioButton)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[0])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[1])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[2])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[3])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[4])
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotVar[5])
                .addComponent(plotStock)
                .addGap(104, 104, 104))
        );
        plotPanelLayout.setVerticalGroup(
            plotPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(plotPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                .addComponent(deleteHistButton)
                .addComponent(plotHist[0])
                .addComponent(plotHist[1])
                .addGroup(plotPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(plotHist[2])
                    .addComponent(plotHist[3])
                    .addGroup(plotPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(plotHist[4])
                        .addComponent(plotHist[5])
                        .addComponent(plotHist[6])
                        .addComponent(plotHist[7])
                        .addComponent(plotHist[8])
                        .addComponent(plotHist[9])
                        .addComponent(plotHist[10])
                        .addComponent(plotHist[11])
                        .addComponent(plotHist[12])
                        .addComponent(plotHist[13])
                        .addComponent(plotHist[14])
                        .addComponent(plotHist[15])
                        .addComponent(launchTradeAnalysis)))
                .addGroup(plotPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(plotVar[0])
                    .addComponent(plotVar[1])
                    .addComponent(plotVar[2])
                    .addComponent(plotVar[3])
                    .addComponent(plotVar[4])
                    .addComponent(plotVar[5])
                    .addComponent(plotStock)
                    .addComponent(deleteVarRatioButton)))
        );

        GroupLayout layout = new GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(plotPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(mdfaAnalysisCanvas, GroupLayout.Alignment.TRAILING, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(mdfaSettingsPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(varRatioPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(mdfaAnalysisCanvas, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addComponent(varRatioPanel, GroupLayout.PREFERRED_SIZE, 212, GroupLayout.PREFERRED_SIZE)
                    .addComponent(mdfaSettingsPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        
        
       //Listeners ------
       nobsBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setTradingCost(((JScrollBar)e.getSource()).getValue());
               nobsText.setText(""+((JScrollBar)e.getSource()).getValue());
            }
         });
         
       stopLossBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setStopLoss((((JScrollBar)e.getSource()).getValue())*.0001);
               stopLossText.setText(""+((JScrollBar)e.getSource()).getValue());
            }
         });  
         
       recompDayBar.addAdjustmentListener(new AdjustmentListener()  {
            public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setRecompDay((((JScrollBar)e.getSource()).getValue()));
               recompDayText.setText(""+((JScrollBar)e.getSource()).getValue());
            }
         }); 
         
         
       nobsText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_nobs(nobsText.getText());}} );  
       
       
       paramFileText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_parameterFile(paramFileText.getText());}} );
       
       historicalFileText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_historicalFile(historicalFileText.getText());}});       
       
       ncoresBar.addAdjustmentListener(new AdjustmentListener() {
           public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               numberOfThreads(((JScrollBar)e.getSource()).getValue());
               ncoresText.setText(""+n_threads);
              
            }
       });
       
       kdiffBar.addAdjustmentListener(new AdjustmentListener() {
           public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setKDifferences(((JScrollBar)e.getSource()).getValue());
               kdiffText.setText(""+Kdiff);
               if(auto) {computeVarRatio();}
            }
       });
       
       samplesBar.addAdjustmentListener(new AdjustmentListener() {
           public void adjustmentValueChanged(AdjustmentEvent e) 
            {
               setNSamples(((JScrollBar)e.getSource()).getValue());
               samplesText.setText(""+n_samples);
               if(auto) {computeVarRatio();}
               computeRollingIndicators();
            }
       });
       
       
       
     //computeButton,deleteHistButton,deleteVarRatioButton
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
          if(event.getSource() == paramFileButton1)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               setMetaFilterFile(file);
             }
             else {System.out.println("Open command cancelled by user.");}                    
          }          
          else if(event.getSource() == historicalFileButton)
          {
             val = fc.showOpenDialog(parent);
             if(val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               setHistoricalDataFile(file);
             }
             else {System.out.println("Open command cancelled by user.");}                
          }
          else if(event.getSource() == launchTradeAnalysis)
          {
            tradeAnalysisDialog.setModal(false);
            tradeAnalysisDialog.setVisible(true);
          }
          else if(event.getSource() == computeButton) 
          {
            tradeAnalysis.keepFrameFixed(false);
            
            dailyStrat = false;
            if(!meta_mode) 
            {
             try{computeMDFA_EvolutionStrategy();} catch(InterruptedException ie){} catch(ExecutionException ee){}
            }
            else
            {
              for(int i = 0; i < filter_files.length; i++)
              {
                System.out.println("Processing performance for filter " + filter_files[i]);
                filter_file = new File(filter_files[i]);
                try{computeMDFA_EvolutionStrategy();} catch(InterruptedException ie){} catch(ExecutionException ee){}
              }               
            }
          }  
          else if(event.getSource() == deleteHistButton) {clearReturnSeries();}
          else if(event.getSource() == deleteVarRatioButton) {clearVarSeries();}
          else if(event.getSource() == varRatioButton) {computeVarRatio();}       
       }
      };
      
      launchTradeAnalysis.addActionListener(buttonActionListener);
      computeButton.addActionListener(buttonActionListener);
      paramFileButton.addActionListener(buttonActionListener);
      paramFileButton1.addActionListener(buttonActionListener);
      historicalFileButton.addActionListener(buttonActionListener);
      deleteHistButton.addActionListener(buttonActionListener);
      deleteVarRatioButton.addActionListener(buttonActionListener);
      varRatioButton.addActionListener(buttonActionListener);
      
      ItemListener MyItemListener = new ItemListener() 
      {
        public void itemStateChanged(ItemEvent e)
        {         
         int i; boolean sel; Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}
 
         for(i=0;i<max_series;i++)
         {     
           if(source == plotHist[i]) 
           {mdfaAnalysisCanvas.setplot(i,sel);}
           else if(source == plotVar[i])
           {mdfaAnalysisCanvas.setvarplot(i,sel);}
         }
         if(source == plotStock) {mdfaAnalysisCanvas.plotStock(sel);}
         if(source == autoVar) {auto = sel;}
         if(source == multivarBox) {multivarMDFA(sel);}
         if(source == shortBox) {ib_data = sel; System.out.println("Using IB data"); /*shortSell(sel);*/}
         
         for(i=0;i<4;i++) 
         {
           if(source == rollingCheck[i])
           {mdfaAnalysisCanvas.plotRollIndic(i,sel);}
         }
       }  
      };      

      plotStock.addItemListener(MyItemListener);
      shortBox.addItemListener(MyItemListener);
      autoVar.addItemListener(MyItemListener);
      multivarBox.addItemListener(MyItemListener);
      for(int i = 0; i < 6; i++)
      {plotHist[i].addItemListener(MyItemListener); plotVar[i].addItemListener(MyItemListener);}
      for(int i = 0; i < 4; i++)
      {rollingCheck[i].addItemListener(MyItemListener);}
      
      
  }
  
  public void test_parameterFile(String s)
  {
     File file = new File(s);
     setFilterFile(file);  
  }

  public void test_historicalFile(String s)
  {
     File file = new File(s);
     setHistoricalDataFile(file); 
  }  
  
  
  public void test_nobs(String s)
  {   
       int i;
       try{i = Integer.parseInt(s.trim()); if(i >=150 && i < 500) {nobsBar.setValue(i); nobs = i; nobsText.setText(""+nobs);}
                                          else nobsText.setText(""+nobs);}
       catch (NumberFormatException nfe)
       {System.out.println("NumberFormatException: " + nfe.getMessage()); nobsText.setText(""+nobs);}
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
  
  public double rankCoefficient(double[] account, int from, int to)
  {

   int i; int n = to - from;
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double sum = 0.0;
   double spear;
   double[] rets = new double[n];
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i; rets[i] = account[from + i];} 
      
   double[] _acc = cumsum(rets, n);
   
   sort_sims(n, _acc, index);
   
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
     double[] cmaxx = new double[ret.length];
     
     for(i=0;i<ret.length;i++)
     {
      cmaxx[i] = cmax[i] - ret[i]; 
      if(cmaxx[i] > max) {max = cmaxx[i];}
     }
     
     return max; 
  }
  
  public double ulcerIndex(double[] rets)
  {
  
    double SumSq = 0; 
    double MaxValue = -100;
    double mean = 0;
// for T = 1 to NumOfPeriods do
// if Value[T] > MaxValue then MaxValue = Value[T] 
// else SumSq = SumSq + sqr(100 * ((Value[T] / MaxValue) - 1))
// UI = sqrt(SumSq / NumOfPeriods)
//   
    double[] tempval = cumsum(rets,rets.length);
    
    MaxValue = tempval[9];
    for(int i = 0; i < rets.length; i++)
    {
      mean = mean + rets[i]; 
      if(tempval[i] > MaxValue) {MaxValue = tempval[i];}
      else {SumSq = SumSq + (((tempval[i]/MaxValue) - 1.0))*(((tempval[i]/MaxValue) - 1.0));}
      //System.out.println(MaxValue + " " + mean);
    }
    //System.out.println((mean/rets.length) + " " + Math.sqrt(SumSq/rets.length)
    return (mean/rets.length)/Math.sqrt(SumSq/rets.length);
    
  }  
  
  public double ulcerIndex(double[] re, int from, int to)
  {
  
    double SumSq = 0; 
    double MaxValue = -100;
    double mean = 0;
    int n = to - from;

    double[] rets = new double[n]; 
    System.arraycopy(re, from, rets, 0, n);
    
    double[] tempval = cumsum(rets,rets.length);
    
    MaxValue = tempval[0];
    for(int i = 0; i < rets.length; i++)
    {
      mean = mean + rets[i]; 
      if(tempval[i] > MaxValue) {MaxValue = tempval[i];}
      else {SumSq = SumSq + (((tempval[i]/MaxValue) - 1.0))*(((tempval[i]/MaxValue) - 1.0));}
      //System.out.println(MaxValue + " " + mean);
    }
    //System.out.println((mean/rets.length) + " " + Math.sqrt(SumSq/rets.length)
    return (mean/rets.length)/Math.sqrt(SumSq/rets.length);
    
  }    
  
  public double kellyIndex(double[] rets, int from, int to)
  {
    int n = to - from;
    int wins = 0; 
    double winAvg,lossAvg;
    
    winAvg = 0; lossAvg = 0;
    for(int i = from; i < to; i++)
    {
      if(rets[i] > 0) {wins++; winAvg = winAvg + rets[i];}
      else
      {lossAvg = lossAvg - rets[i];}
    }
    
    if(wins>0) {winAvg = winAvg/wins;} 
    if((n - wins) > 0) {lossAvg = lossAvg/(n-wins);}
    
    double winRatio = (double)wins/n;
    
    return (winRatio - (1.0 - winRatio)*(lossAvg/winAvg));
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
  
  public static double[] mean_std( double[] data, int from, int to ) 
  { 

       double mean = 0; 
       final int n = to - from; 
       
       for ( int i=from; i<to; i++ )  {  mean += data[i]; } 
       mean /= n; 

       double sum = 0; 
       for ( int i=from; i<to; i++ ) { final double v = data[i] - mean;  sum += v * v; 
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
    
    
}    












class mdfaStrategyCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int n_obs; //---total number of observations
    double[] tseries;  //raw time series n_rep x n_obs

    int min_obs;
    double dataMax, dataMin, dataNorm, dataMaxVar, dataMinVar, dataNormVar, sdataMax, sdataMin,sdataNorm;
    boolean plot_me; 
    int n_rep; 
    int n_sym = 0;
    int width; int height;
    Graphics2D g2d;
    String position = ""; 
    String[] asset_names;
    boolean names_set = false;
    double[] target_series;
    ArrayList<double[]> sim_series;  //----- plots historical saved plots
    ArrayList<Double[]> var_series;
    int n_hist; 
    boolean[] plots,var_plots;
    boolean plot_tar;
    Color[] colorHist;
    int track_pos_t,track_pos_xtarg;
    int[] track_pos_x, track_pos_xvar,track_pos_indic;  
    double[] acf_series_1;
    double[] acf_series_2;
    double[] cross_12; 
    double[] cross_21;  
    boolean acf1,acf2,cross1,cross2;
    boolean plot_acf;
    int pressed_t;
    String valuetarg;
    boolean plot_tracker = true;
    //------- canvas stuff ----------------------------
    String date; 
    String[] value; 
    DecimalFormat df,df2;
    BasicStroke dashed;
    BasicStroke orig;
    float[] dash1;
    Color myGray;
    Color myGray2;
    Color plotGray;
    Color foreRed;
    Color[] indicColor, indicColorlight;
    int max_series;
    ArrayList<String> varPlotNames, seriesPlotNames;
    Color highlight; 
    int hlight_indicator; 
    String[] info; 
    boolean display;
    Font mono;
    int canvas_shift = 0;
    Ellipse2D ellipse;
    Rectangle2D rectangle;
    DecimalFormat df3;
    boolean dates_set = false; 
    String[] dates; 
    double[] stock_perf; 
    boolean stock_set = false;
    boolean rolling_ind_set = false;
    
    String[] rollingPerf;
    String[] indicValue;
    
    boolean plot_stock = false;
    boolean[] plot_rolling;
    double[][] rollingIndics;
    double[] indicMax,indicMin,indicNorm;
    
    double rankNorm = 0; double sharpeNorm = 0; double ulcerNorm = 0; double kellyNorm = 0;
    double sharpeMax = 10000; double sharpeMin = -19999; double rankMax = 10000; double rankMin = -19999;
    double ulcerMax = 10000; double ulcerMin = -19999; double kellyMax = 10000; double kellyMin = -19999;
    
    public mdfaStrategyCanvas(int w, int h, int nobs, int nrep)
    {
      // Initilize everything-------------------
      int i; 
      n_obs = nobs; n_rep = nrep; 
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;
      max_series = 20;
      n_hist = 0; 
      sim_series = new ArrayList<double[]>(); 
      var_series = new ArrayList<Double[]>();
      
      target_series = new double[n_obs]; plot_tar = false;
  
      colorHist = new Color[21];
      //mono = new Font(parent.getDisplay(), "Monospaced", 10, SWT.NONE);
      mono  = new Font("Monospaced", Font.PLAIN, 12);
      setBackground(Color.BLACK);
      //setBackground(new Color(255,250,250));
      setPreferredSize(new Dimension(w, h));    
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
      myGray = new Color(67,71,73);
      myGray2 = new Color(20,21,19);
      plotGray = new Color(145,63,63);
      df = new DecimalFormat("##.##"); 
      df2 = new DecimalFormat("##.###");
      df3 = new DecimalFormat("##.####");
      foreRed = new Color(232, 138 , 187);
      highlight = new Color(255,255,255);
      plot_acf = false;
      acf1 = false; acf2 = false; cross1 = false; cross2 = false;

      track_pos_indic = new int[4];
      track_pos_x = new int[max_series]; track_pos_xvar = new int[max_series];
      valuetarg = "";
      value = new String[max_series];
      for(i=0;i<max_series;i++) {value[i] = new String(""); track_pos_x[i] = 0; track_pos_xvar[i] = 0;} 
      
      varPlotNames = new ArrayList<String>();
      seriesPlotNames = new ArrayList<String>();
      
      indicValue = new String[4];
      hlight_indicator = -1; display = false;
      colorHist[0] = new Color(103,177,246); colorHist[1] = new Color(255,177,246); 
      colorHist[2] = new Color(103,100,246); colorHist[3] = new Color(103,255,246);  
      colorHist[4] = new Color(103,177,0);   colorHist[5] = new Color(239,177,0);
      colorHist[6] = new Color(239,49,0);    colorHist[7] = new Color(239,77,166); 
      colorHist[8] = new Color(135,73,169);  colorHist[9] = new Color(135,197,167);
      colorHist[10] = new Color(200,41,2);    colorHist[11] = new Color(239,200,166); 
      colorHist[12] = new Color(100,73,169);  colorHist[13] = new Color(135,100,167);
      colorHist[14] = new Color(135,197,10);  colorHist[15] = new Color(135,197,30);
      colorHist[16] = new Color(165,197,10);  colorHist[17] = new Color(135,197,50);
      colorHist[18] = new Color(185,197,10);  colorHist[19] = new Color(135,197,70);
      colorHist[20] = new Color(215,197,10);  

      indicColor = new Color[4]; 
      indicColor[0] =  new Color(67,71,100);
      indicColor[1] =  new Color(100,71,71);
      indicColor[2] =  new Color(67,95,71);
      indicColor[3] =  new Color(97,70,100);
      
      indicColorlight = new Color[4]; 
      indicColorlight[0] =  new Color(113,120,169);
      indicColorlight[1] =  new Color(184,130,130);
      indicColorlight[2] =  new Color(120,169,126);
      indicColorlight[3] =  new Color(173,125,179);
            
      
      indicMax = new double[4]; indicMin = new double[4]; indicNorm = new double[4];

      acf_series_1 = new double[n_obs];
      acf_series_2 = new double[n_obs];
      cross_12 = new double[n_obs]; 
      cross_21 = new double[n_obs];  

      plot_rolling = new boolean[4]; for(i = 0; i < 4; i++) {plot_rolling[i] = false;}
      plots = new boolean[max_series]; for(i=0;i<max_series;i++){plots[i] = false;}
      var_plots = new boolean[max_series]; for(i=0;i<max_series;i++){var_plots[i] = false;}
      
      addMouseMotionListener(new MyArtMouseMotionListener());  
      addMouseListener(new MyArtMouseListener());      
      
      
    } 

    public void setNames(String[] asset_names)
    {
      
      if(asset_names.length <= 12)
      {names_set = true; this.asset_names = asset_names;}
    }  
    
    public void setDates(String[] dates)
    {
       dates_set = true;       
       this.dates = dates;
    }
    
    public void setStockReturns(double[] rets)
    {
     stock_set = true;
     this.stock_perf = rets;
     computeDataMax();
    }
    
    public void plotStock(boolean t) {plot_stock = t; computeDataMax(); go();}
    
    public void setRollingIndicators(String[] rolls)
    {
      rolling_ind_set = true;
      rollingPerf = new String[rolls.length];
      
      System.arraycopy(rolls,0,rollingPerf,0,rolls.length);
      
      if(rolls.length != n_obs)
      {System.out.println("Something is wrong with rolling length");}
      
      rollingIndics = new double[4][rolls.length];
      computeDataMax();
      go();
    }
    
    public void changeHighlight(int c)
    {hlight_indicator = c; computeDataMax(); go(); }


    public void setplot(int i, boolean sel)
    {plots[i] = sel; computeDataMax(); go();}
    
    public void setvarplot(int i, boolean sel)
    {var_plots[i] = sel; computeDataMax(); go();}
  
    public void setTarget(double[] v)
    {

      target_series = new double[v.length];
      System.arraycopy(v, 0, target_series, 0, v.length); 
      computeDataMax();
      go();

    }
    public void plotTarget(boolean sel) {plot_tar = sel; computeDataMax(); go();}

    public void setSeries(int i, double[] v)
    {  
      sim_series.set(i,v); 
      computeDataMax();
    }

    
    public void plotRollIndic(int i, boolean t) {plot_rolling[i] = t; computeDataMax(); go();}
    
    
    
    public void addSeries(double[] v, String name)
    {
      sim_series.add(v); seriesPlotNames.add(name); computeNobs(); setplot(sim_series.size()-1, true);
      computeDataMax();
    }

    public void setSeries(double[] v, String name)
    {
      sim_series.remove(sim_series.size() - 1);
      sim_series.trimToSize();
      
      seriesPlotNames.remove(seriesPlotNames.size() - 1);
      seriesPlotNames.trimToSize();
      
      sim_series.add(v); 
      seriesPlotNames.add(name); 
      computeNobs(); 
      setplot(sim_series.size()-1, true);
      computeDataMax();
    }    
    
    
    public void replaceSeries(int i, double[] v)
    {
      sim_series.remove(i); 
      sim_series.set(i,v);
    }

    public void addVarSeries(Double[] var, String name)
    {
      var_series.add(var); varPlotNames.add(name); computeDataMax(); computeNobs(); setvarplot(var_series.size()-1, true);
    }
    
    public void setVarSeries(int i, Double[] v, String name)
    {var_series.set(i,v);  varPlotNames.set(i,name); computeDataMax(); computeNobs(); setvarplot(i, true);}    
    
    
    public void computeDataMax()
    {
      int i,k; dataMax = -1000000; dataMin = 10000000; sdataMax = -1000000; sdataMin = 10000000;
      double[] t_series; Double[] T_series;
      double ind0,ind1,ind2,ind3;
      
      sharpeMax = -10000; sharpeMin = 19999; rankMax = -10000; rankMin = 19999;
      ulcerMax = -10000; ulcerMin = 19999; kellyMax = -10000; kellyMin = 19999;
      rankNorm = 0; sharpeNorm = 0; ulcerNorm = 0; kellyNorm = 0;
      dataMaxVar = -100000; dataMinVar = 100000; 
      
      
      for(k=0; k < sim_series.size(); k++)
      {
        if(plots[k])
        {
         t_series = sim_series.get(k);
         n_obs = t_series.length;
         for(i=0; i < n_obs; i++)
         {
            if(t_series[i] > dataMax) dataMax = t_series[i];
            else if(t_series[i] < dataMin) dataMin = t_series[i];   
         }
        }
      } 
      
      for(k=0; k < var_series.size(); k++)
      {
        if(var_plots[k])
        {
         T_series = var_series.get(k);
         n_obs = T_series.length;
         for(i=0; i < n_obs; i++)
         {
            if(T_series[i] > dataMaxVar) dataMaxVar = T_series[i];
            else if(T_series[i] < dataMinVar) dataMinVar = T_series[i];   
         }
        }
      }       
      
      if(stock_set)
      {
        n_obs = stock_perf.length;
        for(i=0; i < n_obs; i++)
         {  
            //System.out.println(stock_perf[i]);
            if(stock_perf[i] > sdataMax) sdataMax = stock_perf[i];
            else if(stock_perf[i] < sdataMin) sdataMin = stock_perf[i];   
         }
      }
      
      if(rolling_ind_set)
      {
        for(i=0;i<rollingPerf.length;i++)
        {
          
          String[] inds = rollingPerf[i].split("[ ]+");
          ind0 = (new Double(inds[0])).doubleValue(); 
          ind1 = (new Double(inds[1])).doubleValue();
          ind2 = (new Double(inds[2])).doubleValue(); 
          ind3 = (new Double(inds[3])).doubleValue();
          
          if(Double.isNaN(ind0)) {ind0 = 0;}
          if(Double.isNaN(ind1)) {ind1 = 0;}
          if(Double.isNaN(ind2)) {ind2 = 0;}
          if(Double.isNaN(ind3)) {ind3 = 0;}
          
          if(ind0 > sharpeMax) {sharpeMax = ind0; indicMax[0] = ind0;}
          if(ind1 > rankMax) {rankMax = ind1; indicMax[1] = ind1;}
          if(ind2 > ulcerMax) {ulcerMax = ind2; indicMax[2] = ind2;}
          if(ind3 > kellyMax) {kellyMax = ind3; indicMax[3] = ind3;}
 
          if(ind0 < sharpeMin) {sharpeMin = ind0; indicMin[0] = ind0;}
          if(ind1 < rankMin) {rankMin = ind1;   indicMin[1] = ind1;}
          if(ind2 < ulcerMin) {ulcerMin = ind2; indicMin[2] = ind2;}
          if(ind3 < kellyMin) {kellyMin = ind3; indicMin[3] = ind3;} 
        
          rollingIndics[0][i] = ind0;
          rollingIndics[1][i] = ind1;
          rollingIndics[2][i] = ind2;
          rollingIndics[3][i] = ind3;
          //System.out.println(rollingIndics[0][i] + " " + rollingIndics[1][i] + " " + rollingIndics[2][i] + " " + rollingIndics[3][i]);
        }
      }
      
      sharpeNorm = Math.abs(sharpeMax - sharpeMin);
      kellyNorm = Math.abs(kellyMax - kellyMin);
      ulcerNorm = Math.abs(ulcerMax - ulcerMin);
      rankNorm = Math.abs(rankMax - rankMin);
      
      indicNorm[0] = Math.abs(indicMax[0] - indicMin[0]);
      indicNorm[1] = Math.abs(indicMax[1] - indicMin[1]);
      indicNorm[2] = Math.abs(indicMax[2] - indicMin[2]);
      indicNorm[3] = Math.abs(indicMax[3] - indicMin[3]);
      
      //System.out.println(indicMax[0] + " " + indicMin[0]);
      sdataNorm =  Math.abs(sdataMax - sdataMin);
      dataNorm = Math.abs(dataMax - dataMin);
      dataNormVar = Math.abs(dataMaxVar - dataMinVar);
    }

    public void setNobs(int n) {n_obs = n;}


    
    public void clearSeries()
    {
      sim_series.clear(); seriesPlotNames.clear(); 
      for(int i=0;i<max_series;i++) {plots[i]=false;} 
      if(var_series.size() == 0) {dates_set = false;}
      
      rolling_ind_set = false;     
      
      go();
    }
   
    public void clearVarSeries()
    {
      var_series.clear();  varPlotNames.clear();
      for(int i=0;i<max_series;i++) {var_plots[i]=false;} 
      if(sim_series.size() == 0) {dates_set = false;}
      go();
    }    
   
    public void computeNobs()
    {
      min_obs = 1000000;
      
      for(int i=0;i<sim_series.size();i++)
      { 
       double[] s = sim_series.get(i);
       if(s.length < min_obs) {min_obs = s.length;}
      }
 
      for(int i=0;i<var_series.size();i++)
      { 
       Double[] s = var_series.get(i);
       if(s.length < min_obs) {min_obs = s.length;}
      }
      
      //System.out.println("Min obs = " + min_obs);
      
    }
   

    public void plotACFCross(boolean t)
    {
      int k;
      plot_acf = t; 
      for(k=0;k<sim_series.size();k++)  {plots[k] = false;}
      go();
    }

    public void plot_acf_series(boolean t1, boolean t2, boolean t3, boolean t4)
    {acf1 = t1; acf2 = t2; cross1 = t3; cross2 = t4; go();}


    public void setACFCross(double[] ts1, double[] ts2, double[] ts3, double[] ts4)
    {
      acf_series_1 = new double[n_obs];
      acf_series_2 = new double[n_obs];
      cross_12 = new double[n_obs]; 
      cross_21 = new double[n_obs];  

      System.arraycopy(ts1, 0, acf_series_1, 0, n_obs);
      System.arraycopy(ts2, 0, acf_series_2, 0, n_obs);
      System.arraycopy(ts3, 0, cross_12, 0, n_obs);
      System.arraycopy(ts4, 0, cross_21, 0, n_obs);        
    }
    

    public void deleteSeries(int i)
    {
       plots[sim_series.size()-1] = false;
       sim_series.remove(i);
       sim_series.trimToSize();  
       go();
    }


    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N,k;
     int t0, t1, x0, x1;
     double[] tseries;
     Double[] vseries;
     String titleBar = ""; 
     //----- computeDataMax -------
     //computeDataMax();
     computeNobs();
          
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;

     int diff_n; 
     BasicStroke orig = new BasicStroke((float)1.3);
     BasicStroke thick = new BasicStroke((float)1.6);

   
       //--------------------- Draw grid with or without forecast points ------------
     float[] dash1 = {10.0f};
     new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);
                                          
     n_obs = min_obs; N = n_obs;
 
     g2d.setStroke(orig);


     for(k=0;k<sim_series.size();k++)
     {
       if(plots[k])
       {
         tseries = sim_series.get(k); diff_n = tseries.length - min_obs;
         g2d.setPaint(colorHist[k]); if(hlight_indicator == k){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((tseries[diff_n + j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((tseries[diff_n + j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
       }
     }
         
     for(k=0;k<var_series.size();k++)
     {
       if(var_plots[k])
       {
         vseries = var_series.get(k); diff_n = vseries.length - min_obs;
         g2d.setPaint(colorHist[6+k]); if(hlight_indicator == k){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((vseries[diff_n + j] - dataMinVar)/(2*dataNormVar))*(double)height);
          x1 = (int)(((vseries[diff_n + j + 1] - dataMinVar)/(2*dataNormVar))*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
       }
     }
   
   g2d.setPaint(myGray);
   if(stock_set && (sim_series.size() > 0 || var_series.size() > 0))
   {
     if(plot_stock)
     {
         diff_n = stock_perf.length - min_obs;
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((stock_perf[diff_n + j] - sdataMin)/(sdataNorm))*(double)height);
          x1 = (int)(((stock_perf[diff_n + j + 1] - sdataMin)/(sdataNorm))*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
     }    
   }
   
   
   if(rolling_ind_set)
   {
     diff_n = rollingPerf.length - min_obs;
     
     for(k=0;k<4;k++)
     {
      if(plot_rolling[k])
      {
       g2d.setPaint(indicColor[k]);
       for(j = 0; j < N-1; j++)
       {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((rollingIndics[k][diff_n + j] - indicMin[k])/(indicNorm[k]))*(double)height);
          x1 = (int)(((rollingIndics[k][diff_n + j + 1] - indicMin[k])/(indicNorm[k]))*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
       }
      }
     } 
   }
   
   
   
   
   if(plot_tracker)
   {
     
    int n_count = 0; int w = 4; int h = 4; 
    if(!plot_acf)
    { 
    
     if(dates_set)
     {
      g.setFont(mono);
      g2d.setPaint(Color.GREEN);    
      g.drawString(date, 5, 15);
     }
     
     
     if(sim_series.size() > 0)
     {
      titleBar = new String("sharpe drwDn rnkRtio sL  tradeIntval  nobs filter_name       value");  
      g.setFont(mono);
      g2d.setPaint(Color.GREEN);    
      g.drawString(titleBar, 5, 30);     
     }
     
     n_count = 1;
     for(k=0;k<sim_series.size();k++)
     {
       if(plots[k])
       {
        n_count++;
        
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t, (height) - track_pos_x[k], w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);     
     
         
        g.setFont(mono);
        //g2d.setPaint(Color.GREEN);
        g2d.setPaint(colorHist[k]);
        if(seriesPlotNames.size()>0) {g.drawString(seriesPlotNames.get(k) + " " + value[k], 5, 15 + n_count*15);}
        else {g.drawString(value[k], 5, 20 + n_count*13);}
 
       }
     } 
     
     if(rolling_ind_set)
     {
       g.setFont(mono);
         
       if(plot_rolling[0]) {  
       g2d.setPaint(indicColorlight[0]);
       g.drawString(" Sharpe = " + indicValue[0] + " ", 90, 15);}
         
       if(plot_rolling[1]) {  
       g2d.setPaint(indicColorlight[1]);
       g.drawString(" Rank = " + indicValue[1] + " ", 190, 15);}
         
       if(plot_rolling[2]) {  
       g2d.setPaint(indicColorlight[2]);
       g.drawString(" Ulcer = " + indicValue[2] + " ", 290, 15);}
         
       if(plot_rolling[3]) {  
       g2d.setPaint(indicColorlight[3]);
       g.drawString(" Kelly = " + indicValue[3] + " ", 390, 15);}
       
       for(k=0;k<4;k++)
       {
        if(plot_rolling[k]) {
        g2d.setPaint(Color.magenta);
        ellipse = new Ellipse2D.Double(track_pos_t, (height) - track_pos_indic[k], w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);}            
       }
     }
     
     for(k=0;k<var_series.size();k++)
     {
       if(var_plots[k])
       {
        n_count++;
        
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t, (height) - track_pos_xvar[k], w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);     
     
         
        g.setFont(mono);
        //g2d.setPaint(Color.GREEN);
        g2d.setPaint(colorHist[6+k]);
        if(varPlotNames.size()>0) {g.drawString(varPlotNames.get(k) + " " + value[k], 5, 15 + n_count*15);}
        else {g.drawString(value[k], 5, 20 + n_count*13);}
 
       }
     }     
     
     
    } 
   }
 
  }
  
     class MyArtMouseMotionListener extends MouseMotionAdapter 
     {
      
      public void mouseDragged(MouseEvent e) 
      { }

      public void mouseMoved(MouseEvent e) 
      {

        int j,k; double[] tseries; Double[] vseries;
        int N = min_obs; date = ""; int diff_n; 
        //System.out.println(e.getX() + " " + width); 
        if(plot_tracker)
        {
            
           j = (int)(((double)(N))*e.getX()/(double)width);  
        
           track_pos_t = (int)(((double)j/(double)(N-1))*(double)width) - 2;
  
           if(rolling_ind_set)
           {
             for(k=0;k<4;k++)
             {
              diff_n = rollingIndics[k].length - min_obs;          
              track_pos_indic[k] = (int)(((rollingIndics[k][diff_n+j] - indicMin[k])/indicNorm[k])*(double)height) + 2; 
              indicValue[k] = df2.format(rollingIndics[k][diff_n+j]);                         
             }
           }
  
  
           for(k=0;k<sim_series.size();k++)
           {
            if(plots[k])
            {
              tseries = sim_series.get(k);  diff_n = tseries.length - min_obs;          
              track_pos_x[k] = (int)(((tseries[diff_n+j] - dataMin)/dataNorm)*(double)height) + 2; 
              value[k] = df2.format(tseries[diff_n+j]);
            }   
           }    
           
           for(k=0;k<var_series.size();k++)
           {
            if(var_plots[k])
            {
              vseries = var_series.get(k); diff_n = vseries.length - min_obs;           
              track_pos_xvar[k] = (int)(((vseries[diff_n+j] - dataMinVar)/(2*dataNormVar))*(double)height) + 2; 
              value[k] = df2.format(vseries[diff_n+j]);
            }   
           }            
           
           if(dates_set) {diff_n = dates.length - min_obs; date = dates[diff_n+j];}// System.out.println(diff_n + " " + dates.length + " " + min_obs + " " + j);}
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
