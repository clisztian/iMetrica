package ch.imetrica.mdfaTradingStrategies;

/*!
Copyright (C) 2016 Christian D. Blakely
This file is part of iMetrica, a free-software/open-source application
for interactive graphical econometric analysis - http://imetricablog.com/

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

import java.io.*;
import java.util.*;
import java.text.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import org.joda.time.IllegalInstantException;
import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;

import ch.imetrica.bayesCronos.Cronos;
import ch.imetrica.mdfa.IMDFA;


public class MDFAStrategyEvolution
{


    int n_obs, n_rep, K, L, lag,S;
    int i1, i2, dd, DD, K1, output, flength;  
    boolean mdfa_iter, dfa_iter, iter, imdfaReg, short_sell;
    boolean long_buy = true;
    public double[] Gamma; 
    
    String sampMin = "00";
    //---- extra Filter controls ----------------
    double lambda, expweight, cutoff, cutoff0, cutoff1, lambda_3;   
    double smooth, decay, cross, decay2;
    ArrayList<String> expvar;
    int[] histo_stat;
    int bad_starts = 0;
    int trade_succ_ratio = 0;
    boolean diff_sig_trading = false;
    double mean_perf, rank_perf; 
    int num_full_positive_returns = 0;
    //--- adaptive parameters
    int minute_df;
    int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
    double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;    
    boolean exp_var; 
    double intra_delta = 0.019;
    boolean forex24 = false; 
    boolean log_prices = false;
    boolean futures = false;
    boolean printall = false;
    boolean ronin_print = false;
    double phi1,phi2,mu, sigma_2t, sigma_3t, sigma_t;
    DecimalFormat df4 = new DecimalFormat("##.#####");
     //rt = new double[le];
    double[] vol;
    double[] vol2;
    double[] vol3;    
    boolean first_trade_loss = false;
    IMDFA mdfa;
    IMDFA forecast;
    IMDFA mdfa_h0;
    
    
    double bandCutoff = 0;
    boolean allPass = false;
    double kurtosis; 
    double skewness; 
    boolean simulated_data = false;
    int max_delay = 0;
    boolean filter_computed = false;
    boolean print_coeffs = false;
    boolean gaussianized=false;
    double ibspread = .00015;
    boolean pip_take_mode = false;
    int L_white = 100;
    double spreadplus = 0.0;
    double mean_rev_amnt = .0000;
    boolean ln_trans = true;
    public ArrayList<String> forkedFiles; 
    public String[] dataFiles;
    public DecimalFormat df2 = new DecimalFormat("##.##");
    public double[] tseries;     
    public double[] b_coeffs;
    public double[] b_old;
    public double[] b_update; 
    public double[] b_copy;
    public double[] b_lag;
    public double[] h0b0;
    public double[] b_avg;
    public double[] xt;
    public double[] gauss_rt1;
    public double[] gauss_rt2;
    public double[] gauss_rt3;
    public boolean futures_data = false;
    public double win_ratio;
    public String final_trade_time;
    public int day_count = 0; 
    public int recompute_day = 0;
    public double[] ret_dist;
    public int[] pos_ret_dist,neg_trades_started;
    public int[] neg_ret_dist,pos_trades_started;
    public double[] diff_account,pos_ret_mean_time,neg_ret_mean_time; 
    public double[] neg_trades_started_mean,pos_trades_started_mean;
    public int n_files = 0;
    public int cur_hour, start_hour;
    public boolean trading_closed = false;
    public int trade_started;
    public double min_pnl_dd, max_pnl_uu;
    public int n_pos_ret,n_neg_ret;
    public double neg_ret_mean,pos_ret_mean;
    public boolean sig_inverse,classical_data=false;
    public ArrayList<StrategyParameters> filter_strategies;
    public String[] asset_name; 
    public ArrayList<Double> crit_0,ask,mid,bid;
    public ArrayList<Double> crit_1;
    public ArrayList<double[]> full_returns_array;
    public ArrayList<double[]> morning_returns;
    public ArrayList<Double> deg_0;
    public ArrayList<Double> deg_1; 
    public String[] dates_price;
    public boolean its_closing_time, friday_closing;
    public ArrayList<StrategyParameters> strategies = new ArrayList<StrategyParameters>();
    public ArrayList<String> allsignal = new ArrayList<String>();
    public ArrayList<Double> perf_returns = new ArrayList<Double>();
    public ArrayList<Double> vol_0;
    public ArrayList<Double> vol_1;
    public ArrayList<Filter> filters; 
    public ArrayList<Double> max_ranks;
    public double[] fore_gamma;
    public ArrayList<Double> interp;
    public ArrayList<Double> maxIntValue;
    public double sharpeRatio;
    public boolean useH0;
    public double daily_return;
    public double log_price;
    public double total_ROI = 0.0;
    public int total_succ = 0;
    public int total = 0;  
    public double[] fore_coeffs;
    public double final_price;
    public double diff_band;
    public int avg_count = 0;
    public int update_i1=0;
    public int update_i2=0;
    public int L_update; //--- the update filter length
    public int N_update; //--- the total length of new data (N = number of out-of_sample points + L)
    public double[] update_signal;
    public int K_update; //--- the update filter K

    public boolean time_filter = true;
     
    public String start_year = "2011-"; 
    public double kellyPerc = 0;
    public double ulcer_index = 0;
    public boolean reset_signal = true;
    public int last_trade_index;
    public int end_time_index;
    public boolean morning_buy = false;
    public double[] sub_price;
    public double[] old_tseries,out_series;
    public double[] xf; 
    public double[] xf_turn_val;
    public double[] pnl,png;
    public double[] sub_account;
    public double[] sub_signal;
    
    double tradecost = 0;
  
    int n_stages=0;
    
    double execution_delay = 0;
    StrategyParameters strategy;
    boolean ib_data; 
    boolean asian_close = false;
    int min_freq,n_freq;
    int fore_obs;      
    double high_cut = 1.46;
    double fore_lambda = 0;
    double fore_exp = 12.9;
    double fore_smooth = 0.545;
    double fore_decay = 0.198;
    double fore_decay2 = 0.644;
    double fore_cross = 0;
    int Lfore = 126;  
    int fore_lag = -1;    
    
    double criteria;
    double degrees;     
    double MDFAmin;
    double avg_n_trades;
    double avg_rate; 
    int close_time = 0;
    int start_seq;
    int start_index;
    int end_index;
    int total_obs;
    int recomp_pulse;
    int pulse_count = 0;
    int adaptive_recomp_pulse = 0;
    int adapt_pulse_count = 0;
    double avg_vol;
    double vola_thresh;
    boolean volatility_thresh;
    ArrayList<Double> b0_trend;
    ArrayList<Double> final_trades;
    ArrayList<Double> dailyoutret; 
    ArrayList<String> crits,svm;
    ArrayList<String> date_returns; 
    String[] string_trades;
    String adaptfilterParamFile;
    String dataFile; 
    int n_toks; int n_series; int count; 
    String delims = "[,]+";
    String[] tokens; 
    String[] date_tokens;
    String date_delims = "[ ]+";
    String startingTime = "09:30";   
    String endingTime = "15:45";
    String insampStart; 
    String permEndTrading;
    String optimizeTime = "12:30";  //time to check for optimal parameter
    FileInputStream fin; DataInputStream din; BufferedReader br; String names; 

    boolean min30trade = false;
    int spread_method = 0;
    double maxdraw;
    int profit_baby = 0;
    boolean morning_optimize = false;
    double ROI,longROI,ratio,shortROI;
    double max_drop;
    int succ_trades;
    int total_trades;
    double drop_down;
    double[] mnmx;
    double volume, open_ask, open_bid, end_ask, end_bid, high_ask, high_bid, low_ask, low_bid, close;
    int trade_obs;
    double volatility;
    double criteria_h0;
    double degrees_h0;      
    double MDFAmin_h0;      
    boolean stop_loss = false;
    double stop_loss_thresh = 0;
    double yesterday_vol;
    double max_interp_value=0;
    double max_ret_interp = 0.0;
    boolean H0set = false;
    int h0_outsamp; 
    double localdecay;
    double h0_delta;
    double interp_start=0.0;
    ArrayList<Double> interp_vals;
    double trading_cost = 0.0;
    double tradeDelta = 0.0;;
    double[] lo_prix,hi_prix;
    double fridayROI = 0; 
    int fridayROI_pos = 0; 
    int fridays = 0;
    boolean jpy = false;
    String final_time; 
    ArrayList<String> aux_explanatory;
    double t1=0,t2=0;
    boolean gaussianize = false;
    double max_peak = 0; double peak_loc = 0;
    int last_trade;
    double final_trade; 
    int losses_in_arow = 0;
    double current_signal=0;
    double lag_signal=0;
    boolean recomputeH0;
    int num_trade_days;
    double mean = 0.0;
    int num_gains = 0; int num_losses = 0;
    
    boolean spread_on = false;
    boolean binary_rule, signal_strength_rule,downtick_strategy,signal_profit;
    boolean take_profit;
    double take_profit_thresh;
    String finalCall;
    double tradingCost = 0.0000;
    boolean reg_trading_hours = false;
    boolean file_valid;
    boolean cont_lookback = false;
    boolean adaptive_filtering = false;
    boolean compute_lag = false;
    boolean save_volatility = true;
    double[] yesterday_series;
    double[] actual_signal;
    double[] actual_price,actual_loprice,actual_hiprice;
    double risk,reward;
    boolean set_expseries = false;
    String explanatorySeries;    
    boolean let_shop = false;
    
    int filtered_data_method = 0;
    int out_transaction; 
    int in_transaction;
    boolean red_zone; 
    double global_stop_loss;
    double profitable_stop;
    
    DateTime liquidateTimeDT, startTimeDT, endTimeDT;
    
    
  
    ArrayList<String> imetrica_report = new ArrayList<String>(); 
    ArrayList<Double> daily_signal = new ArrayList<Double>();
    ArrayList<Double> daily_price = new ArrayList<Double>();
    ArrayList<Double> daily_returns = new ArrayList<Double>();
    ArrayList<String> daily_dates = new ArrayList<String>(); 
    
    ArrayList<Double> gauss_target = new ArrayList<Double>();
    ArrayList<Double> gauss_1 = new ArrayList<Double>();
    ArrayList<Double> gauss_2 = new ArrayList<Double>();

    ArrayList<Integer> last_trades = new ArrayList<Integer>();
    ArrayList<Double> avg_volatility = new ArrayList<Double>();
    ArrayList<Double> close_series = new ArrayList<Double>();
    ArrayList<Double> highlow_series = new ArrayList<Double>();    
    ArrayList<Double> exp_series_1 = new ArrayList<Double>();  
    ArrayList<Double> exp_series_2 = new ArrayList<Double>();
    ArrayList<Double> price = new ArrayList<Double>();    
    ArrayList<Double> lo_price = new ArrayList<Double>();   
    ArrayList<Double> hi_price = new ArrayList<Double>();   
    ArrayList<String> dates_series = new ArrayList<String>();       
    ArrayList<Double> sub_returns = new ArrayList<Double>();
    ArrayList<Double> live_series = new ArrayList<Double>();
    ArrayList<String> trade_days = new ArrayList<String>();
    ArrayList<Double> returns = new ArrayList<Double>();
    ArrayList<Double> longreturns = new ArrayList<Double>();
    ArrayList<Double> shortreturns = new ArrayList<Double>();
    ArrayList<Double> dropdowns = new ArrayList<Double>();
    ArrayList<Double> success = new ArrayList<Double>();
    ArrayList<String> dates_low_high = new ArrayList<String>();  
    ArrayList<String> dailyReport = new ArrayList<String>(); 
    PrintWriter spread;
    NumberFormat formatter,formatter2,formatter3;
    int ndays;
//     int rolling_length = 30;
//     ArrayList<String> rolling_ind;
    
    boolean daily_strategy = false;
    ArrayList<MDFATrade> mdfaTrades;
    boolean change_time_zone = false;
    double[] account;
    double[] signal;
    double[] prix;
    double[] lag_signals;
    double[] dreturns;
    double[] trades;
    double[] morningSeries;
    int insamp_start_int; 
    ArrayList<Double> lookback_returns = new ArrayList<Double>();
    int num_pos_returns;
    boolean print_debug = false;
    boolean friday = false;
    double standard_deviation;
    double sum_sd;
    double zero_ret;
    boolean lookback_ready = false;
    int nfiles = 0;
    boolean trading_hours;
    double rank_coeff;
    boolean recompute_filter = true;
    boolean bid_ask_data = true;
    int trade_count;
    String date;
    String filterParamFile;
    
    int n_out_samp = 0;
    class ibHash extends HashMap<String,String> {
	private static final long serialVersionUID = 1L; }
    ibHash ib_data_hash;
    String ib_data_file;
    DateTimeFormatter fmt;
    
    
    public MDFAStrategyEvolution(String file, String dataf)
    {
      filterParamFile = new String(file); dataFiles = new String[1];        
      dataFiles[0] = new String(dataf);   
      nfiles = 1;
    }
    
    public MDFAStrategyEvolution(int _nfiles, String file)
    {
      
      filterParamFile = new String(file);
      dataFiles = new String[6];
      

      dataFiles[0] = "lastSeriesData_0.dat";
      dataFiles[1] = "lastSeriesData_1.dat";
      dataFiles[2] = "lastSeriesData_2.dat";      
      dataFiles[3] = "lastSeriesData_3.dat";
      dataFiles[4] = "lastSeriesData_4.dat";      
      dataFiles[5] = "lastSeriesData_5.dat";    
    
   
      nfiles = _nfiles;
    }
    
    
    public MDFAStrategyEvolution(String file, FilterParameters param)
    {
      dataFiles = new String[1];
      dataFiles[0] = file;       
      nfiles = 0;  
      int n; int k;
 
      
      min_freq = 15; n_freq = 4; 
      n_obs = param.n_obs;
      n_rep = param.n_rep;
      n = n_obs;
      K = (int)(n/2); double om;
      K1 = K+1;  
      trade_obs = 26;
      L = param.L; 
      lag = param.lag;
      smooth = param.smooth; 
      decay = param.decay;
      decay2 = param.decay2;
      cross = param.cross;
      expweight = param.expweight;
      lambda = param.lambda;
      i1 = 0; 
      cutoff0 = 0; 
      cutoff1 = param.cutoff;
     
      mdfa = new IMDFA(n_obs, n_rep, L, 0, 0, 0, .52, 0, 1);
      mdfa.set_L(L); 
      mdfa.set_nobs(n_obs);  
      mdfa.set_nreps(n_rep);
      mdfa.set_lag(lag);
  
       mdfa.setRegularization(smooth, decay, decay2, cross);
       mdfa.set_dd(0);
       mdfa.set_DD(0);
       mdfa.set_bconstraints(i1, 0); 
       mdfa.set_lambda(lambda); 
       mdfa.set_exp(expweight); 
       mdfa.set_cutoff0(cutoff0);
       mdfa.set_cutoff(cutoff1);
       mdfa.set_mdfa(true);
  
       Gamma = new double[K1];       
       for(k=0; k<=K;k++)
       {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
       }       
       mdfa.set_Gamma(Gamma);              
    }    
    
    public void togglePrint(boolean t) {print_debug = t;}
    
    public void setCutoff(double cut)
    {
       double om;
       cutoff1 = cut; 
       mdfa.set_cutoff(cutoff1);
       mdfa.set_mdfa(true);
  
       Gamma = new double[K1];       
       for(int k=0; k<=K;k++)
       {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
       }       
       mdfa.set_Gamma(Gamma);       
    }    
    
    
    public void uploadInterpParams(String fname)
    {
       
       
       
       interp = new ArrayList<Double>();
       Double D; 
       String strline; String[] tokens; String delims = "[ ]+"; int n_toks; 
       try
       {
          
         FileInputStream fin = new FileInputStream(new File(fname));
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 

         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
        
           //System.out.println(tokens[0] + " " + tokens[1]);
         
           n_toks = tokens.length; 
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           D = new Double(tokens[0]); //take only the value, no dates
           //val = D.doubleValue();
           interp.add(D);
          
          }
          System.out.println("Past interp values set");
          br.close();
       } 
       catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
       catch(IOException ioe){System.out.println("IO out error..." + ioe);}   
  
    }
    
    public void setDailyStrategy(int nda)
    {
      daily_strategy = true;
      ndays = nda;
    }
    
    
    public void setFinalTradeTime(String f) {final_trade_time = new String(f);}
    
    public void saveVolatility(boolean v) {save_volatility = v;}
    
    public void setAdaptiveFiltering(boolean t) {adaptive_filtering = t;}
    
    public void setStartIndex(int c) {start_index = c;}
    
    public void setParameterFile(String file)
    {filterParamFile = new String(file);}
    
    public void setAdaptiveParameterFile(String file)
    {adaptfilterParamFile = new String(file); readAdaptiveFilterParams();}
    
    public void setTradingParameters(boolean shortse, boolean longt, double trans)
    {short_sell = shortse; long_buy = longt; tradingCost = trans;}
        
    public void setLagLookback(boolean l, boolean comp)
    {cont_lookback = l; compute_lag = comp;}
    
    public void turnOnH0Filter(boolean f) {useH0 = f;}
    
    public void setH0Variables(int _h0_outsamp, double _localdecay, double _h0_delta, double int_start, double intra_del)
    {
      h0_outsamp = _h0_outsamp;
      localdecay = _localdecay;
      h0_delta = _h0_delta;
      interp_start = int_start;
      intra_delta = intra_del;
    }
    
    public void turnOnGaussianize(boolean f) {gaussianize = f;}
    
    public void turnOnVolatilityThresh(boolean f) {volatility_thresh = f;}
    
    public void setNobs(int n)
    {
      int k;
      n_obs = n; 
      mdfa.set_nobs(n_obs);
    
             
      K = (int)(n/2); double om;
      K1 = K+1;  
      trade_obs = 26;

      Gamma = new double[K1];       
      for(k=0; k<=K;k++)
      {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
         
      }       
      mdfa.set_Gamma(Gamma);
//       for(k=0;k<Gamma.length;k++)
//       {
//        System.out.println(Gamma[k]);
//       }
    }
    

    
    public void setStopLoss(double stp) 
    {
       stop_loss = false;
       stop_loss_thresh = stp;
       if(stop_loss_thresh > 0) {stop_loss = true;}
    }       
    
    
    public void setOneStepOptimization(double v) {mdfa.setOnestep(v);}
    public void setOnestepDiff(double v) {mdfa.setDiffOnestep(v);}
    public void seti2Phase(int i, double ph) {mdfa.set_shift(ph); mdfa.set_bconstraints(0,i);}
    
    
    public void setStartingTradeTime(int c)
    {start_seq = c + 4 - insamp_start_int;}    
    
    public void setStrategyParameters(StrategyParameters sp, String data_file, double tradcost)
    {
        //the basic crap    
        turnOnVolatilityThresh(false); saveVolatility(true);
        setAdaptiveFiltering(false); //setAdaptiveParameterFile("child_filter.params");
        setAdaptiveRecompPulse(31); setLagLookback(false,false);
        setTradingParameters(true, true, tradcost); setStartIndex(0); setRecompPulse(0);     
        
        //set filter parameters 
        n_obs = sp.n_obs;  n_rep = sp.n_rep; L = sp.L;
        
        K = (int)(n_obs/2); double om; K1 = K+1; int k; 
        mdfa = new IMDFA(n_obs, n_rep, sp.L, 0, 0, 0, sp.cutoff, 0, 1);
        mdfa.set_L(sp.L);  
        mdfa.set_nobs(n_obs);  
        mdfa.set_nreps(n_rep);
        mdfa.set_lag(sp.lag);
  
        mdfa.setRegularization(sp.smooth, sp.decay, sp.decay2, sp.cross);
        mdfa.set_dd(0);
        mdfa.set_DD(0);
        mdfa.set_bconstraints(sp.i1, sp.i2); 
        mdfa.set_lambda(sp.lambda); 
        mdfa.set_exp(sp.expweight); 
        mdfa.set_cutoff0(sp.cutoff0);
        mdfa.set_cutoff(sp.cutoff);
        mdfa.set_mdfa(true); 
        diff_sig_trading = false;
        sig_inverse = false;    
        stop_loss = false;
        
        if(sp.b0trend == 1) {mdfa.useH0Set(true);}
        else if(sp.b0trend == 2) {reset_signal = true; classical_data = false;}
        else if(sp.b0trend == 3) {reset_signal = true; classical_data = true;}
        else {mdfa.useH0Set(false); reset_signal = false; classical_data = false;}
       
        //if(sp.sig_diff > 0) {diff_sig_trading = true;} 
        if(sp.sig_inv == 1) {sig_inverse = true; diff_sig_trading = false;}
        else if(sp.sig_inv == 2) {sig_inverse = false; diff_sig_trading = true;}
        else if(sp.sig_inv == 3) {sig_inverse = true; diff_sig_trading = true;}
        else {sig_inverse = false; diff_sig_trading = false;}
        
        
//         if(sp.stop_loss < 3)
//         {
//           stop_loss_thresh = 0;
//           spread_method = sp.stop_loss;
//         }
//         else
//         {
//           stop_loss_thresh = .0003*sp.stop_loss;
//           if(stop_loss_thresh > 0) {stop_loss = true;}

          take_profit_thresh = .0001*sp.stop_loss;
          if(take_profit_thresh > 0) {take_profit = true;}
          
          if(sp.take_profit > 0)
          {
            setStopLossInt(sp.take_profit); 
          }
          
          
          if(n_rep <= 2) //use cross instead for spread_
          {spread_method = (int)sp.cross;}
          else
          {spread_method = 0;}
          
       // }
        
        Gamma = new double[K1];       
        for(k=0; k<=K;k++)
        {       
         om = (k*Math.PI/K);
         if(om < sp.cutoff0) {Gamma[k] = 0.0;}
         else if(om >= sp.cutoff0 && om <= sp.cutoff)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
        }       
        mdfa.set_Gamma(Gamma);      
                  
        setOneStepOptimization(sp.hybrid_weight);
        setOnestepDiff(sp.hybrid_weight_diff);
        seti2Phase(sp.i2,sp.time_shift);
/*        setStartingTime(sp.insampStart);
        setStartingTradeTime(sp.startTradingInt);
        setEndingTime(sp.endTrading); */     
        

        
        
        forex24 = false;
        setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time);
        
        if(sp.insampStart.equals("03:30") || sp.insampStart.equals("04:00"))
        {
           futures = true;
           setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time); 
        }    
        else if(sp.insampStart.equals("00:00:00") || sp.insampStart.equals("01:00"))
        {
         forex24 = true;
         setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time);
        }
        else if(sp.insampStart.equals("00:15:00") || sp.insampStart.equals("00:30:00"))
        {
         forex24 = true;
         setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time);
        }        
        if(sp.endTrading.equals("16:00") || sp.endTrading.equals("16:15")) 
        {
         forex24 = false;
         setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time); 
        }   
        forex24 = true;
        setTimeInterval(sp.insampStart, sp.endTrading, final_trade_time);
        
         
        setFirstTradeTime(sp.startTrading);       
        permEndTrading = new String(sp.endTrading);
        
        
     
        dataFiles[0] = data_file; //"lastSeriesData_13.dat"; //dataFiles[i];
 
    }
    
    public void setL(int l)
    {
      mdfa.set_L(l);            
    }    
    
    public void setRecomputeDay(int r) {recompute_day = r;}
    
    public void setMDFAParameters()
    {
 
       int n; int k;

       n_obs = 408;

       //n_rep = 3;       
       n_rep = 4;
       n = n_obs;
       K = (int)(n/2); double om;
       K1 = K+1;  
       trade_obs = 26;
       
     
       mdfa = new IMDFA(n_obs, n_rep, L, 0, 0, 0, .52, 0, 1);
    
       mdfa.set_L(L); 
       mdfa.set_nobs(n_obs);  
       mdfa.set_nreps(n_rep);
       mdfa.set_lag(lag);
  
       mdfa.setRegularization(smooth, decay, decay2, cross);
       mdfa.set_dd(0);
       mdfa.set_DD(0);
       mdfa.set_bconstraints(i1, 0); 
       mdfa.set_lambda(lambda); 
       mdfa.set_exp(expweight); 
       mdfa.set_cutoff0(cutoff0);
       mdfa.set_cutoff(cutoff1);
       mdfa.set_mdfa(true);
  
       Gamma = new double[K1];       
       for(k=0; k<=K;k++)
       {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
       }       
       mdfa.set_Gamma(Gamma);      
          
          
          
    }
    
    public void setStopLossInt(int stp) 
    {
       stop_loss = false;
       stop_loss_thresh = .0001*stp;
       if(stop_loss_thresh > 0) {stop_loss = true;}
    }    
    
    public void setStopLoss(int stp) 
    {
       stop_loss = false;
       stop_loss_thresh = .0001*stp;
       if(stop_loss_thresh > 0) {stop_loss = true;}
    }     
    
    
    public void setH0Parameters(double _lambda, double _expweight, double _smooth, double _decay1, double _decay2, double _cross)
    {
       lambda = _lambda; expweight = _expweight; smooth = _smooth; cross = _cross;
       
       mdfa.set_lambda(lambda); 
       mdfa.set_exp(expweight); 
       mdfa.setRegularization(smooth, decay, decay2, cross);
       
    }
    
    public void setTrendH0()
    {
      mdfa.setTrendH0();  //sets the trend H0 in the mdfa object
    }
    
    public void setRecomputeH0(boolean recomp2)
    {recomputeH0 = recomp2;}
    
    public void setRecompPulse(int p)
    {recomp_pulse = p;}
    
    public void setAdaptiveRecompPulse(int p)
    {adaptive_recomp_pulse = p;}
    
    
    public void setStartingTime(String time)
    {
       startingTime = new String(time);
       if(time.equals("08:30")) {trade_obs = 30; insamp_start_int = 0;}
       else if(time.equals("08:45")) {trade_obs = 29; insamp_start_int = 1;}
       else if(time.equals("09:00")) {trade_obs = 28; insamp_start_int = 2;}
       else if(time.equals("09:15")) {trade_obs = 27; insamp_start_int = 3;}
       else if(time.equals("09:30")) {trade_obs = 26; insamp_start_int = 4;}
       else if(time.equals("09:45")) {trade_obs = 25; insamp_start_int = 5;}
       else if(time.equals("10:00")) {trade_obs = 24; insamp_start_int = 6;}
       else if(time.equals("10:15")) {trade_obs = 23; insamp_start_int = 7;}
    }
    
    public void setEndingTime(String time) //relative to normal ending time at 15:45
    {
       endingTime = new String(time);
       if(time.equals("16:00")) {trade_obs = trade_obs + 1;}
       else if(time.equals("15:30")) {trade_obs = trade_obs-1;}
       else if(time.equals("15:15")) {trade_obs = trade_obs-2;}
       else if(time.equals("15:00")) {trade_obs = trade_obs-3;}
       else if(time.equals("14:45")) {trade_obs = trade_obs-4;}
       else if(time.equals("14:30")) {trade_obs = trade_obs-5;}
       else if(time.equals("14:15")) {trade_obs = trade_obs-6;}
       else if(time.equals("14:00")) {trade_obs = trade_obs-7;}
       else if(time.equals("13:45")) {trade_obs = trade_obs-8;}
       else if(time.equals("13:30")) {trade_obs = trade_obs-9;}
       else if(time.equals("13:15")) {trade_obs = trade_obs-10;}
       else if(time.equals("13:00")) {trade_obs = trade_obs-11;}
    }    
    
    public void setObservationFrequency(int min)
    {
     if(min > 0)
     {
      min_freq = min;
      n_freq = (int)60/min;
     }
    }    
    
    
    //--------------------------------------------------------------------------------------------------------
    public void setTimeInterval(String starttime, String endtime, String finalC)
    {
    
        insampStart = new String(starttime);
        endingTime = new String(endtime);
        finalCall = new String(finalC);
        
        String delimsint = "[:]+";
        
        String[] start_tokens = starttime.split(delimsint);
        String[] end_tokens = finalCall.split(delimsint);
    
        int hourDiff = (new Integer(end_tokens[0])).intValue() - (new Integer(start_tokens[0])).intValue(); //6
        
        
        if(hourDiff > 0)
        {
         //the number of trade obs accounting for hours, assumes at least one hour trading 
         trade_obs = n_freq*hourDiff+1; //25
         
         //now add any additional minutes from end         
         trade_obs = trade_obs + (int)(new Integer(end_tokens[1])).intValue()/min_freq; //28
         
         //now subtract any additional minutes from start
         trade_obs = trade_obs - (int)(new Integer(start_tokens[1])).intValue()/min_freq; //26
         
        }
        //System.out.println("TRADE OBS = " + trade_obs + " " + finalCall);
    }
    

    
    
    //first trade of day, should be equal to or greater than startingTime
    public void setFirstTradeTime(String first)
    {
       String delimsint = "[:]+";
       startingTime = new String(first);
       
       String[] first_tokens = first.split(delimsint);
       String[] start_tokens = insampStart.split(delimsint);
      
      
       int hourDiff = (new Integer(first_tokens[0])).intValue() - (new Integer(start_tokens[0])).intValue(); //0
       
       start_seq = n_freq*hourDiff;
       start_seq = start_seq + (int)(new Integer(first_tokens[1])).intValue()/min_freq;
       start_seq = start_seq - (int)(new Integer(start_tokens[1])).intValue()/min_freq;
     
    }
    
    //------------------------------------------------------------------------------------------------------------
  
    
    
    public static ArrayList<String> addExplanatorySeries(String name)
    { 

      String strline; 
      ArrayList<String> aux_explanatory = new ArrayList<String>();
      
      try
      {
          
         FileInputStream fin = new FileInputStream(new File(name));
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           aux_explanatory.add(strline);
         }  
         br.close();     
      }
      catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
      catch(IOException ioe){System.out.println("IO out error..." + ioe);}   
      
      return aux_explanatory; 
    }
    
    
    public void setTimeStandards(File file)
    {
    
          String strline; 
       int n_toks;   
       String delims = "[,]+";
       String[] tokens; 
       int n_intervals = 0; 
       new ArrayList<String>();
 
      new String(" ");  
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
       
       
       //String[] datetime = tokens[0].split(date_delims);
       //System.out.println(tokens[0] + " " + insampStart);
       
       if(tokens[0].indexOf(startingTime) != -1) //date_stamp.indexOf(endingTime) != -1
       {n_intervals = 0;}
       else if(tokens[0].indexOf(endingTime) != -1)
       {last_trade_index = n_intervals + 1;}
       else
       {n_intervals++;} //System.out.println(tokens[0] + " " + insampStart);}
         
       }  
       br.close();
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);} 
     
      if(print_debug) System.out.println("ENDING TIME and final index = " + endingTime + " " + last_trade_index);
    }
     
    
    
    public boolean startStrategyDaily()
    {
      int j,i,l,N,file_count,k;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      int nobs_count=0;
      String date_stamp,strline; 
      date_returns = new ArrayList<String>();
      double mean = 0; int dtotal_trades = 0;
      trade_obs = 100;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0;
      dailyoutret = new ArrayList<Double>();
      mdfaTrades = new ArrayList<MDFATrade>();
      
      ArrayList<String> datesAll = new ArrayList<String>();
            last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;

      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      date_stamp = "";
      ndays = 30;
      
      stop_loss = true;
      stop_loss_thresh = 50; 
      
      
      ArrayList<String> dailyData = new ArrayList<String>();
      start_seq = 0;
      
        for(file_count=0;file_count<1;file_count++)
        {
         
         dailyData = new ArrayList<String>();
         
         if(dataFiles[file_count].indexOf("JPY") != -1 && dataFiles[file_count].indexOf("SEKJPY") == -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
         }
         
         //setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         
         try{
         
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         //spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; nobs_count = 0; computed = false;
         
         trade_count = 0;
         recompute_filter = true;
         
         while((strline = br.readLine()) != null)
         {dailyData.add(strline);}
         
         trade_obs = dailyData.size();
         
         
         
         signal = new double[trade_obs];
         xt = new double[trade_obs];
         lag_signals = new double[trade_obs];
         prix = new double[trade_obs];
         lo_prix = new double[trade_obs];
         hi_prix = new double[trade_obs];
         ret_dist = new double[trade_obs];
         pos_ret_dist = new int[trade_obs];
         neg_ret_dist = new int[trade_obs];
         neg_trades_started = new int[trade_obs];
         pos_trades_started = new int[trade_obs];
         neg_trades_started_mean = new double[trade_obs];
         pos_trades_started_mean = new double[trade_obs];      
         diff_account = new double[trade_obs];
         pos_ret_mean_time = new double[trade_obs];
         neg_ret_mean_time = new double[trade_obs];         
      
         
         
         for(int day = 0; day < dailyData.size(); day++)
         {
          
          strline = dailyData.get(day);
         
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           

           
           
           
          date_stamp = tokens[0];              
          String[] intdates = date_stamp.split(ddelims);     
          new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
         
          
      
      
      
      
      
      
          if(bid_ask_data)
          {
            
            D = new Double(tokens[4]); close_series.add(D);
           
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
            }
            else
            {
             bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
            }

             if(current_signal > 0)
             {
              D = new Double(tokens[2]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(bid.get(bid.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(bid.get(bid.size()-1)) - log(ask.get(ask.size()-2)));} //high urgency
             }
             else if(current_signal < 0)
             {
              D = new Double(tokens[3]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(ask.get(ask.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(ask.get(ask.size()-1)) - log(bid.get(bid.size()-2)));} //high urgency
             }
             else
             {
              D = new Double(tokens[1]); price.add(D); 
              
              D = new Double(tokens[4]); 
              if(ib_data && ib_data_hash.containsKey(tokens[0]))
              {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
              else
              {live_series.add(new Double(tokens[4]));}
             }
             lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));
             
             //System.out.println(date_stamp + " " + live_series.get(live_series.size() - 1));
//              if(live_series.size() > 2)
//              {
//               System.out.println("\n" + date_stamp + " " + live_series.get(live_series.size() - 1) + ", ask = " + ask.get(ask.size()-1) + ", bid = " + bid.get(bid.size()-1) + 
//               ", ask_t-1 " + ask.get(ask.size()-2) + ", bid_t-1 " + bid.get(bid.size()-2));
//              }
             //spread.println(date_stamp + " " + (new Double(tokens[3]) - new Double(tokens[2])));
                   
          }
          else
          {          
          
           D = new Double(tokens[1]); price.add(D);
           D = new Double(tokens[2]); close_series.add(D);
          
          //if(print_debug)System.out.println(date_stamp + " " + price.get(price.size()-1) + " " + close_series.get(close_series.size()-1));
          //System.out.println(date_stamp + " " + price.get(price.size()-1) + " " + close_series.get(close_series.size()-1));
           if(n_rep > 2 && tokens.length > 3)
           {
            D = new Double(tokens[3]); exp_series_1.add(D); 
           }
           if(n_rep > 3 && tokens.length > 4)
           {
            D = new Double(tokens[4]); exp_series_2.add(D);               
           }           
          }
          
          
          
          nobs_count++;
          
          //System.out.println(nobs_count + " " + n_obs);  
          if(nobs_count >= n_obs) 
          {
          
             datesAll.add(date_stamp);
             tseries = new double[n_rep*n_obs];
             out_series = new double[n_obs];  
                    

                   if(classical_data)
                   {
                    for(i=0;i<n_obs;i++)
                    {
                     tseries[n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }      
                    
                    for(i=0;i<n_obs;i++)
                    {
                     out_series[n_obs-1-i] = tseries[n_obs-1-i];
                    }                    
                    
                   }
                   else
                   {
                    for(i=0;i<n_obs;i++)
                    {
                     tseries[n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }
                    
                    for(i=0;i<n_out_samp;i++)
                    {
                     out_series[n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                    }
                    for(i=n_out_samp;i<n_obs;i++)
                    {
                     out_series[n_obs-1-i] = tseries[n_obs-1-i];
                    }                    
                   }

  
                          

                                 
                  mdfa.set_tseries(tseries,n_obs,n_rep);
              
                  
                  if(recompute_filter) 
                  {   
                   
//                      System.out.println(date_stamp + " " + tseries[n_obs-1] + " " + tseries[n_obs - 2] + " " + tseries[n_obs-3] + " " + tseries[n_obs-4] + " " 
//                         + price.get(price.size()-1));
                   
                    mdfa.computeFilterGeneral(true, false);        
                    
                    
                    
                    b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                    for(l=0;l<L;l++)
                    {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                    }
                    //if(print_debug) System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                    pulse_count = 0;
                    b_copy = new double[mdfa.b.length];
                    System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                    //b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                    
                  }
               
                  
                  
                  //---- apply the mother filter------------
                  sum = 0.0;
                  if(!bid_ask_data)
                  {
                   for(j=1;j<n_rep;j++)
                   {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                   }
                  }
                  else 
                  {
                   for(l=0;l<L;l++) {sum = sum + b_coeffs[l]*out_series[n_obs-1-l];}// System.out.println(out_series[n_obs-1-l]);}
                  } 
                  //System.out.println("");
                  
                  current_signal = sum; 
                                    
                  if(sig_inverse) {current_signal = -current_signal;}                  
                                    
 
                  dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + formatter3.format(lo_price.get(lo_price.size()-1)) + ", " + formatter3.format(hi_price.get(hi_price.size()-1)) + ", " + formatter.format(out_series[n_obs-1]) + ", " + formatter.format(current_signal));

                  
                  xt[trade_count] = close_series.get(close_series.size() - 1).doubleValue();
                  signal[trade_count] = current_signal;    
                  
                  if(diff_sig_trading && trade_count > 0) {signal[trade_count] = signal[trade_count] - signal[trade_count-1];} 
                  
                  prix[trade_count] = price.get(price.size()-1).doubleValue();
                  lo_prix[trade_count] = lo_price.get(lo_price.size()-1).doubleValue();
                  hi_prix[trade_count] = hi_price.get(hi_price.size()-1).doubleValue();
              
              
                  trade_count++;       
                  
                  if(trade_count%ndays == 0) {recompute_filter = true;}
                  else {recompute_filter = false;}
                      
             }
      
          } //finish with running the data
      
      
         //---- now run through trading rules and compute stats------
      
         last_trade_index = prix.length;
         insampleTradingDiff_Cust_SL(prix, lo_prix, hi_prix, signal, prix.length);
         ROI = account[end_time_index];
         if(jpy) {ROI = .01*ROI;}
                   
                    
         mean = 0; dtotal_trades = 0;          
         diff_account = new double[account.length]; diff_account[0] = 0; 
                    

                    //System.out.println(date_tokens[0]);
                    for(k=1;k<account.length;k++)
                    {
                     diff_account[k] = account[k] - account[k-1];
                     
                     ret_dist[k] = ret_dist[k] + diff_account[k];
                     
                     if(diff_account[k] > 0) 
                     {
                      pos_ret_dist[k] = pos_ret_dist[k] + 1;
                      pos_ret_mean_time[k] = pos_ret_mean_time[k] + diff_account[k];
                      pos_ret_mean = pos_ret_mean + diff_account[k];
                      n_pos_ret++;
                     }
                     else if(diff_account[k] < 0) 
                     {
                      neg_ret_dist[k] = neg_ret_dist[k] + 1;
                      neg_ret_mean_time[k] = neg_ret_mean_time[k] + diff_account[k];
                      neg_ret_mean = neg_ret_mean + diff_account[k];
                      n_neg_ret++;
                     }
                    }
                    
                    //----- now do the trade analysis ---------
                    trade_started = 0;  
                    for(k=1;k<account.length;k++)
                    {
                     
                      if(k > end_time_index) {break;}
                      
                      if(diff_account[k] != 0) {}
                      
                      if(diff_account[k] > 0)
                      {
                        
                        mean = mean + diff_account[k];
                        dtotal_trades = dtotal_trades+1;
                        
                        min_pnl_dd = 0; 
                        max_pnl_uu = 0;
                        
                        l=k;
                        while(l > 0 )
                        {   
                          if(pnl[l] < min_pnl_dd) {min_pnl_dd = pnl[l];}  
                          if(png[l] > max_pnl_uu) {max_pnl_uu = png[l];}
                          
                          l--;
                          
                          if(diff_account[l] != 0) {break;}
                        }
                        //System.out.println("trade, " + k + " " + trade_started + " " + end_time_index);
                        if(trade_started != end_time_index) 
                        {
                         pos_trades_started[trade_started] = pos_trades_started[trade_started] + 1;
                         pos_trades_started_mean[trade_started] = pos_trades_started_mean[trade_started] + diff_account[k];
                         if(jpy) {mdfaTrades.add(new MDFATrade(date_stamp, trade_started, k, .01*min_pnl_dd, .01*max_pnl_uu, .01*diff_account[k]));}
                         else {mdfaTrades.add(new MDFATrade(date_stamp, trade_started, k, min_pnl_dd, max_pnl_uu, diff_account[k]));}
                         //System.out.println("trade, " + trade_started + " " + k + " " + min_pnl_dd + " " + diff_account[k]);
                        } 
                        trade_started = k;
                      }
                      if(diff_account[k] < 0)
                      {
                        
                        mean = mean + diff_account[k]; dtotal_trades = dtotal_trades+1;
                        
                        min_pnl_dd = 0; 
                        max_pnl_uu = 0;
                        
                        l=k;
                        while(l > 0 )
                        {   
                          if(pnl[l] < min_pnl_dd) {min_pnl_dd = pnl[l];}
                          if(png[l] > max_pnl_uu) {max_pnl_uu = png[l];}
                          
                          l--;
                          
                          if(diff_account[l] != 0) {break;}
                        }
                        //System.out.println("trade, " + k + " " + trade_started + " " + end_time_index);
                        if(trade_started != end_time_index)
                        {
                         neg_trades_started[trade_started] = neg_trades_started[trade_started] + 1;
                         neg_trades_started_mean[trade_started] = neg_trades_started_mean[trade_started] - diff_account[k];
                         if(jpy) {mdfaTrades.add(new MDFATrade(date_stamp, trade_started, k, .01*min_pnl_dd, .01*max_pnl_uu, .01*diff_account[k]));}
                         else {mdfaTrades.add(new MDFATrade(date_stamp, trade_started, k, min_pnl_dd, max_pnl_uu, diff_account[k]));}
                         //System.out.println("trade, " + trade_started + " " + k + " " + min_pnl_dd + " " + diff_account[k]);
                        }

                        trade_started = k;
                      }  
                     }  
      
      
          mean = mean/(double)dtotal_trades; standard_deviation = 0;
          for(k=0;k<account.length;k++)
          {
            if(diff_account[k] != 0)
            {
              standard_deviation = standard_deviation + (diff_account[k] - mean)*(diff_account[k] - mean);
            }
          }
          standard_deviation = Math.sqrt(standard_deviation)/(double)dtotal_trades;
          
      
          System.out.println("size = " + dailyReport.size());
      
          for(k=0;k<dailyReport.size();k++)
          {
            System.out.println(dailyReport.get(k+start_seq) + ", " + formatter.format(account[k]) + ", " + formatter.format(pnl[k]) + ", " + formatter.format(png[k]));
          }                  
          longROI = account[end_time_index];                   
          shortROI = account[end_time_index];

                    
                    
          for(i = 0; i < datesAll.size(); i++)
          {
            returns.add(diff_account[i]);
            dailyoutret.add(prix[i]);
            date_returns.add(datesAll.get(i)+ ", " + diff_account[i]);
            System.out.println(datesAll.get(i)+ ", " + diff_account[i]);
          }         

         total_ROI = account[account.length-1];

         
        
        mean_perf = mean;
        double risk = -neg_ret_mean/(double)n_neg_ret;
        System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        double reward = pos_ret_mean/(double)n_pos_ret;
        System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
        kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
        ulcer_index = ulcerIndex(returns.toArray(new Double[0])); 
        
        System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
        System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        

      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      

      
        
        
        
     }   
     return computed; 
      
  }      
      
      
      
      
      
      
      
      
      
      
      
     
    
    
    
    
    
    public boolean startStrategy()
    {
      
      int j,i,l,N,file_count,k;
      double sum = 0;
      Double D; 
      double spreadnow;
      String ddelims = "[-]";
      boolean computed = false;
      double b0_sum = 0;
      int nobs_count=0;
      boolean print_filter = false;
      boolean no_return = true;
      String date_stamp,strline; 
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0;
      int advance = 8; // number of observations afterwhich we optimize
      double prev_signal;
      int changed_signs = 0;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      
      int all_trades, all_succ_trades, tot_all_trades, tot_succ_trades;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      tot_all_trades = 0; tot_succ_trades = 0;
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      
      
      for(i=0;i<trade_obs;i++) 
      {
       neg_ret_dist[i] = 0; pos_ret_dist[i] = 0; ret_dist[i] = 0; pos_ret_mean_time[i] = 0; neg_ret_mean_time[i] = 0;
       neg_trades_started[i] = 0; pos_trades_started[i] = 0; 
       neg_trades_started_mean[i] = 0; pos_trades_started_mean[i] = 0;
      }
      pos_ret_mean = 0; neg_ret_mean = 0; n_pos_ret = 0; n_neg_ret = 0;
      
      try{
        //PrintWriter overall = new PrintWriter(new FileWriter("performance_total.dat"));
        PrintWriter max_int = new PrintWriter(new FileWriter("max_interpolation.dat"));
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance.dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
        PrintWriter svmout = new PrintWriter(new FileWriter("neural.dat"));
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         
         if(dataFiles[file_count].indexOf("JPY") != -1 && dataFiles[file_count].indexOf("NOKJPY") == -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
         }
         
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; nobs_count = 0; computed = false;
         while((strline = br.readLine()) != null)
         {

          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
          if(change_time_zone) 
          {
            String[] dts = tokens[0].split("[ ]+");
            String[] ymd = dts[0].split("[-]+");
            String[] hmss = dts[1].split("[:]+");
            
            DateTime localTime = new DateTime((new Integer(ymd[0])).intValue(), (new Integer(ymd[1])).intValue(), (new Integer(ymd[2])).intValue(), 
                                            (new Integer(hmss[0])).intValue(), (new Integer(hmss[1])).intValue(), (new Integer(hmss[2])).intValue(), 0);
                                            
            localTime = localTime.plusHours(13); 
            tokens[0] = localTime.toString(fmt);         
          }
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
         
          //System.out.println(strline + " " + final_trade_time);
         
          if(let_shop) {setTimeInterval(insampStart, permEndTrading, final_trade_time);}
          else {setTimeInterval(insampStart, permEndTrading, permEndTrading);}
          friday = false;
          //System.out.println(strline + " " + forex24);
          if(forex24) //need to change friday trading end time 
          {         
           //System.out.println(weekend.dayOfWeek().getAsText());
           if(weekend.dayOfWeek().getAsText().equals("Friday") && !change_time_zone) //change number of trading days
           {
            if(min30trade) 
            {setTimeInterval(insampStart, "16:30:00", "16:30:00");}
            else 
            {setTimeInterval(insampStart, "16:00:00", "16:00:00");}
            friday = true;            
           }
           else if(weekend.dayOfWeek().getAsText().equals("Saturday") && change_time_zone) 
           {
            if(min30trade) 
            {setTimeInterval(insampStart, "05:30:00", "05:30:00");}
            else 
            {setTimeInterval(insampStart, "05:00:00", "05:00:00");}
            friday = true;     
           }
           else
           {
            if(let_shop) {setTimeInterval(insampStart, permEndTrading, final_trade_time);}
            else {setTimeInterval(insampStart, permEndTrading, permEndTrading);}
           }
           
           //System.out.println(startingTime + " " + insampStart + " " + endingTime + " " + final_trade_time + " " + start_seq + " " + trade_obs);           
          }
          //System.out.println(startingTime + " " + insampStart + " " + endingTime + " " + start_seq + " " + trade_obs); 
          
          
           if(date_stamp.indexOf(insampStart) != -1)
           {
            
            signal = new double[trade_obs];
            xt = new double[trade_obs];
            lag_signals = new double[trade_obs];
            prix = new double[trade_obs];
            lo_prix = new double[trade_obs];
            hi_prix = new double[trade_obs];
           }          
          
          
          //print_filter = true;
          latestDates.add(date_stamp);
          
          
          //use only iqfeed data for insample, use ib for out-of-sample
          
          //strategy is as follows 
          
          
          if(bid_ask_data)
          {
            
            D = new Double(tokens[4]); close_series.add(D);
           
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
            }
            else
            {
             bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
            }
            
            spreadnow = ask.get(ask.size() - 1) - bid.get(bid.size() - 1); 
            //System.out.println(date_stamp + " spread = " + spreadnow);

             if((date_stamp.indexOf(startingTime) != -1) && current_signal > 0)
             {
              D = new Double(tokens[3]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(bid.get(bid.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(bid.get(bid.size()-1)) - log(ask.get(ask.size()-2)));} //high urgency
             }
             else if((date_stamp.indexOf(startingTime) != -1) && current_signal < 0)
             {
              D = new Double(tokens[2]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(ask.get(ask.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(ask.get(ask.size()-1)) - log(bid.get(bid.size()-2)));} //high urgency
             }
             else if((date_stamp.indexOf(startingTime) == -1) && current_signal > 0)
             {
              D = new Double(tokens[2]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(bid.get(bid.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(bid.get(bid.size()-1)) - log(ask.get(ask.size()-2)));} //high urgency
             }
             else if((date_stamp.indexOf(startingTime) == -1) && current_signal < 0)
             {
              D = new Double(tokens[3]); price.add(D);
              if(spread_method == 0)  
              {live_series.add(log(bid.get(bid.size()-1)) - log(bid.get(bid.size()-2)));} //low urgency
              else if(spread_method == 1) 
              {live_series.add(log(ask.get(ask.size()-1)) - log(mid.get(mid.size()-2)));} //med urgency
              else
              {live_series.add(log(ask.get(ask.size()-1)) - log(bid.get(bid.size()-2)));} //high urgency
             }             
             else
             {
              D = new Double(tokens[1]); price.add(D); 
              
              D = new Double(tokens[4]); 
              if(ib_data && ib_data_hash.containsKey(tokens[0]))
              {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
              else
              {live_series.add(new Double(tokens[4]));}
             }
             //lo_price.add((new Double(tokens[7])) + spreadnow/2.0); hi_price.add((new Double(tokens[8])) - spreadnow/2.0);
             
             if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
             {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
             else
             {lo_price.add((new Double(tokens[7])) + spreadnow); hi_price.add((new Double(tokens[8])));}
             
             
             
             //System.out.println(date_stamp + " " + live_series.get(live_series.size() - 1));
//              if(live_series.size() > 2)
//              {
//               System.out.println("\n" + date_stamp + " " + live_series.get(live_series.size() - 1) + ", ask = " + ask.get(ask.size()-1) + ", bid = " + bid.get(bid.size()-1) + 
//               ", ask_t-1 " + ask.get(ask.size()-2) + ", bid_t-1 " + bid.get(bid.size()-2));
//              }
             //spread.println(date_stamp + " " + (new Double(tokens[3]) - new Double(tokens[2])));
                   
          }
          else
          {          
          
           D = new Double(tokens[1]); price.add(D);
           D = new Double(tokens[2]); close_series.add(D);
          
          //if(print_debug)System.out.println(date_stamp + " " + price.get(price.size()-1) + " " + close_series.get(close_series.size()-1));
          //System.out.println(date_stamp + " " + price.get(price.size()-1) + " " + close_series.get(close_series.size()-1));
           if(n_rep > 2 && tokens.length > 3)
           {
            D = new Double(tokens[3]); exp_series_1.add(D); 
           }
           if(n_rep > 3 && tokens.length > 4)
           {
            D = new Double(tokens[4]); exp_series_2.add(D);               
           }           
          }
          
          //System.out.println(nobs_count + " " + n_obs);  
          if(nobs_count >= n_obs) 
          {
                computed = true;
                 //---- collect data and put in time series                 
                if(daily_strategy) {trading_hours = true;}
                
                if(date_stamp.indexOf(insampStart) != -1 || trading_hours) //start at 9:30
                {
//                  if(print_debug)System.out.println("Start the day trading");
//                  System.out.println("Starting TIME = " + startingTime);
                 //System.out.println("ENDING TIME = " + endingTime + " " + finalCall);
                 if(trade_count < trade_obs)
                 {
                  //--- new series 
                  
                  //System.out.println("trade count = " + trade_count + " " + trade_obs);
                  n_out_samp++;
                  trading_hours = true;
                  tseries = new double[n_rep*n_obs];
                  out_series = new double[n_obs];  
                    
                  if(!bid_ask_data)
                  {
                   if(n_rep == 3)
                   {
                    for(i=0;i<n_obs;i++)
                    {
                      tseries[n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                      if(n_rep > 2 && exp_series_1.size() > 0)
                      {tseries[n_obs + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                      if(exp_series_2.size() > 0)
                      {tseries[n_obs*2 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }                 
                   }
                   else        
                   {        
                    for(i=0;i<n_obs;i++)
                    {
                     tseries[n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }
                   }
                  }
                  else
                  {
                  
                   if(classical_data)
                   {
                    for(i=0;i<n_obs;i++)
                    {
                     tseries[n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }      
                    
                    for(i=0;i<n_obs;i++)
                    {
                     out_series[n_obs-1-i] = tseries[n_obs-1-i];
                    }                    
                    
                   }
                   else
                   {
                    for(i=0;i<n_obs;i++)
                    {
                     tseries[n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = close_series.get(close_series.size() - 1 - i);
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                    }
                    
                    for(i=0;i<n_out_samp;i++)
                    {
                     out_series[n_obs-1-i] = live_series.get(live_series.size() - 1 - i);
                    }
                    for(i=n_out_samp;i<n_obs;i++)
                    {
                     out_series[n_obs-1-i] = tseries[n_obs-1-i];
                    }                    
                   }

                  
                  }
                          
                  //---- set new time series------------      
//                   if((gaussianize || volatility_thresh) || save_volatility)
//                   {
//                     if(date_stamp.indexOf(startingTime) != -1) 
//                     {      
//                      gaussianizeData(.9,.01,2.0);
//                      morningSeries = new double[n_rep*n_obs];
//                      System.arraycopy(tseries,0,morningSeries,0,morningSeries.length);
//                      computePeriodogram();
//                     }
//                     else
//                     {
//                      updateGaussianData();
//                      gauss_target.add(tseries[n_obs-1]);
//                      if(n_rep > 2)
//                      {
//                       gauss_1.add(tseries[n_obs*2 + n_obs-1]);
//                       gauss_2.add(tseries[n_obs*3 + n_obs-1]); 
//                      }
//                      for(i=0;i<trade_count;i++)
//                      {
//                        tseries[n_obs-1-i] = gauss_target.get(gauss_target.size()-1-i);
//                        tseries[n_obs + n_obs - 1 - i] = tseries[n_obs-1-i];
//                        if(n_rep > 2)
//                        {tseries[n_obs*2 + n_obs-1 - i] = gauss_1.get(gauss_1.size()-1-i);}
//                        if(n_rep > 3)
//                        {tseries[n_obs*3 + n_obs-1 - i] = gauss_2.get(gauss_2.size()-1-i);}
//                      }   
//                      
//                      //System.out.println(date_stamp + " " + tseries[n_obs-1] + " " + tseries[n_obs*2 + n_obs-1] + " " + tseries[n_obs*3 + n_obs-1]);
//                      
//                     }
//                     avg_vol = avg_vol + vol[vol.length-1];
//                   }
                
                                 
                  mdfa.set_tseries(tseries,n_obs,n_rep);
              
                  
                  //--- compute new filter --------
                  //if(pulse_count >= recomp_pulse || (date_stamp.indexOf(startingTime) != -1)) 
                  if(date_stamp.indexOf(insampStart) != -1 && day_count == 0) 
                  {   
                   
//                      System.out.println(date_stamp + " " + insampStart + " " + startingTime + " " + tseries[n_obs-1] + " " + tseries[n_obs - 2] + " " + tseries[n_obs-3] + " " + tseries[n_obs-4] + " " 
//                         + price.get(price.size()-1));
                   
                   if(!useH0) //general filtering
                   {
                   

                    //if(print_debug) System.out.println(date_stamp);
                    
//                     if(date_stamp.indexOf("2013-04-03") != -1 || date_stamp.indexOf("2013-04-02") != -1) 
//                     { 
//                       System.out.println(date_stamp);
//                       System.out.println("AT THE END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//                     }
                    
                   
                    mdfa.computeFilterGeneral(true, print_filter);        
                    
                    
                    
                    b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                    for(l=0;l<L;l++)
                    {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                    }
                    //if(print_debug) System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                    pulse_count = 0;
                    b_copy = new double[mdfa.b.length];
                    System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                    b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                    
                    b0_sum = 0;
                    for(l=0;l<L;l++)
                    {b0_sum = b0_sum + b_coeffs[l];}
                    b0_trend.add(b0_sum);
//                     if(filters.size() > 5)
//                     {
//                      System.out.println("Old BO = " + b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[L+L]);     
//                      averageFilters(3,.3);
//                      System.out.println("Avg BO = " + b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[L+L]);
//                     }
                   }
                   else
                   {
                    if(H0set) //H0_Optimization_Lookback2(double delta, double interp_start)
                    {
                     if(!lookback_ready) {H0_Optimization(h0_outsamp, localdecay, h0_delta, interp_start); avg_volatility.add(0.0); System.out.println("First");}
                     //else{H0_Optimization_Lookback(h0_delta, 0.0); max_int.println(max_interp_value + " " + max_ret_interp + " " + ROI + " " + date_stamp);}
                     else
                     {
                       //if(morning_optimize)
                       //{H0_Optimization(h0_outsamp, localdecay, h0_delta, interp_start);}
                       //else
                       //{
                        H0_Optimization_Lookback(h0_delta, 0);  
                        //H0_Optimization_Lookback_Test(yesterday_series, h0_delta, 0.0, optval, cut); 
                        //max_int.println(max_interp_value + " " + optval + " " + max_ret_interp + " " + zero_ret);     
                        max_int.println(max_interp_value + " " + max_ret_interp);     
                        //System.out.println("New BO = " + b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[L+L]);
                       //}
                     }
                     //H0_Optimization_Morning(localdecay, h0_delta, interp_start); 
                     //max_int.println(max_interp_value + " " + max_ret_interp);
                     pulse_count = 0;
                    }
                    else
                    {System.out.println("Problem: Couldn't find H0 filter"); break;}                  
                   }           
                  }
                  //System.out.println("Out of sample " + date_stamp);
                  //else {pulse_count++;}
                              
                  //---- change filter to updated filter -------------------------------- 
                  if((morning_optimize && (date_stamp.indexOf(optimizeTime) != -1)) && lookback_ready) 
                  { 
                     System.out.println("Optimizing at 12:30");
                     advance = 13;   //   !!!! Caution, this value must correspond with the optimizeTime (e.g. 8 = 11:30)  
                     //optval = max_interp_value;        
                     H0_Optimization_Intraday(h0_delta, 0.0, advance); //find optimized parameter at new time
                    // max_int.println(max_interp_value + " " + optval);         
                  } 
                   
                  if(date_stamp.indexOf("09:30") != -1)
                  {reg_trading_hours = true;} 
                   
                  //start_seq = 10; 
                  if(pulse_count < recomp_pulse && reg_trading_hours)
                  {       
                   if(H0set && (pulse_count%start_seq==0) && pulse_count > 0) //5 is best
                   //if(H0set)
                   {H0_Optimization_Update(0, trade_count, true, true);}
                   pulse_count++;
                  }
                  else if(pulse_count < recomp_pulse && date_stamp.indexOf("08:30") != -1)
                  {
                    System.out.println("Starting out for the day");
                    double[] outs = computeInterpolation(tseries, 0);
                       
                    int lag4 = 2*(n_obs-L+1) + 2*(n_rep+1)*K1;
                    b_coeffs = new double[(n_rep-1)*L];
                    for(l=0;l<L;l++)
                    {
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = outs[lag4 + L*(i+1)+l];}
                    }          
                    pulse_count = 0;                  
                  }
                              
                              
                  //---- apply the mother filter------------
                  sum = 0.0;
                  if(!bid_ask_data)
                  {
                   for(j=1;j<n_rep;j++)
                   {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                   }
                  }
                  else 
                  {
                   for(l=0;l<L;l++) {sum = sum + b_coeffs[l]*out_series[n_obs-1-l];}// System.out.println(out_series[n_obs-1-l]);}
                  } 
                  //System.out.println("");
                  
                  prev_signal = current_signal;
                  current_signal = sum; 
                                    
                  if(sig_inverse) {current_signal = -current_signal;}                  
                                    
                  if(morning_optimize && (date_stamp.indexOf(optimizeTime) != -1)) 
                  {
                    if((prev_signal > 0 && current_signal < 0) || (prev_signal < 0 && current_signal > 0))
                    {
                      System.out.println("Sign of signal changed after optimization, forced to enter transaction.");
                      changed_signs++; 
                    }
                  }
                  
                                  
                  
                  if(cont_lookback)  //compute previous lag 
                  {                   
                   if(compute_lag)
                   {
                     mdfa.set_lag(lag+1);
                     mdfa.computeFilterGeneral(true,false);  
                     b_lag = new double[(n_rep-1)*L];
                     for(l=0;l<L;l++)
                     {
                      for(i=0;i<n_rep-1;i++)
                      {b_lag[L*i + l] = mdfa.b[L*(i+1)+l];}
                     }                       
                     sum = 0.0;
                     for(j=1;j<n_rep;j++)
                     {for(l=0;l<L;l++) {sum = sum + b_lag[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}}
                     lag_signal = sum;
                   }
                   else
                   {
                    sum = 0.0;
                    for(j=1;j<n_rep;j++)
                    {for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-2-l];}}                                                        
                    lag_signal = sum;
                   }
                  }
                  
                  
                  
                                   
                  if(adaptive_filtering)
                  {
                   //------------------- Adaptive steps --------------------------------------------
                   if(adapt_pulse_count >= adaptive_recomp_pulse || (date_stamp.indexOf(startingTime) != -1))
                   {
                     adaptiveUpdateFilter();
                     updateAdaptiveSignal();
                     adapt_pulse_count = 0;
                   }
                   else
                   {updateAdaptiveSignal(); adapt_pulse_count++;}
                   
                   current_signal = update_signal[update_signal.length-1];
                   lag_signal = update_signal[update_signal.length-2];
                  }
                  
//                   if(print_debug)
//                   {
//                    System.out.println(date_stamp + ", " + formatter.format(price.get(price.size()-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + formatter.format(current_signal));      
//                   }
//                   if(!forex24) {dailyReport.add(""+date_stamp + ", " + formatter2.format(Math.exp(price.get(price.size()-1))) + ", " + formatter.format(tseries[n_obs-1]) + ", " + formatter.format(current_signal));}
//                   else {dailyReport.add(""+date_stamp + ", " + formatter.format(price.get(price.size()-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + formatter.format(current_signal));}
                  if(!bid_ask_data)
                  {dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + formatter3.format(lo_price.get(lo_price.size()-1)) + ", " + formatter3.format(hi_price.get(hi_price.size()-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + formatter.format(current_signal));}
                  else
                  {dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + formatter3.format(lo_price.get(lo_price.size()-1)) + ", " + formatter3.format(hi_price.get(hi_price.size()-1)) + ", " + formatter.format(out_series[n_obs-1]) + ", " + formatter.format(current_signal));}
                  //                   double binsig = 0;
//                   if(current_signal > 0) binsig = 1; 
//                   else binsig = -1;
//                   
//                   dailyReport.add(""+date_stamp + ", " + formatter.format(price.get(price.size()-1)) + ", " + binsig);
                  
                  //System.out.println(dailyReport.size());
                  
                  xt[trade_count] = close_series.get(close_series.size() - 1).doubleValue();
                  signal[trade_count] = current_signal;    
                  
                  if(diff_sig_trading && trade_count > 0) {signal[trade_count] = signal[trade_count] - signal[trade_count-1];} 
                  
                  prix[trade_count] = price.get(price.size()-1).doubleValue();
                  lo_prix[trade_count] = lo_price.get(lo_price.size()-1).doubleValue();
                  hi_prix[trade_count] = hi_price.get(hi_price.size()-1).doubleValue();
              
                  if(trade_count > 0 && cont_lookback)
                  {lag_signals[trade_count-1] = lag_signal;}
              
                  trade_count++;             
                 }
                 
                 
                 if(date_stamp.indexOf(finalCall) != -1) 
                 //if(date_stamp.indexOf(endingTime) != -1) 
                 {
                    //System.out.println("ENDING TIME Data at " + endingTime);
//                    for(i=0;i<L;i++)
//                    {
//                      System.out.println(tseries[n_obs-1-i] + " " + tseries[n_obs + n_obs-1 - i] + " " + tseries[2*n_obs + n_obs-1 - i] + " " + tseries[3*n_obs + n_obs-1 - i]); 
//                    }
                   
                  if(reset_signal) {current_signal = 0;}
                  current_signal = 0;
                  n_out_samp = 0;
                  reg_trading_hours = false;     
                  daily_return = price.get(price.size()-1).doubleValue() - log_price;
                  log_price = price.get(price.size()-1).doubleValue();
                  trade_count = 0; trading_hours = false;
                  pulse_count = 0;
                  
//                   if(startingTime.equals("09:30"))
//                   {
                   
                    actual_price = new double[trade_obs - start_seq];
                    actual_loprice = new double[trade_obs - start_seq];
                    actual_hiprice = new double[trade_obs - start_seq];
                    actual_signal = new double[trade_obs - start_seq];
                    for(i=0;i<actual_price.length;i++) 
                    {
                     actual_price[i] = prix[i+start_seq]; 
                     actual_signal[i] = signal[i+start_seq];
                     actual_loprice[i] = lo_prix[i+start_seq];
                     actual_hiprice[i] = hi_prix[i+start_seq];
                     
                    }/// System.out.println(actual_price[i] + " " + actual_signal[i]);}                                      
                    
                    short_sell = true; long_buy = true;
                    if(!let_shop)
                    {
                     insampleTradingDiff(actual_price, actual_signal, actual_price.length);
                     ROI = account[account.length-1]; 
                    }
                    else
                    {
                     insampleTradingDiff_Cust_SL(actual_price, actual_loprice, actual_hiprice, actual_signal, actual_price.length);
                     ROI = account[end_time_index];
                    }
                    
                    if(jpy) {ROI = .01*ROI;}
                    //date_tokens = date_stamp.split("[ ]+");
                    
                    
                    diff_account = new double[account.length]; diff_account[0] = 0; 
                    

                    //System.out.println(date_tokens[0]);
                    for(k=1;k<account.length;k++)
                    {
                     diff_account[k] = account[k] - account[k-1];
                     
                     ret_dist[k] = ret_dist[k] + diff_account[k];
                     
                     if(diff_account[k] > 0) 
                     {
                      pos_ret_dist[k] = pos_ret_dist[k] + 1;
                      pos_ret_mean_time[k] = pos_ret_mean_time[k] + diff_account[k];
                      pos_ret_mean = pos_ret_mean + diff_account[k];
                      n_pos_ret++;
                     }
                     else if(diff_account[k] < 0) 
                     {
                      neg_ret_dist[k] = neg_ret_dist[k] + 1;
                      neg_ret_mean_time[k] = neg_ret_mean_time[k] + diff_account[k];
                      neg_ret_mean = neg_ret_mean + diff_account[k];
                      n_neg_ret++;
                     }
                    }
                    
                    //----- now do the trade analysis ---------
                    trade_started = 0;  no_return = true;
                    for(k=1;k<account.length;k++)
                    {
                     
                      if(k > end_time_index) {break;}
                      
                      if(diff_account[k] != 0) {no_return = false;}
                      
                      if(diff_account[k] > 0)
                      {
                        
                        min_pnl_dd = 0; 
                        max_pnl_uu = 0;
                        
                        l=k;
                        while(l > 0 )
                        {   
                          if(pnl[l] < min_pnl_dd) {min_pnl_dd = pnl[l];}  
                          if(png[l] > max_pnl_uu) {max_pnl_uu = png[l];}
                          
                          l--;
                          
                          if(diff_account[l] != 0) {break;}
                        }
                        //System.out.println("trade, " + k + " " + trade_started + " " + end_time_index);
                        if(trade_started != end_time_index) 
                        {
                         pos_trades_started[trade_started] = pos_trades_started[trade_started] + 1;
                         pos_trades_started_mean[trade_started] = pos_trades_started_mean[trade_started] + diff_account[k];
                         if(jpy) {mdfaTrades.add(new MDFATrade(date_tokens[0], trade_started, k, .01*min_pnl_dd, .01*max_pnl_uu, .01*diff_account[k]));}
                         else {mdfaTrades.add(new MDFATrade(date_tokens[0], trade_started, k, min_pnl_dd, max_pnl_uu, diff_account[k]));}
                         //System.out.println("trade, " + trade_started + " " + k + " " + min_pnl_dd + " " + diff_account[k]);
                        } 
                        trade_started = k;
                      }
                      if(diff_account[k] < 0)
                      {
                        
                        min_pnl_dd = 0; 
                        max_pnl_uu = 0;
                        
                        l=k;
                        while(l > 0 )
                        {   
                          if(pnl[l] < min_pnl_dd) {min_pnl_dd = pnl[l];}
                          if(png[l] > max_pnl_uu) {max_pnl_uu = png[l];}
                          
                          l--;
                          
                          if(diff_account[l] != 0) {break;}
                        }
                        //System.out.println("trade, " + k + " " + trade_started + " " + end_time_index);
                        if(trade_started != end_time_index)
                        {
                         neg_trades_started[trade_started] = neg_trades_started[trade_started] + 1;
                         neg_trades_started_mean[trade_started] = neg_trades_started_mean[trade_started] - diff_account[k];
                         if(jpy) {mdfaTrades.add(new MDFATrade(date_tokens[0], trade_started, k, .01*min_pnl_dd, .01*max_pnl_uu, .01*diff_account[k]));}
                         else {mdfaTrades.add(new MDFATrade(date_tokens[0], trade_started, k, min_pnl_dd, max_pnl_uu, diff_account[k]));}
                         //System.out.println("trade, " + trade_started + " " + k + " " + min_pnl_dd + " " + diff_account[k]);
                        }

                        trade_started = k;
                      }

                      
                    
                    }
                    
                    //the day had no action, so place an empty trade
                    if(no_return) 
                    {mdfaTrades.add(new MDFATrade(date_tokens[0], 0, end_time_index, 0, 0, .00001));}
                    
                    
  
 
 
                    if(print_debug) 
                    {
                     System.out.println("");
                     for(k=0;k<account.length;k++)
                     {
                      System.out.println(dailyReport.get(k+start_seq) + ", " + formatter.format(account[k]) + ", " + formatter.format(pnl[k]) + ", " + formatter.format(png[k]));
                     }                  
                    }
                    all_trades = total_trades; all_succ_trades = succ_trades; 
                    tot_all_trades = tot_all_trades + total_trades; 
                    tot_succ_trades = tot_succ_trades + succ_trades; 
                    
                    short_sell = false; long_buy = true;
                    //insampleTradingDiff_Cust_SL(actual_price, actual_signal, actual_price.length); 
                    //longROI = account[account.length-1];
                    longROI = account[end_time_index];
                    
                    short_sell = true; long_buy = false;
                    //insampleTradingDiff_Cust_SL(actual_price, actual_signal, actual_price.length); 
                    //shortROI = account[account.length-1];                    
                    shortROI = account[end_time_index];

                    

                    
                  final_trades.add(final_trade); last_trades.add(last_trade);
                  sub_returns.add(diff_band);
                  //I_MDFA.plotData(xt, signal, trade_obs);
                  
               
                  
                  
                  avg_vol = avg_vol/trade_obs;
                  

                  dailyReport.clear(); 
                  mnmx = minmax(account);             
                  ratio = ((double)succ_trades)/((double)total_trades);              
              
                  if(ROI < 0) {losses_in_arow++;}
                  else {losses_in_arow = 0;}
                  
                  if(useH0)
                  {
                   //System.out.println("using H0");
                   if(lookback_ready)
                   {
                   dailyoutret.add(daily_return);
                   maxIntValue.add(max_interp_value); 
                   avg_volatility.add(avg_vol-0.25);
                   returns.add(ROI - all_trades*tradingCost);
                   longreturns.add(longROI);
                   shortreturns.add(shortROI);
                   dropdowns.add(drop_down);
                   success.add(ratio);
                   trade_days.add(date_tokens[0]);
                   dates_low_high.add(date_stamp + ", " + mnmx[0] + ", " + mnmx[1]);     
                   }
                   else
                   {lookback_ready = true;}
                  }
                  else
                  {
                   //System.out.println("Not using H0");
                   dailyoutret.add(daily_return);
                   maxIntValue.add(max_interp_value); 
                   avg_volatility.add(100000.0*(avg_vol-0.30));
                   returns.add(ROI - all_trades*tradingCost); 
                   longreturns.add(longROI);
                   shortreturns.add(shortROI);
                   dropdowns.add(drop_down);
                   success.add(ratio);
                   trade_days.add(date_tokens[0]);
                   dates_low_high.add(date_stamp + ", " + mnmx[0] + ", " + mnmx[1]);                   
                  }
                 
                  if(ROI < 0 && first_trade_loss) {bad_starts++;}
                  
                  if(mnmx[1] > .005) {profit_baby++;}
                  //out.println("Result for " + date_stamp + ", ROI = " + ROI + ", " + mnmx[0] + ", " + mnmx[1] + ", " + succ_trades + ", " + total_trades);
                  if(print_debug)System.out.println("Result for " + date_stamp + ", ROI = " + (ROI - all_trades*tradingCost) + ", succ_trades = " + all_succ_trades + ", total_trades = " + all_trades);
                  
                  
//                   if(returns.size() > rolling_length) //take latest 30 days and compute a few stats
//                   {
//                     double[] latRets = new double[rolling_length];
//                     Double[] latRetsD = new Double[rolling_length];
//                     for(k=0;k<rolling_length;k++) 
//                     {
//                      latRets[rolling_length - 1 - k] = returns.get(returns.size() - 1 - k);
//                      latRetsD[rolling_length - 1 - k] = returns.get(returns.size() - 1 - k);
//                     } 
//                     //compute sharpe
//                     double[] mstd = mean_std(latRets); 
//                     double sh = Math.sqrt(250)*mstd[0]/mstd[1];
//                     //compute rank
//                     double rc = rankCoefficient(latRets,rolling_length);
//                     //compute Ulcer
//                     double ui = ulcerIndex(latRetsD);
//                     //computer Kelly
//                     
//                     rolling_ind.add(sh + " " + rc + " " + ui + " " + kp);     
//                   }
//                   else
//                   {rolling_ind.add("0 0 0 0");}
                  
                  
                  
                  
                  dailyout.println(daily_return);
                  signal = new double[trade_obs];
                  xt = new double[trade_obs];
                  prix = new double[trade_obs];  
                  lag_signals = new double[trade_obs];
                  
                  avg_vol = 0.0;
                  total_ROI = total_ROI + ROI - all_trades*tradingCost;
                  total_succ = total_succ + succ_trades;
                  total = total + total_trades;
                  //System.out.println("total trades = " + total_succ + "/" + total);
                  
                  if(ROI > 0) {trade_succ_ratio++;}
                  
                  
                  
                   yesterday_vol = 0;//avg_volatility.get(avg_volatility.size() - 2);
                   if(ROI > .006)// && (avg_volatility.get(avg_volatility.size() - 2) < 4.007053e-07 ))
                   { 
                                      
                    avg_count++;
                    for(k=1;k<n_rep;k++)
                    {
                      for(l=0;l<L;l++)
                      {b_avg[L*k + l] = b_avg[L*k + l] + b_coeffs[L*(k-1) + l];}
                    }                    
                  }
                  
                  if(day_count == recompute_day) {day_count = 0;}
                  else {day_count++;}
                  
                  
                }                 
               }                       
            }
            nobs_count++;            
            //System.out.println(nobs_count);
        }
        
        //out.println("Total ROI: " + total_ROI + ", Total successful trades = " + succ_trades + ", total trades = " + total_trades + ", rate = " + (double)succ_trades/(double)total_trades);
        System.out.println("Total ROI: " + total_ROI + ", Total successful trades = " + tot_succ_trades + ", total trades = " + tot_all_trades + ", avg_n_trades = " + (double)tot_all_trades/returns.size() + ", rate = " + (double)tot_succ_trades/(double)tot_all_trades);
        //System.out.println("Avg Friday ROI = " + (fridayROI/fridays) + ", friday success = " + ((double)fridayROI_pos/fridays));
       
        avg_n_trades = (double)tot_all_trades/(double)returns.size();
        avg_rate = (double)tot_succ_trades/(double)tot_all_trades;
        
        dailyoutret.set(0,0.0);
        if(!H0set || recomp_pulse > 0)
        {
        System.out.println("Summing statistics");
        num_trade_days = returns.size(); 
        if(num_trade_days == 0) 
        {for(i=0;i<100;i++) returns.add(0.0);}
        num_trade_days = returns.size();
        mean = 0.0;
        num_gains = 0; num_losses = 0;
        dreturns = new double[num_trade_days];
        //double ret = returns.get(0);
        double ret = 0;
        if(num_trade_days>0) dreturns[0] = 0;
        for(i=0;i<num_trade_days;i++)
        {
         ret = returns.get(i);
         if(ret > 0) {num_gains++;}
         else {num_losses++;}
         
         mean = mean + ret;
         out.println(ret);
         if(i > 0) dreturns[i] = dreturns[i-1] + ret;
        }
        mean = mean/((double)returns.size());
        
        
        double risk = -neg_ret_mean/(double)n_neg_ret;
        System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        double reward = pos_ret_mean/(double)n_pos_ret;
        System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
        kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
        ulcer_index = ulcerIndex(returns.toArray(new Double[0])); 
        
        System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
        System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
        sum_sd = 0;
        for(i=0;i<returns.size();i++)
        {sum_sd = sum_sd + (returns.get(i) - mean)*(returns.get(i) - mean)/((double)returns.size());}
        standard_deviation = Math.sqrt(sum_sd);
        }
        else
        {
         num_trade_days = lookback_returns.size(); 
        mean = 0.0;
        num_gains = 0; num_losses = 0;
        dreturns = new double[num_trade_days];
        double ret = lookback_returns.get(0);
        dreturns[0] = ret;
        for(i=0;i<num_trade_days;i++)
        {
         ret = lookback_returns.get(i);
         if(ret > 0) {num_gains++;}
         else {num_losses++;}
         
         mean = mean + ret;
         out.println(ret);
         if(i > 0) dreturns[i] = dreturns[i-1] + ret;
        }
        mean = mean/((double)lookback_returns.size());
        
        
        sum_sd = 0;
        for(i=0;i<lookback_returns.size();i++)
        {sum_sd = sum_sd + (lookback_returns.get(i) - mean)*(lookback_returns.get(i) - mean)/((double)lookback_returns.size());}
        standard_deviation = Math.sqrt(sum_sd);       
        }
        
        
        out.println("");
/*        for(i=0;i<num_trade_days;i++)
        {
           out.println(dreturns[i] + " " + dailyoutret.get(i));
           System.out.println(dreturns[i] + " " + avg_volatility.get(i));
        }*/       
        maxdraw = computeDrawdown(dreturns);
        
        //rankCoefficient(dreturns, num_trade_days);
       
        rank_coeff = segmentRankCorrelation(15, dreturns);
       
        out.println("Rank Coefficient = " + rank_coeff + ", mean = " + mean);
        System.out.println("Rank Coefficient = " + rank_coeff + ", mean = " + mean);
        mean_perf = mean; rank_perf = rank_coeff;
       }
   
   
       for(i=0;i<svm.size();i++)
       {svmout.println(svm.get(i));}          
   
   
        out.close();  
        dailyout.close();
        perform.close();
        b0_coeff.close();
        max_int.close();
        svmout.close();
        spread.close();
      }     
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}   

      total_obs = returns.size();
      if(trade_days.size() < total_obs) {total_obs = trade_days.size();}

      for(i=0;i<total_obs;i++)
      {
        //System.out.println(trade_days.get(i) + " " + returns.get(i) + " " + avg_volatility.get(i) + " " + maxIntValue.get(i) + " " + interp_vals.get(i) + crits.get(i));      
        if(!morning_optimize)
        {System.out.println(trade_days.get(i) + ", " + returns.get(i)); date_returns.add(trade_days.get(i) + ", " + dailyoutret.get(i));}
        else
        {System.out.println(trade_days.get(i) + ", " + returns.get(i) + ", " + interp_vals.get(i) + ", " + max_ranks.get(i));}  
      }
     


      for(i=0;i<morning_returns.size();i++)
      {
       double[] array = morning_returns.get(i);
       for(j=0;j<array.length-1;j++)
       {System.out.print(df4.format(array[j]) + " ");}
       System.out.println(df4.format(array[array.length-1]));
      }


      if(!morning_optimize)
      {
       for(i=0;i<crits.size();i++)
       {System.out.println(crits.get(i));}
      
       
      }

  
      System.out.println(avg_count);
      for(k=1;k<n_rep;k++)
      {
         for(l=0;l<L;l++)
         {b_avg[L*k + l] = b_avg[L*k + l]/(double)avg_count;}
      }  

      try{
        //PrintWriter overall = new PrintWriter(new FileWriter("performance_total.dat"));
       PrintWriter b0_coeff = new PrintWriter(new FileWriter("h0b0_filter.dat"));      
       b0_coeff.println(L + " " + n_rep);
       for(l=0;l<L;l++)
       {
        for(k=0;k<n_rep-1;k++)
        {b0_coeff.print(b_avg[L*k + l] + " ");}
        b0_coeff.println(b_avg[L*(n_rep-1) + l]);
       }
       b0_coeff.close();
     }
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
      
     try{
     
       PrintWriter ret_dists = new PrintWriter(new FileWriter("return_dists_"+dataFiles[0]));   
       for(l=1;l<trade_obs;l++) 
       {
        ret_dists.println((ret_dist[l]/(double)total_obs) + " " + (pos_ret_dist[l]/(double)(pos_ret_dist[l] + neg_ret_dist[l])) + " " +
        (pos_ret_mean_time[l]/(double)pos_ret_dist[l]) + " " + (neg_ret_mean_time[l]/(double)neg_ret_dist[l]));
      
       if(pos_trades_started[l-1] == 0 && neg_trades_started[l-1] == 0)
       {
         System.out.println("0    0     0");
       }
       else
       {
        System.out.println(formatter.format((pos_trades_started[l-1]/(double)(pos_trades_started[l-1] + neg_trades_started[l-1]))) + " " + formatter.format(pos_trades_started_mean[l-1]/(double)pos_trades_started[l-1])      
          + " " + formatter.format(neg_trades_started_mean[l-1]/(double)neg_trades_started[l-1]));       
       }   
       }
       ret_dists.close();   
     }
     catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}  
      

     
      for(i=0;i<full_returns_array.size();i++)
      {
       double[] array = full_returns_array.get(i);
       for(j=0;j<array.length-1;j++)
       {System.out.print(df4.format(array[j]) + " ");}
       System.out.println(df4.format(array[array.length-1]));
      }
     
     

     
     
     if(morning_optimize)
     {System.out.println("Signal changed signs at " + optimizeTime + " " + changed_signs + " times");}
      
     return computed; 
      
  }

  
 
  
    public boolean startStrategyDailyIntraday()
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      
      double pnl; 
      String[] dates; 
      int cur_month = 0;
      double convers = 1.0;
      ArrayList<String> allsignal = new ArrayList<String>();
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      perf_returns = new ArrayList<Double>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      daily_returns = new ArrayList<Double>();
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      daily_price = new ArrayList<Double>();
      
      daily_dates = new ArrayList<String>();
      ArrayList<String> daily_data = new ArrayList<String>();
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      int line; 
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      

      
      try{
        
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
        PrintWriter svmout = new PrintWriter(new FileWriter("neural_"+n_files+".dat"));
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         convers = 1;  
         
         if(dataFiles[file_count].indexOf("JPY") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(dataFiles[file_count].indexOf("NOK") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*8;
          take_profit_thresh = take_profit_thresh*8;
          global_stop_loss = global_stop_loss*8;
          stop_loss = stop_loss_thresh;
         }   
         else if(dataFiles[file_count].indexOf("XAU") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*1000;
          take_profit_thresh = take_profit_thresh*1000;
          global_stop_loss = global_stop_loss*1000;
          stop_loss = stop_loss_thresh;
         }                 
         
//          if(dataFiles[file_count].indexOf("JPY") != -1 && (dataFiles[file_count].indexOf("SEKJPY") == -1 && dataFiles[file_count].indexOf("NOKJPY") == -1))
//          {
//           System.out.println("Changed time zone to Tokyo");
//           jpy = true; 
//           //stop_loss_thresh = stop_loss_thresh*100;
//           take_profit_thresh = take_profit_thresh*100;
//           //global_stop_loss = global_stop_loss*100;
//           stop_loss = stop_loss_thresh;
//           convers = .01;
//          }
//          else if(dataFiles[file_count].indexOf("SEKJPY") != -1 || dataFiles[file_count].indexOf("NOKJPY") != -1)
//          {
//            scandi_jpy = true;
//            convers = .01;
//          }
//          else if(dataFiles[file_count].indexOf("NOK") != -1 || dataFiles[file_count].indexOf("SEK") != -1)
//          {
//            scandi = true;
//            convers = 1/6.0; 
//          }
         
         
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         while((strline = br.readLine()) != null)
         {
           daily_data.add(strline);
         }  
         
         for(line=0;line<daily_data.size()-1;line++)
         {
         
          strline = daily_data.get(line);
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           

           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
         
          String[] hours = date_tokens[1].split("[:]+");
         
         
         //--- TIME FILTER----------------------
         if((!time_filter && !hours[0].equals("17")) || ((date_tokens[1].indexOf(":00:00") != -1 || date_tokens[1].indexOf(":30:00")  != -1) && !hours[0].equals("17")))
         //if(!time_filter || ((date_tokens[1].indexOf(":00:00") != -1 || date_tokens[1].indexOf(":30:00")  != -1) && !hours[0].equals("17")))
         {
          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(date_stamp.indexOf(insampStart) != -1)// && weekend.dayOfWeek().getAsText().equals("Wednesday"))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(new Double(tokens[1]));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
          else
          {
           bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
            

          D = new Double(tokens[1]); 
          price.add(D); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
          if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
          {lo_price.add(new Double(tokens[2])); hi_price.add(new Double(tokens[2]));}
          else
          {lo_price.add((new Double(tokens[2]))); hi_price.add((new Double(tokens[2])));}
             
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && (date_stamp.indexOf(insampStart) != -1))// && weekend.dayOfWeek().getAsText().equals("Tuesday"))) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if((new Integer(intdates[1])).intValue() != cur_month || day_count == 0)
                {   
                   cur_month = (new Integer(intdates[1])).intValue();
                   //System.out.println("Recomputing filter..." + intdates[0] + "-" + intdates[1] + "-" + intdates[2]);
                   if(printall) System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, false);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall) System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                   
                   //System.out.println("\n");
                   //for(l=0;l<L;l++) {System.out.println(b_coeffs[l]);}
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                
                
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                svmout.println(date_stamp + ", " + current_signal);
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal);
                if(binary_rule)
                {
                
                
                 if(signal_profit)
                 {
                 
                   if((current_signal > 0 && in_transaction == 1) && (daily_price.get(daily_size-1) > last_price))
                   {
                   
                   
                     price_sold = daily_price.get(daily_size-1);
	             profit = price_sold - price_bought;
	             if(profit > 0) {succ_trades=succ_trades+1;}
	             total_trades=total_trades+1; 
	 
	             amount = getAmount(account.get(account.size()-1));
	             account.add(date_stamp + " " + (amount + profit));   
	             amount = getAmount(account.get(account.size()-1));
	            
	             made_trade = true;
	             dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                     log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	             in_transaction = 0;
	            
	             if(printall)
	             {
                     if(profit>0) System.out.println("Sold for a profit of " + profit);
                     else if(profit<0) System.out.println("Sold for a loss of " + profit);
                     System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                     }
                   
                     price_bought = daily_price.get(daily_size-1);
	             in_transaction = 1; 
	             last_price = daily_price.get(daily_size-1); 
	             if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	              
                   }
                   else if((current_signal < 0 && out_transaction == 1) && (daily_price.get(daily_size-1) < last_price))
                   {
                   
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    if(printall)
	             {
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);                   
                    }
                   
                    price_borrowed = daily_price.get(daily_size-1);
                    out_transaction = 1;	 
                  
                    if(printall)
	             {System.out.println("Entered short transaction at " + price_borrowed);}
                   
                   }
                 
                 }
                
                
                
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  last_price = daily_price.get(daily_size-1); 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    if(printall)
	             {
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);}
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           
	           if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall)
	             {
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                  
                   if(printall)
	           {
                   System.out.println("Entered short transaction at " + price_borrowed);
                   }
	          }
                 }
                }
                
                if(signal_strength_rule && (in_transaction == 0 && out_transaction == 0))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  last_price = daily_price.get(daily_size-1); 
                  
                  
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  
                    if(printall)
	             { 
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           
	           if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	          } 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
	           in_transaction = 0;
	           
	           if(printall){
                   System.out.println("Bought for a profit of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);}	           
	           
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                  
                   if(printall){System.out.println("Entered short transaction at " + price_borrowed);}
	          }
                 }              
                }
              
                if(!made_trade)
                {
                 account.add(date_stamp + " " + amount);
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {day_count=0;}
                
                allsignal.add(date_stamp + " " + current_signal);
                 
           }
           else if(trading_hours)
           {
                
            

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && price.get(price.size()-1) > last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
              
               cur_pnl = price.get(price.size()-1) - last_price;    
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall){System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);}
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 price_sold = price.get(price.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         //account.add(date_stamp + " " + (amount - stop_loss));   
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 
                 
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall){System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);}
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall){System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);}
                 
                 price_sold = price.get(price.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall){System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);}
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && price.get(price.size()-1) < last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
               cur_pnl =  last_price - price.get(price.size()-1);               
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall){System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);}
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         account.add(date_stamp + " " + (amount + profit));  
	         amount = getAmount(account.get(account.size()-1));

                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall){System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);}
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall){System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);}
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall){System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);}


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               //System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall){ System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
               + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl)); }
               
             
             allsignal.add(date_stamp + " " + current_signal);
               
           }      
                    
          }             
        }              
                    
      }
 
      
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i)) - getAmount(account.get(i-1));
        
        dreturns[i] = convers*dreturns[i]; //approximate pips in dollars
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0]) && !isSunday(dates[0])) //first date entry
        {
         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
         }
         
         pnl = dreturns[i];
         perf_dates.add(dates[0]);
    
        }
        else //already contains the date, so add on pnl
        {pnl = pnl + dreturns[i];}
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       //System.out.println(perf_dates.get(i) + " " + perf_returns.get(i));
       date_returns.add(perf_dates.get(i) + " " + perf_returns.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
      }

      perform.close();
      // load in intraday information Date log-diff price signal - EXCLUDES LATEST OBSERVATION!!! 
      dates_price = new String[n_obs];
      for(i=0;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2]);
//        }
      }
      
      
      n_files++;
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff);
 
      out.close(); dailyout.close(); svmout.close(); b0_coeff.close();
      
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
 
      return computed; 
      
  }   
  
  
 
 
    public boolean startStrategyDailyIntraday_MR()
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      
      boolean waiting_meanrev_down = false;
      boolean waiting_meanrev_up = false;
      double mean_rev_amnt = 0;
      
      double pnl; 
      String[] dates; 
      int cur_month = 0;
      double convers = 1.0;
      ArrayList<String> allsignal = new ArrayList<String>();
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      perf_returns = new ArrayList<Double>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      daily_returns = new ArrayList<Double>();
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = false;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      daily_price = new ArrayList<Double>();
      
      daily_dates = new ArrayList<String>();
      ArrayList<String> daily_data = new ArrayList<String>();
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      int line; 
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      

      
      try{
        
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
        PrintWriter svmout = new PrintWriter(new FileWriter("neural_"+n_files+".dat"));
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         convers = 1;  
         
         if(dataFiles[file_count].indexOf("JPY") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(dataFiles[file_count].indexOf("NOK") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*8;
          take_profit_thresh = take_profit_thresh*8;
          global_stop_loss = global_stop_loss*8;
          stop_loss = stop_loss_thresh;
         }   
         else if(dataFiles[file_count].indexOf("XAU") != -1)
         {
          jpy = true;
          stop_loss_thresh = stop_loss_thresh*1000;
          take_profit_thresh = take_profit_thresh*1000;
          global_stop_loss = global_stop_loss*1000;
          stop_loss = stop_loss_thresh;
         }                 
         
//          if(dataFiles[file_count].indexOf("JPY") != -1 && (dataFiles[file_count].indexOf("SEKJPY") == -1 && dataFiles[file_count].indexOf("NOKJPY") == -1))
//          {
//           System.out.println("Changed time zone to Tokyo");
//           jpy = true; 
//           //stop_loss_thresh = stop_loss_thresh*100;
//           take_profit_thresh = take_profit_thresh*100;
//           //global_stop_loss = global_stop_loss*100;
//           stop_loss = stop_loss_thresh;
//           convers = .01;
//          }
//          else if(dataFiles[file_count].indexOf("SEKJPY") != -1 || dataFiles[file_count].indexOf("NOKJPY") != -1)
//          {
//            scandi_jpy = true;
//            convers = .01;
//          }
//          else if(dataFiles[file_count].indexOf("NOK") != -1 || dataFiles[file_count].indexOf("SEK") != -1)
//          {
//            scandi = true;
//            convers = 1/6.0; 
//          }
         
         
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         while((strline = br.readLine()) != null)
         {
           daily_data.add(strline);
         }  
         
         for(line=0;line<daily_data.size()-1;line++)
         {
         
          strline = daily_data.get(line);
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
         

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(date_stamp.indexOf(insampStart) != -1)// && weekend.dayOfWeek().getAsText().equals("Wednesday"))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(log(new Double(tokens[1])));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(log(new Double(tokens[2]))); ask.add(log(new Double(tokens[3]))); mid.add(log(new Double(tokens[1])));
          }
          else
          {
           bid.add(log(new Double(tokens[2]))); ask.add(log(new Double(tokens[3]))); mid.add(log(new Double(tokens[1])));
          }
            
     
          D = new Double(tokens[1]); 
          price.add(log(D)); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
          if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
          {lo_price.add(new Double(tokens[2])); hi_price.add(new Double(tokens[2]));}
          else
          {lo_price.add((new Double(tokens[2]))); hi_price.add((new Double(tokens[2])));}
             
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && (date_stamp.indexOf(insampStart) != -1))// && weekend.dayOfWeek().getAsText().equals("Tuesday"))) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if((new Integer(intdates[1])).intValue() != cur_month || day_count == 0)
                {   
                   cur_month = (new Integer(intdates[1])).intValue();
                   //System.out.println("Recomputing filter..." + intdates[0] + "-" + intdates[1] + "-" + intdates[2]);
                   if(printall) System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, false);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall) System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                   
                   //System.out.println("\n");
                   //for(l=0;l<L;l++) {System.out.println(b_coeffs[l]);}
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                
                
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                svmout.println(date_stamp + ", " + current_signal);
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal);
                if(binary_rule)
                {
                
                
                 if(signal_profit)
                 {
                 
                   if((current_signal > 0 && in_transaction == 1) && (daily_price.get(daily_size-1) > last_price))
                   {
                   
                   
                     price_sold = daily_price.get(daily_size-1);
	             profit = price_sold - price_bought;
	             if(profit > 0) {succ_trades=succ_trades+1;}
	             total_trades=total_trades+1; 
	 
	             amount = getAmount(account.get(account.size()-1));
	             account.add(date_stamp + " " + (amount + profit));   
	             amount = getAmount(account.get(account.size()-1));
	            
	             made_trade = true;
	             dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                     log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	             in_transaction = 0;
	            
	             if(printall)
	             {
                     if(profit>0) System.out.println("Sold for a profit of " + profit);
                     else if(profit<0) System.out.println("Sold for a loss of " + profit);
                     System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                     }
                   
                     price_bought = daily_price.get(daily_size-1);
	             in_transaction = 1; 
	             last_price = daily_price.get(daily_size-1); 
	             if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	              
                   }
                   else if((current_signal < 0 && out_transaction == 1) && (daily_price.get(daily_size-1) < last_price))
                   {
                   
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    if(printall)
	             {
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);                   
                    }
                   
                    price_borrowed = daily_price.get(daily_size-1);
                    out_transaction = 1;	 
                  
                    if(printall)
	             {System.out.println("Entered short transaction at " + price_borrowed);}
                   
                   }
                 
                 }
                
                
                
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  last_price = daily_price.get(daily_size-1); 
                  waiting_meanrev_down = false;
                  waiting_meanrev_up = false;
                                 
                                 
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    if(printall)
	             {
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);}
                  }
                  
                  waiting_meanrev_down = true;
                  if(printall) {System.out.println("Waiting for better price to buy to enter long at " + last_price);}   
                  
                  
/*                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           
	           if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	          }*/ 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  last_price = daily_price.get(daily_size-1); 
                  waiting_meanrev_down = false;
                  waiting_meanrev_up = false;
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall)
	             {
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
                  waiting_meanrev_up = true;
                  if(printall) {System.out.println("Waiting for better price to buy to enter long at " + last_price);}   
               	         
	         
	         
// 	          if(short_sell && out_transaction == 0)
// 	          {
// 	           price_borrowed = daily_price.get(daily_size-1);
//                    out_transaction = 1;	 
//                   
//                    if(printall)
// 	           {
//                    System.out.println("Entered short transaction at " + price_borrowed);
//                    }
// 	          }
                 }
                }
                
                if(signal_strength_rule && (in_transaction == 0 && out_transaction == 0))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  last_price = daily_price.get(daily_size-1); 
                  waiting_meanrev_down = true;
                  
/*                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  
                    if(printall)
	             { 
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           
	           if(printall)
	             {System.out.println("Entered long transaction at " + price_bought);}
	          }*/ 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = daily_price.get(daily_size-1); 
                  waiting_meanrev_up = true;
                 
//                   if(long_buy && in_transaction == 1)
//                   {
// 	           price_sold = daily_price.get(daily_size-1);
// 	           profit = price_sold - price_bought;
// 	           if(profit > 0) {succ_trades=succ_trades+1;}
// 	           total_trades=total_trades+1; 
// 	 
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
// 
// 	            made_trade = true;
//                     dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                     log_ret = daily_price.get(daily_size-1).doubleValue();	            
// 	            
// 	           in_transaction = 0;
// 	           
// 	           if(printall){
//                    System.out.println("Bought for a profit of " + profit);
//                    System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);}	           
// 	           
// 	           
// 	          }
// 	         
// 	          if(short_sell && out_transaction == 0)
// 	          {
// 	           price_borrowed = daily_price.get(daily_size-1);
//                    out_transaction = 1;	 
//                   
//                    if(printall){System.out.println("Entered short transaction at " + price_borrowed);}
// 	          }
                 }              
                }
              
                if(!made_trade)
                {
                 account.add(date_stamp + " " + amount);
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {day_count=0;}
                
                allsignal.add(date_stamp + " " + current_signal);
                 
           }
           else if(trading_hours)
           {
                
             if(waiting_meanrev_up) //waiting for tick price to go up in order to sell the bid
             {
             
               if((bid.get(bid.size()-1) - last_price) > mean_rev_amnt)
               {
               
               	  if(short_sell && out_transaction == 0)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   last_price = price_borrowed;
                   if(printall)System.out.println("Entered short transaction at " + price_borrowed + " after " + (bid.get(bid.size()-1) - last_price) + " reversion");
	          }
                  waiting_meanrev_up = false; waiting_meanrev_down = false;
               }
             }
             else if(waiting_meanrev_down)
             {
               if((last_price - ask.get(ask.size() - 1)) > mean_rev_amnt)
               {
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           last_price = price_bought;
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                  waiting_meanrev_down = false; waiting_meanrev_up = false;
               }
             }            

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && price.get(price.size()-1) > last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
              
               cur_pnl = price.get(price.size()-1) - last_price;    
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall){System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);}
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 price_sold = price.get(price.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         //account.add(date_stamp + " " + (amount - stop_loss));   
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 
                 
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall){System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);}
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall){System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);}
                 
                 price_sold = price.get(price.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall){System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);}
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && price.get(price.size()-1) < last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
               cur_pnl =  last_price - price.get(price.size()-1);               
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall){System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);}
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         account.add(date_stamp + " " + (amount + profit));  
	         amount = getAmount(account.get(account.size()-1));

                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall){System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);}
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall){System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);}
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall){System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);}


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               //System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall){ System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
               + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl)); }
               
             
             allsignal.add(date_stamp + " " + current_signal);
               
           }      
                    
          }             
        }              
                    

 
      
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i)) - getAmount(account.get(i-1));
        
        dreturns[i] = convers*dreturns[i]; //approximate pips in dollars
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0]) && !isSunday(dates[0])) //first date entry
        {
         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
         }
         
         pnl = dreturns[i];
         perf_dates.add(dates[0]);
    
        }
        else //already contains the date, so add on pnl
        {pnl = pnl + dreturns[i];}
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       //System.out.println(perf_dates.get(i) + " " + perf_returns.get(i));
       date_returns.add(perf_dates.get(i) + " " + perf_returns.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
      }

      perform.close();
      // load in intraday information Date log-diff price signal - EXCLUDES LATEST OBSERVATION!!! 
      dates_price = new String[n_obs];
      for(i=0;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2]);
//        }
      }
      
      
      n_files++;
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff);
 
      out.close(); dailyout.close(); svmout.close(); b0_coeff.close();
      
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
 
      return computed; 
      
  }   
  

 
   
 
 
 
 
 
 
 
 
 
  
  
  
  
  
    public boolean startStrategyDailyIntradayOneHour()
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean print_filter = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      ArrayList<Double> perf_returns = new ArrayList<Double>();
      double pnl; 
      String[] dates; 
      boolean inverse_hours = false;
      String time; 
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> sunday = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      int end_hour;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      ArrayList<String> trade_times = new ArrayList<String>();
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      friday_closing = false;
      //asian_close = true;
      
      
      //----- Fill up with times here -------------------------
      for(i=0;i<10;i++) {trade_times.add("0"+i+":00:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":00:00");}
//       trade_times.add("00:00:00");
//       trade_times.add("06:00:00");
//       trade_times.add("12:00:00");
//       trade_times.add("18:00:00");     
      
      String[] hourToks = startingTime.split("[:]+");
      start_hour = (new Integer(hourToks[0])).intValue();
      
      hourToks = endingTime.split("[:]+");
      end_hour = (new Integer(hourToks[0])).intValue();
      
      inverse_hours = false;
      if(start_hour > end_hour)  //must switch hours here
      {
        int temphour = end_hour;
        
        inverse_hours = true;
        end_hour = start_hour;
        start_hour = temphour;
      }
      
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      
      String[] myname = dataFiles[0].split("[.]+");
      
      try{

        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results_"+myname[0]+".dat"));
        PrintWriter svmout = new PrintWriter(new FileWriter("neural.dat"));
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         
         if(dataFiles[file_count].indexOf("JPY") != -1 || dataFiles[file_count].indexOf("Y") != -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(futures_data)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          
          stop_loss_thresh = stop_loss_thresh*10000;
          stop_loss = stop_loss_thresh;
          take_profit_thresh = take_profit_thresh*10000;
          global_stop_loss = global_stop_loss*10000;
         }         
         else if(dataFiles[file_count].indexOf("MXN") != -1)
         {
          stop_loss_thresh = stop_loss_thresh*10;
          take_profit_thresh = take_profit_thresh*10;
          global_stop_loss = global_stop_loss*10;
          stop_loss = stop_loss_thresh;
         }
         else if(dataFiles[file_count].indexOf("HKD") != -1)
         {
          stop_loss_thresh = stop_loss_thresh*10;
          take_profit_thresh = take_profit_thresh*10;
          global_stop_loss = global_stop_loss*10;
          stop_loss = stop_loss_thresh;
         }
         else if(dataFiles[file_count].indexOf("XAU") != -1)
         {
          stop_loss_thresh = stop_loss_thresh*1000;
          take_profit_thresh = take_profit_thresh*1000;
          global_stop_loss = global_stop_loss*1000;
          stop_loss = stop_loss_thresh;
         }
         else if(dataFiles[file_count].indexOf("GC") != -1)
         {
          stop_loss_thresh = stop_loss_thresh*1000;
          take_profit_thresh = take_profit_thresh*1000;
          global_stop_loss = global_stop_loss*1000;
          stop_loss = stop_loss_thresh;
         }         
         else if(dataFiles[file_count].indexOf("JNJ") != -1)
         {
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }               
         
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         while((strline = br.readLine()) != null)
         {

          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
          time = date_tokens[1];

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(trade_times.contains(time))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(new Double(tokens[1]));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          print_filter = false;
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
          else
          {
           bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
            
 
          D = new Double(tokens[1]); 
          price.add(D); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
          if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
          {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
          else
          {lo_price.add((new Double(tokens[7]))); hi_price.add((new Double(tokens[8])));}
             
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          
          String[] hours = time.split("[:]+");
          cur_hour = (new Integer(hours[0])).intValue();
          (new Integer(hours[1])).intValue();
          
          trading_closed = false;
          
          //if currently not in a transaction and between the hours of midnight and start-hour, then no new 
          //positions will be opened
          
          
          if(asian_close)  //only closed if most recent transaction was closed artificially through SL or TP after end hours
          {
           if((in_transaction == 0 && out_transaction == 0) && cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true; if(printall) System.out.println(cur_hour + " " + start_hour);}
          }
          else 
          {
           if(cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true; if(printall) {System.out.println(cur_hour + " " + start_hour);}}          
          }
          
//           if(cur_hour == 14) {release_first = true;}
//           else {elease_first = false;}
          
          if(inverse_hours) {trading_closed = !trading_closed;}
          
          its_closing_time = (weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1);
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && trade_times.contains(time)) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if(day_count == 0 || weekend.dayOfWeek().getAsText().equals("Sunday") && date_stamp.indexOf("18:00:00") != -1)
                {   
                   
                   if(printall)System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, print_filter);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall)System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                if(date_stamp.indexOf("16:00:00") != -1) {svmout.println(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)));}
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                //if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                //if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal + ", trading_closed = " + trading_closed);
               if(friday_closing && its_closing_time)
               {
                 if(printall) System.out.println("\nIt's Friday at 5pm, time to close shop for week");
                 if(current_signal > 0 && in_transaction == 1) //in a long transaction
                 {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));
	            account.add(date_stamp + " " + daily_price.get(daily_price.size()-1) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit)); 
	            amount = getAmount(account.get(account.size()-1),4);
	            
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                  }
                 }
                 else if(current_signal < 0 && out_transaction == 1) //in a short transaction
                 {
                 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " " + daily_price.get(daily_price.size()-1) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                 }
                 
               }
               else 
               { 
                //if(binary_rule && !trading_closed)
                if(binary_rule)
                {
                   
                
                 made_trade = false;
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  last_price = daily_price.get(daily_size-1); 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));  
	            account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    
                    if(printall){
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if((long_buy && in_transaction == 0) && !trading_closed)
	          {
	           profit = 0;
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	           dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	          if((short_sell && out_transaction == 0) && !trading_closed)
	          {
	           profit = 0;
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                   account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                   
                   dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue(); 
                   
                   if(printall) System.out.println("Entered short transaction at " + price_borrowed);
                  
	          }
                 }
                }
                
                if(signal_strength_rule && ((in_transaction == 0 && out_transaction == 0) && !trading_closed))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  last_price = daily_price.get(daily_size-1); 
                  
                  
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
	           profit = 0;
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	           dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();   
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
	           in_transaction = 0;
	           
                   if(printall)System.out.println("Bought for a profit of " + profit);
                   if(printall)System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           profit = 0;
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                   account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                  
                   dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();   
                  
                   if(printall)System.out.println("Entered short transaction at " + price_borrowed);
	          }
                 }              
                }
              
                if(!made_trade)
                {
                 profit=0;
                 //account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                 //dailyoutret.add(price.get(price.size()-1) - log_ret);
                 //log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {y_count=0;}
                if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));             
                allsignal.add(date_stamp + " " + current_signal);
            }
           
           }
           else if(trading_hours)// && !trading_closed)// && cur_min != 30)
           {
                
            

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && price.get(price.size()-1) > last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
              
               //cur_pnl = price.get(price.size()-1) - last_price;    
               cur_pnl = bid.get(bid.size()-1) - last_price;
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 //price_sold = price.get(price.size()-1);
	         price_sold = bid.get(bid.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount - stop_loss));   
	         //account.add(date_stamp + " " + (amount + profit));   
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " -" + formatter3.format(bid.get(bid.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
                 if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 
                 dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                 log_ret = bid.get(bid.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = bid.get(bid.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount + profit));   
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " -" + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && price.get(price.size()-1) < last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
               //cur_pnl =  last_price - price.get(price.size()-1);               
               cur_pnl =  last_price - ask.get(ask.size()-1); 
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         //account.add(date_stamp + " " + (amount + profit));  
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " " + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                if(printall) System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size() - 1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount + profit));   
	         account.add(date_stamp + " " + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	         
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall)System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 if(printall)System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 if(printall)System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 if(printall)System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               if(printall)System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl)); 
           
             allsignal.add(date_stamp + " " + current_signal);
           }      
                    
          }             
        }              
                    

 
      double mean_ntrades = 0;
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      int n_trades = 0;
      ArrayList<Integer> n_trades_day = new ArrayList<Integer>();
      
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i),4) - getAmount(account.get(i-1),4);
      
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0])) //first date entry
        {
         
         

         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
           n_trades_day.add(n_trades);
         }
         
         
         perf_dates.add(dates[0]);
         pnl = dreturns[i];
         if(dreturns[i] != 0) {n_trades = 1;} 
         else {n_trades = 0;}
        }
        else //already contains the date, so add on pnl
        {
         pnl = pnl + dreturns[i]; 
         //System.out.println(dreturns[i]);
         if(dreturns[i] != 0) {n_trades++;}// System.out.println(n_trades);}
        }
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      n_trades_day.add(n_trades);
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       if(printall) System.out.println(perf_dates.get(i) + " " + perf_returns.get(i) + " " + n_trades_day.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
       mean_ntrades = mean_ntrades + n_trades_day.get(i);
      }
      mean_ntrades = mean_ntrades/n_trades_day.size();
      
  
      dates_price = new String[n_obs];
      for(i=0;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2]);
//        }
      }  
  
      double sunday_roi  = 0; double sunday_trade = 0; 
      double sunday_pos_mean=0; double sunday_neg_mean=0; 
      int nsunday_pos_mean=0; int nsunday_neg_mean=0;
      
      System.out.println("\nSunday Performance");
      for(i=0;i<sunday.size();i++)
      {
        //System.out.println(sunday.get(i));
        sunday_trade = getAmount(sunday.get(i),4);
        sunday_roi = sunday_roi + sunday_trade;
        
        
        
        if(sunday_trade > 0)
        {sunday_pos_mean += sunday_trade; nsunday_pos_mean++;}
        if(sunday_trade < 0)
        {sunday_neg_mean += sunday_trade; nsunday_neg_mean++;}
      }
      
      System.out.println("\nSunday Stats");
      System.out.println("Sunday ROI = " + sunday_roi);
      sunday_roi = sunday_roi/(nsunday_pos_mean+nsunday_neg_mean);
      System.out.println("Sunday meantrade = " + sunday_roi);
      System.out.println("Sunday pos mean = " + (sunday_pos_mean/nsunday_pos_mean) + ", npos_trades = " + nsunday_pos_mean);
      System.out.println("Sunday neg mean = " + (sunday_neg_mean/nsunday_neg_mean) + ", nneg_trades = " + nsunday_neg_mean);
  
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      double risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      double reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff + ", avg_n_trades = " + mean_ntrades);
 
      out.close(); dailyout.close(); perform.close(); svmout.close(); b0_coeff.close();
      
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
      n_files++;
      return computed; 
      
  }    
  
  
  
  
    public boolean startStrategyDailyIntradayStocks()
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean print_filter = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      ArrayList<Double> perf_returns = new ArrayList<Double>();
      double pnl; 
      String[] dates; 
      boolean inverse_hours = false;
      String time; 
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> sunday = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      int end_hour;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      ArrayList<String> trade_times = new ArrayList<String>();
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      friday_closing = true;
      //asian_close = true;
      
      
      //----- Fill up with times here -------------------------
//       for(i=0;i<10;i++) 
//       {
//        if(i%4 == 0) {trade_times.add("0"+i+":00:00");}
//       }
//       for(i=10;i<24;i++) 
//       {
//       trade_times.add(i+":00:00");}
//       for(i=0;i<10;i++) {trade_times.add("0"+i+":00:00");}
//       for(i=10;i<24;i++) {trade_times.add(i+":00:00");} 

      trade_times.add("10:00:00");
      trade_times.add("11:00:00");
      trade_times.add("12:00:00");
      trade_times.add("13:00:00");
      trade_times.add("14:00:00");
      trade_times.add("15:00:00");
      
     
     
      String[] hourToks = startingTime.split("[:]+");
      start_hour = (new Integer(hourToks[0])).intValue();
      
      hourToks = endingTime.split("[:]+");
      end_hour = (new Integer(hourToks[0])).intValue();
      
      inverse_hours = false;
      if(start_hour > end_hour)  //must switch hours here
      {
        int temphour = end_hour;
        
        inverse_hours = true;
        end_hour = start_hour;
        start_hour = temphour;
      }
      
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      
      String[] myname = dataFiles[0].split("[.]+");
      
      try{
    
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results_"+myname[0]+".dat"));
        
        
        
        for(file_count=0;file_count<1;file_count++)
        {
          
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         
         
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         while((strline = br.readLine()) != null)
         {

          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
          time = date_tokens[1];

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(trade_times.contains(time))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(new Double(tokens[1]));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          print_filter = false;
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); 
          close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
          else
          {
           bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
            
         
          D = new Double(tokens[1]); 
          price.add(D); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
          if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
          {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
          else
          {lo_price.add((new Double(tokens[7]))); hi_price.add((new Double(tokens[8])));}
             
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          
          String[] hours = time.split("[:]+");
          cur_hour = (new Integer(hours[0])).intValue();
          (new Integer(hours[1])).intValue();
          
          trading_closed = false;
          
          //if currently not in a transaction and between the hours of midnight and start-hour, then no new 
          //positions will be opened
          
          
          if(asian_close)  //only closed if most recent transaction was closed artificially through SL or TP after end hours
          {
           if((in_transaction == 0 && out_transaction == 0) && cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true; if(printall) System.out.println(cur_hour + " " + start_hour);}
          }
          else 
          {
           if(cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true; if(printall) {System.out.println(cur_hour + " " + start_hour);}}          
          }
          
//           if(cur_hour == 14) {release_first = true;}
//           else {elease_first = false;}
          
          if(inverse_hours) {trading_closed = !trading_closed;}
          
          its_closing_time = false;
          //its_closing_time = (date_stamp.indexOf("16:00:00") != -1);
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && trade_times.contains(time)) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if(day_count == 0 || weekend.dayOfWeek().getAsText().equals("Sunday") && date_stamp.indexOf("18:00:00") != -1)
                {   
                   
                   if(printall)System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, print_filter);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall)System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                //if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                //if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal + ", trading_closed = " + trading_closed);
               if(friday_closing && its_closing_time)
               {
                 if(printall) System.out.println("\nIt's Friday at 5pm, time to close shop for week");
                 if(current_signal > 0 && in_transaction == 1) //in a long transaction
                 {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));
	            account.add(date_stamp + " " + daily_price.get(daily_price.size()-1) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit)); 
	            amount = getAmount(account.get(account.size()-1),4);
	            
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                  }
                 }
                 else if(current_signal < 0 && out_transaction == 1) //in a short transaction
                 {
                 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " " + daily_price.get(daily_price.size()-1) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                 }
                 
               }
               else 
               { 
                //if(binary_rule && !trading_closed)
                if(binary_rule && (cur_hour <= 16 && cur_hour >= 10))
                {
                
                
//                  if(signal_profit)
//                  {
//                    //if(printall) System.out.println("Trading_closed = " + trading_closed);
//                    if((current_signal > 0 && in_transaction == 1) && (daily_price.get(daily_size-1) > last_price))
//                    {
//                    
//                    
//                      price_sold = daily_price.get(daily_size-1);
// 	             profit = price_sold - price_bought;
// 	             if(profit > 0) {succ_trades=succ_trades+1;}
// 	             total_trades=total_trades+1; 
// 	 
// 	             amount = getAmount(account.get(account.size()-1));
// 	             account.add(date_stamp + " " + (amount + profit));   
// 	             amount = getAmount(account.get(account.size()-1));
// 	            
// 	             made_trade = true;
// 	             dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                      log_ret = daily_price.get(daily_size-1).doubleValue();
// 	 
// 	             in_transaction = 0;
// 	           
// 	             if(printall){
//                      if(profit>0) System.out.println("Sold for a profit of " + profit);
//                      else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                       System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
//                      }
//                    
//                      if(!trading_closed)
//                      {price_bought = daily_price.get(daily_size-1);
// 	             //in_transaction = 1; 
// 	             out_transaction = 1;
// 	             last_price = daily_price.get(daily_size-1); 
// 	             if(printall) System.out.println("Entered long transaction at " + price_bought);
// 	             } 
//                    }
//                    else if((current_signal < 0 && out_transaction == 1) && (daily_price.get(daily_size-1) < last_price))
//                    {
//                    
//                     price_sold = daily_price.get(daily_size-1);
//                     profit = price_borrowed - price_sold;
// 
// 	            if(profit > 0) {succ_trades=succ_trades+1;}
// 	            total_trades=total_trades+1; 
// 	  
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
//                     
//                     made_trade = true;
//                     dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                     log_ret = daily_price.get(daily_size-1).doubleValue();
//                     
//                     
//                     out_transaction = 0;
//                     
//                     if(printall)
//                     {
//                     if(profit>0) System.out.println("Sold for a profit of " + profit);
//                     else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                     }
//                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);                   
//                    
//                     if(!trading_closed)
//                     {
//                     price_borrowed = daily_price.get(daily_size-1);
//                     //out_transaction = 1;	 
//                     in_transaction = 1;
//                     if(printall)System.out.println("Entered short transaction at " + price_borrowed);
//                     }
//                    }
//                  
//                  }
                
                
                 made_trade = false;
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  last_price = daily_price.get(daily_size-1); 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));  
	            account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                    
                    made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    
                    if(printall){
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if((long_buy && in_transaction == 0) && !trading_closed)
	          {
	           profit = 0;
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	           dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	          if((short_sell && out_transaction == 0) && !trading_closed)
	          {
	           profit = 0;
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                   account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                   
                   dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue(); 
                   
                   if(printall) System.out.println("Entered short transaction at " + price_borrowed);
                  
	          }
                 }
                }
                
                if(signal_strength_rule && ((in_transaction == 0 && out_transaction == 0) && !trading_closed && (cur_hour <= 16 && cur_hour >= 10)))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  last_price = daily_price.get(daily_size-1); 
                  
                  
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
	           profit = 0;
                   price_bought = daily_price.get(daily_size-1);
	           in_transaction = 1; 
	           account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	           dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();   
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = daily_price.get(daily_size-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1),4);
	            //account.add(date_stamp + " " + (amount + profit));   
	            account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	            amount = getAmount(account.get(account.size()-1),4);
	            if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
	           in_transaction = 0;
	           
                   if(printall)System.out.println("Bought for a profit of " + profit);
                   if(printall)System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           profit = 0;
	           price_borrowed = daily_price.get(daily_size-1);
                   out_transaction = 1;	 
                   account.add(date_stamp + " -" + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                  
                   dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                   log_ret = daily_price.get(daily_size-1).doubleValue();   
                  
                   if(printall)System.out.println("Entered short transaction at " + price_borrowed);
	          }
                 }              
                }
              
                if(!made_trade)
                {
                 profit=0;
                 //account.add(date_stamp + " " + formatter3.format(daily_price.get(daily_price.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
                 //dailyoutret.add(price.get(price.size()-1) - log_ret);
                 //log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {y_count=0;}
                if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));             
                allsignal.add(date_stamp + " " + current_signal);
            }
           
           }
           else if(trading_hours && (cur_hour <= 16 && cur_hour >= 10))// && !trading_closed)// && cur_min != 30)
           {
                
            

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && price.get(price.size()-1) > last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
              
               //cur_pnl = price.get(price.size()-1) - last_price;    
               cur_pnl = bid.get(bid.size()-1) - last_price;
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 //price_sold = price.get(price.size()-1);
	         price_sold = bid.get(bid.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount - stop_loss));   
	         //account.add(date_stamp + " " + (amount + profit));   
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " -" + formatter3.format(bid.get(bid.size()-1)) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
                 if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 
                 dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                 log_ret = bid.get(bid.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = bid.get(bid.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount + profit));   
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " -" + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && price.get(price.size()-1) < last_price) //check if new high price
               {last_price = price.get(price.size()-1);}
                
               //cur_pnl =  last_price - price.get(price.size()-1);               
               cur_pnl =  last_price - ask.get(ask.size()-1); 
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         //account.add(date_stamp + " " + (amount + profit));  
	         //account.add(date_stamp + " " + price.get(price.size()-1) + " " + profit + " " + (amount + profit));
	         account.add(date_stamp + " " + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                if(printall) System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size() - 1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count),4);
	         //account.add(date_stamp + " " + (amount + profit));   
	         account.add(date_stamp + " " + formatter3.format(price_sold) + " " + formatter3.format(profit) + " " + formatter3.format(amount + profit));
	         amount = getAmount(account.get(account.size()-1),4);
	         if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	         
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall)System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 if(printall)System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 if(printall)System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 if(printall)System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               if(printall)System.out.println("End of week");
            
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl)); 
           
             allsignal.add(date_stamp + " " + current_signal);
           }      
                    
          }             
        }              
                    

 
      double mean_ntrades = 0;
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      int n_trades = 0;
      ArrayList<Integer> n_trades_day = new ArrayList<Integer>();
      
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i),4) - getAmount(account.get(i-1),4);
      
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0])) //first date entry
        {
         
         

         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
           n_trades_day.add(n_trades);
         }
         
         
         perf_dates.add(dates[0]);
         pnl = dreturns[i];
         if(dreturns[i] != 0) {n_trades = 1;} 
         else {n_trades = 0;}
        }
        else //already contains the date, so add on pnl
        {
         pnl = pnl + dreturns[i]; 
         //System.out.println(dreturns[i]);
         if(dreturns[i] != 0) {n_trades++;}// System.out.println(n_trades);}
        }
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      n_trades_day.add(n_trades);
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       if(printall) System.out.println(perf_dates.get(i) + " " + perf_returns.get(i) + " " + n_trades_day.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
       mean_ntrades = mean_ntrades + n_trades_day.get(i);
      }
      mean_ntrades = mean_ntrades/n_trades_day.size();
      
  
      dates_price = new String[n_obs];
      for(i=0;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2]);
//        }
      }  
  
      double sunday_roi  = 0; double sunday_trade = 0; 
      double sunday_pos_mean=0; double sunday_neg_mean=0; 
      int nsunday_pos_mean=0; int nsunday_neg_mean=0;
      
      System.out.println("\nSunday Performance");
      for(i=0;i<sunday.size();i++)
      {
        //System.out.println(sunday.get(i));
        sunday_trade = getAmount(sunday.get(i),4);
        sunday_roi = sunday_roi + sunday_trade;
        
        
        
        if(sunday_trade > 0)
        {sunday_pos_mean += sunday_trade; nsunday_pos_mean++;}
        if(sunday_trade < 0)
        {sunday_neg_mean += sunday_trade; nsunday_neg_mean++;}
      }
      
      System.out.println("\nSunday Stats");
      System.out.println("Sunday ROI = " + sunday_roi);
      sunday_roi = sunday_roi/(nsunday_pos_mean+nsunday_neg_mean);
      System.out.println("Sunday meantrade = " + sunday_roi);
      System.out.println("Sunday pos mean = " + (sunday_pos_mean/nsunday_pos_mean) + ", npos_trades = " + nsunday_pos_mean);
      System.out.println("Sunday neg mean = " + (sunday_neg_mean/nsunday_neg_mean) + ", nneg_trades = " + nsunday_neg_mean);
  
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      double risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      double reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff + ", avg_n_trades = " + mean_ntrades);
 
      out.close(); dailyout.close(); perform.close(); b0_coeff.close();
      
      }
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
      n_files++;
      return computed; 
      
  }    
  
  
   public boolean startStrategyDailyIntradayOneHourBidAsk()
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean print_filter = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      ArrayList<Double> perf_returns = new ArrayList<Double>();
      double pnl; 
      String[] dates; 
      boolean inverse_hours = false;
      String time; 
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      int end_hour;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      ArrayList<String> trade_times = new ArrayList<String>();
      ArrayList<String> fxcm_dates = new ArrayList<String>();
      
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      friday_closing = false;
      //asian_close = true;
      
      
      //----- Fill up with times here -------------------------
      for(i=0;i<10;i++) {trade_times.add("0"+i+":00:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":00:00");}
//       trade_times.add("00:00:00");
//       trade_times.add("06:00:00");
//       trade_times.add("12:00:00");
//       trade_times.add("18:00:00");     
      
      String[] hourToks = startingTime.split("[:]+");
      start_hour = (new Integer(hourToks[0])).intValue();
      
      hourToks = endingTime.split("[:]+");
      end_hour = (new Integer(hourToks[0])).intValue();
      
      inverse_hours = false;
      if(start_hour > end_hour)  //must switch hours here
      {
        int temphour = end_hour;
        
        inverse_hours = true;
        end_hour = start_hour;
        start_hour = temphour;
      }
      
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      

      jpy = false;
      try{
        
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
      
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         
         if(dataFiles[file_count].indexOf("JPY") != -1 && dataFiles[file_count].indexOf("NOKJPY") == -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(futures_data)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          
          stop_loss_thresh = stop_loss_thresh*10000;
          stop_loss = stop_loss_thresh;
          take_profit_thresh = take_profit_thresh*10000;
          global_stop_loss = global_stop_loss*10000;
         }         
         
         ArrayList<String> dayData = new ArrayList<String>();
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         
         
         //---- add the latest FXCM data -
         while((strline = br.readLine()) != null)
         {
           dayData.add(strline);
           
           //add the date
           tokens = strline.split("[,]+");
           fxcm_dates.add(tokens[0]);
           
         }
         
         //----- now add the lates IB data ------------
         System.out.println("opening " + dataFiles[file_count] + " IB data");
         String[] fxpair = dataFiles[file_count].split("[.]+");
         
         if((new File(fxpair[0] + ".IB.dat")).exists())
         {
          fin = new FileInputStream(fxpair[0] + ".IB.dat");
          din = new DataInputStream(fin);
          br = new BufferedReader(new InputStreamReader(din)); 
         
          String lastFXCMdate = fxcm_dates.get(fxcm_dates.size()-1);
         
          SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-DD HH:mm:ss");
          Date date1 = sdf.parse(lastFXCMdate);

         
         
          while((strline = br.readLine()) != null)
          {
           tokens = strline.split("[,]+");
           
           Date date2 = sdf.parse(tokens[0]);
           
           if(!fxcm_dates.contains(tokens[0]) && date1.before(date2))
           {
             System.out.println("New time and observation " + strline);
             dayData.add(strline);
           }
          }
         }
         
         for(int ts = 0; ts < dayData.size(); ts++)
         {
          strline = dayData.get(ts);
          
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
          time = date_tokens[1];

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(trade_times.contains(time))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(new Double(tokens[1]));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          print_filter = false;
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
          else
          {
           bid.add(new Double(tokens[2])); ask.add(new Double(tokens[3])); mid.add(new Double(tokens[1]));
          }
            

          D = new Double(tokens[1]); 
          price.add(D); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
//           if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
//           {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
//           else
//           {lo_price.add((new Double(tokens[7]))); hi_price.add((new Double(tokens[8])));}
          lo_price.add((new Double(tokens[2]))); hi_price.add((new Double(tokens[3])));   
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          
          String[] hours = time.split("[:]+");
          cur_hour = (new Integer(hours[0])).intValue();
          (new Integer(hours[1])).intValue();
          
          trading_closed = false;
          
          //if currently not in a transaction and between the hours of midnight and start-hour, then no new 
          //positions will be opened
          
          
          if(asian_close)  //only closed if most recent transaction was closed artificially through SL or TP after end hours
          {
           if((in_transaction == 0 && out_transaction == 0) && cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;}// if(printall) System.out.println(cur_hour + " " + start_hour);}
          }
          else 
          {
           if(cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;} // if(printall) {System.out.println(cur_hour + " " + start_hour);}}          
          }
          
          if(inverse_hours) {trading_closed = !trading_closed;}
          
          its_closing_time = (weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1);
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && trade_times.contains(time)) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if(day_count == 0 || weekend.dayOfWeek().getAsText().equals("Sunday") && date_stamp.indexOf("18:00:00") != -1)
                {   
                   
                   if(printall)System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, print_filter);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall)System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                //if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                //if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal + ", trading_closed = " + trading_closed);
               if(friday_closing && its_closing_time)
               {
                 if(printall) System.out.println("\nIt's Friday at 5pm, time to close shop for week");
                 if(current_signal > 0 && in_transaction == 1) //in a long transaction
                 {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                  }
                 }
                 else if(current_signal < 0 && out_transaction == 1) //in a short transaction
                 {
                 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                 }
                 
               }
               else 
               { 
                //if(binary_rule && !trading_closed)
                if(binary_rule)
                {
                
                
//                  if(signal_profit)
//                  {
//                    //if(printall) System.out.println("Trading_closed = " + trading_closed);
//                    if((current_signal > 0 && in_transaction == 1) && (daily_price.get(daily_size-1) > last_price))
//                    {
//                    
//                    
//                      price_sold = daily_price.get(daily_size-1);
// 	             profit = price_sold - price_bought;
// 	             if(profit > 0) {succ_trades=succ_trades+1;}
// 	             total_trades=total_trades+1; 
// 	 
// 	             amount = getAmount(account.get(account.size()-1));
// 	             account.add(date_stamp + " " + (amount + profit));   
// 	             amount = getAmount(account.get(account.size()-1));
// 	            
// 	             made_trade = true;
// 	             dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                      log_ret = daily_price.get(daily_size-1).doubleValue();
// 	 
// 	             in_transaction = 0;
// 	           
// 	             if(printall){
//                      if(profit>0) System.out.println("Sold for a profit of " + profit);
//                      else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                       System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
//                      }
//                    
//                      if(!trading_closed)
//                      {price_bought = daily_price.get(daily_size-1);
// 	             //in_transaction = 1; 
// 	             out_transaction = 1;
// 	             last_price = daily_price.get(daily_size-1); 
// 	             if(printall) System.out.println("Entered long transaction at " + price_bought);
// 	             } 
//                    }
//                    else if((current_signal < 0 && out_transaction == 1) && (daily_price.get(daily_size-1) < last_price))
//                    {
//                    
//                     price_sold = daily_price.get(daily_size-1);
//                     profit = price_borrowed - price_sold;
// 
// 	            if(profit > 0) {succ_trades=succ_trades+1;}
// 	            total_trades=total_trades+1; 
// 	  
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
//                     
//                     made_trade = true;
//                     dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                     log_ret = daily_price.get(daily_size-1).doubleValue();
//                     
//                     
//                     out_transaction = 0;
//                     
//                     if(printall)
//                     {
//                     if(profit>0) System.out.println("Sold for a profit of " + profit);
//                     else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                     }
//                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);                   
//                    
//                     if(!trading_closed)
//                     {
//                     price_borrowed = daily_price.get(daily_size-1);
//                     //out_transaction = 1;	 
//                     in_transaction = 1;
//                     if(printall)System.out.println("Entered short transaction at " + price_borrowed);
//                     }
//                    }
//                  
//                  }
                
                
                 made_trade = false;
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);                 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    //price_sold = daily_price.get(daily_size-1);
                    price_sold = ask.get(ask.size()-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(ask.get(ask.size()-1) - log_ret);
                    log_ret = ask.get(ask.size()-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    
                    if(printall){
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if((long_buy && in_transaction == 0) && !trading_closed)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = bid.get(bid.size()-1);
                  
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = bid.get(bid.size()-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                    log_ret = bid.get(bid.size()-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	          if((short_sell && out_transaction == 0) && !trading_closed)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   if(printall) System.out.println("Entered short transaction at " + price_borrowed);
                  
	          }
                 }
                }
                
                if(signal_strength_rule && ((in_transaction == 0 && out_transaction == 0) && !trading_closed))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);
                  
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = ask.get(ask.size()-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(ask.get(ask.size()-1) - log_ret);
                    log_ret = ask.get(ask.size()-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = bid.get(bid.size()-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = bid.get(bid.size()-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                    log_ret = bid.get(bid.size()-1).doubleValue();	            
	            
	           in_transaction = 0;
	           
                   if(printall)System.out.println("Bought for a profit of " + profit);
                   if(printall)System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   if(printall)System.out.println("Entered short transaction at " + price_borrowed);
	          }
                 }              
                }
              
                if(!made_trade)
                {
                 account.add(date_stamp + " " + amount);
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {day_count=0;}
                //if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));             
                allsignal.add(date_stamp + " " + current_signal);
            }
           
           }
           if(trading_hours)// && !trading_closed)// && cur_min != 30)
           {
                
             if(cur_pnl != 0 && date_stamp.indexOf("22:30:00") != -1) //close out position
             {
               if(in_transaction == 1)
               {
               
                 //System.out.println("Ending trading for next startup");

                 
                 price_sold = price.get(price.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));

	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 
                 
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);                 
 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         //System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;           
               
               }
               else if(out_transaction == 1)
               {
               
                 //System.out.println("Ending trading for next startup");
               
               
                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         account.add(date_stamp + " " + (amount + profit));  
	         amount = getAmount(account.get(account.size()-1));

                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
               
               }
             }

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && bid.get(bid.size()-1) > last_price) //check if new high price
               {last_price = bid.get(bid.size()-1);}
                
              
               //cur_pnl = price.get(price.size()-1) - last_price;    
               cur_pnl = bid.get(bid.size()-1) - last_price;
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 //price_sold = price.get(price.size()-1);
	         price_sold = bid.get(bid.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 
                 dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                 log_ret = bid.get(bid.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = bid.get(bid.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && ask.get(ask.size()-1) < last_price) //check if new high price
               {last_price = ask.get(ask.size()-1);}
                
               //cur_pnl =  last_price - price.get(price.size()-1);               
               cur_pnl =  last_price - ask.get(ask.size()-1); 
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                if(printall) System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size() - 1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  	         
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));	         
	         
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	         
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall)System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 if(printall)System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 if(printall)System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 if(printall)System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               if(printall)System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall) 
             {
               if(current_signal > 0) {System.out.println(""+date_stamp + ", " + formatter3.format(bid.get(bid.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
               else {System.out.println(""+date_stamp + ", " + formatter3.format(ask.get(ask.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
             }
           
             allsignal.add(date_stamp + " " + current_signal);
           }      
                    
          }             
        }                        
          
          
           

      double mean_ntrades = 0;
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      int n_trades = 0;
      ArrayList<Integer> n_trades_day = new ArrayList<Integer>();
      double pret;
      ArrayList<Double> perf_rets = new ArrayList<Double>();
      ArrayList<String> perf_datestimes = new ArrayList<String>();
      
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i)) - getAmount(account.get(i-1));
      
        if(jpy) 
        {dreturns[i] = dreturns[i]*.01;}
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0])) //first date entry
        {
         
         

         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
           n_trades_day.add(n_trades);
         }
         
         
         perf_dates.add(dates[0]);
         if(dreturns[i] != 0)
         {
          perf_datestimes.add(dates[0] + " " + dates[1]);
          perf_rets.add(dreturns[i]);
         }
         pnl = dreturns[i];
         if(dreturns[i] != 0) {n_trades = 1;} 
         else {n_trades = 0;}
        }
        else //already contains the date, so add on pnl
        {
         pnl = pnl + dreturns[i]; 
         //System.out.println(dreturns[i]);
         if(dreturns[i] != 0) 
         {
           n_trades++;
           perf_datestimes.add(dates[0] + " " + dates[1]);
           perf_rets.add(dreturns[i]);                   
         }// System.out.println(n_trades);}
        }
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      n_trades_day.add(n_trades);
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       if(printall) System.out.println(perf_dates.get(i) + " " + perf_returns.get(i) + " " + n_trades_day.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
       date_returns.add(perf_dates.get(i) + " " + perf_returns.get(i));
       mean_ntrades = mean_ntrades + n_trades_day.get(i);
      }
      mean_ntrades = mean_ntrades/n_trades_day.size();
      
      
       
      
      
      dates_price = new String[n_obs];
      
      for(i=1;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        
        if(perf_datestimes.contains(latestDates.get(latestDates.size() - n_obs + i)))
        {
          pret = perf_rets.get(perf_datestimes.indexOf(latestDates.get(latestDates.size() - n_obs + i)));
          System.out.println(latestDates.get(latestDates.size() - n_obs + i) + " " + pret);
        }
        else
        {pret = 0;}
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2] + " " + pret);
//        }
      }  
      dates_price[0] = dates_price[1];
  
  
  
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      double risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      double reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff + ", avg_n_trades = " + mean_ntrades);
 
      out.close(); dailyout.close(); perform.close(); b0_coeff.close();
      
      }
      catch(ParseException fe){System.out.println("ParseException");}
      catch(NullPointerException fe){System.out.println("Null pointer");}
      catch(IllegalArgumentException fe){System.out.println("IllegalArgument");}
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
      n_files++;
     
      return computed; 
      
  }                   
  
    
  
  
  
   public boolean startStrategyDailyIntradayOneHourBidAsk15Min(double sprd, boolean spreadon)
   {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean print_filter = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      ArrayList<Double> perf_returns = new ArrayList<Double>();
      double pnl; 
      String[] dates; 
      boolean inverse_hours = false;
      String time; 
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      int end_hour;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      ArrayList<String> trade_times = new ArrayList<String>();
      ArrayList<String> fxcm_dates = new ArrayList<String>();
      
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      friday_closing = false;
      //asian_close = true;
      
      boolean spread_on = spreadon;
      
      //----- Fill up with times here -------------------------
      for(i=0;i<10;i++) {trade_times.add("0"+i+":00:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":00:00");}
      for(i=0;i<10;i++) {trade_times.add("0"+i+":15:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":15:00");}
      for(i=0;i<10;i++) {trade_times.add("0"+i+":30:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":30:00");}
      for(i=0;i<10;i++) {trade_times.add("0"+i+":45:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":45:00");}      
      
//       trade_times.add("00:00:00");
//       trade_times.add("06:00:00");
//       trade_times.add("12:00:00");
//       trade_times.add("18:00:00");     
      
      String[] hourToks = startingTime.split("[:]+");
      start_hour = (new Integer(hourToks[0])).intValue();
      
      hourToks = endingTime.split("[:]+");
      end_hour = (new Integer(hourToks[0])).intValue();
      
      inverse_hours = false;
      if(start_hour > end_hour)  //must switch hours here
      {
        int temphour = end_hour;
        
        inverse_hours = true;
        end_hour = start_hour;
        start_hour = temphour;
      }
      
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      

      jpy = false;
      try{
      
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
       
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         
         if(dataFiles[file_count].indexOf("JPY") != -1 && dataFiles[file_count].indexOf("NOKJPY") == -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(futures_data)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          
          stop_loss_thresh = stop_loss_thresh*10000;
          stop_loss = stop_loss_thresh;
          take_profit_thresh = take_profit_thresh*10000;
          global_stop_loss = global_stop_loss*10000;
         }         
         
         ArrayList<String> dayData = new ArrayList<String>();
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         
         
         //---- add the latest FXCM data -
         while((strline = br.readLine()) != null)
         {
           dayData.add(strline);
           
           //add the date
           tokens = strline.split("[,]+");
           fxcm_dates.add(tokens[0]);
           
         }
         
         //----- now add the lates IB data ------------
         System.out.println("opening " + dataFiles[file_count] + " IB data");
         String[] fxpair = dataFiles[file_count].split("[.]+");
         
         if((new File(fxpair[0] + ".IB.dat")).exists())
         {
          fin = new FileInputStream(fxpair[0] + ".IB.dat");
          din = new DataInputStream(fin);
          br = new BufferedReader(new InputStreamReader(din)); 
         
          String lastFXCMdate = fxcm_dates.get(fxcm_dates.size()-1);
         
          SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-DD HH:mm:ss");
          Date date1 = sdf.parse(lastFXCMdate);

         
         
          while((strline = br.readLine()) != null)
          {
           tokens = strline.split("[,]+");
           
           Date date2 = sdf.parse(tokens[0]);
           
           if(!fxcm_dates.contains(tokens[0]) && date1.before(date2))
           {
             System.out.println("New time and observation " + strline);
             dayData.add(strline);
           }
          }
         }
         
         for(int ts = 0; ts < dayData.size(); ts++)
         {
          strline = dayData.get(ts);
          
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
          time = date_tokens[1];

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(trade_times.contains(time))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(new Double(tokens[1]));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          print_filter = false;
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add((new Double(tokens[2])) - sprd); ask.add((new Double(tokens[3])) + sprd); mid.add(new Double(tokens[1]));
          }
          else
          {
           bid.add((new Double(tokens[2])) - sprd); ask.add((new Double(tokens[3])) + sprd); mid.add(new Double(tokens[1]));
          }
            
     
          D = new Double(tokens[1]); 
          //price.add(D); 
              
              
            if(spread_on)
            {
             if(current_signal > 0)
             {price.add(ask.get(ask.size()-1));}
             else
             {price.add(bid.get(bid.size()-1));}

//              if(current_signal > 0)
//              {price.add(mid.get(mid.size()-1) - sprd );}
//              else
//              {price.add(mid.get(mid.size()-1) + sprd );}


            }
            else
            {
              price.add(mid.get(mid.size()-1));
            }               
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
//           if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
//           {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
//           else
//           {lo_price.add((new Double(tokens[7]))); hi_price.add((new Double(tokens[8])));}
          lo_price.add((new Double(tokens[2]))); hi_price.add((new Double(tokens[3])));   
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          
          String[] hours = time.split("[:]+");
          cur_hour = (new Integer(hours[0])).intValue();
          (new Integer(hours[1])).intValue();
          
          trading_closed = false;
          
          //if currently not in a transaction and between the hours of midnight and start-hour, then no new 
          //positions will be opened
          
          
          if(asian_close)  //only closed if most recent transaction was closed artificially through SL or TP after end hours
          {
           if((in_transaction == 0 && out_transaction == 0) && cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;}// if(printall) System.out.println(cur_hour + " " + start_hour);}
          }
          else 
          {
           if(cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;} // if(printall) {System.out.println(cur_hour + " " + start_hour);}}          
          }
          
          if(inverse_hours) {trading_closed = !trading_closed;}
          
          its_closing_time = (weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1);
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && trade_times.contains(time)) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if(day_count == 0 || weekend.dayOfWeek().getAsText().equals("Sunday") && date_stamp.indexOf("18:00:00") != -1)
                {   
                   
                   if(printall)System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, print_filter);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall)System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                //if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                //if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal + ", trading_closed = " + trading_closed);
               if(friday_closing && its_closing_time)
               {
                 if(printall) System.out.println("\nIt's Friday at 5pm, time to close shop for week");
                 if(current_signal > 0 && in_transaction == 1) //in a long transaction
                 {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                  }
                 }
                 else if(current_signal < 0 && out_transaction == 1) //in a short transaction
                 {
                 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                 }
                 
               }
               else 
               { 
                //if(binary_rule && !trading_closed)
                if(binary_rule)
                {
                
                
//                  if(signal_profit)
//                  {
//                    //if(printall) System.out.println("Trading_closed = " + trading_closed);
//                    if((current_signal > 0 && in_transaction == 1) && (daily_price.get(daily_size-1) > last_price))
//                    {
//                    
//                    
//                      price_sold = daily_price.get(daily_size-1);
// 	             profit = price_sold - price_bought;
// 	             if(profit > 0) {succ_trades=succ_trades+1;}
// 	             total_trades=total_trades+1; 
// 	 
// 	             amount = getAmount(account.get(account.size()-1));
// 	             account.add(date_stamp + " " + (amount + profit));   
// 	             amount = getAmount(account.get(account.size()-1));
// 	            
// 	             made_trade = true;
// 	             dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                      log_ret = daily_price.get(daily_size-1).doubleValue();
// 	 
// 	             in_transaction = 0;
// 	           
// 	             if(printall){
//                      if(profit>0) System.out.println("Sold for a profit of " + profit);
//                      else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                       System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
//                      }
//                    
//                      if(!trading_closed)
//                      {price_bought = daily_price.get(daily_size-1);
// 	             //in_transaction = 1; 
// 	             out_transaction = 1;
// 	             last_price = daily_price.get(daily_size-1); 
// 	             if(printall) System.out.println("Entered long transaction at " + price_bought);
// 	             } 
//                    }
//                    else if((current_signal < 0 && out_transaction == 1) && (daily_price.get(daily_size-1) < last_price))
//                    {
//                    
//                     price_sold = daily_price.get(daily_size-1);
//                     profit = price_borrowed - price_sold;
// 
// 	            if(profit > 0) {succ_trades=succ_trades+1;}
// 	            total_trades=total_trades+1; 
// 	  
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
//                     
//                     made_trade = true;
//                     dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
//                     log_ret = daily_price.get(daily_size-1).doubleValue();
//                     
//                     
//                     out_transaction = 0;
//                     
//                     if(printall)
//                     {
//                     if(profit>0) System.out.println("Sold for a profit of " + profit);
//                     else if(profit<0) System.out.println("Sold for a loss of " + profit);
//                     }
//                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);                   
//                    
//                     if(!trading_closed)
//                     {
//                     price_borrowed = daily_price.get(daily_size-1);
//                     //out_transaction = 1;	 
//                     in_transaction = 1;
//                     if(printall)System.out.println("Entered short transaction at " + price_borrowed);
//                     }
//                    }
//                  
//                  }
                
                
                 made_trade = false;
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);                 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, buy
                  { 
                    //price_sold = daily_price.get(daily_size-1);
                    price_sold = ask.get(ask.size()-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(ask.get(ask.size()-1) - log_ret);
                    log_ret = ask.get(ask.size()-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    
                    if(printall){
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  
                  
                  if((long_buy && in_transaction == 0) && !trading_closed)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = bid.get(bid.size()-1);
                  
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = bid.get(bid.size()-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                    log_ret = bid.get(bid.size()-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	          if((short_sell && out_transaction == 0) && !trading_closed)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   if(printall) System.out.println("Entered short transaction at " + price_borrowed);
                  
	          }
                 }
                }
                
                if(signal_strength_rule && ((in_transaction == 0 && out_transaction == 0) && !trading_closed))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);
                  
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    price_sold = ask.get(ask.size()-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(ask.get(ask.size()-1) - log_ret);
                    log_ret = ask.get(ask.size()-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                  }
                  
                  
                  
                  
                  if(long_buy && in_transaction == 0)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           if(printall) System.out.println("Entered long transaction at " + price_bought);
	          } 
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = bid.get(bid.size()-1); 
                 
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = bid.get(bid.size()-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                    log_ret = bid.get(bid.size()-1).doubleValue();	            
	            
	           in_transaction = 0;
	           
                   if(printall)System.out.println("Bought for a profit of " + profit);
                   if(printall)System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           
	           
	          }
	         
	          if(short_sell && out_transaction == 0)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   if(printall)System.out.println("Entered short transaction at " + price_borrowed);
	          }
                 }              
                }
              
                if(!made_trade)
                {
                 account.add(date_stamp + " " + amount);
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {day_count=0;}
                //if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));             
                allsignal.add(date_stamp + " " + current_signal);
            }
           
           }
           if(trading_hours)// && !trading_closed)// && cur_min != 30)
           {
                
//              if(cur_pnl != 0 && date_stamp.indexOf("22:30:00") != -1) //close out position
//              {
//                if(in_transaction == 1)
//                {
//                
//                  //System.out.println("Ending trading for next startup");
// 
//                  
//                  price_sold = price.get(price.size()-1);
// 	         profit = price_sold - price_bought;                 
//                  
//                  if(profit > 0) {succ_trades=succ_trades+1;} 
//                  total_trades=total_trades+1;
//                  
//                  count = account.size()-1;
// 	         amount = getAmount(account.get(count));
// 
// 	         account.add(date_stamp + " " + (amount + profit));   
// 	         amount = getAmount(account.get(account.size()-1));
//                  
//                  
//                  dailyoutret.add(price.get(price.size()-1) - log_ret);
//                  log_ret = price.get(price.size()-1);                 
//  
//                  in_transaction = 0;
// 	         
// 	         //-- return stop loss to original setting --------
// 	         //System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
//                  stop_loss = global_stop_loss;
//                  red_zone = false;           
//                
//                }
//                else if(out_transaction == 1)
//                {
//                
//                  //System.out.println("Ending trading for next startup");
//                
//                
//                  price_sold = price.get(price.size()-1);
//                  profit = price_borrowed - price_sold;                 
//                  
//                  
// //                  price_sold = price_borrowed + stop_loss;
// //                  profit = -stop_loss;
//                  
// 	         if(profit > 0) {succ_trades=succ_trades+1;}
// 	         total_trades=total_trades+1; 
// 	  
//                  count = account.size()-1;
// 	         amount = getAmount(account.get(count));
// 	           
// 	         //account.add(date_stamp + " " + (amount - stop_loss));  
// 	         account.add(date_stamp + " " + (amount + profit));  
// 	         amount = getAmount(account.get(account.size()-1));
// 
//                  dailyoutret.add(price.get(price.size()-1) - log_ret);
//                  log_ret = price.get(price.size()-1);    	         
// 	         
// 	         
// 	         out_transaction = 0;
// 
// 	         //-- return stop loss to original setting --------
//                  stop_loss = global_stop_loss;
//                  red_zone = false;
//                
//                }
//              }

                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && bid.get(bid.size()-1) > last_price) //check if new high price
               {last_price = bid.get(bid.size()-1);}
                
              
               //cur_pnl = price.get(price.size()-1) - last_price;    
               cur_pnl = bid.get(bid.size()-1) - last_price;
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 //price_sold = price.get(price.size()-1);
	         price_sold = bid.get(bid.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 
                 dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                 log_ret = bid.get(bid.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = bid.get(bid.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && ask.get(ask.size()-1) < last_price) //check if new high price
               {last_price = ask.get(ask.size()-1);}
                
               //cur_pnl =  last_price - price.get(price.size()-1);               
               cur_pnl =  last_price - ask.get(ask.size()-1); 
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                if(printall) System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size() - 1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  	         
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));	         
	         
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	         
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall)System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 if(printall)System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 if(printall)System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 if(printall)System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               if(printall)System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall) 
             {
               if(current_signal > 0) {System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
               else {System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
             }
           
             allsignal.add(date_stamp + " " + current_signal);
           }      
                    
          }             
        }                        
          
          
           

      double mean_ntrades = 0;
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      int n_trades = 0;
      ArrayList<Integer> n_trades_day = new ArrayList<Integer>();
      double pret;
      ArrayList<Double> perf_rets = new ArrayList<Double>();
      ArrayList<String> perf_datestimes = new ArrayList<String>();
      
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i)) - getAmount(account.get(i-1));
      
        if(jpy) 
        {dreturns[i] = dreturns[i]*.01;}
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0])) //first date entry
        {
         
         

         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
           n_trades_day.add(n_trades);
         }
         
         
         perf_dates.add(dates[0]);
         if(dreturns[i] != 0)
         {
          perf_datestimes.add(dates[0] + " " + dates[1]);
          perf_rets.add(dreturns[i]);
         }
         pnl = dreturns[i];
         if(dreturns[i] != 0) {n_trades = 1;} 
         else {n_trades = 0;}
        }
        else //already contains the date, so add on pnl
        {
         pnl = pnl + dreturns[i]; 
         //System.out.println(dreturns[i]);
         if(dreturns[i] != 0) 
         {
           n_trades++;
           perf_datestimes.add(dates[0] + " " + dates[1]);
           perf_rets.add(dreturns[i]);                   
         }// System.out.println(n_trades);}
        }
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      n_trades_day.add(n_trades);
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       if(printall) System.out.println(perf_dates.get(i) + " " + perf_returns.get(i) + " " + n_trades_day.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
       date_returns.add(perf_dates.get(i) + " " + perf_returns.get(i));
       mean_ntrades = mean_ntrades + n_trades_day.get(i);
      }
      mean_ntrades = mean_ntrades/n_trades_day.size();
      
      
       
      
      
      dates_price = new String[n_obs];
      
      for(i=1;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        
        if(perf_datestimes.contains(latestDates.get(latestDates.size() - n_obs + i)))
        {
          pret = perf_rets.get(perf_datestimes.indexOf(latestDates.get(latestDates.size() - n_obs + i)));
          System.out.println(latestDates.get(latestDates.size() - n_obs + i) + " " + pret);
        }
        else
        {pret = 0;}
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2] + " " + pret);
//        }
      }  
      dates_price[0] = dates_price[1];
  
  
  
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      double risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      double reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff + ", avg_n_trades = " + mean_ntrades);
 
      out.close(); dailyout.close(); perform.close(); b0_coeff.close();
      
      }
      catch(ParseException fe){System.out.println("ParseException");}
      catch(NullPointerException fe){System.out.println("Null pointer");}
      catch(IllegalArgumentException fe){System.out.println("IllegalArgument");}
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
      n_files++;
     
      return computed; 
      
  }            
  
  
    public boolean startStrategyDailyIntradayOneHourBidAsk_MR(StrategyParameters strat)
    {
      
      int j,i,l,N,file_count;
      double sum = 0;
      Double D; 
      String ddelims = "[-]";
      boolean computed = false;
      boolean print_filter = false;
      boolean made_trade = true;
      String date_stamp,strline; 
      int daily_size;
      double profit,price_borrowed,price_sold,price_bought;
      double current_price,prev_price;
      double last_price,cur_pnl,stop_loss,lo_pnl,hi_pnl;
      double log_ret = 0;
      signal = new double[trade_obs];
      xt = new double[trade_obs];
      lag_signals = new double[trade_obs];
      prix = new double[trade_obs];
      lo_prix = new double[trade_obs];
      hi_prix = new double[trade_obs];
      total_succ = 0; total = 0;
      log_price = 0;
      N = n_obs; avg_vol = 0.0;
      b_avg = new double[L*n_rep];
      count=0; 
      trade_succ_ratio = 0; 
      double amount = 0;
      double prev_signal;
      reg_trading_hours = false;
      String[] intdates; 
      //make sure arraylists empty
      ArrayList<String> perf_dates = new ArrayList<String>();
      ArrayList<Double> perf_returns = new ArrayList<Double>();
      double pnl; 
      String[] dates; 
      boolean inverse_hours = false;
      String time; 
      ArrayList<String> account = new ArrayList<String>();
      ArrayList<String> latestDates = new ArrayList<String>();
      last_trades = new ArrayList<Integer>();
      final_trades = new ArrayList<Double>();
      dailyoutret = new ArrayList<Double>();
      maxIntValue = new ArrayList<Double>();
      avg_volatility = new ArrayList<Double>();
      close_series = new ArrayList<Double>();
      highlow_series = new ArrayList<Double>();    
      exp_series_1 = new ArrayList<Double>();  
      exp_series_2 = new ArrayList<Double>();
      price = new ArrayList<Double>();    
      lo_price = new ArrayList<Double>();    
      hi_price = new ArrayList<Double>();    
      mid = new ArrayList<Double>();
      bid = new ArrayList<Double>();
      ask = new ArrayList<Double>();
      dates_series = new ArrayList<String>();       
      dailyReport = new ArrayList<String>();
      b0_trend = new ArrayList<Double>();
      vol_0 = new ArrayList<Double>();
      vol_1 = new ArrayList<Double>();
      sub_returns = new ArrayList<Double>();
      trade_days = new ArrayList<String>();
      returns = new ArrayList<Double>();
      longreturns = new ArrayList<Double>();
      shortreturns = new ArrayList<Double>();
      dropdowns = new ArrayList<Double>();
      success = new ArrayList<Double>();
      dates_low_high = new ArrayList<String>();       
      crits = new ArrayList<String>();
      svm = new ArrayList<String>();
      filters = new ArrayList<Filter>();
      date_returns = new ArrayList<String>();
      
      live_series = new ArrayList<Double>(); //the data to be applied out of sample
      ib_data_hash = new ibHash();
     
      fridayROI = 0; fridayROI_pos = 0; fridays = 0;
      int end_hour;
      lookback_returns = new ArrayList<Double>();
      num_pos_returns=0;
      deg_0 = new ArrayList<Double>();
      deg_1 = new ArrayList<Double>();
      crit_0 = new ArrayList<Double>();
      crit_1 = new ArrayList<Double>();
      full_returns_array = new ArrayList<double[]>();
      morning_returns = new ArrayList<double[]>();
      
      boolean waiting_meanrev_down = false;
      boolean waiting_meanrev_up = false;
      double mean_rev_amnt = .0005;
      
      morning_buy = true;        //enter transaction at morning open
      morning_optimize = false;   //optimize in the morning trading hours
      num_full_positive_returns = 0;
      //--- Now get historical interp values ------
      //uploadInterpParams("max_int.dat"); 
      //-------------------------------------------
      forex24 = true;
      ret_dist = new double[trade_obs];
      pos_ret_dist = new int[trade_obs];
      neg_ret_dist = new int[trade_obs];
      neg_trades_started = new int[trade_obs];
      pos_trades_started = new int[trade_obs];
      neg_trades_started_mean = new double[trade_obs];
      pos_trades_started_mean = new double[trade_obs];      
      diff_account = new double[trade_obs];
      pos_ret_mean_time = new double[trade_obs];
      neg_ret_mean_time = new double[trade_obs];
      
      
      mdfaTrades = new ArrayList<MDFATrade>();
      fmt = DateTimeFormat.forPattern("y-MM-dd HH:mm:ss");
      formatter = new DecimalFormat("#0.000000");   
      formatter3 = new DecimalFormat("#0.00000");   
      formatter2 = new DecimalFormat("#0.00");
      histo_stat = new int[100];
      interp_vals = new ArrayList<Double>();
      max_ranks = new ArrayList<Double>();
      profit_baby = 0;
      //setForecastDFAParameters();
      bad_starts = 0; 
      n_out_samp = 0;
      
      //take_profit = true;
      //take_profit_thresh = .0020;
      current_signal = 0;
      prev_price = 0;
      cur_pnl = 0;
      stop_loss = stop_loss_thresh;
      out_transaction = 0; 
      in_transaction = 0;
      red_zone = false; 
      global_stop_loss = stop_loss_thresh;
      profitable_stop = .0005;
      count = 0;
      short_sell = true; long_buy = true;
      day_count = 0;
      ArrayList<String> trade_times = new ArrayList<String>();
      ArrayList<String> fxcm_dates = new ArrayList<String>();
      
      
      lo_pnl = 0; hi_pnl = 0;
      price_borrowed = 0; price_sold = 0; price_bought = 0; last_price = 0;
      
      binary_rule = true; 
      signal_strength_rule = true;
      downtick_strategy = false;
      signal_profit = false;
      friday_closing = false;
      //asian_close = true;
      
      
      //----- Fill up with times here -------------------------
      for(i=0;i<10;i++) {trade_times.add("0"+i+":00:00");}
      for(i=10;i<24;i++) {trade_times.add(i+":00:00");}
//       trade_times.add("00:00:00");
//       trade_times.add("06:00:00");
//       trade_times.add("12:00:00");
//       trade_times.add("18:00:00");     
      
      String[] hourToks = startingTime.split("[:]+");
      start_hour = (new Integer(hourToks[0])).intValue();
      
      hourToks = endingTime.split("[:]+");
      end_hour = (new Integer(hourToks[0])).intValue();
      
      inverse_hours = false;
      if(start_hour > end_hour)  //must switch hours here
      {
        int temphour = end_hour;
        
        inverse_hours = true;
        end_hour = start_hour;
        start_hour = temphour;
      }
      
      
      if(ib_data && ib_data_file != null )
      {
        try{
        
         fin = new FileInputStream(ib_data_file);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));

         while((strline = br.readLine()) != null)
         {
           String[] sp = strline.split("[,]+");
           ib_data_hash.put(sp[0], new String(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]));
           //System.out.println(sp[1] + " " + sp[2] + " " + sp[3] + " " + sp[4] + " " + sp[5]  + " " + sp[6]);
         }
        }
        catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
        catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}          
      }
      
      

      jpy = false;
      try{
       
        PrintWriter b0_coeff = new PrintWriter(new FileWriter("b0_coeff.dat"));
        PrintWriter perform = new PrintWriter(new FileWriter("intraday_performance_"+n_files+".dat"));
        PrintWriter dailyout = new PrintWriter(new FileWriter("daily_nasdaq.dat"));
        PrintWriter out = new PrintWriter(new FileWriter("strategy_results.dat"));
       
        
        
        for(file_count=0;file_count<1;file_count++)
        {
         
         if(dataFiles[file_count].indexOf("JPY") != -1 && dataFiles[file_count].indexOf("NOKJPY") == -1)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          jpy = true; 
          stop_loss_thresh = stop_loss_thresh*100;
          take_profit_thresh = take_profit_thresh*100;
          global_stop_loss = global_stop_loss*100;
          stop_loss = stop_loss_thresh;
         }
         else if(futures_data)
         {
          //change_time_zone = true; System.out.println("Changed time zone to Tokyo");
          
          stop_loss_thresh = stop_loss_thresh*10000;
          stop_loss = stop_loss_thresh;
          take_profit_thresh = take_profit_thresh*10000;
          global_stop_loss = global_stop_loss*10000;
         }         
         
         ArrayList<String> dayData = new ArrayList<String>();
         setTimeStandards(new File(dataFiles[file_count]));
         System.out.println("opening " + dataFiles[file_count]);
         fin = new FileInputStream(dataFiles[file_count]);
         din = new DataInputStream(fin);
         br = new BufferedReader(new InputStreamReader(din));     
         lookback_ready = false;
         spread = new PrintWriter(new FileWriter("spread_" + dataFiles[file_count] + ".dat"));
         //if(print_debug)System.out.println("Entering loop...");
         trading_hours = false; computed = false;
         
         
         //---- add the latest FXCM data -
         while((strline = br.readLine()) != null)
         {
           dayData.add(strline);
           
           //add the date
           tokens = strline.split("[,]+");
           fxcm_dates.add(tokens[0]);
           
         }
         
         //----- now add the lates IB data ------------
         System.out.println("opening " + dataFiles[file_count] + " IB data");
         String[] fxpair = dataFiles[file_count].split("[.]+");
         
         if((new File(fxpair[0] + ".IB.dat")).exists())
         {
          fin = new FileInputStream(fxpair[0] + ".IB.dat");
          din = new DataInputStream(fin);
          br = new BufferedReader(new InputStreamReader(din)); 
         
          String lastFXCMdate = fxcm_dates.get(fxcm_dates.size()-1);
         
          SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-DD HH:mm:ss");
          Date date1 = sdf.parse(lastFXCMdate);

         
         
          while((strline = br.readLine()) != null)
          {
           tokens = strline.split("[,]+");
           
           Date date2 = sdf.parse(tokens[0]);
           
           if(!fxcm_dates.contains(tokens[0]) && date1.before(date2))
           {
             System.out.println("New time and observation " + strline);
             dayData.add(strline);
           }
          }
         }
         
         for(int ts = 0; ts < dayData.size(); ts++)
         {
          strline = dayData.get(ts);
          
          //System.out.println(strline);
          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          if(n_toks >= 6) {bid_ask_data = true;}
          else {bid_ask_data = false;}
           
           
  
           
           
           
          date_stamp = tokens[0];             
          date_tokens = date_stamp.split(date_delims); 
          intdates = date_tokens[0].split(ddelims);     
          DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 14, 0); 
          time = date_tokens[1];

          
          //insampStart is the time we collect daily data 
          
          
          //if(date_stamp.indexOf(insampStart) != -1)
          if(trade_times.contains(time))
          {
            
            //get bid/mid/ask data
            if(ib_data && ib_data_hash.containsKey(tokens[0]))
            {
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
            }          
            
            daily_price.add(log(new Double(tokens[1])));
            current_price = (new Double(tokens[1])).doubleValue();
            
            if(daily_price.size() == 1) {daily_returns.add(new Double(0.0)); prev_price = current_price;}
            else
            {
             daily_returns.add(log(current_price) - log(prev_price));
             prev_price = current_price;
            }
            
            daily_dates.add(date_stamp);
            
          }
          
          
          print_filter = false;
          latestDates.add(date_stamp);
          D = new Double(tokens[4]); close_series.add(D);
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {
              
              String[] hashed = ib_data_hash.get(tokens[0]).split("[ ]+");
              //System.out.println("Contains " + tokens[0] + ", lengths = " + hashed.length + ", " + tokens.length);
              for(i=1;i<hashed.length;i++) {tokens[i] = hashed[i];}// System.out.print(tokens[i] + " ");} 
              bid.add(log(new Double(tokens[2]))); ask.add(log(new Double(tokens[3]))); mid.add(log(new Double(tokens[1])));
          }
          else
          {
           bid.add(log(new Double(tokens[2]))); ask.add(log(new Double(tokens[3]))); mid.add(log(new Double(tokens[1])));
          }
            
 
          D = new Double(tokens[1]); 
          price.add(log(D)); 
              
          D = new Double(tokens[4]); 
          if(ib_data && ib_data_hash.containsKey(tokens[0]))
          {live_series.add(log(mid.get(mid.size()-1)) - log(mid.get(mid.size()-2)));}
          else
          {live_series.add(D);}
             
//           if(ib_data && ib_data_hash.containsKey(tokens[0])) //use as is
//           {lo_price.add(new Double(tokens[7])); hi_price.add(new Double(tokens[8]));}
//           else
//           {lo_price.add((new Double(tokens[7]))); hi_price.add((new Double(tokens[8])));}
          lo_price.add((new Double(tokens[2]))); hi_price.add((new Double(tokens[3])));   
          
          //---- start the account ------
          if(account.size() == 0) {account.add(date_stamp + " " + 0); dailyoutret.add(0.0);}  
          
          
          String[] hours = time.split("[:]+");
          cur_hour = (new Integer(hours[0])).intValue();
          (new Integer(hours[1])).intValue();
          
          trading_closed = false;
          
          //if currently not in a transaction and between the hours of midnight and start-hour, then no new 
          //positions will be opened
          
          
          if(asian_close)  //only closed if most recent transaction was closed artificially through SL or TP after end hours
          {
           if((in_transaction == 0 && out_transaction == 0) && cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;}// if(printall) System.out.println(cur_hour + " " + start_hour);}
          }
          else 
          {
           if(cur_hour >= start_hour && cur_hour <= end_hour)
           {trading_closed = true;} //if(printall) {System.out.println(cur_hour + " " + start_hour);}}          
          }
          
          if(inverse_hours) {trading_closed = !trading_closed;}
          
          its_closing_time = (weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1);
          
          made_trade = false;
          if(daily_returns.size() >= n_obs && trade_times.contains(time)) //a new day begineth
          {
 
                computed = true;
                trading_hours = true;
                
                tseries = new double[n_rep*n_obs];
                
                
                for(i=0;i<n_obs;i++)
                {
                     tseries[n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     tseries[n_obs + n_obs-1-i] = daily_returns.get(daily_returns.size() - 1 - i);
                     
                     if(n_rep > 2 && exp_series_1.size() > 0)
                     {tseries[n_obs*2 + n_obs-1-i] = exp_series_1.get(exp_series_1.size() - 1 - i);}
                     if(n_rep > 3 && exp_series_2.size() > 0)
                     {tseries[n_obs*3 + n_obs-1-i] = exp_series_2.get(exp_series_2.size() - 1 - i);}
                }
               
                mdfa.set_tseries(tseries,n_obs,n_rep);
              
                //if(day_count == 0)  //recompute filter coefficients
                if(day_count == 0 || weekend.dayOfWeek().getAsText().equals("Sunday") && date_stamp.indexOf("18:00:00") != -1)
                {   
                   
                   if(printall)System.out.println("Recomputing filter...");
                   mdfa.computeFilterGeneral(true, print_filter);        
                   b_coeffs = new double[(n_rep-1)*L]; //System.out.println(b_coeffs.length + " " + L + n_rep); 
                   for(l=0;l<L;l++)
                   {
                     
                     for(i=0;i<n_rep-1;i++)
                     {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}// System.out.println(b_coeffs[l]);}               
                     //if(date_stamp.indexOf("2013-12-17") != -1) {System.out.println(b_coeffs[l]);}
                   }
                   if(printall)System.out.println(date_stamp + " b_coeffs = " + b_coeffs[0] + " " + b_coeffs[1] + " " + b_coeffs[2]);
                   
                   b_copy = new double[mdfa.b.length];
                   System.arraycopy(mdfa.b, 0, b_copy, 0, b_copy.length);
                   b0_coeff.println(b_coeffs[0]); // + ", " + b_coeffs[L] + ", " + b_coeffs[2*L]);
                
                }

                sum = 0.0;
                for(j=1;j<n_rep;j++)
                {
                    for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + n_obs-1-l];}
                }
                 
                  
                prev_signal = current_signal;
                current_signal = sum; 
                if(sig_inverse) {current_signal = -current_signal;}   
                
                
                //----final signal ---
                daily_signal.add(current_signal);             
                daily_size = daily_price.size();
                dailyReport.add("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
                //if(printall) System.out.println("New day "+date_stamp + ", " + formatter3.format(daily_price.get(daily_size-1)) + ", " + formatter.format(tseries[n_obs-1]) + ", " + current_signal);
           
                //--compute binary trading rule ---
                
                //if(printall) System.out.println("Current signal = " + current_signal + " Prev signal = " + prev_signal + ", trading_closed = " + trading_closed);
               if(friday_closing && its_closing_time)
               {
                 if(printall) System.out.println("\nIt's Friday at 5pm, time to close shop for week");
                 if(current_signal > 0 && in_transaction == 1) //in a long transaction
                 {
	           price_sold = daily_price.get(daily_size-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);
                  }
                 }
                 else if(current_signal < 0 && out_transaction == 1) //in a short transaction
                 {
                 
                    price_sold = daily_price.get(daily_size-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));

	            made_trade = true;
                    dailyoutret.add(daily_price.get(daily_size-1) - log_ret);
                    log_ret = daily_price.get(daily_size-1).doubleValue();	            
	            
                    out_transaction = 0;
                  

                     if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                 }
                 
               }
               else 
               { 
                //if(binary_rule && !trading_closed)
                if(binary_rule)
                {
                
    
                
                 made_trade = false;
                 if(current_signal > 0 && prev_signal <= 0) //new point positive, we see momentum, buy
                 {
                  
                  waiting_meanrev_down = false;
                  waiting_meanrev_up = false;
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);                 
                                   
                  if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
                  { 
                    //price_sold = daily_price.get(daily_size-1);
                    price_sold = ask.get(ask.size()-1);
                    profit = price_borrowed - price_sold;

	            if(profit > 0) {succ_trades=succ_trades+1;}
	            total_trades=total_trades+1; 
	  
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
                    
                    made_trade = true;
                    dailyoutret.add(ask.get(ask.size()-1) - log_ret);
                    log_ret = ask.get(ask.size()-1).doubleValue();
                    
                    
                    out_transaction = 0;
                    
                    if(printall){
                    if(profit>0) System.out.println("Sold for a profit of " + profit);
                    else if(profit<0) System.out.println("Sold for a loss of " + profit);
                    
                    System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
                    }
                  }
                  
                  
                  waiting_meanrev_down = true;
                  
//                   if((long_buy && in_transaction == 0) && !trading_closed)
// 	          {
//                    price_bought = ask.get(ask.size()-1);
// 	           in_transaction = 1; 
// 	           
// 	           if(printall) System.out.println("Entered long transaction at " + price_bought);
// 	          } 
                
                 }
                 else if(current_signal < 0 && prev_signal >= 0) //if in transaction and signal goes below, sell
                 {
                  //last_price = daily_price.get(daily_size-1); 
       
                  waiting_meanrev_down = false;
                  waiting_meanrev_up = false;       
                  
                  last_price = bid.get(bid.size()-1);
                  
                  if(long_buy && in_transaction == 1)
                  {
	           price_sold = bid.get(bid.size()-1);
	           profit = price_sold - price_bought;
	           if(profit > 0) {succ_trades=succ_trades+1;}
	           total_trades=total_trades+1; 
	 
	            amount = getAmount(account.get(account.size()-1));
	            account.add(date_stamp + " " + (amount + profit));   
	            amount = getAmount(account.get(account.size()-1));
	            
	            made_trade = true;
	            dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                    log_ret = bid.get(bid.size()-1).doubleValue();
	 
	           in_transaction = 0;
	           
	           if(printall){
                   if(profit>0) System.out.println("Sold for a profit of " + profit);
                   else if(profit<0) System.out.println("Sold for a loss of " + profit);
                   System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
	           }
	           
	          }
	         
	         
	          waiting_meanrev_up = true;
// 	          if((short_sell && out_transaction == 0) && !trading_closed)
// 	          {
// 	           price_borrowed = bid.get(bid.size()-1);
//                    out_transaction = 1;	 
//                   
//                    if(printall) System.out.println("Entered short transaction at " + price_borrowed);
//                   
// 	          }
                 }
                }
                
                if(signal_strength_rule && ((in_transaction == 0 && out_transaction == 0) && !trading_closed))
                {
                 
                 if(current_signal > 0) //new point positive, we see momentum, buy
                 {
                  //last_price = daily_price.get(daily_size-1); 
                  last_price = ask.get(ask.size()-1);
                  waiting_meanrev_down = true;
//                   if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
//                   { 
//                     price_sold = ask.get(ask.size()-1);
//                     profit = price_borrowed - price_sold;
// 
// 	            if(profit > 0) {succ_trades=succ_trades+1;}
// 	            total_trades=total_trades+1; 
// 	  
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
// 
// 	            made_trade = true;
//                     dailyoutret.add(ask.get(ask.size()-1) - log_ret);
//                     log_ret = ask.get(ask.size()-1).doubleValue();	            
// 	            
//                     out_transaction = 0;
//                   
// 
//                      if(printall) System.out.println("profit, price_borrowed, price_sold = " + profit + ", " + price_borrowed + ", " + price_sold);
//                   }
                  
                  
                  
                  
// //                   if(long_buy && in_transaction == 0)
// // 	          {
// //                    price_bought = ask.get(ask.size()-1);
// // 	           in_transaction = 1; 
// // 	           
// // 	           if(printall) System.out.println("Entered long transaction at " + price_bought);
// // 	          } 
                
                
                 }
                 else if(current_signal < 0) //if in transaction and signal goes below, sell
                 {
                 
                  last_price = bid.get(bid.size()-1); 
                  waiting_meanrev_up = true;
//                   if(long_buy && in_transaction == 1)
//                   {
// 	           price_sold = bid.get(bid.size()-1);
// 	           profit = price_sold - price_bought;
// 	           if(profit > 0) {succ_trades=succ_trades+1;}
// 	           total_trades=total_trades+1; 
// 	 
// 	            amount = getAmount(account.get(account.size()-1));
// 	            account.add(date_stamp + " " + (amount + profit));   
// 	            amount = getAmount(account.get(account.size()-1));
// 
// 	            made_trade = true;
//                     dailyoutret.add(bid.get(bid.size()-1) - log_ret);
//                     log_ret = bid.get(bid.size()-1).doubleValue();	            
// 	            
// 	           in_transaction = 0;
// 	           
//                    if(printall)System.out.println("Bought for a profit of " + profit);
//                    if(printall)System.out.println("profit, price_bought, price_sold = " + profit + ", " + price_bought + ", " + price_sold);	           
// 	           
// 	           
// 	          }
	         
// 	          if(short_sell && out_transaction == 0)
// 	          {
// 	           price_borrowed = bid.get(bid.size()-1);
//                    out_transaction = 1;	 
//                   
//                    if(printall)System.out.println("Entered short transaction at " + price_borrowed);
// 	          }
                 }              
                }
              
                if(!made_trade)
                {
                 account.add(date_stamp + " " + amount);
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);   
                }
              
     
                
                day_count++;
              
                //if(recompute_day == day_count) {day_count=0;}
                if(printall) System.out.println(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));             
                allsignal.add(date_stamp + " " + current_signal);
            }
           
           }
           if(trading_hours)// && !trading_closed)// && cur_min != 30)
           {
                
             if(cur_pnl != 0 && date_stamp.indexOf("22:30:00") != -1) //close out position
             {
               if(in_transaction == 1)
               {
               
                 //System.out.println("Ending trading for next startup");

                 
                 price_sold = price.get(price.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));

	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 
                 
                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);                 
 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         //System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;           
               
               }
               else if(out_transaction == 1)
               {
               
                 //System.out.println("Ending trading for next startup");
               
               
                 price_sold = price.get(price.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	           
	         //account.add(date_stamp + " " + (amount - stop_loss));  
	         account.add(date_stamp + " " + (amount + profit));  
	         amount = getAmount(account.get(account.size()-1));

                 dailyoutret.add(price.get(price.size()-1) - log_ret);
                 log_ret = price.get(price.size()-1);    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 waiting_meanrev_up = false;
                 waiting_meanrev_down = false;
               
               }
             }

             
             if(waiting_meanrev_up) //waiting for tick price to go up in order to sell the bid
             {
             
               if((bid.get(bid.size()-1) - last_price) > mean_rev_amnt)
               {
               
               	  if((short_sell && out_transaction == 0)  && !trading_closed)
	          {
	           price_borrowed = bid.get(bid.size()-1);
                   out_transaction = 1;	 
                  
                   last_price = price_borrowed;
                   if(printall)System.out.println("Waiting over Entered short transaction at " + price_borrowed + " after " + (bid.get(bid.size()-1) - last_price) + " reversion");
	          }
                  waiting_meanrev_up = false;
               }
             }
             else if(waiting_meanrev_down)
             {
               if((last_price - ask.get(ask.size() - 1)) > mean_rev_amnt)
               {
                  if((long_buy && in_transaction == 0)  && !trading_closed)
	          {
                   price_bought = ask.get(ask.size()-1);
	           in_transaction = 1; 
	           
	           last_price = price_bought;
	           if(printall) System.out.println("Waiting over, Entered long transaction at " + price_bought + " after " + (last_price - ask.get(ask.size() - 1)) + " reversion");
	          } 
                  waiting_meanrev_down = false;
               }
             }
             
                          
             if(in_transaction == 1) //in a long transaction 
             {
                
               if(red_zone && bid.get(bid.size()-1) > last_price) //check if new high price
               {last_price = bid.get(bid.size()-1);}
                
              
               //cur_pnl = price.get(price.size()-1) - last_price;    
               cur_pnl = bid.get(bid.size()-1) - last_price;
               lo_pnl = lo_price.get(lo_price.size()-1) - last_price; 
               hi_pnl = hi_price.get(hi_price.size()-1) - last_price;
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but lowest price in bar was " + lo_price.get(lo_price.size()-1));
                 //--------------sell---------- 
               
                 //price_sold = price.get(price.size()-1);
                 
                 //price_sold = price.get(price.size()-1);
	         price_sold = bid.get(bid.size()-1);
	         profit = price_sold - price_bought;                 
                 
                 
                 
//                  price_sold = price_bought - stop_loss;
//                  profit = -stop_loss;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
                 //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 
                 dailyoutret.add(bid.get(bid.size()-1) - log_ret);
                 log_ret = bid.get(bid.size()-1);                 
                 
                 
                 
                 
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;
                 red_zone = false;
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = bid.get(bid.size()-1);
                 profit = price_sold - price_bought;
                 
                 if(profit > 0) {succ_trades=succ_trades+1;} 
                 total_trades=total_trades+1;
                 
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
                 in_transaction = 0;
	         
	         //-- return stop loss to original setting --------
	         if(printall)System.out.println("Sold at " + date_stamp + " for a profit/loss of " + profit);
                 stop_loss = global_stop_loss;                 
                 
                 
                 
                 
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;
               }
             }
             else if(out_transaction == 1)
             {
               
               if(red_zone && ask.get(ask.size()-1) < last_price) //check if new high price
               {last_price = ask.get(ask.size()-1);}
                
               //cur_pnl =  last_price - price.get(price.size()-1);               
               cur_pnl =  last_price - ask.get(ask.size()-1); 
               lo_pnl = last_price - hi_price.get(hi_price.size()-1);
               hi_pnl = last_price - lo_price.get(lo_price.size()-1);
               
               if(cur_pnl < -stop_loss)
               {
                 if(printall)System.out.println("Stop-loss Activated since lo_pnl = " + cur_pnl + " < -" + stop_loss);
                 //System.out.println("Closing price at bar was " + price.get(price.size()-1) + " but highest price in bar was " + hi_price.get(hi_price.size()-1));
                 //--------------sell---------- 
                 
                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size()-1);
                 profit = price_borrowed - price_sold;                 
                 
                 
//                  price_sold = price_borrowed + stop_loss;
//                  profit = -stop_loss;
                 
	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}

                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
	         out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                if(printall) System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);
                 
               } 
               else if(cur_pnl >= take_profit_thresh)
               {               
               
                 if(printall)System.out.println("Profit zone activated since cur_pnl = " + cur_pnl + " > " + take_profit_thresh);
//                  stop_loss = profitable_stop;
//                  last_price = price.get(price.size()-1);
//                  red_zone = true;

                 //price_sold = price.get(price.size()-1);
                 price_sold = ask.get(ask.size() - 1);
                 profit = price_borrowed - price_sold;

	         if(profit > 0) {succ_trades=succ_trades+1;}
	         total_trades=total_trades+1; 
	  	         
                 count = account.size()-1;
	         amount = getAmount(account.get(count));
	         account.add(date_stamp + " " + (amount + profit));   
	         amount = getAmount(account.get(account.size()-1));	         
	         
	         //if(weekend.dayOfWeek().getAsText().equals("Sunday")) {sunday.add(date_stamp + " " + profit);}
	         
                 dailyoutret.add(price_sold - log_ret);
                 log_ret = price_sold;    	         
	         
	         
                 out_transaction = 0;

	         //-- return stop loss to original setting --------
                 stop_loss = global_stop_loss;
                 red_zone = false;
                 
                 if(printall)System.out.println("Bought at " + date_stamp + " for a profit/loss of " + profit);


               }
             }  
             else if(downtick_strategy)
             {
             
                //strategy here is to buy/sell according to signal iff downtick has occurred
             
               if(current_signal > 0 && (price_sold > price.get(price.size()-1)))
               {
                
                 //let's buy some more 
                 if(printall)System.out.println("Buying at " + date_stamp + " since last price sold = " + price_sold + " > " + price.get(price.size()-1));
                 price_bought = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
	         in_transaction = 1; 
            
               }
               else if(current_signal < 0  && (price_sold < price.get(price.size()-1)))
               {
               
                 //let's short some more
                 if(printall)System.out.println("Shorting at " + date_stamp + " since last price bought back at = " + price_sold + " < " + price.get(price.size()-1));
               
                 price_borrowed = price.get(price.size()-1);
                 last_price = price.get(price.size()-1);
                 out_transaction = 1;	 
                  
                 if(printall)System.out.println("Entered short transaction at " + price_borrowed);
               }
               cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;
             }
             else
             {cur_pnl = 0; lo_pnl = 0; hi_pnl = 0;}
             
             
             
             if(weekend.dayOfWeek().getAsText().equals("Friday") && date_stamp.indexOf("17:00:00") != -1)
             {
               if(printall)System.out.println("End of week");
             
             }
             
             
             
             
             dailyReport.add(""+date_stamp + ", " + formatter3.format(price.get(price.size()-1)) 
              + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));     
             
             if(printall) 
             {
               if(current_signal > 0) {System.out.println(""+date_stamp + ", " + formatter3.format(bid.get(bid.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
               else {System.out.println(""+date_stamp + ", " + formatter3.format(ask.get(ask.size()-1)) + ", " + current_signal + ", " + formatter3.format(cur_pnl) + ", " + formatter3.format(lo_pnl) + ", " + formatter3.format(hi_pnl));} 
             }
           
             allsignal.add(date_stamp + " " + current_signal);
           }      
                    
          }             
        }                        
          
          
           

      double mean_ntrades = 0;
      computed = true;
      dailyoutret.set(1,0.0);
      double[] dreturns = new double[account.size()]; 
      dreturns[0] = 0;
      
      double mean=0; double sd=0;
      n_neg_ret = 0; 
      n_pos_ret = 0;
      neg_ret_mean = 0; 
      pos_ret_mean = 0;
      pnl = 0;
      
      int n_trades = 0;
      ArrayList<Integer> n_trades_day = new ArrayList<Integer>();
      double pret;
      ArrayList<Double> perf_rets = new ArrayList<Double>();
      ArrayList<String> perf_datestimes = new ArrayList<String>();
      
      
      for(i=1;i<account.size();i++)
      {
        out.println(account.get(i));
        dailyout.println(account.get(i) + " " + dailyoutret.get(i)); 
        //System.out.println(account.get(i));
        dreturns[i] = getAmount(account.get(i)) - getAmount(account.get(i-1));
      
        if(jpy) 
        {dreturns[i] = dreturns[i]*.01;}
      
      
        dates = account.get(i).split("[ ]+");
        if(!perf_dates.contains(dates[0])) //first date entry
        {
         
         

         
         if(perf_dates.size() != 0)
         {
           perf_returns.add(pnl);
           n_trades_day.add(n_trades);
         }
         
         
         perf_dates.add(dates[0]);
         if(dreturns[i] != 0)
         {
          perf_datestimes.add(dates[0] + " " + dates[1]);
          perf_rets.add(dreturns[i]);
         }
         pnl = dreturns[i];
         if(dreturns[i] != 0) {n_trades = 1;} 
         else {n_trades = 0;}
        }
        else //already contains the date, so add on pnl
        {
         pnl = pnl + dreturns[i]; 
         //System.out.println(dreturns[i]);
         if(dreturns[i] != 0) 
         {
           n_trades++;
           perf_datestimes.add(dates[0] + " " + dates[1]);
           perf_rets.add(dreturns[i]);                   
         }// System.out.println(n_trades);}
        }
       

      
        if(dreturns[i] > 0)
        {
          n_pos_ret++; 
          pos_ret_mean = pos_ret_mean + dreturns[i];
          mean = mean + dreturns[i];
        }
        else if(dreturns[i] < 0)
        {
          n_neg_ret++; 
          neg_ret_mean = neg_ret_mean - dreturns[i];
          mean = mean + dreturns[i];
        }
        
      }
      perf_returns.add(pnl); 
      n_trades_day.add(n_trades);
      
      
      for(i=0;i<perf_dates.size();i++)
      {
       if(printall) System.out.println(perf_dates.get(i) + " " + perf_returns.get(i) + " " + n_trades_day.get(i));
       perform.println(perf_dates.get(i) + " " + perf_returns.get(i));
       date_returns.add(perf_dates.get(i) + " " + perf_returns.get(i));
       mean_ntrades = mean_ntrades + n_trades_day.get(i);
      }
      mean_ntrades = mean_ntrades/n_trades_day.size();
      
      
       
      
      
      dates_price = new String[n_obs];
      
      for(i=1;i<n_obs;i++)
      {
      
        tokens = allsignal.get(allsignal.size() - n_obs + i).split("[ ]+");
        
        
        if(perf_datestimes.contains(latestDates.get(latestDates.size() - n_obs + i)))
        {
          pret = perf_rets.get(perf_datestimes.indexOf(latestDates.get(latestDates.size() - n_obs + i)));
          System.out.println(latestDates.get(latestDates.size() - n_obs + i) + " " + pret);
        }
        else
        {pret = 0;}
        
        //System.out.println(tokens[0] + " " + latestDates.get(latestDates.size() - n_obs + i - 1));
        
//         if(tokens[0].equals(latestDates.get(latestDates.size() - n_obs + i - 1)))
//         {
         dates_price[i] = new String(latestDates.get(latestDates.size() - n_obs + i) + " " + (mid.get(mid.size() - n_obs + i) - mid.get(mid.size() - n_obs + i-1)) 
           + " " + mid.get(mid.size() - n_obs + i) + " " + tokens[2] + " " + pret);
//        }
      }  
      dates_price[0] = dates_price[1];
  
  
  
      mean = mean/(n_pos_ret + n_neg_ret);
  
      System.out.println("ROI = " + account.get(account.size()-1));
      //--- compute stats---------------
      double risk = neg_ret_mean/(double)n_neg_ret;
      System.out.println("neg_ret_mean = " + (-neg_ret_mean) + ", " + n_neg_ret);
        
      double reward = pos_ret_mean/(double)n_pos_ret;
      System.out.println("pos_ret_mean = " + pos_ret_mean + ", " + n_pos_ret);
        
      double win_ratio = (double)(n_pos_ret)/(n_pos_ret + n_neg_ret);
        
      kellyPerc = win_ratio - (1.0 - win_ratio)*(risk/reward);
      ulcer_index = ulcerIndex(dreturns); 
        
      System.out.println("win ratio = " + win_ratio + ", risk = " + risk + ", reward = " + reward);
      System.out.println("kelly and ulcer = " + kellyPerc + " " + ulcer_index);
        
      for(i=0;i<dreturns.length;i++)
      {sd = sd + (dreturns[i] - mean)*(dreturns[i] - mean)/((double)dreturns.length);}
        
      standard_deviation = Math.sqrt(sd);
 
      sharpeRatio = Math.sqrt(250)*mean/standard_deviation;
      maxdraw = computeDrawdown(dreturns);        
      rank_coeff = segmentRankCorrelation(30, dreturns);      
 
      
      System.out.println("MeanRet = " + mean + ", Sharpe = " + sharpeRatio + ", MaxDD = " + maxdraw + ", Rank = " + rank_coeff + ", avg_n_trades = " + mean_ntrades);
 
      out.close(); dailyout.close(); perform.close(); b0_coeff.close();
      
      }
      catch(ParseException fe){System.out.println("ParseException");}
      catch(NullPointerException fe){System.out.println("Null pointer");}
      catch(IllegalArgumentException fe){System.out.println("IllegalArgument");}
      catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
      catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
 
      n_files++;
     
      return computed; 
      
  }                    
  
  
  
  
  
  
  
    
  
  



  
 


  
  
  
  
  
  
  
  


    public void updateEndTimeAsia(String time, int start_hour, int close_hour, int liquide_hour, int liquide_min)
    {
        String date_delims = "[ ]+";
        String ddelims = "[-]"; 
        
        String[] date_tokens = time.split(date_delims); 
        String[] intdates = date_tokens[0].split(ddelims);     
  
        int signal_min = liquide_min;
        
        try
        {
            startTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),start_hour,signal_min);
        }
        catch(IllegalInstantException ie)
        {
            startTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),start_hour+1,signal_min);
        }

        endTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),close_hour,signal_min);
        liquidateTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), liquide_hour, liquide_min);
        
        if(liquide_hour > start_hour && liquidateTimeDT.dayOfWeek().getAsText().equals("Friday"))
        {liquidateTimeDT = liquidateTimeDT.plusDays(2);}
        
        //--- assume starting after 10am HST (16.00 EST or 4pm EST) produces new system setup for next day
        //System.out.println(cur_hour);
  
         if(close_hour < start_hour)
         {
             endTimeDT = endTimeDT.plusDays(1); 
             
             if(endTimeDT.dayOfWeek().getAsText().equals("Saturday"))
             {
                 endTimeDT = endTimeDT.plusDays(2);
             }   
         }
         else if(endTimeDT.dayOfWeek().getAsText().equals("Friday"))
         {
             endTimeDT = endTimeDT.plusDays(2);
         }
         
         if(liquide_hour < start_hour)
         {liquidateTimeDT = liquidateTimeDT.plusDays(1);} 
         
         if(liquidateTimeDT.dayOfWeek().getAsText().equals("Saturday"))
         {liquidateTimeDT = liquidateTimeDT.plusDays(2);}
         
         //System.out.println("Updated start/close at " + startTimeDT + " " + endTimeDT + " " + liquidateTimeDT);
         
         
    }         
  
  
    public void updateEndTime(String time, int start_hour, int close_hour, int liquide_hour, int liquide_min)
    {
        String date_delims = "[ ]+";
        String ddelims = "[-]"; 
        
        String[] date_tokens = time.split(date_delims); 
        String[] intdates = date_tokens[0].split(ddelims);     
  
        int signal_min = liquide_min;
        
        
        
        try
        {
            startTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),start_hour,signal_min);
        }
        catch(IllegalArgumentException ie)
        {
            
            startTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),start_hour+1,signal_min);
        }

       
        
        try
        {
         endTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),close_hour,signal_min);
        }
        catch(IllegalArgumentException ie)
        {          
          endTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(),close_hour+1,signal_min);
        }
        
        liquidateTimeDT = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), liquide_hour, liquide_min);
        //System.out.println("here2??");
        
         if(close_hour < start_hour)
         {
             endTimeDT = endTimeDT.plusDays(1); 
             
             if(endTimeDT.dayOfWeek().getAsText().equals("Saturday"))
             {
                 endTimeDT = endTimeDT.plusDays(2);
             }   
         }

         //System.out.println("here3??");
         if(liquide_hour < start_hour)
         {liquidateTimeDT = liquidateTimeDT.plusDays(1);}          
         
         if(liquidateTimeDT.dayOfWeek().getAsText().equals("Saturday"))
         {liquidateTimeDT = liquidateTimeDT.plusDays(2);}        
        
        
        
        if(start_hour > 17 && startTimeDT.dayOfWeek().getAsText().equals("Friday"))
        {
          startTimeDT = startTimeDT.plusDays(2);
        }        
        
        
        if(liquide_hour > 17 && liquidateTimeDT.dayOfWeek().getAsText().equals("Friday"))
        {
          liquidateTimeDT = liquidateTimeDT.plusDays(2);
        }
        
        if(close_hour > 17 && endTimeDT.dayOfWeek().getAsText().equals("Friday"))
        {
          endTimeDT = endTimeDT.plusDays(2);
        }
        
        
         
         
    }   
  
  
 
  
    
  
  
  
  

  
  
  
  
  
  
  
  
  
  public double getAmount(String d, int r)
  {
    String[] as = d.split("[ ]+");
    
    if(r < as.length)
    {
     return (new Double(as[r])).doubleValue();
    }
    else
    {return 0;}
  }
  
  
  public double getAmount(String d)
  {
    String[] as = d.split("[ ]+");
    

     return (new Double(as[2])).doubleValue();

  }
  
  
  
  
  
  
  
  
  
  
  
  
  public void setStrategies(ArrayList<StrategyParameters> strats)
  {filter_strategies = strats;}
  
  public void setAssetUniverse(String[] assetUniv)
  {asset_name = assetUniv;}
  
  
  
  private static class Result 
  {
        public Result(String stats, ArrayList<String> per, ArrayList<Double> r, ArrayList<Double> r2,ArrayList<Double> r3) 
        { 
        }
  }  
  

  public static Result compute(MDFAStrategyEvolution strategy) throws InterruptedException 
  {
        
     boolean output = strategy.startStrategy();
     if(output)
     {
      double ratio_trades = strategy.trade_succ_ratio/(double)strategy.returns.size();
      strategy.returns.size();
      double sharpe = Math.sqrt(250)*strategy.mean_perf/strategy.standard_deviation;
      String perf = new String(" n_obs = " + strategy.n_obs + ", avg_rank = " + strategy.rank_coeff + ", maxdraw = " + strategy.maxdraw + ", mean = " + strategy.mean_perf + ", success_ratio = " + ratio_trades + ", sharpe " + sharpe + ", rate = " + strategy.avg_rate + ", avg_n_trades = " + strategy.avg_n_trades);        
      
      return new Result(perf,strategy.date_returns,strategy.returns,strategy.longreturns,strategy.shortreturns);
     }
     else
     {
       String perf = new String("");
       return new Result(perf,strategy.date_returns,strategy.returns,strategy.longreturns,strategy.shortreturns);
     }
  }  
  
  
  //----------- Splits data file into m_days seperate components to such that each processor can compute each seperately
  
  
  
  
  public static String[] forkDataFile(File file, int m_days, int nobs, int cushion)
  {
     
     int i,m;
          
     ArrayList<String> data_line = new ArrayList<String>();
     //int cushion = 31;
    
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
     
     //System.out.println("local_lines = " + local_lines + "mdays = " + m_days + ", nlines = " + n_lines);
     int start = 0; m = 0;
     while(start < n_lines)
     { 
       File file1 = new File(""+file.getName()+"_"+m);
       PrintWriter out = new PrintWriter(new FileWriter(file1));    
       
       for(i = start; i < start + local_lines; i++)
       {
        if(i < n_lines && i >= 0) 
        {out.println(data_line.get(i));}
        else {break;}        
       }
       out.close();
       start = start + local_lines - nobs - cushion;
       if(start < 0) {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       forkedFiles.add(""+file.getName()+"_"+m);
       m++;
     } 
    }
    catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
    
    return forkedFiles.toArray(new String[0]);
  }
  
  
  public Double log(Double g) 
  {
    if(forex24) {return g;}
    else {return Math.log(g);}
  }
  
  public Double Log(Double g) 
  {
    if(forex24) {return g;}
    else {return Math.log(g);}
  }  
  
  public static String[] forkDataFileExp(File file, int m_days, int nobs, ArrayList<String> expvar)
  {
     
     int i,m;

     ArrayList<String> data_line = new ArrayList<String>();
     int cushion = 27;
    
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
       
       for(i = start; i < start + local_lines; i++)
       {
        if(i < n_lines && i >= 0) 
        { 
         if(n_lines == expvar.size())
         {
          tokens = (expvar.get(i)).split(delims);
          out.println(data_line.get(i) + ", " + tokens[2]);
         }
         else
         {out.println(data_line.get(i));}
        }
        else {break;}        
       }
       out.close();
       start = start + local_lines - nobs - cushion;
       if(start < 0) {if(file1.delete()) {System.out.println(file1.getName() + " has been deleted");} break;}
       forkedFiles.add(""+file.getName()+"_"+m);
       m++;
     } 
    }
    catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
    catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}       
    
    return forkedFiles.toArray(new String[0]);
  }
    
  public void setFilterFile(File file)
  {
     String[] tokens; String delims = "[ ]+";
     String strline;  
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
              strategy = new StrategyParameters("Filter_0", tokens[0], tokens[1], (new Integer(tokens[2])).intValue(),
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
             }           
           }
    
           din.close();
   }               
   catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
   catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}      
   
  }    

  
  
  public void setFilterFileConfig(File file)
  {
     String[] tokens; String delims = "[ ]+";
     String strline;  
     int n_files = 0;
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
//              System.out.println(strline);
//              System.out.println(tokens.length);
                          
             String LiquidateTime = tokens[0];
             String StopTime = tokens[1];
             String StartTime = tokens[3];
             
      
             String[] liquidHours = LiquidateTime.split("[:]+");
             int start_min = (new Integer(liquidHours[1])).intValue();
      
             String[] starthours = StartTime.split("[:]+");
             int start_hour = (new Integer(starthours[0])).intValue();
             
             String[] endhours = StopTime.split("[:]+");
             int close_hour = (new Integer(endhours[0])).intValue();
        
             (new Integer(tokens[4])).intValue(); 
             int L = (new Integer(tokens[5])).intValue();
             (new Integer(tokens[6])).intValue();
             double omega0 = (new Double(tokens[7])).doubleValue();
   
             double sigma = (new Double(tokens[8])).doubleValue();
             double delta = (new Double(tokens[9])).doubleValue(); 
             double delta2 =(new Double(tokens[10])).doubleValue(); 
             (new Double(tokens[11])).doubleValue();
             (new Double(tokens[12])).doubleValue();
             (new Double(tokens[13])).doubleValue();
      
             int fc = (new Integer(tokens[14])).intValue();  
             (new Integer(tokens[15])).intValue();     
             int i2 = (new Integer(tokens[16])).intValue();
      
             (new Double(tokens[17])).doubleValue();
             (new Double(tokens[18])).doubleValue();  
             (new Double(tokens[19])).doubleValue();
      

             (new Integer(tokens[20])).intValue();  
             int sig_diff = (new Integer(tokens[21])).intValue();  
             int sig_inv = (new Integer(tokens[22])).intValue();  
             int stop_loss = (new Integer(tokens[23])).intValue(); 
             int tp = 0;
             if(tokens.length > 24)
             {
              tp = (new Integer(tokens[24])).intValue();              
             }             
             // "Filter-0" -> signalConfig( 4, 0,  8, 0, filterConfig(5, 0.28, 0.3, 0.1, 0.1,  1, 0).build),             
                          
             String printme = "\"Filter-" + n_files + "\" -> signalConfig(" + start_hour + ", " + start_min +  ", " + close_hour + ", " + start_min + ", " + (0.0001*stop_loss) + ", 0,  " + (0.0001*tp) 
               + ", " + tokens[2] + ")(filterConfig(" 
               + L + ", " + omega0 + ", " + sigma + ", " +delta + ", " +delta2 + ", " +fc + ", " +sig_diff + ", "; //build),";
             
             if(i2 == 0)
             {printme = printme + "useI2 = false,";}
             else 
             {printme = printme + "useI2 = true,";}
             
             if(sig_inv == 0) {printme = printme + "coeffScale = 1.0";}
             else {printme = printme + "coeffScale = -1.0)),";}

             
             

             System.out.println(printme);
             n_files++;
           
           }
           
   
           din.close();
   }               
   catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
   catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}      
   
  }   
  
  
  
  public void setFilterFilesiMetrica(File file, String name, int sl, int sigFreq, int pnlFreq, double spreadFilter, int shift)
  {
  
      String[] tokens; 
     String strline; int i,m;
     int n_files = 0;
     
     ArrayList<String> params = new ArrayList<String>();
     String header;  
     try
     {
          
           FileInputStream fin = new FileInputStream(file);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           //names = new ArrayList<String>();
           while((strline = br.readLine()) != null)
           {
           
             tokens = strline.split("[ ]+");
             
             String[] liquidHour = tokens[0].split("[:]+");
             String[] endHour = tokens[1].split("[:]+");
             String[] startHour = tokens[3].split("[:]+");
           
             DateTime liquidTime = new DateTime(2015, 11, 17, (new Integer(liquidHour[0])).intValue(), (new Integer(liquidHour[1])).intValue());
             DateTime endTime = new DateTime(2015, 11, 17, (new Integer(endHour[0])).intValue(), (new Integer(liquidHour[1])).intValue());
             DateTime startTime =  new DateTime(2015, 11, 17, (new Integer(startHour[0])).intValue(), (new Integer(liquidHour[1])).intValue());
           
             liquidTime = liquidTime.minusHours(shift);
             endTime = endTime.minusHours(shift);
             startTime = startTime.minusHours(shift);
           
             int newLiquidhour = liquidTime.getHourOfDay();
             int newEndhour = endTime.getHourOfDay();
             int newStarthour = startTime.getHourOfDay();
           
             String newLiquid = "" + newLiquidhour + ":" + liquidHour[1] + ":00";
             String newEnd = "" + newEndhour + ":" + liquidHour[1];
             String newStart = "" + newStarthour + ":" + liquidHour[1];
           
             strline = newLiquid + " " + newEnd + " " + tokens[2] + " " + newStart + " "; 
             
             for(i = 4; i < tokens.length; i++)
             {
              strline = strline + tokens[i] + " ";
             }
           
             params.add(strline);
             
             
             
             n_files++;
           }
           
           header = "String[] " + name + "_parameters = new String[" + n_files + "];";
           System.out.println(header);
  
           for( m = 0; m < n_files; m++)
           {
             System.out.println(name + "_parameters[" + m + "] = \"" + params.get(m) + pnlFreq + " " + sigFreq + " " + spreadFilter + " " + name + "-signal-" + m + "\";");
           }
  
           din.close();
   }               
   catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
   catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}      
  
  }
  
  
  
  
  

  
  public void setFilterFiles(File file)
  {
      
     String[] tokens; String delims = "[ ]+";
     String strline;  
     strategies = new ArrayList<StrategyParameters>();
     String[] timeTokens,startTokens,endTokens; 
     String startTime, endTime;
     
     sampMin = "00";
   
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
                          
             timeTokens = tokens[0].split("[:]+");
             
             sampMin = timeTokens[1];
             
             endTokens = tokens[1].split("[:]+");
             startTokens = tokens[3].split("[:]+");
             
             endTime = endTokens[0] + ":" + sampMin;
             startTime = startTokens[0] + ":" + sampMin;
                          
             if(tokens.length == 24)
             {
              StrategyParameters sp = new StrategyParameters("Filter_"+strategies.size(), tokens[0], endTime, (new Double(tokens[2])).doubleValue(),
                          startTime, (new Integer(tokens[4])).intValue(), (new Integer(tokens[6])).intValue(),
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
             
             }    
             else if(tokens.length == 25)
             {
              StrategyParameters sp = new StrategyParameters("Filter_"+strategies.size(), tokens[0], endTime, (new Double(tokens[2])).doubleValue(),
                          startTime, (new Integer(tokens[4])).intValue(), (new Integer(tokens[6])).intValue(),
                          (new Integer(tokens[5])).intValue(), (new Integer(tokens[14])).intValue(),                           
                          (new Double(tokens[7])).doubleValue(), (new Double(tokens[12])).doubleValue(), 
                          (new Double(tokens[13])).doubleValue(),                          
                          (new Double(tokens[8])).doubleValue(),  (new Double(tokens[9])).doubleValue(), 
                          (new Double(tokens[10])).doubleValue(),  (new Double(tokens[11])).doubleValue(), 
                          (new Integer(tokens[15])).intValue(), (new Integer(tokens[16])).intValue(),
                          (new Double(tokens[17])).doubleValue(),  (new Double(tokens[18])).doubleValue(),
                          (new Double(tokens[19])).doubleValue(), (new Integer(tokens[20]).intValue()), (new Integer(tokens[21]).intValue()),
                          (new Integer(tokens[22])).intValue(), (new Integer(tokens[23]).intValue()), (new Integer(tokens[24]).intValue()));
             strategies.add(sp);  
             
             }
           }
    
           din.close();
   }               
   catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
   catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}      
   
  }      
  
  
   
  public boolean isSunday(String d)
  {
     
     String[] intdates = d.split("[-]+");     
     DateTime weekend = new DateTime((new Integer(intdates[0])).intValue(), (new Integer(intdates[1])).intValue(), (new Integer(intdates[2])).intValue(), 16, 0); 
   
     return (weekend.dayOfWeek().getAsText().equals("Sunday"));
     
  }
   
  
  

  
  
  
  
  
      
  
  public double getNextIV(double t1, double t2)
  {
     double t0=0;
     
     if(t1 > .60 )
     {t0 = 0;}
     else if(t1 <= .061 && t2 <= .061)
     {t0 = .42;}
     else if(t1 >= .30 && t2 >= .15)
     {t0 = .88;}
     else 
     {t0 = t1 + .24;}
  
     return t0;  
  }
      
      
      
  //--------- Additional functions -----------------------------------------    
  public void insampleTradingDiff(double[] _price, double[] sig, int n)
  {
    //System.out.println("Short sell is " + short_sell);
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    double price_diff;
    xf_turn_val = new double[n];
   
    succ_trades = 0; total_trades = 0;
    
    first_trade_loss = false;
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    pnl = new double[n];
    pnl[0] = 0.0;
    
    int trade_obs_loc = n;

    for(i=0;i<trade_obs_loc;i++)
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

    double last_price;
    
//     while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = 0;
      
    last_price = prix[0];
    
    if(morning_buy)
    {
      //System.out.println("Signal at beginning of day = " + t_signal[0]);
      if(t_signal[0] > 0 && long_buy) //positive signal to start day, buy
      {
         price_bought = prix[0];
         xf_turn_val[0] = 1;
	 in_transaction = 1;
	 last_trade = 0;
      }
      else if(t_signal[0] < 0.0 && short_sell)
      {
      	 price_borrowed = prix[0];
         out_transaction = 1;	
         xf_turn_val[0] = 0;
	 last_trade = 0;
      }
      
    }
    
    
      
    amount = 0.0; 
    for(i = start; i < trade_obs_loc-2; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {
        last_price = prix[i+1]; last_trade = i+1;
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

        if(long_buy && in_transaction == 0)
	{
	 xf_turn_val[i+1] = 1;
	 price_bought = prix[i+1];
	 
	 in_transaction = 1; 
	} 
       
      }
      else if(t_signal[i+1] < -tradeDelta && t_signal[i] > 0) //if in transaction and signal goes below, sell
      {
        last_price = prix[i+1]; last_trade = i+1;
        if(long_buy && in_transaction == 1)
        {
	 xf_turn_val[i+1] = -1;
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
	if(short_sell && out_transaction == 0)
	{
	 
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;	 
	}
       	
     }
     
     if(in_transaction == 1 || out_transaction == 1)
     {
      price_diff = prix[i+1] - last_price;
      if(t_signal[i+1] > 0)
      {pnl[i+1] = price_diff;} //pnl[i+1] = pnl[i] + price_diff;}
      else
      {pnl[i+1] = -price_diff;} //pnl[i+1] = pnl[i] - price_diff;}
     }
     else
     {pnl[i+1] = 0;}
     
     //System.out.println("stop loss = " + stop_loss + " " + stop_loss_thresh);
     if(stop_loss)
     { 
       //if unrealized loss less than specified amout get rid of position
       //System.out.println("unrealized loss of " + (pnl[last_trade] - pnl[i+1])); 
       if((pnl[last_trade] - pnl[i+1]) > stop_loss_thresh)
       {
         
        if(out_transaction == 1) //in a short-sell transaction, sell
        { 
        
          if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from short position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
          price_sold = prix[i+1]; 
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
          xf_turn_val[i+1] = 1;
        }
        else if(in_transaction == 1)
        {
        
         if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from long position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
	 price_sold = prix[i+1];
	 profit = price_sold - price_bought;
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 

	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 in_transaction = 0;
	 xf_turn_val[i+1] = -1; 
                  
	}          
       }     
     }
     
     
   }
   //---- close any transaction-------
   price_diff = prix[trade_obs_loc-1] - last_price;
   if(t_signal[trade_obs_loc-1] > 0)
   {pnl[trade_obs_loc-1] = price_diff;}
   else
   {pnl[trade_obs_loc-1] = -price_diff;}   
   
   
   account[trade_obs_loc-2] = amount;
   if(in_transaction == 1)
   {
      
      price_sold = prix[trade_obs_loc-1]; 
      profit = price_sold - price_bought;
	 
      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
      amount = account[trade_obs_loc-1]; 
      in_transaction = 0; 
      final_trade = profit;
   }
   else if(out_transaction == 1)
   {

      price_sold = prix[trade_obs_loc-1];
      profit = price_borrowed - price_sold;

      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      //printf("total = %u\n", *total_trades);
      account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
      amount = account[trade_obs_loc-1];
      out_transaction = 0;
      final_trade = profit;
      
   }
   else //nothing to close, already in neutral 
   {account[trade_obs_loc-1] = account[trade_obs_loc-2];}
   //else{System.out.println("What the fuck????");}  //else  //not in transaction, last worth is previous observation worth

//    System.out.println("");
//    if(print_debug && (short_sell && long_buy)){
//    for(i = 0; i < n; i++)
//    {
//      System.out.println(i + " " + prix[i] + " " + t_signal[i] + " " + account[i] + " " + pnl[i] + " " + xf_turn_val[i]); 
//    }}

  
 }  

 
  //--------- Additional functions -----------------------------------------    
  public void insampleTradingDiff_Cust(double[] _price, double[] sig, int n)
  {
    //System.out.println("Short sell is " + short_sell);
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    double price_diff;
    xf_turn_val = new double[n];
    
    end_time_index = n-1;
    boolean store_closed = false;
    succ_trades = 0; total_trades = 0;
    
    boolean finished = false;
    first_trade_loss = false;
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    pnl = new double[n];
    pnl[0] = 0.0;
    
    int trade_obs_loc = n;

    for(i=0;i<trade_obs_loc;i++)
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

    double last_price;
    
//     while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = 0;
      
    last_price = prix[0];
    
    if(morning_buy)
    {
      //System.out.println("Signal at beginning of day = " + t_signal[0]);
      if(t_signal[0] > 0 && long_buy) //positive signal to start day, buy
      {
         price_bought = prix[0];
         xf_turn_val[0] = 1;
	 in_transaction = 1;
	 last_trade = 0;
      }
      else if(t_signal[0] < 0.0 && short_sell)
      {
      	 price_borrowed = prix[0];
         out_transaction = 1;	
         xf_turn_val[0] = 0;
	 last_trade = 0;
      }
      
    }
    
    
      
    amount = 0.0; 
    for(i = start; i < trade_obs_loc-2; i++)
    {

      if((i+1) == last_trade_index) //last trade index, store closing
      {
       store_closed = true;
       //if not currently in a transaction, we're done
       if(out_transaction == 0 && in_transaction == 0) {end_time_index = i+1; break;} 
      }
      
      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {
        last_price = prix[i+1]; last_trade = i+1;
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
   
          if(store_closed) {end_time_index = i+1; finished = true; break;}
        }

        if(long_buy && in_transaction == 0)
	{
	 xf_turn_val[i+1] = 1;
	 price_bought = prix[i+1];
	 
	 in_transaction = 1; 
	} 
       
      }
      else if(t_signal[i+1] < -tradeDelta && t_signal[i] > 0) //if in transaction and signal goes below, sell
      {
        last_price = prix[i+1]; last_trade = i+1;
        if(long_buy && in_transaction == 1)
        {
	 xf_turn_val[i+1] = -1;
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
	 
	 if(store_closed) {end_time_index = i+1; finished = true; break;}

	} 
	if(short_sell && out_transaction == 0)
	{
	 
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;	 
	}
       	
     }
     
     if(in_transaction == 1 || out_transaction == 1)
     {
      price_diff = prix[i+1] - last_price;
      if(t_signal[i+1] > 0)
      {pnl[i+1] = price_diff;} //pnl[i+1] = pnl[i] + price_diff;}
      else
      {pnl[i+1] = -price_diff;} //pnl[i+1] = pnl[i] - price_diff;}
     }
     else
     {pnl[i+1] = 0;}
     
     //System.out.println("stop loss = " + stop_loss + " " + stop_loss_thresh);
     if(stop_loss)
     { 
       //if unrealized loss less than specified amout get rid of position
       //System.out.println("unrealized loss of " + (pnl[last_trade] - pnl[i+1])); 
       if((pnl[last_trade] - pnl[i+1]) > stop_loss_thresh)
       {
         
        if(out_transaction == 1) //in a short-sell transaction, sell
        { 
        
          if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from short position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
          price_sold = prix[i+1]; 
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
          xf_turn_val[i+1] = 1;
        }
        else if(in_transaction == 1)
        {
        
         if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from long position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
	 price_sold = prix[i+1];
	 profit = price_sold - price_bought;
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 

	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 in_transaction = 0;
	 xf_turn_val[i+1] = -1; 
                  
	}          
       }     
     }
     
     
   }
   
   if(!finished)
   {
    //---- close any transaction-------
    price_diff = prix[trade_obs_loc-1] - last_price;
    if(t_signal[trade_obs_loc-1] > 0)
    {pnl[trade_obs_loc-1] = price_diff;}
    else
    {pnl[trade_obs_loc-1] = -price_diff;}   
   
   
    account[trade_obs_loc-2] = amount;
    if(in_transaction == 1)
    {
      
      price_sold = prix[trade_obs_loc-1]; 
      profit = price_sold - price_bought;
	 
      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
      amount = account[trade_obs_loc-1]; 
      in_transaction = 0; 
      final_trade = profit;
    }
    else if(out_transaction == 1)
    {

      price_sold = prix[trade_obs_loc-1];
      profit = price_borrowed - price_sold;

      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      //printf("total = %u\n", *total_trades);
      account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
      amount = account[trade_obs_loc-1];
      out_transaction = 0;
      final_trade = profit;
      
    }
    else //nothing to close, already in neutral 
    {account[trade_obs_loc-1] = account[trade_obs_loc-2];}
   //else{System.out.println("What the fuck????");}  //else  //not in transaction, last worth is previous observation worth
   
    end_time_index = trade_obs_loc - 1;
   }
   
   
   //System.out.println("");
//    if(print_debug && (short_sell && long_buy)){
//    System.out.println("");
//    for(i = 0; i < n; i++)
//    {
//      System.out.println(i + " " + prix[i] + " " + t_signal[i] + " " + account[i] + " " + pnl[i] + " " + xf_turn_val[i]); 
//    }}

  
 }   
 
 
 
 
 
  public void insampleTradingDiff_Cust_SL(double[] _price, double[] lo_price, double[] hi_price, double[] sig, int n)
  {
    //System.out.println("Short sell is " + short_sell);
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
    
    end_time_index = n-1;
    boolean store_closed = false;
    succ_trades = 0; total_trades = 0;
    boolean apply_stop = false;
    boolean end_trading = false;
    boolean finished = false;
    first_trade_loss = false;
    double[] prix = new double[n];
    double[] lo_prix = new double[n];
    double[] hi_prix = new double[n];
    double[] t_signal = new double[n];
    account = new double[n];
    pnl = new double[n];
    png = new double[n];
    png[0] = 0.0;
    pnl[0] = 0.0;
    
    int trade_obs_loc = n;

    for(i=0;i<trade_obs_loc;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i]; lo_prix[i] = lo_price[i]; hi_prix[i] = hi_price[i];}// System.out.println(prix[i] + "  " + t_signal[i]);}

       


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

    double last_price;
    
//     while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = 0;
      
    last_price = prix[0];
    
    if(morning_buy)
    {
      //System.out.println("Signal at beginning of day = " + t_signal[0]);
      if(t_signal[0] > 0 && long_buy) //positive signal to start day, buy
      {
         price_bought = prix[0];
         xf_turn_val[0] = 1;
	 in_transaction = 1;
	 last_trade = 0;
	 last_price = prix[0];
      }
      else if(t_signal[0] < 0.0 && short_sell)
      {
      	 price_borrowed = prix[0];
         out_transaction = 1;	
         xf_turn_val[0] = 0;
	 last_trade = 0;
	 last_price = prix[0];
      }
      
    }
    
    
      
    amount = 0.0; 
    for(i = start; i < trade_obs_loc-1; i++)
    {
      apply_stop = false;
      account[i] = amount;
      
      if((i+1) == last_trade_index) //last trade index, store closing
      {
       
       store_closed = true;
       
       if(out_transaction == 0 && in_transaction == 0) //if not currently in a transaction, we're done
       {
       
        account[i+1] = account[i];
        //System.out.println(i + " " + prix[i+1] + ", " + account[i] + " " + account[i+1]);
        end_time_index = i+1; 
        finished = true; break;        
       } 
      }
      
      

        
      
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {
      

        
        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
        
          pnl[i+1] = last_price - hi_prix[i+1];
          png[i+1] = last_price - lo_prix[i+1];        
        
        
          price_sold = prix[i+1];
          
          profit = price_borrowed - price_sold;


		  
	  if(stop_loss && (-pnl[i+1]) > stop_loss_thresh) {profit = -stop_loss_thresh;}
	  if(take_profit && png[i+1] >= take_profit_thresh) {profit = take_profit_thresh;}
	  
	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  
	  
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
          apply_stop = true;
   
          if(store_closed) {end_trading = true;} 
        }

        if(long_buy && in_transaction == 0)
	{
	 if(!apply_stop) {pnl[i+1] = 0; png[i+1] = 0;}
	 xf_turn_val[i+1] = 1;
	 price_bought = prix[i+1];
	 
         last_price = prix[i+1]; last_trade = i+1;	 
	 apply_stop = true;
	 in_transaction = 1; 
	} 
       
      }
      else if(t_signal[i+1] < -tradeDelta && t_signal[i] > 0) //if in transaction and signal goes below, sell
      {
      

        if(long_buy && in_transaction == 1)
        {

         pnl[i+1] = lo_prix[i+1] - last_price;
         png[i+1] = hi_prix[i+1] - last_price;       
       
         xf_turn_val[i+1] = -1;
	 price_sold = prix[i+1];
	 profit = price_sold - price_bought;
	 
	 //printf("profit = %lf\n", profit);
	 

	 
	 if(stop_loss && (-pnl[i+1]) > stop_loss_thresh) {profit = -stop_loss_thresh;}
	 if(take_profit && png[i+1] >= take_profit_thresh) {profit = take_profit_thresh;}
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 
	 //printf("total = %u\n", *total_trades);	 
	 
	 
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 //printf("account = %lf\n", account[i+1]);
	 in_transaction = 0;
	 apply_stop = true;  
	   
	 if(store_closed) {end_trading = true;} 

	} 
	if(short_sell && out_transaction == 0)
	{
	 
	  if(!apply_stop) {pnl[i+1] = 0; png[i+1] = 0;}
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;	 
          last_price = prix[i+1]; last_trade = i+1;
          apply_stop = true;
	}
       	
     }
     else
     {account[i+1] = account[i];} // if(store_closed) {end_trading = true;} }
     

     if(!apply_stop)
     {
     
      //-----------check take_profit or stop_loss first --------------- 
      if(in_transaction == 1 || out_transaction == 1)
      {
      
       if(t_signal[i+1] > 0)
       {
        pnl[i+1] = lo_prix[i+1] - last_price;
        png[i+1] = hi_prix[i+1] - last_price;
       } 
       else
       {
        pnl[i+1] = last_price - hi_prix[i+1];
        png[i+1] = last_price - lo_prix[i+1];
       } 
      }
      else
      {pnl[i+1] = 0; png[i+1] = 0;}
     
     //System.out.println("stop loss = " + stop_loss + " " + stop_loss_thresh);
      if(stop_loss)
      { 
       //if unrealized loss less than specified amout get rid of position
       //System.out.println("unrealized loss of " + (pnl[last_trade] - pnl[i+1])); 
       //if((pnl[last_trade] - pnl[i+1]) > stop_loss_thresh)
       if((-pnl[i+1]) > stop_loss_thresh)
       {
         
        if(out_transaction == 1) //in a short-sell transaction, sell
        { 
        
          if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from short position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
          price_sold = prix[i+1]; 
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  
	  profit = -stop_loss_thresh;
	  
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
          xf_turn_val[i+1] = 1;
          pnl[i+1] = -stop_loss_thresh;
        }
        else if(in_transaction == 1)
        {
        
         if(print_debug) System.out.println("Activated Stop loss at " + (i+1) + " from long position since last_trade " + last_trade + ": unrealized loss of " + (pnl[last_trade] - pnl[i+1]) + ", " + pnl[last_trade] + " " + pnl[i+1]);
	 price_sold = prix[i+1];
	 profit = price_sold - price_bought;
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 

	 profit = -stop_loss_thresh;
	 
	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 in_transaction = 0;
	 xf_turn_val[i+1] = -1; 
         pnl[i+1] = -stop_loss_thresh;         
	}  
	
	if(store_closed) {end_trading=true;}
	
        }     
      }
      
      if(take_profit)
      {
       if((png[i+1]) >= take_profit_thresh)
       {
         
        if(out_transaction == 1) //in a short-sell transaction, sell
        { 
        
          if(print_debug) System.out.println("Activated Take Profit at " + (i+1) + " from short position since last_trade " + last_trade + ": unrealized gain of " + (png[last_trade] - png[i+1]) + ", " + png[last_trade] + " " + png[i+1]);
          price_sold = prix[i+1]; 
          profit = price_borrowed - price_sold;

          profit = take_profit_thresh;
          
          
	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  
	  
	  
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
          xf_turn_val[i+1] = 1;
          png[i+1] = take_profit_thresh;
        }
        else if(in_transaction == 1)
        {
        
         if(print_debug) System.out.println("Activated Take Profit at " + (i+1) + " from long position since last_trade " + last_trade + ": unrealized gain of " + (png[last_trade] - png[i+1]) + ", " + png[last_trade] + " " + png[i+1]);
	 price_sold = prix[i+1];
	 profit = price_sold - price_bought;
	 
	 
	 profit = take_profit_thresh;
	 
	 if(profit > 0) {succ_trades=succ_trades+1;}
	 total_trades=total_trades+1; 

	 account[i+1] = account[i] + profit - trading_cost;   
	 amount = account[i+1]; 
	 in_transaction = 0;
	 xf_turn_val[i+1] = -1; 
         png[i+1] = take_profit_thresh;         
	}  
	
	if(store_closed) {end_trading=true;}
	
       }     
      }      
      
     }
     
     
     
          
     
      //System.out.println(i + " " + prix[i+1] + ", " + account[i] + " " + account[i+1]);

     if(end_trading) {end_time_index = i+1; finished = true; break;} 
     
     
   }
   
   if(!finished)
   {

    //System.out.println("Not finished, last call...");  
    account[trade_obs_loc-2] = amount;
    if(in_transaction == 1)
    {
      
      
      if(t_signal[trade_obs_loc-1] > 0)
      {
       //price_diff = prix[i+1] - last_price;
       //price_diff = lo_prix[trade_obs_loc-1] - last_price;
       pnl[trade_obs_loc-1] = lo_prix[trade_obs_loc-1] - last_price;
       png[trade_obs_loc-1] = hi_prix[trade_obs_loc-1] - last_price;
       
      } 
      else
      {
       //price_diff = prix[i+1] - last_price;
       //price_diff = hi_prix[trade_obs_loc-1] - last_price;
       pnl[trade_obs_loc-1] = last_price - hi_prix[trade_obs_loc-1];
       png[trade_obs_loc-1] = last_price - lo_prix[trade_obs_loc-1];
      }       
             
      
      if(-pnl[trade_obs_loc-1] > stop_loss_thresh && stop_loss)
      {
         pnl[trade_obs_loc-1] = -stop_loss_thresh;
         profit = -stop_loss_thresh;
         account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;  
      }
      else if(png[trade_obs_loc-1] > take_profit_thresh && take_profit)
      {
         png[trade_obs_loc-1] = take_profit_thresh;
         profit = take_profit_thresh;
         account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;  
      
      }
      else
      {
      
       price_sold = prix[trade_obs_loc-1]; 
       profit = price_sold - price_bought;
	 
       if(profit > 0) {succ_trades=succ_trades+1;}
       total_trades=total_trades+1; 
       
       if(stop_loss && (-profit) > stop_loss_thresh) {System.out.println("LOSS = " + profit); profit = -stop_loss_thresh;}
       if(take_profit && profit > take_profit_thresh) {System.out.println("PROFIT = " + profit); profit = take_profit_thresh;}
       
       account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
       amount = account[trade_obs_loc-1]; 
       in_transaction = 0; 
       final_trade = profit;      
       //pnl[trade_obs_loc-1] = 0;
      }  
      
    }
    else if(out_transaction == 1)
    {
    
      if(t_signal[trade_obs_loc-1] > 0)
      {
       //price_diff = prix[i+1] - last_price;
       //price_diff = lo_prix[trade_obs_loc-1] - last_price;
       pnl[trade_obs_loc-1] = lo_prix[trade_obs_loc-1] - last_price;
       png[trade_obs_loc-1] = hi_prix[trade_obs_loc-1] - last_price;
       
      } 
      else
      {
       //price_diff = prix[i+1] - last_price;
       //price_diff = hi_prix[trade_obs_loc-1] - last_price;
       pnl[trade_obs_loc-1] = last_price - hi_prix[trade_obs_loc-1];
       png[trade_obs_loc-1] = last_price - lo_prix[trade_obs_loc-1];
      }       
      
      if(-pnl[trade_obs_loc-1] > stop_loss_thresh && stop_loss)
      {
         pnl[trade_obs_loc-1] = -stop_loss_thresh;
         profit = -stop_loss_thresh;
         account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;  
      }
      else if(png[trade_obs_loc-1] > take_profit_thresh && take_profit)
      {
         png[trade_obs_loc-1] = take_profit_thresh;
         profit = take_profit_thresh;
         account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;  
      
      }      
      else
      {

       price_sold = prix[trade_obs_loc-1];
       profit = price_borrowed - price_sold;

       if(profit > 0) {succ_trades=succ_trades+1;}
       total_trades=total_trades+1; 
       //printf("total = %u\n", *total_trades);
       
       if(stop_loss && (-profit) > stop_loss_thresh) {System.out.println("LOSS = " + profit); profit = -stop_loss_thresh;}
       if(take_profit && profit > take_profit_thresh) {System.out.println("PROFIT = " + profit); profit = take_profit_thresh;}
       
       account[trade_obs_loc-1] = account[trade_obs_loc-2] + profit - trading_cost;   
       amount = account[trade_obs_loc-1];
       out_transaction = 0;
       final_trade = profit;
       //pnl[trade_obs_loc-1] = 0;
      } 
    }
    else //nothing to close, already in neutral 
    {account[trade_obs_loc-1] = account[trade_obs_loc-2]; pnl[trade_obs_loc-1] = 0; png[trade_obs_loc-1] = 0;}
   
   //else{System.out.println("What the fuck????");}  //else  //not in transaction, last worth is previous observation worth
    
    
    end_time_index = trade_obs_loc - 1;
   }
   

  
 }    
 
 
 
 
 
 
 
 
 
 
  public void insampleTradingDiffVolatility(double[] _price, double[] sig, double[] vola, double vola_thresh, int n)
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

    System.out.println("volatility_thresholding on - " + vola_thresh);


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

//     while(t_signal[i] >= 0 && i < t_signal.length-1) {i++;}  
    start = 0;
      
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && t_signal[i] < 0.0) //new point positive, we see momentum, buy
      {
       if((vola[i+1] >= vola_thresh) || (out_transaction == 1 && price_borrowed > prix[i+1]))
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

        if(in_transaction == 0)
	{
	 xf_turn_val[i+1] = 1;
	 price_bought = prix[i+1];
	 
	 in_transaction = 1;
	 
	} 
       }
      }
      else if(t_signal[i+1] < -tradeDelta && t_signal[i] > 0) //if in transaction and signal goes below, sell
      {
       if((vola[i+1] >= vola_thresh) || (in_transaction == 1 && prix[i+1] > price_bought))
       {
      
        if(in_transaction == 1)
        {
	 xf_turn_val[i+1] = -1;
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
	if(short_sell && out_transaction == 0)
	{
	 
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
          out_transaction = 1;	 
	}
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
   max_drop = dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
  
 }    
 
 
//   public void insampleTradingDouble(double[] _price, double[] sig, double[] lagsig, int n)
//   {
//   
//     int i,start;
//     double dXf = 0.0;
//     int n_xf_turn = 0;
//     int xf_pos; 
//     double price_bought, price_sold, price_borrowed;
//     double profit,amount;
//     double price;
//     xf_turn_val = new double[n];
//    
//     succ_trades = 0; total_trades = 0;
//    
//     string_trades = new String[n];
//     trades = new double[n];
//     double[] prix = new double[n];
//     double[] t_signal = new double[n];
//     double[] lag_signal = new double[n];
//     account = new double[n];
//     
//     trade_obs = n;
// 
//     for(i=0;i<trade_obs;i++)
//     {prix[i] = _price[i]; t_signal[i] = sig[i]; lag_signal[i] = lagsig[i]; string_trades[i] = "";}// System.out.println(prix[i] + "  " + t_signal[i]);}
// 
//        
// 
//     boolean positive = false;
//     // ------------boolean value, if currently owning or short_sell an asset in_transaction = 1-------
//     int in_transaction = 0;  
//     // ------------boolean value, if currently owning borrowed asset to sell cheaper (short sell------
//     int out_transaction = 0; 
//     // ------------boolean value, if currently waiting to buy asset waiting = 1-----------------------
//     int waiting = 0;    
//     // ------------number of days passed -------------------------------------------------------------
//     int days_passed = 0;
//    
//     account[0] = 0.0;   
//     double[] transaction = new double[trade_obs];
//      
//     n_xf_turn = 0; profit = 0; 
//     //--- find where to start, first negative value
//     i=0;
//     price_bought = 0.0;
//     price_borrowed = 0.0;
// 
// //     if(t_signal[i] > 0) {positive = true;}
// //     else {positive = false;}
// //     
// //     if(positive)
// //     {
// //       while(t_signal[i] >= 0 && i < t_signal.length) {i++;}  
// //       start = i;
// //     }
// //     else
// //     {
// //       while(t_signal[i] < 0 && i < t_signal.length) {i++;}  
// //       start = i;
// //     }
//     
//     start = 0;    
//     amount = 0.0; 
//     for(i = start; i < trade_obs-1; i++)
//     {
// 
//       account[i] = amount;
//       if(t_signal[i+1] > tradeDelta && lag_signal[i] < 0.0) //new point positive, we see momentum, buy
//       {
// 
//         if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
//         { 
//           price_sold = prix[i+1];
//           trades[i+1] = price_sold;
//           profit = price_borrowed - price_sold;
// 
//           string_trades[i+1] = "sold short position ";
// 	  if(profit > 0) {succ_trades=succ_trades+1;}
// 	  total_trades=total_trades+1; 
// 	  //printf("total = %u\n", *total_trades);
// 	  account[i+1] = account[i] + profit - trading_cost;   
// 	  amount = account[i+1];
//           out_transaction = 0;
//         }
// 
//       
// 	   xf_turn_val[i+1] = 1;
// 	   xf_pos = 1;
// 	   n_xf_turn++;
// 	   price_bought = prix[i+1];
// 	   trades[i+1] = price_bought;
// 	   in_transaction = 1;
// 	   string_trades[i+1] = string_trades[i+1] + "long";
// 	 
// 
//       }
//       else if(in_transaction == 1 && t_signal[i+1] < -tradeDelta) //if in transaction and signal goes below, sell
//       {
// 
// 	 xf_pos = 0;
// 	 n_xf_turn++;
// 	 xf_turn_val[i+1] = -1;
// 	 price_sold = prix[i+1];
// 	 trades[i+1] = -price_sold;
// 	 string_trades[i+1] = new String("sold");
// 	 profit = price_sold - price_bought;
// 	 
// 	 //printf("profit = %lf\n", profit);
// 	 
// 	 if(profit > 0) {succ_trades=succ_trades+1;}
// 	 total_trades=total_trades+1; 
// 	 //printf("total = %u\n", *total_trades);
// 	 account[i+1] = account[i] + profit - trading_cost;   
// 	 amount = account[i+1]; 
// 	 //printf("account = %lf\n", account[i+1]);
// 	 in_transaction = 0;
//       }
//       if(t_signal[i+1] < -tradeDelta && lag_signal[i] > 0.0) //went from positive to negative slope, short sell
//       {
//          //if(short_sell && out_transaction == 0)
//          if(short_sell)
//          {   
// 	  xf_turn_val[i+1] = -1;
// 	  xf_pos = 1;
// 	  n_xf_turn++;
// 	  price_borrowed = prix[i+1];
// 	  string_trades[i+1] = string_trades[i+1] + " and short";
//           out_transaction = 1;
//          }
//       }
//    }
//    //---- close any transaction-------
//    if(in_transaction == 1)
//    {
//       price_sold = prix[trade_obs-1]; 
//       profit = price_sold - price_bought;
//       trades[trade_obs-1] = -price_sold;	 
//       if(profit > 0) {succ_trades=succ_trades+1;}
//       total_trades=total_trades+1; 
//       account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
//       amount = account[trade_obs-1]; 
//       in_transaction = 0; 
//    }
//    else if(out_transaction == 1)
//    {
// 
//       price_sold = prix[trade_obs-1];
//       profit = price_borrowed - price_sold;
//       trades[trade_obs-1] = price_sold;
//       if(profit > 0) {succ_trades=succ_trades+1;}
//       total_trades=total_trades+1; 
//       //printf("total = %u\n", *total_trades);
//       account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
//       amount = account[trade_obs-1];
//       out_transaction = 0;
//    }
//    else{account[trade_obs-1] = account[trade_obs-2];}  //else  //not in transaction, last worth is previous observation worth
// 
//    //----- with new account data, compute dropDown();
//    max_drop = dropDown();
//   
//    double[] temp_acc = new double[trade_obs];
//    System.arraycopy(account, 0, temp_acc, 0, trade_obs);
//   
//  }    
//  
 

  public void insampleTradingDouble(double[] _price, double[] sig, double[] lagsig, int n)
  {
  
    int i,start;
    double price_bought, price_sold, price_borrowed;
    double profit,amount;
    xf_turn_val = new double[n];
   
    succ_trades = 0; total_trades = 0;
   
    string_trades = new String[n];
    trades = new double[n];
    double[] prix = new double[n];
    double[] t_signal = new double[n];
    double[] lag_signal = new double[n];
    account = new double[n];
    
    trade_obs = n;

    for(i=0;i<trade_obs;i++)
    {prix[i] = _price[i]; t_signal[i] = sig[i]; lag_signal[i] = lagsig[i]; string_trades[i] = "not in market";}// System.out.println(prix[i] + "  " + t_signal[i]);}

       

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
    
    start = 0;    
    amount = 0.0; 
    for(i = start; i < trade_obs-1; i++)
    {

      account[i] = amount;
      if(t_signal[i+1] > tradeDelta && lag_signal[i] < 0.0) //new point positive, we see momentum, buy
      {

        if(short_sell && out_transaction == 1) //in a short-sell transaction, sell
        { 
          price_sold = prix[i+1];
          trades[i+1] = price_sold;
          profit = price_borrowed - price_sold;

	  if(profit > 0) {succ_trades=succ_trades+1;}
	  total_trades=total_trades+1; 
	  //printf("total = %u\n", *total_trades);
	  account[i+1] = account[i] + profit - trading_cost;   
	  amount = account[i+1];
          out_transaction = 0;
        }

        if(in_transaction == 0)
	{
	 xf_turn_val[i+1] = 1;
	 price_bought = prix[i+1];
	 trades[i+1] = price_bought;
	 in_transaction = 1;
	 string_trades[i+1] = new String("long");
	} 

      }
      else if(t_signal[i+1] < -tradeDelta && lag_signal[i] > 0) //if in transaction and signal goes below, sell
      {
        if(in_transaction == 1)
        {
	 xf_turn_val[i+1] = -1;
	 price_sold = prix[i+1];
	 trades[i+1] = -price_sold;
	 string_trades[i+1] = new String("short");
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
	if(short_sell && out_transaction == 0)
	{
	 
	  xf_turn_val[i+1] = -1;
	  price_borrowed = prix[i+1];
	  string_trades[i+1] = new String("short");
          out_transaction = 1;	 
	}
     }
   }  
   //---- close any transaction-------
   if(in_transaction == 1)
   {
      price_sold = prix[trade_obs-1]; 
      profit = price_sold - price_bought;
      trades[trade_obs-1] = -price_sold;	 
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
      trades[trade_obs-1] = price_sold;
      if(profit > 0) {succ_trades=succ_trades+1;}
      total_trades=total_trades+1; 
      //printf("total = %u\n", *total_trades);
      account[trade_obs-1] = account[trade_obs-2] + profit - trading_cost;   
      amount = account[trade_obs-1];
      out_transaction = 0;
   }
   else{account[trade_obs-1] = account[trade_obs-2];}  //else  //not in transaction, last worth is previous observation worth

   //----- with new account data, compute dropDown();
   max_drop = dropDown();
  
   double[] temp_acc = new double[trade_obs];
   System.arraycopy(account, 0, temp_acc, 0, trade_obs);
  
 }     
 
 
 
 
 public void adaptiveUpdateFilter()
 {
   
    int j,l;    
    //--- copy the current filter
    mdfa.set_tseries(tseries,n_obs,n_rep);
    //mdfa.setCopyFilter(b_copy);       
       
    if(univariate == 1)
    {
         //mdfa.setUpdateConstraints(adapt_i1,adapt_i2);
         mdfa.updateSignalOut_Univ(adapt_n, adapt_L, adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross);
         b_update = new double[adapt_L];
         System.arraycopy(mdfa.b_update,0,b_update,0,adapt_L);
         
    }
    else
    {
        //mdfa.setUpdateConstraints(adapt_i1,adapt_i2);
        mdfa.updateSignalOut(adapt_n, adapt_L, adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross); 
        b_update = new double[(n_rep-1)*adapt_L];
      
        for(j=0;j<n_rep-1;j++)
        {
         for(l=0;l<adapt_L;l++)
         {b_update[adapt_L*j + l] = mdfa.b_update[adapt_L*(j+1) + l];}
        }
    }       
  } 
 

   public void updateAdaptiveSignal()
   {
   
     int i,l,j; double sum;
   
     int N = n_obs;
        
     int end = n_obs;
     int flength_update = adapt_n - adapt_L + 1;
     double[] tsdata = new double[(n_rep-1)*adapt_n];
     update_signal = new double[flength_update];        
     int start = end - adapt_n;
     
     if(univariate == 0)
     {
     //------ extract sequences ------     
      for(i=start;i<end;i++)
      {       
       for(j=0;j<n_rep-1;j++)
       {
         sum = 0.0;
         for(l=0;l<L;l++) {sum = sum + b_coeffs[L*j + l]*tseries[N*(j+1) + i-l];}         
         tsdata[adapt_n*j + i-start] = sum;
       }
      }              
           
        //---- with new tsdata of length n, apply update filter of length l
       
      for(i=adapt_L-1;i<adapt_n;i++)
      {
       sum = 0.0; 
       for(j=0;j<n_rep-1;j++)
       {     
         for(l=0;l<adapt_L;l++)
         {sum = sum + b_update[adapt_L*j + l]*tsdata[adapt_n*j + i-l]; }   
       }    
       update_signal[i-(adapt_L-1)] = sum; 
      }
     }
     else
     {
     
         tsdata = new double[adapt_n];
         for(i=start;i<end;i++)
         { 
          sum = 0.0;
          for(j=0;j<n_rep-1;j++)
          {
           for(l=0;l<L;l++) {sum = sum + b_coeffs[L*j + l]*tseries[N*(j+1) + i-l];}                     
          }
          tsdata[i-start] = sum;  //System.out.println(sum);
         }   
         
         for(i=adapt_L-1;i<adapt_n;i++)
         { 
           sum = 0.0;
           for(l=0;l<adapt_L;l++)
           {sum = sum + b_update[l]*tsdata[i-l];}
           update_signal[i-(adapt_L-1)] = sum; 
         }    
     } 
   } 
 
 

  public void rankCoefficient(double[] _account, int n)
  {

   int i; 
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double[] account = new double[n];
   System.arraycopy(_account,0,account,0,n);
   double sum = 0.0;
   double spear;
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i;} 
      
   sort_sims(n, account, index);
   
   for (i=0 ; i<n ; ++i)
   {d[i] = Math.abs(index[i] - rank[i]); d[i] = d[i]*d[i]; sum = sum + 1.0*d[i];} 
   
   spear = 1.0 - (6.0*sum)/(1.0*n*(n*n-1.0));
      
   rank_coeff = spear;   
  }  

  public double rankCoefficientSeg(double[] _account, int n)
  {

   int i; 
   int[] index = new int[n];
   int[] rank = new int[n]; 
   int[] d = new int[n];
   double[] account = new double[n];
   System.arraycopy(_account,0,account,0,n);
   double sum = 0.0;
   double spear;
   
   for (i=0 ; i<n ; ++i) {index[i] = i; rank[i] = i;} 
      
   sort_sims(n, account, index);
   
   for (i=0 ; i<n ; ++i)
   {d[i] = Math.abs(index[i] - rank[i]); d[i] = d[i]*d[i]; sum = sum + 1.0*d[i];} 
   
   spear = 1.0 - (6.0*sum)/(1.0*n*(n*n-1.0));
      
  
   
   return spear; 
   
  }  
 
 
//--------- computes largest drop_down given an account, also computes mean and std of returns -----
  public double dropDown()
  {
   int i;
   double max_drop = 0.0; 
   double drop;
   double mean, std;   
   int ticks = 0;
   int ltrade_obs = account.length;

   int n_drops = 0; mean = 0.0; std = 0.0;
   for(i=0;i<ltrade_obs-1;i++)
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
   for(i=0;i<ltrade_obs-1;i++)
   {
     drop = account[i+1] - account[i];
     if(drop != 0)
     {std = std + (drop - mean)*(drop - mean);}
   }
   std = Math.sqrt(std);
 
   return max_drop;
  } 
  
  public double[] minmax(double[] a)
  {
    double[] mm = new double[2];
  
    double min = 100; double max = -100;
        
    for(int i=0;i<a.length;i++)
    {
      if(a[i] < min) {min = a[i];}
      else if(a[i] > max) {max = a[i];}    
    }
    mm[0] = min;
    mm[1] = max;
    
    return mm;    
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
  
  
  public void averageFilters(int depth, double a)
  {
     int j,l,n;
     double weight;
     double[] avg_filter = new double[b_coeffs.length];
     
     weight = 1.0+a;
     for(j=0;j<n_rep-1;j++)
     {
       for(l=0;l<L;l++)
       {
         avg_filter[j*L + l] = b_coeffs[j*L + l]/weight;
       }
     }     
          
     for(n=0;n<depth;n++)
     {
      weight = Math.pow(1.0+a,n+2);
      Filter f = filters.get(filters.size()-1-n);
      for(j=0;j<n_rep-1;j++)
      {
       for(l=0;l<L;l++)
       {
         avg_filter[j*L + l] = avg_filter[j*L + l] + f.b_coeffs[(j+1)*L + l]/weight;
       }
      }
     }
     
     System.arraycopy(avg_filter,0,b_coeffs,0,avg_filter.length);
  }
  
  
   
    public void readFilterParamFile()
    {
         String strline; 
         file_valid = true;
         double[] parameters = new double[10];
         Double D; Integer I;
         //default parameters
         parameters[1] = 0.52; parameters[8] = 28;
         String[] tokens; String delims = "[ ]+";
         int n_toks;
         
         try
         {
          
           FileInputStream fin = new FileInputStream(filterParamFile);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           strline = br.readLine();
           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks <= 2)
           {
             D = new Double(tokens[0]); cutoff0 = D.doubleValue();
             D = new Double(tokens[1]); cutoff1 = D.doubleValue();                 
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

//            strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
//            I = new Integer(tokens[0]); i1 = I.intValue();
                      
           br.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}
      
        //System.out.println(cutoff0 + " " + cutoff1 + " " + lambda + " " + expweight + " " + smooth + " " + decay +
        //                   decay2 + " " + cross + " " + L + " " + lag);
    }  
  
  
  
    public void addFilter(String filterParamFile2)
    {
         int k;
         String strline; 
         file_valid = true;
         double[] parameters = new double[10];
         Double D; Integer I;
         //default parameters
         parameters[1] = 0.52; parameters[8] = 28;
         String[] tokens; String delims = "[ ]+";
         int K = (int)(n_obs/2); double om;
         int K1 = K+1;  
         double _cutoff0,_cutoff1,_lambda,_expweight,_smooth,_decay,_decay2,_cross; int _L,_lag;
         double[] _Gamma;
         
         System.out.println("Attempting to add filter " + K1);
         try
         {
          
           FileInputStream fin = new FileInputStream(filterParamFile2);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           strline = br.readLine();
           tokens = strline.split(delims); 
           D = new Double(tokens[0]); _cutoff0 = D.doubleValue();
           D = new Double(tokens[1]); _cutoff1 = D.doubleValue();                 
           
          
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _lambda = D.doubleValue();

           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _expweight = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _smooth = D.doubleValue();           
    
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _decay = D.doubleValue();    
   
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _decay2 = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); _cross = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); _L = I.intValue();           
           
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); _lag = I.intValue();

//            strline = br.readLine(); tokens = strline.split(delims); n_toks = tokens.length; 
//            I = new Integer(tokens[0]); i1 = I.intValue();
               
           mdfa_h0 =  new IMDFA(n_obs, n_rep, _L, 0, 0, 0, .52, 0, 1);
    
           mdfa_h0.set_L(_L); 
           mdfa_h0.set_nobs(n_obs);  
           mdfa_h0.set_nreps(n_rep);
           mdfa_h0.set_lag(_lag);
  
           mdfa_h0.setRegularization(_smooth, _decay, _decay2, _cross);
           mdfa_h0.set_dd(0);
           mdfa_h0.set_DD(0);
           mdfa_h0.set_bconstraints(i1, 0); 
           mdfa_h0.set_lambda(_lambda); 
           mdfa_h0.set_exp(_expweight); 
           mdfa_h0.set_cutoff0(_cutoff0);
           mdfa_h0.set_cutoff(_cutoff1);
           mdfa_h0.set_mdfa(true);
  
           _Gamma = new double[K1];       
            for(k=0; k<=K;k++)
            {       
               om = (k*Math.PI/K);
               if(om < _cutoff0) {_Gamma[k] = 0.0;}
               else if(om >= _cutoff0 && om <= _cutoff1)
               {_Gamma[k] = 1.0;}
               else
               {_Gamma[k] = 0.0;}
            }       
            mdfa_h0.set_Gamma(_Gamma);          
  
            System.out.println("New filter set");
            H0set = true;
            
            br.close();
       }
       catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
       catch(IOException ioe){System.out.println("IO out error..." + ioe);}
      
        //System.out.println(cutoff0 + " " + cutoff1 + " " + lambda + " " + expweight + " " + smooth + " " + decay +
        //                   decay2 + " " + cross + " " + L + " " + lag);
    }  
    
  
  
  
  
   public void readAdaptiveFilterParams()
   {
       //int adapt_L, adapt_n;
    //double adapt_lambda, adapt_alpha;
    //double adapt_sm, adapt_dec, adapt_dec2, adapt_cross;

         String strline; 
         file_valid = false; 
         Double D; Integer I;
         //default parameters

         String[] tokens; String delims = "[ ]+";
         try
         {
          
           FileInputStream fin = new FileInputStream(adaptfilterParamFile);
           DataInputStream din = new DataInputStream(fin);
           BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
           //----- first get the frequencies --------------
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); adapt_n = I.intValue();           
           
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); adapt_L = I.intValue(); 
          
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_lambda = D.doubleValue();

           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_alpha = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_sm = D.doubleValue();           
    
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_dec = D.doubleValue();    
   
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_dec2 = D.doubleValue();
           
           strline = br.readLine(); tokens = strline.split(delims);  
           D = new Double(tokens[0]); adapt_cross = D.doubleValue();
 
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); adapt_i1 = I.intValue(); 
           
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); adapt_i2 = I.intValue();    
           
           strline = br.readLine(); tokens = strline.split(delims);  
           I = new Integer(tokens[0]); univariate = I.intValue();            
           
           br.close();           
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}
      
        System.out.println(adapt_n + " " + adapt_L + " " + adapt_lambda + " " + adapt_alpha + " " + adapt_sm + " " + adapt_dec +
                           adapt_dec2 + " " + adapt_cross + " " + adapt_i1 + " " + adapt_i2 + " " + univariate);
                                 
              
   }  
  
  
  public void gaussianizeData(double p1, double p2, double m)
  {
     int i;    
     //diff_vol = true;
     int[] index = new int[n_obs+1];
     double[] vola = new double[n_obs+1];
     phi1 = p1; phi2 = p2; mu = m;
     int burnin = 81; 
    
     int le = n_obs; 
     int ex_length = n_obs + burnin;
     double[] extend_burn = new double[ex_length];
     double[] extend_burn2 = new double[ex_length];
     double[] extend_burn3 = new double[ex_length];

     double[] sigma = new double[ex_length+1];
     double[] sigma2 = new double[ex_length+1];
     double[] sigma3 = new double[ex_length+1];
     
     gauss_rt1 = new double[le];
     gauss_rt2 = new double[le];
     gauss_rt3 = new double[le];
   

     //rt = new double[le];
     vol = new double[le+1];
     vol2 = new double[le+1];
     vol3 = new double[le+1];
    
         
     //System.out.println("Gaussianizing Data - Erst Mal");
     
     vol[0] = 0; vol2[0] = 0; vol3[0] = 0;
      //-- mirror the data
     for(i=0;i<burnin;i++)
     {
       extend_burn[burnin-1-i] = tseries[i+1];
       if(n_rep > 2)
       {
       extend_burn2[burnin-1-i] = tseries[2*n_obs + i + 1];
       extend_burn3[burnin-1-i] = tseries[3*n_obs + i + 1];
       }
     }
     for(i=0;i<le;i++)
     {
       extend_burn[burnin+i] = tseries[i]; 
       if(n_rep > 2)
       {
       extend_burn2[burnin+i] = tseries[2*n_obs + i];
       extend_burn3[burnin+i] = tseries[3*n_obs + i];       
       }
     } 
    
     sigma[0] = 0; sigma2[0] = 0; sigma3[0] = 0;
      
     for(i=0;i<ex_length;i++)
     {
       sigma[i+1] = mu*(1.0 - phi1 - phi2) + phi1*sigma[i] + phi2*extend_burn[i]*extend_burn[i];
       if(n_rep > 2)
       {
       sigma2[i+1] = mu*(1.0 - phi1 - phi2) + phi1*sigma2[i] + phi2*extend_burn2[i]*extend_burn2[i];
       sigma3[i+1] = mu*(1.0 - phi1 - phi2) + phi1*sigma3[i] + phi2*extend_burn3[i]*extend_burn3[i];
       }
     } 
   
      
     for(i=0;i<le;i++)
     {
     
      if(gaussianize)
      {
       tseries[i] = tseries[i]/Math.sqrt(sigma[burnin+i]);
       tseries[n_obs+i] = tseries[n_obs+i]/Math.sqrt(sigma[burnin+i]);
       if(n_rep > 2)
       {
       tseries[2*n_obs+i] = tseries[2*n_obs+i]/Math.sqrt(sigma2[burnin+i]);
       tseries[3*n_obs+i] = tseries[3*n_obs+i]/Math.sqrt(sigma3[burnin+i]);
       }
       gauss_rt1[i] = tseries[i];
       if(n_rep > 2)
       {
       gauss_rt2[i] = tseries[2*n_obs+i]; 
       gauss_rt3[i] = tseries[3*n_obs+i]; 
       }
      }
      
      vol[i] = sigma[burnin+i]-1.5; 
      if(n_rep > 2)
      {
      vol2[i] = sigma2[burnin+i]-1.5; 
      vol3[i] = sigma3[burnin+i]-1.5; 
      }       
       
      index[i] = i;
     }
     vol[le] = sigma[ex_length]-1.5;
     if(n_rep > 2)
     {
     vol2[le] = sigma2[ex_length]-1.5;
     vol3[le] = sigma3[ex_length]-1.5;
     }
     index[n_obs] = n_obs;
  
     System.arraycopy(vol,0,vola,0,vol.length);
     sort_sims(n_obs+1, vola, index);     
  
     vola_thresh = vola[10];
  
     sigma_t = sigma[ex_length];
     if(n_rep > 2)
     {
     sigma_2t = sigma2[ex_length];
     sigma_3t = sigma3[ex_length];
     } 
  }
    
    
    
  public void gaussianizeDataUniv(double p1, double p2, double m)
  {
     int i;    
     //diff_vol = true;
     int[] index = new int[n_obs+1];
     phi1 = p1; phi2 = p2; mu = m;
     int burnin = 81; 
    
     int le = n_obs; 
     int ex_length = n_obs + burnin;
     double[] extend_burn = new double[ex_length];
     double[] sigma = new double[ex_length+1];

     gauss_rt1 = new double[le];
     vol = new double[le+1];
     vol[0] = 0;
     
      //-- mirror the data
     for(i=0;i<burnin;i++)
     {
       extend_burn[burnin-1-i] = tseries[i+1];
     }
     for(i=0;i<le;i++)
     {
       extend_burn[burnin+i] = tseries[i]; 
     } 
    
     sigma[0] = 0; 
      
     for(i=0;i<ex_length;i++)
     {
       sigma[i+1] = mu*(1.0 - phi1 - phi2) + phi1*sigma[i] + phi2*extend_burn[i]*extend_burn[i];

     } 
        
     for(i=0;i<le;i++)
     {

      tseries[i] = tseries[i]/Math.sqrt(Math.log(sigma[burnin+i]));
      gauss_rt1[i] = tseries[i];

      vol[i] = sigma[burnin+i];     
      index[i] = i;
     
     }
     vol[le] = sigma[ex_length]; 
     sigma_t = sigma[ex_length];

  }    
    
    
    


  public void updateGaussianData()
  {
  
     int i;
     sigma_t = mu*(1.0 - phi1 - phi2) + phi1*sigma_t + phi2*tseries[n_obs-1]*tseries[n_obs-1];
     
     //shift over 
     for(i=0;i<n_obs;i++) {vol[i] = vol[i+1]; tseries[i] = tseries[i]/Math.sqrt(Math.log(vol[i]));} 
     vol[n_obs] = sigma_t;
     
   }      
    
    
   /*------------------------------------------------------
     Function the optimizes over the interpolation between current filter 
     and H0. The interpolation is goverrned by the decay1,decay2, with decay1 determining 
     how the coefficients interpolate (0 for full, 1 for none) and decay2 determines how 
     strong the interpolation is (0 for none, 1 for full H0)
     Here the decay1 and all other parameters are kept fixed, and we do a grid search on decay2
     
     Input: n_outsamp 
      The number of out-of-sample points to include in the optimization
     Input: interp_delta, the sampling grid of the interpolation, eg .10 
   -------------------------------------------------------*/
    
   public void H0_Optimization(int n_outsamp, double localdecay, double delta, double interp_start)
   {
   
     int i,j,k,l;
     int trade_length = 100;
     int K,K1;
      
   
     int n = n_obs - n_outsamp;
     
     K = (int)(n/2); double om; K1 = K+1;  
     double[] _Gamma = new double[K1];       
     for(k=0; k<=K;k++)
     {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {_Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {_Gamma[k] = 1.0;}
         else
         {_Gamma[k] = 0.0;}
     }         
     
     double[] logprice = new double[trade_length];
     double[] out;
     double[] series = new double[n*n_rep];
   
     for(i=0;i<n_rep;i++)
     {
       for(j=0;j<n;j++)
       {series[n*i + j] = tseries[n_obs*i + j];}  
     }   
     
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-1-i] = price.get(price.size()-1-i).doubleValue();}
          
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;

     
     if(recomputeH0 && mdfa_h0 != null)
     {computeH0Filter(series);}
     
     //---- now use optimal interpolation value
     out = mdfa.H0_IMDFAreg(series, n, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, localdecay, 0, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }
     
     old_tseries = new double[tseries.length];
     System.arraycopy(tseries,0,old_tseries,0,tseries.length);
     
   } 
   
   
   
   
   
   public void H0_Optimization_Lookback(double delta, double interp_start)
   {
   
     int i,j,k,l; int count;
     int trade_length = trade_obs;  
     double interpolation = interp_start;
     int K1;
     int shift = 3;
    
     //start_trades 9 - n_obs = 385, rank = 0.5360576923076923, mean = 0.011154870280718323, success_ratio = 0.9538461538461539, sharpe 21.249762580002255
    
     double price_move;
     double full_rank; double full_rank_max = -1.0;
     int start_trades = 8;
     double prev_return = -100;
     double full_return;
     double standard_dev=0;
     double full_return_int = 0;
     double max_ret = -100; double max_rank = -1.5; 
     double max_interp = 0.0; 
     double full_return_max = -199;
     int n = n_obs; 
     int n_total = n_obs + trade_length-1;
     yesterday_series = new double[n_rep*n_total];
     double opt_rank_coeff;
     int actual_n = trade_length - start_trades + 1;
     double u;
     double[] short_sig = new double[start_trades];
     double[] short_price = new double[start_trades];
     
     double crite_0=0;
     double avg_vol = 0;
     
     double[] hist_sig = new double[n_obs-L+1];
     double[] hist_price = new double[n_obs-L+1];
     
     double[] full_returns = new double[100];
     
     for(i=0;i<n_obs;i++)
     {
       yesterday_series[i] = old_tseries[i];
       yesterday_series[n_total + i] = old_tseries[n_obs+i];
       if(n_rep>2){
       yesterday_series[2*n_total + i] = old_tseries[2*n_obs+i];
       yesterday_series[3*n_total + i] = old_tseries[3*n_obs+i];
       }
     }
     
     
     //--- create yesterday_series
     for(i=0;i<trade_length-1;i++)
     {
       yesterday_series[n_total - 1 - i] = gauss_target.get(gauss_target.size() - 1 - i); 
       yesterday_series[n_total + n_total - 1 - i] = gauss_target.get(gauss_target.size() - 1 - i);
       if(n_rep>2){
       yesterday_series[n_total*2 + n_total - 1 - i] = gauss_1.get(gauss_1.size() - 1 - i);
       yesterday_series[n_total*3 + n_total - 1 - i] = gauss_2.get(gauss_2.size() - 1 - i);       
       }
     }
          

     double opt_ret;
     K1 = (int)n_obs/2+1; 
     //--- create new Gamma
     double sum; double total_return;       
     int num_postive=0;
     int start; int flength = n_obs-L+1;
     double[] x = new double[actual_n];
     double[] logprice = new double[actual_n];
     double[] localsignal = new double[actual_n];
     double[] b_h0 = new double[n_rep*L];
     double[] out;
     double[] x_morn = new double[start_trades];
     double morn_vol;
     x_morn[0] = 0.0;
     double zero_ret = 0;
     double one_ret = 0; 
     
     
     //-------- Now get previous day's time series data 
     for(i=0;i<actual_n;i++)
     {logprice[actual_n-1-i] = price.get(price.size()-shift-i).doubleValue();}
     //logprice[0] = logprice[1];
     
     for(i=0;i<n_obs-L+1;i++)
     {hist_price[hist_price.length-1-i] = price.get(price.size()-shift-31-i).doubleValue();}
     
     for(i=0;i<start_trades;i++) {short_price[start_trades-1-i] = price.get(price.size()-shift-(actual_n-1)-i).doubleValue();}
     short_price[0] = short_price[1];
     //System.out.println(Math.exp(short_price[0]) + " " + Math.exp(short_price[1]) + Math.exp(short_price[2]));
     
     for(i=1;i<start_trades;i++) {x_morn[i] = (short_price[i] - short_price[i-1])*(short_price[i] - short_price[i-1]);}     
     
     sum = 0;
     for(i=1;i<start_trades;i++) {sum = sum + x_morn[i]/(double)start_trades;} 
     morn_vol = sum;
     
     computePeriodogram(old_tseries);
     
     price_move = short_price[0] - short_price[start_trades-1];         
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
     count = 0; 
     //System.out.println("price move = " + price_move);
     int n_interps = 24;
     double[] interpolations = new double[n_interps];
     double[] interp_max = new double[n_interps];
     interpolations[0] = 0;
     //interpolations[1]
     //interpolations[1] = .6;
     //interpolations[2] = .999;
     //interpolations[2] = .92;
     for(k=1;k<n_interps;k++) {interpolations[k] = interpolations[k-1] + 1.0/24.0;}
     interpolations[n_interps-1] = 1;
     //{0,0.1328}
     if(interp_vals.size() > 1) {interp_vals.get(interp_vals.size()-1);}     
     
     double[] insamp_stat = new double[n_interps];
//      if(prev_val == 0.0)
//      {fore_val = 0.0;}
//      else if(prev_val < .20)
//      {fore_val = .10;}
//      else if(prev_val > .20 && prev_val < .50)
//      {fore_val = .30;}
//      else if(prev_val > .80)
//      {fore_val = .04;}
     
     double[] avg_filt = new double[b_coeffs.length];
     
     if(recomputeH0 && mdfa_h0 != null)
     {computeH0Filter(old_tseries);}
     
     //System.out.println("H0 = " + h0b0[0] + " " + h0b0[1] + " " + h0b0[2]);
     
     //.133 .893
     //while(interpolation <= .999)
     for(k=0;k<n_interps;k++)
     {
     
       interpolation = interpolations[k];
       System.out.println("k = " + k + ", " + interpolation);
       //mdfa.setRegularization(mdfa.smooth, localdecay, interpolation, mdfa.cross);
       //mdfa.computeFilterGeneral(true);

       out = computeInterpolationFrequency(old_tseries, interpolation);
       
       if(k==0)
       {
        crite_0 = out[out.length-1];
        criteria = out[out.length-3];
        degrees = out[out.length-2];      
        MDFAmin = out[out.length-1];  
        //System.out.println("Crist = " + criteria + " " + degrees + " " + MDFAmin);
       }
       else if(k==n_interps-1)
       {
       }
       
       //System.out.println("interpolation value = " + interpolation);
       System.arraycopy(out,lag4,b_h0,0,n_rep*L);
       
       //System.arraycopy(mdfa.b,0,b_h0,0,mdfa.b.length);
       for(j=1;j<n_rep;j++)
       {
        for(l=0;l<L;l++) {avg_filt[(j-1)*L+l] = avg_filt[(j-1)*L+l] + b_h0[L*j + l]/(double)n_interps;} // System.out.println(yesterday_series[n_total*j + i-l]);}
       }       

       //apply on time series
       
       start = n_total-actual_n;
       for(i=start;i<n_total;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
           //for(l=0;l<L;l++) {sum = sum + avg_filt[(j-1)*L+l]*yesterday_series[n_total*j + i-l];}
        }
        x[i-start] = yesterday_series[i];
        localsignal[i-start] = sum;
       }
       standard_dev = stand_dev(x);
   
       if(k==n_interps-1)
       {
       start = n_total-actual_n;
       for(i=start;i<n_total;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + avg_filt[L*(j-1) + l]*yesterday_series[n_total*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
        }
        localsignal[i-start] = sum;
       }       
       
       
       }
   
       //if(k==0 || k==30)
       //{I_MDFA.plotData(x,localsignal,x.length);}
       start = n_total-actual_n;
       for(i=flength;i<n_obs;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*old_tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
        }
        hist_sig[i-flength] = sum;
       }
       insampleTradingDiff(hist_price, hist_sig, flength); 
       rankCoefficient(account,flength);
       insamp_stat[k] = rank_coeff;
 
       start = n_total - trade_length;
       for(i=start;i<start+start_trades;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
        }
        x_morn[i-start] = yesterday_series[i];
        short_sig[i-start] = sum;
        //System.out.println(short_sig[i-start]);
       } 
       

       
//        if(k==15)
//        {
//         
//         System.out.println("H0 = " + b_h0[L] + " " + b_h0[L+L] + " " + b_h0[L+L+L]);
//         for(i=0;i<actual_n;i++)
//         {
//           System.out.println(Math.exp(logprice[i]) + " " + localsignal[i]);        
//         }
//        
//        }
              
       
       
       
       //compute estimated optimal parameter
       insampleTradingDiff(short_price, short_sig, start_trades);       
       rankCoefficient(pnl,start_trades);       
       opt_ret = account[account.length-1];
       opt_rank_coeff = rank_coeff;
       interp_max[k] = opt_ret; //*Math.abs(opt_rank_coeff);
       //rankCoefficient(pnl,start_trades); 
       full_rank = rank_coeff;
       //System.out.println("opt_rank_coeff = " + opt_rank_coeff);
//        I_MDFA.plotData(pnl, pnl, pnl.length);
//        I_MDFA.plotData(x_morn, short_sig, x_morn.length);
       if(k==0) {zero_ret = opt_ret;}//*Math.abs(opt_rank_coeff);}
       if(k==n_interps-1) {one_ret = opt_ret;}//*Math.abs(opt_rank_coeff);}
       
       //---- compute the the real optimal parameter
       insampleTradingDiff(logprice, localsignal, actual_n);
       
//        if(k==0)
//        {  
//          System.out.println("At optimization k=0");
//          for(i=0;i<actual_n;i++)
//          {System.out.println(logprice[i] + " " + localsignal[i]);}
//        }
       
       
       //rankCoefficient(pnl,actual_n);
       //full_rank = rank_coeff;
       full_return = account[account.length - 1]; 
       total_return = full_return;
       
       System.out.println("max_ret = " + opt_rank_coeff + " at " + interpolation + ", return = " + total_return + ", stand_dev = " + standard_dev);
       
       //if(total_return > max_ret) 
       if(opt_ret >= max_rank)
       //if(k==0)
       {
         max_rank = opt_rank_coeff; max_ret = opt_ret; max_interp = interpolation;
         System.out.println("NEW MAX FOUND: max_ret = " + opt_rank_coeff + " at " + max_interp + ", return = " + total_return);         
         //best_account = new double[trade_length]; System.arraycopy(account,0,best_account,0,trade_length);
         //best_signal = new double[trade_length]; System.arraycopy(localsignal,0,best_signal,0,trade_length); 
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;           
       
       
       //total_return = account[trade_length-1];
       //I_MDFA.plotData(account, account, account.length);
       

       //I_MDFA.plotData(account, account, account.length-1);
       //rankCoefficient(pnl, pnl.length);
              
       //if(full_return >= 0)
      // {histo_stat[count] = histo_stat[count] + 1;}

       if(full_return > 0) {num_postive++;}
       
       full_returns[count] = full_return;
       
       if(full_return > 0) {num_full_positive_returns++;}
       
       if(full_return > full_return_max)
       //if(full_rank >= full_rank_max)
       {
          full_rank_max = full_rank;
          full_return_max = full_return;
          full_return_int = interpolation;
          System.out.println("REAL MAX FOUND at " + full_return_int + ", rank = " +  ",return = " + full_return);
    
       }
     
     
       //if(interpolation == 0)
       //{crit_0.add(mdfa.criteria); deg_0.add(mdfa.diff_band);}
       //else
       //{crit_1.add(mdfa.criteria); deg_1.add(mdfa.diff_band);}
       
       //interpolation = interpolation + delta;
       count++;      
     }
     
     lookback_returns.add(full_return_max);
     if(num_postive>0)
     {num_pos_returns++;}
    
     full_returns_array.add(full_returns);
     morning_returns.add(insamp_stat);   
     
     avg_vol = 0;
     if(avg_volatility.size() > 5)
     {avg_vol = (avg_volatility.get(avg_volatility.size()-1) + avg_volatility.get(avg_volatility.size()-2))/2.0;}
     double current_vol = avg_volatility.get(avg_volatility.size()-1);
     //svm.add(full_return_int + " " + MDFAmin + " " + criteria + " " + max_peak + " " + morn_vol + " " + zero_ret + " " + one_ret + " " + avg_vol);     
     
     prev_return = 0; 
     if(returns.size() > 2)
     {prev_return = returns.get(returns.size() - 2);}
     
     svm.add(crite_0 + " " + peak_loc + " " + max_peak + " " + morn_vol + " " + avg_vol + " " + zero_ret + " " + one_ret + " " + current_vol + " " + price_move + " " + b_h0[L] + " " + full_return_int + " " + prev_return);
     
//      out = computeInterpolation(old_tseries,full_return_int);
//      
// //       out = mdfa.H0_IMDFAreg(old_tseries, n, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
// //                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, full_return_int, mdfa.cross, 
// //                        mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
// //                        mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);    
//      
//       System.arraycopy(out,lag4,b_h0,0,n_rep*L);
//        
//        //System.arraycopy(mdfa.b,0,b_h0,0,mdfa.b.length);
//        
//        //apply on time series
//        start = n_total - trade_length;
//        for(i=start;i<start+start_trades;i++)
//        {
//         sum = 0.0;
//         for(j=1;j<n_rep;j++)
//         {
//            for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
//         }
//         best_signal[i-start] = sum;
//        }
//      
//        insampleTradingDiff(short_price, best_signal, start_trades);
//        full_rank_max = account[account.length-1];
       //rankCoefficient(pnl,start_trades);     
       //full_rank_max = rank_coeff;
     
//      System.out.println(b_h0[L+0] + " " + b_h0[L+L] + " " + b_h0[L+2*L]);
//      System.out.println(b_h0[L+L-1] + " " + b_h0[L+L+L-1] + " " + b_h0[L+2*L+L-1]);
//          
     
//      for(i=0;i<12;i++)
//      {System.out.print(df4.format(best_signal[i]) + " ");}
//      System.out.println(df4.format(best_signal[12]));

//      for(i=0;i<trade_length;i++)
//      {System.out.println(best_signal[i] + " " + best_account[i]);}
 
     /*for(i=0;i<12;i++)
     {System.out.print(df4.format(x[i]) + " ");}
     System.out.println(df4.format(x[12]));   */   

     
/*     for(i=0;i<actual_n;i++)
     {
       actual_price[actual_n-1-i] = logprice[logprice.length-1-i];
       actual_sig[actual_n-1-i] = best_signal[best_signal.length - 1 - i];
       actual_best[actual_n-1-i] = localsignal[localsignal.length-1-i];
       //localsignal[localsignal.length - 1 - i] = best_signal[best_signal.length - 1 - i];
     }    
     //if(best_account[9] == 0) {morning_buy = true;}
     //else {morning_buy = false;}     
     
     for(i=0;i<actual_price.length;i++)
     {
       System.out.println(actual_sig[i] + " " + actual_best[i] + " " + Math.exp(actual_price[i])); 
     } */    
     
     //morning_buy = true;
// 
     //insampleTradingDiff(logprice, localsignal, trade_length);
//      //System.out.println(account[trade_length-1]);
//      insampleTradingDiff(actual_price, actual_best, actual_n);
//      full_return_max = account[account.length-1];
//      
//      insampleTradingDiff(actual_price, actual_sig, actual_n);
//      System.out.println("Total = " + account[account.length-1] + ", " + account.length);
     
     //I_MDFA.plotData(account,account,account.length);
     
     
     //crits.add(max_interp_value + " " + account[account.length-1] + " " + full_return_int  + " " + full_return_max  + " " + max_ret + " " + best_account[start_trades] + " " + (best_account[trade_length-1] - best_account[start_trades]));     
//      if(max_ret_interp > 0.013)
//      {
//      
//        out = mdfa.H0_IMDFAreg(in_sample, n_obs, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
//                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, localdecay, max_interp_value, mdfa.cross, 
//                        mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
//                        mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
// 
//        //System.out.println("interpolation value = " + interpolation);
//        System.arraycopy(out,lag4,b_h0,0,n_rep*L);     
//        Filter f = new Filter(L,n_rep);
//        f.setCoeffs(b_h0, L, n_rep);
//        filters.add(f);
//      
//      }
     
//      if(max_interp == 0)
//      {max_interp = .999;}
//      else
//      {max_interp = 0;}
     
     //------ Get next interp value based on past values 
//      t2 = t1;
//      t1 = max_interp_value;
//      double interp_val = getNextIV(t1,t2);
     //-------------------------------------------------
     
//      interp.add(max_interp);
//      p0 = interp.get(interp.size()-1);
//      p1 = interp.get(interp.size()-2);
//      p2 = interp.get(interp.size()-3);
//      if(p0 > 0) p0 = 1; 
//      if(p1 > 0) p1 = 1; 
//      if(p2 > 0) p2 = 1;
//      
//      
//      if((p0 == 0) && (p1 == 0) && (p2 == 0)) {interp_val = 1;}
//      else if((p0 == 0) && (p1 == 0) && (p2 == 1)) {interp_val = 1;}
//      else if((p0 == 0) && (p1 == 1) && (p2 == 0)) {interp_val = 0;}
//      else if((p0 == 0) && (p1 == 1) && (p2 == 1)) {interp_val = 1;}
//      else if((p0 == 1) && (p1 == 0) && (p2 == 0)) {interp_val = 1;}
//      else if((p0 == 1) && (p1 == 0) && (p2 == 1)) {interp_val = 0;}
//      else if((p0 == 1) && (p1 == 1) && (p2 == 0)) {interp_val = 0;}
//      else if((p0 == 1) && (p1 == 1) && (p2 == 1)) {interp_val = 1;}
//      else
//      {interp_val = max_interp;}
     
     //interp_val = forecastInterpolation(1);
     //interp_val = forecastInterpolationDFA();
     //if(losses_in_arow >= 2)
     //{interp_val = .00;}
     
//      System.out.println("Computing filter with new data"); 
//      for(i=0;i<n_obs;i++)
//      {
//        System.out.println(tseries[i] + " " + tseries[n_obs + i] + " " + tseries[2*n_obs + i] + " " + tseries[3*n_obs + i]);
//      }
/*     Random generator = new Random();
     double r = generator.nextGaussian();
     if(Math.abs(r) < .40) {r = r/2;}
     else if(r > 1.0) {r = .90;}
     else if(r < 1.0) {r = -.90;}
//      u = Math.random()/1.001;
//      xo = Math.log(1.0-u)/(-6.0);
xo = full_return_int + r; 
     if(xo < 0) {xo = 0;}
     else if(xo > .999) {xo = .999 - Math.random()/10.0;}*/ 
     crits.add(max_interp_value + " " + max_ret_interp + " " + full_return_int  + " " + full_return_max  +  " " + full_rank_max + " " + MDFAmin + " " + criteria + " " + max_peak + " " + morn_vol);     
      //---- First get current filter, use 0 for simplicity
     
     
     
     if(returns.size() > 0) 
     {
      if(returns.get(returns.size() - 1) < 0.0) {}
     }
     
     u = Math.random(); if(u<.50) {} else{} 
     
     if(full_return_int == .7)
     {}
     else
     { 
      if(u > .60) {} else {}
     }
     //if(k==0) {interpol = interpolations[k_opt+1];}
     //else if(k==1) {interpol = interpolations[k_opt+1];}
     //else if(k==2) {interpol = interpolations[k_opt-add];}
     System.out.println("Used interpolation = " + u);
     //out = computeInterpolationFrequency(tseries,interpolations[15]);
     out = computeInterpolation(tseries, .99999);
     
     
//      out = mdfa.H0_IMDFAreg(tseries, n, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
//                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, full_return_int, mdfa.cross, 
//                        mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
//                        mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     
     interp_vals.add(full_return_int);
     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }
     
     //System.out.println("New BO = " + b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[L+L]);     
     
     
     b_copy = new double[L*n_rep];
     System.arraycopy(out,lag4,b_copy,0,n_rep*L);      
     
     
     old_tseries = new double[tseries.length];
     System.arraycopy(tseries,0,old_tseries,0,tseries.length);
     
     

     
       
   }    
   
   //--- compute the H0 filter with respect to the given filter
   public void computeH0Filter(double[] series)
   {
     int l,i,lag4,K1;
     double[] out;
     
     System.out.println("Recomputing H0 Filter");
     mdfa_h0.set_tseries(series,n_obs,n_rep);
     
     out = mdfa_h0.GEN_IMDFAreg(series, n_obs, mdfa_h0.L, n_rep, 0, 0, mdfa_h0.Gamma, 0, mdfa_h0.lambda, mdfa_h0.expweight, 0, mdfa_h0.cutoff0,
                       mdfa_h0.cutoff, mdfa_h0.Lag, mdfa_h0.i1, mdfa_h0.i2, mdfa_h0.w_const, mdfa_h0.smooth, 0, 0, mdfa_h0.cross, mdfa_h0.Q_cdev, mdfa_h0.iter, mdfa_h0.spec_dens, mdfa_h0.spec,
                       mdfa_h0.ar_p, mdfa_h0.ma_q, mdfa_h0.ar_params, mdfa_h0.ma_params, mdfa_h0.innvar,mdfa_h0.arma);
  
     h0b0 = new double[(n_rep-1)*L];
     K1 = (int)n_obs/2+1;
     lag4 = 2*(n_obs-L+1) + 2*(n_rep+1)*K1;

     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {h0b0[L*i + l] = out[lag4 + L*(i+1)+l];}
     } 
     
     H0set = true;
     criteria_h0 = out[out.length-3];
     degrees_h0 = out[out.length-2];      
     MDFAmin_h0 = out[out.length-1];     
 
   }
   

   
   public double[] computeInterpolation(double[] series, double interp)
   {
     
      double[] out;
      if(interp > 0 && interp < 1)
      {
      
       out = mdfa.H0_IMDFAreg(series, n_obs, L, n_rep, 0, 1, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, mdfa.Lag, mdfa.i1, mdfa.i2, mdfa.w_const, mdfa.smooth, 0, interp, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);    
      }
      else if(interp == 0)
      {
         System.out.println(series.length);
         System.out.println(mdfa.lambda);
         System.out.println("length = " + series.length + "nobs = " + mdfa.n_obs + ", lambda = " + mdfa.lambda + ", tseries[0] = " + series[n_obs-2] + ", n_rep = " + mdfa.n_rep);          
         //for(i=0;i<n_obs;i++) {System.out.println(tseries[i]);}
         System.out.println("cut0 = " + mdfa.cutoff0 + ", cut = " + mdfa.cutoff);
         System.out.println("length Gamma = " + Gamma.length + ", expweight = " + mdfa.expweight + ", regs = " + mdfa.smooth + " " + mdfa.decay + " " + mdfa.decay2 + " " + mdfa.cross);
         System.out.println("L = " + mdfa.L + ", lag = " + mdfa.Lag + ", i1 = " + mdfa.i1 + ",i2 = " + mdfa.i2);        
               
      
      
        out = mdfa.GEN_IMDFAreg(series, n_obs, L, n_rep, 0, 1, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, mdfa.Lag, mdfa.i1, mdfa.i2, mdfa.w_const, mdfa.smooth, mdfa.decay, mdfa.decay2, mdfa.cross, mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma);     
      }
      else 
      {
        if(H0set)
        {
        out = mdfa_h0.GEN_IMDFAreg(series, n_obs, mdfa_h0.L, n_rep, 0, 0, mdfa_h0.Gamma, 0, mdfa_h0.lambda, mdfa_h0.expweight, 0, mdfa_h0.cutoff0,
                       mdfa_h0.cutoff, mdfa_h0.Lag, mdfa_h0.i1, mdfa_h0.i2, mdfa_h0.w_const, mdfa_h0.smooth, 0, 0, mdfa_h0.cross, mdfa_h0.Q_cdev, mdfa_h0.iter, mdfa_h0.spec_dens, mdfa_h0.spec,
                       mdfa_h0.ar_p, mdfa_h0.ma_q, mdfa_h0.ar_params, mdfa_h0.ma_params, mdfa_h0.innvar,mdfa_h0.arma);   
        }
        else
        {out = new double[040];}
      }
      return out; 
   }
  
  
   public double[] computeInterpolationFrequency(double[] series, double interp)
   {
     
      double[] out; int k; double om;
      double lcutoff1;
      
      lcutoff1 = .16*(1.0 - interp) + .64*interp;
      System.out.println("Frequency = " + lcutoff1);  
      double[] _Gamma = new double[K1];
      for(k=0; k<=K;k++)
      {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {_Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= lcutoff1)
         {_Gamma[k] = 1.0;}
         else
         {_Gamma[k] = 0.0;}
      }         
      
      out = mdfa.GEN_IMDFAreg(series, n_obs, L, n_rep, 0, 1, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       lcutoff1, mdfa.Lag, mdfa.i1, mdfa.i2, mdfa.w_const, mdfa.smooth, mdfa.decay, mdfa.decay2, mdfa.cross, mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma);     

      return out; 
   }  
  
  
   public double[] computeInterpolationSmooth(double[] series, double interp)
   {
     
      double[] out;  
      double smooth;
      
      smooth = 10.0*(1.0-interp) + 35.0*interp;

      out = mdfa.GEN_IMDFAreg(series, n_obs, L, n_rep, 0, 1, Gamma, 0, mdfa.lambda, smooth, 0, cutoff0,
                       cutoff, mdfa.Lag, mdfa.i1, mdfa.i2, mdfa.w_const, mdfa.smooth, mdfa.decay, mdfa.decay2, mdfa.cross, mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                      mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma);     

      return out; 
   }   
  
  
  public static double regularity(Double[] returns)
  {
    
     int n_changes = 0;
     
     for(int i = 0; i < returns.length-1; i++)
     {
       if((returns[i+1] > 0 && returns[i] < 0) || (returns[i+1] < 0 && returns[i] > 0)) {n_changes++;}
     }   
     return (double)(returns.length - n_changes)/(double)returns.length;
  }
  
  
  
  
  public static double[] ComputeACF(Double[] data, int maxLag)
  {
   int i, j, n; 
   double val,val2,tx;
   
   n = data.length; 
   double total = 0;
   double[] acf = new double[maxLag+1];
   double[] std = mean_std(data); 
  
   
 
   for (i = 0; i <= maxLag; i++)
   {      
   
      for (j = i; j < n; j++)
      {
	 val = data[j] - std[0];   
	 val2 = data[j-i] - std[0];
	 tx = val*val2;
	 total += tx;       
      }
      //---- copy to acf and set to zero--------
      
      acf[i] = total;
      total = 0.0;      
   }
   
   for(i=0;i<=maxLag;i++) {acf[i] = acf[i]/acf[0];}
   
   return acf; 
   
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
 

  public static double mean_std_avg( Double[] data, int b ) 
  { 

       double mean = 0; 
       int buckets = b;
       final int n = (int)data.length/buckets; 
       double[] bucket_vols = new double[buckets];
       double su = 0; int k;
       
       for(k = 0; k < buckets; k++)
       {
       
        for ( int i=0; i<n; i++ )  {  mean += data[i + k*b]; } 
        mean /= n; 

        double sum = 0; 
        for ( int i=0; i<n; i++ ) { final double v = data[i + k*b] - mean;  sum += v * v; } 
        double[] co = {mean,  Math.sqrt( sum / n )};
        
        bucket_vols[k] = co[1];       
       }  
       
       for(k = 1; k < buckets; k++)
       {su = su + Math.abs(bucket_vols[k] - bucket_vols[k-1]);}
        
       return su; 
  }     
  
  
  
   
   public int forecastOptimize()                    
   {
   
      int N = n_obs; int l,j,i;
      int max_lag = 0; 
      double max_ret = -1.0;
      double max_rank = -11.0;
      double sum=0;
      int lag;
      double ret;
      
      for(lag = -4; lag < -3 ; lag++)
      {
   
       mdfa.set_lag(lag);
       mdfa.computeFilterGeneral(true,false);              
       b_coeffs = new double[(n_rep-1)*L];
      
       for(l=0;l<L;l++)
       {
        for(i=0;i<n_rep-1;i++) {b_coeffs[L*i + l] = mdfa.b[L*(i+1)+l];}
       }
   
       int flength = n_obs - L + 1;
       double[] local_price = new double[flength];
       double[] local_sig = new double[flength];
      
      
       for(i=0;i<flength;i++)
       {local_price[flength - i - 1] = price.get(price.size()-i-1).doubleValue();}
 
       for(i=L-1;i<n_obs;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
         for(l=0;l<L;l++) {sum = sum + b_coeffs[L*(j-1) + l]*tseries[N*j + i - l];}
        } 
        local_sig[i - L + 1] = sum;
       }
 
       insampleTradingDiff(local_price, local_sig, local_sig.length);
       ret = account[account.length-1];
       rankCoefficient(account,account.length);
              
       //if(rank_coeff > max_rank)
       if(max_ret < ret)
       {max_lag = lag; max_rank = rank_coeff; max_ret = ret; System.out.println("Rank coeff = " + max_rank + ", at lag = " + lag);}
   
      }
      
      return max_lag;
   }
   
   public void H0_Optimization_Morning(double localdecay, double delta, double interp_start)
   {
   
     int i,j,k,l;
     int trade_length = 12; 
     double interpolation = interp_start;
     int K,K1;
     
   
     double max_ret = -100.0;  
     double max_interp = 0.0;
     int n = n_obs;
     int n_total = n_obs + trade_length;
     double[] yesterday_series = new double[n_rep*n_total];
     double[] in_sample = new double[n_rep*n_obs];
     int shift = 1;
     
     //--- create yesterday_series
     for(i=0;i<n_total;i++)
     {
       yesterday_series[n_total-1-i] = gauss_target.get(gauss_target.size() - shift - i);
       yesterday_series[n_total + n_total-1-i] = gauss_target.get(gauss_target.size() - shift - i);
       yesterday_series[n_total*2 + n_total-1-i] = gauss_1.get(gauss_1.size() - shift - i);
       yesterday_series[n_total*3 + n_total-1-i] = gauss_2.get(gauss_2.size() - shift - i);         
     }
          
     for(i=0;i<n_obs;i++)
     {
       in_sample[i] = yesterday_series[i];
       in_sample[n_obs   + i] = yesterday_series[n_total + i];
       in_sample[n_obs*2 + i] = yesterday_series[n_total*2 + i];
       in_sample[n_obs*3 + i] = yesterday_series[n_total*3 + i];
     }     
     
     
     //--- create new Gamma
     double sum; double total_return;
     K = (int)(n/2); double om; K1 = K+1;  
     double[] _Gamma = new double[K1];       
     for(k=0; k<=K;k++)
     {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {_Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {_Gamma[k] = 1.0;}
         else
         {_Gamma[k] = 0.0;}
     }         
     
     int start; 
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] out;
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-1-i] = price.get(price.size()-1-i).doubleValue();}
          
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
  
     while(interpolation <= .9999)
     {
     
       
       out = mdfa.H0_IMDFAreg(in_sample, n_obs, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, localdecay, interpolation, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     

       //System.out.println("interpolation value = " + interpolation);
       System.arraycopy(out,lag4,b_h0,0,n_rep*L);
       
       //apply on time series
       start = n_total-trade_length;
       for(i=start;i<n_total;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];}
        }
        x[i-start] = yesterday_series[i];
        localsignal[i-start] = sum;
       } 
     
       //if(interpolation == 0.0)
       //{I_MDFA.plotData(x,localsignal,trade_length);}
       
       //---- compute the insample-outsample trading
       insampleTradingDiff(logprice, localsignal, trade_length);
       total_return = account[trade_length-1];
       
       //I_MDFA.plotData(pnl, account, account.length-1);
       rankCoefficient(pnl, pnl.length);
       System.out.println("max_rank = " + rank_coeff + " at " + interpolation + ", return = " + total_return);
       
       
       if(total_return > max_ret) 
       {
         max_ret = total_return; max_interp = interpolation;
         max_interp = interpolation;
         System.out.println("NEW MAX FOUND: max_rank = " + rank_coeff + " at " + max_interp + ", return = " + total_return); 
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;
     
       interpolation = interpolation + delta;      
     }
     
 
 
      //---- First get current filter, use 0 for simplicity
     out = mdfa.H0_IMDFAreg(tseries, n, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, 0, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }  
  
   }       
   
   
   
   public void H0_Optimization_Lookback2(double delta, double interp_start)
   {
   
     int i,j,k,l;
     int trade_length = 13; double localdecay = 0.0; 
     double interpolation = interp_start;
     int K,K1;
     
   
     double max_ret = -100; double max_rank = -1.0; 
     double max_interp = 0.0;
     int n = n_obs;
     int n_total = n_obs + trade_length;
     double[] yesterday_series = new double[n_rep*n_total];
     double[] in_sample = new double[n_rep*n_obs];
     
     int shift = 14;
     //--- create yesterday_series
     for(i=0;i<n_total;i++)
     {
       yesterday_series[n_total-1-i] = gauss_target.get(gauss_target.size() - shift - i);
       yesterday_series[n_total + n_total-1-i] = gauss_target.get(gauss_target.size() - shift - i);
       yesterday_series[n_total*2 + n_total-1-i] = gauss_1.get(gauss_1.size() - shift - i);
       yesterday_series[n_total*3 + n_total-1-i] = gauss_2.get(gauss_2.size() - shift - i);             
     }
          
     for(i=0;i<n_obs;i++)
     {
       in_sample[i] = yesterday_series[i];
       in_sample[n_obs   + i] = yesterday_series[n_total + i];
       in_sample[n_obs*2 + i] = yesterday_series[n_total*2 + i];
       in_sample[n_obs*3 + i] = yesterday_series[n_total*3 + i];
     }     
     
     
     //--- create new Gamma
     double sum; double total_return;
     K = (int)(n/2); double om; K1 = K+1;  
     double[] _Gamma = new double[K1];       
     for(k=0; k<=K;k++)
     {       
         om = (k*Math.PI/K);
         if(om < cutoff0) {_Gamma[k] = 0.0;}
         else if(om >= cutoff0 && om <= cutoff1)
         {_Gamma[k] = 1.0;}
         else
         {_Gamma[k] = 0.0;}
     }         
     
     int start; 
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] out;
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-1-i] = price.get(price.size()-14-i).doubleValue();}
          
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
  
     while(interpolation <= .9999)
     {
 
       out = mdfa.H0_IMDFAreg(in_sample, n_obs, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, localdecay, interpolation, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     

       //System.out.println("interpolation value = " + interpolation);
       System.arraycopy(out,lag4,b_h0,0,n_rep*L);
       
       //apply on time series
       start = n_obs;
       for(i=start;i<n_total;i++)
       {
        sum = 0.0;
        for(j=1;j<n_rep;j++)
        {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];}
        }
        x[i-start] = in_sample[i];
        localsignal[i-start] = sum;
       } 
     
       //if(interpolation == 0.0)
       //{I_MDFA.plotData(x,localsignal,flength);}
       
       //---- compute the insample-outsample trading
       insampleTradingDiff(logprice, localsignal, trade_length);
       total_return = account[trade_length-1];
       
       //I_MDFA.plotData(account, account, account.length-1);
       rankCoefficient(pnl, pnl.length-1);
       System.out.println("max_rank = " + rank_coeff + " at " + interpolation + ", return = " + total_return);
       
       
       if(rank_coeff > max_rank) 
       {
         max_interp = interpolation; max_rank = rank_coeff;
         System.out.println("NEW MAX FOUND: max_rank = " + rank_coeff + " at " + max_interp + ", return = " + total_return);  
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;
     
       interpolation = interpolation + delta;      
     }
     
     if(max_interp == 0) {max_interp = .999;}
     else {max_interp = 0;}
 
 
      //---- First get current filter, use 0 for simplicity
     out = mdfa.H0_IMDFAreg(tseries, n, L, n_rep, 0, 0, _Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, max_interp, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }  
 
 
 
        
   }       
   
   
   
   
   public void H0_Optimization_Lookback_Test(double[] yesterday_series, double delta, double interp_start, double optval, int cut)
   {
   
     int i,j,l;
     int trade_length = trade_obs;  
     int K,K1;

     double max_rank = -1.0; 
     double max_interp = 0.0;
     int n = n_obs;
     int n_total = n_obs + trade_length-1;
     //double[] yesterday_series = new double[n_rep*n_total];
     double[] in_sample = new double[n_rep*n_obs];
     int shift = 6;
     
     
     
     for(i=0;i<n_obs;i++)
     {
       in_sample[i] = yesterday_series[i];
       in_sample[n_obs+i] = yesterday_series[n_total + i];
       in_sample[2*n_obs+i] = yesterday_series[2*n_total + i];
       in_sample[3*n_obs+i] = yesterday_series[3*n_total + i];
     }     
     
     
   
   
   
   
     //--- create new Gamma
     double sum; double total_return, total_return2;
     K = (int)(n/2); K1 = K+1;  
        
     
     int start; 
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] b_h1 = new double[n_rep*L];
     double[] out0, out1;
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-1-i] = price.get(price.size()-shift-i).doubleValue();}
          
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
  
     out0 = mdfa.H0_IMDFAreg(in_sample, n_obs, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, 0, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     System.arraycopy(out0,lag4,b_h0,0,n_rep*L);
     
     out1 = mdfa.H0_IMDFAreg(in_sample, n_obs, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, .999, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                       mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     System.arraycopy(out1,lag4,b_h1,0,n_rep*L);     
     
       
       //apply on time series
     start = n_total-trade_length;
     for(i=start;i<n_total-cut;i++)
     {
       sum = 0.0;
       for(j=1;j<n_rep;j++)
       {
         for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*yesterday_series[n_total*j + i-l];}
       }
       x[i-start] = yesterday_series[i];
       localsignal[i-start] = sum;
     }  
     
     insampleTradingDiff(logprice, localsignal, trade_length-cut);
     
     System.out.println("Morning price: ");
     for(i=0;i<trade_length-cut-1;i++) {System.out.print(df4.format(logprice[i]) + " ");}
     System.out.println(""+df4.format(logprice[trade_length-cut-1]));      
     
     System.out.println("Morning signal: ");
     for(i=0;i<trade_length-cut-1;i++) {System.out.print(df4.format(localsignal[i]) + " ");}
     System.out.println(""+df4.format(localsignal[trade_length-cut-1]));     
     
     System.out.println("Morning returns: ");
     for(i=0;i<account.length-1;i++) {System.out.print(df4.format(account[i]) + " ");}
     System.out.println(""+df4.format(account[account.length-1]));
     

     
     total_return = account[trade_length-1-cut];
     rankCoefficient(pnl, trade_length-cut);      
     double rank1 = rank_coeff;
     //System.out.println("Total return at this cut " + total_return);
     zero_ret = final_trade;
   
     for(i=start;i<n_total-cut;i++)
     {
       sum = 0.0;
       for(j=1;j<n_rep;j++)
       {
         for(l=0;l<L;l++) {sum = sum + b_h1[L*j + l]*yesterday_series[n_total*j + i-l];}
       }
       x[i-start] = yesterday_series[i];
       localsignal[i-start] = sum;
     }  
     
     insampleTradingDiff(logprice, localsignal, trade_length-cut);
  
     System.out.println("Morning signal: ");
     for(i=0;i<trade_length-cut-1;i++) {System.out.print(df4.format(localsignal[i]) + " ");}
     System.out.println(""+df4.format(localsignal[trade_length-cut-1]));     
  
     System.out.println("Morning returns: ");
     for(i=0;i<account.length-1;i++) {System.out.print(df4.format(account[i]) + " ");}
     System.out.println(""+df4.format(account[account.length-1]));     
     
     total_return2 = account[trade_length-1-cut];
     rankCoefficient(pnl, trade_length-cut);      
     double rank2 = rank_coeff;   
   
     if(final_trade > zero_ret)
     {zero_ret = final_trade;}
     
   
     //if(total_return >= total_return2)
     if(rank1 >= rank2)
     {
       max_rank = rank1;
       max_interp = 0;

       max_interp_value = max_interp;
       max_ret_interp = total_return;
     }
     else
     {
       max_rank = rank2;
       max_interp = .999;

       max_interp_value = max_interp;
       max_ret_interp = total_return2;
     }        
       
     System.out.println("max_rank = " + max_rank + " " + max_interp + " " + optval); 
   
 
      //---- First get current filter, use 0 for simplicity
     out1 = mdfa.H0_IMDFAreg(tseries, n, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, max_interp_value, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out1[lag4 + L*(i+1)+l];}
     }  
     
      
     
     criteria = out1[out1.length-3];
     degrees = out1[out1.length-2];      
     MDFAmin = out1[out1.length-1];
     
     crits.add(""+criteria + " " + degrees + " " + MDFAmin);
    
   }       
   
   
   public double computeVarRatio(ArrayList<Double> recent, ArrayList<Double> total)
   {
   
      int i;
      
      double sum = 0;
      double meanv0 = 0;
      double meanv0base = 0;
      double var0;
      
      Double[] v0 = recent.toArray(new Double[0]);
      Double[] vbase = total.toArray(new Double[0]);
      
      double[] logRetV0 = new double[v0.length];
      double[] logRetVbase = new double[vbase.length];
      
      //System.out.println("logRetV0 length = " + logRetV0.length + ", " + logRetVbase.length);
      
      logRetV0[0] = 0; logRetVbase[0] = 0;
      
      for(i = 1; i < v0.length; i++)
      {
        logRetV0[i] = Math.log(v0[i]) - Math.log(v0[i-1]);
        meanv0 = meanv0 + logRetV0[i];
      }
      meanv0 = meanv0/(1.0*v0.length);
      
      
      for(i = 1; i < vbase.length; i++)
      {
        logRetVbase[i] = Math.log(vbase[i]) - Math.log(vbase[i-1]);
        meanv0base = meanv0base + logRetVbase[i];
      }      
      meanv0base = meanv0base/(1.0*logRetVbase.length);
    
      //System.out.println("means = " + meanv0 + ", " + meanv0base);
      //System.out.println("last vals = " + logRetV0[logRetV0.length-1] + ", " + logRetVbase[logRetVbase.length-1]);
    
      sum = 0;
      for(i=0; i < v0.length; i++)
      {sum = sum + (logRetV0[i] - meanv0)*(logRetV0[i] - meanv0);}     
      var0 = Math.sqrt(sum)/(1.0*v0.length);    
    
      sum = 0;
      for(i=0; i < logRetVbase.length; i++)
      {sum = sum + (logRetVbase[i] - meanv0base)*(logRetVbase[i] - meanv0base);}     

      
      return (var0);
      //return var1;
    
   }

   
   public void H0_Optimization_Intraday(double delta, double interp_start, int advance)
   {
   
     int i,j,l;
     int trade_length = advance;  
     double interpolation = interp_start;
     int K,K1;
     double max_ret = -100; double max_rank = -1.0; 
     double max_interp = 0.0;
     
     int n = n_obs;
     int shift = 1;
     
     //--- create new Gamma
     double sum; double total_return;
     K = (int)(n/2); K1 = K+1;  
         
     
     int start; 
     
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] out;
     
     
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-i-1] = price.get(price.size()-shift-i).doubleValue();}
     
     
//      for(i=0;i<12;i++)
//      {System.out.print(df4.format(logprice[i]) + " ");}
//      System.out.println(df4.format(logprice[12]));       
//  
//      for(i=0;i<12;i++)
//      {System.out.print(df4.format(morningSeries[n_obs - 1 - i]) + " ");}
//      System.out.print(df4.format(morningSeries[n_obs - 1 - 12]));     
//      System.out.println("series lenght = " + morningSeries.length);
     
     
//      System.out.println("09:30 data at revision: " + in_sample[n_obs-1] + " " + in_sample[n_obs + n_obs-1] + " " + in_sample[n_obs*2 + n_obs-1] + " " + in_sample[n_obs*3 + n_obs-1] + " " 
//       + logprice[0]);
                 
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
  
  
      while(interpolation <= .9999)
      {
//      
//        //mdfa.setRegularization(mdfa.smooth, localdecay, interpolation, mdfa.cross);
//        //mdfa.computeFilterGeneral(true);
//      
        out = mdfa.H0_IMDFAreg(morningSeries, n_obs, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, interpolation, mdfa.cross, 
                        mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                        mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
 
 
 
//        //System.out.println("interpolation value = " + interpolation);
        System.arraycopy(out,lag4,b_h0,0,n_rep*L);
//        
//        //System.arraycopy(mdfa.b,0,b_h0,0,mdfa.b.length);
//        
//        //apply on time series
         start = n_obs-trade_length;
         for(i=start;i<n_obs;i++)
         {
          sum = 0.0;
          for(j=1;j<n_rep;j++)
          {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
          }
          x[i-start] = tseries[i];
          localsignal[i-start] = sum;
         }

//        //{I_MDFA.plotData(x,localsignal,trade_length);}
//        
       //---- compute the insample-outsample trading
       insampleTradingDiff(logprice, localsignal, trade_length);
       total_return = account[trade_length-1];
         
       //I_MDFA.plotData(account, account, account.length-1);
       rankCoefficient(pnl, pnl.length);
       System.out.println("max_rank = " + rank_coeff + " at " + interpolation + ", return = " + total_return);
       
       //if(total_return > max_ret) 
       if(rank_coeff > max_rank)
       {
         max_rank = rank_coeff; max_ret = total_return; max_interp = interpolation;
         System.out.println("NEW MAX FOUND: max_rank = " + rank_coeff + " at " + max_interp + ", return = " + total_return);
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;
       
       interpolation = interpolation + intra_delta;   
     }  
  
     //System.out.println("pnl length = " + pnl.length);
    
  
      //---- First get current filter, use 0 for simplicity
     out = mdfa.H0_IMDFAreg(morningSeries, n, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, max_interp, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);    
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }  
     
     System.out.println(b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[2*L]);
     System.out.println(b_coeffs[L-1] + " " + b_coeffs[L+L-1] + " " + b_coeffs[2*L+L-1]);

      System.arraycopy(out,lag4,b_h0,0,n_rep*L);
      start = n_obs-trade_length;
      for(i=start;i<n_obs;i++)
      {
          sum = 0.0;
          for(j=1;j<n_rep;j++)
          {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
           //System.out.println("sum = " + sum);
          }
          x[i-start] = tseries[i];
          localsignal[i-start] = sum;
     }    
     
     for(i=0;i<localsignal.length-1;i++)
     {System.out.print(df4.format(localsignal[i]) + " ");}
     System.out.println(df4.format(localsignal[localsignal.length-1]));
 
     for(i=0;i<x.length-1;i++)
     {System.out.print(df4.format(x[i]) + " ");}
     System.out.println(df4.format(x[x.length-1])); 
     
     
     interp_vals.add(max_interp); max_ranks.add(max_rank);
     System.out.println("Interp value choses = " + max_interp);
        
   }         
   
   
   public void H0_Optimization_Update(double delta, int advance, boolean new_day, boolean new_signal)
   {

     int i,j,k,l;
     
     //--- trade length is 8:30 - 9:30 plus additional times 
     int trade_length = 1 + advance;  
     double interpolation = interp_start;
     int K,K1;
     int[] order; double[] morn_rets;
     double max_ret = 100;  
     double max_interp = 0.0;   
   
     int n = n_obs;
     int shift = 1;   

     double sum; double total_return;
     K = (int)(n/2); K1 = K+1;  
         
     
     int start; 
     
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] out;     
     
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     
     //if(recomputeH0 && mdfa_h0 != null)
     //{computeH0Filter(morningSeries);}     
     
     
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-i-1] = price.get(price.size()-shift-i).doubleValue();}
     logprice[0] = logprice[1];     
          
     int n_interps = 4;
     double[] interpolations = new double[n_interps];
     //interpolations[0] = 0; 
     //interpolations[1] = 0.15;
//      interpolations[1]
     interpolations[0] = .10;
     interpolations[1] = .30;
     interpolations[2] = .66;
     //for(k=1;k<n_interps;k++) {interpolations[k] = interpolations[k-1] + 1.0/24.0;}
     interpolations[n_interps-1] = 1;
     order = new int[n_interps];
     morn_rets = new double[n_interps];
   
     for(k=0;k<n_interps;k++)
     {
       order[k]=k;
       interpolation = interpolations[k];
       
       if(new_day)
       {out = computeInterpolationFrequency(tseries, interpolation);}
       else
       {out = computeInterpolationFrequency(morningSeries, interpolation);} 
       
       System.arraycopy(out,lag4,b_h0,0,n_rep*L);
     
     
       start = n_obs-trade_length;
       for(i=start;i<n_obs;i++)
       {
         sum = 0.0;
         for(j=1;j<n_rep;j++)
         {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
         }
         x[i-start] = tseries[i];
         localsignal[i-start] = sum;
       }     
     
       if(!new_signal) //replace first trade_length-1
       {
        for(j=0;j<trade_length-1;j++)
        {
          //System.out.println("localsig = " + localsignal[j] + ", global = " + this.signal[j]);
          localsignal[j] = this.signal[j];
        }
       }
       System.out.println("new sigval = " + localsignal[localsignal.length-1] + ", price = " + Math.exp(logprice[logprice.length-1])); 
    
       System.out.println("signal length = " + localsignal.length + ", " + logprice.length + ", " + trade_length);
       insampleTradingDiff(logprice, localsignal, trade_length);
       total_return = account[trade_length-1];
       morn_rets[k] = total_return;   
         
       //I_MDFA.plotData(account, account, account.length-1);
       rankCoefficient(pnl, pnl.length);
       System.out.println("max_rank = " + rank_coeff + " at " + interpolation + ", return = " + total_return);
       
       if(total_return > max_ret) 
       //if(rank_coeff > max_rank)
       {
         max_ret = total_return; max_interp = interpolation;
         System.out.println("NEW MAX FOUND: max_rank = " + rank_coeff + " at " + max_interp + ", return = " + total_return);
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;
       
       count++;   
     }       
         
     //sort_sims(n_interps,morn_rets, order);
     
     //for(k=0;k<n_interps;k++)
     //{System.out.println(order[k]);}
     
     //max_interp_value = morn_rets[order[5]];
     if(new_day)
     {out = computeInterpolationFrequency(tseries, max_interp_value);}
     else
     {out = computeInterpolationFrequency(morningSeries, max_interp_value);}     
     
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }       
    
   }   
   
   
 /*  
   public void H0_Optimization_Update(double delta, double interp_start)
   {
   
     int i,j,k,l;
     int trade_length = advance; double localdecay = 0.0; 
     double interpolation = interp_start;
     int K,K1;
     double prev_choice = max_interp_value;
     
     
     double max_ret = -100; double max_rank = -1.0; 
     double max_interp = 0.0;
     
     int n = n_obs;
     int n_total = n_obs + advance;

     int shift = 1;
     
     //--- create new Gamma
     double sum; double total_return, total_return2;
     K = (int)(n/2); double om; K1 = K+1;  
         
     
     int start; int flength = n_obs-L+1;
     
     double[] x = new double[trade_length];
     double[] logprice = new double[trade_length];
     double[] localsignal = new double[trade_length];
     double[] b_h0 = new double[n_rep*L];
     double[] out;
     
     
     //-------- Now get previous day's time series data 
     for(i=0;i<trade_length;i++)
     {logprice[trade_length-i-1] = price.get(price.size()-shift-i).doubleValue();}
     
     
//      for(i=0;i<12;i++)
//      {System.out.print(df4.format(logprice[i]) + " ");}
//      System.out.println(df4.format(logprice[12]));       
//  
//      for(i=0;i<12;i++)
//      {System.out.print(df4.format(morningSeries[n_obs - 1 - i]) + " ");}
//      System.out.print(df4.format(morningSeries[n_obs - 1 - 12]));     
//      System.out.println("series lenght = " + morningSeries.length);
     
     
//      System.out.println("09:30 data at revision: " + in_sample[n_obs-1] + " " + in_sample[n_obs + n_obs-1] + " " + in_sample[n_obs*2 + n_obs-1] + " " + in_sample[n_obs*3 + n_obs-1] + " " 
//       + logprice[0]);
                 
     int lag4 = 2*(n-L+1) + 2*(n_rep+1)*K1;
     //int adapt_n, adapt_L,  adapt_i1, adapt_i2, univariate; 
     //double adapt_lambda, adapt_alpha, adapt_sm, adapt_dec, adapt_dec2, adapt_cross;  
  
  
      while(interpolation <= .9999)
      {
//      
//        //mdfa.setRegularization(mdfa.smooth, localdecay, interpolation, mdfa.cross);
//        //mdfa.computeFilterGeneral(true);
//      
        out = mdfa.H0_IMDFAreg(morningSeries, n_obs, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                        cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, interpolation, mdfa.cross, 
                        mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec,
                        mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);     
 
 
 
//        //System.out.println("interpolation value = " + interpolation);
        System.arraycopy(out,lag4,b_h0,0,n_rep*L);
//        
//        //System.arraycopy(mdfa.b,0,b_h0,0,mdfa.b.length);
//        
//        //apply on time series
         start = n_obs-trade_length;
         for(i=start;i<n_obs;i++)
         {
          sum = 0.0;
          for(j=1;j<n_rep;j++)
          {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
          }
          x[i-start] = tseries[i];
          localsignal[i-start] = sum;
         }

//        //{I_MDFA.plotData(x,localsignal,trade_length);}
//        
       //---- compute the insample-outsample trading
       insampleTradingDiff(logprice, localsignal, trade_length);
       total_return = account[trade_length-1];
         
       //I_MDFA.plotData(account, account, account.length-1);
       rankCoefficient(pnl, pnl.length);
       System.out.println("max_rank = " + rank_coeff + " at " + interpolation + ", return = " + total_return);
       
       //if(total_return > max_ret) 
       if(rank_coeff > max_rank)
       {
         max_rank = rank_coeff; max_ret = total_return; max_interp = interpolation;
         System.out.println("NEW MAX FOUND: max_rank = " + rank_coeff + " at " + max_interp + ", return = " + total_return);
       }
       //System.out.println("MAX RETURN AT = " + max_ret);
       max_interp_value = max_interp;
       max_ret_interp = max_ret;
       
       interpolation = interpolation + intra_delta;   
     }  
  
     //System.out.println("pnl length = " + pnl.length);
    
  
      //---- First get current filter, use 0 for simplicity
     out = mdfa.H0_IMDFAreg(morningSeries, n, L, n_rep, 0, 0, Gamma, 0, mdfa.lambda, mdfa.expweight, 0, cutoff0,
                       cutoff, lag, i1, i2, mdfa.w_const, mdfa.smooth, 0, max_interp, mdfa.cross, 
                       mdfa.Q_cdev, mdfa.iter, mdfa.spec_dens, mdfa.spec, mdfa.ar_p, mdfa.ma_q, mdfa.ar_params, 
                       mdfa.ma_params, mdfa.innvar, mdfa.arma, h0b0);    
     
     b_coeffs = new double[(n_rep-1)*L];
     for(l=0;l<L;l++)
     {
        for(i=0;i<n_rep-1;i++)
        {b_coeffs[L*i + l] = out[lag4 + L*(i+1)+l];}
     }  
     
     System.out.println(b_coeffs[0] + " " + b_coeffs[L] + " " + b_coeffs[2*L]);
     System.out.println(b_coeffs[L-1] + " " + b_coeffs[L+L-1] + " " + b_coeffs[2*L+L-1]);

      System.arraycopy(out,lag4,b_h0,0,n_rep*L);
      start = n_obs-trade_length;
      for(i=start;i<n_obs;i++)
      {
          sum = 0.0;
          for(j=1;j<n_rep;j++)
          {
           for(l=0;l<L;l++) {sum = sum + b_h0[L*j + l]*tseries[n_obs*j + i-l];} // System.out.println(yesterday_series[n_total*j + i-l]);}
           //System.out.println("sum = " + sum);
          }
          x[i-start] = tseries[i];
          localsignal[i-start] = sum;
     }    
     
     for(i=0;i<localsignal.length-1;i++)
     {System.out.print(df4.format(localsignal[i]) + " ");}
     System.out.println(df4.format(localsignal[localsignal.length-1]));
 
     for(i=0;i<x.length-1;i++)
     {System.out.print(df4.format(x[i]) + " ");}
     System.out.println(df4.format(x[x.length-1])); 
     
     
     interp_vals.add(max_interp); max_ranks.add(max_rank);
     System.out.println("Interp value choses = " + max_interp);
        
   }            
   
   */
   
   
   
   
   public double forecastInterpolation(int nsims)
   {
     
     double val; double pval; int i; int n_fsteps = 1;
     double sum=0; double[] iseries; double[] diseries;
     int obs = interp.size();
     Double[] temp = interp.toArray(new Double[0]);
     iseries = new double[obs]; diseries = new double[obs]; diseries[0] = 0.0;
     for(i=0;i<obs;i++) {iseries[i] = temp[i].doubleValue();}
     for(i=1;i<obs;i++) {diseries[i] = iseries[i] - iseries[i-1];}
     //for(i=0;i<obs;i++) {diseries[i] = temp[i].doubleValue();}
     
     int model = 1; int method = 1; n_fsteps = 1; 

     Cronos arma = new Cronos(obs, model, method);
     arma.setData(diseries);
     arma.setARMA_Params(2, 0, 0);
     arma.setNForecastSteps(n_fsteps);
     arma.setNPredictiveSims(n_fsteps, nsims);        
     arma.computeARIMAModel(1);

     sum = 0;
     for(i=0;i<nsims;i++)
     {sum = sum + arma.predictive[i];}   
   
     val = sum/(double)nsims;
     val = arma.predictive[0]; 
     //--- spread out ----------------------
     //if(val < 0) {val = val - .30;}
     //else if(val >= 0) {val = val + .30;}
     //-------------------------------------
     
     pval = iseries[iseries.length-1]; 
     System.out.println("previous = " + pval + "val = " + val); 
     val = pval + val;
     
     if(val < .20)
     {val = val - val/2.0;}
     else if(val > .20 && val < .30)
     {val = val - val/3.0;}
     else if(val > .40 && val < .50)
     {val = val + val/3.0;}
     
     
     if(val >= 1.0) {val = .99;}
     else if(val < 0) {val = 0;}
   
     System.out.println("forecasted interp val = " + val);
     return val; 
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
   
   
   public double forecastInterpolationDFA()
   {
   
     double[] fore;  int l;
     double val; double pval; int i; 
     double sum=0; double[] iseries; int k; double[] diseries;
     fore_obs = interp.size();
    
     
     int Kn = (int)(fore_obs/2); double om;
     int K1n = Kn+1;     
     
     Double[] temp = interp.toArray(new Double[0]);
     iseries = new double[fore_obs]; diseries = new double[fore_obs]; diseries[0] = 0.0;
     for(i=0;i<fore_obs;i++) {iseries[i] = temp[i].doubleValue();}
     for(i=1;i<fore_obs;i++) {diseries[i] = iseries[i] - iseries[i-1];}
     
     forecast.set_nobs(fore_obs); 
     int fle = fore_obs - (Lfore - 1); 
     int lag4 = 2*fle+2*K1n;    
    
    
    
     fore_gamma = new double[K1n];       
     for(k=0; k<=Kn;k++)
     {       
         om = (k*Math.PI/Kn);
         if(om >= high_cut) {fore_gamma[k] = 1.0;}
         else {fore_gamma[k] = 0.0;}
     }       
     forecast.set_Gamma(fore_gamma);    
     
     
     fore = forecast.GEN_IDFAreg(diseries, fore_obs, Lfore, dd, 0, fore_gamma, 0, fore_lambda, fore_exp, 0, 
                        0, high_cut, -1, 0, 0, fore_smooth, fore_decay, fore_decay2, forecast.iter,forecast.spec_dens,forecast.spec,
                        forecast.ar_p, forecast.ma_q, forecast.ar_params, forecast.ma_params, forecast.innvar, forecast.arma);     
   
   
   
   
     fore_coeffs = new double[Lfore];
     for(k=0;k<Lfore;k++) 
     {fore_coeffs[k] = fore[lag4 + k];}
    

     sum = 0.0;
     for(l=0;l<Lfore;l++)
     {sum = sum + fore_coeffs[l]*diseries[diseries.length-1-l];}    
     val = sum;
 
     val = val + .4*val;
 
     pval = iseries[iseries.length-1]; 
     System.out.println("previous = " + pval + "val = " + val); 
     val = pval + val;
     
     
     if(val >= 1.0) {val = .99;}
     else if(val < 0) {val = 0;} 
 
     System.out.println("forecasted interp val = " + val);
     return val; 

   }
   
   
   public void setForecastDFAParameters()
   {
   
   
       int k;
       
       fore_obs = interp.size(); int rep = 1; 
       int Kn = (int)(fore_obs/2); double om;
       int K1n = Kn+1;
       
       high_cut = 1.46;
       fore_lambda = 0;
       fore_exp = 12.9;
       fore_smooth = 0.545;
       fore_decay = 0.198;
       fore_decay2 = 0.644;
       fore_cross = 0;
       Lfore = 126;  
       fore_lag = -1;
       forecast = new IMDFA(fore_obs, 2, Lfore, 0, 0, 0, 1.26, 0, 1);
    
       forecast.set_L(Lfore); 
       forecast.set_nobs(fore_obs);  
       forecast.set_nreps(rep);
       forecast.set_lag(fore_lag);
  
       forecast.setRegularization(fore_smooth, fore_decay, fore_decay2, fore_cross);
       forecast.set_dd(0);
       forecast.set_DD(0);
       forecast.set_bconstraints(0, 0); 
       forecast.set_lambda(fore_lambda); 
       forecast.set_exp(fore_exp); 
       forecast.set_cutoff0(0);
       forecast.set_cutoff(high_cut);
       forecast.set_mdfa(false);
  
       fore_gamma = new double[K1n];       
       for(k=0; k<=Kn;k++)
       {       
         om = (k*Math.PI/Kn);
         if(om >= high_cut) {fore_gamma[k] = 1.0;}
         else {fore_gamma[k] = 0.0;}
       }       
       forecast.set_Gamma(fore_gamma);      
   
   }
   
   
   
    public void computeSharpeRatio(int ann)
    {
      int i; double sum; double mean; double stand;
      int length = account.length;
      double[] diff_account = new double[length];
      
      diff_account[0] = 0.0;
      
      for(i=0;i<length-1;i++)
      {diff_account[i+1] = account[i+1] - account[i];}
      
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
     
    }      
   
  public void setLogTrans(boolean t) {ln_trans = t;}
  
  public double exp(double p) 
  {
    if(ln_trans) {return Math.exp(p);}
    else {return p;}
  }
   
  public void setTrendH0(double cut) //------- sets the prior to the ideal trend
  {
    
    int i; 
    double sum;
    h0b0 = new double[(n_rep-1)*L]; 
    
    h0b0[0] = cut/Math.PI; sum= h0b0[0];
    for(i=1;i<L;i++)
    {h0b0[i] = (1/Math.PI)*Math.sin(cut*i)/((double)i); sum = sum+h0b0[i];} 
    
    for(i=0;i<L;i++)
    {h0b0[i] = h0b0[i]/(sum+(sum-cut/Math.PI));}   
    H0set = true;
  }
  


   public void setRampTrend(double _w0, double _w1)
   {
      int i; 
      double sum,coeff0;
      double w0 = 0; double w1 = 0;
      if(_w0 < _w1)
      {
       w0 = _w0; w1 = _w1; 
       h0b0 = new double[(n_rep-1)*L]; 
              
       //-----compute symmetric filter----
       coeff0 = .5*(w1 + w0)/Math.PI; h0b0[0] = coeff0;  sum = h0b0[0];
       for(i=1;i<L;i++)
       {
          h0b0[i] = (Math.cos(w1*i) - Math.cos(w0*i))/((double)i*i);
          h0b0[i] = -1.0/(Math.PI*(w1 - w0))*h0b0[i]; 
          sum= sum+h0b0[i];
       } 
       for(i=0;i<L;i++)
       {
         h0b0[i] = h0b0[i]/(sum+(sum-coeff0));
       } 
       H0set = true;
      }
      else 
      {System.out.println("omega_0 must be less than omega_1");}
   }


   public void setForecastTrend(int lag)
   {

       int i,k,n; double sum = 0.0; double sum2 = 0.0; double sumi = 0.0;
       h0b0 = new double[(n_rep-1)*L];  
   
       for(k=0;k<L;k++)
       {
         sum=0.0;
         for(n=0;n<=K;n++)
         {
          sum = sum + Gamma[n]*Math.cos(-lag*Math.PI*n/K)*Math.cos(Math.PI*n*k/K);
          sumi = sumi + Gamma[n]*Math.sin(-lag*Math.PI*n/K)*Math.sin(Math.PI*n*k/K);
         }     
         if(k==0) {h0b0[0] = sum*sum;}
         else
         {h0b0[k] = 2.0*(sum*sum + sumi*sumi);} 
         sum2 = sum2 + h0b0[k];
       }
       for(i=0;i<L;i++)
       {h0b0[i] = h0b0[i]/(sum2-h0b0[0]/2.0);}    
   }
  
   
  
  
  
  
 
 
 
  public double stand_dev(double[] x)
  {
    int i; double sd;
    double sum=0;
    for(i=0;i<x.length;i++)
    {sum = sum + x[i];}
    
    double mean = sum/x.length;
    sum = 0;
    for(i=0;i<x.length;i++)
    {sum = sum + (x[i] - mean)*(x[i] - mean);}
    
    sd = Math.sqrt(sum/(double)x.length);
    
    return sd;
  }
   
   
  public double segmentRankCorrelation(int length_rank, double[] rets)
  {
     int i; int obs = rets.length;
     int n_sections = (int)(obs/length_rank);
      
     double min_rank = 1.0; double avg_rank = 0; double rank_coeff = 0;
     for(i=0;i<n_sections;i++)
     {
        double[] section = new double[length_rank];
        System.arraycopy(rets,i*length_rank,section,0,length_rank);
      
        rank_coeff = rankCoefficientSeg(section,length_rank);
        avg_rank = avg_rank + rank_coeff;
        
        if(rank_coeff < min_rank)
        {min_rank = rank_coeff;}       
      }
      avg_rank = avg_rank/n_sections;
      return avg_rank;
  } 
   
   
  public double ulcerIndex(Double[] rets)
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
  
  public double ulcerIndex(double[] rets)
  {
  
    double SumSq = 0; 
    double MaxValue = -100;
    double mean = 0;
  
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
  
   
    double[] cumsum(Double[] data, int n)
    {

      double[] cs = new double[n]; double sum; int k;
      //sum=Math.abs(data[0]); double min = 1000000;
      sum=0; double min = 1000000;

      for(k=0;k<n;k++)
      {
        sum = sum+data[k]; cs[k] = sum; 
        if(cs[k] < min) {min = cs[k];}
      }


      return cs;  
    }   
   

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


      return cs;  
    }     
   
   
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
           n_rep = temp_nreps;
          
             System.out.println("Setting H0 Filter");
             L = tempL;
             
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
             setH0B0(vals, n_rep, L);   
           
           System.out.println("H0 Filter set");
           br.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
         catch(IOException ioe){System.out.println("IO out error..." + ioe);}
         
  }     
   
   
  public void setH0B0(double[] b0, int _nr, int _L)
  {
    int i,l;
    System.out.println("nr = " + _nr + ", L = " + _L);
    if(_nr == n_rep && _L == L)
    {
      h0b0 = new double[(n_rep-1)*L];
      //--- scrape off first column,, only last n_rep-1 columns are nonzero
      for(i=1;i<n_rep;i++)
      {
        for(l=0;l<L;l++)
        {
          h0b0[(i-1)*L + l] = b0[L*i + l];
        }
      }   
      H0set = true;
      System.out.println("Filter set -" + H0set);
    }
    else
    {System.out.println("Dimensions of the H0 filter must match the current filter");}
  }   
   
  public void computeH0B0(int num)
  {
     int l,k,count;
     
     h0b0 = new double[L*(n_rep-1)];
     for(count = 0; count < num; count++)
     {
       Filter coeffs = filters.get(filters.size()-1-count);
       
       for(k=1;k<n_rep;k++)
       {
         for(l=0;l<L;l++)
         {h0b0[L*(k-1) + l] = h0b0[L*(k-1) + l] + coeffs.b_coeffs[L*k + l];}            
       }       
     } 
    
     for(k=1;k<n_rep;k++)
     {
         for(l=0;l<L;l++)
         {h0b0[L*(k-1) + l] = h0b0[L*(k-1) + l]/(double)num;}
     }     
  }
   
  //assumes morning data already set 
  public void computePeriodogram()
  {
    int i,j;
    max_peak = 0; peak_loc = 0;
    diff_band = 0.0;
    int K1 = (int)(n_obs/2)+1; int K = K1-1;
    double[] data = new double[n_obs];
    double[] periodogram = new double[K1];    
    System.arraycopy(morningSeries,0,data,0,n_obs);
    double sumr,sumi;
    for(j=0;j<K1;j++)
    {  
      sumr = 0.0; sumi=0.0;  
      for(i=0;i<n_obs;i++)
      {
        sumr = sumr + data[i]*Math.cos(Math.PI*(i+1.0)*j/K); 
        sumi = sumi + data[i]*Math.sin(Math.PI*(i+1.0)*j/K);
      }
      periodogram[j] = (sumr*sumr + sumi*sumi)/Math.sqrt(2*Math.PI*n_obs);
      
      
      if(j <= (int)(K/2))      
      {
       if(periodogram[j] > max_peak) {max_peak = periodogram[j]; peak_loc = Math.PI*j/(double)K;}            
       diff_band = diff_band + periodogram[j];
      }      
    }
    //System.out.println("max peak = " + max_peak + ", peak_loc = " + peak_loc);
    diff_band = diff_band/K;   
  }
   
   
  public void computePeriodogram(double[] series)
  {
    int i,j;
    max_peak = 0; peak_loc = 0;
    diff_band = 0.0;
    int K1 = (int)(n_obs/2)+1; int K = K1-1;
    double[] data = new double[n_obs];
    double[] periodogram = new double[K1];    
    System.arraycopy(series,0,data,0,n_obs);
    double sumr,sumi;
    int w1 = (int)(cutoff1*K/Math.PI);
    
    
    for(j=0;j<K1;j++)
    {  
      sumr = 0.0; sumi=0.0;  
      for(i=0;i<n_obs;i++)
      {
        sumr = sumr + data[i]*Math.cos(Math.PI*(i+1.0)*j/K); 
        sumi = sumi + data[i]*Math.sin(Math.PI*(i+1.0)*j/K);
      }
      periodogram[j] = (sumr*sumr + sumi*sumi)/Math.sqrt(2*Math.PI*n_obs);
      
      
      if(j <= w1)      
      {
       if(periodogram[j] > max_peak) {max_peak = periodogram[j]; peak_loc = Math.PI*j/K;}            
       diff_band = diff_band + periodogram[j];
      }      
    }
    diff_band = diff_band/K;   
  }   
   
   

   
   
   
   
}







