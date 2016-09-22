package ch.imetrica.mdfaTradingStrategies;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.swing.JPanel;

public class EvolutionCanvas extends JPanel
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
    int[] track_pos_x, track_pos_xvar;  
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
    String date, aggregateName; 
    String[] value; 
    DecimalFormat df,df2,df3,df4;
    BasicStroke dashed;
    BasicStroke orig;
    float[] dash1;
    Color myGray;
    Color myGray2;
    Color plotGray;
    Color foreRed;
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
    int track_pos_tar; 
    
    String[] saved_asset; 
    ArrayList<String[]> asset_info_array;
    
    boolean meta_set = false; 
    boolean asset_set = false;
    boolean dates_set = false; 
    String[] dates,meta_info,asset_info; 
    ArrayList<double[]> stock_perf; 
    String asset,meta;
    boolean stock_set = false;
    boolean invest_in_all = false;
    String[] tokens, tokens2; String delims = "[ ]+";
    
    public EvolutionCanvas(int w, int h, int nobs, int nrep)
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
      
      
      
      asset_info_array = new ArrayList<String[]>();
      foreRed = new Color(232, 138 , 187);
      highlight = new Color(255,255,255);
      plot_acf = false;
      acf1 = false; acf2 = false; cross1 = false; cross2 = false;

      track_pos_tar = 0;
      target_series = new double[0];
      track_pos_x = new int[max_series]; track_pos_xvar = new int[max_series];
      valuetarg = ""; 
      value = new String[max_series];
      for(i=0;i<max_series;i++) {value[i] = new String(""); track_pos_x[i] = 0; track_pos_xvar[i] = 0;} 
      
      varPlotNames = new ArrayList<String>();
      seriesPlotNames = new ArrayList<String>();
      
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

      stock_perf = new ArrayList<double[]>();
      acf_series_1 = new double[n_obs];
      acf_series_2 = new double[n_obs];
      cross_12 = new double[n_obs]; 
      cross_21 = new double[n_obs];  

      
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
    
    public boolean setDates(String[] d)
    {
       boolean match = true;
       if(dates_set)
       {
          tokens = d[d.length-1].split(delims); 
          tokens2 = dates[dates.length-1].split(delims);
          match = tokens[0].equals(tokens2[0]);
       }
      
       if(dates_set && match)
       {
          this.dates = new String[d.length];
          System.arraycopy(d,0,this.dates,0,d.length); 
       }
       else if(!dates_set)
       {
        dates_set = true;       
        this.dates = new String[d.length];
        System.arraycopy(d,0,this.dates,0,d.length); 
       }
       else
       { 
         System.out.println("Problem setting dates, there is a discrepency as the dates don't match");
         System.out.println(d[d.length-1] + " does not equal " + dates[dates.length-1]);
         dates_set = false;
       }
       return dates_set; 
    }
    
    public void setAsset(String[] dates)
    {
       asset_set = true;       
       asset_info = new String[dates.length];
       System.arraycopy(dates,0,asset_info,0,dates.length); 
    }
    
    public void saveAsset()
    {
      if(asset_set)
      { 
        asset_info_array.add(asset_info);
        saved_asset = new String[asset_info_array.size()];
        for(int i=0;i<saved_asset.length;i++)
        {saved_asset[i] = new String(" ");}
      }
    }
    
    public void deleteSavedAssets()
    {
      asset_info_array.clear();
    }
    
    public void setMeta(String[] dates)
    {
       meta_set = true;       
       meta_info = new String[dates.length];
       System.arraycopy(dates,0,meta_info,0,dates.length); 
    }    
    
    
    public void addAggregate(double[] v, String name)
    {
      target_series = new double[v.length];
      System.arraycopy(v, 0, target_series, 0, v.length);
      aggregateName = new String(name);
      go();
    }
    

    
    public void setStockReturns(double[] rets)
    {
     stock_set = true;
     stock_perf.add(rets);
     go();
    }
    
    public void clearAssetPlots()
    {
      stock_set = false;
      stock_perf.clear();
      go();
    }
    
    public void changeHighlight(int c)
    {hlight_indicator = c; go();}


    public void setplot(int i, boolean sel)
    {plots[i] = sel; go();}
    
    public void setvarplot(int i, boolean sel)
    {var_plots[i] = sel; go();}
  
    public void setTarget(double[] v)
    {

      target_series = new double[v.length];
      System.arraycopy(v, 0, target_series, 0, v.length);
      go();

    }
    public void plotTarget(boolean sel) {plot_tar = sel; go();}

    public void setSeries(int i, double[] v, String name)
    {  
      sim_series.set(i,v); seriesPlotNames.set(i,name); computeNobs(); go(); //setplot(sim_series.size()-1, true);
    }

    public void addSeries(double[] v, String name)
    {
      sim_series.add(v); seriesPlotNames.add(name); computeNobs(); go();//setplot(sim_series.size()-1, true);
    }

    public void replaceSeries(int i, double[] v)
    {
      sim_series.remove(i); 
      sim_series.set(i,v);
    }

    public void addVarSeries(Double[] var, String name)
    {
      var_series.add(var); varPlotNames.add(name); computeNobs(); setvarplot(var_series.size()-1, true);
    }
    
    public void setVarSeries(int i, Double[] v, String name)
    {var_series.set(i,v);  varPlotNames.set(i,name); computeNobs(); setvarplot(i, true);}    
    
    
    public void computeDataMax()
    {
      int i,k; dataMax = -1000000; dataMin = 10000000; sdataMax = -1000000; sdataMin = 10000000;
      double[] t_series; 
      
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
      
      if(plot_tar)
      {
         n_obs = target_series.length;
         for(i=0; i < n_obs; i++)
         {
            if(target_series[i] > dataMax) dataMax = target_series[i];
            else if(target_series[i] < dataMin) dataMin = target_series[i];   
         }
      }
      
      
      
//       if(stock_set)
//       {
//         n_obs = stock_perf.length;
//         for(i=0; i < n_obs; i++)
//          {  
//             //System.out.println(stock_perf[i]);
//             if(stock_perf[i] > sdataMax) sdataMax = stock_perf[i];
//             else if(stock_perf[i] < sdataMin) sdataMin = stock_perf[i];   
//          }
//       }
      
      //sdataNorm =  Math.abs(sdataMax - sdataMin);
      dataNorm = Math.abs(dataMax - dataMin);
      dataNormVar = Math.abs(dataMaxVar - dataMinVar);
    }

    public void setNobs(int n) {n_obs = n;}

    public void clearSeries()
    {
      sim_series.clear(); seriesPlotNames.clear(); 
      for(int i=0;i<max_series;i++) {plots[i]=false;} 
      if(var_series.size() == 0) {dates_set = false;}
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

    public double min(double[] s)
    {
      double mins = 100000; 
      for(int i=0;i<s.length;i++) {if(s[i] < mins) {mins = s[i];}}
      return mins;  
    }
    
    public double max(double[] s)
    {
      double mins = -100000; 
      for(int i=0;i<s.length;i++) {if(s[i] > mins) {mins = s[i];}}
      return mins;  
    }    

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N,k;
     int t0, t1, x0, x1;
     double[] tseries;
     //----- computeDataMax -------
     computeDataMax();
     computeNobs();
     //valuetarg = "";
     
     double stmin,stmax,stnorm;     
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;

     int diff_n; 
     BasicStroke orig = new BasicStroke((float)1.3);
     BasicStroke thick = new BasicStroke((float)1.6);
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
         
   
   if(plot_tar && target_series.length > 0) //plot the aggregate strategy
   {
         diff_n = target_series.length - min_obs;
         g2d.setPaint(new Color(0,204,204)); if(hlight_indicator == k){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((target_series[diff_n + j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((target_series[diff_n + j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }   
   }
   
   
   
   g2d.setPaint(myGray);
   //if(stock_set && (sim_series.size() > 0))
   if(stock_set)
   {
     for(k = 0; k < stock_perf.size(); k++)
     {
       g2d.setPaint(new Color(5*k,69,75));
       
       double[] stock = stock_perf.get(k);
       stmin = min(stock);
       stmax = max(stock); 
       stnorm = Math.abs(stmax - stmin);
       
       diff_n = stock.length - min_obs;
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((stock[diff_n + j] - stmin)/(stnorm))*(double)height);
          x1 = (int)(((stock[diff_n + j + 1] - stmin)/(stnorm))*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
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
     
     if(asset_set)
     {
      g.setFont(mono);
      g2d.setPaint(Color.PINK);    
      g.drawString(asset, 5, height - 15);
      
      if(asset_info_array.size()>0)
      {
        for(j=0;j<asset_info_array.size();j++)
        {
         g2d.setPaint(new Color(100,100,115));    
         g.drawString(saved_asset[j], 5, height - 15*(j+2));         
        }
      }
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
        else {g.drawString(value[k], 5, 15 + n_count*15);}
 
       }
     } 
     
     if(plot_tar && target_series.length > 0)
     {
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t, (height) - track_pos_tar, w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);     
     
         
        g.setFont(mono);
        //g2d.setPaint(Color.GREEN);
        g2d.setPaint(new Color(0,204,204));
        g.drawString("Aggregate: " + aggregateName + " " + valuetarg, 5, 30);
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

        int j,k; double[] tseries; 
        int N = min_obs; date = ""; int diff_n; asset = ""; meta = ""; 
      
        if(plot_tracker)
        {
            
           j = (int)(((double)(N))*e.getX()/(double)width);  
        
           track_pos_t = (int)(((double)j/(double)(N-1))*(double)width) - 2;
  
           for(k=0;k<sim_series.size();k++)
           {
            if(plots[k])
            {
              tseries = sim_series.get(k);  diff_n = tseries.length - min_obs;          
              track_pos_x[k] = (int)(((tseries[diff_n+j] - dataMin)/dataNorm)*(double)height) + 2; 
              value[k] = df2.format(tseries[diff_n+j]);
            }   
           }  
           
           if(plot_tar && target_series.length > 0)
           { 
              diff_n = target_series.length - min_obs;          
              track_pos_tar = (int)(((target_series[diff_n+j] - dataMin)/dataNorm)*(double)height) + 2; 
              valuetarg = df2.format(target_series[diff_n+j]);
           
           }
           
           if(dates_set) {diff_n = dates.length - min_obs; date = dates[diff_n+j];}// System.out.println(diff_n + " " + dates.length + " " + min_obs + " " + j);}
           if(asset_set) {diff_n = asset_info.length - min_obs; asset = asset_info[diff_n+j];}
           
           
           for(k=0;k<asset_info_array.size();k++)
           {
             String[] asset_infotemp = asset_info_array.get(k);   
             diff_n = asset_infotemp.length - min_obs; saved_asset[k] = asset_infotemp[diff_n+j];
           } 
  
           if(meta_set) {diff_n = meta_info.length - min_obs; meta = meta_info[diff_n+j];}
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

