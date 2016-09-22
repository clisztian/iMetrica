package ch.imetrica.dataControl;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;


/*------------------------------------------------
  Canvas for plotting series
--------------------------------------------------*/
public class SimCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int n_obs; //---total number of observations
    double[] tseries;  //raw time series n_rep x n_obs

    double dataMax, dataMin, dataNorm;
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
    int n_hist; 
    boolean[] plots;
    boolean plot_tar;
    Color[] colorHist;
    int track_pos_t,track_pos_xtarg;
    int[] track_pos_x;  
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

    Color highlight; 
    int hlight_indicator; 
    String[] info; 
    boolean display;
    Font mono;
    int canvas_shift = 0;
    Ellipse2D ellipse;
    Rectangle2D rectangle;
    
    boolean dates_set = false; 
    String[] dates; 
    
    public SimCanvas(int w, int h, int nobs, int nrep)
    {
      // Initilize everything-------------------
      int i; 
      n_obs = nobs; n_rep = nrep; 
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      n_hist = 0; 
      sim_series = new ArrayList<double[]>(); 
      target_series = new double[n_obs]; plot_tar = false;
  
      colorHist = new Color[15];
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
      df2 = new DecimalFormat("##.#####");
      foreRed = new Color(232, 138 , 187);
      highlight = new Color(255,255,255);
      plot_acf = false;
      acf1 = false; acf2 = false; cross1 = false; cross2 = false;

      track_pos_x = new int[12]; 
      valuetarg = "";
      value = new String[12];
      for(i=0;i<12;i++) {value[i] = new String(""); track_pos_x[i] = 0;}
      
      
      
      hlight_indicator = -1; display = false;
      colorHist[0] = new Color(103,177,246); colorHist[1] = new Color(255,177,246); 
      colorHist[2] = new Color(103,100,246); colorHist[3] = new Color(103,255,246);  
      colorHist[4] = new Color(103,177,0);   colorHist[5] = new Color(239,177,0);
      colorHist[6] = new Color(239,49,0);    colorHist[7] = new Color(239,77,166); 
      colorHist[8] = new Color(135,73,169);  colorHist[9] = new Color(135,197,167);
      colorHist[10] = new Color(200,41,2);    colorHist[11] = new Color(239,200,166); 
      colorHist[12] = new Color(100,73,169);  colorHist[13] = new Color(135,100,167);
      colorHist[14] = new Color(135,197,10);


      acf_series_1 = new double[n_obs];
      acf_series_2 = new double[n_obs];
      cross_12 = new double[n_obs]; 
      cross_21 = new double[n_obs];  

      plots = new boolean[15]; for(i=0;i<15;i++){plots[i] = false;}

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
       
       this.dates = new String[n_obs];
       
       if(dates.length == n_obs)
       {this.dates = dates;}
       else if(dates.length > n_obs)
       {
         for(int i = 0; i < n_obs;i++)
         {this.dates[n_obs - 1 - i] = dates[dates.length-1-i];}
       }
       else if(dates.length < n_obs)
       {
         for(int i = 0; i < dates.length;i++)
         {this.dates[n_obs - 1 - i] = dates[dates.length-1-i];}
       }  
    }
    
    public void changeHighlight(int c)
    {hlight_indicator = c; go();}


    public void setplot(int i, boolean sel)
    {plots[i] = sel; go();}
  
    public void setTarget(double[] v)
    {

      target_series = new double[v.length];
      System.arraycopy(v, 0, target_series, 0, v.length);
      go();

    }
    public void plotTarget(boolean sel) {plot_tar = sel; go();}

    public void setSeries(int i, double[] v)
    {  
      sim_series.set(i,v); 
    }

    public void addSeries(double[] v)
    {
      sim_series.add(v); 
    }

    public void replaceSeries(int i, double[] v)
    {
      sim_series.remove(i); 
      sim_series.set(i,v);
    }

    public void computeDataMax()
    {
      int i,k; dataMax = -1000000; dataMin = 10000000;
      double[] t_series; 
      if(!plot_acf)
      {
       for(k=0; k < sim_series.size(); k++)
       {
        if(plots[k])
        {
         t_series = sim_series.get(k);
         for(i=0; i < n_obs; i++)
         {
            if(t_series[i] > dataMax) dataMax = t_series[i];
            else if(t_series[i] < dataMin) dataMin = t_series[i];   
         }
        }
       }
       if(plot_tar)
       {
         for(i=0; i < n_obs; i++)
         {
            if(target_series[i] > dataMax) dataMax = target_series[i];
            else if(target_series[i] < dataMin) dataMin = target_series[i];   
         }

       }     
      dataNorm = Math.abs(dataMax - dataMin);
     }
     else
     {dataMax = 1.0; dataMin = -1.0; dataNorm = 2.0;}
    }

    public void setNobs(int n) {n_obs = n;}

    public void clearSeries()
    {
      sim_series.clear();  
      for(int i=0;i<15;i++) {plots[i]=false;} plot_tar = false;
      go();
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

     int i,j,N,k;
     
     int t0, t1, x0, x1;
     int nobsP;  int p = 10;
     double[] tseries;
     //----- computeDataMax -------
     computeDataMax();
     
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;


     BasicStroke orig = new BasicStroke((float)1.3);
     BasicStroke thick = new BasicStroke((float)1.6); 
     BasicStroke big = new BasicStroke((float)1.5);

        //Stroke def = g2d.getStroke();

       //--------------------- Draw grid with or without forecast points ------------
     float[] dash1 = {10.0f};
     BasicStroke dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);

     N = n_obs;
     g2d.setStroke(dashed);
     g2d.setPaint(myGray2);
     for(i=0; i < 9; i++)
     {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
     }
     g2d.setPaint(myGray);
     g.drawString((String)df.format(dataMax), 5, 15);
     g.drawString((String)df.format(dataMin), 5, height - 24);
     nobsP = (int)Math.floor((double)N/10); 
     for(i=1; i <= nobsP; i++) 
     {
        t0 = (int)(((double)(i*10)/N)*(double)width);
        g2d.setPaint(myGray2);
	g2d.drawLine(t0, 0, t0, height-20);
        g2d.setPaint(myGray);
        g.drawString((String)"" + p, t0, height - 5); 
        p = p + 10;
     }
 
     g2d.setStroke(orig);

    if(!plot_acf)
    {

     for(k=0;k<sim_series.size();k++)
     {
       if(plots[k])
       {
         tseries = sim_series.get(k);
         g2d.setPaint(colorHist[k]); if(hlight_indicator == k){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((tseries[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((tseries[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
       }
     }
     if(plot_tar)
     {
         g2d.setPaint(new Color(155,204,224)); g2d.setStroke(big);
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((target_series[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((target_series[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
     }          
   }
   else
   {
      if(acf1)
      {
         g2d.setPaint(colorHist[0]); 
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((acf_series_1[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((acf_series_1[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }

      }
      if(acf2)
      {
         g2d.setPaint(colorHist[0]); 
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((acf_series_2[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((acf_series_2[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
      }
      if(cross1)
      {
         g2d.setPaint(colorHist[1]); 
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((cross_12[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((cross_12[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }

      }
      if(cross2)
      {
         g2d.setPaint(colorHist[2]); 
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((cross_21[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((cross_21[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
         }
      }
   }
   
   
   if(plot_tracker)
   {
     
    int n_count = 0; int w = 4; int h = 4; 
    if(!plot_acf)
    {

     if(plot_tar)
     {    
        
        g2d.setPaint(Color.white);
        ellipse = new Ellipse2D.Double(track_pos_t, (height) - track_pos_xtarg, w, h);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);     
     
         
        g.setFont(mono);
        g2d.setPaint(Color.GREEN);
        g.drawString("portfolio = " + valuetarg, 5, 28 + n_count*10);
     }   
    
     if(dates_set)
     {
      g.setFont(mono);
      g2d.setPaint(Color.GREEN);    
      g.drawString(date, 5, 15);
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
        if(names_set) {g.drawString(asset_names[k] + " = " + value[k], 5, 15 + n_count*15);}
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

        int j,k; double[] tseries; 
        int N = n_obs; date = "";
        if(plot_tracker)
        {
            
           j = (int)(((double)(N))*e.getX()/(double)width);  
        
           track_pos_t = (int)(((double)j/(double)(N-1))*(double)width) - 2;
           //t1 = (int)(((double)(j+1)/(double)trade_obs)*(double)width);
           if(plot_tar) //account
           { 
             track_pos_xtarg = (int)(((target_series[j] - dataMin)/dataNorm)*(double)height) + 2; 
             valuetarg = df.format(target_series[j]);
           }
           
           for(k=0;k<sim_series.size();k++)
           {
            if(plots[k])
            {
              tseries = sim_series.get(k);           
              track_pos_x[k] = (int)(((tseries[j] - dataMin)/dataNorm)*(double)height) + 2; 
              value[k] = df2.format(tseries[j]);
            }   
           }       
           if(dates_set) {date = dates[j];}
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
