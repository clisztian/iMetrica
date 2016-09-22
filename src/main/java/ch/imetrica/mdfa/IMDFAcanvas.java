package ch.imetrica.mdfa;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;
import java.awt.geom.Rectangle2D;
/*------------------------------------------------
  Canvas for plotting mdfa tseries data on time domain 
  - Only accepts data of form  tseries by row
--------------------------------------------------*/
public class IMDFAcanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int n_obs; //---total number of observations
    int n_rep; //---number of series (if n_rep = 1, then dfa) 
    int L;     //---number of coefficients in realTime filter
    int xf_obs; //--- Filtered length, ya'ni xf_obs = n_obs - L + 1   
    int lag; 

    boolean zeroLine; 

    double[] tseries;  //raw time series n_rep x n_obs
    double[] xf;       //filtered series n_rep x (n_obs - L + 1)
    double[] xf_hybrid;
    double[] xf_orig;  //original data
    double[] xf_zpc;   //zpc only filtered data
    
    boolean plot_hybrid;
    boolean mdfa;      //mdfa or just dfa

    boolean[] ts_plots;   //---- original time series plots n_rep array   
    boolean xf_plot;      //---- filtered data 
    double[] x_sym;       //---- symmetric filtered data, if available
    boolean sym_plot;     //---- sym plot
    boolean plotForecast; //---- plot extended forecast data if simulated
    double[] extend;

    ArrayList<double[]> historical;  //----- plots historical saved plots
    ArrayList<String> histParams;
    int n_hist; 
    boolean plotHist;
    Color[] colorHist;

    boolean shade = false; 
    int shade_points = 0;
    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;
    Rectangle2D rectangle;
    DecimalFormat df;
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
    int sym_shift = 0;
    int outSamp = 0;
    int win_shift = 0;

    public IMDFAcanvas(int w, int h, int nobs, int nrep, int _L)
    {
      // Initilize everything-------------------
      int i; 
      n_obs = nobs; n_rep = nrep; L = _L; lag = 0;
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      n_hist = 0; zeroLine = true;
      historical = new ArrayList<double[]>(); 
      histParams = new ArrayList<String>();
      plotHist=false; colorHist = new Color[10];
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
      foreRed = new Color(232, 138 , 187);
      highlight = new Color(255,255,255);

      hlight_indicator = -1; display = false;
      colorHist[0] = new Color(103,177,246); colorHist[1] = new Color(255,177,246); 
      colorHist[2] = new Color(103,100,246); colorHist[3] = new Color(103,255,246);  
      colorHist[4] = new Color(103,177,0);   colorHist[5] = new Color(239,177,0);
      colorHist[6] = new Color(239,49,0);    colorHist[7] = new Color(239,77,166); 
      colorHist[8] = new Color(135,73,169);  colorHist[9] = new Color(135,197,167);
 
      //---- initialize plots, set all to false------
      ts_plots = new boolean[11]; 
      for(i=0;i<10;i++) {ts_plots[i] = false;}
      xf_plot = false;  plotForecast = false; extend = new double[36];

      tseries = new double[n_rep*n_obs];
      xf = new double[n_obs-L+1]; 
      x_sym = new double[n_obs-L+1];
      xf_obs = n_obs-L+1;
      mdfa = true; outSamp = 0; plot_hybrid = false;
    } 

    public void changeHighlight(int c)
    {hlight_indicator = c; displayInfo(c); go();}

    public void setPlotDim(double[] series, int N)
    {
        int i;
        dataMax = -1000000.0; dataMin = 1000000.0;
        for(i=0;i<n_obs;i++)
        {
         if(series[i] < dataMin) dataMin = series[i];
         else if(series[i] > dataMax) dataMax = series[i]; 
        } 
        dataNorm = Math.abs(dataMax - dataMin);    
    }
  
    public void setNRep(int n) {n_rep = n;}
    public void setNobs(int n) {n_obs = n; xf_obs = n_obs - L + 1;}
    public void setL(int n) {L = n; xf_obs = n_obs - L + 1;}
    public void setLag(int l) {lag = l;}

    public void setOutSamp(int o) {outSamp = o;}
 
    public void setSymShift(int _L) {sym_shift = _L;}

    public void shadeRegion(boolean s, int p)
    {shade = s; shade_points = p; go();}
    
    //----- set raw data and filtered data
    public void setTseries(double[] ts, int N, int rep)
    {   
       int le = ts.length; this.n_obs = N; this.n_rep = rep; plotForecast = false;
       //if(le != N*rep) {System.out.println("Dimensions not equal, truble tuble\n");}
       //else{this.n_obs = N; this.n_rep = rep;}
       this.tseries = new double[le];  
       System.arraycopy(ts, 0, this.tseries, 0, le);  
       this.setPlotDim(ts,le);
    }

    public void setRTSignal(double[] _xf, int N)
    {
       int le = _xf.length;
       xf_obs = le; 
       this.xf = new double[le]; this.xf_orig = new double[le];
       System.arraycopy(_xf, 0, this.xf, 0, le); 
       System.arraycopy(_xf, 0, this.xf_orig, 0, le);
       xf_hybrid = new double[xf_obs];
       //--- if filter fine, no redim needed 
    }

    public void setHybridZPC(double[] _xf)
    {
       xf_obs = _xf.length;
       this.xf_hybrid = new double[xf_obs];
       System.arraycopy(_xf, 0, this.xf_hybrid, 0, xf_obs); 
       go();
    }

    public void setZPC(double[] _xf)
    {
       xf_obs = _xf.length;
       this.xf_zpc = new double[xf_obs];
       System.arraycopy(_xf, 0, this.xf_zpc, 0, xf_obs); 
       go();
    }

 
    public void plotHybrid(int c)
    {
       int le = xf_obs;
       if(c==2) {System.arraycopy(this.xf_hybrid, 0, this.xf, 0, le);}
       else if(c==0) {System.arraycopy(this.xf_orig, 0, this.xf, 0, le);}
       else if(c==1) {System.arraycopy(this.xf_zpc, 0, this.xf, 0, le);}
       go();
    }
    

    public void setSymSignal(double[] _xf)
    {
       int le = _xf.length; 
       this.x_sym = new double[le];
       System.arraycopy(_xf, 0, this.x_sym, 0, le);        
    }
    

    public void setPlots(int i, boolean sel)
    {
       if(i==0) 
       {xf_plot = sel;} 
       else if(i == 12)
       {sym_plot = sel;}
       else
       {ts_plots[i] = sel;}
       go();
    }

   public void saveToHistorial()
   { 
     if(n_hist < 10)
     {
      double[] temp = new double[xf.length+1];
      for(int i=0; i < xf.length; i++) {temp[i+1] = xf[i];}
      temp[0] = 1.0*lag;
      historical.add(temp);
      n_hist++;
     }
     else
     {System.out.println("Only ten plots can be saved.");}
   }

   public void setParamSnapshot(String par) 
   {
     if(n_hist < 10)  
     {
       histParams.add(par);
     }    
   }    

   public void clearHistorical()
   {historical.clear(); histParams.clear(); n_hist=0; go();}

   public void plotHistorical(boolean p) 
   {plotHist = p; go();}
  
   public void saveFromQueue(int k)
   {
      String file = "iMetric-filter-"+k+".dat";
      double[] temp = historical.get(k);
      try{  
            
           PrintWriter out = new PrintWriter(new FileWriter(file));
           
           out.println(temp.length);
           for(int i =0; i < temp.length; i++) {out.println(temp[i]);}
           out.println(histParams.get(k));

           out.close(); System.out.println("Data successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();} 
   } 

   /*
       paramString = Integer.toString(L)+"\n"+Integer.toString(Lag)+"\n"+Integer.toString(i1)+"\n"+Integer.toString(i2)+"\n"
                    +Double.toString(lambda)+"\n"+Double.toString(expweight)+"\n"+Double.toString(smooth)+"\n"
                    +Double.toString(decay)+"\n"+Double.toString(cross)+"\n"
                    +Double.toString(mdfa.criteria)+"\n"+Double.toString(mdfa.degrees);
   */

   public void setCanvasShift(int _shift)
   {canvas_shift = _shift;}

   public void displayInfo(int k)
   {

      if(k >= 0)
      {
       String par = histParams.get(k); String[] tokens; String delims = "[\n]+";
      
       tokens = par.split(delims); 
       info = new String[4];      

       info[0] = new String("L = " + tokens[0] + ", lag = " + tokens[1] + ", i1 = " + tokens[2] + ", i2 = " + tokens[3]);
       info[1] = new String("Customization: \u03BB = " + tokens[4] + ", \u03B1 = " + tokens[5]); 
       info[2] = new String("Regularization: \u03BB_S = " + tokens[6] + ", \u03BB_D = " + tokens[7] + " \u03BB_C = " + tokens[8]);
       info[3] = new String("Pseudo AIC: " + tokens[9] + ", Effective DOF = " + tokens[10] + ", ZPC-Gene = " + plot_hybrid);    

       display = true;
      }
      else
      {display = false;}
   }

   public void setForeExt(double[] ex, boolean forecast)
   {System.arraycopy(ex, 0, this.extend, 0, 36); plotForecast = forecast;}

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,j,N,k;
     
     int t0, t1, x0, x1;
     int nobsP;  
     int p = 12+L-1;
     N = xf_obs;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;


     orig = new BasicStroke((float)1.7);
     BasicStroke thick = new BasicStroke((float)2.0);
     BasicStroke skinny = new BasicStroke((float)1.3);
 
     g2d.setStroke(skinny);
     g2d.setPaint(new Color(135,110,203));

     
     if(shade)
     {
        t0 = (int)((double)((double)(N-shade_points)/(double)N)*(double)width);
        t1 = width;
     
        g2d.setPaint(new Color(7,9,34));
        rectangle = new Rectangle(t0, 0, t1-t0, height);
        g2d.draw(rectangle);  
        g2d.fill(rectangle);
           
        g2d.setPaint(new Color(14,38,128));
        g2d.setStroke(orig);
        g2d.drawLine(t0, 0, t0, height);
        g2d.drawLine(t1, 0, t1, height);
         
     }
     
     
     
     
     
     if(zeroLine)
     {
      g2d.setStroke(dashed);
      x0 = (int)(((0.0 - dataMin)/dataNorm)*(double)height);
      g2d.drawLine(0, height - x0, width, height - x0);
     }

     g2d.setStroke(orig);
     g2d.setPaint(myGray2);
     g2d.setStroke(dashed);

     for(i=0; i < 9; i++)
     {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
     }
     g2d.setPaint(myGray);
     g.drawString((String)df.format(dataMax), 5, 15);
     g.drawString((String)df.format(dataMin), 5, height - 24);
     nobsP = (int)Math.floor((double)N/12); 
     for(i=1; i <= nobsP; i++) 
     {
        t0 = (int)(((double)(i*12)/N)*(double)width);
        g2d.setPaint(myGray2);
	g2d.drawLine(t0, 0, t0, height-20);
        g2d.setPaint(myGray);
        g.drawString((String)"" + p, t0, height - 5); 
        p = p + 12;
     } 



     g2d.setStroke(orig);
     g2d.setPaint(plotGray);  

     //System.out.println("N-rep " + n_rep + "n_obs = " + N);

     for(k=0;k<n_rep;k++)
     {
       if(ts_plots[k+1])
       {
 
        if(lag < 0)
        {
         for(j = (-lag); j < N-1; j++)
         {
	  t0 = (int)(((double)(j+lag)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1+lag)/(double)N)*(double)width);
	  x0 = (int)(((tseries[n_obs*k + j + L-1] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((tseries[n_obs*k + j + L] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }
        }
        else
        {
         for(j = 0; j < N-1; j++)
         {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((tseries[n_obs*k + j + L-1] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((tseries[n_obs*k + j + L] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }
        }
        g2d.setPaint(new Color((k+1)*25,(k+1)*10+111,102));
       }
     }
     if(xf_plot)
     {

        if(lag >= 0)  //shift data to left
        {
         g2d.setPaint(new Color(64,224,208));
         for(j = lag; j < N-1; j++)
         {
	  t0 = (int)(((double)(j-lag)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1-lag)/(double)N)*(double)width);
	  x0 = (int)(((xf[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((xf[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }
        }
        else
        {

         g2d.setPaint(new Color(64,224,208));
         for(j = 0; j < N-1; j++)
         {
	  t0 = (int)(((double)(j)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((xf[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((xf[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
          
         }
        }


     }
     if(sym_plot) //---- if selected and available plots x_sym ---
     {
        int N_sym = x_sym.length; int _lag = Math.min(lag,0); 
        int shift_r = Math.max(L-sym_shift,0); int shift_l = Math.min(L-sym_shift,0);

        //System.out.println("Symmetric Concurrent");
        

        int ex = Math.min(0,win_shift);
        int win = Math.max(0,win_shift);
        //System.out.println("New series\n");
        g2d.setPaint(new Color(139,156,173));
        for(j = shift_r-_lag; j < N_sym-1+ex; j++)
        {
	  t0 = (int)(((double)(j-shift_r-shift_l+_lag+win)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1-shift_r-shift_l+_lag+win)/(double)N)*(double)width);
	  x0 = (int)(((x_sym[j-ex] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((x_sym[j+1-ex] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);

          //System.out.println(x_sym[j] + "  " + xf[j]);
        }
     }
     //System.out.println("");
 
     //plotForecast = true; 
     if(lag < 0 && plotForecast)
     {
         g2d.setPaint(new Color(128,140,140));         
         for(j = 0; j < (-lag); j++)
         {
         
	  t0 = (int)(((double)(N-1 +lag + j)/(double)N)*(double)width);
	  t1 = (int)(((double)(N   +lag + j)/(double)N)*(double)width);
	  x0 = (int)(((extend[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((extend[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }

     }

     //---------- Now plot historical--------------------------------

     if(plotHist)
     {
 
      for(k=0;k<n_hist;k++)
      { 
        g2d.setStroke(orig);
        double[] temp = historical.get(k);

        int _lag = (int)temp[0];
        g2d.setPaint(colorHist[k]); if(hlight_indicator == k) {g2d.setPaint(highlight); g2d.setStroke(thick);}
        if(_lag >= 0)  //shift data to left
        {         
         for(j = _lag; j < temp.length-2; j++)
         {
	  t0 = (int)(((double)(j   -_lag)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1 -_lag)/(double)N)*(double)width);
	  x0 = (int)(((temp[j+1] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((temp[j+2] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }
        }
        else
        {
         _lag = -_lag;
         for(j = 0; j < temp.length-2-_lag; j++)
         {
	  t0 = (int)(((double)(j+_lag)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1 +_lag)/(double)N)*(double)width);
	  x0 = (int)(((temp[j+1] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((temp[j+2] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
         }
        }
      }
     }

    

    //-------------- If set display info about filter---------
    if(display)
    {
     
     g.setFont(mono);
     g2d.setPaint(Color.GREEN);
     g.drawString(info[0], canvas_shift, 15);
     g.drawString(info[1], canvas_shift, 30);
     g.drawString(info[2], canvas_shift, 45);
     g.drawString(info[3], canvas_shift, 60);

    }
   }

 


}
