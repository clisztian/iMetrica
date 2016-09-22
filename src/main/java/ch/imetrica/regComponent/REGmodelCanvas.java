package ch.imetrica.regComponent;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;

/*------------------------------------------------
  Canvas for plotting mdfa tseries data on time domain 
  - Only accepts data of form  tseries by row
--------------------------------------------------*/
public class REGmodelCanvas extends JPanel
{

   /**
	 * 
	 */
   private static final long serialVersionUID = 1L;
   public int n_obs,NFore;                 // number of observations, observations+forecasts
   public int S;
   public int lag;
   public int n_comps;                     // number of components
   public String date_string;              // starting date of data 
   public double[] series;                 // raw data 
   public boolean plotForecasts;           // plot forecasts

   double[] usimCmpnts;                    // all the sig ext components + forecasts
   double[] usimH;                         // all the scaling factors + forecasts
   double[] usimXbeta;                     // all the regression coponents + forecasts
   double[] usimXbetaH;                    // all the regression coponents*H + forecasts
   double[] usimQ;                         // mean squared errors  + forecasts
   double[] usimQbar;                      // mean squared error variances  + forecasts 

   double[] y_frcstm;                      // aggregate data forecasts - mid
   double[] y_frcstl;                      // aggregate data forecasts - low
   double[] y_frcsth;                      // aggregate data forecasts - high
   double[] y_xB_frcst;                    // aggregate data forecasts - regression
   boolean[] comp_plots;                  
   boolean series_plot;
 
   Color highlight; 
   int hlight_indicator;
   Color seasonalColor;
   Color trendColor; 
   Color saColor; 
   Color cycleColor;  
   Color seriesColor;
   Color seriesColor2;
   Color myGray,myGray2;
   Color forecol; 
   boolean reg_computed; 
    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;

    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color plotGray;



    public REGmodelCanvas(int w, int h, int _N)
    {
      // Initilize everything-------------------
      int i;
      n_obs = _N; NFore=n_obs+24;
      n_comps = 0;
      hlight_indicator = -1;
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      reg_computed = false;

      setBackground(Color.BLACK);
      setPreferredSize(new Dimension(w, h));    
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
      myGray = new Color(67,71,73);
      myGray2 = new Color(20,21,19);
      plotGray = new Color(92,172,238);
      df = new DecimalFormat("##.##"); 

       seasonalColor = new Color(255,47,47);
       trendColor = new Color(201,105,168); 
       saColor = new Color(97,155,193); 
       cycleColor = new Color(185,109,249);  
       seriesColor = new Color(142,238,255);
       seriesColor2 = new Color(135,228,232);
       forecol = new Color(149,107,167);
  
      //---- initialize plots, set all to false------
      comp_plots = new boolean[10]; 
      for(i=0;i<10;i++) {comp_plots[i] = false;}
      series_plot = true; plotForecasts = false;

      y_frcstl = new double[24];
      y_frcstm = new double[24];
      y_frcsth = new double[24];
      y_xB_frcst = new double[24];

      hlight_indicator = -1;

    } 

    public void setPlotDim()
    {
        int i,k,N; 
        dataMax = -100000.0; dataMin = 100000.0; 
   
        N = n_obs+24; NFore=n_obs+24;

        if(reg_computed)
        {
         for(k=0;k<n_comps;k++)
         {
          if(comp_plots[k])
          {
            for(i=0;i<N;i++)
            {
             if(usimXbetaH[k*N+i] < dataMin) dataMin = usimXbetaH[k*N+i];
             else if(usimXbetaH[k*N+i] > dataMax) dataMax = usimXbetaH[k*N+i]; 
            } 
           } 
         }
        }
        if(series_plot)
        {
            for(i=0;i<n_obs;i++)
            {
             if(series[i] < dataMin) dataMin = series[i];
             else if(series[i] > dataMax) dataMax = series[i]; 
            } 
        } 
        dataNorm = Math.abs(dataMax - dataMin)*1.20;                
    }
 
    public void setForecasts(boolean fr) {plotForecasts = fr; go();}

    public void setPlots(int i, boolean sel)
    {
       if(i==0) 
       {series_plot = sel;} 
       else if(i <= n_comps+1)
       {comp_plots[i-1] = sel;}
       setPlotDim(); go();
    }   

    public void setN(int n) {n_obs = n; NFore=n_obs+24;}
    public void setNumberComps(int n) {n_comps = n;}

    public void setData(double[] d)
    { 
       reg_computed = false;
       n_obs = d.length; NFore=n_obs+24; series = new double[n_obs];
       System.arraycopy(d, 0, series, 0, n_obs);
       setPlotDim(); go();
    }

    public void setComponents(double[] d)
    {  
       reg_computed = true; usimXbetaH = new double[d.length];
       System.arraycopy(d, 0, usimXbetaH, 0, d.length); 
    }


    public void setForecasts(double[] d, double[] d1 , double[] d2, double[] d3)
    { 
      reg_computed = true;
      System.arraycopy(d,0,y_frcstl,0,24);
      System.arraycopy(d1,0,y_frcstm,0,24);
      System.arraycopy(d2,0,y_frcsth,0,24);
      System.arraycopy(d3,0,y_xB_frcst,0,24);

    }
  
    public void setQvariances(double[] d, double[] d1)   
    {
       reg_computed = true; usimQ = new double[d.length]; usimQbar = new double[d1.length];
       System.arraycopy(d, 0, usimQ, 0, d.length); 
       System.arraycopy(d1, 0, usimQbar, 0, d1.length);
    }
    
    public void changeHighlight(int c)
    {hlight_indicator = c;}

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,k,p;
     int t0, t1, x0, x1;
     int nObs = n_obs;  

     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
 
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     BasicStroke bold = new BasicStroke((float)1.8);
     BasicStroke def = new BasicStroke((float)1.5);
     BasicStroke dashed = new BasicStroke(1.0f, 
                                          BasicStroke.CAP_BUTT, 
                                          BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);

     g2d.setStroke(dashed);
     g2d.setPaint(myGray);

     //-------------------- Plot dashed lines first -----------------


     int nObsex = nObs;
     if(!plotForecasts) 
     {   

        int nobsP = (int)Math.floor((double)nObs/12);  p = 12; 
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
        nObsex = nObs+24;       
        int nobsP = (int)Math.floor((double)nObsex/12); p = 12; 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*12)/nObsex)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-10);
          g.drawString((String)"" + p, t0, height - 5);
          p = p + 12;
        } 
     }


      g2d.setStroke(def);
			
      if(reg_computed){
      if(series_plot)
      {
	  g2d.setPaint(seriesColor); g2d.setStroke(def);
	  if(!plotForecasts) 
          { 
           for(i = 0; i < nObs-1; i++)
	   {
	    t0 = (int)(((double)i/(double)nObs)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObs)*(double)width);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	   }
          }
          else if(plotForecasts)
          {
           g2d.setStroke(def);

           //if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}

           for(i = 0; i < nObs-1; i++)
	   {            	
	    t0 = (int)(((double)i/(double)NFore)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)NFore)*(double)width);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	   }	

            i=nObs-1; g2d.setPaint(seriesColor2); 
            g2d.setStroke(def); 
            //if(hlight_indicator == 1){g2d.setPaint(highlight);  g2d.setStroke(bold);}

	    t0 = (int)(((double)i/(double)NFore)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)NFore)*(double)width);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((y_frcstl[0] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            g2d.setPaint(forecol);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((y_frcstm[0] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            g2d.setPaint(forecol);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
 	    x1 = (int)(((y_frcsth[0] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

           for(i = 0; i < 23; i++)
	   {            	
             g2d.setPaint(forecol); g2d.setStroke(def);
	     t0 = (int)(((double)(i+nObs)/(double)NFore)*(double)width);
	     t1 = (int)(((double)(i+nObs+1)/(double)NFore)*(double)width);
	     x0 = (int)(((y_frcstm[i] - dataMin)/(dataNorm*1.20))*(double)height);
	     x1 = (int)(((y_frcstm[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	     g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	     g2d.setPaint(seriesColor2); g2d.setStroke(def); //if(hlight_indicator == 1){g2d.setPaint(highlight);}
             x0 = (int)(((y_frcstl[i] - dataMin)/(dataNorm*1.20))*(double)height);
	     x1 = (int)(((y_frcstl[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	     g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1); 
             g2d.setPaint(forecol); g2d.setStroke(def);
	     x0 = (int)(((y_frcsth[i] - dataMin)/(dataNorm*1.20))*(double)height);
	     x1 = (int)(((y_frcsth[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	     g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
           }
          }	
        }


        if(plotForecasts){nObs = NFore;}
        else {nObs = n_obs;} 

        for(k=0; k < n_comps; k++)
        {
          g2d.setPaint(new Color(15+(k+1)*40,102,255));
          if(comp_plots[k])
          {
           if(hlight_indicator == k){g2d.setPaint(highlight);  g2d.setStroke(bold);}
           for(i = 0; i < nObs-1; i++)
    	   {

	    t0 = (int)(((double)i/(double)nObsex)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObsex)*(double)width);

            x0 = (int)(((usimXbetaH[NFore*k + i]- dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((usimXbetaH[NFore*k + i+1]- dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);      
           }
          }
        }
      }
      else  // just plot the series
      {

	  g2d.setPaint(seriesColor); g2d.setStroke(def);
          for(i = 0; i < n_obs-1; i++)
	  {
	    t0 = (int)(((double)i/(double)nObs)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)nObs)*(double)width);
	    x0 = (int)(((series[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((series[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	  }
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

        

        hlight_indicator = -1;
    }
}
