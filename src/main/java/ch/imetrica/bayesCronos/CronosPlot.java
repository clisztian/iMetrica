package ch.imetrica.bayesCronos;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;

/*------------------------------------------------
  Canvas for plotting series
--------------------------------------------------*/
public class CronosPlot extends JPanel
{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int n_obs; //---total number of observations
    double[] t_series;  //raw time series n_rep x n_obs
    double[][] forecasts; 
    double[][] alpha; 
    double[][] factors;
    double[] residual;
    double[] predictive;

    
    int n_fcsts;
    int n_steps;
    int n_factors;
 
    ArrayList<double[]> complete_data;
    double dataMax, dataMin, dataNorm;
    boolean plot_me; 
    int n_rep; 
    int n_sym = 0;
    int width; int height;
    Graphics2D g2d;

    boolean plot_residuals, plot_series, plot_predictive;
    boolean[] plots; 
    boolean[] plot_fore;
    boolean[] plot_alpha;
    boolean[] plot_factors;
    boolean plot_tar;
    Color[] colorHist;

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
    PrintWriter out;
    String file;
    boolean printer;
    FileWriter fileWriter;

    public CronosPlot(int w, int h, int nobs, int nrep)
    {

      int i; 
      n_obs = nobs; n_rep = nrep; t_series = new double[144];
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      colorHist = new Color[15];
      
      mono  = new Font("Monospaced", Font.PLAIN, 12);
      setBackground(Color.BLACK);
    
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
      colorHist[10] = new Color(200,41,2);    colorHist[11] = new Color(239,200,166); 
      colorHist[12] = new Color(100,73,169);  colorHist[13] = new Color(135,100,167);
      colorHist[14] = new Color(135,197,10);

      printer = false;
      plots = new boolean[10]; for(i=0;i<10;i++) {plots[i] = false;} plots[0] = true;
 
      plot_fore = new boolean[10]; for(i=0;i<10;i++) {plot_fore[i] = false;} 
   
      plot_alpha = new boolean[10]; for(i=0;i<10;i++) {plot_alpha[i] = false;} 

      plot_factors = new boolean[10]; for(i=0;i<10;i++) {plot_factors[i] = false;} 

      plot_residuals = false; plot_series = true; plot_predictive = false;
      
    }


    public void changeHighlight(int c)
    {hlight_indicator = c; go();}

    public void setplot(int i, boolean sel)
    {plots[i] = sel; computeDataMax(); go();}
  
    public void setplot_fore(int i, boolean sel)
    {plot_fore[i] = sel; go();}

    public void setplot_alpha(int i, boolean sel)
    {plot_alpha[i] = sel; computeDataMax(); go();}

    public void setplot_factors(int i, boolean sel)
    {plot_factors[i] = sel; computeDataMax(); go();}

    public void computeDataMax()
    {
      int i,k; dataMax = -1000000; dataMin = 10000000;

      //------- the time series--------- 
      for(k = 0; k < n_rep; k++)
      {
       for(i=0; i < n_obs; i++)
       {
          if(plots[k])
          {
            if(t_series[k*n_obs+i] > dataMax) dataMax = t_series[k*n_obs+i];
            else if(t_series[k*n_obs+i] < dataMin) dataMin = t_series[k*n_obs+i];   
          }
       }
      }
 
      //----- the alphas -------------
      for(k = 0; k < n_factors; k++)
      {
       for(i=0; i < n_obs; i++)
       {
          if(plot_alpha[k])
          {
            if(alpha[k][i] > dataMax) dataMax = alpha[k][i];
            else if(alpha[k][i] < dataMin) dataMin = alpha[k][i];   
          }

          if(plot_factors[k])
          {
            if(factors[k][i] > dataMax) dataMax = factors[k][i];
            else if(factors[k][i] < dataMin) dataMin = factors[k][i];   
          }

       }
      }      

      if(plot_residuals)
      {
        for(i=0; i < n_obs; i++)
        {
            if(residual[i] > dataMax) dataMax = residual[i];
            else if(residual[i] < dataMin) dataMin = residual[i];   
        }
      }

      if( plot_predictive )
      {
        for(i=0; i < n_obs; i++)
        {
            if(predictive[i] > dataMax) dataMax = predictive[i];
            else if(predictive[i] < dataMin) dataMin = predictive[i];   
        }
      }

      dataNorm = Math.abs(dataMax - dataMin) + .3*Math.abs(dataMax - dataMin);    
    }

    
    

    public void setData(double[] ts, int nrep, int nobs)
    {
       int i;
       n_rep = nrep; 
       n_obs = nobs; 
       //System.out.println("NREP data Changed = " + n_rep);
 
       t_series = new double[n_obs*n_rep]; 
       
       System.arraycopy(ts, 0, t_series, 0, n_obs*n_rep);
       computeDataMax();
 
       for(i=0;i<n_rep;i++) plots[i] = true;
 
       go();
    }
  
    public void setResidual(double[] ts, int nrep, int nobs)
    {
       n_obs = nobs; 
       residual = new double[n_obs];        
       System.arraycopy(ts, 0, residual, 0, n_obs);      
    }  
    
    public void setAlpha(double[] alphas, int n_fact, int nobs)
    {
       int i,k;
       n_factors = n_fact; 
       n_obs = nobs; 

       alpha = new double[n_fact][n_obs];

       for(i=0;i<n_factors;i++)
       {
         for(k=0;k<n_obs;k++)
         {
           alpha[i][k] = alphas[i*n_obs + k];
         }
       }
    }

    public void setFactors(double[] fact, int n_fact, int nobs)
    {
       int i,k;
       n_factors = n_fact; 
       n_obs = nobs; 

       factors = new double[n_fact][n_obs];

       for(i=0;i<n_factors;i++)
       {
         for(k=0;k<n_obs;k++)
         {
           factors[i][k] = fact[i*n_obs + k];
         }
       }
    }
   
   
    public void setPredictive(double[] fcasts, int nobs)
    {
      //System.out.println("SETS PREDICTIVE = " + fcasts.length);
      predictive = new double[nobs];
      System.arraycopy(fcasts, 0, predictive, 0, nobs);
    }

    public void plotPredictive(boolean pred)
    {
      plot_predictive = pred; computeDataMax(); go();
    }

 
    public void turnOnPrinter()
    {printer = true; go();}

  
    public void setForecasts(double[] fcts, int steps, int n_fore)
    {
      int i,j;
      n_steps = steps;
      n_fcsts = n_fore;
      forecasts = new double[n_fcsts][n_steps];
 
      for(i = 0; i < n_fcsts; i++)
      {
        for(j = 0; j < n_steps; j++)
        {
          forecasts[i][j] = fcts[i*n_steps + j];
        }
      }   
      go();  
    }
    
    void plotResiduals(boolean t) 
    {
      
       plot_residuals = t; computeDataMax();
       go(); 
    }

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,j,N,k;
     
     int t0, t1, x0, x1;
     int nobsP;  int p = 10;

     
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;


     BasicStroke orig = new BasicStroke((float)1.5);
     BasicStroke thick = new BasicStroke((float)1.7);
     new BasicStroke((float)1.1);
     new BasicStroke((float)1.5);

        //Stroke def = g2d.getStroke();

       //--------------------- Draw grid with or without forecast points ------------
     float[] dash1 = {10.0f};
     BasicStroke dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);

     N = n_obs+n_steps;
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
     nobsP = (int)Math.floor((double)n_obs/10); 
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
     double[] savedData = new double[n_obs];
 

     

    //System.out.println("plots, nrep  = " + plots[0] + "  " + plots[1] + " " + n_rep); 

    if(plot_series)
    {
     for(k=0;k<n_rep;k++)
     {
       if(plots[k])
       {
         //System.out.println("plots = " + k);
         g2d.setPaint(colorHist[k]); if(hlight_indicator == (k+1)){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < n_obs-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((t_series[k*n_obs+j] - dataMin)/dataNorm)*(double)height*1.0);
          x1 = (int)(((t_series[k*n_obs+j + 1] - dataMin)/dataNorm)*(double)height*1.0);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);

          if(hlight_indicator == (k+1) && printer)
          {
            savedData[j] = t_series[k*n_obs+j];
          }
         }
         if(hlight_indicator == (k+1) && printer)
         {
            savedData[j] = t_series[k*n_obs+n_obs-1];
         }
       }
     }
     
     g2d.setStroke(orig);
     
     //System.out.println(n_fcsts + " " + n_steps);
     for(k=0;k<n_fcsts;k++)
     {
        if(plot_fore[k])
        {

          g2d.setPaint(colorHist[colorHist.length-1-k]); //if(hlight_indicator == k){g2d.setPaint(highlight); g2d.setStroke(thick);}

          t0 = (int)(((double)(n_obs-1)/(double)(N-1))*(double)width);
          t1 = (int)(((double)(n_obs)/(double)(N-1))*(double)width); 
          x0 = (int)(((t_series[n_obs-1]  - dataMin)/dataNorm)*(double)height*1.0);
          x1 = (int)(((forecasts[k][0] - dataMin)/dataNorm)*(double)height*1.0);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);


       for(j = 0; j < n_steps-1; j++)
       {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)(n_obs+j)/(double)(N-1))*(double)width);
          t1 = (int)(((double)(n_obs+j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((forecasts[k][j] - dataMin)/dataNorm)*(double)height*1.0);
          x1 = (int)(((forecasts[k][j+1] - dataMin)/dataNorm)*(double)height*1.0);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          
         // System.out.println(forecasts[k][j]);
       }
      }
     }
     
          
        
      for(k=0;k<n_factors;k++)
      {
        if(plot_alpha[k])
        {
         
         g2d.setPaint(colorHist[k+1]); 
         if(hlight_indicator == (6+k)){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < n_obs-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j])
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((alpha[k][j] - dataMin)/dataNorm)*(double)height*1.0);
          x1 = (int)(((alpha[k][j + 1] - dataMin)/dataNorm)*(double)height*1.0);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          //System.out.println(alpha[j]);
          if(hlight_indicator == (k+6) && printer)
          {
            savedData[j] = alpha[k][j];
          }
         }
         if(hlight_indicator == (k+6) && printer)
         {
            savedData[n_obs-1] = alpha[k][n_obs-1];
         }         
        }
      }


      for(k=0;k<n_factors;k++)
      {
        if(plot_factors[k])
        {
         
         g2d.setPaint(colorHist[k+7]); 
         if(hlight_indicator == (9+k)){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < n_obs-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j])
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((factors[k][j] - dataMin)/dataNorm)*(double)height*1.0);
          x1 = (int)(((factors[k][j + 1] - dataMin)/dataNorm)*(double)height*1.0);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          //System.out.println(alpha[j]);
          if(hlight_indicator == (k+9) && printer)
          {
            savedData[j] = factors[k][j];
          }
         }
         if(hlight_indicator == (k+9) && printer)
         {
            savedData[n_obs-1] = factors[k][n_obs-1];
         }          
        }
      }  

     
     
     }
     
     if(plot_residuals)
     {
       
         k=0;  
         g2d.setPaint(colorHist[1]); 
         if(hlight_indicator == 12){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < n_obs-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((residual[k*n_obs+j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((residual[k*n_obs+j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          if(hlight_indicator == 12 && printer)
          {
            savedData[j] = residual[j];
          }
         }
         if(hlight_indicator == 12 && printer)
         {
            savedData[n_obs-1] = residual[n_obs-1];
         }  
     } 


     if(plot_predictive)
     {
         
         g2d.setPaint(colorHist[7]);
         if(hlight_indicator == 13){g2d.setPaint(highlight); g2d.setStroke(thick);}
         for(j = 0; j < predictive.length-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
          t0 = (int)(((double)j/(double)(N-1))*(double)width);
          t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width); 
          x0 = (int)(((predictive[j] - dataMin)/dataNorm)*(double)height);
          x1 = (int)(((predictive[j + 1] - dataMin)/dataNorm)*(double)height);
          g2d.drawLine(t0, (height) - x0, t1, (height) - x1);
          if(hlight_indicator == 13 && printer)
          { 
            if(j >= n_steps)
            {savedData[j-n_steps] = predictive[j];}
          }
         }
         if(hlight_indicator == 13 && printer)
         {
            savedData[n_obs-1] = predictive[predictive.length-1];
         }          
     }


     if(printer && hlight_indicator > 0)
     {

       try
       {
         file = "iMetrica-BayesCronos-ts_"+hlight_indicator+".dat";       
         fileWriter = new FileWriter(file);
         out = new PrintWriter(fileWriter);
   
         for(i=0;i<n_obs;i++) {out.println(savedData[i]);}
         out.close();
         System.out.println("Data successfully saved in " + file);
       }
       catch (IOException e) {e.printStackTrace();} 
     }
     
     printer = false;
    }
       
}
