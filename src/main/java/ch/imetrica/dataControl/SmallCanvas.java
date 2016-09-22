package ch.imetrica.dataControl;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

/*------------------------------------------------
  Canvas for plotting miniature series
--------------------------------------------------*/
public class SmallCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int n_obs; //---total number of observations
    double[] tseries;  //raw time series n_rep x n_obs

    double dataMax, dataMin, dataNorm;
    boolean plot_me; 
    int n_series;
    int width; int height;
    Graphics2D g2d;
 
    Color[] colors; 
    int colorx;


    public SmallCanvas(int w, int h, int nobs, int x)
    {
      // Initilize everything-------------------
      n_obs = nobs; colorx = x;
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;
      plot_me = false; n_series = 1;
      colors = new Color[8]; 

      colors[0] = new Color(117,159,159); colors[1] = new Color(117,159,159);
      colors[2] = new Color(117,159,159); colors[3] = new Color(117,159,159);
      colors[4] = new Color(180,115,208); colors[5] = new Color(4,237,208);
      colors[6] = new Color(4,237,230); colors[7] = new Color(4,210,208);

      tseries = new double[n_obs];
      setBackground(Color.BLACK);
      setPreferredSize(new Dimension(w, h));  
    }
     
    public void setPlotDim(double[] series, int N)
    {
        int i;
        dataMax = -1000000.0; dataMin = 1000000.0;
        for(i=0;i<N;i++)
        {
         if(series[i] < dataMin) dataMin = series[i];
         else if(series[i] > dataMax) dataMax = series[i]; 
        } 
        dataNorm = Math.abs(dataMax - dataMin)*1.20;    
    }
  
    public void clear()
    {plot_me = false; go();}

    public void setNSeries(int n) {n_series = n;}
    public void setNobs(int n) {n_obs = n;}
    

    //----- set raw data and filtered data
    public void setTseries(double[] ts, int N, int rep)
    {   
       plot_me = true;  
       n_obs = N; n_series = rep;
       tseries = new double[rep*N];

       if(rep==1){System.arraycopy(ts, 0, tseries, 0, n_obs);}
       else if(rep>1)
       {
         System.arraycopy(ts, 0, tseries, 0, ts.length);
       }  
       this.setPlotDim(tseries,rep*N);
    }
     


    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

       int j,N,k;  int t0, t1, x0, x1;

       N = n_obs;
       super.paintComponent(g); g2d = (Graphics2D)g; 
   
       Dimension ds = this.getSize();
       width = ds.width; height = ds.height;

       BasicStroke orig = new BasicStroke((float)1.0);
     
       g2d.setStroke(orig);
       g2d.setPaint(colors[colorx]);     
   
       if(plot_me)
       {     
        for(j = 0; j < N-1; j++)
        {
         t0 = (int)(((double)j/(double)N)*(double)width); t1 = (int)(((double)(j+1)/(double)N)*(double)width);
         x0 = (int)(((tseries[j] - dataMin)/dataNorm)*(double)height); x1 = (int)(((tseries[j+1] - dataMin)/dataNorm)*(double)height);
         g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);           
        }
 
        g2d.setPaint(colors[4]); 
        if(n_series > 1)
        {
         for(k=1;k<n_series;k++)
         {
          g2d.setPaint(colors[k]);
          for(j = 0; j < N-1; j++)
          {
           t0 = (int)(((double)j/(double)N)*(double)width); t1 = (int)(((double)(j+1)/(double)N)*(double)width);
           x0 = (int)(((tseries[N*k+j] - dataMin)/dataNorm)*(double)height); 
           x1 = (int)(((tseries[N*k+j+1] - dataMin)/dataNorm)*(double)height);
           g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);       
          }
         }
        }
       }
    }    
}