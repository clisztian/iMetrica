package ch.imetrica.mdfa;

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
public class CrystalBallcanvas extends JPanel
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
    double[] in_signal;       //filtered series n_rep x (n_obs - L + 1)
    double[] out_signal;
    
 
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

    public CrystalBallcanvas(int w, int h, int nobs, int nrep, int _L)
    {
      n_obs = nobs; n_rep = nrep; L = _L; lag = 0;
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      //mono = new Font(parent.getDisplay(), "Monospaced", 10, SWT.NONE);
      mono  = new Font("Monospaced", Font.PLAIN, 12);
      setBackground(Color.BLACK);
      //setBackground(new Color(255,250,250));
      setPreferredSize(new Dimension(w, h));    
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
      myGray = new Color(67,71,73);
      myGray2 = new Color(50,51,59);
      plotGray = new Color(145,63,63);
      df = new DecimalFormat("##.##"); 
      foreRed = new Color(232, 138 , 187);
      highlight = new Color(255,255,255);

 
      tseries = new double[n_obs];
      in_signal = new double[n_obs];
      out_signal = new double[n_obs];
      
    } 


    public void setPlotDim()
    {
        int i; 
        dataMax = -1000000.0; dataMin = 1000000.0;
        for(i=0;i<n_obs;i++)
        {
         if(tseries[i] < dataMin) dataMin = tseries[i];
         else if(tseries[i] > dataMax) dataMax = tseries[i]; 
        } 
        dataNorm = Math.abs(dataMax - dataMin);    
    }
  

    public void setSeries(double[] ts, double[] sig1, double[] sig2)
    {   
      if((ts.length == sig1.length) && (sig1.length == sig2.length))
      {
       n_obs = ts.length;
       tseries = new double[ts.length];
       in_signal = new double[ts.length];
       out_signal = new double[ts.length];
      
       System.arraycopy(ts,0,tseries,0,ts.length);
       System.arraycopy(sig1,0,in_signal,0,sig1.length);
       System.arraycopy(sig2,0,out_signal,0,sig2.length);
      }
    }
    
    
    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N;
     
     int t0, t1, x0, x1;
     N = n_obs;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;


     orig = new BasicStroke((float)1.5);
     new BasicStroke((float)2.0);
     BasicStroke skinny = new BasicStroke((float)1.2);
 
     g2d.setStroke(skinny);
     g2d.setPaint(new Color(135,110,203));

     setPlotDim();
     
     
    
      g2d.setStroke(dashed); g2d.setPaint(myGray);
      x0 = (int)(((0.0 - dataMin)/dataNorm)*(double)height);
      g2d.drawLine(0, height - x0, width, height - x0);
     


      g2d.setStroke(orig);
      g2d.setPaint(new Color(181,44,44));  

      for(j = 0; j < N-1; j++)
      {
	  t0 = (int)(((double)(j)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((tseries[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((tseries[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
      }

      g2d.setStroke(orig);
      g2d.setPaint(new Color(213,132,255));       

      
      for(j = 0; j < N-1; j++)
      {
	  t0 = (int)(((double)(j)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((in_signal[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((in_signal[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
      }

      g2d.setStroke(orig);
      g2d.setPaint(new Color(144,174,231));           
      
      for(j = 0; j < N-1; j++)
      {
	  t0 = (int)(((double)(j)/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((out_signal[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((out_signal[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
      }      
        
 
   }


}
