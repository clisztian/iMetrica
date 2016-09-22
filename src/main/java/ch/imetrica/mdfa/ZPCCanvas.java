package ch.imetrica.mdfa;

/*!
Copyright (C) 2016 Christian D. Blakely
This file is part of iMetrica, a free-software/open-source application
for interactive graphical econometric analysis - http://imetricablog.com/

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;

/*------------------------------------------------
  Canvas for plotting mdfa tseries data on time domain 
  - Only accepts data of form  tseries by row
--------------------------------------------------*/
public class ZPCCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int K,K1; //---total number of observations
    double[] Gamma; //gamma K+1
    double[] amp; 
    double[] time_delay; 
    double[] phase;
    boolean[] plot;

    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;

    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;

    public ZPCCanvas()
    {
      // Initilize everything-------------------
      dataMax = -1000000.0; dataMin = 1000000.0;

      setBackground(Color.BLACK);
      //setPreferredSize(new Dimension(w, h));    

      myGray = new Color(67,71,73);
      myGray2 = new Color(20,21,19);
      plotGray = new Color(92,172,238);
      df = new DecimalFormat("##.##"); 
      plot = new boolean[4];
      plot[0] = true; plot[1] = true; plot[2] = false; plot[3] = false; 

    } 

    public void setPlotDim()
    {
        int i; int K1 = K+1;
        dataMax = -100000.0; dataMin = 100000.0; 

        for(i=0;i<K1;i++)
        {
          
           if(Gamma[i] < dataMin) dataMin = Gamma[i];
           if(amp[i] < dataMin) dataMin = amp[i]; 
           if(phase[i] < dataMin) dataMin = phase[i];
           if(time_delay[i] < dataMin) dataMin = time_delay[i];
       
           if(Gamma[i] > dataMax) dataMax = Gamma[i];
           if(amp[i] > dataMax) dataMax = amp[i]; 
           if(phase[i] > dataMax) dataMax = phase[i];
           if(time_delay[i] > dataMax) dataMax = time_delay[i];
        }
        dataNorm = Math.abs(dataMax - dataMin);     
    
        dataMax = 1.5; dataMin = -2.0; dataNorm = Math.abs(dataMax - dataMin);
           
    }

    public void setK(int _k)
    {K1 = _k+1; Gamma = new double[K1];}

    public void setPlotDraw(boolean s, int i)
    {plot[i] = s; go();}


    public void setGamma(double[] _gamma)
    {
       K1 = _gamma.length;
       Gamma = new double[_gamma.length];
       System.arraycopy(_gamma,0,Gamma,0,K1);
       amp = new double[K1];
       phase = new double[K1];
       time_delay = new double[K1];
       //setPlotDim();
    }

    public void setZPCFilter(double[] a, double[] b, double[] c)
    {
        K1 = a.length;
        amp = new double[K1];
        phase = new double[K1];
        time_delay = new double[K1];

        System.arraycopy(a,0,amp,0,K1);
        System.arraycopy(b,0,phase,0,K1);
        System.arraycopy(c,0,time_delay,0,K1);
        setPlotDim(); go();
    }

    
    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N;
     int t0, t1, x0, x1;
     N = K1;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 

     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;

     Font font = new Font("Arial", Font.PLAIN, 10);
    
     g2d.setFont(font);
     g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);


     g2d.setPaint(Color.PINK);
     //g.setFont(sanSerifFont);
     g.drawString("0", 5, height-5);
     g.drawString("\u03C0/4", (int)width/4, height-5);
     g.drawString("\u03C0/2", (int)width/2, height-5);
     g.drawString("3*\u03C0/4", (int)3*width/4, height-5);
     g.drawString("\u03C0", width-10, height-5);
     
     g2d.setPaint(myGray);

     //System.out.println(dataMax + "  " + dataMin);
     g.drawString((String)df.format(dataMax), 5, 15);
     g.drawString((String)df.format(dataMin), 5, height - 24);


     orig = new BasicStroke((float)1.5);
     g2d.setStroke(orig);

   if(plot[0])
   {
     g2d.setPaint(new Color(239,177,177));
     for(j = 0; j < K1-1; j++)
     {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((Gamma[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((Gamma[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
     }
   }
   if(plot[1])
   {
     g2d.setPaint(new Color(18,177,239));
     for(j = 0; j < K1-1; j++)
     {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((amp[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((amp[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height- x1);
     }
   }
   if(plot[2])
   {
     g2d.setPaint(new Color(28,203,10));
     for(j = 0; j < K1-1; j++)
     {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((phase[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((phase[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
     }
   }
   if(plot[3])
   {
     g2d.setPaint(new Color(246,102,214));
     for(j = 0; j < K1-1; j++)
     {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);         
	  x0 = (int)(((time_delay[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((time_delay[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, height - x0, t1, height - x1);
     }
   }           
  }

}
