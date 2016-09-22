package ch.imetrica.mdfa;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;

             //tcimfPanel.setIMFData(mdfa.amMap, mdfa.fmMap, n_imfs, flength);
             //tcimfPanel.setTrend(mdfa.trend_cycle);
             //tcimfPanel.setIMF(0,K+1);


public class IMFPlotPanel extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//------ Other stuff
    Graphics2D g2d;
    int height, width;
    double dataMax, dataMin, dataNorm;
    double tdataMax, tdataMin, tdataNorm;
    BasicStroke dashed;
    float[] dash1;
    Color myGray;
    DecimalFormat df;
    //---- EMD stuff -------
    boolean plotok;
    int plot_imf;
    int N,n_imfs,L;
    double[][] imfs;
    double[] trend; 
    double[] xf;
    Font font;
    boolean plot_data = true;
    int plot_type = 0;

    public IMFPlotPanel()
    {

     dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
     tdataMax = -100000000.0; tdataMin = 10000000.0; tdataNorm = 0.0;
     dash1 = new float[1]; dash1[0] = 5.0f;
     dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          5.0f, dash1, 0.0f); 
     myGray = new Color(37,41,43);
     plotok = false; n_imfs = 11;  font = new Font("serif", Font.PLAIN, 7);
     setBackground(Color.BLACK); df = new DecimalFormat("##.##"); 
     setPreferredSize(new Dimension(775, 90));
    }


    public void setIMFData(double[][] _am, double[][] _fm, int nmfs, int _n)
    {
       int i,m;
       N = _n; n_imfs = nmfs;
       imfs = new double[n_imfs][N];
   
      if(plot_type == 0) //----------------------------------------------------
      {     
       for(m=0; m < nmfs; m++) 
       {
         for(i=0;i<N;i++)
         {imfs[m][i] = _am[m][i]*_fm[m][i];}
       }       
      
       for(m=0; m < nmfs; m++) 
       {
         for(i=0;i<N;i++)
         {
           if(imfs[m][i] < dataMin) dataMin = imfs[m][i];
           else if(imfs[m][i] > dataMax) dataMax = imfs[m][i];
         }
       }
      }
      else  //----------------------------------------------------
      {
       for(m=0; m < nmfs; m++) 
       {
         for(i=0;i<N;i++)
         {imfs[m][i] = _am[m][i];}
       }       
      
       for(m=0; m < nmfs; m++) 
       {
         for(i=0;i<N;i++)
         {
           if(imfs[m][i] < dataMin) dataMin = imfs[m][i];
           else if(imfs[m][i] > dataMax) dataMax = imfs[m][i];
         }
       }
      }

       dataNorm = Math.abs(dataMax - dataMin);
       go();
    }
        
    public void setL(int l) {L = l;}
    public void setData(double[] _xf)
    {
       int i;
       xf = new double[_xf.length];
       dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
       System.arraycopy(_xf, 0, xf, 0, _xf.length);
       N = _xf.length;
         for(i=0;i<N;i++)
         {
           if(xf[i] < dataMin) dataMin = xf[i];
           else if(xf[i] > dataMax) dataMax = xf[i];
         }
       dataNorm = Math.abs(dataMax - dataMin);
    } 
    public void plotData()
    {
       plot_data = !plot_data;   
       go();
    }

    public void plotType(int x)
    {plot_type = x;}


    public void setTrend(double[] _trend)
    {
       int i;
       trend = new double[_trend.length];
       tdataMax = -100000000.0; tdataMin = 10000000.0; tdataNorm = 0.0;
       System.arraycopy(_trend, 0, trend, 0, _trend.length);
       N = _trend.length;
         for(i=0;i<N;i++)
         {
           if(trend[i] < tdataMin) tdataMin = trend[i];
           else if(trend[i] > tdataMax) tdataMax = trend[i];
         }
       tdataNorm = Math.abs(tdataMax - tdataMin);
    }

    public void setIMF(int val, int K, double[] freq_ints)
    {
       plotok = true;
       double pi = Math.PI/K;
       plot_imf = 0;
       //if(val*pi < freq_ints[2]) {plot_imf = 0;}
       for(int k = 2; k < freq_ints.length; k++)
       {
          if((val*pi <= freq_ints[k]) && (val*pi > freq_ints[k-1]))
          {
             plot_imf = k-1; break;
          }
       }
       if(val*pi > freq_ints[freq_ints.length-1]) {plot_imf = freq_ints.length - 2;}
       go();
    }
  
    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {     
      int i,j,colx,red;
      super.paintComponent(g);
      int t0, t1, x0, x1;
      g2d = (Graphics2D)g;
      
      Dimension ds = this.getSize();
      width = ds.width; height = ds.height;  

      

        //Draw dashed lines 
        g2d.setStroke(dashed);
        g2d.setPaint(myGray);

        for(i=0; i < 8; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }
        
        int nobsP = (int)Math.floor((double)300/20); 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*20)/300)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-20);
        }

    if(plotok)
    { 
     
      g2d.setStroke(new BasicStroke(1.0f));
      colx = (int)(140/n_imfs);      
      red = 253; int p;
      Color grad = new Color(54,182,219);
     
      if(plot_data)
      {
          g2d.setPaint(new Color(144,88,137));
          for(j = 0; j < N-1; j++)
	  {  
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((xf[j] - dataMin)/(dataNorm*1.10))*(double)height);
	    x1 = (int)(((xf[j+1] - dataMin)/(dataNorm*1.10))*(double)height);
	    g2d.drawLine(t0, height - x0, t1, height - x1);
	  }
      }
      

       
      if(plot_imf > 0)
      { 
          g2d.setPaint(myGray);
          g.drawString((String)df.format(dataMax), 5, 15);
          g.drawString((String)df.format(dataMin), 5, height - 5);
          red = red - colx;
          //grad = new Color(165,red,211);
          g2d.setPaint(grad);    
          for(j = 0; j < N-1; j++)
	  {
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((imfs[plot_imf-1][j] - dataMin)/(dataNorm*1.10))*(double)height);
	    x1 = (int)(((imfs[plot_imf-1][j+1] - dataMin)/(dataNorm*1.10))*(double)height);
	    g2d.drawLine(t0, height - x0, t1, height - x1);
	  }

     }
     else 
     {    
           g2d.setPaint(myGray); g2d.setFont(font);
           g.drawString((String)df.format(dataMax), 5, 15);
           g.drawString((String)df.format(dataMin), 5, height - 5);
          //grad = new Color(165,red,211);
          g2d.setPaint(grad);     
          for(j = 0; j < N-1; j++)
	  {  
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((trend[j] - dataMin)/(dataNorm*1.10))*(double)height);
	    x1 = (int)(((trend[j+1] - dataMin)/(dataNorm*1.10))*(double)height);
	    g2d.drawLine(t0, height - x0, t1, height - x1);
	  }

     }

        g2d.setFont(font); p = 28;        
        nobsP = (int)Math.floor((double)N/20); 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*20)/N)*(double)width);
          //g2d.setPaint(myGray);
	  //g2d.drawLine(t0, 0, t0, height-20);
          g2d.setPaint(myGray);
          g.drawString((String)"" + p, t0, height - 2); 
          p = p + 20;
        } 



    }




   }


}