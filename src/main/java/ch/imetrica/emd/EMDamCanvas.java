package ch.imetrica.emd;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;

public class EMDamCanvas extends JPanel
{


	private static final long serialVersionUID = 1L;
	//------ Other stuff
    Graphics2D g2d;
    int height, width;
    double dataMax, dataMin, dataNorm;
    BasicStroke dashed;
    float[] dash1;
    Color myGray;
    DecimalFormat df;
    //---- EMD stuff -------
    int N,n_imfs;
    double[][] am;
    boolean[][] plot_am;

    public EMDamCanvas(int w, int h, int _N)
    {
     int i;
     this.height = h; this.width = w; 
     dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
     dash1 = new float[1]; dash1[0] = 10.0f;
     dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
     myGray = new Color(67,71,73);
    
     N = _N; n_imfs = 6;
     am = new double[6][N]; plot_am = new boolean[6][3];
     for(i=0; i < 6; i++) plot_am[i][0] = false;

     setBackground(Color.BLACK); df = new DecimalFormat("##.##"); 
     //setPreferredSize(new Dimension(w, h));
    }

    public void setNObs(int n) {N = n;}

    public void setAM(double[][] _am, int _n, int nmfs)
    {
       int i,m;
       N = _n; n_imfs = nmfs;

       am = _am;
       dataMax = -100000000.0; dataMin = 10000000.0; dataNorm = 0.0;
       for(m=0; m < nmfs; m++) 
       {
         for(i=0;i<N;i++)
         {
           if(am[m][i] < dataMin) dataMin = am[m][i];
           else if(am[m][i] > dataMax) dataMax = am[m][i];
         }
       }
       dataNorm = Math.abs(dataMax - dataMin);
       go();
    }
        
    public void updateAM(int i, boolean w) {plot_am[i][0] = w; go();}

    public void setPlotOptionsAM(boolean[][] _plot_am) {plot_am = _plot_am;}
  
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
        g.drawString((String)df.format(dataMax), 5, 15);
        g.drawString((String)df.format(dataMin), 5, height - 5);
        int nobsP = (int)Math.floor((double)N/12); 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*12)/N)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-20);
        }
 
      g2d.setStroke(new BasicStroke(1.0f));
      colx = (int)(240/n_imfs);      
      red = 253;
      Color grad = new Color(165,red,211);
      g2d.setPaint(grad);

    
      for(i=0; i < n_imfs; i++)
      {
        if(plot_am[i][0])
        {
          for(j = 0; j < N-1; j++)
	  {
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((am[i][j] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((am[i][j+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	  }
          red = red - colx;
          grad = new Color(165,red,211);
          g2d.setPaint(grad);    
        }
      }

    }


}

