package ch.imetrica.mdfa;
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
public class SpectralDensityCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int K; //---total number of observations
    int K1;
    int n_rep; 

    double[] per,mod,arg;

    boolean[] period_plots;   //---- periodograms plot   
 
    boolean plot_per, plot_mod, plot_arg;  


    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;
    
    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;
   

    public SpectralDensityCanvas(int w, int h, int _K, int _nrep)
    {
      // Initilize everything-------------------
      int i;
      K = _K; K1 = _K+1; 
      n_rep = _nrep;
      
      this.width = w; this.height = h; 
      dataMax = -1000000.0; dataMin = 1000000.0;

      setBackground(Color.BLACK);
      setPreferredSize(new Dimension(w, h));    
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
      myGray = new Color(67,71,73);
      myGray2 = new Color(20,21,19);
      plotGray = new Color(92,172,238);
      df = new DecimalFormat("##.##"); 

      //---- initialize plots, set all to false------
      period_plots = new boolean[9]; 
      for(i=0;i<9;i++) {period_plots[i] = false;}
     
      per = new double[n_rep*K1];
      mod = new double[n_rep*K1];
      arg = new double[n_rep*K1];

      plot_per = false; plot_mod = false; plot_arg = false;

    } 

    public void setDimensions(int _K, int _nr) 
    {
       K = _K; n_rep = _nr; K1 = K+1;
       per = new double[n_rep*K1];
       mod = new double[n_rep*K1];
       arg = new double[n_rep*K1]; 

    }


    public void setPeriodogram(double[] _per, int nr, int _K)
    {
      if(nr == n_rep && _K == K)
      {
        per = new double[_per.length];
        System.arraycopy(_per, 0, per, 0, _per.length);
        go();
      }       
    }

    public void setMod(double[] _mod, int nr, int _K)
    {
      if(nr == n_rep && _K == K)
      {
        mod = new double[_mod.length];
        System.arraycopy(_mod, 0, mod, 0, _mod.length);
        go();
      }       
    }

    public void setArg(double[] _arg, int nr, int _K)
    {
      if(nr == n_rep && _K == K)
      {
        arg = new double[_arg.length];
        System.arraycopy(_arg, 0, arg, 0, _arg.length);
        go();
      }   
    }

    public void setPlotPer(boolean t) {plot_per = t; go();} 
    public void setPlotMod(boolean t) {plot_mod = t; go();}
    public void setPlotArg(boolean t) {plot_arg = t; go();}    
    public void setPlot(boolean t, int c) {period_plots[c] = t; go();}


    public void dataSetup()
    {
       int i,k;
       dataMax = -100000.0; dataMin = 100000.0;  

       for(k=0;k<n_rep;k++)
       {
         if(period_plots[k])
         {
          for(i=0;i<K1;i++)
          {
            if(plot_per && per[K1*k+i] > dataMax) {dataMax = per[K1*k+i];}
            if(plot_per && per[K1*k+i] < dataMin) {dataMin = per[K1*k+i];}

            if(plot_mod && mod[K1*k+i] > dataMax) {dataMax = mod[K1*k+i];}
            if(plot_mod && mod[K1*k+i] < dataMin) {dataMin = mod[K1*k+i];}

            if(plot_arg && arg[K1*k+i] > dataMax) {dataMax = arg[K1*k+i];}
            if(plot_arg && arg[K1*k+i] < dataMin) {dataMin = arg[K1*k+i];}
          }
         }
       }
       dataNorm = Math.abs(dataMax - dataMin)*1.20;  
    }



    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,j,N,k;
     int t0, t1, x0, x1;

     
     super.paintComponent(g);
     g2d = (Graphics2D)g; 

     dataSetup();
    
     String intt,intd,loc;
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;


     BasicStroke orig = new BasicStroke((float)1.6);
     
        g2d.setStroke(dashed);
        g2d.setPaint(myGray2);
        
        //Draw dashed lines 

        for(i=0; i < 9; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }
 
        intd = Integer.toString(6);
        //int nobsP = (int)Math.floor((double)300/6); 

        g2d.setPaint(myGray);      
        loc = "\u03C0/6";
        t0 = (int)(((double)1/6)*(double)width);
	    g2d.drawLine(t0, 0, t0, height-10);
        g2d.setPaint(Color.PINK);
        g.drawString(loc, t0, height-5);  

        g2d.setPaint(Color.PINK);
        g.drawString("0", 5, height-5);
        g.drawString("\u03C0", width-10, height-5);     
        for(i=2; i < 6; i++) 
        { 
          intt = Integer.toString(i); 
          loc = intt + "\u03C0/" + intd;  

          t0 = (int)(((double)i/6)*(double)width);
          g2d.setPaint(myGray);
	  g2d.drawLine(t0, 0, t0, height-10);
          g2d.setPaint(Color.PINK);
          g.drawString(loc, t0, height-5);
        } 

        
       g2d.setStroke(orig);
        
       N = K; 
       for(k=0;k<n_rep;k++)
       {
         if(period_plots[k])
         {

          for(j=0;j<K1-1;j++)
          {

            if(plot_per)
            {
              g2d.setPaint(new Color(30*k,30,150));
	      t0 = (int)(((double)j/(double)N)*(double)width);
	      t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	      x0 = (int)(((per[K1*k+j] - dataMin)/dataNorm)*(double)height);
	      x1 = (int)(((per[K1*k+j+1] - dataMin)/dataNorm)*(double)height);
	      g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            }

            if(plot_mod)
            {
              g2d.setPaint(new Color(30*k,90,150));
	      t0 = (int)(((double)j/(double)N)*(double)width);
	      t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	      x0 = (int)(((mod[K1*k+j] - dataMin)/dataNorm)*(double)height);
	      x1 = (int)(((mod[K1*k+j+1] - dataMin)/dataNorm)*(double)height);
	      g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            }

            if(plot_arg)
            {
              g2d.setPaint(new Color(30*k,90,210));
	      t0 = (int)(((double)j/(double)N)*(double)width);
	      t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	      x0 = (int)(((arg[K1*k+j] - dataMin)/dataNorm)*(double)height);
	      x1 = (int)(((arg[K1*k+j+1] - dataMin)/dataNorm)*(double)height);
	      g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            }
          }
         }
       }

   }
}
