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
public class FilterCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int K; //---total number of observations
    int n_rep;
    double[] gamma_hat;  //gamma_hat n_rep x (K+1)
    double[] gamma; //gamma K+1
    double[] gamma_zpc;  //gamma_hat n_rep x (K+1)
    double[] gamma_hybrid;  //gamma_hat n_rep x (K+1)
    double[] gamma_orig;

    boolean[] gamma_plots;   //---- gamma functions plot  
    boolean g_plot;      //---- filtered data 

    boolean mdfa;
    

    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;

    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;
   

    public FilterCanvas(int w, int h, int _K, int _nrep)
    {
      // Initilize everything-------------------
      int i;
      K = _K; 
      n_rep = _nrep;
      mdfa = true;

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
      gamma_plots = new boolean[12]; 
      for(i=0;i<12;i++) {gamma_plots[i] = false;}
      g_plot = false; 

      gamma_hat = new double[(n_rep+1)*(K+1)];
      gamma = new double[K+1]; 
    } 

    public void setPlotDim()
    {
        int i,k; int K1 = K+1;
        dataMax = -100000.0; dataMin = 100000.0; 

        //System.out.println("n_reps = " + n_rep);

        if(n_rep > 1)
        {
         for(k=0;k<=n_rep;k++)
         {
          if(gamma_plots[k+1])
          {
            for(i=0;i<K1;i++)
            {
             if(gamma_hat[k*K1+i] < dataMin) dataMin = gamma_hat[k*K1+i];
             else if(gamma_hat[k*K1+i] > dataMax) dataMax = gamma_hat[k*K1+i]; 
            } 
          } 
         }
        }
        else
        {
          if(gamma_hat.length == K1)
          {
            for(i=0;i<K1;i++)
            {
             if(gamma_hat[i] < dataMin) dataMin = gamma_hat[i];
             else if(gamma_hat[i] > dataMax) dataMax = gamma_hat[i]; 
            } 
          }     
        }
        if(g_plot)
        {
            for(i=0;i<K1;i++)
            {
             if(gamma[i] < dataMin) dataMin = gamma[i];
             else if(gamma[i] > dataMax) dataMax = gamma[i]; 
            } 
        } 
        dataNorm = Math.abs(dataMax - dataMin)*1.20;                
    }


    public void setMdfa(boolean m)
    {mdfa = m;}
  
    public void setPlots(int i, boolean sel)
    {
       if(i==0) 
       {g_plot = sel;} 
       else if(i <= n_rep+1)
       {gamma_plots[i] = sel;}
       setPlotDim(); go();
    }   

    public void setNRep(int n) 
    {
      n_rep = n; 
      if(n_rep > 1) 
      {setMdfa(true);}
      else
      {setMdfa(false);}
    }
    public void setK(int n) {K = n; }

    //----- set raw data and filtered data
    public void setGamma_hat(double[] ts, int _K, int nrep)
    {   
       int le = ts.length; this.K = _K; this.n_rep = nrep;  setNRep(nrep);    
       this.gamma_hat = new double[le];
       this.gamma_orig = new double[le];
       System.arraycopy(ts, 0, this.gamma_hat, 0, le);  
       System.arraycopy(ts, 0, this.gamma_orig, 0, le);  
       setPlotDim(); go();
       this.gamma_hybrid = new double[le];
       this.gamma_zpc = new double[le];
    }

    //----- set raw data and filtered data
    public void setGamma_zpc(double[] ts)
    {   
       int le = ts.length;
       this.gamma_zpc = new double[le];
       System.arraycopy(ts, 0, this.gamma_zpc, 0, le);  
    }    
 
    //----- set raw data and filtered data
    public void setGamma_hybrid(double[] ts)
    {   
       int le = ts.length;
       this.gamma_hybrid = new double[le];
       System.arraycopy(ts, 0, this.gamma_hybrid, 0, le);  
    } 

    public void plotHybrid(int c)
    {
       int le = this.gamma_orig.length; 
       if(c==2) {System.arraycopy(this.gamma_hybrid, 0, this.gamma_hat, 0, le);}
       else if(c==0) {System.arraycopy(this.gamma_orig, 0, this.gamma_hat, 0, le);}
       else if(c==1) {System.arraycopy(this.gamma_zpc, 0, this.gamma_hat, 0, le);}
       setPlotDim(); go();
    }    

    public void setGamma(double[] _xf)
    {
       //int i;
       int le = _xf.length;
       K = le-1;
       this.gamma = new double[le];
       System.arraycopy(_xf, 0, this.gamma, 0, le); 
       setPlotDim(); go();

       //for(int i=0;i<le;i++)
       //{System.out.println((i+1) + " "+this.gamma[i]);}
    }


    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,j,N,k;
     int t0, t1, x0, x1;
     N = K+1;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 

     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;

        g2d.setStroke(dashed);
        g2d.setPaint(myGray2);
        
        //Draw dashed lines 

        for(i=0; i < 9; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }

        //int nobsP = (int)Math.floor((double)300/6); 
        for(i=1; i <= 6; i++) 
        {
          t0 = (int)(((double)i/6)*(double)width);
	  g2d.drawLine(t0, 0, t0, height-10);
        } 

        g2d.setPaint(myGray);
        x0 = (int)(((0 - dataMin)/dataNorm)*(double)height);	
        g2d.drawLine(0, (height-20) - x0, width, (height-20) - x0);
        //g.drawString((String)df.format(0), 5, x0+5);

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


     orig = new BasicStroke((float)1.3);
     g2d.setStroke(orig);
     g2d.setPaint(plotGray);
     if(mdfa)
     {     
      for(k=0;k<=n_rep;k++)
      {
       g2d.setPaint(new Color(15+(k+1)*20,102,(k+1)*8+156));
       if(gamma_plots[k+1])
       {
        for(j = 0; j < N-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((gamma_hat[N*k + j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((gamma_hat[N*k + j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }      
       }
      }
     }
     else
     {       
       if(gamma_plots[1])
       {       
        g2d.setPaint(new Color(15+40,102,156));
        for(j = 0; j < N-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((gamma_hat[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((gamma_hat[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }
        
       }
     }
     if(g_plot)
     {
        g2d.setPaint(new Color(15,102,156));
        for(j = 0; j < N-1; j++)
        {
          
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((gamma[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((gamma[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }
     }

   }


}