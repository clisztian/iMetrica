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
public class PeriodCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int K; //---total number of observations
    int Km;
    int nGamma;
    int n_rep; 

    double cutoff0, cutoff, cutoff2, cutoff3; 

    double [] gamma;
    double[] period_hat;  //periodograms of supporting data Y_i(t) Km observations
    double[] period_xf; //peridodograms for Y(t) and X(t) K observations
    double[] period_hybrid;
    double[] period_zpc;
    double[] period_orig;

    boolean[] period_plots;   //---- periodograms plot   
    boolean mdfa; 

    boolean g_plot, line_plot; 
    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;
    
    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;
   

    public PeriodCanvas(int w, int h, int _K, int _nrep)
    {
      // Initilize everything-------------------
      int i;
      K = _K; Km = _K; nGamma = _K;
      n_rep = _nrep;
      mdfa = true;
      
      cutoff = Math.PI/6; cutoff0 = 0.0; 
      g_plot = false; line_plot = true;

      cutoff2 = 0; cutoff3 = 0; 
      
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
      period_plots = new boolean[12]; 
      for(i=0;i<12;i++) {period_plots[i] = false;}
      period_plots[0] = true;
      
      period_hat = new double[(n_rep-1)*(K+1)];
      period_xf = new double[2*(K+1)]; 

    } 

    public void setPlotDim(double[] series, int N)
    {
        int i,j;
        double locMax = -100;
        dataMax = -100000.0; dataMin = 100000.0; 
        for(i=0;i<N;i++)
        {
         if(series[i] < dataMin) dataMin = series[i];
         else if(series[i] > dataMax) dataMax = series[i]; 
        } 
        dataNorm = Math.abs(dataMax - dataMin)*1.20;  
         
         
       //--- now normalize data for periodograms 
       
       for(j=0;j<period_hat.length;j++)
       {if(period_hat[j] > locMax) {locMax = period_hat[j];} }
       
       for(j=0;j<period_hat.length;j++)
       {period_hat[j]= period_hat[j]*dataNorm/locMax;}       
       
                     
    }

    public void setCutoffs(double c, double c1) 
    {cutoff0 = c; cutoff = c1; cutoff2 = 0; cutoff3 = 0; go();}

    public void setCutoffs2(double c, double c1, double c2, double c3) 
    {cutoff0 = c; cutoff = c1; cutoff2 = c2; cutoff3 = c3; go();}    
    
    public void setGammaPlotType(boolean graph, boolean line)
    {g_plot = graph; line_plot = line;}

    public void setMdfa(boolean m)
    {mdfa = m;}
  
    public void setPlots(int i, boolean sel)
    { period_plots[i] = sel; go(); }   

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
    public void setPeriodograms(double[] xf, double[] hat, int N, int rep)
    {   
       int le = hat.length; int le2 = xf.length;
       n_rep = rep; 
       if(n_rep > 1){setMdfa(true); this.Km = N; }
       else {setMdfa(false); this.K = N; }

       this.period_hybrid = new double[le2/2];
       this.period_zpc = new double[le2/2];
       this.period_hat = new double[le];
       this.period_xf = new double[le2]; 
       this.K = (int)(le2-2)/2;    

       System.arraycopy(hat, 0, this.period_hat, 0, le);  
       System.arraycopy(xf, 0, this.period_xf, 0, le2); 

       this.setPlotDim(xf,le2);
       this.go();
    }

    public void setPeriodogramXf(double[] xf, int N)
    {
       int le2 = xf.length; this.K=N;
       this.period_xf = new double[le2];
       this.period_orig = new double[le2];
       System.arraycopy(xf, 0, this.period_xf, 0, le2);
       System.arraycopy(xf, 0, this.period_orig, 0, le2);
       this.setPlotDim(xf,le2);
       this.go(); 
    }


    //----- set raw data and filtered data
    public void setPeriodogram_zpc(double[] ts)
    {   
       int le = ts.length/2;      
       this.period_zpc = new double[le];
       System.arraycopy(ts, 0, this.period_zpc, 0, le);  
    }    
 
    //----- set raw data and filtered data
    public void setPeriodogram_hybrid(double[] ts)
    {   
       int le = ts.length/2; 
       this.period_hybrid = new double[le];
       System.arraycopy(ts, le, this.period_hybrid, 0, le);  
    } 

    public void plotHybrid(int c)
    {
       int le = this.period_xf.length; int le2 = le/2;
       if(c==2) 
       {
          System.arraycopy(this.period_orig, 0, this.period_xf, 0,le2);
          //System.arraycopy(this.period_hybrid, 0, this.period_xf, le/2, le-1);
          for(int k = 0;k<this.period_hybrid.length;k++) {this.period_xf[le2+k] = this.period_hybrid[k];}
  
       }
       else if(c==0) {System.arraycopy(this.period_orig, 0, this.period_xf, 0, le-1);}
       else if(c==1) 
       {
          System.arraycopy(this.period_orig, 0, this.period_xf, 0, le2);
          //System.arraycopy(this.period_zpc, 0, this.period_xf, le2, le-1); 
          for(int k = 0;k<this.period_zpc.length;k++) {this.period_xf[le2+k] = this.period_zpc[k];}    
       }
       this.setPlotDim(this.period_xf,le); go();
    }    



    public void setGamma(double[] _xf)
    {
       int le = _xf.length;
       nGamma = le-1;
       this.gamma = new double[le];
       System.arraycopy(_xf, 0, this.gamma, 0, le); 
       //this.setPlotDim(_xf,le);
       this.go();
       //--- if filter fine, no redim needed 
    }



    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int i,j,N,k;
     int t0, t1, x0, x1, x2, x3;
     int cutoffInt0, cutoffInt;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
    
     double[] agg_period; int count = 0;
     GradientPaint whitetopurp;
     String intt,intd,loc;
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
    
     cutoffInt = (int)(cutoff*K/Math.PI);
     cutoffInt0 = (int)(cutoff0*K/Math.PI);

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


     g2d.setPaint(new Color(118,0,52));
     if(g_plot)
     {
        N = nGamma;
        for(j = 0; j < N-1; j++)
        {
          
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)((gamma[j]/1.5)*(double)height);
	  x1 = (int)((gamma[j+1]/1.5)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }
     }

     if(line_plot)
     { 
        t0 = (int)((cutoff/Math.PI)*(double)width);
        g2d.drawLine(t0-1, 0, t0-1, height-10);
        if(cutoff0 > 0.0)
        {
          t0 = (int)((cutoff0/Math.PI)*(double)width);  
          g2d.drawLine(t0-1, 0, t0-1, height-10); 
        }
        g2d.setPaint(new Color(100,0,92));
        if(cutoff2 > 0.0)
        {
          t0 = (int)((cutoff2/Math.PI)*(double)width);  
          g2d.drawLine(t0-1, 0, t0-1, height-10); 
        }
        if(cutoff3 > 0.0)
        {
          t0 = (int)((cutoff3/Math.PI)*(double)width);  
          g2d.drawLine(t0-1, 0, t0-1, height-10); 
        }   
     }
     
     g2d.setPaint(myGray);
     g.drawString((String)df.format(dataMax), 5, 15);
     orig = new BasicStroke((float)1.3); g2d.setStroke(orig);
     x0 = (int)(((0 - dataMin)/dataNorm)*(double)height);	
     g2d.drawLine(0, (height-20) - x0, width, (height-20) - x0);


     N = K+1; //System.out.println("length of periodogram = " + period_xf.length + "N = " + N);
     if(period_plots[0] && period_plots[1])
     {        //new Color(230,230,250)
        whitetopurp = new GradientPaint(0, height, new Color(72,61,139), 0, 0, new Color(149,137,219));  
        whitetopurp = new GradientPaint(0, height, new Color(0,132,201), 0, 0, new Color(255,132,201));      
        //orig = new BasicStroke((float)1.3);
       
        for(j = cutoffInt0; j < cutoffInt; j++)
	{   
            g2d.setPaint(whitetopurp);
            //g2d.setStroke(orig);
            //g2d.setPaint(Color.GRAY);
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j] - dataMin)/(dataNorm))*(double)height);
            x1 = (int)(((period_xf[j+1] - dataMin)/(dataNorm))*(double)height);
	    x2 = (int)(((period_xf[j+N+1] - dataMin)/(dataNorm))*(double)height);
            x3 = (int)(((period_xf[j+N] - dataMin)/(dataNorm))*(double)height);
            int[] ts1 = {(height-20) - x0, (height-20) - x1, (height-20) -  x2, (height-20) - x3}; 
            int[] xs1 = {t0,t1,t1,t0};
	    //g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //poly = new Polygon({t0,t0,t1,t1}, {x0,x1,x0,x1}, 4) 
            g.fillPolygon(xs1, ts1, 4);

            g2d.setPaint(new Color(181,99,138)); //g2d.setPaint(new Color(132,83,155));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

            g2d.setPaint(new Color(98,130,201));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j+N] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1+N] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

       }
       for(j = 0; j < cutoffInt0; j++)
       {   
            g2d.setPaint(whitetopurp);
            //g2d.setStroke(orig);
            //g2d.setPaint(Color.GRAY);
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((0.0 - dataMin)/(dataNorm))*(double)height);
            x1 = (int)(((0.0 - dataMin)/(dataNorm))*(double)height);
	    x2 = (int)(((period_xf[j+N+1] - dataMin)/(dataNorm))*(double)height);
            x3 = (int)(((period_xf[j+N] - dataMin)/(dataNorm))*(double)height);
            int[] ts1 = {(height-20) - x0, (height-20) - x1, (height-20) -  x2, (height-20) - x3}; 
            int[] xs1 = {t0,t1,t1,t0};
	    //g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //poly = new Polygon({t0,t0,t1,t1}, {x0,x1,x0,x1}, 4) 
            g.fillPolygon(xs1, ts1, 4);

            g2d.setPaint(new Color(181,99,138)); //g2d.setPaint(new Color(132,83,155));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

            g2d.setPaint(new Color(98,130,201)); //g2d.setPaint(new Color(0,130,255));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j+N] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1+N] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

       }
       for(j = cutoffInt; j < N-1; j++)
       {   
            g2d.setPaint(whitetopurp);
            //g2d.setStroke(orig);
            //g2d.setPaint(Color.GRAY);
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((0.0 - dataMin)/(dataNorm))*(double)height);
            x1 = (int)(((0.0 - dataMin)/(dataNorm))*(double)height);
	    x2 = (int)(((period_xf[j+N+1] - dataMin)/(dataNorm))*(double)height);
            x3 = (int)(((period_xf[j+N] - dataMin)/(dataNorm))*(double)height);
            int[] ts1 = {(height-20) - x0, (height-20) - x1, (height-20) -  x2, (height-20) - x3}; 
            int[] xs1 = {t0,t1,t1,t0};
	    //g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //poly = new Polygon({t0,t0,t1,t1}, {x0,x1,x0,x1}, 4) 
            g.fillPolygon(xs1, ts1, 4);

            g2d.setPaint(new Color(181,99,138)); //g2d.setPaint(new Color(132,83,155));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

            g2d.setPaint(new Color(98,130,201)); //g2d.setPaint(new Color(0,130,255));
	    t0 = (int)(((double)j/(double)N)*(double)width);
	    t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	    x0 = (int)(((period_xf[j+N] - dataMin)/dataNorm)*(double)height);
	    x1 = (int)(((period_xf[j+1+N] - dataMin)/dataNorm)*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);

       }
    }
    else
    {  
      g2d.setStroke(orig);  
      if(period_plots[0])
      {       
        g2d.setPaint(new Color(181,99,138));
        for(j = 0; j < N-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((period_xf[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((period_xf[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }
        
      }     
      if(period_plots[1])
      { 
           
        g2d.setPaint(new Color(164,97,255));
        for(j = 0; j < N-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((period_xf[j+N] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((period_xf[j+N+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }     
      } 
     }
  

     //----- plot aggregate peridodogram
     N = Km+1;
     orig = new BasicStroke((float)1.5); 
     g2d.setStroke(orig);    
     agg_period = new double[N]; count = 0;
     if(mdfa)
     {     
      for(k=0;k<n_rep-1;k++)
      {
       if(period_plots[k+2])
       {for(j = 0; j < N; j++) {agg_period[j] = agg_period[j] + period_hat[N*k + j];} count++;}
      }
      if(count > 0)
      {for(j = 0; j < N; j++) {agg_period[j] = agg_period[j]/count;}}
     }
      
     if(count > 0)
     {
      g2d.setPaint(new Color(255-2*20,55+2*20,255));
      for(j = 0; j < N-1; j++)
      {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((agg_period[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((agg_period[j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
      }
     }
     
   }


}
