package ch.imetrica.emd;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;

public class EMDamfm extends JPanel
{
  
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int nObs;    
    double[] tseries; 

    //------- canvas stuff -------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;

    //----------- EMD stuff ------------------
    private int n_imfs;
    private double[][] am; private double[][] fm;
    private double[][] inst_f; private double[][] phase;
    double[] res_trend;
    double[] agg;
    //-------------------------------------------

    //---------- Plotting options --------------     
    boolean plot_ts;
    boolean plot_res;
    boolean[] aggregate;   

    boolean[][] plots_amfm;
    //------------------------------------------

    DecimalFormat df;
    BasicStroke dashed;
    float[] dash1;
    Color myGray;  
    Color myred,mygreen;
    //boolean IMFs;

   public EMDamfm(int w, int h, int _nObs)
   {
    //IMFs = _IMFs;
    System.loadLibrary("emd_util");
    this.width = w; this.height = h; nObs = _nObs;
    plot_ts = true; plot_res = false; aggregate = new boolean[6];  
    this.initPlotOptions(6);  
    dataMax = -1000000.0; dataMin = 1000000.0;

    setBackground(Color.BLACK);
    setPreferredSize(new Dimension(w, h));

     dash1 = new float[1]; dash1[0] = 10.0f;
     dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f); 
     myGray = new Color(67,71,73);

     myred = new Color(142,238,255);
     mygreen = new Color(228,87,137);

     df = new DecimalFormat("##.##"); 
     aggregate = new boolean[6];
   }
 
   public void setNObs(int n) {nObs = n;}

   public void initPlotOptions(int nimfs)
   {
     int i;
     plots_amfm = new boolean[4][nimfs];
     for(i=0;i<nimfs;i++)
     {plots_amfm[0][i] = false; plots_amfm[1][i] = false; plots_amfm[2][i] = false; plots_amfm[3][i] = false;} 
   }

   public void plotResidual(boolean p) {plot_res = p;}
   public void plotSeries(boolean p) {plot_ts = p;}

   public void AM_FM_decomp(double[] series, int _N)
   {
    int len,m,i,N; 
    
    nObs = _N; N=_N;
    tseries = new double[_N];
    tseries = series;

    double[] output = getAM_FM(series, _N);
    len = output.length;  
    setN_imfs((int)output[len-1]);

    setAm(new double[getN_imfs()][_N]);
    setFm(new double[getN_imfs()][_N]);
    setInst_f(new double[getN_imfs()][_N]);
    setPhase(new double[getN_imfs()][_N]);
    res_trend = new double[_N];
    agg = new double[nObs];
       
    for(m=0;m<getN_imfs();m++)
    {
     for(i=0;i<N;i++)
     {
        getAm()[m][i] = output[4*_N*m+i]; 
        getFm()[m][i] = output[4*_N*m+_N+i];     
        getInst_f()[m][i] = output[4*_N*m+2*_N+i];
        getPhase()[m][i] = output[4*_N*m+3*_N+i];
     }
    }
    dataMax = -1000000.0; dataMin = 1000000.0;
    for(i=0;i<N;i++)
    {
      res_trend[i] = output[4*_N*getN_imfs()+i]; 
      agg[i] = res_trend[i];
      
      if(series[i] < dataMin) dataMin = series[i];
      else if(series[i] > dataMax) dataMax = series[i]; 
    } 
    compAgg();  
    dataNorm = Math.abs(dataMax - dataMin)*1.20;    
    go();
  }

  public void setPlotOptions(boolean[][] _plot) 
  {plots_amfm = _plot;}

  public void setAggregateBool(int i, boolean s)
  {if(i < getN_imfs()) aggregate[i] = s;}

  public void addAggregate(int w)
  {
    int j;
    if(w < getN_imfs())
    {
     for(j=0;j<nObs;j++)
     {agg[j] = agg[j] + getAm()[w][j]*getFm()[w][j];}    
    }
    go();
  }

  public void subAggregate(int w)
  {
    int j;
    if(w < getN_imfs())
    {
     for(j=0;j<nObs;j++)
     {agg[j] = agg[j] - getAm()[w][j]*getFm()[w][j];}    
    }
    go();
  }

  public void compAgg()
  {
     int m,i;
     for(m=0;m<getN_imfs();m++)
     {
       if(aggregate[m])
       {
        for(i=0;i<nObs;i++)
        {agg[i] = agg[i] + getAm()[m][i]*getFm()[m][i];}
       }
     }
  }
 
 
  public native double[] getAM_FM(double[] data, int N);
  static {System.loadLibrary("emd_util");}
 
  public void go() {repaint();}

  //--------- Plots the aggregate series ---------
  public void paintComponent(Graphics g)
  {  
     int i,j,N;
     int t0, t1, x0, x1;
     N = nObs;  

     super.paintComponent(g);
     g2d = (Graphics2D)g; 

     BasicStroke def = new BasicStroke((float)1.3);

     g2d.setStroke(def);
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
     //System.out.println(width + "  " + height);

     g2d.setPaint(myred);     
     if(plot_ts)
     {
       for(j = 0; j < N-1; j++)
       {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((tseries[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((tseries[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
       }
     } 
     if(plot_res)
     { 
       g2d.setPaint(Color.BLUE);  
       for(j = 0; j < N-1; j++)
       {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((res_trend[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((res_trend[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
       }
     }      

     // ------- Plot aggregate series ---------
     //g2d.setPaint(new Color(24,152,105));
     g2d.setPaint(mygreen);
     for(j = 0; j < N-1; j++)
     {
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((agg[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((agg[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
          
     }

    
        //Draw dashed lines 
        g2d.setStroke(dashed);
        g2d.setPaint(myGray);

        for(i=0; i < 9; i++)
        {
            x0 = (int)(((double)i/(double)8)*(double)height);
            g2d.drawLine(0, x0, width, x0);
        }
        g.drawString((String)df.format(dataMax), 5, 15);
        g.drawString((String)df.format(dataMin), 5, height - 5);
        int nobsP = (int)Math.floor((double)N/12); int p = 12;
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*12)/N)*(double)width);
	      g2d.drawLine(t0, 0, t0, height-20);
          g.drawString((String)"" + p, t0, height - 5); p = p + 12;
        } 
   
  }

public double[][] getFm() {
	return fm;
}

public void setFm(double[][] fm) {
	this.fm = fm;
}

public double[][] getInst_f() {
	return inst_f;
}

public void setInst_f(double[][] inst_f) {
	this.inst_f = inst_f;
}

public int getN_imfs() {
	return n_imfs;
}

public void setN_imfs(int n_imfs) {
	this.n_imfs = n_imfs;
}

public double[][] getAm() {
	return am;
}

public void setAm(double[][] am) {
	this.am = am;
}

public double[][] getPhase() {
	return phase;
}

public void setPhase(double[][] phase) {
	this.phase = phase;
}
}

