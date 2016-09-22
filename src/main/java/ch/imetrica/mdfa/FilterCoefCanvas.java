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
public class FilterCoefCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int L,M; //---total number of observations
    int n_rep;

    double[] b;     //----- coefficients 
    double[] sym_b;
    boolean[] b_plots;   //---- periodograms plot   
    boolean mdfa; 
    boolean update_coeff = false;
    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;
    
    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;
   

    public FilterCoefCanvas(int w, int h, int _L, int _nrep)
    {
      // Initilize everything-------------------
      int i;
      L = _L;  n_rep = _nrep; mdfa = true;
      update_coeff = false;
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
      b_plots = new boolean[11]; 
      for(i=0;i<11;i++) {b_plots[i] = false;}
      b_plots[1] = true;    
      b = new double[n_rep*L];

    } 

    public void setUpdateCoeff(boolean t) 
    {update_coeff = t;}
    
    
    public void setPlotDim()  //only considers dimensions of those selected 
    {
        int i,k;
        dataMax = -100000.0; dataMin = 100000.0; 

        for(k=0;k<n_rep;k++)
        {
          if(b_plots[k+1])
          {
            for(i=0;i<L;i++)
            {
             if(b[k*L+i] < dataMin) dataMin = b[k*L+i];
             else if(b[k*L+i] > dataMax) dataMax = b[k*L+i]; 
            } 
          } 
        }
        dataNorm = Math.abs(dataMax - dataMin)*1.20;           
    }


    public void setMdfa(boolean m)
    {mdfa = m;}
  
    public void setPlots(int i, boolean sel)
    { b_plots[i] = sel; setPlotDim(); go(); }   

    public void setNRep(int n) 
    {
      n_rep = n; 
      if(n_rep > 1) 
      {setMdfa(true);}
      else
      {setMdfa(false);}
    }
    public void setL(int n) {L = n;}

    public void setSymB(double[] _b, int _M)
    {
        
        int le = _b.length; if(_M != le) {System.out.println("Problems with length");} 
        sym_b = new double[le];
        M = _M; System.arraycopy(_b, 0, sym_b, 0, le);        
    }


    public void setBCoeffs(double[] _b, int _L, int _nreps)
    {
        L = _L; n_rep = _nreps; int le = _b.length; int i; 
        if(n_rep > 1){setMdfa(true);}
        else {setMdfa(false);}
  
        b = new double[n_rep*L];
        if(le != L*n_rep)
        {System.out.println("You got length problems " + le + "  " + L*n_rep);}
        System.arraycopy(_b, 0, this.b, 0, le);
        setPlotDim(); go(); 
        
        if(update_coeff) //---- subtract out the 1 from the first index
        {
          for(i=0;i<n_rep;i++)
          {
            b[L*i] = b[L*i] - 1.0; 
          }
        }
        
    }

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N,k;
     int t0, t1, x0, x1;
     
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
    
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
    
     Stroke orig = new BasicStroke((float)1.5); 
   
     g2d.setPaint(myGray);
     g.drawString((String)df.format(dataMax), 5, 15);
     g.drawString((String)df.format(dataMin), 5, height - 24);

     g2d.setStroke(orig);
     
     N = L;
     if(mdfa)
     {     
       for(k=0;k<n_rep;k++)
       {
        g2d.setPaint(new Color(15+(k+1)*20,102-(k+1)*5,186));
        if(b_plots[k+1])
        {
         for(j = 0; j < N-1; j++)
         {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)(N-1))*(double)width);
	  t1 = (int)(((double)(j+1)/(double)(N-1))*(double)width);
	  x0 = (int)(((b[N*k + j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((b[N*k + j + 1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
         }      
        }
       }
     }
     else
     {       
       if(b_plots[1])
       {       
        g2d.setPaint(new Color(15+40,102,156));
        for(j = 0; j < N-1; j++)
        {
          //System.out.println(gamma_hat[N*k+j]);
	  t0 = (int)(((double)j/(double)N)*(double)width);
	  t1 = (int)(((double)(j+1)/(double)N)*(double)width);
	  x0 = (int)(((b[j] - dataMin)/dataNorm)*(double)height);
	  x1 = (int)(((b[j+1] - dataMin)/dataNorm)*(double)height);
	  g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
        }
        
       }
     }

   
   }
}
