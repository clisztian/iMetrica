package ch.imetrica.usimx13;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;



public class SARIMAspectrumCanvas extends JPanel
{
    //-------SARMA Model stuff------

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int npts; // npts/2 per graph
    double[] sampleQIf;
    double[] sampleIn, sampleFw, sampleDn;

    //-------Graphics stuff
    Graphics2D g2d;
    int height, width;
    double  T0, T1, X0, X1, Y0, Y1;
    double dataMax, dataMin, dataNorm, dataMax2, dataMin2, dataNorm2;
    boolean plot_diff, dataSetup;
    Color myGray;
    boolean pseudoPlot; DecimalFormat df;
   
    boolean plotIn,plotFw,plotDn; 
    boolean plotGF,plotGI,plotG;      


    public SARIMAspectrumCanvas(int w, int h, int _nObs)
    {
       this.height = h; this.width = w; 
       this.npts = _nObs; 
       sampleQIf = new double[900];
       sampleIn = new double[_nObs]; sampleFw = new double[_nObs]; 
       sampleDn = new double[_nObs];
       setBackground(Color.BLACK);
       setPreferredSize(new Dimension(w, h));
       dataSetup = false;
       plot_diff = true;
       myGray = new Color(67,71,73);
       pseudoPlot = false;

       plotIn = false; plotFw = false; plotDn = false; 
       plotGF = true; plotGI = true; plotG = false; 

       df = new DecimalFormat("##.##"); 
    }

    public void setPlotDiff(boolean p)
    {plot_diff = p;}
  
    public void setSpectrumData(double[] sdata, int n)
    {
      int i; pseudoPlot = false;
      if(n != npts) {System.out.println("Error: Problem with data points");}
      for(i=0; i < n; i++) sampleQIf[i] = sdata[i];
         
      dataMax = -1000000;
      dataMin = 10000000;
  
      for(i=0; i < n; i++)
      {
        if(sdata[i] > dataMax) dataMax = sdata[i];
        else if(sdata[i] < dataMin) dataMin = sdata[i];   
      }
      dataNorm = Math.abs(dataMax - dataMin);
    }

    public void setSpectrumData2(double[] sdata1, double[] sdata2)
    {
      int i;
      dataMax2 = -1000000;
      dataMin2 = 10000000;
  
      for(i=0; i < 300; i++)
      {
        if(sdata1[i] > dataMax2) dataMax2 = sdata1[i];
        else if(sdata1[i] < dataMin2) dataMin2 = sdata1[i];   
      }
      for(i=0; i < 300; i++)
      {
        if(sdata2[i] > dataMax2) dataMax2 = sdata2[i];
        else if(sdata2[i] < dataMin2) dataMin2 = sdata2[i];   
      }

      dataNorm2 = Math.abs(dataMax2 - dataMin2);
    }


    public void setSpectralData(double[] sdata, double[] sdata2, double[] sdata3)
    {
      int i; 
      
      pseudoPlot = true;
      for(i=0; i < 300; i++) {sampleQIf[i] = sdata[i]; sampleQIf[i+300] = sdata2[i];
                              sampleQIf[i+600] = sdata3[i];}
         
      dataMax = -1000000;
      dataMin = 10000000;
  
      for(i=0; i < 900; i++)
      {
        if(sampleQIf[i] > dataMax) dataMax = sampleQIf[i];
        else if(sampleQIf[i] < dataMin) dataMin = sampleQIf[i];   
      }
      dataNorm = Math.abs(dataMax - dataMin);
    }

    public void setSpectralPlots(boolean _plotIn, boolean _plotFw, boolean _plotDn, 
                                 boolean _plotGF, boolean _plotGI, boolean _plotG)
    {
       plotIn = _plotIn;  plotFw = _plotFw; plotDn = _plotDn; 
       plotGF = _plotGF;  plotGI = _plotGI; plotG = _plotG;      
    }

    public void setSpectralFuncs(double[] In, double[] Fw, double[] Dn)
    {
      sampleIn = In; sampleFw = Fw; sampleDn = Dn;
      setSpectrumData2(sampleIn, sampleFw);
    }

    public void go()
    {repaint();}

    public void paintComponent(Graphics g)
    {
        int i;
        int t0,t1,x0,x1;
        float[] dash1 = {10.0f};
        BasicStroke dashed = new BasicStroke(1.0f, 
                                          BasicStroke.CAP_BUTT, 
                                          BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);
    
        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;

	super.paintComponent(g);
	g2d = (Graphics2D)g;	
        
        if(plotGF)
        {  
          
          g2d.setPaint(Color.GREEN);
          for(i = 0; i < 299; i++)
	  {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleQIf[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((sampleQIf[i+1] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //System.out.println(sampleQIf[i]);
	  }
        }
        if(plotGI)
        {
          //System.out.println("PRINTING GI....");
          //System.out.println("");
          g2d.setPaint(Color.MAGENTA);
          for(i = 0; i < 299; i++)
	  {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleQIf[300+i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((sampleQIf[301+i] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //System.out.println(sampleQIf[300+i]);
	  }
        }
        if(plotIn)
        {
          g2d.setPaint(Color.RED);
          for(i = 0; i < 299; i++)
	  {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleIn[i] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    x1 = (int)(((sampleIn[i+1] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	  }
        }
        if(plotFw)
        {
          g2d.setPaint(Color.BLUE);
          for(i = 0; i < 299; i++)
	  {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleFw[i] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    x1 = (int)(((sampleFw[i+1] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            
	  }
        }
        if(plotDn)
        {


          g2d.setPaint(Color.LIGHT_GRAY);
          for(i = 0; i < 299; i++)
	  {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleDn[i] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    x1 = (int)(((sampleDn[i+1] - dataMin2)/(dataNorm2*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
            //System.out.println(sampleDn[i]);
	  }
        }
        if(plotGF && plotGI)
        {        
         g2d.setPaint(Color.GRAY);
         for(i = 0; i < 300; i++)
	 {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)i/(double)300)*(double)width);
	    x0 = (int)(((sampleQIf[i] - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((sampleQIf[i+300] - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	 }
        }
        if(plotG)
        {      
         //System.out.println("Gets to print G");  
         //System.out.println("");
         g2d.setPaint(Color.YELLOW);
         for(i = 0; i < 299; i++)
	 {
	    t0 = (int)(((double)i/(double)300)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)300)*(double)width);
	    x0 = (int)(((sampleQIf[i+600]*10.0 - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((sampleQIf[i+601]*10.0 - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	    //System.out.println(sampleQIf[i+600]);
         }
        }

        g2d.setStroke(dashed);
        g2d.setPaint(myGray);
        
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
        g.drawString((String)df.format(dataMin), 5, height - 5);

     }	
}
