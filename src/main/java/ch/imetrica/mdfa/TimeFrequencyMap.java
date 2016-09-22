package ch.imetrica.mdfa;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;

public class TimeFrequencyMap extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//------ Plot stuff
    Graphics2D g2d;
    int height, width;
    int K1,line;
    BasicStroke dashed;
    DecimalFormat df;
    float[] dash1;
    Color myGray;
    Color[] colorArray,colorAM,colorFM,colorIF,colorPh;
    int colorCount;

    int xCanvasPanelCount;  // determines the resolution of the canvas
    int yCanvasPanelCount;

    int xDataCount;
    int yDataCount;
    int dataType;
    boolean plotok; 

    int colorType; 

    double[][] contourData;
    double dataRangeMin;
    double dataRangeMax;


    double[] freq_ints;
    int n_div,L;

    //---- EMD stuff -------
    int N,n_imfs;
    double[][] fm; //

    public TimeFrequencyMap()
    {
     
     dataRangeMax = -100000000.0; dataRangeMin = 10000000.0;
     //setPreferredSize(new Dimension(w, h));
     plotok = false;
     int co = 150;
     N=300;
     n_imfs = 11; colorType = 0;
     fm = new double[n_imfs][N];
     colorCount = co;
     colorArray = new Color[co];

     colorAM = new Color[co]; 
     colorFM = new Color[co]; 
     colorIF = new Color[co];
     colorPh = new Color[co];
     setupColors(); 
     dataType = 0;
     dash1 = new float[1]; dash1[0] = 10.0f;
     dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER,10.0f, dash1, 0.0f);    
     setBackground(Color.BLACK);
    }

    public void setNDiv(double[] f)
    {
      freq_ints = f; n_div = f.length;
    }
 
    public void setL(int l) {L = l;}

    public void setLine(int val, int _K1)
    {line = val; K1 = _K1; go(); }
 
    public void changeColor(int col) 
    {
       colorType = col;  
       setupColors(); 

      if(dataType == 0) System.arraycopy(colorAM, 0, colorArray, 0, colorAM.length);  
      else if(dataType == 1) System.arraycopy(colorFM, 0, colorArray, 0, colorFM.length);
      else if(dataType == 3) System.arraycopy(colorIF, 0, colorArray, 0, colorIF.length);
      else if(dataType == 2) System.arraycopy(colorPh, 0, colorArray, 0, colorPh.length);    
      go();
    }

    public void setResolution(int PanelCountx, int PanelCounty)
    {xCanvasPanelCount = PanelCountx; yCanvasPanelCount = PanelCounty;}

    public void updateFM(double[][] _fm, int _N, int _n_imfs, int sel)
    {
      int i,m; 
      N = _N; n_imfs = _n_imfs;      
      fm = _fm;
  
      dataRangeMax = -100000000.0; dataRangeMin = 10000000.0;
      for(m=0;m<n_imfs;m++)
      {
         for(i=0; i < N; i++)
         {
          if(fm[m][i] < dataRangeMin) dataRangeMin = fm[m][i];
          else if(fm[m][i] > dataRangeMax) dataRangeMax = fm[m][i];
         } 
      }

      Dimension ds = this.getSize();
      width = ds.width; height = ds.height; 
      computeColorInterp(height/2);
 
      if(sel == 0) System.arraycopy(colorAM, 0, colorArray, 0, colorAM.length);  
      else if(sel == 1) System.arraycopy(colorFM, 0, colorArray, 0, colorFM.length);
      else if(sel == 3) System.arraycopy(colorIF, 0, colorArray, 0, colorIF.length);
      else if(sel == 2) System.arraycopy(colorPh, 0, colorArray, 0, colorPh.length);    
      dataType = sel; plotok = true;
      go();
    }
    
    /*--------------------------------------------
      Compute the color map by interpolation between
      each FM mode
    ----------------------------------------------*/
  
    public void computeColorInterp(int hres)
    {
       int m,i,j,start; 
       int[] hress;
       int count;
       //yDataCount = (n_imfs-1)*hres+fade;
       yDataCount = hres;
       xDataCount = 2*N;

       contourData = new double[yDataCount][xDataCount];
       double p1,p2,pp1,pp2,t; 
  
       hress = new int[n_div-1];
  
       Dimension ds = this.getSize();
       width = ds.width; height = ds.height; 
 

       hress[0] = (int)(yDataCount*((freq_ints[0])/Math.PI));  //System.out.println(freq_ints[0] + " " + hress[0]);
       for(i=2;i<n_div;i++)
       {hress[i-1] = (int)(yDataCount*((freq_ints[i] - freq_ints[i-1])/Math.PI)); } //System.out.println(freq_ints[i] + " " + hress[i-1]);}
    
       for(i=0;i<N-1;i++) {contourData[0][2*i] = 0.0; contourData[0][2*i+1] = 0.0;}

 
        start=0; count=0;
       for(i=0;i<N-1;i++)
       {       
           count=0;
           p1 = 0; p2 = fm[0][i];  
           pp1 = 0; 
           pp2 = fm[0][i] + .5*(fm[0][i+1] - fm[0][i]);
          
           for(j=0;j<hress[0]+hress[1];j++)
           {
             t = (double)j/(hress[0]+hress[1]-1); 
             contourData[count][2*i] =  p1 + t*(p2 - p1); 
             contourData[count][2*i+1] =  pp1 + t*(pp2 - pp1);
             
             count++;             
           }
           start = count;
    
           
        }
        
  
       for(i=0;i<N-1;i++)
       {       
         count=start;
         for(m=1;m<n_div-2;m++)
         {
        
           p1 = fm[m-1][i]; p2 = fm[m][i];  
           pp1 = fm[m-1][i] + .5*(fm[m-1][i+1] - fm[m-1][i]); 
           pp2 = fm[m][i] + .5*(fm[m][i+1] - fm[m][i]);

           
           for(j=0;j<hress[m+1];j++)
           {
             t = (double)j/(hress[m+1]-1); 
             contourData[count][2*i] =  p1 + t*(p2 - p1); 
             contourData[count][2*i+1] =  pp1 + t*(pp2 - pp1);
             count++;
           }
          }
        }
       

       /*for(i=0;i<N-1;i++)
       {      
         //-------------------------- Vertical interpolation first ------------
         for(m=0;m<n_imfs-1;m++)
         { 
           p1 = fm[m][i]; p2 = fm[m+1][i];  
           pp1 = fm[m][i] + .5*(fm[m][i+1] - fm[m][i]); pp2 = fm[m+1][i] + .5*(fm[m+1][i+1] - fm[m+1][i]);  
           for(j=0;j<hres;j++)
           {t = (double)j/(hres-1); contourData[m*hres + j][2*i] =  p1 + t*(p2 - p1); contourData[m*hres + j][2*i+1] =  pp1 + t*(pp2 - pp1);}
         }
         p1 = fm[n_imfs-1][i]; p2 = base; pp1 = fm[m][i] + .5*(fm[m][i+1] - fm[m][i]); pp2 = base;  
         for(j=0;j<fade;j++)
         {t = (double)j/(fade-1); contourData[(n_imfs-1)*hres + j][2*i] =  p1 + t*(p2 - p1); contourData[(n_imfs-1)*hres + j][2*i+1] =  pp1 + t*(pp2 - pp1);} 
       }

         for(m=0;m<n_imfs-1;m++)
         {  
           p1 = fm[m][N-1]; p2 = fm[m+1][N-1];
           for(j=0;j<hres;j++)
           {t = (double)j/(hres-1); contourData[m*hres + j][xDataCount-1] = p1 + t*(p2 - p1);}
         }
         p1 = fm[n_imfs-1][N-1]; p2 = base;
         for(j=0;j<fade;j++)
         {t = (double)j/(fade-1); contourData[(n_imfs-1)*hres + j][xDataCount-1] =  p1 + t*(p2 - p1);}
       */
  

    }
          

    public void setupColors()
    {      
        int i,color;
   
        float hinc;
        float h;
        float saturation;
        float intensity;
      
        float hred = (float)2.0/colorCount;   
        float red = (float)0.0;                
        Color color2;
        

        if(colorType == 0)
        {     
         //-------------IF Colors--------------- 
         hinc = (float)(1.0/150); h=(float)0;
         saturation = (float)(0.7);
         intensity  = (float)(0.7);
         for(i = 0; i < 150; i++)
         {
           color = Color.HSBtoRGB(h, saturation, intensity); 

           colorIF[i] = new Color(color); colorAM[i] = new Color(color); 
           colorPh[i] = new Color(color); colorFM[i] = new Color(color);       
           h +=  hinc;
         
         } 
        }

       else
       {
       
        h = (float)225/(float)360;
        red = (float)0;
        hred = (float)(1.0/50.0);
        //----------- Sahar Colors ----------------------------
        for(i=1;i<=50;i++)
        {
           color = Color.HSBtoRGB(h,(float)1.0,(float)i/(float)50); //for 1/4
           colorIF[i-1] = new Color(color); colorAM[i-1] = new Color(color); 
           colorPh[i-1] = new Color(color); colorFM[i-1] = new Color(color);
        }
        for(i=0;i<50;i++)
        {          
          color2 = new Color(red,(float)65/255,(float)1.0); //for 1/2
          colorIF[50+i] = color2; colorAM[50+i] = color2; 
          colorPh[50+i] = color2; colorFM[50+i] = color2;
          red = red + hred; 
        } 
        red = (float)0;
        for(i=0;i<50;i++)
        {
           color2 = new Color((float)1.0, (float)65/255, (float)1.0 - red); //for 1/4
           colorIF[100+i] = color2; colorAM[100+i] = color2; 
           colorPh[100+i] = color2; colorFM[100+i] = color2;
           red = red + hred; 
        }
        
      }
      


    }   



   /*------------------------------------------------------------
     Using bilinear interpolation of the input data; Data is
     assumed to be at the nodes of a grid. The rectangle color value
     is taken to be the value at the center of the rectangle of
     a "canvas" grid decomposed into xCanvasPanelCount X yCanvasPanelCount 
     panels 
   ------------------------------------------------------------------------*/

  public void go() {repaint();}

  public void paintComponent(Graphics g)
  {
    int x0; super.paintComponent(g);
    if(contourData == null) return;
    if(Math.abs(dataRangeMax- dataRangeMin) < 1.0e-10) return;
    
    g2d = (Graphics2D)g;
    Stroke orig = g2d.getStroke();
    Dimension D =  this.getSize(); 
    setResolution(3*D.width/4,1*D.height/2);
    //setResolution(2*N,1*D.height);  
  
   if(plotok)
   {
    //System.out.println(D.width + " " + D.height);
    int boxWidth  = D.width/xCanvasPanelCount;
    int boxHeight = D.height/yCanvasPanelCount;
    if(boxWidth  <= 1)
    {
     xCanvasPanelCount = D.width/2;
     boxWidth  =  D.width/xCanvasPanelCount;
    }
    //boxHeight = 1;
    if(boxHeight <= 1)
    {
    yCanvasPanelCount = D.height/2;
    boxHeight = D.height/yCanvasPanelCount;
    }
    int xOffset = (D.width  - xCanvasPanelCount*boxWidth)/2;
    int xCanvasLocation;
    int yCanvasLocation;

    double canvasHx; double canvasHy;
    double dataHx;   double dataHy;
    canvasHx = 1.0/((double)(xCanvasPanelCount));
    canvasHy = 1.0/((double)(yCanvasPanelCount));
    dataHx   = 1.0/((double)(xDataCount-1));
    dataHy   = 1.0/((double)(yDataCount-1));

    double xLocation;
    double yLocation;

    int xDataIndex;
    int yDataIndex;

    int i; int j;
    double colorLocation;
    double dataValue;
    int colorIndex;

    double p; double q;
    double v1; double v2; double v3; double v4;
    double phi0; double phi1;

    for(i = 0; i < xCanvasPanelCount; i++)
    {
       xCanvasLocation  = xOffset + i*boxWidth;
       xLocation        = canvasHx/2.0 + ((double)i)*canvasHx;
       xDataIndex = (int)(xLocation/dataHx);
       for(j = 0; j < yCanvasPanelCount; j++)
       {
         yCanvasLocation  = j*boxHeight;
         yLocation        = canvasHy/2.0 + ((double)j)*canvasHy;
         yDataIndex = (int)(yLocation/dataHy);
     
         //do bilinear interpolation
  
         p = (xLocation - xDataIndex*dataHx)/dataHx;
         q = (yLocation - yDataIndex*dataHy)/dataHy;
         v1 = contourData[(yDataCount-1) - yDataIndex][xDataIndex]; //*yDataCount];
         v2 = contourData[(yDataCount-1) - yDataIndex][xDataIndex+1];//*yDataCount];
         v3 = contourData[(yDataCount-1)  - (yDataIndex + 1)][(xDataIndex+1)];//*yDataCount];
         v4 = contourData[(yDataCount-1) - (yDataIndex + 1)][xDataIndex];//*yDataCount];
         phi0 = v1 + p*(v2 - v1);
         phi1 = v4 + p*(v3 - v4);
         dataValue = phi0 + q*(phi1 - phi0);

         colorLocation =
           ((colorCount-1.0)*(dataValue - dataRangeMin))/(dataRangeMax - dataRangeMin);
         colorIndex = (int)colorLocation;
         if(colorIndex >= colorCount) colorIndex = colorCount-1;
         if(colorIndex < 0) colorIndex = 0;
         g.setColor(colorArray[colorIndex]);
         g.fillRect(xCanvasLocation,yCanvasLocation,boxWidth,boxHeight);
       }
    }


    g2d.setColor(Color.black);


     g2d.setStroke(dashed);
     g2d.setPaint(myGray);
     x0 = (int)(((double)line/K1)*(double)D.height);
     g2d.drawLine(0, D.height-x0, D.width, D.height-x0);
    
   }

   g2d.setStroke(orig);
   g2d.setColor(Color.black);
  }

}

