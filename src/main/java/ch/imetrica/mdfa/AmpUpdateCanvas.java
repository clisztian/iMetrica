package ch.imetrica.mdfa;

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
public class AmpUpdateCanvas extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//--- I-MDFA object attributes
    int K; //---total number of observations
    int n_rep;
    double[] gamma_hat;  //gamma_hat n_rep x (K+1)
   
    //------- canvas stuff ----------------------------
    int height, width;
    double dataMax, dataMin, dataNorm;
    Graphics2D g2d;

    DecimalFormat df;
    BasicStroke dashed,orig;
    float[] dash1;
    Color myGray, myGray2;  
    Color plotGray;
    double mse_val=0.0;
    public Font mono;

    public AmpUpdateCanvas(int w, int h, int _K, int _nrep)
    {
      // Initilize everything-------------------
      @SuppressWarnings("unused")
	  int i,k;
      K = _K; 
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
      df = new DecimalFormat("##.#######"); 

      mono  = new Font("Monospaced", Font.PLAIN, 12);
      gamma_hat = new double[(n_rep+1)*(K+1)];
      
    } 

    public void setPlotDim()
    {
        int i,k; int K1 = K+1;
        dataMax = -100000.0; dataMin = 100000.0; 

        for(i=0;i<=n_rep;i++)
        {

            for(k=0;k<K1;k++)
            {
             if(gamma_hat[i*K1+k] < dataMin) dataMin = gamma_hat[i*K1+k];
             else if(gamma_hat[i*K1+k] > dataMax) dataMax = gamma_hat[i*K1+k]; 
            } 
        }
        dataNorm = Math.abs(dataMax - dataMin);                
    }

    public void setMSEValue(double v)
    {
      mse_val = v;
      //System.out.println("new mse value = " + mse_val);
    }
    
    //----- set raw data and filtered data
    public void setGammas(double[] ts, int _K, int nrep)
    {   
       int i,k;
       K = _K; n_rep = nrep; int K1=K+1;
       gamma_hat = new double[(n_rep+1)*K1];
       
       for(i=0;i<=n_rep;i++)
       {
         for(k=0;k<K1;k++)
         {gamma_hat[i*K1+k] = ts[i*K1+k];} 
       }       
       go();
    }

    
    
 

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {  

     int j,N,k;
     int t0, t1, x0, x1;
     N = K+1;
     super.paintComponent(g);
     g2d = (Graphics2D)g; 

     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;

     setPlotDim();
     orig = new BasicStroke((float)1.3);
     g2d.setStroke(orig);

     
     for(k=0;k<=n_rep;k++)
     {
       g2d.setPaint(new Color(15+(k+1)*20,102,(k+1)*8+156));
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
     
     g.setFont(mono);
     g2d.setPaint(Color.GREEN);
     g.drawString("MSE - "+df.format(mse_val), 5, 15); 

   }


}

