package ch.imetrica.bayesCronos;


import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.text.*;

/*------------------------------------------------
  Canvas for plotting series
--------------------------------------------------*/
public class BayesPlot extends JPanel
{


    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int width; int height;
    Graphics2D g2d;
    Color boxColor;
    Color highlightColor, myGray;
    DecimalFormat df;
    int box_width;
    int n_bins;
    int n_samps; 
    int[] boxes;               //  int array of length n_bins, gives height of each box
    double min_data, max_data;
    int max_height;            // max height of boxes
    boolean plot_ready;
    int shade_bin_start, shade_bin_end;
    
    
    public BayesPlot()
    {
       
       //width = w; height = h; 
 
       n_bins = 150; 
       plot_ready = false;
       boxColor = new Color(103,177,246);
       highlightColor = new Color(135,100,167);
       myGray = new Color(187,191,193);
       
       shade_bin_start = 0; shade_bin_end = 0;
       df = new DecimalFormat("##.##"); 
    
    }

    public void plotReady(boolean a) {plot_ready = a;}

    //---------- Assumes the data is sorted from min to max -------------
    public void createBoxes(double[] data)
    {
       plot_ready = true;
       int i; int next = 0; double bin_val; double bin_length;
       n_samps = data.length; 
       
       boxes = new int[n_bins];
       
       //--- Find sample interval, leave out first and last two samples, outliers-------
       min_data = data[2]; max_data = data[n_samps-3]; 
       
       bin_length = (max_data - min_data)/n_bins; 
       bin_val = min_data + bin_length;
       
        for(i=2;i<n_samps-2;i++)
        {
          
          if((data[i] >= bin_val - bin_length) && data[i] <= bin_val)
          {boxes[next] = boxes[next] + 1;}  //add to bin bucket         
          else
          {next++; bin_val = bin_val + bin_length;} 
          if(next == 150) {break;}
          
          //System.out.println(data[i]);
        }
       
        max_height = 0;
        //--- Find max height of boxes --------
        for(i=0; i < n_bins; i++)
        {
         if(max_height < boxes[i]) {max_height = boxes[i]; }        
        }
    }
    
    
    public void colorShadeBoxes(int start, int end, int res)
    {
        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;

        shade_bin_start = (int)(((double)start/(double)res)*width);
        shade_bin_end = (int)(((double)end/(double)res)*width);   
        go();
    }
    
    public void go() {repaint();}

    public void paintComponent(Graphics g)
    { 
    
     int box_start = 0; int box_height = 0; int j;
        
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
     
     setBackground(new Color(0, 0, 0));
     
     if(plot_ready){
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;   
     //System.out.println(height + "  " + width);

     box_width = (int)(width/n_bins);
    
     
    
     double mid_data = (max_data + min_data)/2.0;
     g2d.setPaint(myGray);
     g.drawString((String)df.format(min_data), 5, 10);
     g.drawString((String)df.format(max_data), width - 30, 10);
     g.drawString((String)df.format(mid_data), width/2-5, 10);
     
     box_start = 0;  
     g2d.setPaint(boxColor);
     for(j = 0; j < n_bins; j++)
     {
        g2d.setPaint(boxColor);
       
        box_height = (int)(((double)boxes[j]/max_height)*height*.9);
        
       
        //System.out.println(box_height + "  " + max_height + "  "  + boxes[j] + "  " + height);
       
        g2d.draw(new Rectangle(box_start, height-box_height, box_width, box_height));
       
        //System.out.println(shade_bin_start + "  " + shade_bin_end);
        //System.out.println(box_start + "  " + (height-15) + "  " + box_width  + "  " + (height-box_height));
        if((box_start >= shade_bin_start) && (box_start <= shade_bin_end))
        {
          g2d.setPaint(highlightColor);
          g2d.fill(new Rectangle(box_start, height-box_height, box_width, box_height));
          
        }
        box_start = box_start + box_width;
     }
     
     }
    }
    
}   
