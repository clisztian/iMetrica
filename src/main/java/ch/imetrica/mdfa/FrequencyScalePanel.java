package ch.imetrica.mdfa;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;

public class FrequencyScalePanel extends JPanel
{

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//------ Plot stuff
    Graphics2D g2d;
    int height, width;
    Color myGray;
    Font font;
    DecimalFormat df;
    double[] freq_ints; 
    int n_div; 
    
    public FrequencyScalePanel() // 12 294
    { 
  
      setPreferredSize(new Dimension(30, 294));
      setBackground(Color.BLACK);  df = new DecimalFormat("##.##");        
      myGray = new Color(87,91,93); n_div = 15;
      font = new Font("serif", Font.PLAIN, 8);
      freq_ints = new double[15];
    }

    public void setNDiv(double[] _ints)
    {
      freq_ints = _ints;
      n_div = _ints.length;
      go();
    }

    public void go() {repaint();}
    public void paintComponent(Graphics g)
    {
        int i;  int h;
        super.paintComponent(g);
        g2d = (Graphics2D)g; 
        Dimension D =  this.getSize(); 

        ///System.out.println("FrequencyScalePanel size = "+D.width + " " + D.height);


        h = D.height;
        g2d.setPaint(myGray);
        g2d.setFont(font);
 
        for(i=0; i < n_div; i++) 
        {g.drawString(df.format(freq_ints[i]), 5, D.height-(int)((freq_ints[i]/Math.PI)*h)); }
   
    }

}