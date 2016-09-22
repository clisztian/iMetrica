package ch.imetrica.mdfa;

import java.awt.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.text.*;

public class TimeScalePanel extends JPanel
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
    int n_obs; int L;
    
    public TimeScalePanel()
    { 
  
      //setPreferredSize(new Dimension(w, h));
      setBackground(Color.BLACK);  df = new DecimalFormat("##.##");        
      myGray = new Color(67,71,73); n_obs = 300;
      font = new Font("serif", Font.PLAIN, 9);
      
    }

    public void setL(int l) {L=l;}

    public void setNObs(int n)
    {
      n_obs = n; go();
    }

    public void go() {repaint();}
    public void paintComponent(Graphics g)
    {
        int i; int t0; int p=L;  
        super.paintComponent(g);
        g2d = (Graphics2D)g; 
        Dimension D =  this.getSize(); 
        width = D.width; height = D.height;
        g2d.setPaint(myGray);
        g2d.setFont(font);
 
        int nobsP = (int)Math.floor((double)n_obs/10); 
        for(i=1; i <= nobsP; i++) 
        {
          t0 = (int)(((double)(i*10)/n_obs)*(double)width);
          g2d.setPaint(myGray);
	  g2d.drawLine(t0, 0, t0, height-20);
          g2d.setPaint(myGray);
          g.drawString((String)"" + p, t0, height - 2); 
          p = p + 10;
        } 
     }

}