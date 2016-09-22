package ch.imetrica.mdfa;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

public class ScalePanel extends JPanel implements MouseMotionListener
{


    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//------ Plot stuff
    Graphics2D g2d;
    int height, width;
    int co;
    Color[] colorArray;
    boolean plotok;
    double min, max;

    public ScalePanel()  //62 294
    { 
      co = 60; setBackground(Color.BLACK); 
      colorArray = new Color[co];  plotok =false;          
    }

    public void setExtremum(double _min, double _max) {min = _min; max = _max;}

    public void updateColor(Color[] c)
    {co = c.length; colorArray = new Color[c.length]; System.arraycopy(c, 0, colorArray, 0, c.length); plotok = true; go();}

    public void go() {repaint();}

    public void paintComponent(Graphics g)
    {
        int i; super.paintComponent(g);
        Dimension D =  this.getSize(); 
        if(plotok)
      {
        for(i=0; i<co; i++)
        {
          g.setColor(colorArray[i]);
          g.fillRect(0,2*i,D.width,2*(i+1));  
        }
      }
    }


      public void mouseDragged(MouseEvent e) 
      { }

      public void mouseMoved(MouseEvent e) 
      {this.setToolTipText("here");}
     


}
