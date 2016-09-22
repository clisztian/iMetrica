package ch.imetrica.mdfaTradingStrategies;

import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;


public class StrategyCanvas extends JTextPane
{
  
  /**
	 * 
	 */
  private static final long serialVersionUID = 1L;
  StyledDocument doc;
  Style main_style;
  Style alt_style;

  public StrategyCanvas(int w, int h)
  { 
     
       setBackground(Color.BLACK);
       setPreferredSize(new Dimension(w,h));

       doc = getStyledDocument();
     
       main_style = addStyle("Main Style",null);     
  }
  
  public void clearText()
  { 
    this.setText("");
  }
  
  public void addTextNewLine(String text, Color col)
  {
     StyleConstants.setForeground(main_style, col);
     try { doc.insertString(doc.getLength(), text + "\n", main_style); }
     catch (BadLocationException e){}
  }
     
  public void addSameLineText(String text, Color col)
  {
     StyleConstants.setForeground(main_style, col);
     try { doc.insertString(doc.getLength(), " " + text + " ", main_style); }
     catch (BadLocationException e){}
  }
     
  public void addNewLineText(String text, Color col)
  {
     StyleConstants.setForeground(main_style, col);
     try { doc.insertString(doc.getLength(), "\n" + text, main_style); }
     catch (BadLocationException e){}
  }
     
  public void setIntroduction()
  {  

     Style def = StyleContext.getDefaultStyleContext().
                        getStyle(StyleContext.DEFAULT_STYLE);  

     StyleConstants.setFontFamily(def, "SansSerif");
  
     StyleConstants.setAlignment(main_style, StyleConstants.ALIGN_LEFT);
     StyleConstants.setForeground(main_style, Color.green);
     try { doc.insertString(doc.getLength(), "Welcome to Evolution\n", main_style); }
     catch (BadLocationException e){}      
     
     
     
     
  }
      
}

