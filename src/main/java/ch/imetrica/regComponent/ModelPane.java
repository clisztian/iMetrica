package ch.imetrica.regComponent;

import javax.swing.*; 
import javax.swing.text.*; 
import java.text.*;
import java.awt.*; 
 

public class ModelPane extends JPanel
{ 

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int serial = 0; 
    JTextPane pane = null;
    ARIMAModel myModel = null;

    StyleContext sc; // = new StyleContext();
    DefaultStyledDocument doc;// = new DefaultStyledDocument(sc);
    Style mainStyle;
    Style cwStyle;
    Style heading2Style;
    DecimalFormat df;
    // Style names
    public static final String mainStyleName = "MainStyle";
    public static final String heading2StyleName = "Heading2";
    public static final String charStyleName = "ConstantWidth";


    public ModelPane(int ser, ARIMAModel my)
    { 
          
       myModel = my; 
       serial = ser; 
       setBackground(Color.BLACK);
       setPreferredSize(new Dimension(300, 170));
       df = new DecimalFormat("##.##");
       sc = new StyleContext(); 
       doc = new DefaultStyledDocument(sc);
       pane = new JTextPane(doc);
       pane.setBackground(Color.BLACK);

       createDocumentStyles(sc);
       postModelOnWall();
       add(pane, BorderLayout.CENTER);
    }


    public void postModelOnWall()
    {

      boolean n_sea,sea; 
      n_sea = myModel.n_seas; sea = myModel.seas;
   
      String vari = new String("Variance = "); 
      String arimaS = new String(" (");
      String type = new String("Model type: "); 
      String coeffs = new String("Uploaded Coefficients: ");

      type = type + myModel.model_string; 
      if(n_sea) {arimaS = arimaS+Integer.toString(myModel.dims[0])+", "+Integer.toString(myModel.dims[1])+", "+Integer.toString(myModel.dims[2]) + ")";} 
      if(sea) {arimaS = arimaS+"("+Integer.toString(myModel.dims[3])+ ", "
                        +Integer.toString(myModel.dims[4]) + ", " + Integer.toString(myModel.dims[5]) + ")_" + myModel.S;} 
   
      //vari = vari + Double.toString(df.format(myModel.var)); 
      vari = vari + df.format(myModel.var);      

      if(myModel.n_ma > 0) 
      {
        coeffs = coeffs + "MA "; if(myModel.mafixed){coeffs = coeffs + "(fixed) ";} 
      }
      if(myModel.n_ar > 0) 
      {
         coeffs = coeffs + "AR "; if(myModel.arfixed){coeffs = coeffs + "(fixed) ";}    
      } 
      if(myModel.n_sma > 0) 
      {
         coeffs = coeffs + "SMA "; if(myModel.smafixed){coeffs = coeffs + "(fixed) ";}    
      }
      if(myModel.n_sar > 0) 
      {
         coeffs = coeffs + "SAR "; if(myModel.sarfixed){coeffs = coeffs + "(fixed) ";}    
      } 
      if(myModel.fix) {coeffs = coeffs + "(all fixed) ";}
      if(myModel.hfile) {coeffs = coeffs + "\nh_t file = " + "'" + myModel.ht_file+ "'" + " \n";}

      //------ Lets write this shit----------

      final Paragraph[] content = new Paragraph[] 
      {
          new Paragraph(heading2StyleName, new Run[] {new Run(heading2StyleName, type+"\n")}), // print out model type
          new Paragraph(null, new Run[] { new Run(charStyleName, arimaS+"\n")}),
          new Paragraph(null, new Run[] { new Run(charStyleName, coeffs+"\n")}),
          new Paragraph(null, new Run[] { new Run(charStyleName, vari+"\n")})
      };

    
     addText(pane, sc, sc.getStyle(mainStyleName), content);

    }

    
    public void changeModelOnWall(ARIMAModel newModel)
    {
      myModel = newModel;
      postModelOnWall(); 
    }
         


 public static void addText(JTextPane pane, StyleContext sc, Style logicalStyle, Paragraph[] content) 
 {
   // The outer loop adds paragraphs, while the
   // inner loop adds character runs.
   int paragraphs = content.length;
   for (int i = 0; i < paragraphs; i++) {
     Run[] runs = content[i].content;
     for (int j = 0; j < runs.length; j++) {
       pane.setCharacterAttributes(
         runs[j].styleName == null ? SimpleAttributeSet.EMPTY :
                 sc.getStyle(runs[j].styleName), true);
       pane.replaceSelection(runs[j].content);
     }

     // At the end of the paragraph, add the logical style and
     // any overriding paragraph style and then terminate the 
     // paragraph with a newline.
     pane.setParagraphAttributes(SimpleAttributeSet.EMPTY, true);
     
     if (logicalStyle != null) {
       pane.setLogicalStyle(logicalStyle);
     }

     if (content[i].styleName != null) {
       pane.setParagraphAttributes(sc.getStyle(content[i].styleName), false);
     }
     
     pane.replaceSelection("\n");
   }    
 }


 public static void createDocumentStyles(StyleContext sc) 
 {  

   Style defaultStyle = sc.getStyle(StyleContext.DEFAULT_STYLE);

   // Create and add the main document style
   Style mainStyle = sc.addStyle(mainStyleName, defaultStyle);
   StyleConstants.setLeftIndent(mainStyle, 16);
   StyleConstants.setRightIndent(mainStyle, 16);
   StyleConstants.setFirstLineIndent(mainStyle, 16);
   StyleConstants.setFontFamily(mainStyle, "serif");
   StyleConstants.setFontSize(mainStyle, 12);

   // Create and add the constant width style
   Style cwStyle = sc.addStyle(charStyleName, null);
   StyleConstants.setFontFamily(cwStyle, "monospaced");
   StyleConstants.setForeground(cwStyle, Color.green);

   // Create and add the heading style
   Style heading2Style = sc.addStyle(heading2StyleName, null);
   StyleConstants.setForeground(heading2Style, Color.green);
   StyleConstants.setFontSize(heading2Style, 14);
   StyleConstants.setFontFamily(heading2Style, "sansserif");
   StyleConstants.setBold(heading2Style, true);
   StyleConstants.setLeftIndent(heading2Style, 8);
   StyleConstants.setFirstLineIndent(heading2Style, 0);
 }

 // ----------------Inner classes used to define paragraph structure
 public static class Run 
 {
   public Run(String styleName, String content) {
     this.styleName = styleName;
     this.content = content;
   }

   public String styleName;
   public String content;
 }

 //------------------------------------------------------------------
 public static class Paragraph 
 {
   public Paragraph(String styleName, Run[] content) {
     this.styleName = styleName;
     this.content = content;
   }

   public String styleName;
   public Run[] content;
 }

} 

