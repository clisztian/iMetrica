package ch.imetrica.regComponent;

import java.io.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.JCheckBox;
import javax.swing.border.*;
import javax.swing.border.LineBorder;
import java.awt.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.border.TitledBorder;
import java.text.*;
import javax.swing.filechooser.FileFilter;

/*----------------------------------------------------------------
 
   Panel for REGmodel Controls - 

----------------------------------------------------------------*/


public class REGmodelPanel extends JPanel
{



      /**
	 * 
	 */
	 private static final long serialVersionUID = 1L;
	 public REGmodelJava reg;                          // --- Main model engine
      public REGmodelCanvas reg_canvas;                 // --- Canvas for time series stuff
      //public REGcheckCanvas chk_canvas;                 // --- Canvas for diagnostic stuff

      public FileNameExtensionFilter filter;
      private Component parent;
      public boolean change;
      public int change_index, n_tabs, n_models;
      public double[] series; 
      public int n_obs; 
      public DecimalFormat df;
      public JPanel interfacePane;   
   
      public JTabbedPane modelLiszt;			// --- TextPane for showing models
      
      public JScrollBar yearscroll, monthscroll;        // --- year,month start date
      public JTextField date;                           // --- text field start date

      //--------- Panel for creating model ----------        
      private JComboBox<String> p,d,q,P,D,Q;            // --- combo boxes for SARIMA choise
      private JScrollBar innvar;                        // --- scrollbar for var 
      private double innv;
      private int i_innv;
    

      private JButton addModel; 
      private JButton deleteModel;
      private JButton compute;
      private JButton reset; 
      private JButton changeModel;

      private JButton hfileBut;
      private JCheckBox hfile;
      
      private ButtonGroup model_group;
      private JRadioButton TrendModel;                  // --- preset Trend Model
      private JRadioButton SeasonalModel;               // --- preset Seasonal Model
      private JRadioButton IrregularModel;              // --- present Irregular model
      private JRadioButton RandomWalkModel;             // --- present random walk model
      private JRadioButton TimeVarModel;                // --- present random walk time varying
      private JRadioButton airlineModel;    

      private JCheckBox predefinedSwitch;
      private JCheckBox allfix;                            // --- all of the params fixed? 
      private JCheckBox arfix;                        // --- ar params fixed? 
      private JCheckBox mafix;                        // --- ma params fixed? 
      private JCheckBox sarfix;                       // --- sar params fixed? 
      private JCheckBox smafix;                       // --- sma params fixed? 
      public int n_ar, n_ma, n_sar, n_sma,predefined;

      // -------- Attributes for the coefficients ----
      public boolean lfix;          	
      public boolean larfixed;
      public boolean lmafixed;
      public boolean lsarfixed;
      public boolean lsmafixed;  
      public boolean larima,lpre; 

      public JTextField arfield; 
      public JTextField mafield; 
      public JTextField sarfield; 
      public JTextField smafield; 

      public JButton arbut; 
      public JButton mabut; 
      public JButton sarbut; 
      public JButton smabut; 

      public double[] arcoefs;				// --- array of fixed ar params for model
      public double[] macoefs;     			// --- array of fixed ma params for model
      public double[] sarcoefs;    			// --- array of fixed ar params for model
      public double[] smacoefs;    			// --- array of fixed ma params for model 

      // ------- Panel for regression effects ------------
      public boolean trans, td, easter, outlier, tv, constant,cmpntreg;
      public JCheckBox tdBox, easterBox, outlierBox;
      public JCheckBox transBox, constantBox;
      public JComboBox<String> easterDay; 
      public JComboBox<String> tdComboBox; 
      public JComboBox<String> easterComboBox;
      public int easterDays;
      public JCheckBox tvreg; 
      public boolean input_ar,input_ma,input_sar,input_sma;
      public int easter_model;
      public int td_model;
      public Color myBlue;

      public String curDir;
      public String hfile_name; 
      public JTextField f_field;
      public JFileChooser fc;
      public boolean hfile_x;
  
      public JTextField innvarText, easterText; 
      

      public JLabel modelLabel,pLabel,qLabel,QLabel,PLabel, dLabel, DLabel;
      public JLabel premodelLabel, easterLabel, innvText;
     
      public JPanel pX,qX,dX,PX,QX,DX;
      FileFilter dfaFiles;
  

      //------------ Plotting options -----------------------
      JCheckBox[] timePlot; 
      JCheckBox plotFore; 
      boolean plotForecasts;


      public REGmodelPanel(double[] _data, Component _p, int _w, int _h)
      {

        int width = _w; int height = _h; n_models=0;
        parent = _p; plotForecasts = false;
        reg = new REGmodelJava(_data.length,100);  //start your engines
        reg_canvas = new REGmodelCanvas(width, height, _data.length);
        //reg.defaultModel(_data);
        setData(_data);
        
          dfaFiles = new ExtensionFilter("Data files", new String[] {".sig", ".dat", ".freq", ".coeff"});
          curDir=System.getProperty("user.dir") + File.separator;
          fc = new JFileChooser(curDir);
          fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
          fc.addChoosableFileFilter(dfaFiles);  

        initializePrimitives();
        setupModelBuild();
        setupRegression();
        setupInterface();
        setData(_data); 
      }

     public int getNcmpnts()
     {return reg.Ncmpnt;}
    
     public double[] GetComponent(int i, boolean fore)
     {
       int t; int nn; 
       if(fore) {nn = reg.NFore;} else {nn = reg.N;}

       double[] temp = new double[nn];
       for(t = 0; t < nn; t++)
       {temp[t] = reg.usimXbetaH[reg.NFore*i + t];}
       return temp;
     }

     public boolean isForecasting() 
     {
       return reg_canvas.plotForecasts;
     }

     public void reInitializePanel()
     {
          for(int i=1;i<11;i++)
          {timePlot[i].setSelected(false); timePlot[i].setEnabled(false);}            
     }

     //------------------- Test data --------------------------------
     public void setData(double[] _data)
     {
        n_obs = _data.length; 
        reg.setSeries(_data, n_obs);  
        series = new double[n_obs];     
        System.arraycopy(_data, 0, series, 0, n_obs);
        reg_canvas.setData(_data);

     }

     //------------------ Put any primitives in here you want to initialize ---- 
     public void initializePrimitives()
     {
        change = false; change_index = -1; myBlue = new Color(146,196,210);
        predefined = 6; n_tabs = 0; df = new DecimalFormat("##.##"); 
        n_ar=0; n_ma=1; n_sar=0; n_sma=1; larfixed = false; lmafixed = false; 
        lsmafixed = false; lsarfixed = false;  cmpntreg = false; n_tabs = 0;
        innv = 1.0; i_innv = 10; arcoefs = null;   macoefs = null;  sarcoefs = null;  smacoefs = null;
    
        input_ar = false; input_ma = false; input_sar = false; input_sma = false; hfile_x = false;
     }
 
     public void resetInputs()
     {
        String name = new String("");
        arfield.setText(name); mafield.setText(name); smafield.setText(name); sarfield.setText(name);
        input_ar = false; input_ma = false; input_sar = false; input_sma = false; hfile_x = false;
        hfile_name = name; f_field.setText(name); hfile.setSelected(false);

        larfixed = false; lmafixed = false; lsmafixed = false; lsarfixed = false; 
        arfix.setSelected(false);  mafix.setSelected(false); sarfix.setSelected(false); smafix.setSelected(false);

     } 
 

    //------------------- Test the innvar input --------------------------------
    public void test_innvar(String s)
    {   

      double e; int expint;
      try
      {
          e = Double.parseDouble(s.trim()); 
          if(e >=0.0 && e <= 4000.0) {expint=(int)(10*e); i_innv=expint; innv = e; innvarText.setText(""+df.format(e));}
          else innvarText.setText(""+df.format(innv));
      }
      catch (NumberFormatException nfe)
      {System.out.println("NumberFormatException: " + nfe.getMessage()); innvarText.setText(""+df.format(innv));}
    }


    //------------------- Test the coefficient input field --------------------------------

    public void test_field(String s, int modelx)
    {
       int i;
       String delims = "[ ]+";
       String[] tokens = s.split(delims);

       if((modelx == 0) && (p.getSelectedIndex() == tokens.length))
       {
          arcoefs = new double[tokens.length]; n_ar = tokens.length;
          for(i=0;i<tokens.length;i++) {arcoefs[i] = Double.valueOf(tokens[i]).doubleValue();}
          input_ar = true; System.out.println("Uploaded AR coefficients");
       }
       else if((modelx == 1) && (q.getSelectedIndex() == tokens.length))
       {
          macoefs = new double[tokens.length]; n_ma = tokens.length;
          for(i=0;i<tokens.length;i++) {macoefs[i] = Double.valueOf(tokens[i]).doubleValue();} 
          input_ma = true; System.out.println("Uploaded MA coefficients");  
       }
       else if((modelx == 2) && (P.getSelectedIndex() == tokens.length))
       {
          sarcoefs = new double[tokens.length]; n_sar = tokens.length;
          for(i=0;i<tokens.length;i++) {sarcoefs[i] = Double.valueOf(tokens[i]).doubleValue();}
          input_sar = true; System.out.println("Uploaded SAR coefficients");
       }
       else if((modelx == 3) && (Q.getSelectedIndex() == tokens.length))
       {
          smacoefs = new double[tokens.length]; n_sma = tokens.length;
          for(i=0;i<tokens.length;i++) {smacoefs[i] = Double.valueOf(tokens[i]).doubleValue();}
          input_sma = true;  System.out.println("Uploaded SMA coefficients");
       }           
       else
       {
        System.out.println("Trouble reading parameters. Number of entered parameters not sufficient or incorrect format.");
        System.out.println("Separate them with space, like .45 .81 .10 ... "); 
       }
    }

 
    //------------- Gets the H_file and then adds ones for the forecasts -----------

    public boolean checkHFile(File file)
    {
       /*
       String strline; Double D;
       String[] tokens; String delims = "[ ]+";
       int n_toks; int len;
       double[] values = new double[500];      
       double val = 0;
       int i = 0; len = n_obs;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));

         
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           D = new Double(tokens[n_toks-1]); //take only the value, no dates
           val = D.doubleValue();

           if(i < 500) {values[i] = val; i++;}
           else 
           {System.out.println("Maximum times series length is 500"); break;} 
           
         } 
         din.close(); len = i;
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}

        //---- now re-output with 24 ones
        try{  
           PrintWriter out = new PrintWriter(new FileWriter(file));

           for(i=0; i < len; i++) {out.println(values[i]);}
           for(i=0; i < 24; i++) {out.println(1.0);}

           out.close(); System.out.println("Series successfully saved in " + file);
        } catch (IOException e) {e.printStackTrace();}
       */
       return true;
    }


   public void setComboBoxs()
   {
         int i;
         for(i=0;i<15;i++) {p.addItem(Integer.toString(i));} p.setSelectedIndex(0);
         for(i=0;i<15;i++) {q.addItem(Integer.toString(i));} q.setSelectedIndex(1);
         for(i=0;i<5;i++) {P.addItem(Integer.toString(i));} P.setSelectedIndex(0);
         for(i=0;i<5;i++) {Q.addItem(Integer.toString(i));} Q.setSelectedIndex(1);
       
          D.addItem("0"); D.addItem("1"); D.addItem("2"); D.setSelectedIndex(1);     
          d.addItem("0"); d.addItem("1"); d.addItem("2"); d.setSelectedIndex(1);

   }

    //------------------- Read in coefficient data from file ----------------------------
 
    public boolean readCoeffData(File file, int modelx)
    {
       int n;
       String name = file.getName(); String strline; Double D;  String[] tokens; String delims = "[ ]+";
       int n_toks; double[] values = new double[50]; double val = 0; int i = 0; boolean suc = false;

       try{
          
         FileInputStream fin = new FileInputStream(file);  DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {
           tokens = strline.split(delims);  n_toks = tokens.length; 
           if(n_toks == 0) {System.out.println("End of file"); break;}
  
           D = new Double(tokens[n_toks-1]); //take only the value, no dates
           val = D.doubleValue();

           if(i < 50) {values[i] = val; i++;}
           else  {System.out.println("Maximum coefficient length is 50"); break;}           
         } 
         din.close(); suc = true;
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}        

        if(suc)
        {
         if(modelx == 0)
         { 
          n = (int)values[0]; if(n != n_ar) {System.out.println("Changing number of AR parameters to " + n + ".\n"); n_ar = n;}
          arcoefs = new double[n_ar]; for(i=0;i<n_ar;i++) {arcoefs[i] = values[i+1];} arfield.setText(name); input_ar = true;
         }
         else if(modelx == 1)
         {  
          n = (int)values[0]; if(n != n_ma) {System.out.println("Changing number of MA parameters to " + n + ".\n"); n_ma = n;}
          macoefs = new double[n_ma]; for(i=0;i<n_ma;i++) {macoefs[i] = values[i+1];} mafield.setText(name); input_ma = true;
         }
         else if(modelx == 2)
         { 
          sarcoefs = new double[n_sar]; for(i=0;i<n_sar;i++) {sarcoefs[i] = values[i];} sarfield.setText(name); input_sar = true;
         }
         else if(modelx == 3)
         { 
          smacoefs = new double[n_sma]; for(i=0;i<n_sma;i++) {smacoefs[i] = values[i];} smafield.setText(name); input_sma = true;
         }
        }
        else {System.out.println("Couldn't access coefficient names. Bad format perhaps.");}

        return suc;
    }


    public void setupModelBuild()
    { 
         
          int i;
          addModel = new JButton("Add Model");     addModel.setToolTipText("Add defined model to the modeling list.");
          deleteModel = new JButton("Delete"); deleteModel.setToolTipText("Delete all models to the modeling list.");
          addModel.addActionListener(new MyActionListener()); deleteModel.addActionListener(new MyActionListener());
          compute = new JButton("Compute");  compute.setToolTipText("Estimate and forecast the data using the given regression and model components.");
          compute.addActionListener(new MyActionListener()); deleteModel.setEnabled(false); 
          reset = new JButton("Reset"); reset.setToolTipText("Reset input coefficients to null"); reset.addActionListener(new MyActionListener());
          changeModel = new JButton("Modify"); 
          changeModel.setToolTipText("Modify the choice current model"); changeModel.addActionListener(new MyActionListener());
          changeModel.setEnabled(false);

          //getData = new JButton("Get Data"); getData.setToolTipText("Import time series data from file.");

          //---------------- Radio Buttons to switch between model selection ------------
          //model_switch = new ButtonGroup(); 
          larima=true;lpre=false;
          //arimaSwitch  = new JRadioButton("SARIMA:",true);  
          //arimaSwitch.setToolTipText("Choose SARIMA model using dimension controls.");               
          predefinedSwitch = new JCheckBox("Predefined:"); predefinedSwitch.setToolTipText("Select from an assortment of pre-defined models.");
          predefinedSwitch.setHorizontalTextPosition(JMenuItem.LEFT);           
          
          predefinedSwitch.addItemListener(new MyItemListener());


          // ------ General ARIMA model build --------------------------------------
          premodelLabel = new JLabel("Predefined:"); 
          modelLabel = new JLabel("SARIMA Models:");  
          pLabel = new JLabel("  p:"); pLabel.setToolTipText("Nonseasonal autoregressive order");
          qLabel = new JLabel("  q:"); qLabel.setToolTipText("Nonseasonal moving-average order");
          QLabel = new JLabel("  Q:"); QLabel.setToolTipText("Seasonal moving-average order");
          PLabel = new JLabel("  P:"); PLabel.setToolTipText("Seasonal autoregressive order");
          dLabel = new JLabel("  d:"); dLabel.setToolTipText("Nonseasonal differencing order");
          DLabel = new JLabel("  D:"); DLabel.setToolTipText("Seasonal differencing order");
          
          p = new JComboBox<String>();          
          q = new JComboBox<String>();  
          P = new JComboBox<String>(); 
          Q = new JComboBox<String>();          
          D = new JComboBox<String>();
          d = new JComboBox<String>();
          setComboBoxs();
          p.addActionListener(new MyActionListener()); P.addActionListener(new MyActionListener());
          q.addActionListener(new MyActionListener()); Q.addActionListener(new MyActionListener());
          

          //---------------- Button groupe for predefined ----------------------------------------------------------
          model_group = new ButtonGroup();
          TrendModel = new JRadioButton("Trend:");           TrendModel.setHorizontalTextPosition(JMenuItem.LEFT);               
          SeasonalModel = new JRadioButton("Seas.:");     SeasonalModel.setHorizontalTextPosition(JMenuItem.LEFT);            
          IrregularModel = new JRadioButton("Irreg.:");   IrregularModel.setHorizontalTextPosition(JMenuItem.LEFT);            
          RandomWalkModel = new JRadioButton("Walk:");  RandomWalkModel.setHorizontalTextPosition(JMenuItem.LEFT);           
          TimeVarModel = new JRadioButton("T.Var.:");     TimeVarModel.setHorizontalTextPosition(JMenuItem.LEFT);           
          airlineModel = new JRadioButton("Airline:",true);  airlineModel.setHorizontalTextPosition(JMenuItem.LEFT);
          TrendModel.setToolTipText("Add an MA(2) trend model"); 
          SeasonalModel.setToolTipText("Add an MA(11) seasonal model"); 
          IrregularModel.setToolTipText("Add an irregular/error model"); 
          RandomWalkModel.setToolTipText("Add an Random Walk model"); 
          TimeVarModel.setToolTipText("Add a time-varying trading day (6 component) model"); 
          airlineModel.setToolTipText("Add an standard (0,1,1)(0,1,1)_S Airline model"); 

          TrendModel.setActionCommand("0"); SeasonalModel.setActionCommand("1");
          IrregularModel.setActionCommand("2"); RandomWalkModel.setActionCommand("3");
          TimeVarModel.setActionCommand("4"); airlineModel.setActionCommand("6");

          model_group.add(TrendModel); model_group.add(SeasonalModel); model_group.add(IrregularModel); 
          model_group.add(TimeVarModel); model_group.add(RandomWalkModel); model_group.add(airlineModel); 

          TrendModel.addActionListener(new MyActionListener());  SeasonalModel.addActionListener(new MyActionListener());
          IrregularModel.addActionListener(new MyActionListener());  RandomWalkModel.addActionListener(new MyActionListener());
          TimeVarModel.addActionListener(new MyActionListener());  airlineModel.addActionListener(new MyActionListener());
 
          TrendModel.setEnabled(false); SeasonalModel.setEnabled(false); IrregularModel.setEnabled(false);
          RandomWalkModel.setEnabled(false); TimeVarModel.setEnabled(false); airlineModel.setEnabled(false);
 

          //---------------- Fixed Buttons -------------------------------------------------
          arfix = new JCheckBox("AR:");    arfix.setHorizontalTextPosition(JMenuItem.LEFT);  arfix.addItemListener(new MyItemListener());
          mafix = new JCheckBox("MA:");    mafix.setHorizontalTextPosition(JMenuItem.LEFT);  mafix.addItemListener(new MyItemListener());
          sarfix = new JCheckBox("SAR:");  sarfix.setHorizontalTextPosition(JMenuItem.LEFT); sarfix.addItemListener(new MyItemListener());
          smafix = new JCheckBox("SMA:");   smafix.setHorizontalTextPosition(JMenuItem.LEFT); smafix.addItemListener(new MyItemListener());
          allfix = new JCheckBox("Fix All:");  allfix.setHorizontalTextPosition(JMenuItem.LEFT); allfix.addItemListener(new MyItemListener());
           
          arfix.setToolTipText("Fix the nonseasonal coefficients"); sarfix.setToolTipText("Fix the seasonal coefficients"); 
          mafix.setToolTipText("Fix the nonseasonal coefficients"); smafix.setToolTipText("Fix the seasonal coefficients"); 
          allfix.setToolTipText("Fix all the parameters");


 
          //------------------ Add Fields for coefficients -------------------
          arfield = new JTextField(12);  
          arfield.setToolTipText("File name for AR coefficients or enter them in individually separated by space. Press Enter when done.");
          mafield = new JTextField(12); 
          mafield.setToolTipText("File name for for MA coefficients or enter them in individually separated by space. Press Enter when done.");
          sarfield = new JTextField(12); 
          sarfield.setToolTipText("File name for Seasonal AR coefficients or enter them in individually separated by space. Press Enter when done.");
          smafield = new JTextField(12); 
          smafield.setToolTipText("File name for Seasonal MA coefficients or enter them in individually separated by space. Press Enter when done.");
          arfield.setText(""); 
          arfield.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_field(arfield.getText(),0);}} );  
          mafield.setText(""); 
          mafield.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_field(mafield.getText(),1);}} );  
          sarfield.setText(""); 
          sarfield.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
           if(e.getKeyCode() == KeyEvent.VK_ENTER) test_field(sarfield.getText(),2);}} );            
          smafield.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
          smafield.setText(""); 
          if(e.getKeyCode() == KeyEvent.VK_ENTER) test_field(smafield.getText(),3);}} );  

          arbut = new JButton("AR file"); arbut.setToolTipText("Choose file for AR coefficients.");
          mabut = new JButton("MA file"); mabut.setToolTipText("Choose file for MA coefficients.");
          sarbut = new JButton("SAR file");  sarbut.setToolTipText("Choose file for SAR coefficients.");
          smabut = new JButton("SMA file");  smabut.setToolTipText("Choose file for SMA coefficients.");

          arbut.addActionListener(new MyActionListener()); arbut.setEnabled(false);
          mabut.addActionListener(new MyActionListener()); mabut.setEnabled(true);
          sarbut.addActionListener(new MyActionListener()); sarbut.setEnabled(false);
          smabut.addActionListener(new MyActionListener()); smabut.setEnabled(true);

          //--------------- Innovation variance ---------------------------------------------
          innvar = new JScrollBar(JScrollBar.HORIZONTAL,0,10,0,4000); innvar.setValue(i_innv);
          innvar.addAdjustmentListener(new MyAdjustmentListener()); innvar.setUnitIncrement(1); 
          innvarText = new JTextField(3); innvarText.setColumns(4);
          innvarText.addKeyListener(new KeyAdapter() {public void keyPressed(KeyEvent e) {
            if(e.getKeyCode() == KeyEvent.VK_ENTER) test_innvar(innvarText.getText());}} );  
              innvarText.setText("" + innv);
          innvText = new JLabel("Variance:");
          //----------------------------------------------------------------------------------

          hfileBut = new JButton("h-File"); hfileBut.addActionListener(new MyActionListener());
          hfileBut.setToolTipText("Choose file for the h_t sequence."); hfileBut.setEnabled(false);
          hfile = new JCheckBox("h_t:"); hfile.setToolTipText("Provide an h_t sequence for the model");
          hfile.addItemListener(new MyItemListener()); hfile.setHorizontalTextPosition(JMenuItem.LEFT); 
          f_field = new JTextField(5); 

          pX = new JPanel(); pX.setLayout(new FlowLayout()); pX.add(pLabel); pX.add(p); 
          dX = new JPanel(); dX.setLayout(new FlowLayout()); dX.add(dLabel); dX.add(d);    
          qX = new JPanel(); qX.setLayout(new FlowLayout()); qX.add(qLabel); qX.add(q); 
          PX = new JPanel(); PX.setLayout(new FlowLayout()); PX.add(PLabel); PX.add(P);    
          DX = new JPanel(); DX.setLayout(new FlowLayout()); DX.add(DLabel); DX.add(D); 
          QX = new JPanel(); QX.setLayout(new FlowLayout()); QX.add(QLabel); QX.add(Q);     


          //------- setup the model panel ------------------
          modelLiszt = new JTabbedPane(JTabbedPane.TOP);
          modelLiszt.setPreferredSize(new Dimension(310, 180));
	  modelLiszt.setBorder(new BevelBorder(BevelBorder.RAISED));
 
          plotFore = new JCheckBox("Forecasts:"); 
          plotFore.setSelected(false); plotFore.setEnabled(false);
          plotFore.setHorizontalTextPosition(JMenuItem.LEFT);
          plotFore.addItemListener(new MyItemListener());     

          timePlot = new JCheckBox[11];
          timePlot[0] = new JCheckBox("Y(t):"); timePlot[0].setSelected(true); 
          timePlot[0].setHorizontalTextPosition(JMenuItem.LEFT);
          timePlot[0].addItemListener(new MyItemListener());
          // ------ set up time plots -------------------------------
          for(i=1;i<11;i++)
          {
           timePlot[i] = new JCheckBox("  X_" + Integer.toString(i) + "(t):"); 
           timePlot[i].setSelected(false); timePlot[i].setEnabled(false);
           timePlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
           timePlot[i].addItemListener(new MyItemListener());
          }            



     }


     //------------------- Setup Regression stuff  --------------------------------
     public void setupRegression()
     {

       int i;

       easterDay = new JComboBox<String>(); easterDay.setToolTipText("Check and adjust for Easter-Day effects in data using x-day Easter regressors."); 

       for(i=1;i<16;i++) {easterDay.addItem(Integer.toString(i));}  easterDay.setSelectedIndex(0); 
       easterDay.addActionListener(new MyActionListener());
   
       tdComboBox = new JComboBox<String>(); 
       tdComboBox.setToolTipText("Choose model to associate trading day regression with. Model 0 is implies separate regression term.");
       for(i = 0; i < 1; i++) {tdComboBox.addItem("Separate");}  
       tdComboBox.setSelectedIndex(0);

       easterComboBox = new JComboBox<String>(); 
       easterComboBox.setToolTipText("Choose model to associate easter regression with. Model 0 is implies separate regression term.");
       for(i = 0; i < 1; i++) {easterComboBox.addItem("Separate");}  
       easterComboBox.setSelectedIndex(0);

       tdBox = new JCheckBox("TD:"); tdBox.setHorizontalTextPosition(JLabel.LEFT);
       outlierBox = new JCheckBox("Outliers:");  outlierBox.setHorizontalTextPosition(JLabel.LEFT);

       constantBox = new JCheckBox("Constant:"); constantBox.setHorizontalTextPosition(JLabel.LEFT); 
       easterBox = new JCheckBox("Easter:"); easterBox.setHorizontalTextPosition(JLabel.LEFT);  
       transBox = new JCheckBox("BC-Trans:"); transBox.setHorizontalTextPosition(JLabel.LEFT);  
       transBox.setToolTipText("Set Box-Cox Transformation"); 
       easterBox.setToolTipText("Check and adjust for Easter-Day effects in data using x-day Eastor regressors."); 
       tdBox.setToolTipText("Check and adjust for Trading-Day effects in data. Select models to apply trading day or select time varying trading day.");
       outlierBox.setToolTipText("Check and adjust for outliers in data.");
       constantBox.setToolTipText("Adjust for a constant regression term");
       constantBox.addItemListener(new MyItemListener());

       tdBox.addItemListener(new MyItemListener()); outlierBox.addItemListener(new MyItemListener()); 
       transBox.addItemListener(new MyItemListener()); easterBox.addItemListener(new MyItemListener()); 
       tvreg = new JCheckBox("TV:"); tvreg.setSelected(false); tvreg.addItemListener(new MyItemListener());
       tvreg.setToolTipText("Employ time varying trading day"); tvreg.setHorizontalTextPosition(JLabel.LEFT); 

     }


 
    // ------ Gathers up all the current information and sends it the the reg engine---- 

    public void addModel()
    {
 
     
      ARIMAModel model;
      if(n_tabs < 10)
      {
       if(larima)
       {        
         model = new ARIMAModel();          //---- an ARIMA (0,2,1) model
         model.setDimensions(p.getSelectedIndex(),d.getSelectedIndex(),q.getSelectedIndex(),
                             P.getSelectedIndex(),D.getSelectedIndex(),Q.getSelectedIndex()); 
         model.setVar(innv);
       }
       else
       {model = new ARIMAModel(predefined, innv);}

       //------ If h_file, set h_file -----------------------
       if(hfile_x) {System.out.println("attaching file " + hfile_name); model.setHtFile(hfile_x, hfile_name);} 

       //------ If any parameters fixed -----------------------
       model.setAllFixed(lfix); 
       model.fixAR(larfixed); 
       model.fixMA(lmafixed); 
       model.fixSAR(lsarfixed); 
       model.fixSMA(lsmafixed); 
 
       //------ Set any coefficients -----------------------
       if(input_ar) {model.setAR(arcoefs);}        
       if(input_ma) {model.setMA(macoefs);}  
       if(input_sar) {model.setSAR(sarcoefs);}        
       if(input_sma) {model.setSMA(smacoefs);}  
    
       reg.addModel(model);
       n_models++;

       //----- add it to combo boxes
       tdComboBox.addItem(Integer.toString(reg.models.size()));
       easterComboBox.addItem(Integer.toString(reg.models.size()));

       //Now create new panel and place on model liszt 
       ModelPane newMod = new ModelPane(n_models-1,model);
       modelLiszt.addTab("M"+Integer.toString(n_tabs+1), newMod);
  
       n_tabs = modelLiszt.getTabCount();
       changeModel.setEnabled(true);  
       deleteModel.setEnabled(true);
       resetInputs();
       
      }
  
    }

    public void makeTDModels()
    {
       boolean make = false;
       if(!tv)
       {
        tv=true;
        ARIMAModel model = new ARIMAModel(4,innv);   
        make = reg.addTV(model);
        if(make)
        {
         reg.setTradingDay(true, Integer.toString(n_models+1), Integer.toString(n_models+2), 
                                Integer.toString(n_models+3), Integer.toString(n_models+4), 
                                Integer.toString(n_models+5), Integer.toString(n_models+6));
         n_models++;
  
         //Now create new panel and place on model liszt 
         ModelPane newMod = new ModelPane(n_models-1,model);
         modelLiszt.addTab("M"+Integer.toString(n_tabs+1), newMod);  
         n_tabs = modelLiszt.getTabCount();
         TimeVarModel.setEnabled(false); 
        }
       }
    }

    public void deleteTDModels()
    { 
       tv = false;
       reg.deleteTV();
       TimeVarModel.setEnabled(true); 
       n_models--;
    }


    public void changeState(int i)
    {

       ARIMAModel model = reg.getModel(i);
       int mode = model.model_type;
 
       if((mode >=0) && (mode <=6))  //--- predefined model interface
       {
         predefinedSwitch.setSelected(true); lpre=true; larima=false;      
         if(mode == 0) {TrendModel.setEnabled(true);} 
         else if(mode == 1) { SeasonalModel.setEnabled(true);}
         else if(mode == 2) { IrregularModel.setEnabled(true);}
         else if(mode == 3) { RandomWalkModel.setEnabled(true);}     
         else if(mode == 4) { TimeVarModel.setEnabled(true);}        
         else if(mode == 6) { airlineModel.setEnabled(true);}

         p.setEnabled(false);  d.setEnabled(false); q.setEnabled(false);
         P.setEnabled(false);  D.setEnabled(false); Q.setEnabled(false);
 
       }
       else
       {
          //arimaSwitch.setSelected(true); 
          larima = true; lpre = false;
          p.setEnabled(true);  d.setEnabled(true); q.setEnabled(true);
          P.setEnabled(true);  D.setEnabled(true); Q.setEnabled(true);
 
          p.setSelectedIndex(model.dims[0]); d.setSelectedIndex(model.dims[1]); q.setSelectedIndex(model.dims[2]);
          P.setSelectedIndex(model.dims[3]); D.setSelectedIndex(model.dims[4]); Q.setSelectedIndex(model.dims[5]);
       }
 
       arfix.setSelected(model.arfixed); mafix.setSelected(model.mafixed);
       sarfix.setSelected(model.sarfixed); mafix.setSelected(model.smafixed);

       //------ If h_file, set h_file -----------------------
       hfile_x = model.hfile;
       if(hfile_x) {f_field.setText(model.ht_file);} 

       if(model.n_ar > 0)
       {arfield.setText(""+df.format(model.arcoefs[0]));}
       if(model.n_ma > 0)
       {mafield.setText(""+df.format(model.macoefs[0]));}
       if(model.n_sar > 0)
       {sarfield.setText(""+df.format(model.sarcoefs[0]));}
       if(model.n_sma > 0)
       {smafield.setText(""+df.format(model.smacoefs[0]));}

       innv = model.var; innvarText.setText(""+ df.format(innv));
       i_innv = (int)innv*10; innvar.setValue(i_innv);

    }

    public void updateModel(int i)
    {
       
       ARIMAModel model = reg.getModel(i);

       if(larima)
       {                 
         model.setDimensions(p.getSelectedIndex(),d.getSelectedIndex(),q.getSelectedIndex(),
                             P.getSelectedIndex(),D.getSelectedIndex(),Q.getSelectedIndex()); 
         model.setVar(innv);
       }
       else
       {model = new ARIMAModel(predefined, innv);}

       //------ If h_file, set h_file -----------------------
       if(hfile_x) {model.setHtFile(hfile_x, hfile_name);} 

       //------ If any parameters fixed -----------------------
       model.setAllFixed(lfix); 
       model.fixAR(larfixed); 
       model.fixMA(lmafixed); 
       model.fixSAR(lsarfixed); 
       model.fixSMA(lsmafixed); 
 
       //------ Set any coefficients -----------------------
       if(input_ar) {model.setAR(arcoefs);}        
       if(input_ma) {model.setMA(macoefs);}  
       if(input_sar) {model.setSAR(sarcoefs);}        
       if(input_sma) {model.setSMA(smacoefs);}  

       reg.replaceModel(i, model);  

    }
  

    //--------------- Compute the reg Component model -----------------
    public void computeRegComp()
    {
        int ind1, i; 
        ind1 = easterComboBox.getSelectedIndex(); 
        String tdm;

        //---- set up regcomponent stuff -------
        if(constant)    {reg.setConstant(true, "t"); reg.setRegression(true);}      
        if(trans) {reg.setPowerTransform(true, 0.0);}

        if(easter) {reg.setEaster(true, easterDay.getSelectedIndex()+1, (String)easterComboBox.getItemAt(ind1)); reg.setRegression(true);}
        if(cmpntreg) {reg.setCmpntLom(true, "0"); reg.setRegression(true);}
        if(td && !tv)  
        {
           reg.setRegression(true);
           td_model = tdComboBox.getSelectedIndex();  
           tdm = (String)tdComboBox.getItemAt(td_model);
           reg.setTradingDay(true, tdm, tdm, tdm, tdm, tdm, tdm);
        }
       
        //---- print nml file ------------------ 
        reg.setNMLFile();
        //---- compute -------------------------
        System.out.println("-------- Computing regCmpnt model--------");
        reg.computeRegComponent();    
        reg_canvas.setNumberComps(reg.Ncmpnt);

        //------- copy data to canvas for painting ----------------------
        reg_canvas.setComponents(reg.usimXbetaH);
        reg_canvas.setForecasts(reg.y_frcstl,reg.y_frcstm,reg.y_frcsth,reg.y_xB_frcst);
        reg_canvas.setQvariances(reg.usimQ, reg.usimQbar);   

        //------- repaint the canvas ------------------------------------
        reg_canvas.setPlotDim();
        reg_canvas.go();

        //--------- update plot checkboxes------------------------------
        for(i=1; i <= reg.Ncmpnt; i++)
        {timePlot[i].setEnabled(true);}
        for(i = reg.Ncmpnt+1; i < 11; i++)
        {timePlot[i].setEnabled(false);}
        plotFore.setEnabled(true);

    }
 
   //-------------------------- Setup up interface ------------------------------


   public void changeHighlight(int i)
   {     
     reg_canvas.changeHighlight(i);  reg_canvas.go(); 
   }

    class MyActionListener implements ActionListener {
        public void actionPerformed(ActionEvent e)
        {
           int val; boolean check; File file; 
           
           if(e.getSource() == reset)
           {
              resetInputs();
           } 
           else if(e.getSource() == addModel) 
           {                     
              if(change) 
              {
                 updateModel(change_index); 
                 change = false; addModel.setText("Add Model");
                 deleteModel.setEnabled(true);
                 reset.setEnabled(true);
              }
              else
              {
                if(predefined == 4)
                {makeTDModels();}                 
                else
                {addModel();}                           
              }
           }
           else if(e.getSource() == deleteModel) 
           {            
             if(modelLiszt.getTabCount() > 0) 
             {
               //--- get selected model
               val = modelLiszt.getSelectedIndex();

               //--- find out if TD ------
               if(reg.getModel(val).model_type == 4)      
               {
                deleteTDModels();         
                modelLiszt.remove(val);
                n_tabs = modelLiszt.getTabCount();
               }
               else
               {
                //--- erase in combo boxes
                tdComboBox.removeItemAt(val+1);
                easterComboBox.removeItemAt(val+1);
                 
                n_models=n_models-1;
                reg.deleteModel(val);
                modelLiszt.remove(val);
                n_tabs = modelLiszt.getTabCount();
               }
             }
           }
           else if(e.getSource() == changeModel) // ---- in change 
           {
             change_index = modelLiszt.getSelectedIndex();
             changeState(change_index); change = true;  

              addModel.setText("Done"); 
              deleteModel.setEnabled(false);
              reset.setEnabled(false);

           } 
           else if(e.getSource() == compute) 
           {
              computeRegComp();
           }           
           else if(e.getSource() == arbut)
           {
                         
             val = fc.showOpenDialog(parent);
 
             if (val == JFileChooser.APPROVE_OPTION) {check = readCoeffData(fc.getSelectedFile(), 0);} 
             else {System.out.println("Open command cancelled by user.");}
           }
           else if(e.getSource() == mabut)
           {
               
             val = fc.showOpenDialog(parent);
             if (val == JFileChooser.APPROVE_OPTION)  {check = readCoeffData(fc.getSelectedFile(), 1);} 
             else {System.out.println("Open command cancelled by user.");}
           }             
           else if(e.getSource() == sarbut)
           {
             
             val = fc.showOpenDialog(parent);
             if (val == JFileChooser.APPROVE_OPTION)  {check = readCoeffData(fc.getSelectedFile(), 2);} 
             else {System.out.println("Open command cancelled by user.");}
           }
           else if(e.getSource() == smabut)
           {
             val = fc.showOpenDialog(parent);
             if (val == JFileChooser.APPROVE_OPTION)  {check = readCoeffData(fc.getSelectedFile(), 3);} 
             else {System.out.println("Open command cancelled by user.");}
           }
           else if(e.getSource() == easterDay) 
           {
             easterDays = easterDay.getSelectedIndex() + 1;
           }
           else if(e.getSource() == hfileBut)
           {
              
             val = fc.showOpenDialog(parent);
             if (val == JFileChooser.APPROVE_OPTION) 
             {
               file = fc.getSelectedFile();
               check = checkHFile(file);     //make sure its legit
               if(check) 
               {hfile_name = file.getName();f_field.setText(hfile_name);}  
             } 
             else {System.out.println("Open command cancelled by user.");}   
           }
           else if(e.getSource() == p)
           {
             
             n_ar = p.getSelectedIndex(); System.out.println(""+n_ar);  if(n_ar > 0) {arbut.setEnabled(true);} else {arbut.setEnabled(false);}
             resetInputs();
           }
           else if(e.getSource() == q)
           {
             n_ma = q.getSelectedIndex(); if(n_ma > 0) {mabut.setEnabled(true);} else {mabut.setEnabled(false);} 
             resetInputs();
           }
           else if(e.getSource() == Q)
           {
             n_sma = Q.getSelectedIndex(); if(n_sma > 0) {smabut.setEnabled(true);} else {smabut.setEnabled(false);}
             resetInputs();
           }
           else if(e.getSource() == P)
           { 
             n_sar = P.getSelectedIndex();  if(n_sar > 0) {sarbut.setEnabled(true);} else {sarbut.setEnabled(false);}
             resetInputs();
           }  

           else if(e.getSource() == TrendModel)
           {
             predefined = 0;   n_ar = 0; n_ma = 2; n_sar = 0; n_sma = 0;
             mabut.setEnabled(true); arbut.setEnabled(true);
             smabut.setEnabled(false); sarbut.setEnabled(false);
             mafix.setEnabled(true); arfix.setEnabled(true);
             smafix.setEnabled(false); sarfix.setEnabled(false);
             resetInputs(); 
           }
           else if(e.getSource() == SeasonalModel)
           {
            predefined = 1;  n_ar = 0; n_ma = 11; n_sar = 0; n_sma = 0;
            mabut.setEnabled(true); arbut.setEnabled(true); 
            smabut.setEnabled(false); sarbut.setEnabled(false);
            mafix.setEnabled(true); arfix.setEnabled(true);
            smafix.setEnabled(false); sarfix.setEnabled(false); resetInputs();
             
           }
           else if(e.getSource() == RandomWalkModel)
           {
            predefined = 3;   n_ar = 0; n_ma = 0; n_sar = 0; n_sma = 0;
            mabut.setEnabled(false); arbut.setEnabled(false);
            smabut.setEnabled(false); sarbut.setEnabled(false);
            mafix.setEnabled(false); arfix.setEnabled(false);
            smafix.setEnabled(false); sarfix.setEnabled(false); resetInputs();
            
           } 
           else if(e.getSource() == IrregularModel)
           {
            predefined = 2;  n_ar = 0; n_ma = 0; n_sar = 0; n_sma = 0;
            mabut.setEnabled(true); arbut.setEnabled(true);
            smabut.setEnabled(false); sarbut.setEnabled(false); 
            mafix.setEnabled(true); arfix.setEnabled(true);
            smafix.setEnabled(false); sarfix.setEnabled(false); resetInputs();
            
           }
           else if(e.getSource() == airlineModel)
           {
            predefined = 6;   n_ar = 0; n_ma = 1; n_sar = 0; n_sma = 1;
            mabut.setEnabled(true); arbut.setEnabled(false);
            smabut.setEnabled(true); sarbut.setEnabled(false);
            mafix.setEnabled(true); arfix.setEnabled(false);
            smafix.setEnabled(true); sarfix.setEnabled(false); resetInputs();
            
           } 
           else if(e.getSource() == TimeVarModel)
           {
            
            predefined = 4;  n_ar = 0; n_ma = 0; n_sar = 0; n_sma = 0;
            mabut.setEnabled(false); arbut.setEnabled(false);
            smabut.setEnabled(false); sarbut.setEnabled(false);
            mafix.setEnabled(false); arfix.setEnabled(false); 
            smafix.setEnabled(false); sarfix.setEnabled(false); resetInputs();
             
           }

 
       }
    }


    class MyAdjustmentListener implements AdjustmentListener {
        public void adjustmentValueChanged(AdjustmentEvent e)
        {
           if(e.getAdjustable() == innvar)
           {
             i_innv = innvar.getValue(); innv = i_innv*.1; innvarText.setText(""+ df.format(innv));
           }
       }      
    }

      class MyItemListener implements ItemListener {
        public void itemStateChanged(ItemEvent e)
        {
         int i; boolean sel; 
         Object source = e.getItemSelectable();
         if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
         else{sel = true;}

         if(source == hfile)
         {hfile_x = sel; if(sel) {hfileBut.setEnabled(true);} else{hfileBut.setEnabled(false);}}
         else if(source == easterBox)
         {easter = sel;}
         else if(source == transBox)
         {trans = sel;}
	 else if(source == tdBox)
         {td = sel; TimeVarModel.setEnabled(sel);} 
         else if(source == constantBox)
         {constant = sel;} 
         else if(source == outlierBox)
         {outlier = sel;}
         else if(source == arfix)
         {
           larfixed=sel;   
         }
         else if(source == mafix)
         {  
           lmafixed=sel;
         } 
         else if(source == sarfix)
         {  
           lsarfixed=sel;  
         } 
         else if(source == smafix)
         {  
           lsmafixed=sel;  
         } 
         else if(source == allfix)
         {
           lfix = sel;   
         } 
         else if(source == predefinedSwitch)
         { 
            lpre = sel; larima = !sel;
            p.setEnabled(larima); q.setEnabled(larima); P.setEnabled(larima); Q.setEnabled(larima); d.setEnabled(larima); D.setEnabled(larima);
            if(td){TimeVarModel.setEnabled(lpre);}  airlineModel.setEnabled(lpre);
            TrendModel.setEnabled(lpre);  SeasonalModel.setEnabled(lpre); 
            IrregularModel.setEnabled(lpre);  RandomWalkModel.setEnabled(lpre); 
            resetInputs();               
         }

         else if(source == plotFore)
         {plotForecasts = sel; reg_canvas.setForecasts(sel);}
 
         for(i=0; i <= reg.Ncmpnt; i++)
         {          
           if(source == timePlot[i]) 
           {updateTime(i,sel);}         
         }
        }
      }

      public void updateTime(int i, boolean sel)
      {reg_canvas.setPlots(i,sel);}


      public void setupInterface()
      {


       GroupLayout paramLayout;
 

       //---------- TimeCheckBoxes layout -------------------------------------------
       BevelBorder timeSelBorder = new BevelBorder(BevelBorder.RAISED);
       JLabel tdp = new JLabel("Component Plots   ");
       JPanel timeGrid = new JPanel(); timeGrid.setSize(new Dimension(900, 50));      
        paramLayout = new GroupLayout(timeGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(tdp).addComponent(plotFore).addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
          .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6]) 
          .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[10]));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
           .addComponent(tdp).addComponent(plotFore).addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
           .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(timePlot[5]).addComponent(timePlot[6])
           .addComponent(timePlot[7]).addComponent(timePlot[8]).addComponent(timePlot[9]).addComponent(timePlot[10])));    
        timeGrid.setLayout(paramLayout); 
        timeGrid.setBorder(timeSelBorder);  
        //-----------------------------------------------------------------------
       Box timePane = Box.createVerticalBox();
       

       //---------- SARIMA model layout -------------------------------------------
       //--------------------------------------------------------------------------
       JPanel prePanel = new JPanel();
       JPanel modelXPanel = new JPanel();
       JPanel varPanel = new JPanel();
       new JPanel();
       
 
       //modelXPanel.setLayout(new GridLayout(2, 3, 0, 0));       
       //modelXPanel.add(pX); modelXPanel.add(dX); modelXPanel.add(qX); 
       //modelXPanel.add(PX); modelXPanel.add(DX); modelXPanel.add(QX);

        paramLayout = new GroupLayout(modelXPanel);                
        paramLayout.setAutoCreateGaps(false);
        paramLayout.setAutoCreateContainerGaps(false);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(pX).addComponent(PX))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(dX).addComponent(DX))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(qX).addComponent(QX)));

         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(pX).addComponent(dX).addComponent(qX))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(PX).addComponent(DX).addComponent(QX)));
       modelXPanel.setLayout(paramLayout); //modelXPanel.setSize(new Dimension(20,50));
       

       //varPanel.add(innvText); varPanel.add(innvar); varPanel.add(innvarText); varPanel.add(new JLabel("    ")); 
       //varPanel.add(addModel);
       //varPanel.setLayout(new FlowLayout(FlowLayout.LEFT)); 

        paramLayout = new GroupLayout(varPanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup() 
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)  
             .addComponent(innvText))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(innvar).addComponent(addModel))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
               .addComponent(innvarText)));

         paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(innvText).addComponent(innvar).addComponent(innvarText))
          .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(addModel)));
         varPanel.setLayout(paramLayout); 
       




        /*paramLayout = new GroupLayout(sarimaPanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(arimaSwitch))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(modelXPanel)));

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
              .addComponent(arimaSwitch).addComponent(modelXPanel)));
        sarimaPanel.setLayout(paramLayout);*/

        paramLayout = new GroupLayout(prePanel);                
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(predefinedSwitch))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(TrendModel)
              .addComponent(RandomWalkModel)
              .addComponent(IrregularModel))  
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(SeasonalModel)
              .addComponent(TimeVarModel)
              .addComponent(airlineModel))); 

        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(predefinedSwitch).addComponent(TrendModel).addComponent(SeasonalModel))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(RandomWalkModel).addComponent(TimeVarModel))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)     
              .addComponent(IrregularModel).addComponent(airlineModel)));
        prePanel.setLayout(paramLayout);

       new TitledBorder(new LineBorder(myBlue),"Coefficient and Data Upload");
       JTabbedPane coefs = new JTabbedPane();
       //coefs.setBorder(coeffBorder);       

       JPanel coeffPanel1 = new JPanel();  JPanel coeffPanel2 = new JPanel();  JPanel coeffPanel3= new JPanel(); 

       paramLayout = new GroupLayout(coeffPanel1); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(arbut).addComponent(sarbut))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(arfield) .addComponent(sarfield))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(arfix) .addComponent(sarfix)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(arbut).addComponent(arfield).addComponent(arfix))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(sarbut).addComponent(sarfield).addComponent(sarfix)));
        coeffPanel1.setLayout(paramLayout);

       paramLayout = new GroupLayout(coeffPanel2); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(mabut).addComponent(smabut))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(mafield) .addComponent(smafield))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(mafix) .addComponent(smafix)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(mabut).addComponent(mafield).addComponent(mafix))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(smabut).addComponent(smafield).addComponent(smafix)));
        coeffPanel2.setLayout(paramLayout);

       paramLayout = new GroupLayout(coeffPanel3); paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
       paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(hfileBut).addComponent(reset))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(f_field))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(hfile) .addComponent(allfix)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(hfileBut).addComponent(f_field).addComponent(hfile))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(reset).addComponent(allfix)));
        coeffPanel3.setLayout(paramLayout);
       
        coefs.addTab("AR",coeffPanel1);  coefs.addTab("MA",coeffPanel2); coefs.addTab("Other",coeffPanel3);



    
       JPanel regPanel = new JPanel();
       TitledBorder regBorder = new TitledBorder(new LineBorder(myBlue),"Regression Models");
       regPanel.setBorder(regBorder);     
        
        paramLayout = new GroupLayout(regPanel);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()        
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(tdBox)
              .addComponent(easterBox)
              .addComponent(outlierBox))
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              .addComponent(tdComboBox)
              .addComponent(easterComboBox)
              .addComponent(constantBox))  
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
              //.addComponent(tvreg)
              .addComponent(easterDay)
              .addComponent(transBox)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(tdBox).addComponent(tdComboBox))//.addComponent(tvreg))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(easterBox).addComponent(easterComboBox).addComponent(easterDay))
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
             .addComponent(outlierBox).addComponent(constantBox).addComponent(transBox)));  
        regPanel.setLayout(paramLayout); 


        //--------- Start your engine panel --------------------------

        //JButton addModel,deleteModel,compute,reset; 

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new FlowLayout()); 
        buttonPanel.add(deleteModel); buttonPanel.add(changeModel); 
        buttonPanel.add(compute); 
        
        //-------- Now put everything togehter -----------------------------------------------

        JPanel leftPanel = new JPanel();  //---- Left Panel ------
        paramLayout = new GroupLayout(leftPanel);
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(regPanel).addComponent(coefs)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
            .addComponent(regPanel).addComponent(coefs)) ;
        leftPanel.setLayout(paramLayout);

        JPanel centerPanel = new JPanel();  //---- Left Panel ------
        paramLayout = new GroupLayout(centerPanel);
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(modelXPanel).addComponent(prePanel))); //.addComponent(varPanel)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
            .addComponent(modelXPanel).addComponent(prePanel)); //.addComponent(varPanel)) ;
        centerPanel.setLayout(paramLayout);

       //Box centerPanel1 = Box.createVerticalBox();
       //centerPanel1.add(modelXPanel); centerPanel1.add(prePanel); //centerPanel.add(varPanel);
       TitledBorder sarimaBorder = new TitledBorder(new LineBorder(myBlue),"SARIMA Models");
       centerPanel.setBorder(sarimaBorder);

        Box centerPanel1 = Box.createVerticalBox(); centerPanel1.add(centerPanel); centerPanel1.add(varPanel);

        JPanel rightPanel = new JPanel(); //---- right Panel ------
        paramLayout = new GroupLayout(rightPanel);
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(modelLiszt).addComponent(buttonPanel)));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
            .addComponent(modelLiszt).addComponent(buttonPanel)) ;
        rightPanel.setLayout(paramLayout);        
        

        interfacePane = new JPanel();
        paramLayout = new GroupLayout(interfacePane);
        paramLayout.setAutoCreateGaps(true); paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()   
           .addComponent(leftPanel).addComponent(centerPanel1).addComponent(rightPanel));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(leftPanel).addComponent(centerPanel1).addComponent(rightPanel)));
        interfacePane.setLayout(paramLayout);         
        

        /*-------------------------------------------------------------
            Puts everything together in one final stucke
        --------------------------------------------------------------*/
        timePane.add(reg_canvas);
        timePane.add(timeGrid);

        paramLayout = new GroupLayout(this);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()           
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(timePane)
             .addComponent(interfacePane)));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
           .addComponent(timePane)
           .addComponent(interfacePane));
        this.setLayout(paramLayout);


   } 


  public class ExtensionFilter extends FileFilter {
    private String extensions[];

    private String description;

    public ExtensionFilter(String description, String extension) {
      this(description, new String[] { extension });
    }

    public ExtensionFilter(String description, String extensions[]) {
      this.description = description;
      this.extensions = (String[]) extensions.clone();
    }

    public boolean accept(File file) {
      if (file.isDirectory()) {
        return true;
      }
      int count = extensions.length;
      String path = file.getAbsolutePath();
      for (int i = 0; i < count; i++) {
        String ext = extensions[i];
        if (path.endsWith(ext)
            && (path.charAt(path.length() - ext.length()) == '.')) {
          return true;
        }
      }
      return false;
    }

    public String getDescription() {
      return (description == null ? extensions[0] : description);
    }
  } 


    public void readData(File file)
    {
          
       String strline; Double D;
       String[] tokens; String delims = "[ ]+";
       int n_toks;
       double[] values = new double[500];  double[] real_series;    
       double val = 0;
       int i = 0;
       try{
          
         FileInputStream fin = new FileInputStream(file);
         DataInputStream din = new DataInputStream(fin);
         BufferedReader br = new BufferedReader(new InputStreamReader(din));
 
         while((strline = br.readLine()) != null)
         {

           tokens = strline.split(delims); 
           n_toks = tokens.length; 
           if(n_toks == 0)
           {System.out.println("End of file"); break;}
  
           D = new Double(tokens[n_toks-1]); //take only the value, no dates
           val = D.doubleValue();

           if(i < 500) {values[i] = val; i++;}
           else 
           {System.out.println("Maximum times series length is 500"); break;} 
           
         } 
         din.close();
        }
        catch(FileNotFoundException fe){System.out.println("File not found" + fe);}
        catch(IOException ioe){System.out.println("IO out error..." + ioe);}

        // ------ Update the inputed data--------------------------------------
        n_obs = i; 
        real_series = new double[n_obs];
        for(i=0; i < n_obs; i++)
        {real_series[i] = values[i]; }//System.out.println(real_series[i]);}
        setData(real_series);
    }   


    public void inputX13data(double[] series)
    {
      setData(series);
    }
}
