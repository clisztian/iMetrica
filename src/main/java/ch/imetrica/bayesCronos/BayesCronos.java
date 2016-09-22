package ch.imetrica.bayesCronos;

import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.JCheckBox;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.BorderFactory; 
import java.text.*;


public class BayesCronos extends JPanel
{


   /**
	 * 
	 */
   private static final long serialVersionUID = 1L;
   private JFrame frame;
   JDialog parameterDialog;
   
   JPanel armagarchPanel;
   JPanel svparameterPanel;
   JPanel heavyparameterPanel;

   int n_obs;
   int n_steps;
   int n_sims;
   private int n_factors;
   int n_samps;      //-- number of Bayesian samples
   public int n_rep;
   int p,q,d;        //-- params for ARIMA GARCH
   int paramPlot;
   int n_iters;
 
   boolean f_ready,series_ready,bayes_ready;

   JCheckBox avg_forecast;   

   DecimalFormat df;
   double AICC;
   int model;        //0 = ARMA, 1 = GARCH, 2 = EGARCH, 3 = SV, 4 = FSV
   int methodtype;   //1 = MLE,  2 = BAYESIAN

   double[] tseries;
   double[] residuals;
   double[] fcasts; 
   double[] predictive;
   double[] fmse;

   public JCheckBox[] timePlot;
   public JCheckBox[] factorPlot;
   public JCheckBox[] seriesPlot;
   public JCheckBox[] alphasPlot; 
   public JCheckBox residualPlot; 
   public JCheckBox predictivePlot;
   

   int start, end;

     public Cronos cr;
     public BayesPlot histogram; 
     public CronosPlot cronos_canvas; 
 

    // Variables declaration 
     JLabel arLabel;
     JTextField arText;
     JSlider arimaHistSlider;
     JCheckBox arimaLowerCheck;
     JPanel arimaParamPanel;
     JComboBox<String> arimaPlotCombo; //JComboBox<String>
     JLabel arimaPlotHistLabel;
     
     JCheckBox arimaUpperCheck;
     JButton bayesButton;
     DisabledItemsComboBox chooseModelBox;
     JLabel chooseModelLabel;
     JLabel computeLabel;
     JComboBox<String> ddimBox;
     JLabel ddimLabel;
     JLabel diffLabel;
     JTextField diffText;
     JComboBox<String> factorsBox;
     JLabel factorsLabel;
     JComboBox<String> forecastSampBox;
     JLabel forecastSampLabel;
     JTextField forecastStepsBox;
     JLabel forecastStepsLabel;
     JScrollBar forecastStepsSlider;
     JLabel fphiLabel;
     JTextField fphiText;
     JLabel fsigmaLabel;
     JTextField fsigmaText;
     JLabel gsLabel;
     JTextField gsText;
     JLabel maLabel;
     JTextField maText;
     JButton mleButton;
     JPanel modelDimensionPanel;
     JLabel muLabel;
     JTextField muText;
     JComboBox<String> pdimBox;
     JLabel pdimLabel;
     JLabel phisLabel;
     JTextField phisText;
     JComboBox<String> qdimBox;
     JLabel qdimLabel;
     JLabel sigmaLabel;
     JTextField sigmaText;
     JComboBox<String> sampsCombo;
     JLabel sampsLabel;
     JScrollPane timeScrollPane;
     JScrollBar timeScrollBar;
     JLabel AICCLabel; 
     JTextField AICCText;

     JPanel mcmcPanel = new JPanel();
     JLabel burnLabel = new JLabel();
     @SuppressWarnings("rawtypes")
	 JComboBox burnCombo = new JComboBox();
     JButton parameterButton;
     JPanel modelParamPanel = new JPanel();
     JPanel forecastPanel = new JPanel();
     JPanel computePanel = new JPanel();


    


     //--------- HEAVY Stuff ----------------------
    JLabel alpha2Label, alpha1Label;
    JTextField alpha2Text, alpha1Text,  alpha3Text;
    JLabel alpha3Label, beta1Label, beta2Label, beta3Label;
    JTextField beta3Text, beta1Text, musText;
    JLabel heavyLabel, infoLabel;
    JCheckBox varcheckBox, r2checkBox;
    JLabel w2Label, w1Label;
    JTextField infoText, beta2Text, w1Text, w3Text, w2Text, scoreText, lambda2Text, lambda1Text;
    JLabel fmuLabel, lambda2Label, lambda1Label, w3Label, scoreLabel, fsvLabel, musLabel, arima1Label;

    
         
     public BayesCronos(JFrame fra)
     {
        frame = fra;
        n_obs = 144; n_steps = 12; n_sims = 1;
        setN_factors(1); n_rep = 1;
   
        n_samps = 1000;      //-- number of Bayesian samples
        n_iters = 5000; 
        paramPlot = 0;
        methodtype = 1;
        
        f_ready = false; // ready for forecasts
        series_ready = false; //ready for models
        bayes_ready = false; //ready for bayes plots
        
        AICC = 0.0; model = 0;        
        df = new DecimalFormat("##.##"); 
        tseries = new double[n_obs];
        //System.arraycopy(ts, 0, tseries, 0, n_obs);
        
        p = 1; q = 1; d = 0; 
        start = 0; end = 500;     
    
        //---- Initialize Cronos engine -------------
        cr = new Cronos(n_obs, model, methodtype);
        cr.setData(tseries);
        cr.setNForecastSteps(n_steps);
        cr.setNPredictiveSims(n_steps, n_sims);
        cr.setARMA_Params(p, q, d); 
        cr.setGARCH_Pararms(p, q);
        cr.setSVM();
        cr.setNIterations(n_iters, n_samps);        
        cr.setInitialParamsISV(0.8, 0.1, 1.0);
//        cr.adjustSampleStorage();    
//        cr.computeBayesianISVModel();    
        
        //---- Initialize the cronos time domain ------------------------
        int width = 600; int height = 350;  
        cronos_canvas = new CronosPlot(width+600, height, n_obs, n_rep);
        timeScrollPane = new JScrollPane();
        timeScrollPane.setPreferredSize(new Dimension(width, height)); 
        timeScrollPane.setViewportView(cronos_canvas);
        timeScrollBar = timeScrollPane.getHorizontalScrollBar();     
        
             
        parameterDialog = new JDialog(frame,true);
        parameterButton = new JButton("Parameters");
        parameterButton.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        parameterButton.setText("Parameter Values");
        //----Initialize the histogram plot -----------------------------
        histogram = new BayesPlot();
        histogram.setPreferredSize(new Dimension(450,181));
       
        initControlPanel();
        initParameterControlPanel();
        
        


        readData(new File("tseries2.dat"));
        setupButtons();
        
     }
    
     //--------- Input time series in the form n_rep x n_obs
    
    
     public void setTimeSeries(double[] ts, int nrep, int nobs)
     {
        n_rep = nrep; n_obs = nobs; int i;         
        tseries = new double[n_rep*n_obs]; 
        System.arraycopy(ts, 0, tseries, 0, n_rep*n_obs); 
        cr.setData(tseries);
        cronos_canvas.setData(tseries, nrep, nobs);
        f_ready = false; series_ready = true;

        for(i=0;i<5;i++) 
        { 
          timePlot[i].setSelected(false); timePlot[i].setEnabled(false); 
        }         
     }
    
    

   
    //--------- Data in the form of an ArrayList ---------
    
    @SuppressWarnings("unchecked")
	public void inputMultipleData(ArrayList<double[]> list, boolean[] rep, int nobs)
    {
        
         int i,k,j; int _n_rep = 0; double[] temp; 
         for(i=0;i<rep.length;i++)  {if(rep[i]) {_n_rep++;}}
   
         n_rep = _n_rep; n_obs = nobs;
         tseries = new double[n_rep*n_obs];
        
         //System.out.println("N_REP = " + n_rep + ", N_OBS = " + n_obs);

         _n_rep = 0;
         for(j=0;j<rep.length;j++)
         {
           if(rep[j])
           {
            temp = list.get(j); 
            for(i=0;i<n_obs;i++)
            {
              tseries[_n_rep*n_obs + i] = temp[i];           
            }
            _n_rep++;
           }
         }
      
         System.out.println("n_reps = " + n_rep);
   
         if(n_rep == 1)
         {cr.setData(tseries);}
         else if(n_rep > 1 )
         {cr.setMultiData(tseries, n_rep, n_obs);}
  
         cronos_canvas.setData(tseries, n_rep, n_obs);    
         f_ready = false; 

         for(k=0;k<5;k++)
         {
          seriesPlot[k].setSelected(false); seriesPlot[k].setEnabled(false);
         }
         for(k=0;k<n_rep;k++) {seriesPlot[k].setEnabled(true); seriesPlot[k].setSelected(true);}

         for(k=0;k<3;k++) 
         {
           alphasPlot[k].setSelected(false); alphasPlot[k].setEnabled(false);
           factorPlot[k].setEnabled(false); factorPlot[k].setSelected(false);
         }   

         //if(rep > 1) disable index 1-4
         //chooseModelBox.set
      
         if(n_rep > 1)
         {
           chooseModelBox.removeAllItems();
           chooseModelBox.addItem("ARIMA",true);
           chooseModelBox.addItem("GARCH",true);
           chooseModelBox.addItem("EGARCH",true);
           chooseModelBox.addItem("Stochastic Vol.",true);
           chooseModelBox.addItem("Factor Stochastic Vol.");
           chooseModelBox.addItem("HEAVY");
           chooseModelBox.setSelectedIndex(4);
         }
         else
         {
           chooseModelBox.removeAllItems();
           chooseModelBox.addItem("ARIMA");
           chooseModelBox.addItem("GARCH");
           chooseModelBox.addItem("EGARCH");
           chooseModelBox.addItem("Stochastic Vol.");
           chooseModelBox.addItem("Factor Stochastic Vol.",true);
           chooseModelBox.addItem("HEAVY",true);
           chooseModelBox.setSelectedIndex(0);
         }

    }    
    
    
    
    
    
    
     class MyActionListener implements ActionListener
     {  
      
      @SuppressWarnings("static-access")
	  public void actionPerformed(ActionEvent e)
      {
        int i;
        double[] temp; int[] index; 
        if(e.getSource() == pdimBox)
        {            
          p = pdimBox.getSelectedIndex();
          if(model == 0)
          {cr.setARMA_Params(p, q, d);}          
          else if(model == 1 || model == 2)
          {cr.setGARCH_Pararms(p, q);}
        }
        else if(e.getSource() == qdimBox)
        {            
          q = qdimBox.getSelectedIndex();
          if(model == 0)
          {cr.setARMA_Params(p, q, d);}          
          else if(model == 1 || model == 2)
          {cr.setGARCH_Pararms(p, q);}
        }
        else if(e.getSource() == ddimBox)
        {            
          d = ddimBox.getSelectedIndex();
          if(model == 0)
          {cr.setARMA_Params(p, q, d);}          
        }
        else if(e.getSource() == chooseModelBox)
        {
           model = chooseModelBox.getSelectedIndex();
           setupButtons(); 
        }
        else if(e.getSource() == arimaPlotCombo)
        {
          if(bayes_ready)
          {

            paramPlot = arimaPlotCombo.getSelectedIndex();        
                      
            if(model == 0 || model == 1)
            {
               //System.out.println(n_samps);
               temp = new double[n_samps]; index = new int[n_samps];
               System.arraycopy(cr.bayes_sims[paramPlot], 0, temp, 0, n_samps);
               cr.sortsims(n_samps, temp, index);         
               histogram.createBoxes(temp);
               histogram.colorShadeBoxes(start, end, 1000);
            }
            
            if(model == 3)
            {
               temp = new double[n_samps]; index = new int[n_samps];
               if(paramPlot == 0) {System.arraycopy(cr.phi_sim, 0, temp, 0, n_samps); }
               else if(paramPlot == 1) {System.arraycopy(cr.mu_sim, 0, temp, 0, n_samps); cr.sortsims(n_samps, temp, index); }
               else if(paramPlot == 2) {System.arraycopy(cr.sig_sim, 0, temp, 0, n_samps); cr.sortsims(n_samps, temp, index); }
               histogram.createBoxes(temp);
               histogram.colorShadeBoxes(start, end, 1000);
            }
            
            if(model == 4)
            {

              temp = new double[n_samps]; index = new int[n_samps];             

              for(i=0;i<getN_factors();i++) 
              {
                if(paramPlot == i) {System.arraycopy(cr.fphi_sim, i*n_samps, temp, 0, n_samps);}
                else if(paramPlot == (i+getN_factors())) {System.arraycopy(cr.fsig_sim, i*n_samps, temp, 0, n_samps);}
              }
              
              for(i=0;i<n_rep;i++)
              {
                if(paramPlot == (2*getN_factors() + i)) {System.arraycopy(cr.phi_sim, i*n_samps, temp, 0, n_samps);}
                else if(paramPlot == (2*getN_factors() + n_rep + i)) {System.arraycopy(cr.mu_sim, i*n_samps, temp, 0, n_samps);}
                else if(paramPlot == (2*getN_factors() + 2*n_rep + i)) {System.arraycopy(cr.sig_sim, i*n_samps, temp, 0, n_samps);}
              }

              cr.sortsims(n_samps, temp, index);
              histogram.createBoxes(temp);
              histogram.colorShadeBoxes(start, end, 1000);
            
            }
          }
        }
        else if(e.getSource() == factorsBox)
        {
          if(n_rep > 1)
          {
            if(factorsBox.getSelectedIndex() < n_rep) 
            {
             setN_factors(factorsBox.getSelectedIndex()+1);        
             cr.setNFactors(getN_factors()); 
             cr.setMSVM();
            } 
          }
        }
        else if(e.getSource() == forecastSampBox)
        {
          if(f_ready)
          {
            n_sims = forecastSampBox.getSelectedIndex() + 1;        
            cr.setNPredictiveSims(n_steps,n_sims); 
            computeForecasts();

            for(i=0;i<n_sims;i++) 
            { 
              timePlot[i].setEnabled(true); 
            }
          }
        }
        else if(e.getSource() == sampsCombo)
        {
            int indx = sampsCombo.getSelectedIndex();        
            if(indx == 0)      {n_samps = 1000;}
            else if(indx == 1) {n_samps = 1500;}
            else if(indx == 2) {n_samps = 3000;}
            else if(indx == 3) {n_samps = 5000;}
            else if(indx == 4) {n_samps = 10000;}   
            cr.setNIterations(n_iters, n_samps);
        }
       }
    }





    //---------------------------------------------------------
    // Changes the set of buttons reflection model choice 
    //---------------------------------------------------------
    public void setupButtons()
    {
       String[] options; int i;
 
       bayes_ready = false;
       if(model == 0)        
       {
         pdimBox.setEnabled(true); qdimBox.setEnabled(true); ddimBox.setEnabled(true);
         mleButton.setEnabled(true); bayesButton.setEnabled(true);
         factorsBox.setEnabled(false); gsText.setEnabled(false);
         diffText.setEnabled(true);
         
         maText.setEnabled(true); arText.setEnabled(true); 
         muText.setEnabled(true); sigmaText.setEnabled(true);
           
         phisText.setEnabled(false); fphiText.setEnabled(false); 
         fsigmaText.setEnabled(false);
 
         p = pdimBox.getSelectedIndex(); q = qdimBox.getSelectedIndex();         
         d = ddimBox.getSelectedIndex(); 

         options = new String[p+q+1]; arimaPlotCombo.removeAllItems();
         options[0] = "\u03C3"; arimaPlotCombo.addItem(options[0]);
         for(i=0;i<p;i++) {options[1+i] = "\u03C6" + "_" + (i+1); arimaPlotCombo.addItem(options[1+i]); }  
         for(i=0;i<q;i++) {options[1+p+i] = "\u03B8" + "_" + (i+1); arimaPlotCombo.addItem(options[1+p+i]); } 
                
         //arimaPlotCombo.setModel(new DefaultComboBoxModel(options));
          
       }
       else if(model == 1)        
       {
         pdimBox.setEnabled(true); qdimBox.setEnabled(true); ddimBox.setEnabled(false);
         mleButton.setEnabled(true); bayesButton.setEnabled(true);
         factorsBox.setEnabled(false); gsText.setEnabled(false);
         diffText.setEnabled(false);

         maText.setEnabled(true); arText.setEnabled(true); 
         muText.setEnabled(true); sigmaText.setEnabled(false);
           
         phisText.setEnabled(false); fphiText.setEnabled(false); 
         fsigmaText.setEnabled(false); 
 
         p = pdimBox.getSelectedIndex(); q = qdimBox.getSelectedIndex();      

         arimaPlotCombo.removeAllItems();
         options = new String[p+q+1+1];
         options[0] = "\u03BC"; arimaPlotCombo.addItem(options[0]);
         for(i=0;i<=p;i++) {options[1+i] = "\u03C6" + "_" + i; arimaPlotCombo.addItem(options[1+i]);} 
         for(i=0;i<q;i++) {options[2+p+i] = "\u03B8" + "_" + (i+1); arimaPlotCombo.addItem(options[2+p+i]);} 
         

         //arimaPlotCombo.setModel(new DefaultComboBoxModel(options));


       }
       if(model == 2)        
       {
         pdimBox.setEnabled(true); qdimBox.setEnabled(true); ddimBox.setEnabled(false);
         mleButton.setEnabled(true); bayesButton.setEnabled(false);
         factorsBox.setEnabled(false); gsText.setEnabled(true);
         diffText.setEnabled(false);

         maText.setEnabled(true); arText.setEnabled(true); 
         muText.setEnabled(true); sigmaText.setEnabled(false);
           
         phisText.setEnabled(false); fphiText.setEnabled(false); 
         fsigmaText.setEnabled(false); 

         options = new String[0]; //options[0] = "-"; 
         arimaPlotCombo.removeAllItems();
         //arimaPlotCombo.addItem(options);
         //arimaPlotCombo.setModel(new DefaultComboBoxModel(options));

       }
       if(model == 3)        
       {
         pdimBox.setEnabled(false); qdimBox.setEnabled(false); ddimBox.setEnabled(false);
         mleButton.setEnabled(true); bayesButton.setEnabled(true);
         factorsBox.setEnabled(false); gsText.setEnabled(false);
         diffText.setEnabled(false);

         maText.setEnabled(false); arText.setEnabled(false); 
         muText.setEnabled(true); sigmaText.setEnabled(true);
           
         phisText.setEnabled(true); fphiText.setEnabled(false); 
         fsigmaText.setEnabled(false); 
 
         options = new String[3];
         options[0] = "\u03C6"; options[1] = "\u03BC"; options[2] = "\u03C3"; 
         arimaPlotCombo.removeAllItems();
         arimaPlotCombo.addItem(options[0]);
         arimaPlotCombo.addItem(options[1]);
         arimaPlotCombo.addItem(options[2]);
         //arimaPlotCombo.setModel(new DefaultComboBoxModel(options));

       }
       if(model == 4)        
       {
         pdimBox.setEnabled(false); qdimBox.setEnabled(false); ddimBox.setEnabled(false);
         mleButton.setEnabled(false); bayesButton.setEnabled(true);
         factorsBox.setEnabled(true); gsText.setEnabled(false);
         diffText.setEnabled(false);

         maText.setEnabled(false); arText.setEnabled(false); 
         muText.setEnabled(true); sigmaText.setEnabled(true);
           
         phisText.setEnabled(true); fphiText.setEnabled(true); 
         fsigmaText.setEnabled(true); 
 
         arimaPlotCombo.removeAllItems();

         options = new String[getN_factors()*2 + n_rep*3];
      
         for(i=0;i<getN_factors();i++) {options[i] = "f_\u03C6_" + (i+1); arimaPlotCombo.addItem(options[i]);}
         for(i=0;i<getN_factors();i++) {options[getN_factors() + i] = "f_\u03C3_" + (i+1); arimaPlotCombo.addItem(options[getN_factors()+i]);}
          
         for(i=0;i<n_rep;i++) {options[2*getN_factors() + i] = "\u03C6_" + (i+1); arimaPlotCombo.addItem(options[2*getN_factors()+i]); }
         for(i=0;i<n_rep;i++) {options[2*getN_factors() + n_rep + i] = "\u03C5_" + (i+1); arimaPlotCombo.addItem(options[2*getN_factors() + n_rep + i]);}
         for(i=0;i<n_rep;i++) {options[2*getN_factors() + 2*n_rep + i] = "\u03C3_" + (i+1); arimaPlotCombo.addItem(options[2*getN_factors() + 2*n_rep + i]);}
        
         //arimaPlotCombo.setModel(new DefaultComboBoxModel(options)); 
       }
       if(model == 5)        
       {
         pdimBox.setEnabled(false); qdimBox.setEnabled(false); ddimBox.setEnabled(false);
         mleButton.setEnabled(true); bayesButton.setEnabled(false);
         factorsBox.setEnabled(false); gsText.setEnabled(false);
         diffText.setEnabled(false);

         maText.setEnabled(false); arText.setEnabled(false); 
         muText.setEnabled(false); sigmaText.setEnabled(false);
           
         phisText.setEnabled(false); fphiText.setEnabled(false); 
         fsigmaText.setEnabled(false); 
 
         arimaPlotCombo.removeAllItems();
       }




       bayes_ready = true;       
    }


    public void computeModel(int method)
    {
      int k,i;   
 
      if(series_ready)
      {
        arimaHistSlider.setEnabled(true);
        for(k=0;k<5;k++)
        {          
          seriesPlot[k].setSelected(false); seriesPlot[k].setEnabled(false);
        }
        for(k=0;k<3;k++) 
        {
          alphasPlot[k].setSelected(false); alphasPlot[k].setEnabled(false);
          factorPlot[k].setEnabled(false); factorPlot[k].setSelected(false);
        } 

        alphasPlot[0].setEnabled(true); seriesPlot[0].setEnabled(true);


        if(model == 0)
        {
          cr.computeARIMAModel(method);
        }
        else if(model == 1)
        {
          cr.computeGARCHModel(method);
        }
        else if(model == 2)
        {
          cr.computeEGARCHModel();
        }
        else if(model == 3)
        {
          if(method == 1)
          {cr.computeSVModel(); cronos_canvas.setplot_alpha(0,false); }
          else
          {
            cr.adjustSampleStorage();    
            cr.computeBayesianISVModel();
            cr.computeBayesianAvgISV(0, n_samps); 
            cronos_canvas.setAlpha(cr.alpha, 1, n_obs);
            cronos_canvas.setplot_alpha(0,true);         
          }      
        }
        else if(model == 4) 
        {

           setN_factors(factorsBox.getSelectedIndex()+1);        
           cr.setNFactors(getN_factors()); 
           cr.setMSVM();          

           double[] temp2 = new double[getN_factors()*n_obs];
           double[] temp = new double[getN_factors()*n_obs];
           cr.setDefaultParameterValues();
           cr.computeFSVModel(); 
           //cr.computes the averages already
 
        
           for(k=0;k<getN_factors();k++)
           {
            for(i=0;i<n_obs;i++)
            {
             temp[k*n_obs + i] = cr.falpha_avg[k][i]; 
             temp2[k*n_obs + i] = cr.factors_avg[k][i];
            }
            factorPlot[k].setEnabled(true); factorPlot[k].setSelected(false); 
            alphasPlot[k].setSelected(false); alphasPlot[k].setEnabled(true);
           }

           cronos_canvas.setAlpha(temp, getN_factors(), n_obs);
           cronos_canvas.setFactors(temp2,getN_factors(),n_obs);

           for(k=0;k<n_rep;k++)
           {seriesPlot[k].setSelected(true); seriesPlot[k].setEnabled(true);}
                     
           //---- residuals and predictive off------------
           residualPlot.setSelected(false); residualPlot.setEnabled(false);
           predictivePlot.setSelected(false); predictivePlot.setEnabled(false);

        } 
        if(model == 5) // -- HEAVY Model (factor model with n_rep = 2/3 n_factors = n_rep
        {
  
           int retrack = 0;
           setN_factors(n_rep); 
           cr.setNFactors(getN_factors());
 
                     
           double[] residual = new double[n_obs];
           double[] temp2 = new double[getN_factors()*n_obs];
           double[] temp = new double[getN_factors()*n_obs];

           n_samps = forecastSampBox.getSelectedIndex() + 1;

           //System.out.println("n_samps = " + n_samps);
           cr.setNForecastSteps(n_steps);
           cr.setNPredictiveSims(n_steps, n_samps);
       

           //--- get any initial parameters---------------
           getHeavyParametersText();

           if(varcheckBox.isSelected()) {retrack = 1;}
           cr.setTrackReparameter(retrack);
           cr.computeHeavyModel();
 
           for(k=0;k<getN_factors();k++)
           {
            for(i=0;i<n_obs;i++)
            {
             temp[k*n_obs + i] = cr.falpha_avg[k][i]; 
             temp2[k*n_obs + i] = cr.factors_avg[k][i];
            }
            factorPlot[k].setEnabled(true); factorPlot[k].setSelected(false); 
            alphasPlot[k].setSelected(false); alphasPlot[k].setEnabled(true);
           }

           for(i=0;i<n_obs;i++) {residual[i] = cr.heavy.llht[i];}

           cronos_canvas.setAlpha(temp, getN_factors(), n_obs);
           cronos_canvas.setFactors(temp2,getN_factors(),n_obs);
           cronos_canvas.setResidual(residual, 1, n_obs);
           cronos_canvas.setForecasts(cr.predictive, n_steps, n_sims);

           cronos_canvas.setPredictive(cr.fcasts, cr.fcasts.length);
           cronos_canvas.plotPredictive(false);
           predictivePlot.setEnabled(true);

 
           for(k=0;k<n_rep;k++)
           {seriesPlot[k].setSelected(true); seriesPlot[k].setEnabled(true);}

            for(i=0;i<n_sims;i++) 
            { 
              timePlot[i].setEnabled(true); 
            }
                
           //---- residuals and predictive off------------
           residualPlot.setSelected(false); residualPlot.setEnabled(true);
           predictivePlot.setSelected(false); predictivePlot.setEnabled(true);

           arimaHistSlider.setEnabled(false);
           arimaPlotCombo.setEnabled(false);
     
           setHeavyParametersText();
        }
        //----- plot predictive volatility/MSE-------------- 
        if(model < 4)
        {
         cronos_canvas.setPredictive(cr.fcasts, cr.fcasts.length);
         cronos_canvas.plotPredictive(false);
         predictivePlot.setEnabled(true);

         cronos_canvas.setResidual(cr.residuals,1,n_obs);
         residualPlot.setEnabled(true);

        }

        setText(); f_ready = true;
      }
    }
       
    //--------------------------------------------------
    //  Assumes parameters have been chosen 
    //--------------------------------------------------   
       
    public void computeForecasts()
    {
      if(f_ready)
      {
       if(n_steps > 0 && model != 5)
       {
        cr.setNForecastSteps(n_steps);
        cr.setNPredictiveSims(n_steps, n_samps);

        if(model == 0)
        {cr.forecastARIMA();}
        else if(model == 1)
        {cr.forecastGARCH();}
        else if(model == 2)
        {cr.forecastEGARCH();}
        else if(model == 3)
        {cr.computeBayesianForecastISVModel();}

        if(model < 4)
        {
          cronos_canvas.setForecasts(cr.predictive, n_steps, n_sims);
          cronos_canvas.setPredictive(cr.fcasts, cr.fcasts.length);
          cronos_canvas.plotPredictive(true);
        }        
       }
      }
    }

    public void setAlpha()
    {
      int i,k;
      if(model == 3)
      {
        cronos_canvas.setAlpha(cr.alpha, 1, n_obs); 
        alphasPlot[0].setEnabled(true); alphasPlot[0].setSelected(true);

      }
      else if(model == 4)
      {
          double[] temp = new double[getN_factors()*n_obs];
          double[] temp2 = new double[getN_factors()*n_obs];

          for(k=0;k<getN_factors();k++)
          {
           for(i=0;i<n_obs;i++)
           {
            temp[k*n_obs + i] = cr.falpha_avg[k][i]; 
            temp2[k*n_obs + i] = cr.factors_avg[k][i];
           }
          }

          cronos_canvas.setAlpha(temp, getN_factors(), n_obs);
          cronos_canvas.setFactors(temp2, getN_factors(), n_obs);
          cronos_canvas.go();
          //for(k=0;k<n_factors;k++)
          //{alphasPlot[k].setEnabled(true); alphasPlot[k].setSelected(true);}          

      }
    }

    public void changeIntervalSample()
    {
       if(model == 0)
       {cr.computeBayesianAvgARIMA(start, end);}
       else if(model == 1)
       {cr.computeBayesianAvgGARCH(start, end);}
       else if(model == 3)
       {cr.computeBayesianAvgISV(start, end);}
       else if(model == 4) 
       {cr.computeBayesianAvgFISV(start, end);}
  
       //---- set Text with new parameters -------
       setText();
       //---- compute new forecasts --------------
       computeForecasts();
       //---- set new alpha ----------------------
       setAlpha(); 
    }    

    public void highlight_plot(int sel)
    {
      if(sel == 0) {sel = -1;}
      cronos_canvas.changeHighlight(sel);
    }


    public void setText()
    {
      String ar_string = ""; String ma_string = ""; String gs_string = ""; int i;

      if(model == 0) //ARIMA
      {
        ar_string = ""; ma_string = "";
        if(cr.ar_params.length > 0) ar_string = ""+df.format(cr.ar_params[0]); 
        if(cr.ma_params.length > 0) ma_string = ""+df.format(cr.ma_params[0]); 
        
        for(i=1;i<cr.ar_params.length;i++)
        {ar_string = ar_string + " " + df.format(cr.ar_params[i]);} 
  
        for(i=1;i<cr.ma_params.length;i++)
        {ma_string = ma_string + " " + df.format(cr.ma_params[i]);}               
 
     
        arText.setText(ar_string); 
        maText.setText(ma_string);
        muText.setText(df.format(cr.mu));
        sigmaText.setText(df.format(cr.innvar));
        
      }
      else if(model == 1) //Garch
      {
        ar_string = ""; ma_string = "";
        if(cr.garch_ar.length > 0) ar_string = ""+df.format(cr.garch_ar[0]); 
        if(cr.garch_ma.length > 0) ma_string = ""+df.format(cr.garch_ma[0]); 
        
        for(i=1;i<cr.garch_ar.length;i++)
        {ar_string = ar_string + " " + df.format(cr.garch_ar[i]);} 
  
        for(i=1;i<cr.garch_ma.length;i++)
        {ma_string = ma_string + " " + df.format(cr.garch_ma[i]);}               
 
        arText.setText(ar_string); 
        maText.setText(ma_string);
        muText.setText(df.format(cr.mu));
      }
      else if(model == 2) //Garch
      {
        ar_string = ""; ma_string = ""; gs_string = "";
        if(cr.garch_ar.length > 0) ar_string = ""+df.format(cr.garch_ar[0]); 
        if(cr.garch_ma.length > 0) ma_string = ""+df.format(cr.garch_ma[0]); 
        if(cr.egarch_gs.length > 0) ma_string = ""+df.format(cr.egarch_gs[0]);         

        for(i=1;i<cr.garch_ar.length;i++)
        {ar_string = ar_string + " " + df.format(cr.garch_ar[i]);} 
  
        for(i=1;i<cr.garch_ma.length;i++)
        {ma_string = ma_string + " " + df.format(cr.garch_ma[i]);}               
 
        for(i=1;i<cr.egarch_gs.length;i++)
        {gs_string = gs_string + " " + df.format(cr.egarch_gs[i]);}      

        arText.setText(ar_string); 
        maText.setText(ma_string);
        gsText.setText(gs_string);
        muText.setText(df.format(cr.mu));
      }
      else if(model == 3) //SVM
      {
        sigmaText.setText(df.format(cr.mux));
        phisText.setText(df.format(cr.phi));
        muText.setText(df.format(cr.mean));
      }
      else if(model == 4) //FSVM
      {

        if(cr.fphi_avg.length > 0) ar_string = ""+df.format(cr.fphi_avg[0]); 
        if(cr.fsig_avg.length > 0) ma_string = ""+df.format(cr.fsig_avg[0]); 
        if(cr.wphi_avg.length > 0) gs_string = ""+df.format(cr.wphi_avg[0]);         

        for(i=1;i<cr.fphi_avg.length;i++)
        {ar_string = ar_string + " " + df.format(cr.fphi_avg[i]);} 
  
        for(i=1;i<cr.fsig_avg.length;i++)
        {ma_string = ma_string + " " + df.format(cr.fsig_avg[i]);}    
         
        for(i=1;i<cr.wphi_avg.length;i++)
        {gs_string = gs_string + " " + df.format(cr.wphi_avg[i]);}

        fphiText.setText(ar_string);
        fsigmaText.setText(ma_string);
        phisText.setText(gs_string);
      }
      AICCText.setText(df.format(cr.AICC));

    }

    public void plotTimeSeries()
    {
      int k,i;
      double temp[] = new double[getN_factors()*n_obs];
      cronos_canvas.setForecasts(cr.predictive, n_steps, n_sims);
      if(model == 0 || model == 1)
      {
        cronos_canvas.setResidual(cr.residuals, n_rep, n_obs);
      }
      if(model == 3)
      {
        if(methodtype == 0)
        {cronos_canvas.setAlpha(cr.alpha, n_rep, n_obs);}
      }
      if(model == 4)
      {        
        for(k=0;k<getN_factors();k++)
        {
         for(i=0;i<n_obs;i++)
         {temp[k*n_obs + i] = cr.alpha_sim[n_obs*n_samps*k + n_obs*(n_samps/2) + i];}
        }
        cronos_canvas.setAlpha(temp, getN_factors(), n_obs);
      }
      cronos_canvas.go();
   
    }
    

    void getHeavyParametersText()
    {
      double val;
  
      if(!w1Text.getText().equals(""))
      {
       val = Double.valueOf(w1Text.getText()).doubleValue();
       if(val > 0.0 && val < 10.0) {cr.w[0] = val;}     
      }

      if(!w2Text.getText().equals(""))
      {
       val = Double.valueOf(w2Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.w[1] = val;}     
      } 

      if(!w3Text.getText().equals(""))
      {
       val = Double.valueOf(w3Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.w[2] = val;}     
      } 

      if(!alpha1Text.getText().equals(""))
      {
       val = Double.valueOf(alpha1Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_alpha[0] = val;}     
      }
 
      if(!alpha2Text.getText().equals(""))
      {
       val = Double.valueOf(alpha2Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_alpha[1] = val;}     
      }

      if(!alpha3Text.getText().equals(""))
      {
       val = Double.valueOf(alpha3Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_alpha[2] = val;}     
      }


      if(!beta1Text.getText().equals(""))
      {
       val = Double.valueOf(beta1Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_beta[0] = val;}     
      }
 
      if(!beta2Text.getText().equals(""))
      {
       val = Double.valueOf(beta2Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_beta[1] = val;}     
      }

      if(!beta3Text.getText().equals(""))
      {
       val = Double.valueOf(beta3Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_beta[2] = val;}     
      }

      if(!lambda1Text.getText().equals(""))
      {
       val = Double.valueOf(lambda1Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_lambda[0] = val;}     
      }
 
      if(!lambda2Text.getText().equals(""))
      {
       val = Double.valueOf(lambda2Text.getText()).doubleValue();
       if(val > 0.0 && val < 1.0) {cr.h_lambda[1] = val;}     
      }

    }

    void setHeavyParametersText()
    {

      w1Text.setText(""+df.format(cr.heavy.W_mean[0]));
      w2Text.setText(""+df.format(cr.heavy.W_mean[1]));
     
      alpha1Text.setText(""+df.format(cr.heavy.alphas[0]));
      alpha2Text.setText(""+df.format(cr.heavy.alphas[1]));

      beta1Text.setText(""+df.format(cr.heavy.betas[0]));
      beta2Text.setText(""+df.format(cr.heavy.betas[1]));

      if(cr.heavy.rt == 1) 
      {lambda1Text.setText(""+df.format(cr.heavy.lambdas[0]));}     

      if(n_rep == 3)
      {
        w3Text.setText(""+df.format(cr.heavy.W_mean[2]));
        alpha3Text.setText(""+df.format(cr.heavy.alphas[2]));
        beta3Text.setText(""+df.format(cr.heavy.betas[2]));

        if(cr.heavy.rt == 1) 
        {lambda2Text.setText(""+df.format(cr.heavy.lambdas[1]));}
      }
    }    


    public ArrayList<double[]> getPlottedData()  
    {
       int i; int k;
       int count = 0;
       ArrayList<double[]> list = new ArrayList<double[]>();


       for(k=0;k<n_rep;k++)
       {
         if(seriesPlot[k].isSelected())
         { 
           double[] temp = new double[n_obs];
           for(i=0;i<n_obs;i++)
           {temp[i] = tseries[k*n_obs + i];}
           list.add(count,temp); count++;
         }
       }


     if(model == 4 || model == 5)
     {    
       for(k=0;k<getN_factors();k++)
       {
         if(factorPlot[k].isSelected())
         {
           double[] temp = new double[n_obs];
           for(i=0;i<n_obs;i++)
           {temp[i] = cr.factors_avg[k][i];}
           list.add(count,temp); count++;

         }
       }

       for(k=0;k<getN_factors();k++)
       {
         if(alphasPlot[k].isSelected())
         {
           double[] temp = new double[n_obs];
           for(i=0;i<n_obs;i++)
           {temp[i] = cr.falpha_avg[k][i];}
           list.add(count,temp); count++;
         }
       }
     }
     else if(model == 3)
     {
       if(alphasPlot[0].isSelected())
       {
           double[] temp = new double[n_obs];
           for(i=0;i<n_obs;i++)
           {temp[i] = cr.alpha[i];}
           list.add(count,temp); count++;
       }

       if(predictivePlot.isSelected())
       {
        double[] temp = new double[n_obs];
        for(i=0;i<n_obs;i++)
        {temp[i] = cr.fcasts[i];}
        list.add(count,temp); count++;
       }
     }
     else
     {

       if(residualPlot.isSelected())
       {
        double[] temp = new double[n_obs];
        for(i=0;i<n_obs;i++)
        {temp[i] = cr.residuals[i];}
        list.add(count,temp); count++;
       }       

       if(predictivePlot.isSelected())
       {
        double[] temp = new double[n_obs];
        for(i=0;i<n_obs;i++)
        {temp[i] = cr.fcasts[i];}
        list.add(count,temp); count++;
       }
     }

     return list;
    }

    public void print_to_file()
    {cronos_canvas.turnOnPrinter();}

      

  







    public void initParameterControlPanel()
    {


      JLabel sigmaLabel2;
      JLabel sigmaLabel3;




        w1Text = new JTextField();
        w1Label = new JLabel();
        heavyLabel = new JLabel();
        w2Text = new JTextField();
        w2Label = new JLabel();
        w3Label = new JLabel();
        w3Text = new JTextField();
        alpha3Label = new JLabel();
        alpha3Text = new JTextField();
        alpha1Label = new JLabel();
        alpha1Text = new JTextField();
        alpha2Label = new JLabel();
        alpha2Text = new JTextField();
        beta1Text = new JTextField();
        beta3Text = new JTextField();
        beta2Text = new JTextField();
        beta2Label = new JLabel();
        beta1Label = new JLabel();
        beta3Label = new JLabel();
        lambda1Label = new JLabel();
        lambda2Label = new JLabel();
        lambda2Text = new JTextField();
        lambda1Text = new JTextField();
        r2checkBox = new JCheckBox();
        varcheckBox = new JCheckBox();
        scoreLabel = new JLabel();
        scoreText = new JTextField();
        infoText = new JTextField();
        infoLabel = new JLabel();
        arima1Label = new JLabel(); 
 

        musLabel = new javax.swing.JLabel();
        sigmaText = new javax.swing.JTextField();
        fphiText = new javax.swing.JTextField();
        musText = new javax.swing.JTextField();
        fsvLabel = new javax.swing.JLabel();
        fmuLabel = new javax.swing.JLabel();
        phisLabel = new javax.swing.JLabel();
        phisText = new javax.swing.JTextField();
        sigmaLabel2 = new javax.swing.JLabel();
        fsigmaText = new javax.swing.JTextField();
        sigmaLabel3 = new javax.swing.JLabel();


      sigmaLabel2.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
      sigmaLabel2.setText("\"\\u03C3:\"");

      sigmaLabel3.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
      sigmaLabel3.setText("\"f_\\u03C3:\"");


        arima1Label.setText("ARIMA/(E)GARCH Model - Estimated Parameters");


        musLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        musLabel.setText("\"\\u03BC:\"");

        sigmaText.setBackground(new java.awt.Color(1, 1, 1));
        sigmaText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sigmaText.setForeground(new java.awt.Color(1, 249, 11));

        fphiText.setBackground(new java.awt.Color(1, 1, 1));
        fphiText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        fphiText.setForeground(new java.awt.Color(1, 249, 11));

        musText.setBackground(new java.awt.Color(1, 1, 1));
        musText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        musText.setForeground(new java.awt.Color(1, 249, 11));

        fsvLabel.setText("Factor /Stochastic Volatility Model - Estimated Parameters");

        fmuLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        fmuLabel.setText("\"f_\\u03C6:\"");

        phisLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        phisLabel.setText("\"\\u03C6:\"");

        phisText.setBackground(new java.awt.Color(1, 1, 1));
        phisText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        phisText.setForeground(new java.awt.Color(1, 249, 11));



        fsigmaText.setBackground(new java.awt.Color(1, 1, 1));
        fsigmaText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        fsigmaText.setForeground(new java.awt.Color(1, 249, 11));


        w1Text.setBackground(new java.awt.Color(1, 1, 1));
        w1Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w1Text.setForeground(new java.awt.Color(1, 249, 11));

        w1Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w1Label.setText("\"\\u03C9\"");

        heavyLabel.setText("High-Frequency Volatility Model - Estimated Parameters");

        w2Text.setBackground(new java.awt.Color(1, 1, 1));
        w2Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w2Text.setForeground(new java.awt.Color(1, 249, 11));

        w2Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w2Label.setText("\"\\u03C9\"_RM");

        w3Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w3Label.setText("\"\\u03C9\"_SM");

        w3Text.setBackground(new java.awt.Color(1, 1, 1));
        w3Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        w3Text.setForeground(new java.awt.Color(1, 249, 11));

        alpha3Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha3Label.setText("\"\\u03B1\"_SM");

        alpha3Text.setBackground(new java.awt.Color(1, 1, 1));
        alpha3Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha3Text.setForeground(new java.awt.Color(1, 249, 11));

        alpha1Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha1Label.setText("\"\\u03B1\"");

        alpha1Text.setBackground(new java.awt.Color(1, 1, 1));
        alpha1Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha1Text.setForeground(new java.awt.Color(1, 249, 11));

        alpha2Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha2Label.setText("\"\\u03B1\"_RM");

        alpha2Text.setBackground(new java.awt.Color(1, 1, 1));
        alpha2Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        alpha2Text.setForeground(new java.awt.Color(1, 249, 11));

        beta1Text.setBackground(new java.awt.Color(1, 1, 1));
        beta1Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta1Text.setForeground(new java.awt.Color(1, 249, 11));

        beta3Text.setBackground(new java.awt.Color(1, 1, 1));
        beta3Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta3Text.setForeground(new java.awt.Color(1, 249, 11));

        beta2Text.setBackground(new java.awt.Color(1, 1, 1));
        beta2Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta2Text.setForeground(new java.awt.Color(1, 249, 11));

        beta2Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta2Label.setText("\"\\u03D0\"_RM");

        beta1Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta1Label.setText("\"\\u03D0\"");

        beta3Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        beta3Label.setText("\"\\u03D0\"_SM");

        lambda1Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lambda1Label.setText("\"\\u03BB\"");

        lambda2Label.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lambda2Label.setText("\"\\u03BB\"_RM");

        lambda2Text.setBackground(new java.awt.Color(1, 1, 1));
        lambda2Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lambda2Text.setForeground(new java.awt.Color(1, 249, 11));

        lambda1Text.setBackground(new java.awt.Color(1, 1, 1));
        lambda1Text.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        lambda1Text.setForeground(new java.awt.Color(1, 249, 11));

        r2checkBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        r2checkBox.setText("Include r^2_t in model estimation");
        r2checkBox.addItemListener(new ItemListener() {
          public void itemStateChanged(ItemEvent e)
          {

             boolean sel; //computeFilter = true;
             if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
             else{sel = true;}
             if(sel)
             {lambda1Text.setText(""+0.2);}
             else
             {lambda1Text.setText(""+0.0);}
          }
         });
 
        varcheckBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        varcheckBox.setText("Variance Targeting");

        scoreLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        scoreLabel.setText("Score");

        scoreText.setBackground(new java.awt.Color(1, 1, 1));
        scoreText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N

        infoText.setBackground(new java.awt.Color(1, 1, 1));
        infoText.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N

        infoLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        infoLabel.setText("Information Matrix");
 
        fsvLabel = new JLabel();
        fsvLabel.setText("Factor /Stochastic Volatility Model - Estimated Parameters");



       //---- ARIMA/GARCH PANEL------
       armagarchPanel = new JPanel(); 
 
        GroupLayout jPanel1Layout = new GroupLayout(modelDimensionPanel);
        modelDimensionPanel.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(pdimLabel)
                .addGap(18, 18, 18)
                .addComponent(pdimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(qdimLabel)
                .addGap(18, 18, 18)
                .addComponent(qdimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ddimLabel)
                .addGap(18, 18, 18)
                .addComponent(ddimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(90, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(ddimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(ddimLabel))
                    .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(pdimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(pdimLabel)
                        .addComponent(qdimBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(qdimLabel)))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );


        GroupLayout layout = new GroupLayout(armagarchPanel);
        armagarchPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(arima1Label))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(13, 13, 13)
                                .addComponent(arLabel))
                            .addGroup(layout.createSequentialGroup()
                                .addContainerGap()
                                .addComponent(maLabel))
                            .addGroup(layout.createSequentialGroup()
                                .addContainerGap()
                                .addComponent(gsLabel)))
                        .addGap(31, 31, 31)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                                    .addComponent(arText, GroupLayout.PREFERRED_SIZE, 255, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(maText, GroupLayout.PREFERRED_SIZE, 255, GroupLayout.PREFERRED_SIZE)
                                    .addComponent(gsText, GroupLayout.PREFERRED_SIZE, 255, GroupLayout.PREFERRED_SIZE))
                                .addGap(18, 18, 18)
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(muLabel)
                                        .addGap(18, 18, 18)
                                        .addComponent(muText, GroupLayout.DEFAULT_SIZE, 50, Short.MAX_VALUE))
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                            .addComponent(sigmaLabel)
                                            .addComponent(diffLabel))
                                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                                            .addComponent(sigmaText, GroupLayout.DEFAULT_SIZE, 52, Short.MAX_VALUE)
                                            .addComponent(diffText)))))
                            .addComponent(modelDimensionPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(141, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(arima1Label, GroupLayout.PREFERRED_SIZE, 27, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(arText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(arLabel)
                    .addComponent(sigmaLabel)
                    .addComponent(sigmaText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(diffText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(maLabel)
                    .addComponent(diffLabel)
                    .addComponent(maText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(muText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(gsLabel)
                    .addComponent(muLabel)
                    .addComponent(gsText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(modelDimensionPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );



       //--------SV Panel ------------------------------------------------------------

       svparameterPanel = new JPanel();
       layout = new GroupLayout(svparameterPanel);
        svparameterPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(13, 13, 13)
                                .addComponent(phisLabel))
                            .addGroup(layout.createSequentialGroup()
                                .addContainerGap()
                                .addComponent(fmuLabel)))
                        .addGap(31, 31, 31)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING, false)
                            .addComponent(phisText, GroupLayout.Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 207, Short.MAX_VALUE)
                            .addComponent(fphiText, GroupLayout.Alignment.LEADING)
                            .addComponent(musText))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(sigmaLabel2)
                            .addComponent(sigmaLabel3))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.TRAILING, false)
                            .addComponent(fsigmaText, GroupLayout.Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 164, Short.MAX_VALUE)
                            .addComponent(sigmaText, GroupLayout.Alignment.LEADING)))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(fsvLabel)
                            .addComponent(musLabel))))
                .addContainerGap(89, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(fsvLabel, GroupLayout.PREFERRED_SIZE, 27, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(phisText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(phisLabel)
                    .addComponent(sigmaLabel2)
                    .addComponent(sigmaText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(sigmaLabel3)
                    .addComponent(fsigmaText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(fphiText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(fmuLabel))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(musLabel)
                    .addComponent(musText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(92, Short.MAX_VALUE))
        );
    




        //-------------------- HEAVY parameter panel ------------------

        heavyparameterPanel = new JPanel();

        layout = new GroupLayout(heavyparameterPanel);
        heavyparameterPanel.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(13, 13, 13)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(beta1Label)
                            .addComponent(lambda1Label)
                            .addComponent(w1Label)
                            .addComponent(alpha1Label))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(w1Text, GroupLayout.PREFERRED_SIZE, 43, GroupLayout.PREFERRED_SIZE)
                            .addComponent(lambda1Text, GroupLayout.PREFERRED_SIZE, 43, GroupLayout.PREFERRED_SIZE)
                            .addComponent(alpha1Text, GroupLayout.PREFERRED_SIZE, 43, GroupLayout.PREFERRED_SIZE)
                            .addComponent(beta1Text, GroupLayout.PREFERRED_SIZE, 43, GroupLayout.PREFERRED_SIZE))
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(lambda2Label))
                            .addGroup(layout.createSequentialGroup()
                                .addGap(43, 43, 43)
                                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                                    .addComponent(beta2Label)
                                    .addComponent(w2Label)
                                    .addComponent(alpha2Label))))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(beta2Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)
                            .addComponent(lambda2Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)
                            .addComponent(w2Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)
                            .addComponent(alpha2Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE))
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(28, 28, 28)
                                .addComponent(w3Label))
                            .addGroup(layout.createSequentialGroup()
                                .addGap(30, 30, 30)
                                .addComponent(alpha3Label))
                            .addGroup(layout.createSequentialGroup()
                                .addGap(32, 32, 32)
                                .addComponent(beta3Label)))
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING)
                            .addComponent(w3Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)
                            .addComponent(alpha3Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)
                            .addComponent(beta3Text, GroupLayout.PREFERRED_SIZE, 39, GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(264, 264, 264)
                        .addComponent(varcheckBox))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(heavyLabel))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(r2checkBox)))
                .addContainerGap(141, Short.MAX_VALUE))
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(scoreLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(scoreText, GroupLayout.PREFERRED_SIZE, 244, GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(infoLabel)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(infoText, GroupLayout.PREFERRED_SIZE, 244, GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(heavyLabel, GroupLayout.PREFERRED_SIZE, 27, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(w1Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(w1Label)
                    .addComponent(w2Label)
                    .addComponent(w2Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(w3Label)
                    .addComponent(w3Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(alpha1Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(alpha1Label)
                    .addComponent(alpha2Label)
                    .addComponent(alpha2Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(alpha3Label)
                    .addComponent(alpha3Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(beta1Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(beta1Label)
                    .addComponent(beta2Label)
                    .addComponent(beta2Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(beta3Label)
                    .addComponent(beta3Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(lambda1Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(lambda1Label)
                    .addComponent(lambda2Label)
                    .addComponent(lambda2Text, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(r2checkBox)
                    .addComponent(varcheckBox))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(scoreLabel)
                    .addComponent(scoreText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(infoLabel)
                    .addComponent(infoText, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(20, Short.MAX_VALUE))
        );

        JTabbedPane tabParameters; 
        tabParameters = new JTabbedPane(JTabbedPane.TOP);
        tabParameters.addTab("ARIMA/(E)GARCH Model",armagarchPanel);
        tabParameters.addTab("(Factor) Stochastic Volatility Model",svparameterPanel);
        tabParameters.addTab("HEAVY Model",heavyparameterPanel);
        


        parameterDialog.getContentPane().add(tabParameters);  
        parameterDialog.pack(); 
        parameterDialog.setLocationRelativeTo(frame); 
        parameterDialog.setModal(false); 
        parameterDialog.setVisible(false);


        parameterButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
            parameterDialog.setModal(false); parameterDialog.setVisible(true);



            }
        });

    }





    @SuppressWarnings({ "rawtypes", "unchecked" })
	public void initControlPanel() 
    {
        int i;
        JLabel space = new JLabel(" | ");
        JLabel space1 = new JLabel(" | ");
        JLabel space2 = new JLabel(" | ");
        JLabel space3 = new JLabel(" | ");
        arimaParamPanel = new JPanel();
        arimaPlotCombo = new JComboBox<String>();
        arimaPlotHistLabel = new JLabel();
        bayesButton = new JButton();
        arText = new JTextField();
        maText = new JTextField();
        muText = new JTextField();
        chooseModelBox = new DisabledItemsComboBox();// = new JComboBox();
        chooseModelLabel = new JLabel();
        arLabel = new JLabel();
        maLabel = new JLabel();
        sigmaText = new JTextField();
        phisText = new JTextField();
        diffText = new JTextField();
        muLabel = new JLabel();
        diffLabel = new JLabel();
        sigmaLabel = new JLabel();
        gsText = new JTextField();
        gsLabel = new JLabel();
        phisLabel = new JLabel();
        fphiLabel = new JLabel();
        fsigmaLabel = new JLabel();
        fphiText = new JTextField();
        fsigmaText = new JTextField();
        mleButton = new JButton();
        modelDimensionPanel = new JPanel();
        pdimBox = new JComboBox<String>();
        pdimLabel = new JLabel();
        qdimLabel = new JLabel();
        qdimBox = new JComboBox<String>();
        ddimLabel = new JLabel();
        ddimBox = new JComboBox<String>();
        factorsLabel = new JLabel();
        factorsBox = new JComboBox<String>();
        forecastStepsLabel = new JLabel();
        forecastSampBox = new JComboBox<String>();
        forecastSampLabel = new JLabel();
        //forecastStepsSlider = new JScrollBar();
        forecastStepsBox = new JTextField();
        sampsCombo = new JComboBox(new String[] {"1000", "1500", "3000", "5000", "10000"}); 
        arimaHistSlider = new JSlider();
        arimaUpperCheck = new JCheckBox();
        arimaLowerCheck = new JCheckBox();
        computeLabel = new JLabel();
        sampsLabel = new JLabel();
        avg_forecast = new JCheckBox();
        arimaParamPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

        mcmcPanel = new JPanel();
        burnLabel = new JLabel();
        burnCombo = new JComboBox();
        
        modelParamPanel = new JPanel();
        forecastPanel = new JPanel();
        computePanel = new JPanel();
 
        arimaPlotCombo.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        //arimaPlotCombo.setModel(new DefaultComboBoxModel(new String[] { "sigma" }));

        arimaPlotHistLabel.setText("Parameter Histogram");

        bayesButton.setText("Bayes");
        bayesButton.setFont(new Font("Ubuntu", 0, 12));
        
        
        arText.setBackground(new Color(0, 0, 0));
        arText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        arText.setForeground(new Color(46, 239, 14));
        arText.setText("0.0 ");

        maText.setBackground(new Color(0, 0, 0));
        maText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        maText.setForeground(new Color(7, 247, 11));
        maText.setText("0.0 ");

        muText.setBackground(new Color(0, 0, 0));
        muText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        muText.setForeground(new Color(11, 251, 17));
        muText.setText("0.0 ");

        AICCLabel = new JLabel(""); 
        AICCText = new JTextField(); 
        AICCText.setBackground(new Color(0, 0, 0));
        AICCText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        AICCText.setForeground(new Color(150, 150, 150));
        AICCText.setText("0.0 ");        
        
        //chooseModelBox.setModel(new DefaultComboBoxModel(new String[] { "ARIMA", "GARCH", "EGARCH", "SVM", "Factor SVM", "Dynamic Factor" }));
        chooseModelBox.addItem("ARIMA");
        chooseModelBox.addItem("GARCH");
        chooseModelBox.addItem("EGARCH");
        chooseModelBox.addItem("Stochastic Vol.");
        chooseModelBox.addItem("Factor Stochastic Vol.",true);
        chooseModelBox.addItem("HEAVY",true);
       
        //chooseModelBox.addItem("HEAVY Model");
 
        ActionListener modelActionListener = new ActionListener() {
        public void actionPerformed(ActionEvent event) 
        {
          if(event.getSource() == chooseModelBox)
          {
                  
            if(chooseModelBox.getSelectedIndex() == 5)
            {arimaPlotCombo.setEnabled(false); bayesButton.setEnabled(false);}         
            else
            {arimaPlotCombo.setEnabled(true); bayesButton.setEnabled(true);}
          }
         }
        };

        chooseModelLabel.setText("Choose Model:");
        chooseModelLabel.setToolTipText("Choose model to apply to time series data and forecasting. Factor models only available for multivariate data.");
        chooseModelBox.addActionListener(modelActionListener);

        arLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        arLabel.setText("AR:");

        maLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        maLabel.setText("MA:");

        sigmaText.setBackground(new Color(0, 0, 0));
        sigmaText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        sigmaText.setForeground(new Color(11, 251, 17));
        sigmaText.setText("0.0 ");

        phisText.setBackground(new Color(0, 0, 0));
        phisText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        phisText.setForeground(new Color(11, 251, 17));
        phisText.setText("0.0 ");

        diffText.setBackground(new Color(0, 0, 0));
        diffText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        diffText.setForeground(new Color(11, 251, 17));
        diffText.setText("0.0 ");

        muLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        muLabel.setText("    \u03BC:");
        
        sampsLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        sampsLabel.setText("samps:");

        diffLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        diffLabel.setText("d:");

        sigmaLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        sigmaLabel.setText("\u03C3:");

        gsText.setBackground(new Color(0, 0, 0));
        gsText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        gsText.setForeground(new Color(7, 247, 11));
        gsText.setText("0.0 ");

        gsLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        gsLabel.setText("GS:");

        phisLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        phisLabel.setText("\u03C6:");

        fphiLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        fphiLabel.setText("f \u03C6:");

        fsigmaLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        fsigmaLabel.setText("f \u03C3:");

        fphiText.setBackground(new Color(0, 0, 0));
        fphiText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        fphiText.setForeground(new Color(11, 251, 17));
        fphiText.setText("0.0 ");

        fsigmaText.setBackground(new Color(0, 0, 0));
        fsigmaText.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        fsigmaText.setForeground(new Color(11, 251, 17));
        fsigmaText.setText("0.0 ");

        mleButton.setText("MLE");
        mleButton.setFont(new Font("Ubuntu", 0, 12));
        
        modelDimensionPanel.setBorder(BorderFactory.createTitledBorder("Model Dimensions"));

        pdimBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        //pdimBox.setModel(new DefaultComboBoxModel(new String[] { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13" }));
        for(i=0;i<13;i++)
        {pdimBox.addItem(""+(i+1));}

 

        pdimLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        pdimLabel.setText("p:");
        pdimBox.setSelectedIndex(1);
       
        qdimLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        qdimLabel.setText("q:");
        
        
        qdimBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        //qdimBox.setModel(new DefaultComboBoxModel(new String[] { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13" }));
        for(i=0;i<13;i++)
        {qdimBox.addItem(""+(i+1));}
        qdimBox.setSelectedIndex(1);
        
        ddimLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        ddimLabel.setText("d:");



        ddimBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        //ddimBox.setModel(new DefaultComboBoxModel(new String[] { "0", "1"}));
        for(i=0;i<2;i++)
        {ddimBox.addItem(""+(i+1));}
        ddimBox.setSelectedIndex(0);
 

        
        timePlot = new JCheckBox[5];
        for(i=0;i<5;i++) 
        { 
          timePlot[i] = new JCheckBox("F_"+(i+1)+":"); 
          timePlot[i].setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
          timePlot[i].setSelected(false); 
          timePlot[i].setEnabled(false);
          timePlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
        }

        seriesPlot = new JCheckBox[5];
        for(i=0;i<5;i++) 
        { 
          seriesPlot[i] = new JCheckBox("y_"+(i+1)+":"); 
          seriesPlot[i].setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
          seriesPlot[i].setSelected(false); 
          seriesPlot[i].setEnabled(false);
          seriesPlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
        }       
 
        alphasPlot = new JCheckBox[3]; factorPlot = new JCheckBox[3];
        for(i=0;i<3;i++)
        {
          alphasPlot[i] = new JCheckBox("h_"+(i+1)+":");
          alphasPlot[i].setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
          alphasPlot[i].setSelected(false); 
          alphasPlot[i].setEnabled(false);
          alphasPlot[i].setHorizontalTextPosition(JMenuItem.LEFT); 

          factorPlot[i] = new JCheckBox("f_"+(i+1)+":");
          factorPlot[i].setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
          factorPlot[i].setSelected(false); 
          factorPlot[i].setEnabled(false);
          factorPlot[i].setHorizontalTextPosition(JMenuItem.LEFT);
 
        }

        residualPlot = new JCheckBox("r:");
        residualPlot.setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
        residualPlot.setSelected(false); 
        residualPlot.setEnabled(false);
        residualPlot.setHorizontalTextPosition(JMenuItem.LEFT);

        predictivePlot = new JCheckBox("e:");
        predictivePlot.setFont(new Font("Ubuntu", 0, 10)); // NOI18N 
        predictivePlot.setSelected(false); 
        predictivePlot.setEnabled(false);
        predictivePlot.setHorizontalTextPosition(JMenuItem.LEFT);


       BevelBorder timeSelBorder = new BevelBorder(BevelBorder.RAISED);
       JPanel timeGrid = new JPanel(); 
       timeGrid.setSize(new Dimension(600, 50));      
       GroupLayout paramLayout = new GroupLayout(timeGrid);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);        
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()
          .addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
          .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(space)
          .addComponent(seriesPlot[0]).addComponent(seriesPlot[1]).addComponent(seriesPlot[2])
          .addComponent(seriesPlot[3]).addComponent(seriesPlot[4]).addComponent(space1)
          .addComponent(alphasPlot[0]).addComponent(alphasPlot[1]).addComponent(alphasPlot[2])
          .addComponent(space2).addComponent(factorPlot[0]).addComponent(factorPlot[1]).addComponent(factorPlot[2])
          .addComponent(space3).addComponent(residualPlot).addComponent(predictivePlot));
        paramLayout.setVerticalGroup(
         paramLayout.createSequentialGroup()
         .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
          .addComponent(timePlot[0]).addComponent(timePlot[1]).addComponent(timePlot[2])
          .addComponent(timePlot[3]).addComponent(timePlot[4]).addComponent(space)
          .addComponent(seriesPlot[0]).addComponent(seriesPlot[1]).addComponent(seriesPlot[2])
          .addComponent(seriesPlot[3]).addComponent(seriesPlot[4]).addComponent(space1)
          .addComponent(alphasPlot[0]).addComponent(alphasPlot[1]).addComponent(alphasPlot[2])
          .addComponent(space2).addComponent(factorPlot[0]).addComponent(factorPlot[1]).addComponent(factorPlot[2])
          .addComponent(space3).addComponent(residualPlot).addComponent(predictivePlot)));    
        timeGrid.setLayout(paramLayout); 
        timeGrid.setBorder(timeSelBorder);  



//==========================================================



        arimaHistSlider.setFont(new java.awt.Font("Ubuntu", 0, 3)); // NOI18N

        arimaUpperCheck.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        arimaUpperCheck.setSelected(true);
        arimaUpperCheck.setText("upper");
        arimaUpperCheck.setHorizontalTextPosition(SwingConstants.LEFT);

        arimaLowerCheck.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        arimaLowerCheck.setText("lower");

        mcmcPanel.setBorder(BorderFactory.createTitledBorder("MCMC Sampling"));

        sampsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sampsLabel.setText("Samples:");

        sampsCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        sampsCombo.setModel(new DefaultComboBoxModel(new String[] { "1000", "1500", "2000", "3000", "5000", "10000" }));

        burnLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        burnLabel.setText("Burnin:");

        burnCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        burnCombo.setModel(new DefaultComboBoxModel(new String[] { "1000", "1500", "2000", "3000", "5000", "10000" }));

        GroupLayout mcmcPanelLayout = new GroupLayout(mcmcPanel);
        mcmcPanelLayout.setAutoCreateGaps(true);
        mcmcPanelLayout.setAutoCreateContainerGaps(true); 
        mcmcPanel.setLayout(mcmcPanelLayout);
        mcmcPanelLayout.setHorizontalGroup(
            mcmcPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(mcmcPanelLayout.createSequentialGroup()
                .addComponent(sampsLabel)
                .addGap(6, 6, 6)
                .addComponent(sampsCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(burnLabel)
                .addGap(3, 3, 3)
                .addComponent(burnCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(13, Short.MAX_VALUE))
        );
        mcmcPanelLayout.setVerticalGroup(
            mcmcPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(mcmcPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(mcmcPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(sampsLabel)
                    .addComponent(sampsCombo, GroupLayout.PREFERRED_SIZE, 23, GroupLayout.PREFERRED_SIZE)
                    .addComponent(burnLabel)
                    .addComponent(burnCombo, GroupLayout.PREFERRED_SIZE, 23, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        forecastPanel.setBorder(BorderFactory.createTitledBorder("Forecasting"));

        forecastStepsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsLabel.setText("Forecast Steps:");

        forecastSampLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        forecastSampLabel.setText("Forecast Samples:");
        forecastSampLabel.setToolTipText("Number of forecasting samples");

        forecastStepsSlider = new JScrollBar(JScrollBar.HORIZONTAL,0,2,0,48);
        forecastStepsSlider.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsSlider.setUnitIncrement(1);       
        forecastStepsSlider.setMinimum(0);
        forecastStepsSlider.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsSlider.setOrientation(JScrollBar.HORIZONTAL);

        forecastSampBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        forecastSampBox.setModel(new DefaultComboBoxModel(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" }));

        forecastStepsBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsBox.setText("12");






        factorsLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        factorsLabel.setText("Factors:");
        factorsLabel.setToolTipText("Choose number of factors to use in stochastic factor models");

        factorsBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        for(i=0;i<3;i++)
        {factorsBox.addItem(""+(i+1));}
        //factorsBox.setModel(new DefaultComboBoxModel(new String[] { "1", "2", "3" }));




        forecastStepsLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsLabel.setText("Forecast Steps:");

        forecastSampBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        for(i=0;i<5;i++)
        {forecastSampBox.addItem(""+(i+1));} 
        //forecastSampBox.setModel(new DefaultComboBoxModel(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" }));




        forecastSampLabel.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        forecastSampLabel.setText("Forecast Samples:");
        forecastSampLabel.setToolTipText("Number of forecasting samples");


  
    
        forecastStepsBox.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        forecastStepsBox.setText("0");

        histogram.setBackground(new Color(0, 0, 0));

        GroupLayout arimaTheatreLayout2 = new GroupLayout(histogram);
        histogram.setLayout(arimaTheatreLayout2);
        arimaTheatreLayout2.setHorizontalGroup(
            arimaTheatreLayout2.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        arimaTheatreLayout2.setVerticalGroup(
            arimaTheatreLayout2.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGap(0, 181, Short.MAX_VALUE)
        );

        arimaHistSlider.setFont(new Font("Ubuntu", 0, 3)); // NOI18N
        arimaHistSlider.setMaximum(1000); arimaHistSlider.setMinimum(0); arimaHistSlider.setValue(500);
        
        arimaUpperCheck.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        arimaUpperCheck.setSelected(false);
        arimaUpperCheck.setText("upper");
        arimaUpperCheck.setHorizontalTextPosition(SwingConstants.LEFT);

        arimaLowerCheck.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        arimaLowerCheck.setText("lower");
        arimaUpperCheck.setSelected(true);
 
        avg_forecast.setFont(new Font("Ubuntu", 0, 12)); // NOI18N
        avg_forecast.setText("     Avg:");
        avg_forecast.setSelected(false);
        avg_forecast.setHorizontalTextPosition(SwingConstants.LEFT);








        GroupLayout forecastPanelLayout = new GroupLayout(forecastPanel);
        forecastPanelLayout.setAutoCreateGaps(true);
        forecastPanelLayout.setAutoCreateContainerGaps(true); 

        forecastPanel.setLayout(forecastPanelLayout);
        forecastPanelLayout.setHorizontalGroup(
            forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(forecastPanelLayout.createSequentialGroup()
                .addGroup(forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(forecastSampLabel)
                    .addComponent(forecastStepsLabel))
                //.addGap(18, 18, 18)
                .addGroup(forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(forecastPanelLayout.createSequentialGroup()
                        .addComponent(forecastStepsSlider, GroupLayout.PREFERRED_SIZE, 91, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(forecastStepsBox, GroupLayout.PREFERRED_SIZE, 40, GroupLayout.PREFERRED_SIZE))
                    .addGroup(forecastPanelLayout.createSequentialGroup()
                        .addComponent(forecastSampBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        //.addGap(0, 0, Short.MAX_VALUE)
                      ))
                .addContainerGap())
        );
        forecastPanelLayout.setVerticalGroup(
            forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(forecastPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(forecastStepsLabel)
                    .addComponent(forecastStepsSlider, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(forecastPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(forecastSampLabel)
                    .addComponent(forecastSampBox, GroupLayout.PREFERRED_SIZE, 15, GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
            .addGroup(forecastPanelLayout.createSequentialGroup()
                .addComponent(forecastStepsBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(0, 0, Short.MAX_VALUE))
        );

        computePanel.setBorder(BorderFactory.createTitledBorder("Compute Model:"));
        computePanel.setToolTipText("");
        computePanel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N

        bayesButton.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        bayesButton.setText("Bayesian");

        mleButton.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        mleButton.setText("MLE");

        GroupLayout jPanel1Layout = new GroupLayout(computePanel);
        jPanel1Layout.setAutoCreateGaps(true);
        jPanel1Layout.setAutoCreateContainerGaps(true); 

        computePanel.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(mleButton, GroupLayout.PREFERRED_SIZE, 100, GroupLayout.PREFERRED_SIZE)
                //.addGap(18, 18, 18)
                .addComponent(bayesButton, GroupLayout.PREFERRED_SIZE, 100, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                //.addGap(0, 12, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(mleButton)
                    .addComponent(bayesButton)))
        );

        modelParamPanel.setBorder(BorderFactory.createTitledBorder("Model and Parameters"));

        chooseModelLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        chooseModelLabel.setText("Model:");
        chooseModelLabel.setToolTipText("Choose model to apply to  time series data and forecasting. Factor models only available for multivariate data.");

        chooseModelBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        chooseModelBox.setModel(new DefaultComboBoxModel(new String[] { "ARIMA", "GARCH", "EGARCH", "SVM", "Factor SVM", "Dynamic Factor" }));
        chooseModelBox.setToolTipText("");

        factorsLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        factorsLabel.setText("Factors:");
        factorsLabel.setToolTipText("Choose number of factors to use in stochastic factor models");

        factorsBox.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        factorsBox.setModel(new DefaultComboBoxModel(new String[] { "1" }));

        arimaPlotHistLabel.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        arimaPlotHistLabel.setText("Histograms:");

        arimaPlotCombo.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        arimaPlotCombo.setModel(new DefaultComboBoxModel(new String[] { "sigma" }));

        /*parameterButton.setFont(new java.awt.Font("Ubuntu", 0, 12)); // NOI18N
        parameterButton.setText("Parameter Values");
        parameterButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                parameterButtonActionPerformed(evt);
            }
        });*/


        GroupLayout modelParamPanelLayout = new javax.swing.GroupLayout(modelParamPanel);
        modelParamPanelLayout.setAutoCreateGaps(true);
        modelParamPanelLayout.setAutoCreateContainerGaps(true);
        modelParamPanel.setLayout(modelParamPanelLayout);
        modelParamPanelLayout.setHorizontalGroup(
            modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(modelParamPanelLayout.createSequentialGroup()
                .addGroup(modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(parameterButton, javax.swing.GroupLayout.PREFERRED_SIZE, 132, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(modelParamPanelLayout.createSequentialGroup()
                        .addGroup(modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addGroup(modelParamPanelLayout.createSequentialGroup()
                                .addComponent(arimaPlotHistLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(arimaPlotCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(modelParamPanelLayout.createSequentialGroup()
                                .addComponent(chooseModelLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(chooseModelBox, javax.swing.GroupLayout.PREFERRED_SIZE, 107, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(factorsLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(factorsBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                /*.addGap(0, 23, Short.MAX_VALUE)*/)
        );
        modelParamPanelLayout.setVerticalGroup(
            modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(modelParamPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(chooseModelLabel)
                    .addComponent(chooseModelBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(factorsLabel)
                    .addComponent(factorsBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                //.addGap(18, 18, 18)
                .addGroup(modelParamPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(arimaPlotHistLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(arimaPlotCombo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(parameterButton))
        );








        /*
        GroupLayout modelParamPanelLayout = new GroupLayout(modelParamPanel);
        modelParamPanelLayout.setAutoCreateGaps(true);
        modelParamPanelLayout.setAutoCreateContainerGaps(true);       
        modelParamPanel.setLayout(modelParamPanelLayout);
        modelParamPanelLayout.setHorizontalGroup(
            modelParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(modelParamPanelLayout.createSequentialGroup()
                .addComponent(chooseModelLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(chooseModelBox, GroupLayout.PREFERRED_SIZE, 100, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(factorsLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(factorsBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addGap(0, 0, Short.MAX_VALUE))
            .addGroup(modelParamPanelLayout.createSequentialGroup()
                .addComponent(parameterButton, GroupLayout.PREFERRED_SIZE, 100, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(arimaPlotHistLabel)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(arimaPlotCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
        );
        modelParamPanelLayout.setVerticalGroup(
            modelParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(modelParamPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(modelParamPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(chooseModelLabel)
                    .addComponent(chooseModelBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(factorsLabel)
                    .addComponent(factorsBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                //.addGap(18, 18, 18)
                .addGroup(modelParamPanelLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
                    .addComponent(parameterButton)
                    .addComponent(arimaPlotHistLabel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(arimaPlotCombo, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );*/

        GroupLayout arimaParamPanelLayout = new GroupLayout(arimaParamPanel);
        arimaParamPanel.setLayout(arimaParamPanelLayout);
        arimaParamPanelLayout.setHorizontalGroup(
            arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(arimaParamPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                    .addComponent(modelParamPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(computePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(forecastPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                    .addComponent(mcmcPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addGroup(arimaParamPanelLayout.createSequentialGroup()
                        .addComponent(arimaLowerCheck)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(arimaHistSlider, GroupLayout.PREFERRED_SIZE, 321, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(arimaUpperCheck))
                    .addComponent(histogram, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        arimaParamPanelLayout.setVerticalGroup(
            arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(arimaParamPanelLayout.createSequentialGroup()
                .addComponent(histogram, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                    .addComponent(arimaUpperCheck, GroupLayout.Alignment.TRAILING)
                    .addGroup(GroupLayout.Alignment.TRAILING, arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(arimaLowerCheck)
                        .addComponent(arimaHistSlider, GroupLayout.PREFERRED_SIZE, 19, GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
            .addGroup(arimaParamPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(arimaParamPanelLayout.createParallelGroup(GroupLayout.Alignment.TRAILING)
                    .addGroup(arimaParamPanelLayout.createSequentialGroup()
                        .addComponent(forecastPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mcmcPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                    .addGroup(arimaParamPanelLayout.createSequentialGroup()
                        .addComponent(modelParamPanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(computePanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
               ) //.addGap(82, 82, 82))
        );

        GroupLayout layout = new GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(arimaParamPanel, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(arimaParamPanel, GroupLayout.PREFERRED_SIZE, 224, GroupLayout.PREFERRED_SIZE)
                .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
         );       
  




































        //========================================================================

        pdimBox.addActionListener(new MyActionListener()); 
        qdimBox.addActionListener(new MyActionListener()); 
        ddimBox.addActionListener(new MyActionListener()); 
        chooseModelBox.addActionListener(new MyActionListener()); 
        arimaPlotCombo.addActionListener(new MyActionListener()); 
        factorsBox.addActionListener(new MyActionListener()); 
        forecastSampBox.addActionListener(new MyActionListener()); 
        sampsCombo.addActionListener(new MyActionListener());
        
 
        AdjustmentListener al = new AdjustmentListener()  {
         public void adjustmentValueChanged(AdjustmentEvent e) {
          
           if(f_ready)
           {
             //n_sims = forecastSampBox.getSelectedIndex() + 1; 
             if(e.getAdjustable() == forecastStepsSlider)
             { 
                n_steps = forecastStepsSlider.getValue(); forecastStepsBox.setText(""+n_steps);
                cr.setNPredictiveSims(n_steps,n_sims);
                cr.setNForecastSteps(n_steps);
                if(n_steps > 0) {computeForecasts();}
             }
           }
           else
           {
             if(e.getAdjustable() == forecastStepsSlider)
             { 
                n_sims = forecastSampBox.getSelectedIndex() + 1;
                n_steps = forecastStepsSlider.getValue(); forecastStepsBox.setText(""+n_steps);
                cr.setNPredictiveSims(n_steps,n_sims);
                cr.setNForecastSteps(n_steps);
             }

           }
         }
        };
        forecastStepsSlider.addAdjustmentListener(al);
 
 
       ChangeListener cl = new ChangeListener() 
       {
         public void stateChanged(ChangeEvent e) 
         {   
 
          //---- gives number in range of 0 n_samp-10 
          int newval = (int)(((double)(n_samps - 10)/(double)n_samps)*arimaHistSlider.getValue());
 
          if(arimaUpperCheck.isSelected() && !arimaLowerCheck.isSelected()) 
          {
           end = n_samps; start = newval;
          }        
          else if(!arimaUpperCheck.isSelected() && arimaLowerCheck.isSelected())
          {
            start = 0; end = 10 + newval;
          }   
              
          histogram.colorShadeBoxes(start, end, 1000);
 
          if(!arimaHistSlider.getValueIsAdjusting()) {changeIntervalSample();}
         }
       };
       arimaHistSlider.addChangeListener(cl);
 
 
       ItemListener il = new ItemListener()
       {  
         public void itemStateChanged(ItemEvent e)
         {

          int i; boolean sel; //computeFilter = true;
          Object source = e.getItemSelectable();
          if(e.getStateChange() == ItemEvent.DESELECTED){sel = false;}
          else{sel = true;}

          if(source == avg_forecast)
          { 
           if(f_ready)
           {
             if(sel) {cronos_canvas.setForecasts(cr.predictiveAvg, n_steps, 1);}         
             else {cronos_canvas.setForecasts(cr.predictive, n_steps, n_sims);} 
           }
          }
 
          for(i=0;i<5;i++) 
          { 
            if(source == timePlot[i])
            {cronos_canvas.setplot_fore(i, sel);}

            else if(source == seriesPlot[i])
            {cronos_canvas.setplot(i, sel);}
       
          }

          for(i=0;i<3;i++) 
          {
            if(source == alphasPlot[i])
            {cronos_canvas.setplot_alpha(i,sel);}    
          
            else if(source == factorPlot[i])
            {cronos_canvas.setplot_factors(i,sel);} 
          }

          if(source == residualPlot)
          {cronos_canvas.plotResiduals(sel);}
 
          if(source == predictivePlot)
          {cronos_canvas.plotPredictive(sel);}          

        } 
       };
       avg_forecast.addItemListener(il);
       for(i=0;i<5;i++) {timePlot[i].addItemListener(il);}
       for(i=0;i<3;i++) {alphasPlot[i].addItemListener(il); factorPlot[i].addItemListener(il);}
       for(i=0;i<5;i++) {seriesPlot[i].addItemListener(il);}
       residualPlot.addItemListener(il);
       predictivePlot.addItemListener(il);

       ActionListener all = new ActionListener()
       {
         public void actionPerformed(ActionEvent e)
         {
          if(series_ready)
          {
           if(e.getSource() == mleButton)
           {computeModel(1); methodtype = 1; arimaPlotCombo.setEnabled(false);}
           else if(e.getSource() == bayesButton)
           {computeModel(0); methodtype = 0; arimaPlotCombo.setEnabled(true); setupButtons();}
          }
         }
       };
       mleButton.addActionListener(all);
       bayesButton.addActionListener(all);
 
 

     Box timePane2 = Box.createVerticalBox();
     timePane2.add(timeGrid);
     timePane2.add(arimaParamPanel);

 
      Box timePane = Box.createVerticalBox();
      timePane.add(timeScrollPane,BorderLayout.NORTH);
      //timePane.add(timeGrid,BorderLayout.SOUTH);
 
      paramLayout = new GroupLayout(this);
        paramLayout.setAutoCreateGaps(true);
        paramLayout.setAutoCreateContainerGaps(true);
        paramLayout.setHorizontalGroup(paramLayout.createSequentialGroup()           
           .addGroup(paramLayout.createParallelGroup(GroupLayout.Alignment.LEADING)
             .addComponent(timePane)
             .addComponent(timePane2)));

        paramLayout.setVerticalGroup(
        paramLayout.createSequentialGroup()
           .addComponent(timePane)
           .addComponent(timePane2));
        this.setLayout(paramLayout);
 
    }


    
    
    public void readData(File file)
    {
          
       String strline; Double D; 
       String[] tokens; String delims = "[ ]+";
       int n_toks;
       double[] values = new double[3000];      
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

           if(i < 3000) {values[i] = val; i++;}
           else 
           {System.out.println("Maximum times series length is 3000"); break;} 
           
         } 
         
         n_rep = 1;
         setTimeSeries(values, n_rep, i);   
         din.close();
         }
         catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
         catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}         

   }

	public int getN_factors() {
		return n_factors;
	}

	public void setN_factors(int n_factors) {
		this.n_factors = n_factors;
	}

}



