package ch.imetrica.usimx13;

import java.awt.*;
import javax.swing.*;

public class SARIMApseudoCanvas extends JPanel
{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//-------SARMA Model stuff------
    SARIMAmodelJava pseudo;
    int nObs;    
    double[] tseries; 
      
    //------ true model stuff
    double[] true_params;
    int[] true_dim;
    int n_true;
    int m_p, m_d, m_q, m_P, m_D, m_Q;

    //------ pseudo-true model stuff
    double[] pseudo_params;
    int[] pseudo_dim;
    int n_pseudo;

    //------ results
    double[] sigex_diags;
    double[] efficacies;
    double[] results;
    double minKL;  
    double[] sampQf_pseudo;
    double[] sampQI_pseudo;
    double[] sampQf_true;
    double[] sampG;
    
    //--------- 5-23-11 ---------------------- 
    // Added likelihood and LB states for estimated model 
    double[] lk;
    double[] lbv;

    //------ Other stuff
    Graphics2D g2d;
    int height, width,S;
    double T0, T1, X0, X1, Y0, Y1;
    double dataMax, dataMin, dataNorm;

    //------ even more stuff
    boolean timeDomain;
    boolean simData;
    boolean realData;
    int model; int seed; int burnin;
    double Qftrue_max, Qftrue_min;
    double Qfpseudo_max, Qfpseudo_min;
    double QIpseudo_max, QIpseudo_min;   
    double Qfpseudonorm, QIpseudonorm, Qftruenorm;
  

    public SARIMApseudoCanvas(int w, int h, int _nObs, int _burnin)
    {
       nObs = _nObs; burnin = _burnin; S = 12;
       pseudo = new SARIMAmodelJava(_nObs, _burnin);
       tseries = new double[nObs];
      

       true_dim = new int[6];
       pseudo_dim = new int[6];
       
       sigex_diags = new double[4];
       efficacies = new double[4];
       sampQf_pseudo = new double[300];
       sampQf_true = new double[300];
       sampQI_pseudo = new double[300];
       sampG = new double[300];

       //--------- 5-23-11 ---------------------- 
       // Added likelihood and LB states for estimated model 
       lk = new double[4]; lbv = new double[2];
       this.height = h; this.width = w;  
       timeDomain = true; 
       model = 0; seed = 500; minKL =0.0;
    }

    public void setTrueParameters(int nparams, int[] dim, double[] params, double _innvar)
    { 
        int j=0;
        n_true = nparams;
        true_params = new double[nparams];
        m_p = dim[0]; m_d = dim[1]; m_q = dim[2]; 
        m_P = dim[3]; m_D = dim[4]; m_Q = dim[5]; 
 
        if(m_p >= 1) {true_params[j] = params[0]; j++;}      
        if(m_p == 2) {true_params[j] = params[1]; j++;}      
        if(m_q >= 1) {true_params[j] = params[2]; j++;}      
        if(m_q == 2) {true_params[j] = params[3]; j++;}      
        if(m_P == 1) {true_params[j] = params[4]; j++;}            
        if(m_Q == 1) {true_params[j] = params[5]; j++;}       
       
        true_params[nparams-1] = _innvar;
        true_dim = dim;
        pseudo.SetSARIMAparams(true_params, true_dim, n_true, S);
    }

 


    public void setPseudoParameters(int[] dim)
    { 
        n_pseudo = dim[0] + dim[2] + dim[3] + dim[5] + 1;;
        pseudo_params = new double[n_pseudo];
        this.pseudo_dim = dim; 
        pseudo.changeModelDimension(dim);
    }

    public void changeModelDimension()
    {pseudo.changeModelDimension(this.pseudo_dim);}

    public void setNobs(int _nObs) {nObs = _nObs; pseudo.initializeUseless(_nObs);}
    public void setBurnin(int _burnin) {burnin = _burnin;}
    public void setModel(int _model) {model = _model;} 
    public void setSeed(int _seed) {seed = _seed;}
    public void setLag(int _lag) {pseudo.lag = _lag;}

    //--------  Simulate data from true dgp ------     
    public void simulateData(int seed)
    { 
     pseudo.SetSARIMAparams(true_params, true_dim, n_true, S);    
     tseries = pseudo.sampleSARIMAmodel(seed); simData = true;
    }

    //--------  Simulate data from true dgp ------     
    public void setData(int N, double[] data)
    {nObs = N; tseries = data; pseudo.initializeUseless(nObs);}

    public void computeEfficacies()
    {
       int i; //System.out.println("Computing Efficacies");
       results = pseudo.computeEfficacies(nObs, tseries, true_dim, true_params, n_true, pseudo_dim, n_pseudo, model);     
       for(i=0; i < n_pseudo; i++)
       {pseudo_params[i] = results[8+i];}
       for(i=0; i < 4; i++)
       {sigex_diags[i] = results[i];}  
       
       for(i=0; i < 4; i++)
       {efficacies[i] = results[4+i]; }//System.out.print(efficacies[i] + " ");}
       minKL = results[8+n_pseudo];
  
       for(i=0; i < 300; i++)
       {sampQf_pseudo[i] = results[9+n_pseudo+i];}
       for(i=0; i < 300; i++)
       {sampQI_pseudo[i] = results[9+n_pseudo+300+i];}      
       for(i=0; i < 300; i++)
       {sampQf_true[i] = results[9+n_pseudo+600+i];}
       for(i=0; i < 300; i++)
       {sampG[i] = results[9+n_pseudo+900+i];}

       //------ Added 5-23-11-----------------------------------------
       lbv[0] = results[8+n_pseudo+1+1200]; lbv[1] = results[8+n_pseudo+1+1201];
       lk[0] = results[8+n_pseudo+1+1202];  lk[1] = results[8+n_pseudo+1+1203];
       lk[2] = results[8+n_pseudo+1+1204];  lk[3] = results[8+n_pseudo+1+1205]; 

       //System.out.println(lk[0] + " " + lk[1] + " " + lk[2] + " " + lk[3]);

    }

    public void setSpectrumDataLimits()
    {
      int i; 
      
      Qfpseudo_max = -1000000;
      Qfpseudo_min = 10000000; 
      QIpseudo_max = -1000000;
      QIpseudo_min = 10000000; 
      Qftrue_max = -1000000;
      Qftrue_min = 10000000; 


      for(i=0; i < 300; i++)
      {
        if(sampQf_pseudo[i] > Qfpseudo_max) Qfpseudo_max = sampQf_pseudo[i];
        else if(sampQf_pseudo[i] < Qfpseudo_min) Qfpseudo_min = sampQf_pseudo[i];  

        if(sampQI_pseudo[i] > QIpseudo_max) QIpseudo_max = sampQI_pseudo[i];
        else if(sampQI_pseudo[i] < QIpseudo_min) QIpseudo_min = sampQI_pseudo[i]; 

        if(sampQf_true[i] > Qftrue_max) Qftrue_max = sampQf_true[i];
        else if(sampQf_true[i] < Qftrue_min) Qftrue_min = sampQf_true[i]; 
 
      }

      dataMin = min(Qfpseudo_min, QIpseudo_min, Qftrue_min);
      dataMax = max(Qfpseudo_max, QIpseudo_max, Qftrue_max); 
      dataNorm = Math.abs(dataMax - dataMin);      
    }

    public static double min(double d1, double d2, double d3)
    {double dt = Math.min(d1,d2); return Math.min(dt,d3);}
    public static double max(double d1, double d2, double d3)
    {double dt = Math.max(d1,d2); return Math.max(dt,d3);}

    /* public void go() {repaint();}
    public void paintComponent(Graphics g)
    {
        int i,j;
        int t0,t1,x0,x1;
        int n2 = 300; 		
        double point1, point2;
        Dimension ds = this.getSize();
        width = ds.width; height = ds.height;
 



        float[] dash = {44,46};
        BasicStroke bs = new BasicStroke(10f, BasicStroke.CAP_BUTT, 
                     BasicStroke.JOIN_BEVEL,10f,dash,dash[1]);

	super.paintComponent(g);
	g2d = (Graphics2D)g;
			
	//Draw Q(ftrue) Data
	g2d.setPaint(Color.GRAY);
        for(i = 0; i < 300; i++)
	{
	    t0 = (int)(((double)i/(double)n2)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)n2)*(double)width);
	    x0 = (int)(((sampQf_pseudo[i] - Qfpseudo_min)/(Qfpseudonorm*1.20))*(double)height);
	    x1 = (int)(((sampQf_pseudo[i+1] - Qfpseudo_min)/(Qfpseudonorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	}
        g2d.setPaint(Color.BLUE);
        for(i = 0; i < 300; i++)
	{
	    t0 = (int)(((double)i/(double)n2)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)n2)*(double)width);
	    x0 = (int)(((sampQI_pseudo[i] - QIpseudo_min)/(QIpseudonorm*1.20))*(double)height);
	    x1 = (int)(((sampQI_pseudo[i+1] - QIpseudo_min)/(QIpseudonorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	}

        g2d.setPaint(Color.BLUE);
        for(i = 0; i < 300; i++)
	{
	    t0 = (int)(((double)i/(double)n2)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)n2)*(double)width);
	    x0 = (int)(((sampQf_true[i] - Qftrue_min)/(Qftruenorm*1.20))*(double)height);
	    x1 = (int)(((sampQf_true[i+1] - Qftrue_min)/(Qftruenorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	}
      }
        /*g2d.setPaint(Color.BLUE);
        if(plot_diff)
        {
          for(i = 0; i < n2-1; i++)
	  {
            point1 = Math.abs(sampleQIf[i] - sampleQIf[i+n2]);
            point2 = Math.abs(sampleQIf[i+1] - sampleQIf[i+1+n2]);
	    t0 = (int)(((double)i/(double)n2)*(double)width);
	    t1 = (int)(((double)(i+1)/(double)n2)*(double)width);
	    x0 = (int)(((point1 - dataMin)/(dataNorm*1.20))*(double)height);
	    x1 = (int)(((point2 - dataMin)/(dataNorm*1.20))*(double)height);
	    g2d.drawLine(t0, (height-20) - x0, t1, (height-20) - x1);
	  }          
        }
	//g2d.setPaint(Color.GRAY);
        //g2d.setStroke(bs);
    }
    */
}
	
      

    

