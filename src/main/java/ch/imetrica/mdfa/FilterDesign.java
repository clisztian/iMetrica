package ch.imetrica.mdfa;

/*----------------------------------------------------------------
 
   Filter design ---
   One can design three different types - lowpass, bandpass, and smoothed lowpass (ramp) 
   Cutoffs can be defined directly with omega0, omega1 or by defining periods for the bandpass p1,p2

----------------------------------------------------------------*/

public class FilterDesign
{


   double[] Gamma; //------ main output for symmetric filter 
   double[] bsym;  //------ symmetric filter coefficients   

   int K,K1;       //------ K length of filter on K1 points
   int L1;         //------ resolution for approximations
   int p1,p2;      //------ period parameters  2pi/p1, 2pi/p2
   double w0,w1;   //------ cutoff omega1, omega2
   int num1, num2;  //----- numerators for cutoffs
   int den1, den2;  //----- denominators for cutoffs
   double w2,w3;    //----- second set of cutoffs for multipass

   boolean low,high,ramp,band;

   public FilterDesign(int _K, int _L1, double _w0, double _w1)
   {  
     
     K=_K; K1 = _K+1; L1 = _L1;
     Gamma = new double[K1]; bsym = new double[L1+1];
  
     p1=300;p2=12;  num1 = 0; num2 = 1; den1 = 20; den2 = 6; 
     low=true; high=false; ramp=false; band=false;
     w0 = _w0; w1 = _w1; w2 = .70; w3 = 1.00;
     setBand(w0, w1);

   }

   public void setK(int _K)
   {K=_K; K1 = K+1; Gamma = new double[K1];}
  
   public void recomputeGamma()
   {setBand(w0,w1);}
     
   public void setL1(int _L)
   {L1 = _L; bsym = new double[L1+1];}

   public void setFilterType(int which)
   {
     if(which == 0) {low = true; high = false; ramp = false; band = false;}
     else if(which == 1) {low = false; high = true; ramp = false; band = false;}
     else if(which == 2) {low = false; high = false; ramp = true; band = false;}
     else if(which == 3) {low = false; high = false; ramp = false; band = true;}
     setBand(w0,w1);
   }

   public void setBand(double _w0, double _w1)
   {
       w0 = _w0; w1 = _w1; 
    
       if(low)   
       {setLowPassTrend(w1);}
       else if(high)
       {setHighPass(w0,w1,w2,w3);}
       else if(ramp)
       {setRampTrend(w0, w1);}
       else
       {setBandPassTrend(w0, w1);} 
   }
   
   public void setBand2(double _w2, double _w3)
   {
     w2 = _w2; w3 = _w3;      
     if(high) {setHighPass(w0,w1,w2,w3);}   
   }

  //--------set ramp trend -------------------------
   public void setRampTrend(double _w0, double _w1)
   {
      int i,k; 
      double sum,om,coeff0;
 
      if(_w0 < _w1)
      {
       w0 = _w0; w1 = _w1; 
       bsym = new double[L1+1]; Gamma = new double[K1];
              
       //-----compute symmetric filter----
       coeff0 = .5*(w1 + w0)/Math.PI; bsym[0] = coeff0;  sum = bsym[0];
       for(i=1;i<=L1;i++)
       {
          bsym[i] = (Math.cos(w1*i) - Math.cos(w0*i))/(i*i);
          bsym[i] = -1.0/(Math.PI*(w1 - w0))*bsym[i]; 
          sum= sum+bsym[i];
       } 
       for(i=0;i<=L1;i++)
       {
         bsym[i] = bsym[i]/(sum+(sum-coeff0));
       }     
       //---- now compute Gamma -------
       for(k=0; k<=K;k++)
       {       
         om = (k*Math.PI/K);
         if(om <= w0) {Gamma[k] = 1.0;}
         else if(om > w0 && om <= w1)
         {Gamma[k] = (w1 - om)/(w1 - w0);}
         else
         {Gamma[k] = 0.0;}
       }           
      }
      else 
      {System.out.println("omega_0 must be less than omega_1");}
   }


   public void setLowPassTrend(double _w1)
   {
      int i,k; 
      double sum,cutoff;
 
      if(_w1 < Math.PI)
      {
        w0 = 0.0; w1 = _w1; cutoff = w1;
        bsym = new double[L1+1]; Gamma = new double[K1];

        bsym[0] = cutoff/Math.PI; sum= bsym[0];
        for(i=1;i<=L1;i++)
        {bsym[i] = (1/Math.PI)*Math.sin(cutoff*i)/i; sum= sum+bsym[i];} 
        for(i=0;i<=L1;i++)
        {bsym[i] = bsym[i]/(sum+(sum-cutoff/Math.PI));}     
        for(k=0; k<=K;k++)
        {       
         if((k*Math.PI/K) <= cutoff) {Gamma[k] = 1.0;}
         else {Gamma[k] = 0.0;}
        }           
      }
      else 
      {System.out.println("omega_1 must be less than PI");}
   }


   public void setBandPassTrend(double _w0, double _w1)
   {
      int i,k; 
      double om,coeff0;
 
      if(_w0 < _w1)
      {
       w0 = _w0; w1 = _w1; 
       bsym = new double[L1+1]; Gamma = new double[K1];
 
       //-----compute symmetric filter----
       coeff0 = (w1 - w0)/Math.PI; bsym[0] = coeff0;  
       for(i=1;i<=L1;i++)
       {
          bsym[i] = (Math.sin(w1*i) - Math.sin(w0*i))/(Math.PI*i);
       } 

       //---- now compute Gamma -------
       for(k=0; k<=K;k++)
       {       
         om = (k*Math.PI/K);
         if(om < w0) {Gamma[k] = 0.0;}
         else if(om >= w0 && om <= w1)
         {Gamma[k] = 1.0;}
         else
         {Gamma[k] = 0.0;}
       }
      }
      else 
      {System.out.println("omega_0 must be less than omega_1");}       

   }


   public void setHighPass(double _w0, double _w1, double _w2, double _w3)
   {
      int k; 
      double om;
      w0 = _w0; w1 = _w1; w2 = _w2; w3 = _w3; 

      if(w1 >= 0 && w0 >= 0)
      {
        
        Gamma = new double[K1];
        for(k=0; k<=K;k++) 
        {   
          
          Gamma[k] = 0.0;
          om = (k*Math.PI/K);
          
          if(om >= w1)
          {Gamma[k] = 1.0;}         
        }  
      
//       if(w0 >= 0 && w2 >= 0)
//       {      
//        if((w1 > w0) && (w3 > w2))
//        {
//          //System.out.println("w0 = " + w0 + ", w1 = " + w1 + ", w2 = " + w2 + ", w3 = " + w3);
//          Gamma = new double[K1];
//          //---- now compute Gamma -------
//          for(k=0; k<=K;k++)
//          {       
//           Gamma[k] = 0.0;
//           om = (k*Math.PI/K);
//          
//           
//          
//           //if(om >= w0 && om <= w1) {Gamma[k] = 1.0;}// System.out.println("Gamma[" + k + "] = 1.0 at " + om);}
//          
//           //if(om >= w2 && om <= w3) {Gamma[k] = 1.0;}// System.out.println("Gamma[" + k + "] = 1.0 at " + om);}         
//          }
//          setGeneralSymmetric();
//        }
//        
       }
       else 
       {System.out.println("omega_0 must be greater than 0");}       
   }

   public void setGeneralSymmetric()
   {

       int i,k,n; double sum = 0.0; double sum2 = 0.0;
       bsym = new double[L1+1]; 
       //System.out.println(K + "  " + Gamma.length);

   
       for(k=0;k<=L1;k++)
       {
         sum=0.0;
         for(n=0;n<=K;n++)
         {
          sum = sum + Gamma[n]*Math.cos(Math.PI*n*k/K);
         }     
         if(k==0) {bsym[0] = sum;}
         else
         {bsym[k] = 2.0*sum;} 
         sum2 = sum2 + bsym[k];
       }
       for(i=0;i<=L1;i++)
       {bsym[i] = bsym[i]/(sum2-bsym[0]/2.0);}    
   }

   public void setGeneralSymmetric(double[] _Gamma)
   {

       int i,k,n; double sum = 0.0; double sum2 = 0.0;
       bsym = new double[L1+1]; setK(_Gamma.length); 
       Gamma = new double[_Gamma.length];
       
       //System.out.println(K + "  " + Gamma.length);
       System.arraycopy(_Gamma,0,Gamma,0,_Gamma.length);
   
       for(k=0;k<=L1;k++)
       {
         sum=0.0;
         for(n=0;n<=K;n++)
         {
          sum = sum + Gamma[n]*Math.cos(Math.PI*n*k/K);
         }     
         if(k==0) {bsym[0] = sum;}
         else
         {bsym[k] = 2.0*sum;} 
         sum2 = sum2 + bsym[k];
       }
       for(i=0;i<=L1;i++)
       {bsym[i] = bsym[i]/(sum2-bsym[0]/2.0);}    
   }

  public void setReproducingKernelFilter(int choice)
  {
     int i,k; 
     double sum,norm; 
     double x1,x;
     norm = 0.0;
     bsym = new double[L1+1]; Gamma = new double[K1];
     

     sum = 0; w0 = 0.0;
     for(i=0;i<=L1;i++)
     {
       x = 1.0*i/(L1+1.0);

       //System.out.println(x);
       switch(choice)
       {
         //Beta(2,3) 3
	case 1:  
	x1 = (Math.pow(Math.abs(x),2) - 1);
	bsym[i] = (35.0*((99.0*Math.pow(x,2.0))/16.0 - 27.0/16.0)*Math.pow(x1,3.0))/32.0;
	break;

	//Beta(3,3) 5
	case 2:  
	x1 = (Math.pow(Math.abs(x),3) - 1);
	bsym[i] = (70.0*Math.pow(x1,3.0)*((2778817527263606855434240.0*Math.pow(x,2))/476710297726655559415563.0 - 39697393246622955077632.0/21579478497914449191651.0))/81.0;
	break;

	//Beta(3,3) 4
	case 3:  
	x1 = (Math.pow(Math.abs(x),3) - 1);
	bsym[i] = (70.0*((282465768628677509120.0*Math.pow(x,2.0))/48457424548190474559.0 - 4035225266123964416.0/2193545967202037943.0)*Math.pow(x1,3.0))/81.0;
	break;

	//Beta(4,4) 3
	case 4:  
	x1 = (Math.pow(Math.abs(x),4.0) - 1.0);
	bsym[i] = -(3315.0*Math.pow(x1,4.0)*((112619678413006026309632.0*Math.pow(x,2.0))/18489259006422606298281.0 - 509591305036226363392.0/265396062293147458827.0))/4096.0;
	break;

	//Beta(2,2) 3
  
	case 5:  
	x1 = Math.pow(Math.abs(x),2);
	bsym[i] = -(15.0*((21.0*Math.pow(x,2.0))/4.0 - 7.0/4.0)*(x1 - 1.0)*(x1 - 1.0))/16.0;
	break;


	//Beta(2,2) 4
	case 6:  
	bsym[i] = -(15.0*((590295810358705651712.0*Math.pow(x,2.0))/112437297211182026391.0 - 590295810358705651712.0/337311891633546079173.0)*(Math.pow(Math.abs(x),2.0) - 1.0)*(Math.pow(Math.abs(x),2.0) - 1.0))/16.0;
	break;

	//Henderson Order 2
	case 7:  
	bsym[i] = -(15.0*(49.0*Math.pow(x,2.0) - 49.0)*(49.0*Math.pow(x,2.0) - 64.0)*(49.0*Math.pow(x,2.0) - 81.0))/3889424.0;
	break;

	//Henderson Order 3
	case 8:  
	bsym[i] = (15.0*((1191281797138502140297216.0*Math.pow(x,2.0))/230094499825316391106941.0 - 737739194587985339219968.0/432128694793886880859377.0)*(49.0*Math.pow(x,2.0) - 49.0)*(49.0*Math.pow(x,2.0) - 64.0)*(49.0*Math.pow(x,2.0) - 81.0))/3889424.0;
	break;

	//Henderson  Order 4
	case 9:  
	bsym[i] = (15.0*(49.0*Math.pow(x,2.0) - 49.0)*(49.0*Math.pow(x,2.0) - 64.0)*(49.0*Math.pow(x,2.0) - 81.0)*((15199586199233736933150816957440.0*Math.pow(x,2.0))/2935779924166721551069856602281.0 - 9412827852845683619693496581120.0/5513537906361891693472657521357.0))/3889424.0;
	break;

	//Gaussian Order 3
	case 10:  
	bsym[i] = -(7186705221432913.0*(Math.pow(x,2.0)/2.0 - 3.0/2.0))/(18014398509481984.0*Math.exp(Math.pow(x,2.0)/2.0));
	break;
	//Bandwidth = 2.8000
       
        default:
	x1 = (Math.pow(Math.abs(x),2.0) - 1.0);
	bsym[i] = (35.0*((99.0*Math.pow(x,2.0))/16.0 - 27.0/16.0)*Math.pow(x1,3.0))/32.0;
	break;
 
       }

       sum = sum + 2.0*bsym[i];
     }
       
     //----- Normalize coefficients -------------------
     for(i=0;i<=L1;i++)
     {bsym[i] = bsym[i]/sum; norm = norm + bsym[i];}   
     //(sum-bsym[0]/2.0)
   
     //System.out.println(norm);

     //----------- Now compute Gamma ---------------------
     boolean set = false;
     for(k=0;k<K1;k++)
     {

        sum = 0;

        for(i=0; i <= L1; i++)
        {sum = sum + bsym[i]*Math.cos(Math.PI*i*k/(double)K);}         
        Gamma[k] = (2.0*sum - bsym[0]);        
        //System.out.println(Gamma[k]);

        if(!set && Gamma[k] <= 0.0)
        {w1 = Math.PI*k/K; set = true;}
       
        /*
        if(k>0) 
          Gamma[k] = 2.0*sum - bsym[0];
        else 
          Gamma[0] = 2.0*sum;
        */  
     } 
     for(k=0;k<K1;k++)
     {Gamma[k] =  Gamma[k]/Gamma[0];} 
       
     //System.out.println(w1);
  }


  
}

