package ch.imetrica.mdfaTradingStrategies;

import java.io.*;
import java.util.ArrayList;

public class VarRatio
{

  

  public VarRatio() {}

  ArrayList<Double> varianceRatio; 
  ArrayList<String> dates;
  
  public double var(Double[] rets)
  {
    int i; int n = rets.length;
    
    double mean = 0;
    for (i=0; i<n; i++)
    {mean = mean + (double)rets[i];}
    mean = mean/n;
    
    double sum = 0; 
    for (i=0; i<n; i++) 
    {final double v = (double)rets[i] - mean;  sum += v * v;} 
  
    return sum;
  }

  public double var_k(Double[] rets, int k)
  {
   int i; double val = 0; 
   //difference k
   Double[] diff_k = new Double[rets.length - k];
   
   for(i=k; i < rets.length; i++)
   {diff_k[i-k] = rets[i] - rets[i-k];}

   val = var(diff_k);
   return val;
  }


/*-----------------------------------

  file = price data file from 8:30-16:15
  k = differences 
  
------------------------------------*/

  public void computeVarRatio(File file, int k, String begin_time, String end_time, int span)
  {

   
       int i = 0; 
       String strline; 
       Double D; String Dt; 
       double var_ratio, delta, delta_k;
      
       String prev_date; 
       boolean get_data = false; boolean same_day = false;
       int n_toks; int count, tot_length; 
       String delims = "[,]+";
       String[] tokens; 
       String date_delims = "[ ]+"; 
       String[] date_tokens;
       int n_intervals = 0; 
       FileInputStream fin; DataInputStream din; BufferedReader br;  
       int sample_no=0; 
       
       int l = 0; int n_intervals_k = 0;
       ArrayList<String> dates = new ArrayList<String>();
       ArrayList<Double> close_series = new ArrayList<Double>();
       new ArrayList<Double>(); 
       ArrayList<Double> price = new ArrayList<Double>();    
       ArrayList<Double> local_close = new ArrayList<Double>();
       new ArrayList<Double>();
       ArrayList<String> local_dates = new ArrayList<String>();
       new ArrayList<Double>(); 
       new ArrayList<Double>(); 
       ArrayList<Double> price_diff_k = new ArrayList<Double>();
       ArrayList<String> date_diff_k = new ArrayList<String>(); 
       ArrayList<Double> variance_close = new ArrayList<Double>();
       ArrayList<Double> variance_diff_k = new ArrayList<Double>();
        
       varianceRatio = new ArrayList<Double>(); 
       int count2 = 0;
       try{
       fin = new FileInputStream(file);
       din = new DataInputStream(fin);
       br = new BufferedReader(new InputStreamReader(din));  

       //--- break up by time
       while((strline = br.readLine()) != null)
       {

          tokens = strline.split(delims); 
          n_toks = tokens.length; //System.out.println("Number of toks = "+n_toks);
          if(n_toks == 0)
          {System.out.println("End of file"); break;}
  
          Dt = new String(tokens[0]);
          
          //filter data for time 

          
          if(Dt.indexOf(begin_time) != -1 || get_data)
          {      
           dates.add(Dt);     
           D = new Double(tokens[1]); price.add(D);
           D = new Double(tokens[2]); close_series.add(D);
           get_data = true;  count2++;
          }      
          
          if(Dt.indexOf(end_time) != -1) //end of day, take differencing now
          {
            n_intervals = count2; count2 = 0;
            get_data = false; same_day = true; l = 0;
            //--- now go back until begin time for K-differencing
            while(same_day && dates.size() > 0)
            {
              //System.out.println(dates.get(dates.size()-1-l) + " - " + dates.get(dates.size()-1-l-k));
              price_diff_k.add(price.get(price.size() - 1 - l) - price.get(price.size() - 1 - l - k)); 
              date_diff_k.add(dates.get(dates.size()-1-l));
              
              if(dates.get(dates.size()-1-l-k).indexOf(begin_time) != -1)
              {same_day = false;}
              l++;
            }  
          }          
  
       }
       din.close();
       }catch(FileNotFoundException fe){System.out.println("File not found..." + fe);}
       catch(IOException ioe){System.out.println("IO procedure faulty..." + ioe);}    
   
   
       tot_length = close_series.size();
       count = 0; 
       int n_day_sample = 0;
      // double delta = var(close_series.toArray(new Double[0]));
       //filter data for time
       String date_stamp = dates.get(0);
       date_tokens = date_stamp.split(date_delims);
       prev_date = date_tokens[0];

       //System.out.println("total length = " + tot_length);
        
       //--compute the variance for lag 1 on span days 
       for(i = 0; i < tot_length; i++)
       {
         date_stamp = dates.get(i);
         date_tokens = date_stamp.split(date_delims);
         
         local_close.add(close_series.get(i));        
         local_dates.add(dates.get(i));
         //System.out.println(date_tokens[0]);
         
         if(prev_date.equals(date_tokens[0])) //no date change
         {n_day_sample++;}
         else //change in day, increase count
         {count++;  n_day_sample = 0;}
         
         //System.out.println(n_day_sample + "  " + n_intervals);
         
         if(n_day_sample == n_intervals-1)
         {
          if(count >= span)
          {         
            
            Double[] close = local_close.toArray(new Double[0]);
            String[] lesdates = local_dates.toArray(new String[0]);
            
            Double[] temp_close = new Double[n_intervals*span];
            String[] temp_dates = new String[n_intervals*span];
            
            System.arraycopy(close, n_intervals*sample_no, temp_close, 0, n_intervals*span);
            System.arraycopy(lesdates, n_intervals*sample_no, temp_dates, 0, n_intervals*span);
            
//             System.out.println("Printing sample no = " + sample_no); 
//             System.out.println();
//             
//             for(int m=0;m<temp_close.length;m++)
//             {
//               System.out.println(temp_dates[m] + " " + temp_close[m]);
//             }
                        
            delta = var(temp_close);          
            variance_close.add(delta);
            sample_no++;
           }
          }
          prev_date = date_tokens[0];  
        }
          

       ///---- now do close-k series and empty arrays 
       n_day_sample = 0; count = 0;  
       n_intervals_k = n_intervals - k; sample_no=0;
       local_close.clear(); local_dates.clear();
       tot_length = price_diff_k.size();

       date_stamp = date_diff_k.get(0);
       date_tokens = date_stamp.split(date_delims);
       prev_date = date_tokens[0];       
            
       for(i = 0; i < tot_length; i++)
       {
         date_stamp = date_diff_k.get(i);
         date_tokens = date_stamp.split(date_delims);
         
         local_close.add(price_diff_k.get(i));        
         local_dates.add(date_diff_k.get(i));
         //System.out.println(date_diff_k.get(i) + " " + price_diff_k.get(i));
         
         if(prev_date.equals(date_tokens[0])) //no date change
         {n_day_sample++;}
         else //change in day, increase count
         {count++;  n_day_sample = 0;}
         
         //System.out.println(n_day_sample + " " + n_intervals_k);
         
         if(n_day_sample == n_intervals_k-1)
         {
          if(count >= span)
          {         
            
            Double[] close = local_close.toArray(new Double[0]);
            String[] lesdates = local_dates.toArray(new String[0]);
            
            Double[] temp_close = new Double[n_intervals_k*span];
            String[] temp_dates = new String[n_intervals_k*span];
           
            
            System.arraycopy(close, n_intervals_k*sample_no, temp_close, 0, n_intervals_k*span);
            System.arraycopy(lesdates, n_intervals_k*sample_no, temp_dates, 0, n_intervals_k*span);
            
//             System.out.println("Printing sample no = " + sample_no); 
//             System.out.println();
//             
//             for(int m=0;m<temp_close.length;m++)
//             {
//               System.out.println(temp_dates[m] + " " + temp_close[m]);
//             }
                        
            delta_k = var(temp_close);          
            variance_diff_k.add(delta_k);
            sample_no++;
           }
          }
          prev_date = date_tokens[0];  
        }          
          
          
        //System.out.println("Computing variance ratio"); 
        //System.out.println("Number in var_close = " + variance_close.size() + ", number in var_k = " + variance_diff_k.size());
        
        for(i=0;i<variance_close.size();i++)
        {
          var_ratio = variance_diff_k.get(i)/(k*variance_close.get(i));
          varianceRatio.add(var_ratio);
        }        
  }


  public static void main(String args[])
  { 
     VarRatio vr = new VarRatio();
    
     File file = new File("./AUDSGD.FXCM.dat");
     
     int K = 10; //differencing 
     int span = 20;
     String begin_time = "14:00"; 
     String end_time = "23:00";
     
     vr.computeVarRatio(file, K, begin_time, end_time, span);
      
     for(int i=0;i<vr.varianceRatio.size();i++)
     {
       System.out.println(vr.varianceRatio.get(i));
     }
      
  }
  
}  
