package ch.imetrica.mdfaTradingStrategies;



  public class MDFATrade
  {
   int startTransTime; //which time trade started
   int endTransTime;   //which time trade ended
 
   String date;
   double drawdown;    //lowest drawdown
   double drawup;      //highest drawup
   double pl;          //profit > 0, or loss < 0

   public MDFATrade(String d, int start, int end, double dd, double ma, double p)
   {
     date = new String(d);
     startTransTime = start; 
     endTransTime = end;
     drawdown = dd; 
     drawup = ma;
     pl = p;
   }
   
   public void printTrade()
   {
     System.out.println(date + " " + startTransTime + " " + endTransTime + " " + drawdown + " " + drawup + " " + pl);
   }
   
   
  }