package ch.imetrica.mdfaTradingStrategies;

import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.text.*;
import java.awt.geom.Ellipse2D;




public class MDFAStrategyTradeCanvas extends JPanel
{


   /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
   MDFATrade[] trades;
   HashMap<Integer,String> trade_info;
   DecimalFormat df3;
   int height, width;
    
   boolean keep_fixed = false;
   NumberFormat formatter;
   double take_profit; 
   Graphics2D g2d; 
   Ellipse2D ellipse; 
   double max_profloss;
   double max_dd;
   String info_bitch; 
   boolean data_set = false;
   boolean trades_set = false;
   Font mono;
   int max_x,max_y;
   BasicStroke dashed,orig;
   //--- conditional means------
   double avg_max_pl_loss;
   double avg_min_pl_loss;
   double avg_max_pl_win;
   double avg_min_pl_win;

   double avg_profit,avg_loss;
   double win_ratio,kelly_crit;
   double risk_reward;
   int total_trades;
   int n_wins;
   boolean paint_sl_line = true;
   double stop_loss_line = 0;
   float[] dash1;
   
   
    public MDFAStrategyTradeCanvas(int w, int h)
    {
      // Initilize everything-------------------
    
      
      this.width = w; this.height = h; 
      mono  = new Font("Monospaced", Font.PLAIN, 12);
      setBackground(Color.BLACK);      
      setPreferredSize(new Dimension(w, h));
      
      formatter = new DecimalFormat("#0.00000"); 
      trade_info = new HashMap<Integer,String>();
      data_set = false;
      info_bitch = new String(" ");
      trades_set = false;
      take_profit = 0;
      keep_fixed = false;
      
      addMouseMotionListener(new MyArtMouseMotionListener());  
      addMouseListener(new MyArtMouseListener());   
      
      avg_max_pl_loss=0; 
      avg_min_pl_loss=0;
      avg_max_pl_win=0;
      avg_min_pl_win=0;

      win_ratio = 0;
      risk_reward = 0;
      total_trades = 0;
      n_wins = 0;
      avg_profit = 0;
      avg_loss = 0;
      
      dash1 = new float[1]; dash1[0] = 10.0f;
      dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER, 
                                          10.0f, dash1, 0.0f);  
      
    }
    
    
    public void setTrades(MDFATrade[] t)
    {
      trades = new MDFATrade[t.length];
      System.arraycopy(t,0,trades,0,t.length);
      trades_set = true;
      setDataNormalization();
    }
    
    public void setTakeProfit(double t) {take_profit = t; go();}
    public void keepFrameFixed(boolean t) {keep_fixed = t;}
    
    //----- for normalization of data canvas assumes .00001 increments-------------
    public void setDataNormalization()
    {
      int i;
      int len = trades.length;
      double log_pl,log_dd;
      int hashmarkpl,hashmarkdd;
      avg_max_pl_loss=0; 
      avg_min_pl_loss=0;
      avg_max_pl_win=0;
      avg_min_pl_win=0;
      win_ratio = 0;
      risk_reward = 0;
      total_trades = 0;
      n_wins = 0;      
      avg_profit = 0;
      avg_loss = 0;      
      kelly_crit = 0;
      
      Dimension ds = this.getSize();
      width = ds.width; height = ds.height;      
      
      if(!keep_fixed)
      {
       max_profloss = 0;
       max_dd = 0;
      }
      
      
      total_trades = len;
      
      //--- find max_min and compute trade statistics
      for(i = 0; i < len; i++)
      {
      
        if(!keep_fixed)
        {
         if(max_dd < Math.abs(trades[i].drawdown)) {max_dd = Math.abs(trades[i].drawdown);}
         if(max_profloss < Math.abs(trades[i].pl)) {max_profloss = Math.abs(trades[i].pl);}
        }
        
        if(trades[i].pl > 0) //profit
        {
          avg_max_pl_win = avg_max_pl_win + trades[i].drawup;
          avg_min_pl_win = avg_min_pl_win - trades[i].drawdown;
          avg_profit = avg_profit + trades[i].pl;
          n_wins++;
        }
        else
        {
          avg_max_pl_loss = avg_max_pl_loss + trades[i].drawup;
          avg_min_pl_loss = avg_min_pl_loss - trades[i].drawdown;
          avg_loss = avg_loss - trades[i].pl;
        }
      
      }
 
      avg_profit = avg_profit/n_wins;
      avg_loss = avg_loss/(len - n_wins);
      
      risk_reward = avg_loss/avg_profit;
      
      avg_max_pl_win = avg_max_pl_win/n_wins;
      avg_min_pl_win = avg_min_pl_win/n_wins;
      
      avg_max_pl_loss = avg_max_pl_loss/(len - n_wins);
      avg_min_pl_loss = avg_min_pl_loss/(len - n_wins);
      win_ratio = (double)n_wins/len;
 
      kelly_crit = win_ratio - (1.0 - win_ratio)*risk_reward;
 
      if(!keep_fixed)
      {
       stop_loss_line = max_dd;
       max_x = (int)(100*Math.log(100000*max_dd));
       max_y = (int)(100*Math.log(100000*max_profloss)); 
      }
      
      //second round, after finding the extremes 
      for(i = 0; i < len; i++)
      {
      
        log_dd = 0;
        if(Math.abs(trades[i].drawdown) > 0)
        {
          log_dd = Math.log(100000.0*Math.abs(trades[i].drawdown));
        }
        
        log_pl = 0;        
        if(Math.abs(trades[i].pl) > 0)
        {
          log_pl = Math.log(100000.0*Math.abs(trades[i].pl));
        }
       
       
        hashmarkpl = (int)(100*log_pl);
        hashmarkdd = (int)(100*log_dd);
        
        pips2canvas(hashmarkdd,hashmarkpl);
        
        //System.out.println(hashmark + " " + trades[i].startTransTime + " " + trades[i].endTransTime + " " + trades[i].drawdown + " " + trades[i].pl);
        //trade_info.put(new Integer(hashmark), new String(trades[i].startTransTime + " " + trades[i].endTransTime + " " + trades[i].drawdown + " " + trades[i].pl));
        
        
      }
      data_set = true;
      
    }
    
    // maps the pips in log form to current screen dimensions
    public int pips2canvas(int x, int y)
    {
     
     if(max_x > 0 && max_y > 0)
     {
    
      int wi = (int)(((double)x/max_x)*width);      
      int hi = (int)(((double)y/max_y)*height);
     
      return (1000*wi + hi); 
     }
     else
     {return 0;}
     
    }
    
    
  public void setStopLossLine(double v) {stop_loss_line = v; go();}  
    
  public void go() {repaint();}

  public void paintComponent(Graphics g)
  {  
   if(trades_set)
   {
     int i;
     double log_pl,log_dd;
     
     double pips_pl;
     double pips_dd;
     double pips_max;
     
     new BasicStroke((float)1.7); 
     Stroke thin = new BasicStroke((float)1.2);   
     if(!data_set)
     {setDataNormalization();}
     //setDataNormalization();
          
     super.paintComponent(g);
     g2d = (Graphics2D)g; 
   
     Dimension ds = this.getSize();
     width = ds.width; height = ds.height;
    
    
    
     //--- paint the canvas with dots
     for(i = 0; i < trades.length; i++)
     {
       // --- get the info-----
       pips_pl = trades[i].pl;
       pips_dd = trades[i].drawdown;
       pips_max = trades[i].drawup;
      
      
       // paint trade either green or red
       if(pips_pl > 0)
       {g2d.setPaint(Color.green);}
       else 
       {
        g2d.setPaint(Color.red);
        
        if(pips_max > take_profit)
        {g2d.setPaint(Color.yellow);} 
       }
       
       
       
        log_dd = 0;
        if(Math.abs(pips_dd) > 0)
        {
          log_dd = Math.log(100000.0*Math.abs(pips_dd));
        }
        
        log_pl = 0;        
        if(Math.abs(pips_pl) > 0)
        {
          log_pl = Math.log(100000.0*Math.abs(pips_pl));
        }
       

        //1000*hashmarkdd + hashmarkpl;
       
        int wi = (int)(((100*log_dd)/max_x)*width);      
        int hi = (int)(((100*log_pl)/max_y)*height);
       
       
        //where do I paint the trade 
        ellipse = new Ellipse2D.Double(wi, (height) - hi, 3, 3);
        g2d.draw(ellipse);  
        g2d.fill(ellipse);         
        
        if(!trade_info.containsKey(new Integer(1000*wi + ((height) - hi))))
        {
         trade_info.put(new Integer(1000*wi + ((height) - hi)), new String(trades[i].date + " " + trades[i].startTransTime + " " + trades[i].endTransTime + " " + formatter.format(trades[i].drawdown) + " " + 
          formatter.format(trades[i].drawup) + " " + formatter.format(trades[i].pl)));
        }        
        
      }
     
      // if touching a dot, put info up in left corner  
      g.setFont(mono);
      g2d.setPaint(Color.magenta);    
      g.drawString(info_bitch, 5, 15);       
        
        
      //put trade statistics in bottom left  
      g2d.setPaint(Color.green);  
      g.drawString("Avg unrlzd profit: " + formatter.format(avg_max_pl_win), 5, height - 15*5);
      g.drawString("Avg unrlzd loss  : " + formatter.format(avg_min_pl_win), 5, height - 15*4);

      g2d.setPaint(Color.red);  
      g.drawString("Avg unrlzd profit: " + formatter.format(avg_max_pl_loss), 5, height - 15*3);
      g.drawString("Avg unrlzd loss  : " + formatter.format(avg_min_pl_loss), 5, height - 15*2);
      
      g2d.setPaint(Color.magenta); 
      g.drawString("Avg profit/loss,win_ratio,risk/reward,kelly: " + formatter.format(avg_profit) + ", " + formatter.format(avg_loss) + ", " + formatter.format(win_ratio) + ", " + formatter.format(risk_reward) + ", " + formatter.format(kelly_crit), 5, height - 15);  
    
    
      g2d.setPaint(Color.gray);
      g.drawString(""+formatter.format(max_dd), width - 60, height - 15*2);
      g.drawString(""+formatter.format(max_profloss), 5, 30);
      
      
      if(paint_sl_line)
      {
        double sl = Math.log(100000.0*Math.abs(stop_loss_line));
        int wi = (int)(((100*sl)/max_x)*width);
        g2d.setStroke(thin); 
        g2d.setStroke(dashed);
        g2d.setPaint(new Color(100,140,170));
        g2d.drawLine(wi, 0, wi, height);
      }
      
      
    }
  }
    
    
     class MyArtMouseMotionListener extends MouseMotionAdapter 
     {
      
      public void mouseDragged(MouseEvent e) 
      { }

      public void mouseMoved(MouseEvent e) 
      {

         //j = (int)(((double)(N))*e.getX()/(double)width);  
         info_bitch = new String(" ");
         if(trade_info.containsKey(new Integer(1000*e.getX() + e.getY())))
         {info_bitch = trade_info.get(1000*e.getX() + e.getY());}

         //System.out.println(1000*e.getX() + e.getY());
         go();           
      }    
       
     }
   
     
    class MyArtMouseListener extends MouseAdapter 
    {       
       public void mousePressed(MouseEvent e) 
       {}

       public void mouseReleased(MouseEvent e) { }
       public void mouseClicked(MouseEvent e) { }
        
    }      
 
 
} 
