package ch.imetrica.dataControl;



public class JAsset
{

  //--- asset information -----------
  String name;              //name of the asset
  int n_obs;                //total number of observations in price
  public double[] price;           //the price at the given frequency
  double[] volume;          //the volume at the given frequency
  public double[] log_return;
  
  int period;               //periods of frequency (ie. 1,5,30, etc.)
  String frequency;         //frequency of data (sec, min, hour, daily, etc.) 
  
  
  //--- realized volatilty
  int n_days;               //number of days 
  public double[] realized_vol;    //daily realized volat

   
  public JAsset(int n)
  {n_obs = n;}
  
  public void setName(String n) {name = n;}
  
  public void setFrequency(String f) {frequency = f;}
  
  public void setPeriod(int p) {period = p;}
  
  public void setPrice(double[] _p)
  {
    n_obs = _p.length;
    price = new double[n_obs]; 
    System.arraycopy(_p, 0, price, 0, n_obs);
  }
  
  public void setPriceVolume(double[] _p, double[] _v)
  {
    n_obs = _p.length;
    price = new double[n_obs]; volume = new double[n_obs];
    System.arraycopy(_p, 0, price, 0, n_obs); System.arraycopy(_v, 0, volume, 0, n_obs);
  }

  public void setVolume(double[] _v)
  {
    n_obs = _v.length;
    volume = new double[n_obs];
    System.arraycopy(_v, 0, volume, 0, n_obs);
  }  
  
  
  public void setLogReturn(double[] _p)
  {
    n_obs = _p.length;
    log_return = new double[n_obs]; 
    System.arraycopy(_p, 0, log_return, 0, n_obs); 
  }
  
  public void setRV(double[] rv)
  {
    n_days = rv.length;
    realized_vol = new double[n_days];
    System.arraycopy(rv, 0, realized_vol, 0, n_days);
  }
  
   
}
