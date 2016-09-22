# iMetrica
 A fast, interactive, GUI-oriented software suite for predictive modeling, multivariate time series analysis, real-time signal extraction, Bayesian financial econometrics, and much more.

This build is for Linux Ubuntu 64 bit. Windows 64 version of iMetrica coming soon. 

Dependencies: 
GNU Scientific Library (libgsl, libgslblas and gsl header files)  
libgfortran
libf2c

To run:
The iMetrica software project was built using the open source build automation system called Gradle (https://gradle.org/).  

In main imetrica directory type: 
./gradlew build
then 
./gradlew run


The principal use of iMetrica is to provide an interactive environment for the analysis of (multivariate) time series modeling and real-time filtering and signal extraction. The interactive features in iMetrica are meant to provide a modeling and graphics environment for users and students of econometrics, finance, and real-time data analysis that have no coding experience, as no coding is required to produce any models or filtering designs. 

The design is interactive, meaning one can change modeling data/parameter inputs and see the effects in both graphical and numerical form automatically. This feature is designed to help understand the underlying of the modeling or filtering process. One can test many attributes of the modeling or filtering process this way both visually and numerically such as sensitivity, nonlinearity, goodness-of-fit, overfitting issues, stability, etc.     

All the computational libraries were written in Fortran and/or GNU C and provided as Native libraries to the Java platform, where Java provides the user-interface, control, graphics, and several other components. All the libraries 

1) Data simulation for several econometric models
   - (S)ARIMA, (E)GARCH, (Multivariate) Stochastic Volatility, HEAVY models, Cycles/Trends, and more
   - Random number generators from several distributions to create shocks, outliers, regression    components, etc.  

2) An interactive GUI for X-13-ARIMA-SEATS called uSimX13
   - Perform automatic seasonal adjustment on thousands of economic time series
   - Compare SARIMA model choices using several different novel signal extraction diagnostics and tools available only in iMetrica
   - Visualize in real-time several components of modeling process
   - Analyze forecasts and compare with other models
   - All of the most important features of X-13-ARIMA-SEATS included  
   
3) An interactive GUI for RegComponent (State Space and Unobserved Component Models) 
4) An interactive GUI for multivariate real-time signal extraction using the multivariate direct filter approach (MDFA)
5) Empirical Mode Decomposition
6) Bayesian Time Series Modeling
      
      
Tutorials on how to use iMetrica can be found at imetricablog.com and will be added on a weekly basis. New tools, features, modules are constantly being added.  

Please send any bug reports, comments, complaints, to clisztian@gmail.com. 
      

    
  