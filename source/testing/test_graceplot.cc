#include "grace_plot.h"
#include "xml_handler.h"
#include "task_utils.h"

using namespace std;


void testGracePlot(XMLHandler& xml_in)
{

 if (xml_tag_count(xml_in,"TestGracePlot")==0)
 return;

 bool keep=true;
 {XMLHandler xmlr(xml_in,"TestGracePlot");
 if (xml_child_tag_count(xmlr,"PlotPrompt")>0) keep=false;}

 UserInterface *ui=new UserInterface;

 GracePlot gplot;
 gplot.setLabels("time value","m\\s\\f{0}eff\\f{}\\N(t)","Effective mass");
 gplot.addText("Label 1",0.3,0.56);
 gplot.addText("Label 2",0.4,0.7,true,2.0,"red");
 gplot.addXYDataSet("circle","open","solid","blue","Legend 1");
 gplot.addXYDataPoint(0.2,0.4);
 gplot.addXYDataPoint(0.45,0.5);
 gplot.addXYDataPoint(0.58,0.6);
 XYPoint pt(0.8,1.12);
 gplot.addXYDataPoint(pt);
 vector<XYPoint> pts(3);
 pts[0].xval=1.2; pts[1].xval=1.3; pts[2].xval=1.5;
 pts[0].yval=0.45; pts[1].yval=0.55; pts[2].yval=0.62;
 gplot.addXYDataPoints(pts);

 gplot.addXYDYDataSet("square","solid","none","red","Legend 2");
 gplot.addXYDYDataPoint(0.3, 0.3, 0.04);
 gplot.addXYDYDataPoint(0.44, 0.44, 0.12);
 gplot.addXYDYDataPoint(0.55, 0.50, 0.02);
 gplot.addXYDYDataPoint(0.66, 0.48, 0.22);

 gplot.addXYDYDYDataSet("circle","solid","none","green","Legend 3");
 gplot.addXYDYDYDataPoint(0.63, 0.53, 0.04, 0.08);
 gplot.addXYDYDYDataPoint(0.74, 0.64, 0.12, 0.24);
 gplot.addXYDYDYDataPoint(0.85, 0.60, 0.09, 0.02);
 gplot.addXYDYDYDataPoint(0.96, 0.38, 0.44, 0.22);

 gplot.addXYDXDataSet("diamond","solid","dash","cyan","Legend 4");
 gplot.addXYDXDataPoint(0.5,0.1,0.04);
 gplot.addXYDXDataPoint(0.6,0.03,0.08);
 gplot.addXYDXDataPoint(0.8,0.1,0.04);

 gplot.setFonts("helvetica-boldoblique","times-bold","helvetica-oblique","times-bolditalic");
 gplot.setLegend(0.25,0.85); 
 gplot.setLimits(0.0,0.7,0.2,0.9);
 gplot.autoScale(0.4);
 gplot.drawToScreenAndSave("crap1.agr",ui,keep);

 gplot.drawToScreen(ui,keep);
 gplot.saveToFile("crapA.agr",0.1);

 GracePlot gplot2;
 gplot2.setLabels("u","v");
 gplot2.addXYDXDXDataSet("triangleup","solid","dot","blue","Legend A");
 gplot2.addXYDXDXDataPoint(1.0, 2.0, 0.1, 0.3);
 gplot2.addXYDXDXDataPoint(2.0, 3.0, 0.5, 0.2);
 gplot2.addXYDXDXDataPoint(3.0, 5.0, 0.2, 0.4);

 gplot2.addXYDXDYDataSet("square","open","dash","red","Legend B");
 gplot2.addXYDXDYDataPoint(1.5, 6.0, 0.1, 0.3);
 gplot2.addXYDXDYDataPoint(2.5, 4.1, 0.5, 0.2);
 gplot2.addXYDXDYDataPoint(3.5, 2.3, 0.2, 0.4);

 gplot2.addXYDXDXDYDYDataSet("star","open","dashdot","orange","Legend C");
 gplot2.addXYDXDXDYDYDataPoint(2.0, 2.0, 0.1, 0.3, 0.4, 0.9);
 gplot2.addXYDXDXDYDYDataPoint(2.7, 3.1, 0.5, 0.2, 0.8, 0.1);
 gplot2.addXYDXDXDYDYDataPoint(3.9, 4.3, 0.2, 0.4, 0.7, 0.3);

 gplot2.autoScale(0.4);
 gplot2.drawToScreen(ui,keep);
 gplot2.saveToFile("crapB.agr",0.1);


 GracePlot gplot3;
 gplot3.setLabels("n","Z","Histogram");
 gplot3.addBarDataSet("cyan","red",0.6,"Legend A");
 gplot3.addBarDataPoint(1.0,1.0);
 gplot3.addBarDataPoint(2.0,1.5);
 gplot3.addBarDataPoint(3.0,2.1);
 gplot3.addBarDataPoint(4.0,0.3);
 gplot3.addBarDataPoint(5.0,1.6);
 gplot3.addBarDataPoint(6.0,2.7);
 gplot3.saveToFile("crapC.agr",0.1);

 GracePlot gplot4;
 gplot4.setLabels("n","Z","Histogram");
 gplot4.addBarDYDataSet("cyan","black",0.6,"Legend A");
 gplot4.addBarDYDataPoint(1.0,1.0,0.3);
 gplot4.addBarDYDataPoint(2.0,1.5,0.5);
 gplot4.addBarDYDataPoint(3.0,2.1,0.1);
 gplot4.addBarDYDataPoint(4.0,0.3,0.22);
 gplot4.addBarDYDataPoint(5.0,1.6,0.28);
 gplot4.addBarDYDataPoint(6.0,2.7,0.11);
 gplot4.autoScale(0.1,0.3,0.5,0.0);
 gplot4.saveToFile("crapD.agr");


 GracePlot gplot5;
 gplot5.setLabels("n","Z","Histogram");
 gplot5.addBarDYDYDataSet("green","blue",0.6,"Legend A");
 gplot5.addBarDYDYDataPoint(1.0,1.0,0.3,0.1);
 gplot5.addBarDYDYDataPoint(2.0,1.5,0.5,0.1);
 gplot5.addBarDYDYDataPoint(3.0,2.1,0.1,0.6);
 gplot5.addBarDYDYDataPoint(4.0,0.3,0.22,0.1);
 gplot5.addBarDYDYDataPoint(5.0,1.6,0.28,0.5);
 gplot5.addBarDYDYDataPoint(6.0,2.7,0.11,0.7);
 gplot5.saveToFile("crapE.agr",0.1);

 delete ui;
}

// ************************************************************

void testEffEnergy(XMLHandler& xml_in)
{

 if (xml_tag_count(xml_in,"TestEffEnergy")==0)
 return;

 XMLHandler xmlf(xml_in,"FuncDerivTest");
 uint step,Kvalue;
 double rvalue;
 xmlreadchild(xmlf,"Rvalue",rvalue,"FuncDerivTest");
 xmlreadchild(xmlf,"Kvalue",Kvalue,"FuncDerivTest");
 xmlreadchild(xmlf,"Stepvalue",step,"FuncDerivTest");

 cout.precision(16);
 cout << "PeriodicExpFuncDeriv funcd1("<<rvalue<<","<<step<<","<<Kvalue<<")"<<endl;
 PeriodicExpFuncDeriv funcd1(rvalue,step,Kvalue);
 double f,df,b;
 b=0.3; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.4; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.5; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.6; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.7; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.8; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.9; funcd1(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;

 cout <<endl<<endl;
 cout << "PeriodicExp2FuncDeriv funcd1("<<rvalue<<","<<step<<","<<Kvalue<<")"<<endl;
 PeriodicExp2FuncDeriv funcd2(rvalue,step,Kvalue);
 b=0.3; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.4; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.5; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.6; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.7; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.8; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;
 b=0.9; funcd2(b,f,df); cout <<"b="<< b<<" f="<<f<<" df="<<df<<endl;


 XMLHandler xmlc(xml_in,"Correlator");
 uint Textent;
 double A,m,B0;
 xmlreadchild(xmlc,"Avalue",A,"Correlator");
 xmlreadchild(xmlc,"Textent",Textent,"Correlator");
 xmlreadchild(xmlc,"Energy",m,"Correlator");
 xmlreadchild(xmlc,"AddConstant",B0,"Correlator");
 xmlreadchild(xmlc,"Step",step,"Correlator");

 cout <<endl<<endl;
 uint fittype=0; cout << "Fit type = "<<fittype<<endl;
 EffectiveEnergyCalculator meff0(step,Textent,fittype);
 for (unsigned int t=0;t<Textent;t++){
    double corrt=A*exp(-m*t);
    double corrtstep=A*exp(-m*(t+step));
    double msolve;
    if (meff0.calculate(msolve,t,corrt,corrtstep))
       cout << "t="<<t<<" msolve = "<<msolve<<endl;
    else
       cout << "t="<<" could not solve"<<endl;}

 cout <<endl<<endl;
 fittype=1;cout << "Fit type = "<<fittype<<endl;
 EffectiveEnergyCalculator meff1(step,Textent,fittype);
 for (unsigned int t=0;t<Textent;t++){
    double corrt=A*(exp(-m*t)+exp(-m*(Textent-t)));
    double corrtstep=A*(exp(-m*(t+step))+exp(-m*(Textent-t-step)));
    double msolve;
    if (meff1.calculate(msolve,t,corrt,corrtstep))
       cout << "t="<<t<<" msolve = "<<msolve<<endl;
    else
       cout << "t="<<" could not solve"<<endl;}

 fittype=2;cout << "Fit type = "<<fittype<<endl;
 EffectiveEnergyCalculator meff2(step,Textent,fittype);
 for (unsigned int t=0;t<Textent;t++){
    double corrt=A*exp(-m*t)+B0;
    double corrtstep=A*exp(-m*(t+step))+B0;
    double corrtback=A*exp(-m*(t-step))+B0;
    double msolve;
    if (meff2.calculate(msolve,t,corrt,corrtstep,corrtback))
       cout << "t="<<t<<" msolve = "<<msolve<<endl;
    else
       cout << "t="<<" could not solve"<<endl;}

 cout <<endl<<endl;
 fittype=3;cout << "Fit type = "<<fittype<<endl;
 EffectiveEnergyCalculator meff3(step,Textent,fittype);
 for (unsigned int t=0;t<Textent;t++){
    double corrt=A*(exp(-m*t)+exp(-m*(Textent-t)))+B0;
    double corrtstep=A*(exp(-m*(t+step))+exp(-m*(Textent-t-step)))+B0;
    double corrtback=A*exp(-m*(t-step))+B0;
    double msolve;
    if (meff3.calculate(msolve,t,corrt,corrtstep,corrtback))
       cout << "t="<<t<<" msolve = "<<msolve<<endl;
    else
       cout << "t="<<" could not solve"<<endl;}

}

// ************************************************************
