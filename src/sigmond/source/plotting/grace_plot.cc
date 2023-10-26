#include "grace_plot.h"
#include <cstring>
#include <cstdlib>
#include <grace_np.h>
#include <limits> 
#include <fstream>
#include <stdexcept>

using namespace std;

// *************************************************************


GracePlot::GracePlot() : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                         m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
}

GracePlot::GracePlot(double xmin, double xmax, double ymin, double ymax) 
               : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                 m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
 setLimits(xmin,xmax,ymin,ymax);
}

GracePlot::GracePlot(double xmin, double xmax, double ymin, double ymax,
                     const string& xlabel, const string& ylabel)
               : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                 m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
 setLimits(xmin,xmax,ymin,ymax);
 setLabels(xlabel,ylabel);
}

GracePlot::GracePlot(double xmin, double xmax, double ymin, double ymax,
                     const string& xlabel, const string& ylabel,
                     const string& title)
               : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                 m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
 setLimits(xmin,xmax,ymin,ymax);
 setLabels(xlabel,ylabel,title);
}

GracePlot::GracePlot(const string& xlabel, const string& ylabel)
               : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                 m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
 setLabels(xlabel,ylabel);
}

GracePlot::GracePlot(const string& xlabel, const string& ylabel,
                     const string& title)
               : m_dparams(dparamsize), m_curr_dataset_index(-1), 
                 m_curr_dataset_type(-1), m_dcount(0)
{
 setDefaults();
 setLabels(xlabel,ylabel,title);
}




GracePlot::GracePlot(const GracePlot& gplot) 
    :  m_dparams(gplot.m_dparams), m_title(gplot.m_title), m_xlabel(gplot.m_xlabel),
       m_ylabel(gplot.m_ylabel), m_curr_dataset_index(gplot.m_curr_dataset_index),
       m_curr_dataset_type(gplot.m_curr_dataset_type), m_data(gplot.m_data),
       m_text(gplot.m_text), m_dcount(gplot.m_dcount), m_dset(gplot.m_dset)

{}

GracePlot::GracePlot(GracePlot& gplot)
    :  m_dparams(gplot.m_dparams), m_title(gplot.m_title), m_xlabel(gplot.m_xlabel),
       m_ylabel(gplot.m_ylabel), m_curr_dataset_index(gplot.m_curr_dataset_index),
       m_curr_dataset_type(gplot.m_curr_dataset_type), m_data(gplot.m_data),
       m_text(gplot.m_text), m_dcount(gplot.m_dcount), m_dset(gplot.m_dset)
{}

GracePlot& GracePlot::operator=(const GracePlot& gplot)
{
 m_dparams=gplot.m_dparams;
 m_title=gplot.m_title;
 m_xlabel=gplot.m_xlabel;
 m_ylabel=gplot.m_ylabel;
 m_curr_dataset_index=gplot.m_curr_dataset_index;
 m_curr_dataset_type=gplot.m_curr_dataset_type;
 m_data=gplot.m_data;
 m_text=gplot.m_text;
 m_dcount=gplot.m_dcount;
 m_dset=gplot.m_dset;
 return *this;
}

GracePlot& GracePlot::operator=(GracePlot& gplot)
{
 m_dparams=gplot.m_dparams;
 m_title=gplot.m_title;
 m_xlabel=gplot.m_xlabel;
 m_ylabel=gplot.m_ylabel;
 m_curr_dataset_index=gplot.m_curr_dataset_index;
 m_curr_dataset_type=gplot.m_curr_dataset_type;
 m_data=gplot.m_data;
 m_text=gplot.m_text;
 m_dcount=gplot.m_dcount;
 m_dset=gplot.m_dset;
 return *this;
}


GracePlot::~GracePlot()
{
 if (GraceIsOpen()) GraceClose();
}



void GracePlot::setDefaults()
{
 //setLabels("x","y","Plot Title");
 setFontsizes(1.5,2.0,2.0,1.0);
 setFonts("times-italics","times-italics","times-roman","times-roman");
 setView(0.15,0.95,0.15,0.85);
 setLimits(0.0,1.0,0.0,1.0);
 m_dparams[i_legend_flag]=-1.0;
 m_dparams[i_val_xmin]=std::numeric_limits<double>::infinity();
 m_dparams[i_val_xmax]=-std::numeric_limits<double>::infinity();
 m_dparams[i_val_ymin]=std::numeric_limits<double>::infinity();
 m_dparams[i_val_ymax]=-std::numeric_limits<double>::infinity();
 m_dcount=0;
 m_dset="s"+make_string(0);
}


void GracePlot::setView(double view_xmin, double view_xmax, double view_ymin, double view_ymax)
{
 m_dparams[i_view_xmin]=view_xmin;
 m_dparams[i_view_xmax]=view_xmax;
 m_dparams[i_view_ymin]=view_ymin;
 m_dparams[i_view_ymax]=view_ymax;
}

void GracePlot::setLabels(const string& xlabel, const string& ylabel)
{
 m_xlabel=xlabel;
 m_ylabel=ylabel;
}

void GracePlot::setLabels(const string& xlabel, const string& ylabel,
                          const string& title)
{
 m_xlabel=xlabel;
 m_ylabel=ylabel;
 m_title=title;
}

void GracePlot::setFontsizes(double ticklabelcharsize, double axislabelcharsize, 
                             double titlecharsize, double textcharsize)
{
 m_dparams[i_tick_label_char_size]=ticklabelcharsize;
 m_dparams[i_axis_label_char_size]=axislabelcharsize;
 m_dparams[i_title_char_size]=titlecharsize;
 m_dparams[i_text_char_size]=textcharsize;
}

void GracePlot::setFonts(const std::string& xlabelfont, const std::string& ylabelfont,
                         const std::string& titlefont, const std::string& legendfont)
{
 m_dparams[i_xlabel_font]=encode_font(xlabelfont);
 m_dparams[i_ylabel_font]=encode_font(ylabelfont);
 m_dparams[i_title_font]=encode_font(titlefont);
 m_dparams[i_legend_font]=encode_font(legendfont);
}


void GracePlot::setLimits(double xmin, double xmax, double ymin, double ymax)
{
 m_dparams[i_xmin]=xmin;
 m_dparams[i_xmax]=xmax;
 m_dparams[i_ymin]=ymin;
 m_dparams[i_ymax]=ymax;
}


void GracePlot::setVerticalLimits(double ymin, double ymax)
{
 m_dparams[i_ymin]=ymin;
 m_dparams[i_ymax]=ymax;
}


void GracePlot::setHorizontalLimits(double xmin, double xmax)
{
 m_dparams[i_xmin]=xmin;
 m_dparams[i_xmax]=xmax;
}


void GracePlot::setLegend(double xviewpos, double yviewpos, bool viewport)
{
 m_dparams[i_legend_flag]=1.0;
 m_dparams[i_legend_xpos]=xviewpos;
 m_dparams[i_legend_ypos]=yviewpos;
 m_dparams[i_legend_loctype]=(viewport)?1.0:-1.0;
}


string GracePlot::make_string(int data) 
{
 try{
    ostringstream oss;
    oss.setf(ios_base::boolalpha);
    oss.precision(12);
    oss << data;
    return oss.str();}
 catch(const std::exception& msg){ 
    ostringstream err;
    err << "Error: Failed to convert numerical to string: data = "<<data;
    throw(std::invalid_argument(err.str()));}
}

string GracePlot::make_string(double data) 
{
 try{
    ostringstream oss;
    oss.setf(ios_base::boolalpha);
    oss.precision(12);
    oss << data;
    return oss.str();}
 catch(const std::exception& msg){ 
    ostringstream err;
    err << "Error: Failed to convert numerical to string: data = "<<data;
    throw(std::invalid_argument(err.str()));}
}



  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing blanks.

string GracePlot::tidy_string(const string& str)   
{
 string tmp;
 for (size_t i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 size_t start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 size_t len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}


void GracePlot::addXYDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(0,false,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDataPoint(const XYPoint& pt)
{
 if (m_curr_dataset_type!=0){
    throw(std::invalid_argument("Current data set is not XY: could not add XY point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDataPoint(double x, double y)
{
 addXYDataPoint(XYPoint(x,y));
}

void GracePlot::addXYDataPoints(const vector<XYPoint>& pts)
{
 if (m_curr_dataset_type!=0){
    throw(std::invalid_argument("Current data set is not XY: could not add XY point"));}
 string prefix="g0."+m_dset;
 for (vector<XYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}




void GracePlot::addXYDYDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(1,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDYDataPoint(const XYDYPoint& pt)
{
 if (m_curr_dataset_type!=1){
    throw(std::invalid_argument("Current data set is not XYDY: could not add XYDY point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDYDataPoint(double x, double y, double dy)
{
 addXYDYDataPoint(XYDYPoint(x,y,dy));
}

void GracePlot::addXYDYDataPoints(const vector<XYDYPoint>& pts)
{
 if (m_curr_dataset_type!=1){
    throw(std::invalid_argument("Current data set is not XYDY: could not add XYDY point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}





void GracePlot::addXYDXDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(2,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDXDataPoint(const XYDXPoint& pt)
{
 if (m_curr_dataset_type!=2){
    throw(std::invalid_argument("Current data set is not XYDX: could not add XYDX point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDXDataPoint(double x, double y, double dx)
{
 addXYDXDataPoint(XYDXPoint(x,y,dx));
}

void GracePlot::addXYDXDataPoints(const vector<XYDXPoint>& pts)
{
 if (m_curr_dataset_type!=2){
    throw(std::invalid_argument("Current data set is not XYDX: could not add XYDX point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDXPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}





void GracePlot::addXYDXDYDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(3,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDXDYDataPoint(const XYDXDYPoint& pt)
{
 if (m_curr_dataset_type!=3){
    throw(std::invalid_argument("Current data set is not XYDXDY: could not add XYDXDY point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDXDYDataPoint(double x, double y, double dx, double dy)
{
 addXYDXDYDataPoint(XYDXDYPoint(x,y,dx,dy));
}

void GracePlot::addXYDXDYDataPoints(const vector<XYDXDYPoint>& pts)
{
 if (m_curr_dataset_type!=3){
    throw(std::invalid_argument("Current data set is not XYDXDY: could not add XYDXDY point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDXDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}




void GracePlot::addXYDXDXDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(4,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDXDXDataPoint(const XYDXDXPoint& pt)
{
 if (m_curr_dataset_type!=4){
    throw(std::invalid_argument("Current data set is not XYDXDX: could not add XYDXDX point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDXDXDataPoint(double x, double y, double dxup, double dxdn)
{
 addXYDXDXDataPoint(XYDXDXPoint(x,y,dxup,dxdn));
}

void GracePlot::addXYDXDXDataPoints(const vector<XYDXDXPoint>& pts)
{
 if (m_curr_dataset_type!=4){
    throw(std::invalid_argument("Current data set is not XYDXDX: could not add XYDXDX point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDXDXPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}





void GracePlot::addXYDYDYDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(5,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDYDYDataPoint(const XYDYDYPoint& pt)
{
 if (m_curr_dataset_type!=5){
    throw(std::invalid_argument("Current data set is not XYDYDY: could not add XYDYDY point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDYDYDataPoint(double x, double y, double dyup, double dydn)
{
 addXYDYDYDataPoint(XYDYDYPoint(x,y,dyup,dydn));
}

void GracePlot::addXYDYDYDataPoints(const vector<XYDYDYPoint>& pts)
{
 if (m_curr_dataset_type!=5){
    throw(std::invalid_argument("Current data set is not XYDYDY: could not add XYDYDY point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDYDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}





void GracePlot::addXYDXDXDYDYDataSet(const string& symbol, 
             const string& symbolfill, const string& linestyle, 
             const string& color, const string& legendtext)
{
 addDataSet(6,true,symbol,symbolfill,linestyle,color,legendtext);
}

void GracePlot::addXYDXDXDYDYDataPoint(const XYDXDXDYDYPoint& pt)
{
 if (m_curr_dataset_type!=6){
    throw(std::invalid_argument("Current data set is not XYDXDXDYDY: could not add XYDXDXDYDY point"));}
 addPoint(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addXYDXDXDYDYDataPoint(double x, double y, double dxup, double dxdn,
                                       double dyup, double dydn)
{
 addXYDXDXDYDYDataPoint(XYDXDXDYDYPoint(x,y,dxup,dxdn,dyup,dydn));
}

void GracePlot::addXYDXDXDYDYDataPoints(const vector<XYDXDXDYDYPoint>& pts)
{
 if (m_curr_dataset_type!=6){
    throw(std::invalid_argument("Current data set is not XYDXDXDYDY: could not add XYDXDXDYDY point"));}
 string prefix="g0."+m_dset;
 for (vector<XYDXDXDYDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addPoint(prefix, *it, m_dcount);}
}





void GracePlot::addDataSet(int dataset_type, bool error_bars, const string& symbol, 
                           const string& symbolfill, const string& linestyle, 
                           const string& color, const string& legendtext)
{
 m_curr_dataset_type=dataset_type;
 m_curr_dataset_index++;
 m_dset="s"+make_string(m_curr_dataset_index);
 m_dcount=0;
 string colcode(encode_color(color));
 string cmd;
 cmd=m_dset+" on"; m_data.push_back(cmd);
 cmd=m_dset+" type "+encode_datatype(dataset_type); m_data.push_back(cmd);
 cmd=m_dset+" symbol "+encode_symbol(symbol); m_data.push_back(cmd);
 cmd=m_dset+" symbol size 0.7"; m_data.push_back(cmd);
 cmd=m_dset+" symbol color "+colcode; m_data.push_back(cmd);
 cmd=m_dset+" symbol fill color "+colcode; m_data.push_back(cmd);
 cmd=m_dset+" symbol fill pattern "+encode_symbolfill(symbolfill); m_data.push_back(cmd);
 cmd=m_dset+" symbol linewidth 2.0"; m_data.push_back(cmd);
 cmd=m_dset+" line type 1"; m_data.push_back(cmd);
 cmd=m_dset+" line linestyle "+encode_linestyle(linestyle); m_data.push_back(cmd);
 cmd=m_dset+" line linewidth 2.0"; m_data.push_back(cmd);
 cmd=m_dset+" line color "+colcode; m_data.push_back(cmd);
 if (error_bars){
    cmd=m_dset+" errorbar on"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar place both"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar color "+colcode; m_data.push_back(cmd);
    cmd=m_dset+" errorbar pattern 1"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar size 1.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar linewidth 2.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar linestyle 1"  ; m_data.push_back(cmd);  
    cmd=m_dset+" errorbar riser linewidth 2.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser linestyle 1"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser clip off"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser clip length 0.10"; m_data.push_back(cmd);}

 if (!legendtext.empty()){
    cmd=m_dset+" legend \""+legendtext+"\""; m_data.push_back(cmd);}
}


void GracePlot::addPoint(const string& prefix, const XYPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_horizontal_span(pt.xval);
 update_vertical_span(pt.yval);
 count++;
}


void GracePlot::addPoint(const string& prefix, const XYDYPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_horizontal_span(pt.xval);
 update_vertical_span(pt.yval+pt.yerr);
 update_vertical_span(pt.yval-pt.yerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.yerr);
 m_data.push_back(cmd);
 count++;
}




void GracePlot::addPoint(const string& prefix, const XYDXPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_vertical_span(pt.yval);
 update_horizontal_span(pt.xval+pt.xerr);
 update_horizontal_span(pt.xval-pt.xerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.xerr);
 m_data.push_back(cmd);
 count++;
}

void GracePlot::addPoint(const string& prefix, const XYDXDXPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_vertical_span(pt.yval);
 update_horizontal_span(pt.xval+pt.xuperr);
 update_horizontal_span(pt.xval-pt.xdnerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.xuperr);
 m_data.push_back(cmd);
 cmd=prefix+".y2["+make_string(count)+"]="+make_string(pt.xdnerr);
 m_data.push_back(cmd);
 count++;
}

void GracePlot::addPoint(const string& prefix, const XYDXDYPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_horizontal_span(pt.xval+pt.xerr);
 update_horizontal_span(pt.xval-pt.xerr);
 update_vertical_span(pt.yval+pt.yerr);
 update_vertical_span(pt.yval-pt.yerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.xerr);
 m_data.push_back(cmd);
 cmd=prefix+".y2["+make_string(count)+"]="+make_string(pt.yerr);
 m_data.push_back(cmd);
 count++;
}

void GracePlot::addPoint(const string& prefix, const XYDYDYPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_horizontal_span(pt.xval);
 update_vertical_span(pt.yval+pt.yuperr);
 update_vertical_span(pt.yval-pt.ydnerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.yuperr);
 m_data.push_back(cmd);
 cmd=prefix+".y2["+make_string(count)+"]="+make_string(pt.ydnerr);
 m_data.push_back(cmd);
 count++;
}

void GracePlot::addPoint(const string& prefix, const XYDXDXDYDYPoint& pt, int& count)
{
 string cmd=prefix+" point "+make_string(pt.xval)+", "+make_string(pt.yval);
 m_data.push_back(cmd);
 update_horizontal_span(pt.xval+pt.xuperr);
 update_horizontal_span(pt.xval-pt.xdnerr);
 update_vertical_span(pt.yval+pt.yuperr);
 update_vertical_span(pt.yval-pt.ydnerr);
 cmd=prefix+".y1["+make_string(count)+"]="+make_string(pt.xuperr);
 m_data.push_back(cmd);
 cmd=prefix+".y2["+make_string(count)+"]="+make_string(pt.xdnerr);
 m_data.push_back(cmd);
 cmd=prefix+".y3["+make_string(count)+"]="+make_string(pt.yuperr);
 m_data.push_back(cmd);
 cmd=prefix+".y4["+make_string(count)+"]="+make_string(pt.ydnerr);
 m_data.push_back(cmd);
 count++;
}



void GracePlot::addBarDataSet(const string& color, const string& bordercolor,
                              double barwidth, const string& legendtext)
{
 addBarSet(7,false,color,bordercolor,barwidth,legendtext);
}

void GracePlot::addBarDataPoint(const XYPoint& pt)
{
 if (m_curr_dataset_type!=7){
    throw(std::invalid_argument("Current data set is not BAR: could not add BAR data"));}
 addBar(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addBarDataPoint(double x, double y)
{
 addBarDataPoint(XYPoint(x,y));
}

void GracePlot::addBarDataPoints(const vector<XYPoint>& pts)
{
 if (m_curr_dataset_type!=7){
    throw(std::invalid_argument("Current data set is not BAR: could not add BAR data"));}
 string prefix="g0."+m_dset;
 for (vector<XYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addBar(prefix, *it, m_dcount);}
}



void GracePlot::addBarDYDataSet(const string& color, const string& bordercolor,
                                double barwidth, const string& legendtext)
{
 addBarSet(8,true,color,bordercolor,barwidth,legendtext);
}

void GracePlot::addBarDYDataPoint(const XYDYPoint& pt)
{
 if (m_curr_dataset_type!=8){
    throw(std::invalid_argument("Current data set is not BARDY: could not add BARDY data"));}
 addBar(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addBarDYDataPoint(double x, double y, double yerr)
{
 addBarDYDataPoint(XYDYPoint(x,y,yerr));
}

void GracePlot::addBarDYDataPoints(const vector<XYDYPoint>& pts)
{
 if (m_curr_dataset_type!=8){
    throw(std::invalid_argument("Current data set is not BARDY: could not add BARDY data"));}
 string prefix="g0."+m_dset;
 for (vector<XYDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addBar(prefix, *it, m_dcount);}
}


void GracePlot::addBarDYDYDataSet(const string& color, const string& bordercolor,
                                  double barwidth, const string& legendtext)
{
 addBarSet(9,true,color,bordercolor,barwidth,legendtext);
}

void GracePlot::addBarDYDYDataPoint(const XYDYDYPoint& pt)
{
 if (m_curr_dataset_type!=9){
    throw(std::invalid_argument("Current data set is not BARDYDY: could not add BARDYDY data"));}
 addBar(string("g0.")+m_dset, pt, m_dcount);
}

void GracePlot::addBarDYDYDataPoint(double x, double y, double yuperr, double ydnerr)
{
 addBarDYDYDataPoint(XYDYDYPoint(x,y,yuperr,ydnerr));
}

void GracePlot::addBarDYDYDataPoints(const vector<XYDYDYPoint>& pts)
{
 if (m_curr_dataset_type!=9){
    throw(std::invalid_argument("Current data set is not BARDYDY: could not add BARDYDY data"));}
 string prefix="g0."+m_dset;
 for (vector<XYDYDYPoint>::const_iterator it=pts.begin();it!=pts.end();++it){
    addBar(prefix, *it, m_dcount);}
}





void GracePlot::addBarSet(int dataset_type, bool error_bars, const string& color, 
                          const string& bordercolor, double barwidth, 
                          const string& legendtext)
{
 m_dparams[i_temp]=(barwidth>0)?barwidth:0.0;
 m_curr_dataset_type=dataset_type;
 m_curr_dataset_index++;
 m_dset="s"+make_string(m_curr_dataset_index);
 m_dcount=0;
 string colcode(encode_color(color));
 string bcolcode(encode_color(bordercolor));
 string cmd;
 cmd=m_dset+" on"; m_data.push_back(cmd);
 cmd=m_dset+" type "+encode_datatype(dataset_type); m_data.push_back(cmd);
 cmd=m_dset+" symbol "+encode_symbol("none"); m_data.push_back(cmd);
 cmd=m_dset+" line type 1"; m_data.push_back(cmd);
 cmd=m_dset+" line linestyle "+encode_linestyle("solid"); m_data.push_back(cmd);
 cmd=m_dset+" line linewidth 2.0"; m_data.push_back(cmd);
 cmd=m_dset+" line color "+bcolcode; m_data.push_back(cmd);
 cmd=m_dset+" fill type 2"; m_data.push_back(cmd);
 cmd=m_dset+" fill rule 0"; m_data.push_back(cmd);
 cmd=m_dset+" fill pattern 1"; m_data.push_back(cmd);
 cmd=m_dset+" fill color "+colcode; m_data.push_back(cmd);
 if (error_bars){
    cmd=m_dset+" errorbar on"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar place both"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar color "+bcolcode; m_data.push_back(cmd);
    cmd=m_dset+" errorbar pattern 1"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar size 1.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar linewidth 2.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar linestyle 1"  ; m_data.push_back(cmd);  
    cmd=m_dset+" errorbar riser linewidth 2.0"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser linestyle 1"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser clip off"; m_data.push_back(cmd);
    cmd=m_dset+" errorbar riser clip length 0.10"; m_data.push_back(cmd);}

 if (!legendtext.empty()){
    cmd=m_dset+" legend \""+legendtext+"\""; m_data.push_back(cmd);}
}


void GracePlot::addBar(const string& prefix, const XYPoint& pt, int& count)
{
 double barhalfwidth=0.5*m_dparams[i_temp];
 XYPoint xy(pt.xval-barhalfwidth,0.0);
 addPoint(prefix,xy,count);
 xy.yval=pt.yval;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval+barhalfwidth;
 addPoint(prefix,xy,count);
 xy.yval=0.0;
 addPoint(prefix,xy,count);
}

void GracePlot::addBar(const string& prefix, const XYDYPoint& pt, int& count)
{
 double barhalfwidth=0.5*m_dparams[i_temp];
 XYDYPoint xy(pt.xval-barhalfwidth,0.0,0.0);
 addPoint(prefix,xy,count);
 xy.yval=pt.yval;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval; xy.yerr=pt.yerr;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval+barhalfwidth; xy.yerr=0.0;
 addPoint(prefix,xy,count);
 xy.yval=0.0;
 addPoint(prefix,xy,count);
}

void GracePlot::addBar(const string& prefix, const XYDYDYPoint& pt, int& count)
{
 double barhalfwidth=0.5*m_dparams[i_temp];
 XYDYDYPoint xy(pt.xval-barhalfwidth,0.0,0.0,0.0);
 addPoint(prefix,xy,count);
 xy.yval=pt.yval;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval; xy.yuperr=pt.yuperr; xy.ydnerr=pt.ydnerr;
 addPoint(prefix,xy,count);
 xy.xval=pt.xval+barhalfwidth; xy.yuperr=xy.ydnerr=0.0;
 addPoint(prefix,xy,count);
 xy.yval=0.0;
 addPoint(prefix,xy,count);
}









void GracePlot::renderCommands(list<string>& commands)
{
 commands.clear();
 string labelcharsize(make_string(m_dparams[i_axis_label_char_size]));
 string tickcharsize(make_string(m_dparams[i_tick_label_char_size]));
 string titlecharsize(make_string(m_dparams[i_title_char_size]));
 string cmd;
 cmd="with g0"; commands.push_back(cmd);
 cmd="view xmin "+make_string(m_dparams[i_view_xmin]); commands.push_back(cmd);
 cmd="view xmax "+make_string(m_dparams[i_view_xmax]); commands.push_back(cmd);
 cmd="view ymin "+make_string(m_dparams[i_view_ymin]); commands.push_back(cmd);
 cmd="view ymax "+make_string(m_dparams[i_view_ymax]); commands.push_back(cmd);
 cmd="world xmin "+make_string(m_dparams[i_xmin]); commands.push_back(cmd);
 cmd="world xmax "+make_string(m_dparams[i_xmax]); commands.push_back(cmd);
 cmd="world ymin "+make_string(m_dparams[i_ymin]); commands.push_back(cmd);
 cmd="world ymax "+make_string(m_dparams[i_ymax]);  commands.push_back(cmd);
 cmd="title \""+m_title+"\""; commands.push_back(cmd);
 cmd="title font "+make_string(int(m_dparams[i_title_font])); commands.push_back(cmd);
 cmd="title size "+titlecharsize; commands.push_back(cmd);
 cmd="frame linewidth 2.0"; commands.push_back(cmd);
 cmd="autoticks"; commands.push_back(cmd);
 cmd="page background fill off"; commands.push_back(cmd);

 cmd="xaxis bar linewidth 2.0"; commands.push_back(cmd);
 cmd="xaxis label \""+m_xlabel+"\""; commands.push_back(cmd);
 cmd="xaxis label char size "+labelcharsize; commands.push_back(cmd);
 cmd="xaxis tick major size 1.50"; commands.push_back(cmd);
 cmd="xaxis tick major linewidth 2.0"; commands.push_back(cmd);
 cmd="xaxis tick minor linewidth 2.0"; commands.push_back(cmd);
 cmd="xaxis tick minor size 0.75"; commands.push_back(cmd);
 cmd="xaxis ticklabel char size "+tickcharsize; commands.push_back(cmd);
 cmd="xaxis label font "+make_string(int(m_dparams[i_xlabel_font])); commands.push_back(cmd);
 cmd="yaxis bar linewidth 2.0"; commands.push_back(cmd);
 cmd="yaxis label \""+m_ylabel+"\""; commands.push_back(cmd);
 cmd="yaxis label char size "+labelcharsize; commands.push_back(cmd);
 cmd="yaxis tick major size 1.50"; commands.push_back(cmd);
 cmd="yaxis tick major linewidth 2.0"; commands.push_back(cmd);
 cmd="yaxis tick minor linewidth 2.0"; commands.push_back(cmd);
 cmd="yaxis tick minor size 0.75"; commands.push_back(cmd);
 cmd="yaxis ticklabel char size "+tickcharsize; commands.push_back(cmd);
 cmd="yaxis label font "+make_string(int(m_dparams[i_ylabel_font])); commands.push_back(cmd);

 if (m_dparams[i_legend_flag]>0.0){
   cmd="legend on"; commands.push_back(cmd);
   cmd="legend loctype "; cmd+=(m_dparams[i_legend_loctype]>0.0)?"view":"world"; 
        commands.push_back(cmd);
   cmd="legend "+make_string(m_dparams[i_legend_xpos])+", "
                +make_string(m_dparams[i_legend_ypos]); commands.push_back(cmd);
   cmd="legend box color 1"; commands.push_back(cmd);
   cmd="legend box pattern 1"; commands.push_back(cmd);
   cmd="legend box linewidth 1.0"; commands.push_back(cmd);
   cmd="legend box linestyle 1"; commands.push_back(cmd);
   cmd="legend box fill color 0"; commands.push_back(cmd);
   cmd="legend box fill pattern 1"; commands.push_back(cmd);
   cmd="legend font "+make_string(int(m_dparams[i_legend_font])); commands.push_back(cmd);
   cmd="legend char size 1.000000"; commands.push_back(cmd);
   cmd="legend color 1"; commands.push_back(cmd);
   cmd="legend length 4"; commands.push_back(cmd);
   cmd="legend vgap 1"; commands.push_back(cmd);
   cmd="legend hgap 1"; commands.push_back(cmd);
   cmd="legend invert false"; commands.push_back(cmd);}
 else{
   cmd="legend off"; commands.push_back(cmd);}
}


void GracePlot::addText(const string& text, double xloc, double yloc,
                        bool viewport, double fontsize, const string& color,
                        const string& justify, const string& font)
{
 double fsize=(fontsize>0.0)?fontsize:m_dparams[i_text_char_size];
 string textcharsize(make_string(fsize));
 string cmd;
 cmd="with string"; m_text.push_back(cmd);
 cmd="string on"; m_text.push_back(cmd);
 cmd="string loctype "; cmd+=(viewport)?"view":"world"; m_text.push_back(cmd);
 cmd="string g0"; m_text.push_back(cmd);
 cmd="string "+make_string(xloc)+", "+make_string(yloc); m_text.push_back(cmd);
 cmd="string color "+encode_color(color); m_text.push_back(cmd);
 cmd="string rot 0"; m_text.push_back(cmd);
 cmd="string font "+make_string(encode_font(font)); m_text.push_back(cmd);
 cmd="string just "+make_string(encode_justification(justify)); m_text.push_back(cmd);
 cmd="string char size "+textcharsize; m_text.push_back(cmd);
 cmd="string def \""+text+"\""; m_text.push_back(cmd);
}


string GracePlot::encode_color(const string& color)
{
 if ((color=="black")||(color.empty()))
    return "1";
 else if (color=="red")
    return "2";
 else if (color=="blue")
    return "4";
 else if (color=="green")
    return "15";
 else if (color=="yellow")
    return "5";
 else if (color=="magenta")
    return "10";
 else if (color=="cyan")
    return "9";
 else if (color=="orange")
    return "11";
 else if (color=="violet")
    return "8";
 else if (color=="maroon")
    return "13";
 else
    return "1";
}


string GracePlot::encode_symbol(const string& symbol)
{
 if ((symbol=="none")||(symbol.empty()))
    return "0";
 else if (symbol=="circle")
    return "1";
 else if (symbol=="square")
    return "2";
 else if (symbol=="diamond")
    return "3";
 else if (symbol=="triangleup")
    return "4";
 else if (symbol=="triangleleft")
    return "5";
 else if (symbol=="triangledown")
    return "6";
 else if (symbol=="triangleright")
    return "7";
 else if (symbol=="plus")
    return "8";
 else if (symbol=="X")
    return "9";
 else if (symbol=="star")
    return "10";
 else
    return "0";
}

string GracePlot::encode_linestyle(const string& linestyle)
{
 if ((linestyle=="none")||(linestyle.empty()))
    return "0";
 else if (linestyle=="solid")
    return "1";
 else if (linestyle=="dot")
    return "2";
 else if (linestyle=="dash")
    return "3";
 else if (linestyle=="dash2")
    return "4";
 else if (linestyle=="dashdot")
    return "5";
 else if (linestyle=="dashdot2")
    return "6";
 else if (linestyle=="dashdot3")
    return "7";
 else if (linestyle=="dashdot4")
    return "8";
 else
    return "0";
}


string GracePlot::encode_symbolfill(const string& symbolfill)
{
 if ((symbolfill=="none")||(symbolfill=="open")||(symbolfill.empty()))
    return "0";
 else if (symbolfill=="solid")
    return "1";
 else if (symbolfill=="checkered")
    return "3";
 else
    return "0";
}


string GracePlot::encode_datatype(int dataset_type)
{
 if ((dataset_type==0)||(dataset_type==7))
    return "xy";
 else if ((dataset_type==1)||(dataset_type==8))
    return "xydy";
 else if (dataset_type==2)
    return "xydx";
 else if (dataset_type==3)
    return "xydxdy";
 else if (dataset_type==4)
    return "xydxdx";
 else if ((dataset_type==5)||(dataset_type==9))
    return "xydydy";
 else if (dataset_type==6)
    return "xydxdxdydy";
 else
    return "";
}


int GracePlot::encode_font(const std::string& fontname)
{
 if (fontname.empty())
    return 1;
 else if (fontname=="times-roman")
    return 0;
 else if (fontname=="times-italic")
    return 1;
 else if (fontname=="times-bold")
    return 2;
 else if (fontname=="times-bolditalic")
    return 3;
 else if (fontname=="helvetica")
    return 4;
 else if (fontname=="helvetica-oblique")
    return 5;
 else if (fontname=="helvetica-bold")
    return 6;
 else if (fontname=="helvetica-boldoblique")
    return 7;
 else
    return 1;
}


int GracePlot::encode_justification(const std::string& justify)
{
 if (justify.empty())
    return 4;
 else if (justify=="bottom-left")
    return 4;
 else if (justify=="bottom-center")
    return 6;
 else if (justify=="bottom-right")
    return 5;
 else if (justify=="center-left")
    return 12;
 else if (justify=="center-center")
    return 14;
 else if (justify=="center-right")
    return 13;
 else if (justify=="top-left")
    return 8;
 else if (justify=="top-center")
    return 10;
 else if (justify=="top-right")
    return 9;
 else
    return 4;
}



void GracePlot::update_horizontal_span(double x)
{
 if (x<m_dparams[i_val_xmin])
    m_dparams[i_val_xmin]=x;
 if (x>m_dparams[i_val_xmax])
    m_dparams[i_val_xmax]=x;
}

void GracePlot::update_vertical_span(double y)
{
 if (y<m_dparams[i_val_ymin])
    m_dparams[i_val_ymin]=y;
 if (y>m_dparams[i_val_ymax])
    m_dparams[i_val_ymax]=y;
}


void GracePlot::autoScale(double borderfrac)
{
 double sc=borderfrac;
 if (sc<0.0) sc=0.0;
 else if (sc>1.0) sc=1.0;
 double xmin,xmax,ymin,ymax;
 double range=m_dparams[i_val_xmax]-m_dparams[i_val_xmin];
 if (range>0){
    xmin=m_dparams[i_val_xmin]-sc*range;
    xmax=m_dparams[i_val_xmax]+sc*range;
    range=m_dparams[i_val_ymax]-m_dparams[i_val_ymin];
    if (range>0){
       ymin=m_dparams[i_val_ymin]-sc*range;
       ymax=m_dparams[i_val_ymax]+sc*range;}
    else{
       ymin=m_dparams[i_val_ymax]-0.1;
       ymax=m_dparams[i_val_ymax]+0.1;}
    setLimits(xmin,xmax,ymin,ymax);}
}


void GracePlot::autoScale(double leftborderfrac, double rightborderfrac,
                          double topborderfrac, double bottomborderfrac)
{
 double leftsc=leftborderfrac;
 if (leftsc<0.0) leftsc=0.0;
 else if (leftsc>1.0) leftsc=1.0;
 double rightsc=rightborderfrac;
 if (rightsc<0.0) rightsc=0.0;
 else if (rightsc>1.0) rightsc=1.0;
 double topsc=topborderfrac;
 if (topsc<0.0) topsc=0.0;
 else if (topsc>1.0) topsc=1.0;
 double bottomsc=bottomborderfrac;
 if (bottomsc<0.0) bottomsc=0.0;
 else if (bottomsc>1.0) bottomsc=1.0;
 double xmin,xmax,ymin,ymax;
 double range=m_dparams[i_val_xmax]-m_dparams[i_val_xmin];
 if (range>0){
    xmin=m_dparams[i_val_xmin]-leftsc*range;
    xmax=m_dparams[i_val_xmax]+rightsc*range;
    range=m_dparams[i_val_ymax]-m_dparams[i_val_ymin];
    if (range>0){
       ymin=m_dparams[i_val_ymin]-bottomsc*range;
       ymax=m_dparams[i_val_ymax]+topsc*range;}
    else{
       ymin=m_dparams[i_val_ymax]-0.1;
       ymax=m_dparams[i_val_ymax]+0.1;}
    setLimits(xmin,xmax,ymin,ymax);}
}


/*
void GracePlot::drawToScreen(UserInterface *ui, bool keep)
{
 list<string> commands;
 renderCommands(commands);
 if (GraceIsOpen()) GraceClose();
 if (GraceOpen(2048) == -1) {
    throw(std::invalid_argument("Can't run Grace"));}
 doDraw(true);
 if ((ui==0)||(keep)) GraceClosePipe();
 else{
    ui->pressEnterToContinue();
    if (GraceIsOpen()) GraceClose();}
}

void GracePlot::drawToScreenAndSave(const string& filename, UserInterface *ui, 
                                    bool keep)
{
 string fname(tidy_string(filename));
 if (fname.empty()) throw(std::invalid_argument("Invalid file name in GracePlot::drawToScreenAndSave"));
 list<string> commands;
 renderCommands(commands);
 if (GraceIsOpen()) GraceClose();
 if (GraceOpen(2048) == -1) {
    throw(std::invalid_argument("Can't run Grace"));}
 doDraw(true);
 string cmd=string("saveall \"")+fname+"\"";
 GraceCommand(cmd.c_str());
 if ((ui==0)||(keep)) GraceClosePipe();
 else{
    ui->pressEnterToContinue();
    if (GraceIsOpen()) GraceClose();}
}
*/

void GracePlot::saveToFile(const string& filename)
{
 string fname(tidy_string(filename));
 if (fname.empty()) throw(std::invalid_argument("Invalid file name in GracePlot::printToFile"));
 if (GraceIsOpen()) GraceClose();
 char prog[9]="gracebat";
 if (GraceOpenVA(prog,2048,"-hardcopy","-noprint","-nosafe", "-noask", NULL)==-1) {
    throw(std::invalid_argument("Can't run Grace"));}
 doDraw(false);
 string cmd=string("saveall \"")+fname+"\"";
 GraceCommand(cmd.c_str());
 GraceClose();
}

void GracePlot::doDraw(bool redraw)
{
 GraceCommand("page size 440, 440");
 list<string> commands;
 renderCommands(commands);
 if (redraw) GraceCommand("redraw");
 for (list<string>::const_iterator it=commands.begin();it!=commands.end();it++)
    GraceCommand(it->c_str());
 for (list<string>::const_iterator dt=m_data.begin();dt!=m_data.end();dt++)
    GraceCommand(dt->c_str());
 for (list<string>::const_iterator tt=m_text.begin();tt!=m_text.end();tt++)
    GraceCommand(tt->c_str());
 if (redraw) GraceCommand("redraw");
}

/*
void GracePlot::drawToScreen(double borderfrac, UserInterface *ui, 
                             bool keep)
{
 autoScale(borderfrac);
 drawToScreen(ui,keep);
}


void GracePlot::drawToScreenAndSave(const string& filename, double borderfrac, 
                                    UserInterface *ui, bool keep)
{
 autoScale(borderfrac);
 drawToScreenAndSave(filename,ui,keep);
}
*/

void GracePlot::saveToFile(const string& filename, double borderfrac)
{
 autoScale(borderfrac);
 saveToFile(filename);
}


// *************************************************************
