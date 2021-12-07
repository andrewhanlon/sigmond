#ifndef GRACE_PLOT_H
#define GRACE_PLOT_H

#include <list>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
//#include "user_interface.h"


        // forward declarations

class XYPoint;
class XYDYPoint;
class XYDXPoint;
class XYDXDXPoint;
class XYDXDYPoint;
class XYDYDYPoint;
class XYDXDXDYDYPoint;


 // *****************************************************************************
 // *                                                                           *
 // *   This class creates an xmgrace plot, then draws it on the screen         *
 // *   or sends to a file.  The end user creates an object of this class,      *
 // *   then issues commands to add labels, add points, set limits, etc.        *
 // *   Once finished, the command "drawToScreen" or "saveToFile"               *
 // *   member is used.  "drawToScreen" is only useful in an interactive        *
 // *   mode, whereas "saveToFile" is useful in batch and interactive mode.     *
 // *   This is mainly a wrapper to the grace_np library,                       *
 // *   but much more user friendly, since colors, symbols, etc. are            *
 // *   specified in strings.  This class does not use the full functionality   *
 // *   of the grace library; it can create only those kinds of plots that      *
 // *   are usually needed in Monte Carlo data analysis.  This class is         *
 // *   meant to be used by objects of the class "PlotHandler".                 *
 // *                                                                           *
 // *   To add XY data, one first uses an "addXYDataSet" which specifies        *
 // *   the symbol type, color, line style, legend information, etc.            *
 // *   One then uses "addXYDataPoint" or "addXYDataPoints" to add              *
 // *   data points to the data.  You can only add points to the current        *
 // *   data set.  For other kinds of data points, the                          *
 // *                                                                           *
 // *   Color choices:                                                          *
 // *    "black" "red" "blue" "green" "yellow" "magenta" "cyan"                 *
 // *    "orange" "violet" "maroon"                                             *
 // *                                                                           *
 // *   Symbol choices:                                                         *
 // *    "none" "circle" "square" "diamond" "triangleup" "triangleleft"         *
 // *    "triangledown" "triangleright" "plus" "X" "star"                       *
 // *                                                                           *
 // *   Line style choices:                                                     *
 // *    "none" "solid" "dot" "dash" "dash2" "dashdot" "dashdot2"               *
 // *    "dashdot3" "dashdot4"                                                  *
 // *                                                                           *
 // *   Symbol fill choices:                                                    *
 // *    "solid", "open", "checkered"                                           *
 // *                                                                           *
 // *   fill choices:                                                           *
 // *    "solid", "open", "checkered", "dotted"                                 *
 // *                                                                           *
 // *   Font choices:                                                           *
 // *     "times-roman", "times-italic", "times-bold", "times-bolditalic",      *
 // *     "helvetica", "helvetica-oblique", "helvetica-bold",                   *
 // *     "helvetica-boldoblique"                                               *
 // *                                                                           *
 // *   Justification choices:                                                  *
 // *     "bottom-left", "bottom-center", "bottom-right",                       *
 // *     "center-left", "center-center", "center-right",                       *
 // *     "top-left", "top-center", "top-right"                                 *
 // *                                                                           *
 // *   When using "drawToScreen", a boolean "keep" is used.  If true,          *
 // *   the plot is drawn on the screen, and the program continues,             *
 // *   leaving the plot on the screen.  No further interaction with            *
 // *   the plot is possible; the user must close the xmgrace window.           *
 // *   If "keep" is false, the plot is drawn on the screen for the             *
 // *   user to see, but the program waits for the user to press <enter>        *
 // *   to continue.  Upon pressing <enter>, the plot is erased from            *
 // *   the screen.                                                             *
 // *                                                                           *
 // *****************************************************************************



class GracePlot
{

   std::vector<double> m_dparams;
   std::string m_title, m_xlabel, m_ylabel;
   int m_curr_dataset_index, m_curr_dataset_type;
   std::list<std::string> m_data;
   std::list<std::string> m_text;
   int m_dcount;
   std::string m_dset;

 public:

   GracePlot();
   GracePlot(double xmin, double xmax, double ymin, double ymax);
   GracePlot(double xmin, double xmax, double ymin, double ymax,
             const std::string& xlabel, const std::string& ylabel);
   GracePlot(double xmin, double xmax, double ymin, double ymax,
             const std::string& xlabel, const std::string& ylabel,
             const std::string& title);
   GracePlot(const std::string& xlabel, const std::string& ylabel);
   GracePlot(const std::string& xlabel, const std::string& ylabel,
             const std::string& title);
   GracePlot(const GracePlot& gplot);
   GracePlot(GracePlot& gplot);
   GracePlot& operator=(const GracePlot& gplot);
   GracePlot& operator=(GracePlot& gplot);
   ~GracePlot();

   void setDefaults();
   void setView(double view_xmin, double view_xmax, double view_ymin, double view_ymax);
   void setLabels(const std::string& xlabel, const std::string& ylabel);
   void setLabels(const std::string& xlabel, const std::string& ylabel,
                  const std::string& title);
   void setFontsizes(double ticklabelcharsize, double axislabelcharsize, 
                     double titlecharsize, double textcharsize);
   void setFonts(const std::string& xlabelfont, const std::string& ylabelfont,
                 const std::string& titlefont, const std::string& legendfont);
   void setLimits(double xmin, double xmax, double ymin, double ymax);
   void setVerticalLimits(double ymin, double ymax);
   void setHorizontalLimits(double xmin, double xmax);
   void setLegend(double xviewpos, double yviewpos, bool viewport=true);

   void addText(const std::string& text, double xloc, double yloc, bool viewport=true,
                double fontsize=0, const std::string& color="black", 
                const std::string& justify="bottom-left",
                const std::string& font="times-roman");

   void autoScale(double borderfrac);
   void autoScale(double leftborderfrac, double rightborderfrac,
                  double topborderfrac, double bottomborderfrac);


   void addXYDataSet(const std::string& symbol, const std::string& symbolfill,
                     const std::string& linestyle, const std::string& color, 
                     const std::string& legendtext="");
   void addXYDataPoint(double xval, double yval);
   void addXYDataPoint(const XYPoint& pt);
   void addXYDataPoints(const std::vector<XYPoint>& pts);

   void addXYDYDataSet(const std::string& symbol, const std::string& symbolfill,
                       const std::string& linestyle, const std::string& color, 
                       const std::string& legendtext="");
   void addXYDYDataPoint(double x, double y, double dy);
   void addXYDYDataPoint(const XYDYPoint& pt);
   void addXYDYDataPoints(const std::vector<XYDYPoint>& pts);

   void addXYDXDataSet(const std::string& symbol, const std::string& symbolfill,
                       const std::string& linestyle, const std::string& color, 
                       const std::string& legendtext="");
   void addXYDXDataPoint(double x, double y, double dx);
   void addXYDXDataPoint(const XYDXPoint& pt);
   void addXYDXDataPoints(const std::vector<XYDXPoint>& pts);

   void addXYDXDYDataSet(const std::string& symbol, const std::string& symbolfill,
                         const std::string& linestyle, const std::string& color, 
                         const std::string& legendtext="");
   void addXYDXDYDataPoint(double x, double y, double dx, double dy);
   void addXYDXDYDataPoint(const XYDXDYPoint& pt);
   void addXYDXDYDataPoints(const std::vector<XYDXDYPoint>& pts);

   void addXYDXDXDataSet(const std::string& symbol, const std::string& symbolfill,
                         const std::string& linestyle, const std::string& color, 
                         const std::string& legendtext="");
   void addXYDXDXDataPoint(double x, double y, double dxup, double dxdn);
   void addXYDXDXDataPoint(const XYDXDXPoint& pt);
   void addXYDXDXDataPoints(const std::vector<XYDXDXPoint>& pts);

   void addXYDYDYDataSet(const std::string& symbol, const std::string& symbolfill,
                         const std::string& linestyle, const std::string& color, 
                         const std::string& legendtext="");
   void addXYDYDYDataPoint(double x, double y, double dyup, double dydn);
   void addXYDYDYDataPoint(const XYDYDYPoint& pt);
   void addXYDYDYDataPoints(const std::vector<XYDYDYPoint>& pts);

   void addXYDXDXDYDYDataSet(const std::string& symbol, const std::string& symbolfill,
                             const std::string& linestyle, const std::string& color, 
                             const std::string& legendtext="");
   void addXYDXDXDYDYDataPoint(double x, double y, double dxup, double dxdn,
                               double dyup, double dydn);
   void addXYDXDXDYDYDataPoint(const XYDXDXDYDYPoint& pt);
   void addXYDXDXDYDYDataPoints(const std::vector<XYDXDXDYDYPoint>& pts);


   void addBarDataSet(const std::string& color, const std::string& bordercolor,
                      double barwidth, const std::string& legendtext="");
   void addBarDataPoint(double x, double y);
   void addBarDataPoint(const XYPoint& pt);
   void addBarDataPoints(const std::vector<XYPoint>& pts);


   void addBarDYDataSet(const std::string& color, const std::string& bordercolor,
                        double barwidth, const std::string& legendtext="");
   void addBarDYDataPoint(double x, double y, double yerr);
   void addBarDYDataPoint(const XYDYPoint& pt);
   void addBarDYDataPoints(const std::vector<XYDYPoint>& pts);


   void addBarDYDYDataSet(const std::string& color, const std::string& bordercolor,
                          double barwidth, const std::string& legendtext="");
   void addBarDYDYDataPoint(double x, double y, double dyup, double dydn);
   void addBarDYDYDataPoint(const XYDYDYPoint& pt);
   void addBarDYDYDataPoints(const std::vector<XYDYDYPoint>& pts);


   void addFillDataSet(const std::string& color, const std::string& fillpattern);
   void addFillDataPoint(double x, double y);
   void addFillDataPoint(const XYPoint& pt);
   void addFillDataPoints(const std::vector<XYPoint>& pts);


   //void drawToScreen(UserInterface *ui=0, bool keep=true);
   //void drawToScreenAndSave(const std::string& filename, UserInterface *ui=0, 
   //                         bool keep=true);
   void saveToFile(const std::string& filename);

   //void drawToScreen(double borderfrac, UserInterface *ui=0, bool keep=true);
   //void drawToScreenAndSave(const std::string& filename, double borderfrac, 
   //                         UserInterface *ui=0, bool keep=true);
   void saveToFile(const std::string& filename, double borderfrac);


 private:

   void renderCommands(std::list<std::string>& commands);
   void doDraw(bool redraw);

   static const int i_axis_label_char_size=0;
   static const int i_tick_label_char_size=1;
   static const int i_title_char_size=2;
   static const int i_text_char_size=3;
   static const int i_view_xmin=4;
   static const int i_view_xmax=5; 
   static const int i_view_ymin=6;
   static const int i_view_ymax=7;
   static const int i_xmin=8;
   static const int i_xmax=9; 
   static const int i_ymin=10;
   static const int i_ymax=11;
   static const int i_legend_flag=12;
   static const int i_legend_xpos=13;
   static const int i_legend_ypos=14;
   static const int i_legend_loctype=15;
   static const int i_val_xmin=16;
   static const int i_val_xmax=17; 
   static const int i_val_ymin=18;
   static const int i_val_ymax=19;
   static const int i_temp=20;
   static const int i_xlabel_font=21;
   static const int i_ylabel_font=22;
   static const int i_title_font=23;
   static const int i_legend_font=24;
   static const int dparamsize=25;

   std::string make_string(int ival);
   std::string make_string(double dval);
   std::string tidy_string(const std::string& str);

      // encode input strings into string codes (cycle if empty string)
   std::string encode_color(const std::string& color);
   std::string encode_symbol(const std::string& symbol);
   std::string encode_symbolfill(const std::string& symbolfill);
   std::string encode_fill(const std::string& fill);
   std::string encode_linestyle(const std::string& linestyle);
   std::string encode_datatype(int dataset_type);
   int encode_font(const std::string& fontname);
   int encode_justification(const std::string& justify);

   void addDataSet(int dataset_type, bool errorbars, const std::string& symbol, 
                   const std::string& symbolfill, const std::string& linestyle, 
                   const std::string& color, const std::string& legendtext);

   void addPoint(const std::string& prefix, const XYPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDYPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDXPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDXDXPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDXDYPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDYDYPoint& pt, int& count);
   void addPoint(const std::string& prefix, const XYDXDXDYDYPoint& pt, int& count);

   void addBarSet(int dataset_type, bool errorbars, const std::string& color, 
                  const std::string& bordercolor, double barwidth, 
                  const std::string& legendtext);

   void addBar(const std::string& prefix, const XYPoint& pt, int& count);
   void addBar(const std::string& prefix, const XYDYPoint& pt, int& count);
   void addBar(const std::string& prefix, const XYDYDYPoint& pt, int& count);

   void addFillSet(int dataset_type, const std::string& color, const std::string& fillpattern);

   void update_horizontal_span(double x);
   void update_vertical_span(double y);

};


  // *******************************************************************


class XYPoint
{
 public:
   double xval;
   double yval;

   XYPoint(double x, double y) : xval(x), yval(y) {}
   XYPoint& operator=(const XYPoint& indata) 
    {xval=indata.xval; yval=indata.yval; return *this;}
   XYPoint() : xval(0.0), yval(0.0) {}
};

class XYDYPoint
{
 public:
   double xval;
   double yval;
   double yerr;

   XYDYPoint(double x, double y, double dy) : xval(x), yval(y), yerr(dy) {}
   XYDYPoint& operator=(const XYDYPoint& indata) 
    {xval=indata.xval; yval=indata.yval; yerr=indata.yerr; return *this;}
   XYDYPoint() : xval(0.0), yval(0.0), yerr(0.0) {}
};


class XYDXPoint
{
 public:
   double xval;
   double yval;
   double xerr;

   XYDXPoint(double x, double y, double dx) : xval(x), yval(y), xerr(dx) {}
   XYDXPoint& operator=(const XYDXPoint& indata) 
    {xval=indata.xval; yval=indata.yval; xerr=indata.xerr; return *this;}
   XYDXPoint() : xval(0.0), yval(0.0), xerr(0.0) {}
};


class XYDXDXPoint
{
 public:
   double xval;
   double yval;
   double xuperr;
   double xdnerr;

   XYDXDXPoint(double x, double y, double dxup, double dxdn) 
              : xval(x), yval(y), xuperr(dxup), xdnerr(dxdn) {}
   XYDXDXPoint& operator=(const XYDXDXPoint& indata) 
    {xval=indata.xval; yval=indata.yval; 
     xuperr=indata.xuperr; xdnerr=indata.xdnerr; return *this;}
   XYDXDXPoint() : xval(0.0), yval(0.0), xuperr(0.0), xdnerr(0.0){}
};


class XYDXDYPoint
{
 public:
   double xval;
   double yval;
   double xerr;
   double yerr;

   XYDXDYPoint(double x, double y, double dx, double dy) 
              : xval(x), yval(y), xerr(dx), yerr(dy) {}
   XYDXDYPoint& operator=(const XYDXDYPoint& indata) 
    {xval=indata.xval; yval=indata.yval; 
     xerr=indata.xerr; yerr=indata.yerr; return *this;}
   XYDXDYPoint() : xval(0.0), yval(0.0), xerr(0.0), yerr(0.0){}
};


class XYDYDYPoint
{
 public:
   double xval;
   double yval;
   double yuperr;
   double ydnerr;

   XYDYDYPoint(double x, double y, double dyup, double dydn) 
                : xval(x), yval(y), yuperr(dyup), ydnerr(dydn) {}
   XYDYDYPoint& operator=(const XYDYDYPoint& indata) 
    {xval=indata.xval; yval=indata.yval; 
     yuperr=indata.yuperr; ydnerr=indata.ydnerr; return *this;}
   XYDYDYPoint() : xval(0.0), yval(0.0), yuperr(0.0), ydnerr(0.0) {}
};


class XYDXDXDYDYPoint
{
 public:
   double xval;
   double yval;
   double xuperr;
   double xdnerr;
   double yuperr;
   double ydnerr;

   XYDXDXDYDYPoint(double x, double y, double dxup, double dxdn,
                   double dyup, double dydn) 
                : xval(x), yval(y), xuperr(dxup), xdnerr(dxdn),
                  yuperr(dyup), ydnerr(dydn) {}
   XYDXDXDYDYPoint& operator=(const XYDXDXDYDYPoint& indata) 
    {xval=indata.xval; yval=indata.yval; 
     xuperr=indata.xuperr; xdnerr=indata.xdnerr;
     yuperr=indata.yuperr; ydnerr=indata.ydnerr; return *this;}
   XYDXDXDYDYPoint() : xval(0.0), yval(0.0), xuperr(0.0), xdnerr(0.0),
                       yuperr(0.0), ydnerr(0.0) {}
};


// **************************************************************
#endif
