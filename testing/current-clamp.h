/*
 Copyright (C) 2011 Georgia Institute of Technology

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <default_gui_model.h>
#include "include/scatterplot.h"
#include "include/runningstat.h"
#include "include/RTXIprintfilter.h"
#include <QtGui>

class Clamp : public DefaultGUIModel
{

Q_OBJECT

public:

  Clamp(void);
  virtual
  ~Clamp(void);

  void
  execute(void);
  void
  createGUI(DefaultGUIModel::variable_t *, int);

  //Custom flags for clamp mode
  enum mode_t
  {
    RAMP, STEP,
  };

public slots:

  signals: // custom signals

  void setFIRange(double newminamp, double newmaxamp, double minHz, double maxHz);
  void
  newDataPoint(double newx, double newy);
  void
  setStepMode(bool);
  void
  setPlotMode(bool);
  void
  saveImage(QString fileName);
  void
  drawFit(double* x, double* y, int size);
  void
  setEqnMsg(const QString &);

protected:

  virtual void
  update(DefaultGUIModel::update_flags_t);

private:

  // inputs, states, calculated values
  double systime;
  double Vm;
  double dt;
  // parameters
  double minamp;
  double maxamp;
  double stepsize;
  double delay;
  double width;
  double repeat;
  // options
  bool randomize;
  bool plotFI;
  mode_t mode;
  // bookkeeping
  double spktime;
  double prevspktime;
  double newminamp;
  double newmaxamp;
  int stepcount;
  int cyclecount;
  int nstep;
  double delay0;
  bool random;
  long long count;
  int spikecount;
  double spikestate;
  RunningStat ISI;

  double* arrstep;
  double* arrrate;
  double* arrFIamp; // internally store FI data points
  double* arrFIHz;
  double yrangemin;
  double yrangemax;
  // QT components
  QString eqnmsg;
  ScatterPlot *splot;

  // Clamp functions
  void
  initParameters();
  void
  initStepArray();
  void
  bookkeep();
  void
  countspikes();
  bool
  OpenFile(QString);
  QFile dataFile;
  QDataStream stream;

private slots:

  void togglerandom(bool);
  void
  updateClampMode(int);
  void
  toggleFIplot(bool);
  void
  saveFIData();
  void
  fitData();
  void
  print();
  void
  exportSVG();
  void
  clearData();
};
