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

/*
 * This module allows you to deliver current clamp stimuli. With the SpikeDetect module it can also
 * plot an F-I curve.
 */

#include <current-clamp.h>
#include <algorithm>
#include <main_window.h>
#include <QtGui>

#if QT_VERSION >= 0x040300
#ifdef QT_SVG_LIB
#endif
#endif
#if QT_VERSION >= 0x040000
#else
#endif

#include <time.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new Clamp();
}

static DefaultGUIModel::variable_t vars[] =
  {
    { "Spike State", "Spike State", DefaultGUIModel::INPUT, },
    { "Command", "Command", DefaultGUIModel::OUTPUT, },
    { "Min Amplitude (pA)", "Minimum amplitude (pA) of command",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Max Amplitude (pA)", "Maximum amplitude (pA) of command",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Step Size (pA)",
        "Step size between minimum and maximum amplitude (pA) of command",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Delay (s)", "Delay (s) after last command", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "Width (s)", "Width (s) of command", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "Repeat", "Number of times to repeat cycle", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, };

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

Clamp::Clamp(void) :
  DefaultGUIModel("Current Clamp", ::vars, ::num_vars)
{

  setWhatsThis(
      "<p><b>Current Clamp</b></p><p>This plugin requires the SpikeDetect plugin, which indicates when a spike has occurred using a threshold method. This plugin allows you to deliver current steps or ramps. Choose the <b>Clamp Mode</b>, then set the parameters for the input amplitudes and step size. This plugin uses an absolute delay between the end of a command input and the beginning of the next command. You can repeat the series of command inputs as many times as you like. Use the <b>Randomize</b> checkbox to randomize amplitudes within each cycle. To run a current clamp protocol, toggle the <b>Pause</b> button. You can edit parameter values directly in the textboxes, but you must click <b>Modify</b> to commit the changes. In <b>Step Mode</b>, you can choose to plot the FI Curve as it is generated and save a screenshot of the plot or save the data to a plain text file. The frequency is computed by averaging all the ISIs detected during a command step and taking the reciprocal. A linear fit is performed on all data points currently displayed in the scatterplot, not simply the data points acquired in the last run. Saving data points to a file will only save the data points acquired in the last run. You can completely overwrite or append the data to an existing file.</p>");
  initParameters();
  initStepArray();
  createGUI(vars, num_vars); // this is required to create the GUI
  update(INIT);
  refresh(); // this is required to update the GUI with parameter and state values
  printf("Loaded Current Clamp:\n");

}

Clamp::~Clamp(void)
{
}

void
Clamp::execute(void)
{
  if (plotFI == true)
    {
      spikestate = input(0);
    }
  systime = count * dt; // current time, s

  if (cyclecount < repeat)
    { // as long as there are cycles remaining
      if (stepcount == 0 && randomize == true && random == false)
        { // shuffle amplitudes once at beginning of each cycle
          std::random_shuffle(arrstep, arrstep + nstep);
          random = true;
        }
      if (stepcount < nstep)
        { // as long as there are steps remaining in cycle
          if (systime >= delay && systime < (delay + width))
            { // during command window, set output
              if (mode == STEP)
                {
                  output(0) = arrstep[stepcount];
                }
              else if (mode == RAMP)
                {
                  if (systime >= delay && systime < (delay + width / 2))
                    {
                      output(0) += arrrate[stepcount];
                    }
                  else
                    {
                      output(0) -= arrrate[stepcount];
                    }
                }
              if ((spikestate == 1) && (plotFI == true))
                countspikes();
            }
          else if (systime >= delay + width)
            { // after command window, shift delay, clean up
              delay = delay + delay0 + width;
              if (plotFI == true)
                {
                  arrFIamp[static_cast<int> (cyclecount * nstep + stepcount)]
                      = arrstep[stepcount] * 1e12; // in pA
                  if (spikecount >= 2)
                    { // at least 2 spikes (1 ISI) occurred
                      arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)]
                          = 1 / ISI.mean();
                      if (ISI.mean() != 0 && 1 / ISI.mean() > yrangemax)
                        { // resize window
                          yrangemax = 1 / ISI.mean();
                          emit setFIRange(minamp * 1e12, maxamp * 1e12,
                              yrangemin, yrangemax);
                        }
                      if (ISI.mean() != 0)
                        {
                          printf("%f         %f\n", arrstep[stepcount] * 1e12,
                              1 / ISI.mean());
                          emit newDataPoint(arrFIamp[cyclecount * nstep
                              + stepcount], arrFIHz[cyclecount * nstep
                              + stepcount]);
                        }
                    }
                  else if (spikecount == 1)
                    { // only 1 spike occurred
                      arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)]
                          = GSL_NAN;
                    }
                  else
                    { // no spikes occurred
                      arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)]
                          = 0;
                    }
                  spikecount = 0;
                  ISI.clear();
                }
              stepcount++;
              output(0) = 0;
            }
          else
            { // outside command window, set output to zero
              output(0) = 0;
            }
        }
      if (stepcount == nstep)
        { // increment cycle count, reset flag for random shuffle
          cyclecount++;
          splot->nextSymbol();
          stepcount = 0;
          if (randomize == true)
            random = false;
        }
    }
  else
    { // all cycles are done
    } // end if cyclecount
  count++; // increment count to measure time
  return;
}

void
Clamp::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
    {
  case INIT:
    setParameter("Min Amplitude (pA)", QString::number(minamp * 1e12)); // initialized in A, display in pA
    setParameter("Max Amplitude (pA)", QString::number(maxamp * 1e12)); // initialized in A
    setParameter("Step Size (pA)", QString::number(stepsize * 1e12)); // initialized in A
    setParameter("Delay (s)", QString::number(delay)); //
    setParameter("Width (s)", QString::number(width)); //
    setParameter("Repeat", QString::number(repeat)); // initially 1
    setState("Time (s)", systime);
    //emit setFIRange(minamp * 1e12, maxamp * 1e12, yrangemin, yrangemax);
    break;
  case MODIFY:
    delay = getParameter("Delay (s)").toDouble();
    delay0 = delay;
    width = getParameter("Width (s)").toDouble();
    minamp = getParameter("Min Amplitude (pA)").toDouble() * 1e-12; // set by user in pA
    maxamp = getParameter("Max Amplitude (pA)").toDouble() * 1e-12; // set by user in pA
    if (minamp == maxamp)
      {
        QMessageBox::information(this, "Current Clamp",
            "The minimum and maximum amplitude cannot\n"
              "be the same. The defaults will be set.\n");
        minamp = -100;
        maxamp = 400;
        setParameter("Min Amplitude (pA)", QString::number(minamp * 1e12)); // initialized in A, display in pA
        setParameter("Max Amplitude (pA)", QString::number(maxamp * 1e12)); // initialized in A, display in pA
      }
    stepsize = getParameter("Step Size (pA)").toDouble() * 1e-12; // set by user in pA
    repeat = getParameter("Repeat").toDouble();
    nstep = int((maxamp + 1e-12 - minamp) / stepsize) + 1; // recalculate number of steps and array
    initStepArray();
    emit setFIRange(minamp * 1e12, maxamp * 1e12, yrangemin, yrangemax);
    eqnmsg = "Y = c0 + c1 * X";
    //emit setEqnMsg(eqnmsg);
    bookkeep();
    break;
  case UNPAUSE:
    bookkeep();
    delay = getParameter("Delay (s)").toDouble(); // set by user in s, convert to ms
    delay0 = delay;
    initStepArray();
    printf(
        "Starting protocol: minamp: %f maxamp: %f stepsize %f delay %f width %f\n",
        minamp * 1e12, maxamp * 1e12, stepsize * 1e12, delay, width);
    printf("Amplitude (pA)  Rate (Hz)\n");
    break;
  case PAUSE:
    output(0) = 0; // stop command in case pause occurs in the middle of command
    printf("Protocol paused.\n");
    break;
  case PERIOD:
    dt = RT::System::getInstance()->getPeriod() * 1e-9;
    initStepArray();
    break;
  default:
    break;

    }

}

void
Clamp::createGUI(DefaultGUIModel::variable_t *var, int size)
{

  setMinimumSize(200, 300); // Qt API for setting window size

  //overall GUI layout with a "horizontal box" copied from DefaultGUIModel
  QBoxLayout *layout = new QHBoxLayout(this);

  // left and right panels
  // Right side GUI with buttons and FI plot
  QBoxLayout *rightlayout = new QVBoxLayout();
  QHButtonGroup *plotBox = new QHButtonGroup("FI Plot:", this);
  QPushButton *clearButton = new QPushButton("&Clear", plotBox);
  QPushButton *linearfitButton = new QPushButton("Linear &Fit", plotBox);
  QPushButton *savePlotButton = new QPushButton("Save Screenshot", plotBox);
  QPushButton *printButton = new QPushButton("Print", plotBox);
  QPushButton *saveDataButton = new QPushButton("Save FI Data", plotBox);
  QLineEdit *eqnLine = new QLineEdit(this, "Linear Equation");
  eqnLine->setText("Y = c0 + c1 * X");
  eqnLine->setFrame(false);
  splot = new ScatterPlot(this);
  // Connect buttons to functions
  QObject::connect(clearButton, SIGNAL(clicked()), splot, SLOT(clear()));
  QObject::connect(clearButton, SIGNAL(clicked()), this, SLOT(clearData()));
  QObject::connect(savePlotButton, SIGNAL(clicked()), this, SLOT(exportSVG()));
  QObject::connect(printButton, SIGNAL(clicked()), this, SLOT(print()));
  QObject::connect(saveDataButton, SIGNAL(clicked()), this, SLOT(saveFIData()));
  QObject::connect(linearfitButton, SIGNAL(clicked()), this, SLOT(fitData()));
  QToolTip::add(clearButton, "Clear");
  QToolTip::add(savePlotButton, "Save screenshot");
  QToolTip::add(saveDataButton, "Save data");
  QToolTip::add(linearfitButton, "Perform linear least-squares regression");
  QToolTip::add(printButton, "Print plot");

  rightlayout->addWidget(plotBox);
  rightlayout->addWidget(eqnLine);
  rightlayout->addWidget(splot);
  plotBox->hide();
  eqnLine->hide();
  splot->setFixedSize(540, 300);
  splot->hide();

  QBoxLayout *leftlayout = new QVBoxLayout();

  QHButtonGroup *modeBox = new QHButtonGroup("Clamp Mode", this);
  modeBox->setRadioButtonExclusive(true);
  QRadioButton *stepButton = new QRadioButton("Step", modeBox);
  stepButton->setChecked(true);
  QRadioButton *rampButton = new QRadioButton("Ramp", modeBox);
  QObject::connect(modeBox,SIGNAL(clicked(int)),this,SLOT(updateClampMode(int)));
  QToolTip::add(stepButton, "Set mode to current steps");
  QToolTip::add(rampButton, "Set mode to triangular current ramps");

  QVBox *optionBox = new QVBox(this);
  QHBox *optionRow1 = new QHBox(optionBox);
  QCheckBox *randomCheckBox = new QCheckBox("Randomize", optionRow1);
  QHBox *optionRow2 = new QHBox(optionBox);
  QCheckBox *plotFICheckBox = new QCheckBox("Plot FI Curve", optionRow2);
  QObject::connect(randomCheckBox,SIGNAL(toggled(bool)),this,SLOT(togglerandom(bool)));
  QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),eqnLine,SLOT(setShown(bool)));
  QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),splot,SLOT(setShown(bool)));
  QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),plotBox,SLOT(setShown(bool)));
  QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleFIplot(bool)));
  QToolTip::add(randomCheckBox, "Randomize input amplitudes within a cycle");
  QToolTip::add(plotFICheckBox, "Show/Hide FI plot area");

  QHBox *utilityBox = new QHBox(this);
  pauseButton = new QPushButton("Pause", utilityBox);
  pauseButton->setToggleButton(true);
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),this,SLOT(pause(bool)));
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),savePlotButton,SLOT(setEnabled(bool)));
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),printButton,SLOT(setEnabled(bool)));
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),saveDataButton,SLOT(setEnabled(bool)));
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),linearfitButton,SLOT(setEnabled(bool)));
  QPushButton *modifyButton = new QPushButton("Modify", utilityBox);
  QObject::connect(modifyButton,SIGNAL(clicked(void)),this,SLOT(modify(void)));
  QPushButton *unloadButton = new QPushButton("Unload", utilityBox);
  QObject::connect(unloadButton,SIGNAL(clicked(void)),this,SLOT(exit(void)));
  QObject::connect(pauseButton,SIGNAL(toggled(bool)),modifyButton,SLOT(setEnabled(bool)));
  QToolTip::add(pauseButton, "Start/Stop current clamp protocol");
  QToolTip::add(modifyButton, "Commit changes to parameter values");
  QToolTip::add(unloadButton, "Close module");

  QObject::connect(this,SIGNAL(newDataPoint(double,double)),splot,SLOT(appendPoint(double,double)));
  QObject::connect(this,SIGNAL(setFIRange(double, double, double, double)),splot,SLOT(setAxes(double, double, double, double)));
  QObject::connect(this,SIGNAL(setPlotMode(bool)),plotFICheckBox,SLOT(setChecked(bool)));
  QObject::connect(this,SIGNAL(setStepMode(bool)),plotFICheckBox,SLOT(setEnabled(bool)));
  QObject::connect(this,SIGNAL(setStepMode(bool)),plotBox,SLOT(setEnabled(bool)));
  QObject::connect(this,SIGNAL(drawFit(double*, double*, int)),splot,SLOT(appendLine(double*, double*, int)));
  QObject::connect(this,SIGNAL(setEqnMsg(const QString &)), eqnLine,SLOT(setText(const QString &)));

  // add custom button group at the top of the layout
  leftlayout->addWidget(modeBox);

  // create default_gui_model GUI DO NOT EDIT
  QScrollView *sv = new QScrollView(this);
  sv->setResizePolicy(QScrollView::AutoOneFit);
  leftlayout->addWidget(sv);

  QWidget *viewport = new QWidget(sv->viewport());
  sv->addChild(viewport);
  QGridLayout *scrollLayout = new QGridLayout(viewport, 1, 2);

  size_t nstate = 0, nparam = 0, nevent = 0, ncomment = 0;
  for (size_t i = 0; i < num_vars; i++)
    {
      if (vars[i].flags & (PARAMETER | STATE | EVENT | COMMENT))
        {
          param_t param;

          param.label = new QLabel(vars[i].name, viewport);
          scrollLayout->addWidget(param.label, parameter.size(), 0);
          param.edit = new DefaultGUILineEdit(viewport);
          scrollLayout->addWidget(param.edit, parameter.size(), 1);

          QToolTip::add(param.label, vars[i].description);
          QToolTip::add(param.edit, vars[i].description);

          if (vars[i].flags & PARAMETER)
            {
              if (vars[i].flags & DOUBLE)
                {
                  param.edit->setValidator(new QDoubleValidator(param.edit));
                  param.type = PARAMETER | DOUBLE;
                }
              else if (vars[i].flags & UINTEGER)
                {
                  QIntValidator *validator = new QIntValidator(param.edit);
                  param.edit->setValidator(validator);
                  validator->setBottom(0);
                  param.type = PARAMETER | UINTEGER;
                }
              else if (vars[i].flags & INTEGER)
                {
                  param.edit->setValidator(new QIntValidator(param.edit));
                  param.type = PARAMETER | INTEGER;
                }
              else
                param.type = PARAMETER;
              param.index = nparam++;
              param.str_value = new QString;
            }
          else if (vars[i].flags & STATE)
            {
              param.edit->setReadOnly(true);
              param.edit->setPaletteForegroundColor(Qt::darkGray);
              param.type = STATE;
              param.index = nstate++;
            }
          else if (vars[i].flags & EVENT)
            {
              param.edit->setReadOnly(true);
              param.type = EVENT;
              param.index = nevent++;
            }
          else if (vars[i].flags & COMMENT)
            {
              param.type = COMMENT;
              param.index = ncomment++;
            }

          parameter[vars[i].name] = param;
        }
    }

  // end default_gui_model GUI DO NOT EDIT

  leftlayout->addWidget(optionBox);
  leftlayout->addWidget(utilityBox);
  // Add left and right side layouts to the overall layout
  layout->addLayout(leftlayout);
  layout->addLayout(rightlayout);
  layout->setResizeMode(QLayout::Fixed);

  show(); // this line is required to render the GUI
}

void
Clamp::initParameters()
{
  minamp = -100e-12; // A
  maxamp = 400e-12; // A
  stepsize = 50e-12; // A
  delay = 2; // s
  width = 1; // s
  bookkeep();
  randomize = false;
  plotFI = false;
  repeat = 1;
  dt = RT::System::getInstance()->getPeriod() * 1e-9; // s
  mode = STEP;
  srand(time(NULL));
  delay0 = delay; // delay0 is the incremental delay, delay is the absolute total delay from time zero
  nstep = int((maxamp + 1e-12 - minamp) / stepsize) + 1; // calculate the number of amplitude steps
  yrangemin = 0;
  yrangemax = 50;

}

void
Clamp::initStepArray()
{
  arrstep = new double[nstep];
  arrrate = new double[nstep];
  arrFIamp = new double[static_cast<int> (nstep * repeat)];
  arrFIHz = new double[static_cast<int> (nstep * repeat)];
  for (int i = 0; i < nstep * repeat; i++)
    {
      arrFIamp[i] = 0;
      arrFIHz[i] = 0;
    }
  for (int i = 0; i < nstep; i++)
    {
      arrstep[i] = minamp + i * stepsize;
      arrrate[i] = arrstep[i] / width * 2 * dt;
    }
  random = false;
}

void
Clamp::updateClampMode(int index)
{
  if (index == 0)
    { // STEP
      mode = STEP;
      update(MODIFY);
      //emit setStepMode(true);
      //emit setPlotMode(false);
      printf("Entering STEP mode\n");
    }
  else if (index == 1)
    { // RAMP
      mode = RAMP;
      minamp = 0;
      setParameter("Min Amplitude (pA)", QString::number(minamp * 1e12)); // initialized in A, display in pA
      update(MODIFY);
      //emit setStepMode(false);
      plotFI = false;
      printf("Entering RAMP mode\n");
    }
}

void
Clamp::bookkeep()
{
  stepcount = 0;
  cyclecount = 0;
  count = 0;
  systime = 0;
  spikestate = 0;
  spikecount = 0;
  spktime = 0;
  prevspktime = 0;
  ISI.clear();
}

void
Clamp::clearData()
{
  yrangemax = 50;
  //emit setFIRange(minamp * 1e12, maxamp * 1e12, yrangemin, yrangemax);
  eqnmsg = "Y = c0 + c1 * X";
  //emit setEqnMsg(eqnmsg);
}

void
Clamp::saveFIData()
{
  QFileDialog* fd = new QFileDialog(this, "Save File As", TRUE);
  fd->setMode(QFileDialog::AnyFile);
  fd->setViewMode(QFileDialog::Detail);
  QString fileName;
  if (fd->exec() == QDialog::Accepted)
    {
      fileName = fd->selectedFile();

      if (OpenFile(fileName))
        {
          stream.setPrintableData(true);
          for (int i = 0; i < cyclecount * nstep + stepcount; i++)
            {
              stream << (double) arrFIamp[i] << (double) arrFIHz[i];
            }
          dataFile.close();
          printf("File closed.\n");
        }
      else
        {
          QMessageBox::information(this, "Current Clamp: Save FI Curve",
              "There was an error writing to this file. You can view\n"
                "the values that should be plotted in the terminal.\n");
        }
    }
}

void
Clamp::fitData()
{
  if (splot->dataExists())
    {
      int n = splot->dataSize();
      if (n >= 2)
        {
          const double* tempx = splot->xData();
          const double* tempy = splot->yData();
          double x[n];
          double y[n];
          int j = 0;
          double maxVal = 0;
          double minVal = 0;
          for (int i = 0; i < n; i++)
            { // exclude any inf or NaN or zero rates from the FI linear fit
              if (gsl_finite(tempy[i]) == 1 || tempy[i] == 0)
                {
                  x[j] = tempx[i];
                  y[j] = tempy[i];
                  if (x[j] > maxVal)
                    maxVal = x[j];
                  if (x[j] < minVal)
                    minVal = x[j];
                  j++;
                }
            }
          double c0, c1, cov00, cov01, cov11, sumsq;
          gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11,
              &sumsq);
          // gsl_stats_tss does not exist in GSL 1.10 included in Ubuntu 9.10
          // use this function if you manually upgrade GSL, otherwise use the next section of code
          //double SST = gsl_stats_tss(y,1,j);
          double SST = 0; //calculating total sum of squares around the mean
          double ymean = gsl_stats_mean(y, 1, j);
          for (int i = 0; i < j; i++)
            {
              double delta = y[i] - ymean;
              SST += delta * delta;
            }
          if (gsl_finite(c0) == 1 && gsl_finite(c1) == 1)
            {
              printf("Best fit: Y = %g + %g X\n", c0, c1);
              printf("SSE = %g\n", sumsq);
              printf("SST = %g\n", SST);
              printf("R^2 = %g\n", 1 - sumsq / SST);
              eqnmsg = "Y = " + QString::number(c0) + " + " + QString::number(
                  c1) + " X, R^2 = " + QString::number(1 - sumsq / SST);
              emit setEqnMsg(eqnmsg);
            }
          else
            {
              eqnmsg = "Error.";
              emit setEqnMsg(eqnmsg);
            }
          double fitx[2] =
            { minVal, maxVal };
          double fity[2];
          fity[0] = c0 + c1 * fitx[0];
          fity[1] = c0 + c1 * fitx[1];
          emit drawFit(fitx, fity, 2);
        }
      else
        {
          eqnmsg = "No data for a linear fit.";
          emit setEqnMsg(eqnmsg);
        }
    }
  else
    {
      eqnmsg = "Not enough data for a linear fit.";
      emit setEqnMsg(eqnmsg);
    }
}

void
Clamp::countspikes()
{
  spikecount++;
  prevspktime = spktime;
  spktime = systime;
  if (spikecount >= 2)
    {
      ISI.push(spktime - prevspktime);
    }
}

void
Clamp::print()
{
#if 1
  QPrinter printer;
#else
  QPrinter printer(QPrinter::HighResolution);
#if QT_VERSION < 0x040000
  printer.setOutputToFile(true);
  printer.setOutputFileName("/tmp/FI.ps");
  printer.setColorMode(QPrinter::Color);
#else
  printer.setOutputFileName("/tmp/FI.pdf");
#endif
#endif

  QString docName = splot->title().text();
  if (!docName.isEmpty())
    {
      docName.replace(QRegExp(QString::fromLatin1("\n")), tr(" -- "));
      printer.setDocName(docName);
    }

  printer.setCreator("RTXI");
  printer.setOrientation(QPrinter::Landscape);

#if QT_VERSION >= 0x040000
  QPrintDialog dialog(&printer);
  if ( dialog.exec() )
    {
#else
  if (printer.setup())
    {
#endif
      RTXIPrintFilter filter;
      if (printer.colorMode() == QPrinter::GrayScale)
        {
          int options = QwtPlotPrintFilter::PrintAll;
          filter.setOptions(options);
          filter.color(QColor(29, 100, 141),
              QwtPlotPrintFilter::CanvasBackground);
          filter.color(Qt::white, QwtPlotPrintFilter::CurveSymbol);
        }
      splot->print(printer, filter);
    }
}

void
Clamp::exportSVG()
{
  QString fileName = "FI.svg";

#if QT_VERSION < 0x040000

#ifndef QT_NO_FILEDIALOG
  fileName = QFileDialog::getSaveFileName("FI.svg", "SVG Documents (*.svg)",
      this);
#endif
  if (!fileName.isEmpty())
    {
      // enable workaround for Qt3 misalignments
      QwtPainter::setSVGMode(true);
      QPicture picture;
      QPainter p(&picture);
      splot->print(&p, QRect(0, 0, 800, 600));
      p.end();
      picture.save(fileName, "svg");
    }

#elif QT_VERSION >= 0x040300

#ifdef QT_SVG_LIB
#ifndef QT_NO_FILEDIALOG
  fileName = QFileDialog::getSaveFileName(
      this, "Export File Name", QString(),
      "SVG Documents (*.svg)");
#endif
  if ( !fileName.isEmpty() )
    {
      QSvgGenerator generator;
      generator.setFileName(fileName);
      generator.setSize(QSize(800, 600));
      splot->print(generator);
    }
#endif
#endif
}

void
Clamp::togglerandom(bool on)
{
  randomize = on;
}

void
Clamp::toggleFIplot(bool on)
{
  plotFI = on;
}

bool
Clamp::OpenFile(QString FName)
{
  dataFile.setName(FName);
  if (dataFile.exists())
    {
      switch (QMessageBox::warning(this, "Current Clamp", tr(
          "This file already exists: %1.\n").arg(FName), "Overwrite", "Append",
          "Cancel", 0, 2))
        {
      case 0: // overwrite
        dataFile.remove();
        if (!dataFile.open(IO_Raw | IO_WriteOnly))
          {
            return false;
          }
        break;
      case 1: // append
        if (!dataFile.open(IO_Raw | IO_WriteOnly | IO_Append))
          {
            return false;
          }
        break;
      case 2: // cancel
        return false;
        break;
        }
    }
  else
    {
      if (!dataFile.open(IO_Raw | IO_WriteOnly))
        return false;
    }
  stream.setDevice(&dataFile);
  stream.setPrintableData(false); // write binary
  printf("File opened: %s\n", FName.latin1());
  return true;
}
