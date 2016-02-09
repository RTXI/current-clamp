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
 * This module allows you to deliver current clamp stimuli. With the Spike 
 * Detector module it can also plot an F-I curve.
 */

#include "current-clamp.h"
#include <algorithm>
#include <main_window.h>

#include <QSvgGenerator>
#include <QFileInfo>
#include <QtGlobal>

#if QT_VERSION >= 0x050000
	#include <QtPrintSupport/QPrintDialog>
	#include <QtPrintSupport/QPrinter>
#else
	#include <QPrintDialog>
	#include <QPrinter>
#endif

#include <time.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include <iostream>

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new Clamp();
}

static DefaultGUIModel::variable_t vars[] = {
	{ "Spike State", "Spike State", DefaultGUIModel::INPUT, },
	{ "Command", "Command", DefaultGUIModel::OUTPUT, },
	{ "Min Amplitude (pA)", "Minimum amplitude (pA) of command", 
	  DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Max Amplitude (pA)", "Maximum amplitude (pA) of command", 
	  DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Step Size (pA)", "Step size between minimum and maximum amplitude (pA) of command",
	  DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Delay (s)", "Delay (s) after last command", DefaultGUIModel::PARAMETER
	  | DefaultGUIModel::DOUBLE, },
	{ "Width (s)", "Width (s) of command", DefaultGUIModel::PARAMETER
	  | DefaultGUIModel::DOUBLE, },
	{ "Repeat", "Number of times to repeat cycle", DefaultGUIModel::PARAMETER
	  | DefaultGUIModel::DOUBLE, },
	{ "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, 
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

Clamp::Clamp(void) :  DefaultGUIModel("Current Clamp", ::vars, ::num_vars) {
	setWhatsThis(
		"<p><b>Current Clamp</b></p><p>This plugin requires the SpikeDetect plugin, which indicates when a spike has occurred using a threshold method. This plugin allows you to deliver current steps or ramps. Choose the <b>Clamp Mode</b>, then set the parameters for the input amplitudes and step size. This plugin uses an absolute delay between the end of a command input and the beginning of the next command. You can repeat the series of command inputs as many times as you like. Use the <b>Randomize</b> checkbox to randomize amplitudes within each cycle. To run a current clamp protocol, toggle the <b>Pause</b> button. You can edit parameter values directly in the textboxes, but you must click <b>Modify</b> to commit the changes. In <b>Step Mode</b>, you can choose to plot the FI Curve as it is generated and save a screenshot of the plot or save the data to a plain text file. The frequency is computed by averaging all the ISIs detected during a command step and taking the reciprocal. A linear fit is performed on all data points currently displayed in the scatterplot, not simply the data points acquired in the last run. Saving data points to a file will only save the data points acquired in the last run. You can completely overwrite or append the data to an existing file.</p>");
	initParameters();
	initStepArray();
	DefaultGUIModel::createGUI(vars, num_vars); // this is required to create the GUI
	customizeGUI();
	update(INIT);
	refresh(); // this is required to update the GUI with parameter and state values
	QTimer::singleShot(0, this, SLOT(resizeMe()));
	printf("Loaded Current Clamp:\n");
}

Clamp::~Clamp(void) {}

void Clamp::execute(void) {
	if (plotFI == true) spikestate = input(0);
	systime = count * dt; // current time, s

	if (cyclecount < repeat) { // as long as there are cycles remaining
		if (stepcount == 0 && randomize == true && random == false) { // shuffle amplitudes once at beginning of each cycle
			std::random_shuffle(arrstep, arrstep + nstep);
			random = true;
		}

		if (stepcount < nstep) { // as long as there are steps remaining in cycle
			if (systime >= delay && systime < (delay + width)) { // during command window, set output
				if (mode == STEP)	{
					output(0) = arrstep[stepcount];
				}
				else if (mode == RAMP) {
					if (systime >= delay && systime < (delay + width / 2)) {
						output(0) += arrrate[stepcount];
	  				}
					else {
						output(0) -= arrrate[stepcount];
					}
				}
				if ((spikestate == 1) && (plotFI == true)) countspikes();
			}

			else if (systime >= delay + width) { // after command window, shift delay, clean up
				delay = delay + delay0 + width;
				if (plotFI == true) {
					arrFIamp[static_cast<int> (cyclecount * nstep + stepcount)] = arrstep[stepcount] * 1e12; // in pA
					if (spikecount >= 2) { // at least 2 spikes (1 ISI) occurred
						arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)] = 1 / ISI.mean();
						if (ISI.mean() != 0 && 1 / ISI.mean() > yrangemax) { // resize window
							yrangemax = 1 / ISI.mean();
							emit setFIRange(minamp * 1e12, maxamp * 1e12,
							yrangemin, yrangemax);
						}
						if (ISI.mean() != 0) {
							printf("%f%f\n", arrstep[stepcount] * 1e12, 1 / ISI.mean());
							emit newDataPoint(arrFIamp[cyclecount * nstep + stepcount], arrFIHz[cyclecount * nstep + stepcount]);
						}
					}
					else if (spikecount == 1) { // only 1 spike occurred
						arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)] = GSL_NAN;
					}
					else { // no spikes occurred
						arrFIHz[static_cast<int> (cyclecount * nstep + stepcount)] = 0;
					}
					spikecount = 0;
					ISI.clear();
				}
				stepcount++;
				output(0) = 0;
			}
			else { // outside command window, set output to zero
			  output(0) = 0;
			}
		}
		if (stepcount == nstep) { // increment cycle count, reset flag for random shuffle
		 cyclecount++;
		 splot->nextSymbol();
		 stepcount = 0;
		 if (randomize == true) random = false;
		}
	}
	else { // all cycles are done
	} // end if cyclecount
	count++; // increment count to measure time
	return;
}

void Clamp::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
	case INIT:
		setParameter("Min Amplitude (pA)", QString::number(minamp * 1e12)); // initialized in A, display in pA
		setParameter("Max Amplitude (pA)", QString::number(maxamp * 1e12)); // initialized in A
		setParameter("Step Size (pA)", QString::number(stepsize * 1e12)); // initialized in A
		setParameter("Delay (s)", QString::number(delay)); //
		setParameter("Width (s)", QString::number(width)); //
		setParameter("Repeat", QString::number(repeat)); // initially 1
		setState("Time (s)", systime);
		emit setFIRange(minamp * 1e12, maxamp * 1e12, yrangemin, yrangemax);
		break;

	case MODIFY:
		delay = getParameter("Delay (s)").toDouble();
		delay0 = delay;
		width = getParameter("Width (s)").toDouble();
		minamp = getParameter("Min Amplitude (pA)").toDouble() * 1e-12; // set by user in pA
		maxamp = getParameter("Max Amplitude (pA)").toDouble() * 1e-12; // set by user in pA
		if (minamp == maxamp) {
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
		emit setEqnMsg(eqnmsg);
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
};

void Clamp::customizeGUI(void) {

	QGridLayout *customlayout = DefaultGUIModel::getLayout(); 
	customlayout->setColumnStretch(1,1);

	//overall GUI layout with a "horizontal box" copied from DefaultGUIModel
	QGroupBox *plotBox = new QGroupBox("FI Plot");
	QHBoxLayout *plotBoxLayout = new QHBoxLayout;
	plotBox->setLayout(plotBoxLayout);

	QPushButton *clearButton = new QPushButton("&Clear");
	QPushButton *linearfitButton = new QPushButton("Linear &Fit");
	QPushButton *savePlotButton = new QPushButton("Screenshot");
	QPushButton *printButton = new QPushButton("Print");
	QPushButton *saveDataButton = new QPushButton("Save Data");
	plotBoxLayout->addWidget(clearButton);
	plotBoxLayout->addWidget(linearfitButton);
	plotBoxLayout->addWidget(printButton);
	plotBoxLayout->addWidget(savePlotButton);
	plotBoxLayout->addWidget(saveDataButton);
	QLineEdit *eqnLine = new QLineEdit("Linear Equation");
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
	clearButton->setToolTip("Clear");
	savePlotButton->setToolTip("Save screenshot");
	saveDataButton->setToolTip("Save data");
	linearfitButton->setToolTip("Perform linear least-squares regression");
	printButton->setToolTip("Print plot");

	plotBox->hide();
	eqnLine->hide();
//	splot->setMinimumSize(450, 270);
	splot->hide();
	customlayout->addWidget(plotBox, 0, 1, 1, 1);
	customlayout->addWidget(eqnLine, 10, 1, 1, 1);
	customlayout->addWidget(splot, 1, 1, 3, 1);

	QGroupBox *modeBox = new QGroupBox("Clamp Mode");
	QHBoxLayout *modeBoxLayout = new QHBoxLayout;
	modeBox->setLayout(modeBoxLayout);
	QButtonGroup *modeButtons = new QButtonGroup;
	modeButtons->setExclusive(true);
	QRadioButton *stepButton = new QRadioButton("Step");
	modeBoxLayout->addWidget(stepButton);
	modeButtons->addButton(stepButton);
	stepButton->setEnabled(true);
	QRadioButton *rampButton = new QRadioButton("Ramp"); 
	modeBoxLayout->addWidget(rampButton);
	modeButtons->addButton(rampButton);
	stepButton->setChecked(true);
	QObject::connect(modeButtons,SIGNAL(buttonClicked(int)),this,SLOT(updateClampMode(int)));
	stepButton->setToolTip("Set mode to current steps");
	rampButton->setToolTip("Set mode to triangular current ramps");
	customlayout->addWidget(modeBox, 0, 0);

	QHBoxLayout *optionBoxLayout = new QHBoxLayout;
	QGroupBox *optionBoxGroup = new QGroupBox;
	QCheckBox *randomCheckBox = new QCheckBox("Randomize");
	optionBoxLayout->addWidget(randomCheckBox);
	QCheckBox *plotFICheckBox = new QCheckBox("Plot FI Curve");
	optionBoxLayout->addWidget(plotFICheckBox);
	QObject::connect(randomCheckBox,SIGNAL(toggled(bool)),this,SLOT(togglerandom(bool)));
	QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),eqnLine,SLOT(setVisible(bool)));
	QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),splot,SLOT(setVisible(bool)));
	QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),plotBox,SLOT(setVisible(bool)));
	QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleFIplot(bool)));
	QObject::connect(plotFICheckBox,SIGNAL(toggled(bool)),this,SLOT(resizeMe()));
	randomCheckBox->setToolTip("Randomize input amplitudes within a cycle");
	plotFICheckBox->setToolTip("Show/Hide FI plot area");
	optionBoxGroup->setLayout(optionBoxLayout);
	customlayout->addWidget(optionBoxGroup, 3, 0);

	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),savePlotButton,SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),printButton,SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),saveDataButton,SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),linearfitButton,SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),DefaultGUIModel::modifyButton,SLOT(setEnabled(bool)));
	DefaultGUIModel::pauseButton->setToolTip("Start/Stop current clamp protocol");
	DefaultGUIModel::modifyButton->setToolTip("Commit changes to parameter values");
	DefaultGUIModel::unloadButton->setToolTip("Close module");

	QObject::connect(this,SIGNAL(newDataPoint(double,double)),splot,SLOT(appendPoint(double,double)));
	QObject::connect(this,SIGNAL(setFIRange(double, double, double, double)),splot,SLOT(setAxes(double, double, double, double)));
	QObject::connect(this,SIGNAL(setPlotMode(bool)),plotFICheckBox,SLOT(setChecked(bool)));
	QObject::connect(this,SIGNAL(setStepMode(bool)),plotFICheckBox,SLOT(setEnabled(bool)));
	QObject::connect(this,SIGNAL(setStepMode(bool)),plotBox,SLOT(setEnabled(bool)));
	QObject::connect(this,SIGNAL(drawFit(double*, double*, int)),splot,SLOT(appendLine(double*, double*, int)));
	QObject::connect(this,SIGNAL(setEqnMsg(const QString &)), eqnLine,SLOT(setText(const QString &)));

	setLayout(customlayout);
}

void Clamp::initParameters() {
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

void Clamp::initStepArray() {
	arrstep = new double[nstep];
	arrrate = new double[nstep];
	arrFIamp = new double[static_cast<int> (nstep * repeat)];
	arrFIHz = new double[static_cast<int> (nstep * repeat)];
	for (int i = 0; i < nstep * repeat; i++) {
		arrFIamp[i] = 0;
		arrFIHz[i] = 0;
	}
	for (int i = 0; i < nstep; i++) {
		arrstep[i] = minamp + i * stepsize;
		arrrate[i] = arrstep[i] / width * 2 * dt;
	}
	random = false;
}

void Clamp::updateClampMode(int index) {
	if (index == 0) { // STEP
		mode = STEP;
		update(MODIFY);
		emit setStepMode(true);
		emit setPlotMode(false);
		printf("Entering STEP mode\n");
	}
	else if (index == 1)	{ // RAMP
		mode = RAMP;
		minamp = 0;
		setParameter("Min Amplitude (pA)", QString::number(minamp * 1e12)); // initialized in A, display in pA
		update(MODIFY);
		emit setStepMode(false);
		plotFI = false;
		printf("Entering RAMP mode\n");
	}
}

void Clamp::bookkeep() {
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

void Clamp::clearData() {
	yrangemax = 50;
	emit setFIRange(minamp * 1e12, maxamp * 1e12, yrangemin, yrangemax);
	eqnmsg = "Y = c0 + c1 * X";
	emit setEqnMsg(eqnmsg);
}

void Clamp::saveFIData() {
	QFileDialog* fd = new QFileDialog(this, "Save File As");//, TRUE);
	fd->setFileMode(QFileDialog::AnyFile);
	fd->setViewMode(QFileDialog::Detail);
	QString fileName;
	if (fd->exec() == QDialog::Accepted) {
		QStringList fileNames = fd->selectedFiles();
		if (!fileNames.isEmpty()) { fileName = fileNames.takeFirst(); };

		if (OpenFile(fileName)) {
//			stream.setPrintableData(true);
			for (int i = 0; i < cyclecount * nstep + stepcount; i++) {
				textStream << (double) arrFIamp[i] << (double) arrFIHz[i];
			}
			dataFile.close();
			printf("File closed.\n");
		}
		else {
			QMessageBox::information(this, "Current Clamp: Save FI Curve",
			                         "There was an error writing to this file. You can view\n"
			                         "the values that should be plotted in the terminal.\n");
		}
	}
}

void Clamp::fitData() {
	if (splot->dataExists()) {
		int n = splot->dataSize();
		if (n >= 2) {
			const double* tempx = splot->xData();
			const double* tempy = splot->yData();
			double x[n];
			double y[n];
			int j = 0;
			double maxVal = 0;
			double minVal = 0;
			for (int i = 0; i < n; i++) { // exclude any inf or NaN or zero rates from the FI linear fit
				if (gsl_finite(tempy[i]) == 1 || tempy[i] == 0) {
					x[j] = tempx[i];
					y[j] = tempy[i];
					if (x[j] > maxVal) maxVal = x[j];
					if (x[j] < minVal) minVal = x[j];
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
			for (int i = 0; i < j; i++) {
				double delta = y[i] - ymean;
				SST += delta * delta;
			}

			if (gsl_finite(c0) == 1 && gsl_finite(c1) == 1) {
				printf("Best fit: Y = %g + %g X\n", c0, c1);
				printf("SSE = %g\n", sumsq);
				printf("SST = %g\n", SST);
				printf("R^2 = %g\n", 1 - sumsq / SST);
				eqnmsg = "Y = " + QString::number(c0) + " + " + QString::number(c1) + " X, R^2 = " + QString::number(1 - sumsq / SST);
				emit setEqnMsg(eqnmsg);
			}
			else {
				eqnmsg = "Error.";
				emit setEqnMsg(eqnmsg);
			}
			double fitx[2] = { minVal, maxVal };
			double fity[2];
			fity[0] = c0 + c1 * fitx[0];
			fity[1] = c0 + c1 * fitx[1];
			emit drawFit(fitx, fity, 2);
		}
		else {
			eqnmsg = "No data for a linear fit.";
			emit setEqnMsg(eqnmsg);
		}
	}
	else {
		eqnmsg = "Not enough data for a linear fit.";
		emit setEqnMsg(eqnmsg);
	}
}

void Clamp::countspikes() {
	spikecount++;
	prevspktime = spktime;
	spktime = systime;
	if (spikecount >= 2) {
		ISI.push(spktime - prevspktime);
	}
}

void Clamp::print() {
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
	if (!docName.isEmpty()) {
		docName.replace(QRegExp(QString::fromLatin1("\n")), tr(" -- "));
		printer.setDocName(docName);
	}

	printer.setCreator("RTXI");
	printer.setOrientation(QPrinter::Landscape);

#if QT_VERSION >= 0x040000
	QPrintDialog dialog(&printer);
	if ( dialog.exec() ) {
#else
		if (printer.setup())	{
#endif 
/*
		RTXIPrintFilter filter;
		if (printer.colorMode() == QPrinter::GrayScale) {
		 int options = QwtPlotPrintFilter::PrintAll;
		 filter.setOptions(options);
		 filter.color(QColor(29, 100, 141),
			  QwtPlotPrintFilter::CanvasBackground);
		 filter.color(Qt::white, QwtPlotPrintFilter::CurveSymbol);
		}
*/
//	splot->print(printer, filter);
	QwtPlotRenderer *renderer = new QwtPlotRenderer;
	renderer->renderTo(splot, printer);
	}
}

void Clamp::exportSVG() {
	QString fileName = "FI.svg";

std::cout<<"flag 0"<<std::endl;

#if QT_VERSION < 0x040000
std::cout<<"flag 1"<<std::endl;

#ifndef QT_NO_FILEDIALOG
std::cout<<"flag 2"<<std::endl;
	fileName = QFileDialog::getSaveFileName("FI.svg", "SVG Documents (*.svg)",	this);
#endif
std::cout<<"flag 3"<<std::endl;
	if (!fileName.isEmpty()) {
		// enable workaround for Qt3 misalignments
		QwtPainter::setSVGMode(true);
		QPicture picture;
		QPainter p(&picture);
//		splot->print(&p, QRect(0, 0, 800, 600));
		QwtPlotRenderer *renderer = new QwtPlotRenderer;
		renderer->renderTo(splot, p, QRect(0, 0, 800, 600));
		p.end();
		picture.save(fileName, "svg");
	}

#elif QT_VERSION >= 0x040300
std::cout<<"flag 4"<<std::endl;
#ifdef QT_SVG_LIB
std::cout<<"flag 5"<<std::endl;
#ifndef QT_NO_FILEDIALOG
std::cout<<"flag 6"<<std::endl;
		fileName = QFileDialog::getSaveFileName(this, "Export File Name", 
		                                        QString(), "SVG Documents (*.svg)");
#endif
std::cout<<"flag 7"<<std::endl;
	if ( !fileName.isEmpty() ) {
		QSvgGenerator generator;
		generator.setFileName(fileName);
		generator.setSize(QSize(800, 600));
//		splot->print(generator);
	}
#endif
#endif
}

void Clamp::togglerandom(bool on) {
	randomize = on;
}

void Clamp::toggleFIplot(bool on) {
	plotFI = on;
}

bool Clamp::OpenFile(QString FName) {
	dataFile.setFileName(FName);
	if (dataFile.exists()) {
		switch (QMessageBox::warning(this, "Current Clamp", 
		                             tr("This file already exists: %1.\n").arg(FName), 
		                             "Overwrite", "Append","Cancel", 0, 2)) {
		case 0: // overwrite
			dataFile.remove();
			if (!dataFile.open( QIODevice::Unbuffered | QIODevice::WriteOnly)) return false;
			break;

		case 1: // append
			if (!dataFile.open( QIODevice::Unbuffered | QIODevice::WriteOnly | QIODevice::Append)) {
				return false;
			}
			break;

		case 2: // cancel
			return false;
			break;
		}
	}
	else {
		if (!dataFile.open( QIODevice::Unbuffered | QIODevice::WriteOnly )) return false;
	}
	binaryStream.setDevice(&dataFile);
//	stream.setPrintableData(false); // write binary
	printf("File opened: %s\n", FName.toStdString().data());
	return true;
}
