/*******************************************************************************
* Copyright (c) 2012-2013, The Microsystems Design Labratory (MDL)
* Department of Computer Science and Engineering, The Pennsylvania State University
* Exascale Computing Lab, Hewlett-Packard Company
* All rights reserved.
* 
* This source code is part of NVSim - An area, timing and power model for both 
* volatile (e.g., SRAM, DRAM) and non-volatile memory (e.g., PCRAM, STT-RAM, ReRAM, 
* SLC NAND Flash). The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Author list: 
*   Cong Xu	    ( Email: czx102 at psu dot edu 
*                     Website: http://www.cse.psu.edu/~czx102/ )
*   Xiangyu Dong    ( Email: xydong at cse dot psu dot edu
*                     Website: http://www.cse.psu.edu/~xydong/ )
*******************************************************************************/


#include "MemCell.h"
#include "formula.h"
#include "global.h"
#include "macros.h"
#include <math.h>

MemCell::MemCell() {
	// TODO Auto-generated constructor stub
	area                = 0;
	aspectRatio         = 0;
	resistanceOn        = 0;
	resistanceOff       = 0;
	readMode            = true;
	readVoltage         = 0;
	readCurrent         = 0;
	readPower           = 0;
    wordlineBoostRatio  = 1.0;
	resetMode           = true;
	resetVoltage        = 0;
	resetCurrent        = 0;
	minSenseVoltage     = 0.08;
	resetPulse          = 0;
	resetEnergy         = 0;
	setMode             = true;
	setVoltage          = 0;
	setCurrent          = 0;
	setPulse            = 0;
	accessType          = CMOS_access;
	processNode         = 0;
	setEnergy           = 0;

	/* Optional */
	stitching         = 0;
	gateOxThicknessFactor = 2;
	widthSOIDevice = 0;
	widthAccessCMOS   = 0;
	voltageDropAccessDevice = 0;
	leakageCurrentAccessDevice = 0;
}

MemCell::~MemCell() {
	// TODO Auto-generated destructor stub
}

void MemCell::ReadCellFromFile(const string & inputFile)
{
	FILE *fp = fopen(inputFile.c_str(), "r");
	char line[5000];
	char tmp[5000];

	if (!fp) {
		cout << inputFile << " cannot be found!\n";
		exit(-1);
	}

	while (fscanf(fp, "%[^\n]\n", line) != EOF) {
		if (!strncmp("-ProcessNode", line, strlen("-ProcessNode"))) {
			sscanf(line, "-ProcessNode: %d", &processNode);
			continue;
		}
		if (!strncmp("-CellArea", line, strlen("-CellArea"))) {
			sscanf(line, "-CellArea (F^2): %lf", &area);
			continue;
		}
		if (!strncmp("-CellAspectRatio", line, strlen("-CellAspectRatio"))) {
			sscanf(line, "-CellAspectRatio: %lf", &aspectRatio);
			heightInFeatureSize = sqrt(area * aspectRatio);
			widthInFeatureSize = sqrt(area / aspectRatio);
			continue;
		}

		if (!strncmp("-ResistanceOn", line, strlen("-ResistanceOn"))) {
			sscanf(line, "-ResistanceOn (ohm): %lf", &resistanceOn);
			continue;
		}
		if (!strncmp("-ResistanceOff", line, strlen("-ResistanceOff"))) {
			sscanf(line, "-ResistanceOff (ohm): %lf", &resistanceOff);
			continue;
		}
		if (!strncmp("-CapacitanceOn", line, strlen("-CapacitanceOn"))) {
			sscanf(line, "-CapacitanceOn (F): %lf", &capacitanceOn);
			continue;
		}
		if (!strncmp("-CapacitanceOff", line, strlen("-CapacitanceOff"))) {
			sscanf(line, "-CapacitanceOff (F): %lf", &capacitanceOff);
			continue;
		}

		if (!strncmp("-GateOxThicknessFactor", line, strlen("-GateOxThicknessFactor"))) {
			sscanf(line, "-GateOxThicknessFactor: %lf", &gateOxThicknessFactor);
			continue;
		}

		if (!strncmp("-SOIDeviceWidth (F)", line, strlen("-SOIDeviceWidth (F)"))) {
			sscanf(line, "-SOIDeviceWidth (F): %lf", &widthSOIDevice);
			continue;
		}

		if (!strncmp("-ReadMode", line, strlen("-ReadMode"))) {
			sscanf(line, "-ReadMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				readMode = true;
			else
				readMode = false;
			continue;
		}
		if (!strncmp("-ReadVoltage", line, strlen("-ReadVoltage"))) {
			sscanf(line, "-ReadVoltage (V): %lf", &readVoltage);
			continue;
		}
		if (!strncmp("-ReadCurrent", line, strlen("-ReadCurrent"))) {
			sscanf(line, "-ReadCurrent (uA): %lf", &readCurrent);
			readCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-ReadPower", line, strlen("-ReadPower"))) {
			sscanf(line, "-ReadPower (uW): %lf", &readPower);
			readPower /= 1e6;
			continue;
		}
		if (!strncmp("-WordlineBoostRatio", line, strlen("-WordlineBoostRatio"))) {
			sscanf(line, "-WordlineBoostRatio: %lf", &wordlineBoostRatio);
			continue;
		}
		if (!strncmp("-MinSenseVoltage", line, strlen("-MinSenseVoltage"))) {
			sscanf(line, "-MinSenseVoltage (mV): %lf", &minSenseVoltage);
			minSenseVoltage /= 1e3;
			continue;
		}


		if (!strncmp("-ResetMode", line, strlen("-ResetMode"))) {
			sscanf(line, "-ResetMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				resetMode = true;
			else
				resetMode = false;
			continue;
		}
		if (!strncmp("-ResetVoltage", line, strlen("-ResetVoltage"))) {
			sscanf(line, "-ResetVoltage (V): %lf", &resetVoltage);
			continue;
		}
		if (!strncmp("-ResetCurrent", line, strlen("-ResetCurrent"))) {
			sscanf(line, "-ResetCurrent (uA): %lf", &resetCurrent);
			resetCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-ResetVoltage", line, strlen("-ResetVoltage"))) {
			sscanf(line, "-ResetVoltage (V): %lf", &resetVoltage);
			continue;
		}
		if (!strncmp("-ResetPulse", line, strlen("-ResetPulse"))) {
			sscanf(line, "-ResetPulse (ns): %lf", &resetPulse);
			resetPulse /= 1e9;
			continue;
		}
		if (!strncmp("-ResetEnergy", line, strlen("-ResetEnergy"))) {
			sscanf(line, "-ResetEnergy (pJ): %lf", &resetEnergy);
			resetEnergy /= 1e12;
			continue;
		}

		if (!strncmp("-SetMode", line, strlen("-SetMode"))) {
			sscanf(line, "-SetMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				setMode = true;
			else
				setMode = false;
			continue;
		}
		if (!strncmp("-SetVoltage", line, strlen("-SetVoltage"))) {
			sscanf(line, "-SetVoltage (V): %lf", &setVoltage);
			continue;
		}
		if (!strncmp("-SetCurrent", line, strlen("-SetCurrent"))) {
			sscanf(line, "-SetCurrent (uA): %lf", &setCurrent);
			setCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-SetVoltage", line, strlen("-SetVoltage"))) {
			sscanf(line, "-SetVoltage (V): %lf", &setVoltage);
			continue;
		}
		if (!strncmp("-SetPulse", line, strlen("-SetPulse"))) {
			sscanf(line, "-SetPulse (ns): %lf", &setPulse);
			setPulse /= 1e9;
			continue;
		}
		if (!strncmp("-SetEnergy", line, strlen("-SetEnergy"))) {
			sscanf(line, "-SetEnergy (pJ): %lf", &setEnergy);
			setEnergy /= 1e12;
			continue;
		}

		if (!strncmp("-AccessType", line, strlen("-AccessType"))) {
			sscanf(line, "-AccessType: %s", tmp);
			if (!strcmp(tmp, "CMOS"))
				accessType = CMOS_access;
			else if (!strcmp(tmp, "BJT"))
				accessType = BJT_access;
			else if (!strcmp(tmp, "diode"))
				accessType = diode_access;
			else
				accessType = none_access;
			continue;
		}

		if (!strncmp("-AccessCMOSWidth", line, strlen("-AccessCMOSWidth"))) {
			if (accessType != CMOS_access)
				cout << "Warning: The input of CMOS access transistor width is ignored because the cell is not CMOS-accessed." << endl;
			else
				sscanf(line, "-AccessCMOSWidth (F): %lf", &widthAccessCMOS);
			continue;
		}

		if (!strncmp("-VoltageDropAccessDevice", line, strlen("-VoltageDropAccessDevice"))) {
			sscanf(line, "-VoltageDropAccessDevice (V): %lf", &voltageDropAccessDevice);
			continue;
		}

		if (!strncmp("-LeakageCurrentAccessDevice", line, strlen("-LeakageCurrentAccessDevice"))) {
			sscanf(line, "-LeakageCurrentAccessDevice (uA): %lf", &leakageCurrentAccessDevice);
			leakageCurrentAccessDevice /= 1e6;
			continue;
		}
	}

	fclose(fp);
}


void MemCell::CellScaling(int _targetProcessNode) {
	if ((processNode > 0) && (processNode != _targetProcessNode)) {		
	double scalingFactor = (double)processNode / _targetProcessNode;
		resistanceOn *= scalingFactor * scalingFactor;
		resistanceOff *= scalingFactor * scalingFactor;
		if (!setMode) {
			setCurrent /= scalingFactor;
		} else {
			setVoltage *= scalingFactor;
		}
		if (!resetMode) {
			resetCurrent /= scalingFactor;
		} else {
			resetVoltage *= scalingFactor;
		}
		if (accessType == diode_access) {
			capacitanceOn /= scalingFactor; //TO-DO
			capacitanceOff /= scalingFactor; //TO-DO
		}
	processNode = _targetProcessNode;
	}
}

void MemCell::CalculateWriteEnergy() {
	if (resetEnergy == 0) {
		if (resetMode) {
			resetEnergy = fabs(resetVoltage) * (fabs(resetVoltage) - voltageDropAccessDevice) / resistanceOn * resetPulse;
		} else {
		if (resetVoltage == 0) {
			resetEnergy = tech->vdd * fabs(resetCurrent) * resetPulse; /* TO-DO consider charge pump*/
		} else {
			resetEnergy = fabs(resetVoltage) * fabs(resetCurrent) * resetPulse;
		}
		/* previous model seems to be problematic
		resetEnergy = resetCurrent * (resetCurrent * resistanceOff + voltageDropAccessDevice) * resetPulse;
		*/
		}
	}
	if (setEnergy == 0) {
		if (setMode) {
			setEnergy = fabs(setVoltage) * (fabs(setVoltage) - voltageDropAccessDevice) / resistanceOn * setPulse;
		} else {
		if (resetVoltage == 0){
			setEnergy = tech->vdd * fabs(setCurrent) * setPulse; /*TO-DO consider charge pump*/
		} else {
			setEnergy = fabs(setVoltage) * fabs(setCurrent) * setPulse;
		}
		/* previous model seems to be problematic
		setEnergy = setCurrent * (setCurrent * resistanceOff + voltageDropAccessDevice) * setPulse;
		*/
		}
	}
}

double MemCell::CalculateReadPower() { /* TO-DO consider charge pumped read voltage */
	if (readPower == 0) {
		if (cell->readMode) {	/* voltage-sensing */
			if (readVoltage == 0) { /* Current-in voltage sensing */
				return tech->vdd * readCurrent;
			}
			if (readCurrent == 0) { /* Voltage-divider sensing */
				double resInSerialForSenseAmp, maxBitlineCurrent;
				resInSerialForSenseAmp = sqrt(resistanceOn * resistanceOff);
				maxBitlineCurrent = (readVoltage - voltageDropAccessDevice) / (resistanceOn + resInSerialForSenseAmp);
				return tech->vdd * maxBitlineCurrent;
			}
		} else { /* current-sensing */
			double maxBitlineCurrent = (readVoltage - voltageDropAccessDevice) / resistanceOn;
			return tech->vdd * maxBitlineCurrent;
		}
	} else {
		return -1.0; /* should not call the function if read energy exists */
	}
	return -1.0;
}

void MemCell::PrintCell()
{
	cout << "Memory Cell: MRAM (Magnetoresistive)" << endl;
		
	cout << "Cell Area (F^2)    : " << area << " (" << heightInFeatureSize << "Fx" << widthInFeatureSize << "F)" << endl;
	cout << "Cell Aspect Ratio  : " << aspectRatio << endl;

	if (resistanceOn < 1e3 )
		cout << "Cell Turned-On Resistance : " << resistanceOn << "ohm" << endl;
	else if (resistanceOn < 1e6)
		cout << "Cell Turned-On Resistance : " << resistanceOn / 1e3 << "Kohm" << endl;
	else
		cout << "Cell Turned-On Resistance : " << resistanceOn / 1e6 << "Mohm" << endl;
	if (resistanceOff < 1e3 )
		cout << "Cell Turned-Off Resistance: "<< resistanceOff << "ohm" << endl;
	else if (resistanceOff < 1e6)
		cout << "Cell Turned-Off Resistance: "<< resistanceOff / 1e3 << "Kohm" << endl;
	else
		cout << "Cell Turned-Off Resistance: "<< resistanceOff / 1e6 << "Mohm" << endl;

	if (readMode) {
		cout << "Read Mode: Voltage-Sensing" << endl;
		if (readCurrent > 0)
			cout << "  - Read Current: " << readCurrent * 1e6 << "uA" << endl;
		if (readVoltage > 0)
			cout << "  - Read Voltage: " << readVoltage << "V" << endl;
	} else {
		cout << "Read Mode: Current-Sensing" << endl;
		if (readCurrent > 0)
			cout << "  - Read Current: " << readCurrent * 1e6 << "uA" << endl;
		if (readVoltage > 0)
			cout << "  - Read Voltage: " << readVoltage << "V" << endl;
	}

	if (resetMode) {
		cout << "Reset Mode: Voltage" << endl;
		cout << "  - Reset Voltage: " << resetVoltage << "V" << endl;
	} else {
		cout << "Reset Mode: Current" << endl;
		cout << "  - Reset Current: " << resetCurrent * 1e6 << "uA" << endl;
	}
	cout << "  - Reset Pulse: " << TO_SECOND(resetPulse) << endl;

	if (setMode) {
		cout << "Set Mode: Voltage" << endl;
		cout << "  - Set Voltage: " << setVoltage << "V" << endl;
	} else {
		cout << "Set Mode: Current" << endl;
		cout << "  - Set Current: " << setCurrent * 1e6 << "uA" << endl;
	}
	cout << "  - Set Pulse: " << TO_SECOND(setPulse) << endl;

	switch (accessType) {
	case CMOS_access:
		cout << "Access Type: CMOS" << endl;
		break;
	case BJT_access:
		cout << "Access Type: BJT" << endl;
		break;
	case diode_access:
		cout << "Access Type: Diode" << endl;
		break;
	default:
		cout << "Access Type: None Access Device" << endl;
	}
}
