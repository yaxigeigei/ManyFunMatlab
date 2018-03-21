/*
ManyRig.h - Library for rig specific utilities.
Created by Duo Xu, November 23, 2017.
*/

#ifndef ManyRig_h
#define ManyRig_h

#include "Arduino.h"

class ManyRig
{
public:
	ManyRig();

	// Pin-out (values defined in constructor)
	byte lickDetectorPin;
	byte pawDetectorPin;

	byte waterValvePin;
	byte speakerPin;
	byte airPuffPin;
	byte odorEmptyPin;
	byte odorPins[2];
	byte syncPin;
	byte camPin;
	byte randPin;

	// IO
	bool isPawOn();
	bool isLickOn();
	void playTone(unsigned int freq, unsigned long durInMs);
	void playSweep(unsigned int freqStart, unsigned int freqEnd, unsigned long durInMs);
	void playNoise(unsigned long durInMs);
	void deliverWater(unsigned long durInMs);
	void deliverOdor(byte odorId, unsigned long durInMs);
	void deliverAirPuff(unsigned long durInMs);
	void sendNumTTL(unsigned long num);
	void sendNumTTL(byte pin, unsigned long num);

	// Computing
	byte choose(byte* probVector, byte numChoices);

private:

};

#endif