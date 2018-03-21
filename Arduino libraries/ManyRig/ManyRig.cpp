/*
ManyRig.cpp - Library for rig specific utilities.
Created by Duo Xu, November 23, 2017.
*/

#include "Arduino.h"
#include "ManyRig.h"

ManyRig::ManyRig()
{
	// Define pins in array
	lickDetectorPin = 14;
	pawDetectorPin = 15;

	waterValvePin = 4;
	speakerPin = 6;
	airPuffPin = 16;
	odorEmptyPin = 17;
	odorPins[0] = 18;
	odorPins[1] = 19;
	syncPin = 22;
	camPin = 23;
	randPin = 28;

	// Setup pin mode
	pinMode(lickDetectorPin, INPUT);
	pinMode(pawDetectorPin, INPUT);

	pinMode(waterValvePin, OUTPUT);
	pinMode(speakerPin, OUTPUT);
	pinMode(airPuffPin, OUTPUT);
	pinMode(odorEmptyPin, OUTPUT);
	pinMode(odorPins[0], OUTPUT);
	pinMode(odorPins[1], OUTPUT);
	pinMode(syncPin, OUTPUT);
	pinMode(camPin, OUTPUT);

	// Set random seed from an unused analog pin
	randomSeed(randPin);
}

bool ManyRig::isPawOn()
{
	return digitalRead(pawDetectorPin) == HIGH;
}

bool ManyRig::isLickOn()
{
	return digitalRead(lickDetectorPin) == HIGH;
}

void ManyRig::playTone(unsigned int freq, unsigned long durInMs)
{
	if (durInMs > 0) {
		tone(speakerPin, freq, durInMs);
		delay(durInMs);
	}
}

void ManyRig::playSweep(unsigned int freqStart, unsigned int freqEnd, unsigned long durInMs)
{
	unsigned int dT = 10;
	unsigned int numSteps = durInMs / dT;
	unsigned int dFreq = (freqEnd - freqStart) / numSteps;

	unsigned long t0 = millis();
	unsigned int currentStep = 0;

	while (currentStep < numSteps)
	{
		if (millis() - t0 >= dT)
		{
			currentStep++;
			t0 = millis();
			tone(speakerPin, freqStart + currentStep * dFreq);
		}
	}
	noTone(speakerPin);
}

void ManyRig::playNoise(unsigned long durInMs)
{
	unsigned long dur = durInMs * 1e3;
	unsigned long t = micros();

	while (micros() - t < dur)
		digitalWrite(speakerPin, random(2));

	digitalWrite(speakerPin, LOW);
}

void ManyRig::deliverWater(unsigned long durInMs)
{
	digitalWrite(waterValvePin, HIGH);
	delay(durInMs);
	digitalWrite(waterValvePin, LOW);
}

void ManyRig::deliverOdor(byte odorId, unsigned long durInMs)
{
	digitalWrite(odorPins[odorId], HIGH);
	digitalWrite(odorEmptyPin, HIGH);
	delay(durInMs);
	digitalWrite(odorEmptyPin, LOW);
	digitalWrite(odorPins[odorId], LOW);
}

void ManyRig::deliverAirPuff(unsigned long durInMs)
{
	digitalWrite(airPuffPin, HIGH);
	delay(durInMs);
	digitalWrite(airPuffPin, LOW);
}

void ManyRig::sendNumTTL(unsigned long num)
{
	// Send integer via TTL. The duration equals the number times 1000us.
	sendNumTTL(syncPin, num);
}

void ManyRig::sendNumTTL(byte pin, unsigned long num)
{
	// Send integer via TTL. The duration equals the number times 1000us.

	digitalWrite(pin, HIGH);
	delayMicroseconds(num * 1000);
	digitalWrite(pin, LOW);
}

byte ManyRig::choose(byte* probVector, byte numChoices)
{
	byte randNum = random(0, 100);
	int sumProb = 0;

	for (byte i = 0; i < numChoices; i++)
	{
		sumProb = sumProb + probVector[i];

		if (i == numChoices - 1)
			sumProb = 100;

		if (randNum < sumProb)
			return i;
	}

	return 0;
}


