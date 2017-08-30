/*
ManyRig.h - Library facilitating the control of behaviral experiments.
Created by Duo Xu, September 24, 2016.
Released into the public domain.
*/

#include "Arduino.h"
#include "ManyRig.h"



ManyRig::ManyRig() {
	// do nothing
}

void ManyRig::setDelimiter(char d) {
	_delimiter = d;
}

char ManyRig::getDelimiter() {
	return _delimiter;
}



void ManyRig::attachParser(void(*f)(void)) {
	_parserFunc = f;
}

void ManyRig::detachParser() {
	_parserFunc = NULL;
}

unsigned int ManyRig::getInputIndex() {
	return _numDelimiter;
}

unsigned long ManyRig::getInputValue() {
	return _inputVal;
}

String ManyRig::getCmdString() {
	return _cmdString;
}

void ManyRig::serialRead() {
	// Handle command identification and dispatching

	if (Serial.available())
	{
		// Read byte as char
		char ch = Serial.read();

		if (isDigit(ch))
		{
			// Accumulate digits to assemble the incoming value
			_inputVal = _inputVal * 10 + ch - '0';
		}
		else if (isAlpha(ch))
		{
			// Accumulate letters to assemble the incoming string
			_cmdString += ch;
		}
		else if (ch == _delimiter || isControl(ch))
		{
			// Parse command and incoming value
			if (_parserFunc != NULL && _cmdString.length() > 0)
				_parserFunc();

			// Clear incoming value of the current input
			_inputVal = 0;
		}

		// Keep track of the number of delimiters for indexing inputs
		if (ch == _delimiter)
			_numDelimiter++;

		// Reset reader state for identification
		if (isControl(ch))
		{
			_inputVal = 0;
			_cmdString = "";
			_numDelimiter = 0;
		}
	}
}



void ManyRig::delay(unsigned long dur) {
	// Delay and read serial command

	unsigned long t0 = millis();
	while (millis() - t0 < dur)
		serialRead();
}

bool ManyRig::delayUntil(bool(*f)(void)) {
	// Delay and read serial command until function returns true

	bool b = false;

	while (!b) {
		serialRead();
		b = f();
	}

	return b;
}

bool ManyRig::delayUntil(bool(*f)(void), unsigned long timeout) {
	// Delay and read serial command until function returns true or timeout

	bool b = false;

	unsigned long t0 = millis();
	while (millis() - t0 < timeout && !b) {
		serialRead();
		b = f();
	}

	return b;
}

bool ManyRig::delayContinue(bool(*f)(void), unsigned long unitTime) {
	// Delay and read serial command reiterate when function returns true

	bool b = false;
	unsigned long t0 = millis();
	unsigned long t = t0;

	while (millis() - t < unitTime) {
		serialRead();
		if (f()) {
			t = millis();
			b = true;
		}
	}

	return b;
}



void ManyRig::attachSender(void(*f)(String)) {
	if (f != NULL)
		_senderFunc = f;
}

void ManyRig::serialSend(String msg) {
	Serial.println(msg);
}

unsigned long ManyRig::sendData(char* tag, unsigned long t) {
	// Send data message with the header

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile byte num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile int num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile unsigned int num) {
	// Send data message by event type, time, and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile unsigned long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile float num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile byte* dataArray, byte numData) {
	// Send data message with the header and an array of bytes, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile int* dataArray, byte numData) {
	// Send data message with the header and an array of ints, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile unsigned int* dataArray, byte numData) {
	// Send data message with the header and an array of unsigned ints, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile long* dataArray, byte numData) {
	// Send data message with the header and an array of long integers, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile unsigned long* dataArray, byte numData) {
	// Send data message with the header and an array of unsigned long integers, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}

unsigned long ManyRig::sendData(char* tag, unsigned long t, volatile float* dataArray, byte numData) {
	// Send data message with the header and an array of floating point numbers, separated by commas

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	for (int i = 0; i < numData; i++) {
		msg += _delimiter;
		msg += dataArray[i];
	}

	_senderFunc(msg);

	return micros() - tStart;
}



unsigned long ManyRig::sendNumTTL(byte pin, unsigned long num) {
	// Send integer via TTL. The duration equals the number times 100us.

	unsigned long durInUs = num * 100;

	digitalWrite(pin, HIGH);

	delayMicroseconds(durInUs);

	digitalWrite(pin, LOW);

	return durInUs;
}