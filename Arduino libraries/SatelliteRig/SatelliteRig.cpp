/*
SatelliteRig.cpp - Library for the embedded part of Project Satellites.
Created by Duo Xu, September 18, 2016.
Released into the public domain.
*/

#include "SatelliteRig.h"

#include "Arduino.h"
#include "WString.h"
#include "Printers.h"
#include "XBee.h"

#include "ManyRig.h"



SatelliteRig::SatelliteRig(byte id)
{
	_rigId = id;
}

byte SatelliteRig::getRigId()
{
	return _rigId;
}

void SatelliteRig::setCommunicationMode(byte modeId)
{
	modeId = min(modeId, API_MODE);
	_communicationMode = modeId;
}



unsigned int SatelliteRig::getInputIndex() {
	return _numDelimiter - 1;
}

void SatelliteRig::serialRead() {
	// Handle command identification and dispatching

	if (_communicationMode == TRANSPARENT_MODE && Serial.available())
	{
		// Read byte as char
		char ch = Serial.read();

		// Basic parseing of the incoming char
		parseChar(ch);
	}
	else if (_communicationMode == API_MODE)
	{
		// Read the message in a XBee frame
		String msg = readXBeeFrame();

		// Parse message one char at a time
		for (byte i = 0; i < msg.length(); i++)
			parseChar(msg.charAt(i));
	}
}

void SatelliteRig::parseChar(char ch) {
	// Handle command identification and dispatching

	if (_isWaiting && _numDelimiter == 0)
	{
		if (isDigit(ch))
		{
			// Accumulate digits to assemble the incoming value
			_inputVal = _inputVal * 10 + ch - '0';
		}
		else if (ch == _delimiter && _inputVal == getRigId())
		{
			// Rig ID matches that in the header of the incoming command
			_isWaiting = false;
			_inputVal = 0;
			_cmdString = "";
		}
	}
	else if (!_isWaiting)
	{
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
	}


	// Keep track of the number of delimiters for indexing inputs
	if (ch == _delimiter)
		_numDelimiter++;

	// Reset reader state for identification
	if (isControl(ch))
	{
		_inputVal = 0;
		_cmdString = "";

		_isWaiting = true;
		_numDelimiter = 0;
	}
}

String SatelliteRig::readXBeeFrame() {
	// Read an API frame

	xbee.readPacket();

	String msgString = String();

	// Got something
	if (xbee.getResponse().isAvailable())
	{
		// Got a zb rx packet
		byte apiId = xbee.getResponse().getApiId();

		if (apiId == ZB_RX_RESPONSE)
		{
			// now fill our zb rx class
			xbee.getResponse().getZBRxResponse(rx);

			// Put message into the an array object
			byte msgLength = rx.getDataLength();
			char msgArray[msgLength + 2];
			for (byte i = 0; i < msgLength; i++)
				msgArray[i] = rx.getData(i);

			// Append return chars for later parsing
			msgArray[msgLength] = '\r';
			msgArray[msgLength + 1] = '\n';

			msgString = msgArray;
		}
	}
	else if (xbee.getResponse().isError())
	{
		//nss.print("Error reading packet.  Error code: ");
		//nss.println(xbee.getResponse().getErrorCode());
	}

	return msgString;
}





void SatelliteRig::sendString(String msg) {
	String fullMsg = String();

	// Add rig ID (with a delimiter) in the front of the message
	if (_rigId != 255) {
		fullMsg += _rigId;
		fullMsg += _delimiter;
	}
	fullMsg += msg;
	fullMsg += "\r\n";

	// Send out message via the selected physical path
	if (_communicationMode == TRANSPARENT_MODE)
	{
		Serial.print(fullMsg);
	}
	else if (_communicationMode == API_MODE)
	{
		sendXBeeFrame(fullMsg);
	}
}

void SatelliteRig::sendXBeeFrame(String msg) {
	// Send an API frame from a formatted string

	unsigned int msgLength = msg.length();
	byte byteArray[msgLength+1];
	msg.getBytes(byteArray, msgLength+1);

	ZBTxRequest zbTx = ZBTxRequest(addr64, byteArray, msgLength);
	xbee.send(zbTx);

}





unsigned long SatelliteRig::sendData(const char* tag, unsigned long t) {
	// Send data message with the header

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile byte num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile int num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile unsigned int num) {
	// Send data message by event type, time, and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile unsigned long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile float num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile byte* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile int* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile unsigned int* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile long* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile unsigned long* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const char* tag, unsigned long t, volatile float* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}



unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t) {
	// Send data message with the header

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile byte num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile int num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned int num) {
	// Send data message by event type, time, and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned long num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile float num) {
	// Send data message with the header and value

	unsigned long tStart = micros();

	String msg = String();

	msg += tag;
	msg += _delimiter;
	msg += t;
	msg += _delimiter;
	msg += num;

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile byte* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile int* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned int* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile long* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned long* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}

unsigned long SatelliteRig::sendData(const __FlashStringHelper* tag, unsigned long t, volatile float* dataArray, byte numData) {
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

	sendString(msg);

	return micros() - tStart;
}