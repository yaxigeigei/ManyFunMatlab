/*
SatelliteRig.h - Library for the embedded part of Project Satellites.
Created by Duo Xu, September 18, 2016.
Released into the public domain.
*/

#ifndef SatelliteRig_h
#define SatelliteRig_h

#include "Arduino.h"
#include "WString.h"
#include "Printers.h"
#include "XBee.h"

#include "ManyRig.h"



const byte TRANSPARENT_MODE = 0;
const byte API_MODE = 1;



class SatelliteRig: public ManyRig
{
public:
	// Constructor
	SatelliteRig(byte id);

	// Rig info
	byte getRigId();
	void setCommunicationMode(byte modeId);

	// Handling incoming messages
	unsigned int getInputIndex();
	void serialRead();

	// Format and send data message
	unsigned long sendData(const char* tag, unsigned long t = millis());
	unsigned long sendData(const char* tag, unsigned long t, volatile byte num);
	unsigned long sendData(const char* tag, unsigned long t, volatile int num);
	unsigned long sendData(const char* tag, unsigned long t, volatile unsigned int num);
	unsigned long sendData(const char* tag, unsigned long t, volatile long num);
	unsigned long sendData(const char* tag, unsigned long t, volatile unsigned long num);
	unsigned long sendData(const char* tag, unsigned long t, volatile float num);
	unsigned long sendData(const char* tag, unsigned long t, volatile byte* dataArray, byte numData);
	unsigned long sendData(const char* tag, unsigned long t, volatile int* dataArray, byte numData);
	unsigned long sendData(const char* tag, unsigned long t, volatile unsigned int* dataArray, byte numData);
	unsigned long sendData(const char* tag, unsigned long t, volatile long* dataArray, byte numData);
	unsigned long sendData(const char* tag, unsigned long t, volatile unsigned long* dataArray, byte numData);
	unsigned long sendData(const char* tag, unsigned long t, volatile float* dataArray, byte numData);

	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t = millis());
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile byte num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile int num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned int num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile long num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned long num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile float num);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile byte* dataArray, byte numData);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile int* dataArray, byte numData);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned int* dataArray, byte numData);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile long* dataArray, byte numData);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile unsigned long* dataArray, byte numData);
	unsigned long sendData(const __FlashStringHelper* tag, unsigned long t, volatile float* dataArray, byte numData);


private:
	// Rig variables
	byte _rigId = 255;
	byte _communicationMode = TRANSPARENT_MODE;
	bool _isWaiting = true;			// a flag used for parsing incoming address (i.e. rig ID) but not wasting time on irrelevant stuff

	// Create the XBee object
	XBee xbee = XBee();
	XBeeResponse response = XBeeResponse();
	XBeeAddress64 addr64 = XBeeAddress64(0x00000000, 0x00000000);	// SH + SL address of the coordinator

	// Create reusable response objects for responses we expect to handle
	ZBRxResponse rx = ZBRxResponse();
	ZBTxStatusResponse txStatus = ZBTxStatusResponse();
	ModemStatusResponse msr = ModemStatusResponse();



	// General character parser
	void parseChar(char ch);

	// Read a XBee frame
	String readXBeeFrame();

	// Sending formatted data message
	void sendString(String msg);
	void sendXBeeFrame(String msg);
};

#endif